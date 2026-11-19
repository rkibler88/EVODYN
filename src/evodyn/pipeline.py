"""The run: stages, in order, with their failures recorded rather than hidden.

A stage that cannot run is skipped and its reason is written into the manifest
and into the limits section of the report. A stage that produces a *wrong*
answer - a numbering mismatch, a degenerate score column - is fatal. The
distinction is the whole design: missing data narrows what a run concludes,
whereas a broken mapping invalidates it.
"""
from __future__ import annotations

import pathlib
import time
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from evodyn import __version__
from evodyn.config import EvodynConfig
from evodyn.exceptions import EvodynError, MappingError, StageSkipped
from evodyn.numbering import CanonicalMap
from evodyn.provenance import Manifest


@dataclass
class RunState:
    """Everything the stages hand to each other."""

    cfg: EvodynConfig
    manifest: Manifest
    cmap: CanonicalMap | None = None
    alignment: pathlib.Path | None = None
    couplings: pd.DataFrame | None = None
    network: Any = None
    distances: pd.DataFrame | None = None
    dccm: pd.DataFrame | None = None
    dccm_replicates: dict[int, pd.DataFrame] = field(default_factory=dict)
    rmsf: pd.DataFrame | None = None
    secondary: pd.DataFrame | None = None
    conservation: pd.DataFrame | None = None
    residues: pd.DataFrame | None = None
    pairs: pd.DataFrame | None = None
    enrichment: pd.DataFrame | None = None
    summary: dict[str, Any] = field(default_factory=dict)


class Pipeline:
    """Runs the configured stages against one target."""

    def __init__(self, cfg: EvodynConfig):
        cfg.validate()
        self.cfg = cfg
        self.out = cfg.out
        self.out.mkdir(parents=True, exist_ok=True)
        self.manifest = Manifest(config=cfg.to_dict(), outdir=self.out)
        self.state = RunState(cfg=cfg, manifest=self.manifest)
        self.state.summary["evodyn_version"] = __version__

    # ------------------------------------------------------------------ run
    def run(self) -> RunState:
        order = ["sequences", "coevolution", "network", "dynamics",
                 "structure", "controls", "export"]
        for name in [s for s in order if s in self.cfg.stages]:
            rec = self.manifest.stage(name)
            t0 = time.time()
            try:
                msg = getattr(self, f"_stage_{name}")()
                rec.status, rec.message = "ok", msg or ""
            except StageSkipped as exc:
                rec.status, rec.message = "skipped", str(exc)
                self.manifest.limit(f"stage '{name}' did not run: {exc}")
            except MappingError:
                rec.status = "failed"
                rec.seconds = round(time.time() - t0, 2)
                self.manifest.write()
                raise
            except EvodynError as exc:
                rec.status, rec.message = "failed", str(exc)
                rec.seconds = round(time.time() - t0, 2)
                self.manifest.write()
                raise
            rec.seconds = round(time.time() - t0, 2)
        # the report is written last so that it can describe every stage,
        # including the one that wrote the tables
        from evodyn.io.report import write_report
        write_report(self.manifest, self.state.summary, self.out / "report.md",
                     self.state.enrichment)
        self.manifest.write()
        return self.state

    # --------------------------------------------------------------- stages
    def _stage_sequences(self) -> str:
        from evodyn.sequences import (alignment_composition, cdhit_filter,
                                      fetch_pfam_alignment, residue_conservation,
                                      rescue_reference, stockholm_to_fasta)
        cfg, st = self.cfg, self.state

        # 1. the numbering map, which every later stage depends on
        t = cfg.target
        if t.reference_row:
            st.cmap = CanonicalMap.from_reference_row(
                t.reference_row, t.reference_start, t.reference_id)
        elif t.sequence:
            st.cmap = CanonicalMap.from_sequence(
                t.sequence, t.reference_start or (t.region[0] if t.region else 1),
                t.reference_id)
        if st.cmap and t.region:
            st.cmap = st.cmap.subset(*t.region)

        # 2. the alignment itself
        if cfg.sequences.alignment:
            st.alignment = pathlib.Path(cfg.sequences.alignment)
        elif cfg.sequences.source == "pfam" and t.pfam:
            sto = fetch_pfam_alignment(t.pfam, self.out / "sequences")
            fasta = self.out / "sequences" / f"{t.pfam}.full.fasta"
            n = stockholm_to_fasta(sto, fasta)
            clustered = self.out / "sequences" / f"{t.pfam}.clustered.fasta"
            cdhit_filter(fasta, clustered,
                         identity=cfg.sequences.cdhit_identity,
                         word_size=cfg.sequences.cdhit_word_size)
            final = self.out / "sequences" / f"{t.pfam}.final.fasta"
            rescued = (rescue_reference(clustered, t.reference_id, fasta, final)
                       if cfg.sequences.rescue_reference and t.reference_id
                       else bool(final.write_bytes(clustered.read_bytes())))
            st.alignment = final
            self.manifest.limit(
                f"alignment is family-wide ({t.pfam}); it is not restricted to the "
                "target's own subfamily unless the configured family is")
            if rescued:
                self.manifest.limit(
                    f"the reference sequence {t.reference_id} was not selected as a "
                    "cluster representative and was re-inserted into the alignment")
        elif cfg.sequences.source == "none":
            self.manifest.limit(
                "no alignment was read in this run: alignment composition and "
                "conservation-matched controls are unavailable")
        else:
            raise StageSkipped("no alignment supplied and no Pfam accession to fetch")

        # 3. composition and conservation, always
        msg = []
        if st.alignment and st.alignment.exists():
            comp = alignment_composition(st.alignment)
            comp.to_csv(self.out / "alignment_composition.csv", index=False)
            top = comp.iloc[0] if len(comp) else None
            if top is not None:
                self.manifest.limit(
                    f"alignment composition: {top.group} is the largest group at "
                    f"{top.percent:.1f}% of {comp.attrs['total_sequences']} sequences")
            if st.cmap is not None:
                st.conservation = residue_conservation(
                    st.alignment, st.cmap,
                    max_gap_fraction=cfg.sequences.max_gap_fraction)
                msg.append(f"conservation for {len(st.conservation)} residues")
            self.manifest.record_file(st.alignment)
        if st.cmap is not None:
            msg.append(f"{st.cmap.n_scored} scored positions "
                       f"(focus region spans {st.cmap.n_focus})")
            st.summary["n_scored_positions"] = st.cmap.n_scored
            st.summary["n_focus_residues"] = st.cmap.n_focus
        return "; ".join(msg)

    def _stage_coevolution(self) -> str:
        from evodyn.coevolution import read_couplings, run_plmc, translate_couplings
        cfg, st = self.cfg, self.state
        c = cfg.coevolution

        if c.couplings and c.couplings_format == "translated":
            st.couplings = pd.read_csv(c.couplings)
            missing = {"res_i", "res_j", "score"} - set(st.couplings.columns)
            if missing:
                raise EvodynError(f"{c.couplings}: missing column(s) {sorted(missing)}")
            st.couplings["sep"] = (st.couplings.res_i - st.couplings.res_j).abs()
        else:
            if c.couplings:
                raw = pathlib.Path(c.couplings)
            else:
                if c.method == "none":
                    raise StageSkipped("coevolution.method is 'none'")
                if st.alignment is None:
                    raise StageSkipped("no alignment available for DCA")
                self.manifest.record_binary(c.plmc_binary)
                raw = run_plmc(
                    st.alignment, self.out / "couplings.txt",
                    focus=cfg.target.reference_id,
                    out_params=self.out / "model.plmc",
                    lambda_e=c.lambda_e, lambda_h=c.lambda_h, max_iter=c.max_iter,
                    reweight=c.reweight_threshold, threads=c.threads,
                    binary=c.plmc_binary)
            table = read_couplings(raw, score_column=c.score_column)
            if st.cmap is None:
                raise EvodynError(
                    "translating couplings needs a numbering map; supply "
                    "target.reference_row or target.sequence")
            st.couplings = translate_couplings(table, st.cmap)
            self.manifest.record_file(raw)

        frac_neg = float((st.couplings.score < 0).mean())
        st.summary["n_coupled_pairs"] = len(st.couplings)
        st.summary["fraction_negative_scores"] = round(frac_neg, 3)
        if frac_neg > 0.2:
            self.manifest.limit(
                f"{100 * frac_neg:.1f}% of coupling scores are negative: the score is "
                "signed and approximately mean-zero, so it must not be summed over "
                "arbitrary residue-pair sets as though it were a magnitude")
        return f"{len(st.couplings)} pairs, {100 * frac_neg:.1f}% negative"

    def _stage_network(self) -> str:
        from evodyn.network import ResidueNetwork, knockout_table
        cfg, st = self.cfg, self.state
        if st.couplings is None:
            raise StageSkipped("no couplings")
        n = cfg.network
        st.network = ResidueNetwork.from_couplings(
            st.couplings,
            n_nodes=st.cmap.n_scored if st.cmap else None,
            rule=n.threshold, multiplier=n.multiplier,
            n_focus=st.cmap.n_focus if st.cmap else None,
            fixed_edges=n.fixed_edges, seeds=n.louvain_seeds)
        s = st.network.summary()
        s["louvain_seeds_tried"] = n.louvain_seeds
        st.summary.update({f"network_{k}": v for k, v in s.items()})

        st.network.node_table().to_csv(self.out / "network_nodes.csv", index=False)
        st.network.edge_table().to_csv(self.out / "network_edges.csv", index=False)
        st.network.intermodule_coupling().to_csv(
            self.out / "intermodule_coupling.csv", index=False)
        if n.knockout_targets:
            knockout_table(st.network.graph, n.knockout_targets).to_csv(
                self.out / "network_knockout.csv", index=False)
        if not st.network.connected:
            self.manifest.limit(
                "the network is not connected at the chosen threshold; path-based "
                "statistics are computed within components only")
        return (f"{s['nodes']} nodes, {s['edges']} edges, "
                f"{s['communities']} communities, Q = {s['weighted_modularity']}")

    def _stage_dynamics(self) -> str:
        from evodyn.dynamics import (dccm_from_trajectory, mean_matrix, read_xvg,
                                     replicate_agreement, rmsf_from_trajectory,
                                     run_equilibrator)
        cfg, st = self.cfg, self.state
        d = cfg.dynamics
        if d.mode == "none":
            self.manifest.limit(
                "no dynamic layer: this run reports evolutionary and structural "
                "information only")
            raise StageSkipped("dynamics.mode is 'none'")

        residues = st.cmap.residues if st.cmap else None
        offset = -(st.cmap.structure_offset or 0) if st.cmap else 0

        if d.mode == "matrix":
            m = pd.read_csv(d.dccm_matrix, index_col=0)
            m.index = m.index.astype(int); m.columns = m.columns.astype(int)
            st.dccm = m
            if d.rmsf:
                st.rmsf = read_xvg(d.rmsf)
            self.manifest.record_file(d.dccm_matrix)
            self.manifest.limit(
                "the dynamic layer was supplied as a precomputed matrix; EVODYN "
                "cannot verify the trajectory that produced it")
            return f"precomputed matrix, {len(m)} residues"

        if d.mode == "equilibrate":
            self.manifest.record_binary(d.equilibrator_binary)
            self.manifest.record_binary("gmx")
            structure = cfg.structure.pdb or cfg.target.structure
            runs = []
            for rep in range(1, d.replicates + 1):
                runs.append(run_equilibrator(
                    structure, self.out / "md", replicate=rep,
                    binary=d.equilibrator_binary, ligands=d.ligands,
                    force_field=d.force_field, water_model=d.water_model,
                    box_type=d.box_type, box_distance=d.box_distance,
                    pos_ion=d.pos_ion, neg_ion=d.neg_ion, extra_args=d.extra_args))
            for r in runs:
                st.dccm_replicates[r.replicate] = dccm_from_trajectory(
                    r.trajectory, r.structure, residues=residues, stride=d.stride)
            st.rmsf = read_xvg(runs[0].rmsf) if runs[0].rmsf else rmsf_from_trajectory(
                runs[0].trajectory, runs[0].structure, stride=d.stride)
            self.manifest.limit(
                f"molecular dynamics: {d.replicates} replicate(s) of the supplied "
                f"starting structure with {d.force_field}/{d.water_model}; no "
                "alternative conformational state, mutant or oligomer was simulated "
                "unless one was supplied as a separate run")
            note = f"{d.replicates} EquilibraTor replicate(s)"
        else:  # trajectory
            trajs = ([d.trajectory] if isinstance(d.trajectory, str) else list(d.trajectory))
            for rep, traj in enumerate(trajs, start=1):
                st.dccm_replicates[rep] = dccm_from_trajectory(
                    traj, d.topology, residues=residues, stride=d.stride, offset=offset)
                self.manifest.record_file(traj)
            st.rmsf = (read_xvg(d.rmsf) if d.rmsf else
                       rmsf_from_trajectory(trajs[0], d.topology, stride=d.stride,
                                            offset=offset))
            self.manifest.limit(
                "the trajectory was supplied by the user; its force field, length "
                "and starting state are outside EVODYN's record unless given in the "
                "configuration")
            note = f"{len(trajs)} supplied trajectory/ies"

        mats = st.dccm_replicates
        st.dccm = mean_matrix(list(mats.values())) if len(mats) > 1 else next(iter(mats.values()))
        if len(mats) > 1:
            agree = replicate_agreement(mats)
            agree.to_csv(self.out / "replicate_agreement.csv", index=False)
            st.summary["replicate_dccm_r_min"] = round(float(agree.pearson_r.min()), 3)
            st.summary["replicate_dccm_r_median"] = round(float(agree.pearson_r.median()), 3)
            self.manifest.limit(
                "replicates share a force field, a water model and a starting "
                "structure: their agreement measures precision, not accuracy")
        return note

    def _stage_structure(self) -> str:
        from evodyn.structure import ca_distance_matrix, contact_fraction, secondary_structure
        cfg, st = self.cfg, self.state
        pdb = cfg.structure.pdb or cfg.target.structure
        if not pdb:
            raise StageSkipped("no structure supplied; spatial controls unavailable")

        offset = 0
        if st.cmap is not None:
            offset = st.cmap.check_structure(pdb)
            if offset:
                self.manifest.limit(
                    f"structure {pathlib.Path(pdb).name} is numbered with a {offset:+d} "
                    "offset from canonical numbering; EVODYN applied it after checking "
                    "residue identities")
        residues = ([r for r in st.cmap.residues] if st.cmap else None)
        have = set(__import__("evodyn.numbering", fromlist=["x"])
                   .read_structure_residues(pdb))
        if residues:
            residues = [r for r in residues if (r - offset) in have]
        st.distances = ca_distance_matrix(pdb, residues=residues, offset=offset)
        self.manifest.record_file(pdb)

        msg = [f"{len(st.distances)} residues"]
        if st.network is not None:
            edges = st.network.edge_table()
            frac = contact_fraction(st.distances, edges, cfg.structure.contact_cutoff)
            med = float(np.median([st.distances.at[a, b] for a, b in
                                   zip(edges.res_i, edges.res_j)
                                   if a in st.distances.index and b in st.distances.index]))
            st.summary["network_edge_contact_fraction"] = round(frac, 3)
            st.summary["network_edge_median_distance_A"] = round(med, 2)
            msg.append(f"{100 * frac:.0f}% of network edges are contacts "
                       f"(median {med:.1f} A)")
            if frac > 0.8:
                self.manifest.limit(
                    f"{100 * frac:.0f}% of network edges are spatial contacts at "
                    f"{cfg.structure.contact_cutoff:.0f} A with a median length of "
                    f"{med:.1f} A: the network is close to a contact map, which is "
                    "expected behaviour for direct coupling analysis and bounds any "
                    "claim of long-range coupling")
        if cfg.structure.compute_dssp:
            try:
                st.secondary = secondary_structure(pdb, offset=offset)
                st.secondary.to_csv(self.out / "secondary_structure.csv", index=False)
                msg.append("DSSP assigned")
            except Exception as exc:  # mdtraj or dssp missing
                self.manifest.limit(f"secondary structure not computed: {exc}")
        return "; ".join(msg)

    def _stage_controls(self) -> str:
        from evodyn.controls import enrichment_panel
        from evodyn.io.tables import build_pair_table, build_residue_table
        cfg, st = self.cfg, self.state
        if st.couplings is None:
            raise StageSkipped("no couplings to test")

        st.residues = build_residue_table(
            cfg.target.name, st.cmap, network=st.network,
            conservation=st.conservation, rmsf=st.rmsf, secondary=st.secondary,
            sites=cfg.structure.sites)
        st.pairs = build_pair_table(
            cfg.target.name, st.couplings, distances=st.distances, dccm=st.dccm,
            network=st.network, residue_table=st.residues,
            contact_cutoff=cfg.structure.contact_cutoff)

        if "dccm" not in st.pairs:
            self.manifest.limit(
                "no dynamic layer, so no co-evolution/dynamics enrichment was tested")
            return "pair and residue tables assembled; no dynamic signal to test"

        c = cfg.controls
        st.enrichment = enrichment_panel(
            st.pairs, signal="dccm", score="score",
            sequence_separation=c.sequence_separation,
            distance_column="ca_distance" if c.distance_matched else None,
            conservation_column="conservation_min" if c.conservation_matched else None,
            ss_column=None, top_n=c.top_n_pairs, bin_width=c.distance_bin_width,
            n_draws=c.n_draws, n_bootstrap=c.n_bootstrap, seed=c.seed)
        st.enrichment.to_csv(self.out / "enrichment.csv", index=False)

        for _, row in st.enrichment.iterrows():
            st.summary[f"enrichment_{row['match']}"] = (
                f"{row['ratio']:.2f} (95% CI {row['ci_low']:.2f}-{row['ci_high']:.2f})")
        unmatched = st.enrichment[st.enrichment.match == "none"]
        matched = st.enrichment[st.enrichment.match != "none"]
        if len(unmatched) and len(matched):
            u = float(unmatched.ratio.iloc[0])
            for _, m in matched.iterrows():
                if m["ci_spans_one"] and u > 1.2:
                    self.manifest.limit(
                        f"the {u:.2f}-fold enrichment against an unmatched null falls to "
                        f"{m['ratio']:.2f} (95% CI {m['ci_low']:.2f}-{m['ci_high']:.2f}) "
                        f"when matched on {m['match']}: that covariate, not "
                        "co-evolution, accounts for the signal")
        return "; ".join(f"{r['match']}: {r['ratio']:.2f}" for _, r in st.enrichment.iterrows())

    def _stage_export(self) -> str:
        from evodyn.io.tables import write_table
        cfg, st = self.cfg, self.state
        written = []
        if cfg.export.per_residue and st.residues is not None:
            written.append(write_table(st.residues, self.out / "residues.csv",
                                       cfg.export.compress))
        if cfg.export.per_pair and st.pairs is not None:
            written.append(write_table(st.pairs, self.out / "pairs.csv",
                                       cfg.export.compress))
        if cfg.export.spheronizator:
            from evodyn.adapters import voxelize_structure
            pdb = cfg.structure.pdb or cfg.target.structure
            got = voxelize_structure(pdb, self.out / "ml") if pdb else None
            if got is None:
                self.manifest.limit("Spheronizator export requested but unavailable")
            else:
                written.append(got)
        if cfg.export.interactor:
            from evodyn.adapters import interaction_features
            pdb = cfg.structure.pdb or cfg.target.structure
            got = (interaction_features(pathlib.Path(pdb).parent,
                                        self.out / "ml" / "interactor.csv")
                   if pdb else None)
            if got is None:
                self.manifest.limit("InteracTor features requested but unavailable")
        for p in written:
            self.manifest.record_file(p)
        return ", ".join(p.name for p in written) or "report only"


def run_pipeline(cfg: EvodynConfig) -> RunState:
    return Pipeline(cfg).run()
