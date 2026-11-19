"""Regression test against the EGFR case study.

EVODYN generalises a published single-protein analysis, so the first thing it
has to demonstrate is that the generalisation did not change the answer. This
module runs the whole pipeline over the committed EGFR intermediates and
compares what comes out against the values in the manuscript and in
``EGFRAllostery/analysis/verify_numbers.py``.

Two classes of check are distinguished:

**exact**
    Properties of the data and of deterministic algorithms: the number of
    scored positions, the Louvain partition and its modularity, the
    percolation threshold, residue identities at named motifs.

**tolerance**
    Statistics that depend on a random draw. EVODYN's null sampler is not
    bit-identical to the one in the original analysis - it sorts the pair table
    deterministically before drawing, where the original consumed the table in
    file order - so these are compared within a stated tolerance. The
    scientific content is unchanged: the unmatched enrichment is large, the
    distance-matched enrichment is not distinguishable from 1.0.
"""
from __future__ import annotations

import pathlib
from dataclasses import dataclass

import numpy as np
import pandas as pd

DATA = pathlib.Path(__file__).resolve().parents[2] / "examples" / "egfr" / "data"
CONFIG = pathlib.Path(__file__).resolve().parents[2] / "examples" / "egfr" / "config.yaml"

#: published values, EGFRAllostery analysis/verify_numbers.py and the manuscript
PUBLISHED = {
    "scored_positions": 253,
    "focus_residues": 257,
    "coupled_pairs": 31878,
    "network_edges": 462,
    "communities": 15,
    "weighted_modularity": 0.8036,
    "percolation_edges": 458,
    "residue_853": "I",
    "motif_855": "DFG",
    "gatekeeper_790": "T",
    "enrichment_unmatched": 1.83,
    "enrichment_distance_matched": 1.02,
    "edge_contact_fraction": 0.92,
    "edge_median_distance_A": 3.9,
}


@dataclass
class Check:
    name: str
    expected: object
    got: object
    tolerance: float | None = None

    @property
    def ok(self) -> bool:
        if self.tolerance is None:
            return self.expected == self.got
        try:
            return abs(float(self.expected) - float(self.got)) <= self.tolerance
        except (TypeError, ValueError):
            return False

    def __str__(self) -> str:
        tol = f" +/- {self.tolerance}" if self.tolerance else ""
        return (f"{'PASS' if self.ok else 'FAIL'}  {self.name:<34} "
                f"expected {self.expected}{tol}  got {self.got}")


def run_benchmark(data_dir=None, verbose: bool = True) -> bool:
    """Run the EGFR configuration and compare against the published values."""
    from evodyn.config import load_config
    from evodyn.network.build import percolation_threshold
    from evodyn.pipeline import Pipeline

    cfg = load_config(CONFIG)
    if data_dir:
        d = pathlib.Path(data_dir)
        cfg.target.structure = str(d / "egfr_clean.pdb")
        cfg.target.reference_row = str(d / "reference_aligned.fasta")
        cfg.structure.pdb = str(d / "egfr_clean.pdb")
        cfg.coevolution.couplings = str(d / "couplings_translated.csv.gz")
        cfg.dynamics.dccm_matrix = str(d / "dccm_matrix.csv.gz")
        cfg.dynamics.rmsf = str(d / "rmsf.xvg")
    pipe = Pipeline(cfg)
    st = pipe.run()
    s = st.summary
    cmap = st.cmap

    enr = st.enrichment.set_index("match")
    checks = [
        Check("scored positions", PUBLISHED["scored_positions"], s["n_scored_positions"]),
        Check("focus residues", PUBLISHED["focus_residues"], s["n_focus_residues"]),
        Check("coupled pairs", PUBLISHED["coupled_pairs"], s["n_coupled_pairs"]),
        Check("network nodes", PUBLISHED["scored_positions"], s["network_nodes"]),
        Check("network edges", PUBLISHED["network_edges"], s["network_edges"]),
        Check("communities", PUBLISHED["communities"], s["network_communities"]),
        Check("weighted modularity", PUBLISHED["weighted_modularity"],
              s["network_weighted_modularity"], 0.0005),
        Check("percolation threshold", PUBLISHED["percolation_edges"],
              percolation_threshold(st.couplings, cmap.n_scored)),
        Check("residue 853 identity", PUBLISHED["residue_853"], cmap.res2aa[853]),
        Check("DFG motif at 855-857", PUBLISHED["motif_855"],
              "".join(cmap.res2aa[r] for r in (855, 856, 857))),
        Check("gatekeeper 790 identity", PUBLISHED["gatekeeper_790"], cmap.res2aa[790]),
        Check("edge contact fraction", PUBLISHED["edge_contact_fraction"],
              s["network_edge_contact_fraction"], 0.02),
        Check("edge median distance", PUBLISHED["edge_median_distance_A"],
              s["network_edge_median_distance_A"], 0.2),
        Check("enrichment, unmatched null", PUBLISHED["enrichment_unmatched"],
              round(float(enr.loc["none", "ratio"]), 2), 0.10),
        Check("enrichment, distance-matched", PUBLISHED["enrichment_distance_matched"],
              round(float(enr.loc["ca_distance", "ratio"]), 2), 0.08),
        Check("distance-matched CI spans 1.0", True,
              bool(enr.loc["ca_distance", "ci_spans_one"])),
    ]
    if verbose:
        print(f"EVODYN benchmark: EGFR kinase domain "
              f"({PUBLISHED['scored_positions']} positions)\n")
        for c in checks:
            print(c)
        n_ok = sum(c.ok for c in checks)
        print(f"\n{n_ok}/{len(checks)} checks passed")
        if n_ok < len(checks):
            print("FAILED: " + ", ".join(c.name for c in checks if not c.ok))
    return all(c.ok for c in checks)
