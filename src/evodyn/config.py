"""Run configuration.

A run is fully described by one YAML file. Nothing that affects a number is
read from the environment, from argv beyond ``--config``, or from a hardcoded
constant inside a stage - the EGFR working directory accumulated scripts with
edge counts and file names baked in, and reconciling them afterwards was most
of the provenance audit.
"""
from __future__ import annotations

import pathlib
from dataclasses import dataclass, field, asdict
from typing import Any, Literal

import yaml

from evodyn.exceptions import ConfigError


@dataclass
class TargetConfig:
    """The protein the run is about."""

    name: str
    uniprot: str | None = None
    pfam: str | None = None
    #: inclusive canonical residue range analysed, e.g. [714, 966]
    region: list[int] | None = None
    #: reference sequence identifier inside the family alignment
    reference_id: str | None = None
    structure: str | None = None
    #: aligned reference row (FASTA). Uppercase = match state, lowercase =
    #: insert state. This is the column -> residue map the whole run uses.
    reference_row: str | None = None
    #: canonical number of the first residue of ``reference_row``
    reference_start: int | None = None
    #: ungapped reference sequence, as an alternative to ``reference_row``
    sequence: str | None = None


@dataclass
class SequencesConfig:
    alignment: str | None = None          # pre-built alignment, skips retrieval
    #: "pfam" downloads from InterPro, "local" uses ``alignment``, "none"
    #: skips sequence retrieval entirely (numbering map only - used when the
    #: couplings are already computed)
    source: Literal["pfam", "local", "none"] = "pfam"
    cdhit_identity: float = 0.62
    cdhit_word_size: int = 4
    rescue_reference: bool = True
    max_gap_fraction: float = 0.5


@dataclass
class CoevolutionConfig:
    method: Literal["plmc", "none"] = "plmc"
    plmc_binary: str = "plmc"
    lambda_e: float = 1.0
    lambda_h: float = 0.01
    max_iter: int = 500
    reweight_threshold: float = 0.2
    threads: int = 8
    #: column of the plmc output to use. Column 5 (Frobenius norm) is zero in
    #: some plmc builds; column 6 is the average-product-corrected score.
    score_column: int = 6
    run_sca: bool = False
    #: a coupling file that already exists. Set this to analyse published
    #: couplings, or to re-run downstream stages without repeating DCA.
    couplings: str | None = None
    #: "plmc" for raw plmc output, "translated" for a CSV already carrying
    #: res_i, res_j and score in canonical numbering.
    couplings_format: Literal["plmc", "translated"] = "plmc"


@dataclass
class NetworkConfig:
    #: "percolation" locates the smallest edge count at which the graph spans
    #: every scored position and is connected. "multiplier" uses int(L * value).
    threshold: Literal["percolation", "multiplier", "fixed"] = "percolation"
    multiplier: float = 1.80
    fixed_edges: int | None = None
    louvain_seeds: int = 100
    knockout_targets: list[int] = field(default_factory=list)


@dataclass
class DynamicsConfig:
    """How the dynamic layer is obtained.

    ``mode`` is the decision that bounds every dynamic claim a run can make:

    ``equilibrate``
        Run molecular dynamics inside EVODYN by delegating to EquilibraTor
        (Dias Lab), which drives GROMACS from topology generation through
        energy minimisation, NVT, NPT, convergence detection and production.
    ``trajectory``
        Consume a trajectory the user already has. Intended for testing, for
        reusing published simulations, and for anyone whose MD is run by a
        facility rather than by this pipeline.
    ``matrix``
        Consume a precomputed cross-correlation matrix.
    ``none``
        Skip the dynamic layer. The run reports evolutionary and structural
        layers only, and says so in its limits section.
    """

    mode: Literal["equilibrate", "trajectory", "matrix", "none"] = "trajectory"
    replicates: int = 1
    # mode="equilibrate"
    force_field: str = "amber99sb-ildn"
    water_model: str = "tip3p"
    box_type: str = "cubic"
    box_distance: float = 1.2
    pos_ion: str = "NA"
    neg_ion: str = "CL"
    ligands: list[str] = field(default_factory=list)
    equilibrator_binary: str = "EquilibraTor"
    extra_args: list[str] = field(default_factory=list)
    # mode="trajectory"
    trajectory: str | None = None
    topology: str | None = None
    stride: int = 1
    # mode="matrix"
    dccm_matrix: str | None = None
    rmsf: str | None = None


@dataclass
class StructureConfig:
    pdb: str | None = None
    contact_cutoff: float = 8.0
    compute_dssp: bool = True
    #: named residue sets (orthosteric site, allosteric pocket, interface ...)
    sites: dict[str, list[int]] = field(default_factory=dict)


@dataclass
class ControlsConfig:
    """Null models. These are not optional switches on an analysis; they are
    the analysis. An enrichment without its matched null is not reported."""

    distance_matched: bool = True
    conservation_matched: bool = True
    ss_matched: bool = False
    sequence_separation: int = 5
    top_n_pairs: int = 500
    distance_bin_width: float = 2.0
    n_draws: int = 200
    n_bootstrap: int = 2000
    seed: int = 0


@dataclass
class ExportConfig:
    per_residue: bool = True
    per_pair: bool = True
    compress: bool = True
    spheronizator: bool = False
    interactor: bool = False


@dataclass
class EvodynConfig:
    target: TargetConfig
    outdir: str = "evodyn_run"
    sequences: SequencesConfig = field(default_factory=SequencesConfig)
    coevolution: CoevolutionConfig = field(default_factory=CoevolutionConfig)
    network: NetworkConfig = field(default_factory=NetworkConfig)
    dynamics: DynamicsConfig = field(default_factory=DynamicsConfig)
    structure: StructureConfig = field(default_factory=StructureConfig)
    controls: ControlsConfig = field(default_factory=ControlsConfig)
    export: ExportConfig = field(default_factory=ExportConfig)
    stages: list[str] = field(default_factory=lambda: [
        "sequences", "coevolution", "network", "dynamics",
        "structure", "controls", "export",
    ])

    # ---------------------------------------------------------------- checks
    def validate(self) -> None:
        if self.target.region is not None and len(self.target.region) != 2:
            raise ConfigError("target.region must be [first, last] canonical residue numbers")
        if self.target.region and self.target.region[0] > self.target.region[1]:
            raise ConfigError("target.region is reversed")
        d = self.dynamics
        if d.mode == "trajectory" and not d.trajectory:
            raise ConfigError("dynamics.mode='trajectory' requires dynamics.trajectory")
        if d.mode == "matrix" and not d.dccm_matrix:
            raise ConfigError("dynamics.mode='matrix' requires dynamics.dccm_matrix")
        if d.mode == "equilibrate" and not (self.structure.pdb or self.target.structure):
            raise ConfigError(
                "dynamics.mode='equilibrate' requires a starting structure "
                "(structure.pdb or target.structure)"
            )
        if d.replicates < 1:
            raise ConfigError("dynamics.replicates must be >= 1")
        if self.network.threshold == "fixed" and not self.network.fixed_edges:
            raise ConfigError("network.threshold='fixed' requires network.fixed_edges")
        if (self.coevolution.method == "none" and not self.coevolution.couplings
                and "network" in self.stages):
            raise ConfigError(
                "the network stage needs either coevolution.method='plmc' or a "
                "precomputed coevolution.couplings file")
        if self.target.reference_row and self.target.reference_start is None:
            raise ConfigError("target.reference_row requires target.reference_start")
        if not 0 < self.sequences.cdhit_identity <= 1:
            raise ConfigError("sequences.cdhit_identity must be in (0, 1]")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def out(self) -> pathlib.Path:
        return pathlib.Path(self.outdir)


_SECTIONS = {
    "sequences": SequencesConfig,
    "coevolution": CoevolutionConfig,
    "network": NetworkConfig,
    "dynamics": DynamicsConfig,
    "structure": StructureConfig,
    "controls": ControlsConfig,
    "export": ExportConfig,
}


#: fields holding input paths, resolved relative to the configuration file so
#: that a config is portable and can be run from any working directory
_PATH_FIELDS = {
    "target": ("structure", "reference_row"),
    "sequences": ("alignment",),
    "coevolution": ("couplings",),
    "dynamics": ("trajectory", "topology", "dccm_matrix", "rmsf"),
    "structure": ("pdb",),
}


def _resolve_paths(raw: dict, base: pathlib.Path) -> dict:
    for section, fields in _PATH_FIELDS.items():
        block = raw.get(section)
        if not isinstance(block, dict):
            continue
        for f in fields:
            v = block.get(f)
            if isinstance(v, str):
                q = pathlib.Path(v)
                block[f] = str(q if q.is_absolute() else (base / q).resolve())
            elif isinstance(v, list):
                block[f] = [str(pathlib.Path(x) if pathlib.Path(x).is_absolute()
                                else (base / x).resolve()) for x in v]
    return raw


def load_config(path: str | pathlib.Path) -> EvodynConfig:
    """Read and validate a run configuration.

    Input paths are interpreted relative to the configuration file, so a run
    directory can be moved or cloned without editing every path in it.
    """
    path = pathlib.Path(path).resolve()
    raw = yaml.safe_load(path.read_text()) or {}
    raw = _resolve_paths(raw, path.parent)
    if "target" not in raw:
        raise ConfigError(f"{path}: no 'target' section")

    def build(cls, data, label):
        known = {f for f in cls.__dataclass_fields__}
        unknown = set(data) - known
        if unknown:
            raise ConfigError(f"{label}: unknown key(s) {sorted(unknown)}")
        return cls(**data)

    cfg = EvodynConfig(
        target=build(TargetConfig, raw["target"], "target"),
        outdir=raw.get("outdir", "evodyn_run"),
        stages=raw.get("stages", EvodynConfig.__dataclass_fields__["stages"].default_factory()),
        **{k: build(c, raw.get(k, {}) or {}, k) for k, c in _SECTIONS.items()},
    )
    cfg.validate()
    return cfg
