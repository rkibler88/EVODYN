"""Molecular dynamics, delegated to EquilibraTor.

EVODYN runs MD rather than asking the user for a trajectory, but it does not
reimplement the setup. The simulation workflow is EquilibraTor's, by Jose D. D.
Cediel-Becerra, Jose Cleydson F. Silva and Raquel Dias (Dias Lab, University of
Florida), *Computational and Structural Biotechnology Journal* **27**,
5318-5325 (2025), doi:10.1016/j.csbj.2025.11.034,
github.com/Dias-Lab/EquilibraTor. It drives GROMACS from topology generation
through box construction, solvation, ion placement, energy minimisation, NVT,
NPT, convergence detection and the production stage, for proteins and
protein-ligand complexes.

EVODYN supplies the structure and the parameters, runs one invocation per
replicate in its own working directory, and collects the trajectory, the RMSF
profile and the convergence diagnostics. **If you use this stage, cite
EquilibraTor.**

Three things this adapter deliberately does:

* **one working directory per replicate.** EquilibraTor names its project
  directory after the input structure, so concurrent replicates would collide.
* **no silent fallback.** If GROMACS or EquilibraTor is missing, the stage
  fails with an instruction rather than quietly producing a run with no
  dynamic layer. Users who want the evolutionary layers alone set
  ``dynamics.mode: none`` and the limits section records it.
* **replicates are independent runs, and the manifest says what they share.**
  Agreement between replicates that share a force field and a starting model
  is precision, not accuracy. The EGFR study made that explicit after review
  and EVODYN writes the same caveat automatically.
"""
from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
from dataclasses import dataclass, field

from evodyn.exceptions import ExternalToolError


@dataclass
class EquilibratorRun:
    """Outputs of one EquilibraTor invocation."""

    replicate: int
    workdir: pathlib.Path
    project_dir: pathlib.Path
    trajectory: pathlib.Path | None = None
    structure: pathlib.Path | None = None
    rmsf: pathlib.Path | None = None
    log: pathlib.Path | None = None
    returncode: int = 0
    command: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.returncode == 0 and self.trajectory is not None


def equilibrator_available(binary: str = "EquilibraTor") -> bool:
    return shutil.which(binary) is not None


def build_command(
    structure,
    binary: str = "EquilibraTor",
    ligands: list[str] | None = None,
    force_field: str = "amber99sb-ildn",
    water_model: str = "tip3p",
    box_type: str = "cubic",
    box_distance: float = 1.2,
    pos_ion: str = "NA",
    neg_ion: str = "CL",
    first_step: int | None = None,
    last_step: int | None = None,
    extra_args: list[str] | None = None,
) -> list[str]:
    """Assemble the EquilibraTor argv.

    Flag names follow EquilibraTor v1.0.0 (``-p``, ``-l``, ``-gff``, ``-gwm``,
    ``-gbt``, ``-gd``, ``-gpi``, ``-gni``, ``-fs``, ``-ls``).
    """
    cmd = [binary, "-p", str(structure)]
    if ligands:
        cmd += ["-l", *[str(l) for l in ligands]]
    cmd += ["-gff", force_field, "-gwm", water_model,
            "-gbt", box_type, "-gd", str(box_distance),
            "-gpi", pos_ion, "-gni", neg_ion]
    if first_step is not None:
        cmd += ["-fs", str(first_step)]
    if last_step is not None:
        cmd += ["-ls", str(last_step)]
    if extra_args:
        cmd += list(extra_args)
    return cmd


def run_equilibrator(
    structure,
    outdir,
    replicate: int = 1,
    binary: str = "EquilibraTor",
    dry_run: bool = False,
    timeout: float | None = None,
    **kwargs,
) -> EquilibratorRun:
    """Run one replicate and locate its production outputs.

    Parameters
    ----------
    dry_run
        Assemble the command, create the directory and return without
        executing. Used by the test suite and by ``evodyn plan``, so that a
        user can see exactly what would be submitted before spending GPU hours.
    """
    structure = pathlib.Path(structure).resolve()
    workdir = pathlib.Path(outdir).resolve() / f"replicate_{replicate:02d}"
    workdir.mkdir(parents=True, exist_ok=True)
    stem = structure.stem
    ligands = kwargs.get("ligands") or []
    project = workdir / (f"{stem}_" + "_".join(pathlib.Path(l).stem for l in ligands)
                         if ligands else stem)

    cmd = build_command(structure, binary=binary, **kwargs)
    run = EquilibratorRun(replicate, workdir, project, command=cmd)
    if dry_run:
        return run

    if not equilibrator_available(binary):
        raise ExternalToolError(
            f"{binary} not found on PATH. Install it with `pip install EquilibraTor` "
            "and make sure GROMACS is loaded, or set dynamics.mode to 'trajectory' "
            "to supply your own simulation."
        )
    log = workdir / "equilibrator.log"
    with open(log, "w") as fh:
        proc = subprocess.run(cmd, cwd=workdir, stdout=fh,
                              stderr=subprocess.STDOUT, timeout=timeout)
    run.returncode, run.log = proc.returncode, log
    if proc.returncode != 0:
        raise ExternalToolError(
            f"{binary} exited {proc.returncode} for replicate {replicate}; see {log}")

    prod = project / "production_output"
    run.trajectory = _first_existing(prod / "production_stage.xtc",
                                     prod / "production_stage.trr")
    run.structure = _first_existing(prod / "production.pdb",
                                    prod / "production_last_frame.pdb")
    run.rmsf = _first_existing(prod / "production_rmsf.xvg")
    if run.trajectory is None:
        raise ExternalToolError(
            f"{binary} finished but no production trajectory was found under {prod}")
    return run


def _first_existing(*paths) -> pathlib.Path | None:
    for p in paths:
        if pathlib.Path(p).exists():
            return pathlib.Path(p)
    return None
