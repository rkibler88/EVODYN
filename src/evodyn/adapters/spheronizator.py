"""Voxel export, via Spheronizator.

Spheronizator is by Matthew Richardson, Jose Cleydson F. Silva, Jose D. D.
Cediel-Becerra and Raquel Dias (Dias Lab, University of Florida),
github.com/Dias-Lab/spheronizator. Cite it if you use this export.

EVODYN's output is tabular: per-residue and per-pair features in canonical
numbering. A model that learns from local geometry wants the structure itself
in a fixed-size grid. Spheronizator voxelises a PDB for exactly that purpose,
so the export stage can emit a voxel tensor beside the tables, keyed on the
same residue numbering, and a downstream study can join them.
"""
from __future__ import annotations

import importlib.util
import pathlib


def spheronizator_available() -> bool:
    return importlib.util.find_spec("spheronizator") is not None


def voxelize_structure(pdb, outdir, **kwargs):
    """Voxelise a structure and write the grid next to the run's tables.

    Returns the output path, or ``None`` when Spheronizator is not installed -
    the caller records a skip rather than failing the run.
    """
    if not spheronizator_available():
        return None
    import spheronizator  # type: ignore

    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dest = outdir / (pathlib.Path(pdb).stem + "_voxels.npz")
    grid = spheronizator.voxelize(str(pdb), **kwargs) if hasattr(spheronizator, "voxelize") else None
    if grid is None:
        return None
    import numpy as np
    np.savez_compressed(dest, voxels=grid)
    return dest
