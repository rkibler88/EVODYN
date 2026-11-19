"""Physics-based interaction features, via InteracTor.

InteracTor is by Jose Cleydson F. Silva, Layla Schuster, Nick Sexson, Melissa
Erdem, Ryan Hulk, Matias Kirst, Marcio F. R. Resende and Raquel Dias (Dias Lab,
University of Florida), doi:10.1101/2025.04.10.648139 and
doi:10.1101/2024.10.07.616705. Cite it if you use this export.

EVODYN's structural layer is geometric: Ca-Ca distances, contacts, secondary
structure. InteracTor extracts interaction terms - hydrogen bonds, hydrophobic
contacts, van der Waals and electrostatic terms, solvent accessibility - from
a structure and its MOL2 conversion. Where it is installed, those terms are
attached to the run as an additional structural feature block, which is the
natural complement to the evolutionary and dynamic layers for a
structure-function model.
"""
from __future__ import annotations

import importlib.util
import pathlib
import shutil
import subprocess

import pandas as pd


def interactor_available() -> bool:
    return (importlib.util.find_spec("InteracTor") is not None
            or shutil.which("InteracTor") is not None)


def interaction_features(structure_dir, out_csv, binary: str = "InteracTor"):
    """Run InteracTor over a directory of structures.

    InteracTor expects PDB files with matching MOL2 conversions in the same
    directory (``obabel -ipdb x.pdb -omol2 > x.pdb.mol2``). Returns the parsed
    feature table, or ``None`` if the tool is absent.
    """
    if shutil.which(binary) is None:
        return None
    out_csv = pathlib.Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run([binary, str(structure_dir), str(out_csv)],
                          capture_output=True, text=True)
    if proc.returncode != 0 or not out_csv.exists():
        return None
    return pd.read_csv(out_csv)
