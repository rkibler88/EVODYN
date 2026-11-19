"""Direct coupling analysis with plmc.

plmc writes one line per scored pair::

    col_i  aa_i  col_j  aa_j  frobenius_norm  apc_score

Two properties of that output cause recurring errors and are handled here
rather than left to the user:

* the Frobenius-norm column is **zero throughout** in some builds. EVODYN
  checks the chosen score column for degeneracy and raises instead of
  building a network out of zeros.
* the APC-corrected score is **signed and approximately mean-zero**. It is a
  rank statistic, not a magnitude. :func:`read_couplings` records the fraction
  of negative scores so the run report can state it.

The columns plmc reports are *match-state* positions of the alignment, not
residue numbers. :func:`translate_couplings` converts them through the
numbering map, and refuses to emit a pair whose columns do not both map.
"""
from __future__ import annotations

import pathlib
import shutil
import subprocess

import numpy as np
import pandas as pd

from evodyn.exceptions import ExternalToolError, MappingError
from evodyn.numbering import CanonicalMap


def run_plmc(
    alignment, out_couplings, focus: str, out_params=None,
    lambda_e: float = 1.0, lambda_h: float = 0.01, max_iter: int = 500,
    reweight: float = 0.2, threads: int = 8, binary: str = "plmc",
) -> pathlib.Path:
    """Run pseudo-likelihood maximisation DCA.

    ``focus`` is the identifier of the reference row; plmc restricts scoring to
    the columns where that sequence has a residue, which is what makes the
    output interpretable in the reference protein's own numbering.
    """
    if shutil.which(binary) is None:
        raise ExternalToolError(
            f"{binary} not found on PATH. Build it from github.com/debbiemarkslab/plmc "
            "or set coevolution.method to 'none'.")
    cmd = [binary, "-c", str(out_couplings), "-f", focus,
           "-le", str(lambda_e), "-lh", str(lambda_h), "-m", str(max_iter),
           "-t", str(reweight), "-n", str(threads)]
    if out_params:
        cmd += ["-o", str(out_params)]
    cmd.append(str(alignment))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise ExternalToolError(f"plmc failed: {(proc.stderr or proc.stdout)[-2000:]}")
    return pathlib.Path(out_couplings)


def read_couplings(path, score_column: int = 6) -> pd.DataFrame:
    """Read plmc output into a table of columns, residue identities and score."""
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=["col_i", "aa_i", "col_j", "aa_j", "fn", "apc"])
    col = {5: "fn", 6: "apc"}.get(score_column)
    if col is None:
        raise ValueError("score_column must be 5 (Frobenius norm) or 6 (APC)")
    if np.allclose(df[col].to_numpy(), 0.0):
        other = "apc" if col == "fn" else "fn"
        raise ValueError(
            f"{path}: column {score_column} is zero for every pair. "
            f"This plmc build reports its scores in the '{other}' column; "
            "set coevolution.score_column accordingly.")
    out = df.rename(columns={col: "score"})[
        ["col_i", "aa_i", "col_j", "aa_j", "score"]].copy()
    out.attrs["fraction_negative"] = float((out.score < 0).mean())
    out.attrs["score_column"] = score_column
    return out


def translate_couplings(couplings: pd.DataFrame, cmap: CanonicalMap) -> pd.DataFrame:
    """Map alignment columns to canonical residue numbers, checking identities.

    Both the column mapping and the amino-acid identity are verified. A pair
    whose reported residue letter disagrees with the reference row means the
    alignment and the numbering map describe different sequences, and the run
    stops.
    """
    res_i = couplings["col_i"].map(cmap.col2res)
    res_j = couplings["col_j"].map(cmap.col2res)
    unmapped = int(res_i.isna().sum() + res_j.isna().sum())
    df = couplings.assign(res_i=res_i, res_j=res_j).dropna(subset=["res_i", "res_j"])
    df = df.astype({"res_i": int, "res_j": int})

    bad = 0
    for col_res, col_aa in (("res_i", "aa_i"), ("res_j", "aa_j")):
        ref = df[col_res].map(cmap.res2aa)
        bad += int((ref.str.upper() != df[col_aa].str.upper()).sum())
    if bad:
        raise MappingError(
            f"{bad} residue identities in the coupling file disagree with the "
            "alignment reference row. The numbering map and the alignment are "
            "not describing the same sequence.")

    out = df[["res_i", "res_j", "score", "aa_i", "aa_j", "col_i", "col_j"]].copy()
    out["sep"] = (out.res_i - out.res_j).abs()
    out.attrs.update(couplings.attrs)
    out.attrs["unmapped_column_references"] = unmapped
    out.attrs["n_scored_positions"] = len(set(out.res_i) | set(out.res_j))
    return out.reset_index(drop=True)
