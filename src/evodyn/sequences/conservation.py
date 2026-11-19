"""Per-column conservation, for the conservation-matched null.

The EGFR study could not recompute its conservation-matched null from a clone
because the pySCA positional conservation vector was never committed. EVODYN
computes conservation itself, from the same alignment the couplings came from,
so the control is always reproducible.

Conservation is relative entropy to a background amino-acid distribution,
which is the quantity pySCA's ``Di`` also expresses:

.. math::
    D_i = \\sum_a f_{ia} \\ln \\frac{f_{ia}}{q_a}
"""
from __future__ import annotations

import numpy as np
import pandas as pd

AA = "ACDEFGHIKLMNPQRSTVWY"

#: background frequencies, UniProt/Swiss-Prot composition
BACKGROUND = np.array([
    0.0825, 0.0137, 0.0545, 0.0675, 0.0386, 0.0707, 0.0227, 0.0596, 0.0584,
    0.0966, 0.0242, 0.0406, 0.0470, 0.0393, 0.0553, 0.0656, 0.0534, 0.0687,
    0.0108, 0.0292,
])


def column_conservation(
    fasta, columns: list[int] | None = None, max_gap_fraction: float = 0.5,
    pseudocount: float = 1e-6,
) -> pd.DataFrame:
    """Relative entropy and gap fraction for each alignment column.

    Columns are 1-based, matching the convention used by the numbering map.
    """
    counts: np.ndarray | None = None
    gaps: np.ndarray | None = None
    n_seq = 0
    index = {a: k for k, a in enumerate(AA)}
    for seq in _iter_sequences(fasta):
        if counts is None:
            counts = np.zeros((len(seq), 20))
            gaps = np.zeros(len(seq))
        n_seq += 1
        arr = np.frombuffer(seq.upper().encode(), dtype="S1")
        for k, ch in enumerate(seq.upper()):
            j = index.get(ch)
            if j is None:
                gaps[k] += 1
            else:
                counts[k, j] += 1
    if counts is None:
        raise ValueError(f"{fasta}: no sequences")

    tot = counts.sum(axis=1, keepdims=True)
    freq = np.divide(counts, np.maximum(tot, 1.0))
    freq = np.clip(freq, pseudocount, None)
    freq /= freq.sum(axis=1, keepdims=True)
    rel_ent = (freq * np.log(freq / BACKGROUND)).sum(axis=1)

    df = pd.DataFrame({
        "column": np.arange(1, counts.shape[0] + 1),
        "conservation": rel_ent,
        "gap_fraction": gaps / max(n_seq, 1),
        "n_sequences": n_seq,
    })
    df.loc[df.gap_fraction > max_gap_fraction, "conservation"] = np.nan
    return df.loc[df.column.isin(columns)].reset_index(drop=True) if columns else df


def residue_conservation(fasta, cmap, **kwargs) -> pd.DataFrame:
    """Column conservation mapped onto canonical residue numbers."""
    cols = column_conservation(fasta, **kwargs)
    cols["residue"] = cols["column"].map(cmap.col2res)
    return (cols.dropna(subset=["residue"])
            .astype({"residue": int})
            [["residue", "conservation", "gap_fraction"]]
            .reset_index(drop=True))


def _iter_sequences(fasta):
    buf = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if buf:
                    yield "".join(buf)
                    buf = []
            else:
                buf.append(line.strip())
    if buf:
        yield "".join(buf)
