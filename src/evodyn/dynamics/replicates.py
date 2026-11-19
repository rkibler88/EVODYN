"""Agreement between replicate trajectories.

Reported as a bound, not as a validation. Replicates that share a force field,
a water model and a starting model can only tell you how reproducible the
sampling is; they cannot tell you the ensemble is right. EVODYN writes that
sentence into the limits section whenever more than one replicate is present.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def mean_matrix(matrices: list[pd.DataFrame]) -> pd.DataFrame:
    """Element-wise mean over replicate matrices sharing an index."""
    idx = matrices[0].index
    for m in matrices[1:]:
        if not m.index.equals(idx):
            raise ValueError("replicate matrices do not share a residue index")
    return pd.DataFrame(np.mean([m.to_numpy() for m in matrices], axis=0),
                        index=idx, columns=matrices[0].columns)


def replicate_agreement(matrices: dict[int, pd.DataFrame]) -> pd.DataFrame:
    """Pearson r between the off-diagonal elements of every replicate pair."""
    keys = sorted(matrices)
    n = len(matrices[keys[0]])
    iu = np.triu_indices(n, k=1)
    vec = {k: matrices[k].to_numpy()[iu] for k in keys}
    rows = [{"replicate_i": a, "replicate_j": b, "pearson_r": float(pearsonr(vec[a], vec[b])[0])}
            for a, b in combinations(keys, 2)]
    return pd.DataFrame(rows)


def profile_agreement(profiles: dict[int, pd.Series]) -> pd.DataFrame:
    """Pearson r between replicate per-residue profiles (RMSF, for example)."""
    keys = sorted(profiles)
    rows = []
    for a, b in combinations(keys, 2):
        x, y = profiles[a].align(profiles[b], join="inner")
        rows.append({"replicate_i": a, "replicate_j": b,
                     "pearson_r": float(pearsonr(x, y)[0])})
    return pd.DataFrame(rows)
