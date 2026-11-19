"""Statistical coupling analysis, optional.

SCA answers a different question from DCA - groups of positions that co-vary
collectively, rather than pairs with direct statistical dependence - and the
EGFR study used both. pySCA is not a hard dependency because it is not pip
installable in the usual way; when it is absent the stage is skipped and the
run manifest records the skip instead of the run silently losing a layer.

Sectors are reported with their spatial compactness and the nulls that test
it, because a sector that is compact only because its positions are
consecutive in sequence, or only because they are conserved, is not evidence
of a functional unit. In the EGFR analysis exactly one of six sectors survived
both nulls.
"""
from __future__ import annotations

import importlib.util

import numpy as np
import pandas as pd


def sca_available() -> bool:
    return importlib.util.find_spec("pysca") is not None


def run_sca(alignment, cmap, kpos: int = 6, p_cut: float = 0.95, seed: int = 0):
    """Run pySCA and return ``(sector_table, Csca)``.

    Raises
    ------
    ImportError
        If pySCA is not installed. The caller is expected to record a skip.
    """
    if not sca_available():
        raise ImportError(
            "pySCA not installed; see github.com/ranganathanlab/pySCA. "
            "Set coevolution.run_sca: false to skip this stage.")
    from pysca import scaTools as sca  # type: ignore

    np.random.seed(seed)
    headers, seqs = sca.readAlg(str(alignment))
    alg = sca.lett2num(seqs)
    _, _, _, _, Csca, _, Di = _sca_matrix(sca, alg)
    Vsca, Lsca = np.linalg.eigh(Csca)[::-1][0], None
    Vpica = sca.rotICA(Vsca, kmax=kpos)[0]
    ics, _ = sca.icList(Vpica, kpos, Csca, p_cut=p_cut)

    rows = []
    for k, ic in enumerate(ics, start=1):
        for pos in ic.items:
            rows.append({"sector": k, "column": int(pos) + 1,
                         "residue": cmap.col2res.get(int(pos) + 1)})
    return pd.DataFrame(rows), Csca


def _sca_matrix(sca, alg):
    """Thin indirection so the pySCA call site is in one place."""
    return sca.scaMat(alg)


def sector_compactness(sector_residues, distance_matrix, n_draws: int = 2000,
                       conservation: pd.Series | None = None, seed: int = 0) -> dict:
    """Mean pairwise Ca-Ca distance of a sector, against two nulls.

    The run-length-matched null places runs of the same sequence lengths at
    random, which controls for sectors that are compact because they are
    contiguous. The conservation-matched null resamples within deciles of
    conservation, which controls for sectors that are compact because
    conserved positions cluster in the core.
    """
    rng = np.random.default_rng(seed)
    res = [r for r in sector_residues if r in distance_matrix.index]
    if len(res) < 3:
        return {"n": len(res), "mean_distance": float("nan")}
    sub = distance_matrix.loc[res, res].to_numpy()
    iu = np.triu_indices(len(res), k=1)
    observed = float(sub[iu].mean())
    pool = list(distance_matrix.index)

    def _draw(sample):
        s = distance_matrix.loc[sample, sample].to_numpy()
        return float(s[np.triu_indices(len(sample), k=1)].mean())

    rand = [_draw(list(rng.choice(pool, len(res), replace=False))) for _ in range(n_draws)]
    out = {"n": len(res), "mean_distance": observed,
           "p_random": float((np.asarray(rand) <= observed).mean())}

    if conservation is not None:
        cons = conservation.reindex(pool).dropna()
        deciles = pd.qcut(cons, 10, labels=False, duplicates="drop")
        want = deciles.reindex(res).dropna().value_counts()
        draws = []
        for _ in range(n_draws):
            pick = []
            for d, k in want.items():
                cand = deciles[deciles == d].index.to_numpy()
                pick.extend(rng.choice(cand, size=min(k, len(cand)), replace=False))
            if len(pick) >= 3:
                draws.append(_draw(list(pick)))
        if draws:
            out["p_conservation_matched"] = float((np.asarray(draws) <= observed).mean())
    return out
