"""Matched null models.

This module exists because of a specific result. In the EGFR study, the top
co-evolving distal pairs showed 1.83-fold enrichment in correlated motion over
a null matched on *sequence* separation - a large and publishable effect. The
same comparison against a null matched on Ca-Ca *distance* gives 1.02-fold
(bootstrap 95% CI 0.94-1.11), in every one of five replicates. The enrichment
was spatial proximity, and the sequence-separation filter did not control for
it.

The lesson generalises past that one protein, so in EVODYN a matched null is
not an option a user can forget to switch on. :func:`enrichment_panel` reports
every available null side by side, and the exporter refuses to write an
enrichment statistic that has no matched comparison.

Matching is stratified resampling: the test set is binned on the covariate,
and the null draw takes the same number of pairs from each bin. Confidence
intervals are bootstrap percentiles over both the test set and the pooled null.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Literal, Sequence

import numpy as np
import pandas as pd


@dataclass
class EnrichmentResult:
    """One enrichment statistic and the null it was measured against."""

    match: str                 # "none" | the covariate matched on
    ratio: float               # mean(signal | top pairs) / mean(signal | null)
    ci_low: float
    ci_high: float
    n_test: int
    n_null_draws: int
    test_mean: float
    null_mean: float
    note: str = ""

    @property
    def spans_one(self) -> bool:
        """True when the interval includes no enrichment."""
        return self.ci_low <= 1.0 <= self.ci_high

    def to_dict(self) -> dict:
        d = asdict(self)
        d["ci_spans_one"] = self.spans_one
        return d


def matched_enrichment(
    pairs: pd.DataFrame,
    signal: str,
    score: str = "score",
    match_on: str | None = None,
    top_n: int = 500,
    bin_width: float | None = 2.0,
    bins: Sequence[float] | None = None,
    categorical: bool = False,
    n_draws: int = 200,
    n_bootstrap: int = 2000,
    seed: int = 0,
    absolute: bool = True,
) -> EnrichmentResult:
    """Enrichment of ``signal`` among top-``top_n`` pairs by ``score``.

    Parameters
    ----------
    pairs
        One row per residue pair. Must already be filtered to the comparison
        set (for example, pairs above a sequence-separation cutoff).
    signal
        Column whose mean is compared, e.g. the dynamic cross-correlation.
    match_on
        Column the null is matched on. ``None`` gives the unmatched null -
        a plain resample of the remaining pairs - which is reported for
        contrast, never alone.
    categorical
        Match on exact values of ``match_on`` rather than numeric bins. Used
        for secondary-structure class.
    absolute
        Compare ``|signal|``. Correlated and anti-correlated motion are both
        dynamic coupling; sign-preserving means cancel.

    Returns
    -------
    EnrichmentResult
    """
    if signal not in pairs or score not in pairs:
        raise KeyError(f"pairs needs columns {score!r} and {signal!r}")
    df = pairs.dropna(subset=[score, signal]).copy()
    if match_on:
        df = df.dropna(subset=[match_on])
    # Deterministic row order: the draw must not depend on how the table was
    # read. The score is included as a tiebreak so that the order is total even
    # if a caller passes a table with repeated pairs.
    sort_cols = [c for c in ("res_i", "res_j") if c in df.columns] + [score]
    df = df.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)

    vals = df[signal].abs() if absolute else df[signal]
    df = df.assign(_sig=vals.to_numpy())

    top = df.nlargest(top_n, score)
    rest = df.drop(top.index)
    if rest.empty:
        raise ValueError("no pairs left for the null after taking the top set")
    test = top["_sig"].to_numpy()

    rng = np.random.default_rng(seed)
    draws: list[np.ndarray] = []

    if match_on is None:
        for _ in range(n_draws):
            draws.append(rng.choice(rest["_sig"].to_numpy(),
                                    size=min(top_n, len(rest)), replace=False))
        note = "unmatched: pairs resampled from the remainder"
    else:
        if categorical:
            tb = top[match_on].astype(str).to_numpy()
            rb = rest[match_on].astype(str).to_numpy()
        else:
            edges = (np.asarray(bins) if bins is not None
                     else np.arange(0.0,
                                    float(np.nanmax(df[match_on])) + bin_width,
                                    bin_width))
            tb = np.digitize(top[match_on].to_numpy(), edges)
            rb = np.digitize(rest[match_on].to_numpy(), edges)
        counts = pd.Series(tb).value_counts()
        pools = {b: np.flatnonzero(rb == b) for b in counts.index}
        usable = {b: k for b, k in counts.items() if len(pools.get(b, ())) > 0}
        covered = sum(usable.values()) / len(tb)
        rest_sig = rest["_sig"].to_numpy()
        for _ in range(n_draws):
            picks = [rng.choice(pools[b], size=k, replace=len(pools[b]) < k)
                     for b, k in usable.items()]
            draws.append(rest_sig[np.concatenate(picks)])
        note = (f"matched on {match_on}"
                + (" (exact)" if categorical else f" (bin width {bin_width})")
                + f"; {100 * covered:.1f}% of test pairs had a non-empty donor stratum")

    pooled = np.concatenate(draws)
    ratio = float(test.mean() / pooled.mean())
    boot = np.empty(n_bootstrap)
    for b in range(n_bootstrap):
        boot[b] = (rng.choice(test, test.size, replace=True).mean()
                   / rng.choice(pooled, test.size, replace=True).mean())
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return EnrichmentResult(
        match=match_on or "none", ratio=ratio, ci_low=float(lo), ci_high=float(hi),
        n_test=int(test.size), n_null_draws=n_draws,
        test_mean=float(test.mean()), null_mean=float(pooled.mean()), note=note,
    )


def enrichment_panel(
    pairs: pd.DataFrame,
    signal: str,
    score: str = "score",
    sequence_separation: int = 5,
    distance_column: str | None = "ca_distance",
    conservation_column: str | None = "conservation_min",
    ss_column: str | None = None,
    top_n: int = 500,
    bin_width: float = 2.0,
    n_draws: int = 200,
    n_bootstrap: int = 2000,
    seed: int = 0,
) -> pd.DataFrame:
    """Every available null for one signal, side by side.

    The unmatched null is always included, because the contrast between it and
    the matched nulls is the result: a large unmatched enrichment that
    collapses under matching says the covariate explains the signal.
    """
    df = pairs
    if {"res_i", "res_j"} <= set(df.columns) and sequence_separation:
        df = df[(df.res_i - df.res_j).abs() > sequence_separation]

    rows = [matched_enrichment(df, signal, score, None, top_n,
                               n_draws=n_draws, n_bootstrap=n_bootstrap, seed=seed)]
    if distance_column and distance_column in df:
        rows.append(matched_enrichment(df, signal, score, distance_column, top_n,
                                       bin_width=bin_width, n_draws=n_draws,
                                       n_bootstrap=n_bootstrap, seed=seed))
    if conservation_column and conservation_column in df:
        q = df[conservation_column].quantile(np.linspace(0, 1, 11)).to_numpy()
        rows.append(matched_enrichment(df, signal, score, conservation_column, top_n,
                                       bins=np.unique(q), n_draws=n_draws,
                                       n_bootstrap=n_bootstrap, seed=seed))
    if ss_column and ss_column in df:
        rows.append(matched_enrichment(df, signal, score, ss_column, top_n,
                                       categorical=True, n_draws=n_draws,
                                       n_bootstrap=n_bootstrap, seed=seed))
    return pd.DataFrame([r.to_dict() for r in rows])
