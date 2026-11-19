"""The null models, tested against data with a known answer."""
import numpy as np
import pandas as pd
import pytest

from evodyn.controls import matched_enrichment, enrichment_panel


@pytest.fixture
def confounded_pairs():
    """A signal driven entirely by distance, not by score.

    ``score`` and ``signal`` are both functions of ``ca_distance`` and are
    conditionally independent given it. An unmatched null must show
    enrichment; a distance-matched null must not.
    """
    rng = np.random.default_rng(0)
    n = 6000
    # distinct pairs, as a real pair table has
    ii, jj = np.triu_indices(130, k=6)
    pick = rng.choice(len(ii), size=n, replace=False)
    i, j = ii[pick] + 1, jj[pick] + 1
    dist = rng.uniform(4, 40, n)
    score = -dist + rng.normal(0, 3, n)
    signal = np.exp(-dist / 12) + rng.normal(0, 0.05, n)
    return pd.DataFrame({"res_i": i, "res_j": j, "score": score,
                         "ca_distance": dist, "signal": signal})


def test_unmatched_null_reports_the_confounded_effect(confounded_pairs):
    r = matched_enrichment(confounded_pairs, signal="signal", match_on=None,
                           top_n=500, n_bootstrap=400, seed=1)
    assert r.ratio > 1.5
    assert not r.spans_one


def test_distance_matching_removes_it(confounded_pairs):
    r = matched_enrichment(confounded_pairs, signal="signal", match_on="ca_distance",
                           top_n=500, bin_width=2.0, n_bootstrap=400, seed=1)
    assert 0.9 < r.ratio < 1.1
    assert r.spans_one


def test_panel_reports_every_null_together(confounded_pairs):
    panel = enrichment_panel(confounded_pairs, signal="signal",
                             conservation_column=None, n_bootstrap=300, seed=1)
    assert list(panel["match"]) == ["none", "ca_distance"]
    assert panel.loc[panel.match == "none", "ratio"].iat[0] > \
           panel.loc[panel.match == "ca_distance", "ratio"].iat[0]


def test_draw_is_independent_of_row_order(confounded_pairs):
    """The original analysis's null depended on the order the table was read."""
    a = matched_enrichment(confounded_pairs, signal="signal", match_on="ca_distance",
                           top_n=300, n_bootstrap=200, seed=7)
    shuffled = confounded_pairs.sample(frac=1.0, random_state=99)
    b = matched_enrichment(shuffled, signal="signal", match_on="ca_distance",
                           top_n=300, n_bootstrap=200, seed=7)
    assert a.ratio == pytest.approx(b.ratio, rel=1e-12)
