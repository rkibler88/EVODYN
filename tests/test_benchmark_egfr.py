"""The EGFR regression: does the generalisation reproduce the case study?"""
import pytest

from evodyn.benchmark import PUBLISHED, run_benchmark


@pytest.mark.slow
def test_egfr_case_study_reproduces(egfr_data):
    assert run_benchmark(data_dir=egfr_data, verbose=True)


@pytest.mark.slow
def test_enrichment_collapses_under_distance_matching(egfr_data):
    """The result that motivated the controls being built in."""
    from evodyn.benchmark import CONFIG
    from evodyn.config import load_config
    from evodyn.pipeline import Pipeline

    cfg = load_config(CONFIG)
    enr = Pipeline(cfg).run().enrichment.set_index("match")
    assert enr.loc["none", "ratio"] > 1.5
    assert not enr.loc["none", "ci_spans_one"]
    assert enr.loc["ca_distance", "ci_spans_one"]
