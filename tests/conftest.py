import pathlib
import pytest

DATA = pathlib.Path(__file__).resolve().parents[1] / "examples" / "egfr" / "data"


@pytest.fixture(scope="session")
def egfr_data():
    if not (DATA / "couplings_translated.csv.gz").exists():
        pytest.skip("EGFR benchmark data not present")
    return DATA


@pytest.fixture
def reference_row(tmp_path):
    """A small aligned reference row: match states, insert states and gaps."""
    p = tmp_path / "ref.fasta"
    #                  insert   match states            insert
    p.write_text(">TEST_PROT/10-21\nmm-ACDEFGHIK--nn\n")
    return p
