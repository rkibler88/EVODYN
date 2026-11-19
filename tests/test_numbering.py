"""Numbering is the part that must never be wrong silently."""
import pytest

from evodyn.exceptions import MappingError
from evodyn.numbering import CanonicalMap


def test_insert_states_consume_numbers_but_are_not_scored(reference_row):
    cmap = CanonicalMap.from_reference_row(reference_row, start=10)
    # mm = inserts 10,11; ACDEFGHIK = match states 12..20; nn = inserts 21,22
    assert cmap.n_scored == 9
    assert cmap.residues == list(range(12, 21))
    assert cmap.sequence() == "ACDEFGHIK"
    # the focus region spans the inserts too - this distinction drives the
    # L-proportional edge count
    assert cmap.n_focus == 13
    assert cmap.n_focus != cmap.n_scored


def test_subset_preserves_the_focus_region(reference_row):
    cmap = CanonicalMap.from_reference_row(reference_row, start=10)
    sub = cmap.subset(14, 18)
    assert sub.residues == [14, 15, 16, 17, 18]
    assert sub.n_focus == cmap.n_focus          # focus is a property of the row
    assert sub.analysis_region == (14, 18)


def test_motif_assertion_catches_an_off_by_two(reference_row):
    cmap = CanonicalMap.from_reference_row(reference_row, start=10)
    cmap.assert_motif(12, "ACD")
    with pytest.raises(MappingError):
        cmap.assert_motif(13, "ACD")


def test_structure_offset_is_detected(tmp_path):
    seq = "ACDEFGHIKL"
    cmap = CanonicalMap.from_sequence(seq, start=101)
    three = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
             "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU"}
    pdb = tmp_path / "s.pdb"
    pdb.write_text("".join(
        f"ATOM  {i:5d}  CA  {three[a]} A{i:4d}    "
        f"{i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n"
        for i, a in enumerate(seq, start=1)))
    assert cmap.check_structure(pdb) == 100


def test_structure_mismatch_is_fatal(tmp_path):
    cmap = CanonicalMap.from_sequence("ACDEFGHIKL", start=1)
    pdb = tmp_path / "s.pdb"
    pdb.write_text("".join(
        f"ATOM  {i:5d}  CA  TRP A{i:4d}    "
        f"{i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n"
        for i in range(1, 11)))
    with pytest.raises(MappingError):
        cmap.check_structure(pdb)
