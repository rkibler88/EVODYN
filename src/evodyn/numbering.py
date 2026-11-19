"""Canonical residue numbering, and the assertions that keep it honest.

Three coordinate systems are in play in any analysis of this kind:

1. **alignment columns** - the family alignment, including insert states;
2. **canonical residue numbers** - the reference protein's own numbering, which
   is what a reader, a mutation and a structure paper all use;
3. **structure residue numbers** - whatever the PDB or model file happens to
   carry, which is frequently 1-based regardless of the construct.

Everything EVODYN writes is in canonical numbering. The translation happens
once, here, and is checked by comparing residue identities between the
alignment reference row and the structure. A silent off-by-one between these
systems is the single most damaging error this kind of pipeline can make: in
the EGFR study it put the DFG motif at the wrong residue in figures and text
(853 is isoleucine; the motif is D855-F856-G857), and that survived until the
identities were explicitly asserted.

Two lengths are also routinely confused and are therefore kept separate and
named:

``n_focus``
    residues spanned by the requested region, insert states included;
``n_scored``
    match-state positions actually scored by the coevolution stage, which is
    the number of nodes the network has.
"""
from __future__ import annotations

import pathlib
from dataclasses import dataclass, field

from evodyn.exceptions import MappingError

#: three-letter to one-letter, including the residues a prepared model may carry
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K", "SEC": "U",
}

GAP_CHARS = set("-.")


@dataclass
class CanonicalMap:
    """Alignment columns, canonical residues and structure residues.

    Parameters
    ----------
    col2res
        1-based alignment column -> canonical residue number. Match states only.
    res2aa
        canonical residue number -> one-letter amino acid, from the alignment
        reference row.
    region
        ``(first, last)`` canonical residue numbers spanned by the reference row
        including insert states - the *focus* region, which is what a
        coevolution run was asked to score.
    analysis_region
        the sub-range actually analysed, when :meth:`subset` has been applied.
        Defaults to ``region``.

    The two are kept apart on purpose. In the EGFR study the focus region was
    712-968 (257 residues) while the analysed domain was 714-966 (253 scored
    positions), because four columns were insert states; an L-proportional edge
    count computed from the wrong one of those numbers silently changes the
    network.
    """

    col2res: dict[int, int]
    res2aa: dict[int, str]
    region: tuple[int, int]
    reference_id: str | None = None
    structure_offset: int | None = field(default=None)
    analysis_region: tuple[int, int] | None = None

    def __post_init__(self):
        if self.analysis_region is None:
            self.analysis_region = self.region

    # ------------------------------------------------------------- builders
    @classmethod
    def from_reference_row(
        cls,
        fasta: str | pathlib.Path,
        start: int,
        reference_id: str | None = None,
    ) -> "CanonicalMap":
        """Build from an aligned reference row.

        The row follows the Pfam/HMMER convention: uppercase letters are match
        states, lowercase letters are insert states, ``-`` and ``.`` are gaps.
        Insert states consume a residue number but are not scored, which is
        exactly why ``n_focus`` and ``n_scored`` differ.

        Parameters
        ----------
        fasta
            FASTA file whose first record is the aligned reference row.
        start
            Canonical number of the first residue in that row.
        """
        header, seq = _read_first_record(fasta)
        col2res: dict[int, int] = {}
        res2aa: dict[int, str] = {}
        cur = start
        for i, ch in enumerate(seq, start=1):
            if ch in GAP_CHARS:
                continue
            if ch.isupper():
                col2res[i] = cur
                res2aa[cur] = ch
            cur += 1
        if not col2res:
            raise MappingError(f"{fasta}: reference row has no match-state columns")
        return cls(col2res, res2aa, (start, cur - 1), reference_id or header)

    @classmethod
    def from_sequence(cls, sequence: str, start: int, reference_id: str | None = None):
        """Build from an ungapped sequence; every position is a match state."""
        col2res = {i: start + i - 1 for i in range(1, len(sequence) + 1)}
        res2aa = {start + i: ch.upper() for i, ch in enumerate(sequence)}
        return cls(col2res, res2aa, (start, start + len(sequence) - 1), reference_id)

    # ------------------------------------------------------------ accessors
    @property
    def residues(self) -> list[int]:
        """Canonical residue numbers of the scored (match-state) positions."""
        return sorted(self.col2res.values())

    @property
    def n_scored(self) -> int:
        """Match-state positions - the number of nodes the network will have."""
        return len(self.col2res)

    @property
    def n_focus(self) -> int:
        """Residues spanned by the reference row, insert states included.

        This is the ``L`` of an L-proportional edge count, not the number of
        network nodes.
        """
        return self.region[1] - self.region[0] + 1

    def residue(self, column: int) -> int | None:
        return self.col2res.get(column)

    def subset(self, first: int, last: int) -> "CanonicalMap":
        """Restrict to a canonical residue range."""
        col2res = {c: r for c, r in self.col2res.items() if first <= r <= last}
        if not col2res:
            raise MappingError(f"no match-state positions in range {first}-{last}")
        res2aa = {r: a for r, a in self.res2aa.items() if first <= r <= last}
        return CanonicalMap(col2res, res2aa, self.region, self.reference_id,
                            self.structure_offset, (first, last))

    def sequence(self) -> str:
        return "".join(self.res2aa[r] for r in self.residues)

    # ----------------------------------------------------------- assertions
    def check_structure(
        self,
        pdb: str | pathlib.Path,
        chain: str | None = None,
        max_mismatch_fraction: float = 0.02,
        allow_offset: bool = True,
    ) -> int:
        """Assert that a structure agrees with the alignment reference row.

        Residue identities are compared position by position. If the structure
        is numbered on a different origin (a model of a construct numbered from
        1 is the common case) the offset is detected and returned; it is then
        applied to every coordinate EVODYN reads from that file.

        Returns
        -------
        int
            The offset added to structure numbering to reach canonical
            numbering. Zero when the file is already canonical.

        Raises
        ------
        MappingError
            If no offset brings identity agreement within
            ``max_mismatch_fraction``. The run stops; it is not annotated and
            continued.
        """
        struct = read_structure_residues(pdb, chain=chain)
        if not struct:
            raise MappingError(f"{pdb}: no standard residues found"
                               + (f" in chain {chain}" if chain else ""))

        offsets = [0]
        if allow_offset:
            lo = min(self.residues) - max(struct)
            hi = max(self.residues) - min(struct)
            # candidate offsets that at least align the two ranges
            offsets = sorted({0, min(self.residues) - min(struct)}
                             | set(range(lo, hi + 1)) if (hi - lo) <= 4000 else {0, min(self.residues) - min(struct)})

        best, best_score, best_overlap = None, -1.0, 0
        for off in offsets:
            matched = total = 0
            for rnum, aa in struct.items():
                canon = rnum + off
                ref = self.res2aa.get(canon)
                if ref is None:
                    continue
                total += 1
                matched += (ref == aa)
            if total < max(10, 0.2 * self.n_scored):
                continue
            score = matched / total
            if score > best_score:
                best, best_score, best_overlap = off, score, total

        if best is None:
            raise MappingError(
                f"{pdb}: no numbering offset gives meaningful overlap with the "
                f"reference row ({self.n_scored} scored positions, "
                f"structure spans {min(struct)}-{max(struct)})"
            )
        if 1.0 - best_score > max_mismatch_fraction:
            raise MappingError(
                f"{pdb}: best offset {best:+d} still disagrees with the alignment "
                f"reference at {100 * (1 - best_score):.1f}% of {best_overlap} "
                f"compared residues (tolerance {100 * max_mismatch_fraction:.1f}%). "
                "Check that the structure, the reference row and target.region "
                "describe the same construct."
            )
        self.structure_offset = best
        return best

    def assert_motif(self, first_residue: int, motif: str) -> None:
        """Assert that a named motif sits where the run claims it does.

        Intended for the handful of landmarks a manuscript names explicitly
        (DFG, HRD, a gatekeeper, a catalytic lysine). Cheap, and it would have
        caught the EGFR DFG mislabelling before it reached a figure.
        """
        got = "".join(self.res2aa.get(first_residue + i, "?") for i in range(len(motif)))
        if got != motif.upper():
            raise MappingError(
                f"motif assertion failed: expected {motif.upper()} at "
                f"{first_residue}-{first_residue + len(motif) - 1}, found {got}"
            )


# --------------------------------------------------------------------------
def _read_first_record(fasta: str | pathlib.Path) -> tuple[str, str]:
    header, seq = None, []
    with open(fasta) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:].strip()
            elif header is not None:
                seq.append(line.strip())
    if header is None:
        raise MappingError(f"{fasta}: not a FASTA file")
    return header, "".join(seq)


def read_structure_residues(
    pdb: str | pathlib.Path, chain: str | None = None
) -> dict[int, str]:
    """Residue number -> one-letter code, from a PDB file.

    Insertion codes are ignored (a residue with one is skipped rather than
    silently merged with its neighbour) and only the first model is read.
    """
    out: dict[int, str] = {}
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            if chain and line[21] != chain:
                continue
            if line[26] != " ":
                continue
            aa = THREE_TO_ONE.get(line[17:20].strip().upper())
            if aa is None:
                continue
            out[int(line[22:26])] = aa
    return out


def read_ca_coordinates(
    pdb: str | pathlib.Path, chain: str | None = None, offset: int = 0
) -> tuple[list[int], "object"]:
    """CA coordinates in canonical numbering.

    Returns ``(residues, xyz)`` with ``xyz`` an ``(N, 3)`` array ordered to
    match ``residues``.
    """
    import numpy as np

    res, xyz = [], []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[12:16].strip() != "CA" or line[26] != " ":
                continue
            if chain and line[21] != chain:
                continue
            if line[17:20].strip().upper() not in THREE_TO_ONE:
                continue
            res.append(int(line[22:26]) + offset)
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return res, np.asarray(xyz, dtype=float)
