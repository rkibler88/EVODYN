"""Family alignment retrieval.

The family alignment is the input that most strongly determines what a run can
conclude, and it is the one users are least likely to interrogate. PF07714
spans tyrosine *and* serine/threonine kinases; after redundancy filtering the
EGFR study's retained set was 34.8% Viridiplantae, 26.2% Metazoa, 24.0% Fungi
and contained exactly one ErbB-family sequence. Claims of receptor-tyrosine-
kinase specificity had to be withdrawn. :func:`evodyn.sequences.filter.
alignment_composition` runs by default so that this is visible in the run
report instead of being discovered at review.
"""
from __future__ import annotations

import gzip
import pathlib
import urllib.request

INTERPRO_ALIGNMENT = (
    "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{accession}/"
    "?annotation=alignment:{kind}&download"
)


def fetch_pfam_alignment(
    accession: str, outdir, kind: str = "full", timeout: float = 600.0
) -> pathlib.Path:
    """Download a Pfam alignment from InterPro, in Stockholm format.

    ``kind`` is ``full``, ``seed`` or ``uniprot``. Full alignments of large
    families are large - PF07714 is 2.5 GB uncompressed, 400,321 sequences by
    3,839 columns - so the gzip is kept on disk and decompressed on demand.
    """
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dest = outdir / f"{accession}.{kind}.sto.gz"
    if dest.exists():
        return dest
    url = INTERPRO_ALIGNMENT.format(accession=accession, kind=kind)
    with urllib.request.urlopen(url, timeout=timeout) as resp, open(dest, "wb") as fh:
        while chunk := resp.read(1 << 20):
            fh.write(chunk)
    return dest


def _open(path):
    path = pathlib.Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def stockholm_to_fasta(stockholm, fasta, keep_gaps: bool = True) -> int:
    """Stream a Stockholm alignment to FASTA, preserving alignment columns.

    Streaming rather than parsing into memory is the point: a Biopython
    ``AlignIO.read`` of a full Pfam alignment will not fit on a laptop, and
    the conversion is the first stage every run performs.
    """
    rows: dict[str, list[str]] = {}
    order: list[str] = []
    with _open(stockholm) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "//")):
                continue
            name, _, seq = line.partition(" ")
            seq = seq.strip()
            if not seq:
                continue
            if name not in rows:
                rows[name] = []
                order.append(name)
            rows[name].append(seq if keep_gaps else seq.replace(".", "").replace("-", ""))
    with open(fasta, "w") as out:
        for name in order:
            out.write(f">{name}\n{''.join(rows[name])}\n")
    return len(order)
