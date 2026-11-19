"""Redundancy filtering, reference rescue and alignment composition."""
from __future__ import annotations

import collections
import pathlib
import shutil
import subprocess

import pandas as pd

from evodyn.exceptions import ExternalToolError


def cdhit_filter(
    fasta, out_fasta, identity: float = 0.62, word_size: int = 4,
    memory_mb: int = 16000, threads: int = 8, binary: str = "cd-hit",
) -> pathlib.Path:
    """Cluster an alignment at a sequence-identity threshold.

    The word size must be compatible with the identity threshold or CD-HIT
    refuses to run; the usual table is n=5 for 0.7-1.0, n=4 for 0.6-0.7,
    n=3 for 0.5-0.6, n=2 for 0.4-0.5. EVODYN checks this rather than letting
    the user discover it from CD-HIT's error, and it records the command that
    actually ran - a mismatch between a script's comment and its command
    (90% / n=5 documented, 62% / n=4 executed) was one of the discrepancies
    the EGFR provenance audit turned up.
    """
    expected = 5 if identity >= 0.7 else 4 if identity >= 0.6 else 3 if identity >= 0.5 else 2
    if word_size != expected:
        raise ValueError(
            f"cd-hit word size {word_size} is not valid at identity {identity}; "
            f"use {expected}")
    if shutil.which(binary) is None:
        raise ExternalToolError(
            f"{binary} not found on PATH (conda install -c bioconda cd-hit)")
    cmd = [binary, "-i", str(fasta), "-o", str(out_fasta), "-c", str(identity),
           "-n", str(word_size), "-M", str(memory_mb), "-T", str(threads)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise ExternalToolError(f"cd-hit failed: {proc.stderr[-2000:]}")
    return pathlib.Path(out_fasta)


def rescue_reference(aligned_fasta, reference_id: str, source_fasta, out_fasta) -> bool:
    """Re-insert the reference sequence if clustering dropped it.

    The reference protein is frequently not chosen as a cluster
    representative. If it is missing, every downstream coordinate mapping is
    against a sequence that is not in the alignment. EVODYN restores it and
    records that it did; the EGFR alignment needed exactly this (32,911
    representatives plus one rescued row).
    """
    have = {l[1:].split()[0] for l in open(aligned_fasta) if l.startswith(">")}
    if reference_id in have:
        if pathlib.Path(aligned_fasta) != pathlib.Path(out_fasta):
            shutil.copy(aligned_fasta, out_fasta)
        return False
    row = None
    with open(source_fasta) as fh:
        keep, buf = False, []
        for line in fh:
            if line.startswith(">"):
                if keep:
                    break
                keep = line[1:].split()[0] == reference_id
                if keep:
                    buf.append(line)
            elif keep:
                buf.append(line)
        row = "".join(buf) if buf else None
    if row is None:
        raise KeyError(f"reference {reference_id!r} not found in {source_fasta}")
    shutil.copy(aligned_fasta, out_fasta)
    with open(out_fasta, "a") as fh:
        fh.write(row if row.endswith("\n") else row + "\n")
    return True


def alignment_composition(fasta, rank_from_header=None, top: int = 10) -> pd.DataFrame:
    """Who is actually in the alignment.

    Default grouping is by the species suffix of a UniProt-style identifier
    (``EGFR_HUMAN/712-968`` -> ``HUMAN``). Supply ``rank_from_header`` to group
    by taxonomy instead, if a mapping is available.
    """
    counter = collections.Counter()
    n = 0
    for line in open(fasta):
        if not line.startswith(">"):
            continue
        n += 1
        head = line[1:].split()[0]
        if rank_from_header:
            key = rank_from_header(head)
        else:
            name = head.split("/")[0]
            key = name.split("_")[-1] if "_" in name else "unknown"
        counter[key] += 1
    df = pd.DataFrame(counter.most_common(top), columns=["group", "sequences"])
    df["percent"] = 100.0 * df.sequences / max(n, 1)
    df.attrs["total_sequences"] = n
    return df
