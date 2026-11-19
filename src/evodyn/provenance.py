"""Run manifest.

Every run writes ``manifest.json``: the configuration as resolved, the version
of every package and external binary that produced a number, the SHA-256 of
every input and output file, per-stage timings and status, and the accumulated
limits. A result without its manifest is not reproducible, and the EGFR study
spent a revision cycle reconstructing exactly this information by hand.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import platform
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from typing import Any


def sha256(path: str | pathlib.Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def binary_version(binary: str) -> str | None:
    """Best-effort version string for an external executable."""
    where = shutil.which(binary)
    if where is None:
        return None
    for flag in ("--version", "-version", "-h"):
        try:
            p = subprocess.run([binary, flag], capture_output=True, text=True, timeout=30)
        except (OSError, subprocess.SubprocessError):
            continue
        text = (p.stdout or "") + (p.stderr or "")
        if text.strip():
            return text.strip().splitlines()[0][:200]
    return where


def python_packages() -> dict[str, str]:
    import importlib.metadata as md

    out = {}
    for name in ("numpy", "pandas", "scipy", "networkx", "python-louvain",
                 "biopython", "mdtraj", "evodyn"):
        try:
            out[name] = md.version(name)
        except md.PackageNotFoundError:
            continue
    return out


@dataclass
class StageRecord:
    name: str
    status: str = "pending"          # pending | ok | skipped | failed
    seconds: float = 0.0
    message: str = ""
    outputs: list[str] = field(default_factory=list)


@dataclass
class Manifest:
    """Accumulates everything needed to audit a run."""

    config: dict[str, Any]
    outdir: pathlib.Path
    started: float = field(default_factory=time.time)
    stages: list[StageRecord] = field(default_factory=list)
    limits: list[str] = field(default_factory=list)
    files: dict[str, str] = field(default_factory=dict)
    binaries: dict[str, str | None] = field(default_factory=dict)

    def stage(self, name: str) -> StageRecord:
        rec = StageRecord(name)
        self.stages.append(rec)
        return rec

    def limit(self, text: str) -> None:
        """Record a bound on what this run's output supports."""
        if text not in self.limits:
            self.limits.append(text)

    def record_file(self, path: str | pathlib.Path) -> None:
        p = pathlib.Path(path)
        if p.is_file():
            self.files[str(p)] = sha256(p)

    def record_binary(self, binary: str) -> None:
        self.binaries[binary] = binary_version(binary)

    def write(self) -> pathlib.Path:
        self.outdir.mkdir(parents=True, exist_ok=True)
        payload = {
            "evodyn_version": __import__("evodyn").__version__,
            "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(self.started)),
            "elapsed_seconds": round(time.time() - self.started, 2),
            "platform": {"python": sys.version.split()[0], "system": platform.platform()},
            "packages": python_packages(),
            "external_binaries": self.binaries,
            "config": self.config,
            "stages": [vars(s) for s in self.stages],
            "limits": self.limits,
            "file_sha256": self.files,
        }
        path = self.outdir / "manifest.json"
        path.write_text(json.dumps(payload, indent=2, sort_keys=False))
        return path
