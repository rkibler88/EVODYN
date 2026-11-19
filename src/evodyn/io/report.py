"""The run report, and its limits section.

Every EVODYN run writes what its own output does not support. This is not
decoration: the EGFR study's headline co-evolution/dynamics result survived
review only as a spatial effect, and the difference between the published
claim and the defensible one was entirely a matter of which null was used and
whether the scope of the simulation was stated. A tool that generates these
datasets at scale will propagate that error at scale unless the limits travel
with the numbers.
"""
from __future__ import annotations

import pathlib
import textwrap

import pandas as pd


def _markdown_table(df: pd.DataFrame, floatfmt: str = ".3f") -> str:
    """Markdown table without pulling in an optional dependency."""
    def cell(v):
        return format(v, floatfmt) if isinstance(v, float) else str(v)

    head = "| " + " | ".join(df.columns) + " |"
    rule = "|" + "|".join("---" for _ in df.columns) + "|"
    rows = ["| " + " | ".join(cell(v) for v in r) + " |"
            for r in df.itertuples(index=False)]
    return "\n".join([head, rule, *rows])


def write_report(manifest, summary: dict, path, enrichment: pd.DataFrame | None = None) -> pathlib.Path:
    """Write ``report.md`` - what ran, what came out, what it does not support."""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cfg = manifest.config
    target = cfg.get("target", {})

    lines = [
        f"# EVODYN run report: {target.get('name', 'unnamed')}",
        "",
        f"EVODYN {summary.get('evodyn_version') or __import__('evodyn').__version__} | "
        f"{len(manifest.stages)} stages | "
        f"{sum(1 for s in manifest.stages if s.status == 'ok')} completed",
        "",
        "## Target",
        "",
        f"- name: `{target.get('name')}`",
        f"- UniProt: `{target.get('uniprot')}`",
        f"- family: `{target.get('pfam')}`",
        f"- canonical region: {target.get('region')}",
        "",
        "## Stages",
        "",
        "| stage | status | seconds | note |",
        "|---|---|---|---|",
    ]
    for s in manifest.stages:
        lines.append(f"| {s.name} | {s.status} | {s.seconds:.1f} | {s.message} |")

    if summary:
        lines += ["", "## Results", "", "| quantity | value |", "|---|---|"]
        for k, v in summary.items():
            if k == "evodyn_version":
                continue
            lines.append(f"| {k} | {v} |")

    if enrichment is not None and len(enrichment):
        lines += ["", "## Enrichment, with matched nulls", "",
                  "An enrichment measured against an unmatched null is reported "
                  "here only alongside its matched comparisons. Where the matched "
                  "interval spans 1.0, the covariate explains the signal.", "",
                  _markdown_table(enrichment)]

    lines += ["", "## Limits of this run", ""]
    if manifest.limits:
        lines += [f"- {textwrap.fill(l, 100, subsequent_indent='  ')}"
                  for l in manifest.limits]
    else:
        lines.append("- none recorded")

    lines += ["", "## Reproducing", "",
              "```bash", "evodyn run --config config.yaml", "```", "",
              "Package versions, external binary versions and SHA-256 of every "
              "input and output are in `manifest.json`.", ""]
    path.write_text("\n".join(lines))
    return path
