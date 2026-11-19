"""Command-line interface.

::

    evodyn run      --config config.yaml    # the pipeline
    evodyn plan     --config config.yaml    # what would run, including the MD command
    evodyn validate --config config.yaml    # configuration only
    evodyn benchmark                        # the EGFR regression test

``plan`` exists because the dynamics stage can cost GPU-days. It resolves the
configuration, reports which stages would run, which external tools are
present, and prints the exact EquilibraTor invocation for each replicate -
without executing anything.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import shutil
import sys

from evodyn import __version__
from evodyn.config import load_config
from evodyn.exceptions import EvodynError


def _cmd_validate(args) -> int:
    cfg = load_config(args.config)
    print(f"configuration valid: target {cfg.target.name}")
    print(f"  stages      : {', '.join(cfg.stages)}")
    print(f"  dynamics    : {cfg.dynamics.mode}"
          + (f" x{cfg.dynamics.replicates}" if cfg.dynamics.mode == "equilibrate" else ""))
    print(f"  network     : {cfg.network.threshold}")
    print(f"  controls    : "
          + ", ".join(n for n, on in (("distance", cfg.controls.distance_matched),
                                      ("conservation", cfg.controls.conservation_matched),
                                      ("secondary structure", cfg.controls.ss_matched)) if on))
    return 0


def _cmd_plan(args) -> int:
    from evodyn.dynamics.equilibrator import build_command

    cfg = load_config(args.config)
    print(f"EVODYN {__version__} plan for {cfg.target.name}\n")
    print("stages:")
    for s in cfg.stages:
        print(f"  - {s}")

    print("\nexternal tools:")
    for binary, why in (("plmc", "direct coupling analysis"),
                        ("cd-hit", "redundancy filtering"),
                        ("EquilibraTor", "molecular dynamics"),
                        ("gmx", "GROMACS, used by EquilibraTor"),
                        ("mkdssp", "secondary structure")):
        where = shutil.which(binary)
        print(f"  {'found   ' if where else 'MISSING '} {binary:<14} {why}"
              + (f"  [{where}]" if where else ""))

    if cfg.dynamics.mode == "equilibrate":
        structure = cfg.structure.pdb or cfg.target.structure
        cmd = build_command(
            structure, binary=cfg.dynamics.equilibrator_binary,
            ligands=cfg.dynamics.ligands, force_field=cfg.dynamics.force_field,
            water_model=cfg.dynamics.water_model, box_type=cfg.dynamics.box_type,
            box_distance=cfg.dynamics.box_distance, pos_ion=cfg.dynamics.pos_ion,
            neg_ion=cfg.dynamics.neg_ion, extra_args=cfg.dynamics.extra_args)
        print(f"\nmolecular dynamics: {cfg.dynamics.replicates} replicate(s), "
              "each run in its own directory under "
              f"{cfg.out / 'md'}/replicate_NN/\n")
        print("  " + " ".join(cmd))
        print("\n  This is the expensive stage. Each replicate is a full GROMACS "
              "\n  setup-to-production run; budget accordingly.")
    else:
        print(f"\nmolecular dynamics: not run (dynamics.mode = {cfg.dynamics.mode})")
    return 0


def _cmd_run(args) -> int:
    from evodyn.pipeline import run_pipeline

    cfg = load_config(args.config)
    if args.outdir:
        cfg.outdir = args.outdir
    state = run_pipeline(cfg)
    print(f"\nEVODYN {__version__} - {cfg.target.name}")
    for rec in state.manifest.stages:
        mark = {"ok": "+", "skipped": "~", "failed": "!"}.get(rec.status, "?")
        print(f"  {mark} {rec.name:<12} {rec.seconds:>6.1f}s  {rec.message}")
    print(f"\nwritten to {cfg.out}/ (report.md, manifest.json)")
    if state.manifest.limits:
        print("\nlimits of this run:")
        for l in state.manifest.limits:
            print(f"  - {l}")
    return 0


def _cmd_benchmark(args) -> int:
    from evodyn.benchmark import run_benchmark

    ok = run_benchmark(data_dir=args.data, verbose=True)
    return 0 if ok else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="evodyn",
        description="Evolutionary, structural and dynamic datasets of protein allostery")
    p.add_argument("--version", action="version", version=f"evodyn {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    for name, fn, helptext in (
        ("run", _cmd_run, "run the pipeline"),
        ("plan", _cmd_plan, "show what would run, without running it"),
        ("validate", _cmd_validate, "check a configuration file"),
    ):
        sp = sub.add_parser(name, help=helptext)
        sp.add_argument("-c", "--config", required=True, help="run configuration (YAML)")
        if name == "run":
            sp.add_argument("-o", "--outdir", help="override outdir")
        sp.set_defaults(func=fn)

    sp = sub.add_parser("benchmark", help="regression test against the EGFR case study")
    sp.add_argument("--data", default=None, help="benchmark data directory")
    sp.set_defaults(func=_cmd_benchmark)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except EvodynError as exc:
        print(f"evodyn: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
