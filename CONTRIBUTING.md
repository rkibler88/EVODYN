# Contributing

EVODYN is at v0.1.0 and the parts that need work are listed in
[docs/ROADMAP.md](docs/ROADMAP.md). The most useful contribution right now is a
second protein taken end to end, from Pfam accession to tables.

## Ground rules

1. **A new analysis arrives with its null.** A statistic that compares a
   selected set against a background is not merged without a matched null and
   an interval. This is the project's one non-negotiable rule.
2. **Numbering assertions are not softened.** If `check_structure` or
   `assert_motif` fails, the input is wrong. Do not widen the tolerance to make
   a run pass.
3. **A stage that cannot run is skipped with a recorded reason.** A stage that
   would produce a wrong answer raises. Keep the distinction.
4. **New outputs go into the existing schema.** A column added to
   `residues.csv` or `pairs.csv` must be meaningful for every protein, not just
   the one that motivated it.

## Development

```bash
conda env create -f environment.yml
conda activate evodyn
pip install -e ".[dev]"
make test          # fast unit tests
make all-tests     # with the EGFR regression
```

The EGFR regression is the acceptance test: 16 values from the published case
study, recomputed through the generalised code. If a change moves one of them,
either the change is wrong or the change is a correction — and a correction
needs a note in `docs/ROADMAP.md` saying which value moved and why.
