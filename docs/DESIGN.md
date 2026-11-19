# Design

EVODYN exists because one protein's analysis was worth generalising, and
because generalising it is risky in a specific way: a pipeline that runs over
many proteins propagates a methodological error at the same rate it produces
results. Every design decision below is a response to something that went wrong
once, in a case where a human was checking each step by hand.

## 1. The controls are stages

The EGFR study reported that top-coupled distal residue pairs were enriched
1.83-fold in correlated motion. The null was matched on sequence separation,
which is the field-standard filter. Matched on Cα–Cα distance instead, the same
comparison gives 1.02-fold with a bootstrap 95% interval of 0.94–1.11, in every
one of five replicates. The top-coupled pairs sit at a median 11.4 Å against
24.3 Å for the rest: the effect was spatial proximity.

`evodyn.controls.nulls.enrichment_panel` therefore reports every available null
together, and the unmatched null is included precisely so that the contrast is
visible. When a matched interval spans 1.0 and the unmatched ratio does not,
the pipeline writes a sentence saying which covariate accounts for the signal
into the run's limits section. It is harder to publish the unmatched number
than to publish the matched one.

Three matched nulls are implemented. All are stratified resampling: bin the
test set on the covariate, draw the same number of null pairs from each bin.

| null | matched on | controls for |
|---|---|---|
| distance | Cα–Cα distance, 2 Å bins | the network being mostly a contact map |
| conservation | relative entropy, deciles | conserved positions clustering in the core |
| secondary structure | DSSP class pair, exact | β-pairing and helix geometry |

Conservation is computed by EVODYN from the same alignment the couplings came
from. In the original study the conservation vector came from pySCA and was
never committed, so that control could not be recomputed from a clone.

## 2. Numbering is asserted

Three coordinate systems, all in play at once:

* alignment columns, including insert states;
* canonical residue numbers, which is what a reader and a mutation use;
* structure numbering, frequently 1-based regardless of the construct.

`CanonicalMap` holds the translation and checks it. `check_structure` compares
residue identities between the alignment reference row and the structure across
candidate offsets, returns the offset that agrees, and raises `MappingError` if
none agrees within tolerance. `assert_motif` checks a named landmark. These are
fatal, not advisory: in the EGFR study residue 853 was labelled as the DFG motif
in figures and text, and 853 is isoleucine — the motif is D855-F856-G857.

Two lengths are also kept separate and named, because confusing them silently
changes the network: `n_focus` is the residues spanned by the reference row
including insert states (257 for EGFR), `n_scored` is the match-state positions
actually scored (253). An L-proportional edge count uses the first; the node
count is the second.

## 3. Every run states its limits

`report.md` ends with a limits section that the stages write as they run. The
EGFR run produces, unprompted:

* 55.5% of coupling scores are negative, so the score must not be summed over
  arbitrary residue-pair sets as though it were a magnitude;
* 92% of network edges are spatial contacts with a median length of 3.9 Å, so
  the network is close to a contact map;
* the 1.83-fold enrichment falls to 1.02 under distance matching, so distance
  accounts for the signal.

A run with `dynamics.mode: equilibrate` adds what the simulations did and did
not sample. A run with replicates adds that replicates sharing a force field
and a starting structure measure precision, not accuracy. A run whose alignment
is family-wide adds its taxonomic composition.

## 4. Thresholds are derived where they can be

The only free parameter in the network is how many top-scoring pairs to keep as
edges. `threshold: percolation` locates the smallest edge count at which the
graph contains every scored position and is connected — a property of the
coupling distribution. `threshold: multiplier` reproduces the L-proportional
convention for continuity with published partitions. Louvain is stochastic, so
it is run over `louvain_seeds` seeds, the best weighted modularity is kept, and
the winning seed is written into the manifest.

## 5. Failure modes are distinguished

Three kinds of thing can go wrong, and they are not the same:

* **a stage cannot run** — a tool is missing, an input was not supplied. The
  stage is skipped, the run continues, and the limits section records what the
  output therefore cannot support.
* **a stage produces a wrong answer** — a numbering mismatch, a degenerate
  score column. Fatal. The manifest is written and the run stops.
* **a configuration is wrong** — rejected before anything runs, by name.

## 6. Reusing rather than reimplementing

The MD stage delegates to EquilibraTor (Cediel-Becerra, Silva and Dias, CSBJ
2025), a published and maintained GROMACS driver, rather than shipping another
set of `.mdp` files. The
coevolution stage calls plmc. Sector decomposition calls pySCA. EVODYN's
contribution is the numbering spine, the controls, the schema and the
provenance — not another implementation of DCA.
