# Roadmap

v0.1.0 is an alpha. This file is the honest state of it.

## Implemented and tested

| area | state |
|---|---|
| configuration, validation, config-relative paths | tested |
| canonical numbering, offset detection, motif assertions | tested |
| network construction, percolation threshold, Louvain seed sweep | tested |
| centrality, knockout, backbone-excluded shortest paths | tested |
| matched nulls (distance, conservation, categorical) with bootstrap CIs | tested |
| covariance → DCCM conversion | tested |
| per-residue and per-pair tables, run report, manifest | tested |
| EGFR regression: 16/16 published values reproduced | tested |

## Implemented, exercised only in dry run or on the example

| area | what is untested |
|---|---|
| InterPro alignment retrieval | not run against the live endpoint in CI |
| CD-HIT and reference rescue | no end-to-end run in the test suite |
| plmc invocation | argument assembly only; no binary in CI |
| EquilibraTor adapter | command assembly and directory layout only; no GROMACS in CI |
| trajectory → DCCM / RMSF | needs an mdtraj-readable fixture |
| pySCA sectors | written against the pySCA 6.1 API, not executed |
| Spheronizator and InteracTor adapters | shape of the call, not the call |

## Next

1. **A second protein, end to end.** The claim that this is a dataset generator
   is not demonstrated until a target other than EGFR has run from Pfam
   accession to tables. Until then it is a generalised pipeline with one
   validated case.
2. **A small MD fixture.** A few hundred frames of a small system, committed,
   so that the trajectory path is covered in CI rather than asserted.
3. **A multi-protein driver.** `evodyn batch` over a target list, with the
   per-target tables concatenated and the per-target limits preserved.
4. **Secondary-structure-matched nulls wired into the panel.** The sampler
   supports exact categorical matching; the pipeline does not yet populate the
   pair-level SS class column.
5. **Subfamily alignments.** The EGFR analysis left open whether an alignment
   restricted to receptor tyrosine kinases, rather than the whole kinase
   superfamily, changes the couplings. EVODYN should make that a configuration
   change rather than a new project.
6. **Replicate handling in the equilibrate path.** Replicates currently differ
   only by working directory; they should differ by initial velocity seed, and
   the seed should be in the manifest.

## Known limitations

* `check_structure` searches candidate offsets by brute force and is only
  practical for single-domain constructs.
* The conservation-matched null uses deciles of relative entropy with a fixed
  background composition; a family-specific background would be better.
* Secondary structure needs mdtraj and a DSSP executable; it is skipped with a
  recorded limit when either is absent.
* Nothing yet handles multi-chain targets or inter-chain couplings.
