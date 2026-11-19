# EGFR benchmark data

Committed intermediates from the case study EVODYN generalises,
[rkibler88/EGFRAllostery](https://github.com/rkibler88/EGFRAllostery)
(`data/processed/`). They are here so that `evodyn benchmark` and the EGFR
regression test run from a clone with no external tools and no large downloads.

| file | origin | what it is |
|---|---|---|
| `couplings_translated.csv.gz` | `egfr_couplings_translated.csv` | plmc couplings over PF07714, mapped to UniProt P00533 numbering; 31,878 pairs, APC-corrected score |
| `dccm_matrix.csv.gz` | `egfr_dccm_matrix.csv` | 253 × 253 dynamic cross-correlation, replicate 1 of five 100 ns trajectories |
| `rmsf.xvg` | `rmsf.xvg` | per-residue Cα RMSF, replicate 1 |
| `egfr_clean.pdb` | `egfr_clean.pdb` | AlphaFold2 model of the P00533 kinase domain; the source of Cα–Cα distances |
| `reference_aligned.fasta` | `egfr_reference_aligned.fasta` | the aligned `EGFR_HUMAN/712-968` row — uppercase match states, lowercase inserts. This is the column → residue map |
| `site_definitions.json` | `site_definitions.json` | residue sets for the ATP site, the allosteric pocket and the dimer interface |

Numbering is UniProt P00533 throughout. The reference row spans 712–968 (257
residues); four of those columns are insert states, so 253 positions are
scored and the network has 253 nodes.

Not shipped, and not needed by the benchmark: the 2.5 GB PF07714 full
alignment, the filtered 122 MB alignment, the plmc parameter file, and the
8.2 GB-per-replicate trajectories. A run that needs them fetches or regenerates
them — see `examples/new_protein/config.yaml`.

## What the benchmark checks

Sixteen values from the published analysis, recomputed through EVODYN: the
number of scored positions and focus residues, the coupling count, the network
at 462 edges (253 nodes, 15 communities, weighted *Q* = 0.8036), the
percolation threshold at 458 edges, three residue identities including the DFG
motif, the fraction of network edges that are spatial contacts and their median
length, and both enrichment statistics with the distance-matched interval
spanning 1.0.
