---
title: "EVODYN: reproducible generation of evolutionary, structural and dynamic datasets of protein allostery"
tags:
  - Python
  - protein allostery
  - coevolution
  - direct coupling analysis
  - molecular dynamics
  - residue interaction networks
authors:
  - name: Ryan Kibler
    corresponding: true
    affiliation: 1
  - name: Raquel Dias
    affiliation: 1
affiliations:
  - name: Department of Microbiology and Cell Science, University of Florida, United States
    index: 1
date: 19 November 2026
bibliography: paper.bib
---

<!-- DRAFT. ORCIDs still to be added for both authors; JOSS requires at least
     one. Affiliations to be confirmed before submission. -->

# Summary

Allostery is studied one protein at a time, with analyses assembled from
several tools that do not share a residue-numbering convention, a null model,
or a record of what was run. EVODYN is a configuration-driven pipeline that
produces the same multi-layer dataset for any protein family: co-evolutionary
couplings from a family alignment, a residue interaction network and its
community structure, molecular dynamics and the correlated motions it produces,
and the geometry and secondary structure that both layers must be controlled
against. A run writes two tables — one row per residue, one row per residue
pair — with identical column names across proteins and in the reference
protein's own numbering, so that a study across many targets is a
concatenation rather than a reconciliation.

EVODYN calls established tools rather than reimplementing them: CD-HIT
[@cdhit] for redundancy filtering, plmc [@plmc] for direct coupling analysis,
pySCA [@pysca] for sector decomposition, and EquilibraTor [@equilibrator] for
the GROMACS [@gromacs] setup-to-production molecular dynamics workflow.
Optional export adapters target Spheronizator [@spheronizator] for voxelised
structural input to machine-learning models and InteracTor
[@interactor2024; @interactor2025] for physics-based interaction features. Its
own contribution is the layer those tools do not supply: a canonical numbering
spine that is asserted rather than assumed, matched null models that are
executed as pipeline stages rather than offered as options, a fixed output
schema, and a per-run manifest and limits report.

# Statement of need

EVODYN generalises a published analysis of the EGFR kinase domain
[@egfrallostery], and the specific failure that motivated it is instructive.
That study reported a 1.83-fold enrichment in correlated motion among
top-coupled distal residue pairs, measured against a null matched on sequence
separation — the field-standard proximity filter. Matched instead on Cα–Cα
distance, the same comparison yields 1.02-fold with a bootstrap 95% interval of
0.94–1.11, reproduced in each of five independent trajectories: the top-coupled
pairs sit at a median 11.4 Å against 24.3 Å for the remainder, and the
apparent evolutionary signal was spatial proximity. Separately, a residue
labelled as the catalytic DFG motif in figures and text was an isoleucine two
positions away, a consequence of three coordinate systems — alignment columns,
UniProt numbering, and a model numbered from one — being mapped by hand.

Both errors are of a kind that a pipeline multiplies. EVODYN is built so that
neither can recur silently. Every enrichment statistic is reported alongside
its matched nulls, and where a matched interval spans unity the run states in
its own report which covariate accounts for the signal. Numbering is
established once, the structure's offset is recovered by comparing residue
identities, and a disagreement halts the run rather than annotating the output.
Every run ends with a limits section written by the stages themselves —
alignment composition, how much of the network is simply a contact map, what
the simulations did and did not sample.

A regression test reproduces sixteen published values of the EGFR case study
from the generalised code, including the Louvain partition and its weighted
modularity ($Q = 0.8036$ at 462 edges), the percolation threshold, the
proportion of network edges that are spatial contacts, and both enrichment
statistics with their intervals. The intended users are groups performing
comparative studies of allosteric mechanism, and groups assembling training
data for machine-learning models of protein function, for whom consistency and
an auditable record across many proteins matter more than any single number.

# Acknowledgements

EVODYN grew out of the EGFR kinase-domain study and a suggestion to turn that
one-protein workflow into a reusable tool.

The molecular dynamics stage runs on **EquilibraTor**, by José D. D.
Cediel-Becerra, Jose Cleydson F. Silva and Raquel Dias [@equilibrator]; EVODYN
supplies parameters and collects outputs, and the simulation workflow itself is
theirs. The optional voxel export targets **Spheronizator**, by Matthew
Richardson, Jose Cleydson Silva, Jose Cediel-Becerra and Raquel Dias
[@spheronizator], and the optional interaction-feature export targets
**InteracTor**, by Jose Cleydson F. Silva, Layla Schuster, Nick Sexson, Melissa
Erdem, Ryan Hulk, Matias Kirst, Marcio F. R. Resende and Raquel Dias
[@interactor2024; @interactor2025].

EVODYN also depends on plmc [@plmc], pySCA [@pysca], CD-HIT [@cdhit], GROMACS
[@gromacs], MDTraj [@mdtraj], NetworkX and the Louvain method [@louvain].

# References
