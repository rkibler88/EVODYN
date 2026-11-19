# Scope: what EVODYN is, and what it is not

This file exists to answer one question precisely, because the answer changes
what the tool is, how a reviewer should judge it and how much compute a user
needs.

## Does EVODYN run molecular dynamics?

**Yes, by default, and it also accepts a trajectory you already have.**

With `dynamics.mode: equilibrate`, EVODYN runs the simulations. It does so by
invoking [EquilibraTor](https://github.com/Dias-Lab/EquilibraTor) once per
replicate — the GROMACS driver of José D. D. Cediel-Becerra, Jose Cleydson F.
Silva and Raquel Dias ([CSBJ **27**, 5318–5325,
2025](https://doi.org/10.1016/j.csbj.2025.11.034)), which handles topology
generation, ligand
parameterisation through acpype, box construction, solvation, ion placement,
energy minimisation, NVT, NPT, convergence detection and the production stage.
EVODYN supplies the structure and the parameters, isolates each replicate in
its own working directory, and collects the trajectory, the RMSF profile and
the convergence diagnostics.

This is the expensive path. `evodyn plan` prints the exact command and checks
for GROMACS before anything starts.

With `dynamics.mode: trajectory`, EVODYN analyses a trajectory you supply. This
is the streamlined path for testing, for reusing published simulations, and for
users whose MD is run by a facility or a cluster workflow outside EVODYN.

`matrix` consumes a precomputed correlation matrix; `none` skips the dynamic
layer and says so in the report.

## What EVODYN does not do

* It does not implement direct coupling analysis, statistical coupling
  analysis, or molecular dynamics. It calls plmc, pySCA and
  GROMACS-via-EquilibraTor, and records their versions.
* It does not choose a conformational state for you. A run describes the
  starting structure it was given.
* It does not simulate mutants, oligomers or transitions unless you configure
  them as separate runs with their own structures.
* It does not decide whether an alignment is appropriate for a question. It
  reports the alignment's composition so you can.
* It does not produce publication figures. It produces tables, a report and a
  manifest; the figures in the EGFR case study live in that repository.

## What a run's output does and does not support

Every run writes this itself, in the limits section of `report.md`, from what
actually happened rather than from a template.
