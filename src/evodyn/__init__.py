"""EVODYN - evolutionary, structural and dynamic datasets of protein allostery.

Authors: Ryan Kibler and Raquel Dias, PhD (Department of Microbiology and Cell
Science, University of Florida).

EVODYN generalises the single-protein workflow of the EGFR kinase-domain study
(github.com/rkibler88/EGFRAllostery) into a configuration-driven pipeline that
runs over arbitrary protein families.

EVODYN invokes rather than reimplements: plmc for direct coupling analysis,
pySCA for sectors, CD-HIT for redundancy filtering, and - for the dynamic layer
and the machine-learning exports - EquilibraTor, Spheronizator and InteracTor
from the Dias Lab. See CITATION.cff for the authors of each, and cite the tools
whose stages you use.

The organising principle is that *the controls are part of the pipeline*, not
something added afterwards. Every enrichment statistic EVODYN reports is
accompanied by a matched null (distance-matched, conservation-matched or
secondary-structure-matched), every coordinate mapping is asserted rather than
assumed, and every run writes a limits section describing what its own output
does not support.
"""

__version__ = "0.1.0"

from evodyn.config import EvodynConfig, load_config
from evodyn.numbering import CanonicalMap

__all__ = ["EvodynConfig", "load_config", "CanonicalMap", "__version__"]
