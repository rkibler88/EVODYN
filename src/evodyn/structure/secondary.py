"""Secondary structure per residue, over a trajectory or a single structure."""
from __future__ import annotations

import pandas as pd

from evodyn.exceptions import ExternalToolError

#: DSSP eight-state to the three classes used for matched nulls
SS_CLASS = {"H": "helix", "G": "helix", "I": "helix",
            "E": "beta", "B": "beta",
            "T": "loop", "S": "loop", "C": "loop", " ": "loop", "NA": "loop"}


def secondary_structure(
    topology, trajectory=None, stride: int = 1, offset: int = 0, simplified: bool = False
) -> pd.DataFrame:
    """Per-residue DSSP assignment and the fraction of frames in each class.

    With a trajectory, assignment is the modal class over frames and the
    fractions describe how stable it is - a residue that is beta in 55% of
    frames is not the same object as one that is beta in 99%, and only the
    second belongs in a "beta-sheet pairing" analysis.
    """
    try:
        import mdtraj as md
    except ImportError as exc:  # pragma: no cover
        raise ExternalToolError(
            "secondary structure needs mdtraj: pip install 'evodyn[dynamics]'") from exc

    traj = (md.load(trajectory, top=topology, stride=stride) if trajectory
            else md.load(topology))
    ss = md.compute_dssp(traj, simplified=simplified)
    resseq = [r.resSeq + offset for r in traj.topology.residues if r.is_protein]
    ss = ss[:, : len(resseq)]

    rows = []
    for k, res in enumerate(resseq):
        col = pd.Series(ss[:, k]).map(lambda c: SS_CLASS.get(c, "loop"))
        frac = col.value_counts(normalize=True)
        rows.append({
            "residue": res,
            "ss_class": col.mode().iat[0],
            "frac_helix": float(frac.get("helix", 0.0)),
            "frac_beta": float(frac.get("beta", 0.0)),
            "frac_loop": float(frac.get("loop", 0.0)),
            "n_frames": int(traj.n_frames),
        })
    return pd.DataFrame(rows)
