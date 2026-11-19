"""Dynamic cross-correlation and fluctuation profiles.

The correlation is computed on Ca displacements from the mean position:

.. math::
    C_{ij} = \\frac{\\langle \\Delta r_i \\cdot \\Delta r_j \\rangle}
                   {\\sqrt{\\langle \\Delta r_i^2 \\rangle
                            \\langle \\Delta r_j^2 \\rangle}}

after superposition onto the mean structure, which removes global rotation and
translation. Both a trajectory and a GROMACS ``gmx covar -ascii`` matrix are
accepted; the second path exists so that a correlation matrix computed on a
cluster can be brought in without moving a multi-gigabyte trajectory.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from evodyn.exceptions import ExternalToolError


def _require_mdtraj():
    try:
        import mdtraj as md
        return md
    except ImportError as exc:  # pragma: no cover
        raise ExternalToolError(
            "trajectory analysis needs mdtraj: pip install 'evodyn[dynamics]'") from exc


def dccm_from_trajectory(
    trajectory, topology, residues: list[int] | None = None,
    stride: int = 1, offset: int = 0, superpose: bool = True,
) -> pd.DataFrame:
    """Ca-Ca dynamic cross-correlation matrix, indexed by canonical residue."""
    md = _require_mdtraj()
    traj = md.load(trajectory, top=topology, stride=stride)
    ca = traj.topology.select("name CA and protein")
    traj = traj.atom_slice(ca)
    if superpose:
        traj.superpose(traj, 0)
        traj.superpose(traj, 0, atom_indices=None, ref_atom_indices=None)
    xyz = traj.xyz.astype(np.float64)               # (frames, atoms, 3)
    resseq = [r.resSeq + offset for r in traj.topology.residues]

    delta = xyz - xyz.mean(axis=0, keepdims=True)
    cov = np.einsum("fid,fjd->ij", delta, delta) / xyz.shape[0]
    sd = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sd, sd)
    out = pd.DataFrame(corr, index=resseq, columns=resseq)
    if residues is not None:
        keep = [r for r in residues if r in out.index]
        out = out.loc[keep, keep]
    return out


def dccm_from_covariance(
    ascii_path, n_atoms: int | None = None, residues: list[int] | None = None
) -> pd.DataFrame:
    """Convert a ``gmx covar -ascii`` 3N x 3N covariance matrix to an N x N DCCM."""
    flat = np.loadtxt(ascii_path).ravel()
    n3 = int(round(np.sqrt(flat.size)))
    if n3 * n3 != flat.size:
        raise ValueError(f"{ascii_path}: {flat.size} values is not a square matrix")
    cov3 = flat.reshape(n3, n3)
    n = n3 // 3
    if n_atoms and n_atoms != n:
        raise ValueError(f"{ascii_path} holds {n} atoms, expected {n_atoms}")
    blocks = cov3.reshape(n, 3, n, 3)
    cov = np.einsum("idjd->ij", blocks)
    sd = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sd, sd)
    idx = residues if residues is not None and len(residues) == n else list(range(1, n + 1))
    return pd.DataFrame(corr, index=idx, columns=idx)


def rmsf_from_trajectory(
    trajectory, topology, stride: int = 1, offset: int = 0
) -> pd.DataFrame:
    """Per-residue Ca RMSF in nanometres, to match GROMACS units."""
    md = _require_mdtraj()
    traj = md.load(trajectory, top=topology, stride=stride)
    traj = traj.atom_slice(traj.topology.select("name CA and protein"))
    traj.superpose(traj, 0)
    rmsf = md.rmsf(traj, traj, 0)
    return pd.DataFrame({
        "residue": [r.resSeq + offset for r in traj.topology.residues],
        "rmsf": rmsf,
    })


def read_xvg(path, offset: int = 0) -> pd.DataFrame:
    """Read a two-column GROMACS ``.xvg`` (``gmx rmsf -res`` and friends)."""
    rows = [l.split() for l in open(path) if l[:1] not in ("#", "@", "")]
    arr = np.asarray(rows, dtype=float)
    return pd.DataFrame({"residue": arr[:, 0].astype(int) + offset, "rmsf": arr[:, 1]})
