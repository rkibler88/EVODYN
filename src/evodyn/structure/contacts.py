"""Spatial geometry: the covariate that every evolutionary claim is matched on.

Proximity filtering by sequence index (``|i-j| > 5``) is near-universal in
coevolution work and does not control for spatial proximity. In the EGFR
alignment the top-coupled "distal" pairs sat at a median 11.4 A Ca-Ca against
24.3 A for the rest. EVODYN therefore computes the Ca-Ca distance for every
pair it scores, and carries it into the pair table whether or not the user
asked for a distance-matched null.
"""
from __future__ import annotations

from itertools import combinations

import networkx as nx
import numpy as np
import pandas as pd

from evodyn.numbering import read_ca_coordinates


def ca_distance_matrix(
    pdb, residues: list[int] | None = None, chain: str | None = None, offset: int = 0
) -> pd.DataFrame:
    """Pairwise Ca-Ca distances in canonical numbering, in angstrom."""
    res, xyz = read_ca_coordinates(pdb, chain=chain, offset=offset)
    if residues is not None:
        idx = {r: k for k, r in enumerate(res)}
        missing = [r for r in residues if r not in idx]
        if missing:
            raise KeyError(
                f"{len(missing)} requested residue(s) absent from {pdb}, "
                f"first few: {missing[:5]}"
            )
        order = [idx[r] for r in residues]
        xyz, res = xyz[order], list(residues)
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    return pd.DataFrame(d, index=res, columns=res)


def pair_distance_table(dist: pd.DataFrame) -> pd.DataFrame:
    """Upper-triangle of a distance matrix as a long table."""
    res = list(dist.index)
    i, j = np.triu_indices(len(res), k=1)
    arr = dist.to_numpy()
    return pd.DataFrame({
        "res_i": np.asarray(res)[i],
        "res_j": np.asarray(res)[j],
        "ca_distance": arr[i, j],
    })


def contact_graph(dist: pd.DataFrame, cutoff: float = 8.0) -> nx.Graph:
    """Contact graph at a Ca-Ca cutoff, used for sector compactness tests."""
    G = nx.Graph()
    G.add_nodes_from(dist.index)
    arr, res = dist.to_numpy(), list(dist.index)
    i, j = np.triu_indices(len(res), k=1)
    for a, b in zip(i[arr[i, j] <= cutoff], j[arr[i, j] <= cutoff]):
        G.add_edge(res[a], res[b], distance=float(arr[a, b]))
    return G


def contact_fraction(dist: pd.DataFrame, pairs: pd.DataFrame, cutoff: float = 8.0) -> float:
    """Fraction of the supplied pairs that are in spatial contact.

    Reported next to any "long-range coupling" claim. In the EGFR network 92%
    of edges were spatial contacts with a median length of 3.9 A, which is the
    expected behaviour of direct coupling analysis and bounds what the network
    can be said to show.
    """
    d = [dist.at[int(a), int(b)] for a, b in zip(pairs.res_i, pairs.res_j)
         if int(a) in dist.index and int(b) in dist.index]
    return float(np.mean(np.asarray(d) <= cutoff)) if d else float("nan")
