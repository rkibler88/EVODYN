"""The two tables that are EVODYN's actual product.

Everything a run computes lands in one of two shapes, both in canonical
numbering, both one file per run, both with a column naming convention that
does not change between proteins:

``residues.csv``
    one row per scored residue - amino acid, conservation, community,
    centrality, RMSF, secondary structure, site membership.
``pairs.csv``
    one row per scored residue pair - coupling, sequence separation, Ca-Ca
    distance, dynamic correlation, community pair, contact flags.

A comparative or machine-learning study over many proteins is then a
concatenation, which is the whole purpose: one case study becomes a dataset
generator.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

#: columns that must be present and identically named in every run's output
RESIDUE_SCHEMA = ["protein", "residue", "aa"]
PAIR_SCHEMA = ["protein", "res_i", "res_j", "sep"]


def build_residue_table(
    protein: str,
    cmap,
    network=None,
    conservation: pd.DataFrame | None = None,
    rmsf: pd.DataFrame | None = None,
    secondary: pd.DataFrame | None = None,
    sites: dict[str, list[int]] | None = None,
) -> pd.DataFrame:
    """Assemble the per-residue table from whichever layers a run produced."""
    df = pd.DataFrame({"residue": cmap.residues})
    df["protein"] = protein
    df["aa"] = df.residue.map(cmap.res2aa)

    if network is not None:
        df = df.merge(network.node_table(), on="residue", how="left")
        from evodyn.network.metrics import centrality_table
        cen = centrality_table(network.graph)
        df = df.merge(cen[["residue", "closeness", "node_strength_rank",
                           "betweenness_rank"]], on="residue", how="left")
    if conservation is not None:
        df = df.merge(conservation[["residue", "conservation"]], on="residue", how="left")
    if rmsf is not None:
        df = df.merge(rmsf[["residue", "rmsf"]], on="residue", how="left")
    if secondary is not None:
        keep = [c for c in ("residue", "ss_class", "frac_helix", "frac_beta", "frac_loop")
                if c in secondary]
        df = df.merge(secondary[keep], on="residue", how="left")
    for name, members in (sites or {}).items():
        df[f"site_{name}"] = df.residue.isin(members)

    front = [c for c in RESIDUE_SCHEMA if c in df]
    return df[front + [c for c in df.columns if c not in front]]


def build_pair_table(
    protein: str,
    couplings: pd.DataFrame,
    distances: pd.DataFrame | None = None,
    dccm: pd.DataFrame | None = None,
    network=None,
    residue_table: pd.DataFrame | None = None,
    contact_cutoff: float = 8.0,
) -> pd.DataFrame:
    """Assemble the per-pair table.

    The Ca-Ca distance is attached whether or not the user asked for a
    distance-matched null, because an enrichment reported without it cannot be
    interpreted.
    """
    df = couplings[["res_i", "res_j", "score"]].copy()
    df.insert(0, "protein", protein)
    df["sep"] = (df.res_i - df.res_j).abs()

    if distances is not None:
        idx = set(distances.index)
        df["ca_distance"] = [
            distances.at[a, b] if a in idx and b in idx else np.nan
            for a, b in zip(df.res_i, df.res_j)]
        df["contact"] = df.ca_distance <= contact_cutoff
    if dccm is not None:
        idx = set(dccm.index)
        df["dccm"] = [dccm.at[a, b] if a in idx and b in idx else np.nan
                      for a, b in zip(df.res_i, df.res_j)]
        df["abs_dccm"] = df.dccm.abs()
    if network is not None:
        part = network.partition
        df["community_i"] = df.res_i.map(part)
        df["community_j"] = df.res_j.map(part)
        df["intermodule"] = df.community_i.notna() & (df.community_i != df.community_j)
        edges = {tuple(sorted(e)) for e in network.graph.edges}
        df["network_edge"] = [tuple(sorted((a, b))) in edges
                              for a, b in zip(df.res_i, df.res_j)]
    if residue_table is not None and "conservation" in residue_table:
        cons = residue_table.set_index("residue")["conservation"]
        ci, cj = df.res_i.map(cons), df.res_j.map(cons)
        df["conservation_min"] = np.minimum(ci, cj)
        df["conservation_mean"] = (ci + cj) / 2.0

    front = [c for c in PAIR_SCHEMA if c in df]
    return df[front + [c for c in df.columns if c not in front]]


def write_table(df: pd.DataFrame, path, compress: bool = True) -> pathlib.Path:
    path = pathlib.Path(path)
    if compress and path.suffix != ".gz":
        path = path.with_suffix(path.suffix + ".gz")
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path
