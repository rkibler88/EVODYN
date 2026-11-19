"""Centrality, knockout and path analysis on a residue network."""
from __future__ import annotations

from typing import Iterable, Sequence

import networkx as nx
import pandas as pd


def centrality_table(G: nx.Graph) -> pd.DataFrame:
    """Weighted degree, betweenness, closeness and their ranks.

    Ranks matter more than values when comparing proteins of different sizes,
    so both are written out.
    """
    strength = pd.Series(dict(G.degree(weight="weight")))
    between = pd.Series(nx.betweenness_centrality(G, weight="weight"))
    # shortest-path centrality needs distance, not affinity
    dist = {(u, v): 1.0 / d["weight"] for u, v, d in G.edges(data=True) if d["weight"] > 0}
    H = nx.Graph()
    H.add_nodes_from(G.nodes)
    for (u, v), w in dist.items():
        H.add_edge(u, v, weight=w)
    close = pd.Series(nx.closeness_centrality(H, distance="weight"))
    df = pd.DataFrame({
        "node_strength": strength,
        "betweenness": between,
        "closeness": close,
    }).rename_axis("residue").reset_index()
    for col in ("node_strength", "betweenness", "closeness"):
        df[f"{col}_rank"] = df[col].rank(ascending=False, method="min").astype(int)
    return df.sort_values("residue").reset_index(drop=True)


def knockout_table(G: nx.Graph, targets: Iterable[int] | None = None) -> pd.DataFrame:
    """Loss of global efficiency when each residue is removed.

    With ``targets=None`` every node is knocked out, which is O(N) global
    efficiency evaluations and is the expensive option for large networks.
    """
    eff_wt = nx.global_efficiency(G)
    nodes = list(targets) if targets else list(G.nodes)
    rows = []
    for n in nodes:
        if n not in G:
            rows.append({"residue": n, "in_network": False,
                         "efficiency_drop_pct": float("nan")})
            continue
        H = G.copy()
        H.remove_node(n)
        rows.append({
            "residue": n, "in_network": True,
            "efficiency_drop_pct": 100.0 * (eff_wt - nx.global_efficiency(H)) / eff_wt,
        })
    out = pd.DataFrame(rows).sort_values("efficiency_drop_pct", ascending=False)
    out.attrs["wild_type_efficiency"] = eff_wt
    return out.reset_index(drop=True)


def shortest_path(
    G: nx.Graph, source: int, target: int, exclude_backbone: bool = False
) -> list[int] | None:
    """Minimum-cost path, with edge cost the reciprocal of coupling strength.

    ``exclude_backbone=True`` deletes ``|i-j| == 1`` edges first. A path that
    runs along the backbone is not evidence of an allosteric connection, and a
    residue whose only edges are backbone edges has no non-trivial path at all
    - which is itself a result worth reporting rather than hiding.
    """
    H = nx.Graph()
    for u, v, d in G.edges(data=True):
        if exclude_backbone and abs(u - v) == 1:
            continue
        if d["weight"] <= 0:
            continue
        H.add_edge(u, v, weight=1.0 / d["weight"])
    if source not in H or target not in H:
        return None
    try:
        return nx.shortest_path(H, source, target, weight="weight")
    except nx.NetworkXNoPath:
        return None
