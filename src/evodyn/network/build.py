"""Residue interaction network and its community partition.

The network is the top-scoring coupled pairs, taken as edges. The only free
parameter is how many pairs to take, and EVODYN prefers to derive it rather
than set it: ``percolation`` finds the smallest edge count at which the graph
contains every scored position and is connected, which is a property of the
coupling distribution rather than a choice. ``multiplier`` reproduces the
L-proportional convention used in the EGFR study and in much of the DCA
literature.

A warning that the EGFR analysis earned the hard way: average-product-corrected
coupling scores are signed and approximately mean-zero. Summing them over
arbitrary residue-pair sets accumulates background and scales with set size.
:func:`ResidueNetwork.intermodule_coupling` therefore sums over *network edges*
only, and also reports a rank statistic.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal

import community as community_louvain
import networkx as nx
import numpy as np
import pandas as pd


def percolation_threshold(
    couplings: pd.DataFrame,
    n_nodes: int,
    lo: int | None = None,
    hi: int | None = None,
) -> int:
    """Smallest edge count at which the network spans and connects every node.

    Parameters
    ----------
    couplings
        Columns ``res_i``, ``res_j``, ``score``, sorted or not.
    n_nodes
        Number of scored positions the network should contain.

    Returns
    -------
    int
        The edge count. If the graph never connects within the search range the
        largest count tried is returned, and the caller is expected to record
        that in the run limits.
    """
    df = couplings.sort_values("score", ascending=False)
    lo = lo or max(n_nodes - 1, 1)
    hi = hi or min(len(df), 10 * n_nodes)
    for n in range(lo, hi + 1):
        G = _graph(df.head(n))
        if G.number_of_nodes() == n_nodes and nx.is_connected(G):
            return n
    return hi


def _graph(df: pd.DataFrame) -> nx.Graph:
    G = nx.Graph()
    for i, j, s in zip(df["res_i"].to_numpy(), df["res_j"].to_numpy(), df["score"].to_numpy()):
        G.add_edge(int(i), int(j), weight=float(s))
    return G


@dataclass
class ResidueNetwork:
    """A coupling network with its best Louvain partition."""

    graph: nx.Graph
    n_edges: int
    partition: dict[int, int]
    modularity: float
    seed: int
    threshold_rule: str

    # ---------------------------------------------------------------- build
    @classmethod
    def from_couplings(
        cls,
        couplings: pd.DataFrame,
        n_nodes: int | None = None,
        rule: Literal["percolation", "multiplier", "fixed"] = "percolation",
        multiplier: float = 1.80,
        n_focus: int | None = None,
        fixed_edges: int | None = None,
        seeds: int = 100,
    ) -> "ResidueNetwork":
        """Build the network and choose its partition.

        ``seeds`` Louvain runs are performed and the partition with the highest
        weighted modularity is kept, with the winning seed recorded. Louvain is
        stochastic; reporting one run without its seed is not reproducible.
        """
        df = couplings.sort_values("score", ascending=False)
        nodes = n_nodes or len(set(df["res_i"]) | set(df["res_j"]))
        if rule == "fixed":
            if not fixed_edges:
                raise ValueError("rule='fixed' needs fixed_edges")
            n_edges, label = int(fixed_edges), f"fixed({fixed_edges})"
        elif rule == "multiplier":
            L = n_focus or nodes
            n_edges = int(L * multiplier)
            label = f"multiplier({multiplier} x L={L})"
        else:
            n_edges = percolation_threshold(df, nodes)
            label = f"percolation(n_nodes={nodes})"

        G = _graph(df.head(n_edges))
        best_q, best_p, best_seed = -np.inf, None, -1
        for s in range(seeds):
            p = community_louvain.best_partition(G, weight="weight", random_state=s)
            q = community_louvain.modularity(p, G, weight="weight")
            if q > best_q:
                best_q, best_p, best_seed = q, p, s
        return cls(G, n_edges, best_p, float(best_q), best_seed, label)

    # ------------------------------------------------------------- summary
    @property
    def n_nodes(self) -> int:
        return self.graph.number_of_nodes()

    @property
    def n_communities(self) -> int:
        return len(set(self.partition.values()))

    @property
    def connected(self) -> bool:
        return nx.is_connected(self.graph)

    def summary(self) -> dict:
        return {
            "nodes": self.n_nodes,
            "edges": self.n_edges,
            "threshold_rule": self.threshold_rule,
            "communities": self.n_communities,
            "weighted_modularity": round(self.modularity, 4),
            "louvain_seed": self.seed,
            "louvain_seeds_tried": None,
            "connected": self.connected,
        }

    # -------------------------------------------------------------- tables
    def node_table(self) -> pd.DataFrame:
        """Per-residue community, weighted degree and betweenness."""
        return (
            pd.DataFrame({
                "community": pd.Series(self.partition),
                "weighted_degree": pd.Series(dict(self.graph.degree(weight="weight"))),
                "betweenness": pd.Series(
                    nx.betweenness_centrality(self.graph, weight="weight")),
            })
            .rename_axis("residue")
            .reset_index()
            .sort_values("residue")
            .reset_index(drop=True)
        )

    def edge_table(self) -> pd.DataFrame:
        rows = [
            {"res_i": min(u, v), "res_j": max(u, v), "score": d["weight"],
             "sep": abs(u - v),
             "community_i": self.partition[u], "community_j": self.partition[v],
             "intermodule": self.partition[u] != self.partition[v]}
            for u, v, d in self.graph.edges(data=True)
        ]
        return pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)

    def intermodule_coupling(self) -> pd.DataFrame:
        """Coupling between community pairs, over network edges only.

        Both a sum over edges and the median edge rank are reported. The sum is
        interpretable here because it runs over edges, which are by
        construction the strongest pairs; it would not be over an arbitrary
        residue set.
        """
        e = self.edge_table()
        inter = e[e.intermodule].copy()
        if inter.empty:
            return pd.DataFrame(columns=["community_i", "community_j", "coupling",
                                         "edges", "median_rank"])
        inter["rank"] = e["score"].rank(ascending=False, method="min").loc[inter.index]
        inter["pair"] = [tuple(sorted(p)) for p in zip(inter.community_i, inter.community_j)]
        agg = (inter.groupby("pair")
               .agg(coupling=("score", "sum"), edges=("score", "count"),
                    median_rank=("rank", "median"))
               .sort_values("coupling", ascending=False)
               .reset_index())
        agg[["community_i", "community_j"]] = pd.DataFrame(agg["pair"].tolist(), index=agg.index)
        return agg.drop(columns="pair")[
            ["community_i", "community_j", "coupling", "edges", "median_rank"]]
