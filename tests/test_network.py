"""Network construction, on a graph whose answer is known by hand."""
import numpy as np
import pandas as pd
import pytest

from evodyn.network import ResidueNetwork, knockout_table, percolation_threshold, shortest_path


@pytest.fixture
def two_cliques():
    """Two 5-cliques joined by a single bridge edge."""
    rows = []
    for block in (range(1, 6), range(11, 16)):
        nodes = list(block)
        for a in range(len(nodes)):
            for b in range(a + 1, len(nodes)):
                rows.append((nodes[a], nodes[b], 10.0))
    rows.append((5, 11, 1.0))                      # the bridge
    rows.append((1, 15, 0.01))                     # a weak distractor
    return pd.DataFrame(rows, columns=["res_i", "res_j", "score"])


def test_percolation_finds_the_bridge(two_cliques):
    # 20 clique edges connect nothing across; the 21st (the bridge) does
    assert percolation_threshold(two_cliques, n_nodes=10) == 21


def test_partition_recovers_the_two_cliques(two_cliques):
    net = ResidueNetwork.from_couplings(two_cliques, n_nodes=10, rule="percolation", seeds=20)
    assert net.n_nodes == 10
    assert net.n_communities == 2
    assert net.modularity > 0.3
    labels = net.partition
    assert len({labels[n] for n in range(1, 6)}) == 1
    assert labels[1] != labels[11]


def test_multiplier_rule_uses_the_focus_length(two_cliques):
    net = ResidueNetwork.from_couplings(two_cliques, n_nodes=10, rule="multiplier",
                                        multiplier=2.0, n_focus=10, seeds=5)
    assert net.n_edges == 20


def test_knockout_ranks_the_bridge_first(two_cliques):
    net = ResidueNetwork.from_couplings(two_cliques, n_nodes=10, rule="percolation", seeds=5)
    out = knockout_table(net.graph)
    assert out.iloc[0].residue in (5, 11)
    assert out.iloc[0].efficiency_drop_pct > 0


def test_backbone_exclusion_can_disconnect(two_cliques):
    """A residue whose only edges are backbone edges has no non-trivial path."""
    df = pd.DataFrame([(1, 2, 5.0), (2, 3, 5.0), (3, 4, 5.0)],
                      columns=["res_i", "res_j", "score"])
    net = ResidueNetwork.from_couplings(df, n_nodes=4, rule="fixed", fixed_edges=3, seeds=3)
    assert shortest_path(net.graph, 1, 4) == [1, 2, 3, 4]
    assert shortest_path(net.graph, 1, 4, exclude_backbone=True) is None


def test_intermodule_coupling_runs_over_edges_only(two_cliques):
    net = ResidueNetwork.from_couplings(two_cliques, n_nodes=10, rule="percolation", seeds=20)
    inter = net.intermodule_coupling()
    assert len(inter) == 1
    assert inter.iloc[0].edges == 1
