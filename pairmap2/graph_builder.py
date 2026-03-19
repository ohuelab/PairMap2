"""Incremental graph pruner for PairMap2."""
import logging

import networkx as nx

logger = logging.getLogger(__name__)


def _find_bridges(G: nx.Graph) -> set:
    """Return the set of bridge edges in *G* as frozenset pairs.

    Uses iterative Tarjan's bridge-finding algorithm.
    """
    bridges: set = set()
    visited = {u: False for u in G.nodes()}
    low = {u: 0 for u in G.nodes()}
    pre = {u: 0 for u in G.nodes()}
    parent = {u: None for u in G.nodes()}
    count = 0
    for u in G.nodes():
        if not visited[u]:
            stack = [(u, iter(G.neighbors(u)))]
            visited[u] = True
            count += 1
            pre[u] = count
            low[u] = count
            while stack:
                u, children = stack[-1]
                for v in children:
                    if not visited[v]:
                        parent[v] = u
                        stack.append((v, iter(G.neighbors(v))))
                        visited[v] = True
                        count += 1
                        pre[v] = count
                        low[v] = count
                        break
                    elif v != parent[u]:
                        low[u] = min(low[u], pre[v])
                else:
                    stack.pop()
                    if parent[u] is not None:
                        low[parent[u]] = min(low[parent[u]], low[u])
                        if low[u] > pre[parent[u]]:
                            bridges.add(frozenset([u, parent[u]]))
    return bridges
