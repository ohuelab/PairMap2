"""Path-finding utilities using BFS and Yen's K-shortest paths."""
from typing import List, Optional, Tuple

import networkx as nx
import numpy as np


def find_optimal_path(
    graph: nx.Graph,
    source: int,
    target: int,
    max_path_length: int = 4,
    squared_sum: bool = True,
    k_max: int = 100,
    force_path_length: Optional[int] = None,
) -> List[int]:
    """Find the optimal path from *source* to *target* using Yen's algorithm.

    Uses ``nx.shortest_simple_paths`` (lazy evaluation) and stops early once
    ``k_max`` candidates or the path length limit is exceeded.
    Path cost = sum(1/score^2) when ``squared_sum=True``, else -sum(score).
    Raises ``ValueError`` if no suitable path exists.
    """

    def path_cost(path: List[int]) -> float:
        scores = [
            graph[path[k]][path[k + 1]]["score"] for k in range(len(path) - 1)
        ]
        if squared_sum:
            return float(np.sum(1.0 / (np.array(scores) ** 2 + 1e-5)))
        else:
            return float(-np.sum(scores))

    best_path: Optional[List[int]] = None
    best_cost = float("inf")
    count = 0

    try:
        for path in nx.shortest_simple_paths(graph, source, target):
            n_edges = len(path) - 1
            if n_edges > max_path_length:
                break
            count += 1
            if count > k_max:
                break
            if force_path_length is not None and n_edges != force_path_length:
                continue
            cost = path_cost(path)
            if cost < best_cost:
                best_cost = cost
                best_path = path
    except nx.NetworkXNoPath:
        raise ValueError(f"No path found from {source} to {target}")

    if best_path is None:
        raise ValueError(
            f"No path found within max_path_length={max_path_length}"
            + (f" and force_path_length={force_path_length}" if force_path_length is not None else "")
        )
    return best_path


def get_reachable_subgraph(
    graph: nx.Graph,
    source: int,
    target: int,
    max_path_length: int,
    required_nodes: Optional[List[int]] = None,
) -> nx.Graph:
    """Return the subgraph of nodes that lie on some source→target path.

    A node *v* is included iff:
    ``dist(source, v) + dist(v, target) <= max_path_length``

    Uses two BFS passes (O(V+E)).  ``required_nodes`` are always included.
    """
    dist_from_source = nx.single_source_shortest_path_length(
        graph, source, cutoff=max_path_length
    )
    dist_from_target = nx.single_source_shortest_path_length(
        graph, target, cutoff=max_path_length
    )

    reachable = {
        v
        for v in graph.nodes
        if v in dist_from_source
        and v in dist_from_target
        and dist_from_source[v] + dist_from_target[v] <= max_path_length
    }

    if required_nodes:
        reachable.update(required_nodes)

    return graph.subgraph(reachable).copy()


def get_cycled_edges(
    graph: nx.Graph,
    cycle_links: List[Tuple[int, int]],
    cycle_length: int = 3,
) -> set:
    """Return the subset of *cycle_links* already covered by an alternative path.

    For each edge (u, v), temporarily removes it and checks via BFS whether an
    alternative path of length ≤ ``cycle_length - 1`` exists.
    """
    cycled_edges: set = set()
    for u, v in cycle_links:
        if not graph.has_edge(u, v):
            continue
        edge_data = dict(graph[u][v])
        graph.remove_edge(u, v)
        try:
            if nx.has_path(graph, u, v):
                dist = nx.shortest_path_length(graph, u, v)
                if dist <= cycle_length - 1:
                    cycled_edges.add((u, v))
        except nx.NetworkXNoPath:
            pass
        graph.add_edge(u, v, **edge_data)
    return cycled_edges
