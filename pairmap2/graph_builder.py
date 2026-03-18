"""Incremental graph pruner for PairMap2."""
import itertools
import logging
from operator import itemgetter
from typing import List, Optional, Tuple

import networkx as nx
import numpy as np

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



class GraphBuilder:
    """Build and prune a PairMap graph via incremental edge removal.

    Uses Tarjan's bridge-finding algorithm (pre-computed, O(1) lookup) and
    BFS-based reachability for constraint checks.  Supports chunk-based edge
    removal for efficiency.

    Parameters
    ----------
    graph:
        Initial graph to prune.  Edges must carry a ``'similarity'`` (or
        ``'score'``) attribute.
    max_path_length:
        Maximum shortest-path distance allowed between any two essential nodes.
    max_dist_from_actives:
        Maximum distance from any non-active essential node to the nearest
        active essential node.
    require_cycle_covering:
        If ``True``, the set of non-cyclic nodes / bridge edges must not
        increase during pruning.
    cycle_length:
        Unused in GraphBuilder; retained for API compatibility with
        ``MapGenerator2``.
    chunk_scale:
        Base for chunk-size computation (default 10).
    chunk_terminate_factor:
        Factor controlling the boundary between fast-chunk and slow-chunk
        phases (default 2.0).
    chunk_mode:
        If ``False``, fall back to simple iterative edge removal.
    ignore_intermediates:
        If ``True``, nodes with ``intermediate=True`` are treated as optional
        and constraints are applied only to essential nodes.
    verbose:
        Extra logging.
    """

    def __init__(
        self,
        graph: nx.Graph,
        max_path_length: int = 4,
        max_dist_from_actives: int = 6,
        require_cycle_covering: bool = True,
        cycle_length: int = 3,
        chunk_scale: int = 10,
        chunk_terminate_factor: float = 2.0,
        chunk_mode: bool = True,
        ignore_intermediates: bool = True,
        verbose: bool = False,
    ):
        self.max_path_length = max_path_length
        self.max_dist_from_actives = max_dist_from_actives
        self.require_cycle_covering = require_cycle_covering
        self.cycle_length = cycle_length
        self.chunk_scale = chunk_scale
        self.chunk_terminate_factor = chunk_terminate_factor
        self.chunk_mode = chunk_mode
        self.verbose = verbose

        if ignore_intermediates:
            self.intermediate_nodes = [
                n for n in graph.nodes if graph.nodes[n].get("intermediate")
            ]
            self.essential_nodes = [
                n for n in graph.nodes if not graph.nodes[n].get("intermediate")
            ]
        else:
            self.intermediate_nodes = []
            self.essential_nodes = list(graph.nodes)

        self._subgraph = graph.copy()
        self._precompute()
        self._minimize_edges()

    def _precompute(self):
        g = self._subgraph
        self._bridges = _find_bridges(g) if len(g.edges) > 0 else set()
        self._non_cycle_nodes = self._find_non_cyclic_nodes(g)
        self._non_cycle_edges = self._find_non_cyclic_edges(g)
        self._dist_to_active_failures = self._count_dist_to_active_failures(g)

    def _is_bridge(self, u: int, v: int) -> bool:
        return frozenset([u, v]) in self._bridges

    def _find_non_cyclic_nodes(self, g: nx.Graph) -> set:
        cycle_nodes = set(itertools.chain.from_iterable(nx.cycle_basis(g)))
        return {n for n in self.essential_nodes if n not in cycle_nodes}

    def _find_non_cyclic_edges(self, g: nx.Graph) -> set:
        all_bridges = _find_bridges(g)
        # Remove bridges that touch intermediate nodes
        return {
            b
            for b in all_bridges
            if not any(n in b for n in self.intermediate_nodes)
        }

    def _count_dist_to_active_failures(self, g: nx.Graph) -> int:
        active_nodes = [
            n for n in self.essential_nodes if g.nodes[n].get("active", False)
        ]
        if not active_nodes:
            return 0
        failures = 0
        for n in self.essential_nodes:
            if g.nodes[n].get("active", False):
                continue
            ok = any(
                nx.has_path(g, n, a)
                and nx.shortest_path_length(g, n, a) <= self.max_dist_from_actives
                for a in active_nodes
            )
            if not ok:
                failures += 1
        return failures

    def _remains_connected(self, g: nx.Graph) -> bool:
        components = list(nx.connected_components(g))
        if len(components) == 1:
            return True
        essential_set = set(self.essential_nodes)
        return any(essential_set.issubset(c) for c in components)

    def _check_cycle_covering(self, g: nx.Graph) -> bool:
        new_non_cycle_nodes = self._find_non_cyclic_nodes(g)
        if new_non_cycle_nodes.difference(self._non_cycle_nodes):
            return False
        new_non_cycle_edges = self._find_non_cyclic_edges(g)
        if new_non_cycle_edges.difference(self._non_cycle_edges):
            return False
        return True

    def _check_max_distance(self, g: nx.Graph) -> bool:
        essential = self.essential_nodes
        for n1 in essential:
            try:
                lengths = nx.single_source_shortest_path_length(
                    g, n1, cutoff=self.max_path_length
                )
            except Exception:
                return False
            for n2 in essential:
                if n1 == n2:
                    continue
                if n2 not in lengths:
                    return False
        return True

    def _check_dist_to_active(self, g: nx.Graph) -> bool:
        return (
            self._count_dist_to_active_failures(g)
            <= self._dist_to_active_failures
        )

    def _check_constraints(self, g: nx.Graph) -> bool:
        if not self._remains_connected(g):
            logger.debug("Constraint fail: connectivity")
            return False
        if self.require_cycle_covering and not self._check_cycle_covering(g):
            logger.debug("Constraint fail: cycle covering")
            return False
        if not self._check_max_distance(g):
            logger.debug("Constraint fail: max distance")
            return False
        if not self._check_dist_to_active(g):
            logger.debug("Constraint fail: dist-to-active")
            return False
        return True

    def _minimize_edges(self):
        subgraph = self._subgraph

        weights = sorted(
            [
                (u, v, d.get("similarity", d.get("score", 0)))
                for u, v, d in subgraph.edges(data=True)
            ],
            key=itemgetter(2),
        )

        if len(subgraph.edges()) <= 2:
            self.result_graph = subgraph
            return

        if not self.chunk_mode:
            for u, v, w in weights:
                if w >= 1.0:
                    continue
                if self._is_bridge(u, v):
                    continue
                edge_data = subgraph.get_edge_data(u, v)
                if edge_data is None:
                    continue
                subgraph.remove_edge(u, v)
                if self._check_constraints(subgraph):
                    self._bridges = (
                        _find_bridges(subgraph) if len(subgraph.edges) > 0 else set()
                    )
                    logger.debug(f"Removed edge {u}--{v} (sim={w:.3f})")
                else:
                    subgraph.add_edge(u, v, **edge_data)
        else:
            N = len(subgraph)
            M = len(weights)
            edges = [(u, v) for u, v, _ in weights]
            data = [subgraph.get_edge_data(u, v) for u, v, _ in weights]

            if M == 0:
                self.result_graph = subgraph
                return

            chunk_size = self.chunk_scale ** int(
                np.log(max(M, 1)) / np.log(self.chunk_scale)
            )
            terminate_n = int(self.chunk_terminate_factor * N)

            fast_idx = list(range(0, max(0, M - terminate_n), chunk_size))
            slow_idx = list(range(max(0, M - terminate_n), M))
            chunk_list = fast_idx + slow_idx

            self._tmp_subgraph = subgraph

            for i, idx_i in enumerate(chunk_list):
                idx_j = chunk_list[i + 1] if i < len(chunk_list) - 1 else M
                edge_chunk = edges[idx_i:idx_j]
                data_chunk = data[idx_i:idx_j]
                if len(edge_chunk) > 1:
                    self._chunk_process(edge_chunk, data_chunk, len(edge_chunk), idx_i)
                else:
                    self._check_chunk(edge_chunk, data_chunk)

            subgraph = self._tmp_subgraph

        # keep only the connected component that contains all essential nodes
        components = list(nx.connected_components(subgraph))
        essential_set = set(self.essential_nodes)
        for comp in components:
            if essential_set.issubset(comp):
                subgraph = subgraph.subgraph(comp).copy()
                break

        self.result_graph = subgraph

    def _check_chunk(self, edge_chunk, data_chunk) -> bool:
        subgraph = self._tmp_subgraph
        removables = [
            d.get("similarity", d.get("score", 0)) < 1.0 for d in data_chunk
        ]
        if not all(removables):
            if not any(removables):
                logger.debug(
                    f"Skip chunk (all similarity=1.0), size={len(edge_chunk)}"
                )
                return True
            return False

        subgraph.remove_edges_from(edge_chunk)
        if self._check_constraints(subgraph):
            logger.debug(f"Removed chunk size={len(edge_chunk)}")
            self._bridges = (
                _find_bridges(subgraph) if len(subgraph.edges) > 0 else set()
            )
            return True
        # Revert
        for (u, v), d in zip(edge_chunk, data_chunk):
            subgraph.add_edge(u, v, **d)
        return False

    def _chunk_process(self, edge_chunk, data_chunk, chunk_size, idx) -> bool:
        if self._check_chunk(edge_chunk, data_chunk):
            return True
        if chunk_size == 1:
            return False
        # Recursively split
        chunk_size = max(chunk_size // self.chunk_scale, 1)
        for i in range(0, len(edge_chunk), chunk_size):
            self._chunk_process(
                edge_chunk[i: i + chunk_size],
                data_chunk[i: i + chunk_size],
                chunk_size,
                idx + i,
            )
        return True

    def get_graph(self) -> nx.Graph:
        return self.result_graph
