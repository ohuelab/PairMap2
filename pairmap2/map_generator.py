"""MapGenerator2 – map generation using BFS-based path/subgraph utilities."""
import itertools
import logging
from typing import List, Optional, Tuple

import networkx as nx
import numpy as np

from .path_solver import (
    find_optimal_path,
    get_reachable_subgraph,
    get_cycled_edges,
)
from .score_engine import ScoreEngine
from .score_cache import ScoreCache

logger = logging.getLogger(__name__)


class MapGenerator2:
    """Generate a pairmap between source and target through intermediates.

    Parameters
    ----------
    intermediate_list:
        List of RDKit molecules.  Index 0 = source, index 1 = target.
    optimal_path_mode:
        Return only the optimal path graph (no full pruning).
    max_path_length:
        Maximum number of hops in the final map.
    cycle_length:
        Maximum cycle length for cycle-covering check.
    max_optimal_path_length:
        Maximum hops considered for the optimal path search.
    rough_score_threshold:
        Score threshold for the rough path check (informational only).
    min_score_threshold:
        Minimum score for an edge to appear in the initial graph.
    cycle_link_threshold:
        Edges on the optimal path below this score must be cycle-covered.
    force_optimal_path_length:
        Restrict optimal-path search to exactly ``max_optimal_path_length``.
    chunk_scale:
        Base for the chunk-size computation.
    squared_sum:
        Use ``sum(1/score^2)`` as path cost if ``True``, else ``-sum(score)``.
    source_node_index:
        Index of the source molecule in ``intermediate_list``.
    target_node_index:
        Index of the target molecule in ``intermediate_list``.
    score_engine:
        Pre-configured ``ScoreEngine``.  Created fresh (in-memory cache) if
        not supplied.
    custom_score_matrix:
        Pre-computed N×N score matrix.  Takes precedence over ``score_engine``.
    verbose:
        Extra logging.
    lomap_options:
        Options forwarded to the LOMAP MCS scorer.
    jobs:
        Parallel workers for score-matrix computation.
    """

    def __init__(
        self,
        intermediate_list: list,
        optimal_path_mode: bool = False,
        max_path_length: int = 4,
        cycle_length: int = 3,
        max_optimal_path_length: int = 4,
        rough_score_threshold: float = 0.5,
        min_score_threshold: float = 0.2,
        cycle_link_threshold: float = 0.6,
        force_optimal_path_length: bool = False,
        chunk_scale: int = 10,
        squared_sum: bool = True,
        source_node_index: int = 0,
        target_node_index: int = 1,
        score_engine: Optional[ScoreEngine] = None,
        custom_score_matrix: Optional[np.ndarray] = None,
        verbose: bool = False,
        lomap_options: Optional[dict] = None,
        jobs: int = -1,
    ):
        self.intermediate_list = intermediate_list
        self.intermediate_names = [
            m.GetProp("_Name") if m.HasProp("_Name") else f"intermediate-{i:04d}"
            for i, m in enumerate(intermediate_list)
        ]

        if custom_score_matrix is not None:
            n = len(intermediate_list)
            if len(custom_score_matrix) != n or len(custom_score_matrix[0]) != n:
                raise ValueError(
                    "custom_score_matrix shape does not match intermediate_list length."
                )
            self.score_matrix: Optional[np.ndarray] = np.asarray(custom_score_matrix)
        else:
            self.score_matrix = None

        self.N = len(intermediate_list)
        self.source_node_index = source_node_index
        self.target_node_index = target_node_index

        self.optimal_path_mode = optimal_path_mode
        self.max_path_length = max_path_length
        self.cycle_length = cycle_length
        self.max_optimal_path_length = max_optimal_path_length
        self.rough_score_threshold = rough_score_threshold
        self.min_score_threshold = min_score_threshold
        self.cycle_link_threshold = cycle_link_threshold
        self.force_optimal_path_length = force_optimal_path_length
        self.chunk_scale = chunk_scale
        self.squared_sum = squared_sum
        self.verbose = verbose
        self.lomap_options = lomap_options or {}

        # Score engine (used only when score_matrix is not pre-supplied)
        self._score_engine = score_engine or ScoreEngine(
            cache=ScoreCache(), jobs=jobs
        )

        self.found_path: List[int] = [source_node_index, target_node_index]
        self.found_links: List[Tuple[int, int]] = [
            (source_node_index, target_node_index)
        ]
        self.cycle_links: List[Tuple[int, int]] = []

    def _get_score_matrix(self) -> np.ndarray:
        if self.score_matrix is None:
            self.score_matrix = self._score_engine.get_score_matrix(
                self.intermediate_list, self.lomap_options
            )
        return self.score_matrix

    def _make_graph(
        self,
        min_score: Optional[float] = None,
        forced_links: Optional[List[Tuple[int, int]]] = None,
    ) -> nx.Graph:
        """Build graph from score matrix with edges above *min_score*."""
        if min_score is None:
            min_score = self.min_score_threshold
        forced_links = forced_links or []
        sm = self._get_score_matrix()

        g = nx.Graph()
        for i, name in enumerate(self.intermediate_names):
            g.add_node(i, label=name)

        for u, v in itertools.combinations(range(self.N), 2):
            score = round(float(sm[u][v]), 2)
            is_forced = (u, v) in forced_links or (v, u) in forced_links
            if score >= min_score or is_forced:
                g.add_edge(u, v, score=score)
        return g

    def _make_optimal_path_graph(self) -> nx.Graph:
        sm = self._get_score_matrix()
        g = nx.Graph()
        for i in self.found_path:
            g.add_node(i, label=self.intermediate_names[i])
        for k in range(len(self.found_path) - 1):
            u, v = self.found_path[k], self.found_path[k + 1]
            g.add_edge(u, v, score=float(sm[u][v]))
        return g

    def _find_optimal_path(self) -> List[int]:
        """Determine the optimal path using Yen's K-shortest paths."""
        sm = self._get_score_matrix()

        # Rough-score check (informational)
        g_rough = self._make_graph(self.rough_score_threshold)
        src, tgt = self.source_node_index, self.target_node_index
        if nx.has_path(g_rough, src, tgt):
            dist = nx.shortest_path_length(g_rough, src, tgt)
            if dist <= 2:
                logger.warning(
                    "Short high-quality path already exists; "
                    "introducing intermediates may be unnecessary."
                )

        g = self._make_graph()
        force_len = (
            self.max_optimal_path_length if self.force_optimal_path_length else None
        )
        found_path = find_optimal_path(
            g,
            src,
            tgt,
            max_path_length=self.max_optimal_path_length,
            squared_sum=self.squared_sum,
            force_path_length=force_len,
        )

        self.found_path = found_path
        self.found_links = [
            (found_path[k], found_path[k + 1])
            if found_path[k] < found_path[k + 1]
            else (found_path[k + 1], found_path[k])
            for k in range(len(found_path) - 1)
        ]
        self.cycle_links = [
            (u, v)
            for u, v in self.found_links
            if round(float(sm[u][v]), 2) < self.cycle_link_threshold
        ]
        return self.found_path

    def _check_optimal_path(self, g: nx.Graph) -> bool:
        return all(g.has_edge(u, v) for u, v in self.found_links)

    def _check_cycle_covering(self, g: nx.Graph) -> bool:
        current_cycled = get_cycled_edges(g, self.cycle_links, self.cycle_length)
        return not self._initial_cycled_edges.difference(current_cycled)

    def _check_constraints(self, g: nx.Graph) -> bool:
        if not self._check_optimal_path(g):
            return False
        if not self._check_cycle_covering(g):
            return False
        return True

    def _get_main_subgraph(self, g: nx.Graph) -> nx.Graph:
        for nodes in nx.connected_components(g):
            if all(n in nodes for n in self.found_path):
                return g.subgraph(nodes).copy()
        raise ValueError("found_path nodes are not in the same connected component")

    def _get_reachable_subgraph(self, g: nx.Graph) -> nx.Graph:
        sub = get_reachable_subgraph(
            g,
            self.source_node_index,
            self.target_node_index,
            self.max_path_length,
            required_nodes=list(self.found_path),
        )
        if not all(n in sub.nodes for n in self.found_path):
            raise ValueError("found_path nodes missing after get_reachable_subgraph")
        return sub

    def _check_chunk(self, edge_chunk, data_chunk) -> bool:
        subgraph = self._tmp_subgraph
        removables = [
            d["score"] < 1.0 and not d.get("found_path", False)
            for d in data_chunk
        ]
        if not all(removables):
            return not any(removables)  # skip-all if none removable, else False

        subgraph.remove_edges_from(edge_chunk)

        try:
            ex = self._get_reachable_subgraph(subgraph)
            ex = self._get_main_subgraph(ex).copy()
        except ValueError:
            for (u, v), d in zip(edge_chunk, data_chunk):
                subgraph.add_edge(u, v, **d)
            return False

        if not self._check_constraints(ex):
            for (u, v), d in zip(edge_chunk, data_chunk):
                subgraph.add_edge(u, v, **d)
            return False

        self._tmp_subgraph = ex
        return True

    def _chunk_process(self, edge_chunk, data_chunk, chunk_size, idx) -> bool:
        if self._check_chunk(edge_chunk, data_chunk):
            return True
        if chunk_size == 1:
            return False
        chunk_size = max(chunk_size // self.chunk_scale, 1)
        crt = 0
        while crt < len(edge_chunk):
            subgraph = self._tmp_subgraph
            edge_in, data_in = [], []
            while len(edge_in) < chunk_size and crt < len(edge_chunk):
                u, v = edge_chunk[crt]
                if subgraph.get_edge_data(u, v):
                    edge_in.append(edge_chunk[crt])
                    data_in.append(data_chunk[crt])
                crt += 1
            ret = self._chunk_process(edge_in, data_in, chunk_size, idx + crt)
            if not ret:
                rest_e = [
                    (u, v)
                    for u, v in edge_chunk[crt:]
                    if self._tmp_subgraph.get_edge_data(u, v) is not None
                ]
                rest_d = [
                    d
                    for (u, v), d in zip(edge_chunk[crt:], data_chunk[crt:])
                    if self._tmp_subgraph.get_edge_data(u, v) is not None
                ]
                if self._check_chunk(rest_e, rest_d):
                    break
        return True

    def build_map(self) -> nx.Graph:
        """Generate and return the final pairmap graph."""
        _ = self._get_score_matrix()
        self._find_optimal_path()

        if self.verbose:
            logger.info(f"Found path: {self.found_path}")
            logger.info(f"Found links: {self.found_links}")

        self.optimal_path_graph = self._make_optimal_path_graph()
        if self.optimal_path_mode:
            self.final_graph = self.optimal_path_graph
            return self.final_graph

        subgraph = self._make_graph(forced_links=self.found_links)
        for u, v in subgraph.edges:
            subgraph[u][v]["found_path"] = False
        for k in range(len(self.found_path) - 1):
            u, v = self.found_path[k], self.found_path[k + 1]
            subgraph[u][v]["found_path"] = True

        scores_list = sorted(subgraph.edges(data="score"), key=lambda e: e[2])
        edges = [(u, v) for u, v, _ in scores_list]
        data = [subgraph[u][v] for u, v, _ in scores_list]

        M = len(scores_list)
        chunk_size = self.chunk_scale ** int(
            np.log(max(M, 1)) / np.log(self.chunk_scale)
        )

        self._initial_cycled_edges = get_cycled_edges(
            subgraph, self.cycle_links, self.cycle_length
        )
        if self.verbose:
            logger.info(f"Initial cycled edges: {self._initial_cycled_edges}")

        subgraph = self._get_main_subgraph(subgraph).copy()
        self._tmp_subgraph = subgraph

        crt = 0
        while crt < len(data):
            edge_chunk, data_chunk = [], []
            while len(edge_chunk) < chunk_size and crt < len(data):
                u, v = edges[crt]
                if self._tmp_subgraph.get_edge_data(u, v):
                    edge_chunk.append(edges[crt])
                    data_chunk.append(data[crt])
                crt += 1
            self._chunk_process(edge_chunk, data_chunk, chunk_size, crt)
            self._tmp_subgraph = self._get_main_subgraph(self._tmp_subgraph).copy()

        subgraph = self._tmp_subgraph.copy()
        final = self._get_reachable_subgraph(subgraph)
        final = self._get_main_subgraph(final).copy()

        self.final_graph = final
        return self.final_graph
