"""Data classes for PairMap2 pipeline configuration and results."""
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class PipelineConfig:
    # similarity thresholds
    similarity_threshold: float = 0.6
    min_score_threshold: float = 0.2
    rough_score_threshold: float = 0.5

    # path parameters
    max_path_length: int = 4
    max_optimal_path_length: int = 4

    # graph parameters
    max_intermediate: int = -1
    allow_tree: bool = False
    max_dist_from_actives: int = 6

    # score engine
    jobs: int = -1
    atom_count_diff_threshold: int = 10
    tanimoto_prefilter: float = 0.0

    # cache
    cache_db_path: Optional[str] = None  # None = in-memory only

    # lomap options
    lomap_options: dict = field(default_factory=dict)

    # map generator
    cycle_length: int = 3
    chunk_scale: int = 10
    squared_sum: bool = True
    optimal_path_mode: bool = True

    # search
    is_atom_modification_enabled: bool = True
    cap_ring_with_carbon: bool = True
    cap_ring_with_hydrogen: bool = True
    no_backward_search: bool = False
    use_seed: bool = True
    search_mode: str = "random"  # "bfs" or "random"
    search_random_seed: int = 42
    ionize: bool = True

    # output
    input_dir: str = "./input"
    output_dir: str = "./output"
    save_output: bool = True
    verbose: bool = False


@dataclass
class StageTimings:
    stage: str
    wall_time: float
    cpu_time: float
    peak_rss_mb: float


@dataclass
class PipelineResult:
    graphs: list  # list of nx.Graph (history)
    node_mols: dict  # node_id -> RDKit Mol
    timings: list  # list of StageTimings
