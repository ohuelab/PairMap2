"""Pydantic models for PairMap2 WebUI API.

Dependencies required (all present in pyproject.toml):
  fastapi>=0.100, uvicorn[standard]>=0.20, python-multipart>=0.0.5
"""
from __future__ import annotations

from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel


class JobStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    DONE = "done"
    FAILED = "failed"
    CANCELLED = "cancelled"


class JobConfig(BaseModel):
    similarity_threshold: float = 0.6
    max_path_length: int = 4
    max_intermediate: int = -1
    jobs: int = -1
    verbose: bool = False


class JobCreate(BaseModel):
    config: JobConfig = JobConfig()


class StageTimingModel(BaseModel):
    stage: str
    wall_time: float
    cpu_time: float
    peak_rss_mb: float


class JobResponse(BaseModel):
    job_id: str
    status: JobStatus
    created_at: float
    updated_at: float
    config: Dict[str, Any]
    timings: List[StageTimingModel] = []
    error: Optional[str] = None
    n_nodes: Optional[int] = None
    n_edges: Optional[int] = None


class GraphNode(BaseModel):
    id: str
    label: str
    is_intermediate: bool
    smiles: Optional[str] = None


class GraphEdge(BaseModel):
    source: str
    target: str
    score: float


class GraphResponse(BaseModel):
    nodes: List[GraphNode]
    edges: List[GraphEdge]


# ── Map Mode v1 job models ────────────────────────────────────────────────────

class MapJobStatus(BaseModel):
    id: str
    status: str
    engine: str
    config: Dict[str, Any]
    created_at: str
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    error: Optional[str] = None
    progress: Optional[str] = None


class MapJobList(BaseModel):
    jobs: List[MapJobStatus]


# ── Pair Mode models ──────────────────────────────────────────────────────────

from pydantic import Field  # noqa: E402


class SearchConfig(BaseModel):
    is_atom_modfication_enabled: bool = True
    cap_ring_with_carbon: bool = True
    cap_ring_with_hydrogen: bool = True
    no_backward_search: bool = False
    use_seed: bool = True
    max_intermediate: int = 100
    ionize: bool = False


class MapGenConfig(BaseModel):
    maxOptimalPathLength: int = 3
    roughScoreThreshold: float = Field(0.5, ge=0.0, le=1.0)
    minScoreThreshold: float = Field(0.2, ge=0.0, le=1.0)
    optimal_path_mode: bool = True
    CycleLinkThreshold: float = Field(0.6, ge=0.0, le=1.0)
    squared_sum: bool = True


class PairRequest(BaseModel):
    smiles_a: str
    smiles_b: str
    name_a: str = "Molecule A"
    name_b: str = "Molecule B"
    search: SearchConfig = Field(default_factory=SearchConfig)
    mapgen: MapGenConfig = Field(default_factory=MapGenConfig)


# ── Map Mode v1 config ────────────────────────────────────────────────────────

class MapV1Config(BaseModel):
    similarity_threshold: float = Field(0.6, ge=0.0, le=1.0)
    max_intermediate: int = -1
    max: int = 6
    max_dist_from_actives: int = 6
    allow_tree: bool = False
    radial: bool = False
    max_path_length: int = 4
    maxOptimalPathLength: int = 4
    roughScoreThreshold: float = Field(0.5, ge=0.0, le=1.0)
    optimal_path_mode: bool = True
    minScoreThreshold: float = Field(0.2, ge=0.0, le=1.0)
