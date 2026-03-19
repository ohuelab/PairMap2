"""PairMap2 – fast FEP intermediate-insertion engine.

Public API::

    from pairmap2 import Pipeline, PipelineConfig

    result = Pipeline(PipelineConfig(input_dir="./input")).run()
"""
from .pipeline import Pipeline
from .types import PipelineConfig, PipelineResult, StageTimings
from .intermediate_search import SearchIntermediates
from .map_generator import MapGenerator

__all__ = [
    "Pipeline",
    "PipelineConfig",
    "PipelineResult",
    "StageTimings",
    "SearchIntermediates",
    "MapGenerator",
]
