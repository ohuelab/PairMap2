"""PairMap2 engine — wraps pairmap2.Pipeline."""
from __future__ import annotations

import os

from .base import EngineResult, PairMapEngine


class PairMapDefaultEngine(PairMapEngine):
    """Wraps pairmap2.Pipeline."""

    def run(self, input_dir: str, config: dict) -> EngineResult:
        from pairmap2 import Pipeline, PipelineConfig

        default_jobs = int(os.environ.get("PAIRMAP_SCORE_JOBS", "-1"))
        cfg = PipelineConfig(
            input_dir=input_dir,
            output_dir=config.get("output_dir", "./output"),
            save_output=config.get("save_output", False),
            similarity_threshold=config.get("similarity_threshold", 0.6),
            max_path_length=config.get("max_path_length", 4),
            max_intermediate=config.get("max_intermediate", -1),
            jobs=config.get("jobs", default_jobs),
            verbose=config.get("verbose", False),
            ionize=config.get("ionize", True),
        )
        pipeline = Pipeline(cfg)
        result = pipeline.run(input_dir=input_dir)
        return EngineResult(
            graphs=result.graphs,
            node_mols=result.node_mols,
            timings=result.timings,
        )

    def run_from_moldf(self, mols: list, df, config: dict) -> EngineResult:
        from pairmap2 import Pipeline, PipelineConfig

        default_jobs = int(os.environ.get("PAIRMAP_SCORE_JOBS", "-1"))
        cfg = PipelineConfig(
            save_output=False,
            similarity_threshold=config.get("similarity_threshold", 0.6),
            max_path_length=config.get("max_path_length", 4),
            max_intermediate=config.get("max_intermediate", -1),
            jobs=config.get("jobs", default_jobs),
            verbose=config.get("verbose", False),
            ionize=config.get("ionize", True),
        )
        pipeline = Pipeline(cfg)
        result = pipeline.run_from_moldf(mols, df)
        return EngineResult(
            graphs=result.graphs,
            node_mols=result.node_mols,
            timings=result.timings,
        )
