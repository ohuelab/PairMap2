"""Subprocess-based worker for PairMap2 Pipeline jobs."""
from __future__ import annotations

import json
import multiprocessing
import os
from pathlib import Path


def _run_job(job_id: str, input_dir: str, config: dict) -> None:
    """Executed in a child process."""
    import sys
    import dataclasses

    project_root = str(Path(__file__).parent.parent.parent)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)

    from web.backend import job_store
    from web.backend.utils import graph_to_cytoscape

    job_dir = Path(input_dir).parent

    job_store.update_job(job_id, status="running")

    try:
        from pairmap2 import Pipeline, PipelineConfig

        default_jobs = int(os.environ.get("PAIRMAP_SCORE_JOBS", "-1"))
        cfg = PipelineConfig(
            input_dir=input_dir,
            output_dir=str(job_dir / "output"),
            save_output=config.get("save_output", True),
            similarity_threshold=config.get("similarity_threshold", 0.6),
            max_path_length=config.get("max_path_length", 4),
            max_intermediate=config.get("max_intermediate", -1),
            jobs=config.get("jobs", default_jobs),
            verbose=config.get("verbose", False),
        )
        pipeline = Pipeline(cfg)
        result = pipeline.run(input_dir=input_dir)

        cy = graph_to_cytoscape(result.graphs[-1], result.node_mols)
        cy["history_length"] = len(result.graphs)
        with open(job_dir / "graph.json", "w") as f:
            json.dump(cy, f)

        timings_data = [
            dataclasses.asdict(t) if dataclasses.is_dataclass(t) else t
            for t in result.timings
        ]
        with open(job_dir / "timings.json", "w") as f:
            json.dump(timings_data, f)

        final = result.graphs[-1]
        job_store.update_job(
            job_id,
            status="done",
            n_nodes=final.number_of_nodes(),
            n_edges=final.number_of_edges(),
        )

    except Exception:
        import traceback
        job_store.update_job(
            job_id,
            status="failed",
            error=traceback.format_exc(),
        )


def submit_job(
    job_id: str,
    input_dir: Path,
    config: dict,
) -> multiprocessing.Process:
    """Spawn a child process to run the pipeline and return it."""
    p = multiprocessing.Process(
        target=_run_job,
        args=(job_id, str(input_dir), config),
        daemon=False,
    )
    p.start()
    return p
