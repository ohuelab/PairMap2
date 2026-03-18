"""Runs PairMap map-generation jobs in a subprocess and serialises results to disk."""
from __future__ import annotations

import json
import multiprocessing
import shutil
from datetime import datetime
from pathlib import Path

from .utils import graph_to_cytoscape

JOBS_DIR = Path(__file__).parent.parent / "jobs"


def _run_job(job_id: str, engine_name: str, config: dict, input_sdf: str) -> None:
    """Executed in a child process."""
    import sys

    project_root = str(Path(__file__).parent.parent.parent)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)

    from PairMapWeb.backend import job_store
    from PairMapWeb.backend.engine import get_engine

    job_dir = JOBS_DIR / job_id
    input_dir = job_dir / "input"
    output_dir = job_dir / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy(input_sdf, input_dir / Path(input_sdf).name)

    job_store.update_job(
        job_id,
        status="running",
        started_at=datetime.utcnow().isoformat(),
        progress="Starting engine",
    )

    try:
        engine = get_engine(engine_name)
        job_store.update_job(job_id, progress="Running PairMap engine")

        run_config = dict(config)
        run_config["output_dir"] = str(output_dir)
        run_config["save_output"] = True

        result = engine.run(str(input_dir), run_config)

        job_store.update_job(job_id, progress="Saving graph")
        cy = graph_to_cytoscape(result.graphs[-1], result.node_mols)
        cy["history_length"] = len(result.graphs)
        with open(job_dir / "graph.json", "w") as f:
            json.dump(cy, f)

        with open(job_dir / "timings.json", "w") as f:
            json.dump(result.timings, f)

        job_store.update_job(
            job_id,
            status="completed",
            completed_at=datetime.utcnow().isoformat(),
            progress="Done",
        )

    except Exception:
        import traceback
        job_store.update_job(
            job_id,
            status="failed",
            completed_at=datetime.utcnow().isoformat(),
            error=traceback.format_exc(),
        )


def submit_job(job_id: str, engine_name: str, config: dict, input_sdf: str) -> multiprocessing.Process:
    """Spawn a child process and return it."""
    p = multiprocessing.Process(
        target=_run_job,
        args=(job_id, engine_name, config, input_sdf),
        daemon=True,
    )
    p.start()
    return p
