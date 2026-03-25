"""Subprocess-based worker for Map Mode v1 (PairMap engine) jobs."""
from __future__ import annotations

import json
import multiprocessing
import os
from datetime import datetime
from pathlib import Path

_JOBS_DIR = Path(os.environ.get("PAIRMAP_JOBS_DIR", str(Path(__file__).parent.parent / "jobs")))
MAP_JOBS_DIR = _JOBS_DIR / "map"


def _run_job(job_id: str, engine_name: str, config: dict, input_sdf: str) -> None:
    """Executed in a child process."""
    import sys

    project_root = str(Path(__file__).parent.parent.parent)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)

    from web.backend import map_store
    from web.backend.engine import get_engine
    from web.backend.utils import graph_to_cytoscape

    try:
        job_dir = MAP_JOBS_DIR / job_id
        input_dir = job_dir / "input"
        output_dir = job_dir / "output"
        input_dir.mkdir(parents=True, exist_ok=True)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Split multi-molecule SDF into one file per molecule (LOMAP requirement)
        sdf_path = Path(input_sdf)
        raw = sdf_path.read_text()
        records = [r.strip() for r in raw.split("$$$$") if r.strip()]
        if len(records) < 2:
            raise ValueError(f"Input SDF must contain at least 2 molecules, got {len(records)}")
        # Remove the combined SDF so only per-mol files remain in input_dir
        sdf_path.unlink(missing_ok=True)
        for i, molblock in enumerate(records):
            first_line = molblock.splitlines()[0].strip()
            mol_name = first_line if first_line else f"mol_{i}"
            # Sanitise name for use as filename
            safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in mol_name)
            out_path = input_dir / f"{safe_name or f'mol_{i}'}.sdf"
            # Avoid collisions
            if out_path.exists():
                out_path = input_dir / f"{safe_name}_{i}.sdf"
            out_path.write_text(molblock + "\n$$$$\n")

        map_store.update_job(
            job_id,
            status="running",
            started_at=datetime.utcnow().isoformat(),
            progress="Starting engine",
        )

        engine = get_engine(engine_name)
        map_store.update_job(job_id, progress="Running PairMap engine")

        run_config = dict(config)
        run_config["output_dir"] = str(output_dir)
        run_config["save_output"] = True

        result = engine.run(str(input_dir), run_config)

        map_store.update_job(job_id, progress="Saving graph")
        cy = graph_to_cytoscape(result.graphs[-1], result.node_mols)
        cy["history_length"] = len(result.graphs)
        with open(job_dir / "graph.json", "w") as f:
            json.dump(cy, f)

        with open(job_dir / "timings.json", "w") as f:
            import dataclasses
            timings_data = [dataclasses.asdict(t) if dataclasses.is_dataclass(t) else t for t in result.timings]
            json.dump(timings_data, f)

        map_store.update_job(
            job_id,
            status="completed",
            completed_at=datetime.utcnow().isoformat(),
            progress="Done",
        )

    except Exception:
        import traceback
        map_store.update_job(
            job_id,
            status="failed",
            completed_at=datetime.utcnow().isoformat(),
            error=traceback.format_exc(),
        )


def submit_job(
    job_id: str,
    engine_name: str,
    config: dict,
    input_sdf: str,
) -> multiprocessing.Process:
    """Spawn a child process and return it."""
    p = multiprocessing.Process(
        target=_run_job,
        args=(job_id, engine_name, config, input_sdf),
        daemon=False,
    )
    p.start()
    return p
