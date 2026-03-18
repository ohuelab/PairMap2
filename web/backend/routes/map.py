"""Map Mode v1 routes — single SDF → PairMap engine → perturbation map."""
from __future__ import annotations

import json
import uuid
from pathlib import Path

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse

from .. import map_store
from ..map_worker import MAP_JOBS_DIR, submit_job
from ..models import MapJobList, MapJobStatus

router = APIRouter()

_processes: dict[str, object] = {}


@router.post("/jobs", response_model=MapJobStatus, status_code=202)
async def create_map_job(
    file: UploadFile = File(..., description="Input SDF with ligands"),
    engine: str = Form("v1"),
    config: str = Form("{}"),
):
    try:
        cfg = json.loads(config)
    except json.JSONDecodeError as exc:
        raise HTTPException(status_code=422, detail=f"config is not valid JSON: {exc}")

    job_id = str(uuid.uuid4())
    job_dir = MAP_JOBS_DIR / job_id / "input"
    job_dir.mkdir(parents=True, exist_ok=True)

    sdf_path = job_dir / (file.filename or "input.sdf")
    sdf_path.write_bytes(await file.read())

    status = map_store.create_job(job_id, engine, cfg)
    proc = submit_job(job_id, engine, cfg, str(sdf_path))
    _processes[job_id] = proc
    return status


@router.get("/jobs", response_model=MapJobList)
async def list_map_jobs():
    return MapJobList(jobs=map_store.list_jobs())


@router.get("/jobs/{job_id}", response_model=MapJobStatus)
async def get_map_job(job_id: str):
    status = map_store.get_job(job_id)
    if status is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return status


@router.get("/jobs/{job_id}/graph")
async def get_map_graph(job_id: str):
    graph_path = MAP_JOBS_DIR / job_id / "graph.json"
    if not graph_path.exists():
        status = map_store.get_job(job_id)
        if status is None:
            raise HTTPException(status_code=404, detail="Job not found")
        raise HTTPException(
            status_code=404,
            detail=f"Graph not yet available (status: {status.status})",
        )
    with open(graph_path) as f:
        return json.load(f)


@router.get("/jobs/{job_id}/artifacts/{filename}")
async def get_map_artifact(job_id: str, filename: str):
    if ".." in filename or "/" in filename or "\\" in filename:
        raise HTTPException(status_code=400, detail="Invalid filename")
    artifact = MAP_JOBS_DIR / job_id / "output" / filename
    if not artifact.exists():
        raise HTTPException(status_code=404, detail="Artifact not found")
    return FileResponse(str(artifact))


@router.post("/jobs/{job_id}/cancel", response_model=MapJobStatus)
async def cancel_map_job(job_id: str):
    status = map_store.get_job(job_id)
    if status is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if status.status not in ("queued", "running"):
        raise HTTPException(status_code=409, detail=f"Job is already {status.status}")

    proc = _processes.get(job_id)
    if proc is not None:
        try:
            proc.terminate()
        except Exception:
            pass

    from datetime import datetime
    map_store.update_job(
        job_id,
        status="cancelled",
        completed_at=datetime.utcnow().isoformat(),
    )
    return map_store.get_job(job_id)
