"""Job CRUD endpoints for PairMap2 WebUI."""
from __future__ import annotations

import json
import uuid
from pathlib import Path
from typing import List, Optional

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse

from .. import job_store
from ..job_store import JOBS_DIR
from ..models import JobResponse, JobStatus
from ..worker import submit_job

router = APIRouter()


@router.post("/jobs", response_model=JobResponse, status_code=202)
async def create_job(
    source_sdf: UploadFile = File(..., description="Source molecules SDF file"),
    target_sdf: UploadFile = File(..., description="Target molecules SDF file"),
    config: str = Form(
        "{}",
        description=(
            "Optional JSON config (keys: similarity_threshold, max_path_length, "
            "max_intermediate, jobs, verbose)"
        ),
    ),
):
    """Accept source + target SDF uploads, create a job, and queue it."""
    try:
        cfg: dict = json.loads(config)
    except json.JSONDecodeError as exc:
        raise HTTPException(status_code=422, detail=f"config is not valid JSON: {exc}")

    job_id = str(uuid.uuid4())
    input_dir = JOBS_DIR / job_id / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    # Save uploaded files
    source_path = input_dir / "source.sdf"
    target_path = input_dir / "target.sdf"
    source_path.write_bytes(await source_sdf.read())
    target_path.write_bytes(await target_sdf.read())

    status = job_store.create_job(job_id, cfg)
    submit_job(job_id, input_dir, cfg)
    return status


@router.get("/jobs", response_model=List[JobResponse])
async def list_jobs():
    return job_store.list_jobs()


@router.get("/jobs/{job_id}", response_model=JobResponse)
async def get_job(job_id: str):
    status = job_store.get_job(job_id)
    if status is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return status


@router.get("/jobs/{job_id}/graph")
async def get_graph(job_id: str):
    """Return Cytoscape.js-compatible graph data for a completed job."""
    graph_path = JOBS_DIR / job_id / "graph.json"
    if not graph_path.exists():
        job = job_store.get_job(job_id)
        if job is None:
            raise HTTPException(status_code=404, detail="Job not found")
        raise HTTPException(
            status_code=404,
            detail=f"Graph not yet available (job status: {job.status})",
        )
    with open(graph_path) as f:
        return json.load(f)


@router.get("/jobs/{job_id}/artifacts/{filename}")
async def get_artifact(job_id: str, filename: str):
    """Download a result artifact file (e.g. intermediate_mols.sdf)."""
    # Prevent path traversal
    if ".." in filename or "/" in filename or "\\" in filename:
        raise HTTPException(status_code=400, detail="Invalid filename")
    artifact = JOBS_DIR / job_id / "output" / filename
    if not artifact.exists():
        raise HTTPException(status_code=404, detail="Artifact not found")
    return FileResponse(str(artifact))


@router.post("/jobs/{job_id}/cancel", response_model=JobResponse)
async def cancel_job(job_id: str):
    """Cancel a queued or running job."""
    job = job_store.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.status not in (JobStatus.PENDING, JobStatus.RUNNING):
        raise HTTPException(
            status_code=409, detail=f"Job is already {job.status.value}"
        )
    result = job_store.cancel_job(job_id)
    return result
