"""Map Mode v1 routes — single SDF → PairMap engine → perturbation map."""
from __future__ import annotations

import asyncio
import json
import uuid
from functools import partial
from pathlib import Path

from fastapi import APIRouter, File, Form, Header, HTTPException, UploadFile
from fastapi.responses import FileResponse

from .. import map_store
from ..map_worker import MAP_JOBS_DIR, submit_job
from ..models import MapJobList, MapJobStatus

router = APIRouter()

_processes: dict[str, object] = {}


@router.post("/jobs", response_model=MapJobStatus, status_code=202)
async def create_map_job(
    file: UploadFile = File(..., description="Input SDF with ligands"),
    engine: str = Form("v2"),
    config: str = Form("{}"),
    x_session_id: str = Header(...),
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

    status = map_store.create_job(job_id, engine, cfg, session_id=x_session_id)
    proc = submit_job(job_id, engine, cfg, str(sdf_path))
    _processes[job_id] = proc
    return status


@router.get("/jobs", response_model=MapJobList)
async def list_map_jobs(x_session_id: str = Header(...)):
    return MapJobList(jobs=map_store.list_jobs(session_id=x_session_id))


@router.get("/jobs/{job_id}", response_model=MapJobStatus)
async def get_map_job(job_id: str, x_session_id: str = Header(...)):
    status = map_store.get_job(job_id, session_id=x_session_id)
    if status is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return status


@router.get("/jobs/{job_id}/graph")
async def get_map_graph(job_id: str, x_session_id: str = Header(...)):
    if map_store.get_job(job_id, session_id=x_session_id) is None:
        raise HTTPException(status_code=404, detail="Job not found")
    graph_path = MAP_JOBS_DIR / job_id / "graph.json"
    if not graph_path.exists():
        status = map_store.get_job(job_id, session_id=x_session_id)
        raise HTTPException(
            status_code=404,
            detail=f"Graph not yet available (status: {status.status})",
        )
    with open(graph_path) as f:
        return json.load(f)


@router.get("/jobs/{job_id}/artifacts/{filename}")
async def get_map_artifact(job_id: str, filename: str, x_session_id: str = Header(...)):
    if map_store.get_job(job_id, session_id=x_session_id) is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if ".." in filename or "/" in filename or "\\" in filename:
        raise HTTPException(status_code=400, detail="Invalid filename")
    artifact = MAP_JOBS_DIR / job_id / "output" / filename
    if not artifact.exists():
        raise HTTPException(status_code=404, detail="Artifact not found")
    return FileResponse(str(artifact))


@router.post("/jobs/{job_id}/cancel", response_model=MapJobStatus)
async def cancel_map_job(job_id: str, x_session_id: str = Header(...)):
    status = map_store.get_job(job_id, session_id=x_session_id)
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
    return map_store.get_job(job_id, session_id=x_session_id)


def _compute_map_mcs(job_id: str, node_a: str, node_b: str) -> dict:
    """Load molblocks from graph.json and compute MCS-highlighted SVGs."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Geometry import rdGeometry

    graph_path = MAP_JOBS_DIR / job_id / "graph.json"
    if not graph_path.exists():
        raise FileNotFoundError(f"graph.json not found for job {job_id}")

    with open(graph_path) as f:
        graph = json.load(f)

    nodes = {e["data"]["id"]: e["data"] for e in graph["elements"] if e["group"] == "nodes"}
    data_a = nodes.get(node_a)
    data_b = nodes.get(node_b)
    if data_a is None or data_b is None:
        raise ValueError(f"Node not found: {node_a!r} or {node_b!r}")

    def to_2d(molblock_or_smiles):
        mol = Chem.MolFromMolBlock(molblock_or_smiles) if '\n' in molblock_or_smiles else None
        if mol is None:
            mol = Chem.MolFromSmiles(molblock_or_smiles)
        if mol is None:
            raise ValueError("Could not parse molecule")
        m = Chem.RWMol(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
        AllChem.Compute2DCoords(m)
        return m

    mol_a = to_2d(data_a.get("molblock") or data_a.get("smiles", ""))
    mol_b = to_2d(data_b.get("molblock") or data_b.get("smiles", ""))

    mcs_result = rdFMCS.FindMCS(
        [mol_a, mol_b],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        timeout=10,
    )

    mcs_map: list = []
    if mcs_result.numAtoms > 0 and mcs_result.smartsString:
        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        if mcs_mol is not None:
            match_a = mol_a.GetSubstructMatch(mcs_mol)
            match_b = mol_b.GetSubstructMatch(mcs_mol)
            if match_a and match_b:
                mcs_map = list(zip(match_a, match_b))
                if len(mcs_map) >= 3:
                    try:
                        conf_a = mol_a.GetConformer()
                        coord_map = {}
                        for idx_a, idx_b in mcs_map:
                            pos = conf_a.GetAtomPosition(idx_a)
                            coord_map[idx_b] = rdGeometry.Point2D(pos.x, pos.y)
                        AllChem.Compute2DCoords(mol_b, coordMap=coord_map)
                    except Exception:
                        pass

    mcs_set_a = {ia for ia, _ in mcs_map}
    mcs_set_b = {ib for _, ib in mcs_map}
    deleted  = {a.GetIdx() for a in mol_a.GetAtoms() if a.GetAtomicNum() > 0} - mcs_set_a
    inserted = {a.GetIdx() for a in mol_b.GetAtoms() if a.GetAtomicNum() > 0} - mcs_set_b

    COLOR_GREEN = (0.0, 0.7, 0.0)
    COLOR_RED   = (0.8, 0.1, 0.1)
    COLOR_BLUE  = (0.1, 0.4, 0.85)

    def draw_mol(mol, common_set, other_set, other_color):
        atoms  = list(common_set) + list(other_set)
        colors = {idx: COLOR_GREEN for idx in common_set}
        colors.update({idx: other_color for idx in other_set})
        drawer = rdMolDraw2D.MolDraw2DSVG(200, 160)
        drawer.drawOptions().addStereoAnnotation = False
        drawer.DrawMolecule(mol, highlightAtoms=atoms, highlightAtomColors=colors,
                            highlightBonds=[], highlightBondColors={})
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    return {
        "svg_a": draw_mol(mol_a, mcs_set_a, deleted, COLOR_RED),
        "svg_b": draw_mol(mol_b, mcs_set_b, inserted, COLOR_BLUE),
        "label_a": data_a.get("label", node_a),
        "label_b": data_b.get("label", node_b),
        "n_common": len(mcs_map),
        "n_deleted": len(deleted),
        "n_inserted": len(inserted),
    }


@router.get("/jobs/{job_id}/mcs/{node_a}/{node_b}")
async def get_map_mcs(job_id: str, node_a: str, node_b: str, x_session_id: str = Header(...)):
    """Return MCS-highlighted SVGs for two nodes in a map job graph."""
    status = map_store.get_job(job_id, session_id=x_session_id)
    if status is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if status.status != "completed":
        raise HTTPException(status_code=409, detail=f"Job is not completed (status: {status.status})")
    loop = asyncio.get_event_loop()
    try:
        result = await loop.run_in_executor(None, partial(_compute_map_mcs, job_id, node_a, node_b))
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    return result
