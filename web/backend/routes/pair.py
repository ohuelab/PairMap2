from __future__ import annotations

import asyncio
import json
from functools import partial
from io import BytesIO

from fastapi import APIRouter, Form, HTTPException, UploadFile
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field

from ..models import MapGenConfig, PairRequest, SearchConfig
from ..utils import graph_to_cytoscape
from .. import pair_cache as cache
from ..pair_cache import CacheEntry

router = APIRouter()


def _embed_mol(mol):
    """Add a minimal 3D conformer so lomap.mcs.MCS can run (requires conformer)."""
    from rdkit import Chem as _Chem
    from rdkit.Chem import AllChem
    mol = _Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    mol = _Chem.RemoveHs(mol)
    return mol


def _run_search(mol_a, mol_b, name_a: str, name_b: str, search: SearchConfig, engine: str = "v1") -> str:
    """Run intermediate search, embed conformers, cache results, return session_id."""
    from pairmap import SearchIntermediates

    s = search
    si = SearchIntermediates(
        mol_a, mol_b,
        is_atom_modfication_enabled=s.is_atom_modfication_enabled,
        cap_ring_with_carbon=s.cap_ring_with_carbon,
        cap_ring_with_hydrogen=s.cap_ring_with_hydrogen,
        no_backward_search=s.no_backward_search,
        use_seed=s.use_seed,
        max_intermediate=s.max_intermediate,
        ionize=s.ionize,
    )
    si.search()

    intermediates = si.intermediates
    intermediates_all = getattr(si, 'intermediates_all', intermediates)

    # Ensure all intermediates have 3D conformers (for 3D display and scoring)
    for i, mol in enumerate(intermediates):
        if mol is not None and mol.GetNumConformers() == 0:
            intermediates[i] = _embed_mol(mol)

    if len(intermediates) >= 2:
        intermediates[0].SetProp("_Name", name_a)
        intermediates[1].SetProp("_Name", name_b)

    node_mols = {i: mol for i, mol in enumerate(intermediates)}

    entry = CacheEntry(
        intermediates=intermediates,
        intermediates_all=intermediates_all,
        node_mols=node_mols,
        graph=None,
        score_matrix=None,
        name_a=name_a,
        name_b=name_b,
        engine=engine,
    )
    return cache.store(entry)


def _run_map(session_id: str, mapgen: MapGenConfig) -> dict:
    """Build perturbation map from cached intermediates, return cytoscape data."""
    from pairmap import MapGenerator

    entry = cache.get(session_id)
    if entry is None:
        raise ValueError(f"Session {session_id!r} not found or expired")

    m = mapgen
    engine = getattr(entry, 'engine', 'v1')

    if engine == "v2":
        if entry.score_matrix is None:
            from pairmap2.score_engine import ScoreEngine
            from pairmap2.score_cache import ScoreCache
            # jobs=1: avoid multiprocessing Pool spawn overhead in web context (macOS spawn)
            se = ScoreEngine(cache=ScoreCache(None), jobs=1)
            score_matrix = se.get_score_matrix(entry.intermediates, {})
        else:
            score_matrix = entry.score_matrix
        mg = MapGenerator(
            entry.intermediates,
            custom_score_matrix=score_matrix,
            maxOptimalPathLength=m.maxOptimalPathLength,
            roughScoreThreshold=m.roughScoreThreshold,
            minScoreThreshold=m.minScoreThreshold,
            optimal_path_mode=m.optimal_path_mode,
            CycleLinkThreshold=m.CycleLinkThreshold,
            squared_sum=m.squared_sum,
            source_node_index=0,
            target_node_index=1,
        )
    else:
        mg = MapGenerator(
            entry.intermediates,
            custom_score_matrix=entry.score_matrix,
            maxOptimalPathLength=m.maxOptimalPathLength,
            roughScoreThreshold=m.roughScoreThreshold,
            minScoreThreshold=m.minScoreThreshold,
            optimal_path_mode=m.optimal_path_mode,
            CycleLinkThreshold=m.CycleLinkThreshold,
            squared_sum=m.squared_sum,
            source_node_index=0,
            target_node_index=1,
        )

    graph = mg.build_map()
    entry.score_matrix = mg.score_matrix
    entry.graph = graph

    result = graph_to_cytoscape(graph, entry.node_mols, source_idx=0, target_idx=1)
    result["n_intermediates"] = max(0, len(entry.intermediates) - 2)
    result["session_id"] = session_id
    return result


def _run_pair_smiles(req: PairRequest) -> dict:
    from rdkit import Chem

    mol_a = Chem.MolFromSmiles(req.smiles_a)
    mol_b = Chem.MolFromSmiles(req.smiles_b)
    if mol_a is None:
        raise ValueError(f"Invalid SMILES for molecule A: {req.smiles_a!r}")
    if mol_b is None:
        raise ValueError(f"Invalid SMILES for molecule B: {req.smiles_b!r}")

    mol_a = _embed_mol(mol_a)
    mol_b = _embed_mol(mol_b)

    session_id = _run_search(mol_a, mol_b, req.name_a, req.name_b, req.search, engine=req.engine)
    return _run_map(session_id, req.mapgen)


def _run_pair_sdf(
    file_a_bytes: bytes, file_b_bytes: bytes,
    search: SearchConfig, mapgen: MapGenConfig,
    engine: str = "v1",
) -> dict:
    from rdkit.Chem import ForwardSDMolSupplier

    def read_first_mol(data: bytes, label: str):
        suppl = ForwardSDMolSupplier(BytesIO(data), removeHs=False)
        mol = next(suppl, None)
        if mol is None:
            raise ValueError(f"Could not read molecule from SDF {label}")
        if mol.GetNumConformers() == 0:
            mol = _embed_mol(mol)
        return mol

    mol_a = read_first_mol(file_a_bytes, "A")
    mol_b = read_first_mol(file_b_bytes, "B")

    name_a = (mol_a.GetProp("_Name") if mol_a.HasProp("_Name") else "").strip() or "Molecule A"
    name_b = (mol_b.GetProp("_Name") if mol_b.HasProp("_Name") else "").strip() or "Molecule B"

    session_id = _run_search(mol_a, mol_b, name_a, name_b, search, engine=engine)
    return _run_map(session_id, mapgen)


# ── Pydantic model for remap ──────────────────────────────────────────────────

class RemapRequest(BaseModel):
    session_id: str
    mapgen: MapGenConfig = Field(default_factory=MapGenConfig)


# ── Endpoints ─────────────────────────────────────────────────────────────────

@router.post("/pair")
async def run_pair(req: PairRequest):
    """Run PairMap synchronously (SMILES input)."""
    loop = asyncio.get_event_loop()
    try:
        result = await loop.run_in_executor(None, partial(_run_pair_smiles, req))
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except Exception:
        import traceback
        raise HTTPException(status_code=500, detail=traceback.format_exc())
    return result


@router.post("/pair/sdf")
async def run_pair_sdf(
    file_a: UploadFile,
    file_b: UploadFile,
    search: str = Form("{}"),
    mapgen: str = Form("{}"),
    engine: str = Form("v1"),
):
    """Run PairMap using SDF files that already contain 3D coordinates."""
    try:
        bytes_a = await file_a.read()
        bytes_b = await file_b.read()
        search_cfg = SearchConfig(**json.loads(search))
        mapgen_cfg = MapGenConfig(**json.loads(mapgen))
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc))

    loop = asyncio.get_event_loop()
    try:
        result = await loop.run_in_executor(
            None, partial(_run_pair_sdf, bytes_a, bytes_b, search_cfg, mapgen_cfg, engine)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except Exception:
        import traceback
        raise HTTPException(status_code=500, detail=traceback.format_exc())
    return result


@router.post("/pair/remap")
async def remap_pair(req: RemapRequest):
    """Re-run map generation only from cached intermediates (fast, no re-search)."""
    loop = asyncio.get_event_loop()
    try:
        result = await loop.run_in_executor(None, partial(_run_map, req.session_id, req.mapgen))
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    except Exception:
        import traceback
        raise HTTPException(status_code=500, detail=traceback.format_exc())
    return result


# ── MCS highlight endpoint ─────────────────────────────────────────────────────

def _compute_mcs_highlight(session_id: str, node_a: int, node_b: int) -> dict:
    """Compute MCS-based highlighted SVGs for two nodes."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Geometry import rdGeometry

    entry = cache.get(session_id)
    if entry is None:
        raise ValueError(f"Session {session_id!r} not found or expired")

    mol_a_raw = entry.node_mols.get(node_a)
    mol_b_raw = entry.node_mols.get(node_b)
    if mol_a_raw is None or mol_b_raw is None:
        raise ValueError(f"Node not found: {node_a} or {node_b}")

    # Create clean 2D structures
    def to_2d(mol):
        m = Chem.RWMol(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
        AllChem.Compute2DCoords(m)
        return m

    mol_a_2d = to_2d(mol_a_raw)
    mol_b_2d = to_2d(mol_b_raw)

    # Find MCS
    mcs_result = rdFMCS.FindMCS(
        [mol_a_2d, mol_b_2d],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        timeout=10,
    )

    mcs_map: list = []
    if mcs_result.numAtoms > 0 and mcs_result.smartsString:
        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        if mcs_mol is not None:
            match_a = mol_a_2d.GetSubstructMatch(mcs_mol)
            match_b = mol_b_2d.GetSubstructMatch(mcs_mol)
            if match_a and match_b:
                mcs_map = list(zip(match_a, match_b))
                # Align mol_b to mol_a using MCS coordinate map
                if len(mcs_map) >= 3:
                    try:
                        conf_a = mol_a_2d.GetConformer()
                        coord_map = {}
                        for idx_a, idx_b in mcs_map:
                            pos = conf_a.GetAtomPosition(idx_a)
                            coord_map[idx_b] = rdGeometry.Point2D(pos.x, pos.y)
                        AllChem.Compute2DCoords(mol_b_2d, coordMap=coord_map)
                    except Exception:
                        pass

    mcs_set_a = {ia for ia, _ in mcs_map}
    mcs_set_b = {ib for _, ib in mcs_map}
    all_heavy_a = {a.GetIdx() for a in mol_a_2d.GetAtoms() if a.GetAtomicNum() > 0}
    all_heavy_b = {a.GetIdx() for a in mol_b_2d.GetAtoms() if a.GetAtomicNum() > 0}
    deleted = all_heavy_a - mcs_set_a
    inserted = all_heavy_b - mcs_set_b

    COLOR_GREEN = (0.0, 0.7, 0.0)
    COLOR_RED   = (0.8, 0.1, 0.1)
    COLOR_BLUE  = (0.1, 0.4, 0.85)

    def draw_mol(mol, common_set, other_set, other_color):
        atoms = list(common_set) + list(other_set)
        colors = {idx: COLOR_GREEN for idx in common_set}
        colors.update({idx: other_color for idx in other_set})
        drawer = rdMolDraw2D.MolDraw2DSVG(200, 160)
        drawer.drawOptions().addStereoAnnotation = False
        drawer.DrawMolecule(
            mol,
            highlightAtoms=atoms,
            highlightAtomColors=colors,
            highlightBonds=[],
            highlightBondColors={},
        )
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    svg_a = draw_mol(mol_a_2d, mcs_set_a, deleted, COLOR_RED)
    svg_b = draw_mol(mol_b_2d, mcs_set_b, inserted, COLOR_BLUE)

    return {
        "svg_a": svg_a,
        "svg_b": svg_b,
        "mcs_map": mcs_map,
        "n_common": len(mcs_map),
        "n_deleted": len(deleted),
        "n_inserted": len(inserted),
    }


@router.get("/pair/{session_id}/mcs/{node_a}/{node_b}")
async def get_mcs_highlight(session_id: str, node_a: int, node_b: int):
    """Return MCS-highlighted SVGs for two nodes in a pair session."""
    loop = asyncio.get_event_loop()
    try:
        result = await loop.run_in_executor(
            None, partial(_compute_mcs_highlight, session_id, node_a, node_b)
        )
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    except Exception:
        import traceback
        raise HTTPException(status_code=500, detail=traceback.format_exc())
    return result


# ── Download endpoints ────────────────────────────────────────────────────────

def _mols_to_sdf(mols) -> str:
    """Convert a list of RDKit Mols to SDF-format string."""
    from rdkit import Chem
    parts = []
    for mol in mols:
        if mol is None:
            continue
        try:
            parts.append(Chem.MolToMolBlock(mol))
        except Exception:
            pass
    return "".join(f"{mb}$$$$\n" for mb in parts)


@router.get("/pair/{session_id}/download/intermediates.sdf")
async def download_intermediates_sdf(session_id: str):
    """Download SDF of all nodes present in the perturbation map."""
    entry = cache.get(session_id)
    if entry is None:
        raise HTTPException(status_code=404, detail="Session not found or expired")
    if entry.graph is None:
        raise HTTPException(status_code=400, detail="Map has not been generated yet")

    mols = [entry.node_mols.get(nid) for nid in entry.graph.nodes()]
    content = _mols_to_sdf(mols)
    return StreamingResponse(
        iter([content]),
        media_type="chemical/x-mdl-sdfile",
        headers={"Content-Disposition": "attachment; filename=intermediates.sdf"},
    )


@router.get("/pair/{session_id}/download/links.csv")
async def download_links_csv(session_id: str):
    """Download CSV of perturbation map edges (name_a, name_b, score)."""
    entry = cache.get(session_id)
    if entry is None:
        raise HTTPException(status_code=404, detail="Session not found or expired")
    if entry.graph is None:
        raise HTTPException(status_code=400, detail="Map has not been generated yet")

    def _mol_name(mol, fallback: str) -> str:
        if mol and mol.HasProp("_Name"):
            name = mol.GetProp("_Name").strip()
            if name:
                return name
        return fallback

    lines = ["name_a,name_b,score"]
    for u, v, edata in entry.graph.edges(data=True):
        name_u = _mol_name(entry.node_mols.get(u), str(u))
        name_v = _mol_name(entry.node_mols.get(v), str(v))
        score = edata.get("similarity", edata.get("score", 0.0))
        lines.append(f"{name_u},{name_v},{score:.4f}")

    content = "\n".join(lines) + "\n"
    return StreamingResponse(
        iter([content]),
        media_type="text/csv",
        headers={"Content-Disposition": "attachment; filename=links.csv"},
    )


@router.get("/pair/{session_id}/download/all_intermediates.sdf")
async def download_all_intermediates_sdf(session_id: str):
    """Download SDF of all intermediates found during search (including non-map nodes)."""
    entry = cache.get(session_id)
    if entry is None:
        raise HTTPException(status_code=404, detail="Session not found or expired")

    content = _mols_to_sdf(entry.intermediates_all)
    return StreamingResponse(
        iter([content]),
        media_type="chemical/x-mdl-sdfile",
        headers={"Content-Disposition": "attachment; filename=all_intermediates.sdf"},
    )
