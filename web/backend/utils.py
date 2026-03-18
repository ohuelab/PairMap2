from __future__ import annotations

from typing import Any

import networkx as nx


def _compute_aligned_svgs(node_mols: dict) -> dict:
    """Compute global-MCS-aligned 2D SVGs for all node molecules.

    Returns a dict mapping node_id → SVG string.  Uses rdFMCS to find a
    common scaffold and aligns all 2D depictions to that scaffold.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdFMCS
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        return {}

    # Generate clean 2D structures from canonical SMILES
    mols_2d: dict = {}
    for nid, mol in node_mols.items():
        if mol is None:
            continue
        try:
            smi = Chem.MolToSmiles(mol)
            m = Chem.RWMol(Chem.MolFromSmiles(smi))
            AllChem.Compute2DCoords(m)
            mols_2d[nid] = m
        except Exception:
            pass

    mol_list = list(mols_2d.values())
    if len(mol_list) >= 2:
        try:
            mcs_res = rdFMCS.FindMCS(
                mol_list,
                atomCompare=rdFMCS.AtomCompare.CompareAny,
                timeout=10,
            )
            if mcs_res.queryMol is not None and mcs_res.queryMol.GetNumAtoms() >= 3:
                mcs_mol = mcs_res.queryMol
                AllChem.Compute2DCoords(mcs_mol)
                for m in mol_list:
                    try:
                        AllChem.GenerateDepictionMatching2DStructure(m, mcs_mol)
                    except Exception:
                        pass
        except Exception:
            pass

    svgs: dict = {}
    for nid, m in mols_2d.items():
        try:
            drawer = rdMolDraw2D.MolDraw2DSVG(100, 80)
            drawer.DrawMolecule(m)
            drawer.FinishDrawing()
            svgs[nid] = drawer.GetDrawingText()
        except Exception:
            pass
    return svgs


def graph_to_cytoscape(
    graph: nx.Graph,
    node_mols: dict,
    source_idx: int | None = None,
    target_idx: int | None = None,
) -> dict:
    """Convert an nx.Graph + mol dict to a Cytoscape.js elements list.

    Args:
        graph: NetworkX graph from build_map() or IntermediateGraphManager.
        node_mols: Dict mapping node id → RDKit Mol.
        source_idx: Node id of the source molecule (Pair mode).
        target_idx: Node id of the target molecule (Pair mode).
    """
    try:
        from rdkit.Chem import MolToSmiles, MolToMolBlock
    except ImportError:
        MolToSmiles = None
        MolToMolBlock = None

    # Pre-compute MCS-aligned SVGs for all nodes
    aligned_svgs = _compute_aligned_svgs(node_mols)

    elements: list[dict] = []

    for node_id, data in graph.nodes(data=True):
        mol = node_mols.get(node_id)
        label = data.get("label") or data.get("fname") or data.get("NAME") or str(node_id)
        if mol and mol.HasProp("_Name"):
            label = mol.GetProp("_Name")

        if source_idx is not None and target_idx is not None:
            is_source = node_id == source_idx
            is_target = node_id == target_idx
            active = is_source or is_target
        else:
            active = data.get("active", True)
            is_source = False
            is_target = False

        node_data: dict[str, Any] = {
            "id": str(node_id),
            "label": label,
            "active": active,
            "intermediate": not active,
            "is_source": is_source,
            "is_target": is_target,
            "aligned_svg": aligned_svgs.get(node_id, ""),
        }
        if mol and MolToSmiles:
            try:
                node_data["smiles"] = MolToSmiles(mol)
            except Exception:
                pass
        if mol and MolToMolBlock:
            try:
                node_data["molblock"] = MolToMolBlock(mol)
            except Exception:
                pass

        elements.append({"data": node_data, "group": "nodes"})

    for u, v, edata in graph.edges(data=True):
        score = edata.get("similarity", edata.get("score", 0.0))
        elements.append({
            "data": {
                "id": f"{u}-{v}",
                "source": str(u),
                "target": str(v),
                "similarity": round(float(score), 4),
                "bad_edge": edata.get("bad_edge", False),
            },
            "group": "edges",
        })

    return {"elements": elements}
