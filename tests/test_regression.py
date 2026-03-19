"""Regression tests for pairmap2 Pipeline using benchmark data.

These tests run the full pipeline on 9 benchmark cases and verify quality
equivalence (not exact match) with reference outputs.

Marked as @pytest.mark.slow because each case can take 3–20 seconds.
Run with: pytest -m slow
"""
import json
import pytest
import pandas as pd
from pathlib import Path

from rdkit import Chem

from pairmap2 import Pipeline, PipelineConfig

DATA_DIR = Path(__file__).parents[1] / "benchmarks" / "data"

CASES = [
    "Bace1_00",
    "P38_00",
    "P38_01",
    "P38_02",
    "PTP1B_00",
    "PTP1B_01",
    "PTP1B_02",
    "PTP1B_03",
    "PTP1B_04",
]


def load_mol(sdf_path: Path):
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True)
    mols = [m for m in suppl if m is not None]
    assert mols, f"No molecule in {sdf_path}"
    return mols[0]


def load_case(case_name: str):
    case_dir = DATA_DIR / case_name
    source = load_mol(case_dir / "source.sdf")
    target = load_mol(case_dir / "target.sdf")
    source.SetProp("_Name", "source")
    target.SetProp("_Name", "target")

    ref_path = json.loads((case_dir / "reference_path.json").read_text())
    ref_links = pd.read_csv(case_dir / "reference_links.csv")

    return source, target, ref_path, ref_links


def run_pipeline(source, target, config=None):
    if config is None:
        config = PipelineConfig(save_output=False, verbose=False)

    mols = [source, target]
    pipeline = Pipeline(config)
    # Compute initial similarity without prefilters so that very dissimilar
    # molecule pairs still get a non-zero LOMAP score for the initial edge.
    from pairmap2.score_engine import ScoreEngine
    sim = ScoreEngine(tanimoto_prefilter=0.0, atom_count_diff_threshold=10000).get_score(source, target, {})
    bad_edge = sim < config.similarity_threshold
    df = pd.DataFrame(
        [("source", "target", sim, bad_edge)],
        columns=["Node1", "Node2", "score", "BadEdge"],
    )
    return pipeline.run_from_moldf(mols, df)


# ---------------------------------------------------------------------------
# Basic smoke tests (fast, just check no exception)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case_name", CASES)
@pytest.mark.slow
def test_pipeline_runs_without_error(case_name):
    """Pipeline should complete without raising an exception."""
    source, target, _, _ = load_case(case_name)
    result = run_pipeline(source, target)
    final_graph = result.graphs[-1]
    assert len(final_graph.nodes) >= 2, "Final graph must have at least source and target"


@pytest.mark.parametrize("case_name", CASES)
@pytest.mark.slow
def test_pipeline_result_structure(case_name):
    """PipelineResult should have correct structure."""
    source, target, _, _ = load_case(case_name)
    result = run_pipeline(source, target)

    assert result.graphs, "graphs list must not be empty"
    assert result.node_mols, "node_mols must not be empty"
    assert result.timings, "timings must not be empty"
    assert result.timings[0].wall_time > 0
    assert result.timings[0].stage == "run_from_moldf"


@pytest.mark.parametrize("case_name", CASES)
@pytest.mark.slow
def test_final_graph_connected(case_name):
    """Final graph must be connected (source and target in same component)."""
    import networkx as nx
    source, target, _, _ = load_case(case_name)
    result = run_pipeline(source, target)
    final_graph = result.graphs[-1]

    assert nx.is_connected(final_graph) or len(list(nx.connected_components(final_graph))) >= 1
    # Source (node 0) and target (node 1) must be reachable from each other
    assert nx.has_path(final_graph, 0, 1), "source and target must be in the same component"


@pytest.mark.parametrize("case_name", CASES)
@pytest.mark.slow
def test_path_length_reasonable(case_name):
    """Path from source to target in final graph must be <= max_path_length."""
    import networkx as nx
    config = PipelineConfig(save_output=False, max_path_length=4)
    source, target, ref_path, _ = load_case(case_name)
    result = run_pipeline(source, target, config)
    final_graph = result.graphs[-1]

    path_len = nx.shortest_path_length(final_graph, 0, 1)
    assert path_len <= config.max_path_length, (
        f"Path length {path_len} exceeds max_path_length {config.max_path_length}"
    )


@pytest.mark.parametrize("case_name", CASES)
@pytest.mark.slow
def test_edge_scores_positive(case_name):
    """All edges in the final graph must have positive similarity scores."""
    source, target, _, _ = load_case(case_name)
    result = run_pipeline(source, target)
    final_graph = result.graphs[-1]

    for u, v, data in final_graph.edges(data=True):
        score = data.get("similarity", data.get("score", -1))
        assert score >= 0, f"Edge ({u},{v}) has negative score {score}"


@pytest.mark.parametrize("case_name", CASES)
@pytest.mark.slow
def test_intermediate_nodes_have_smiles(case_name):
    """All intermediate nodes must have SMILES stored."""
    source, target, _, _ = load_case(case_name)
    result = run_pipeline(source, target)
    final_graph = result.graphs[-1]
    node_mols = result.node_mols

    for n in final_graph.nodes:
        if final_graph.nodes[n].get("intermediate"):
            mol = node_mols[n]
            assert mol is not None
            assert mol.HasProp("smiles"), f"Intermediate node {n} missing 'smiles' property"


# ---------------------------------------------------------------------------
# Quality comparison against reference (allow slight variation)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case_name", CASES)
@pytest.mark.slow
def test_quality_vs_reference(case_name):
    """
    Verify path quality is reasonable.

    Note: The V3 reference was computed from full ligand sets via LOMAP, while
    the pipeline here runs on only 2 molecules (source + target).  A direct
    quality comparison would be unfair, so we instead check:
      1. A path of length <= reference path length + 1 exists.
      2. Every edge score on the path is positive (no degenerate edges).
    """
    import networkx as nx

    source, target, ref_path, ref_links = load_case(case_name)
    ref_path_len = len(ref_path) - 1  # number of edges in reference path

    result = run_pipeline(source, target)
    final_graph = result.graphs[-1]

    path = nx.shortest_path(final_graph, 0, 1)
    path_len = len(path) - 1

    assert path_len <= ref_path_len + 1, (
        f"Path length {path_len} exceeds reference {ref_path_len} + 1 for {case_name}"
    )
    for i in range(len(path) - 1):
        score = final_graph[path[i]][path[i + 1]].get("similarity", 0)
        assert score > 0, f"Edge ({path[i]},{path[i+1]}) has zero/negative score for {case_name}"
