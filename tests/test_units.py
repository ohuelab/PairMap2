"""Unit tests for pairmap2 components (fast, no LOMAP calls)."""
import numpy as np
import pytest
import networkx as nx

from rdkit import Chem
from rdkit.Chem import AllChem

from pairmap2.score_cache import ScoreCache
from pairmap2.score_engine import ScoreEngine
from pairmap2.path_solver import find_optimal_path, get_reachable_subgraph, get_cycled_edges
from pairmap2.mol_identity import canonical_smiles, normalize_mol


# ---------------------------------------------------------------------------
# ScoreCache
# ---------------------------------------------------------------------------

class TestScoreCache:
    def test_put_and_get(self):
        cache = ScoreCache()
        cache.put("CC", "CCC", 0.75)
        assert cache.get("CC", "CCC") == pytest.approx(0.75)

    def test_canonical_key_order(self):
        cache = ScoreCache()
        cache.put("ZZZ", "AAA", 0.5)
        # reverse order lookup should hit the same entry
        assert cache.get("AAA", "ZZZ") == pytest.approx(0.5)

    def test_missing_returns_none(self):
        cache = ScoreCache()
        assert cache.get("CC", "CCC") is None

    def test_options_hash_isolation(self):
        cache = ScoreCache()
        opts_a = {"time": 20}
        opts_b = {"time": 30}
        cache.put("CC", "CCC", 0.8, opts_a)
        cache.put("CC", "CCC", 0.4, opts_b)
        assert cache.get("CC", "CCC", opts_a) == pytest.approx(0.8)
        assert cache.get("CC", "CCC", opts_b) == pytest.approx(0.4)

    def test_get_matrix_complete(self):
        cache = ScoreCache()
        smiles = ["CC", "CCC", "CCCC"]
        cache.put("CC", "CCC", 0.9)
        cache.put("CC", "CCCC", 0.7)
        cache.put("CCC", "CCCC", 0.8)
        mat = cache.get_matrix(smiles)
        assert mat is not None
        assert mat[0][1] == pytest.approx(0.9)
        assert mat[1][0] == pytest.approx(0.9)
        assert mat[0][2] == pytest.approx(0.7)
        assert mat[1][2] == pytest.approx(0.8)

    def test_get_matrix_incomplete_returns_none(self):
        cache = ScoreCache()
        smiles = ["CC", "CCC", "CCCC"]
        cache.put("CC", "CCC", 0.9)
        # Only 1 of 3 pairs cached
        assert cache.get_matrix(smiles) is None


# ---------------------------------------------------------------------------
# MolIdentity
# ---------------------------------------------------------------------------

class TestMolIdentity:
    def test_canonical_smiles(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        smi = canonical_smiles(mol)
        assert isinstance(smi, str)
        assert len(smi) > 0

    def test_normalize_mol_removes_hs(self):
        mol = Chem.MolFromSmiles("CC")
        mol_h = Chem.AddHs(mol)
        norm = normalize_mol(mol_h)
        assert norm.GetNumAtoms() == 2  # only heavy atoms


# ---------------------------------------------------------------------------
# ScoreEngine pre-filters
# ---------------------------------------------------------------------------

class TestScoreEnginePrefilter:
    def _make_mol(self, smiles):
        return Chem.MolFromSmiles(smiles)

    def test_atom_count_prefilter(self):
        # ethane (2 heavy atoms) vs large molecule (>12 heavy atoms difference)
        mol_small = self._make_mol("CC")
        mol_large = self._make_mol("CCCCCCCCCCCCCCC")  # 15 C
        engine = ScoreEngine(atom_count_diff_threshold=10, tanimoto_prefilter=0.0)
        score = engine.get_score(mol_small, mol_large)
        assert score == 0.0

    def test_tanimoto_prefilter(self):
        # benzene vs a very dissimilar molecule (e.g. long alkane with heteroatoms)
        mol_a = self._make_mol("c1ccccc1")
        mol_b = self._make_mol("CCCCCCCCCC")  # long alkane - dissimilar to benzene
        engine = ScoreEngine(atom_count_diff_threshold=100, tanimoto_prefilter=0.9)
        score = engine.get_score(mol_a, mol_b)
        assert score == 0.0

    def test_cache_hit_skips_lomap(self):
        mol_a = self._make_mol("CC")
        mol_b = self._make_mol("CCC")
        cache = ScoreCache()
        smi_a = canonical_smiles(mol_a)
        smi_b = canonical_smiles(mol_b)
        cache.put(smi_a, smi_b, 0.999)
        engine = ScoreEngine(cache=cache)
        # Should return cached value without calling LOMAP
        assert engine.get_score(mol_a, mol_b) == pytest.approx(0.999)


# ---------------------------------------------------------------------------
# PathSolver
# ---------------------------------------------------------------------------

def _build_test_graph():
    """
    Build a small graph:
      0 --0.8-- 1 --0.7-- 2 --0.6-- 3
      0 --0.3-- 3
      1 --0.9-- 3
    Nodes 0 and 3 are source/target.
    """
    g = nx.Graph()
    g.add_edge(0, 1, score=0.8)
    g.add_edge(1, 2, score=0.7)
    g.add_edge(2, 3, score=0.6)
    g.add_edge(0, 3, score=0.3)
    g.add_edge(1, 3, score=0.9)
    return g


class TestPathSolver:
    def test_find_optimal_path_basic(self):
        g = _build_test_graph()
        path = find_optimal_path(g, source=0, target=3, max_path_length=4)
        assert path[0] == 0
        assert path[-1] == 3
        assert len(path) >= 2

    def test_find_optimal_path_respects_max_length(self):
        g = _build_test_graph()
        path = find_optimal_path(g, source=0, target=3, max_path_length=2)
        assert len(path) - 1 <= 2

    def test_find_optimal_path_no_path_raises(self):
        g = nx.Graph()
        g.add_edge(0, 1, score=0.5)
        g.add_node(2)
        with pytest.raises((ValueError, Exception)):
            find_optimal_path(g, source=0, target=2, max_path_length=4)

    def test_get_reachable_subgraph(self):
        g = _build_test_graph()
        sub = get_reachable_subgraph(g, source=0, target=3, max_path_length=4)
        # All nodes should be reachable within length 4 (0→3 directly or via intermediates)
        assert 0 in sub.nodes
        assert 3 in sub.nodes

    def test_get_reachable_subgraph_excludes_too_distant(self):
        g = nx.path_graph(10)  # 0-1-2-...-9
        for u, v in g.edges():
            g[u][v]['score'] = 0.5
        sub = get_reachable_subgraph(g, source=0, target=9, max_path_length=3)
        # Path 0->...->9 requires 9 hops, which exceeds max_path_length=3.
        # No node satisfies dist(0,v)+dist(v,9)<=3, so subgraph is empty.
        assert len(sub.nodes) == 0

    def test_get_reachable_subgraph_short_path(self):
        g = nx.path_graph(5)  # 0-1-2-3-4
        for u, v in g.edges():
            g[u][v]['score'] = 0.5
        sub = get_reachable_subgraph(g, source=0, target=4, max_path_length=4)
        # The only path is 0-1-2-3-4 (length 4), all nodes reachable
        assert set(sub.nodes) == {0, 1, 2, 3, 4}

    def test_get_cycled_edges(self):
        # Triangle: 0-1-2-0, edge (0,2) should be "cycled" (alt path 0-1-2)
        g = nx.Graph()
        g.add_edge(0, 1, score=0.8)
        g.add_edge(1, 2, score=0.7)
        g.add_edge(0, 2, score=0.6)
        cycled = get_cycled_edges(g, cycle_links=[(0, 2)], cycle_length=3)
        assert (0, 2) in cycled

    def test_get_cycled_edges_no_alternative(self):
        # Linear: 0-1-2, edge (0,1) has no alternative path to 1
        g = nx.Graph()
        g.add_edge(0, 1, score=0.8)
        g.add_edge(1, 2, score=0.7)
        cycled = get_cycled_edges(g, cycle_links=[(0, 1)], cycle_length=3)
        assert (0, 1) not in cycled
