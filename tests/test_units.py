"""Unit tests for pairmap2 components (fast, no LOMAP calls)."""
import numpy as np
import pytest
import networkx as nx

from rdkit import Chem
from rdkit.Chem import AllChem

from pairmap2.score_cache import ScoreCache
from pairmap2.score_engine import ScoreEngine
from pairmap2.mol_identity import canonical_smiles


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


