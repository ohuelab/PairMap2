"""Score engine: pre-filter + LOMAP MCS scoring with persistent caching."""
import copy
import itertools
from multiprocessing import Pool
from typing import Optional

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFMCS
from lomap.dbmol import ecr
from lomap import mcs as lomap_mcs

from .score_cache import ScoreCache
from .mol_identity import canonical_smiles



# Module-level functions so they can be pickled for multiprocessing Pool.
def _tanimoto(mol_a, mol_b) -> float:
    """Morgan fingerprint Tanimoto similarity (radius 2, 1024 bits)."""
    fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, 2, nBits=1024)
    fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, 2, nBits=1024)
    return DataStructs.TanimotoSimilarity(fp_a, fp_b)


def _lomap_score(mol_a, mol_b, options: dict) -> float:
    """Calculate LOMAP MCS score for a molecule pair.

    The score combines ECR with the standard set of LOMAP rules (mncar, mcsr,
    atomic_number_rule, hybridization_rule, sulfonamides_rule,
    heterocycles_rule, transmuting_methyl_into_ring_rule,
    transmuting_ring_sizes_rule).
    """
    mola = copy.deepcopy(mol_a)
    molb = copy.deepcopy(mol_b)
    ecr_score = ecr(mola, molb)
    MC = lomap_mcs.MCS(
        mola,
        molb,
        time=options.get("time", 20),
        verbose=options.get("verbose", "info"),
        threed=options.get("threed", False),
        max3d=options.get("max3d", 1000.0),
        element_change=options.get("element_change", True),
        seed=options.get("seed", ""),
        shift=options.get("shift", True),
    )
    MC.all_atom_match_list()
    score = ecr_score * MC.mncar() * MC.mcsr()
    score *= MC.atomic_number_rule() * MC.hybridization_rule()
    score *= MC.sulfonamides_rule() * MC.heterocycles_rule()
    score *= MC.transmuting_methyl_into_ring_rule() * MC.transmuting_ring_sizes_rule()
    return score


def _compute_score_worker(args):
    """Worker function for multiprocessing Pool.

    Receives mol objects directly (RDKit mols are picklable and preserve
    conformers, which LOMAP requires for its internal position remapping).
    """
    i, j, mol_i, mol_j, options = args
    score = _lomap_score(mol_i, mol_j, options)
    return i, j, score



class ScoreEngine:
    """Computes pairwise similarity scores with pre-filtering and caching.

    Pre-filters (atom-count difference and Tanimoto threshold) quickly rule out
    pairs that will score zero, avoiding expensive LOMAP MCS calls.  Results
    are stored in a ``ScoreCache`` so repeated queries are O(1).
    """

    def __init__(
        self,
        cache: Optional[ScoreCache] = None,
        atom_count_diff_threshold: int = 10,
        tanimoto_prefilter: float = 0.3,
        jobs: int = -1,
    ):
        self.cache = cache if cache is not None else ScoreCache()
        self.atom_count_diff_threshold = atom_count_diff_threshold
        self.tanimoto_prefilter = tanimoto_prefilter
        self.jobs = jobs

    def _prefilter(self, mol_a, mol_b, fps_i=None, fps_j=None) -> bool:
        """Return ``True`` if the pair passes pre-filters (worth computing)."""
        n_a = mol_a.GetNumHeavyAtoms()
        n_b = mol_b.GetNumHeavyAtoms()
        if abs(n_a - n_b) > self.atom_count_diff_threshold:
            return False
        if self.tanimoto_prefilter > 0:
            if fps_i is not None and fps_j is not None:
                if DataStructs.TanimotoSimilarity(fps_i, fps_j) < self.tanimoto_prefilter:
                    return False
            else:
                if _tanimoto(mol_a, mol_b) < self.tanimoto_prefilter:
                    return False
        return True

    def get_score(self, mol_a, mol_b, options: Optional[dict] = None) -> float:
        """Return similarity score for a pair, using cache and pre-filters."""
        options = options or {}
        smi_a = canonical_smiles(mol_a)
        smi_b = canonical_smiles(mol_b)

        cached = self.cache.get(smi_a, smi_b, options)
        if cached is not None:
            return cached

        if not self._prefilter(mol_a, mol_b):
            self.cache.put(smi_a, smi_b, 0.0, options)
            return 0.0

        score = _lomap_score(mol_a, mol_b, options)
        self.cache.put(smi_a, smi_b, score, options)
        return score

    def get_score_matrix(
        self, mols: list, options: Optional[dict] = None
    ) -> np.ndarray:
        """Compute full N×N score matrix with pre-filtering and caching.

        Steps:
        1. Try to serve the entire matrix from cache.
        2. For each uncached pair, apply pre-filters; zero out pairs that fail.
        3. Compute remaining pairs in parallel (or serially when ``jobs==1``).
        4. Store all results in the cache and return the matrix.
        """
        options = options or {}
        n = len(mols)
        smiles_list = [canonical_smiles(m) for m in mols]

        # Fast path: complete cache hit
        cached_matrix = self.cache.get_matrix(smiles_list, options)
        if cached_matrix is not None:
            return cached_matrix

        score_matrix = np.zeros((n, n))
        pairs_to_compute = []

        # Precompute fingerprints once for large N; skip for small N where
        # Tanimoto prefilter rarely fires and computation overhead dominates.
        TANIMOTO_SKIP_THRESHOLD = 10
        if n >= TANIMOTO_SKIP_THRESHOLD and self.tanimoto_prefilter > 0:
            fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024) for m in mols]
        else:
            fps = None

        # Pre-compute MCS seed for all molecules (like Phase1 calc_mcs).
        # Provides a starting substructure so individual LOMAP MCS calls converge faster.
        _mcs_result = rdFMCS.FindMCS(
            mols, timeout=1200,
            atomCompare=rdFMCS.AtomCompare.CompareAny,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            matchValences=False, ringMatchesRingOnly=True,
            completeRingsOnly=True, matchChiralTag=False,
        )
        seed_smarts = _mcs_result.smartsString if _mcs_result.smartsString else ""
        options_with_seed = {**options, "seed": seed_smarts}

        for i, j in itertools.combinations(range(n), 2):
            cached = self.cache.get(smiles_list[i], smiles_list[j], options)
            if cached is not None:
                score_matrix[i][j] = cached
                score_matrix[j][i] = cached
            elif not self._prefilter(
                mols[i], mols[j],
                fps_i=fps[i] if fps is not None else None,
                fps_j=fps[j] if fps is not None else None,
            ):
                self.cache.put(smiles_list[i], smiles_list[j], 0.0, options)
            else:
                pairs_to_compute.append(
                    (i, j, mols[i], mols[j], options_with_seed)
                )

        POOL_THRESHOLD = 4
        if pairs_to_compute:
            n_jobs = self.jobs if self.jobs > 0 else None
            if n_jobs == 1 or len(pairs_to_compute) <= POOL_THRESHOLD:
                results = [_compute_score_worker(p) for p in pairs_to_compute]
            else:
                with Pool(n_jobs) as pool:
                    results = list(pool.imap_unordered(_compute_score_worker, pairs_to_compute))

            for i, j, score in results:
                score_matrix[i][j] = score
                score_matrix[j][i] = score
                self.cache.put(smiles_list[i], smiles_list[j], score, options)

        return score_matrix
