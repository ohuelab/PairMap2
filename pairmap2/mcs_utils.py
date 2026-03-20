"""mcs_utils -- MCS-based scoring utilities."""
import copy
import itertools
from concurrent.futures import as_completed
from multiprocessing import Pool

import numpy as np
from lomap.mcs import MCS
from rdkit.Chem import rdFMCS
from tqdm import tqdm


def formal_charge(mol):
    total_charge_mol = 0.0
    try:
        # Assume mol2
        total_charge_mol = sum([float(a.GetProp('_TriposPartialCharge')) for a in mol.GetAtoms()])
    except Exception:
        # wasn't mol2, so assume SDF with correct formal charge props for mols
        total_charge_mol = sum([a.GetFormalCharge() for a in mol.GetAtoms()])
    return total_charge_mol


def ecr(mol_i, mol_j):
    total_charge_mol_i = formal_charge(mol_i)
    total_charge_mol_j = formal_charge(mol_j)
    if abs(total_charge_mol_j - total_charge_mol_i) < 1e-3:
        scr_ecr = 1.0
    else:
        scr_ecr = 0.0
    return scr_ecr


def score_function(mola, molb, options=None):
    """Calculate the score of two molecules based on various rules."""
    mola = copy.deepcopy(mola)
    molb = copy.deepcopy(molb)
    ecr_score = ecr(mola, molb)
    if options is None:
        options = {'time': 20, 'verbose': 'info', 'max3d': 0, 'threed': False}
    MC = MCS(mola, molb, **options)
    tmp_scr = ecr_score * MC.mncar() * MC.mcsr() * MC.atomic_number_rule() * MC.hybridization_rule()
    tmp_scr *= MC.sulfonamides_rule() * MC.heterocycles_rule()
    tmp_scr *= MC.transmuting_methyl_into_ring_rule()
    tmp_scr *= MC.transmuting_ring_sizes_rule()
    return MC, tmp_scr


def compute_score(pair):
    i, j, mols, options = pair
    try:
        _, score = score_function(mols[i], mols[j], options)
    except Exception:
        _, score = score_function(mols[i], mols[j], options)
    return i, j, score


def compute_score_pair(mol_a, mol_b, options):
    """Worker for ProcessPoolExecutor: receives individual mols, not the full list."""
    try:
        _, score = score_function(mol_a, mol_b, options)
    except Exception:
        _, score = score_function(mol_a, mol_b, options)
    return score


def get_score_matrix(mols, options=None, use_seed=True, jobs=0, executor=None):
    if options is None:
        options = {}
    if use_seed:
        mcs = calc_mcs(mols)
        seedSmarts = mcs.smartsString
    else:
        seedSmarts = ""

    opts = {**options, "seed": seedSmarts}
    N = len(mols)
    score_matrix = np.zeros((N, N))
    pairs = list(itertools.combinations(range(N), 2))

    if executor is not None:
        futures = {
            executor.submit(compute_score_pair, mols[i], mols[j], opts): (i, j)
            for i, j in pairs
        }
        for fut in tqdm(as_completed(futures), total=len(futures)):
            i, j = futures[fut]
            score = fut.result()
            score_matrix[i][j] = score
            score_matrix[j][i] = score
    elif jobs == 1 or jobs == 0:
        serial_pairs = [(i, j, mols, opts) for i, j in pairs]
        for i, j, score in tqdm(map(compute_score, serial_pairs), total=len(serial_pairs)):
            score_matrix[i][j] = score
            score_matrix[j][i] = score
    else:
        if jobs < 0:
            jobs = None
        pool_pairs = [(i, j, mols, opts) for i, j in pairs]
        with Pool(jobs) as pool:
            for i, j, score in tqdm(pool.imap_unordered(compute_score, pool_pairs), total=len(pool_pairs)):
                score_matrix[i][j] = score
                score_matrix[j][i] = score

    return score_matrix


def calc_mcs(mols):
    mcs = rdFMCS.FindMCS(
        mols,
        timeout=30,
        atomCompare=rdFMCS.AtomCompare.CompareAny,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        matchValences=False,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        matchChiralTag=False,
    )
    return mcs
