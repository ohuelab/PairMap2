"""intermediate_search -- BFS-based intermediate molecule search."""
import copy
import logging
import random
from collections import deque

from lomap.mcs import MCS

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem.rdchem import RWMol

from .intermediate_generator import IntermediateGenerator
from .ligand_preparation import execute_ligand_preparation
from .mcs_utils import formal_charge

logger = logging.getLogger(__name__)

__all__ = ["SearchIntermediates"]


class SearchIntermediates:
    def __init__(self, source_ligand, target_ligand, verbose=False, is_atom_modfication_enabled=True, cap_ring_with_carbon=True, cap_ring_with_hydrogen=True, no_backward_search=False, intermediate_name_prefix='intermediate', use_seed=True, score_config=None, ionize=True, max_intermediate=100, search_mode='random', search_random_seed=42):
        # RemoveHs already returns a new molecule object, so deepcopy is not needed here
        self.source_ligand = RWMol(AllChem.RemoveHs(source_ligand))
        self.target_ligand = RWMol(AllChem.RemoveHs(target_ligand))

        self.verbose = verbose
        self.no_backward_search = no_backward_search
        self.use_seed = use_seed
        self.intermediate_name_prefix = intermediate_name_prefix
        self.score_config = score_config if score_config is not None else {}

        self.ionize = ionize
        self.warnings: list[str] = []
        self._charge_mismatch = False
        src_charge = formal_charge(self.source_ligand)
        tgt_charge = formal_charge(self.target_ligand)
        if self.ionize and src_charge != tgt_charge:
            self.warnings.append(
                f"Source (charge={src_charge}) and target (charge={tgt_charge}) have different "
                "formal charges. Charge filtering skipped."
            )
            self._charge_mismatch = True
        self.formal_charge = src_charge

        self.generator = IntermediateGenerator(
            is_atom_modfication_enabled=is_atom_modfication_enabled,
            cap_ring_with_carbon=cap_ring_with_carbon,
            cap_ring_with_hydrogen=cap_ring_with_hydrogen,
            verbose=verbose,
        )

        if self.use_seed:
            _mcs_result = rdFMCS.FindMCS(
                [AllChem.RemoveHs(copy.deepcopy(source_ligand)),
                 AllChem.RemoveHs(copy.deepcopy(target_ligand))],
                timeout=1,
                atomCompare=rdFMCS.AtomCompare.CompareAny,
                bondCompare=rdFMCS.BondCompare.CompareAny,
                matchValences=False, ringMatchesRingOnly=True,
                completeRingsOnly=True, matchChiralTag=False,
            )
            self.seedSmarts = _mcs_result.smartsString if _mcs_result.smartsString else ''
        else:
            self.seedSmarts = ''
        self.max_intermediate = max_intermediate

        self.search_mode = search_mode
        # Single seeded RNG shared across forward and backward searches for determinism
        self._rng = random.Random(search_random_seed)

    def MCS(self, source_ligand, target_ligand):
        try:
            MC = MCS(copy.deepcopy(source_ligand), copy.deepcopy(target_ligand), **self.score_config, seed=self.seedSmarts)
        except Exception:
            MC = MCS(copy.deepcopy(source_ligand), copy.deepcopy(target_ligand), **self.score_config, seed='')
        return MC

    def simplex_search(self, direction='forward'):
        if direction == 'forward':
            source_ligand = copy.deepcopy(self.source_ligand)
            target_ligand = copy.deepcopy(self.target_ligand)
        else:
            source_ligand = copy.deepcopy(self.target_ligand)
            target_ligand = copy.deepcopy(self.source_ligand)

        source_smiles = Chem.MolToSmiles(source_ligand)
        target_smiles = Chem.MolToSmiles(target_ligand)

        intermediate_info_list = [self.get_intermediate_info(source_ligand)]
        intermediate_info_list += [self.get_intermediate_info(target_ligand)]
        source_index = 0
        target_index = 1
        smiles_list = [source_smiles, target_smiles]
        smiles_to_idx = {source_smiles: 0, target_smiles: 1}
        traces = [[source_index, target_index]]
        depth_map = {source_smiles: 0, target_smiles: 0}

        if self.search_mode == 'bfs':
            q = deque([source_ligand])
        else:
            q = [source_ligand]

        while len(q) > 0 and (self.max_intermediate <= 0 or len(intermediate_info_list) < self.max_intermediate):
            if self.search_mode == 'bfs':
                ligand = q.popleft()
            else:
                # Random selection with O(1) swap-with-last removal
                idx = self._rng.randint(0, len(q) - 1)
                q[idx], q[-1] = q[-1], q[idx]
                ligand = q.pop()

            ligand = RWMol(ligand)
            smiles = Chem.MolToSmiles(ligand)
            ligand_index = smiles_to_idx[smiles]
            ligand_depth = depth_map[smiles]
            MC = self.MCS(ligand, target_ligand)
            mcs_map = {a1: a2 for a1, a2 in MC.heavy_atom_mcs_map()}
            intermediates = self.generator.generate_intermediates(ligand, target_ligand, mcs_map)
            for intermediate in intermediates:
                info = self.get_intermediate_info(intermediate)
                if info['smiles'] not in smiles_to_idx:
                    if self.verbose:
                        logger.debug('intermediate: %s', info['smiles'])
                    new_idx = len(smiles_list)
                    smiles_to_idx[info['smiles']] = new_idx
                    smiles_list.append(info['smiles'])
                    intermediate_info_list.append(info)
                    depth_map[info['smiles']] = ligand_depth + 1
                    if self.search_mode == 'bfs':
                        q.append(intermediate)
                    else:
                        q.append(intermediate)
                intermediate_index = smiles_to_idx[info['smiles']]
                traces.append([ligand_index, intermediate_index])
                traces.append([intermediate_index, target_index])
                if self.max_intermediate > 0 and len(intermediate_info_list) >= self.max_intermediate:
                    break
        return intermediate_info_list, traces, depth_map

    def get_intermediate_info(self, ligand):
        if Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(ligand))) != Chem.MolToSmiles(ligand):
            raise ValueError('smiles must be canonical, please report if this error occurs: {} != {}'.format(
                Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(ligand))), Chem.MolToSmiles(ligand)))
        smiles = Chem.MolToSmiles(ligand)
        return {
            'ligand': ligand,
            'smiles': smiles,
        }

    def merge_intermediate_info_list(self, prefix='intermediate'):
        forward_intermediate_info_list = self.forward_intermediate_info_list
        forward_traces = self.forward_traces
        backward_intermediate_info_list = self.backward_intermediate_info_list
        backward_traces = self.backward_traces

        forward_intermediate_smiles = [info['smiles'] for info in forward_intermediate_info_list]
        backward_intermediate_info_list_uniq = [info for info in backward_intermediate_info_list if info['smiles'] not in forward_intermediate_smiles]
        intermediate_info_list = forward_intermediate_info_list + backward_intermediate_info_list_uniq
        intermediate_smiles = [info['smiles'] for info in intermediate_info_list]

        for i, info in enumerate(intermediate_info_list):
            info['name'] = f'{prefix}-{i:04d}'
            info['ligand'].SetProp('_Name', info['name'])
            info['ligand'].SetProp('NAME', info['name'])
            info['ligand'].SetProp('smiles', info['smiles'])

        backward_traces_reindex = []
        reindexmap = {i: intermediate_smiles.index(info['smiles']) for i, info in enumerate(backward_intermediate_info_list)}
        for source_idx, target_idx in backward_traces:
            backward_traces_reindex.append([reindexmap[target_idx], reindexmap[source_idx]])
        intermediate_traces = forward_traces + backward_traces_reindex
        self.intermediate_info_list = intermediate_info_list
        self.intermediate_smiles = intermediate_smiles
        self.intermediate_traces = intermediate_traces

        # Merge depth maps: prefer forward depth; backward depth fills in unknowns
        depth_map = dict(self.forward_depth_map)
        for smi, d in self.backward_depth_map.items():
            if smi not in depth_map:
                depth_map[smi] = d
        self.depth_map = depth_map

        return self.intermediate_info_list

    def show_result(self):
        logger.info('Number of total intermediates: %d', len(self.intermediate_info_list))
        if self.ionize and not self._charge_mismatch:
            logger.info('Number of intermediates after charge filtering: %d', max(0, len(self.intermediates) - 2))
        for w in self.warnings:
            logger.warning(w)

    def search(self):
        self.forward_intermediate_info_list, self.forward_traces, self.forward_depth_map = self.simplex_search('forward')
        if self.no_backward_search:
            self.backward_intermediate_info_list, self.backward_traces, self.backward_depth_map = [], [], {}
        else:
            self.backward_intermediate_info_list, self.backward_traces, self.backward_depth_map = self.simplex_search('backward')

        intermediate_info_list = self.merge_intermediate_info_list(self.intermediate_name_prefix)
        self.intermediates_all = [info['ligand'] for info in intermediate_info_list]
        if self.ionize:
            self.intermediates_ionized = execute_ligand_preparation(self.intermediates_all)
            # Index 0 = source, 1 = target: always kept regardless of charge.
            # Only filter intermediates (index 2+) by charge.
            if self._charge_mismatch:
                self.intermediates = list(self.intermediates_ionized)
            else:
                ionized_src = self.intermediates_ionized[0] if len(self.intermediates_ionized) > 0 else self.intermediates_all[0]
                ionized_tgt = self.intermediates_ionized[1] if len(self.intermediates_ionized) > 1 else self.intermediates_all[1]
                ref_charge = formal_charge(ionized_src)
                filtered = [mol for mol in self.intermediates_ionized[2:] if formal_charge(mol) == ref_charge]
                self.intermediates = [ionized_src, ionized_tgt] + filtered
        else:
            self.intermediates = self.intermediates_all
        if self.verbose:
            self.show_result()
        return self.intermediates
