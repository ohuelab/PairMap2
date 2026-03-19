"""intermediate_generator -- generates intermediate molecules between two ligands."""
import copy
import logging

from rdkit import Chem
from rdkit.Chem import MolStandardize

HYDROGEN_ATOM = Chem.Atom(1)
CARBON_ATOM = Chem.Atom(6)

logger = logging.getLogger(__name__)


class IntermediateGenerator:
    def __init__(self, is_atom_modfication_enabled=True, cap_ring_with_carbon=True, cap_ring_with_hydrogen=True, verbose=False):
        '''
        :param is_atom_modfication_enabled: Whether to enable atom modification.
        :param cap_ring_with_carbon: Whether to cap rings with carbon atoms.
        :param cap_ring_with_hydrogen: Whether to cap rings with hydrogen atoms.
        :param verbose: Whether to print verbose output.
        '''
        if not cap_ring_with_carbon and not cap_ring_with_hydrogen:
            raise ValueError('At least one of the options cap_ring_with_carbon and cap_ring_with_hydrogen must be True.')

        self.is_atom_modfication_enabled = is_atom_modfication_enabled
        self.cap_ring_with_carbon = cap_ring_with_carbon
        self.cap_ring_with_hydrogen = cap_ring_with_hydrogen
        self.verbose = verbose
        self.lfc = MolStandardize.rdMolStandardize.LargestFragmentChooser()

    @staticmethod
    def remove_props(mol):
        '''Remove all properties from a molecule'''
        for prop in mol.GetPropsAsDict():
            mol.ClearProp(prop)

    @staticmethod
    def remove_atom_map(mol):
        '''Remove all atom map from a molecule'''
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return mol

    def postprocess_ligand(self, ligand):
        '''Postprocess a molecule'''
        try:
            self.remove_props(ligand)
            Chem.SanitizeMol(ligand)
            ligand = self.remove_atom_map(ligand)
            ligand = Chem.RemoveHs(ligand)
            Chem.AssignStereochemistryFrom3D(ligand)
            return ligand
        except Exception:
            return None

    @staticmethod
    def get_atom_idx_by_map_num(mol, mapnum):
        '''Get atom index from atom map number'''
        indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() == mapnum]
        return indices[0] if len(indices) == 1 else None

    def get_terminal_rings(self, rings):
        '''Identify rings with exactly one neighboring atom outside the ring.'''
        ligand = self.source_ligand
        terminal_rings = []
        for ring in rings:
            external_neighbor_indices = []
            for idx in ring:
                atom = ligand.GetAtomWithIdx(idx)
                for neighbor_atom in atom.GetNeighbors():
                    neighbor_idx = neighbor_atom.GetIdx()
                    if (neighbor_idx not in ring) and neighbor_idx not in external_neighbor_indices:
                        external_neighbor_indices.append(neighbor_idx)
            if len(external_neighbor_indices) == 1:
                terminal_rings.append(ring)
        return terminal_rings

    def get_fused_rings(self, rings):
        """Get the fused rings from a ligand and a list of rings."""
        ligand = self.source_ligand
        fused_rings = []
        for i in range(len(rings)):
            for j in range(0, len(rings)):
                if i == j:
                    continue
                ring1 = set(rings[i])
                ring2 = set(rings[j])
                if len(ring1.intersection(ring2)) > 0:
                    ring_diff = ring1.difference(ring2)
                    fused_ring = set()
                    for idx in ring_diff:
                        atom = ligand.GetAtomWithIdx(idx)
                        for neighbor_atom in atom.GetNeighbors():
                            neighbor_idx = neighbor_atom.GetIdx()
                            if (neighbor_idx not in ring1):
                                fused_ring.add(idx)
                    fused_rings += [list(ring_diff)]
        return fused_rings

    def handle_deletable_entities(self, atom, terminal_rings, fused_rings):
        '''Identify deletable atoms and rings.'''
        ligand = self.source_ligand
        deletable_atoms = []
        deletable_rings = []
        deletable_fused_rings = []
        atom_idx = atom.GetIdx()

        if not atom.IsInRing():
            if atom.GetDegree() == 1:
                deletable_atoms.append(atom)
        else:
            for ring in terminal_rings:
                if atom_idx in ring:
                    atom_ring = [ligand.GetAtomWithIdx(idx) for idx in ring]
                    deletable_rings.append(atom_ring)
                    terminal_rings.remove(ring)
            for ring in fused_rings:
                if atom_idx in ring:
                    atom_ring = [ligand.GetAtomWithIdx(idx) for idx in ring]
                    deletable_fused_rings.append(atom_ring)
                    fused_rings.remove(ring)

        return deletable_atoms, deletable_rings, deletable_fused_rings

    def extract_atoms_for_modification_and_deletion(self):
        '''Extract atoms from the source ligand that can be modified or deleted.'''
        source_ligand, target_ligand, mcs_map = self.source_ligand, self.target_ligand, self.mcs_map

        source_rings = list(source_ligand.GetRingInfo().AtomRings())
        terminal_rings = self.get_terminal_rings(source_rings)
        fused_rings = self.get_fused_rings(source_rings)

        atoms_for_modification = []
        atoms_for_deletion = []
        rings_for_deletion = []
        fused_rings_for_deletion = []

        for atom in source_ligand.GetAtoms():
            source_atom_idx = atom.GetIdx()
            if source_atom_idx in mcs_map:
                target_atom_idx = mcs_map[source_atom_idx]
                target_atom = target_ligand.GetAtomWithIdx(target_atom_idx)
                if atom.GetSymbol() != target_atom.GetSymbol():
                    atoms_for_modification.append(atom)
            else:
                deletable_atoms, deletable_rings, deletable_fused_rings = self.handle_deletable_entities(
                    atom, terminal_rings, fused_rings)
                atoms_for_deletion.extend(deletable_atoms)
                rings_for_deletion.extend(deletable_rings)
                fused_rings_for_deletion.extend(deletable_fused_rings)
        return (
            atoms_for_modification,
            atoms_for_deletion,
            rings_for_deletion,
            fused_rings_for_deletion,
        )

    def generate_atom_modification_intermediate(self, atom_map_num):
        '''Generate an intermediate by modifying a specific atom.'''
        source_ligand, target_ligand, mcs_map = self.source_ligand, self.target_ligand, self.mcs_map

        intermediate_ligand = copy.deepcopy(source_ligand)
        source_atom_idx = self.get_atom_idx_by_map_num(intermediate_ligand, atom_map_num)
        target_atom_idx = mcs_map[source_atom_idx]
        target_atomic_num = target_ligand.GetAtomWithIdx(target_atom_idx).GetAtomicNum()

        intermediate_ligand.ReplaceAtom(source_atom_idx, Chem.Atom(target_atomic_num))
        return intermediate_ligand

    def generate_atom_deletion_intermediate(self, atom_map_num):
        '''Generate an intermediate by deleting a specific atom.'''
        source_ligand = self.source_ligand

        intermediate_ligand = copy.deepcopy(source_ligand)
        atom_idx = self.get_atom_idx_by_map_num(intermediate_ligand, atom_map_num)
        atom_to_delete = intermediate_ligand.GetAtomWithIdx(atom_idx)
        assert len(atom_to_delete.GetBonds()) == 1
        bond_to_remove = atom_to_delete.GetBonds()[0]

        if bond_to_remove.GetBondType() == Chem.rdchem.BondType.SINGLE:
            intermediate_ligand.ReplaceAtom(atom_idx, HYDROGEN_ATOM)
        else:
            intermediate_ligand.RemoveAtom(atom_idx)

        return intermediate_ligand

    def generate_ring_deletion_intermediate(self, ring_atom_map_nums):
        '''Generate an intermediate by deleting a specific ring structure.'''
        source_ligand = self.source_ligand
        intermediate_ligand = copy.deepcopy(source_ligand)
        atoms_to_cap = []

        for map_num in ring_atom_map_nums:
            atom_idx = self.get_atom_idx_by_map_num(intermediate_ligand, map_num)
            neighbor_indices_outside_ring = self.get_neighbor_indices_outside_ring(atom_idx, ring_atom_map_nums, intermediate_ligand)

            if neighbor_indices_outside_ring:
                atoms_to_cap.append(map_num)
            else:
                intermediate_ligand.RemoveAtom(atom_idx)

        assert len(set(atoms_to_cap)) == 1
        map_num = atoms_to_cap[0]
        atom_idx = self.get_atom_idx_by_map_num(intermediate_ligand, map_num)
        intermediate_ligands = []

        if self.cap_ring_with_carbon:
            new_intermediate_ligand = copy.deepcopy(intermediate_ligand)
            new_intermediate_ligand.ReplaceAtom(atom_idx, CARBON_ATOM)
            intermediate_ligands.append(new_intermediate_ligand)

        if self.cap_ring_with_hydrogen:
            new_intermediate_ligand = copy.deepcopy(intermediate_ligand)
            new_intermediate_ligand.ReplaceAtom(atom_idx, HYDROGEN_ATOM)
            intermediate_ligands.append(new_intermediate_ligand)

        return intermediate_ligands

    def generate_fused_ring_deletion_intermediate(self, ring_atom_map_nums):
        '''Generate intermediates by deleting a fused ring structure.'''
        source_ligand = self.source_ligand
        intermediate_ligand = copy.deepcopy(source_ligand)
        neighbor_map_nums = []

        for map_num in ring_atom_map_nums:
            atom_idx = self.get_atom_idx_by_map_num(intermediate_ligand, map_num)
            neighbor_indices_outside_ring = self.get_neighbor_indices_outside_ring(atom_idx, ring_atom_map_nums, intermediate_ligand)

            if neighbor_indices_outside_ring:
                neighbor_map_nums.extend([intermediate_ligand.GetAtomWithIdx(idx).GetAtomMapNum() for idx in neighbor_indices_outside_ring])
            intermediate_ligand.RemoveAtom(atom_idx)

        neighbor_map_nums = list(set(neighbor_map_nums))
        intermediate_ligands = [intermediate_ligand]

        for i, map_num_i in enumerate(neighbor_map_nums):
            for map_num_j in neighbor_map_nums[i + 1:]:
                idx_i = self.get_atom_idx_by_map_num(intermediate_ligand, map_num_i)
                idx_j = self.get_atom_idx_by_map_num(intermediate_ligand, map_num_j)
                if intermediate_ligand.GetBondBetweenAtoms(idx_i, idx_j):
                    new_intermdiate_ligands = []
                    for intermdiate_ligand in intermediate_ligands:
                        new_intermdiate_ligand = copy.deepcopy(intermdiate_ligand)
                        bond = new_intermdiate_ligand.GetBondBetweenAtoms(idx_i, idx_j)
                        bond.SetBondType(Chem.rdchem.BondType.SINGLE)
                        bond.SetIsAromatic(False)
                        new_intermdiate_ligand.GetAtomWithIdx(idx_i).SetIsAromatic(False)
                        new_intermdiate_ligand.GetAtomWithIdx(idx_j).SetIsAromatic(False)
                        new_intermdiate_ligands.append(new_intermdiate_ligand)
                    intermediate_ligands.extend(new_intermdiate_ligands)
        return intermediate_ligands

    @staticmethod
    def get_neighbor_indices_outside_ring(atom_idx, ring_atom_map_nums, ligand):
        '''Get indices of neighboring atoms not part of the specified ring.'''
        atom_neighbors = ligand.GetAtomWithIdx(atom_idx).GetNeighbors()
        neighbor_indices_outside_ring = []
        for neighbor_atom in atom_neighbors:
            neighbor_idx = neighbor_atom.GetIdx()
            neighbor_map_num = ligand.GetAtomWithIdx(neighbor_idx).GetAtomMapNum()
            if neighbor_map_num not in ring_atom_map_nums:
                neighbor_indices_outside_ring.append(neighbor_idx)
        return neighbor_indices_outside_ring

    def generate_intermediates(self, source_ligand, target_ligand, mcs_map):
        '''Generate intermediates for transforming source into target.'''
        self.source_ligand = source_ligand
        self.target_ligand = target_ligand
        for i, atom in enumerate(source_ligand.GetAtoms(), start=1):
            atom.SetAtomMapNum(i)
        self.mcs_map = mcs_map

        atoms_for_modification, atoms_for_deletion, rings_for_deletion, fused_rings_for_deletion = self.extract_atoms_for_modification_and_deletion()
        if self.verbose:
            logger.debug('atoms_for_modification: %d', len(atoms_for_modification))
            logger.debug('atoms_for_deletion: %d', len(atoms_for_deletion))
            logger.debug('rings_for_deletion: %d', len(rings_for_deletion))
            logger.debug('fused_rings_for_deletion: %d', len(fused_rings_for_deletion))

        atom_map_nums_for_modification = [atom.GetAtomMapNum() for atom in atoms_for_modification]
        atom_map_nums_for_deletion = [atom.GetAtomMapNum() for atom in atoms_for_deletion]
        ring_atom_map_nums_list_for_deletion = [[atom.GetAtomMapNum() for atom in ring] for ring in rings_for_deletion]
        fused_ring_atom_map_nums_list_for_deletion = [[atom.GetAtomMapNum() for atom in ring] for ring in fused_rings_for_deletion]

        intermediates_list = []
        if self.is_atom_modfication_enabled:
            intermediates_list.extend([[self.generate_atom_modification_intermediate(atom_map_num)] for atom_map_num in atom_map_nums_for_modification])

        intermediates_list.extend([[self.generate_atom_deletion_intermediate(atom_map_num)] for atom_map_num in atom_map_nums_for_deletion])
        intermediates_list.extend([self.generate_ring_deletion_intermediate(ring_atom_map_nums) for ring_atom_map_nums in ring_atom_map_nums_list_for_deletion])
        intermediates_list.extend([self.generate_fused_ring_deletion_intermediate(fused_ring_atom_map_nums) for fused_ring_atom_map_nums in fused_ring_atom_map_nums_list_for_deletion])

        intermediates = [self.postprocess_ligand(intermediate) for intermediates in intermediates_list for intermediate in intermediates]
        intermediates = [intermediate for intermediate in intermediates if intermediate is not None]

        return intermediates
