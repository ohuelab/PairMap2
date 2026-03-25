"""ligand_preparation -- OpenBabel-based ligand protonation via Python API."""
from openbabel import pybel

from rdkit import Chem

from .mcs_utils import formal_charge


def execute_ligand_preparation(mols, pH=7.4, extract_same_formal_charge=False, charge_indices=None):
    """Execute ligand preparation using OpenBabel Python API.

    :param mols: A list of RDKit Mol objects.
    :param pH: The pH value for the protonation.
    :param extract_same_formal_charge: Whether to filter to keep only molecules
        with the same formal charge as the reference molecule.
    :param charge_indices: Indices used to determine the reference charge.
        Defaults to [0, 1].
    :return: A list of prepared RDKit Mol objects.
    """
    if charge_indices is None:
        charge_indices = [0, 1]

    prepared = []
    for mol in mols:
        try:
            mol_block = Chem.MolToMolBlock(mol)
            obmol = pybel.readstring("sdf", mol_block)
            obmol.OBMol.AddHydrogens(False, True, pH)
            prepared_block = obmol.write("sdf")
            prepared_mol = Chem.MolFromMolBlock(prepared_block, removeHs=True)
            if prepared_mol is not None:
                # Copy RDKit properties from original mol to the protonated mol
                for prop_name in mol.GetPropsAsDict():
                    prepared_mol.SetProp(prop_name, str(mol.GetPropsAsDict()[prop_name]))
                prepared.append(prepared_mol)
            else:
                prepared.append(mol)
        except Exception:
            prepared.append(mol)

    if extract_same_formal_charge:
        formal_charges = [formal_charge(mol) for mol in prepared]
        if isinstance(charge_indices, int):
            charge_indices = [charge_indices]
        if not isinstance(charge_indices, list):
            raise ValueError('charge_indices must be a list of integers or a single integer.')
        if len(charge_indices) == 0:
            raise ValueError('charge_indices must not be empty.')
        if not all(isinstance(i, int) for i in charge_indices):
            raise ValueError('charge_indices must be a list of integers or a single integer.')

        base_charge = formal_charges[charge_indices[0]]
        mismatched = [i for i in charge_indices if formal_charges[i] != base_charge]
        if mismatched:
            charges_str = ', '.join(f'mol[{i}]={formal_charges[i]}' for i in charge_indices)
            raise ValueError(
                f'Formal charges of the molecules are not the same ({charges_str}). '
                f'Expected all to be {base_charge}.'
            )
        prepared = [mol for mol, charge in zip(prepared, formal_charges) if charge == base_charge]

    return prepared
