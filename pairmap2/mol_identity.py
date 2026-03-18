"""Canonical SMILES normalization utilities."""
from rdkit import Chem


def canonical_smiles(mol) -> str:
    return Chem.MolToSmiles(mol)


def normalize_mol(mol):
    """Return a new mol with canonical atom ordering.

    Removes explicit Hs, assigns stereochemistry from the canonical SMILES
    round-trip so that the molecule is in a reproducible state.
    """
    mol = Chem.RemoveHs(mol)
    smi = Chem.MolToSmiles(mol)
    return Chem.MolFromSmiles(smi)
