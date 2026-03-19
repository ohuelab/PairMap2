"""Canonical SMILES normalization utilities."""
from rdkit import Chem


def canonical_smiles(mol) -> str:
    return Chem.MolToSmiles(mol)
