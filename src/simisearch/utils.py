"""
utils.py

Small shared helpers used across simisearch providers.
"""

from rdkit import Chem


def canonical_smiles(smiles: str):
    """
    Canonicalize a SMILES string via RDKit so structurally identical
    molecules returned by different providers (which use different ID
    systems - PubChem CIDs, DrugBank IDs, ChEMBL IDs) can be matched
    for deduplication.

    Returns the canonical SMILES string, or None if the SMILES cannot
    be parsed by RDKit.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None


def validate_smiles(smiles: str) -> bool:
    """Return True if RDKit can parse the given SMILES string."""
    return canonical_smiles(smiles) is not None
