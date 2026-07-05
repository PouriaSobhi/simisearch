from simisearch.utils import canonical_smiles, validate_smiles


def test_canonical_smiles_valid():
    # Two SMILES for the same molecule (ethanol) should canonicalize
    # to the same string.
    a = canonical_smiles("CCO")
    b = canonical_smiles("OCC")
    assert a is not None
    assert a == b


def test_canonical_smiles_invalid():
    assert canonical_smiles("not_a_smiles(((") is None


def test_validate_smiles():
    assert validate_smiles("CC(=O)Oc1ccccc1C(=O)O") is True  # aspirin
    assert validate_smiles("???") is False
