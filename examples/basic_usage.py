"""
Basic usage examples for simisearch.

Run with:
    python examples/basic_usage.py
"""

from simisearch import search, SearchConfig, SwissSimilarityConfig, PubChemConfig

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


def default_search():
    """Both providers, default thresholds."""
    results = search(ASPIRIN)
    print(f"Found {len(results)} similar compounds (defaults)")
    for r in results[:5]:
        print(f"  {r['similarity']:.2f}  {r['smiles']}")


def custom_thresholds():
    """Tune similarity thresholds and library selection in plain Python."""
    config = SearchConfig(
        swiss=SwissSimilarityConfig(
            similarity_threshold=0.85,
            drug_libraries=["DrugBank"],
            bioactive_libraries=[],  # skip the ChEMBL bioactive screen
        ),
        pubchem=PubChemConfig(threshold=92, max_records=10),
    )
    results = search(ASPIRIN, config=config)
    print(f"Found {len(results)} similar compounds (stricter thresholds)")
    for r in results[:5]:
        print(f"  {r['similarity']:.2f}  {r['smiles']}")


def pubchem_only():
    """Skip SwissSimilarity entirely (faster, no polling wait)."""
    config = SearchConfig(swiss=SwissSimilarityConfig(enabled=False))
    results = search(ASPIRIN, config=config)
    print(f"Found {len(results)} similar compounds (PubChem only)")
    for r in results[:5]:
        print(f"  {r['similarity']:.2f}  {r['smiles']}")


if __name__ == "__main__":
    default_search()
    print()
    custom_thresholds()
    print()
    pubchem_only()
