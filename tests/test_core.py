from simisearch import core
from simisearch.config import SearchConfig, SwissSimilarityConfig, PubChemConfig
from simisearch.clients import swiss_similarity, pubchem

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


def test_search_merges_and_dedupes(monkeypatch):
    # Fake SwissSimilarity hit: caffeine, already-canonical SMILES.
    def fake_swiss(smiles, **kwargs):
        return [
            {
                "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",  # canonical caffeine
                "similarity": 0.75,
                "sources": [{"id": "DB00201", "library": "DrugBank", "similarity": 0.75}],
            }
        ]

    # Fake PubChem hit for the SAME molecule (different SMILES string,
    # same structure) plus a distinct second molecule.
    def fake_pubchem(smiles, **kwargs):
        return [
            {"cid": 2519, "smiles": CAFFEINE, "name": "caffeine", "similarity": 0.90},
            {"cid": 6224, "smiles": "CCO", "name": "ethanol", "similarity": 0.40},
        ]

    monkeypatch.setattr(swiss_similarity, "screen_drugs_and_bioactives", fake_swiss)
    monkeypatch.setattr(pubchem, "similarity_search", fake_pubchem)

    config = SearchConfig(
        swiss=SwissSimilarityConfig(enabled=True),
        pubchem=PubChemConfig(enabled=True),
    )
    results = core.search(ASPIRIN, config=config)

    # Caffeine should appear once (deduped across providers), with the
    # higher of the two similarity scores kept, and both provider
    # sources recorded.
    caffeine_entries = [r for r in results if "n1c" in r["smiles"].lower()]
    assert len(caffeine_entries) == 1
    caffeine_entry = caffeine_entries[0]
    assert caffeine_entry["similarity"] == 0.90
    providers = {s["provider"] for s in caffeine_entry["sources"]}
    assert providers == {"swiss_similarity", "pubchem"}

    # Ethanol should also be present, from PubChem only.
    assert any(r["smiles"] == core.canonical_smiles("CCO") for r in results)

    # Results should be sorted descending by similarity.
    scores = [r["similarity"] for r in results]
    assert scores == sorted(scores, reverse=True)


def test_search_provider_disabled(monkeypatch):
    def fake_pubchem(smiles, **kwargs):
        return [{"cid": 1, "smiles": "CCO", "name": "ethanol", "similarity": 0.5}]

    def boom(*args, **kwargs):
        raise AssertionError("swiss provider should not be called when disabled")

    monkeypatch.setattr(swiss_similarity, "screen_drugs_and_bioactives", boom)
    monkeypatch.setattr(pubchem, "similarity_search", fake_pubchem)

    config = SearchConfig(swiss=SwissSimilarityConfig(enabled=False))
    results = core.search(ASPIRIN, config=config)

    assert len(results) == 1
    assert results[0]["sources"][0]["provider"] == "pubchem"
