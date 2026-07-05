"""
core.py

The main entry point: search(smiles, config) queries whichever
providers are enabled in the config, merges their hits, deduplicates
by canonical SMILES, and returns a single ranked list.
"""

from .config import SearchConfig
from .utils import canonical_smiles
from .clients import swiss_similarity, pubchem


def search(smiles: str, config: SearchConfig = None):
    """
    Find structurally similar compounds for a query SMILES.

    Parameters
    ----------
    smiles : str
        Query SMILES string.
    config : SearchConfig, optional
        Controls which providers run and their parameters (similarity
        thresholds, libraries, result limits, ...). Defaults to
        SearchConfig() if not given.

    Returns
    -------
    list of dicts, sorted by similarity descending (unless
    config.sort_descending is False), each with:
        "smiles"     : str   - canonical SMILES (dedup key)
        "similarity" : float - highest similarity score seen across
                                providers/libraries for this molecule
        "sources"    : list  - one entry per provider/library the
                                molecule was found under, e.g.
                                {"provider": "swiss_similarity",
                                 "similarity": 0.82,
                                 "libraries": [{"id": "DB00945",
                                                "library": "DrugBank",
                                                "similarity": 0.82}]}
                                {"provider": "pubchem", "cid": 2244,
                                 "name": "aspirin", "similarity": 0.9}

    Notes
    -----
    - If neither provider is enabled, returns an empty list.
    - SwissSimilarity screens can take 30s-5min per library (it's a
      polling job against a remote queue); PubChem is faster.
    - Canonicalization/dedup requires RDKit-parseable SMILES; hits with
      unparseable SMILES are dropped.
    """

    if config is None:
        config = SearchConfig()

    canon_query = canonical_smiles(smiles)
    if canon_query is None:
        raise ValueError(f"Could not parse query SMILES: {smiles!r}")

    merged = {}

    def _add(canon, similarity, source_entry):
        if config.dedupe and canon in merged:
            entry = merged[canon]
            entry["similarity"] = max(entry["similarity"], similarity)
            entry["sources"].append(source_entry)
        else:
            key = canon if config.dedupe else f"{canon}#{len(merged)}"
            merged[key] = {
                "smiles": canon,
                "similarity": similarity,
                "sources": [source_entry],
            }

    if config.swiss.enabled:
        swiss_hits = swiss_similarity.screen_drugs_and_bioactives(
            smiles,
            method=config.swiss.method,
            similarity_threshold=config.swiss.similarity_threshold,
            drug_libraries=config.swiss.drug_libraries,
            bioactive_libraries=config.swiss.bioactive_libraries,
            poll_interval=config.swiss.poll_interval_seconds,
            poll_max_tries=config.swiss.poll_max_tries,
        )
        for hit in swiss_hits:
            _add(
                hit["smiles"],
                hit["similarity"],
                {
                    "provider": "swiss_similarity",
                    "similarity": hit["similarity"],
                    "libraries": hit["sources"],
                },
            )

    if config.pubchem.enabled:
        pubchem_hits = pubchem.similarity_search(
            smiles,
            threshold=config.pubchem.threshold,
            max_records=config.pubchem.max_records,
            poll_interval=config.pubchem.poll_interval_seconds,
            poll_max_tries=config.pubchem.poll_max_tries,
        )
        for hit in pubchem_hits:
            canon = canonical_smiles(hit["smiles"])
            if canon is None:
                continue
            _add(
                canon,
                hit["similarity"],
                {
                    "provider": "pubchem",
                    "cid": hit["cid"],
                    "name": hit["name"],
                    "similarity": hit["similarity"],
                },
            )

    results = list(merged.values())

    # Drop the query molecule itself if a provider echoed it back.
    results = [r for r in results if r["smiles"] != canon_query]

    results.sort(key=lambda x: x["similarity"], reverse=config.sort_descending)

    return results
