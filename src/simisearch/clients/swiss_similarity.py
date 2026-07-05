"""
clients/swiss_similarity.py

Client for the SwissSimilarity command-line API
(www.swisssimilarity.ch:1234).

Workflow per library:
    1. POST SMILES file -> startscreen -> get Session Number
    2. Poll checksession until "Calculation is finished"
    3. GET retrievesession -> plain-text "ID score SMILES" lines

Searches real compound libraries (DrugBank, ChEMBL approved drugs,
full ChEMBL bioactive set, and others) for structurally similar
molecules.

Library families (per SwissSimilarity's own grouping - see
swisssimilarity.ch/FAQ.php):
    "Drugs"     group -> DrugBank, CHEMBL_drug, CHEMBL_clinic
    "Bioactive" group -> CHEMBL (full), LigandExpo, CHEBI, Glass, HMDB

Merging and deduplication
--------------------------
Different libraries use different ID systems (DrugBank: "DB00945",
ChEMBL: "CHEMBL25") for what may be the same underlying molecule.
True deduplication is done on canonical SMILES (not on ID), since the
server returns SMILES directly with each hit. When the same molecule
is found in multiple libraries, the highest similarity score is kept
and all source libraries/IDs are recorded.

Reference: https://www.swisssimilarity.ch/command-line.php
"""

import time
import requests

from ..utils import canonical_smiles

SWISS_SIMILARITY_HOST = "https://www.swisssimilarity.ch:1234"

DEFAULT_METHOD = "Combined"
DEFAULT_SIMILARITY_THRESHOLD = 0.7
DEFAULT_POLL_INTERVAL = 5.0
DEFAULT_POLL_MAX_TRIES = 60  # ~5 minutes max wait per library

DRUG_LIBRARIES = ["DrugBank", "CHEMBL_drug"]
BIOACTIVE_LIBRARIES = ["CHEMBL"]

VALID_LIBRARIES = {
    "Asinex", "AsisChem", "CHEBI", "CHEMBL", "CHEMBL_drug", "CHEMBL_clinic",
    "CHEMBL_act", "CHEMBL_GPCR", "CHEMBL_kinase", "CHEMBL_protease",
    "ChemBridge", "ChemDiv", "DrugBank", "Enamine", "Enamine_Tang", "Glass",
    "HMDB", "InnovaPharm", "InnovaPharm_Tang", "LifeChemicals",
    "LifeChemicals_Tang", "LigandExpo", "Maybridge", "Otava", "Otava_Tang",
    "SPECS", "TimTec", "Vitas", "ZINC_drug", "ZINC_frag", "ZINC_lead",
}

VALID_METHODS = {
    "ECFP", "FP2", "MHFP6", "pharmacophore", "ERG", "scaffold",
    "GenScaffold", "ES5D", "E3FP", "Combined",
}

# In-memory cache: (smiles, library, method) -> list of result dicts
_screen_cache = {}


# ==========================================================
# Core single-library workflow
# ==========================================================

def _start_screen(smiles, library, method):

    if library not in VALID_LIBRARIES:
        raise ValueError(
            f"Unknown library '{library}'. Valid options: {sorted(VALID_LIBRARIES)}"
        )

    if method not in VALID_METHODS:
        raise ValueError(
            f"Unknown method '{method}'. Valid options: {sorted(VALID_METHODS)}"
        )

    url = f"{SWISS_SIMILARITY_HOST}/startscreen"

    try:
        response = requests.post(
            url,
            params={"library": library, "method": method},
            files={"mySMILES": ("query.smi", smiles)},
            timeout=30,
        )
        response.raise_for_status()
        session_number = response.text.strip()

        if not session_number or not session_number.isdigit():
            return None

        return session_number

    except Exception:
        return None


def _check_session(session_number):

    url = f"{SWISS_SIMILARITY_HOST}/checksession"

    try:
        response = requests.get(
            url,
            params={"sessionNumber": session_number},
            timeout=30,
        )
        response.raise_for_status()
        return response.text.strip()

    except Exception:
        return None


def _retrieve_results(session_number):

    url = f"{SWISS_SIMILARITY_HOST}/retrievesession"

    try:
        response = requests.get(
            url,
            params={"sessionNumber": session_number},
            timeout=60,
        )
        response.raise_for_status()
        return response.text

    except Exception:
        return None


def _parse_results(raw_text, library):
    """
    Parse "ID score [SMILES]" lines into structured dicts.

    Works for both ChEMBL-style IDs (CHEMBL25) and DrugBank-style IDs
    (DB00945), since the format is identical (ID score SMILES).
    """

    results = []

    for line in raw_text.strip().split("\n"):

        line = line.strip()

        if not line:
            continue

        parts = line.split(maxsplit=2)

        if len(parts) < 2:
            continue

        compound_id = parts[0]
        score_str = parts[1]
        smiles = parts[2] if len(parts) > 2 else ""

        try:
            score = float(score_str)
        except ValueError:
            continue

        results.append({
            "id": compound_id,
            "library": library,
            "similarity": score,
            "smiles": smiles,
        })

    return results


def similarity_screen(
    smiles,
    library,
    method=DEFAULT_METHOD,
    poll_interval=DEFAULT_POLL_INTERVAL,
    poll_max_tries=DEFAULT_POLL_MAX_TRIES,
):
    """
    Screen a query SMILES against a single SwissSimilarity library.

    Returns a list of dicts (unfiltered, unsorted beyond raw order):
        "id"         : str   - compound ID (ChEMBL or DrugBank format)
        "library"    : str   - source library name
        "similarity" : float - similarity score (0-1)
        "smiles"     : str   - SMILES of the hit compound

    Returns empty list on timeout or failure. Cached per
    (smiles, library, method).
    """

    cache_key = (smiles, library, method)
    if cache_key in _screen_cache:
        return _screen_cache[cache_key]

    session_number = _start_screen(smiles, library, method)

    if not session_number:
        _screen_cache[cache_key] = []
        return []

    finished = False

    for _ in range(poll_max_tries):
        time.sleep(poll_interval)

        status = _check_session(session_number)

        if status is None:
            continue

        if "finished" in status.lower():
            finished = True
            break

        if "cancel" in status.lower() or "error" in status.lower():
            break

    if not finished:
        _screen_cache[cache_key] = []
        return []

    raw_text = _retrieve_results(session_number)

    if not raw_text:
        _screen_cache[cache_key] = []
        return []

    results = _parse_results(raw_text, library)

    _screen_cache[cache_key] = results
    return results


# ==========================================================
# Public API: multi-library merge
# ==========================================================

def screen_drugs_and_bioactives(
    smiles,
    method=DEFAULT_METHOD,
    similarity_threshold=DEFAULT_SIMILARITY_THRESHOLD,
    drug_libraries=None,
    bioactive_libraries=None,
    poll_interval=DEFAULT_POLL_INTERVAL,
    poll_max_tries=DEFAULT_POLL_MAX_TRIES,
):
    """
    Screen a query SMILES against both the "Drugs" and "Bioactive"
    library families, merge results, deduplicate by canonical SMILES,
    and filter by similarity threshold.

    Returns
    -------
    list of dicts, sorted by similarity descending, each with:
        "smiles"  : str   - canonical SMILES (dedup key)
        "similarity" : float - highest similarity score across all
                                libraries it was found in
        "sources" : list  - [{"id": str, "library": str,
                              "similarity": float}, ...] for every
                            library/ID this molecule was found under

    Notes
    -----
    - Each library screen takes 30s - 5min depending on size/queue,
      run sequentially.
    - Molecules with unparseable SMILES are dropped.
    - If a library screen fails/times out, it's skipped silently and
      the remaining libraries still contribute results.
    """

    if drug_libraries is None:
        drug_libraries = DRUG_LIBRARIES

    if bioactive_libraries is None:
        bioactive_libraries = BIOACTIVE_LIBRARIES

    all_libraries = list(drug_libraries) + list(bioactive_libraries)

    # canonical_smiles -> merged result dict
    merged = {}

    for library in all_libraries:

        hits = similarity_screen(
            smiles,
            library=library,
            method=method,
            poll_interval=poll_interval,
            poll_max_tries=poll_max_tries,
        )

        for hit in hits:

            if hit["similarity"] < similarity_threshold:
                continue

            canon = canonical_smiles(hit["smiles"])

            if canon is None:
                continue

            source_entry = {
                "id": hit["id"],
                "library": hit["library"],
                "similarity": hit["similarity"],
            }

            if canon not in merged:
                merged[canon] = {
                    "smiles": canon,
                    "similarity": hit["similarity"],
                    "sources": [source_entry],
                }
            else:
                merged[canon]["sources"].append(source_entry)
                if hit["similarity"] > merged[canon]["similarity"]:
                    merged[canon]["similarity"] = hit["similarity"]

    results = list(merged.values())
    results.sort(key=lambda x: x["similarity"], reverse=True)

    return results
