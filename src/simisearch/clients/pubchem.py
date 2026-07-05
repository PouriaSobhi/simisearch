"""
clients/pubchem.py

Client for the PubChem PUG REST similarity search endpoint.

Endpoint:  https://pubchem.ncbi.nlm.nih.gov/rest/pug
No authentication required. Rate limit: max 5 requests/second.

Note: the parent LeadReplacementEngine pipeline also uses PubChem for
target lookup (CID -> human protein targets via GeneID xrefs). That's
a mechanism-of-action concern, not a structural-similarity one, so
it's intentionally left out of this spinoff - simisearch only returns
structurally similar compounds.
"""

import time
import requests


PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

DEFAULT_THRESHOLD = 85       # Tanimoto similarity, 0-100 scale
DEFAULT_MAX_RECORDS = 20
DEFAULT_POLL_INTERVAL = 3.0
DEFAULT_POLL_MAX_TRIES = 20  # ~60 seconds total

_similarity_cache = {}


def _get(url, params=None):
    """GET with error handling; returns parsed JSON or None."""
    try:
        response = requests.get(url, params=params, timeout=30)
        if response.status_code == 404:
            return None
        response.raise_for_status()
        return response.json()
    except Exception:
        return None


def similarity_search(
    smiles,
    threshold=DEFAULT_THRESHOLD,
    max_records=DEFAULT_MAX_RECORDS,
    poll_interval=DEFAULT_POLL_INTERVAL,
    poll_max_tries=DEFAULT_POLL_MAX_TRIES,
):
    """
    Find compounds in PubChem similar to the query SMILES.

    Parameters
    ----------
    smiles      : str - query SMILES (submitted as POST body, not URL,
                        to handle "/" and other special characters safely)
    threshold   : int - Tanimoto similarity threshold 0-100 (default 85)
    max_records : int - maximum results (default 20)

    Returns
    -------
    list of dicts with:
        "cid"        : int   - PubChem CID
        "smiles"     : str   - IsomericSMILES
        "name"       : str   - IUPACName (may be empty)
        "similarity" : float - threshold/100.0

    Returns empty list on timeout or failure.
    """

    cache_key = f"{smiles}|{threshold}|{max_records}"
    if cache_key in _similarity_cache:
        return _similarity_cache[cache_key]

    # Step 1: submit async job
    url = (
        f"{PUG_BASE}/compound/similarity/smiles/JSON"
        f"?Threshold={threshold}&MaxRecords={max_records}"
    )

    try:
        response = requests.post(
            url,
            data={"smiles": smiles},
            timeout=30,
        )
        response.raise_for_status()
        data = response.json()
    except Exception:
        _similarity_cache[cache_key] = []
        return []

    list_key = (data.get("Waiting") or {}).get("ListKey")
    if not list_key:
        _similarity_cache[cache_key] = []
        return []

    # Step 2: poll until done
    cids = []
    for _ in range(poll_max_tries):
        time.sleep(poll_interval)
        poll_data = _get(
            f"{PUG_BASE}/compound/listkey/{list_key}/cids/JSON"
        )
        if poll_data is None:
            continue
        if "Waiting" in poll_data:
            continue
        cids = (poll_data.get("IdentifierList") or {}).get("CID") or []
        break

    if not cids:
        _similarity_cache[cache_key] = []
        return []

    # Step 3: fetch SMILES + name per CID individually
    results = []
    for cid in cids[:max_records]:
        prop_data = _get(
            f"{PUG_BASE}/compound/cid/{cid}"
            f"/property/IsomericSMILES,IUPACName/JSON"
        )
        if prop_data:
            props = (
                prop_data
                .get("PropertyTable", {})
                .get("Properties", [{}])[0]
            )
            results.append({
                "cid": cid,
                "smiles": props.get("IsomericSMILES", props.get("SMILES", "")),
                "name": props.get("IUPACName", ""),
                "similarity": round(threshold / 100.0, 2),
            })
        time.sleep(0.25)

    _similarity_cache[cache_key] = results
    return results
