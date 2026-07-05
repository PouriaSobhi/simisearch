# simisearch

A small Python library that takes a query SMILES and returns a ranked set of
structurally similar compounds, using two independent chemical similarity
sources:

- **[SwissSimilarity](https://www.swisssimilarity.ch)** — screens real
  compound libraries (DrugBank, ChEMBL approved drugs, full ChEMBL bioactive
  set, and others) via its command-line API.
- **[PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)** —
  asynchronous 2D similarity search over the PubChem Compound database.

Results from both sources are deduplicated by canonical SMILES (via RDKit) and
merged into one ranked list, so the same molecule found under a DrugBank ID
and a PubChem CID shows up once, with its best similarity score and all
originating sources recorded.

This is a focused spinoff of a larger drug-discovery pipeline
(`LeadReplacementEngine`) — it only does structural similarity search, nothing
about targets, pathways, or ADME/toxicity filtering.

> **Disclaimer:** This is a personal hobby project, free to use for personal,
> non-commercial purposes only (see [LICENSE](LICENSE)). The author and this
> project are not affiliated with, endorsed by, or sponsored by
> SwissSimilarity, the Swiss Institute of Bioinformatics, PubChem, NCBI,
> DrugBank, ChEMBL/EMBL-EBI, or any other organization whose public API this
> project connects to. All trademarks belong to their respective owners. Use
> of these third-party services is subject to each provider's own terms of
> service and rate limits.

## Install

```bash
git clone https://github.com/PouriaSobhi/simisearch.git
cd simisearch
pip install -e .
```

Requires Python 3.9+. Dependencies (`requests`, `rdkit`) install automatically.

## Quick start

```python
from simisearch import search

aspirin = "CC(=O)Oc1ccccc1C(=O)O"
results = search(aspirin)

for hit in results[:5]:
    print(f"{hit['similarity']:.2f}  {hit['smiles']}")
```

Each result looks like:

```python
{
    "smiles": "Cc1ccc(cc1)...",     # canonical SMILES
    "similarity": 0.87,             # best score seen across sources
    "sources": [
        {"provider": "swiss_similarity", "similarity": 0.87,
         "libraries": [{"id": "DB00945", "library": "DrugBank", "similarity": 0.87}]},
        {"provider": "pubchem", "cid": 2244, "name": "aspirin", "similarity": 0.90},
    ],
}
```

## Setting parameters yourself

All parameters — similarity thresholds, which libraries to screen, how many
PubChem records to pull back, whether to use one provider or both — are set
through plain Python config objects. No config files, no CLI flags to
memorize.

```python
from simisearch import search, SearchConfig, SwissSimilarityConfig, PubChemConfig

config = SearchConfig(
    swiss=SwissSimilarityConfig(
        similarity_threshold=0.8,          # default 0.7
        drug_libraries=["DrugBank"],       # default ["DrugBank", "CHEMBL_drug"]
        bioactive_libraries=[],            # default ["CHEMBL"]; [] to skip
        method="ECFP",                     # default "Combined"
    ),
    pubchem=PubChemConfig(
        threshold=90,                      # 0-100 Tanimoto scale, default 85
        max_records=10,                    # default 20
    ),
)

results = search("CC(=O)Oc1ccccc1C(=O)O", config=config)
```

Use only one provider:

```python
# PubChem only — much faster, since SwissSimilarity is a polling job
# against a remote queue that can take minutes per library.
config = SearchConfig(swiss=SwissSimilarityConfig(enabled=False))
results = search(aspirin, config=config)
```

Reach for a provider directly if you want raw, non-merged output:

```python
from simisearch.clients import pubchem

hits = pubchem.similarity_search(aspirin, threshold=90, max_records=10)
```

See [`examples/basic_usage.py`](examples/basic_usage.py) for more.

## Configuration reference

**`SwissSimilarityConfig`**

| field | default | meaning |
|---|---|---|
| `enabled` | `True` | whether this provider runs |
| `method` | `"Combined"` | screening method (`ECFP`, `FP2`, `MHFP6`, `pharmacophore`, `scaffold`, ...) |
| `similarity_threshold` | `0.7` | minimum similarity (0-1) to keep a hit |
| `drug_libraries` | `["DrugBank", "CHEMBL_drug"]` | libraries in the "Drugs" family to screen |
| `bioactive_libraries` | `["CHEMBL"]` | libraries in the "Bioactive" family to screen |
| `poll_interval_seconds` / `poll_max_tries` | `5.0` / `60` | polling behaviour while SwissSimilarity computes results |

**`PubChemConfig`**

| field | default | meaning |
|---|---|---|
| `enabled` | `True` | whether this provider runs |
| `threshold` | `85` | Tanimoto similarity, 0-100 scale |
| `max_records` | `20` | max number of hits to retrieve |
| `poll_interval_seconds` / `poll_max_tries` | `3.0` / `20` | polling behaviour for the async PubChem job |

**`SearchConfig`**

| field | default | meaning |
|---|---|---|
| `swiss` / `pubchem` | see above | per-provider config |
| `dedupe` | `True` | merge hits from both providers that share a canonical SMILES |
| `sort_descending` | `True` | sort merged results by similarity |

## How it works

1. `search()` sends the query SMILES to every enabled provider.
2. Each provider's hits are canonicalized with RDKit (`Chem.MolToSmiles`).
3. Hits sharing a canonical SMILES are merged into one entry: the highest
   similarity score is kept, and every contributing provider/library/ID is
   recorded under `sources`.
4. The merged list is sorted by similarity, descending.

Canonical-SMILES matching is necessary because DrugBank, ChEMBL, and PubChem
each use their own ID system for what may be the same molecule — matching on
structure, not ID, is what makes cross-provider dedup possible.

## Notes on the providers

- **SwissSimilarity** is a polling job: submit → poll `checksession` until
  finished → retrieve results. Each library takes anywhere from 30 seconds to
  a few minutes depending on load. Screening the default 3 libraries can take
  2–15 minutes total.
- **PubChem** similarity search is also asynchronous (submit → poll `listkey`
  → fetch properties per CID) but is generally faster, and is rate-limited to
  5 requests/second.
- Both clients cache results in memory per process (keyed on the query and
  parameters), so repeated calls with identical parameters are free.

## Running tests

```bash
pip install -e ".[dev]"
pytest
```

The test suite only covers logic that doesn't require live network access
(canonical SMILES handling, merge/dedup behavior via mocked providers). It
does not hit SwissSimilarity or PubChem.

## License

Free for personal, non-commercial use. See [LICENSE](LICENSE) for the full
terms and the third-party affiliation disclaimer.
