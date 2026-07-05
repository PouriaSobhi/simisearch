"""
config.py

Plain-Python configuration objects for simisearch. Every parameter has
a sensible default (matching what was validated in the parent
LeadReplacementEngine pipeline), so `SearchConfig()` works out of the
box, but every value can be overridden.

Example
-------
    from simisearch import SearchConfig, SwissSimilarityConfig, PubChemConfig

    config = SearchConfig(
        swiss=SwissSimilarityConfig(
            similarity_threshold=0.8,
            drug_libraries=["DrugBank"],
            bioactive_libraries=[],
        ),
        pubchem=PubChemConfig(threshold=90, max_records=10),
    )
"""

from dataclasses import dataclass, field
from typing import List


@dataclass
class SwissSimilarityConfig:
    """Settings for the SwissSimilarity command-line API provider."""

    enabled: bool = True

    # One of the SwissSimilarity screening methods, e.g. "Combined",
    # "ECFP", "FP2", "MHFP6", "pharmacophore", "scaffold", ...
    method: str = "Combined"

    # Minimum similarity score (0-1) required to keep a hit.
    similarity_threshold: float = 0.7

    # Library families screened. Defaults are a lighter representative
    # subset (not all 8 possible libraries) to keep runtime manageable.
    drug_libraries: List[str] = field(default_factory=lambda: ["DrugBank", "CHEMBL_drug"])
    bioactive_libraries: List[str] = field(default_factory=lambda: ["CHEMBL"])

    # Polling behaviour while SwissSimilarity computes results.
    poll_interval_seconds: float = 5.0
    poll_max_tries: int = 60  # ~5 minutes max wait per library


@dataclass
class PubChemConfig:
    """Settings for the PubChem PUG REST similarity search provider."""

    enabled: bool = True

    # Tanimoto similarity threshold, 0-100 (PubChem's own scale).
    threshold: int = 85

    # Max number of similar compounds to retrieve.
    max_records: int = 20

    # Polling behaviour while PubChem computes the async similarity job.
    poll_interval_seconds: float = 3.0
    poll_max_tries: int = 20


@dataclass
class SearchConfig:
    """Top-level configuration passed to simisearch.search()."""

    swiss: SwissSimilarityConfig = field(default_factory=SwissSimilarityConfig)
    pubchem: PubChemConfig = field(default_factory=PubChemConfig)

    # Deduplicate hits found by multiple providers via canonical SMILES.
    dedupe: bool = True

    # Sort merged results by similarity, descending.
    sort_descending: bool = True
