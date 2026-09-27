"""Schema constants, evidence classes and provenance (layer L0, stdlib only).

Every value this package outputs carries one of five evidence classes. The
mapping from output field to class lives only in FIELD_CLASSES; output
documents copy their "field_classes" block from it, so a field that is not
registered here cannot be emitted (field_classes() raises KeyError).
"""

import hashlib
import sys
from collections.abc import Iterable, Mapping
from types import MappingProxyType

SCHEMA = "protein-hints/1"
TOOL_VERSION = "0.1.0"

EVIDENCE_CLASSES = (
    "sequence-derived",               # computed directly from the sequence
    "heuristic",                      # rule-based candidates, no physical claim
    "predicted",                      # statistical or propensity model
    "structure-derived-experimental", # measured on an X-ray / cryo-EM / NMR model
    "structure-derived-predicted",    # measured on a predicted model (e.g. AlphaFold)
)

FIELD_CLASSES = MappingProxyType({
    # Whole-sequence features (protein_features).
    "length": "sequence-derived",
    "composition": "sequence-derived",
    "fraction_aromatic": "sequence-derived",
    "fraction_sulfur": "sequence-derived",
    "fraction_pro_gly": "sequence-derived",
    "molecular_weight_avg": "sequence-derived",
    "net_charge_simplified": "sequence-derived",
    "cys_positions": "sequence-derived",
    "pro_positions": "sequence-derived",
    "gly_positions": "sequence-derived",
    "aromatic_positions": "sequence-derived",
    # Per-residue values (protein_profiles).
    "hydropathy": "sequence-derived",
    "charge": "sequence-derived",
    "hydropathy_window": "sequence-derived",
    "charge_window": "sequence-derived",
    "entropy_window": "sequence-derived",
    # Regions of interest (regions_of_interest).
    "regions": "heuristic",
    "summary": "heuristic",
})

# Per-residue keys that locate or label a residue rather than describe it.
ROW_KEYS = ("position", "input_position", "residue", "class")


def sequence_sha256(residues: str) -> str:
    """sha256 of the prepared residue string (upper case, ASCII)."""
    return hashlib.sha256(residues.encode("ascii")).hexdigest()


def field_classes(names: Iterable[str]) -> dict[str, str]:
    """The FIELD_CLASSES entries for `names`; KeyError for an unregistered field."""
    return {name: FIELD_CLASSES[name] for name in sorted(set(names))}


def provenance(*, table_version: str, params: Mapping[str, object]) -> dict:
    """Provenance block: tool, interpreter, property table and parameters."""
    return {
        "tool_version": TOOL_VERSION,
        "python": "{}.{}.{}".format(*sys.version_info[:3]),
        "table_version": table_version,
        "params": dict(params),
        "extras": {},
    }
