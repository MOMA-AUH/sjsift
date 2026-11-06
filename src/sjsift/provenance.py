"""Build a deterministic, bounded-metadata description of one completed run."""

import hashlib
import json
from pathlib import Path

from . import __version__
from .alignment_evidence import AlignmentEvidence
from .catalog import Catalog


COMPATIBILITY_WARNING = {
    "code": "NO_CATALOG_CHROMOSOMES",
    "message": "no catalog chromosome identifiers found; check genome assembly and chromosome naming compatibility",
}


def _fingerprint(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return {"basename": path.name, "sha256": digest.hexdigest(),
            "fingerprint_scope": "exact_file_bytes"}


def render_provenance(
    catalog: Catalog, definitions_path: Path, junctions_path: Path, *,
    sample_name: str,
    sample_source: str = "explicit",
    compatibility_warning: bool = False,
    alignment_evidence: AlignmentEvidence | None = None,
) -> str:
    """Return schema-1 JSON, without output paths, timestamps or record data.

    Inputs must remain unchanged during a run. Fingerprints cover the exact
    file bytes, including the gzip wrapper when STAR input is compressed.
    """
    catalog_input = {"schema_version": catalog.schema_version, **_fingerprint(definitions_path)}
    star_input = _fingerprint(junctions_path)
    with junctions_path.open("rb") as stream:
        star_input["compression"] = "gzip" if stream.read(2) == b"\x1f\x8b" else "none"
    document = {
        "schema_version": 1,
        "sample": {"name": sample_name, "source": sample_source},
        "software": {"name": "sjsift", "version": __version__},
        "genome_assembly": catalog.genome_assembly,
        "catalog": catalog_input,
        "star_input": star_input,
        "warnings": [COMPATIBILITY_WARNING] if compatibility_warning else [],
        "alignment": {"status": "not_requested"},
    }
    if alignment_evidence is not None:
        document["alignment"] = {
            "status": "requested",
            **alignment_evidence.metadata,
            "eligible_junction_memberships": sum(
                sum(group.eligible.values()) for group in alignment_evidence.groups.values()
            ),
        }
    return json.dumps(document, indent=2, ensure_ascii=True, allow_nan=False) + "\n"
