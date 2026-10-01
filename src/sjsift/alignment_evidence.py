"""Extract bounded structural evidence independently of STAR quantification."""

from array import array
from bisect import bisect_left
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import heapq
import os
import tempfile

import pysam

from .catalog import Catalog


MAPPING_CLASSES = ("unique", "multimapping", "unknown")


class AlignmentError(ValueError):
    """A user-correctable alignment input failure."""


@dataclass(frozen=True)
class CigarOperation:
    """One CIGAR operation at a zero-based, half-open reference interval.

    Non-reference-consuming operations have start == end and retain their length.
    """

    op: str
    start: int
    end: int
    length: int


@dataclass(frozen=True)
class AlignmentRecord:
    record_id: str
    read_name: str
    read_group: str
    mate: str
    chromosome: str
    start: int
    end: int
    cigar: str
    mapq: int
    orientation: str
    flag: int
    flags: tuple[str, ...]
    mapping_class: str
    nh: int | None
    strand_status: str
    strand_evidence: tuple[str, ...]
    geometry: tuple[CigarOperation, ...]
    matches: tuple[str, ...]
    anchors: tuple[int, int]


@dataclass
class JunctionEvidence:
    label: str
    eligible: dict[str, int] = field(default_factory=lambda: dict.fromkeys(MAPPING_CLASSES, 0))
    strand_totals: dict[str, int] = field(default_factory=lambda: {"agreeing": 0, "unverified": 0})
    exclusions: dict[str, int] = field(default_factory=lambda: {"opposite": 0, "conflicting": 0})
    records: tuple[AlignmentRecord, ...] = ()


@dataclass(frozen=True)
class AlignmentEvidence:
    groups: dict[str, JunctionEvidence]
    provenance: dict[str, str]
    limit: int


def _geometry(read: pysam.AlignedSegment) -> tuple[CigarOperation, ...]:
    cursor = read.reference_start
    operations = []
    for code, length in read.cigartuples or ():
        if code not in range(9) or length <= 0:
            raise AlignmentError("invalid or unsupported CIGAR operation")
        op = "MIDNSHP=X"[code]
        end = cursor + length if op in "MDN=X" else cursor
        operations.append(CigarOperation(op, cursor, end, length))
        cursor = end
    return tuple(operations)


def _tag(read: pysam.AlignedSegment, name: str):
    return read.get_tag(name, with_value_type=True) if read.has_tag(name) else (None, None)


def _mapping_class(read: pysam.AlignedSegment) -> str:
    value, kind = _tag(read, "NH")
    if kind in ("c", "C", "s", "S", "i", "I") and type(value) is int and value > 0:
        return "unique" if value == 1 else "multimapping"
    return "unknown"


def _program_provenance(read, header) -> set[str]:
    programs = {p["ID"]: p for p in header.get("PG", [])}
    pg, kind = _tag(read, "PG")
    if pg is not None:
        leaves = [pg] if kind == "Z" and pg in programs else []
    else:
        parents = {p.get("PP") for p in programs.values()}
        leaves = [key for key in programs if key not in parents]
    chains = []
    for leaf in leaves:
        names, visited = set(), set()
        while leaf in programs and leaf not in visited:
            visited.add(leaf)
            program = programs[leaf]
            names.add(program.get("PN", program["ID"]).lower())
            leaf = program.get("PP")
        chains.append(names)
    # Without a record PG, only provenance shared by every possible chain is usable.
    return set.intersection(*chains) if chains else set()


def _integer_array(value, kind):
    return kind in ("Bc", "BC", "Bs", "BS", "Bi", "BI") and isinstance(value, array) and value.typecode in "bBhHiIlL"


def _strand(read, introns, intron_index, strand, header):
    sources, notes = [], []
    for name in ("TS", "XS"):
        value, kind = _tag(read, name)
        if kind == "A" and value in ("+", "-"):
            sources.append((f"{name}:A:{value}", value))
        elif value is not None:
            notes.append(f"{name}: unusable type or value")
    programs = _program_provenance(read, header)
    value, kind = _tag(read, "ts")
    if value is not None:
        if "minimap2" in programs and kind == "A" and value in ("+", "-"):
            direction = ("-" if value == "+" else "+") if read.is_reverse else value
            sources.append((f"minimap2 ts:A:{value} interpreted relative to alignment: {direction}", direction))
        else:
            notes.append("ts: unusable value or unverified minimap2 provenance")
    motifs, kind = _tag(read, "jM")
    if motifs is not None:
        coordinates, coordinate_kind = _tag(read, "jI")
        valid = (
            "star" in programs and _integer_array(motifs, kind)
            and len(motifs) == len(introns)
            and all(m in {*range(7), *range(20, 27)} for m in motifs)
            and (coordinates is None or (
                _integer_array(coordinates, coordinate_kind)
                and list(coordinates) == [coordinate for pair in introns for coordinate in pair]
            ))
        )
        if valid:
            motif = motifs[intron_index]
            code = motif - 20 if motif >= 20 else motif
            if code:
                direction = "+" if code % 2 else "-"
                sources.append((f"STAR jM[{intron_index + 1}]={motif}: {direction}; CIGAR order validated", direction))
            else:
                notes.append(f"STAR jM[{intron_index + 1}]={motif}: noncanonical, strand unverified")
        else:
            notes.append("jM: unusable motif/coordinates or unverified STAR provenance")
    values = {value for _, value in sources}
    status = "conflicting" if len(values) > 1 else (
        "opposite" if values and strand not in values else "agreeing" if values else "unverified"
    )
    return status, tuple([source for source, _ in sources] + notes)


def _retain(heap, record, key, limit):
    # Stable bottom-k hash sample, bounded during extraction, independently per class.
    rank = int.from_bytes(hashlib.sha256(f"{key}/{record.record_id}".encode()).digest(), "big")
    item = (-rank, record.record_id, record)
    if len(heap) < limit:
        heapq.heappush(heap, item)
    elif item > heap[0]:
        heapq.heapreplace(heap, item)


def _anchors(geometry, intron_index):
    """Contiguous M/= /X lengths immediately left/right of the selected N."""
    position = [i for i, op in enumerate(geometry) if op.op == "N"][intron_index]
    counts = []
    for step in (-1, 1):
        cursor, count = position + step, 0
        while 0 <= cursor < len(geometry) and geometry[cursor].op in "M=X":
            count += geometry[cursor].length
            cursor += step
        counts.append(count)
    return tuple(counts)


def _serialize(read, record_id, geometry, mapping_class, status, evidence, matches, intron_index):
    flags = tuple(label for mask, label in (
        (1, "paired"), (2, "proper pair"), (8, "mate unmapped"),
        (256, "secondary"), (512, "QC failure"), (1024, "duplicate"), (2048, "supplementary"),
    ) if read.flag & mask)
    mate = "read 1" if read.is_read1 and not read.is_read2 else (
        "read 2" if read.is_read2 and not read.is_read1 else "unspecified"
    )
    rg, kind = _tag(read, "RG")
    nh, nh_kind = _tag(read, "NH")
    nh = nh if mapping_class != "unknown" else None
    return AlignmentRecord(
        record_id, read.query_name or "*", rg if kind == "Z" else "", mate,
        read.reference_name, read.reference_start, read.reference_end,
        read.cigarstring, read.mapping_quality, "reverse" if read.is_reverse else "forward",
        read.flag, flags, mapping_class, nh, status, evidence, geometry, matches,
        _anchors(geometry, intron_index),
    )


def _extract_evidence(
    catalog: Catalog, path: Path, *, index_path: Path, limit: int = 200,
    reference_path: Path | None = None,
) -> AlignmentEvidence:
    """Return selected-junction record evidence; never modify STAR support."""
    targets = []
    groups = {}
    for i, variant in enumerate(catalog.variants):
        for position, junction in enumerate((variant, *variant.reference_junctions)):
            key = f"v{i}-j{position}"
            label = "Defining junction" if position == 0 else junction.role
            targets.append((key, junction))
            groups[key] = JunctionEvidence(f"{variant.identifier} · {label}")
    retained = {key: {kind: [] for kind in MAPPING_CLASSES} for key, _ in targets}
    with pysam.AlignmentFile(str(path), "rc" if reference_path else "rb", index_filename=str(index_path),
                             reference_filename=str(reference_path) if reference_path else None,
                             require_index=True) as stream:
        if not (stream.is_bam or stream.is_cram):
            raise AlignmentError("expected indexed BAM or CRAM")
        stream.check_index()
        header = stream.header.to_dict()
        provenance = _validate_header(catalog, header)
        if reference_path:
            provenance.update(_validate_reference(header, reference_path))
        for _, junction in targets:
            if junction.chromosome not in stream.references:
                raise AlignmentError(f"missing contig {junction.chromosome!r}")
            if junction.intron_end > stream.get_reference_length(junction.chromosome):
                raise AlignmentError(f"junction exceeds contig bounds: {junction.chromosome}:{junction.intron_end}")
        by_chromosome = {}
        for key, junction in targets:
            by_chromosome.setdefault(junction.chromosome, []).append((key, junction))
        for chromosome, selected in by_chromosome.items():
            positions = sorted({j.intron_start - 1 for _, j in selected})
            for region, position in enumerate(positions):
                previous_start = -1
                for ordinal, read in enumerate(stream.fetch(chromosome, position, position + 1)):
                    if read.reference_start < previous_start:
                        raise AlignmentError("records are not coordinate sorted")
                    previous_start = read.reference_start
                    if read.reference_end is not None and read.reference_end > stream.get_reference_length(chromosome):
                        raise AlignmentError("alignment exceeds contig bounds")
                    rg, kind = _tag(read, "RG")
                    if rg is not None and (kind != "Z" or rg not in {r["ID"] for r in header.get("RG", [])}):
                        raise AlignmentError("alignment has an unknown or invalid read group")
                    # Each physical record belongs to the first requested position it overlaps.
                    # This avoids repeat-fetch duplication without collapsing identical records.
                    if read.is_unmapped or bisect_left(positions, read.reference_start) != region:
                        continue
                    geometry = _geometry(read)
                    introns = [(op.start + 1, op.end) for op in geometry if op.op == "N"]
                    memberships = []
                    for key, junction in selected:
                        pair = (junction.intron_start, junction.intron_end)
                        if pair not in introns:
                            continue
                        status, evidence = _strand(read, introns, introns.index(pair), junction.strand, header)
                        group = groups[key]
                        if status in group.exclusions:
                            group.exclusions[status] += 1
                            continue
                        memberships.append((key, status, evidence, introns.index(pair)))
                    mapping_class = _mapping_class(read)
                    matches = tuple(key for key, _, _, _ in memberships)
                    for key, status, evidence, intron_index in memberships:
                        group = groups[key]
                        group.eligible[mapping_class] += 1
                        group.strand_totals[status] += 1
                        record_id = f"r{read.reference_id}-{region}-{ordinal}"
                        record = _serialize(read, record_id, geometry, mapping_class, status, evidence, matches, intron_index)
                        _retain(retained[key][mapping_class], record, key, limit)
    for key, group in groups.items():
        group.records = tuple(
            item[2] for kind in MAPPING_CLASSES
            for item in sorted(retained[key][kind], key=lambda item: item[1])
        )
    provenance.update({"Alignment file": path.name, "Format": "CRAM" if reference_path else "BAM", "pysam": pysam.__version__,
                       "HTSlib": pysam.__samtools_version__, "Index": index_path.name})
    return AlignmentEvidence(groups, provenance, limit)


def _validate_header(catalog, header):
    if header.get("HD", {}).get("SO") != "coordinate":
        raise AlignmentError("alignment header must declare coordinate sort order")
    read_groups = header.get("RG", [])
    identifiers = [r.get("ID") for r in read_groups]
    if None in identifiers or len(set(identifiers)) != len(identifiers):
        raise AlignmentError("missing or duplicate read-group identifiers")
    samples = {r["SM"] for r in read_groups if r.get("SM")}
    if len(samples) > 1:
        raise AlignmentError("multiple samples in alignment header; only one sample is supported")
    required = {j.chromosome for v in catalog.variants for j in (v, *v.reference_junctions)}
    assemblies = {sq["AS"] for sq in header.get("SQ", []) if sq["SN"] in required and sq.get("AS")}
    if assemblies and assemblies != {catalog.genome_assembly}:
        raise AlignmentError("alignment assembly metadata disagrees with catalog assembly")
    return {
        "Alignment sample": next(iter(samples)) if samples else "Unavailable: no sample metadata",
        "Sample assurance": "STAR sample identity cannot be verified; missing read-group/sample metadata may limit alignment assurance",
        "Assembly assurance": "Available assembly labels agree; sequence equivalence and STAR assembly remain unverified" if assemblies else "Assembly unverified: contig names and bounds alone do not establish equivalence",
    }


@contextmanager
def _offline_reference():
    """The single-command process must never inherit a remote refget search path."""
    previous = {key: os.environ.get(key) for key in ("REF_PATH", "REF_CACHE")}
    with tempfile.TemporaryDirectory(prefix="sjsift-reference-") as directory:
        try:
            for key in previous:
                os.environ[key] = str(Path(directory) / "%s")
            yield
        finally:
            for key, value in previous.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


def _check_reference_files(path):
    if not path.is_file():
        raise AlignmentError("--reference must be a local FASTA file")
    index = Path(str(path) + ".fai")
    if not index.is_file():
        raise AlignmentError("reference requires an existing local .fai index; no index is created")
    with path.open("rb") as source:
        compressed = source.read(2) == b"\x1f\x8b"
    if compressed and not Path(str(path) + ".gzi").is_file():
        raise AlignmentError("compressed reference requires an existing local .gzi index")
    # Explicit index paths prevent pysam's automatic index-building fallback.
    with pysam.FastaFile(str(path), filepath_index=str(index),
                         filepath_index_compressed=str(path) + ".gzi" if compressed else None):
        pass


def _validate_reference(header, path):
    checked = 0
    with pysam.FastaFile(str(path), filepath_index=str(path) + ".fai") as reference:
        for sequence in header.get("SQ", []):
            name, length = sequence["SN"], sequence["LN"]
            if name not in reference.references:
                raise AlignmentError(f"reference is missing alignment contig {name!r}")
            if reference.get_reference_length(name) != length:
                raise AlignmentError(f"reference length disagrees with alignment contig {name!r}")
            digest = hashlib.md5()
            # Bound validation memory even for whole human chromosomes.
            for start in range(0, length, 1024 * 1024):
                end = min(length, start + 1024 * 1024)
                bases = reference.fetch(name, start, end)
                if len(bases) != end - start:
                    raise AlignmentError(f"reference sequence/index is truncated or inconsistent for {name!r}")
                digest.update(bases.upper().encode("ascii"))
            if sequence.get("M5"):
                if digest.hexdigest() != sequence["M5"].lower():
                    raise AlignmentError(f"reference M5 checksum disagrees for {name!r}")
                checked += 1
    return {
        "Reference FASTA": path.name,
        "Reference index": path.name + ".fai",
        "Reference assurance": f"Local contig names/lengths checked; {checked} of {len(header.get('SQ', []))} header M5 checksums verified. Missing M5 limits sequence assurance; decoded slice checks remain enabled. STAR assembly remains unverified.",
        "Reference access": "Explicit local FASTA only; inherited REF_PATH/REF_CACHE isolated; no remote reference lookup",
    }


def extract_evidence(
    catalog: Catalog, path: Path, *, index_path: Path | None = None, limit: int = 200,
    reference_path: Path | None = None,
) -> AlignmentEvidence:
    """Validate local indexed BAM/CRAM and return evidence before output creation.

    Regional access counts alignment records, not templates. Every requested
    region is decoded to completion even after the independent sample caps fill.
    """
    try:
        if type(limit) is not int or limit < 1:
            raise AlignmentError("alignment limit must be a positive integer")
        if not path.is_file():
            raise AlignmentError("alignment input must be a local file")
        with path.open("rb") as source:
            is_cram = source.read(4) == b"CRAM"
        if is_cram and reference_path is None:
            raise AlignmentError("CRAM requires --reference with an explicit local indexed FASTA")
        if not is_cram and reference_path is not None:
            raise AlignmentError("--reference applies only to CRAM input")
        if index_path is None:
            candidates = (Path(str(path) + ".crai"), path.with_suffix(".crai")) if is_cram else (
                Path(str(path) + ".bai"), path.with_suffix(".bai"),
                Path(str(path) + ".csi"), path.with_suffix(".csi"),
            )
            index_path = next((candidate for candidate in candidates if candidate.is_file()), None)
        if index_path is None or not index_path.is_file():
            raise AlignmentError("a local alignment index is required; no index is created automatically")
        with _offline_reference() if is_cram else nullcontext():
            if reference_path:
                _check_reference_files(reference_path)
                # AlignmentFile.check_truncation does not reject missing CRAM EOF.
                # Use pysam's bundled header/EOF check; no external executable.
                pysam.quickcheck("--", str(path))
            return _extract_evidence(catalog, path, index_path=index_path, limit=limit, reference_path=reference_path)
    except (OSError, ValueError, KeyError, IndexError, OverflowError, pysam.SamtoolsError) as error:
        raise AlignmentError(f"{path}: cannot extract alignment evidence: {error}") from None
