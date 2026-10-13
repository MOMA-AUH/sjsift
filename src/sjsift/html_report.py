"""Render a self-contained, offline overview and per-variant evidence report."""

from collections.abc import Iterable
from dataclasses import dataclass
from html import escape
from importlib.resources import files
from string import Template
from typing import TextIO

from . import __version__
from .alignment_evidence import AlignmentEvidence, AlignmentRecord, MAPPING_CLASSES
from .quantify import ReferenceJunctionSupport, VariantSupport


_ROLE_LABELS = {"same_donor": "Same donor", "same_acceptor": "Same acceptor"}


def _role_label(role: str) -> str:
    return _ROLE_LABELS.get(role, role)


def _observation(result: VariantSupport) -> tuple[str, str]:
    if not result.reference_junctions:
        label = "Variant support observed" if result.total else "No variant support"
        return (
            f"{label}; no reference context configured",
            "This catalog entry has no reference junctions. Reference evidence is unavailable.",
        )
    reference_seen = any(support.total for support in result.reference_junctions)
    if result.total and reference_seen:
        return (
            "Variant and reference support",
            "The defining junction and at least one selected reference have support. "
            "Compare each reference separately.",
        )
    if result.total:
        return (
            "Variant support only",
            "The defining junction has support; all selected references have zero reported support.",
        )
    if reference_seen:
        return (
            "Reference support only",
            "At least one selected reference has support; the defining junction has zero reported support.",
        )
    return (
        "No selected junction support",
        "The defining junction and all selected references have zero reported support in this STAR file.",
    )


def _counts(support: VariantSupport | ReferenceJunctionSupport) -> str:
    return (
        f'<span class="counts"><strong>{support.unique}</strong> unique'
        f'<span class="multi-count"> + {support.multimapping} multi</span></span>'
    )


def _bar(
    support: VariantSupport | ReferenceJunctionSupport,
    maximum: int,
    reference: bool = False,
) -> str:
    unique_width = 100 * (support.unique / maximum)
    multi_width = 100 * (support.multimapping / maximum)
    kind = "reference" if reference else "variant"
    return (
        '<span class="bar" aria-hidden="true">'
        f'<span class="{kind}" style="width:{unique_width:.6f}%"></span>'
        f'<span class="multimapping" style="width:{multi_width:.6f}%"></span></span>'
    )


def _coordinates(support: VariantSupport | ReferenceJunctionSupport) -> str:
    definition = support.definition
    return escape(
        f"{definition.chromosome}:{definition.intron_start}–{definition.intron_end} "
        f"({definition.strand})"
    )


def _overview_row(
    index: int, result: VariantSupport, roles: list[str], maximum: int
) -> str:
    definition = result.definition
    references = {support.definition.role: support for support in result.reference_junctions}
    cells = []
    for role in roles:
        reference = references.get(role)
        content = (
            _counts(reference) + _bar(reference, maximum, reference=True)
            if reference is not None else '<span class="muted">Not configured</span>'
        )
        cells.append(f"<td>{content}</td>")
    return (
        f'<tr data-index="{index}" data-unique="{result.unique}" data-total="{result.total}">'
        f'<th scope="row"><a href="#variant-{index}">{escape(definition.identifier)}</a>'
        f'<small>{escape(definition.chromosome)} · {escape(definition.strand)}</small></th>'
        f'<td>{_counts(result)}{_bar(result, maximum)}</td>'
        + "".join(cells)
        + f'<td class="observation">{escape(_observation(result)[0])}</td></tr>'
    )


def _schematic(label: str, support: VariantSupport | ReferenceJunctionSupport) -> str:
    junction = support.definition
    annotation = junction.annotation
    donor, acceptor = (junction.intron_start, junction.intron_end)
    if junction.strand == "-":
        donor, acceptor = acceptor, donor
    return (
        '<section class="junction-schematic">'
        f'<h3>{escape(label)}</h3>'
        '<p class="scale-note">Local splice boundaries · transcript direction → · not to scale. '
        'Exon extents and intervening transcript structure are not shown.</p>'
        '<div class="splice-boundaries">'
        f'<div class="exon-boundary donor">Donor exon {escape(annotation.donor_exon)}'
        f'<small>Donor boundary: {donor}</small></div>'
        f'<div class="splice-gap">↗ ··· ↘<small>{junction.intron_end - junction.intron_start + 1} nt intron</small></div>'
        f'<div class="exon-boundary acceptor">Acceptor exon {escape(annotation.acceptor_exon)}'
        f'<small>Acceptor boundary: {acceptor}</small></div></div>'
        f'<p class="coordinates">{_coordinates(support)} · boundaries are intronic bases</p>'
        f'<p>Annotation / comparison basis: {escape(annotation.reference_transcript)}. '
        'This accession supplies labels or a comparator; it need not contain the defining junction.</p>'
        f'<p class="annotation-source">{escape(annotation.annotation_source)}</p></section>'
    )


def _display_interval(record: AlignmentRecord, operation) -> tuple[int, int]:
    """Extend soft clips visually without assigning them reference coordinates."""
    if operation.op == "S":
        if operation.start == record.start:
            return operation.start - operation.length, operation.start
        return operation.start, operation.start + operation.length
    return operation.start, operation.end


@dataclass(frozen=True)
class _PileupScale:
    """One coordinate transform for all embedded records on a figure's chromosome."""

    start: int
    end: int
    breaks: tuple[tuple[int, int], ...]
    left_padding: int
    right_padding: int

    @classmethod
    def for_records(cls, records: tuple[AlignmentRecord, ...], boundaries: tuple[int, ...] = ()):
        start = min([record.start for record in records] + list(boundaries))
        end = max([record.end for record in records] + list(boundaries))
        # Preserve the base scale across aligned blocks and displayed soft clips.
        intervals = sorted([_display_interval(record, op) for record in records
                            for op in record.geometry if op.op not in "ND"]
                           + [(boundary, boundary) for boundary in boundaries])
        left_padding = max(0, start - min(left for left, _ in intervals))
        right_padding = max(0, max(right for _, right in intervals) - end)
        cursor, breaks = start, []
        for left, right in [*intervals, (end, end)]:
            if left - cursor > 80:
                breaks.append((cursor, left))
            cursor = max(cursor, right)
        return cls(start, end, tuple(breaks), left_padding, right_padding)

    def project(self, coordinate: float) -> float:
        return coordinate - self.start + self.left_padding - sum(
            max(0, min(coordinate, right) - left) * (1 - 24 / (right - left))
            for left, right in self.breaks
        )

    def percent(self, coordinate: float) -> float:
        # Leave room at both edges for fixed-width markers.
        return 1 + 98 * self.project(coordinate) / max(1, self.project(self.end) + self.right_padding)

    def ruler(self, chromosome: str) -> str:
        ticks = []
        # Label only the reference span; outer soft-clip space has no coordinates.
        width = self.project(self.end) - self.left_padding
        for index in range(5):
            projected = self.left_padding + width * index / 4
            coordinate = self.start + projected - self.left_padding
            for left, right in self.breaks:
                offset = max(0, min(24, projected - self.project(left)))
                coordinate += offset * ((right - left) / 24 - 1)
            base = min(self.end, int(coordinate) + 1)
            ticks.append(f'<span class="pileup-tick tick-{index}" style="left:{self.percent(coordinate):.6f}%">{base:,}</span>')
        breaks = "".join(
            f'<span class="gap-break" style="left:{self.percent((left + right) / 2):.6f}%" '
            f'title="Compressed gap: {left + 1:,}–{right:,} ({right - left:,} nt)">//</span>'
            for left, right in self.breaks
        )
        return (f'<div class="pileup-chromosome">{escape(chromosome)} · genomic coordinates →</div>'
                '<div class="pileup-ruler" aria-label="Shared genomic coordinate ruler">'
                + "".join(ticks) + breaks + '</div>')


def _read_geometry(record: AlignmentRecord, scale: _PileupScale, junctions=()) -> str:
    operations = []
    for operation in record.geometry:
        coordinates = f"{operation.start + 1}–{operation.end}" if operation.end > operation.start else f"boundary {operation.start}"
        description = f"{operation.op}: {operation.length} nt; {coordinates}"
        if operation.op == "S":
            description += "; unaligned soft clip"
        start, end = _display_interval(record, operation)
        left = scale.percent(start)
        width = scale.percent(end) - left
        selected = next((color for support, color in junctions
                         if operation.op == "N" and operation.start == support.definition.intron_start - 1
                         and operation.end == support.definition.intron_end), None)
        operations.append(
            f'<span class="cigar-op op-{operation.op.replace("=", "equal")}{" selected-junction" if selected else ""}" '
            f'style="left:{left:.6f}%;width:{width:.6f}%'
            + (f';background:{selected}' if selected else '') + f'" title="{description}"></span>'
        )
    return '<span class="read-geometry" aria-hidden="true">' + "".join(operations) + '</span>'


def _read_row(record: AlignmentRecord, position: int, also_matches: tuple[str, ...], scale: _PileupScale,
              memberships, bucket: int, junctions) -> str:
    flags = ", ".join(record.flags) or "primary; no additional flags"
    strand = "strand unverified" if record.strand_status == "unverified" else "strand evidence agrees"
    color = next(color for index, _, _, color in memberships if index == bucket)
    support_labels = "; ".join(label for _, label, _, _ in memberships)
    return (
        f'<li class="read-row {record.strand_status}" data-mapping="{record.mapping_class}"'
        f' data-bucket="{bucket}" data-support="{" ".join(str(index) for index, _, _, _ in memberships)}"'
        f' style="--read-color:{color}"'
        f' data-record-id="{record.record_id}"' + (' hidden' if position >= 10 else '') + '>'
        '<details class="read-inspection"><summary class="read-hit" '
        f'aria-label="Inspect {escape(record.read_name, quote=True)}; {escape(support_labels, quote=True)}; {record.mapping_class}; {strand}">'
        + _read_geometry(record, scale, junctions) + '</summary>'
        '<div class="read-metadata">'
        f'<div class="read-label"><strong>{escape(record.read_name)}</strong> · {record.mate} · '
        f'{record.mapping_class} · {strand}</div>'
        f'<p>Supports: {escape(support_labels)}</p>'
        + (f'<p class="shared-membership">Also matches: {escape("; ".join(also_matches))}</p>' if also_matches else "")
        + '<strong>Read details</strong>'
        f'<p>Record ID: {record.record_id} · Read group: {escape(record.read_group or "not supplied")} · '
        f'{escape(record.chromosome)}:{record.start + 1}–{record.end}</p>'
        f'<p>CIGAR: {escape(record.cigar)} · MAPQ: {record.mapq} · NH: {record.nh if record.nh is not None else "unavailable/unusable"} · Alignment orientation: {record.orientation}</p>'
        f'<p>Flags: {escape(flags)} ({record.flag}) · Transcript {strand}</p>'
        f'<p>Aligned anchors: genomic left {record.anchors[0]} nt; genomic right {record.anchors[1]} nt</p>'
        f'<p>Anchor and transcript-strand details refer to {escape(memberships[0][1])}.</p>'
        f'<p>Strand sources: {escape("; ".join(record.strand_evidence) or "No usable transcript-strand evidence")}</p>'
        '</div></details></li>'
    )


def _alignment_group(evidence: AlignmentEvidence | None, key: str) -> str:
    if evidence is None:
        return '<p class="alignment-unavailable">Alignment evidence not requested for this junction.</p>'
    group = evidence.groups[key]
    eligible = sum(group.eligible.values())
    counts = "".join(
        f'<tr><th scope="row">{kind}</th><td>{group.eligible[kind]}</td>'
        f'<td>{sum(r.mapping_class == kind for r in group.records)}</td></tr>'
        for kind in MAPPING_CLASSES
    )
    return (
        f'<section class="alignment-count-group" id="reads-{key}">'
        '<h4>Alignment records</h4>'
        f'<p>Eligible: {eligible} · Embedded: {len(group.records)}</p>'
        '<table class="alignment-counts"><caption>Alignment counts by mapping class</caption>'
        '<thead><tr><th>Class</th><th>Eligible</th><th>Embedded</th></tr></thead>'
        f'<tbody>{counts}</tbody></table>'
        f'<p>Transcript strand totals: agreeing {group.strand_totals["agreeing"]}; '
        f'unverified {group.strand_totals["unverified"]}. '
        f'Excluded: opposite {group.exclusions["opposite"]}; conflicting {group.exclusions["conflicting"]}.</p>'
        + ('<p>No eligible alignment records (0).</p>' if not eligible else '')
        + '</section>'
    )


def _junction_guide(label, support, color, scale):
    junction = support.definition
    annotation = junction.annotation
    left, right = scale.percent(junction.intron_start - 1), scale.percent(junction.intron_end)
    exons = (annotation.donor_exon, annotation.acceptor_exon)
    if junction.strand == "-":
        exons = exons[::-1]
    return (
        f'<div class="junction-guide" style="--read-color:{color}">'
        f'<div class="junction-guide-label"><span class="support-swatch" style="background:{color}"></span>'
        f'{escape(label)} · exon {escape(annotation.donor_exon)} → {escape(annotation.acceptor_exon)}'
        f' · {support.unique} unique + {support.multimapping} multi (STAR)</div>'
        '<div class="junction-guide-axis">'
        f'<span class="junction-guide-line" style="left:{left:.6f}%;width:{right - left:.6f}%" '
        f'title="{escape(junction.chromosome, quote=True)}:{junction.intron_start}–{junction.intron_end} ({junction.strand})">'
        f'{"→" if junction.strand == "+" else "←"}</span>'
        f'<span class="junction-endpoint left" style="left:{left:.6f}%">{escape(exons[0])}</span>'
        f'<span class="junction-endpoint right" style="left:{right:.6f}%">{escape(exons[1])}</span>'
        '</div></div>'
    )


def _combined_pileup(index: int, result: VariantSupport, evidence: AlignmentEvidence | None) -> str:
    if evidence is None:
        return '<p class="alignment-unavailable">Alignment evidence not requested.</p>'
    # Keys and colors are generated here; catalog text is never used in CSS or selectors.
    junctions = [("Variant", result, "#3159ad"), *(
        (_role_label(support.definition.role), support,
         {"same_donor": "#168375", "same_acceptor": "#9950a5"}.get(support.definition.role, "#677481"))
        for support in result.reference_junctions
    )]
    keys = [f"v{index}-j{position}" for position in range(len(junctions))]
    retained = {}
    for position, key in enumerate(keys):
        for record in evidence.groups[key].records:
            retained.setdefault(record.record_id, (record, position))
    buckets = {}
    for record, owner in retained.values():
        members = [(position, label, support, color) for position, (label, support, color) in enumerate(junctions)
                   if keys[position] in record.matches]
        bucket = min(position for position, _, _, _ in members)
        # Keep the metadata's extraction junction explicit for shared records.
        members.sort(key=lambda member: (member[0] != owner, member[0]))
        buckets.setdefault(bucket, []).append((record, members, owner))
    entries = []
    for bucket, rows in sorted(buckets.items()):
        rows.sort(key=lambda item: (item[0].chromosome, item[0].start, item[0].end, item[0].record_id))
        entries.extend((record, members, owner, bucket, position)
                       for position, (record, members, owner) in enumerate(rows))
    tracks = []
    chromosomes = list(dict.fromkeys(support.definition.chromosome for _, support, _ in junctions))
    for chromosome in chromosomes:
        local_junctions = [(label, support, color) for label, support, color in junctions
                           if support.definition.chromosome == chromosome]
        local_entries = [entry for entry in entries if entry[0].chromosome == chromosome]
        boundaries = tuple(coordinate for _, support, _ in local_junctions
                           for coordinate in (support.definition.intron_start - 1, support.definition.intron_end))
        scale = _PileupScale.for_records(tuple(entry[0] for entry in local_entries), boundaries)
        tracks.append('<div class="pileup-track">' + scale.ruler(chromosome)
                      + "".join(_junction_guide(label, support, color, scale)
                                for label, support, color in local_junctions)
                      + '<ol class="read-rows">' + "".join(
                          _read_row(record, position,
                                    tuple(evidence.groups[match].label for match in record.matches if match != keys[owner]),
                                    scale, members, bucket,
                                    tuple((support, color) for _, _, support, color in members))
                          for record, members, owner, bucket, position in local_entries
                      ) + '</ol></div>')
    legend = "".join(f'<span><i class="support-swatch" style="background:{color}"></i>{escape(label)}</span>'
                     for label, _, color in junctions)
    options = "".join(f'<option value="{position}">{escape(label)}</option>'
                      for position, (label, _, _) in enumerate(junctions))
    return (
        f'<section class="alignment-group" id="pileup-v{index}"><h3>Read support</h3>'
        f'<div class="support-legend">{legend}</div>'
        '<div class="read-controls" hidden><label>Show reads <select class="read-filter" aria-label="Mapping class">'
        '<option value="all" selected>All</option><option value="unique">Unique only</option></select></label>'
        '<label>Support <select class="support-filter" aria-label="Junction support">'
        f'<option value="all">All support</option>{options}</select></label>'
        '<button class="read-expand" type="button" aria-expanded="false">Expand reads</button>'
        '<button class="read-previous" type="button">Previous reads</button>'
        '<button class="read-next" type="button">Next reads</button>'
        '<span class="read-page" role="status" aria-live="polite"></span></div>'
        '<p class="read-filter-empty" hidden>No embedded records match this filter.</p>'
        + ('<p>No eligible alignment records (0).</p>' if not retained else '')
        + '<div class="pileup-scroll">' + "".join(tracks) + '</div>'
        '</section>'
    )


def _detail(index: int, result: VariantSupport, evidence: AlignmentEvidence | None) -> str:
    maximum = max(1, result.total, *(support.total for support in result.reference_junctions))
    junctions = [("Defining junction", result), *(
        (_role_label(support.definition.role), support) for support in result.reference_junctions
    )]
    rows = []
    for position, (label, support) in enumerate(junctions):
        rows.append(
            f'<tr><th scope="row">{escape(label)}</th>'
            f'<td>{_counts(support)}{_bar(support, maximum, reference=position > 0)}</td>'
            f'<td class="number">{support.total}</td>'
            f'<td class="coordinates">{_coordinates(support)}</td></tr>'
        )
    return (
        f'<article id="variant-{index}" class="variant-detail panel">'
        f'<h2 tabindex="-1">{escape(result.definition.identifier)}</h2>'
        '<div class="table-scroll"><table class="detail-table">'
        '<thead><tr><th scope="col">Junction</th><th scope="col">Support</th>'
        '<th scope="col">Total</th><th scope="col">Coordinates (1-based inclusive intron)</th>'
        f'</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'
        + _combined_pileup(index, result, evidence)
        + '<details class="junction-details"><summary>Junction annotations and alignment counts</summary>'
        + "".join(
            _schematic(name, support) + _alignment_group(
                evidence, f"v{index}-j{position}",
            ) for position, (name, support) in enumerate(junctions)
        )
        + '</details><a class="back-link" href="#overview">Back to all variants</a></article>'
    )


def write_html(
    genome_assembly: str,
    results: Iterable[VariantSupport],
    stream: TextIO,
    *,
    junctions_name: str = "",
    definitions_name: str = "",
    sample_name: str = "",
    compatibility_warning: bool = False,
    alignment_evidence: AlignmentEvidence | None = None,
) -> None:
    """Write deterministic HTML with escaped text and no external resources.

    Data is rendered as HTML text, never interpolated into executable JavaScript.
    The static overview and details remain readable without JavaScript.
    """
    rows = tuple(results)
    roles = list(dict.fromkeys(
        support.definition.role for result in rows for support in result.reference_junctions
    ))
    maximum = max(1, *(support.total for result in rows for support in (result, *result.reference_junctions)))
    ordered = sorted(enumerate(rows), key=lambda pair: (-pair[1].unique, -pair[1].total, pair[0]))
    resources = files("sjsift")
    template = Template(resources.joinpath("html_report.html").read_text(encoding="utf-8"))
    stream.write(template.substitute(
        title=escape(f"{sample_name} · SJSIFT junction report" if sample_name else "SJSIFT junction report"),
        sample_subtitle=f'<p class="sample-name">{escape(sample_name)}</p>' if sample_name else "",
        styles=resources.joinpath("html_report.css").read_text(encoding="utf-8"),
        script=resources.joinpath("html_report.js").read_text(encoding="utf-8"),
        version=escape(__version__),
        assembly=escape(genome_assembly),
        junctions=escape(junctions_name or "Not supplied"),
        count=len(rows),
        alignment_file=escape(alignment_evidence.provenance.get("Alignment file", "Not supplied")
                              if alignment_evidence else "Not supplied"),
        observed=sum(bool(result.total) for result in rows),
        warning=(
            '<aside class="warning" role="note"><strong>Check input compatibility.</strong> '
            'No catalog chromosome identifiers were found in the STAR file. '
            'The file may be empty; check genome assembly and chromosome naming.</aside>'
            if compatibility_warning else ""
        ),
        reference_headers="".join(
            f'<th scope="col">Reference · {escape(_role_label(role))}</th>' for role in roles
        ),
        overview_rows="".join(_overview_row(index, result, roles, maximum) for index, result in ordered),
        detail_sections="".join(_detail(index, result, alignment_evidence) for index, result in enumerate(rows)),
    ))
