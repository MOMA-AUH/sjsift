"""Render a self-contained, offline overview and per-variant evidence report."""

from collections.abc import Iterable
from html import escape
from importlib.resources import files
from string import Template
from typing import TextIO

from . import __version__
from .alignment_evidence import AlignmentEvidence, AlignmentRecord, JunctionEvidence, MAPPING_CLASSES
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


def _read_geometry(record: AlignmentRecord) -> str:
    operations = []
    for operation in record.geometry:
        coordinates = f"{operation.start + 1}–{operation.end}" if operation.end > operation.start else f"boundary {operation.start}"
        description = f"{operation.op}: {operation.length} nt; {coordinates}"
        operations.append(
            f'<span class="cigar-op op-{operation.op.replace("=", "equal")}" '
            f'style="flex-basis:{max(12, min(operation.length, 80))}px" title="{description}">'
            f'{operation.op}<small>{operation.length}</small></span>'
        )
    return '<div class="read-geometry" aria-label="CIGAR geometry">' + "".join(operations) + '</div>'


def _read_row(record: AlignmentRecord, position: int) -> str:
    flags = ", ".join(record.flags) or "primary; no additional flags"
    strand = "strand unverified" if record.strand_status == "unverified" else "strand evidence agrees"
    return (
        f'<li class="read-row {record.strand_status}" data-mapping="{record.mapping_class}"'
        f' data-record-id="{record.record_id}"' + (' hidden' if position >= 10 else '') + '>'
        f'<div class="read-label"><strong>{escape(record.read_name)}</strong> · {record.mate} · '
        f'{record.mapping_class} · {strand}</div>'
        + _read_geometry(record)
        + '<details class="read-metadata"><summary>Read details</summary>'
        f'<p>Record ID: {record.record_id} · Read group: {escape(record.read_group or "not supplied")} · '
        f'{escape(record.chromosome)}:{record.start + 1}–{record.end}</p>'
        f'<p>CIGAR: {escape(record.cigar)} · MAPQ: {record.mapq} · NH: {record.nh if record.nh is not None else "unavailable/unusable"} · Alignment orientation: {record.orientation}</p>'
        f'<p>Flags: {escape(flags)} ({record.flag}) · Transcript {strand}</p>'
        f'<p>Strand sources: {escape("; ".join(record.strand_evidence) or "No usable transcript-strand evidence")}</p>'
        '</details></li>'
    )


def _alignment_group(group: JunctionEvidence | None, key: str) -> str:
    if group is None:
        return '<p class="alignment-unavailable">Alignment evidence not requested for this junction.</p>'
    eligible = sum(group.eligible.values())
    counts = "".join(
        f'<tr><th scope="row">{kind}</th><td>{group.eligible[kind]}</td>'
        f'<td>{sum(r.mapping_class == kind for r in group.records)}</td></tr>'
        for kind in MAPPING_CLASSES
    )
    return (
        f'<section class="alignment-group" id="reads-{key}">'
        '<h4>Alignment records</h4>'
        f'<p>Eligible: {eligible} · Embedded: {len(group.records)} · '
        f'<span class="shown-count">Currently shown: {min(10, len(group.records))}</span></p>'
        '<table class="alignment-counts"><caption>Alignment counts by mapping class</caption>'
        '<thead><tr><th>Class</th><th>Eligible</th><th>Embedded</th></tr></thead>'
        f'<tbody>{counts}</tbody></table>'
        f'<p>Transcript strand totals: agreeing {group.strand_totals["agreeing"]}; '
        f'unverified {group.strand_totals["unverified"]}. '
        f'Excluded: opposite {group.exclusions["opposite"]}; conflicting {group.exclusions["conflicting"]}.</p>'
        + ('<p>No eligible alignment records (0).</p>' if not eligible else '')
        + '<ol class="read-rows">' + "".join(_read_row(row, position) for position, row in enumerate(group.records))
        + '</ol></section>'
    )


def _alignment_provenance(evidence: AlignmentEvidence | None) -> str:
    if evidence is None:
        return '<aside class="alignment-provenance"><p>Alignment evidence not requested.</p></aside>'
    fields = "".join(f'<div><dt>{escape(key)}</dt><dd>{escape(value)}</dd></div>' for key, value in evidence.provenance.items())
    return (
        '<aside class="alignment-provenance"><h2>Alignment evidence provenance</h2>'
        f'<dl class="provenance">{fields}</dl>'
        f'<p>Sampling limit: {evidence.limit} records per junction per mapping class. '
        'All eligible records are counted; embedded samples are deterministic. Initial previews show up to 10 records per junction.</p>'
        '<p>STAR support is separate from alignment-record counts. Differences can reflect upstream STAR filtering, '
        'paired-template counting, alignment output choices, and this display policy. '
        'Records are not molecule or transcript counts.</p>'
        '<p>This shareable file retains original read names, read-group/sample identifiers, CIGARs, coordinates, '
        'mapping quality, flags and strand evidence. It contains no read sequence, base qualities, '
        'source alignment file or reference sequence.</p></aside>'
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
    label, explanation = _observation(result)
    return (
        f'<article id="variant-{index}" class="variant-detail panel">'
        f'<h2 tabindex="-1">{escape(result.definition.identifier)}</h2>'
        f'<p class="lead-count"><strong>{result.unique}</strong> unique support</p>'
        f'<p>{result.multimapping} multimapping · {result.total} total defining-junction support</p>'
        '<div class="table-scroll"><table class="detail-table">'
        '<caption>Defining junction and each configured reference, in catalog order.</caption>'
        '<thead><tr><th scope="col">Junction</th><th scope="col">Support</th>'
        '<th scope="col">Total</th><th scope="col">Coordinates (1-based inclusive intron)</th>'
        f'</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'
        '<p class="scale-note">Bars use one linear scale within this variant. '
        'Donor and acceptor roles follow transcript orientation on either strand.</p>'
        f'<div class="evidence"><strong>{escape(label)}</strong><p>{escape(explanation)}</p></div>'
        + "".join(
            _schematic(name, support) + _alignment_group(
                evidence.groups.get(f"v{index}-j{position}") if evidence else None, f"v{index}-j{position}",
            ) for position, (name, support) in enumerate(junctions)
        )
        + '<a class="back-link" href="#overview">Back to all variants</a></article>'
    )


def write_html(
    genome_assembly: str,
    results: Iterable[VariantSupport],
    stream: TextIO,
    *,
    junctions_name: str = "",
    definitions_name: str = "",
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
        title=escape(f"sjsift — {junctions_name}" if junctions_name else "sjsift junction report"),
        styles=resources.joinpath("html_report.css").read_text(encoding="utf-8"),
        script=resources.joinpath("html_report.js").read_text(encoding="utf-8"),
        version=escape(__version__),
        assembly=escape(genome_assembly),
        junctions=escape(junctions_name or "Not supplied"),
        definitions=escape(definitions_name or "Not supplied"),
        count=len(rows),
        alignment_provenance=_alignment_provenance(alignment_evidence),
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
