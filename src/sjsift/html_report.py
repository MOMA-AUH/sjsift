"""Render a self-contained, offline overview and per-variant evidence report."""

from collections.abc import Iterable
from dataclasses import dataclass
from html import escape
from importlib.resources import files
from math import log1p
from string import Template
from typing import TextIO

from . import __version__
from .alignment_evidence import AlignmentEvidence, AlignmentRecord
from .quantify import ReferenceJunctionSupport, VariantSupport


_ROLE_LABELS = {"same_donor": "Same donor", "same_acceptor": "Same acceptor"}


def _support_colors(roles):
    colors = {None: "#3159ad", "same_donor": "#168375", "same_acceptor": "#9950a5"}
    for index, role in enumerate(sorted(set(roles) - colors.keys())):
        colors[role] = f"hsl({(30 + index * 137.508) % 360:.3f}, 60%, 40%)"
    return colors


def _support_color(support, colors):
    return colors[support.definition.role if isinstance(support, ReferenceJunctionSupport) else None]


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
    color: str,
) -> str:
    unique_width = 100 * (support.unique / maximum)
    multi_width = 100 * (support.multimapping / maximum)
    return (
        f'<span class="bar" style="--support-color:{color}" aria-hidden="true">'
        f'<span class="unique-support" style="width:{unique_width:.6f}%"></span>'
        f'<span class="multimapping" style="width:{multi_width:.6f}%"></span></span>'
    )


def _coordinates(support: VariantSupport | ReferenceJunctionSupport) -> str:
    definition = support.definition
    return escape(
        f"{definition.chromosome}:{definition.intron_start}–{definition.intron_end} "
        f"({definition.strand})"
    )


def _overview_row(
    index: int, result: VariantSupport, roles: list[str], maximum: int, colors
) -> str:
    definition = result.definition
    references = {support.definition.role: support for support in result.reference_junctions}
    cells = []
    for role in roles:
        reference = references.get(role)
        content = (
            _counts(reference) + _bar(reference, maximum, colors[role])
            if reference is not None else '<span class="muted">Not configured</span>'
        )
        cells.append(f"<td>{content}</td>")
    return (
        f'<tr data-index="{index}" data-unique="{result.unique}" data-total="{result.total}">'
        f'<th scope="row"><a href="#variant-{index}">{escape(definition.identifier)}</a>'
        f'<small>{escape(definition.chromosome)} · {escape(definition.strand)}</small></th>'
        f'<td>{_counts(result)}{_bar(result, maximum, colors[None])}</td>'
        + "".join(cells)
        + f'<td class="observation">{escape(_observation(result)[0])}</td></tr>'
    )


def _display_interval(record: AlignmentRecord, operation) -> tuple[int, int]:
    """Extend soft clips visually without assigning them reference coordinates."""
    if operation.op == "S":
        if operation.start == record.start:
            return operation.start - operation.length, operation.start
        return operation.start, operation.start + operation.length
    return operation.start, operation.end


def _exon_windows(exon_extents, evidence_intervals):
    """Retain local exon context, without letting uncovered tails set the scale."""
    windows = []
    for start, end in exon_extents:
        local = [(max(start, left - 120), min(end, right + 120))
                 for left, right in evidence_intervals if left <= end and right >= start]
        if not local:
            continue
        first = min(left for left, _ in local)
        last = max(right for _, right in local)
        # Keep short remainders intact, just as we keep short internal gaps.
        if first - start <= 80:
            local.append((start, first))
        if end - last <= 80:
            local.append((last, end))
        windows.extend(local)
    return windows


@dataclass(frozen=True)
class _PileupScale:
    """One coordinate transform for all embedded records on a figure's chromosome."""

    start: int
    end: int
    breaks: tuple[tuple[int, int], ...]
    left_padding: int
    right_padding: int

    @classmethod
    def for_records(cls, records: tuple[AlignmentRecord, ...], boundaries: tuple[int, ...] = (),
                    exon_extents: tuple[tuple[int, int], ...] = ()):
        evidence_intervals = [_display_interval(record, op) for record in records
                              for op in record.geometry if op.op not in "ND"]
        evidence_intervals += [(boundary, boundary) for boundary in boundaries]
        windows = _exon_windows(exon_extents, evidence_intervals)
        start = min([record.start for record in records] + list(boundaries)
                    + [left for left, _ in windows])
        end = max([record.end for record in records] + list(boundaries)
                  + [right for _, right in windows])
        # Preserve the base scale across local exon context, aligned blocks and soft clips.
        intervals = sorted(evidence_intervals + windows)
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

    def ruler(self, junctions) -> str:
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
        # A straight connected backbone, including compressed gaps.
        backbone = []
        cursor = self.percent(self.start) * 10
        for left, right in self.breaks:
            start, end = self.percent(left) * 10, self.percent(right) * 10
            backbone.append(f'<path class="axis-line" d="M {cursor:.6f} 12 H {start:.6f}"/>')
            backbone.append(
                f'<path class="gap-break" d="M {start:.6f} 12 H {end:.6f}"><title>'
                f'Compressed gap: {left + 1:,}–{right:,} ({right - left:,} nt)</title></path>')
            cursor = end
        backbone.append(f'<path class="axis-line" d="M {cursor:.6f} 12 H {self.percent(self.end) * 10:.6f}"/>')
        blocks, labels = _exon_edges(junctions, self)
        return ('<div class="pileup-ruler" aria-label="Shared genomic coordinate ruler">'
                '<svg class="exon-axis" viewBox="0 0 1000 24" preserveAspectRatio="none" aria-hidden="true">'
                + "".join(backbone) + _omitted_exon_marks(junctions, self) + blocks
                + '</svg>' + labels + "".join(ticks) + '</div>')


def _omitted_exon_marks(junctions, scale):
    """Mark missing exon numbers between neighbors from the same transcript."""
    transcripts = {}
    extents = set()
    for _, support, _ in junctions:
        junction = support.definition
        annotation = junction.annotation
        exons = transcripts.setdefault((annotation.reference_transcript, junction.strand), set())
        for exon, extent in ((annotation.donor_exon, annotation.donor_extent),
                             (annotation.acceptor_exon, annotation.acceptor_extent)):
            if extent:
                exons.add((extent, exon))
                extents.add(extent)
    marks = []
    seen = set()
    for (transcript, strand), exons in transcripts.items():
        ordered = sorted(exons)
        for (left_extent, left_exon), (right_extent, right_exon) in zip(ordered, ordered[1:]):
            if not all(exon.isascii() and exon.isdecimal() for exon in (left_exon, right_exon)):
                continue
            left_number, right_number = int(left_exon), int(right_exon)
            if (right_number - left_number) * (1 if strand == "+" else -1) <= 1:
                continue
            left, right = left_extent[1], right_extent[0] - 1
            # Do not mark across another displayed exon or overlapping annotations.
            if left >= right or any(start - 1 < right and end > left for start, end in extents):
                continue
            if (left, right) in seen or left < scale.start or right > scale.end:
                continue
            seen.add((left, right))
            start, end = scale.percent(left) * 10, scale.percent(right) * 10
            center, half = (start + end) / 2, min(9, (end - start) / 4)
            first, last = min(left_number, right_number) + 1, max(left_number, right_number) - 1
            omitted = str(first) if first == last else f"{first}–{last}"
            description = escape(f"Exons omitted from illustration: {omitted} ({transcript})")
            marks.append(
                f'<path class="omitted-exons-mask" d="M {center - half:.6f} 12 H {center + half:.6f}"/>'
                f'<path class="omitted-exons" d="M {center - half:.6f} 12 '
                f'L {center - half / 2:.6f} 8 L {center:.6f} 16 '
                f'L {center + half / 2:.6f} 8 L {center + half:.6f} 12">'
                f'<title>{description}</title></path>')
    return "".join(marks)


def _exon_annotations(junctions):
    """Yield genomic sides and optional inclusive extents in genomic order."""
    for _, support, _ in junctions:
        junction = support.definition
        annotation = junction.annotation
        exons = ((annotation.donor_exon, annotation.donor_extent),
                 (annotation.acceptor_exon, annotation.acceptor_extent))
        if junction.strand == "-":
            exons = exons[::-1]
        for boundary, side, (exon, extent) in (
            (junction.intron_start - 1, "left", exons[0]),
            (junction.intron_end, "right", exons[1]),
        ):
            yield boundary, side, exon, extent, annotation.reference_transcript, annotation.annotation_source


def _exon_edges(junctions, scale):
    """Draw visible exon segments, zig-zag cuts, and unknown-boundary ticks."""
    features = {}
    for boundary, side, exon, extent, transcript, source in _exon_annotations(junctions):
        # Merge a known exon even when different junctions touch its two ends.
        key = (extent, None, None) if extent else (None, boundary, side)
        features.setdefault(key, []).append((exon, transcript, source))
    blocks, labels = [], []
    for (extent, boundary, side), annotations in features.items():
        names = " / ".join(dict.fromkeys(exon for exon, _, _ in annotations))
        sources = "; ".join(dict.fromkeys(
            f"exon {exon} ({transcript}); {source}" for exon, transcript, source in annotations))
        if extent:
            start, end = extent
            description = f"{sources}; {start:,}–{end:,} ({end - start + 1:,} nt); annotated exon extent"
            for visible_start, visible_end in _visible_exon_segments(start - 1, end, scale):
                left = scale.percent(visible_start) * 10
                right = scale.percent(visible_end) * 10
                cut_left, cut_right = visible_start > start - 1, visible_end < end
                segment_description = description
                if cut_left or cut_right:
                    segment_description += (f"; shown segment: {visible_start + 1:,}–{visible_end:,}; "
                                            "zig-zag edge marks omitted exon sequence")
                blocks.append(_exon_segment(left, right, cut_left, cut_right, segment_description))
                labels.append(_exon_label(left, right - left, names, segment_description))
        else:
            x = scale.percent(boundary) * 10
            description = (f"{sources}; exon edge at interbase coordinate {boundary:,}; "
                           "full exon extent unknown")
            blocks.append(
                f'<path class="exon-boundary" d="M {x:.6f} 4 V 20">'
                f'<title>{escape(description)}</title></path>')
            width = 36
            left = max(0, x - width) if side == "left" else min(x, 1000 - width)
            labels.append(_exon_label(left, width, names, description, boundary=True))
    return "".join(blocks), "".join(labels)


def _visible_exon_segments(start, end, scale):
    cursor, stop = max(start, scale.start), min(end, scale.end)
    for left, right in scale.breaks:
        if right <= cursor or left >= stop:
            continue
        if cursor < left:
            yield cursor, left
        cursor = min(stop, right)
    if cursor < stop:
        yield cursor, stop


def _exon_segment(left, right, cut_left, cut_right, description):
    if not (cut_left or cut_right):
        return (f'<rect class="exon-edge" x="{left:.6f}" y="4" width="{right - left:.6f}" '
                f'height="16" rx="2"><title>{escape(description)}</title></rect>')
    depth = min(4, (right - left) / 4)
    outline = f'M {left:.6f} 4 H {right:.6f}'
    if cut_right:
        outline += (f' L {right - depth:.6f} 8 L {right:.6f} 12'
                    f' L {right - depth:.6f} 16 L {right:.6f} 20')
    else:
        outline += ' V 20'
    outline += f' H {left:.6f}'
    if cut_left:
        outline += (f' L {left + depth:.6f} 16 L {left:.6f} 12'
                    f' L {left + depth:.6f} 8 L {left:.6f} 4')
    outline += ' Z'
    return (f'<path class="exon-edge exon-truncated" d="{outline}">'
            f'<title>{escape(description)}</title></path>')


def _exon_label(left, width, names, description, boundary=False):
    kind = "exon-edge-label exon-boundary-label" if boundary else "exon-edge-label"
    return (f'<span class="{kind}" style="left:{left / 10:.6f}%;width:{width / 10:.6f}%" '
            f'title="{escape(description, quote=True)}">{escape(names)}</span>')


def _read_geometry(record: AlignmentRecord, scale: _PileupScale, junctions, read_description: str) -> str:
    operations = []
    aligned = [index for index, operation in enumerate(record.geometry) if operation.op in "M=X"]
    terminal = (aligned[-1] if record.orientation == "forward" else aligned[0]) if aligned else None
    for index, operation in enumerate(record.geometry):
        coordinates = f"{operation.start + 1}–{operation.end}" if operation.end > operation.start else f"boundary {operation.start}"
        description = f"{operation.op}: {operation.length} nt; {coordinates}"
        if operation.op == "S":
            description += "; unaligned soft clip"
        description += f"; {read_description}"
        start, end = _display_interval(record, operation)
        left = scale.percent(start)
        width = scale.percent(end) - left
        selected = next((color for support, color in junctions
                         if operation.op == "N" and operation.start == support.definition.intron_start - 1
                         and operation.end == support.definition.intron_end), None)
        operations.append(
            f'<span class="cigar-op op-{operation.op.replace("=", "equal")}{" selected-junction" if selected else ""}'
            + (f' read-end-{record.orientation}' if index == terminal else '') + '" '
            f'style="left:{left:.6f}%;width:{width:.6f}%'
            + (f';background:{selected}' if selected else '') + f'" title="{escape(description, quote=True)}"></span>'
        )
    return '<span class="read-geometry" aria-hidden="true">' + "".join(operations) + '</span>'


def _read_row(record: AlignmentRecord, position: int, scale: _PileupScale,
              memberships, bucket: int, junctions, pair_id="") -> str:
    flags = ", ".join(record.flags) or "none"
    strand = "strand unverified" if record.strand_status == "unverified" else "strand evidence agrees"
    color = next(color for index, _, _, color in memberships if index == bucket)
    support_labels = "; ".join(label for _, label, _, _ in memberships)
    style = "muted" if record.mapping_class == "unknown" or record.flag & (256 | 512 | 1024 | 2048) else (
        "hatched" if record.mapping_class == "multimapping" else "solid"
    )
    description = (f"{record.read_name}; {support_labels}; {record.mapping_class}; {strand}; "
                   f"{record.mate}, {record.orientation}; Flags: {flags}; "
                   f"NH: {record.nh if record.nh is not None else 'unavailable/unusable'}")
    if record.chimeric_tag is not None:
        description += f"; Chimeric tag: {record.chimeric_tag} (tag only; STAR counting eligibility not inferred)"
    intervals = [_display_interval(record, operation) for operation in record.geometry]
    left = scale.percent(min(start for start, _ in intervals))
    right = scale.percent(max(end for _, end in intervals))
    return (
        f'<li class="read-row {record.strand_status} read-{style}" data-mapping="{record.mapping_class}"'
        f' data-bucket="{bucket}" data-support="{" ".join(str(index) for index, _, _, _ in memberships)}"'
        f' data-pair="{pair_id}"'
        f' data-left="{left:.6f}" data-right="{right:.6f}"'
        f' style="--read-color:{color}"'
        f' data-record-id="{record.record_id}"' + (' hidden' if position >= 10 else '') + '>'
        '<details class="read-inspection"><summary class="read-hit" '
        f'aria-label="Inspect {escape(description, quote=True)}" title="{escape(description, quote=True)}">'
        + _read_geometry(record, scale, junctions, description) + '</summary>'
        '<div class="read-metadata">'
        f'<div class="read-label"><strong>{escape(record.read_name)}</strong></div>'
        f'<p>{escape(record.chromosome)}:{record.start + 1}–{record.end} · '
        f'CIGAR: {escape(record.cigar)} · MAPQ: {record.mapq} · NH: {record.nh if record.nh is not None else "unavailable/unusable"}</p>'
        f'<p>Read: {escape(record.mate.removeprefix("read "))}, {record.orientation}'
        f'<br>Flags: {escape(flags)} ({record.flag})</p>'
        + (f'<p>Chimeric tag: {escape(record.chimeric_tag)} '
           '(tag only; STAR counting eligibility not inferred)</p>' if record.chimeric_tag is not None else '')
        + f'<p>Read group: {escape(record.read_group) if record.read_group else "unassigned"} · '
        f'HI: {record.hit_index if record.hit_index is not None else "unavailable/unusable"}</p>'
        f'<p>Aligned anchors: left {record.anchors[0]} nt, right {record.anchors[1]} nt</p>'
        '</div></details></li>'
    )


def _paired_records(records):
    """Associate reciprocal, unambiguous mates within this figure."""
    candidates = {}
    for record in records:
        if (not record.flag & 1 or record.flag & (8 | 2048) or record.mate == "unspecified"
                or record.read_name == "*" or not record.hit_index_usable):
            continue
        key = (record.read_name, record.read_group, record.hit_index, record.chromosome, record.start)
        candidates.setdefault(key, []).append(record)
    partners = {}
    for records_at_start in candidates.values():
        for record in records_at_start:
            key = (record.read_name, record.read_group, record.hit_index, record.mate_chromosome, record.mate_start)
            matches = [mate for mate in candidates.get(key, ())
                       if mate.mate != record.mate and mate.mate_chromosome == record.chromosome
                       and mate.mate_start == record.start and mate.chromosome == record.chromosome
                       and bool(mate.flag & 16) == bool(record.flag & 32)
                       and bool(record.flag & 16) == bool(mate.flag & 32)
                       and bool(record.flag & 256) == bool(mate.flag & 256)]
            if len(matches) == 1:
                partners[record.record_id] = matches[0].record_id
    return {record_id: min(record_id, mate_id) for record_id, mate_id in partners.items()
            if partners.get(mate_id) == record_id}


def _junction_arcs(junctions, scale, maximum):
    """Sashimi-style connections on the read scale, with details in tooltips."""
    height = 48 + 44 * len(junctions)
    baseline = height - 12
    paths = []
    # Longer connections arch over shorter ones. Separate heights also
    # distinguish coincident junctions.
    ordered = sorted(junctions, key=lambda item: -(
        scale.project(item[1].definition.intron_end)
        - scale.project(item[1].definition.intron_start - 1)))
    for lane, (label, support, color) in enumerate(ordered):
        junction = support.definition
        annotation = junction.annotation
        left = scale.percent(junction.intron_start - 1) * 10
        right = scale.percent(junction.intron_end) * 10
        apex = 16 + lane * 44
        control = baseline - (baseline - apex) * 4 / 3
        width = 1.5 + 4.5 * log1p(support.total) / log1p(maximum)
        description = (f"{label} · exon {annotation.donor_exon} → {annotation.acceptor_exon}; "
                       f"{junction.chromosome}:{junction.intron_start:,}–{junction.intron_end:,} "
                       f"({junction.strand}); {support.unique} unique + {support.multimapping} multi (STAR)")
        paths.append(
            f'<path class="junction-arc{" zero-support" if not support.total else ""}" '
            f'd="M {left:.6f} {baseline} C {left:.6f} {control:.6f}, '
            f'{right:.6f} {control:.6f}, {right:.6f} {baseline}" '
            f'stroke="{color}" stroke-width="{width:.3f}"><title>{escape(description)}</title></path>'
            f'<path class="junction-foot" d="M {left:.6f} {baseline - 4} v 8 '
            f'M {right:.6f} {baseline - 4} v 8" stroke="{color}"/>'
        )
    return (
        f'<div class="junction-arcs" style="height:{height}px">'
        f'<svg viewBox="0 0 1000 {height}" preserveAspectRatio="none" aria-hidden="true">'
        + "".join(paths) + '</svg></div>'
    )


def _combined_pileup(index: int, result: VariantSupport, evidence: AlignmentEvidence | None, colors) -> str:
    if evidence is None:
        return '<p class="alignment-unavailable">Alignment evidence not requested.</p>'
    # Keys and colors are generated here; catalog text is never used in CSS or selectors.
    junctions = [("Variant", result, colors[None]), *(
        (_role_label(support.definition.role), support, _support_color(support, colors))
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
    pairs = _paired_records([record for record, _ in retained.values()])
    chromosomes = list(dict.fromkeys(support.definition.chromosome for _, support, _ in junctions))
    for chromosome in chromosomes:
        local_junctions = [(label, support, color) for label, support, color in junctions
                           if support.definition.chromosome == chromosome]
        local_entries = [entry for entry in entries if entry[0].chromosome == chromosome]
        boundaries = tuple(coordinate for _, support, _ in local_junctions
                           for coordinate in (support.definition.intron_start - 1, support.definition.intron_end))
        exon_extents = tuple(dict.fromkeys(
            (extent[0] - 1, extent[1]) for _, _, _, extent, _, _ in _exon_annotations(local_junctions)
            if extent is not None))
        scale = _PileupScale.for_records(tuple(entry[0] for entry in local_entries), boundaries, exon_extents)
        tracks.append(f'<div class="pileup-track" role="group" aria-label="{escape(chromosome, quote=True)}">'
                      + _junction_arcs(local_junctions, scale, max(1, *(support.total for _, support, _ in junctions)))
                      + scale.ruler(local_junctions)
                      + '<ol class="read-rows">' + "".join(
                          _read_row(record, position, scale, members, bucket,
                                    tuple((support, color) for _, _, support, color in members),
                                    pairs.get(record.record_id, ""))
                          for record, members, owner, bucket, position in local_entries
                      ) + '</ol></div>')
    options = "".join(f'<option value="{position}">{escape(label)}</option>'
                      for position, (label, _, _) in enumerate(junctions))
    return (
        f'<section class="alignment-group" id="pileup-v{index}"><h3>Read support</h3>'
        '<p>Displayed alignments may include records not counted in STAR junction support. '
        'Both mates may be shown separately in record view.</p>'
        '<p class="read-legend"><span><i class="legend-solid" aria-hidden="true"></i>Unique</span> '
        '<span><i class="legend-hatched" aria-hidden="true"></i>Multimapping</span> '
        '<span><i class="legend-muted" aria-hidden="true"></i>Unknown multiplicity or flagged '
        '(secondary, supplementary, duplicate, QC failure)</span></p>'
        '<div class="read-controls" hidden><label>Show reads <select class="read-filter" aria-label="Mapping class">'
        '<option value="all" selected>All</option><option value="unique">Unique only</option></select></label>'
        '<label><input class="read-pairs" type="checkbox"> View as pairs</label>'
        '<label>Support <select class="support-filter" aria-label="Junction support">'
        f'<option value="all">All support</option>{options}</select></label>'
        '<button class="read-expand" type="button" aria-expanded="false">Expand reads</button>'
        '<button class="read-previous" type="button">Previous reads</button>'
        '<button class="read-next" type="button">Next reads</button>'
        '<span class="read-page" role="status" aria-live="polite"></span></div>'
        '<p class="read-pair-help" hidden>Each shaded band groups two mates, with brackets on both sides. '
        'Hover or focus either read to highlight the pair. Records without a visible, '
        'unambiguous mate remain separate.</p>'
        '<p class="read-filter-empty" hidden>No embedded records match this filter.</p>'
        + ('<p>No eligible alignment records (0).</p>' if not retained else '')
        + '<div class="pileup-scroll">' + "".join(tracks) + '</div>'
        '</section>'
    )


def _detail(index: int, result: VariantSupport, evidence: AlignmentEvidence | None, colors) -> str:
    maximum = max(1, result.total, *(support.total for support in result.reference_junctions))
    junctions = [("Defining junction", result), *(
        (_role_label(support.definition.role), support) for support in result.reference_junctions
    )]
    rows = []
    for label, support in junctions:
        color = _support_color(support, colors)
        rows.append(
            f'<tr><th scope="row">{escape(label)}</th>'
            f'<td>{_counts(support)}{_bar(support, maximum, color)}</td>'
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
        + _combined_pileup(index, result, evidence, colors)
        + '</article>'
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
    colors = _support_colors(roles)
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
        support_legend="".join(
            f'<span class="swatch" style="background:{colors[role]}" aria-hidden="true"></span>'
            + escape("Variant unique" if role is None else f"{_role_label(role)} unique")
            for role in (None, *roles)
        ),
        reference_headers="".join(
            f'<th scope="col">Reference · {escape(_role_label(role))}</th>' for role in roles
        ),
        overview_rows="".join(_overview_row(index, result, roles, maximum, colors) for index, result in ordered),
        detail_sections="".join(_detail(index, result, alignment_evidence, colors) for index, result in enumerate(rows)),
    ))
