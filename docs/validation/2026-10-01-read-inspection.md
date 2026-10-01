# Structural inspection validation — 2026-10-01

Issue #32 was checked against real synthetic indexed BAM records. The full
Python suite passed (313 tests), JavaScript syntax passed `node --check`, and
`git diff --check` passed. Independent standards and specification reviews
reported no findings.

Firefox opened the generated report directly from disk. The fixture contained
139 eligible defining-junction records (131 NH=1, five NH>1, three unknown), a
105-per-class cap, a shared reference, a long-gap alignment and an empty variant.
Observed browser behavior:

- The default All preview showed 10 of 113 embedded records. Expansion showed
  50, then 50, then 13 records; the last page said End of embedded reads and
  disabled Next.
- Unique only reset to page one with 105 matching embedded records. Collapse
  restored ten previews. Full eligible/embedded totals remained 139/113.
- Variant Next navigation selected the long-gap and empty variants. The empty
  view separately stated zero eligible records and no reference context.
- The shared reference showed Also matches membership. The long-gap row showed
  distinct H/S/M/I/=/X/N/D operations, a visible break and exact 1100 nt skip.
  Native Read details exposed CIGAR, flags, MAPQ, orientation, identity, strand
  evidence and anchors of 10 and 3 nt for that selected splice.
- A hostile read name containing an image/event-handler-shaped string remained
  literal text. Automated serialization checks also reject embedded read bases,
  qualities and absolute input paths.
- With Firefox Network Monitor recording during reload, the only entry was the
  local `file:` HTML document: one request, 0 B transferred. No HTTP or external
  asset request appeared.
- A copy with document-level CSP `script-src 'none'` prevented all report script
  execution. Its complete static overview, all three variant sections, ten-row
  previews, full counts and provenance remained readable; controls stayed hidden.
  This check changed only the disposable report copy, not browser preferences.

These are software fixtures, not biological or clinical validation. No release
or publication was performed.
