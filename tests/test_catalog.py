"""Tests for loading and validating variant-definition catalogs."""

from pathlib import Path

import pytest

from sjsift.catalog import CatalogError, ReferenceJunction, VariantDefinition, load_catalog


VALID_CATALOG = """\
schema_version = 1
genome_assembly = "GRCh38"

[[variants]]
id = "first"
chromosome = "chr7"
intron_start = 10
intron_end = 20
strand = "+"

[[variants]]
id = "second"
chromosome = "chrX"
intron_start = 30
intron_end = 30
strand = "-"
"""


def load_text(tmp_path: Path, text: str):
    path = tmp_path / "catalog.toml"
    path.write_text(text, encoding="utf-8")
    return path, load_catalog(path)


def assert_invalid(tmp_path: Path, text: str, *messages: str) -> None:
    path = tmp_path / "catalog.toml"
    path.write_text(text, encoding="utf-8")

    with pytest.raises(CatalogError) as caught:
        load_catalog(path)

    assert str(path) in str(caught.value)
    for message in messages:
        assert message in str(caught.value)


def test_valid_catalog_preserves_entry_order(tmp_path: Path) -> None:
    _, catalog = load_text(tmp_path, VALID_CATALOG)

    assert catalog.genome_assembly == "GRCh38"
    assert catalog.variants == (
        VariantDefinition("first", "chr7", 10, 20, "+"),
        VariantDefinition("second", "chrX", 30, 30, "-"),
    )


def test_version_two_catalog_loads_ordered_named_reference_junctions(
    tmp_path: Path,
) -> None:
    text = VALID_CATALOG.replace("schema_version = 1", "schema_version = 2")
    text = text.replace(
        'strand = "+"',
        'strand = "+"\nreference_junctions = []',
        1,
    ).replace(
        'strand = "-"',
        '''strand = "-"

[[variants.reference_junctions]]
role = "same_donor"
chromosome = "chrX"
intron_start = 41
intron_end = 50
strand = "-"

[[variants.reference_junctions]]
role = "same_acceptor"
chromosome = "chrX"
intron_start = 21
intron_end = 29
strand = "-"''',
        1,
    )

    _, catalog = load_text(tmp_path, text)

    assert catalog.variants[0].reference_junctions == ()
    assert catalog.variants[1].reference_junctions == (
        ReferenceJunction("same_donor", "chrX", 41, 50, "-"),
        ReferenceJunction("same_acceptor", "chrX", 21, 29, "-"),
    )


def test_version_two_requires_reference_junctions_for_every_variant(
    tmp_path: Path,
) -> None:
    text = VALID_CATALOG.replace("schema_version = 1", "schema_version = 2")

    assert_invalid(tmp_path, text, "variant 1 ('first')", "'reference_junctions'")


def test_reference_junction_validation_has_variant_and_junction_context(
    tmp_path: Path,
) -> None:
    text = VALID_CATALOG.replace("schema_version = 1", "schema_version = 2")
    text = text.replace(
        'strand = "+"',
        '''strand = "+"

[[variants.reference_junctions]]
role = "same_donor"
chromosome = "chr7"
intron_start = 10
intron_end = 20
strand = "+"''',
        1,
    ).replace(
        'strand = "-"',
        'strand = "-"\nreference_junctions = []',
        1,
    )

    assert_invalid(
        tmp_path,
        text,
        "variant 1 ('first'), reference junction 1 ('same_donor')",
        "must not duplicate the defining junction",
    )


@pytest.mark.parametrize("field", ["schema_version", "genome_assembly", "variants"])
def test_missing_top_level_field_is_rejected(tmp_path: Path, field: str) -> None:
    lines = VALID_CATALOG.splitlines()
    if field == "variants":
        text = "\n".join(lines[:2])
    else:
        text = "\n".join(line for line in lines if not line.startswith(f"{field} ="))

    assert_invalid(tmp_path, text, "catalog", "missing field", repr(field))


def test_unknown_top_level_field_is_rejected(tmp_path: Path) -> None:
    assert_invalid(tmp_path, "extra = true\n" + VALID_CATALOG, "unknown field", "'extra'")


def test_mistyped_top_level_field_reports_missing_and_unknown_names(
    tmp_path: Path,
) -> None:
    text = VALID_CATALOG.replace("genome_assembly", "genome_assemblly", 1)
    assert_invalid(
        tmp_path,
        text,
        "missing field(s) 'genome_assembly'",
        "unknown field(s) 'genome_assemblly'",
    )


@pytest.mark.parametrize("value", ["true", '"1"', "1.0"])
def test_schema_version_must_be_an_integer(tmp_path: Path, value: str) -> None:
    text = VALID_CATALOG.replace("schema_version = 1", f"schema_version = {value}")
    assert_invalid(tmp_path, text, "schema_version", "must be an integer")


def test_unsupported_schema_version_is_explicit(tmp_path: Path) -> None:
    text = VALID_CATALOG.replace("schema_version = 1", "schema_version = 3")
    assert_invalid(
        tmp_path, text, "unsupported catalog schema version 3", "expected 1 or 2"
    )


@pytest.mark.parametrize(
    "value", ['""', '"GR\\tCh38"', '"GR\\nCh38"', '"GR\\rCh38"', "1", "true"]
)
def test_genome_assembly_must_be_a_safe_nonempty_string(
    tmp_path: Path, value: str
) -> None:
    text = VALID_CATALOG.replace('genome_assembly = "GRCh38"', f"genome_assembly = {value}")
    assert_invalid(tmp_path, text, "genome_assembly")


@pytest.mark.parametrize("value", ["[]", '"not an array"'])
def test_variants_must_be_a_nonempty_array(tmp_path: Path, value: str) -> None:
    text = f'schema_version = 1\ngenome_assembly = "GRCh38"\nvariants = {value}\n'
    assert_invalid(tmp_path, text, "variants", "must")


def test_every_variant_must_be_a_table(tmp_path: Path) -> None:
    text = 'schema_version = 1\ngenome_assembly = "GRCh38"\nvariants = [1]\n'
    assert_invalid(tmp_path, text, "variant 1", "must be a table")


@pytest.mark.parametrize(
    "field", ["id", "chromosome", "intron_start", "intron_end", "strand"]
)
def test_missing_variant_field_is_rejected_with_entry_context(
    tmp_path: Path, field: str
) -> None:
    first_variant = VALID_CATALOG.split("[[variants]]", 2)[1]
    text = "\n".join(
        line for line in first_variant.splitlines() if not line.startswith(f"{field} =")
    )
    text = f'schema_version = 1\ngenome_assembly = "GRCh38"\n[[variants]]\n{text}\n'
    assert_invalid(tmp_path, text, "variant 1", "missing field", repr(field))


def test_unknown_variant_field_is_rejected_with_identifier(tmp_path: Path) -> None:
    text = VALID_CATALOG.replace('strand = "+"', 'strand = "+"\nstart = 10', 1)
    assert_invalid(tmp_path, text, "variant 1 ('first')", "unknown field", "'start'")


def test_mistyped_variant_field_reports_missing_and_unknown_names(tmp_path: Path) -> None:
    text = VALID_CATALOG.replace("chromosome", "chromsome", 1)
    assert_invalid(
        tmp_path,
        text,
        "variant 1 ('first')",
        "missing field(s) 'chromosome'",
        "unknown field(s) 'chromsome'",
    )


@pytest.mark.parametrize("field", ["id", "chromosome"])
@pytest.mark.parametrize(
    "value", ['""', '"has\\ttab"', '"has\\nline"', '"has\\rline"', "7", "true"]
)
def test_variant_text_fields_must_be_safe_nonempty_strings(
    tmp_path: Path, field: str, value: str
) -> None:
    old_value = '"first"' if field == "id" else '"chr7"'
    text = VALID_CATALOG.replace(f"{field} = {old_value}", f"{field} = {value}", 1)
    assert_invalid(tmp_path, text, "variant 1", repr(field))


@pytest.mark.parametrize("field", ["intron_start", "intron_end"])
@pytest.mark.parametrize("value", ["0", "-1", "true", '"10"', "10.0"])
def test_coordinates_must_be_positive_integers(
    tmp_path: Path, field: str, value: str
) -> None:
    old_value = "10" if field == "intron_start" else "20"
    text = VALID_CATALOG.replace(f"{field} = {old_value}", f"{field} = {value}", 1)
    assert_invalid(tmp_path, text, "variant 1 ('first')", repr(field))


def test_intron_end_must_not_precede_start(tmp_path: Path) -> None:
    text = VALID_CATALOG.replace("intron_end = 20", "intron_end = 9", 1)
    assert_invalid(tmp_path, text, "variant 1 ('first')", "intron_end", "intron_start")


@pytest.mark.parametrize("value", ['""', '"."', '"1"', "1", "true", "[]"])
def test_strand_must_be_plus_or_minus(tmp_path: Path, value: str) -> None:
    text = VALID_CATALOG.replace('strand = "+"', f"strand = {value}", 1)
    assert_invalid(tmp_path, text, "variant 1 ('first')", "strand", "'+' or '-'")


def test_duplicate_identifiers_are_rejected(tmp_path: Path) -> None:
    text = VALID_CATALOG.replace('id = "second"', 'id = "first"')
    assert_invalid(tmp_path, text, "variant 2 ('first')", "duplicate id", "variant 1")


def test_duplicate_defining_junctions_are_rejected(tmp_path: Path) -> None:
    text = (
        VALID_CATALOG.replace('chromosome = "chrX"', 'chromosome = "chr7"')
        .replace("intron_start = 30", "intron_start = 10")
        .replace("intron_end = 30", "intron_end = 20")
        .replace('strand = "-"', 'strand = "+"')
    )
    assert_invalid(
        tmp_path,
        text,
        "variant 2 ('second')",
        "duplicate defining junction",
        "variant 1 ('first')",
    )


def test_invalid_toml_identifies_source_path(tmp_path: Path) -> None:
    assert_invalid(tmp_path, 'schema_version = "unterminated\n', "invalid TOML")


def test_unreadable_catalog_identifies_source_path(tmp_path: Path) -> None:
    path = tmp_path / "missing.toml"

    with pytest.raises(CatalogError) as caught:
        load_catalog(path)

    assert str(path) in str(caught.value)
    assert "cannot read catalog" in str(caught.value)
