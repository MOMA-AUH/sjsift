"""Tests for validating and quantifying STAR splice junctions."""

from pathlib import Path

import pytest

from sjsift.catalog import Catalog, VariantDefinition
from sjsift.quantify import QuantifyError, VariantSupport, quantify


CATALOG = Catalog(
    genome_assembly="GRCh38",
    variants=(
        VariantDefinition("plus", "chr7", 10, 20, "+"),
        VariantDefinition("minus", "chrX", 30, 40, "-"),
    ),
)


def quantify_text(tmp_path: Path, text: str) -> tuple[VariantSupport, ...]:
    path = tmp_path / "sample.SJ.out.tab"
    path.write_text(text, encoding="utf-8")
    return quantify(CATALOG, path)


def assert_invalid(tmp_path: Path, text: str, *messages: str) -> None:
    path = tmp_path / "sample.SJ.out.tab"
    path.write_text(text, encoding="utf-8")

    with pytest.raises(QuantifyError) as caught:
        quantify(CATALOG, path)

    assert str(path) in str(caught.value)
    for message in messages:
        assert message in str(caught.value)


def test_exact_matches_copy_support_and_calculate_unweighted_total(
    tmp_path: Path,
) -> None:
    results = quantify_text(
        tmp_path,
        "chr7\t10\t20\t1\t0\t0\t23\t4\t0\n"
        "chrX\t30\t40\t2\t6\t1\t7\t11\t99\n",
    )

    assert results == (
        VariantSupport(CATALOG.variants[0], unique=23, multimapping=4),
        VariantSupport(CATALOG.variants[1], unique=7, multimapping=11),
    )
    assert tuple(result.total for result in results) == (27, 18)


@pytest.mark.parametrize(
    "row",
    [
        "chr7\t10\t20\t1\t0\t0\t0\t0\t0",
        "chr7\t10\t20\t1\t6\t1\t1\t2\t3",
        "chrX\t30\t40\t2\t1\t0\t4\t5\t6",
    ],
)
def test_documented_star_field_values_are_accepted(
    tmp_path: Path, row: str
) -> None:
    quantify_text(tmp_path, row + "\n")


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ("chr7\t10\t20\t1\t0\t0\t0\t0", "exactly 9 tab-separated fields"),
        ("chr7\t10\t20\t1\t0\t0\t0\t0\t0\textra", "exactly 9 tab-separated fields"),
        ("\t10\t20\t1\t0\t0\t0\t0\t0", "chromosome must not be empty"),
        ("chr7\tx\t20\t1\t0\t0\t0\t0\t0", "intron start must be an integer"),
        ("chr7\t10\tx\t1\t0\t0\t0\t0\t0", "intron end must be an integer"),
        ("chr7\t10\t20\tx\t0\t0\t0\t0\t0", "strand must be an integer"),
        ("chr7\t10\t20\t1\tx\t0\t0\t0\t0", "motif must be an integer"),
        ("chr7\t10\t20\t1\t0\tx\t0\t0\t0", "annotation status must be an integer"),
        ("chr7\t10\t20\t1\t0\t0\tx\t0\t0", "unique support must be an integer"),
        ("chr7\t10\t20\t1\t0\t0\t0\tx\t0", "multimapping support must be an integer"),
        ("chr7\t10\t20\t1\t0\t0\t0\t0\tx", "overhang must be an integer"),
        ("chr7\t 10\t20\t1\t0\t0\t0\t0\t0", "intron start must be an integer"),
        ("chr7\t0\t20\t1\t0\t0\t0\t0\t0", "intron start must be positive"),
        ("chr7\t10\t0\t1\t0\t0\t0\t0\t0", "intron end must be positive"),
        ("chr7\t20\t10\t1\t0\t0\t0\t0\t0", "intron end must not be less than intron start"),
        ("chr7\t10\t20\t-1\t0\t0\t0\t0\t0", "strand must be one of 0, 1, or 2"),
        ("chr7\t10\t20\t3\t0\t0\t0\t0\t0", "strand must be one of 0, 1, or 2"),
        ("chr7\t10\t20\t1\t-1\t0\t0\t0\t0", "motif must be between 0 and 6"),
        ("chr7\t10\t20\t1\t7\t0\t0\t0\t0", "motif must be between 0 and 6"),
        ("chr7\t10\t20\t1\t0\t-1\t0\t0\t0", "annotation status must be 0 or 1"),
        ("chr7\t10\t20\t1\t0\t2\t0\t0\t0", "annotation status must be 0 or 1"),
        ("chr7\t10\t20\t1\t0\t0\t-1\t0\t0", "unique support must not be negative"),
        ("chr7\t10\t20\t1\t0\t0\t0\t-1\t0", "multimapping support must not be negative"),
        ("chr7\t10\t20\t1\t0\t0\t0\t0\t-1", "overhang must not be negative"),
    ],
)
def test_every_malformed_field_class_is_rejected_with_line_context(
    tmp_path: Path, row: str, message: str
) -> None:
    assert_invalid(
        tmp_path,
        "unsupported\t1\t2\t0\t0\t0\t0\t0\t0\n" + row + "\n",
        "line 2",
        message,
    )


@pytest.mark.parametrize(
    "near_miss",
    [
        "7\t10\t20\t1\t1\t1\t9\t8\t7",
        "chr7\t9\t20\t1\t1\t1\t9\t8\t7",
        "chr7\t10\t19\t1\t1\t1\t9\t8\t7",
        "chr7\t10\t20\t2\t1\t1\t9\t8\t7",
        "chr7\t10\t20\t0\t1\t1\t9\t8\t7",
    ],
)
def test_near_misses_and_undefined_strand_do_not_match(
    tmp_path: Path, near_miss: str
) -> None:
    results = quantify_text(tmp_path, near_miss + "\n")

    assert results[0].unique == 0
    assert results[0].multimapping == 0


def test_motif_annotation_and_overhang_do_not_participate_in_matching(
    tmp_path: Path,
) -> None:
    first = quantify_text(tmp_path, "chr7\t10\t20\t1\t0\t0\t3\t4\t0\n")
    second = quantify_text(tmp_path, "chr7\t10\t20\t1\t6\t1\t3\t4\t999\n")

    assert first == second


def test_star_row_order_does_not_change_catalog_order_or_values(
    tmp_path: Path,
) -> None:
    rows = [
        "chr7\t10\t20\t1\t1\t1\t3\t4\t5",
        "chrX\t30\t40\t2\t2\t0\t6\t7\t8",
    ]

    forward = quantify_text(tmp_path, "\n".join(rows) + "\n")
    reverse = quantify_text(tmp_path, "\n".join(reversed(rows)) + "\n")

    assert forward == reverse
    assert tuple(result.definition.identifier for result in forward) == (
        "plus",
        "minus",
    )


@pytest.mark.parametrize("text", ["", "\n\r\n", "chr2\t1\t2\t0\t0\t0\t5\t6\t7\n"])
def test_empty_or_unmatched_input_produces_ordered_zero_results(
    tmp_path: Path, text: str
) -> None:
    results = quantify_text(tmp_path, text)

    assert tuple(result.definition for result in results) == CATALOG.variants
    assert tuple((result.unique, result.multimapping, result.total) for result in results) == (
        (0, 0, 0),
        (0, 0, 0),
    )


def test_duplicate_rows_matching_one_catalog_target_are_rejected(
    tmp_path: Path,
) -> None:
    row = "chr7\t10\t20\t1\t1\t1\t3\t4\t5\n"

    assert_invalid(tmp_path, row + row, "line 2", "duplicate", "plus")


def test_malformed_unmatched_rows_are_still_rejected(tmp_path: Path) -> None:
    assert_invalid(
        tmp_path,
        "unsupported\tbad\t2\t0\t0\t0\t0\t0\t0\n",
        "line 1",
        "intron start must be an integer",
    )
