"""Tests for the orthologs-to-rbh import path.

Covers:
- BED parsing (4-column, 6-column, CRLF, blank lines, bad rows).
- Orthologs-table parsing (with and without orthogroup-id column,
  auto-detection, header row, --species-order override, multi-copy
  rejection, empty cells, mismatched column counts).
- Join + .rbh emission (column order, sequential rbh ids, midpoint
  position, skipped rows when gene ids aren't in the BED, duplicate
  gene-id handling, full N-species path).
- The CLI subcommand end-to-end via subprocess.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from textwrap import dedent

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
CLI_PATH = REPO_ROOT / "bin" / "odp"


@pytest.fixture(scope="module")
def importer(source_dir):
    """Import the importer module via the source_dir fixture (which inserts
    source/ onto sys.path)."""
    import ortholog_table_importer
    return ortholog_table_importer


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def two_species_bed_pair(tmp_path):
    sp1 = tmp_path / "sp1.bed"
    sp1.write_text(dedent("""\
        chr1\t100\t200\tg1
        chr1\t300\t450\tg2
        chr2\t1000\t1500\tg3
        chr2\t2000\t2400\tg4
    """))
    sp2 = tmp_path / "sp2.bed"
    sp2.write_text(dedent("""\
        scafA\t50\t150\tx1
        scafA\t500\t700\tx2
        scafB\t10\t90\tx3
        scafB\t200\t400\tx4
    """))
    return sp1, sp2


@pytest.fixture
def two_col_orthologs(tmp_path):
    p = tmp_path / "pairs.tsv"
    p.write_text(dedent("""\
        g1\tx1
        g2\tx2
        g3\tx3
        g4\tx4
    """))
    return p


@pytest.fixture
def orthofinder_style(tmp_path):
    """Orthofinder-style: leading orthogroup-id column, then per-species
    gene columns. Header row included."""
    p = tmp_path / "Orthogroups.tsv"
    p.write_text(dedent("""\
        Orthogroup\tsp1\tsp2
        OG0000001\tg1\tx1
        OG0000002\tg2\tx2
        OG0000003\tg3\tx3
        OG0000004\tg4\tx4
    """))
    return p


@pytest.fixture
def orthofinder_no_header(tmp_path):
    """Same as Orthofinder but without the header row."""
    p = tmp_path / "Orthogroups.tsv"
    p.write_text(dedent("""\
        OG0000001\tg1\tx1
        OG0000002\tg2\tx2
        OG0000003\tg3\tx3
        OG0000004\tg4\tx4
    """))
    return p


# ---------------------------------------------------------------------------
# BED parsing
# ---------------------------------------------------------------------------


def test_parse_bed_basic(importer, two_species_bed_pair):
    sp1, _ = two_species_bed_pair
    df = importer.parse_bed(sp1)
    assert len(df) == 4
    assert list(df.columns) == ["chrom", "start", "end", "gene_id", "pos"]
    assert df.loc[0, "gene_id"] == "g1"
    assert df.loc[0, "pos"] == (100 + 200) // 2


def test_parse_bed_extra_columns_ignored(importer, tmp_path):
    p = tmp_path / "extra.bed"
    p.write_text("chr1\t10\t20\tg1\t.\t+\tfoo\n")
    df = importer.parse_bed(p)
    assert len(df) == 1
    assert df.loc[0, "gene_id"] == "g1"


def test_parse_bed_crlf_endings(importer, tmp_path):
    p = tmp_path / "crlf.bed"
    p.write_text("chr1\t10\t20\tg1\r\nchr1\t30\t40\tg2\r\n")
    df = importer.parse_bed(p)
    assert len(df) == 2
    assert df.loc[0, "gene_id"] == "g1"


def test_parse_bed_blank_lines_skipped(importer, tmp_path):
    p = tmp_path / "blanks.bed"
    p.write_text("\n\nchr1\t10\t20\tg1\n\nchr1\t30\t40\tg2\n\n")
    df = importer.parse_bed(p)
    assert len(df) == 2


def test_parse_bed_too_few_columns_errors(importer, tmp_path):
    p = tmp_path / "bad.bed"
    p.write_text("chr1\t10\t20\n")
    with pytest.raises(ValueError, match="at least 4 tab"):
        importer.parse_bed(p)


def test_parse_bed_non_integer_coords_errors(importer, tmp_path):
    p = tmp_path / "bad.bed"
    p.write_text("chr1\tNOT_AN_INT\t20\tg1\n")
    with pytest.raises(ValueError, match="must be integers"):
        importer.parse_bed(p)


def test_parse_bed_empty_file_errors(importer, tmp_path):
    p = tmp_path / "empty.bed"
    p.write_text("")
    with pytest.raises(ValueError, match="empty or has no parseable"):
        importer.parse_bed(p)


def test_detect_coord_format_chrom(importer, tmp_path):
    p = tmp_path / "x.chrom"
    p.write_text("g1\tchr1\t+\t100\t200\n")
    assert importer.detect_coord_format(p) == "chrom"


def test_detect_coord_format_bed(importer, tmp_path):
    p = tmp_path / "x.bed"
    p.write_text("chr1\t100\t200\tg1\n")
    assert importer.detect_coord_format(p) == "bed"


def test_detect_coord_format_skips_blank_then_chrom(importer, tmp_path):
    """Detection should skip blank lines (line 67-68)."""
    p = tmp_path / "x.chrom"
    p.write_text("\n\ng1\tchr1\t+\t100\t200\n")
    assert importer.detect_coord_format(p) == "chrom"


def test_detect_coord_format_unrecognized_raises(importer, tmp_path):
    """Lines that match neither format raise (line 85-89)."""
    p = tmp_path / "junk"
    p.write_text("not\ttab\tseparated\tcoordinates\there\n")
    with pytest.raises(ValueError, match="cannot detect format"):
        importer.detect_coord_format(p)


def test_detect_coord_format_empty_file_raises(importer, tmp_path):
    p = tmp_path / "empty"
    p.write_text("")
    with pytest.raises(ValueError, match="is empty"):
        importer.detect_coord_format(p)


def test_parse_chrom_basic(importer, tmp_path):
    p = tmp_path / "x.chrom"
    p.write_text("g1\tchr1\t+\t100\t200\ng2\tchr2\t-\t300\t450\n")
    df = importer.parse_chrom(p)
    assert len(df) == 2
    assert df.loc[0, "gene_id"] == "g1"
    assert df.loc[0, "strand"] == "+"
    assert df.loc[0, "pos"] == 150


def test_parse_chrom_blank_lines_skipped(importer, tmp_path):
    """Line 107-108: blank lines inside chrom file are skipped."""
    p = tmp_path / "x.chrom"
    p.write_text("\n\ng1\tchr1\t+\t100\t200\n\n")
    df = importer.parse_chrom(p)
    assert len(df) == 1


def test_parse_chrom_too_few_columns_raises(importer, tmp_path):
    """Line 110-115: fewer than 5 columns is rejected."""
    p = tmp_path / "bad.chrom"
    p.write_text("g1\tchr1\t+\t100\n")
    with pytest.raises(ValueError, match="needs at least 5 tab"):
        importer.parse_chrom(p)


def test_parse_chrom_bad_strand_raises(importer, tmp_path):
    """Line 117-121: strand must be +, -, or '.'."""
    p = tmp_path / "bad.chrom"
    p.write_text("g1\tchr1\tX\t100\t200\n")
    with pytest.raises(ValueError, match="strand must be"):
        importer.parse_chrom(p)


def test_parse_chrom_non_integer_coords_raises(importer, tmp_path):
    """Line 125-129: start/stop must be ints."""
    p = tmp_path / "bad.chrom"
    p.write_text("g1\tchr1\t+\tNOT_AN_INT\t200\n")
    with pytest.raises(ValueError, match="start/stop must be integers"):
        importer.parse_chrom(p)


def test_parse_chrom_empty_file_raises(importer, tmp_path):
    """Line 138-139: empty after blank-stripping is rejected."""
    p = tmp_path / "empty.chrom"
    p.write_text("\n\n")
    with pytest.raises(ValueError, match="is empty or has no parseable"):
        importer.parse_chrom(p)


def test_parse_coordinates_unsupported_format_raises(importer, tmp_path):
    """Line 204: explicit fmt that's not 'auto'/'chrom'/'bed' is rejected."""
    p = tmp_path / "x.chrom"
    p.write_text("g1\tchr1\t+\t100\t200\n")
    with pytest.raises(ValueError, match="unsupported coordinate format"):
        importer.parse_coordinates(p, fmt="vcf")


# ---------------------------------------------------------------------------
# Orthologs-table parsing
# ---------------------------------------------------------------------------


def test_parse_two_col_no_header_no_og(importer, two_col_orthologs):
    names, df = importer.parse_ortholog_table(two_col_orthologs)
    # Auto-detect: first column doesn't match orthogroup id pattern, so
    # treated as a gene column.
    assert names == ["sp1", "sp2"]
    assert len(df) == 4
    assert df.loc[0, "sp1"] == "g1"
    assert df.loc[0, "sp2"] == "x1"


def test_parse_orthofinder_style_auto_detects_og_column(importer, orthofinder_no_header):
    names, df = importer.parse_ortholog_table(orthofinder_no_header)
    # First column values match ^[A-Za-z]+\d+$, so auto-detected as OG ids.
    # Species names default to sp1/sp2.
    assert names == ["sp1", "sp2"]
    assert len(df) == 4
    assert df.loc[0, "sp1"] == "g1"
    assert df.loc[0, "sp2"] == "x1"


def test_parse_orthofinder_style_with_header(importer, orthofinder_style):
    names, df = importer.parse_ortholog_table(
        orthofinder_style,
        has_header=True,
    )
    assert names == ["sp1", "sp2"]
    assert len(df) == 4


def test_parse_species_order_overrides_header(importer, orthofinder_style):
    names, df = importer.parse_ortholog_table(
        orthofinder_style,
        species_names=["A", "B"],
        has_header=True,
        has_orthogroup_id_column=True,
    )
    assert names == ["A", "B"]
    assert "A" in df.columns
    assert "B" in df.columns


def test_parse_rejects_multi_copy_row(importer, tmp_path):
    p = tmp_path / "multi.tsv"
    p.write_text("g1\tx1,x1b\n")
    with pytest.raises(ValueError, match="multiple genes"):
        importer.parse_ortholog_table(p, species_names=["sp1", "sp2"])


def test_parse_rejects_empty_cell(importer, tmp_path):
    p = tmp_path / "empty_cell.tsv"
    p.write_text("g1\t\ng2\tx2\n")
    with pytest.raises(ValueError, match="empty cell"):
        importer.parse_ortholog_table(p, species_names=["sp1", "sp2"])


def test_parse_rejects_inconsistent_column_count(importer, tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("g1\tx1\ng2\tx2\tEXTRA\n")
    with pytest.raises(ValueError, match="inconsistent column count"):
        importer.parse_ortholog_table(p, species_names=["sp1", "sp2"])


def test_parse_strips_crlf(importer, tmp_path):
    p = tmp_path / "crlf.tsv"
    p.write_text("g1\tx1\r\ng2\tx2\r\n")
    names, df = importer.parse_ortholog_table(
        p, species_names=["sp1", "sp2"],
    )
    assert len(df) == 2
    # No \r should leak through.
    assert "\r" not in df.loc[0, "sp2"]


def test_parse_species_names_count_mismatch_errors(importer, two_col_orthologs):
    with pytest.raises(ValueError, match="3 entries but the table has 2"):
        importer.parse_ortholog_table(
            two_col_orthologs, species_names=["a", "b", "c"],
        )


def test_parse_empty_table_raises(importer, tmp_path):
    """Line 263-264: empty file is rejected up front."""
    p = tmp_path / "empty.tsv"
    p.write_text("")
    with pytest.raises(ValueError, match="is empty"):
        importer.parse_ortholog_table(p)


def test_parse_only_header_no_data_raises(importer, tmp_path):
    """Line 273-274: header but no data rows."""
    p = tmp_path / "header_only.tsv"
    p.write_text("Orthogroup\tsp1\tsp2\n")
    with pytest.raises(ValueError, match="has no data rows"):
        importer.parse_ortholog_table(p, has_header=True)


def test_parse_species_count_with_og_column_mismatch(importer, tmp_path):
    """Line 301-309: species count doesn't match table cols even with OG."""
    p = tmp_path / "x.tsv"
    p.write_text("OG0000001\tg1\tx1\n")
    with pytest.raises(ValueError, match="Expected either"):
        importer.parse_ortholog_table(p, species_names=["only_one"])


def test_orthologs_to_rbh_warns_on_duplicate_bed_gene_ids(
    importer, tmp_path, capsys
):
    """Line 422-430: duplicate gene_id in a coordinate file emits a warning
    and the first occurrence is kept."""
    sp1 = tmp_path / "sp1.bed"
    sp1.write_text(dedent("""\
        chr1\t100\t200\tg1
        chr1\t300\t400\tg1
        chr2\t500\t600\tg2
    """))
    sp2 = tmp_path / "sp2.bed"
    sp2.write_text(dedent("""\
        scafA\t10\t20\tx1
        scafA\t30\t40\tx2
    """))
    orth = tmp_path / "pairs.tsv"
    orth.write_text("g1\tx1\ng2\tx2\n")
    df = importer.orthologs_to_rbh(
        orthologs_path=orth,
        bed_paths={"sp1": sp1, "sp2": sp2},
        species_order=["sp1", "sp2"],
    )
    assert len(df) == 2
    # First-occurrence wins for the duplicate gene id.
    g1_row = df[df["sp1_gene"] == "g1"].iloc[0]
    assert g1_row["sp1_pos"] == 150
    err = capsys.readouterr().err
    assert "duplicate gene_id entries" in err


# ---------------------------------------------------------------------------
# Join + .rbh emission
# ---------------------------------------------------------------------------


def test_orthologs_to_rbh_two_species(
    importer, two_species_bed_pair, two_col_orthologs
):
    sp1_bed, sp2_bed = two_species_bed_pair
    df = importer.orthologs_to_rbh(
        orthologs_path=two_col_orthologs,
        bed_paths={"sp1": sp1_bed, "sp2": sp2_bed},
        species_order=["sp1", "sp2"],
    )
    assert len(df) == 4
    expected_cols = [
        "rbh", "gene_group",
        "sp1_gene", "sp1_scaf", "sp1_pos",
        "sp2_gene", "sp2_scaf", "sp2_pos",
    ]
    assert list(df.columns) == expected_cols
    # First row: g1 ↔ x1
    row = df.iloc[0]
    assert row["rbh"] == "rbh1"
    assert row["gene_group"] == "None"
    assert row["sp1_gene"] == "g1"
    assert row["sp1_scaf"] == "chr1"
    assert row["sp1_pos"] == (100 + 200) // 2
    assert row["sp2_gene"] == "x1"
    assert row["sp2_scaf"] == "scafA"
    assert row["sp2_pos"] == (50 + 150) // 2


def test_orthologs_to_rbh_skips_missing_bed_genes(
    importer, two_species_bed_pair, tmp_path, capsys
):
    sp1_bed, sp2_bed = two_species_bed_pair
    p = tmp_path / "pairs.tsv"
    # g1 is in sp1.bed, MISSING is not — that row must be skipped.
    p.write_text("g1\tx1\nMISSING\tx2\n")
    df = importer.orthologs_to_rbh(
        orthologs_path=p,
        bed_paths={"sp1": sp1_bed, "sp2": sp2_bed},
        species_order=["sp1", "sp2"],
    )
    assert len(df) == 1
    assert df.iloc[0]["sp1_gene"] == "g1"
    err = capsys.readouterr().err
    assert "skipped 1 ortholog rows" in err


def test_orthologs_to_rbh_missing_bed_for_species_errors(
    importer, two_species_bed_pair, two_col_orthologs
):
    sp1_bed, _ = two_species_bed_pair
    with pytest.raises(ValueError, match="Coordinate files missing for species"):
        importer.orthologs_to_rbh(
            orthologs_path=two_col_orthologs,
            bed_paths={"sp1": sp1_bed},  # sp2 missing
            species_order=["sp1", "sp2"],
        )


def test_orthologs_to_rbh_three_species_orthofinder_style(
    importer, two_species_bed_pair, tmp_path
):
    sp1_bed, sp2_bed = two_species_bed_pair
    # Third species BED
    sp3_bed = tmp_path / "sp3.bed"
    sp3_bed.write_text(dedent("""\
        ctgA\t1\t10\ty1
        ctgA\t20\t30\ty2
        ctgB\t100\t200\ty3
        ctgB\t300\t400\ty4
    """))
    # Orthofinder-style table: OG id + 3 species
    orth = tmp_path / "Single_copy.tsv"
    orth.write_text(dedent("""\
        OG0000001\tg1\tx1\ty1
        OG0000002\tg2\tx2\ty2
        OG0000003\tg3\tx3\ty3
        OG0000004\tg4\tx4\ty4
    """))
    df = importer.orthologs_to_rbh(
        orthologs_path=orth,
        bed_paths={"sp1": sp1_bed, "sp2": sp2_bed, "sp3": sp3_bed},
        species_order=["sp1", "sp2", "sp3"],
        has_orthogroup_id_column=True,
    )
    assert len(df) == 4
    assert "sp3_gene" in df.columns
    assert df.iloc[0]["sp3_gene"] == "y1"
    assert df.iloc[0]["sp3_scaf"] == "ctgA"


def test_orthologs_to_rbh_rejects_zero_rows_after_join(
    importer, two_species_bed_pair, tmp_path
):
    sp1_bed, sp2_bed = two_species_bed_pair
    p = tmp_path / "miss.tsv"
    p.write_text("NONE\tNONE\n")
    with pytest.raises(ValueError, match="zero ortholog rows survived"):
        importer.orthologs_to_rbh(
            orthologs_path=p,
            bed_paths={"sp1": sp1_bed, "sp2": sp2_bed},
            species_order=["sp1", "sp2"],
        )


def test_write_rbh_round_trip(importer, two_species_bed_pair, two_col_orthologs, tmp_path):
    sp1_bed, sp2_bed = two_species_bed_pair
    df = importer.orthologs_to_rbh(
        orthologs_path=two_col_orthologs,
        bed_paths={"sp1": sp1_bed, "sp2": sp2_bed},
        species_order=["sp1", "sp2"],
    )
    out = tmp_path / "rbh.tsv"
    n = importer.write_rbh(df, out)
    assert n == 4
    body = out.read_text()
    assert body.startswith("rbh\tgene_group\tsp1_gene")
    # Reading it back as a plain TSV recovers the same column count.
    lines = [l for l in body.split("\n") if l.strip()]
    assert len(lines) == 5  # 1 header + 4 rows


# ---------------------------------------------------------------------------
# CLI end-to-end
# ---------------------------------------------------------------------------


def test_cli_orthologs_to_rbh_two_species(
    tmp_path, two_species_bed_pair, two_col_orthologs
):
    sp1_bed, sp2_bed = two_species_bed_pair
    out = tmp_path / "rbh.tsv"
    result = subprocess.run(
        [
            sys.executable, str(CLI_PATH),
            "orthologs-to-rbh",
            "--orthologs", str(two_col_orthologs),
            "--bed", f"sp1={sp1_bed}",
            "--bed", f"sp2={sp2_bed}",
            "--output", str(out),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out.is_file()
    body = out.read_text()
    assert "rbh1" in body
    assert "g1" in body


def test_cli_bed_missing_equals_sign(tmp_path):
    result = subprocess.run(
        [
            sys.executable, str(CLI_PATH),
            "orthologs-to-rbh",
            "--orthologs", "/dev/null",
            "--bed", "missing_equals",
            "--output", str(tmp_path / "rbh.tsv"),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert "expected NAME=PATH" in result.stderr


def test_cli_requires_at_least_one_coord(tmp_path, two_col_orthologs):
    result = subprocess.run(
        [
            sys.executable, str(CLI_PATH),
            "orthologs-to-rbh",
            "--orthologs", str(two_col_orthologs),
            "--output", str(tmp_path / "rbh.tsv"),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert "at least one --chrom or --bed" in result.stderr


def test_cli_help_lists_orthologs_to_rbh():
    result = subprocess.run(
        [sys.executable, str(CLI_PATH), "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "orthologs-to-rbh" in result.stdout
