"""Tests that exercise heavy functions using real pipeline outputs.

The fixtures used here live in `tests/test_odp_basic/odp/` and
`tests/testdb/`, products of prior smoke-test runs. The .D.FET.rbh
files carry the full column set (FET, D, ix, breakchrom, break_ix) that
the plotting + statistics functions expect; the .blastp files are real
diamond/blastp pairwise outputs; the mini_LG_db carries a complete
ALG-vs-species set.

These tests target the heavier functions that synthetic fixtures
couldn't cover cleanly:
- ``odp_functions.reciprocal_best_hits_blastp_or_diamond_blastp``
- ``odp_functions.calc_D_for_y_and_x`` (break-mode)
- ``synteny_plot_sheet.synteny_plot_sheet``
- ``odp_ribbon_plot.ribbon_plot``
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
# Frozen pipeline outputs lifted into a dedicated unit-test fixture
# directory. They mustn't live under tests/test_odp_basic/ because
# their presence would short-circuit the snakemake smoke test.
PIPELINE_OUT = REPO_ROOT / "tests" / "testdb" / "test_fixtures"
TESTDB = REPO_ROOT / "tests" / "testdb"


def _skip_if_fixture_missing(*paths: Path) -> None:
    for p in paths:
        if not p.exists():
            pytest.skip(f"pipeline fixture missing: {p}")


@pytest.fixture(scope="module")
def of(source_dir, fasta_parser_dir):
    import odp_functions
    return odp_functions


@pytest.fixture(scope="module")
def sps(source_dir, fasta_parser_dir):
    import synteny_plot_sheet
    return synteny_plot_sheet


@pytest.fixture(scope="module")
def rp(source_dir, fasta_parser_dir):
    import odp_ribbon_plot
    return odp_ribbon_plot


@pytest.fixture(scope="module")
def rt(source_dir):
    import rbh_tools
    return rbh_tools


# ---------------------------------------------------------------------------
# Fixture paths
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def hydra_chrom():
    p = TESTDB / "mini_hydra" / "mini_hydra.chrom.gz"
    if not p.is_file():
        pytest.skip(f"missing fixture: {p}")
    return p


@pytest.fixture(scope="module")
def urchin_chrom():
    p = TESTDB / "mini_urchin" / "mini_urchin.chrom.gz"
    if not p.is_file():
        pytest.skip(f"missing fixture: {p}")
    return p


@pytest.fixture(scope="module")
def hydra_blastp():
    p = PIPELINE_OUT / "Hydra_against_Urchin.blastp"
    if not p.is_file():
        pytest.skip(f"missing fixture: {p}")
    return p


@pytest.fixture(scope="module")
def urchin_blastp():
    p = PIPELINE_OUT / "Urchin_against_Hydra.blastp"
    if not p.is_file():
        pytest.skip(f"missing fixture: {p}")
    return p


@pytest.fixture(scope="module")
def rbh_with_D_FET():
    p = PIPELINE_OUT / "Hydra_Urchin_reciprocal_best_hits.D.FET.rbh"
    if not p.is_file():
        pytest.skip(f"missing fixture: {p}")
    return p


@pytest.fixture(scope="module")
def analysis_chromsize():
    p = PIPELINE_OUT / "Hydra_Urchin.chromsize"
    if not p.is_file():
        pytest.skip(f"missing fixture: {p}")
    return p


@pytest.fixture(scope="module")
def rbh_with_color():
    """Plotted RBH file that includes a `color` column (required for
    odp_ribbon_plot.ribbon_plot)."""
    p = (PIPELINE_OUT /
         "Hydra_Urchin_xy_reciprocal_best_hits.coloredby_mBCnSSimakov2022.plotted.rbh")
    if not p.is_file():
        pytest.skip(f"missing fixture: {p}")
    return p


# ---------------------------------------------------------------------------
# reciprocal_best_hits_blastp_or_diamond_blastp — real blastp pair
# ---------------------------------------------------------------------------


def test_reciprocal_best_permissive_blastp_real_pair(
    of, hydra_blastp, urchin_blastp, tmp_path
):
    """Permissive variant — outputs every row whose evalue ties the
    per-qseqid minimum, not just the single best. Smoke test."""
    out_path = tmp_path / "permissive.tsv"
    of.reciprocal_best_permissive_blastp_or_diamond_blastp(
        str(hydra_blastp), str(urchin_blastp), str(out_path),
    )
    assert out_path.is_file()
    body = out_path.read_text()
    assert len(body.splitlines()) >= 1


def test_reciprocal_best_hits_real_blastp_pair(
    of, hydra_blastp, urchin_blastp, tmp_path
):
    """Run RBH inference on the real Hydra-vs-Urchin blastp output and
    its reverse. Verify the function produces a non-empty TSV of
    reciprocal best pairs."""
    out_path = tmp_path / "recip.tsv"
    of.reciprocal_best_hits_blastp_or_diamond_blastp(
        str(hydra_blastp), str(urchin_blastp), str(out_path),
    )
    assert out_path.is_file()
    body = out_path.read_text()
    lines = [l for l in body.splitlines() if l.strip()]
    # At least one pair survived
    assert len(lines) >= 1
    # Each line has the 12-column blast outfmt 6 shape
    for line in lines[:5]:
        parts = line.split("\t")
        assert len(parts) == 12
        # bitscore on the end is a float
        float(parts[11])


# ---------------------------------------------------------------------------
# rbh_tools.parse_rbh on a real .D.FET.rbh
# ---------------------------------------------------------------------------


def test_parse_rbh_real_DFET_file(rt, rbh_with_D_FET):
    """Real pipeline-produced .D.FET.rbh has every optional column
    (Int64 _pos, float _FET, float _D, object _breakchrom). Confirm
    parse_rbh applies the right dtypes across the board."""
    df = rt.parse_rbh(rbh_with_D_FET)
    assert len(df) > 0
    assert str(df["Hydra_pos"].dtype) == "Int64"
    assert str(df["Urchin_pos"].dtype) == "Int64"
    assert df["whole_FET"].dtype == float
    assert df["Hydra_D"].dtype == float
    assert df["Urchin_D"].dtype == float
    assert df["Hydra_breakchrom"].dtype == object


def test_real_rbh_contains_expected_columns(rt, rbh_with_D_FET):
    df = rt.parse_rbh(rbh_with_D_FET)
    needed = {
        "rbh", "gene_group", "Hydra_gene", "Hydra_scaf", "Hydra_pos",
        "Urchin_gene", "Urchin_scaf", "Urchin_pos",
        "Hydra_breakchrom", "Urchin_breakchrom",
        "Hydra_ix", "Hydra_break_ix",
        "Urchin_ix", "Urchin_break_ix",
        "whole_FET", "break_FET",
        "Hydra_D", "Urchin_D",
    }
    assert needed.issubset(set(df.columns))


# ---------------------------------------------------------------------------
# calc_D_for_y_and_x — break-mode on a real .D.FET.rbh
# ---------------------------------------------------------------------------


def test_calc_D_break_mode_on_real_rbh(of, rt, rbh_with_D_FET):
    """Break-mode of calc_D_for_y_and_x takes a parsed rbh and recomputes
    the D columns. Verify it produces non-empty Hydra_D / Urchin_D
    columns matching the row count of the input."""
    df = rt.parse_rbh(rbh_with_D_FET)
    # Drop the existing D columns so we can verify the function is doing
    # the work rather than passing through.
    df = df.drop(columns=["Hydra_D", "Urchin_D"])
    out = of.calc_D_for_y_and_x(df, xsample="Hydra", ysample="Urchin")
    assert "Hydra_D" in out.columns
    assert "Urchin_D" in out.columns
    assert len(out) == len(df)
    # D values are non-negative (norm of a vector)
    assert (out["Hydra_D"].dropna() >= 0).all()
    assert (out["Urchin_D"].dropna() >= 0).all()


# ---------------------------------------------------------------------------
# chromsize_to_s2c2s with a real analysis chromsize file
# ---------------------------------------------------------------------------


def test_chromsize_to_s2c2s_real_analysis_file(of, analysis_chromsize):
    s2c2s = of.chromsize_to_s2c2s(str(analysis_chromsize))
    # Both Hydra and Urchin should be present
    assert "Hydra" in s2c2s
    assert "Urchin" in s2c2s
    # Some chromosome under each, all-int values
    assert all(isinstance(v, int) for v in s2c2s["Hydra"].values())
    assert all(isinstance(v, int) for v in s2c2s["Urchin"].values())


# ---------------------------------------------------------------------------
# synteny_plot_sheet — end-to-end on real rbh + chromsize
# ---------------------------------------------------------------------------


def test_synteny_plot_sheet_e2e(sps, of, rbh_with_D_FET, analysis_chromsize, tmp_path):
    """Run the main synteny plot on real pipeline outputs. Verifies the
    function (a) doesn't crash and (b) writes both the plotdf TSV and
    the synteny PDF."""
    s2c2s = of.chromsize_to_s2c2s(str(analysis_chromsize))
    out_pdf = tmp_path / "synteny.pdf"
    out_tsv = tmp_path / "synteny_plotdf.tsv"
    # Run with no per-protein color mode and the default sort.
    sps.synteny_plot_sheet(
        df_file=str(rbh_with_D_FET),
        plotdf_file=str(out_tsv),
        synplot=str(out_pdf),
        xsample="Hydra", ysample="Urchin",
        s2c2s=s2c2s,
        xorder=[], yorder=[],
        xbreaks=[], ybreaks=[],
        orthology_method="diamond",
        sort_y_by_x=True,
    )
    assert out_pdf.is_file()
    assert out_pdf.stat().st_size > 0
    assert out_tsv.is_file()
    # Plotted RBH has the same shape as the input RBH plus plot-specific cols.
    plotdf = pd.read_csv(out_tsv, sep="\t")
    assert len(plotdf) > 0
    for col in ("Hydra_gene", "Urchin_gene", "color"):
        assert col in plotdf.columns


def test_synteny_plot_sheet_with_user_provided_order(
    sps, of, rbh_with_D_FET, analysis_chromsize, tmp_path
):
    """User-provided ``xorder``/``yorder`` (chrom-name lists) takes a
    different code branch than the auto-sort path. Run end-to-end."""
    s2c2s = of.chromsize_to_s2c2s(str(analysis_chromsize))
    # Use the first two chromosomes of each species per their s2c2s.
    xorder = list(s2c2s["Hydra"].keys())[:2]
    yorder = list(s2c2s["Urchin"].keys())[:2]
    out_pdf = tmp_path / "synteny_ordered.pdf"
    out_tsv = tmp_path / "synteny_ordered.tsv"
    sps.synteny_plot_sheet(
        df_file=str(rbh_with_D_FET),
        plotdf_file=str(out_tsv),
        synplot=str(out_pdf),
        xsample="Hydra", ysample="Urchin",
        s2c2s=s2c2s,
        xorder=xorder, yorder=yorder,
        xbreaks=[], ybreaks=[],
        orthology_method="diamond",
        sort_y_by_x=False,  # also exercises the no-auto-sort branch
    )
    assert out_pdf.is_file()


# ---------------------------------------------------------------------------
# odp_ribbon_plot.ribbon_plot — end-to-end with real pipeline rbh files
# ---------------------------------------------------------------------------


def test_ribbon_plot_two_species_e2e(rp, of, rbh_with_color, analysis_chromsize, tmp_path):
    """Run the top-level ribbon_plot orchestrator with a colored rbh
    (the only kind the function accepts) and a 2-species species_order.
    Verifies PDF generation end-to-end."""
    s2c2s = of.chromsize_to_s2c2s(str(analysis_chromsize))
    out_pdf = tmp_path / "ribbon.pdf"
    rp.ribbon_plot(
        species_order=["Hydra", "Urchin"],
        rbh_filelist=[str(rbh_with_color)],
        sp_to_chr_to_size=s2c2s,
        sp_min_chr_size={"Hydra": 0, "Urchin": 0},  # keep every scaffold
        outfile=str(out_pdf),
        sp_to_gene_order={},
        chr_sort_order="optimal-size",
        plot_all=False,
    )
    assert out_pdf.is_file()
    assert out_pdf.stat().st_size > 0


def test_ribbon_plot_with_None_gene_order(
    rp, of, rbh_with_color, analysis_chromsize, tmp_path
):
    """The sp_to_gene_order=None code path explicitly normalises to an
    empty dict; ensure that branch doesn't crash."""
    s2c2s = of.chromsize_to_s2c2s(str(analysis_chromsize))
    out_pdf = tmp_path / "ribbon2.pdf"
    rp.ribbon_plot(
        species_order=["Hydra", "Urchin"],
        rbh_filelist=[str(rbh_with_color)],
        sp_to_chr_to_size=s2c2s,
        sp_min_chr_size={"Hydra": 0, "Urchin": 0},
        outfile=str(out_pdf),
        sp_to_gene_order=None,
        chr_sort_order="custom",
    )
    assert out_pdf.is_file()
