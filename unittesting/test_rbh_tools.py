"""Tests for source/rbh_tools.py.

Covers the hex-color validator, the .rbh file parser (with dtype
enforcement), combine_rbh_db (db-mode merge), rbh_to_scafnum, and
guards/error paths in rbhdf_to_alglocdf.
"""
from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import pandas as pd
import pytest


@pytest.fixture(scope="module")
def rt(source_dir):
    import rbh_tools
    return rbh_tools


# ---------------------------------------------------------------------------
# hex_color_legal
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("good", [
    "#000000", "#FFFFFF", "#abcdef", "#012345", "#ABCDEF",
])
def test_hex_color_legal_accepts_well_formed(rt, good):
    assert rt.hex_color_legal(good) is True


@pytest.mark.parametrize("bad", [
    "",
    "#",
    "000000",       # missing #
    "#00000",       # too short
    "#0000000",     # too long
    "#GGGGGG",      # invalid hex char
    "#00 000",      # whitespace inside
    " #000000",     # leading whitespace
    "#000000 ",     # trailing whitespace
])
def test_hex_color_legal_rejects_malformed(rt, bad):
    assert rt.hex_color_legal(bad) is False


# ---------------------------------------------------------------------------
# parse_rbh — happy path, missing files, missing columns, dtype enforcement
# ---------------------------------------------------------------------------


def _write_minimal_rbh(tmp_path: Path, n=3) -> Path:
    p = tmp_path / "x.rbh"
    rows = [
        f"rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tsp2_gene\tsp2_scaf\tsp2_pos",
    ]
    for i in range(n):
        rows.append(f"rbh{i}\tNone\tg{i}\tchr1\t{i*100}\tx{i}\tscafA\t{i*150}")
    p.write_text("\n".join(rows) + "\n")
    return p


def test_parse_rbh_file_missing_raises(rt, tmp_path):
    with pytest.raises(IOError, match="does not exist"):
        rt.parse_rbh(tmp_path / "nope.rbh")


def test_parse_rbh_happy_path(rt, tmp_path):
    p = _write_minimal_rbh(tmp_path, n=3)
    df = rt.parse_rbh(p)
    assert len(df) == 3
    assert "rbh" in df.columns
    assert "sp1_gene" in df.columns
    assert "sp1_scaf" in df.columns
    assert "sp1_pos" in df.columns
    # sp1_pos is enforced to nullable Int64
    assert str(df["sp1_pos"].dtype) == "Int64"


def test_parse_rbh_missing_rbh_column(rt, tmp_path):
    p = tmp_path / "bad.rbh"
    p.write_text("foo\tbar\n1\t2\n")
    with pytest.raises(IOError, match="does not have a column named 'rbh'"):
        rt.parse_rbh(p)


def test_parse_rbh_missing_per_sample_columns(rt, tmp_path):
    """Has sp1_scaf but missing sp1_gene + sp1_pos."""
    p = tmp_path / "bad.rbh"
    p.write_text(
        "rbh\tgene_group\tsp1_scaf\n"
        "rbh1\tNone\tchr1\n"
    )
    with pytest.raises(IOError, match="mandatory columns are missing"):
        rt.parse_rbh(p)


def test_parse_rbh_missing_gene_group(rt, tmp_path):
    p = tmp_path / "bad.rbh"
    p.write_text(
        "rbh\tsp1_gene\tsp1_scaf\tsp1_pos\n"
        "rbh1\tg1\tchr1\t100\n"
    )
    with pytest.raises(IOError, match="mandatory columns are missing"):
        rt.parse_rbh(p)


def test_parse_rbh_legal_colors_pass(rt, tmp_path):
    p = tmp_path / "ok.rbh"
    p.write_text(dedent("""\
        rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos\tsp2_gene\tsp2_scaf\tsp2_pos
        rbh1\tNone\t#aabbcc\tg1\tchr1\t100\tx1\tscafA\t150
        rbh2\tNone\t#FFFFFF\tg2\tchr1\t200\tx2\tscafA\t250
    """))
    df = rt.parse_rbh(p)
    assert len(df) == 2


def test_parse_rbh_illegal_color_raises(rt, tmp_path):
    p = tmp_path / "badcolor.rbh"
    p.write_text(
        "rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos\n"
        "rbh1\tNone\tNOTAHEX\tg1\tchr1\t100\n"
    )
    with pytest.raises(IOError, match="not a legal hex color"):
        rt.parse_rbh(p)


def test_parse_rbh_NaN_in_color_column_is_tolerated(rt, tmp_path):
    """Empty cells in the color column should be treated as NaN, not
    rejected as illegal."""
    p = tmp_path / "ok.rbh"
    p.write_text(dedent("""\
        rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos
        rbh1\tNone\t#aabbcc\tg1\tchr1\t100
        rbh2\tNone\t\tg2\tchr1\t200
    """))
    df = rt.parse_rbh(p)
    assert len(df) == 2


def test_parse_rbh_FET_columns_are_float(rt, tmp_path):
    p = tmp_path / "fet.rbh"
    p.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tsp1_FET\tsp1_D
        rbh1\tNone\tg1\tchr1\t100\t0.001\t0.5
        rbh2\tNone\tg2\tchr1\t200\t0.05\t0.3
    """))
    df = rt.parse_rbh(p)
    assert df["sp1_FET"].dtype == float
    assert df["sp1_D"].dtype == float


def test_parse_rbh_breakchrom_column_is_object(rt, tmp_path):
    p = tmp_path / "br.rbh"
    p.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tsp1_breakchrom
        rbh1\tNone\tg1\tchr1\t100\tchr1
    """))
    df = rt.parse_rbh(p)
    assert df["sp1_breakchrom"].dtype == object


def test_parse_rbh_does_not_require_color_column(rt, tmp_path):
    p = _write_minimal_rbh(tmp_path)
    df = rt.parse_rbh(p)
    # No color column present; parse should still succeed.
    assert "color" not in df.columns or df["color"].notna().any()


# ---------------------------------------------------------------------------
# combine_rbh_db — db-mode merge
# ---------------------------------------------------------------------------


def _write_db_rbh(tmp_path: Path, name: str, rows: list[tuple]) -> Path:
    p = tmp_path / name
    lines = [
        "rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos\tsp2_gene\tsp2_scaf\tsp2_pos",
    ]
    for i, (rbh, group, color) in enumerate(rows):
        lines.append(f"{rbh}\t{group}\t{color}\tg{i}\tchr1\t{i}\tx{i}\tscafA\t{i}")
    p.write_text("\n".join(lines) + "\n")
    return p


def test_combine_rbh_db_basic(rt, tmp_path):
    a = _write_db_rbh(tmp_path, "a.rbh", [
        ("rbh_A1", "groupA", "#FF0000"),
        ("rbh_A2", "groupA", "#00FF00"),
    ])
    b = _write_db_rbh(tmp_path, "b.rbh", [
        ("rbh_B1", "groupB", "#0000FF"),
    ])
    df = rt.combine_rbh_db(a, b)
    assert len(df) == 3
    assert set(df["rbh"]) == {"rbh_A1", "rbh_A2", "rbh_B1"}
    assert set(df.columns) == {"rbh", "gene_group", "color"}


def test_combine_rbh_db_missing_file_raises(rt, tmp_path):
    a = _write_db_rbh(tmp_path, "a.rbh", [("rbh1", "g", "#FF0000")])
    with pytest.raises(IOError, match="does not exist"):
        rt.combine_rbh_db(a, tmp_path / "nonexistent.rbh")


def test_combine_rbh_db_duplicate_rbh_id_raises(rt, tmp_path):
    a = _write_db_rbh(tmp_path, "a.rbh", [("rbh_X", "g", "#FF0000")])
    b = _write_db_rbh(tmp_path, "b.rbh", [("rbh_X", "g", "#00FF00")])
    with pytest.raises(IOError, match="shared entries in the 'rbh' column"):
        rt.combine_rbh_db(a, b)


def test_combine_rbh_db_no_color_column_adds_black(rt, tmp_path):
    """A db-mode file without a color column should get '#000000' filled
    in during the merge."""
    a = tmp_path / "no_color.rbh"
    a.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos
        rbh1\tgA\tg1\tchr1\t100
    """))
    b = _write_db_rbh(tmp_path, "b.rbh", [("rbh2", "gB", "#0000FF")])
    df = rt.combine_rbh_db(a, b)
    assert df[df["rbh"] == "rbh1"]["color"].iloc[0] == "#000000"


# ---------------------------------------------------------------------------
# rbh_to_scafnum
# ---------------------------------------------------------------------------


def test_rbh_to_scafnum_counts_unique_scaffolds(rt, tmp_path):
    p = tmp_path / "many.rbh"
    p.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos
        rbh1\tNone\tg1\tchr1\t1
        rbh2\tNone\tg2\tchr1\t2
        rbh3\tNone\tg3\tchr2\t3
        rbh4\tNone\tg4\tchr3\t4
    """))
    df = rt.parse_rbh(p)
    assert rt.rbh_to_scafnum(df, "sp1") == 3


def test_rbh_to_scafnum_missing_sample_raises(rt, tmp_path):
    p = _write_minimal_rbh(tmp_path)
    df = rt.parse_rbh(p)
    with pytest.raises(IOError, match="does not have a column named"):
        rt.rbh_to_scafnum(df, "no_such_sample")


# ---------------------------------------------------------------------------
# rbhdf_to_alglocdf — guards
# ---------------------------------------------------------------------------


def test_rbhdf_to_alglocdf_rejects_filepath_input(rt, tmp_path):
    """Passing a path instead of a DataFrame is a common user mistake;
    the function detects it and points the user back to parse_rbh()."""
    p = _write_minimal_rbh(tmp_path)
    with pytest.raises(IOError, match="should be a pandas dataframe"):
        rt.rbhdf_to_alglocdf(str(p), 0.05, "sp1")


def test_rbhdf_to_alglocdf_rejects_unknown_ALG_name(rt, tmp_path):
    p = _write_minimal_rbh(tmp_path)
    df = rt.parse_rbh(p)
    df["whole_FET"] = 0.001  # mock the column otherwise filter fails
    with pytest.raises(IOError, match="not in the samples"):
        rt.rbhdf_to_alglocdf(df, 0.05, "BogusALG")


def test_rbhdf_to_alglocdf_rejects_ALG_substring_of_sample(rt, tmp_path):
    """If the ALG name is a substring of another sample name, the
    downstream column-splitting heuristics break — the function detects
    this up front. Sample names containing '_' are forbidden by
    odp_filechecker, so concatenation-style overlap is the case to test
    (e.g. 'BCN' substring of 'BCNHydra')."""
    p = tmp_path / "subname.rbh"
    p.write_text(dedent("""\
        rbh\tgene_group\tBCN_gene\tBCN_scaf\tBCN_pos\tBCNHydra_gene\tBCNHydra_scaf\tBCNHydra_pos
        rbh1\tNone\tg1\tchr1\t100\tx1\tscafA\t100
    """))
    df = rt.parse_rbh(p)
    df["whole_FET"] = 0.001
    with pytest.raises(IOError, match="is a part of the samplename"):
        rt.rbhdf_to_alglocdf(df, 0.05, "BCN")


def test_rbhdf_to_alglocdf_rejects_more_than_one_other_sample(rt, tmp_path):
    """Function expects exactly one species + the ALG. Three columns
    should be rejected with a help message about plot_LGs config."""
    p = tmp_path / "three.rbh"
    p.write_text(dedent("""\
        rbh\tgene_group\tALG_gene\tALG_scaf\tALG_pos\tspA_gene\tspA_scaf\tspA_pos\tspB_gene\tspB_scaf\tspB_pos
        rbh1\tNone\tg1\tchr1\t100\ta1\tscafA\t100\tb1\tscafB\t100
    """))
    df = rt.parse_rbh(p)
    df["whole_FET"] = 0.001
    with pytest.raises(IOError, match="more than one sample in the rbh file"):
        rt.rbhdf_to_alglocdf(df, 0.05, "ALG")


def test_rbhdf_to_alglocdf_missing_gene_group_raises(rt, tmp_path):
    """If the rbh DataFrame is missing the gene_group column, error."""
    df = pd.DataFrame({
        "rbh": ["r1"],
        "ALG_scaf": ["chr1"],
        "sp_scaf": ["scafA"],
    })
    with pytest.raises(IOError, match="does not have a column named 'gene_group'"):
        rt.rbhdf_to_alglocdf(df, 0.05, "ALG")  # noqa: E501


def test_rbhdf_to_alglocdf_happy_path(rt, tmp_path):
    """End-to-end with a fully-formed rbh DataFrame: returns
    (splitsdf, samplename) where splitsdf has rows per
    (gene_group, scaf, pvalue) group filtered by whole_FET."""
    p = tmp_path / "alg.rbh"
    p.write_text(dedent("""\
        rbh\tgene_group\twhole_FET\tALG_gene\tALG_scaf\tALG_pos\tsp_gene\tsp_scaf\tsp_pos
        r1\tA\t0.001\tg1\tchr1\t1\ts1\tscafX\t10
        r2\tA\t0.001\tg2\tchr1\t2\ts2\tscafX\t20
        r3\tA\t0.001\tg3\tchr1\t3\ts3\tscafX\t30
        r4\tB\t0.0001\tg4\tchr2\t4\ts4\tscafY\t40
        r5\tB\t0.0001\tg5\tchr2\t5\ts5\tscafY\t50
        r6\tA\t0.5\tg6\tchr1\t6\ts6\tscafZ\t60
    """))
    df = rt.parse_rbh(p)
    splitsdf, samplename = rt.rbhdf_to_alglocdf(df, 0.05, "ALG")
    assert samplename == "sp"
    # Only A/scafX and B/scafY survive the whole_FET <= 0.05 filter.
    assert len(splitsdf) == 2
    assert set(splitsdf["gene_group"]) == {"A", "B"}
    assert set(splitsdf["scaffold"]) == {"scafX", "scafY"}
    # frac_of_this_ALG_on_this_scaffold for A: 3/4 (3 on scafX of 4 A-rows total).
    a_row = splitsdf[splitsdf["gene_group"] == "A"].iloc[0]
    assert a_row["frac_of_this_ALG_on_this_scaffold"] == pytest.approx(0.75)


def test_rbhdf_to_alglocdf_inconsistent_FET_raises(rt, tmp_path):
    """Same (gene_group, samplescaf) pair with different whole_FET values
    indicates upstream FET ran incorrectly — function detects and
    rejects."""
    df = pd.DataFrame({
        "rbh": ["r1", "r2"],
        "gene_group": ["A", "A"],
        "whole_FET": [0.001, 0.002],  # different on same group/scaf
        "ALG_scaf": ["chr1", "chr1"],
        "sp_scaf": ["scafX", "scafX"],
    })
    with pytest.raises(IOError, match="same value for all rows"):
        rt.rbhdf_to_alglocdf(df, 0.05, "ALG")


def test_combine_rbh_strips_FET_D_columns_from_inputs(rt, tmp_path):
    """Lines 252-255: combine_rbh drops columns ending in _FET, _D,
    _ix etc. from both inputs because the merged frame can't keep stale
    per-input stats."""
    a = tmp_path / "a.rbh"
    a.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tspA_gene\tspA_scaf\tspA_pos\tsp1_FET\tsp1_D
        ra1\tNone\tg1\tchr1\t100\ta1\tscafA\t10\t0.001\t0.5
        ra2\tNone\tg2\tchr1\t200\ta2\tscafA\t20\t0.002\t0.3
    """))
    b = tmp_path / "b.rbh"
    b.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tspB_gene\tspB_scaf\tspB_pos\tsp1_ix
        rb1\tNone\tg3\tchr1\t300\tb1\tscafB\t30\t5
    """))
    df = rt.combine_rbh(a, b)
    # Stat columns get stripped in the merged output.
    for col in df.columns:
        for stale_suffix in ("_FET", "_D", "_ix"):
            assert not col.endswith(stale_suffix), f"stale {col} should have been dropped"


def test_combine_rbh_rejects_shared_rbh_ids_between_files(rt, tmp_path):
    """Line 285-286: same rbh id in both files is rejected — the merge
    is meant to bring in disjoint marker sets."""
    a = tmp_path / "a.rbh"
    a.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tspA_gene\tspA_scaf\tspA_pos
        rXX\tNone\tg1\tchr1\t100\ta1\tscafA\t10
    """))
    b = tmp_path / "b.rbh"
    b.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tspB_gene\tspB_scaf\tspB_pos
        rXX\tNone\tg2\tchr1\t200\tb1\tscafB\t30
    """))
    with pytest.raises(IOError, match="shared entries"):
        rt.combine_rbh(a, b)


def test_combine_rbh_error_paths_do_not_crash_in_format_string(rt, tmp_path):
    """Regression: the error path for `df1_unique_list != 1` previously
    referenced the unbound name `df1_unique`. We can't directly trigger
    that branch through `combine_rbh` because the surrounding shared-sample
    check fires first, but the fix makes sure the f-string is well-formed."""
    a = tmp_path / "a.rbh"
    a.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos
        r1\tNone\tg1\tchr1\t1
    """))
    b = tmp_path / "b.rbh"
    b.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos
        r2\tNone\tg2\tchr1\t2
    """))
    # No non-shared sample at all — fails on shared check, not on _unique check.
    with pytest.raises(IOError):
        rt.combine_rbh(a, b)


# ---------------------------------------------------------------------------
# parse_ALG_rbh_to_colordf — ALG-color/size summary
# ---------------------------------------------------------------------------


def test_parse_ALG_rbh_to_colordf_basic(rt, tmp_path):
    p = tmp_path / "alg.rbh"
    p.write_text(dedent("""\
        rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos
        r1\tA1\t#FF0000\tg1\tchr1\t1
        r2\tA1\t#FF0000\tg2\tchr1\t2
        r3\tA1\t#FF0000\tg3\tchr1\t3
        r4\tB2\t#00FF00\tg4\tchr2\t4
        r5\tB2\t#00FF00\tg5\tchr2\t5
    """))
    df = rt.parse_ALG_rbh_to_colordf(p)
    assert list(df.columns) == ["ALGname", "Color", "Size"]
    # Sorted by increasing Size — B2 (size 2) before A1 (size 3).
    assert df.iloc[0]["ALGname"] == "B2"
    assert df.iloc[0]["Size"] == 2
    assert df.iloc[1]["ALGname"] == "A1"
    assert df.iloc[1]["Size"] == 3


def test_parse_ALG_rbh_to_colordf_picks_most_common_color(rt, tmp_path):
    """If an ALG has rows with multiple colors, picks the most common."""
    p = tmp_path / "alg.rbh"
    p.write_text(dedent("""\
        rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos
        r1\tA1\t#FF0000\tg1\tchr1\t1
        r2\tA1\t#FF0000\tg2\tchr1\t2
        r3\tA1\t#00FF00\tg3\tchr1\t3
    """))
    df = rt.parse_ALG_rbh_to_colordf(p)
    assert df.iloc[0]["Color"] == "#FF0000"


def test_parse_ALG_rbh_to_colordf_missing_columns_raises(rt, tmp_path):
    """If the file is missing `gene_group` or `color`, error."""
    p = tmp_path / "bad.rbh"
    p.write_text("rbh\tsp1_gene\tsp1_scaf\tsp1_pos\nr1\tg1\tchr1\t1\n")
    with pytest.raises(IOError, match="does not have the correct columns"):
        rt.parse_ALG_rbh_to_colordf(p)


# ---------------------------------------------------------------------------
# combine_rbh — full RBH merge (more than just rbh/gene_group/color)
# ---------------------------------------------------------------------------


def test_combine_rbh_basic(rt, tmp_path):
    """Combine two rbh files with a shared species (sp1) sharing
    scaffold names — the merged output contains all rbh entries plus
    merged species columns."""
    a = tmp_path / "a.rbh"
    a.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tspA_gene\tspA_scaf\tspA_pos
        ra1\tNone\tg1\tchr1\t100\ta1\tscafA\t10
        ra2\tNone\tg2\tchr1\t200\ta2\tscafA\t20
    """))
    b = tmp_path / "b.rbh"
    # Shares sp1 scaffold name `chr1` so combine_rbh's "shared scaf"
    # consistency check passes.
    b.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tspB_gene\tspB_scaf\tspB_pos
        rb1\tNone\tg3\tchr1\t300\tb1\tscafB\t30
    """))
    df = rt.combine_rbh(a, b)
    assert len(df) >= 2
    # Output retains the shared species (sp1) and unifies the two
    # non-shared species into a single `MergedCol_*` column group.
    assert "sp1_gene" in df.columns
    assert any(c.startswith("MergedCol") for c in df.columns)


def test_combine_rbh_no_shared_scaffolds_raises(rt, tmp_path):
    """If the two rbh files don't share scaffold names on the common
    species, combine_rbh refuses to silently merge incompatible data."""
    a = tmp_path / "a.rbh"
    a.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tspA_gene\tspA_scaf\tspA_pos
        ra1\tNone\tg1\tchr1\t100\ta1\tscafA\t10
    """))
    b = tmp_path / "b.rbh"
    # Different sp1 scaffold name.
    b.write_text(dedent("""\
        rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\tspB_gene\tspB_scaf\tspB_pos
        rb1\tNone\tg3\tchr_other\t300\tb1\tscafB\t30
    """))
    with pytest.raises(IOError, match="no shared entries"):
        rt.combine_rbh(a, b)
