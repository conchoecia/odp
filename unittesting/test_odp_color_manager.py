"""Tests for source/odp_color_manager.py.

Focus on the pure-function color utilities — hex↔rgb conversion, hex
validation, contrasting-color picker, brightness / saturation math,
random / palette color generators. These run without any I/O, so they
exercise the helpers used everywhere downstream (synteny plots, ribbon
plots, ALG colorings).
"""
from __future__ import annotations

import pytest


@pytest.fixture(scope="module")
def cm(source_dir):
    """Import source/odp_color_manager.py via the source_dir fixture."""
    import odp_color_manager
    return odp_color_manager


# ---------------------------------------------------------------------------
# h2r — hex string → normalized RGB triple (each in [0,1])
# ---------------------------------------------------------------------------


def test_h2r_black(cm):
    r, g, b = cm.h2r("#000000")
    assert (r, g, b) == (0.0, 0.0, 0.0)


def test_h2r_white(cm):
    r, g, b = cm.h2r("#FFFFFF")
    assert r == pytest.approx(1.0)
    assert g == pytest.approx(1.0)
    assert b == pytest.approx(1.0)


def test_h2r_red(cm):
    r, g, b = cm.h2r("#FF0000")
    assert r == pytest.approx(1.0)
    assert g == 0.0
    assert b == 0.0


def test_h2r_handles_lowercase(cm):
    assert cm.h2r("#ff0000") == cm.h2r("#FF0000")


def test_h2r_handles_no_hash_prefix(cm):
    assert cm.h2r("ff0000") == cm.h2r("#FF0000")


def test_h2r_known_value(cm):
    # 0x33aa66 -> (51, 170, 102) -> /255 each
    r, g, b = cm.h2r("#33aa66")
    assert r == pytest.approx(51 / 255)
    assert g == pytest.approx(170 / 255)
    assert b == pytest.approx(102 / 255)


# ---------------------------------------------------------------------------
# inverse_color — given a background hex, returns "#000000" or "#ffffff"
# ---------------------------------------------------------------------------


def test_inverse_color_black_returns_black(cm):
    """Black-on-black is the documented short-circuit; lets callers know
    no text was rendered."""
    assert cm.inverse_color("#000000") == "#000000"


def test_inverse_color_white_returns_black(cm):
    """Bright background → dark text."""
    assert cm.inverse_color("#FFFFFF") == "#000000"


def test_inverse_color_dark_returns_white(cm):
    """Dark background → light text."""
    assert cm.inverse_color("#222222") == "#ffffff"


def test_inverse_color_red_returns_white(cm):
    # red has weighted brightness 255*0.299 ≈ 76 → < 186 → white
    assert cm.inverse_color("#FF0000") == "#ffffff"


def test_inverse_color_yellow_returns_black(cm):
    # yellow has weighted brightness ≈ 226 → > 186 → black
    assert cm.inverse_color("#FFFF00") == "#000000"


# ---------------------------------------------------------------------------
# is_valid_hex_code — strict 7-char `#XXXXXX` form
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("good", [
    "#000000", "#FFFFFF", "#abcdef", "#012345", "#ABCDEF",
])
def test_is_valid_hex_code_accepts_well_formed(cm, good):
    assert cm.is_valid_hex_code(good) is True


@pytest.mark.parametrize("bad", [
    "",            # empty
    "#",           # just hash
    "000000",      # missing hash
    "#00000",      # 6 chars total — too short
    "#0000000",    # 8 chars total — too long
    "#GGGGGG",     # non-hex chars
    "#00 000",     # whitespace inside
    "  #000000",   # leading whitespace
    "#000000 ",    # trailing whitespace
    None,          # not a string at all
])
def test_is_valid_hex_code_rejects_malformed(cm, bad):
    if bad is None:
        # Function uses len(); None has no len → TypeError. Accept that
        # too as "rejection".
        with pytest.raises(TypeError):
            cm.is_valid_hex_code(bad)
    else:
        assert cm.is_valid_hex_code(bad) is False


# ---------------------------------------------------------------------------
# calculate_brightness — int color → weighted brightness
# ---------------------------------------------------------------------------


def test_calculate_brightness_black_is_zero(cm):
    assert cm.calculate_brightness(0x000000) == 0


def test_calculate_brightness_white_at_max(cm):
    # 255 * (0.299 + 0.587 + 0.114) = 255 → weighted sum / 1000 = ...
    # Actually formula: (r*299 + g*587 + b*114) / 1000
    expected = (255 * 299 + 255 * 587 + 255 * 114) / 1000
    assert cm.calculate_brightness(0xFFFFFF) == pytest.approx(expected)


def test_calculate_brightness_red(cm):
    expected = (255 * 299) / 1000
    assert cm.calculate_brightness(0xFF0000) == pytest.approx(expected)


def test_calculate_brightness_returns_positive(cm):
    """Any non-zero color should give a positive brightness."""
    assert cm.calculate_brightness(0x010101) > 0


def test_calculate_brightness_is_monotone_in_each_channel(cm):
    """Raising any single channel must not decrease brightness."""
    base = cm.calculate_brightness(0x808080)
    redder = cm.calculate_brightness(0xFF8080)
    greener = cm.calculate_brightness(0x80FF80)
    bluer = cm.calculate_brightness(0x8080FF)
    assert redder >= base
    assert greener >= base
    assert bluer >= base


# ---------------------------------------------------------------------------
# calculate_saturation — int color → [0, 1]
# ---------------------------------------------------------------------------


def test_calculate_saturation_grayscale_is_zero(cm):
    """Gray colors have zero saturation (max == min)."""
    assert cm.calculate_saturation(0x000000) == 0
    assert cm.calculate_saturation(0x808080) == 0
    assert cm.calculate_saturation(0xFFFFFF) == 0


def test_calculate_saturation_pure_red_is_one(cm):
    """Pure red: max=255, min=0 → saturation = 1."""
    assert cm.calculate_saturation(0xFF0000) == pytest.approx(1.0)


def test_calculate_saturation_in_unit_range(cm):
    for c in [0x123456, 0xABCDEF, 0x010101, 0xFEDCBA]:
        s = cm.calculate_saturation(c)
        assert 0.0 <= s <= 1.0


def test_calculate_saturation_black_handles_zero_max(cm):
    """The implementation special-cases max==0 to avoid /0."""
    assert cm.calculate_saturation(0x000000) == 0


# ---------------------------------------------------------------------------
# return_random_color — hex string with leading "#"
# ---------------------------------------------------------------------------


def test_return_random_color_format(cm):
    import re
    for _ in range(20):
        c = cm.return_random_color()
        assert isinstance(c, str)
        assert re.match(r"^#[0-9a-f]{6}$", c), f"bad format: {c!r}"


# ---------------------------------------------------------------------------
# generate_random_color — bounded brightness + saturation
# ---------------------------------------------------------------------------


def test_generate_random_color_is_within_bounds(cm):
    """Sample 30 colors and verify each falls inside the brightness and
    saturation windows the function advertises."""
    for _ in range(30):
        hex_c = cm.generate_random_color()
        assert hex_c.startswith("#") and len(hex_c) == 7
        n = int(hex_c.lstrip("#"), 16)
        b = cm.calculate_brightness(n)
        s = cm.calculate_saturation(n)
        assert 127 < b < 300
        assert 0.3 < s < 0.8


# ---------------------------------------------------------------------------
# yield_color — deterministic palette generator
# ---------------------------------------------------------------------------


def test_yield_color_mypals25_cycles_through_all_25(cm):
    gen = cm.yield_color("mypals25")
    palette = [next(gen) for _ in range(25)]
    assert len(set(palette)) == 25  # all distinct on first cycle
    # 26th value loops back to first
    assert next(gen) == palette[0]


def test_yield_color_tableau20_default_length(cm):
    gen = cm.yield_color("tableau20")
    palette = [next(gen) for _ in range(20)]
    assert len(set(palette)) == 20
    assert next(gen) == palette[0]


def test_yield_color_unknown_palette_raises(cm):
    with pytest.raises(KeyError):
        # KeyError comes from the dict lookup before the generator yields.
        next(cm.yield_color("does-not-exist"))


def test_yield_color_each_entry_is_valid_hex(cm):
    gen = cm.yield_color("mypals25")
    import re
    for _ in range(25):
        c = next(gen)
        assert re.match(r"^#[0-9A-Fa-f]{6}$", c), f"bad palette entry: {c!r}"


# ---------------------------------------------------------------------------
# LG_db class — directory validation + duplicate-species detection
# ---------------------------------------------------------------------------


def _write_min_lg_db(tmp_path, rbh_content, hmm_content):
    """Build a fake LG-database directory with one .rbh + one .hmm."""
    d = tmp_path / "minLG"
    d.mkdir()
    (d / "minLG.rbh").write_text(rbh_content)
    (d / "minLG.hmm").write_text(hmm_content)
    return d


def test_LG_db_missing_directory_raises(cm, tmp_path):
    with pytest.raises(IOError, match="does not exist"):
        cm.LG_db("missing", str(tmp_path / "no_such_dir"), [])


def test_LG_db_no_rbh_file_raises(cm, tmp_path):
    d = tmp_path / "norbh"
    d.mkdir()
    (d / "x.hmm").write_text("")
    with pytest.raises(IOError, match="There must be a single .rbh"):
        cm.LG_db("name", str(d), [])


def test_LG_db_no_hmm_file_raises(cm, tmp_path):
    d = tmp_path / "nohmm"
    d.mkdir()
    (d / "x.rbh").write_text("rbh\tgene_group\tcolor\n")
    with pytest.raises(IOError, match="There must be a single .hmm"):
        cm.LG_db("name", str(d), [])


def test_LG_db_two_rbh_files_raises(cm, tmp_path):
    d = tmp_path / "tworbh"
    d.mkdir()
    (d / "a.rbh").write_text("rbh\tgene_group\tcolor\n")
    (d / "b.rbh").write_text("rbh\tgene_group\tcolor\n")
    (d / "x.hmm").write_text("")
    with pytest.raises(IOError, match="There must be a single .rbh"):
        cm.LG_db("name", str(d), [])


def test_LG_db_basic_instantiation(cm, tmp_path):
    """Build a minimal but valid LG_db with one ALG and one species, and
    confirm the lookups are populated."""
    rbh = (
        "rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos\n"
        "rbh1\tA1\t#FF0000\tg1\tchrA\t100\n"
        "rbh2\tA1\t#FF0000\tg2\tchrA\t200\n"
        "rbh3\tB2\t#00FF00\tg3\tchrB\t300\n"
    )
    # Single ALG model named `rbh1` in the HMM file (hmmer NAME line).
    hmm = "NAME  rbh1\nNAME  rbh2\nNAME  rbh3\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)
    obj = cm.LG_db("minLG", str(d), [])
    assert obj.name == "minLG"
    assert obj.rbhfile.endswith("minLG.rbh")
    assert obj.hmmfile.endswith("minLG.hmm")
    # group_to_color picks the most-common color per group
    assert obj.group_to_color["A1"] == "#FF0000"
    assert obj.group_to_color["B2"] == "#00FF00"
    # None/NA defaults to black for downstream plotting
    assert obj.group_to_color["None"] == "#000000"
    assert obj.group_to_color["NA"] == "#000000"


def test_LG_db_rejects_mismatched_rbh_and_hmm(cm, tmp_path):
    """The .hmm file must contain a NAME entry for every rbh in the
    .rbh file. Missing entries raise."""
    rbh = (
        "rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos\n"
        "rbh1\tA1\t#FF0000\tg1\tchrA\t100\n"
        "rbh2\tA1\t#FF0000\tg2\tchrA\t200\n"
    )
    # Only rbh1 has a NAME entry — rbh2 is missing.
    hmm = "NAME  rbh1\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)
    with pytest.raises(IOError, match="contains some orthologs not in the hmm"):
        cm.LG_db("minLG", str(d), [])


def test_LG_db_rejects_missing_required_column(cm, tmp_path):
    """rbh file missing one of rbh/gene_group/color raises."""
    rbh = (
        # color column missing
        "rbh\tgene_group\tsp1_gene\tsp1_scaf\tsp1_pos\n"
        "rbh1\tA1\tg1\tchrA\t100\n"
    )
    hmm = "NAME  rbh1\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)
    with pytest.raises(IOError, match="must have the following column"):
        cm.LG_db("minLG", str(d), [])


def test_LG_db_parse_rbhdf_removes_duplicate_species_columns(cm, tmp_path):
    """When two species columns have identical gene IDs (e.g., one was
    annotated using the other's proteome), the method should drop one."""
    rbh = (
        "rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos\t"
        "sp2_gene\tsp2_scaf\tsp2_pos\n"
        "rbh1\tA\t#FF0000\tg1\tchrA\t100\tg1\tscafX\t100\n"
        "rbh2\tA\t#FF0000\tg2\tchrA\t200\tg2\tscafX\t200\n"
    )
    hmm = "NAME  rbh1\nNAME  rbh2\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)
    obj = cm.LG_db("minLG", str(d), [])
    # After dedup, only one of sp1/sp2 should survive (alphabetically
    # earlier wins per the docstring).
    survived = [c for c in obj.rbhdf.columns if c.endswith("_gene")]
    assert len(survived) == 1
    assert survived[0] == "sp1_gene"


# ---------------------------------------------------------------------------
# LG_db._gen_sp_to_gene_to_color + _gen_sp_to_gene_to_group
# ---------------------------------------------------------------------------


def test_LG_db_sp_to_gene_to_color_built_from_rbh(cm, tmp_path):
    rbh = (
        "rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos\n"
        "rbh1\tA1\t#FF0000\tg1\tchrA\t100\n"
        "rbh2\tA1\t#FF0000\tg2\tchrA\t200\n"
    )
    hmm = "NAME  rbh1\nNAME  rbh2\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)
    obj = cm.LG_db("minLG", str(d), [])
    assert obj.sp_to_gene_to_color["sp1"]["g1"] == "#FF0000"
    assert obj.sp_to_gene_to_color["sp1"]["g2"] == "#FF0000"


def test_LG_db_sp_to_gene_to_group_built_from_rbh(cm, tmp_path):
    rbh = (
        "rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos\n"
        "rbh1\tA1\t#FF0000\tg1\tchrA\t100\n"
        "rbh2\tB2\t#00FF00\tg2\tchrA\t200\n"
    )
    hmm = "NAME  rbh1\nNAME  rbh2\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)
    obj = cm.LG_db("minLG", str(d), [])
    assert obj.sp_to_gene_to_group["sp1"]["g1"] == "A1"
    assert obj.sp_to_gene_to_group["sp1"]["g2"] == "B2"


def test_LG_db_invalid_color_in_row_falls_back_to_black(cm, tmp_path):
    """When the rbh `color` cell doesn't look like a hex code (no
    '#'), the per-gene lookup uses black as a fallback."""
    rbh = (
        "rbh\tgene_group\tcolor\tsp1_gene\tsp1_scaf\tsp1_pos\n"
        "rbh1\tA1\tnotacolor\tg1\tchrA\t100\n"
    )
    # This should fail because hex_color_legal would reject "notacolor"
    # — but actually LG_db doesn't use hex_color_legal; it just checks
    # for "#" in the cell, so "notacolor" becomes #000000.
    hmm = "NAME  rbh1\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)
    obj = cm.LG_db("minLG", str(d), [])
    assert obj.sp_to_gene_to_color["sp1"]["g1"] == "#000000"


# ---------------------------------------------------------------------------
# LG_db._sp_matches_which_db_species
# ---------------------------------------------------------------------------


def test_LG_db_sp_matches_returns_match_when_majority_overlap(cm, tmp_path):
    """When 10:1 majority of input protein IDs match one DB species,
    that species should be returned."""
    rbh = (
        "rbh\tgene_group\tcolor\tdbsp_gene\tdbsp_scaf\tdbsp_pos\n"
        + "\n".join(
            f"rbh{i}\tA1\t#FF0000\tprot{i}\tchrA\t{i}"
            for i in range(20)
        ) + "\n"
    )
    hmm = "\n".join(f"NAME  rbh{i}" for i in range(20)) + "\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)
    obj = cm.LG_db("minLG", str(d), [])
    # Pass in IDs matching the DB species' proteins.
    match = obj._sp_matches_which_db_species(
        [f"prot{i}" for i in range(20)]
    )
    assert match == "dbsp"


def test_LG_db_sp_matches_returns_none_when_no_match(cm, tmp_path):
    rbh = (
        "rbh\tgene_group\tcolor\tdbsp_gene\tdbsp_scaf\tdbsp_pos\n"
        "rbh1\tA1\t#FF0000\tprot1\tchrA\t1\n"
    )
    hmm = "NAME  rbh1\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)
    obj = cm.LG_db("minLG", str(d), [])
    assert obj._sp_matches_which_db_species(["totally_unrelated"]) is None


# ---------------------------------------------------------------------------
# LG_db.color_dataframe via the HMM-results path
# ---------------------------------------------------------------------------


def test_LG_db_color_dataframe_via_hmm_results(cm, tmp_path):
    """If the plot species' protein IDs don't match any DB species, the
    function falls back to using HMM results to identify orthologs."""
    rbh = (
        "rbh\tgene_group\tcolor\tdbsp_gene\tdbsp_scaf\tdbsp_pos\n"
        "rbh1\tA1\t#FF0000\torth_a\tchrA\t1\n"
        "rbh2\tB2\t#00FF00\torth_b\tchrA\t2\n"
    )
    hmm = "NAME  rbh1\nNAME  rbh2\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)

    # Build an HMM results file (12-column blastp-like format).
    # column 1 = rbh id (qseqid); column 2 = species' protein id (sseqid).
    hmm_results = tmp_path / "results.tsv"
    hmm_results.write_text(
        "rbh1\tspecies_prot_1\t99\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200\n"
        "rbh2\tspecies_prot_2\t99\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200\n"
    )
    obj = cm.LG_db("minLG", str(d), [str(hmm_results)])
    # The HMM-prot-to-rbh map should be populated.
    assert obj.hmm_prot_to_rbh == {
        "species_prot_1": "rbh1",
        "species_prot_2": "rbh2",
    }

    import pandas as pd
    plotdf = pd.DataFrame({
        "plotsp_gene": ["species_prot_1", "species_prot_2", "unknown_prot"],
        "plotsp_scaf": ["chrZ", "chrZ", "chrZ"],
        "plotsp_pos":  [1, 2, 3],
        "color":       ["", "", ""],
        "gene_group":  ["", "", ""],
    })
    result = obj.color_dataframe(plotdf)
    # gene_group + color come from the rbh file via the hmm lookup
    assert result.iloc[0]["gene_group"] == "A1"
    assert result.iloc[0]["color"] == "#FF0000"
    assert result.iloc[1]["gene_group"] == "B2"
    assert result.iloc[1]["color"] == "#00FF00"
    # Unknown protein stays black/None
    assert result.iloc[2]["color"] == "#000000"
    assert result.iloc[2]["gene_group"] == "None"


def test_LG_db_color_dataframe_via_direct_species_match(cm, tmp_path):
    """If the plot species' protein IDs match a DB species' IDs
    directly, the function uses the direct path (no HMM lookup)."""
    rbh = (
        "rbh\tgene_group\tcolor\tdbsp_gene\tdbsp_scaf\tdbsp_pos\n"
        + "\n".join(
            f"rbh{i}\tA1\t#FF0000\tprot{i}\tchrA\t{i}"
            for i in range(20)
        ) + "\n"
    )
    hmm = "\n".join(f"NAME  rbh{i}" for i in range(20)) + "\n"
    d = _write_min_lg_db(tmp_path, rbh, hmm)
    obj = cm.LG_db("minLG", str(d), [])

    import pandas as pd
    plotdf = pd.DataFrame({
        "plotsp_gene": ["prot0", "prot1", "prot2"],
        "plotsp_scaf": ["chrZ", "chrZ", "chrZ"],
        "plotsp_pos":  [1, 2, 3],
        "color":       ["", "", ""],
        "gene_group":  ["", "", ""],
    })
    result = obj.color_dataframe(plotdf)
    assert result.iloc[0]["gene_group"] == "A1"
    assert result.iloc[0]["color"] == "#FF0000"
    # color_method records which path was taken
    assert "protein ids" in obj.color_method.lower()
