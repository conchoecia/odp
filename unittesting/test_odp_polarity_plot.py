"""Tests for source/odp_polarity_plot.py.

The polarity plotter is mostly rendering code, so these tests use tiny
in-memory-style TSV fixtures and inspect the matplotlib objects it adds.
"""
from __future__ import annotations

import math

import pandas as pd
import pytest


@pytest.fixture(scope="module")
def pp(source_dir):
    import odp_polarity_plot

    return odp_polarity_plot


def _write_tsv(path, rows):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def _rbh_paths(tmp_path, tm_rows, bm_rows):
    top_mid = _write_tsv(tmp_path / "top_mid.rbh", tm_rows)
    bot_mid = _write_tsv(tmp_path / "bot_mid.rbh", bm_rows)
    return top_mid, bot_mid


def _base_tm_rows():
    return [
        {
            "Top_gene": "t1",
            "Top_scaf": "eupHapAv0.3_chr1",
            "Top_pos": 100,
            "Mid_gene": "m1",
            "Mid_scaf": "Mid_chrA",
            "Mid_pos": 20,
        },
        {
            "Top_gene": "t2",
            "Top_scaf": "HiC_scaffold_22",
            "Top_pos": 200,
            "Mid_gene": "m2",
            "Mid_scaf": "Mid_chrB",
            "Mid_pos": 30,
        },
        {
            "Top_gene": "t3",
            "Top_scaf": "HiC_scaffold_22",
            "Top_pos": 50,
            "Mid_gene": "m3",
            "Mid_scaf": "Mid_chrA",
            "Mid_pos": 10,
        },
    ]


def _base_bm_rows():
    return [
        {
            "Bot_gene": "b1",
            "Bot_scaf": "Bot_fused",
            "Bot_pos": 1000,
            "Mid_gene": "m1",
            "Mid_scaf": "Mid_chrA",
            "Mid_pos": 20,
        },
        {
            "Bot_gene": "b2",
            "Bot_scaf": "Bot_fused",
            "Bot_pos": 2000,
            "Mid_gene": "m2",
            "Mid_scaf": "Mid_chrB",
            "Mid_pos": 30,
        },
        {
            "Bot_gene": "b3",
            "Bot_scaf": "Bot_fused",
            "Bot_pos": 3000,
            "Mid_gene": "m3",
            "Mid_scaf": "Mid_chrA",
            "Mid_pos": 10,
        },
    ]


def _draw_panel(pp, tmp_path, tm_rows=None, bm_rows=None, **kwargs):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top_mid, bot_mid = _rbh_paths(
        tmp_path,
        tm_rows if tm_rows is not None else _base_tm_rows(),
        bm_rows if bm_rows is not None else _base_bm_rows(),
    )
    fig, ax = plt.subplots()
    stats = pp.draw_polarity_panel(
        ax,
        top_sp="Top",
        top_chrA="eupHapAv0.3_chr1",
        top_chrB="HiC_scaffold_22",
        mid_sp="Mid",
        mid_chrA="Mid_chrA",
        mid_chrB="Mid_chrB",
        bot_sp="Bot",
        bot_chr_fused="Bot_fused",
        rbh_top_mid=top_mid,
        rbh_bot_mid=bot_mid,
        **kwargs,
    )
    return fig, ax, stats


@pytest.mark.parametrize(
    ("chrom", "expected"),
    [
        ("HiC_scaffold_12.1", "HiC12"),
        ("NW_123456789.2", "NW789"),
        ("NC_000001.11", "NC01"),
        ("OZ12345678.1", "OZ78"),
        ("CM98765432.1", "CM32"),
        ("eupHapAv0.3_chr14", "chr14"),
        ("asm_chrom5", "chr5"),
        ("plain_scaffold", "plain_scaffold"),
    ],
)
def test_chrom_acronym_common_accession_patterns(pp, chrom, expected):
    assert pp.chrom_acronym(chrom) == expected


def test_layout_chroms_handles_zero_width_chromosomes(pp):
    layout = pp.layout_chroms(["chrA", "chrB"], {"chrA": 0, "chrB": 0}, 0, 1)

    assert set(layout) == {"chrA", "chrB"}
    for x0, x1 in layout.values():
        assert math.isfinite(x0)
        assert math.isfinite(x1)
        assert x0 <= x1


def test_gene_x_single_gene_chromosome_uses_midpoint(pp):
    layout = {"chr1": (0.2, 0.8)}

    assert pp.gene_x("chr1", pos_index=99, n_on_chrom=1, layout=layout) == 0.5


def test_load_bcns_map_missing_malformed_and_valid_files(pp, tmp_path):
    missing = tmp_path / "missing.rbh"
    assert pp.load_bcns_map(missing, "Top_gene") == {}

    malformed = _write_tsv(
        tmp_path / "malformed.rbh",
        [{"Top_gene": "t1", "wrong_col": "D"}],
    )
    assert pp.load_bcns_map(malformed, "Top_gene") == {}

    valid = _write_tsv(
        tmp_path / "valid.rbh",
        [{"Top_gene": "t1", "BCnSSimakov2022_scaf": "D"}],
    )
    assert pp.load_bcns_map(valid, "Top_gene") == {"t1": "D"}


def test_draw_panel_empty_triplet_join_returns_no_data_stats(pp, tmp_path):
    tm_rows = [_base_tm_rows()[0]]
    bm_rows = [{**_base_bm_rows()[0], "Mid_gene": "not_shared"}]

    fig, ax, stats = _draw_panel(pp, tmp_path, tm_rows, bm_rows, color_by="raw")
    try:
        assert stats == {"n_triplets": 0, "n_unmapped": 0}
        assert len(ax.patches) == 0
        assert any("no 3-way triplets" in t.get_text() for t in ax.texts)
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)


def test_draw_panel_lg2_singletons_and_chrom_label_overrides(pp, tmp_path):
    fig, ax, stats = _draw_panel(
        pp,
        tmp_path,
        color_by="lg2",
        chrom_acronyms=True,
        chrom_label_map={
            "eupHapAv0.3_chr1": "mapped-A",
            "Bot_fused": "mapped-fused",
        },
        lg_a_label="side A",
        lg_b_label="side B",
    )
    try:
        text = {t.get_text() for t in ax.texts}
        legend_text = {t.get_text() for t in ax.get_legend().get_texts()}

        assert stats == {"n_triplets": 3, "n_unmapped": 0}
        assert len(ax.patches) == 6
        assert "mapped-A" in text
        assert "mapped-fused" in text
        assert "HiC22" in text
        assert legend_text == {"side A", "side B"}
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)


def test_draw_panel_bcns_requires_palette_for_nonempty_panel(pp, tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top_mid, bot_mid = _rbh_paths(tmp_path, _base_tm_rows(), _base_bm_rows())
    fig, ax = plt.subplots()
    try:
        with pytest.raises(ValueError, match="bcns_palette path is required"):
            pp.draw_polarity_panel(
                ax,
                top_sp="Top",
                top_chrA="eupHapAv0.3_chr1",
                top_chrB="HiC_scaffold_22",
                mid_sp="Mid",
                mid_chrA="Mid_chrA",
                mid_chrB="Mid_chrB",
                bot_sp="Bot",
                bot_chr_fused="Bot_fused",
                rbh_top_mid=top_mid,
                rbh_bot_mid=bot_mid,
                color_by="bcns",
            )
    finally:
        plt.close(fig)


def test_draw_panel_bcns_precedence_fallback_and_unmapped_count(pp, tmp_path):
    import matplotlib.colors as mcolors

    tm_rows = _base_tm_rows() + [
        {
            "Top_gene": "t4",
            "Top_scaf": "eupHapAv0.3_chr1",
            "Top_pos": 400,
            "Mid_gene": "m4",
            "Mid_scaf": "Mid_chrB",
            "Mid_pos": 40,
        },
    ]
    bm_rows = _base_bm_rows() + [
        {
            "Bot_gene": "b4",
            "Bot_scaf": "Bot_fused",
            "Bot_pos": 4000,
            "Mid_gene": "m4",
            "Mid_scaf": "Mid_chrB",
            "Mid_pos": 40,
        },
    ]
    palette = _write_tsv(
        tmp_path / "palette.tsv",
        [
            {"gene_group": "TOP_ALG", "color": "#111111"},
            {"gene_group": "MID_ALG", "color": "#222222"},
            {"gene_group": "BOT_ALG", "color": "#333333"},
        ],
    )
    hmm_top = _write_tsv(
        tmp_path / "top.hmm.rbh",
        [{"Top_gene": "t1", "BCnSSimakov2022_scaf": "TOP_ALG"}],
    )
    hmm_mid = _write_tsv(
        tmp_path / "mid.hmm.rbh",
        [
            {"Mid_gene": "m1", "BCnSSimakov2022_scaf": "MID_ALG"},
            {"Mid_gene": "m2", "BCnSSimakov2022_scaf": "MID_ALG"},
        ],
    )
    hmm_bot = _write_tsv(
        tmp_path / "bot.hmm.rbh",
        [
            {"Bot_gene": "b1", "BCnSSimakov2022_scaf": "BOT_ALG"},
            {"Bot_gene": "b3", "BCnSSimakov2022_scaf": "BOT_ALG"},
        ],
    )

    fig, ax, stats = _draw_panel(
        pp,
        tmp_path,
        tm_rows=tm_rows,
        bm_rows=bm_rows,
        color_by="bcns",
        bcns_hmm_top=hmm_top,
        bcns_hmm_mid=hmm_mid,
        bcns_hmm_bot=hmm_bot,
        bcns_palette=palette,
    )
    try:
        edgecolors = [mcolors.to_hex(p.get_edgecolor()) for p in ax.patches]
        legend_text = {t.get_text() for t in ax.get_legend().get_texts()}

        assert stats == {"n_triplets": 4, "n_unmapped": 1}
        assert edgecolors.count("#888888") == 2
        assert edgecolors.count("#111111") == 2
        assert edgecolors.count("#222222") == 2
        assert edgecolors.count("#333333") == 2
        assert "no BCnS call (n=1)" in legend_text
        assert {"BCnS-TOP_ALG", "BCnS-MID_ALG", "BCnS-BOT_ALG"} <= legend_text
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)
