"""Tests for source/odp_ribbon_options.py, the optional ribbon-plot features of odp_rbh_to_ribbon."""
from __future__ import annotations

import pandas as pd
import pytest


@pytest.fixture(scope="module")
def RO(source_dir):
    import odp_ribbon_options
    return odp_ribbon_options


def test_short_chrom_names_drops_shared_text_prefix(RO):
    assert RO.short_chrom_names(["PMA3", "PMA12"]) == {"PMA3": "3", "PMA12": "12"}
    assert RO.short_chrom_names(["Ch_6", "Ch_10", "Ch_1"]) == {"Ch_6": "6", "Ch_10": "10", "Ch_1": "1"}
    assert RO.short_chrom_names(["chr01", "chr10"]) == {"chr01": "1", "chr10": "10"}


def test_short_chrom_names_keeps_accessions_whole(RO):
    names = ["CAIIXF020000008.1", "CAIIXF020000012.1"]
    assert RO.short_chrom_names(names) == {n: n for n in names}
    names = ["NC_047017.1", "NC_047015.1"]
    assert RO.short_chrom_names(names) == {n: n for n in names}


def test_short_chrom_names_single_name(RO):
    assert RO.short_chrom_names(["chr1"]) == {"chr1": "chr1"}


def test_relative_min_size(RO):
    assert RO.relative_min_size([10, 20, 30], 0.2) == pytest.approx(4.0)
    assert RO.relative_min_size([10, 20, 30, 40], 0.5) == pytest.approx(12.5)
    assert RO.relative_min_size([], 0.2) == 0
    assert RO.relative_min_size([10, 20], 0) == 0


def test_shade_runs_light_to_dark(RO):
    import matplotlib.colors as mc
    base = "#336699"
    assert RO.shade(base, 0.5).lower() == base
    light, dark = mc.to_rgb(RO.shade(base, 0.0)), mc.to_rgb(RO.shade(base, 1.0))
    assert sum(light) > sum(mc.to_rgb(base)) > sum(dark)


def test_text_color_for(RO):
    assert RO.text_color_for("#000080") == "white"
    assert RO.text_color_for("#FFFF66") == "black"
    assert RO.text_color_for("#000080", {"K": "black"}, "K") == "black"


def test_is_unassigned(RO):
    assert RO.is_unassigned("#000000", "A1a")
    assert RO.is_unassigned("#123456", "None")
    assert not RO.is_unassigned("#123456", "A1a")


def _pair(a_pos, b_scaf, b_pos, groups=None, colors=None):
    n = len(a_pos)
    return pd.DataFrame({
        "X_gene": ["x%d" % i for i in range(n)], "X_scaf": ["c1"] * n, "X_pos": a_pos,
        "Y_gene": ["y%d" % i for i in range(n)], "Y_scaf": b_scaf, "Y_pos": b_pos,
        "gene_group": groups or ["A"] * n, "color": colors or ["#ff0000"] * n,
    })


def test_collinear_runs_forward_and_reversed(RO):
    # neighbors 0-1-2 stay neighbors forward; 3-4-5 stay neighbors reversed
    df = _pair([1, 2, 3, 4, 5, 6], ["d1"] * 6, [10, 20, 30, 90, 80, 70])
    runs = RO.collinear_runs(df, "X", "Y", min_run=3)
    as_sets = sorted((sorted(r), rev) for r, rev in runs)
    assert as_sets == [([0, 1, 2], False), ([3, 4, 5], True)]


def test_collinear_runs_break_on_crossing(RO):
    # ortholog 1 jumps elsewhere on Y: nothing crossing is allowed inside a run
    df = _pair([1, 2, 3, 4], ["d1"] * 4, [10, 50, 20, 30])
    assert RO.collinear_runs(df, "X", "Y", min_run=3) == []


def test_microsynteny_rows_window(RO):
    df = _pair([1, 2, 3, 4], ["d1"] * 4, [10, 30, 20, 40])
    assert RO.microsynteny_rows(df, "X", "Y", window=3, min_run=3) == {0, 1, 2, 3}
    assert RO.microsynteny_rows(df, "X", "Y", window=1, min_run=3) == set()


def test_parse_follow_and_mask(RO):
    genes, regions = RO.parse_follow(["X:x1", "Y:d1:25-35"])
    assert genes["X"] == {"x1"} and regions == [("Y", "d1", 25.0, 35.0)]
    df = _pair([1, 2, 3], ["d1"] * 3, [10, 20, 30])
    assert RO.follow_mask(df, "X", "Y", genes, regions) == [False, True, True]
    with pytest.raises(ValueError):
        RO.parse_follow(["nonsense"])


def test_order_shades_carry_through_pairs(RO):
    # reference X: three genes of group A on c1; the shades travel to Y and then to Z by gene ID
    xy = _pair([1, 2, 3], ["d1"] * 3, [5, 6, 7])
    yz = pd.DataFrame({"Y_gene": ["y0", "y1", "y2"], "Y_scaf": ["d1"] * 3, "Y_pos": [5, 6, 7],
                       "Z_gene": ["z0", "z1", "z2"], "Z_scaf": ["e1"] * 3, "Z_pos": [9, 8, 7],
                       "gene_group": ["A"] * 3, "color": ["#ff0000"] * 3})
    t = RO.order_shades([xy, yz], ["X", "Y", "Z"], "X")
    assert t[("X", "x0")] == 0.0 and t[("X", "x2")] == 1.0
    assert t[("Z", "z1")] == pytest.approx(0.5)


def test_place_labels_stay_over_their_chromosome(RO):
    homes = RO.label_homes([("A", "c1", 0.1), ("B", "c1", 0.2), ("C", "c2", 0.7)])
    spans = {"c1": (0.0, 0.3), "c2": (0.5, 0.9)}
    widths = {"A": 0.1, "B": 0.1, "C": 0.1}
    placed, ntiers = RO.place_labels(homes, spans, widths, 0.0, 1.0, 0.01)
    for g, cx, bw, tier in placed:
        x0, x1 = spans["c2"] if g == "C" else spans["c1"]
        assert x0 - 1e-9 <= cx - bw / 2 and cx + bw / 2 <= x1 + 1e-9
    assert ntiers == 1


def test_place_labels_stack_when_crowded(RO):
    homes = RO.label_homes([("A", "c1", 0.1), ("B", "c1", 0.15), ("C", "c1", 0.2)])
    placed, ntiers = RO.place_labels(homes, {"c1": (0.0, 0.2)}, {"A": 0.15, "B": 0.15, "C": 0.15}, 0.0, 1.0, 0.01)
    assert ntiers == 3 and sorted(t for *_, t in placed) == [0, 1, 2]


def test_band_patch_and_shuffle(RO):
    from matplotlib.patches import PathPatch
    assert isinstance(RO.band_patch(0, 1, 2, 3, 0, 1, reversed_=True), PathPatch)
    order = RO.shuffled(10, seed=1)
    assert sorted(order) == list(range(10)) and order == RO.shuffled(10, seed=1)


def test_micro_color(RO):
    assert RO.micro_color("gray", "#ff0000") == "#4d4d4d"
    assert RO.micro_color("#123456", None) == "#123456"
    assert RO.micro_color("dark", "#ff0000") != "#ff0000"
