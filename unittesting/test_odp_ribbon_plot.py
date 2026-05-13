"""Tests for source/odp_ribbon_plot.py helper functions.

The top-level `ribbon_plot` function is a 400-line plot orchestrator
that needs heavy fixtures (multiple .rbh files, per-species size dicts,
chromosome orderings, etc.) to exercise end-to-end. Tested here are
the importable helpers and module-level setup that don't require
that scaffolding.
"""
from __future__ import annotations

import pytest


@pytest.fixture(scope="module")
def rp(source_dir, fasta_parser_dir):
    import odp_ribbon_plot
    return odp_ribbon_plot


# ---------------------------------------------------------------------------
# Module sets the Agg matplotlib backend at import time
# ---------------------------------------------------------------------------


def test_module_imports_without_display(rp):
    """Module-level `import matplotlib; matplotlib.use("Agg")` runs on
    first import. After that, get_backend() returns the canonical name."""
    import matplotlib
    assert matplotlib.get_backend().lower() == "agg"


# ---------------------------------------------------------------------------
# plot_bezier_lines — draws cubic Bezier curves on a matplotlib axes
# ---------------------------------------------------------------------------


def test_plot_bezier_lines_adds_one_patch_per_line(rp):
    """One PathPatch added per (topx, bottomx) pair."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    n_before = len(ax.patches)
    out = rp.plot_bezier_lines(
        panel=ax,
        topxL=[1.0, 5.0, 10.0],
        bottomxL=[2.0, 4.0, 11.0],
        colors=["#FF0000", "#00FF00", "#0000FF"],
        alpha=[1.0, 0.5, 0.2],
        topy=10.0,
        bottomy=0.0,
    )
    plt.close(fig)
    assert out is ax
    assert len(ax.patches) - n_before == 3


def test_plot_bezier_lines_black_lines_get_low_zorder(rp):
    """Black-colored ribbons are pushed back behind colored ones
    (zorder = -50)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    rp.plot_bezier_lines(
        panel=ax,
        topxL=[1.0, 2.0],
        bottomxL=[3.0, 4.0],
        colors=["#000000", "#FF0000"],
        alpha=[1.0, 1.0],
        topy=10.0,
        bottomy=0.0,
    )
    # First patch is black, must have zorder -50.
    assert ax.patches[-2].get_zorder() == -50
    # Second patch is red with alpha 1.0, must have zorder 1.
    assert ax.patches[-1].get_zorder() == 1
    plt.close(fig)


def test_plot_bezier_lines_low_alpha_pushed_to_back(rp):
    """alpha < 0.5 sends the ribbon to zorder -99."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    rp.plot_bezier_lines(
        panel=ax,
        topxL=[1.0],
        bottomxL=[2.0],
        colors=["#FF0000"],
        alpha=[0.2],
        topy=10.0,
        bottomy=0.0,
    )
    assert ax.patches[-1].get_zorder() == -99
    plt.close(fig)


def test_plot_bezier_lines_empty_inputs(rp):
    """Empty inputs add no patches and return the panel unchanged."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    n_before = len(ax.patches)
    out = rp.plot_bezier_lines(
        panel=ax,
        topxL=[], bottomxL=[], colors=[], alpha=[],
        topy=10.0, bottomy=0.0,
    )
    assert out is ax
    assert len(ax.patches) == n_before
    plt.close(fig)


def test_plot_bezier_lines_color_alpha_attached(rp):
    """The patch picks up the supplied color and alpha verbatim."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    rp.plot_bezier_lines(
        panel=ax,
        topxL=[1.0], bottomxL=[2.0],
        colors=["#0099CC"], alpha=[0.8],
        topy=10.0, bottomy=0.0,
    )
    patch = ax.patches[-1]
    # PathPatch.get_alpha() returns the supplied alpha.
    assert patch.get_alpha() == pytest.approx(0.8)
    plt.close(fig)


# ---------------------------------------------------------------------------
# _quality_check_chromosome_list — internal sanity check
# ---------------------------------------------------------------------------


def test_quality_check_rejects_duplicate_chromosomes(rp):
    """Line 40: duplicate entries in the per-species chromosome list
    raise (catches config typos / pipeline bugs)."""
    sp_to_chr_to_size = {"sp": {"chr1": 100, "chr2": 200}}
    with pytest.raises(IOError, match="duplicate chromosome entries"):
        rp._quality_check_chromosome_list(
            "sp", ["chr1", "chr1"], sp_to_chr_to_size, {}, {"sp": 0},
        )


def test_quality_check_rejects_chrom_missing_from_fasta(rp):
    """Line 44: chromosome in the order list but not in the species'
    genome FASTA is rejected."""
    sp_to_chr_to_size = {"sp": {"chr1": 100}}
    with pytest.raises(IOError, match="not in the fasta file"):
        rp._quality_check_chromosome_list(
            "sp", ["chr1", "chr_phantom"], sp_to_chr_to_size, {}, {"sp": 0},
        )


def test_quality_check_filters_to_gene_order_when_provided(rp):
    """Lines 47-49: if a per-species gene order is supplied, the output
    is filtered to only those chromosomes that appear in it."""
    sp_to_chr_to_size = {"sp": {"chr1": 100, "chr2": 200, "chr3": 300}}
    # Order has only chr1 and chr3 — chr2 should be dropped from the result.
    sp_to_gene_order = {"sp": {"chr1": 0, "chr3": 1}}
    out = rp._quality_check_chromosome_list(
        "sp", ["chr1", "chr2", "chr3"], sp_to_chr_to_size,
        sp_to_gene_order, {"sp": 0},
    )
    assert out == ["chr1", "chr3"]


def test_quality_check_drops_small_chromosomes(rp):
    """Lines 51-53: with no per-species order, chromosomes smaller than
    min size are dropped."""
    sp_to_chr_to_size = {"sp": {"chr1": 1000, "chr_tiny": 50}}
    out = rp._quality_check_chromosome_list(
        "sp", ["chr1", "chr_tiny"], sp_to_chr_to_size,
        sp_to_gene_order={},  # no per-species order
        sp_min_chr_size={"sp": 500},
    )
    assert out == ["chr1"]
