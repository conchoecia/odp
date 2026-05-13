"""Tests for source/odp_plotting_functions.py.

Two public functions: `format_matplotlib` (just configures rcParams)
and `plot_decay` (writes a PDF + a TSV from an ALG-decay datastruct).
Both run under the `Agg` backend that the module forces at import time.
"""
from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def pf(source_dir):
    import odp_plotting_functions
    return odp_plotting_functions


# ---------------------------------------------------------------------------
# Backend is non-interactive
# ---------------------------------------------------------------------------


def test_matplotlib_backend_is_agg(pf):
    """matplotlib.get_backend() returns the canonical name in lowercase
    after `use('Agg')`. We just need to confirm it's non-interactive."""
    import matplotlib
    assert matplotlib.get_backend().lower() == "agg"


# ---------------------------------------------------------------------------
# format_matplotlib — sets rcParams
# ---------------------------------------------------------------------------


def test_format_matplotlib_sets_pdf_fonttype_to_42(pf):
    import matplotlib
    pf.format_matplotlib()
    assert matplotlib.rcParams["pdf.fonttype"] == 42
    assert matplotlib.rcParams["ps.fonttype"] == 42


def test_format_matplotlib_sets_image_composite_false(pf):
    import matplotlib
    pf.format_matplotlib()
    assert matplotlib.rcParams["image.composite_image"] is False


def test_format_matplotlib_idempotent(pf):
    """Calling twice produces the same final state."""
    import matplotlib
    pf.format_matplotlib()
    snap1 = (
        matplotlib.rcParams["pdf.fonttype"],
        matplotlib.rcParams["ps.fonttype"],
        matplotlib.rcParams["image.composite_image"],
    )
    pf.format_matplotlib()
    snap2 = (
        matplotlib.rcParams["pdf.fonttype"],
        matplotlib.rcParams["ps.fonttype"],
        matplotlib.rcParams["image.composite_image"],
    )
    assert snap1 == snap2


# ---------------------------------------------------------------------------
# plot_decay — writes PDF + TSV
# ---------------------------------------------------------------------------


def test_plot_decay_writes_outputs(pf, tmp_path):
    """Synthetic datastruct: 3 ALGs, each with a "conserved" dict and a
    "scattered" dict. The function should produce a PDF and a TSV with
    the right columns."""
    datastruct = {
        "ALG_A1": [{"chr1": 30, "chr2": 5}, {"chr3": 2}],
        "ALG_B1": [{"chr1": 20}, {"chr2": 8}],
        "ALG_C1": [{"chr1": 50}, {}],
    }
    out_pdf = tmp_path / "decay.pdf"
    out_tsv = tmp_path / "decay.tsv"
    # The function uses `pd` as a module-global; need it in scope inside
    # the module. odp_plotting_functions doesn't import pandas at the
    # top; it relies on the snakemake-rule body importing it. Replicate
    # that behavior here.
    import pandas as pd
    pf.pd = pd

    pf.plot_decay(datastruct, str(out_pdf), str(out_tsv))
    assert out_pdf.is_file()
    assert out_pdf.stat().st_size > 0
    assert out_tsv.is_file()
    # TSV should have columns ALG, conserved, scattered, total
    body = out_tsv.read_text()
    header = body.splitlines()[0]
    assert "ALG" in header
    assert "conserved" in header
    assert "scattered" in header
    assert "total" in header
    # 3 ALGs + header
    assert len(body.splitlines()) == 4
