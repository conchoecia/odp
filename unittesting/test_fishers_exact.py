"""Tests for source/fishers_exact.py.

Covers the graph-building / connected-components helpers (``FET_graph``,
``FET_bfs``, ``FET_connected_components``), the chrom-to-CC translator
(``FET_sample_to_chrom_to_CC``), and the Fisher's-exact-test wrapper
(``calculate_FET``) on small synthetic RBH tables.
"""
from __future__ import annotations

import pandas as pd
import pytest


@pytest.fixture(scope="module")
def fet(source_dir):
    import fishers_exact
    return fishers_exact


# ---------------------------------------------------------------------------
# FET_graph + FET_bfs + FET_connected_components
# ---------------------------------------------------------------------------


def test_FET_graph_builds_undirected_adjacency(fet):
    g = fet.FET_graph([("a", "b"), ("b", "c")])
    # Edges are undirected — both endpoints get the neighbour.
    assert g["a"] == {"b"}
    assert g["b"] == {"a", "c"}
    assert g["c"] == {"b"}


def test_FET_graph_extending_existing_graph(fet):
    g1 = fet.FET_graph([("a", "b")])
    g2 = fet.FET_graph([("b", "c")], G=g1)
    # Same dict mutated; "b" now connects to both a and c.
    assert g1 is g2
    assert g2["b"] == {"a", "c"}


def test_FET_graph_self_loop_is_recorded(fet):
    g = fet.FET_graph([("a", "a")])
    assert g["a"] == {"a"}


def test_FET_graph_handles_duplicate_edges(fet):
    g = fet.FET_graph([("a", "b"), ("a", "b"), ("b", "a")])
    # Set-based — duplicates collapse.
    assert g["a"] == {"b"}
    assert g["b"] == {"a"}


def test_FET_bfs_single_component(fet):
    g = fet.FET_graph([("a", "b"), ("b", "c"), ("c", "d")])
    visited = fet.FET_bfs(g, "a")
    assert visited == {"a", "b", "c", "d"}


def test_FET_bfs_isolated_start_node(fet):
    g = fet.FET_graph([])
    g.setdefault("solo", set())
    visited = fet.FET_bfs(g, "solo")
    assert visited == {"solo"}


def test_FET_connected_components_two_components(fet):
    g = fet.FET_graph([
        ("a", "b"), ("b", "c"),  # component {a, b, c}
        ("x", "y"),              # component {x, y}
    ])
    comps = list(fet.FET_connected_components(g))
    by_size = sorted(comps, key=len, reverse=True)
    assert by_size[0] == {"a", "b", "c"}
    assert by_size[1] == {"x", "y"}


def test_FET_connected_components_single_node_components(fet):
    g = fet.FET_graph([("a", "b")])
    # Add an isolated node manually
    g.setdefault("lonely", set())
    comps = list(fet.FET_connected_components(g))
    # 2 components: {a, b} and {lonely}
    found = sorted([sorted(c) for c in comps])
    assert ["a", "b"] in [sorted(c) for c in comps]
    assert {"lonely"} in comps


def test_FET_connected_components_empty_graph(fet):
    g = fet.FET_graph([])
    assert list(fet.FET_connected_components(g)) == []


# ---------------------------------------------------------------------------
# FET_sample_to_chrom_to_CC — assigns chromosomes to connected components
# ---------------------------------------------------------------------------


def test_sample_to_chrom_to_CC_two_components_two_species(fet):
    """Two species, two scaffold pairs, sharing connected components."""
    components = [
        # Component 0: sp1 chr1 ↔ sp2 scafA
        {"sp1|chr1", "sp2|scafA"},
        # Component 1: sp1 chr2 ↔ sp2 scafB
        {"sp1|chr2", "sp2|scafB"},
    ]
    out = fet.FET_sample_to_chrom_to_CC(components, "sp1", "sp2", "|")
    assert out == {
        "sp1": {"chr1": 0, "chr2": 1},
        "sp2": {"scafA": 0, "scafB": 1},
    }


def test_sample_to_chrom_to_CC_rejects_chrom_in_two_components(fet):
    """A chromosome appearing in more than one component is a bug —
    callers expect an error, not silent overwrite."""
    components = [
        {"sp1|chr1", "sp2|scafA"},
        {"sp1|chr1", "sp2|scafB"},  # chr1 again
    ]
    with pytest.raises(ValueError, match="we have seen this chrom before"):
        fet.FET_sample_to_chrom_to_CC(components, "sp1", "sp2", "|")


def test_sample_to_chrom_to_CC_uses_correct_splitter(fet):
    """Splitter is configurable — the function must respect it."""
    components = [{"sp1__chr1", "sp2__scafA"}]
    out = fet.FET_sample_to_chrom_to_CC(components, "sp1", "sp2", "__")
    assert out["sp1"] == {"chr1": 0}
    assert out["sp2"] == {"scafA": 0}


def test_sample_to_chrom_to_CC_empty_components(fet):
    out = fet.FET_sample_to_chrom_to_CC([], "sp1", "sp2", "|")
    assert out == {"sp1": {}, "sp2": {}}


# ---------------------------------------------------------------------------
# calculate_FET — Fisher's exact test on a small RBH dataframe
# ---------------------------------------------------------------------------


def _make_rbh_df(rows):
    """Build a minimal RBH dataframe with the columns calculate_FET expects:
    sp1_gene/sp2_gene + sp1_scaf/sp2_scaf + sp1_breakchrom/sp2_breakchrom.
    Each row in ``rows`` is a tuple
    (sp1_scaf, sp2_scaf, sp1_breakchrom, sp2_breakchrom)."""
    data = []
    for i, (s1, s2, b1, b2) in enumerate(rows):
        data.append({
            "sp1_gene": f"g1_{i}", "sp2_gene": f"g2_{i}",
            "sp1_scaf": s1, "sp2_scaf": s2,
            "sp1_breakchrom": b1, "sp2_breakchrom": b2,
        })
    return pd.DataFrame(data)


def test_calculate_FET_rejects_more_than_two_species(fet):
    df = pd.DataFrame({"sp1_gene": ["a"], "sp2_gene": ["b"], "sp3_gene": ["c"]})
    with pytest.raises(IOError, match="only be performed with two species"):
        fet.calculate_FET(df, {"sp1": [], "sp2": [], "sp3": []})


def test_calculate_FET_rejects_missing_species_in_scafdict(fet):
    df = _make_rbh_df([("chrA", "scafX", "chrA", "scafX")])
    with pytest.raises(IOError, match="not in the scafdict"):
        fet.calculate_FET(df, {"sp1": ["chrA"]})  # sp2 missing


def test_calculate_FET_filters_rows_outside_allowed_scaffolds(fet):
    """Rows on scaffolds absent from scafdict are dropped before testing."""
    df = _make_rbh_df([
        ("chrA", "scafX", "chrA", "scafX"),
        ("chrZ", "scafY", "chrZ", "scafY"),  # chrZ/scafY not allowed
    ])
    scafdict = {"sp1": ["chrA"], "sp2": ["scafX"]}
    result = fet.calculate_FET(df, scafdict)
    # Only one row should survive.
    assert len(result) == 1
    assert result.iloc[0]["sp1_scaf"] == "chrA"


def test_calculate_FET_produces_fet_columns(fet):
    """The function must add `whole_FET` and `break_FET` columns with
    p-values (any non-negative numbers; 1.0 for a contingency that has
    no signal)."""
    rows = [
        ("chrA", "scafX", "chrA", "scafX"),
        ("chrA", "scafX", "chrA", "scafX"),
        ("chrA", "scafX", "chrA", "scafX"),
        ("chrB", "scafY", "chrB", "scafY"),
        ("chrB", "scafY", "chrB", "scafY"),
    ]
    df = _make_rbh_df(rows)
    scafdict = {"sp1": ["chrA", "chrB"], "sp2": ["scafX", "scafY"]}
    out = fet.calculate_FET(df, scafdict)
    assert "whole_FET" in out.columns
    assert "break_FET" in out.columns
    # Every row should have a non-negative FET value (was -1 sentinel,
    # then overwritten by stats.fisher_exact * num-combos).
    assert all(v >= 0 for v in out["whole_FET"])
    assert all(v >= 0 for v in out["break_FET"])


def test_calculate_FET_signals_strong_association(fet):
    """A clean diagonal contingency (chrA always partnered with scafX)
    should produce a small (significant) p-value for the (chrA, scafX)
    combination."""
    rows = [("chrA", "scafX", "chrA", "scafX")] * 20 + \
           [("chrB", "scafY", "chrB", "scafY")] * 20
    df = _make_rbh_df(rows)
    scafdict = {"sp1": ["chrA", "chrB"], "sp2": ["scafX", "scafY"]}
    out = fet.calculate_FET(df, scafdict)
    chrA_scafX_rows = out[(out["sp1_scaf"] == "chrA") & (out["sp2_scaf"] == "scafX")]
    # p-value (Bonferroni-corrected by num-combos = 4) should be small.
    assert chrA_scafX_rows["whole_FET"].iloc[0] < 0.05
