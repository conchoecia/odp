"""Tests for the Orthofinder importer + multi-copy disambiguation.

Covers:
- Orthofinder run-dir discovery (Results_*/, Orthogroups/, missing files).
- Orthogroups.tsv parsing (header, comma+space delimiters, blank cells,
  ragged rows).
- partition_single_vs_multi_copy with and without an explicit single-
  copy list.
- Anchor-graph construction from a synthetic single-copy backbone.
- Each multi-copy strategy: skip, longest, most-common-chrom, synteny.
- min-confidence threshold drops low-score OGs.
- Top-level orthofinder_to_rbh() against a synthetic 3-species run.
- CLI subcommand end-to-end.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from textwrap import dedent

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
CLI_PATH = REPO_ROOT / "bin" / "odp"


@pytest.fixture(scope="module")
def ofi(source_dir):
    import orthofinder_importer
    return orthofinder_importer


@pytest.fixture(scope="module")
def oti(source_dir):
    import ortholog_table_importer
    return ortholog_table_importer


# ---------------------------------------------------------------------------
# Synthetic Orthofinder run
# ---------------------------------------------------------------------------


@pytest.fixture
def orthofinder_synthetic(tmp_path):
    """Build a synthetic 3-species Orthofinder run plus matching .chrom
    files. The single-copy backbone establishes that:

      sp1 chr1 partners with sp2 scafA and sp3 ctgX
      sp1 chr2 partners with sp2 scafB and sp3 ctgY

    Two multi-copy orthogroups are added:
      OG0000005: sp1=[g5a, g5b], sp2=[x5], sp3=[y5]
        g5a is on chr1; g5b is on chr2. sp2.x5 is on scafA → expects chr1.
        So synteny should pick g5a.
      OG0000006: sp1=[g6], sp2=[x6a, x6b], sp3=[y6]
        x6a is on scafA; x6b is on scafB. sp1.g6 is on chr2 → expects
        scafB. Synteny should pick x6b.
    """
    root = tmp_path / "Results_Test"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)

    # Single-copy orthogroups (4 OGs anchoring two chromosome-pair edges).
    # plus the multi-copy ones at the end.
    og_tsv = og_dir / "Orthogroups.tsv"
    og_tsv.write_text(dedent("""\
        Orthogroup\tsp1\tsp2\tsp3
        OG0000001\tg1\tx1\ty1
        OG0000002\tg2\tx2\ty2
        OG0000003\tg3\tx3\ty3
        OG0000004\tg4\tx4\ty4
        OG0000005\tg5a, g5b\tx5\ty5
        OG0000006\tg6\tx6a, x6b\ty6
    """))

    # Single-copy list (explicit).
    (og_dir / "Orthogroups_SingleCopyOrthologues.txt").write_text(
        "OG0000001\nOG0000002\nOG0000003\nOG0000004\n"
    )

    # .chrom files for each species.
    sp1_chrom = tmp_path / "sp1.chrom"
    sp1_chrom.write_text(dedent("""\
        g1\tchr1\t+\t100\t200
        g2\tchr1\t+\t300\t450
        g3\tchr2\t-\t1000\t1500
        g4\tchr2\t+\t2000\t2400
        g5a\tchr1\t+\t500\t600
        g5b\tchr2\t-\t3000\t3200
        g6\tchr2\t+\t2500\t2700
    """))
    sp2_chrom = tmp_path / "sp2.chrom"
    sp2_chrom.write_text(dedent("""\
        x1\tscafA\t+\t50\t150
        x2\tscafA\t+\t500\t700
        x3\tscafB\t-\t10\t90
        x4\tscafB\t+\t200\t400
        x5\tscafA\t+\t800\t900
        x6a\tscafA\t-\t1000\t1200
        x6b\tscafB\t+\t500\t700
    """))
    sp3_chrom = tmp_path / "sp3.chrom"
    sp3_chrom.write_text(dedent("""\
        y1\tctgX\t+\t1\t10
        y2\tctgX\t+\t20\t30
        y3\tctgY\t+\t100\t200
        y4\tctgY\t+\t300\t400
        y5\tctgX\t+\t40\t50
        y6\tctgY\t-\t500\t600
    """))

    return {
        "root": root,
        "og_tsv": og_tsv,
        "sp1_chrom": sp1_chrom,
        "sp2_chrom": sp2_chrom,
        "sp3_chrom": sp3_chrom,
    }


# ---------------------------------------------------------------------------
# discover_orthofinder_paths
# ---------------------------------------------------------------------------


def test_discover_paths_not_a_directory_errors(ofi, tmp_path):
    """Line 75: non-existent path / non-directory is rejected up front."""
    f = tmp_path / "not_a_dir.txt"
    f.write_text("hi")
    with pytest.raises(FileNotFoundError, match="does not exist"):
        ofi.discover_orthofinder_paths(f)


def test_discover_paths_finds_optional_dirs(ofi, tmp_path):
    """Lines 99, 103, 107: discover_orthofinder_paths sets the optional
    output directory attributes when those dirs exist."""
    root = tmp_path / "Results_X"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text("Orthogroup\tsp1\nOG1\tg1\n")
    (root / "Single_Copy_Orthologue_Sequences").mkdir()
    (root / "Phylogenetic_Hierarchical_Orthogroups").mkdir()
    (root / "Resolved_Gene_Trees").mkdir()
    paths = ofi.discover_orthofinder_paths(root)
    assert paths.single_copy_sequences_dir is not None
    assert paths.hierarchical_orthogroups_dir is not None
    assert paths.resolved_gene_trees_dir is not None


def test_parse_orthogroups_gzipped(ofi, tmp_path):
    """Line 131: .gz extension triggers gzip opener."""
    import gzip
    p = tmp_path / "Orthogroups.tsv.gz"
    body = "Orthogroup\tsp1\tsp2\nOG1\tg1\tx1\nOG2\tg2\tx2\n"
    with gzip.open(p, "wt") as fh:
        fh.write(body)
    species, df = ofi.parse_orthogroups_tsv(p)
    assert species == ["sp1", "sp2"]
    assert len(df) == 2


def test_parse_orthogroups_skips_blank_lines(ofi, tmp_path):
    """Line 139-140: blank lines inside the file are skipped."""
    p = tmp_path / "Orthogroups.tsv"
    p.write_text("Orthogroup\tsp1\tsp2\n\nOG1\tg1\tx1\n\nOG2\tg2\tx2\n")
    species, df = ofi.parse_orthogroups_tsv(p)
    assert len(df) == 2


def test_parse_orthogroups_header_only_raises(ofi, tmp_path):
    """Line 142-143: only-header file has no data rows; reject."""
    p = tmp_path / "Orthogroups.tsv"
    p.write_text("Orthogroup\tsp1\tsp2\n")
    with pytest.raises(ValueError, match="no data rows"):
        ofi.parse_orthogroups_tsv(p)


def test_parse_orthogroups_warns_on_non_orthogroup_header_col0(ofi, capsys, tmp_path):
    """Lines 146-150: warning when first header cell isn't 'Orthogroup'."""
    p = tmp_path / "Orthogroups.tsv"
    p.write_text("OG_ID\tsp1\tsp2\nOG1\tg1\tx1\n")
    ofi.parse_orthogroups_tsv(p)
    err = capsys.readouterr().err
    assert "expected 'Orthogroup'" in err


def test_parse_orthogroups_no_species_columns_errors(ofi, tmp_path):
    """Lines 152-155: header with only the OG column and no species."""
    p = tmp_path / "Orthogroups.tsv"
    p.write_text("Orthogroup\nOG1\n")
    with pytest.raises(ValueError, match="no species columns"):
        ofi.parse_orthogroups_tsv(p)


def test_parse_orthogroups_pads_ragged_row(ofi, tmp_path):
    """Line 159-161: row with fewer cells than the header is padded with
    empty strings (treated as no genes for the missing species)."""
    p = tmp_path / "Orthogroups.tsv"
    p.write_text("Orthogroup\tsp1\tsp2\tsp3\nOG1\tg1\tx1\n")  # row missing sp3
    species, df = ofi.parse_orthogroups_tsv(p)
    assert species == ["sp1", "sp2", "sp3"]
    assert df.iloc[0]["sp3"] == []


def test_disambiguate_unknown_strategy_raises(ofi, oti, orthofinder_synthetic):
    """Line 378-382: unknown strategy is rejected with the four-option list."""
    fx = orthofinder_synthetic
    species_names, og_df = ofi.parse_orthogroups_tsv(fx["og_tsv"])
    _, multi = ofi.partition_single_vs_multi_copy(og_df, species_names)
    coords = {
        sp: oti.parse_coordinates(fx[f"{sp}_chrom"]) for sp in species_names
    }
    with pytest.raises(ValueError, match="unknown multi-copy strategy"):
        ofi.disambiguate_multicopy(
            multi, species_names, coords, None,
            strategy="not-a-real-strategy",
        )


def test_disambiguate_synteny_without_anchor_graph_raises(
    ofi, oti, orthofinder_synthetic
):
    """Line 384-388: synteny + most-common-chrom need the anchor_graph
    arg; passing None is rejected."""
    fx = orthofinder_synthetic
    species_names, og_df = ofi.parse_orthogroups_tsv(fx["og_tsv"])
    _, multi = ofi.partition_single_vs_multi_copy(og_df, species_names)
    coords = {
        sp: oti.parse_coordinates(fx[f"{sp}_chrom"]) for sp in species_names
    }
    with pytest.raises(ValueError, match="needs an anchor graph"):
        ofi.disambiguate_multicopy(
            multi, species_names, coords, None,
            strategy="synteny",
        )


def test_emit_rbh_skips_rows_with_empty_gene_list(ofi, oti, tmp_path, capsys):
    """Line 571-573: a row with an empty list for any species is dropped
    from the .rbh output."""
    import pandas as pd
    bed = tmp_path / "sp.bed"
    bed.write_text(
        "chr1\t100\t200\tg1\n"
        "chr1\t300\t400\tg2\n"
    )
    coords = {"sp1": oti.parse_coordinates(bed), "sp2": oti.parse_coordinates(bed)}
    df = pd.DataFrame({
        "OG": ["OG1", "OG2"],
        "sp1": [["g1"], []],     # OG2 has empty list for sp1 — drop
        "sp2": [["g1"], ["g2"]],
    })
    rbh = ofi.emit_rbh_from_orthogroups(df, ["sp1", "sp2"], coords)
    assert len(rbh) == 1
    err = capsys.readouterr().err
    assert "dropped 1 OG" in err


def test_disambiguate_skips_og_with_zero_genes_in_species(
    ofi, oti, orthofinder_synthetic, tmp_path
):
    """Lines 414-424: an OG whose species column is empty triggers an
    audit entry with reason 'no candidates in this species'."""
    import pandas as pd
    fx = orthofinder_synthetic
    species_names = ["sp1", "sp2", "sp3"]
    coords = {
        sp: oti.parse_coordinates(fx[f"{sp}_chrom"]) for sp in species_names
    }
    # An OG with one species empty; the disambiguator should record a
    # "no candidates" audit row and drop the OG.
    multi = pd.DataFrame([{
        "OG": "OG9",
        "sp1": ["g1", "g2"],  # multi-copy in sp1
        "sp2": [],            # empty — triggers "no candidates"
        "sp3": ["y1"],
    }])
    resolved, audit = ofi.disambiguate_multicopy(
        multi, species_names, coords, anchor_graph=None,
        strategy="longest",
    )
    assert len(resolved) == 0
    reasons = [d.reason for d in audit]
    assert any("no candidates in this species" in r for r in reasons)


def test_disambiguate_when_no_candidate_in_coord_file(
    ofi, oti, orthofinder_synthetic, tmp_path
):
    """Lines 433-444: when a multi-copy species' candidate genes are
    none in the coord lookup, the OG is dropped with a specific audit
    reason."""
    import pandas as pd
    fx = orthofinder_synthetic
    species_names = ["sp1", "sp2", "sp3"]
    coords = {
        sp: oti.parse_coordinates(fx[f"{sp}_chrom"]) for sp in species_names
    }
    multi = pd.DataFrame([{
        "OG": "OGZ",
        "sp1": ["ghost1", "ghost2"],  # neither in sp1.chrom
        "sp2": ["x5"],
        "sp3": ["y5"],
    }])
    resolved, audit = ofi.disambiguate_multicopy(
        multi, species_names, coords, anchor_graph=None,
        strategy="longest",
    )
    assert len(resolved) == 0
    reasons = [d.reason for d in audit]
    assert any("no candidate gene found in coord file" in r for r in reasons)


def test_discover_paths_with_results_root(ofi, orthofinder_synthetic):
    paths = ofi.discover_orthofinder_paths(orthofinder_synthetic["root"])
    assert paths.orthogroups_tsv.name == "Orthogroups.tsv"
    assert paths.single_copy_list is not None
    assert paths.single_copy_list.name == "Orthogroups_SingleCopyOrthologues.txt"


def test_discover_paths_with_orthogroups_subdir(ofi, orthofinder_synthetic):
    paths = ofi.discover_orthofinder_paths(orthofinder_synthetic["og_tsv"].parent)
    assert paths.orthogroups_tsv.name == "Orthogroups.tsv"


def test_discover_paths_missing_orthogroups_errors(ofi, tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError, match="Orthogroups.tsv"):
        ofi.discover_orthofinder_paths(empty)


# ---------------------------------------------------------------------------
# parse_orthogroups_tsv
# ---------------------------------------------------------------------------


def test_parse_orthogroups_basic(ofi, orthofinder_synthetic):
    species, df = ofi.parse_orthogroups_tsv(orthofinder_synthetic["og_tsv"])
    assert species == ["sp1", "sp2", "sp3"]
    assert len(df) == 6
    # OG0000005 multi-copy in sp1
    og5 = df[df["OG"] == "OG0000005"].iloc[0]
    assert og5["sp1"] == ["g5a", "g5b"]
    assert og5["sp2"] == ["x5"]


def test_parse_orthogroups_handles_comma_only(ofi, tmp_path):
    p = tmp_path / "Orthogroups.tsv"
    p.write_text("Orthogroup\tsp1\tsp2\nOG1\tg1,g1b\tx1\n")
    species, df = ofi.parse_orthogroups_tsv(p)
    assert df.iloc[0]["sp1"] == ["g1", "g1b"]


def test_parse_orthogroups_handles_empty_cells(ofi, tmp_path):
    p = tmp_path / "Orthogroups.tsv"
    p.write_text("Orthogroup\tsp1\tsp2\nOG1\t\tx1\n")
    species, df = ofi.parse_orthogroups_tsv(p)
    assert df.iloc[0]["sp1"] == []


# ---------------------------------------------------------------------------
# partition_single_vs_multi_copy
# ---------------------------------------------------------------------------


def test_partition_with_explicit_list(ofi, orthofinder_synthetic):
    species, df = ofi.parse_orthogroups_tsv(orthofinder_synthetic["og_tsv"])
    sc_list = ofi.load_single_copy_list(
        orthofinder_synthetic["root"] / "Orthogroups" /
        "Orthogroups_SingleCopyOrthologues.txt"
    )
    single, multi = ofi.partition_single_vs_multi_copy(df, species, sc_list)
    assert set(single["OG"]) == {"OG0000001", "OG0000002", "OG0000003", "OG0000004"}
    assert set(multi["OG"]) == {"OG0000005", "OG0000006"}


def test_partition_without_explicit_list(ofi, orthofinder_synthetic):
    species, df = ofi.parse_orthogroups_tsv(orthofinder_synthetic["og_tsv"])
    single, multi = ofi.partition_single_vs_multi_copy(df, species, None)
    # Auto-partition by cell content should match the explicit list.
    assert set(single["OG"]) == {"OG0000001", "OG0000002", "OG0000003", "OG0000004"}


# ---------------------------------------------------------------------------
# build_anchor_graph
# ---------------------------------------------------------------------------


def test_build_anchor_graph(ofi, oti, orthofinder_synthetic):
    species, df = ofi.parse_orthogroups_tsv(orthofinder_synthetic["og_tsv"])
    single, _ = ofi.partition_single_vs_multi_copy(df, species, None)
    coords = {
        "sp1": oti.parse_chrom(orthofinder_synthetic["sp1_chrom"]),
        "sp2": oti.parse_chrom(orthofinder_synthetic["sp2_chrom"]),
        "sp3": oti.parse_chrom(orthofinder_synthetic["sp3_chrom"]),
    }
    g = ofi.build_anchor_graph(single, species, coords)
    # sp1 chr1 ↔ sp2 scafA should appear twice (OG1 + OG2)
    assert g.edge[("sp1", "sp2")].get(("chr1", "scafA"), 0) == 2
    # sp1 chr2 ↔ sp2 scafB twice (OG3 + OG4)
    assert g.edge[("sp1", "sp2")].get(("chr2", "scafB"), 0) == 2
    # sp1 chr1 ↔ sp3 ctgX twice
    assert g.edge[("sp1", "sp3")].get(("chr1", "ctgX"), 0) == 2
    # Reverse direction populated
    assert g.edge[("sp2", "sp1")].get(("scafA", "chr1"), 0) == 2
    # Positions recorded
    assert len(g.positions["sp1"]["chr1"]) == 2


# ---------------------------------------------------------------------------
# Disambiguation strategies
# ---------------------------------------------------------------------------


def test_strategy_skip_drops_all_multi(ofi, oti, orthofinder_synthetic):
    species, df = ofi.parse_orthogroups_tsv(orthofinder_synthetic["og_tsv"])
    single, multi = ofi.partition_single_vs_multi_copy(df, species, None)
    coords = {
        sp: oti.parse_chrom(orthofinder_synthetic[f"{sp}_chrom"])
        for sp in species
    }
    resolved, audit = ofi.disambiguate_multicopy(
        multi, species, coords, None, strategy="skip",
    )
    assert len(resolved) == 0
    assert len(audit) == 0


def test_strategy_longest_picks_longest_gene(ofi, oti, orthofinder_synthetic):
    species, df = ofi.parse_orthogroups_tsv(orthofinder_synthetic["og_tsv"])
    single, multi = ofi.partition_single_vs_multi_copy(df, species, None)
    coords = {
        sp: oti.parse_chrom(orthofinder_synthetic[f"{sp}_chrom"])
        for sp in species
    }
    resolved, audit = ofi.disambiguate_multicopy(
        multi, species, coords, None, strategy="longest",
    )
    # OG0000005: sp1=[g5a (len 100), g5b (len 200)]. Longest = g5b.
    og5 = resolved[resolved["OG"] == "OG0000005"].iloc[0]
    assert og5["sp1"] == ["g5b"]
    # OG0000006: sp2=[x6a (len 200), x6b (len 200)] — tie. Tie-break by
    # alphabetic gene id ascending → x6a wins on tie.
    og6 = resolved[resolved["OG"] == "OG0000006"].iloc[0]
    assert og6["sp2"][0] in ("x6a", "x6b")


def test_strategy_synteny_picks_chromosome_anchored_paralog(
    ofi, oti, orthofinder_synthetic
):
    species, df = ofi.parse_orthogroups_tsv(orthofinder_synthetic["og_tsv"])
    single, multi = ofi.partition_single_vs_multi_copy(df, species, None)
    coords = {
        sp: oti.parse_chrom(orthofinder_synthetic[f"{sp}_chrom"])
        for sp in species
    }
    g = ofi.build_anchor_graph(single, species, coords)
    resolved, audit = ofi.disambiguate_multicopy(
        multi, species, coords, g, strategy="synteny", min_confidence=1.0,
    )
    # OG0000005: sp1's other species are sp2.x5 (scafA) and sp3.y5 (ctgX).
    # scafA pairs with chr1; ctgX pairs with chr1. So sp1 should pick g5a (chr1).
    og5 = resolved[resolved["OG"] == "OG0000005"].iloc[0]
    assert og5["sp1"] == ["g5a"]
    # OG0000006: sp2 candidates are x6a (scafA) vs x6b (scafB). sp1.g6
    # is on chr2 and sp3.y6 is on ctgY. chr2 pairs with scafB, ctgY pairs
    # with scafB. So sp2 should pick x6b.
    og6 = resolved[resolved["OG"] == "OG0000006"].iloc[0]
    assert og6["sp2"] == ["x6b"]


def test_strategy_most_common_chrom_same_picks_synteny(
    ofi, oti, orthofinder_synthetic
):
    species, df = ofi.parse_orthogroups_tsv(orthofinder_synthetic["og_tsv"])
    single, multi = ofi.partition_single_vs_multi_copy(df, species, None)
    coords = {
        sp: oti.parse_chrom(orthofinder_synthetic[f"{sp}_chrom"])
        for sp in species
    }
    g = ofi.build_anchor_graph(single, species, coords)
    resolved, _audit = ofi.disambiguate_multicopy(
        multi, species, coords, g,
        strategy="most-common-chrom", min_confidence=1.0,
    )
    # In this fixture the chromosome vote alone is enough to resolve
    # both OGs, identical to the synteny strategy.
    og5 = resolved[resolved["OG"] == "OG0000005"].iloc[0]
    og6 = resolved[resolved["OG"] == "OG0000006"].iloc[0]
    assert og5["sp1"] == ["g5a"]
    assert og6["sp2"] == ["x6b"]


def test_min_confidence_drops_low_score(ofi, oti, orthofinder_synthetic):
    """A very high threshold should drop both multi-copy OGs."""
    species, df = ofi.parse_orthogroups_tsv(orthofinder_synthetic["og_tsv"])
    single, multi = ofi.partition_single_vs_multi_copy(df, species, None)
    coords = {
        sp: oti.parse_chrom(orthofinder_synthetic[f"{sp}_chrom"])
        for sp in species
    }
    g = ofi.build_anchor_graph(single, species, coords)
    resolved, _ = ofi.disambiguate_multicopy(
        multi, species, coords, g,
        strategy="synteny", min_confidence=1000.0,
    )
    assert len(resolved) == 0


# ---------------------------------------------------------------------------
# Top-level orthofinder_to_rbh
# ---------------------------------------------------------------------------


def test_orthofinder_to_rbh_synteny_e2e(ofi, orthofinder_synthetic):
    rbh, audit = ofi.orthofinder_to_rbh(
        orthofinder_dir=orthofinder_synthetic["root"],
        coord_paths={
            "sp1": orthofinder_synthetic["sp1_chrom"],
            "sp2": orthofinder_synthetic["sp2_chrom"],
            "sp3": orthofinder_synthetic["sp3_chrom"],
        },
        strategy="synteny",
        min_confidence=1.0,
    )
    # 4 single-copy + 2 rescued = 6 total
    assert len(rbh) == 6
    # Confirm column shape
    expected_cols = [
        "rbh", "gene_group",
        "sp1_gene", "sp1_scaf", "sp1_pos",
        "sp2_gene", "sp2_scaf", "sp2_pos",
        "sp3_gene", "sp3_scaf", "sp3_pos",
    ]
    assert list(rbh.columns) == expected_cols
    # Audit covers both multi-copy decisions (sp1 in OG5, sp2 in OG6)
    assert len(audit) == 2


def test_orthofinder_to_rbh_skip_strategy(ofi, orthofinder_synthetic):
    rbh, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=orthofinder_synthetic["root"],
        coord_paths={
            "sp1": orthofinder_synthetic["sp1_chrom"],
            "sp2": orthofinder_synthetic["sp2_chrom"],
            "sp3": orthofinder_synthetic["sp3_chrom"],
        },
        strategy="skip",
    )
    # Only the 4 single-copy OGs.
    assert len(rbh) == 4


def test_orthofinder_to_rbh_missing_coord_for_species_errors(ofi, orthofinder_synthetic):
    with pytest.raises(ValueError, match="Coord file missing"):
        ofi.orthofinder_to_rbh(
            orthofinder_dir=orthofinder_synthetic["root"],
            coord_paths={
                "sp1": orthofinder_synthetic["sp1_chrom"],
                # sp2, sp3 missing
            },
        )


# ---------------------------------------------------------------------------
# CLI end-to-end
# ---------------------------------------------------------------------------


def test_cli_orthofinder_import_synteny(tmp_path, orthofinder_synthetic):
    out = tmp_path / "rbh.tsv"
    report = tmp_path / "audit.tsv"
    result = subprocess.run(
        [
            sys.executable, str(CLI_PATH),
            "orthofinder-import",
            "--orthofinder-dir", str(orthofinder_synthetic["root"]),
            "--chrom", f"sp1={orthofinder_synthetic['sp1_chrom']}",
            "--chrom", f"sp2={orthofinder_synthetic['sp2_chrom']}",
            "--chrom", f"sp3={orthofinder_synthetic['sp3_chrom']}",
            "--multi-copy-strategy", "synteny",
            "--output", str(out),
            "--report", str(report),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    body = out.read_text()
    assert "rbh1" in body
    # The synteny strategy should rescue g5a (not g5b) and x6b (not x6a).
    rbh_df = pd.read_csv(out, sep="\t")
    assert "g5a" in set(rbh_df["sp1_gene"])
    assert "g5b" not in set(rbh_df["sp1_gene"])
    assert "x6b" in set(rbh_df["sp2_gene"])
    assert "x6a" not in set(rbh_df["sp2_gene"])
    # Audit report exists and has the two disambiguation decisions
    audit_df = pd.read_csv(report, sep="\t")
    assert len(audit_df) == 2


def test_cli_orthofinder_import_requires_coords(tmp_path, orthofinder_synthetic):
    out = tmp_path / "rbh.tsv"
    result = subprocess.run(
        [
            sys.executable, str(CLI_PATH),
            "orthofinder-import",
            "--orthofinder-dir", str(orthofinder_synthetic["root"]),
            "--output", str(out),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert "at least one --chrom or --bed" in result.stderr


def test_cli_help_lists_orthofinder_import():
    result = subprocess.run(
        [sys.executable, str(CLI_PATH), "--help"],
        capture_output=True, text=True,
    )
    assert "orthofinder-import" in result.stdout
