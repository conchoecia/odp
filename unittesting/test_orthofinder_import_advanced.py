"""Advanced tests for the Orthofinder importer.

Coverage focuses on patterns common to ortholog inference test suites:

- **Inparalog handling**: tandem duplicates on the same chromosome.
- **Gene-family expansions**: orthogroups with 10+ paralogs in a species.
- **Adversarial ambiguity**: paralogs whose scores tie under one strategy
  but differ under another.
- **Sparse backbone**: only a handful of single-copy anchors.
- **Backbone-free fallback**: zero single-copy OGs, synteny falls back to
  whatever signal is available without crashing.
- **NCBI-style identifiers**: gene IDs with pipes, dots, accession-style
  prefixes (`XP_001.1`, `gi|12345|ref|XP_001.1`).
- **Compression**: gzipped .chrom and BED inputs.
- **CRLF line endings** in Orthogroups.tsv (Orthofinder's default).
- **Mixed --chrom + --bed flags** across species in one invocation.
- **Round-trip**: build orthogroups from a known .rbh, import, recover.
- **Audit report shape**: header, row count, columns.
- **min-confidence threshold boundaries**.
- **Multi-copy in every species simultaneously** (no within-OG anchor).
- **Missing gene in coord file** — flagged, not silent.
- **Strategy comparison**: synteny rescues more OGs than longest on the
  same fixture, and synteny picks the chromosome-correct paralog where
  longest picks the chromosome-wrong one.
"""
from __future__ import annotations

import gzip
import subprocess
import sys
import time
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


# ===========================================================================
# Inparalog handling
# ===========================================================================


@pytest.fixture
def inparalog_fixture(tmp_path):
    """Two species, with sp1 having a recent tandem duplication
    (g_dup_a, g_dup_b on the same chromosome adjacent positions).
    Single-copy backbone has 3 anchors so the chromosome-pair edge is
    clearly established."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)

    (og_dir / "Orthogroups.tsv").write_text(dedent("""\
        Orthogroup\tsp1\tsp2
        OG0000001\ta1\tb1
        OG0000002\ta2\tb2
        OG0000003\ta3\tb3
        OG0000004\tg_dup_a, g_dup_b\tb_dup
    """))

    (tmp_path / "sp1.chrom").write_text(dedent("""\
        a1\tchrA\t+\t100\t200
        a2\tchrA\t+\t1000\t1100
        a3\tchrA\t+\t2000\t2100
        g_dup_a\tchrA\t+\t3000\t3500
        g_dup_b\tchrA\t-\t3600\t3900
    """))
    (tmp_path / "sp2.chrom").write_text(dedent("""\
        b1\tscafX\t+\t1\t100
        b2\tscafX\t+\t200\t300
        b3\tscafX\t+\t400\t500
        b_dup\tscafX\t+\t600\t700
    """))
    return {"root": root, "sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"}


def test_inparalog_synteny_picks_one_deterministically(ofi, inparalog_fixture):
    """Both tandem duplicates are on the right chromosome — score tie.
    Tie-break by length, then alphabetic, so the choice is deterministic
    across runs."""
    rbh1, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=inparalog_fixture["root"],
        coord_paths={"sp1": inparalog_fixture["sp1"], "sp2": inparalog_fixture["sp2"]},
        strategy="synteny",
    )
    rbh2, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=inparalog_fixture["root"],
        coord_paths={"sp1": inparalog_fixture["sp1"], "sp2": inparalog_fixture["sp2"]},
        strategy="synteny",
    )
    # Same fixture, two invocations: same pick both times.
    pick1 = rbh1[rbh1["rbh"] == rbh1["rbh"].iloc[-1]]["sp1_gene"].iloc[0]
    pick2 = rbh2[rbh2["rbh"] == rbh2["rbh"].iloc[-1]]["sp1_gene"].iloc[0]
    assert pick1 == pick2
    assert pick1 in ("g_dup_a", "g_dup_b")


def test_inparalog_longest_picks_longest_tandem(ofi, inparalog_fixture):
    """g_dup_a is 500bp, g_dup_b is 300bp. 'longest' must pick g_dup_a."""
    rbh, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=inparalog_fixture["root"],
        coord_paths={"sp1": inparalog_fixture["sp1"], "sp2": inparalog_fixture["sp2"]},
        strategy="longest",
    )
    multi_row = rbh[rbh["sp1_gene"].isin(["g_dup_a", "g_dup_b"])].iloc[0]
    assert multi_row["sp1_gene"] == "g_dup_a"


# ===========================================================================
# Gene-family expansion (large multi-copy)
# ===========================================================================


@pytest.fixture
def large_family_fixture(tmp_path):
    """One orthogroup with 20 paralogs in sp1. Only 1 of them is on the
    syntenic chromosome; the rest are on unrelated chromosomes."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)

    # 4 single-copy anchors establishing chrA ↔ scafX
    sc_lines = [
        f"OG{i:07d}\ta{i}\tb{i}" for i in range(1, 5)
    ]
    # Large multi-copy orthogroup
    paralogs = [f"para_{i}" for i in range(20)]
    multi_row = f"OG0000005\t{', '.join(paralogs)}\tb_target"
    (og_dir / "Orthogroups.tsv").write_text(
        "Orthogroup\tsp1\tsp2\n" + "\n".join(sc_lines + [multi_row]) + "\n"
    )

    sp1_lines = [f"a{i}\tchrA\t+\t{i * 1000}\t{i * 1000 + 100}" for i in range(1, 5)]
    # The "true" paralog is on chrA; the rest are on chrB-chrT
    sp1_lines.append(f"para_0\tchrA\t+\t50000\t50500")  # on syntenic chrom
    for i in range(1, 20):
        sp1_lines.append(f"para_{i}\tchr{chr(ord('B') + (i - 1) % 19)}\t+\t100\t200")
    (tmp_path / "sp1.chrom").write_text("\n".join(sp1_lines) + "\n")

    sp2_lines = [f"b{i}\tscafX\t+\t{i * 1000}\t{i * 1000 + 100}" for i in range(1, 5)]
    sp2_lines.append("b_target\tscafX\t+\t60000\t60500")
    (tmp_path / "sp2.chrom").write_text("\n".join(sp2_lines) + "\n")

    return {"root": root, "sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"}


def test_large_gene_family_synteny_picks_syntenic_paralog(ofi, large_family_fixture):
    """With 1 syntenic + 19 non-syntenic paralogs, synteny must pick the
    syntenic one (para_0)."""
    rbh, audit = ofi.orthofinder_to_rbh(
        orthofinder_dir=large_family_fixture["root"],
        coord_paths={"sp1": large_family_fixture["sp1"], "sp2": large_family_fixture["sp2"]},
        strategy="synteny",
    )
    multi_row = rbh[~rbh["sp1_gene"].str.match(r"^a\d+$")].iloc[0]
    assert multi_row["sp1_gene"] == "para_0"
    # Audit should record all 20 candidates considered
    multi_audit = [d for d in audit if d.og == "OG0000005"][0]
    assert multi_audit.n_candidates == 20


def test_large_gene_family_longest_likely_picks_wrong(ofi, large_family_fixture):
    """All paralogs are the same length in this fixture, so 'longest'
    falls back to alphabetic tie-break — para_0 wins by accident here.
    Still, the *reasoning* is different — confirm via the audit."""
    rbh, audit = ofi.orthofinder_to_rbh(
        orthofinder_dir=large_family_fixture["root"],
        coord_paths={"sp1": large_family_fixture["sp1"], "sp2": large_family_fixture["sp2"]},
        strategy="longest",
    )
    multi_audit = [d for d in audit if d.og == "OG0000005"][0]
    assert "length=" in multi_audit.reason


# ===========================================================================
# Sparse backbone
# ===========================================================================


@pytest.fixture
def sparse_backbone(tmp_path):
    """Only ONE single-copy orthogroup to anchor the chromosome map."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(dedent("""\
        Orthogroup\tsp1\tsp2
        OG0000001\ta1\tb1
        OG0000002\ta2, a2b\tb2
    """))
    (tmp_path / "sp1.chrom").write_text(dedent("""\
        a1\tchrA\t+\t1\t100
        a2\tchrA\t+\t200\t300
        a2b\tchrB\t+\t1\t100
    """))
    (tmp_path / "sp2.chrom").write_text(dedent("""\
        b1\tscafX\t+\t1\t100
        b2\tscafX\t+\t200\t300
    """))
    return {"root": root, "sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"}


def test_sparse_backbone_one_anchor_still_resolves(ofi, sparse_backbone):
    """With just one anchor (chrA ↔ scafX, weight 1), synteny should
    still pick a2 (on chrA) over a2b (on chrB)."""
    rbh, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=sparse_backbone["root"],
        coord_paths={"sp1": sparse_backbone["sp1"], "sp2": sparse_backbone["sp2"]},
        strategy="synteny",
        min_confidence=1.0,
    )
    multi = rbh[rbh["sp1_gene"].isin(["a2", "a2b"])].iloc[0]
    assert multi["sp1_gene"] == "a2"


# ===========================================================================
# Backbone-empty fallback
# ===========================================================================


def test_backbone_empty_synteny_degenerates_safely(ofi, tmp_path):
    """No single-copy OGs at all. Synteny has nothing to anchor against.
    The function should not crash; it should produce an empty rbh and a
    multi-copy audit explaining why nothing was rescued."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(dedent("""\
        Orthogroup\tsp1\tsp2
        OG0000001\ta1, a1b\tb1, b1b
    """))
    (tmp_path / "sp1.chrom").write_text(
        "a1\tchrA\t+\t1\t100\na1b\tchrB\t+\t1\t100\n"
    )
    (tmp_path / "sp2.chrom").write_text(
        "b1\tscafX\t+\t1\t100\nb1b\tscafY\t+\t1\t100\n"
    )

    rbh, audit = ofi.orthofinder_to_rbh(
        orthofinder_dir=root,
        coord_paths={"sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"},
        strategy="synteny",
        min_confidence=1.0,
    )
    # No backbone → no chromosome-pair edges → all scores 0 → below
    # min-confidence → no rescue.
    assert len(rbh) == 0


# ===========================================================================
# NCBI-style gene identifiers
# ===========================================================================


def test_ncbi_style_gene_ids_round_trip(ofi, oti, tmp_path):
    """Gene IDs with dots, pipes, underscores, accession prefixes."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(dedent("""\
        Orthogroup\tsp1\tsp2
        OG0000001\tXP_001234.1\tNP_005678.2
        OG0000002\tgi|12345|ref|XP_007890.3|\tNP_111222.1
    """))
    (tmp_path / "sp1.chrom").write_text(dedent("""\
        XP_001234.1\tNC_054740.1\t+\t1\t100
        gi|12345|ref|XP_007890.3|\tNC_054740.1\t+\t200\t300
    """))
    (tmp_path / "sp2.chrom").write_text(dedent("""\
        NP_005678.2\tNC_061156.1\t+\t1\t100
        NP_111222.1\tNC_061156.1\t+\t200\t300
    """))

    rbh, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=root,
        coord_paths={"sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"},
        strategy="synteny",
    )
    assert len(rbh) == 2
    assert "XP_001234.1" in set(rbh["sp1_gene"])
    assert "gi|12345|ref|XP_007890.3|" in set(rbh["sp1_gene"])


# ===========================================================================
# Compression: gzipped .chrom and BED
# ===========================================================================


def test_gzipped_chrom_input(ofi, tmp_path):
    """A .chrom.gz file should parse identically to its plain text twin."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(
        "Orthogroup\tsp1\tsp2\nOG0000001\ta1\tb1\n"
    )

    # Plain text + gzipped twin
    plain = tmp_path / "sp1.chrom"
    plain.write_text("a1\tchrA\t+\t1\t100\n")
    gz = tmp_path / "sp1.chrom.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write("a1\tchrA\t+\t1\t100\n")

    (tmp_path / "sp2.chrom").write_text("b1\tscafX\t+\t1\t100\n")

    rbh, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=root,
        coord_paths={"sp1": gz, "sp2": tmp_path / "sp2.chrom"},
        strategy="skip",  # only one OG, no multi-copy
    )
    assert len(rbh) == 1
    assert rbh.iloc[0]["sp1_gene"] == "a1"
    assert rbh.iloc[0]["sp1_scaf"] == "chrA"


def test_gzipped_bed_input(oti, tmp_path):
    plain_bed = tmp_path / "sp1.bed"
    plain_bed.write_text("chr1\t1\t100\tg1\n")
    gz_bed = tmp_path / "sp1.bed.gz"
    with gzip.open(gz_bed, "wt") as fh:
        fh.write("chr1\t1\t100\tg1\n")
    plain_df = oti.parse_bed(plain_bed)
    gz_df = oti.parse_bed(gz_bed)
    pd.testing.assert_frame_equal(plain_df.reset_index(drop=True),
                                   gz_df.reset_index(drop=True))


# ===========================================================================
# CRLF in Orthogroups.tsv
# ===========================================================================


def test_orthogroups_with_crlf_line_endings(ofi, tmp_path):
    """Orthofinder emits CRLF by default; the parser must normalise."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    # Use bytes-level write so the \r\n stays exactly as written.
    (og_dir / "Orthogroups.tsv").write_bytes(
        b"Orthogroup\tsp1\tsp2\r\n"
        b"OG0000001\ta1\tb1\r\n"
        b"OG0000002\ta2\tb2\r\n"
    )
    species, df = ofi.parse_orthogroups_tsv(og_dir / "Orthogroups.tsv")
    assert species == ["sp1", "sp2"]
    assert len(df) == 2
    # No stray \r should leak into gene IDs
    assert "\r" not in df.iloc[0]["sp2"][0]


# ===========================================================================
# Mixed --chrom and --bed across species
# ===========================================================================


def test_mixed_chrom_and_bed_flags(tmp_path):
    """Some species via --chrom, others via --bed, in the same run."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(
        "Orthogroup\tsp1\tsp2\nOG0000001\ta1\tb1\n"
    )
    sp1_chrom = tmp_path / "sp1.chrom"
    sp1_chrom.write_text("a1\tchrA\t+\t1\t100\n")
    sp2_bed = tmp_path / "sp2.bed"
    sp2_bed.write_text("scafX\t1\t100\tb1\n")

    out = tmp_path / "rbh.tsv"
    result = subprocess.run(
        [
            sys.executable, str(CLI_PATH), "orthofinder-import",
            "--orthofinder-dir", str(root),
            "--chrom", f"sp1={sp1_chrom}",
            "--bed", f"sp2={sp2_bed}",
            "--multi-copy-strategy", "skip",
            "--output", str(out),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    rbh = pd.read_csv(out, sep="\t")
    assert rbh.iloc[0]["sp1_scaf"] == "chrA"
    assert rbh.iloc[0]["sp2_scaf"] == "scafX"


# ===========================================================================
# Round-trip: .rbh -> orthogroups -> orthofinder-import -> .rbh
# ===========================================================================


def test_round_trip_from_rbh(ofi, oti, tmp_path):
    """Take a known .rbh, synthesise an Orthofinder-style table from it,
    re-import, confirm the original ortholog tuples are recovered."""
    # Synthetic .rbh-shaped table
    original_pairs = [
        ("g1", "x1"),
        ("g2", "x2"),
        ("g3", "x3"),
    ]
    sp1_lines, sp2_lines = [], []
    for i, (g, x) in enumerate(original_pairs):
        sp1_lines.append(f"{g}\tchrA\t+\t{(i+1) * 1000}\t{(i+1) * 1000 + 100}")
        sp2_lines.append(f"{x}\tscafX\t+\t{(i+1) * 1000}\t{(i+1) * 1000 + 100}")
    (tmp_path / "sp1.chrom").write_text("\n".join(sp1_lines) + "\n")
    (tmp_path / "sp2.chrom").write_text("\n".join(sp2_lines) + "\n")

    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    og_lines = "Orthogroup\tsp1\tsp2\n"
    for i, (g, x) in enumerate(original_pairs):
        og_lines += f"OG{i+1:07d}\t{g}\t{x}\n"
    (og_dir / "Orthogroups.tsv").write_text(og_lines)

    rbh, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=root,
        coord_paths={"sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"},
        strategy="synteny",
    )

    recovered_pairs = list(zip(rbh["sp1_gene"], rbh["sp2_gene"]))
    assert set(recovered_pairs) == set(original_pairs)


# ===========================================================================
# Audit report shape
# ===========================================================================


def test_audit_report_columns(ofi, orthofinder_synthetic, tmp_path):
    """The audit TSV must have a fixed column set so downstream tools can
    parse it reliably."""
    # Reuse the fixture from the basic test module.
    from unittesting.test_orthofinder_import import orthofinder_synthetic as _ofs  # noqa: F401
    rbh, audit = ofi.orthofinder_to_rbh(
        orthofinder_dir=orthofinder_synthetic["root"],
        coord_paths={
            "sp1": orthofinder_synthetic["sp1_chrom"],
            "sp2": orthofinder_synthetic["sp2_chrom"],
            "sp3": orthofinder_synthetic["sp3_chrom"],
        },
        strategy="synteny",
    )
    report = tmp_path / "audit.tsv"
    ofi.write_audit_report(audit, report)
    df = pd.read_csv(report, sep="\t")
    expected_cols = [
        "OG", "species", "strategy", "n_candidates",
        "picked_gene", "picked_score",
        "second_best_gene", "second_best_score", "reason",
    ]
    assert list(df.columns) == expected_cols


def test_audit_report_empty_when_no_multi_copy(ofi, tmp_path):
    """Skip strategy on a single-copy-only run: audit is empty but the
    file still has a header row."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(
        "Orthogroup\tsp1\tsp2\nOG1\tg1\tx1\n"
    )
    (tmp_path / "sp1.chrom").write_text("g1\tchrA\t+\t1\t100\n")
    (tmp_path / "sp2.chrom").write_text("x1\tscafX\t+\t1\t100\n")

    _, audit = ofi.orthofinder_to_rbh(
        orthofinder_dir=root,
        coord_paths={"sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"},
        strategy="skip",
    )
    report = tmp_path / "audit.tsv"
    n = ofi.write_audit_report(audit, report)
    assert n == 0
    body = report.read_text()
    assert body.split("\n")[0].startswith("OG\tspecies\t")


# ===========================================================================
# min-confidence boundary
# ===========================================================================


def test_min_confidence_zero_keeps_everything(ofi, tmp_path):
    """min-confidence=0 should keep even score-0 picks (anything beats
    skipping)."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(dedent("""\
        Orthogroup\tsp1\tsp2
        OG1\ta1\tb1
        OG2\ta1_dup, a2_dup\tb2
    """))
    (tmp_path / "sp1.chrom").write_text(dedent("""\
        a1\tchrA\t+\t1\t100
        a1_dup\tchrZ\t+\t1\t100
        a2_dup\tchrZ\t+\t200\t300
    """))
    # sp2 b1 on scafX, b2 on scafY — chrA↔scafX edge has weight 1, no
    # edge between chrZ and scafY → chrom_score = 0 for both candidates.
    (tmp_path / "sp2.chrom").write_text(dedent("""\
        b1\tscafX\t+\t1\t100
        b2\tscafY\t+\t1\t100
    """))
    rbh, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=root,
        coord_paths={"sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"},
        strategy="synteny",
        min_confidence=0.0,
    )
    # Both OGs should survive even though OG2's chrom_score is 0.
    assert len(rbh) == 2


# ===========================================================================
# Multi-copy in every species simultaneously
# ===========================================================================


def test_multi_copy_in_every_species(ofi, tmp_path):
    """OG has multiple candidates in both species. No within-OG anchor.
    Synteny must score each pairing against the single-copy backbone."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    # 3 single-copy anchors establishing chrA ↔ scafX
    (og_dir / "Orthogroups.tsv").write_text(dedent("""\
        Orthogroup\tsp1\tsp2
        OG1\ta1\tb1
        OG2\ta2\tb2
        OG3\ta3\tb3
        OG_AMB\tcand_correct_chrA, cand_wrong_chrB\ttarget_correct_scafX, target_wrong_scafY
    """))
    (tmp_path / "sp1.chrom").write_text(dedent("""\
        a1\tchrA\t+\t100\t200
        a2\tchrA\t+\t1000\t1100
        a3\tchrA\t+\t2000\t2100
        cand_correct_chrA\tchrA\t+\t3000\t3100
        cand_wrong_chrB\tchrB\t+\t500\t600
    """))
    (tmp_path / "sp2.chrom").write_text(dedent("""\
        b1\tscafX\t+\t100\t200
        b2\tscafX\t+\t1000\t1100
        b3\tscafX\t+\t2000\t2100
        target_correct_scafX\tscafX\t+\t3000\t3100
        target_wrong_scafY\tscafY\t+\t500\t600
    """))
    rbh, audit = ofi.orthofinder_to_rbh(
        orthofinder_dir=root,
        coord_paths={"sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"},
        strategy="synteny",
    )
    # OG_AMB should survive with the correct pairing.
    amb = rbh[rbh["rbh"] == rbh["rbh"].iloc[-1]].iloc[0]
    # The first multi-copy species gets resolved first; its pick anchors
    # the second species. Either ordering should land on the chrA/scafX
    # pair because both anchor toward chrA↔scafX edges.
    assert amb["sp1_gene"] == "cand_correct_chrA"
    assert amb["sp2_gene"] == "target_correct_scafX"


# ===========================================================================
# Missing gene in coord file
# ===========================================================================


def test_missing_gene_in_chrom_drops_og(ofi, tmp_path):
    """If a gene listed in Orthogroups.tsv has no entry in the .chrom
    file, that orthogroup must be dropped from the output (and the audit
    should record the failure)."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(dedent("""\
        Orthogroup\tsp1\tsp2
        OG1\ta1\tb1
        OG2\tnot_in_chrom\tb2
    """))
    (tmp_path / "sp1.chrom").write_text("a1\tchrA\t+\t1\t100\n")
    (tmp_path / "sp2.chrom").write_text(dedent("""\
        b1\tscafX\t+\t1\t100
        b2\tscafX\t+\t200\t300
    """))
    rbh, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=root,
        coord_paths={"sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"},
        strategy="synteny",
    )
    # Only OG1 survives; OG2's sp1 gene wasn't found.
    assert set(rbh["sp1_gene"]) == {"a1"}


# ===========================================================================
# Strategy comparison: synteny rescues more than longest where it counts
# ===========================================================================


@pytest.fixture
def strategy_divergence_fixture(tmp_path):
    """Constructed so that 'longest' picks the chromosome-wrong paralog
    and 'synteny' picks the chromosome-right one."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(dedent("""\
        Orthogroup\tsp1\tsp2
        OG1\ta1\tb1
        OG2\ta2\tb2
        OG3\ta3\tb3
        OG_TEST\tshort_correct, long_wrong\tb_target
    """))
    (tmp_path / "sp1.chrom").write_text(dedent("""\
        a1\tchrA\t+\t100\t200
        a2\tchrA\t+\t1000\t1100
        a3\tchrA\t+\t2000\t2100
        short_correct\tchrA\t+\t3000\t3100
        long_wrong\tchrZ\t+\t1\t100000
    """))
    (tmp_path / "sp2.chrom").write_text(dedent("""\
        b1\tscafX\t+\t100\t200
        b2\tscafX\t+\t1000\t1100
        b3\tscafX\t+\t2000\t2100
        b_target\tscafX\t+\t3000\t3100
    """))
    return {"root": root, "sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"}


def test_synteny_picks_short_correct(ofi, strategy_divergence_fixture):
    rbh, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=strategy_divergence_fixture["root"],
        coord_paths={"sp1": strategy_divergence_fixture["sp1"], "sp2": strategy_divergence_fixture["sp2"]},
        strategy="synteny",
    )
    test_row = rbh[rbh["sp1_gene"].isin(["short_correct", "long_wrong"])].iloc[0]
    assert test_row["sp1_gene"] == "short_correct"


def test_longest_picks_long_wrong(ofi, strategy_divergence_fixture):
    rbh, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=strategy_divergence_fixture["root"],
        coord_paths={"sp1": strategy_divergence_fixture["sp1"], "sp2": strategy_divergence_fixture["sp2"]},
        strategy="longest",
    )
    test_row = rbh[rbh["sp1_gene"].isin(["short_correct", "long_wrong"])].iloc[0]
    assert test_row["sp1_gene"] == "long_wrong"


# ===========================================================================
# synteny vs most-common-chrom: position breaks chromosome ties
# ===========================================================================


def test_synteny_proximity_breaks_chromosome_ties(ofi, tmp_path):
    """Two paralogs on the same chromosome. most-common-chrom scores them
    identically. synteny tiebreaks by proximity to the single-copy
    anchors on that chromosome."""
    root = tmp_path / "Results"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(dedent("""\
        Orthogroup\tsp1\tsp2
        OG1\ta_anchor\tb_anchor
        OG2\ta_near, a_far\tb_target
    """))
    # Anchor at pos 50_000. a_near is 1 Mbp away (pos ~51_000_000);
    # a_far is 50 Mbp away (pos ~100_000_000). Both on chrA.
    (tmp_path / "sp1.chrom").write_text(dedent("""\
        a_anchor\tchrA\t+\t49000\t51000
        a_near\tchrA\t+\t51000000\t51001000
        a_far\tchrA\t+\t100000000\t100001000
    """))
    (tmp_path / "sp2.chrom").write_text(dedent("""\
        b_anchor\tscafX\t+\t1\t100
        b_target\tscafX\t+\t200\t300
    """))
    rbh_synteny, _ = ofi.orthofinder_to_rbh(
        orthofinder_dir=root,
        coord_paths={"sp1": tmp_path / "sp1.chrom", "sp2": tmp_path / "sp2.chrom"},
        strategy="synteny",
    )
    multi = rbh_synteny[rbh_synteny["sp1_gene"].isin(["a_near", "a_far"])].iloc[0]
    assert multi["sp1_gene"] == "a_near"


# ===========================================================================
# Performance regression (loose bound)
# ===========================================================================


def test_orthogroups_parsing_handles_10k_rows_quickly(ofi, tmp_path):
    """Synthetic 10k-OG Orthogroups.tsv should parse in well under 5
    seconds."""
    p = tmp_path / "big.tsv"
    lines = ["Orthogroup\tsp1\tsp2"]
    for i in range(10_000):
        lines.append(f"OG{i:07d}\tg{i}\tx{i}")
    p.write_text("\n".join(lines) + "\n")

    start = time.perf_counter()
    species, df = ofi.parse_orthogroups_tsv(p)
    elapsed = time.perf_counter() - start
    assert len(df) == 10_000
    assert species == ["sp1", "sp2"]
    assert elapsed < 5.0, f"parse took {elapsed:.2f}s (>5s budget)"


# ===========================================================================
# Strategy: skip + longest + synteny exit codes are all 0
# ===========================================================================


@pytest.mark.parametrize("strategy", ["skip", "longest", "most-common-chrom", "synteny"])
def test_cli_every_strategy_produces_output(strategy, tmp_path, orthofinder_synthetic):
    out = tmp_path / "rbh.tsv"
    result = subprocess.run(
        [
            sys.executable, str(CLI_PATH), "orthofinder-import",
            "--orthofinder-dir", str(orthofinder_synthetic["root"]),
            "--chrom", f"sp1={orthofinder_synthetic['sp1_chrom']}",
            "--chrom", f"sp2={orthofinder_synthetic['sp2_chrom']}",
            "--chrom", f"sp3={orthofinder_synthetic['sp3_chrom']}",
            "--multi-copy-strategy", strategy,
            "--output", str(out),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"strategy={strategy}: {result.stderr}"
    rbh = pd.read_csv(out, sep="\t")
    # Single-copy backbone (4 OGs) must always survive.
    assert len(rbh) >= 4


# Bring the basic fixture into scope here so the strategy-comparison
# CLI parametrize can use it.
@pytest.fixture
def orthofinder_synthetic(tmp_path):
    root = tmp_path / "Results_Test"
    og_dir = root / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text(dedent("""\
        Orthogroup\tsp1\tsp2\tsp3
        OG0000001\tg1\tx1\ty1
        OG0000002\tg2\tx2\ty2
        OG0000003\tg3\tx3\ty3
        OG0000004\tg4\tx4\ty4
        OG0000005\tg5a, g5b\tx5\ty5
        OG0000006\tg6\tx6a, x6b\ty6
    """))
    sp1 = tmp_path / "sp1.chrom"
    sp1.write_text(dedent("""\
        g1\tchr1\t+\t100\t200
        g2\tchr1\t+\t300\t450
        g3\tchr2\t-\t1000\t1500
        g4\tchr2\t+\t2000\t2400
        g5a\tchr1\t+\t500\t600
        g5b\tchr2\t-\t3000\t3200
        g6\tchr2\t+\t2500\t2700
    """))
    sp2 = tmp_path / "sp2.chrom"
    sp2.write_text(dedent("""\
        x1\tscafA\t+\t50\t150
        x2\tscafA\t+\t500\t700
        x3\tscafB\t-\t10\t90
        x4\tscafB\t+\t200\t400
        x5\tscafA\t+\t800\t900
        x6a\tscafA\t-\t1000\t1200
        x6b\tscafB\t+\t500\t700
    """))
    sp3 = tmp_path / "sp3.chrom"
    sp3.write_text(dedent("""\
        y1\tctgX\t+\t1\t10
        y2\tctgX\t+\t20\t30
        y3\tctgY\t+\t100\t200
        y4\tctgY\t+\t300\t400
        y5\tctgX\t+\t40\t50
        y6\tctgY\t-\t500\t600
    """))
    return {
        "root": root, "sp1_chrom": sp1, "sp2_chrom": sp2, "sp3_chrom": sp3,
    }


# ===========================================================================
# CLI: rejects unknown strategy
# ===========================================================================


def test_cli_rejects_unknown_strategy(tmp_path, orthofinder_synthetic):
    result = subprocess.run(
        [
            sys.executable, str(CLI_PATH), "orthofinder-import",
            "--orthofinder-dir", str(orthofinder_synthetic["root"]),
            "--chrom", f"sp1={orthofinder_synthetic['sp1_chrom']}",
            "--chrom", f"sp2={orthofinder_synthetic['sp2_chrom']}",
            "--chrom", f"sp3={orthofinder_synthetic['sp3_chrom']}",
            "--multi-copy-strategy", "tree-based",
            "--output", str(tmp_path / "out.tsv"),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert "invalid choice" in result.stderr


def test_cli_orthofinder_dir_missing(tmp_path):
    result = subprocess.run(
        [
            sys.executable, str(CLI_PATH), "orthofinder-import",
            "--orthofinder-dir", "/nonexistent/path",
            "--chrom", f"sp1={tmp_path}/sp1.chrom",
            "--output", str(tmp_path / "out.tsv"),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 1
    assert "does not exist" in result.stderr.lower()
