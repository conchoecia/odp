"""Tests for source/odp_functions.py.

Focuses on the small utility functions that are exercised throughout
the odp pipeline: resource estimators, file-existence guards, chrom
file legality + parsing, expand-style helpers used to build
snakemake input lists, and small text/IO helpers.

Heavier statistical functions (`calc_D_for_y_and_x`,
`reciprocal_best_hits_*`) need much larger fixtures and are left for a
later pass.
"""
from __future__ import annotations

import gzip
from pathlib import Path
from textwrap import dedent

import pytest


@pytest.fixture(scope="module")
def of(source_dir, fasta_parser_dir):
    import odp_functions
    return odp_functions


# ---------------------------------------------------------------------------
# Resource estimator dicts (`*_get_mem_mb`)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("attempt,expected", [
    (1, 2000), (2, 4000), (3, 8000), (4, 16000),
    (5, 32000), (6, 64000), (7, 128000), (8, 256000),
])
def test_hmmsearch_get_mem_mb_scales_with_attempt(of, attempt, expected):
    assert of.hmmsearch_get_mem_mb(None, attempt) == expected


def test_hmmsearch_get_mem_mb_attempt_zero_raises(of):
    """attempt=0 (or anything missing from the dict) raises KeyError —
    callers should never pass a 0-th attempt."""
    with pytest.raises(KeyError):
        of.hmmsearch_get_mem_mb(None, 0)


def test_gzhmm_get_mem_mb_falls_back_to_max_for_unknown_attempt(of):
    """`gzhmm_get_mem_mb` is the one variant that uses .get() with a
    safe maximum fallback, unlike its siblings."""
    mx = of.gzhmm_get_mem_mb(None, 1000)
    assert mx == 640000  # equals attemptdict[max(attemptdict)] == [7]


def test_filtprot_get_mem_mb_default_for_unknown(of):
    mx = of.filtprot_get_mem_mb(None, 1000)
    assert mx == 512000  # max value in the dict


def test_filtprot_get_mem_mb_attempt_1(of):
    assert of.filtprot_get_mem_mb(None, 1) == 2000


def test_tmp_unzip_get_mem_mb_basic(of):
    assert of.tmp_unzip_get_mem_mb(None, 1) == 2000
    assert of.tmp_unzip_get_mem_mb(None, 4) == 16000


def test_hmm_against_prots_get_mem_mb_basic(of):
    assert of.hmm_against_prots_get_mem_mb(None, 1) == 5000
    assert of.hmm_against_prots_get_mem_mb(None, 4) == 256000


def test_filthmm_get_mem_mb_basic(of):
    assert of.filthmm_get_mem_mb(None, 1) == 1000
    assert of.filthmm_get_mem_mb(None, 7) == 640000


# ---------------------------------------------------------------------------
# check_file_exists
# ---------------------------------------------------------------------------


def test_check_file_exists_returns_true_for_real_file(of, tmp_path):
    p = tmp_path / "ok.txt"
    p.write_text("hi")
    assert of.check_file_exists(str(p)) is True


def test_check_file_exists_raises_on_missing_file(of, tmp_path):
    with pytest.raises(IOError, match="This file does not exist"):
        of.check_file_exists(str(tmp_path / "no_such_file"))


def test_check_file_exists_rejects_directories(of, tmp_path):
    """A directory is not a regular file — should raise."""
    sub = tmp_path / "sub"
    sub.mkdir()
    with pytest.raises(IOError):
        of.check_file_exists(str(sub))


# ---------------------------------------------------------------------------
# chrom_file_is_legal
# ---------------------------------------------------------------------------


def test_chrom_file_is_legal_well_formed(of, tmp_path):
    p = tmp_path / "good.chrom"
    p.write_text(dedent("""\
        prot1\tchrA\t+\t100\t200
        prot2\tchrA\t-\t300\t400
        prot3\tchrB\t.\t500\t600
    """))
    assert of.chrom_file_is_legal(str(p)) is True


def test_chrom_file_is_legal_rejects_header_line(of, tmp_path, capsys):
    p = tmp_path / "header.chrom"
    p.write_text(
        "pid\tscaf\tstrand\tstart\tstop\n"
        "prot1\tchrA\t+\t100\t200\n"
    )
    assert of.chrom_file_is_legal(str(p)) is False


def test_chrom_file_is_legal_rejects_bad_strand(of, tmp_path):
    p = tmp_path / "badstrand.chrom"
    p.write_text("prot1\tchrA\tFORWARD\t100\t200\n")
    assert of.chrom_file_is_legal(str(p)) is False


def test_chrom_file_is_legal_rejects_non_integer_start(of, tmp_path):
    p = tmp_path / "badint.chrom"
    p.write_text("prot1\tchrA\t+\tNOT_AN_INT\t200\n")
    assert of.chrom_file_is_legal(str(p)) is False


def test_chrom_file_is_legal_rejects_non_integer_stop(of, tmp_path):
    p = tmp_path / "badint.chrom"
    p.write_text("prot1\tchrA\t+\t100\tNOT_AN_INT\n")
    assert of.chrom_file_is_legal(str(p)) is False


def test_chrom_file_is_legal_handles_gzip(of, tmp_path):
    p = tmp_path / "good.chrom.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("prot1\tchrA\t+\t100\t200\n")
    assert of.chrom_file_is_legal(str(p)) is True


def test_chrom_file_is_legal_raises_on_missing_file(of, tmp_path):
    with pytest.raises(IOError, match="does not exist"):
        of.chrom_file_is_legal(str(tmp_path / "no.chrom"))


# ---------------------------------------------------------------------------
# flatten
# ---------------------------------------------------------------------------


def test_flatten_basic(of):
    out = of.flatten([[1, 2], [3, 4]])
    assert set(out) == {1, 2, 3, 4}


def test_flatten_deduplicates(of):
    """flatten() returns unique items via set()."""
    out = of.flatten([[1, 2], [2, 3], [1]])
    assert sorted(out) == [1, 2, 3]


def test_flatten_empty(of):
    assert of.flatten([]) == []
    assert of.flatten([[]]) == []


def test_flatten_strings(of):
    """flatten works on strings of lists too (set of strings)."""
    out = of.flatten([["a", "b"], ["b", "c"]])
    assert set(out) == {"a", "b", "c"}


# ---------------------------------------------------------------------------
# expand_avoid_matching_x_and_y
# ---------------------------------------------------------------------------


def test_expand_avoid_matching_x_and_y_basic(of):
    out = of.expand_avoid_matching_x_and_y(
        "out/{}-{}.txt",
        xsamples=["A", "B"],
        ysamples=["A", "B"],
    )
    # All combos minus the matched ones (A-A, B-B).
    assert sorted(out) == ["out/A-B.txt", "out/B-A.txt"]


def test_expand_avoid_matching_x_and_y_disjoint_samples(of):
    out = of.expand_avoid_matching_x_and_y(
        "{}_{}",
        xsamples=["A"],
        ysamples=["X", "Y"],
    )
    # No overlap, no skips.
    assert sorted(out) == ["A_X", "A_Y"]


def test_expand_avoid_matching_x_and_y_empty(of):
    assert of.expand_avoid_matching_x_and_y("{}-{}", [], []) == []


# ---------------------------------------------------------------------------
# expand_avoid_matching — variadic version with kwarg dimensions
# ---------------------------------------------------------------------------


def test_expand_avoid_matching_two_dims(of):
    out = of.expand_avoid_matching(
        "{a}_vs_{b}",
        a=["A", "B"], b=["A", "B"],
    )
    # All combos minus the matched ones (A,A) and (B,B).
    assert sorted(out) == ["A_vs_B", "B_vs_A"]


def test_expand_avoid_matching_three_dims(of):
    """All three placeholders must be distinct."""
    out = of.expand_avoid_matching(
        "{a}_{b}_{c}",
        a=["X", "Y"], b=["X", "Y"], c=["X", "Y"],
    )
    # With only 2 values per dim and 3 dims, no triplet has all-distinct
    # values — all combinations should be filtered out.
    assert out == []


def test_expand_avoid_matching_three_dims_three_values(of):
    out = of.expand_avoid_matching(
        "{a}_{b}_{c}",
        a=["X", "Y", "Z"], b=["X", "Y", "Z"], c=["X", "Y", "Z"],
    )
    # 3! = 6 distinct permutations.
    assert len(out) == 6
    assert "X_Y_Z" in out
    assert "X_X_Y" not in out


# ---------------------------------------------------------------------------
# expand_avoid_matching_all_third
# ---------------------------------------------------------------------------


def test_expand_avoid_matching_all_third_basic(of):
    out = of.expand_avoid_matching_all_third(
        "{xsamp}_{ysamp}_{third}",
        xsamp=["A", "B"], ysamp=["A", "B"], third=["t1", "t2"],
    )
    # 2 valid (xsamp, ysamp) pairs × 2 third values = 4
    assert len(out) == 4
    assert set(out) == {"A_B_t1", "A_B_t2", "B_A_t1", "B_A_t2"}


def test_expand_avoid_matching_all_third_missing_third_raises(of):
    with pytest.raises(IOError, match="third not in kwargs"):
        of.expand_avoid_matching_all_third(
            "{xsamp}_{ysamp}", xsamp=["A"], ysamp=["B"],
        )


# ---------------------------------------------------------------------------
# chromsize_to_s2c2s — small TSV parser
# ---------------------------------------------------------------------------


def test_chromsize_to_s2c2s_basic(of, tmp_path):
    p = tmp_path / "sizes.txt"
    p.write_text(dedent("""\
        sp1\tchr1\t12345
        sp1\tchr2\t67890
        sp2\tscafA\t11111
    """))
    d = of.chromsize_to_s2c2s(str(p))
    assert d == {
        "sp1": {"chr1": 12345, "chr2": 67890},
        "sp2": {"scafA": 11111},
    }


def test_chromsize_to_s2c2s_skips_blank_lines(of, tmp_path):
    p = tmp_path / "sizes.txt"
    p.write_text("\nsp1\tchr1\t100\n\nsp1\tchr2\t200\n\n")
    d = of.chromsize_to_s2c2s(str(p))
    assert d["sp1"] == {"chr1": 100, "chr2": 200}


def test_chromsize_to_s2c2s_empty_file(of, tmp_path):
    p = tmp_path / "empty.txt"
    p.write_text("")
    assert of.chromsize_to_s2c2s(str(p)) == {}


# ---------------------------------------------------------------------------
# generate_coord_structs_from_chrom_to_loc
# ---------------------------------------------------------------------------


def test_generate_coord_structs_basic(of, tmp_path):
    p = tmp_path / "test.chrom"
    p.write_text(dedent("""\
        prot1\tchrA\t+\t100\t200
        prot2\tchrB\t-\t300\t400
    """))
    out = of.generate_coord_structs_from_chrom_to_loc(str(p))
    assert out["prot_to_scaf"] == {"prot1": "chrA", "prot2": "chrB"}
    assert out["prot_to_strand"] == {"prot1": "+", "prot2": "-"}
    assert out["prot_to_start"] == {"prot1": 100, "prot2": 300}
    assert out["prot_to_stop"] == {"prot1": 200, "prot2": 400}
    # midpoint = start + (stop - start) // 2 = floor of (start + stop) / 2
    assert out["prot_to_middle"] == {"prot1": 150, "prot2": 350}


def test_generate_coord_structs_skips_blank_lines(of, tmp_path):
    p = tmp_path / "test.chrom"
    p.write_text("prot1\tchrA\t+\t1\t100\n\n\nprot2\tchrB\t-\t200\t300\n")
    out = of.generate_coord_structs_from_chrom_to_loc(str(p))
    assert len(out["prot_to_scaf"]) == 2


# ---------------------------------------------------------------------------
# open_text_maybe_gzip
# ---------------------------------------------------------------------------


def test_open_text_maybe_gzip_plain(of, tmp_path):
    p = tmp_path / "plain.txt"
    p.write_text("hello\nworld\n")
    with of.open_text_maybe_gzip(str(p)) as fh:
        assert fh.read() == "hello\nworld\n"


def test_open_text_maybe_gzip_gzipped(of, tmp_path):
    p = tmp_path / "compressed.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("hello\nworld\n")
    with of.open_text_maybe_gzip(str(p)) as fh:
        assert fh.read() == "hello\nworld\n"


def test_open_text_maybe_gzip_detects_by_magic_not_extension(of, tmp_path):
    """A file misnamed .txt but gzipped inside should still be opened
    correctly because the function looks at the magic bytes."""
    p = tmp_path / "weird.txt"
    with gzip.open(p, "wt") as fh:
        fh.write("payload\n")
    with of.open_text_maybe_gzip(str(p)) as fh:
        assert fh.read() == "payload\n"


# ---------------------------------------------------------------------------
# filter_fasta_chrom — keeps only proteins whose IDs appear in chrom_file
# ---------------------------------------------------------------------------


def test_filter_fasta_chrom_basic(of, tmp_path):
    chrom = tmp_path / "in.chrom"
    chrom.write_text("p1\tchrA\t+\t1\t100\np3\tchrA\t+\t200\t300\n")
    inp = tmp_path / "in.pep"
    inp.write_text(">p1\nMA\n>p2\nMB\n>p3\nMC\n")
    out = tmp_path / "out.pep"
    of.filter_fasta_chrom(str(chrom), str(inp), str(out))
    body = out.read_text()
    assert ">p1" in body
    assert ">p3" in body
    assert ">p2" not in body


def test_filter_fasta_chrom_gzipped_output(of, tmp_path):
    chrom = tmp_path / "in.chrom"
    chrom.write_text("p1\tchrA\t+\t1\t100\n")
    inp = tmp_path / "in.pep"
    inp.write_text(">p1\nMA\n>p2\nMB\n")
    out = tmp_path / "out.pep.gz"
    of.filter_fasta_chrom(str(chrom), str(inp), str(out))
    with gzip.open(out, "rt") as fh:
        body = fh.read()
    assert ">p1" in body
    assert ">p2" not in body


def test_filter_fasta_chrom_handles_gzipped_input(of, tmp_path):
    chrom = tmp_path / "in.chrom"
    chrom.write_text("p1\tchrA\t+\t1\t100\n")
    inp = tmp_path / "in.pep.gz"
    with gzip.open(inp, "wt") as fh:
        fh.write(">p1\nMA\n>p2\nMB\n")
    out = tmp_path / "out.pep"
    of.filter_fasta_chrom(str(chrom), str(inp), str(out))
    body = out.read_text()
    assert ">p1" in body
    assert ">p2" not in body


# ---------------------------------------------------------------------------
# check_legality — top-level config validator
# ---------------------------------------------------------------------------


def test_check_legality_passes_with_legal_config(of):
    cfg = {"species": {"SpeciesOne": {}, "SpeciesTwo": {}}}
    # Returns None (no return statement) on success. Just exercise it.
    assert of.check_legality(cfg) is None


def test_check_legality_rejects_underscore_in_sample_name(of):
    cfg = {"species": {"my_species": {}}}
    with pytest.raises(IOError, match="Sample names can't have '_' char"):
        of.check_legality(cfg)


# ---------------------------------------------------------------------------
# genome_coords_to_plotstart_dict / genome_coords_to_offset_dict
# ---------------------------------------------------------------------------


def _write_genocoords(tmp_path):
    p = tmp_path / "geno.coords"
    p.write_text(dedent("""\
        chr1 1000 1000 0
        chr2 2000 3000 1000
        chr3 1500 4500 3000
    """))
    return p


def test_genome_coords_to_plotstart_dict(of, tmp_path):
    """col4 (the per-scaffold plot starting offset) maps to the scaffold name."""
    p = _write_genocoords(tmp_path)
    d = of.genome_coords_to_plotstart_dict(str(p))
    assert d == {"chr1": 0, "chr2": 1000, "chr3": 3000}


def test_genome_coords_to_offset_dict(of, tmp_path):
    """col3 (cumsum total plot size) maps to the scaffold name."""
    p = _write_genocoords(tmp_path)
    d = of.genome_coords_to_offset_dict(str(p))
    assert d == {"chr1": 1000, "chr2": 3000, "chr3": 4500}


def test_genome_coords_dicts_skip_blank_lines(of, tmp_path):
    p = tmp_path / "geno.coords"
    p.write_text("\n\nchr1 1000 1000 0\n\nchr2 2000 3000 1000\n\n")
    assert of.genome_coords_to_plotstart_dict(str(p)) == {"chr1": 0, "chr2": 1000}


def test_genome_coords_dicts_empty_file(of, tmp_path):
    p = tmp_path / "geno.coords"
    p.write_text("")
    assert of.genome_coords_to_plotstart_dict(str(p)) == {}
    assert of.genome_coords_to_offset_dict(str(p)) == {}


# ---------------------------------------------------------------------------
# generate_coord_structs_from_chrom_to_loc (kwargs-accepting version)
# ---------------------------------------------------------------------------


def test_generate_coord_structs_kwargs_variant(of, tmp_path):
    """The kwargs-accepting variant at line 1072 behaves identically to
    the bare variant at line 935."""
    p = tmp_path / "x.chrom"
    p.write_text(dedent("""\
        prot1\tchrA\t+\t100\t200
        prot2\tchrB\t-\t300\t400
    """))
    out = of.generate_coord_structs_from_chrom_to_loc(str(p))
    assert out["prot_to_scaf"] == {"prot1": "chrA", "prot2": "chrB"}
    assert out["prot_to_start"] == {"prot1": 100, "prot2": 300}
    assert out["prot_to_stop"] == {"prot1": 200, "prot2": 400}


# ---------------------------------------------------------------------------
# general_legal_run — guards against running in install dir
# ---------------------------------------------------------------------------


def test_general_legal_run_outside_install_passes(of, tmp_path, monkeypatch):
    """Running from a directory outside the odp checkout is always legal."""
    monkeypatch.chdir(tmp_path)
    of.general_legal_run()  # should return None, no exception


def test_general_legal_run_inside_install_crashes(of, monkeypatch):
    """Running inside the odp install path (anywhere under odp/ except
    tests/test_odp_basic) raises."""
    import source.odp_functions  # for resolving the install dir
    # Use the repo root as cwd — that's "inside the install dir".
    import sys
    odp_module_dir = None
    for p in sys.path:
        candidate = Path(p) / "odp_functions.py"
        if candidate.is_file():
            odp_module_dir = Path(p)
            break
    if odp_module_dir is None:
        pytest.skip("could not resolve odp install dir for this test")
    install_root = odp_module_dir.parent
    monkeypatch.chdir(install_root)
    with pytest.raises(ValueError, match="odp install directory"):
        of.general_legal_run()
