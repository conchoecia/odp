"""Tests for source/odp_functions.py:check_species_input_legality.

This is the biggest single function in odp_functions and validates the
three input files for one species: genome FASTA, protein FASTA, and
.chrom. Each file has multiple failure modes that produce a multi-line
error message; we cover the happy path and each failure mode.
"""
from __future__ import annotations

from textwrap import dedent

import pytest


@pytest.fixture(scope="module")
def of(source_dir, fasta_parser_dir):
    import odp_functions
    return odp_functions


def _write_minimal_inputs(tmp_path):
    """Write a clean, mutually-consistent (genome, protein, chrom)
    trio that should pass every check."""
    g = tmp_path / "g.fa"
    g.write_text(">chr1\nACGTACGT\n>chr2\nACGTACGT\n")
    p = tmp_path / "p.pep"
    p.write_text(">p1\nMASE\n>p2\nMRTI\n")
    c = tmp_path / "c.chrom"
    c.write_text(
        "p1\tchr1\t+\t100\t200\n"
        "p2\tchr2\t-\t300\t400\n"
    )
    return g, p, c


# ---------------------------------------------------------------------------
# Type check on dup_proteins_allowed
# ---------------------------------------------------------------------------


def test_dup_proteins_allowed_must_be_bool(of, tmp_path):
    g, p, c = _write_minimal_inputs(tmp_path)
    with pytest.raises(IOError, match="must be a boolean"):
        of.check_species_input_legality(str(g), str(p), str(c), dup_proteins_allowed="yes")


# ---------------------------------------------------------------------------
# Genome FASTA checks
# ---------------------------------------------------------------------------


def test_missing_genome_file(of, tmp_path):
    _, p, c = _write_minimal_inputs(tmp_path)
    with pytest.raises(IOError, match="does not exist"):
        of.check_species_input_legality(
            str(tmp_path / "nope.fa"), str(p), str(c),
        )


def test_duplicate_genome_headers_rejected(of, tmp_path):
    """Two records with the same `>chr1` ID is rejected because we can't
    distinguish them downstream."""
    g = tmp_path / "g.fa"
    g.write_text(">chr1\nACGT\n>chr1\nGGGG\n")
    p = tmp_path / "p.pep"
    p.write_text(">p1\nM\n")
    c = tmp_path / "c.chrom"
    c.write_text("p1\tchr1\t+\t1\t10\n")
    with pytest.raises(IOError, match="duplicate sequence headers"):
        of.check_species_input_legality(str(g), str(p), str(c))


# ---------------------------------------------------------------------------
# Protein FASTA checks
# ---------------------------------------------------------------------------


def test_missing_protein_file(of, tmp_path):
    g, _, c = _write_minimal_inputs(tmp_path)
    with pytest.raises(IOError, match="does not exist"):
        of.check_species_input_legality(
            str(g), str(tmp_path / "nope.pep"), str(c),
        )


def test_duplicate_protein_headers_rejected(of, tmp_path):
    g, _, _ = _write_minimal_inputs(tmp_path)
    p = tmp_path / "p.pep"
    p.write_text(">p1\nMASE\n>p1\nMOTHER\n")  # duplicate header
    c = tmp_path / "c.chrom"
    c.write_text("p1\tchr1\t+\t1\t10\n")
    with pytest.raises(IOError, match="duplicate sequence headers|protein"):
        of.check_species_input_legality(str(g), str(p), str(c))


def test_duplicate_protein_sequences_rejected_by_default(of, tmp_path):
    """Two proteins with the same sequence is rejected unless the caller
    opts in via dup_proteins_allowed=True."""
    g, _, _ = _write_minimal_inputs(tmp_path)
    p = tmp_path / "p.pep"
    p.write_text(">p1\nMASE\n>p2\nMASE\n")  # same sequence
    c = tmp_path / "c.chrom"
    c.write_text(
        "p1\tchr1\t+\t1\t10\n"
        "p2\tchr2\t+\t20\t30\n"
    )
    with pytest.raises(IOError, match="duplicate"):
        of.check_species_input_legality(str(g), str(p), str(c))


def test_duplicate_protein_sequences_allowed_when_flag_set(of, tmp_path):
    """With dup_proteins_allowed=True, identical sequences pass."""
    g, _, _ = _write_minimal_inputs(tmp_path)
    p = tmp_path / "p.pep"
    p.write_text(">p1\nMASE\n>p2\nMASE\n")
    c = tmp_path / "c.chrom"
    c.write_text(
        "p1\tchr1\t+\t1\t10\n"
        "p2\tchr2\t+\t20\t30\n"
    )
    assert of.check_species_input_legality(
        str(g), str(p), str(c), dup_proteins_allowed=True,
    ) is True


# ---------------------------------------------------------------------------
# Chrom file checks
# ---------------------------------------------------------------------------


def test_missing_chrom_file(of, tmp_path):
    g, p, _ = _write_minimal_inputs(tmp_path)
    with pytest.raises(IOError, match="does not exist"):
        of.check_species_input_legality(str(g), str(p), str(tmp_path / "nope.chrom"))


def test_chrom_has_protein_not_in_pep(of, tmp_path):
    g, p, _ = _write_minimal_inputs(tmp_path)
    c = tmp_path / "c.chrom"
    c.write_text(
        "p1\tchr1\t+\t1\t10\n"
        "MISSING\tchr2\t+\t20\t30\n"
    )
    with pytest.raises(IOError, match="not seen in the protein"):
        of.check_species_input_legality(str(g), str(p), str(c))


def test_chrom_has_scaffold_not_in_genome(of, tmp_path):
    g, p, _ = _write_minimal_inputs(tmp_path)
    c = tmp_path / "c.chrom"
    c.write_text(
        "p1\tchr1\t+\t1\t10\n"
        "p2\tphantom_scaffold\t+\t20\t30\n"
    )
    with pytest.raises(IOError, match="not seen in the genome"):
        of.check_species_input_legality(str(g), str(p), str(c))


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_happy_path_returns_true(of, tmp_path):
    g, p, c = _write_minimal_inputs(tmp_path)
    assert of.check_species_input_legality(str(g), str(p), str(c)) is True
