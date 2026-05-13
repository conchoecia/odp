import os
import sys
import re
import gzip
from pathlib import Path
import pytest
import hashlib
import itertools


# fasta-parser, odp_filechecker source extraction, etc. — load lazily so
# conftest fixtures get a chance to register the necessary sys.path entries
# (source_dir, scripts_dir, fasta_parser_dir) before we resolve any imports.
@pytest.fixture(scope="module")
def loaded_legality_funcs(scripts_dir, fasta_parser_dir, repo_root):
    import afp as fasta  # noqa: E402

    class _ODPFunctions:
        @staticmethod
        def check_file_exists(path):
            if not Path(path).is_file():
                raise IOError(f"File does not exist: {path}")

        @staticmethod
        def chrom_file_is_legal(path):
            return True

        @staticmethod
        def open_text_maybe_gzip(path, encoding="utf-8"):
            with open(path, "rb") as fh:
                head = fh.read(2)
            if head == b"\x1f\x8b":
                return gzip.open(path, "rt", encoding=encoding)
            return open(path, "rt", encoding=encoding)

    odpf = _ODPFunctions()

    source_text = (repo_root / "scripts" / "odp_filechecker").read_text()

    def _load_func(name):
        start = source_text.index(f"def {name}")
        end = source_text.index("    return True", start) + len("    return True")
        code = source_text[start:end]
        namespace = {
            "fasta": fasta,
            "odpf": odpf,
            "hashlib": hashlib,
            "os": os,
            "itertools": itertools,
        }
        exec(code, namespace)
        return namespace[name]

    return {
        "check_genome_file_legality": _load_func("check_genome_file_legality"),
        "check_protein_file_legality": _load_func("check_protein_file_legality"),
        "check_chrom_file_legality": _load_func("check_chrom_file_legality"),
    }


def write_fasta(tmp_path, text):
    path = tmp_path / "test.fasta"
    path.write_text(text)
    return str(path)


def test_genome_illegal_chars(tmp_path, loaded_legality_funcs):
    fasta_path = write_fasta(
        tmp_path,
        ">seq1\nACGT\n>seq2\nACGTP\n",
    )
    with pytest.raises(IOError):
        loaded_legality_funcs["check_genome_file_legality"](fasta_path)


def test_protein_illegal_chars(tmp_path, loaded_legality_funcs):
    fasta_path = write_fasta(
        tmp_path,
        ">seq1\nACDEFGHIKLMNPQRSTVWY1\n",
    )
    with pytest.raises(IOError):
        loaded_legality_funcs["check_protein_file_legality"](fasta_path)


def test_genome_empty(tmp_path, loaded_legality_funcs):
    fasta_path = write_fasta(tmp_path, "")
    with pytest.raises(IOError):
        loaded_legality_funcs["check_genome_file_legality"](fasta_path)


def test_protein_empty(tmp_path, loaded_legality_funcs):
    fasta_path = write_fasta(tmp_path, "")
    with pytest.raises(IOError):
        loaded_legality_funcs["check_protein_file_legality"](fasta_path)


def test_chrom_empty(tmp_path, loaded_legality_funcs):
    genome_path = tmp_path / "genome.fasta"
    genome_path.write_text(">s1\nACGT\n")
    protein_path = tmp_path / "protein.fasta"
    protein_path.write_text(">p1\nM\n")
    chrom_path = tmp_path / "test.chrom"
    chrom_path.write_text("")
    with pytest.raises(IOError):
        loaded_legality_funcs["check_chrom_file_legality"](
            str(chrom_path), str(genome_path), str(protein_path)
        )
