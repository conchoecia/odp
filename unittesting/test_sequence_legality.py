import os
import sys
import re
from pathlib import Path
import pytest
import hashlib

# set up paths to import dependencies
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(REPO_ROOT / "dependencies" / "fasta-parser"))

import fasta  # noqa: E402
# provide minimal replacement for odp_functions.check_file_exists

class _ODPFunctions:
    @staticmethod
    def check_file_exists(path):
        if not Path(path).is_file():
            raise IOError(f"File does not exist: {path}")

odpf = _ODPFunctions()

# extract function definitions from Snakemake file
SOURCE = (REPO_ROOT / "scripts" / "odp_filechecker").read_text()

def _load_func(name):
    start = SOURCE.index(f"def {name}")
    end = SOURCE.index("    return True", start) + len("    return True")
    code = SOURCE[start:end]
    namespace = {
        'fasta': fasta,
        'odpf': odpf,
        'hashlib': hashlib,
    }
    exec(code, namespace)
    return namespace[name]

check_genome_file_legality = _load_func('check_genome_file_legality')
check_protein_file_legality = _load_func('check_protein_file_legality')


def write_fasta(tmp_path, text):
    path = tmp_path / "test.fasta"
    path.write_text(text)
    return str(path)


def test_genome_illegal_chars(tmp_path):
    fasta_path = write_fasta(
        tmp_path,
        ">seq1\nACGT\n>seq2\nACGTP\n",
    )
    with pytest.raises(IOError):
        check_genome_file_legality(fasta_path)


def test_protein_illegal_chars(tmp_path):
    fasta_path = write_fasta(
        tmp_path,
        ">seq1\nACDEFGHIKLMNPQRSTVWY1\n",
    )
    with pytest.raises(IOError):
        check_protein_file_legality(fasta_path)
