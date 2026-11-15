"""Tests for source/odp_fasta_tools.py, the afp-based steps of the
filter_prots rule that replaced seqtk and fainfo."""
from __future__ import annotations

import bz2
import gzip
import os
import random
import subprocess

import pytest


@pytest.fixture(scope="module")
def FT(source_dir, fasta_parser_dir):
    import odp_fasta_tools
    return odp_fasta_tools


@pytest.fixture(scope="module")
def afp(fasta_parser_dir):
    import afp
    return afp


def _write(path, data: bytes, how="plain"):
    if how == "gzip":
        with gzip.open(path, "wb") as fh:
            fh.write(data)
    elif how == "bzip2":
        with bz2.open(path, "wb") as fh:
            fh.write(data)
    else:
        path.write_bytes(data)
    return path


def _random_fasta(rng, eols=("\n", "\r\n", "\r")) -> bytes:
    """FASTA text with the awkward cases: CRLF and lone-CR line endings, blank
    and whitespace-only lines, spaces and tabs inside sequence lines,
    descriptions, empty names and records, a '>' inside a sequence line, very
    long lines and no final newline."""
    out = [rng.choice(["", "\n", "  \n", "\t" + eols[-1]])]
    for _ in range(rng.randint(1, 6)):
        name = rng.choice(["chr1", "scaf_2.1", "X|y", "", "c%d" % rng.randint(0, 99)])
        desc = rng.choice(["", " len=5", "\tsome description", "   two  spaces "])
        out.append(">" + rng.choice(["", " "]) + name + desc + rng.choice(eols))
        for _ in range(rng.randint(0, 5)):
            kind = rng.random()
            if kind < 0.1:
                line = ""
            elif kind < 0.2:
                line = rng.choice([" ", "\t", "  \t "])
            elif kind < 0.3:
                line = "AC GT\tN" + " " * rng.randint(0, 2)
            elif kind < 0.35:
                line = "AC>GT"
            else:
                line = "".join(rng.choice("ACGTNacgtn") for _ in range(rng.randint(1, 150)))
            out.append(line + rng.choice(eols))
    text = "".join(out)
    if rng.random() < 0.3:
        text = text.rstrip("\r\n")
    return text.encode()


CHUNKS = [1, 2, 3, 5, 8, 13, 64, 1 << 22]


def test_chrom_sizes_match_afp_on_random_files(FT, afp, tmp_path):
    rng = random.Random(20261114)
    for k in range(300):
        # half the files have only \n line endings, which take the fast route
        path = _write(tmp_path / ("r%d.fa" % k), _random_fasta(rng, ("\n",) if k % 2 else ("\n", "\r\n", "\r")))
        expected = [(r.id.encode(), len(r.seq)) for r in afp.parse_fasta(path)]
        for chunk in CHUNKS:
            got = list(FT.chrom_sizes(path, chunk_size=chunk))
            assert got == expected, (k, chunk, path.read_bytes()[:300])


def test_chrom_sizes_compressed(FT, tmp_path):
    data = b">a desc\nACGT\nAC\n>b\n" + b"N" * 1000 + b"\n"
    for how in ("plain", "gzip", "bzip2"):
        path = _write(tmp_path / ("g." + how), data, how)
        assert list(FT.chrom_sizes(path, chunk_size=7)) == [(b"a", 6), (b"b", 1000)]


def test_chrom_sizes_unwrapped_long_line_uses_small_chunks(FT, tmp_path):
    path = _write(tmp_path / "long.fa", b">chr1\n" + b"A" * 300000 + b"\n>chr2\nCC")
    assert list(FT.chrom_sizes(path, chunk_size=4096)) == [(b"chr1", 300000), (b"chr2", 2)]


def test_chrom_sizes_header_at_end_and_empty_file(FT, tmp_path):
    assert list(FT.chrom_sizes(_write(tmp_path / "e.fa", b""))) == []
    path = _write(tmp_path / "h.fa", b">a\nAC\n>last")
    assert list(FT.chrom_sizes(path, chunk_size=2)) == [(b"a", 2), (b"last", 0)]


def test_chrom_sizes_sequence_before_header_raises(FT, afp, tmp_path):
    path = _write(tmp_path / "bad.fa", b"\nACGT\n>a\nAC\n")
    with pytest.raises(ValueError):
        list(afp.parse_fasta(path))
    with pytest.raises(ValueError):
        list(FT.chrom_sizes(path))


def test_write_chrom_sizes_format(FT, tmp_path):
    path = _write(tmp_path / "g.fa", b">c1 x\nACGT\n>c2\nA\n")
    out = tmp_path / "S.chromsize"
    assert FT.write_chrom_sizes(path, "SP", out) == 2
    assert out.read_bytes() == b"SP\tc1\t4\nSP\tc2\t1\n"


def test_chrom_ids(FT, tmp_path):
    text = b"p1\tc1\t+\t10\t20\np2 extra\tc1\t-\t30\t40\n\n   \np1\tc2\t+\t1\t2\n"
    for how in ("plain", "gzip"):
        path = _write(tmp_path / ("x.chrom." + how), text, how)
        assert FT.chrom_ids(path) == {"p1", "p2"}


def test_filter_proteins_keeps_order_duplicates_and_descriptions(FT, tmp_path):
    seq = "M" + "A" * 169
    prots = _write(tmp_path / "p.fa", (
        ">p3 third protein\n" + seq[:60] + "\n" + seq[60:] + "\n"
        ">skip\nMMM\n"
        ">p1\nMK\n"
        ">p3 again\nMV\n"
        ">p2\n\n").encode(), "gzip")
    out = tmp_path / "out.pep.gz"
    written, missing = FT.filter_proteins(prots, {"p1", "p2", "p3", "absent"}, out)
    # p2 has no sequence: seqtk subseq left such records out, and so does odp
    assert written == 3 and missing == {"absent", "p2"}
    with gzip.open(out, "rt") as fh:
        text = fh.read()
    assert text == (">p3 third protein\n" + seq[:80] + "\n" + seq[80:160] + "\n" + seq[160:] + "\n"
                    ">p1\nMK\n>p3 again\nMV\n")


# The old tools, when someone has them: ODP_SEQTK=/path/seqtk ODP_FAINFO=/path/fainfo pytest
SEQTK, FAINFO = os.environ.get("ODP_SEQTK"), os.environ.get("ODP_FAINFO")


@pytest.mark.skipif(not (SEQTK and FAINFO), reason="set ODP_SEQTK and ODP_FAINFO to compare with the old tools")
def test_same_output_as_seqtk_and_fainfo(FT, tmp_path):
    rng = random.Random(7)
    seqs = ["".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(rng.randint(0, 400))) for _ in range(40)]
    names = ["prot%d" % i for i in range(40)]
    fasta = "".join(">%s desc %d\n%s\n" % (n, i, s) for i, (n, s) in enumerate(zip(names, seqs)))
    prots = _write(tmp_path / "p.fa.gz", fasta.encode(), "gzip")
    ids_path = tmp_path / "ids.txt"
    ids = set(names[::3])
    ids_path.write_text("".join(i + "\n" for i in sorted(ids)))
    old = subprocess.run("%s subseq %s %s | %s seq -l 80" % (SEQTK, prots, ids_path, SEQTK),
                         shell=True, check=True, capture_output=True).stdout
    FT.filter_proteins(prots, ids, tmp_path / "new.pep.gz")
    with gzip.open(tmp_path / "new.pep.gz", "rb") as fh:
        assert fh.read() == old

    genome = _write(tmp_path / "g.fa.gz", b">chr1 x\n" + b"ACGT" * 5000 + b"\n>chr2\nAC\nGT\n>chr3\n", "gzip")
    old = subprocess.run([FAINFO, "SP", str(genome)], check=True, capture_output=True).stdout
    FT.write_chrom_sizes(genome, "SP", tmp_path / "new.chromsize")
    assert (tmp_path / "new.chromsize").read_bytes() == old
