"""afp — Another Fastx Parser.

Tiny, dependency-free single-file reader/writer for FASTA and FASTQ files
with transparent gzip / bzip2 / zip / zstandard decompression.

This module is intentionally a single file. To vendor it into another project,
drop ``afp.py`` somewhere on the python path and ``import afp``.

Quick start:

    import afp

    for rec in afp.parse("reads.fq.gz"):
        print(rec.id, len(rec.seq))

    for rec in afp.parse("genome.fa"):
        print(rec.id, rec.desc, rec.seq[:50])

Public API:

    parse(path, format=None)        # auto-detects FASTA vs FASTQ
    parse_fasta(path)
    parse_fastq(path)
    Record                          # id, seq, desc, qual (qual=None for FASTA)
    get_open_func(path)             # the right open() for the file's compression
    detect_compression(path)        # 'gzip' | 'bzip2' | 'zip' | 'zstd' | 'none'
    write(records, path, ...)       # write iterable of Records as FASTA or FASTQ
"""
from __future__ import annotations

import bz2
import gzip
import io
import pathlib
import zipfile
from dataclasses import dataclass
from typing import Callable, Iterable, Iterator, Optional, Union

__version__ = "0.1.2"

PathLike = Union[str, pathlib.Path]

__all__ = [
    "Record",
    "parse",
    "parse_fasta",
    "parse_fastq",
    "get_open_func",
    "detect_compression",
    "write",
    "__version__",
]


# ---------------------------------------------------------------------------
# Record
# ---------------------------------------------------------------------------


@dataclass
class Record:
    """A sequence record from a FASTA or FASTQ file.

    Attributes are mutable so callers can rewrite ids in place.

    Args:
        id: identifier (the token after `>` or `@` up to the first whitespace).
        seq: the sequence as a single string with newlines stripped.
        desc: optional description (everything after the first whitespace on
            the header line). May be ``None``.
        qual: per-base quality string (FASTQ only). ``None`` for FASTA records.
    """

    id: str
    seq: str
    desc: Optional[str] = None
    qual: Optional[str] = None

    @property
    def description(self) -> str:
        """The full header line as it would appear in the file, including the
        leading ``>`` (for FASTA) or ``@`` (for FASTQ). Use ``format()`` for
        the full multi-line record."""
        prefix = "@" if self.qual is not None else ">"
        if self.desc:
            return f"{prefix}{self.id} {self.desc}"
        return f"{prefix}{self.id}"

    def __len__(self) -> int:
        return len(self.seq)

    def __iter__(self):
        return iter(self.seq)

    def __contains__(self, item) -> bool:
        return item in self.seq

    def format(self, wrap: Optional[int] = None) -> str:
        """Return the record as it would appear written to disk.

        For FASTA, wraps sequence lines at ``wrap`` columns if specified
        (default: no wrap). For FASTQ, wrapping is not applied — the record is
        emitted as the standard 4-line form.
        """
        if self.qual is not None:
            return f"@{self._header_tail()}\n{self.seq}\n+\n{self.qual}\n"

        header = f">{self._header_tail()}"
        if wrap is None or wrap <= 0:
            return f"{header}\n{self.seq}\n"

        lines = [self.seq[i : i + wrap] for i in range(0, len(self.seq), wrap)]
        return header + "\n" + "\n".join(lines) + "\n"

    def _header_tail(self) -> str:
        return f"{self.id} {self.desc}" if self.desc else self.id


# ---------------------------------------------------------------------------
# Compression detection + opener selection
# ---------------------------------------------------------------------------


# Magic-byte prefixes. Detection is by content, not filename — a misnamed file
# still works as long as its content is valid.
_MAGIC = {
    b"\x1f\x8b": "gzip",
    b"BZh": "bzip2",
    b"PK\x03\x04": "zip",
    b"PK\x05\x06": "zip",  # empty zip
    b"\x28\xb5\x2f\xfd": "zstd",
}


def detect_compression(path: PathLike) -> str:
    """Return one of ``'gzip'``, ``'bzip2'``, ``'zip'``, ``'zstd'``, or
    ``'none'``. Reads only the first 4 bytes of the file."""
    with open(path, "rb") as fh:
        head = fh.read(4)
    for magic, name in _MAGIC.items():
        if head.startswith(magic):
            return name
    return "none"


def _open_zstd(path: PathLike, mode: str):
    try:
        import zstandard  # type: ignore
    except ImportError as e:
        raise RuntimeError(
            "zstandard package required to read .zst files. "
            "Install with: pip install zstandard"
        ) from e
    dctx = zstandard.ZstdDecompressor()
    if "b" in mode:
        return dctx.stream_reader(open(path, "rb"))
    return io.TextIOWrapper(dctx.stream_reader(open(path, "rb")), encoding="utf-8")


def _open_zip(path: PathLike, mode: str):
    """Open the first file inside a zip archive. Multi-member zips are rejected
    to avoid silent data loss."""
    zf = zipfile.ZipFile(path, "r")
    names = [n for n in zf.namelist() if not n.endswith("/")]
    if len(names) != 1:
        raise ValueError(
            f"zip archive {path!s} must contain exactly one file, "
            f"got {len(names)}: {names!r}"
        )
    stream = zf.open(names[0], "r")
    if "b" in mode:
        return stream
    return io.TextIOWrapper(stream, encoding="utf-8")


def get_open_func(path: PathLike) -> Callable:
    """Return an ``open``-style callable matching ``path``'s compression.

    The returned callable has the signature ``open(path, mode='rt')`` and
    transparently decompresses on read.
    """
    comp = detect_compression(path)
    if comp == "gzip":
        return gzip.open
    if comp == "bzip2":
        return bz2.open
    if comp == "zip":
        return _open_zip
    if comp == "zstd":
        return _open_zstd
    return open


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------


def _split_header(line: str):
    """Split a header line (without its leading '>' or '@') into (id, desc).

    Splits on the first run of any whitespace, so headers separated by a
    tab, a space, or both work the same way."""
    line = line.rstrip("\r\n")
    parts = line.split(None, 1)
    if not parts:
        return "", None
    if len(parts) == 1:
        return parts[0], None
    return parts[0], parts[1]


def parse_fasta(path: PathLike) -> Iterator[Record]:
    """Yield ``Record`` instances from a FASTA file. Multi-line sequences are
    joined into a single string per record."""
    opener = get_open_func(path)
    with opener(path, "rt") as fh:
        rec_id: Optional[str] = None
        rec_desc: Optional[str] = None
        seq_parts: list[str] = []
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                if rec_id is not None:
                    yield Record(id=rec_id, seq="".join(seq_parts), desc=rec_desc)
                rec_id, rec_desc = _split_header(line[1:])
                seq_parts = []
            elif rec_id is None:
                raise ValueError(
                    f"sequence line before any record header in {path!s}: "
                    f"{line[:60]!r}"
                )
            else:
                seq_parts.append(line)
        if rec_id is not None:
            yield Record(id=rec_id, seq="".join(seq_parts), desc=rec_desc)


def parse_fastq(path: PathLike) -> Iterator[Record]:
    """Yield ``Record`` instances from a FASTQ file. Sequence and quality
    strings are joined across continuation lines (rare but legal).

    Records with sequence and quality of differing length raise ``ValueError``.
    """
    opener = get_open_func(path)
    with opener(path, "rt") as fh:
        lineno = 0
        while True:
            header = fh.readline()
            if not header:
                return
            lineno += 1
            header = header.rstrip("\r\n")
            if not header:
                continue
            if not header.startswith("@"):
                raise ValueError(
                    f"expected FASTQ header starting with '@' at line {lineno} "
                    f"of {path!s}, got: {header[:60]!r}"
                )
            rec_id, rec_desc = _split_header(header[1:])

            seq_parts: list[str] = []
            while True:
                line = fh.readline()
                if not line:
                    raise ValueError(
                        f"unexpected EOF inside FASTQ record {rec_id!r} in {path!s}"
                    )
                lineno += 1
                line = line.rstrip("\r\n")
                if line.startswith("+"):
                    break
                seq_parts.append(line)
            seq = "".join(seq_parts)

            qual_parts: list[str] = []
            need = len(seq)
            while sum(len(p) for p in qual_parts) < need:
                line = fh.readline()
                if not line:
                    raise ValueError(
                        f"unexpected EOF in qualities for record {rec_id!r} "
                        f"in {path!s}"
                    )
                lineno += 1
                qual_parts.append(line.rstrip("\r\n"))
            qual = "".join(qual_parts)
            if len(qual) != len(seq):
                raise ValueError(
                    f"record {rec_id!r} in {path!s}: sequence length {len(seq)} "
                    f"!= quality length {len(qual)}"
                )

            yield Record(id=rec_id, seq=seq, desc=rec_desc, qual=qual)


def parse(path: PathLike, format: Optional[str] = None) -> Iterator[Record]:
    """Yield ``Record`` instances from a FASTA or FASTQ file.

    ``format`` may be ``'fasta'``, ``'fastq'``, or ``None`` to auto-detect from
    the file's first non-empty byte (``>`` for FASTA, ``@`` for FASTQ).
    Empty files yield zero records.
    """
    fmt = (format or _autodetect_format(path)).lower()
    if fmt == "fasta":
        yield from parse_fasta(path)
    elif fmt == "fastq":
        yield from parse_fastq(path)
    else:
        raise ValueError(f"unsupported format: {format!r}")


def _autodetect_format(path: PathLike) -> str:
    """Sniff the first non-empty line. Empty files default to ``'fasta'`` so
    callers iterating over ``parse()`` see an empty record stream instead of
    a ``ValueError``."""
    opener = get_open_func(path)
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.lstrip()
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("@"):
                return "fastq"
            raise ValueError(
                f"could not detect FASTA/FASTQ format from first byte of "
                f"{path!s}: {line[:60]!r}"
            )
    return "fasta"


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


def write(
    records: Iterable[Record],
    path: PathLike,
    wrap: Optional[int] = None,
    compression: Optional[str] = None,
) -> int:
    """Write an iterable of ``Record`` to ``path``. Returns the number of
    records written.

    Output format mirrors the records: if the first record has a ``qual``
    field, output is FASTQ; otherwise FASTA. Mixing types in a single output
    is rejected.

    ``compression`` overrides the auto-detection by extension (``.gz``,
    ``.bz2``). If ``None``, output is plain text unless the path ends with
    ``.gz`` / ``.bz2``.
    """
    path_str = str(path)
    if compression is None:
        if path_str.endswith(".gz"):
            compression = "gzip"
        elif path_str.endswith(".bz2"):
            compression = "bzip2"
        else:
            compression = "none"

    if compression == "gzip":
        opener = gzip.open
    elif compression == "bzip2":
        opener = bz2.open
    elif compression == "none":
        opener = open
    else:
        raise ValueError(f"unsupported output compression: {compression!r}")

    count = 0
    first_is_fastq: Optional[bool] = None
    with opener(path, "wt") as fh:
        for rec in records:
            is_fastq = rec.qual is not None
            if first_is_fastq is None:
                first_is_fastq = is_fastq
            elif is_fastq != first_is_fastq:
                raise ValueError(
                    "cannot mix FASTA and FASTQ records in a single output stream"
                )
            fh.write(rec.format(wrap=wrap))
            count += 1
    return count
