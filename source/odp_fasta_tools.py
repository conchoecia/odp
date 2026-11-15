"""FASTA steps of odp's filter_prots rule, in pure Python on afp.

These replace two programs that odp used to compile from Heng Li's kseq.h:
seqtk (``seqtk subseq`` piped into ``seqtk seq -l 80``) and fainfo (one
``sample<TAB>name<TAB>length`` line per sequence). Without them, installing
odp needs no C compiler, no zlib headers and no download of seqtk.

For ordinary FASTA files the outputs are the same as before:

- ``filter_proteins`` keeps the records whose ID is in the .chrom file, in
  file order (an ID that appears twice is written twice, as seqtk did), with
  the header's description, and wraps the sequence at 80 columns. Records
  with an empty sequence are left out, as ``seqtk subseq`` left them out.
- ``chrom_sizes`` gives each sequence's name (the header's first word) and
  length. It reads the genome in fixed-size binary chunks and never holds a
  sequence in memory, so a chromosome written on a single line costs no more
  memory than a short one. fainfo held each whole sequence in memory.

Both read plain, gzip, bzip2, zip and zstandard files, through afp.
"""
from __future__ import annotations

import gzip
import os
import re
import sys

_deps = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dependencies")
if _deps not in sys.path:
    sys.path.insert(1, _deps)
import afp  # noqa: E402

# afp reads text with universal newlines, so \r, \n and \r\n all end a line. Chunks without a \r take
# the fast route: bytes.find and bytes.count, which run at memory speed.
_EOL = re.compile(rb"[\r\n]")
_HEADER_START = re.compile(rb"[\r\n]>")
_SPACE = b" \t\r\n\x0b\x0c"


def _first_eol(buf, pos, has_cr):
    """Index of the first line ending at or after pos, or -1."""
    if has_cr:
        m = _EOL.search(buf, pos)
        return m.start() if m else -1
    return buf.find(b"\n", pos)


def _last_eol(buf, has_cr):
    """Index of the last line ending in buf, or -1."""
    return max(buf.rfind(b"\n"), buf.rfind(b"\r")) if has_cr else buf.rfind(b"\n")


def _next_header(buf, pos, has_cr):
    """Index of the '>' that starts the next header line after pos, or len(buf)."""
    if has_cr:
        m = _HEADER_START.search(buf, pos)
        return m.start() + 1 if m else len(buf)
    i = buf.find(b"\n>", pos)
    return i + 1 if i >= 0 else len(buf)


def _has_inner_space(seg):
    return b" " in seg or b"\t" in seg or b"\x0b" in seg or b"\x0c" in seg


def chrom_ids(chrom_path):
    """The protein IDs of a .chrom file (plain or compressed), as a set:
    the first word of column 1 on every line that has one."""
    ids = set()
    opener = afp.get_open_func(chrom_path)
    with opener(chrom_path, "rt") as fh:
        for line in fh:
            words = line.split("\t", 1)[0].split()
            if words:
                ids.add(words[0])
    return ids


def filter_proteins(prots_path, ids, out_path, wrap=80, compresslevel=6):
    """Write to ``out_path`` (gzip) the records of ``prots_path`` whose ID is
    in ``ids`` and whose sequence is not empty, in file order, as
    ``seqtk subseq | seqtk seq -l 80`` wrote them: the header with its
    description, then the sequence in lines of ``wrap`` characters.
    Returns the number of records written and the set of IDs not written."""
    found = set()
    written = 0
    with gzip.open(out_path, "wt", compresslevel=compresslevel) as out:
        for rec in afp.parse_fasta(prots_path):
            if rec.seq and rec.id in ids:
                out.write(rec.format(wrap=wrap))
                found.add(rec.id)
                written += 1
    return written, set(ids) - found


def chrom_sizes(fasta_path, chunk_size=1 << 22):
    """Yield ``(name, length)`` for each sequence of a FASTA file (plain or
    compressed). Names are bytes. The values are those of ``afp.parse_fasta``:
    the name is the header's first word, and the length counts the characters
    of the sequence lines, leaving out line endings and blank lines.

    The file is read in binary chunks of ``chunk_size`` bytes and sequence
    lines are only counted, so memory stays flat whatever the line lengths.
    A sequence line before the first header raises ``ValueError``, as afp does.
    """
    opener = afp.get_open_func(fasta_path)
    name = None        # current record's name; None until the first header
    length = 0         # characters counted for the current record
    at_line_start = True
    line_text = False  # the current line has a character other than whitespace
    line_space = 0     # whitespace counted on the current line before its first text
    carry = b""        # a header line that continues into the next chunk

    with opener(fasta_path, "rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            buf = carry + chunk if carry else chunk
            carry = b""
            has_cr = b"\r" in buf
            pos, n = 0, len(buf)
            while pos < n:
                if at_line_start and buf[pos] == 62:  # ">": a header line
                    eol = _first_eol(buf, pos, has_cr)
                    if eol < 0:
                        carry = buf[pos:]
                        break
                    if name is not None:
                        yield name, length
                    words = buf[pos + 1:eol].split(None, 1)
                    name, length = (words[0] if words else b""), 0
                    pos = eol + 1
                    line_text, line_space = False, 0
                    continue
                # Sequence text up to the next header (or the end of the chunk).
                end = _next_header(buf, pos, has_cr)
                seg = buf[pos:end]
                pos = end
                if name is None and seg.translate(None, _SPACE):
                    raise ValueError(
                        f"sequence line before any record header in {fasta_path!s}")
                if not _has_inner_space(seg):
                    # No whitespace other than line endings: every other byte is sequence.
                    length += len(seg) - seg.count(b"\n") - (seg.count(b"\r") if has_cr else 0)
                    first = _first_eol(seg, 0, has_cr)
                    if first < 0:
                        line_text, line_space = True, 0
                        at_line_start = False
                        continue
                    if first > 0:
                        line_text = True
                    if not line_text:
                        length -= line_space  # the line that just ended was blank
                    tail = len(seg) - 1 - _last_eol(seg, has_cr)
                    line_text, line_space = tail > 0, 0
                    at_line_start = tail == 0
                    continue
                # Whitespace inside the text: go line by line, as afp does.
                pieces = _EOL.split(seg)
                last = len(pieces) - 1
                for i, piece in enumerate(pieces):
                    if piece:
                        length += len(piece)
                        at_line_start = False
                        if not line_text:
                            if piece.isspace():
                                line_space += len(piece)
                            else:
                                line_text, line_space = True, 0
                    if i < last:  # a line ending follows this piece
                        if not line_text:
                            length -= line_space
                        line_text, line_space = False, 0
                        at_line_start = True

    if carry:  # the file ends inside a header line
        if name is not None:
            yield name, length
        words = carry[1:].split(None, 1)
        name, length = (words[0] if words else b""), 0
    elif not line_text:
        length -= line_space
    if name is not None:
        yield name, length


def write_chrom_sizes(fasta_path, sample, out_path, chunk_size=1 << 22):
    """Write ``sample<TAB>name<TAB>length`` for each sequence of
    ``fasta_path``, one line each, as fainfo did. Returns the number of
    sequences."""
    prefix = sample.encode()
    count = 0
    with open(out_path, "wb") as out:
        for name, length in chrom_sizes(fasta_path, chunk_size):
            out.write(b"%s\t%s\t%d\n" % (prefix, name, length))
            count += 1
    return count
