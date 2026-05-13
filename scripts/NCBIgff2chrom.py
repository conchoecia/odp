#!/usr/bin/env python3
"""
NCBIgff2chrom.py
================

Parse an NCBI GFF annotation and produce a `.chrom` file (ODP format).

This script supports two invocation styles:

1.  **Legacy / odp:** one positional argument (the GFF, possibly gzipped).
    Emits the chrom table to stdout. CDS extents are unioned per protein.

        python NCBIgff2chrom.py annotation.gff(.gz)  > genome.chrom

2.  **Chrombase pipeline:** four flags. Filters proteins to chromosome-scale
    scaffolds (defined by the records in `-f`), writes three files:
        {prefix}.chrom          protein_id <tab> scaf <tab> strand <tab> start <tab> stop
        {prefix}.pep            input protein FASTA filtered to chr-scale proteins
        {prefix}.report.txt     one-line summary of what made it through.

        python NCBIgff2chrom.py \\
            -f genome.chr.fasta.gz -p protein.pep -g annotation.gff \\
            --union -o output_prefix

Notes
-----
- Only CDS lines with `protein_id=` are considered. Pseudogenes / ncRNAs are
  skipped, matching the behaviour of the original odp version.
- `--union` (the default in legacy mode; required in pipeline mode) takes the
  min CDS start and max CDS stop across all CDS lines for the same protein.
- Both GFF and FASTA inputs may be gzipped (detected by suffix).

License: GPL-3.0
"""

from __future__ import annotations

import argparse
import gzip
import sys
from typing import IO


def open_maybe_gz(path: str) -> IO[str]:
    """Open a file in text mode, transparently decompressing gzip suffixes."""
    if path.endswith((".gz", ".gzip", ".GZ", ".GZIP", ".gzipped", ".GZIPPED")):
        return gzip.open(path, "rt")
    return open(path, "r")


def fasta_scaffold_ids(path: str) -> set[str]:
    """Return the set of FASTA record IDs (first whitespace-separated token of each header)."""
    ids: set[str] = set()
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def parse_gff_cds(
    path: str,
    chr_scaffolds: set[str] | None,
    union: bool,
) -> tuple[dict[str, dict], int, int]:
    """
    Scan the GFF for CDS lines with a protein_id attribute.

    If chr_scaffolds is provided, only CDS lines on those scaffolds contribute to
    `proteins`; CDS lines on other scaffolds are counted into `nonchr_cds_count`
    and otherwise ignored.

    Returns (proteins, chr_cds_count, nonchr_cds_count).
    """
    proteins: dict[str, dict] = {}
    chr_cds_count = 0
    nonchr_cds_count = 0
    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parts[8]
            if "protein_id=" not in attrs:
                continue
            scaf = parts[0]
            if chr_scaffolds is not None:
                if scaf in chr_scaffolds:
                    chr_cds_count += 1
                else:
                    nonchr_cds_count += 1
                    continue
            else:
                chr_cds_count += 1
            pid = None
            for kv in attrs.split(";"):
                if kv.startswith("protein_id="):
                    pid = kv.split("=", 1)[1].strip()
                    break
            if not pid:
                continue
            try:
                start = int(parts[3])
                stop = int(parts[4])
            except ValueError:
                continue
            strand = parts[6]
            if pid not in proteins:
                proteins[pid] = {
                    "scaf": scaf,
                    "strand": strand,
                    "start": start,
                    "stop": stop,
                }
            elif union:
                p = proteins[pid]
                if start < p["start"]:
                    p["start"] = start
                if stop > p["stop"]:
                    p["stop"] = stop
    return proteins, chr_cds_count, nonchr_cds_count


def filter_protein_fasta(in_path: str, out_path: str, keep_ids: set[str]) -> tuple[int, int]:
    """Stream a protein FASTA, write a new one keeping only records whose ID is in keep_ids."""
    kept = 0
    total = 0
    keeping = False
    with open_maybe_gz(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                total += 1
                rid = line[1:].split()[0]
                keeping = rid in keep_ids
                if keeping:
                    kept += 1
                    fout.write(line)
            else:
                if keeping:
                    fout.write(line)
    return kept, total


def emit_chrom_lines(proteins: dict[str, dict], out_fh):
    for pid in sorted(proteins):
        p = proteins[pid]
        out_fh.write(f"{pid}\t{p['scaf']}\t{p['strand']}\t{p['start']}\t{p['stop']}\n")


def main_legacy(argv: list[str]) -> None:
    """Single positional GFF argument; chrom -> stdout. Backwards-compatible mode."""
    if len(argv) != 1:
        sys.exit("usage: NCBIgff2chrom.py <annotation.gff[.gz]>")
    proteins, _, _ = parse_gff_cds(argv[0], chr_scaffolds=None, union=True)
    emit_chrom_lines(proteins, sys.stdout)


def main_pipeline() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", "--fasta", required=True, help="Chromosome-scale FASTA (may be gzipped). Used to enumerate the chr-scale scaffold IDs to keep.")
    ap.add_argument("-p", "--protein", required=True, help="Full protein FASTA from the NCBI annotation (may be gzipped).")
    ap.add_argument("-g", "--gff", required=True, help="Full GFF from the NCBI annotation (may be gzipped).")
    ap.add_argument("-o", "--out-prefix", required=True, dest="prefix", help="Output prefix. Creates {prefix}.chrom, {prefix}.pep, {prefix}.report.txt")
    ap.add_argument("--union", action="store_true", help="Take the union of CDS coordinates per protein_id (min start, max stop). Strongly recommended; off by default to match legacy single-CDS behaviour.")
    args = ap.parse_args()

    chr_ids = fasta_scaffold_ids(args.fasta)
    if not chr_ids:
        sys.exit(f"error: no FASTA records found in {args.fasta}")

    proteins, chr_cds, nonchr_cds = parse_gff_cds(args.gff, chr_scaffolds=chr_ids, union=args.union)
    if not proteins:
        sys.exit("error: no CDS lines with protein_id found on the chromosome-scale scaffolds.")

    with open(args.prefix + ".chrom", "w") as fh:
        emit_chrom_lines(proteins, fh)

    kept, total = filter_protein_fasta(args.protein, args.prefix + ".pep", set(proteins))

    with open(args.prefix + ".report.txt", "w") as fh:
        fh.write("NCBIgff2chrom.py report\n")
        fh.write(f"  input genome FASTA       : {args.fasta}\n")
        fh.write(f"  input protein FASTA      : {args.protein}\n")
        fh.write(f"  input GFF                : {args.gff}\n")
        fh.write(f"  output prefix            : {args.prefix}\n")
        fh.write(f"  --union                  : {args.union}\n")
        fh.write("\n")
        fh.write(f"  chr-scale scaffolds      : {len(chr_ids)}\n")
        fh.write(f"  CDS lines on chr-scale   : {chr_cds}\n")
        fh.write(f"  CDS lines on other scaff : {nonchr_cds}\n")
        fh.write(f"  proteins on chr-scale    : {len(proteins)}\n")
        fh.write(f"  proteins kept in .pep    : {kept} of {total} in input\n")

    print(f"[NCBIgff2chrom] {len(proteins)} chr-scale proteins; wrote {args.prefix}.chrom, {args.prefix}.pep, {args.prefix}.report.txt", file=sys.stderr)


def main() -> None:
    # Dispatch: if the first non-program argv item starts with '-' (or argv has --help),
    # use the pipeline argparse parser. Otherwise fall back to the legacy single-arg form.
    rest = sys.argv[1:]
    if any(a in ("-h", "--help") for a in rest) or any(a.startswith("-") for a in rest):
        main_pipeline()
    else:
        main_legacy(rest)


if __name__ == "__main__":
    main()
