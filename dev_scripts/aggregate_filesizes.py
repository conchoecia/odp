#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path

import yaml  # pip install pyyaml

def bytes_to_mb(nbytes: int, base: int = 1024) -> float:
    return nbytes / (base * base)

def size_mb_or_na(path: str | None, base: int = 1024) -> str:
    if not path:
        return "NA"
    p = Path(path)
    try:
        mb = bytes_to_mb(p.stat().st_size, base)
        return f"{mb:.3f}"
    except FileNotFoundError:
        print(f"WARNING: file not found -> {p}", file=sys.stderr)
        return "NA"

def main():
    ap = argparse.ArgumentParser(description="Summarize input file sizes from config.yaml")
    ap.add_argument("-c", "--config", help="Path to config.yaml")
    ap.add_argument("-o", "--out", help="Output TSV path",
                    type=str, default="input_filesizes.tsv")
    ap.add_argument("--base", type=int, choices=(1000, 1024), default=1024,
                    help="MB base (1000 for decimal MB, 1024 for MiB; default 1024)")
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    species = cfg.get("species", {})

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "assembly_accession", "proteins_MB", "chrom_MB", "genome_MB"])
        for sample, info in species.items():
            acc = info.get("assembly_accession", "")
            proteins_mb = size_mb_or_na(info.get("proteins"), args.base)
            chrom_mb    = size_mb_or_na(info.get("chrom"), args.base)
            genome_mb   = size_mb_or_na(info.get("genome"), args.base)
            w.writerow([sample, acc, proteins_mb, chrom_mb, genome_mb])

if __name__ == "__main__":
    main()
