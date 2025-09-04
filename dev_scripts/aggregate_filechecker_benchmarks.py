#!/usr/bin/env python3
"""Aggregate benchmark files from odp_filechecker and append input file sizes.

This script reads Snakemake benchmark files produced by the rules in
`scripts/odp_filechecker`.  Each benchmark file name is expected to contain the
sample (genome assembly accession) and the rule name, e.g.:
```
benchmarks/check_genome_legality/GCF_00000000.1.check_genome_legality.benchmark.txt
```

For every benchmark entry the script looks up the input files for the
corresponding sample in the ODP configuration file and records their sizes in
bytes.  The benchmark metrics together with the computed file sizes are written
to a single TSV file.
"""
import argparse
import glob
import os
import re
from typing import Dict

import pandas as pd
import yaml


def _input_paths(rule: str, sample: str, config: Dict) -> Dict[str, str]:
    """Return a mapping of column name to input file path for a rule/sample."""
    sp_conf = config["species"][sample]
    if rule == "check_genome_legality":
        return {"genome_bytes": sp_conf["genome"]}
    if rule == "check_protein_legality":
        return {"proteins_bytes": sp_conf["proteins"]}
    if rule == "check_chrom_legality":
        return {
            "genome_bytes": sp_conf["genome"],
            "proteins_bytes": sp_conf["proteins"],
            "chrom_bytes": sp_conf["chrom"],
        }
    return {}


def aggregate(config_path: str, benchmarks_dir: str) -> pd.DataFrame:
    """Read benchmark files and append input file sizes."""
    with open(config_path) as fh:
        config = yaml.safe_load(fh)

    records = []
    pattern = re.compile(r"(?P<sample>[^.]+)\.(?P<rule>[^.]+)\.benchmark\.txt$")

    for bfile in glob.glob(os.path.join(benchmarks_dir, "**", "*.benchmark.txt"), recursive=True):
        m = pattern.search(os.path.basename(bfile))
        if not m:
            continue
        sample = m.group("sample")
        rule = m.group("rule")

        bench_df = pd.read_csv(bfile, sep="\t")
        if bench_df.empty:
            continue
        row = bench_df.iloc[0].to_dict()
        row.update({"sample": sample, "rule": rule})

        try:
            inputs = _input_paths(rule, sample, config)
        except KeyError:
            inputs = {}
        for col, path in inputs.items():
            row[col] = os.path.getsize(path) if os.path.exists(path) else pd.NA
        records.append(row)
    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Path to ODP config YAML")
    parser.add_argument(
        "--benchmarks", required=True, help="Directory with benchmark files"
    )
    parser.add_argument(
        "--out",
        default="aggregated_benchmarks_with_sizes.tsv",
        help="Output TSV file",
    )
    args = parser.parse_args()

    df = aggregate(args.config, args.benchmarks)
    df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
