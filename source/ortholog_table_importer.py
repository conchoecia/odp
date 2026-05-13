"""Generic ortholog-table importers and friends.

Provides paths into odp's `.rbh`-based pipeline that don't require the user
to run the full diamond / blastp reciprocal-best-hits step. Two input
shapes are supported:

* **Orthofinder output.** Either the raw `Orthogroups.tsv` filtered to
  single-copy orthologs, or any tab-separated file with one column per
  species and one row per orthogroup.
* **Plain N-column orthologs table.** Tab-separated, no orthogroup ID
  column. One column per species, one row per ortholog tuple.

Both shapes are joined against per-species BED files
(``chrom``/``start``/``end``/``gene_id`` minimum) and emitted as an
odp-compatible `.rbh` table with the columns:

* ``rbh`` — synthetic identifier, ``rbhN`` where N is the row index
* ``gene_group`` — ``"None"`` by default (downstream tools fill it in)
* ``<species>_gene`` — gene identifier for each species
* ``<species>_scaf`` — chromosome / scaffold name for each species
* ``<species>_pos`` — midpoint of the gene's BED coordinates (``(start+end)//2``)

This lets users with pre-computed orthology (Orthofinder runs, in-house
RBH workflows, manually curated ortholog tables) skip the all-vs-all
search and feed orthologs straight into ``odp run`` / ``odp nway-rbh``
/ downstream plotters.

No new dependencies. Standard library + pandas (already a runtime dep).
"""
from __future__ import annotations

import csv
import gzip
import io
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Coordinate file parsing (.chrom and BED)
# ---------------------------------------------------------------------------


def _open_text_maybe_gzip(path: Path | str) -> io.TextIOBase:
    """Open a path as text. Handles gzip by magic bytes (not extension)."""
    path = Path(path)
    with open(path, "rb") as fh:
        head = fh.read(2)
    if head == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8", newline="")


def detect_coord_format(path: Path | str) -> str:
    """Return ``'chrom'``, ``'bed'``, or raise. Inspects the first non-empty
    line:

    - ``.chrom``: 5 columns; col 3 is ``+``/``-``/``.``; cols 4-5 are int.
    - BED:        4+ columns; cols 2-3 are int; col 1 is the chromosome.
    """
    with _open_text_maybe_gzip(path) as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            # chrom format: col[2] in {'+','-','.'} and cols 3, 4 are int
            if (
                len(parts) >= 5
                and parts[2] in ("+", "-", ".")
                and parts[3].lstrip("-").isdigit()
                and parts[4].lstrip("-").isdigit()
            ):
                return "chrom"
            # BED format: cols 1, 2 are int
            if (
                len(parts) >= 4
                and parts[1].lstrip("-").isdigit()
                and parts[2].lstrip("-").isdigit()
            ):
                return "bed"
            raise ValueError(
                f"Coordinate file {path} line 1: cannot detect format. "
                f"Expected either odp .chrom (gene_id, chrom, +/-/., start, "
                f"stop) or BED (chrom, start, end, gene_id). Got: {parts!r}"
            )
    raise ValueError(f"Coordinate file {path} is empty")


def parse_chrom(path: Path | str) -> pd.DataFrame:
    """Read an odp ``.chrom`` file (5 tab-separated fields: ``gene_id``,
    ``scaffold``, ``strand``, ``start``, ``stop``). Transparently handles
    gzip-compressed input (``.chrom.gz``).

    Returns a DataFrame with columns: ``chrom`` (str), ``start`` (int),
    ``end`` (int), ``gene_id`` (str), ``strand`` (str), ``pos`` (int —
    gene midpoint).
    """
    path = Path(path)
    rows: list[dict] = []
    with _open_text_maybe_gzip(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                raise ValueError(
                    f".chrom file {path} line {lineno}: needs at least 5 "
                    f"tab-separated columns (gene_id, scaffold, strand, "
                    f"start, stop), got {len(parts)}: {parts!r}"
                )
            gene_id, chrom, strand, start_s, stop_s = parts[:5]
            if strand not in ("+", "-", "."):
                raise ValueError(
                    f".chrom file {path} line {lineno}: strand must be "
                    f"'+', '-', or '.', got {strand!r}"
                )
            try:
                start = int(start_s)
                end = int(stop_s)
            except ValueError as e:
                raise ValueError(
                    f".chrom file {path} line {lineno}: start/stop must be "
                    f"integers, got {start_s!r}/{stop_s!r}"
                ) from e
            rows.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "gene_id": gene_id,
                "strand": strand,
                "pos": (start + end) // 2,
            })
    if not rows:
        raise ValueError(f".chrom file {path} is empty or has no parseable lines")
    return pd.DataFrame(rows)


def parse_bed(path: Path | str) -> pd.DataFrame:
    """Read a BED-like file with at least four columns
    (``chrom``, ``start``, ``end``, ``gene_id``). Trailing columns are kept
    but unused. Strips CRLF and skips blank lines. Transparently handles
    gzip-compressed input.

    Returns a DataFrame with columns: ``chrom`` (str), ``start`` (int),
    ``end`` (int), ``gene_id`` (str), ``pos`` (int — gene midpoint).
    """
    path = Path(path)
    rows: list[dict] = []
    with _open_text_maybe_gzip(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                raise ValueError(
                    f"BED file {path} line {lineno}: needs at least 4 tab-"
                    f"separated columns (chrom, start, end, gene_id), got "
                    f"{len(parts)}: {parts!r}"
                )
            chrom, start_s, end_s, gene_id = parts[0], parts[1], parts[2], parts[3]
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError as e:
                raise ValueError(
                    f"BED file {path} line {lineno}: start/end must be "
                    f"integers, got {start_s!r}/{end_s!r}"
                ) from e
            rows.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "gene_id": gene_id,
                "pos": (start + end) // 2,
            })
    if not rows:
        raise ValueError(f"BED file {path} is empty or has no parseable lines")
    return pd.DataFrame(rows)


def parse_coordinates(path: Path | str, fmt: str = "auto") -> pd.DataFrame:
    """Parse a coordinate file in either odp .chrom or BED format.

    Args:
      path: file path (gzip-compressed accepted).
      fmt: ``'auto'`` (default), ``'chrom'``, or ``'bed'``. With ``'auto'``,
        the format is detected from the first non-empty line.

    Returns a DataFrame with columns ``chrom``, ``start``, ``end``,
    ``gene_id``, ``pos`` (and ``strand`` if the input was .chrom).
    """
    if fmt == "auto":
        fmt = detect_coord_format(path)
    if fmt == "chrom":
        return parse_chrom(path)
    if fmt == "bed":
        return parse_bed(path)
    raise ValueError(f"unsupported coordinate format: {fmt!r}")


# ---------------------------------------------------------------------------
# Orthologs-table parsing
# ---------------------------------------------------------------------------


def _read_tsv_skipping_blanks(path: Path) -> List[List[str]]:
    """Read a TSV with tolerance for CRLF line endings and blank lines."""
    with open(path, "rb") as fh:
        raw = fh.read()
    # Normalise line endings: Orthofinder output is "ASCII with CRLF
    # terminators" on common installs. csv.reader handles \r\n but not
    # bare \r, so we standardise to \n up front.
    text = raw.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    out: list[list[str]] = []
    for line in text.split("\n"):
        if not line.strip():
            continue
        out.append(line.split("\t"))
    return out


def parse_ortholog_table(
    path: Path | str,
    species_names: Optional[Sequence[str]] = None,
    has_orthogroup_id_column: Optional[bool] = None,
    has_header: bool = False,
) -> Tuple[List[str], pd.DataFrame]:
    """Parse an N-column ortholog table.

    The table is tab-separated with one row per ortholog tuple and one
    column per species (optionally preceded by an orthogroup-id column,
    which is the shape of Orthofinder's ``Orthogroups.tsv``).

    Args:
      path: input path.
      species_names: optional list of species names matching the species
        columns left-to-right. If omitted and the file has a header row,
        the header is used.
      has_orthogroup_id_column: True if the first column is an orthogroup
        identifier rather than a gene. If None, the function auto-detects:
        the column is treated as an orthogroup id if every value matches
        ``^[A-Za-z]+\\d+$`` (e.g. ``OG0000001``).
      has_header: True if the first row is a header.

    Returns:
      ``(species_names, dataframe)`` where the dataframe has one column
      per species named after that species, each cell being a gene id.

    Single-copy orthologs in Orthofinder put exactly one gene per species
    per row. Multi-copy rows (rare but possible if the user feeds the
    full ``Orthogroups.tsv`` instead of the single-copy subset) are
    rejected with a clear error message — the .rbh format assumes one
    gene per species per row.
    """
    path = Path(path)
    rows = _read_tsv_skipping_blanks(path)
    if not rows:
        raise ValueError(f"Orthologs table {path} is empty")

    if has_header:
        header_row = rows[0]
        data_rows = rows[1:]
    else:
        header_row = None
        data_rows = rows

    if not data_rows:
        raise ValueError(f"Orthologs table {path} has no data rows")

    ncols = len(data_rows[0])
    if any(len(r) != ncols for r in data_rows):
        bad = next(i for i, r in enumerate(data_rows, start=1) if len(r) != ncols)
        raise ValueError(
            f"Orthologs table {path}: inconsistent column count. Row {bad} "
            f"has {len(data_rows[bad-1])} columns, first row has {ncols}."
        )

    # Auto-detect orthogroup-id column.
    #
    # Pattern matches real Orthofinder-style identifiers (OG0000001) and
    # similar — at least two letters followed by at least four digits.
    # A bare "g1" is NOT an orthogroup id and should be treated as a gene.
    #
    # If `species_names` is supplied explicitly and the column count is
    # exactly one greater than the species count, we infer the leading
    # column is an OG identifier even if it doesn't match the regex
    # (allows custom-format OG ids).
    if has_orthogroup_id_column is None:
        if species_names is not None:
            # User told us the species; the column count vs. species count
            # disambiguates whether there's a leading OG column.
            if ncols == len(species_names):
                has_orthogroup_id_column = False
            elif ncols == len(species_names) + 1:
                has_orthogroup_id_column = True
            else:
                raise ValueError(
                    f"Orthologs table {path}: --species-order has "
                    f"{len(species_names)} entries but the table has "
                    f"{ncols} columns. Expected either {len(species_names)} "
                    f"(no OG id column) or {len(species_names) + 1} (one "
                    f"OG id column followed by species)."
                )
        else:
            import re
            og_pattern = re.compile(r"^[A-Za-z]{2,}\d{4,}$")
            first_col = [r[0] for r in data_rows]
            has_orthogroup_id_column = all(og_pattern.match(v) for v in first_col)

    if has_orthogroup_id_column:
        species_cols = list(range(1, ncols))
    else:
        species_cols = list(range(ncols))

    # Determine species names.
    if species_names is not None:
        if len(species_names) != len(species_cols):
            raise ValueError(
                f"Orthologs table {path}: --species-order has "
                f"{len(species_names)} entries but the table has "
                f"{len(species_cols)} species columns."
            )
        names = list(species_names)
    elif header_row is not None:
        names = [header_row[i] for i in species_cols]
    else:
        names = [f"sp{i+1}" for i in range(len(species_cols))]

    # Build the dataframe.
    out_rows: list[dict] = []
    for r_idx, r in enumerate(data_rows, start=1):
        row_dict: dict = {}
        for col_idx, sp in zip(species_cols, names):
            cell = r[col_idx].strip()
            if not cell:
                raise ValueError(
                    f"Orthologs table {path}: row {r_idx} has an empty cell "
                    f"for species {sp!r}. Orthologs tables fed to "
                    f"orthologs-to-rbh must have one gene per species per "
                    f"row. Pre-filter to the single-copy subset (e.g. via "
                    f"Orthogroups_SingleCopyOrthologues.txt for Orthofinder)."
                )
            # Multi-copy detection: comma- or space-separated gene lists.
            if "," in cell or (" " in cell and "_" not in cell.split(" ")[0]):
                # Heuristic: if the cell looks like "gene1, gene2" or
                # "gene1 gene2", reject. We allow a single gene id that
                # happens to contain spaces (rare but legal).
                tokens = [t for t in cell.replace(",", " ").split() if t]
                if len(tokens) > 1:
                    raise ValueError(
                        f"Orthologs table {path}: row {r_idx} species {sp!r} "
                        f"has multiple genes ({len(tokens)}): {cell!r}. "
                        f".rbh format requires one gene per species per "
                        f"row. Filter to single-copy orthologs first."
                    )
            row_dict[sp] = cell
        out_rows.append(row_dict)

    return names, pd.DataFrame(out_rows)


# ---------------------------------------------------------------------------
# Join + .rbh emit
# ---------------------------------------------------------------------------


def orthologs_to_rbh(
    orthologs_path: Path | str,
    bed_paths: Dict[str, Path | str],
    species_order: Optional[Sequence[str]] = None,
    has_header: bool = False,
    has_orthogroup_id_column: Optional[bool] = None,
) -> pd.DataFrame:
    """Join an N-column ortholog table against per-species coordinate
    files and return an odp-style `.rbh` DataFrame.

    Args:
      orthologs_path: path to the tab-separated orthologs table.
      bed_paths: mapping of species name → path to that species'
        coordinate file. odp ``.chrom`` and BED formats are both
        accepted; the format is auto-detected from file content.
        gzip-compressed inputs are handled transparently.
      species_order: if supplied, the species columns of the orthologs
        table are interpreted in this order (overrides any header row
        and disables auto-detection of an orthogroup-id column).
      has_header: True if the orthologs table has a header row.
      has_orthogroup_id_column: see ``parse_ortholog_table``.

    Returns:
      A DataFrame with columns
      ``rbh``, ``gene_group``, then for each species
      ``<sp>_gene``, ``<sp>_scaf``, ``<sp>_pos``.
    """
    species_names, ortho_df = parse_ortholog_table(
        orthologs_path,
        species_names=species_order,
        has_orthogroup_id_column=has_orthogroup_id_column,
        has_header=has_header,
    )

    # Every species in the table must have a coordinate file.
    missing_beds = [s for s in species_names if s not in bed_paths]
    if missing_beds:
        raise ValueError(
            f"Coordinate files missing for species {missing_beds!r}. The "
            f"orthologs table has columns {species_names!r}; --bed or "
            f"--chrom must be supplied for every one."
        )

    # Load coordinate files (auto-detect chrom vs BED) and build
    # per-species gene-id → (chrom, pos) lookup.
    sp_lookup: Dict[str, Dict[str, Tuple[str, int]]] = {}
    for sp in species_names:
        bed_df = parse_coordinates(bed_paths[sp])
        # If a gene id is duplicated in the BED, keep the first; warn.
        dup_mask = bed_df["gene_id"].duplicated()
        if dup_mask.any():
            n_dup = dup_mask.sum()
            sys.stderr.write(
                f"orthologs-to-rbh: warning: BED for {sp!r} has {n_dup} "
                f"duplicate gene_id entries; keeping first occurrence of "
                f"each.\n"
            )
            bed_df = bed_df[~dup_mask].copy()
        sp_lookup[sp] = dict(zip(bed_df["gene_id"], zip(bed_df["chrom"], bed_df["pos"])))

    # Build the .rbh rows.
    out_rows: list[dict] = []
    n_skipped = 0
    for idx, row in ortho_df.iterrows():
        new_row: dict = {"rbh": "", "gene_group": "None"}
        skip = False
        for sp in species_names:
            gene_id = row[sp]
            info = sp_lookup[sp].get(gene_id)
            if info is None:
                skip = True
                break
            chrom, pos = info
            new_row[f"{sp}_gene"] = gene_id
            new_row[f"{sp}_scaf"] = chrom
            new_row[f"{sp}_pos"] = pos
        if skip:
            n_skipped += 1
            continue
        out_rows.append(new_row)

    if n_skipped:
        sys.stderr.write(
            f"orthologs-to-rbh: skipped {n_skipped} ortholog rows that "
            f"reference gene ids missing from the BED files.\n"
        )

    if not out_rows:
        raise ValueError(
            "orthologs-to-rbh: zero ortholog rows survived the join with "
            "the BED files. Check that gene ids match between the orthologs "
            "table and the BED's fourth column."
        )

    # Assign sequential rbh ids.
    for i, r in enumerate(out_rows, start=1):
        r["rbh"] = f"rbh{i}"

    df = pd.DataFrame(out_rows)
    # Column order: rbh, gene_group, then triples per species.
    col_order = ["rbh", "gene_group"]
    for sp in species_names:
        col_order += [f"{sp}_gene", f"{sp}_scaf", f"{sp}_pos"]
    return df[col_order]


def write_rbh(df: pd.DataFrame, out_path: Path | str) -> int:
    """Write an `.rbh` DataFrame to a tab-separated file. Returns the
    number of rows written."""
    out_path = Path(out_path)
    df.to_csv(out_path, sep="\t", index=False)
    return len(df)
