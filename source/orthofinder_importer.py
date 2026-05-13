"""Orthofinder-result importer with synteny-aware multi-copy disambiguation.

This module turns a full Orthofinder run directory into an odp-style
`.rbh` table without requiring the user to pre-filter to single-copy
orthogroups. Multi-copy orthogroups, which would normally be discarded
by a naive `fgrep -f Orthogroups_SingleCopyOrthologues.txt` pipeline,
are rescued using one of several strategies — most usefully a
synteny-anchored strategy that picks the candidate gene whose
chromosome and position best fit the syntenic context established by
the single-copy backbone.

Strategies (selected via ``--multi-copy-strategy``):

- ``skip``: drop every multi-copy orthogroup. Matches the behaviour of
  naive Orthofinder → orthologs-table pipelines.
- ``longest``: for each species with N>1 candidate genes in an
  orthogroup, pick the longest CDS (``end - start``). No synteny
  information needed. Fast but doesn't account for paralog placement.
- ``most-common-chrom``: anchor-based, chromosome-only. For each
  species pair (sp_A, sp_B) build a weighted chromosome compatibility
  graph from the single-copy backbone, then for each multi-copy
  ambiguity in sp_A, pick the candidate whose chromosome has the
  highest weight to the chromosomes of the (single-copy) genes of the
  other species in the same orthogroup.
- ``synteny``: same as ``most-common-chrom`` but additionally rewards
  candidates whose position is close to single-copy anchors on the
  picked chromosome. Best signal, slower.

Multi-copy orthogroups whose best candidate score in any species falls
below ``--min-confidence`` are dropped (configurable; default 1).

The module also writes an optional audit report — one row per
(orthogroup, species) decision — recording the strategy used,
candidate count, picked gene, picked score, and the second-best
alternative.
"""
from __future__ import annotations

import gzip
import io
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Orthofinder run directory discovery
# ---------------------------------------------------------------------------


@dataclass
class OrthofinderPaths:
    """Resolved locations of useful files inside an Orthofinder run dir."""
    orthogroups_tsv: Path
    single_copy_list: Optional[Path] = None
    single_copy_sequences_dir: Optional[Path] = None
    hierarchical_orthogroups_dir: Optional[Path] = None
    resolved_gene_trees_dir: Optional[Path] = None


def discover_orthofinder_paths(root: Path | str) -> OrthofinderPaths:
    """Walk an Orthofinder run directory and resolve the useful subpaths.

    Accepts either:
    - the top-level Results_* directory (preferred), or
    - the ``Orthogroups/`` subdirectory directly.

    Raises ``FileNotFoundError`` if ``Orthogroups.tsv`` can't be found.
    """
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f"Orthofinder dir does not exist: {root}")

    # Find Orthogroups/Orthogroups.tsv. Try root/Orthogroups first; if root
    # IS the Orthogroups dir, look for the .tsv directly.
    candidates = [
        root / "Orthogroups" / "Orthogroups.tsv",
        root / "Orthogroups.tsv",
    ]
    orthogroups_tsv = next((p for p in candidates if p.is_file()), None)
    if orthogroups_tsv is None:
        raise FileNotFoundError(
            f"Could not find Orthogroups.tsv under {root}. Expected one of: "
            f"{[str(p) for p in candidates]}"
        )

    out = OrthofinderPaths(orthogroups_tsv=orthogroups_tsv)

    og_dir = orthogroups_tsv.parent
    sc_list = og_dir / "Orthogroups_SingleCopyOrthologues.txt"
    if sc_list.is_file():
        out.single_copy_list = sc_list

    sc_seq = root / "Single_Copy_Orthologue_Sequences"
    if sc_seq.is_dir():
        out.single_copy_sequences_dir = sc_seq

    hho = root / "Phylogenetic_Hierarchical_Orthogroups"
    if hho.is_dir():
        out.hierarchical_orthogroups_dir = hho

    rgt = root / "Resolved_Gene_Trees"
    if rgt.is_dir():
        out.resolved_gene_trees_dir = rgt

    return out


# ---------------------------------------------------------------------------
# Orthogroups.tsv parsing
# ---------------------------------------------------------------------------


def parse_orthogroups_tsv(path: Path | str) -> Tuple[List[str], pd.DataFrame]:
    """Parse Orthofinder ``Orthogroups.tsv``.

    The file has a header row: ``Orthogroup<TAB>sp1<TAB>sp2<...>``. Each
    data row is a single orthogroup. Each species cell may contain zero,
    one, or many comma-separated (Orthofinder default: ``, ``) gene IDs.

    Returns ``(species_names, dataframe)`` where the dataframe has one
    column per species (named after the species) of ``list[str]`` gene
    IDs, plus an ``OG`` column for the orthogroup ID.
    """
    path = Path(path)
    # Tolerate CRLF and gzip.
    if path.suffix == ".gz":
        opener = lambda p: io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8")
    else:
        opener = lambda p: open(p, "r", encoding="utf-8", newline="")

    rows: List[List[str]] = []
    with opener(path) as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            rows.append(line.split("\t"))
    if len(rows) < 2:
        raise ValueError(f"Orthogroups.tsv {path} has no data rows")

    header = rows[0]
    if header[0] != "Orthogroup":
        sys.stderr.write(
            f"orthofinder-import: warning: header[0] of {path} is "
            f"{header[0]!r}, expected 'Orthogroup'. Continuing.\n"
        )
    species = header[1:]
    if not species:
        raise ValueError(
            f"Orthogroups.tsv {path}: header has no species columns"
        )

    data_rows: List[Dict] = []
    for r in rows[1:]:
        if len(r) != len(header):
            # Pad or truncate: missing trailing cells = empty species.
            r = (r + [""] * len(header))[: len(header)]
        og_id = r[0]
        record: Dict = {"OG": og_id}
        for sp, cell in zip(species, r[1:]):
            cell = cell.strip()
            if not cell:
                genes: List[str] = []
            else:
                # Orthofinder's default delimiter is ", ". Tolerate
                # comma-only, tab-only, or whitespace-only separators.
                genes = [t.strip() for t in cell.replace(",", " ").split() if t.strip()]
            record[sp] = genes
        data_rows.append(record)

    return species, pd.DataFrame(data_rows)


def load_single_copy_list(path: Path | str) -> set[str]:
    """Read ``Orthogroups_SingleCopyOrthologues.txt`` (one OG id per line)."""
    out: set[str] = set()
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            v = raw.strip()
            if v:
                out.add(v)
    return out


def partition_single_vs_multi_copy(
    df: pd.DataFrame,
    species: Sequence[str],
    explicit_single_copy: Optional[set[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Split orthogroups into single-copy (exactly one gene per species,
    no zero cells) and everything-else.

    If ``explicit_single_copy`` is supplied (e.g. read from
    ``Orthogroups_SingleCopyOrthologues.txt``), that set defines the
    single-copy partition. Otherwise, the partition is computed from the
    cell contents.
    """
    if explicit_single_copy is not None:
        single_mask = df["OG"].isin(explicit_single_copy)
    else:
        # `.map` was renamed from `.applymap` in pandas 2.1; .applymap still
        # works on 2.x but warns. Use .apply(lambda col: col.map(...)) for
        # cross-version compatibility without triggering the deprecation.
        single_mask = df[species].apply(
            lambda col: col.map(lambda v: len(v) == 1)
        ).all(axis=1)
    single = df[single_mask].copy()
    multi = df[~single_mask].copy()
    return single, multi


# ---------------------------------------------------------------------------
# Anchor graph from single-copy backbone
# ---------------------------------------------------------------------------


@dataclass
class AnchorGraph:
    """Per-species-pair chromosome compatibility graph derived from the
    single-copy orthogroup backbone.

    ``edge[(sp_a, sp_b)][(chr_a, chr_b)] = count`` — the number of
    single-copy orthogroups where sp_a's gene is on ``chr_a`` and sp_b's
    gene is on ``chr_b``.

    ``positions[sp][chr] = sorted list of int positions of single-copy
    anchor genes on that chromosome``. Used for proximity scoring in the
    ``synteny`` strategy.
    """
    edge: Dict[Tuple[str, str], Dict[Tuple[str, str], int]]
    positions: Dict[str, Dict[str, List[int]]]


def build_anchor_graph(
    single_copy_df: pd.DataFrame,
    species: Sequence[str],
    coords: Dict[str, pd.DataFrame],
) -> AnchorGraph:
    """Compute the chromosome-pair edge weights and per-species anchor
    positions from the single-copy orthogroup backbone.

    ``coords[sp]`` must be a DataFrame with at least ``gene_id``,
    ``chrom``, and ``pos`` columns (produced by
    ``ortholog_table_importer.parse_coordinates``).
    """
    edge: Dict[Tuple[str, str], Dict[Tuple[str, str], int]] = {}
    positions: Dict[str, Dict[str, List[int]]] = {sp: {} for sp in species}

    # Pre-build per-species gene -> (chrom, pos) lookups.
    lookup: Dict[str, Dict[str, Tuple[str, int]]] = {}
    for sp in species:
        df = coords[sp]
        lookup[sp] = dict(zip(df["gene_id"], zip(df["chrom"], df["pos"])))

    for _, row in single_copy_df.iterrows():
        gene_per_sp: Dict[str, Tuple[str, int]] = {}
        skip = False
        for sp in species:
            genes = row[sp]
            if len(genes) != 1:
                # Defensive: should never happen if the partition is right.
                skip = True
                break
            info = lookup[sp].get(genes[0])
            if info is None:
                # Gene not in the coord file — drop this row from the
                # anchor backbone. It can't anchor anything.
                skip = True
                break
            gene_per_sp[sp] = info
        if skip:
            continue

        for sp, (chrom, pos) in gene_per_sp.items():
            positions[sp].setdefault(chrom, []).append(pos)

        # Pairwise edges (directed both ways).
        species_list = list(gene_per_sp.keys())
        for i, sp_a in enumerate(species_list):
            chr_a = gene_per_sp[sp_a][0]
            for sp_b in species_list[i + 1:]:
                chr_b = gene_per_sp[sp_b][0]
                key = (sp_a, sp_b)
                rev_key = (sp_b, sp_a)
                edge.setdefault(key, {})
                edge.setdefault(rev_key, {})
                edge[key][(chr_a, chr_b)] = edge[key].get((chr_a, chr_b), 0) + 1
                edge[rev_key][(chr_b, chr_a)] = edge[rev_key].get((chr_b, chr_a), 0) + 1

    for sp in species:
        for chrom in positions[sp]:
            positions[sp][chrom].sort()

    return AnchorGraph(edge=edge, positions=positions)


# ---------------------------------------------------------------------------
# Disambiguation strategies
# ---------------------------------------------------------------------------


@dataclass
class DisambDecision:
    """One species-level pick within an orthogroup, for audit purposes."""
    og: str
    species: str
    strategy: str
    n_candidates: int
    picked_gene: Optional[str]
    picked_score: float
    second_best_gene: Optional[str]
    second_best_score: Optional[float]
    reason: str


def _nearest_anchor_distance(positions: List[int], pos: int) -> int:
    """Return the distance (in bp) from ``pos`` to the nearest entry in
    the sorted ``positions`` list. Empty list returns a sentinel (10 Mbp)."""
    if not positions:
        return 10_000_000
    import bisect
    idx = bisect.bisect_left(positions, pos)
    candidates = []
    if idx > 0:
        candidates.append(abs(pos - positions[idx - 1]))
    if idx < len(positions):
        candidates.append(abs(pos - positions[idx]))
    return min(candidates) if candidates else 10_000_000


def _candidates_for(
    species: str,
    genes: List[str],
    coords_lookup: Dict[str, Dict[str, Tuple[str, int, int, int]]],
) -> List[Tuple[str, str, int, int]]:
    """Return list of ``(gene_id, chrom, pos, length)`` for the given
    gene IDs, dropping any that aren't in the coord file."""
    out: List[Tuple[str, str, int, int]] = []
    sp_lookup = coords_lookup[species]
    for g in genes:
        info = sp_lookup.get(g)
        if info is not None:
            chrom, pos, start, end = info
            out.append((g, chrom, pos, end - start))
    return out


def disambiguate_multicopy(
    multi_copy_df: pd.DataFrame,
    species: Sequence[str],
    coords: Dict[str, pd.DataFrame],
    anchor_graph: Optional[AnchorGraph],
    strategy: str,
    min_confidence: float = 1.0,
) -> Tuple[pd.DataFrame, List[DisambDecision]]:
    """Apply ``strategy`` to each multi-copy orthogroup. Returns

      (resolved_single_copy_df, audit_decisions)

    where ``resolved_single_copy_df`` has the same shape as the
    single-copy partition (one gene per species per row, as a 1-element
    list in each species column) for orthogroups that survive, and
    ``audit_decisions`` records every decision for the report.

    For ``strategy == 'skip'``, ``resolved_single_copy_df`` is empty.
    For ``strategy == 'longest'``, no anchor graph is needed; pass None.
    For ``most-common-chrom`` and ``synteny``, ``anchor_graph`` is
    required and must be built from the single-copy backbone.
    """
    if strategy == "skip":
        return pd.DataFrame(columns=multi_copy_df.columns), []

    valid = {"longest", "most-common-chrom", "synteny"}
    if strategy not in valid:
        raise ValueError(
            f"unknown multi-copy strategy: {strategy!r}. Choose one of "
            f"{sorted({'skip'} | valid)}."
        )

    if strategy in ("most-common-chrom", "synteny") and anchor_graph is None:
        raise ValueError(
            f"strategy {strategy!r} needs an anchor graph built from "
            f"single-copy orthogroups."
        )

    # gene_id -> (chrom, pos, start, end) per species, for length + position.
    coords_lookup: Dict[str, Dict[str, Tuple[str, int, int, int]]] = {}
    for sp in species:
        df = coords[sp]
        coords_lookup[sp] = {
            row["gene_id"]: (row["chrom"], row["pos"], row["start"], row["end"])
            for _, row in df.iterrows()
        }

    resolved_rows: List[Dict] = []
    audit: List[DisambDecision] = []

    for _, og_row in multi_copy_df.iterrows():
        og = og_row["OG"]
        picked: Dict[str, str] = {}
        og_kept = True

        # First pass: species that already have exactly one candidate need
        # no disambiguation — record the single gene and move on. These
        # rows serve as "given" context when scoring multi-copy species
        # in the same OG.
        species_to_disambiguate: List[str] = []
        for sp in species:
            genes = og_row[sp]
            if len(genes) == 0:
                # Missing species — this OG can't be a complete single-
                # copy row across all species. Skip.
                og_kept = False
                audit.append(DisambDecision(
                    og=og, species=sp, strategy=strategy,
                    n_candidates=0, picked_gene=None, picked_score=0.0,
                    second_best_gene=None, second_best_score=None,
                    reason="no candidates in this species",
                ))
                break
            if len(genes) == 1:
                picked[sp] = genes[0]
            else:
                species_to_disambiguate.append(sp)
        if not og_kept:
            continue

        # Second pass: for each multi-copy species, score candidates.
        for sp in species_to_disambiguate:
            cands = _candidates_for(sp, og_row[sp], coords_lookup)
            if not cands:
                audit.append(DisambDecision(
                    og=og, species=sp, strategy=strategy,
                    n_candidates=len(og_row[sp]),
                    picked_gene=None, picked_score=0.0,
                    second_best_gene=None, second_best_score=None,
                    reason="no candidate gene found in coord file",
                ))
                og_kept = False
                break

            scored: List[Tuple[float, str, str]] = []  # (score, gene, reason)

            if strategy == "longest":
                for gene, chrom, pos, length in cands:
                    scored.append((float(length), gene, f"length={length}"))
            else:
                # most-common-chrom or synteny: use anchor graph.
                other_sp_genes: Dict[str, Tuple[str, int]] = {}
                for other_sp in species:
                    if other_sp == sp:
                        continue
                    if other_sp not in picked:
                        # Another multi-copy species we haven't resolved
                        # yet; skip — only score against already-resolved
                        # species (single-copy + earlier picks).
                        continue
                    info = coords_lookup[other_sp].get(picked[other_sp])
                    if info is not None:
                        other_sp_genes[other_sp] = (info[0], info[1])

                for gene, chrom, pos, length in cands:
                    chrom_score = 0.0
                    for other_sp, (other_chr, _other_pos) in other_sp_genes.items():
                        key = (sp, other_sp)
                        weight = anchor_graph.edge.get(key, {}).get((chrom, other_chr), 0)
                        chrom_score += weight

                    if strategy == "synteny":
                        anchor_positions = anchor_graph.positions.get(sp, {}).get(chrom, [])
                        nearest = _nearest_anchor_distance(anchor_positions, pos)
                        # Decay: 1 / (1 + distance_in_Mbp). Anchors within
                        # 1 Mbp give a boost ~0.5; 10 Mbp gives ~0.1; 100
                        # Mbp gives ~0.01. Bounded contribution so a single
                        # chromosome-pair vote always beats positional
                        # tiebreaking.
                        proximity = 1.0 / (1.0 + nearest / 1_000_000.0)
                        score = chrom_score + 0.5 * proximity
                        reason = f"chrom={chrom_score:.0f} prox={proximity:.3f}"
                    else:  # most-common-chrom
                        score = chrom_score
                        reason = f"chrom={chrom_score:.0f}"

                    scored.append((score, gene, reason))

            # Pick highest. Tie-break by length (descending), then alphabetic.
            scored.sort(key=lambda t: (-t[0], -coords_lookup[sp][t[1]][3] + coords_lookup[sp][t[1]][2], t[1]))
            best_score, best_gene, best_reason = scored[0]
            second_gene = scored[1][1] if len(scored) > 1 else None
            second_score = scored[1][0] if len(scored) > 1 else None

            audit.append(DisambDecision(
                og=og, species=sp, strategy=strategy,
                n_candidates=len(cands),
                picked_gene=best_gene, picked_score=best_score,
                second_best_gene=second_gene, second_best_score=second_score,
                reason=best_reason,
            ))

            if best_score < min_confidence:
                og_kept = False
                # Mark in the audit reason for clarity.
                audit[-1] = DisambDecision(
                    og=og, species=sp, strategy=strategy,
                    n_candidates=len(cands),
                    picked_gene=best_gene, picked_score=best_score,
                    second_best_gene=second_gene, second_best_score=second_score,
                    reason=audit[-1].reason + f"; BELOW min-confidence {min_confidence}",
                )
                break

            picked[sp] = best_gene

        if og_kept:
            row: Dict = {"OG": og}
            for sp in species:
                row[sp] = [picked[sp]]
            resolved_rows.append(row)

    if resolved_rows:
        resolved = pd.DataFrame(resolved_rows)
    else:
        resolved = pd.DataFrame(columns=multi_copy_df.columns)
    return resolved, audit


# ---------------------------------------------------------------------------
# Emit .rbh
# ---------------------------------------------------------------------------


def emit_rbh_from_orthogroups(
    df: pd.DataFrame,
    species: Sequence[str],
    coords: Dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Convert a single-copy-per-row orthogroup DataFrame (with list-valued
    species columns of length 1) into an odp-style .rbh DataFrame."""
    # Pre-build coord lookups.
    lookup: Dict[str, Dict[str, Tuple[str, int]]] = {}
    for sp in species:
        cdf = coords[sp]
        lookup[sp] = dict(zip(cdf["gene_id"], zip(cdf["chrom"], cdf["pos"])))

    out_rows: List[Dict] = []
    skipped = 0
    for _, row in df.iterrows():
        new_row: Dict = {"rbh": "", "gene_group": "None"}
        skip = False
        for sp in species:
            gene_list = row[sp]
            gene = gene_list[0] if gene_list else None
            if gene is None:
                skip = True
                break
            info = lookup[sp].get(gene)
            if info is None:
                skip = True
                break
            chrom, pos = info
            new_row[f"{sp}_gene"] = gene
            new_row[f"{sp}_scaf"] = chrom
            new_row[f"{sp}_pos"] = pos
        if skip:
            skipped += 1
            continue
        out_rows.append(new_row)

    if skipped:
        sys.stderr.write(
            f"orthofinder-import: emit_rbh dropped {skipped} OGs with "
            f"genes missing from coord files.\n"
        )
    if not out_rows:
        return pd.DataFrame(columns=["rbh", "gene_group"] + [
            f"{sp}_{kind}" for sp in species for kind in ("gene", "scaf", "pos")
        ])

    for i, r in enumerate(out_rows, start=1):
        r["rbh"] = f"rbh{i}"

    col_order = ["rbh", "gene_group"]
    for sp in species:
        col_order += [f"{sp}_gene", f"{sp}_scaf", f"{sp}_pos"]
    return pd.DataFrame(out_rows)[col_order]


# ---------------------------------------------------------------------------
# Audit report
# ---------------------------------------------------------------------------


def write_audit_report(decisions: List[DisambDecision], path: Path | str) -> int:
    """Write the disambiguation audit as a TSV. Returns row count."""
    if not decisions:
        Path(path).write_text(
            "OG\tspecies\tstrategy\tn_candidates\tpicked_gene\t"
            "picked_score\tsecond_best_gene\tsecond_best_score\treason\n"
        )
        return 0
    rows = [
        {
            "OG": d.og,
            "species": d.species,
            "strategy": d.strategy,
            "n_candidates": d.n_candidates,
            "picked_gene": d.picked_gene if d.picked_gene is not None else "",
            "picked_score": f"{d.picked_score:.4f}",
            "second_best_gene": d.second_best_gene if d.second_best_gene is not None else "",
            "second_best_score": f"{d.second_best_score:.4f}" if d.second_best_score is not None else "",
            "reason": d.reason,
        }
        for d in decisions
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return len(rows)


# ---------------------------------------------------------------------------
# Top-level entry
# ---------------------------------------------------------------------------


def orthofinder_to_rbh(
    orthofinder_dir: Path | str,
    coord_paths: Dict[str, Path | str],
    strategy: str = "synteny",
    min_confidence: float = 1.0,
) -> Tuple[pd.DataFrame, List[DisambDecision]]:
    """Top-level: walk the Orthofinder dir, parse Orthogroups.tsv, build
    the anchor graph from single-copy orthogroups, disambiguate multi-
    copy orthogroups per the chosen strategy, and emit a unified .rbh
    DataFrame along with an audit list.

    Args:
      orthofinder_dir: top-level Orthofinder run directory.
      coord_paths: ``{species: path}``. Species names must match the
        column names in ``Orthogroups.tsv``. Each path may point at a
        .chrom or BED file (gzip-OK; format auto-detected).
      strategy: ``'skip'``, ``'longest'``, ``'most-common-chrom'``, or
        ``'synteny'``.
      min_confidence: orthogroups whose best per-species score in
        any species falls below this threshold are dropped (only
        meaningful for ``most-common-chrom`` and ``synteny``).
    """
    # Lazy import to avoid coupling at module import time.
    from ortholog_table_importer import parse_coordinates

    paths = discover_orthofinder_paths(orthofinder_dir)
    species, og_df = parse_orthogroups_tsv(paths.orthogroups_tsv)

    missing = [s for s in species if s not in coord_paths]
    if missing:
        raise ValueError(
            f"Coord file missing for species {missing!r}. The Orthofinder "
            f"Orthogroups.tsv has columns {species!r}; supply --chrom or "
            f"--bed for every one."
        )

    coords: Dict[str, pd.DataFrame] = {}
    for sp in species:
        coords[sp] = parse_coordinates(coord_paths[sp])

    # Partition into single-copy + multi-copy.
    explicit_sc = None
    if paths.single_copy_list is not None:
        explicit_sc = load_single_copy_list(paths.single_copy_list)
    single_df, multi_df = partition_single_vs_multi_copy(
        og_df, species, explicit_single_copy=explicit_sc
    )

    sys.stderr.write(
        f"orthofinder-import: {len(single_df)} single-copy OGs, "
        f"{len(multi_df)} multi-copy OGs. Strategy: {strategy}.\n"
    )

    if strategy == "skip":
        anchor_graph = None
        resolved_df = pd.DataFrame(columns=multi_df.columns)
        audit: List[DisambDecision] = []
    elif strategy == "longest":
        anchor_graph = None
        resolved_df, audit = disambiguate_multicopy(
            multi_df, species, coords, None, strategy, min_confidence
        )
    else:
        anchor_graph = build_anchor_graph(single_df, species, coords)
        resolved_df, audit = disambiguate_multicopy(
            multi_df, species, coords, anchor_graph, strategy, min_confidence
        )

    sys.stderr.write(
        f"orthofinder-import: rescued {len(resolved_df)}/{len(multi_df)} "
        f"multi-copy OGs.\n"
    )

    # Combine single + resolved-multi into one frame, then emit .rbh.
    combined = pd.concat([single_df, resolved_df], ignore_index=True)
    rbh = emit_rbh_from_orthogroups(combined, species, coords)
    return rbh, audit
