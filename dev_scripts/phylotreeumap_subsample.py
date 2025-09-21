#!/usr/bin/env python3
from ete4 import NCBITaxa
from collections import defaultdict
import pandas as pd
from pathlib import Path
import random
import sys

def return_kingdom_full_sort_order():
    """Return a list of the sort order for the taxonomic rankings."""
    return ["superkingdom",
            "kingdom",
            "subkingdom",
            "infrakingdom",
            "superphylum",
            "phylum",
            "subphylum",
            "infraphylum",
            "superclass",
            "class",
            "subclass",
            "infraclass",
            "parvclass",
            "cohort",
            "subcohort",
            "superorder",
            "order",
            "suborder",
            "infraorder",
            "parvorder",
            "superfamily",
            "family",
            "subfamily",
            "tribe",
            "subtribe",
            "genus",
            "subgenus",
            "section",
            "subsection",
            "series",
            "subseries",
            "species group",
            "species subgroup",
            "species",
            "subspecies"]

def rank_sort_full():
    """Return a dictionary mapping taxonomic rankings to their sort order index."""
    return {rank: index for index, rank in enumerate(return_kingdom_full_sort_order())}

def return_kingdom_limited_order():
    """Return a list of the sort order for the limited set of taxonomic rankings."""
    return ["kingdom",
            "phylum",
            "subphylum",
            "superclass",
            "class",
            "subclass",
            "superorder",
            "order",
            "suborder",
            "infraorder",
            "superfamily",
            "family",
            "subfamily",
            "genus"]

def generate_subsample_priorities(smallest_level    = "family",
                                  largest_level     = "phylum",
                                  output_allsamples = False,
                                  custom_sample_set = None):
    """
    Generate a list of lists of sampling priorities and their fallbacks.
    The list starts with the smallest level and ends with the largest level.

    Inputs:
      - smallest_level: The lowest taxonomic rank to include (default: "genus").
      - largest_level: The highest taxonomic rank to include (default: "class").
      - custom_sample_set: An optional set of samples to prioritize (default: None).
                           Could be a sortable iterable like ["family", "order", "class"].
                           This will limit the outputs and their fallbacks to these ranks,
                           adding the higher rankins than the largest level to avoid
                           incorrectly getting them ranked as None.
    """
    # check that smallest_level is less than largest_level
    full_order = return_kingdom_full_sort_order()
    if full_order.index(smallest_level) <= full_order.index(largest_level):
        print(f"Error: smallest_level '{smallest_level}' must be lower than largest_level '{largest_level}'")
        sys.exit(1)
    sample_order = []
    if custom_sample_set is None:
        sample_order = return_kingdom_limited_order()
    else:
        sample_order = list(custom_sample_set)
        # Add any higher ranks than the largest_level to avoid None rankings
        largest_index = full_order.index(largest_level)
        for rank in full_order[largest_index+1:]:
            if rank not in sample_order:
                sample_order.append(rank)
    # sort sample order from lower to higher rank based on the reverse full order
    sample_order = sorted(sample_order, key=lambda x: full_order.index(x), reverse=True)

    # If the output_allsamples flag is set, add "allsamples" to the front of the list
    if output_allsamples:
        sample_order = ["allsamples"] + sample_order

    # Make a list of lists of priorities and their fallbacks, peeling back the lowest
    #  rank with each iteration.
    output = []
    for i in range(len(sample_order)):
        output.append(sample_order[i:])
        if sample_order[i] == largest_level:
            break
    return output

def subsample_phylogenetically(
    df,
    max_per_bucket   = 10,
    sep             = ";",
    seed            = 42,
    bucket_priority = ("genus", "subfamily", "family", "superfamily", "infraorder", "suborder", "order",
                       "superorder", "subclass", "class", "superclass", "subphylum", "phylum", "kingdom"),
    priority        = False,
    priority_taxids = {9606, 7227, 7739, 6579, 499914},  # H. sapiens, D. melanogaster, B. floridae, P. maximus, R. esculentum
    select_all      = False,
    select_all_rank = "order"
):
    """
    df columns:
      - 'sample'           : ID to return (e.g., assembly accession)
      - 'taxid_list_str'   : semicolon-separated lineage taxids (root > ... > species)

    Args:
        df: pandas DataFrame with 'sample' and 'taxid_list_str' columns
        max_per_bucket: maximum number of samples to select per bucket
        sep: separator used in 'taxid_list_str'
        seed: random seed for reproducibility
        bucket_priority: iterable of taxonomic ranks to use for bucketing (in order of preference, so lower ranks first)
        priority: if True, try to include at most one representative per species of the priority_taxids set
        priority_taxids: set of species-level taxids to prioritize if priority=True
        select_all: If True, also returns an additional bucket with all samples, organized at the specified select_all_rank.
                          This is useful for plotting all the samples in the same pipeline, using the same helpful output functions.
        select_all_rank: The taxonomic rank to use for the "all samples" bucket if select_all=True.

    Returns:
      selected_buckets: {bucket_tid: {
          "rank": str,
          "name": str,
          "chosen": [ {"sample":..., "path":[...], "species_taxid": int} ],
          "original_size": int,
          "final_size": int,
          "priority_count": int,
          "cap_exceeded": bool
      }}
      flat_selected: [sample...]
    """
    rng = random.Random(seed)
    ncbi = NCBITaxa()
    # avoid NameError if a global default isn't defined elsewhere
    priority_set = set(priority_taxids) if priority_taxids is not None else set()

    # -------- Parse lineage paths --------
    recs, all_ids = [], set()
    for _, row in df.iterrows():
        path = [int(x) for x in str(row["taxid_list_str"]).split(sep) if x]
        species_tid = path[-1] if path else None
        recs.append({"sample": row["sample"], "path": path, "species_taxid": species_tid})
        all_ids.update(path)

    rank_map = ncbi.get_rank(list(all_ids))
    name_map = ncbi.get_taxid_translator(list(all_ids))

    def bucket_id(path):
        for rank in bucket_priority:
            hit = next((t for t in path if rank_map.get(t) == rank), None)
            if hit:
                return hit
        return None

    def bucket_id_at_rank(path, rank):
        return next((t for t in path if rank_map.get(t) == rank), None)

    # -------- Bucket --------
    buckets = defaultdict(list)
    if select_all:
        for r in recs:
            buckets[bucket_id_at_rank(r["path"], select_all_rank)].append(r)
    else:
        for r in recs:
            buckets[bucket_id(r["path"])].append(r)

    selected_buckets = {}
    flat_selected = []

    if select_all:
        # take EVERYTHING in each bucket, no caps/priority
        for tid, items in buckets.items():
            chosen = list(items)
            flat_selected.extend(r["sample"] for r in chosen)
            selected_buckets[tid] = {
                "rank": rank_map.get(tid, select_all_rank if tid is not None else "unranked"),
                "name": name_map.get(tid, select_all_rank.capitalize() if tid is not None else "Unranked"),
                "chosen": chosen,
                "original_size": len(items),
                "final_size": len(chosen),
                "priority_count": 0,
                "cap_exceeded": False,
            }
        return selected_buckets, flat_selected

    # -------- Original logic handling --------
    # distance via longest common prefix
    def path_dist(pa, pb):
        i, L = 0, min(len(pa), len(pb))
        while i < L and pa[i] == pb[i]:
            i += 1
        return (len(pa) - i) + (len(pb) - i)

    def farthest_k_with_init(items, k, init_selected):
        """Pick k items maximizing spread given init_selected is already chosen."""
        if k <= 0 or not items:
            return []
        init_ids = {id(x) for x in init_selected}
        remaining = [x for x in items if id(x) not in init_ids]
        if not remaining:
            return []

        picked = list(init_selected)
        if not picked:
            first = rng.choice(remaining)
            picked.append(first)
            remaining.remove(first)

        def nearest_dist(x):
            return min(path_dist(x["path"], p["path"]) for p in picked) if picked else 0

        min_d = {id(x): nearest_dist(x) for x in remaining}
        while remaining and len([p for p in picked if id(p) not in init_ids]) < k:
            nxt = max(remaining, key=lambda x: (min_d[id(x)], x["sample"]))
            picked.append(nxt)
            remaining.remove(nxt)
            for x in remaining:
                d = path_dist(x["path"], nxt["path"])
                if d < min_d[id(x)]:
                    min_d[id(x)] = d

        new_picks = [p for p in picked if id(p) not in init_ids][:k]
        return new_picks

    # Prefer GCF over GCA for the same species (stable fallback = lexicographic sample)
    def pick_one_prefer_gcf(cands):
        if not cands:
            return None
        def tag(rec):
            s = str(rec["sample"]).upper()
            # higher is better
            return (2 if s.startswith("GCF_") else 1 if s.startswith("GCA_") else 0, s)
        return max(cands, key=tag)

    for tid, items in buckets.items():
        # --- Priority: at most ONE representative per priority species in this bucket ---
        chosen = []
        priority_count = 0

        if priority:
            # which priority species are present in this bucket?
            present_priority_species = sorted(set(r["species_taxid"] for r in items) & priority_set)
            for sp_tid in present_priority_species:
                cands = [r for r in items if r["species_taxid"] == sp_tid]
                rep = pick_one_prefer_gcf(cands)
                if rep is not None:
                    chosen.append(rep)
                    priority_count += 1

        # remaining capacity
        remaining_k = max_per_bucket - len(chosen)
        cap_exceeded = remaining_k < 0

        if remaining_k > 0:
            extra = farthest_k_with_init(items, remaining_k, init_selected=chosen)
            chosen.extend(extra)

        flat_selected.extend(r["sample"] for r in chosen)

        selected_buckets[tid] = {
            "rank": rank_map.get(tid, "unranked") if tid is not None else "unranked",
            "name": name_map.get(tid, "Unranked") if tid is not None else "Unranked",
            "chosen": chosen,
            "original_size": len(items),
            "final_size": len(chosen),
            "priority_count": priority_count,
            "cap_exceeded": cap_exceeded,
        }

    return selected_buckets, flat_selected

def _rank_name_maps_full(selected_buckets):
    ncbi = NCBITaxa()
    all_ids = set()
    for tid, info in selected_buckets.items():
        if tid is not None:
            all_ids.add(int(tid))
        for rec in info["chosen"]:
            all_ids.update(rec["path"])
    return ncbi.get_rank(list(all_ids)), ncbi.get_taxid_translator(list(all_ids))

def make_subsampling_report_breadcrumbs(selected_buckets, root_taxid=33208, include_unranked=True):
    """
    One line per selected sample, sorted by the *full* lineage (incl. 'no rank' clades as 'unranked').
    The node corresponding to the bucket taxid is annotated with [rank] <selected/original>.
    """
    # Reuse your helper to collect all taxids and fetch rank/name maps
    rank_map, name_map = _rank_name_maps_full(selected_buckets)

    FULL_RANKS = set(return_kingdom_full_sort_order())
    RANK_POS   = rank_sort_full()  # rank -> position (lower = higher level). Unknown/unranked -> 999.

    def lineage_nodes(path):
        """
        Return [(tid, name, rank)] along the lineage, clipped to root_taxid if present.
        - Ranked nodes: keep first occurrence per rank (prevents duplicates like multiple 'class').
        - Unranked nodes: keep *each taxid* (so Bilateria, Protostomia, etc., are preserved).
        """
        if root_taxid is not None and root_taxid in path:
            path = path[path.index(root_taxid):]

        nodes = []
        seen_ranks = set()
        seen_taxa  = set()
        for t in path:
            rk = rank_map.get(t)
            nm = name_map.get(t, str(t))
            if rk in FULL_RANKS:
                if rk not in seen_ranks:
                    nodes.append((t, nm, rk))
                    seen_ranks.add(rk)
            else:
                if include_unranked and t not in seen_taxa:
                    nodes.append((t, nm, "unranked"))
                    seen_taxa.add(t)
        return nodes

    rows = []
    for tid, info in selected_buckets.items():
        sel, orig = info["final_size"], info["original_size"]

        for rec in info["chosen"]:
            nodes = lineage_nodes(rec["path"])

            # Build the breadcrumb text, marking the bucket node inline
            parts = []
            for t, nm, rk in nodes:
                if t == tid:
                    parts.append(f"{nm} [{rk}] <{sel}/{orig}>")
                else:
                    parts.append(nm)
            line = " › ".join(parts) + f"  -  {rec['sample']}"

            # Uniform sort key: sequence of (rank_pos, name_lower) pairs + tie-break on sample
            key = tuple((RANK_POS.get(rk, 999), nm.lower()) for (_, nm, rk) in nodes) \
                  + ((1000, str(rec["sample"]).lower()),)

            rows.append((key, line))

    rows.sort(key=lambda x: x[0])

    out = []
    out.append("=" * 78)
    out.append("Subsampling Report (breadcrumbs; full-lineage sort, incl. unranked clades)")
    out.append("=" * 78)
    out.extend(l for _, l in rows)
    out.append("-" * 78)
    return "\n".join(out)

def make_subsampling_summary_table(selected_buckets):
    """
    ASCII table sorted by retained % ascending: (final_size / original_size)*100
    """
    # rank sort for grouping (does not affect the main retention sort)
    rank_order = {"infraorder":0, "suborder":1, "order":2, "superorder":3,
                  "subclass":4, "class":5, "unranked":9}

    # collect rows
    rows = []
    for tid, info in selected_buckets.items():
        orig = info["original_size"]
        sel  = info["final_size"]
        retain = (100.0 * sel / orig) if orig else 0.0
        rows.append({
            "rank": info["rank"],
            "name": info["name"],
            "tid": tid,
            "orig": orig,
            "sel": sel,
            "retain": retain,
            "priority": info.get("priority_count", 0),
            "cap_exceeded": info.get("cap_exceeded", False),
            "rank_order": rank_order.get(info["rank"], 99),
        })

    # sort: retained % asc, then bigger buckets first, then name
    rows.sort(key=lambda r: (r["retain"], -r["orig"], r["rank_order"], r["name"].lower()))

    # render
    header = f"{'Bucket (rank)':38} {'TaxID':>9} {'Orig':>6} {'Sel':>6} {'Retain%':>8} {'Prio':>5}"
    line = "-" * len(header)
    out = []
    out.append("=" * len(header))
    out.append("Bucket Summary (sorted by retained % ascending)")
    out.append("=" * len(header))
    out.append(header)
    out.append(line)
    for r in rows:
        label = f"{r['name']} ({r['rank']})"
        prio  = f"{r['priority']}{'!' if r['cap_exceeded'] else ''}"
        out.append(f"{label:38} {str(r['tid']):>9} {r['orig']:>6} {r['sel']:>6} {r['retain']:>7.1f} {prio:>5}")
    out.append(line)
    return "\n".join(out)

def make_subsampling_report_tree(
    selected_buckets,
    root_taxid=33208,
    include_unranked=True,
    suppress_higher_fallback_labels=True
):
    """
    Indented tree from Metazoa (if present) → species, including 'no rank' clades.
    - Internal nodes that are buckets get [rank] <sel/orig>, unless suppressed.
    - Species leaves get:  [bucket: NAME (RANK) <sel/orig>]
    """
    from ete4 import NCBITaxa

    ncbi = NCBITaxa()

    # ---- Collect IDs to translate ----
    all_ids = set()
    for tid, info in selected_buckets.items():
        if tid is not None:
            all_ids.add(int(tid))
        for rec in info["chosen"]:
            all_ids.update(rec["path"])

    rank_map = ncbi.get_rank(list(all_ids))                 # taxid -> rank (e.g., 'phylum', 'no rank', ...)
    name_map = ncbi.get_taxid_translator(list(all_ids))     # taxid -> name

    bucket_info = {tid: info for tid, info in selected_buckets.items()}

    # sample -> its bucket (for leaf labels)
    bucket_for_sample = {}
    for tid, inf in selected_buckets.items():
        for rec in inf["chosen"]:
            bucket_for_sample[rec["sample"]] = {
                "name": inf["name"],
                "rank": inf["rank"],
                "sel":  inf["final_size"],
                "orig": inf["original_size"],
            }

    FULL_RANKS = set(return_kingdom_full_sort_order())
    RANK_POS   = rank_sort_full()  # used only to sort siblings; unknown/unranked -> 999

    def ranked_nodes_from_path(path):
        """Return [(tid,name,rank)] along the lineage, preserving multiple 'no rank' clades."""
        if root_taxid is not None and root_taxid in path:
            path = path[path.index(root_taxid):]

        nodes = []
        seen_ranks = set()  # for ranked nodes
        seen_taxa  = set()  # for unranked nodes so multiple clades survive

        for t in path:
            rk = rank_map.get(t)
            nm = name_map.get(t, str(t))

            if rk in FULL_RANKS:
                if rk not in seen_ranks:
                    nodes.append((t, nm, rk))
                    seen_ranks.add(rk)
            else:
                if include_unranked and t not in seen_taxa:
                    nodes.append((t, nm, "unranked"))
                    seen_taxa.add(t)
        return nodes

    # Build a simple tree
    root = {"tid": None, "name": "ROOT", "rank": "unranked", "children": {}, "samples": []}

    def insert_record(rec):
        nodes = ranked_nodes_from_path(rec["path"])
        cur = root
        for t, nm, rk in nodes:
            if t not in cur["children"]:
                cur["children"][t] = {"tid": t, "name": nm, "rank": rk, "children": {}, "samples": []}
            cur = cur["children"][t]
        cur["samples"].append(rec["sample"])

    for tid, info in selected_buckets.items():
        for rec in info["chosen"]:
            insert_record(rec)

    # Helpers to decide when to hide fallback labels
    def deepest_bucket_rank_pos_in_subtree(node):
        best = None
        if node["tid"] in bucket_info:
            rk = bucket_info[node["tid"]]["rank"]
            best = RANK_POS.get(rk, 999)
        for ch in node["children"].values():
            child_best = deepest_bucket_rank_pos_in_subtree(ch)
            if child_best is not None:
                best = child_best if best is None else max(best, child_best)
        return best

    def has_deeper_bucket_than(node_rank_pos, node):
        best = deepest_bucket_rank_pos_in_subtree(node)
        return best is not None and best > node_rank_pos

    # Emit lines
    def emit(node, depth=0):
        lines = []
        if node["tid"] is not None:
            label = node["name"]
            if node["tid"] in bucket_info:
                inf = bucket_info[node["tid"]]
                node_rank_pos = RANK_POS.get(inf["rank"], 999)
                if suppress_higher_fallback_labels and inf["rank"] in {"class", "subclass"} and has_deeper_bucket_than(node_rank_pos, node):
                    pass
                else:
                    label += f" [{inf['rank']}] <{inf['final_size']}/{inf['original_size']}>"
            lines.append(" " * (2 * depth) + label)

        # Sort children primarily by rank position, then name
        children_sorted = sorted(
            node["children"].values(),
            key=lambda c: (RANK_POS.get(c["rank"], 999), c["name"].lower())
        )
        for ch in children_sorted:
            lines.extend(emit(ch, depth + (0 if node["tid"] is None else 1)))

        # Species leaves with bucket annotation
        if node["samples"]:
            indent = " " * (2 * (depth + 1))
            for s in sorted(node["samples"], key=str):
                b = bucket_for_sample.get(s)
                if b:
                    leaf = f"{s}  [bucket: {b['name']} ({b['rank']}) <{b['sel']}/{b['orig']}>]"
                else:
                    leaf = s
                lines.append(indent + leaf)
        return lines

    lines = []
    lines.append("=" * 78)
    lines.append("Subsampling Report (full tree, incl. unranked clades; anchored at Metazoa if present)")
    lines.append("=" * 78)

    for ch in sorted(root["children"].values(),
                     key=lambda c: (RANK_POS.get(c["rank"], 999), c["name"].lower())):
        lines.extend(emit(ch, 0))
        lines.append("-" * 78)

    return "\n".join(lines)

def main():
    rbhdir    = "/lisc/scratch/molevo/dts/manifold/BCnSSimakov2022_current_rbh_202509"
    sampletsv = "/lisc/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP/sampledf.tsv"

    df = pd.read_csv(sampletsv, sep="\t", index_col=0)
    print(df.columns)

    selected_buckets, flat = subsample_phylogenetically(df, max_per_bucket=10, priority=True)
    summary_txt     = make_subsampling_summary_table(selected_buckets)
    breadcrumbs_txt = make_subsampling_report_breadcrumbs(selected_buckets)
    tree_txt        = make_subsampling_report_tree(selected_buckets)

    outdir = Path("subsampling_reports")
    outdir.mkdir(parents=True, exist_ok=True)

    (outdir / "summary_table.txt").write_text(summary_txt, encoding="utf-8")
    (outdir / "breadcrumbs.txt").write_text(breadcrumbs_txt, encoding="utf-8")
    (outdir / "tree.txt").write_text(tree_txt, encoding="utf-8")

    # Save flat list of selected samples
    (outdir / "selected_samples.txt").write_text("\n".join(sorted(flat)), encoding="utf-8")

    # print out the sort orders for testing
    for sortorder in generate_subsample_priorities():
        print(sortorder)

if __name__ == "__main__":
    main()