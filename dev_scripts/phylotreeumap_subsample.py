#!/usr/bin/env python3
from ete4 import NCBITaxa
from collections import defaultdict
import pandas as pd
from pathlib import Path
import random
import sys


def subsample_phylogenetically(
    df,
    max_per_order   = 10,
    sep             = ";",
    seed            = 42,
    bucket_priority = ("infraorder", "suborder", "order",
                       "superorder", "subclass", "class"),
    priority        = False,
    priority_taxids = {9606, 7227, 7739, 6579, 499914},  # H. sapiens, D. melanogaster, B. floridae, P. maximus, R. esculentum
):
    """
    df columns:
      - 'sample'           : ID to return (e.g., assembly accession)
      - 'taxid_list_str'   : semicolon-separated lineage taxids (root→…→species)
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

    # Parse and collect IDs
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

    # Bucket
    buckets = defaultdict(list)
    for r in recs:
        buckets[bucket_id(r["path"])].append(r)

    selected_buckets = {}
    flat_selected = []

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
        remaining_k = max_per_order - len(chosen)
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

# Full rank ordering + helpers (can live next to your other report fns)
FULL_RANKS = (
    "superkingdom","kingdom","subkingdom",
    "superphylum","phylum","subphylum","infraphylum",
    "superclass","class","subclass","infraclass",
    "cohort","subcohort",
    "superorder","order","suborder","infraorder","parvorder",
    "superfamily","family","subfamily","tribe","subtribe",
    "genus","subgenus","species"
)
RANK_SORT_FULL = {r:i for i,r in enumerate(FULL_RANKS)}
METAZOA_TAXID_DEFAULT = 33208  # anchor if present

def _rank_name_maps_full(selected_buckets):
    ncbi = NCBITaxa()
    all_ids = set()
    for tid, info in selected_buckets.items():
        if tid is not None:
            all_ids.add(int(tid))
        for rec in info["chosen"]:
            all_ids.update(rec["path"])
    return ncbi.get_rank(list(all_ids)), ncbi.get_taxid_translator(list(all_ids))

def _ranked_nodes_from_path(path, rank_map, name_map, root_taxid=METAZOA_TAXID_DEFAULT):
    """Ordered [(tid,name,rank)] along FULL_RANKS; clip to root_taxid if present."""
    if root_taxid is not None and root_taxid in path:
        path = path[path.index(root_taxid):]
    seen = set()
    nodes = []
    for t in path:
        rk = rank_map.get(t)
        if rk in FULL_RANKS and rk not in seen:
            nodes.append((t, name_map.get(t, str(t)), rk))
            seen.add(rk)
    return nodes

def make_subsampling_report_breadcrumbs(selected_buckets, root_taxid=METAZOA_TAXID_DEFAULT):
    """
    One line per selected sample, sorted by full ranked lineage (Metazoa-anchored).
    Bucket node in each line is annotated with [rank] <selected/original>.
    """
    rank_map, name_map = _rank_name_maps_full(selected_buckets)
    bucket_info = {tid: info for tid, info in selected_buckets.items()}

    # Build rows: (sort_key, line)
    rows = []
    for tid, info in selected_buckets.items():
        sel, orig = info["final_size"], info["original_size"]
        for rec in info["chosen"]:
            nodes = _ranked_nodes_from_path(rec["path"], rank_map, name_map, root_taxid=root_taxid)

            # format breadcrumb, marking the bucket node
            parts = []
            for t, nm, rk in nodes:
                if t == tid:
                    parts.append(f"{nm} [{rk}] <{sel}/{orig}>")
                else:
                    parts.append(nm)
            line = " › ".join(parts) + f"  -  {rec['sample']}"

            # full-tree sort key: ((rank_idx, name_lower) ... , sample)
            key = tuple((RANK_SORT_FULL.get(rk, 999), nm.lower()) for _, nm, rk in nodes) + (str(rec["sample"]),)
            rows.append((key, line))

    rows.sort(key=lambda x: x[0])

    out = []
    out.append("=" * 78)
    out.append("Subsampling Report (breadcrumbs, full-tree order; anchored at Metazoa if present)")
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

# A comprehensive ordered rank list (top → bottom). Adjust if you like.
FULL_RANKS = (
    "superkingdom","kingdom","subkingdom",
    "superphylum","phylum","subphylum","infraphylum",
    "superclass","class","subclass","infraclass",
    "cohort","subcohort",
    "superorder","order","suborder","infraorder","parvorder",
    "superfamily","family","subfamily","tribe","subtribe",
    "genus","subgenus","species"
)
RANK_SORT_FULL = {r:i for i,r in enumerate(FULL_RANKS)}
METAZOA_TAXID_DEFAULT = 33208  # Metazoa

def make_subsampling_report_tree(selected_buckets, root_taxid=METAZOA_TAXID_DEFAULT):
    """
    Indented tree from Metazoa (if present) → species, using FULL_RANKS order.
    Bucket nodes are annotated with [rank] <sel/orig>.
    """
    ncbi = NCBITaxa()

    # gather taxids to translate
    all_ids = set()
    for tid, info in selected_buckets.items():
        if tid is not None:
            all_ids.add(int(tid))
        for rec in info["chosen"]:
            all_ids.update(rec["path"])

    rank_map = ncbi.get_rank(list(all_ids))
    name_map = ncbi.get_taxid_translator(list(all_ids))

    bucket_info = {tid: info for tid, info in selected_buckets.items()}

    # Build tree
    root = {"tid": None, "name": "ROOT", "rank": "unranked", "children": {}, "samples": []}

    def ranked_nodes_from_path(path):
        """Return ordered [(tid,name,rank)] along the path, filtered to FULL_RANKS and
           clipped to start at root_taxid if present."""
        # clip to Metazoa root if available
        if root_taxid is not None and root_taxid in path:
            start_idx = path.index(root_taxid)
            path2 = path[start_idx:]
        else:
            # start at the first ranked node we know (top-most in FULL_RANKS)
            path2 = path

        nodes = []
        seen_ranks = set()
        for t in path2:
            rk = rank_map.get(t)
            if rk in FULL_RANKS and rk not in seen_ranks:
                nodes.append((t, name_map.get(t, str(t)), rk))
                seen_ranks.add(rk)
        return nodes

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

    def emit(node, depth=0):
        lines = []
        if node["tid"] is not None:
            label = node["name"]
            if node["tid"] in bucket_info:
                inf = bucket_info[node["tid"]]
                label += f" [{inf['rank']}] <{inf['final_size']}/{inf['original_size']}>"
            lines.append(" " * (2 * depth) + label)

        # sort children by rank, then name (stable)
        children_sorted = sorted(
            node["children"].values(),
            key=lambda c: (RANK_SORT_FULL.get(c["rank"], 999), c["name"].lower())
        )
        for ch in children_sorted:
            lines.extend(emit(ch, depth + (0 if node["tid"] is None else 1)))

        if node["samples"]:
            for s in sorted(node["samples"], key=str):
                lines.append(" " * (2 * (depth + 1)) + s)
        return lines

    lines = []
    lines.append("=" * 78)
    lines.append("Subsampling Report (full tree, anchored at Metazoa if present)")
    lines.append("=" * 78)
    # print each top-level child (e.g., Metazoa) from ROOT
    for ch in sorted(root["children"].values(),
                     key=lambda c: (RANK_SORT_FULL.get(c["rank"], 999), c["name"].lower())):
        lines.extend(emit(ch, 0))
        lines.append("-" * 78)
    return "\n".join(lines)

rbhdir    = "/lisc/scratch/molevo/dts/manifold/BCnSSimakov2022_current_rbh_202509"
sampletsv = "/lisc/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP/sampledf.tsv"

df = pd.read_csv(sampletsv, sep="\t", index_col=0)
print(df.columns)

selected_buckets, flat = subsample_phylogenetically(df, max_per_order=10, priority=True)
summary_txt    = make_subsampling_summary_table(selected_buckets)
breadcrumbs_txt = make_subsampling_report_breadcrumbs(selected_buckets)
tree_txt        = make_subsampling_report_tree(selected_buckets)

outdir = Path("subsampling_reports")
outdir.mkdir(parents=True, exist_ok=True)

(outdir / "summary_table.txt").write_text(summary_txt, encoding="utf-8")
(outdir / "breadcrumbs.txt").write_text(breadcrumbs_txt, encoding="utf-8")
(outdir / "tree.txt").write_text(tree_txt, encoding="utf-8")

# Save flat list of selected samples
(outdir / "selected_samples.txt").write_text("\n".join(sorted(flat)), encoding="utf-8")