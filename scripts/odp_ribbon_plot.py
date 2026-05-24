#!/usr/bin/env python
"""
This plot contains all of the functions
  related to making the ribbon plots.

These functions can be called by the other
  programs.
"""

# ODP-specific imports
import odp_plotting_functions as odp_plot

# Other standard python libraries
import bisect
import logging

# non-standard dependencies
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.path as mpath
import matplotlib.patches as patches
import pandas as pd

import sys


def _count_inversions(seq):
    """Number of pairs (i,j) with i<j and seq[i]>seq[j]. O(n log n) via
    bisect into a sorted side buffer. Used to score ribbon crossings: lay
    out line endpoints sorted by the top species' x, then count inversions
    in the bottom species' x — each inversion is one crossing of bezier
    paths in plot_bezier_lines (paths are monotonic in y, so two paths
    cross iff their endpoint orderings flip)."""
    side = []
    inv = 0
    for x in reversed(list(seq)):
        i = bisect.bisect_left(side, x)
        inv += i
        bisect.insort(side, x)
    return inv


def _optimize_chrom_flips_top_down(species_order, rbh_df_list,
                                   sp_pair_to_rbh_df_list_index,
                                   sp_to_chromorder, sp_to_genesdfs,
                                   max_passes=8, verbose=False,
                                   initial_flip=None):
    """Decide, for every scaffold of every species, whether to plot it in the
    fasta-forward direction or reversed, so that the total number of
    ribbon-line crossings between each adjacent pair of species is
    minimized. Greedy cascade top-down: the top species is fixed
    forward; for each species below, the species above is already pinned
    and only this species' flips are optimized.

    Returns a dict sp -> {scaf: bool}, True means "reverse this scaffold
    relative to the assembly fasta".

    Two-stage optimization per cascade step:
      1. Independent per-scaffold seed: for each bottom-side scaffold,
         compare crossings inside that scaffold (lines on the same B)
         with and without flipping; take the cheaper. Cheap, O(L log L).
      2. Iterative greedy global refinement: try flipping each scaffold
         in turn; keep it flipped if the total crossing count drops.
         Repeats until a full pass makes no improvement.

    Inter-scaffold crossings are dominated by chromosome *order*, not
    rotation, so the greedy global pass usually converges in 1–3 sweeps.
    """
    # If the caller supplied an initial flip state, honor it (esp.
    # the top row, which the rest of this function does NOT re-decide
    # -- it only optimizes rows 1+). Otherwise start everything at
    # False.
    if initial_flip is not None:
        sp_to_flip = {sp: dict(initial_flip.get(sp, {}))
                      for sp in species_order}
        for sp in species_order:
            for scaf in sp_to_chromorder[sp]:
                sp_to_flip[sp].setdefault(scaf, False)
    else:
        sp_to_flip = {sp: {scaf: False for scaf in sp_to_chromorder[sp]}
                      for sp in species_order}

    for i in range(1, len(species_order)):
        topsp = species_order[i - 1]
        botsp = species_order[i]
        top_scaf_col = "{}_scaf".format(topsp)
        bot_scaf_col = "{}_scaf".format(botsp)
        top_gene_col = "{}_gene".format(topsp)
        bot_gene_col = "{}_gene".format(botsp)
        top_order = sp_to_chromorder[topsp]
        bot_order = sp_to_chromorder[botsp]

        lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([topsp, botsp]))]
        df = rbh_df_list[lookup]
        mask = df[top_scaf_col].isin(top_order) & df[bot_scaf_col].isin(bot_order)
        df = df.loc[mask, [top_scaf_col, bot_scaf_col, top_gene_col, bot_gene_col]].copy()

        top_gene_to_rank = dict(zip(sp_to_genesdfs[topsp][top_gene_col],
                                    sp_to_genesdfs[topsp]["{}_chromIx".format(topsp)]))
        bot_gene_to_rank = dict(zip(sp_to_genesdfs[botsp][bot_gene_col],
                                    sp_to_genesdfs[botsp]["{}_chromIx".format(botsp)]))
        top_scaf_size = sp_to_genesdfs[topsp].groupby(top_scaf_col).size().to_dict()
        bot_scaf_size = sp_to_genesdfs[botsp].groupby(bot_scaf_col).size().to_dict()

        df["_top_rank"] = df[top_gene_col].map(top_gene_to_rank)
        df["_bot_rank"] = df[bot_gene_col].map(bot_gene_to_rank)
        df = df.dropna(subset=["_top_rank", "_bot_rank"]).reset_index(drop=True)
        if df.empty:
            continue
        df["_top_sz"]  = df[top_scaf_col].map(top_scaf_size).astype(float)
        df["_bot_sz"]  = df[bot_scaf_col].map(bot_scaf_size).astype(float)
        df["_top_ord"] = df[top_scaf_col].map(top_order).astype(float)
        df["_bot_ord"] = df[bot_scaf_col].map(bot_order).astype(float)

        top_flip = sp_to_flip[topsp]
        top_flipped_mask = df[top_scaf_col].map(lambda s: top_flip.get(s, False))
        df["_top_eff"] = df["_top_rank"].where(~top_flipped_mask,
                                               df["_top_sz"] - 1 - df["_top_rank"])
        df["_top_pos"] = df["_top_ord"] + df["_top_eff"] / df["_top_sz"]
        df = df.sort_values("_top_pos", kind="mergesort").reset_index(drop=True)

        bot_scaf = df[bot_scaf_col].to_numpy()
        bot_rank = df["_bot_rank"].to_numpy()
        bot_sz   = df["_bot_sz"].to_numpy()
        bot_ord  = df["_bot_ord"].to_numpy()

        flips = dict(sp_to_flip[botsp])

        def bot_positions(flips_dict):
            flipped = pd.Series(bot_scaf).map(lambda s: flips_dict.get(s, False)).to_numpy()
            eff = bot_rank.copy()
            eff[flipped] = bot_sz[flipped] - 1 - bot_rank[flipped]
            return bot_ord + eff / bot_sz

        # Stage 1: independent per-chrom seed.
        for scaf, sub in df.groupby(bot_scaf_col, sort=False):
            if len(sub) < 2:
                continue
            ranks = sub["_bot_rank"].tolist()
            sz = bot_scaf_size[scaf] - 1
            if _count_inversions([sz - r for r in ranks]) < _count_inversions(ranks):
                flips[scaf] = True

        # Stage 2: iterated greedy refinement on full crossing count.
        cur = _count_inversions(bot_positions(flips))
        for _ in range(max_passes):
            improved = False
            for scaf in list(flips.keys()):
                flips[scaf] = not flips[scaf]
                new = _count_inversions(bot_positions(flips))
                if new < cur:
                    cur = new
                    improved = True
                else:
                    flips[scaf] = not flips[scaf]
            if not improved:
                break

        sp_to_flip[botsp] = flips
        if verbose:
            print("flip cascade {} -> {}: crossings={}, flipped={}".format(
                topsp, botsp, cur,
                sum(1 for v in flips.values() if v)))

    # Final boundary pass: revisit EVERY chromosome in EVERY row
    # (including row 0, which the cascade above skips) and try
    # flipping it against the full crossing count. Row 0 was
    # otherwise pinned forward forever -- if e.g. Plutella CM/509.1
    # would reduce total crossings when flipped, this is where we
    # catch it. We also re-visit the inner rows because the
    # cascade-greedy order may have missed flips that only pay off
    # once the rows below have settled into their final order.
    def _eval_pair(top_sp, bot_sp):
        lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([top_sp, bot_sp]))]
        df_in = rbh_df_list[lookup]
        tc, bc = "{}_scaf".format(top_sp), "{}_scaf".format(bot_sp)
        tg, bg = "{}_gene".format(top_sp), "{}_gene".format(bot_sp)
        m = df_in[tc].isin(sp_to_chromorder[top_sp]) & df_in[bc].isin(sp_to_chromorder[bot_sp])
        d = df_in.loc[m].copy()
        trk = dict(zip(sp_to_genesdfs[top_sp][tg],
                       sp_to_genesdfs[top_sp]["{}_chromIx".format(top_sp)]))
        brk = dict(zip(sp_to_genesdfs[bot_sp][bg],
                       sp_to_genesdfs[bot_sp]["{}_chromIx".format(bot_sp)]))
        tsz = sp_to_genesdfs[top_sp].groupby(tc).size().to_dict()
        bsz = sp_to_genesdfs[bot_sp].groupby(bc).size().to_dict()
        d["_tr"] = d[tg].map(trk); d["_br"] = d[bg].map(brk)
        d = d.dropna(subset=["_tr", "_br"])
        if d.empty:
            return 0.0
        t_flip = sp_to_flip[top_sp]; b_flip = sp_to_flip[bot_sp]
        d["_tsz"] = d[tc].map(tsz).astype(float)
        d["_bsz"] = d[bc].map(bsz).astype(float)
        d["_te"] = d.apply(lambda r: (r["_tsz"] - 1 - r["_tr"])
                           if t_flip.get(r[tc], False) else r["_tr"], axis=1)
        d["_be"] = d.apply(lambda r: (r["_bsz"] - 1 - r["_br"])
                           if b_flip.get(r[bc], False) else r["_br"], axis=1)
        d["_tp"] = d[tc].map(sp_to_chromorder[top_sp]).astype(float) + d["_te"] / d["_tsz"]
        d["_bp"] = d[bc].map(sp_to_chromorder[bot_sp]).astype(float) + d["_be"] / d["_bsz"]
        d = d.sort_values("_tp", kind="mergesort")
        return _count_inversions(d["_bp"].tolist())

    pair_cache = {k: _eval_pair(species_order[k], species_order[k + 1])
                  for k in range(len(species_order) - 1)}

    def _rescore_neighbors(sp_idx):
        if sp_idx > 0:
            pair_cache[sp_idx - 1] = _eval_pair(species_order[sp_idx - 1],
                                                species_order[sp_idx])
        if sp_idx < len(species_order) - 1:
            pair_cache[sp_idx] = _eval_pair(species_order[sp_idx],
                                            species_order[sp_idx + 1])

    cur_score = sum(pair_cache.values())
    for _ in range(max_passes):
        any_improvement = False
        for sp_idx, sp in enumerate(species_order):
            for scaf in list(sp_to_flip[sp].keys()):
                sp_to_flip[sp][scaf] = not sp_to_flip[sp][scaf]
                _rescore_neighbors(sp_idx)
                new_score = sum(pair_cache.values())
                if new_score < cur_score - 1e-9:
                    cur_score = new_score
                    any_improvement = True
                else:
                    sp_to_flip[sp][scaf] = not sp_to_flip[sp][scaf]
                    _rescore_neighbors(sp_idx)
        if not any_improvement:
            break

    return sp_to_flip

def _quality_check_chromosome_list(sp, templist, sp_to_chr_to_size,
                                   sp_to_gene_order, sp_min_chr_size):
    """
    We perform some quality checks on the list of chromosomes to make sure that they are all valid.

    templist is the temporary list of scaffold IDs
    """
    # Now we perform some checks on the list of chromosomes and add them to the dict
    # check here to make sure that there aren't duplicate entries in the chromosomes
    if len(templist) != len(set(templist)):
        raise IOError("There are some duplicate chromosome entries for {}: {}.".format(sp, templist))
    # Make sure that all the chromosomes in the chromorder are actually in the fasta file
    checklist = [x for x in templist if x not in sp_to_chr_to_size[sp]]
    if len(checklist) > 0:
        raise IOError("These chromosomes in the chromosome order list for sp {} were not in the fasta file: {}".format(
            sp, checklist))
    # if we included a list of chromosomes below, remove those 
    #TODO DELETE
    print("sp_to_gene_order: {}".format(sp_to_gene_order))
    if sp in sp_to_gene_order:
        templist2 = [x for x in templist if x in sp_to_gene_order[sp]]
        templist = templist2
    else:
        # remove things smaller than the minimum size
        templist2 = [x for x in templist if (x in sp_to_chr_to_size[sp]) and (sp_to_chr_to_size[sp][x] >= sp_min_chr_size[sp])]
        templist = templist2
    return templist

def _optimize_spA_based_on_rbh(sp, prevsp, rbhdf, sp_to_chromorder, sp_to_gene_order):
    """
    Takes in a species, and optimizes the order of the chromosomes based on the rbh file and the order of the other species' chromosomes.

    sp is the species for which we want to optimize the order
    prevsp is the other species in the rbh dataframe

    Returns the optimized order of this species' chromosomes
    """
    # get the significance of the interactions
    rbhdf = rbhdf.copy()
    rbhdf = rbhdf[["{}_scaf".format(x) for x in [sp, prevsp]] + ["whole_FET"]]
    countsdf = rbhdf.groupby(by=["{}_scaf".format(sp), "{}_scaf".format(prevsp)]).count().reset_index()
    countsdf.columns = ["{}_scaf".format(sp), "{}_scaf".format(prevsp), "count"]
    countsdf = countsdf.sort_values(by=["count"], ascending=False).reset_index(drop=True)
    rbhdf = rbhdf.drop_duplicates().reset_index(drop=True)
    rbhdf = rbhdf.merge(countsdf, on=["{}_scaf".format(sp), "{}_scaf".format(prevsp)])
    # get rid of scaffolds that weren't in the custom list, if the custom list exists
    if sp in sp_to_gene_order:
        rbhdf = rbhdf.loc[rbhdf["{}_scaf".format(prevsp)].isin(sp_to_chromorder[prevsp])]
    # put the scaffold in the best place based on what is above it
    rbhdf["sort_index"] = rbhdf["{}_scaf".format(prevsp)].map(sp_to_chromorder[prevsp])
    rbhdf.dropna(subset=["sort_index"], inplace=True)
    rbhdf = rbhdf.sort_values(by=["{}_scaf".format(sp), "whole_FET"], ascending=[True, True])
    rbhdf.drop_duplicates(subset=["{}_scaf".format(sp)], inplace=True)
    rbhdf.sort_values(by=["sort_index", "whole_FET"], ascending=[True, True], inplace=True)
    rbhdf.reset_index(drop=True, inplace=True)
    return rbhdf["{}_scaf".format(sp)].tolist()

def _draw_oriented_bar(panel, x1, x2, y, flipped, lw=0.9,
                       head_len_frac=0.15, head_len_cap=0.0045,
                       head_h=0.07, color="black"):
    """Render a chromosome as a chevron-tipped bar: the bar ends in a
    small arrowhead whose tip lies at the fasta-3' end of the scaffold.
    Replaces the separate-triangle marker so the arrowhead IS the bar
    rather than an extra glyph next to it. Forward scaffolds tip right,
    reversed scaffolds tip left."""
    width = x2 - x1
    if width <= 0:
        return
    head_len = min(width * head_len_frac, head_len_cap)
    shaft_w  = max(width - head_len, 0)
    if flipped:
        tip_x   = x1
        base_x  = x1 + head_len
        shaft_x = (base_x, x2)
    else:
        tip_x   = x2
        base_x  = x2 - head_len
        shaft_x = (x1, base_x)
    # shaft
    panel.plot([shaft_x[0], shaft_x[1]], [y, y],
               color=color, lw=lw, solid_capstyle="butt", zorder=3)
    # filled arrowhead
    head = patches.Polygon(
        [(tip_x, y), (base_x, y - head_h), (base_x, y + head_h)],
        closed=True, facecolor=color, edgecolor=color, lw=0.3, zorder=4,
    )
    panel.add_patch(head)


def _shorten_chrom_labels(labels):
    """Compact a scaffold id to "<first 2 chars>/<tail>", where tail is
    chosen to keep the assembly-version suffix visible whenever the id
    ends with a "<digits>.<digits>" pattern (the usual NCBI form):
        NC_051358.1 -> NC/358.1
        NC_051358   -> NC/358
        JAJNFR010000001.1 -> JA/001.1
    Short ids (<=6 chars) are returned unchanged."""
    import re as _re
    decimal_tail = _re.compile(r"\.([^.]+)$")
    out = {}
    for lbl in labels:
        if len(lbl) <= 6:
            out[lbl] = lbl
            continue
        m = decimal_tail.search(lbl)
        if m:
            stem = lbl[:m.start()]
            tail = "{}.{}".format(stem[-3:], m.group(1))
            head = stem[:2]
        else:
            head = lbl[:2]
            tail = lbl[-3:]
        out[lbl] = "{}/{}".format(head, tail)
    return out


_DETANGLE_MODES = {"optimal-barycenter", "optimal-barycenter-iter",
                   "optimal-median", "optimal-swap"}


def _build_fet_lg_clusters(species_order, rbh_df_list, top_k=1):
    """Build a chromosome-level graph from every adjacent-row RBH:
    nodes are (species, scaffold) tuples, edges are kept only when
    they are *mutual best* (each side's top-k FET-strongest partner
    list contains the other) AND have whole_FET <= 0.05. With top_k=1
    each chrom links to at most one chrom per neighboring species,
    which keeps the graph sparse enough that one homology block per
    component falls out cleanly even on data with universal trace
    similarity. top_k>1 accommodates fission events where a single
    ancestral chrom split into k pieces in one lineage.

    Returns
        node_to_cluster : {(sp, scaf): cluster_id}
        cluster_size    : {cluster_id: int}
        cluster_weight  : {cluster_id: total edge weight}
    where clusters are connected components ordered by total weight
    descending.

    Two Cydia chroms whose only sharp link is to the same Bombyx
    chrom (an in-Bombyx fusion of two ancestral LGs) end up in the
    same component via that Bombyx node, even though their dominant
    ALG labels differ. This is the case the per-scaffold
    gene_group-mode approach misses; the FET-graph captures it.
    """
    from collections import defaultdict
    adj = defaultdict(set)
    edge_w = {}
    # Per-pair top-K partner table, kept for the orphan-rescue pass.
    per_pair_data = []
    all_nodes = set()
    for df in rbh_df_list:
        if "whole_FET" not in df.columns:
            continue
        species_scaf_cols = [c for c in df.columns if c.endswith("_scaf")]
        sps = [c[:-5] for c in species_scaf_cols]
        if len(sps) != 2:
            continue
        sa, sb = sps
        sca, scb = "{}_scaf".format(sa), "{}_scaf".format(sb)
        sig = df[df["whole_FET"] <= 0.05]
        if sig.empty:
            continue
        grp = sig.groupby([sca, scb]).agg(
            count=("whole_FET", "size"),
            best_FET=("whole_FET", "min")).reset_index()
        a_top = (grp.sort_values([sca, "best_FET"])
                    .groupby(sca, sort=False).head(top_k))
        b_top = (grp.sort_values([scb, "best_FET"])
                    .groupby(scb, sort=False).head(top_k))
        a_pairs = set(zip(a_top[sca], a_top[scb]))
        b_pairs = set(zip(b_top[sca], b_top[scb]))
        # Strict mutual best (intersection of top-1 lists).
        edges = a_pairs & b_pairs
        for (ca, cb) in edges:
            na = (sa, ca); nb = (sb, cb)
            adj[na].add(nb); adj[nb].add(na)
            cnt = int(grp.loc[(grp[sca] == ca) & (grp[scb] == cb), "count"].iloc[0])
            edge_w[(na, nb)] = edge_w[(nb, na)] = cnt
        for s in set(grp[sca].unique()):
            all_nodes.add((sa, s))
        for s in set(grp[scb].unique()):
            all_nodes.add((sb, s))
        per_pair_data.append((sa, sb, grp, a_top, b_top))

    # Global orphan rescue: chromosomes that received no mutual-best
    # edge anywhere in the graph need one attachment so they don't
    # form singleton "clusters" placed far from their biological
    # siblings. Restricting the rescue to fully-orphan nodes (no
    # edges in any pair) prevents a chrom that already lives in an
    # LG component on one side of the plot from accidentally
    # bridging two distinct LG components on the other side via a
    # fusion hub. Each orphan gets a single edge to its single
    # FET-strongest partner observed across all adjacent-pair RBHs.
    connected = set(adj.keys())
    orphans = all_nodes - connected
    for sp, scaf in orphans:
        best = None  # (count, FET, partner_sp, partner_scaf)
        for (sa, sb, grp, a_top, b_top) in per_pair_data:
            if sp == sa:
                row = a_top[a_top["{}_scaf".format(sa)] == scaf]
                if row.empty:
                    continue
                partner = row["{}_scaf".format(sb)].iloc[0]
                cnt = int(grp.loc[(grp["{}_scaf".format(sa)] == scaf)
                                  & (grp["{}_scaf".format(sb)] == partner),
                                  "count"].iloc[0])
                fet = float(grp.loc[(grp["{}_scaf".format(sa)] == scaf)
                                    & (grp["{}_scaf".format(sb)] == partner),
                                    "best_FET"].iloc[0])
                if best is None or fet < best[1]:
                    best = (cnt, fet, sb, partner)
            elif sp == sb:
                row = b_top[b_top["{}_scaf".format(sb)] == scaf]
                if row.empty:
                    continue
                partner = row["{}_scaf".format(sa)].iloc[0]
                cnt = int(grp.loc[(grp["{}_scaf".format(sa)] == partner)
                                  & (grp["{}_scaf".format(sb)] == scaf),
                                  "count"].iloc[0])
                fet = float(grp.loc[(grp["{}_scaf".format(sa)] == partner)
                                    & (grp["{}_scaf".format(sb)] == scaf),
                                    "best_FET"].iloc[0])
                if best is None or fet < best[1]:
                    best = (cnt, fet, sa, partner)
        if best is None:
            continue
        cnt, _, partner_sp, partner_scaf = best
        na = (sp, scaf); nb = (partner_sp, partner_scaf)
        adj[na].add(nb); adj[nb].add(na)
        edge_w[(na, nb)] = edge_w[(nb, na)] = cnt
    visited = set()
    clusters = []
    for node in adj:
        if node in visited:
            continue
        stack = [node]
        comp = []
        while stack:
            n = stack.pop()
            if n in visited:
                continue
            visited.add(n)
            comp.append(n)
            stack.extend(adj[n])
        clusters.append(comp)

    def cluster_weight(c):
        s = set(c)
        return sum(w for (a, b), w in edge_w.items()
                   if a in s and b in s) // 2

    clusters.sort(key=cluster_weight, reverse=True)
    raw_node_to_cluster = {}
    for i, c in enumerate(clusters):
        for n in c:
            raw_node_to_cluster[n] = i

    # Cross-cluster edge weights, built from EVERY FET-significant
    # pair (not just mutual-best). This is what surfaces fusion /
    # fission relationships: two Cydia chroms in different clusters
    # that both link to the same Bombyx chrom will have a strong
    # cross-cluster signal via that Bombyx node, which gets used by
    # the cluster ordering below to keep them plot-adjacent.
    from collections import Counter
    cross = Counter()
    for df in rbh_df_list:
        if "whole_FET" not in df.columns:
            continue
        scaf_cols = [c for c in df.columns if c.endswith("_scaf")]
        sps = [c[:-5] for c in scaf_cols]
        if len(sps) != 2:
            continue
        sa, sb = sps
        sca, scb = "{}_scaf".format(sa), "{}_scaf".format(sb)
        sig = df[df["whole_FET"] <= 0.05]
        if sig.empty:
            continue
        grp = sig.groupby([sca, scb]).size().reset_index(name="count")
        for _, r in grp.iterrows():
            na = (sa, r[sca]); nb = (sb, r[scb])
            ca = raw_node_to_cluster.get(na)
            cb = raw_node_to_cluster.get(nb)
            if ca is None or cb is None or ca == cb:
                continue
            key = (min(ca, cb), max(ca, cb))
            cross[key] += int(r["count"])

    # Greedy nearest-neighbor cluster ordering: start from the largest
    # cluster, append the unselected cluster with the strongest edge
    # to the most-recently-placed cluster (cluster-weight tiebreak).
    # This puts clusters that share fusion/fission links next to each
    # other in the final plot.
    cluster_ids = list(range(len(clusters)))
    cw = {i: cluster_weight(c) for i, c in enumerate(clusters)}
    if cluster_ids:
        ordered_ids = [max(cluster_ids, key=lambda i: cw[i])]
        remaining = set(cluster_ids) - {ordered_ids[0]}
        while remaining:
            last = ordered_ids[-1]
            def score(j):
                key = (min(last, j), max(last, j))
                return (cross.get(key, 0), cw[j])
            nxt = max(remaining, key=score)
            ordered_ids.append(nxt)
            remaining.discard(nxt)
        # remap to dense ids in plot order
        remap = {old: new for new, old in enumerate(ordered_ids)}
    else:
        remap = {}

    node_to_cluster = {n: remap[c] for n, c in raw_node_to_cluster.items()}
    cluster_size = {remap[i]: len(clusters[i]) for i in range(len(clusters))}
    cluster_wt   = {remap[i]: cw[i] for i in range(len(clusters))}
    return node_to_cluster, cluster_size, cluster_wt


def _order_by_lg(species_order, rbh_df_list, sp_to_chr_to_size,
                 sp_min_chr_size, sp_to_gene_order, sp_pair_to_rbh_df_list_index):
    """Per-species scaffold ordering built from a global FET-graph LG
    assignment (_build_fet_lg_clusters) PLUS a within-cluster cascade
    that uses FET-best-partner so siblings sharing a downstream chrom
    stay adjacent. Cluster id is the primary sort key (so all members
    of one synteny linkage group share an x-region across every row);
    within a cluster, ordering is the same FET-best-partner heuristic
    the row-by-row cascade uses, but restricted to that cluster's
    chromosomes so it can't pull partners from a different LG.

    Scaffolds with no FET-significant edges fall to the end, sorted
    by chrom size descending.
    """
    node_to_cluster, cluster_size, cluster_weight = \
        _build_fet_lg_clusters(species_order, rbh_df_list)
    cluster_ids_sorted = sorted(cluster_size.keys())

    def all_scafs_of(sp):
        out = set()
        for df in rbh_df_list:
            sc = "{}_scaf".format(sp)
            if sc in df.columns:
                out.update(df[sc].unique())
        return out

    sp_to_chromorder = {}
    for i, sp in enumerate(species_order):
        # bucket this species' scaffolds by cluster id
        by_cluster = {cid: [] for cid in cluster_ids_sorted}
        untagged = []
        for s in all_scafs_of(sp):
            cid = node_to_cluster.get((sp, s))
            if cid is None:
                untagged.append(s)
            else:
                by_cluster[cid].append(s)

        # for each cluster, run the FET-best-partner cascade restricted
        # to that cluster's chromosomes. The top row uses row 1 as the
        # anchor; every other row uses the row directly above (already
        # ordered).
        ordered_in_clusters = []
        for cid in cluster_ids_sorted:
            members = set(by_cluster[cid])
            if not members:
                continue
            if i == 0 and len(species_order) > 1:
                ref_sp = species_order[1]
                ref_members = {s for s in all_scafs_of(ref_sp)
                               if node_to_cluster.get((ref_sp, s)) == cid}
                sub = _seed_top_order_within(
                    sp, ref_sp, rbh_df_list[0], members, ref_members)
            elif i == 0:
                sub = sorted(members)
            else:
                prev_sp = species_order[i - 1]
                lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([sp, prev_sp]))]
                prev_members = {s for s in all_scafs_of(prev_sp)
                                if node_to_cluster.get((prev_sp, s)) == cid}
                sub = _cascade_within(
                    sp, prev_sp, rbh_df_list[lookup],
                    members, prev_members, sp_to_chromorder)
            ordered_in_clusters.extend(sub)

        # untagged scaffolds sit at the tail, sorted by chromsize desc
        sized = sp_to_chr_to_size.get(sp, {})
        untagged.sort(key=lambda s: -sized.get(s, 0))
        ordered = ordered_in_clusters + untagged
        ordered = _quality_check_chromosome_list(
            sp, ordered, sp_to_chr_to_size,
            sp_to_gene_order, sp_min_chr_size)
        sp_to_chromorder[sp] = {s: k for k, s in enumerate(ordered)}
    return sp_to_chromorder


def _seed_top_order_within(sp, partner_sp, rbhdf, members, partner_members):
    """Within-cluster top-row seeding. Same idea as
    _seed_top_order_from_partner but restricted to a specific subset
    of `sp` and `partner_sp` chromosomes (one cluster). Groups sp
    members by best partner_sp partner; orders partner groups by total
    shared ortholog count desc; within a group, by FET asc."""
    sa = "{}_scaf".format(sp)
    sb = "{}_scaf".format(partner_sp)
    df = rbhdf[rbhdf[sa].isin(members) & rbhdf[sb].isin(partner_members)]
    if df.empty:
        return sorted(members)
    counts = df.groupby([sa, sb]).size().reset_index(name="_count")
    fets   = df.groupby([sa, sb])["whole_FET"].min().reset_index() \
             if "whole_FET" in df.columns else None
    if fets is None:
        counts["whole_FET"] = 0.0
        pairs = counts
    else:
        pairs = counts.merge(fets, on=[sa, sb])
    best = pairs.sort_values([sa, "whole_FET"]).drop_duplicates(subset=[sa])
    partner_totals = pairs.groupby(sb)["_count"].sum().sort_values(ascending=False)
    partner_pos = {p: i for i, p in enumerate(partner_totals.index)}
    best["_pord"] = best[sb].map(partner_pos)
    best = best.sort_values(["_pord", "whole_FET"])
    ordered = best[sa].tolist()
    leftover = [s for s in members if s not in set(ordered)]
    return ordered + leftover


def _cascade_within(sp, prev_sp, rbhdf, members, prev_members, sp_to_chromorder):
    """Within-cluster cascade for non-top rows. For each sp chromosome
    in `members`, pick its FET-best partner among prev_sp's
    `prev_members`; order sp members by that partner's position in
    sp_to_chromorder[prev_sp], tiebreak by FET. Mirrors
    _optimize_spA_based_on_rbh, but restricted to one cluster on each
    side so it cannot drag siblings across LG boundaries."""
    sa = "{}_scaf".format(sp)
    sb = "{}_scaf".format(prev_sp)
    df = rbhdf[rbhdf[sa].isin(members) & rbhdf[sb].isin(prev_members)]
    if df.empty:
        return sorted(members)
    df = df[[sa, sb, "whole_FET"]].copy()
    df = df.sort_values([sa, "whole_FET"]).drop_duplicates(subset=[sa])
    df["_pos"] = df[sb].map(sp_to_chromorder[prev_sp])
    df = df.dropna(subset=["_pos"])
    df = df.sort_values(["_pos", "whole_FET"]).reset_index(drop=True)
    ordered = df[sa].tolist()
    leftover = [s for s in members if s not in set(ordered)]
    return ordered + leftover


def _seed_top_order_from_partner(sp, partner_sp, rbhdf):
    """Top-row ordering that respects FET-significant chromosome
    relationships with the second species. Chromosomes of the top
    species that share the same best-FET partner in row 2 get grouped
    together; partner groups are ordered by total shared ortholog
    count desc; within a group, members are sorted by FET (most
    significant first).

    Solves the case (odp issue #127 review) where the previous top-row
    ordering, plain value_counts(), would separate two top-row
    chromosomes that both map to the same row-2 chromosome — leading
    to forced crossings the cascade below cannot fix."""
    sa = "{}_scaf".format(sp)
    sb = "{}_scaf".format(partner_sp)
    if rbhdf.empty or sa not in rbhdf.columns or sb not in rbhdf.columns:
        return list(rbhdf[sa].value_counts().index) if sa in rbhdf.columns else []
    counts = rbhdf.groupby([sa, sb]).size().reset_index(name="_count")
    fets   = rbhdf.groupby([sa, sb])["whole_FET"].min().reset_index() \
             if "whole_FET" in rbhdf.columns else None
    if fets is None:
        counts["whole_FET"] = 0.0
        pairs = counts
    else:
        pairs = counts.merge(fets, on=[sa, sb])
    best = pairs.sort_values([sa, "whole_FET"]).drop_duplicates(subset=[sa])
    partner_totals = pairs.groupby(sb)["_count"].sum().sort_values(ascending=False)
    partner_pos = {p: i for i, p in enumerate(partner_totals.index)}
    best["_pord"] = best[sb].map(partner_pos)
    best = best.sort_values(["_pord", "whole_FET"])
    ordered = best[sa].tolist()
    # any scaffold that didn't get a partner row keeps a deterministic
    # tail position so we don't drop it from the plot.
    leftover = [s for s in rbhdf[sa].unique() if s not in set(ordered)]
    return ordered + leftover


def _scaffold_anchor_positions(rbhdf, this_scaf_col, ref_scaf_col,
                               ref_chromorder, this_chromIx_col=None,
                               ref_chromIx_col=None):
    """Return a dict scaf -> sorted list of float positions of this
    scaffold's orthologs as seen on the reference (already-ordered)
    species, where a float position is (ref_chrom_slot + within-chrom
    fraction). If chromIx columns are absent, falls back to using the
    chrom slot only (no within-chrom resolution)."""
    sub = rbhdf[rbhdf[ref_scaf_col].isin(ref_chromorder)].copy()
    if sub.empty:
        return {}
    sub["_ref_slot"] = sub[ref_scaf_col].map(ref_chromorder).astype(float)
    if ref_chromIx_col is not None and ref_chromIx_col in sub.columns:
        sizes = sub.groupby(ref_scaf_col)[ref_chromIx_col].transform("max").replace(0, 1) + 1
        sub["_ref_pos"] = sub["_ref_slot"] + sub[ref_chromIx_col] / sizes
    else:
        sub["_ref_pos"] = sub["_ref_slot"]
    out = {}
    for scaf, grp in sub.groupby(this_scaf_col, sort=False):
        out[scaf] = grp["_ref_pos"].tolist()
    return out


def _detangle_one_side(sp, ref_sp, rbhdf, sp_to_chromorder, sp_to_genesdfs,
                       method, kept_scafs):
    """Decide the order of *sp*'s scaffolds using anchor positions on
    *ref_sp* (already ordered). Returns a list of sp's scaffolds in the
    chosen order. `method` is one of:
      barycenter    - sort by mean ref position
      median        - sort by median ref position
    `kept_scafs` restricts the result to that set (used for parity with
    _quality_check_chromosome_list, which filters by minscafsize)."""
    this_scaf_col = "{}_scaf".format(sp)
    ref_scaf_col  = "{}_scaf".format(ref_sp)
    ref_chromIx_col = "{}_chromIx".format(ref_sp) if "{}_chromIx".format(ref_sp) in rbhdf.columns else None
    df = rbhdf[[this_scaf_col, ref_scaf_col]].copy()
    if ref_chromIx_col is None:
        # Recover ref-side rank from the genesdf if present.
        gdf = sp_to_genesdfs.get(ref_sp)
        if gdf is not None and "{}_chromIx".format(ref_sp) in gdf.columns:
            gene_col = "{}_gene".format(ref_sp)
            if gene_col in rbhdf.columns:
                df[gene_col] = rbhdf[gene_col]
                df = df.merge(gdf[[gene_col, "{}_chromIx".format(ref_sp)]],
                              on=gene_col, how="left")
                ref_chromIx_col = "{}_chromIx".format(ref_sp)
    anchors = _scaffold_anchor_positions(
        df, this_scaf_col, ref_scaf_col, sp_to_chromorder[ref_sp],
        ref_chromIx_col=ref_chromIx_col)

    import statistics
    scores = {}
    for scaf, positions in anchors.items():
        if not positions:
            continue
        if method == "median":
            scores[scaf] = statistics.median(positions)
        else:  # barycenter / mean
            scores[scaf] = sum(positions) / len(positions)
    # scaffolds with no orthologs on the reference get pushed to the end,
    # preserving their relative order among themselves.
    ordered  = sorted(scores.keys(), key=lambda s: (scores[s], s))
    leftover = [s for s in rbhdf[this_scaf_col].unique() if s not in scores]
    candidates = ordered + leftover
    return [s for s in candidates if s in kept_scafs]


def _crossings_between(sp_top, sp_bot, rbhdf, sp_to_chromorder,
                       sp_to_genesdfs):
    """Total bezier-line crossings between two species given their
    current sp_to_chromorder. Inversion count between top-x and bottom-x
    ranks of every line, ordering within a chromosome resolved by the
    gene-rank (chromIx) so two lines on the same top scaffold don't tie."""
    tcol = "{}_scaf".format(sp_top); bcol = "{}_scaf".format(sp_bot)
    tg   = "{}_gene".format(sp_top); bg   = "{}_gene".format(sp_bot)
    df = rbhdf[[tcol, bcol, tg, bg]].copy()
    df = df[df[tcol].isin(sp_to_chromorder[sp_top])
            & df[bcol].isin(sp_to_chromorder[sp_bot])]
    if df.empty:
        return 0
    t_rank = dict(zip(sp_to_genesdfs[sp_top][tg],
                      sp_to_genesdfs[sp_top]["{}_chromIx".format(sp_top)]))
    b_rank = dict(zip(sp_to_genesdfs[sp_bot][bg],
                      sp_to_genesdfs[sp_bot]["{}_chromIx".format(sp_bot)]))
    t_sz = sp_to_genesdfs[sp_top].groupby(tcol).size().to_dict()
    b_sz = sp_to_genesdfs[sp_bot].groupby(bcol).size().to_dict()
    df["_tr"] = df[tg].map(t_rank); df["_br"] = df[bg].map(b_rank)
    df = df.dropna(subset=["_tr", "_br"])
    df["_tp"] = df[tcol].map(sp_to_chromorder[sp_top]).astype(float) \
                + df["_tr"] / df[tcol].map(t_sz)
    df["_bp"] = df[bcol].map(sp_to_chromorder[sp_bot]).astype(float) \
                + df["_br"] / df[bcol].map(b_sz)
    df = df.sort_values("_tp", kind="mergesort")
    return _count_inversions(df["_bp"].tolist())


def _greedy_swap_polish(species_order, rbh_df_list,
                        sp_pair_to_rbh_df_list_index,
                        sp_to_chromorder, sp_to_genesdfs,
                        max_passes=6):
    """Adjacent-swap polish over every species. For each row, try
    swapping every pair of neighboring chromosomes; if the sum of
    crossings against the row above and the row below drops, keep the
    swap. Repeat until a pass yields no improvement or max_passes."""
    def pair_crossings(sp, neighbor):
        if neighbor is None:
            return 0
        lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([sp, neighbor]))]
        return _crossings_between(neighbor, sp, rbh_df_list[lookup],
                                  sp_to_chromorder, sp_to_genesdfs) \
               if neighbor == species_order[species_order.index(sp) - 1] \
               else _crossings_between(sp, neighbor, rbh_df_list[lookup],
                                       sp_to_chromorder, sp_to_genesdfs)

    for _ in range(max_passes):
        improved = False
        for i, sp in enumerate(species_order):
            above = species_order[i - 1] if i > 0 else None
            below = species_order[i + 1] if i + 1 < len(species_order) else None
            order = sorted(sp_to_chromorder[sp].keys(),
                           key=lambda s: sp_to_chromorder[sp][s])
            j = 0
            while j < len(order) - 1:
                a, b = order[j], order[j + 1]
                cur = pair_crossings(sp, above) + pair_crossings(sp, below)
                sp_to_chromorder[sp][a] = j + 1
                sp_to_chromorder[sp][b] = j
                new = pair_crossings(sp, above) + pair_crossings(sp, below)
                if new < cur:
                    order[j], order[j + 1] = b, a
                    improved = True
                else:
                    sp_to_chromorder[sp][a] = j
                    sp_to_chromorder[sp][b] = j + 1
                j += 1
        if not improved:
            break
    return sp_to_chromorder


def _fenwick_weighted_inversion(values, weights):
    """Weighted inversion count. For every pair (i < j) with values[i] >
    values[j], adds weights[i] * weights[j] to the total. O(n log n)
    using a Fenwick tree over the unique values seen in `values`. This
    is the scoring kernel for the crossing-count local search: each
    bezier-line pair is one potential crossing, weighted by the
    importance of the two lines."""
    if not values:
        return 0.0
    sorted_unique = sorted(set(values))
    rank = {v: i + 1 for i, v in enumerate(sorted_unique)}
    n = len(sorted_unique)
    fen = [0.0] * (n + 1)

    def update(i, w):
        while i <= n:
            fen[i] += w
            i += i & -i

    def query(i):
        s = 0.0
        while i > 0:
            s += fen[i]
            i -= i & -i
        return s

    inv = 0.0
    total = 0.0
    for v, w in zip(values, weights):
        r = rank[v]
        gt = total - query(r)
        inv += w * gt
        update(r, w)
        total += w
    return inv


def _build_pair_anchor_data(species_order, rbh_df_list,
                            sp_pair_to_rbh_df_list_index, sp_to_genesdfs,
                            fet_weight=1000.0):
    """Coarsened version of _build_pair_line_data: every distinct
    (top_chrom, bot_chrom) pair becomes ONE record at the
    within-chrom mean rank on each side, weight = (fet_weight if the
    pair's best whole_FET <= 0.05 else 1) * shared-ortholog count.
    Drops the line count from ~2000 to ~200 per pair on Lep-style data
    and brings the local search runtime down by ~10x. Resolution
    loss: we no longer see within-chrom inversions caused by individual
    gene-pair crossings, which the search couldn't fix anyway (the
    moves are chrom flips and chrom swaps, not gene reorderings)."""
    out = []
    for k in range(len(species_order) - 1):
        top_sp = species_order[k]
        bot_sp = species_order[k + 1]
        lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([top_sp, bot_sp]))]
        df = rbh_df_list[lookup]
        tcol = "{}_scaf".format(top_sp); bcol = "{}_scaf".format(bot_sp)
        tg = "{}_gene".format(top_sp); bg = "{}_gene".format(bot_sp)
        t_rank = dict(zip(sp_to_genesdfs[top_sp][tg],
                          sp_to_genesdfs[top_sp]["{}_chromIx".format(top_sp)]))
        b_rank = dict(zip(sp_to_genesdfs[bot_sp][bg],
                          sp_to_genesdfs[bot_sp]["{}_chromIx".format(bot_sp)]))
        t_sz = sp_to_genesdfs[top_sp].groupby(tcol).size().to_dict()
        b_sz = sp_to_genesdfs[bot_sp].groupby(bcol).size().to_dict()
        tmp = df.copy()
        tmp["_trk"] = tmp[tg].map(t_rank)
        tmp["_brk"] = tmp[bg].map(b_rank)
        tmp = tmp.dropna(subset=["_trk", "_brk"])
        if "whole_FET" not in tmp.columns:
            tmp["whole_FET"] = 1.0
        gr = tmp.groupby([tcol, bcol]).agg(
            count=("whole_FET", "size"),
            best_fet=("whole_FET", "min"),
            mean_trk=("_trk", "mean"),
            mean_brk=("_brk", "mean")).reset_index()
        records = []
        for tchrom, bchrom, count, best_fet, mtrk, mbrk in zip(
                gr[tcol].values, gr[bcol].values,
                gr["count"].values, gr["best_fet"].values,
                gr["mean_trk"].values, gr["mean_brk"].values):
            w_unit = fet_weight if best_fet <= 0.05 else 1.0
            records.append((tchrom, bchrom, float(mtrk), float(mbrk),
                            int(t_sz[tchrom]), int(b_sz[bchrom]),
                            float(w_unit * count)))
        out.append(records)
    return out


def _build_pair_line_data(species_order, rbh_df_list,
                          sp_pair_to_rbh_df_list_index, sp_to_genesdfs,
                          fet_weight=1000.0):
    """For every adjacent species pair, precompute per-line records
    (chrom, within-chrom rank, chrom size, weight). Used by the local
    search so we don't re-parse the RBH dataframes for each candidate
    move. fet_weight controls how much FET-significant lines outscore
    non-significant ones — the user-specified 1000x makes the search
    effectively minimize crossings of the bold ribbons first, only
    using the faint ones as a tiebreak."""
    line_data_per_pair = []
    for k in range(len(species_order) - 1):
        top_sp = species_order[k]
        bot_sp = species_order[k + 1]
        lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([top_sp, bot_sp]))]
        df = rbh_df_list[lookup]
        tcol = "{}_scaf".format(top_sp); bcol = "{}_scaf".format(bot_sp)
        tg = "{}_gene".format(top_sp); bg = "{}_gene".format(bot_sp)
        t_rank = dict(zip(sp_to_genesdfs[top_sp][tg],
                          sp_to_genesdfs[top_sp]["{}_chromIx".format(top_sp)]))
        b_rank = dict(zip(sp_to_genesdfs[bot_sp][bg],
                          sp_to_genesdfs[bot_sp]["{}_chromIx".format(bot_sp)]))
        t_sz = sp_to_genesdfs[top_sp].groupby(tcol).size().to_dict()
        b_sz = sp_to_genesdfs[bot_sp].groupby(bcol).size().to_dict()
        records = []
        for tchrom, bchrom, trk, brk, fet in zip(
                df[tcol].values, df[bcol].values,
                df[tg].map(t_rank).values,
                df[bg].map(b_rank).values,
                df.get("whole_FET", pd.Series([1.0] * len(df))).values):
            if trk is None or brk is None or pd.isna(trk) or pd.isna(brk):
                continue
            w = fet_weight if fet <= 0.05 else 1.0
            records.append((tchrom, bchrom, int(trk), int(brk),
                            int(t_sz[tchrom]), int(b_sz[bchrom]), w))
        line_data_per_pair.append(records)
    return line_data_per_pair


def _score_pair_lines(records, top_order, bot_order, top_flip, bot_flip):
    """Weighted crossing count for one adjacent species pair given
    current chrom orders and per-chrom flips. Lines whose endpoints
    sit on chroms not in `top_order`/`bot_order` are ignored
    (typical when a scaffold is filtered out by minscafsize)."""
    vals_top = []
    vals_bot = []
    weights = []
    for tchrom, bchrom, trk, brk, tsz, bsz, w in records:
        if tchrom not in top_order or bchrom not in bot_order:
            continue
        t_eff = (tsz - 1 - trk) if top_flip.get(tchrom, False) else trk
        b_eff = (bsz - 1 - brk) if bot_flip.get(bchrom, False) else brk
        vals_top.append(top_order[tchrom] + t_eff / tsz)
        vals_bot.append(bot_order[bchrom] + b_eff / bsz)
        weights.append(w)
    if not vals_top:
        return 0.0
    order = sorted(range(len(vals_top)), key=lambda i: vals_top[i])
    sorted_bot = [vals_bot[i] for i in order]
    sorted_w   = [weights[i]  for i in order]
    return _fenwick_weighted_inversion(sorted_bot, sorted_w)


def _brushing_sweep(species_order, rbh_df_list,
                    sp_pair_to_rbh_df_list_index,
                    sp_to_chromorder, sp_to_chrom_flip,
                    sp_to_genesdfs, max_iters=50, verbose=False,
                    snapshot_fn=None, patience=4,
                    top_moves_per_iter=15,
                    top_move_offsets=(1, 2, 3, 5),
                    enable_bottom_up=True,
                    bottom_moves_per_iter=15):
    """Modified Sugiyama: FET-weighted barycenter sweeps alternating
    top-down and bottom-up. After each direction the per-scaffold
    flip cascade is re-run so chromosome orientation is re-checked
    in the new positional context. The total bezier-crossing count
    is logged after each direction, with flush so the user can see
    convergence live.

    Only FET-significant orthologs (whole_FET <= 0.05) contribute to
    each barycenter; trace homology is ignored. A chromosome with no
    significant partners keeps its existing position.

    Best-seen state is tracked because the two directions can
    disagree on a chromosome's ideal position (top says near anchor
    X, bottom says near anchor Y). The function returns whichever
    state had the lowest crossing count across all sweeps.

    snapshot_fn, if supplied, is called as
        snapshot_fn(gen, direction, sp_to_chromorder, sp_to_chrom_flip, crossings)
    after each direction. Callers use this to render a PDF per brush
    stroke so the convergence is visually inspectable.
    """
    line_data = _build_pair_anchor_data(
        species_order, rbh_df_list, sp_pair_to_rbh_df_list_index,
        sp_to_genesdfs, fet_weight=1000.0)

    def total_crossings():
        total = 0.0
        for k in range(len(species_order) - 1):
            total += _score_pair_lines(
                line_data[k],
                sp_to_chromorder[species_order[k]],
                sp_to_chromorder[species_order[k + 1]],
                sp_to_chrom_flip[species_order[k]],
                sp_to_chrom_flip[species_order[k + 1]])
        return total

    def barycenter_position(sp, anchor_sp, lookup_idx):
        df = rbh_df_list[lookup_idx]
        if "whole_FET" not in df.columns:
            return {}
        sa = "{}_scaf".format(sp); sb = "{}_scaf".format(anchor_sp)
        ag = "{}_gene".format(anchor_sp)
        sig = df[(df["whole_FET"] <= 0.05)
                 & df[sa].isin(sp_to_chromorder[sp])
                 & df[sb].isin(sp_to_chromorder[anchor_sp])]
        if sig.empty:
            return {}
        anchor_rank = dict(zip(
            sp_to_genesdfs[anchor_sp][ag],
            sp_to_genesdfs[anchor_sp]["{}_chromIx".format(anchor_sp)]))
        anchor_sz = sp_to_genesdfs[anchor_sp].groupby(sb).size().to_dict()
        anchor_flip = sp_to_chrom_flip[anchor_sp]
        anchor_pos_lookup = sp_to_chromorder[anchor_sp]
        rows = sig.assign(
            _ank=sig[ag].map(anchor_rank),
            _asz=sig[sb].map(anchor_sz).astype(float),
            _aslot=sig[sb].map(anchor_pos_lookup).astype(float))
        rows = rows.dropna(subset=["_ank", "_asz", "_aslot"])
        if rows.empty:
            return {}
        rows["_aeff"] = rows.apply(
            lambda r: (r["_asz"] - 1 - r["_ank"])
                      if anchor_flip.get(r[sb], False) else r["_ank"], axis=1)
        rows["_apos"] = rows["_aslot"] + rows["_aeff"] / rows["_asz"]
        return rows.groupby(sa)["_apos"].mean().to_dict()

    initial = total_crossings()
    if verbose:
        print("brushing seed crossings: {:.0f}".format(initial), flush=True)
    if snapshot_fn is not None:
        snapshot_fn(0, "seed", sp_to_chromorder, sp_to_chrom_flip, initial)

    best_score  = initial
    best_orders = {sp: dict(sp_to_chromorder[sp]) for sp in species_order}
    best_flips  = {sp: dict(sp_to_chrom_flip[sp]) for sp in species_order}

    def snapshot_if_best(score):
        nonlocal best_score, best_orders, best_flips
        if score < best_score:
            best_score = score
            best_orders = {sp: dict(sp_to_chromorder[sp]) for sp in species_order}
            best_flips  = {sp: dict(sp_to_chrom_flip[sp]) for sp in species_order}

    prev = initial
    stagnation = 0

    def top_row_tension():
        """Per-top-row-chrom score: sum of weights of lines emerging
        from each top chrom into row 1, multiplied by the number of
        crossings each line participates in. Heavy crossings on
        bezier paths leaving a particular top chrom show up here -
        the "knot" anchored at that chrom."""
        top_sp = species_order[0]
        recs = line_data[0]
        top_order = sp_to_chromorder[top_sp]
        bot_order = sp_to_chromorder[species_order[1]]
        # gather (top_chrom, t_pos, b_pos, w)
        items = []
        for tchrom, bchrom, trk, brk, tsz, bsz, w in recs:
            if tchrom not in top_order or bchrom not in bot_order:
                continue
            t_eff = (tsz - 1 - trk) if sp_to_chrom_flip[top_sp].get(tchrom, False) else trk
            b_eff = (bsz - 1 - brk) if sp_to_chrom_flip[species_order[1]].get(bchrom, False) else brk
            items.append((tchrom, top_order[tchrom] + t_eff/tsz,
                          bot_order[bchrom] + b_eff/bsz, w))
        items.sort(key=lambda x: x[1])
        # for each item, count crossings = pairs (j<i with bot_j > bot_i)
        tension = {tc: 0.0 for tc in top_order}
        # naive O(L^2) is fine here (anchor data ~200 lines per pair)
        for i in range(len(items)):
            ti_c, ti_t, ti_b, ti_w = items[i]
            for j in range(i):
                tj_c, tj_t, tj_b, tj_w = items[j]
                if tj_b > ti_b:
                    tension[ti_c] += ti_w * tj_w
                    tension[tj_c] += ti_w * tj_w
        return tension

    def try_top_moves():
        """Greedy moves on the top row, prioritised by tension. Returns
        True if any move was accepted."""
        if top_moves_per_iter <= 0:
            return False
        top_sp = species_order[0]
        tension = top_row_tension()
        order_by_tension = sorted(tension, key=lambda c: -tension[c])
        top_targets = order_by_tension[:top_moves_per_iter]
        cur_score = total_crossings()
        any_move = False
        for chrom in top_targets:
            # 1) flip
            sp_to_chrom_flip[top_sp][chrom] = not sp_to_chrom_flip[top_sp][chrom]
            new_score = total_crossings()
            if new_score < cur_score - 1e-9:
                cur_score = new_score
                any_move = True
            else:
                sp_to_chrom_flip[top_sp][chrom] = not sp_to_chrom_flip[top_sp][chrom]
            # 2) swap with chrom at offset (1,2,3,5)
            for off in top_move_offsets:
                cur_pos = sp_to_chromorder[top_sp][chrom]
                other = next((c for c, p in sp_to_chromorder[top_sp].items()
                              if p == cur_pos + off), None)
                if other is None:
                    continue
                sp_to_chromorder[top_sp][chrom], sp_to_chromorder[top_sp][other] = \
                    sp_to_chromorder[top_sp][other], sp_to_chromorder[top_sp][chrom]
                new_score = total_crossings()
                if new_score < cur_score - 1e-9:
                    cur_score = new_score
                    any_move = True
                else:
                    sp_to_chromorder[top_sp][chrom], sp_to_chromorder[top_sp][other] = \
                        sp_to_chromorder[top_sp][other], sp_to_chromorder[top_sp][chrom]
        return any_move

    def bottom_row_tension():
        """Symmetric to top_row_tension but for the last species row,
        scoring tension from the last adjacent pair (N-2, N-1)."""
        if len(species_order) < 2:
            return {}
        bot_sp = species_order[-1]
        recs = line_data[-1]
        top_order = sp_to_chromorder[species_order[-2]]
        bot_order = sp_to_chromorder[bot_sp]
        items = []
        for tchrom, bchrom, trk, brk, tsz, bsz, w in recs:
            if tchrom not in top_order or bchrom not in bot_order:
                continue
            t_eff = (tsz - 1 - trk) if sp_to_chrom_flip[species_order[-2]].get(tchrom, False) else trk
            b_eff = (bsz - 1 - brk) if sp_to_chrom_flip[bot_sp].get(bchrom, False) else brk
            items.append((bchrom, top_order[tchrom] + t_eff / tsz,
                          bot_order[bchrom] + b_eff / bsz, w))
        items.sort(key=lambda x: x[1])
        tension = {bc: 0.0 for bc in bot_order}
        for i in range(len(items)):
            bi_c, bi_t, bi_b, bi_w = items[i]
            for j in range(i):
                bj_c, bj_t, bj_b, bj_w = items[j]
                if bj_b > bi_b:
                    tension[bi_c] += bi_w * bj_w
                    tension[bj_c] += bi_w * bj_w
        return tension

    def try_bottom_moves():
        """Symmetric to try_top_moves: greedy moves on the last row
        prioritised by tension from the last pair. Helps after a bu
        sweep the same way try_top_moves helps after a td sweep."""
        if bottom_moves_per_iter <= 0 or len(species_order) < 2:
            return False
        bot_sp = species_order[-1]
        tension = bottom_row_tension()
        order_by_tension = sorted(tension, key=lambda c: -tension[c])
        targets = order_by_tension[:bottom_moves_per_iter]
        cur_score = total_crossings()
        any_move = False
        for chrom in targets:
            sp_to_chrom_flip[bot_sp][chrom] = not sp_to_chrom_flip[bot_sp][chrom]
            new_score = total_crossings()
            if new_score < cur_score - 1e-9:
                cur_score = new_score
                any_move = True
            else:
                sp_to_chrom_flip[bot_sp][chrom] = not sp_to_chrom_flip[bot_sp][chrom]
            for off in top_move_offsets:
                cur_pos = sp_to_chromorder[bot_sp][chrom]
                other = next((c for c, p in sp_to_chromorder[bot_sp].items()
                              if p == cur_pos + off), None)
                if other is None:
                    continue
                sp_to_chromorder[bot_sp][chrom], sp_to_chromorder[bot_sp][other] = \
                    sp_to_chromorder[bot_sp][other], sp_to_chromorder[bot_sp][chrom]
                new_score = total_crossings()
                if new_score < cur_score - 1e-9:
                    cur_score = new_score
                    any_move = True
                else:
                    sp_to_chromorder[bot_sp][chrom], sp_to_chromorder[bot_sp][other] = \
                        sp_to_chromorder[bot_sp][other], sp_to_chromorder[bot_sp][chrom]
        return any_move

    for it in range(max_iters):
        # top-down brush
        for i in range(1, len(species_order)):
            sp = species_order[i]
            anchor = species_order[i - 1]
            lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([sp, anchor]))]
            pos = barycenter_position(sp, anchor, lookup)
            if not pos:
                continue
            chroms = list(sp_to_chromorder[sp].keys())
            chroms.sort(key=lambda c: (pos.get(c, sp_to_chromorder[sp][c]),
                                       sp_to_chromorder[sp][c]))
            sp_to_chromorder[sp] = {c: k for k, c in enumerate(chroms)}
        sp_to_chrom_flip = _optimize_chrom_flips_top_down(
            species_order, rbh_df_list,
            sp_pair_to_rbh_df_list_index,
            sp_to_chromorder, sp_to_genesdfs,
            initial_flip=sp_to_chrom_flip)
        td = total_crossings()
        if verbose:
            print("  gen {} top-down  crossings: {:.0f}  (delta {:+.0f})".format(
                it + 1, td, td - prev), flush=True)
        improved_best = td < best_score
        snapshot_if_best(td)
        if snapshot_fn is not None:
            snapshot_fn(it + 1, "td", sp_to_chromorder, sp_to_chrom_flip, td)
        stagnation = 0 if improved_best else stagnation + 1
        prev = td

        # Downstream-feedback: after the comb stroke settles every row
        # below the top, the remaining tangles reveal where the *top*
        # row was wrong. Try a small set of moves on the highest-
        # tension top-row chromosomes; accept any that drop the global
        # score. The next td sweep starts from this corrected top.
        if try_top_moves():
            ts = total_crossings()
            if verbose:
                print("  gen {} top-moves crossings: {:.0f}  (delta {:+.0f})".format(
                    it + 1, ts, ts - prev), flush=True)
            improved_best = ts < best_score
            snapshot_if_best(ts)
            if snapshot_fn is not None:
                snapshot_fn(it + 1, "tm", sp_to_chromorder, sp_to_chrom_flip, ts)
            stagnation = 0 if improved_best else stagnation + 1
            prev = ts

        if not enable_bottom_up:
            # Patience check based on top-down + top-moves cycles only.
            if stagnation >= patience:
                if verbose:
                    print("  brushing converged at gen {} after {} stagnant sweeps"
                          " (best {:.0f})".format(it + 1, patience, best_score),
                          flush=True)
                break
            continue

        # bottom-up brush
        for i in range(len(species_order) - 2, -1, -1):
            sp = species_order[i]
            anchor = species_order[i + 1]
            lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([sp, anchor]))]
            pos = barycenter_position(sp, anchor, lookup)
            if not pos:
                continue
            chroms = list(sp_to_chromorder[sp].keys())
            chroms.sort(key=lambda c: (pos.get(c, sp_to_chromorder[sp][c]),
                                       sp_to_chromorder[sp][c]))
            sp_to_chromorder[sp] = {c: k for k, c in enumerate(chroms)}
        sp_to_chrom_flip = _optimize_chrom_flips_top_down(
            species_order, rbh_df_list,
            sp_pair_to_rbh_df_list_index,
            sp_to_chromorder, sp_to_genesdfs,
            initial_flip=sp_to_chrom_flip)
        bu = total_crossings()
        if verbose:
            print("  gen {} bottom-up crossings: {:.0f}  (delta {:+.0f})".format(
                it + 1, bu, bu - prev), flush=True)
        improved_best = bu < best_score
        snapshot_if_best(bu)
        if snapshot_fn is not None:
            snapshot_fn(it + 1, "bu", sp_to_chromorder, sp_to_chrom_flip, bu)
        stagnation = 0 if improved_best else stagnation + 1
        prev = bu

        # Symmetric bottom-row feedback after bu sweep. Tries moves on
        # the last species row's chroms ordered by tension from the
        # last adjacent pair. The next bu sweep starts from this
        # corrected bottom.
        if try_bottom_moves():
            bs = total_crossings()
            if verbose:
                print("  gen {} bot-moves crossings: {:.0f}  (delta {:+.0f})".format(
                    it + 1, bs, bs - prev), flush=True)
            improved_best = bs < best_score
            snapshot_if_best(bs)
            if snapshot_fn is not None:
                snapshot_fn(it + 1, "bm", sp_to_chromorder, sp_to_chrom_flip, bs)
            stagnation = 0 if improved_best else stagnation + 1
            prev = bs

        if stagnation >= patience:
            if verbose:
                print("  brushing converged at gen {} after {} stagnant sweeps"
                      " (best {:.0f})".format(it + 1, patience, best_score),
                      flush=True)
            break

    for sp in species_order:
        sp_to_chromorder[sp] = best_orders[sp]
        sp_to_chrom_flip[sp] = best_flips[sp]
    if verbose:
        print("brushing total drop: {:.0f} -> {:.0f}  ({:+.2f}%)".format(
            initial, best_score,
            100.0 * (best_score - initial) / initial), flush=True)
    return sp_to_chromorder, sp_to_chrom_flip


def _crossing_local_search(species_order, rbh_df_list,
                           sp_pair_to_rbh_df_list_index,
                           sp_to_chromorder, sp_to_chrom_flip,
                           sp_to_genesdfs, fet_weight=1000.0,
                           max_passes=50, coarse=True,
                           early_stop_frac=0.0005, verbose=False):
    """Greedy improvement-monotone local search on the *true* total
    weighted crossing count across every adjacent species pair. Moves
    tried per chromosome: F/R flip, swap with the next chrom in the
    same row. Per-pass species traversal order: most-rearranged
    species first ("start on the species that are rearranged, then
    work outward"). Stops on the first pass that yields no
    improvement; any move that reduces the global score is kept.

    Seeds (chrom order and flips) are read from sp_to_chromorder and
    sp_to_chrom_flip and mutated in place; the function returns the
    refined dicts.
    """
    if coarse:
        line_data = _build_pair_anchor_data(
            species_order, rbh_df_list, sp_pair_to_rbh_df_list_index,
            sp_to_genesdfs, fet_weight=fet_weight)
    else:
        line_data = _build_pair_line_data(
            species_order, rbh_df_list, sp_pair_to_rbh_df_list_index,
            sp_to_genesdfs, fet_weight=fet_weight)

    # Build FET-graph cluster membership for the sticky LG-block move:
    # when we move a chromosome in one row we also move every other
    # chromosome that shares its synteny linkage group to the same
    # slot in its own row. This pulls a whole LG block sideways in
    # one step, letting the search escape the local minimum where a
    # single-row move would be uphill but a coordinated multi-row
    # move is downhill.
    node_to_cluster, _, _ = _build_fet_lg_clusters(species_order, rbh_df_list)
    cluster_members = {}
    for (sp, scaf), cid in node_to_cluster.items():
        cluster_members.setdefault(cid, []).append((sp, scaf))

    def pair_score(k):
        return _score_pair_lines(
            line_data[k],
            sp_to_chromorder[species_order[k]],
            sp_to_chromorder[species_order[k + 1]],
            sp_to_chrom_flip[species_order[k]],
            sp_to_chrom_flip[species_order[k + 1]])

    def neighboring_pair_indices(sp_idx):
        out = []
        if sp_idx > 0:
            out.append(sp_idx - 1)
        if sp_idx < len(species_order) - 1:
            out.append(sp_idx)
        return out

    pair_cache = {k: pair_score(k) for k in range(len(species_order) - 1)}

    def total_score():
        return sum(pair_cache.values())

    def rescore_neighbors(sp_idx):
        for k in neighboring_pair_indices(sp_idx):
            pair_cache[k] = pair_score(k)

    def chrom_crossings(sp_idx, chrom):
        """Crossings involving lines incident to `chrom` on the two
        adjacent pairs of species `sp_idx`. Used to sort chroms within
        a species so the worst-offenders are visited first.
        """
        total = 0.0
        for k in neighboring_pair_indices(sp_idx):
            recs = line_data[k]
            top_sp = species_order[k]
            bot_sp = species_order[k + 1]
            top_order = sp_to_chromorder[top_sp]
            bot_order = sp_to_chromorder[bot_sp]
            top_flip  = sp_to_chrom_flip[top_sp]
            bot_flip  = sp_to_chrom_flip[bot_sp]
            this_is_top = (species_order[sp_idx] == top_sp)
            for tchrom, bchrom, trk, brk, tsz, bsz, w in recs:
                if (this_is_top and tchrom != chrom) or \
                   (not this_is_top and bchrom != chrom):
                    continue
                if tchrom not in top_order or bchrom not in bot_order:
                    continue
                # quick proxy: contribution to crossings = weight * line count
                # (rather than exact; we just need ordering of severity)
                total += w
        return total

    cur = total_score()
    prev = cur
    if verbose:
        print("local search seed score = {:.0f}".format(cur))

    for pass_num in range(max_passes):
        # rank species by local crossings, worst first
        sp_local = [0.0] * len(species_order)
        for k, s in pair_cache.items():
            sp_local[k]     += s
            sp_local[k + 1] += s
        order_idx = sorted(range(len(species_order)), key=lambda i: -sp_local[i])

        improved_pass = False
        for sp_idx in order_idx:
            sp = species_order[sp_idx]
            chrom_seq = sorted(sp_to_chromorder[sp].items(), key=lambda x: x[1])
            chroms = [c for c, _ in chrom_seq]
            # Visit chroms in this row worst-first (chrom that has the
            # most weight in lines crossing into the adjacent pairs).
            # The user-cited case -- a single chromosome that crosses
            # many others on the row above -- is the kind we want to
            # try moving first.
            chroms.sort(key=lambda c: -chrom_crossings(sp_idx, c))
            # 1) try a flip on every chrom in this row
            for chrom in chroms:
                sp_to_chrom_flip[sp][chrom] = not sp_to_chrom_flip[sp][chrom]
                rescore_neighbors(sp_idx)
                new = total_score()
                if new < cur - 1e-9:
                    cur = new
                    improved_pass = True
                else:
                    sp_to_chrom_flip[sp][chrom] = not sp_to_chrom_flip[sp][chrom]
                    rescore_neighbors(sp_idx)
            # 2) Swap chrom_j with chrom_{j+offset} for a small set of
            #    offsets. Adjacent swap (offset=1) covers the local case;
            #    larger offsets let the search escape local minima where
            #    moving a chrom one slot only swaps it with a sibling and
            #    raises the score, but moving it two or three slots lands
            #    it next to its real partner.
            for offset in (1, 2, 3, 5):
                j = 0
                while j + offset < len(chroms):
                    a, b = chroms[j], chroms[j + offset]
                    sp_to_chromorder[sp][a], sp_to_chromorder[sp][b] = \
                        sp_to_chromorder[sp][b], sp_to_chromorder[sp][a]
                    rescore_neighbors(sp_idx)
                    new = total_score()
                    if new < cur - 1e-9:
                        cur = new
                        chroms[j], chroms[j + offset] = b, a
                        improved_pass = True
                    else:
                        sp_to_chromorder[sp][a], sp_to_chromorder[sp][b] = \
                            sp_to_chromorder[sp][b], sp_to_chromorder[sp][a]
                        rescore_neighbors(sp_idx)
                    j += 1
            # 3) Insertion move: for each chrom, try moving it to every
            #    other slot in the row. A swap can only exchange a
            #    chrom with a fixed neighbor; insertion lets the worst
            #    offender go directly to its true partner column,
            #    which is the kind of move adjacent-swap cannot make
            #    when the route between the two slots crosses unrelated
            #    chroms that each individual swap raises the score.
            for src_chrom in list(chroms):
                src_pos = sp_to_chromorder[sp][src_chrom]
                best_delta = 0.0
                best_dst = None
                # Try every destination slot
                cur_chroms = sorted(sp_to_chromorder[sp].items(),
                                    key=lambda x: x[1])
                ids_in_order = [c for c, _ in cur_chroms]
                for dst_pos in range(len(ids_in_order)):
                    if dst_pos == src_pos:
                        continue
                    new_order = list(ids_in_order)
                    new_order.pop(src_pos)
                    new_order.insert(dst_pos, src_chrom)
                    new_map = {c: k for k, c in enumerate(new_order)}
                    saved = dict(sp_to_chromorder[sp])
                    sp_to_chromorder[sp] = new_map
                    rescore_neighbors(sp_idx)
                    new = total_score()
                    if new < cur + best_delta - 1e-9:
                        best_delta = new - cur
                        best_dst = dst_pos
                    sp_to_chromorder[sp] = saved
                    rescore_neighbors(sp_idx)
                if best_dst is not None:
                    new_order = list(ids_in_order)
                    new_order.pop(src_pos)
                    new_order.insert(best_dst, src_chrom)
                    sp_to_chromorder[sp] = {c: k for k, c in enumerate(new_order)}
                    rescore_neighbors(sp_idx)
                    cur = cur + best_delta
                    improved_pass = True
                    # Update the local chroms list so the next iter sees
                    # the new order.
                    chroms = list(new_order)
            # 4) Sticky LG-block move: for each chrom whose FET cluster
            #    has members in other species, try moving the entire
            #    cluster as one block to a small set of candidate
            #    slots in each row simultaneously. Pulls a syntenic
            #    block sideways together rather than letting one row
            #    drift away from its partners.
            tried_clusters = set()
            for chrom in list(chroms):
                cid = node_to_cluster.get((sp, chrom))
                if cid is None or cid in tried_clusters:
                    continue
                tried_clusters.add(cid)
                members = cluster_members.get(cid, [])
                if len(members) < 2:
                    continue
                # Candidate destination slots: the current slot of any
                # member chrom in any row. These are the realistic
                # places the block could go.
                cand_slots = sorted({
                    sp_to_chromorder[m_sp].get(m_scaf)
                    for m_sp, m_scaf in members
                    if m_scaf in sp_to_chromorder.get(m_sp, {})
                })
                cand_slots = [s for s in cand_slots if s is not None]
                if not cand_slots:
                    continue
                # save current state for revert
                saved_orders = {sp_n: dict(sp_to_chromorder[sp_n])
                                for sp_n in species_order}
                # Try each candidate slot for the block
                for dst in cand_slots:
                    # restore from saved
                    for sp_n in species_order:
                        sp_to_chromorder[sp_n] = dict(saved_orders[sp_n])
                    # Move every member chrom to slot `dst` in its row
                    for m_sp, m_scaf in members:
                        if m_scaf not in sp_to_chromorder.get(m_sp, {}):
                            continue
                        row_chroms = sorted(sp_to_chromorder[m_sp].items(),
                                            key=lambda x: x[1])
                        ids = [c for c, _ in row_chroms]
                        if m_scaf not in ids:
                            continue
                        ids.remove(m_scaf)
                        dst_clamped = min(max(dst, 0), len(ids))
                        ids.insert(dst_clamped, m_scaf)
                        sp_to_chromorder[m_sp] = {c: k for k, c in enumerate(ids)}
                    # rescore all pairs (block move touches many)
                    for k in range(len(species_order) - 1):
                        pair_cache[k] = pair_score(k)
                    new = total_score()
                    if new < cur - 1e-9:
                        cur = new
                        improved_pass = True
                        # accept this candidate; refresh saved baseline
                        saved_orders = {sp_n: dict(sp_to_chromorder[sp_n])
                                        for sp_n in species_order}
                    else:
                        # revert
                        for sp_n in species_order:
                            sp_to_chromorder[sp_n] = dict(saved_orders[sp_n])
                        for k in range(len(species_order) - 1):
                            pair_cache[k] = pair_score(k)
                # refresh chroms list for next outer iter
                chroms = [c for c, _ in sorted(
                    sp_to_chromorder[sp].items(), key=lambda x: x[1])]
        if verbose:
            print("  pass {} done, score={:.0f}".format(pass_num, cur))
        if not improved_pass:
            break
        # Early stop: each pass costs the same; if the relative drop
        # falls below early_stop_frac we are deep in the long tail
        # and further passes won't visually change the plot.
        if prev > 0 and (prev - cur) / prev < early_stop_frac:
            if verbose:
                print("  early stop: delta {:.0f} < {} * {:.0f}".format(
                    prev - cur, early_stop_frac, prev))
            break
        prev = cur

    # Dedicated final flip-only sweep. After the order has settled, the
    # context in which each chromosome's orientation matters is fully
    # determined; some flips may now be improvable that weren't when
    # the search visited them earlier with a different ordering. User
    # called this out as key: at the end of optimization, re-check
    # every chromosome's LR/RL orientation.
    for _ in range(5):
        improved = False
        # process species worst-rearranged first, same as main loop
        sp_local = [0.0] * len(species_order)
        for k, s in pair_cache.items():
            sp_local[k] += s
            sp_local[k + 1] += s
        order_idx = sorted(range(len(species_order)), key=lambda i: -sp_local[i])
        for sp_idx in order_idx:
            sp = species_order[sp_idx]
            for chrom in list(sp_to_chrom_flip[sp].keys()):
                sp_to_chrom_flip[sp][chrom] = not sp_to_chrom_flip[sp][chrom]
                rescore_neighbors(sp_idx)
                new = total_score()
                if new < cur - 1e-9:
                    cur = new
                    improved = True
                else:
                    sp_to_chrom_flip[sp][chrom] = not sp_to_chrom_flip[sp][chrom]
                    rescore_neighbors(sp_idx)
        if not improved:
            break
    if verbose:
        print("after final flip sweep: score={:.0f}".format(cur))

    return sp_to_chromorder, sp_to_chrom_flip


def _iterative_order_propagation(species_order, rbh_df_list,
                                 sp_pair_to_rbh_df_list_index,
                                 sp_to_chromorder, sp_to_gene_order,
                                 sp_to_chr_to_size, sp_min_chr_size,
                                 sp_to_genesdfs, max_passes=10):
    """Bidirectional FET-best-partner sweep over the whole stack.
    The one-pass top-down cascade fixes only the relationship between
    each row and the row directly above it. A change in row k that
    needs to propagate up to row 0 (or vice versa) is invisible to it
    — exactly the case the user pointed out: CM/518.1 / CM/525.1 line
    up for the top three species but diverge by row 4, and fixing it
    means dragging row-4 scaffolds, flipping one, then propagating
    rows 3 -> 2 -> 1.

    Each pass alternates direction. Within a direction, every interior
    row is re-ordered with `_optimize_spA_based_on_rbh` using the
    currently-pinned neighbor on the other side as anchor; the top
    row is re-seeded with `_seed_top_order_from_partner` using row 1
    as anchor. We score after each pass with the same crossing-count
    metric the flip optimizer uses and keep the best-seen ordering,
    so an oscillation between two equally-good orderings doesn't
    leave us worse off.
    """
    def order_list(sp):
        return [s for s, _ in sorted(sp_to_chromorder[sp].items(),
                                     key=lambda x: x[1])]
    def total_crossings():
        total = 0
        for k in range(len(species_order) - 1):
            a, b = species_order[k], species_order[k + 1]
            lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([a, b]))]
            total += _crossings_between(a, b, rbh_df_list[lookup],
                                        sp_to_chromorder, sp_to_genesdfs)
        return total
    best_orders = {sp: dict(sp_to_chromorder[sp]) for sp in species_order}
    best_score  = total_crossings()
    no_improve  = 0
    for _ in range(max_passes):
        any_change = False
        for direction in (-1, +1):
            rng = (range(len(species_order) - 2, 0, -1) if direction == -1
                   else range(1, len(species_order)))
            for i in rng:
                sp = species_order[i]
                ref_idx = i + direction
                if not (0 <= ref_idx < len(species_order)):
                    continue
                ref = species_order[ref_idx]
                lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([sp, ref]))]
                new = _optimize_spA_based_on_rbh(
                    sp, ref, rbh_df_list[lookup],
                    sp_to_chromorder, sp_to_gene_order)
                new = _quality_check_chromosome_list(
                    sp, new, sp_to_chr_to_size,
                    sp_to_gene_order, sp_min_chr_size)
                if new != order_list(sp):
                    sp_to_chromorder[sp] = {s: k for k, s in enumerate(new)}
                    any_change = True
            # re-seed the top row from row 1 every sub-sweep so a row-1
            # reordering propagates to the very top.
            if len(species_order) > 1:
                top = species_order[0]
                new_top = _seed_top_order_from_partner(
                    top, species_order[1], rbh_df_list[0])
                new_top = _quality_check_chromosome_list(
                    top, new_top, sp_to_chr_to_size,
                    sp_to_gene_order, sp_min_chr_size)
                if new_top != order_list(top):
                    sp_to_chromorder[top] = {s: k for k, s in enumerate(new_top)}
                    any_change = True
        cur = total_crossings()
        if cur < best_score:
            best_score  = cur
            best_orders = {sp: dict(sp_to_chromorder[sp]) for sp in species_order}
            no_improve = 0
        else:
            no_improve += 1
        if not any_change or no_improve >= 2:
            break
    for sp in species_order:
        sp_to_chromorder[sp] = best_orders[sp]
    return sp_to_chromorder


def plot_bezier_lines(panel, topxL, bottomxL, colors, alpha, topy, bottomy):
    """
    Plot bezier curves between chromosome coordinates of different species.
    - Returns the panel, but with the lines now

    Parameters:
      - panel:    The panel on which we want to plot the lines
      - topxL:    X-axis value where the lines start on the top of the plot
      - bottomxL: X-axis value where the lines stop on the bottom of the plot
      - colors:   The colors of each of the lines
      - alpha:    The alpha values of each of the lines
      - topy:     The y-value of the top of the lines
      - bottomy:  The y-value of the bottom of the lines
    """
    indent = 1.0
    # plot the indices
    for i in range(len(topxL)):
        topx     = topxL[i]
        bottomx  = bottomxL[i]
        diff = abs(topx - bottomx)
        middlex = min([topx, bottomx]) + (diff/2)
        leftx  = middlex - ((diff/2)*indent)
        rightx = middlex + ((diff/2)*indent)
        second = (-1,-1)
        third  = (-1,-1)
        if topx <= bottomx:
            second = (leftx,  topy+0.5)
            third  = (rightx, topy+0.5)
        else:
            second = (rightx, topy+0.5)
            third  = (leftx,  topy+0.5)
        path_data = [
            (Path.MOVETO, (topx, topy)),
            (Path.CURVE4, second),
            (Path.CURVE4, third),
            (Path.CURVE4, (bottomx, bottomy)),
            ]
        codes, verts = zip(*path_data)
        path  = mpath.Path(verts, codes)
        if colors[i] == "#000000":
            zord = -50
        elif alpha[i] < 0.5:
            zord = -99
        else:
            zord = 1
        patch = patches.PathPatch(path, fill = False,
                                  facecolor=[0,0,0,0], lw = 0.25,
                                  alpha=alpha[i], edgecolor = colors[i],
                                  zorder = zord)
        panel.add_patch(patch)
    return panel

def _render_ribbon_figure(species_order, sp_to_chromorder, sp_to_chrom_flip,
                          sp_to_chr_to_size, sp_to_genesdfs,
                          sp_pair_to_rbh_df_list_index, rbh_df_list,
                          outfile, plot_all=True,
                          show_orientation_marks=True,
                          species_labels=None):
    """Render the ribbon plot to `outfile`. Extracted from
    ribbon_plot() so brushing can snapshot the layout after every
    sweep."""
    # now construct dataframes describing how to plot the chromosomes based on
    #  gene index or on chromosome coordinate
    sp_to_chrom_to_order = {}
    sp_to_chromdf        = {} # this is for plotting by chromosome or gene ix coordinates
    for i in range(len(species_order) - 1):
        thissp = species_order[i]
        nextsp = species_order[i+1]
        # make a sp_to_chromorder dict entry
        for sp in [thissp, nextsp]:
            ########################################
            # START THE DATAFRAME FOR THIS SPECIES
            ########################################
            chromdf = pd.DataFrame([k for k,v in sorted(sp_to_chromorder[sp].items(),
                                              key=lambda item: item[1])],
                                   columns=['chrom'])

            ########################################
            # MAGIC NUMBERS
            ########################################
            ## We make a special dataframe to figure out how to plot the chromosomes in percent of chroms
            #space_between_chrom = 0.0175
            #percent_chrom_as_spaces = space_between_chrom * (len(sp_to_chromorder[sp]) - 1)
            # This doesn't work well for large numbers of chromosomes, so use a specific percent.
            # - I am finding that it helps to use different percents depending on how many chromosomes there are.
            #   This is an awkward block of code, but it is readable and is easy to change.
            num_chroms = len(sp_to_chromorder[sp])
            min_chr_size = 10
            max_chr_size = 50
            min_gap_percent = 0.12
            max_gap_percent = 0.5
            if num_chroms < min_chr_size:
                percent_chrom_as_spaces = min_gap_percent
            elif num_chroms > max_chr_size:
                percent_chrom_as_spaces = max_gap_percent
            else:
                # scale the percent_chrom_as_spaces between min_gap_percent and max_gap_percent based on the number of chromosomes
                percent_chrom_as_spaces = min_gap_percent + (max_gap_percent - min_gap_percent) * (num_chroms - min_chr_size) / (max_chr_size - min_chr_size)
            print("sp: {}, num_chroms: {}, percent_chrom_as_spaces: {}".format(sp, num_chroms, percent_chrom_as_spaces))
            percent_chrom_as_chroms = 1 - percent_chrom_as_spaces
            space_between_chrom = percent_chrom_as_spaces / (len(sp_to_chromorder[sp]) - 1)

            ######################################################
            # THIS SECTION IS FOR ABSOLUTE CHROMOSOME COORDINATES
            ######################################################
            # We make a special dataframe to figure out how to plot the chromosomes in absolute basepair coordinates.
            #  We determine what percent of the horizontal line will be occupied by gaps.
            #  This number space_between_chrom is a value between 0-1.
            chromdf["chromix"] = chromdf["chrom"].map(sp_to_chromorder[sp])
            chromdf["chrsize"] = chromdf["chrom"].map(sp_to_chr_to_size[sp])
            total_chrom_len = sum(chromdf["chrsize"])
            chromdf["chrPlotPercent"] = (chromdf["chrsize"]/total_chrom_len
                                         ) * percent_chrom_as_chroms
            chromdf["chrPlotOffset"] = chromdf["chrPlotPercent"].cumsum(
                                       ).shift(1).fillna(0) + \
                                       (space_between_chrom * chromdf["chromix"])

            ######################################################
            # THIS SECTION IS FOR RELATIVE CHROMOSOME COORDINATES
            ######################################################
            # This section does the same thing as above, but it plots the
            #  chromosomes as relative coordinates based on number of genes.
            chromdf["ixSize"] = chromdf["chrom"].map(
                sp_to_genesdfs[sp].groupby("{}_scaf".format(sp)).size().to_dict()).fillna(0)
            totalIxSize = sum(chromdf["ixSize"])
            chromdf["ixPlotPercent"] = (chromdf["ixSize"]/totalIxSize) * percent_chrom_as_chroms
            chromdf["ixPlotOffset"] = chromdf["ixPlotPercent"].cumsum(
                                       ).shift(1).fillna(0) + \
                                       (space_between_chrom * chromdf["chromix"])
            #print(chromdf)
            sp_to_chromdf[sp] = chromdf

    # CALL THIS TO GET THE VISUAL STYLE WE NEED
    odp_plot.format_matplotlib()

    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
    # Preserve the vertical order of embedded images:
    matplotlib.rcParams['image.composite_image'] = False
    # text as font in pdf
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    # first we need to figure out the dimensions of the figure.
    #  Just make it the number of samples times the space, plus the buffer
    interspeciesHeight = 0.5
    panelHeight = interspeciesHeight * len(species_order)
    # Panel width leaves room for two-line italic species + grey accession
    # labels on the y-axis. Total figure width is 180 mm (7.087"); the
    # left margin holds the labels.
    panelWidth = 6.3

    #           two panels        top, bottom, middle
    bufferHeight = 1.5
    figHeight = (panelHeight*2) + (bufferHeight * 3)
    figWidth = 7.087  # 180 mm (single-column journal width)

    fig = plt.figure(figsize=(figWidth,figHeight))

    #set the panel dimensions
    leftStart = figWidth - panelWidth - 0.3  # 0.3" right margin
    bottomMargin = bufferHeight
    # pChr will host the chrom coordinate plots
    plt.gcf().text(leftStart/figWidth + (panelWidth/figWidth)/2,
                   (bottomMargin+panelHeight+0.8)/figHeight,
                   "p<=0.05 RBH results (Chr-coords)",
                   fontsize = 10, ha = "center", va = "bottom")
    pChr = plt.axes([leftStart/figWidth, #left
                   bottomMargin/figHeight,    #bottom
                   panelWidth/figWidth,   #width
                   panelHeight/figHeight])     #height
    pChr.tick_params(axis='both',which='both',
                   bottom=False, labelbottom=False,
                   left=False, labelleft=True,
                   right=False, labelright=False,
                   top=False, labeltop=False)
    # pIx title
    plt.gcf().text(leftStart/figWidth + (panelWidth/figWidth)/2,
                   ((bottomMargin*2)+(panelHeight*2)+0.8)/figHeight,
                   "p<=0.05 RBH results (RBH-gene-coords)",
                   fontsize = 10, ha = "center", va = "bottom")
    pIx = plt.axes([leftStart/figWidth, #left
                   ((bottomMargin*2)+panelHeight)/figHeight,    #bottom
                   panelWidth/figWidth,   #width
                   panelHeight/figHeight])     #height
    pIx.tick_params(axis='both',which='both',
                   bottom=False, labelbottom=False,
                   left=False, labelleft=True,
                   right=False, labelright=False,
                   top=False, labeltop=False)
    # make a panel for the legend, too
    panellg = plt.axes([ 0/figWidth, #left
                         0/figHeight,    #bottom
                         figWidth/figWidth,   #width
                         (bottomMargin*0.75)/figHeight])     #height
    panellg.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=False, labelright=False,
                        top=False, labeltop=False)

    color_labels = set()
    # get the gene group column where the filtercol row is less than or equal to 0.05
    for thisdf in rbh_df_list:
        filtercol = "break_FET"
        if filtercol not in thisdf.columns:
            if "alpha" in thisdf.columns:
                filtercol = "alpha"
            elif "pval" in thisdf.columns:
                filtercol = "pval"
            else:
                # there is no acceptable column, raise an error
                raise ValueError("No acceptable column for filtering in colordf")
        # now that we have picked the filter column, get the color labels
        color_labels = color_labels.union( \
                        set(thisdf.loc[thisdf[filtercol] <= 0.05, "gene_group"].tolist())) 
    color_labels = sorted(list(color_labels))

    # Get a dataframe of group to color based on frequency
    colordf = pd.concat(rbh_df_list)[["gene_group", "color"]
               ].groupby(["gene_group", "color"]
               ).size(
               ).reset_index(name='Freq'
               ).sort_values(by = ["Freq"], ascending = False
               ).drop_duplicates(subset = "gene_group"
               ).sort_values(by = ["gene_group"], ascending = True
               ).reset_index(drop=True)

    # get the colors of the labels in color_labels
    # zip color_df["gene_group"] and color_df["color"] to a dict
    color_dict = dict(zip(colordf["gene_group"], colordf["color"]))
    color_colors = [color_dict[x] for x in color_labels]


    # make a legend
    # set a legend if it is prot_to_color_mode
    legend_elements = []
    for ii in range(len(color_labels)):
        legend_elements.append(
            patches.Patch(facecolor=color_colors[ii],
            edgecolor='black', lw = 0,
            label=color_labels[ii])
           )
    panellg.legend(handles=legend_elements,
                   ncol = 10,
                   fontsize = 8, loc='center left')

    # plot all the lines
    for i in range(len(species_order)-1):
        thissp = species_order[i]
        nextsp = species_order[i+1]
        thisspChrom = sp_to_chromdf[thissp]
        nextspChrom = sp_to_chromdf[nextsp]
        thisspGenes = sp_to_genesdfs[thissp]
        nextspGenes = sp_to_genesdfs[nextsp]

        sp_pair_lookup = tuple(sorted([thissp, nextsp])) 
        lookup_index = sp_pair_to_rbh_df_list_index[sp_pair_lookup]

        tr = rbh_df_list[lookup_index].copy()
        # Choose the column used to decide whether a line is "significant".
        filtercol = "break_FET"
        if filtercol not in tr.columns:
            if "alpha" in tr.columns:
                filtercol = "alpha"
            else:
                print("ERROR: no break_FET or alpha column in rbh_df_list[{}]".format(i))
                sys.exit(1)
        if not plot_all:
            tr = tr.loc[tr[filtercol] <= 0.05, ]
            tr["alpha"] = 0.8
        else:
            tr["alpha"] = tr[filtercol].apply(lambda x: 0.8 if x <= 0.05 else 0.15)
        tr = tr[["{}_gene".format(thissp),
                 "{}_scaf".format(thissp),
                  "{}_pos".format(thissp),
                 "{}_gene".format(nextsp),
                 "{}_scaf".format(nextsp),
                  "{}_pos".format(nextsp),
                  "alpha",
                 "color"]]
        # handle the index plotting info
        #print(thisspChrom)
        #print()

        # Per-scaffold flip state from the rotation optimizer (issue #127).
        # When True we mirror the within-chrom position about the chromosome
        # midpoint, so the plotted x is (chrPercent - geneOffset) instead of
        # geneOffset. We apply the same flip to both rank ("Ix") and bp
        # ("Chr") coordinates so the two panels stay consistent.
        top_flip_map = sp_to_chrom_flip[thissp]
        bot_flip_map = sp_to_chrom_flip[nextsp]
        tr["_top_flipped"] = tr["{}_scaf".format(thissp)].map(top_flip_map).fillna(False)
        tr["_bot_flipped"] = tr["{}_scaf".format(nextsp)].map(bot_flip_map).fillna(False)

        # get the index plot position for the top
        chrom_to_ixSize    = dict(zip(thisspChrom["chrom"], thisspChrom["ixSize"]))
        chrom_to_ixPercent = dict(zip(thisspChrom["chrom"], thisspChrom["ixPlotPercent"]))
        chrom_to_ixOffset  = dict(zip(thisspChrom["chrom"], thisspChrom["ixPlotOffset"]))
        gene_to_ix = dict(zip(thisspGenes["{}_gene".format(thissp)],
                              thisspGenes["{}_chromIx".format(thissp)]))
        tr["topIx"] = tr["{}_gene".format(thissp)].map(gene_to_ix)
        tr["topIx_ChromSize"] = tr["{}_scaf".format(thissp)].map(chrom_to_ixSize)
        tr["topIx_ChromPercent"] = tr["{}_scaf".format(thissp)].map(chrom_to_ixPercent)
        _topIx_eff = tr["topIx"].where(~tr["_top_flipped"],
                                       (tr["topIx_ChromSize"] - 1) - tr["topIx"])
        tr["topIx_geneOffset"]   = (_topIx_eff / tr["topIx_ChromSize"]) * tr["topIx_ChromPercent"]
        tr["topIx_ChromOffset"] = tr["{}_scaf".format(thissp)].map(chrom_to_ixOffset)
        tr["topIx_finalOffset"] = tr["topIx_ChromOffset"] + tr["topIx_geneOffset"]
        delete = ["topIx", "topIx_ChromSize", "topIx_ChromPercent",
                       "topIx_geneOffset","topIx_ChromOffset"]
        tr = tr[[x for x in tr.columns if x not in delete]]
        tr = tr.sort_values(by=["topIx_finalOffset"]).reset_index(drop=True)

        # get the index plot position for the bottom
        chrom_to_ixSize = dict(zip(nextspChrom["chrom"], nextspChrom["ixSize"]))
        chrom_to_ixPercent = dict(zip(nextspChrom["chrom"], nextspChrom["ixPlotPercent"]))
        chrom_to_ixOffset = dict(zip(nextspChrom["chrom"], nextspChrom["ixPlotOffset"]))
        gene_to_ix = dict(zip(nextspGenes["{}_gene".format(nextsp)],
                              nextspGenes["{}_chromIx".format(nextsp)]))
        tr["bottomIx"] = tr["{}_gene".format(nextsp)].map(gene_to_ix)
        tr["bottomIx_ChromSize"] = tr["{}_scaf".format(nextsp)].map(chrom_to_ixSize)
        tr["bottomIx_ChromPercent"] = tr["{}_scaf".format(nextsp)].map(chrom_to_ixPercent)
        _botIx_eff = tr["bottomIx"].where(~tr["_bot_flipped"],
                                          (tr["bottomIx_ChromSize"] - 1) - tr["bottomIx"])
        tr["bottomIx_geneOffset"]   = (_botIx_eff / tr["bottomIx_ChromSize"]) * tr["bottomIx_ChromPercent"]
        tr["bottomIx_ChromOffset"] = tr["{}_scaf".format(nextsp)].map(chrom_to_ixOffset)
        tr["bottomIx_finalOffset"] = tr["bottomIx_ChromOffset"] + tr["bottomIx_geneOffset"]
        delete = ["bottomIx", "bottomIx_ChromSize", "bottomIx_ChromPercent",
                       "bottomIx_geneOffset","bottomIx_ChromOffset"]
        tr = tr[[x for x in tr.columns if x not in delete]]
        tr = tr.sort_values(by=["bottomIx_finalOffset"]).reset_index(drop=True)

        # get the chrom plot position for the top
        chrom_to_chrSize =    dict(zip(thisspChrom["chrom"], thisspChrom["chrsize"]))
        chrom_to_chrPercent = dict(zip(thisspChrom["chrom"], thisspChrom["chrPlotPercent"]))
        chrom_to_chrOffset =  dict(zip(thisspChrom["chrom"], thisspChrom["chrPlotOffset"]))
        tr["topChr_ChromSize"] = tr["{}_scaf".format(thissp)].map(chrom_to_chrSize)
        tr["topChr_ChromPercent"] = tr["{}_scaf".format(thissp)].map(chrom_to_chrPercent)
        _topChr_pos = tr["{}_pos".format(thissp)]
        _topChr_eff = _topChr_pos.where(~tr["_top_flipped"],
                                        tr["topChr_ChromSize"] - _topChr_pos)
        tr["topChr_geneOffset"]   = (_topChr_eff / tr["topChr_ChromSize"]) * tr["topChr_ChromPercent"]
        tr["topChr_ChromOffset"] = tr["{}_scaf".format(thissp)].map(chrom_to_chrOffset)
        tr["topChr_finalOffset"] = tr["topChr_ChromOffset"] + tr["topChr_geneOffset"]
        delete = ["topChr", "topChr_ChromSize", "topChr_ChromPercent",
                       "topChr_geneOffset","topChr_ChromOffset"]
        tr = tr[[x for x in tr.columns if x not in delete]]
        tr = tr.sort_values(by=["topChr_finalOffset"]).reset_index(drop=True)

        # get the chrom plot position for the bottom
        chrom_to_chrSize =    dict(zip(nextspChrom["chrom"], nextspChrom["chrsize"]))
        chrom_to_chrPercent = dict(zip(nextspChrom["chrom"], nextspChrom["chrPlotPercent"]))
        chrom_to_chrOffset =  dict(zip(nextspChrom["chrom"], nextspChrom["chrPlotOffset"]))
        tr["bottomChr_ChromSize"] = tr["{}_scaf".format(nextsp)].map(chrom_to_chrSize)
        tr["bottomChr_ChromPercent"] = tr["{}_scaf".format(nextsp)].map(chrom_to_chrPercent)
        _botChr_pos = tr["{}_pos".format(nextsp)]
        _botChr_eff = _botChr_pos.where(~tr["_bot_flipped"],
                                        tr["bottomChr_ChromSize"] - _botChr_pos)
        tr["bottomChr_geneOffset"]   = (_botChr_eff / tr["bottomChr_ChromSize"]) * tr["bottomChr_ChromPercent"]
        tr["bottomChr_ChromOffset"] = tr["{}_scaf".format(nextsp)].map(chrom_to_chrOffset)
        tr["bottomChr_finalOffset"] = tr["bottomChr_ChromOffset"] + tr["bottomChr_geneOffset"]
        delete = ["bottomChr", "bottomChr_ChromSize", "bottomChr_ChromPercent",
                       "bottomChr_geneOffset","bottomChr_ChromOffset"]
        tr = tr[[x for x in tr.columns if x not in delete]]
        tr = tr.sort_values(by=["bottomChr_finalOffset"]).reset_index(drop=True)

        # plot the indices
        pIx = plot_bezier_lines(pIx,
                         tr["topIx_finalOffset"],
                         tr["bottomIx_finalOffset"],
                         tr["color"],
                         tr["alpha"],
                         i, i+1)

        pChr = plot_bezier_lines(pChr,
                         tr["topChr_finalOffset"],
                         tr["bottomChr_finalOffset"],
                         tr["color"],
                         tr["alpha"],
                         i, i+1)

    # Now we plot all the chroms. Each chromosome is drawn as a chevron-
    # tipped bar (when show_orientation_marks is True): the bar IS the
    # arrow rather than a bar + separate triangle, so dense plots stay
    # readable. Labels are rotated 90 degrees and trimmed to drop the
    # long common prefix shared by an assembly's scaffolds (e.g. all
    # "NC_05" prefix), so the distinguishing tail is left.
    for i in range(len(species_order)):
        sp = species_order[i]
        flip_map = sp_to_chrom_flip.get(sp, {})
        all_chroms = list(sp_to_chromdf[sp]["chrom"])
        short_label = _shorten_chrom_labels(all_chroms)
        for index, row in sp_to_chromdf[sp].iterrows():
            flipped = bool(flip_map.get(row["chrom"], False))
            label   = short_label.get(row["chrom"], row["chrom"])
            # chromosome-coordinate panel
            x1 = row["chrPlotOffset"]
            x2 = row["chrPlotOffset"] + row["chrPlotPercent"]
            if show_orientation_marks:
                _draw_oriented_bar(pChr, x1, x2, i, flipped)
            else:
                pChr.plot([x1, x2], [i, i], "k-")
            # Label anchored just left of and above the chromosome
            # bar's left end, rotated 90 degrees CCW.
            pChr.text(x1 - 0.005, i - 0.03, label, ha="left", va="top",
                      fontsize=5, rotation=90, rotation_mode="anchor")

            # rbh-gene-coordinate panel
            x1 = row["ixPlotOffset"]
            x2 = row["ixPlotOffset"] + row["ixPlotPercent"]
            if show_orientation_marks:
                _draw_oriented_bar(pIx, x1, x2, i, flipped)
            else:
                pIx.plot([x1, x2], [i, i], "k-")
            pIx.text(x1 - 0.005, i - 0.03, label, ha="left", va="top",
                     fontsize=4, rotation=90, rotation_mode="anchor")

    # flip the y axes
    for p in [pChr, pIx]:
        # remove some spines
        for side in ["top", "right", "bottom", "left"]:
            p.spines[side].set_visible(False)
        # set the limits
        p.set_xlim([-0.02,1.02])
        # Extra slack at top and bottom for rotated scaffold labels.
        p.set_ylim([-0.03, len(species_order) - 1 + 0.55])
        # flip the axes
        p.set_ylim(p.get_ylim()[::-1])
        # set the axis labels
        p.set_yticks(list(range(len(species_order))))
        # Custom two-line y-axis labels: species name (italic) on top,
        # assembly accession (grey) below. Matplotlib doesn't support
        # mixed styling inside a single tick label, so we suppress the
        # default tick labels and draw each row's label with two
        # Text() artists at the same y. Robust to missing inputs --
        # species_labels may be None or have only some entries, and
        # the species id may or may not contain a GC[AF] accession.
        import re as _re
        p.set_yticklabels([""] * len(species_order))
        for i, sp in enumerate(species_order):
            if species_labels and sp in species_labels:
                primary = species_labels[sp]
            else:
                # primary fallback: strip -taxid-GC...xxx if present,
                # otherwise just use the species id as-is.
                primary = _re.sub(r"-\d+-GC[AF]\d+\.\d+$", "", sp)
                if primary == sp:
                    primary = sp  # unchanged
            acc_m = _re.search(r"(GC[AF]\d+\.\d+)$", sp)
            acc = acc_m.group(1) if acc_m else ""
            # x at the left edge of the panel data area
            x_anchor = p.get_xlim()[0] - 0.005
            p.text(x_anchor, i - 0.06, primary, ha="right", va="bottom",
                   fontsize=7, fontstyle="italic")
            if acc:
                p.text(x_anchor, i + 0.06, acc, ha="right", va="top",
                       fontsize=6, color="#888888")

    plt.savefig(outfile)
    plt.close('all')


def ribbon_plot(species_order, rbh_filelist,
                sp_to_chr_to_size,
                sp_min_chr_size, outfile,
                sp_to_gene_order = {},
                chr_sort_order   = "custom",
                plot_all = False,
                optimize_chrom_rotation = True,
                show_orientation_marks  = True,
                species_labels = None):
    """
    Takes in a list of species as the plotting order,
     a list of rbh files, and a dict of species_to_chr_to_sizes

    In the future for the rbh file parse the headers.

    There are several ways that the chromosomes can be sorted.
      chr_sort_order < custom | optimal-top | optimal-size | optimal-random >
        custom         - use the custom sorting order for EVERY species in chromorder
        optimal-top    - use the custom order for the topmost species, then optimizes everything else
        optimal-size   - sort the top species' chromosomes by number of genes, optimize everything else
        optimal-chr-or - use `chromorder` when possible, optimize everything else
        optimal-random - randomly sort the chromosomes of the top species, optimize everything else

    optimize_chrom_rotation:
        If True (default), every scaffold's plot direction is chosen
        independently of the fasta orientation, cascading top-down, to
        minimize ribbon-line crossings. See _optimize_chrom_flips_top_down.
        Resolves the "twist" artefact described in odp issue #127.
    show_orientation_marks:
        If True (default), each chromosome bar gets a small triangle at
        the fasta-3' end. The triangle points right (▶) for scaffolds
        plotted forward and left (◀) for scaffolds plotted reversed,
        so the viewer can tell the plot direction from the fasta
        direction at a glance.
    """

    # check sp_to_gene_order. Reset to empty dict if None
    if type(sp_to_gene_order) == type(None):
        sp_to_gene_order = {}

    import random

    # make a list of the dataframes to open for these analyses
    rbh_df_list = [pd.read_csv(x, sep = "\t",
                   header = "infer",
                   index_col = None) for x in rbh_filelist]
    
    # This list tells the program later which index of RBH files to look at.
    # Useful for cases where we submit just one .rbh file with all of the entries
    sp_pair_to_rbh_df_list_index = {}
    for i in range(len(rbh_df_list)):
        thisdf = rbh_df_list[i]
        sp_list = [x.split("_")[0] for x in thisdf.columns if "_scaf" in x]
        for j in range(len(sp_list)-1):
            for k in range(j+1, len(sp_list)):
                sp_combo = tuple(sorted([sp_list[j], sp_list[k]]))
                if sp_combo not in sp_pair_to_rbh_df_list_index:
                    sp_pair_to_rbh_df_list_index[sp_combo] = i

    sp_to_genesdfs = {}
    # make composite gene index dataframe for each species
    #  we will use this later to plot by gene index
    #  rather than the chromosome index

    for thisrbh in rbh_df_list:
        thesesp = [x.split("_")[0] for x in thisrbh.columns if "_scaf" in x]
        for thissp in thesesp:
            if thissp not in sp_to_genesdfs:
                sp_to_genesdfs[thissp] = []
                sp_to_genesdfs[thissp].append(thisrbh[[x for x in thisrbh.columns if thissp in x]])

    # we now concat and deduplicate the gene index dfs
    for thissp in sp_to_genesdfs:
        sp_to_genesdfs[thissp] = pd.concat(sp_to_genesdfs[thissp]
                                  ).sort_values(by=["{}_gene".format(thissp)]
                                  ).drop_duplicates(subset=["{}_gene".format(thissp)]
                                  ).sort_values(by=["{}_scaf".format(thissp),
                                                    "{}_pos".format(thissp)],
                                                ascending=[True, True]
                                  ).reset_index(drop=True)
        sp_to_genesdfs[thissp] = sp_to_genesdfs[thissp][["{}_gene".format(thissp),
                                                         "{}_scaf".format(thissp),
                                                         "{}_pos".format(thissp)]]
        # Assign within-scaffold index (chromIx). Replaces a pandas-version-
        # fragile groupby().apply().reset_index() pattern with a direct
        # cumcount, which is the same operation expressed plainly.
        sp_to_genesdfs[thissp]["{}_chromIx".format(thissp)] = sp_to_genesdfs[thissp].groupby(
            "{}_scaf".format(thissp)).cumcount()
        #print(sp_to_genesdfs[thissp])

    # FET-graph LG mode short-circuits the row-by-row cascade. It
    # builds one global synteny linkage graph from every adjacent-row
    # RBH, partitions chromosomes into connected components, and
    # orders each row by component id. Members of the same linkage
    # group share an x-region across every species so ribbons run
    # straight down. Handles in-row fusions/fissions correctly (two
    # chroms of one species linked to the same chrom in another
    # species end up in the same component even when their dominant
    # ALG labels differ).
    if chr_sort_order == "optimal-lg":
        sp_to_chromorder = _order_by_lg(
            species_order, rbh_df_list, sp_to_chr_to_size,
            sp_min_chr_size, sp_to_gene_order,
            sp_pair_to_rbh_df_list_index)
        # iterate / cascade blocks below are bypassed by this early
        # exit from the order-decision stage; pad sp_to_chromorder
        # below for completeness.
    # This block is used for determining the chromosome order for the figure
    sp_to_chromorder     = sp_to_chromorder if chr_sort_order == "optimal-lg" else {}
    print("species_order: {}".format(species_order))
    if chr_sort_order == "optimal-lg":
        # skip the row-by-row cascade entirely
        species_loop_range = []
    else:
        species_loop_range = range(0, len(species_order))
    for i in species_loop_range:
        sp = species_order[i]

        templist = [] # templist is used to hold the chromosome order
        if i == 0:
            # This is the first species, there is a special case for it
            #      chr_sort_order < custom | optimal-top | optimal-size | optimal-random >
            if chr_sort_order in ["custom", "optimal-top", "optimal-chr-or"]:
                # in this case we use the custom order or take the order from the file
                if not sp_to_gene_order or (sp not in sp_to_gene_order):
                    templist = sorted(list(rbh_df_list[i]["{}_scaf".format(sp)].unique()))
                else:
                    templist = sp_to_gene_order[sp]
            elif chr_sort_order in ["optimal-size", "optimal-random"] \
                  or chr_sort_order in _DETANGLE_MODES:
                # Top species seed. Plain value_counts() ignores partner
                # identity, so two top-row chromosomes that both map to
                # the same row-2 chromosome could end up scattered, and
                # the cascade below can't fix that. Seed the top row
                # from FET-significant best-partner clusters in row 2
                # instead (odp issue #127).
                if len(species_order) > 1:
                    templist = _seed_top_order_from_partner(
                        sp, species_order[1], rbh_df_list[i])
                else:
                    templist = list(rbh_df_list[i]["{}_scaf".format(sp)].value_counts().index)
                if sp in sp_to_gene_order:
                    # append the chromosomes that are in the config file but not in the rbh file
                    for x in sp_to_gene_order[sp]:
                        if x not in templist:
                            templist.append(x)
                if chr_sort_order == "optimal-random":
                    random.shuffle(templist)
        else:
            # we construct this variable to lookup the species in sp_pair_to_rbh_df_list_index
            sp_pair_lookup = tuple(sorted([species_order[i-1], sp])) 
            lookup_index = sp_pair_to_rbh_df_list_index[sp_pair_lookup]

            # this is the second or later species, change sort order depending on what was there
            if chr_sort_order == "custom":
                # in this case we use the custom order or take the order from the file
                if not sp_to_gene_order or (sp not in sp_to_gene_order):
                    # This is the case where we don't have a custom order for this species
                    #templist = sorted(list(rbh_df_list[i]["{}_scaf".format(sp)].unique()))
                    templist = sorted(list(rbh_df_list[lookup_index]["{}_scaf".format(sp)].unique()))
                else:
                    templist = sp_to_gene_order[sp]
            elif (chr_sort_order == "optimal-chr-or") and (sp in sp_to_gene_order):
                    templist = sp_to_gene_order[sp]
            elif chr_sort_order in _DETANGLE_MODES:
                # detangle by anchor position on the species directly above
                prevsp = species_order[i - 1]
                method = ("median" if chr_sort_order == "optimal-median"
                          else "barycenter")
                kept = set(rbh_df_list[lookup_index]["{}_scaf".format(sp)].unique())
                if sp in sp_to_gene_order:
                    kept |= set(sp_to_gene_order[sp])
                templist = _detangle_one_side(sp, prevsp, rbh_df_list[lookup_index],
                                              sp_to_chromorder, sp_to_genesdfs,
                                              method, kept)
            else:
                # we optimize every other case
                prevsp = species_order[i-1]
                templist = _optimize_spA_based_on_rbh(sp, prevsp, rbh_df_list[lookup_index],
                                                      sp_to_chromorder, sp_to_gene_order)
        # Now we perform some checks on the list of chromosomes and add them to the dict
        # check here to make sure that there aren't duplicate entries in the chromosomes
        templist = _quality_check_chromosome_list(sp, templist, sp_to_chr_to_size,
                                                  sp_to_gene_order, sp_min_chr_size)
        sp_to_chromorder[sp] = {templist[i]: i for i in range(len(templist))}
        #print(sp, sp_to_chromorder)

    # Iterated barycenter / median sweeps. The top-down pass above already
    # set an initial ordering. Now we sweep bottom-up and top-down again
    # using the *neighbor* row's positions, alternating until no change
    # or max_passes. Each species's first scaffold is held to the row's
    # current top so we don't drift the whole plot sideways.
    if chr_sort_order == "optimal-barycenter-iter":
        for _ in range(8):
            changed = False
            for direction in (-1, +1):
                rng = range(len(species_order) - 2, 0, -1) if direction == -1 \
                      else range(1, len(species_order) - 1)
                for i in rng:
                    sp = species_order[i]
                    ref = species_order[i + direction]
                    lookup = sp_pair_to_rbh_df_list_index[tuple(sorted([sp, ref]))]
                    new = _detangle_one_side(sp, ref, rbh_df_list[lookup],
                                             sp_to_chromorder, sp_to_genesdfs,
                                             "barycenter",
                                             set(sp_to_chromorder[sp].keys()))
                    cur = sorted(sp_to_chromorder[sp].keys(),
                                 key=lambda s: sp_to_chromorder[sp][s])
                    if new != cur:
                        changed = True
                        sp_to_chromorder[sp] = {s: k for k, s in enumerate(new)}
            if not changed:
                break

    # Greedy adjacent-swap polish. Runs on top of whichever method built
    # the ordering. Keeps swaps that reduce the sum of crossings against
    # the row above and the row below.
    if chr_sort_order == "optimal-swap":
        sp_to_chromorder = _greedy_swap_polish(species_order, rbh_df_list,
                                               sp_pair_to_rbh_df_list_index,
                                               sp_to_chromorder, sp_to_genesdfs)

    # There is an edge case where, if we used option "optimal-chr-or", but the 0th through nth species weren't
    #  in the chromosome order, we need to optimize the species in reverse order
    if chr_sort_order == "optimal-chr-or":
        # If someone uses this option when there is no chromosome order
        #  we just skip this step
        if len(sp_to_gene_order) > 0:
            first_optimized = 99999999
            for i in range(len(species_order)):
                sp = species_order[i]
                if sp in sp_to_gene_order:
                    first_optimized = i
                    break
            while first_optimized != 0:
                sp = species_order[first_optimized - 1]
                prevsp = species_order[first_optimized]
                # optimize the chromosome order
                templist = _optimize_spA_based_on_rbh(sp, prevsp, rbh_df_list[first_optimized - 1], sp_to_chromorder, sp_to_gene_order)
                templist = _quality_check_chromosome_list(sp, templist, sp_to_chr_to_size, sp_to_gene_order, sp_min_chr_size)
                sp_to_chromorder[sp] = {templist[i]: i for i in range(len(templist))}
                first_optimized -= 1

    # Iterative bidirectional FET-best-partner propagation. The plain
    # top-down cascade above only relates each row to the row directly
    # above it; chains of dependencies further down can't propagate up.
    # We sweep top-down and bottom-up to fixed point (or best-seen, in
    # case of oscillation between two equally good orderings). Runs
    # only when an "optimal-*" sort mode is selected — the "custom"
    # path is left untouched so user-supplied orderings stick.
    if chr_sort_order not in ("custom", "optimal-lg") and len(species_order) > 1:
        sp_to_chromorder = _iterative_order_propagation(
            species_order, rbh_df_list, sp_pair_to_rbh_df_list_index,
            sp_to_chromorder, sp_to_gene_order, sp_to_chr_to_size,
            sp_min_chr_size, sp_to_genesdfs)

    # Decide per-scaffold plot direction (issue #127). The chromosome
    # ordering is now fixed; we only choose, for each scaffold, whether to
    # display it in the fasta-forward direction or reversed. Cascades
    # top-down so a flip on row i is judged against the already-pinned
    # row i-1.
    if optimize_chrom_rotation:
        sp_to_chrom_flip = _optimize_chrom_flips_top_down(
            species_order, rbh_df_list, sp_pair_to_rbh_df_list_index,
            sp_to_chromorder, sp_to_genesdfs)
    else:
        sp_to_chrom_flip = {sp: {scaf: False for scaf in sp_to_chromorder[sp]}
                            for sp in species_order}

    # Brushing sweep (FET-weighted barycenter) + per-direction flip
    # cascade. Snapshot a PDF after every brush stroke so the user can
    # watch the layout evolve. Local search disabled for this run --
    # brushing is the sole optimizer so its effect is visible without
    # being masked by subsequent refinement passes.
    if optimize_chrom_rotation and len(species_order) > 1:
        import os as _os
        import glob as _glob
        out_dir = _os.path.dirname(_os.path.abspath(outfile))
        base = _os.path.splitext(_os.path.basename(outfile))[0]
        # Pre-run cleanup: delete snapshot PDFs from any prior run so
        # the user only sees this run's progression. The final
        # PDF is overwritten by the last _render_ribbon_figure call.
        for _stale in _glob.glob(_os.path.join(
                out_dir, base + "_brush_gen*.pdf")):
            try:
                _os.remove(_stale)
            except OSError:
                pass

        def _snapshot(gen, direction, order_state, flip_state, crossings):
            # Sweep ordinal within gen for alphabetical sort = run
            # order: seed=00, td=01, tm=02, bu=03, bm=04.
            sub = {"seed": "00", "td": "01", "tm": "02",
                   "bu": "03", "bm": "04"}.get(direction, direction)
            snap_path = _os.path.join(
                out_dir, "{}_brush_gen{:02d}_{}_{}.pdf".format(
                    base, gen, sub, direction))
            _render_ribbon_figure(
                species_order, order_state, flip_state,
                sp_to_chr_to_size, sp_to_genesdfs,
                sp_pair_to_rbh_df_list_index, rbh_df_list,
                snap_path, plot_all=plot_all,
                show_orientation_marks=show_orientation_marks,
                species_labels=species_labels)
            print("wrote snapshot {}  crossings: {:.0f}".format(
                snap_path, crossings), flush=True)

        # Iterated local search with random restarts. Each restart
        # perturbs the best-seen state and re-runs the full brushing
        # sweep; the brushing's own best-seen tracking handles each
        # local descent. Across restarts we keep the global best.
        # The perturbation is a random shuffle of a few chromosomes in
        # randomly chosen rows -- enough to break out of the current
        # local minimum, small enough that the next brushing can
        # recover quickly.
        import random as _random
        import copy as _copy
        n_restarts = 8
        perturb_rows = 3
        perturb_chroms_per_row = 5
        _random.seed(42)

        def _save_state():
            return ({sp: dict(sp_to_chromorder[sp]) for sp in species_order},
                    {sp: dict(sp_to_chrom_flip[sp]) for sp in species_order})

        def _restore_state(saved):
            ords, flps = saved
            for sp in species_order:
                sp_to_chromorder[sp] = dict(ords[sp])
                sp_to_chrom_flip[sp] = dict(flps[sp])

        def _perturb():
            """Randomly shuffle a window of chromosomes within a few
            randomly chosen rows. Each chosen row picks a window of
            perturb_chroms_per_row consecutive slots and shuffles."""
            rows = _random.sample(species_order, k=min(perturb_rows, len(species_order)))
            for sp in rows:
                chroms_by_pos = sorted(sp_to_chromorder[sp].items(),
                                       key=lambda x: x[1])
                ids = [c for c, _ in chroms_by_pos]
                n = len(ids)
                k = min(perturb_chroms_per_row, n)
                if k < 2:
                    continue
                start = _random.randint(0, n - k)
                window = ids[start:start + k]
                _random.shuffle(window)
                ids = ids[:start] + window + ids[start + k:]
                sp_to_chromorder[sp] = {c: i for i, c in enumerate(ids)}

        # initial descent
        sp_to_chromorder, sp_to_chrom_flip = _brushing_sweep(
            species_order, rbh_df_list, sp_pair_to_rbh_df_list_index,
            sp_to_chromorder, sp_to_chrom_flip, sp_to_genesdfs,
            max_iters=50, verbose=True, snapshot_fn=_snapshot)
        # score after first descent
        _ld0 = _build_pair_anchor_data(
            species_order, rbh_df_list, sp_pair_to_rbh_df_list_index,
            sp_to_genesdfs, fet_weight=1000.0)
        def _total_now():
            return sum(_score_pair_lines(
                _ld0[k],
                sp_to_chromorder[species_order[k]],
                sp_to_chromorder[species_order[k + 1]],
                sp_to_chrom_flip[species_order[k]],
                sp_to_chrom_flip[species_order[k + 1]])
                for k in range(len(species_order) - 1))
        global_best = _total_now()
        global_best_state = _save_state()
        print("[restart] initial descent best: {:.0f}".format(global_best),
              flush=True)

        for r in range(1, n_restarts + 1):
            _restore_state(global_best_state)
            _perturb()
            print("[restart {}/{}] perturbed; brushing again..."
                  .format(r, n_restarts), flush=True)

            def _restart_snapshot(gen, direction, order_state, flip_state, crossings):
                sub = {"seed": "00", "td": "01", "tm": "02",
                       "bu": "03", "bm": "04"}.get(direction, direction)
                sp = _os.path.join(
                    out_dir, "{}_brush_rest{:02d}_gen{:02d}_{}_{}.pdf".format(
                        base, r, gen, sub, direction))
                _render_ribbon_figure(
                    species_order, order_state, flip_state,
                    sp_to_chr_to_size, sp_to_genesdfs,
                    sp_pair_to_rbh_df_list_index, rbh_df_list,
                    sp, plot_all=plot_all,
                    show_orientation_marks=show_orientation_marks,
                    species_labels=species_labels)
                print("wrote snapshot {}  crossings: {:.0f}".format(
                    sp, crossings), flush=True)

            sp_to_chromorder, sp_to_chrom_flip = _brushing_sweep(
                species_order, rbh_df_list, sp_pair_to_rbh_df_list_index,
                sp_to_chromorder, sp_to_chrom_flip, sp_to_genesdfs,
                max_iters=50, verbose=True, snapshot_fn=_restart_snapshot)
            score = _total_now()
            if score < global_best:
                print("[restart {}/{}] NEW global best {:.0f} (was {:.0f})".format(
                    r, n_restarts, score, global_best), flush=True)
                global_best = score
                global_best_state = _save_state()
            else:
                print("[restart {}/{}] no improvement ({:.0f} vs best {:.0f})"
                      .format(r, n_restarts, score, global_best), flush=True)

        _restore_state(global_best_state)
        print("[restart] global best across {} restarts: {:.0f}".format(
            n_restarts, global_best), flush=True)

    _render_ribbon_figure(species_order, sp_to_chromorder, sp_to_chrom_flip,
                          sp_to_chr_to_size, sp_to_genesdfs,
                          sp_pair_to_rbh_df_list_index, rbh_df_list,
                          outfile, plot_all=plot_all,
                          show_orientation_marks=show_orientation_marks,
                          species_labels=species_labels)
    return sp_to_chromorder
