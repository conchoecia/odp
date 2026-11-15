"""
Optional ribbon-plot features for odp_rbh_to_ribbon.

These come from the Schultz Lab ribbon standard (2026). Every option is off by
default, so a config without them draws exactly what odp drew before.

Functions here take plain data (positions, colors, panel coordinates) so they can
be tested without Snakemake. odp_rbh_to_ribbon calls them from subway_plot().
"""
import bisect
import collections
import math
import random

import matplotlib.colors as mcolors
from matplotlib.collections import PolyCollection
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path

UNASSIGNED = "#000000"   # odp's color for orthologs in no linkage group


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def is_unassigned(color, group):
    return str(color).upper() == UNASSIGNED or str(group) in ("None", "nan", "")


def shade(color, t):
    """Light (t = 0) through the color (t = 0.5) to dark (t = 1)."""
    c = mcolors.to_rgb(color)
    if t < 0.5:
        k = (0.5 - t) * 1.3   # up to 65% white
        return mcolors.to_hex(tuple(ci + (1 - ci) * k for ci in c))
    k = (t - 0.5) * 1.1       # up to 55% black
    return mcolors.to_hex(tuple(ci * (1 - k) for ci in c))


def _apca_y(rgb):
    y = 0.2126729 * rgb[0] ** 2.4 + 0.7151522 * rgb[1] ** 2.4 + 0.0721750 * rgb[2] ** 2.4
    return y + (0.022 - y) ** 1.414 if y < 0.022 else y


def apca_lc(text_rgb, bg_rgb):
    """APCA lightness contrast (the method drafted for WCAG 3)."""
    yt, yb = _apca_y(text_rgb), _apca_y(bg_rgb)
    if yb > yt:
        sapc = (yb ** 0.56 - yt ** 0.57) * 1.14
        return 0.0 if sapc < 0.1 else (sapc - 0.027) * 100
    sapc = (yb ** 0.65 - yt ** 0.62) * 1.14
    return 0.0 if sapc > -0.1 else (sapc + 0.027) * 100


def text_color_for(fill, overrides=None, name=None):
    """White or black text on a colored label box, by APCA contrast."""
    if overrides and name in overrides:
        return overrides[name]
    rgb = mcolors.to_rgb(fill)
    return "white" if abs(apca_lc((1, 1, 1), rgb)) > abs(apca_lc((0, 0, 0), rgb)) else "black"


def rounded_box(x, y, w, h, rx, ry, **kw):
    """A box with rounded corners in data units; rx and ry are the corner radii along x and y."""
    rx = max(0.0, min(rx, w / 2))
    ry = max(0.0, min(ry, h / 2))
    return FancyBboxPatch((x + rx, y + ry), max(w - 2 * rx, 1e-9), max(h - 2 * ry, 1e-9),
                          boxstyle="round,pad=%f,rounding_size=%f" % (rx, rx),
                          mutation_aspect=(ry / rx) if rx > 0 else 1, **kw)


# ---------------------------------------------------------------------------
# small scaffolds: a size relative to the other chromosomes
# ---------------------------------------------------------------------------
def relative_min_size(sizes, fraction):
    """Absolute minimum size for scaffolds that must reach `fraction` of the median of `sizes`."""
    if not sizes or not fraction:
        return 0
    s = sorted(sizes)
    n = len(s)
    median = s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2
    return fraction * median


# ---------------------------------------------------------------------------
# linkage-group label boxes above the top row
# ---------------------------------------------------------------------------
def label_homes(points):
    """points: [(group, chrom, x)] for the top species. Each group's home is the chromosome holding most of its
    orthologs; groups sharing a home spread along it in the order of their median positions.
    Returns [(group, home chrom, preferred x)] with the home's span looked up by the caller."""
    per = collections.defaultdict(lambda: collections.defaultdict(list))
    for g, c, x in points:
        per[g][c].append(x)
    shared = collections.defaultdict(list)
    for g, by in per.items():
        home = max(by, key=lambda c: len(by[c]))
        xs = sorted(by[home])
        shared[home].append((xs[len(xs) // 2], g))
    out = []
    for home, gs in shared.items():
        gs.sort()
        for i, (_, g) in enumerate(gs):
            out.append((g, home, i, len(gs)))
    return out


def place_labels(homes, spans, widths, lo, hi, gap, max_tiers=8, step=None):
    """Place label boxes over their home chromosomes.
    homes: [(group, chrom, i, n)] from label_homes; spans: {chrom: (x0, x1)}; widths: {group: box width}.
    Each box moves only along its own chromosome; boxes stack in as many tiers as needed.
    Returns [(group, center x, width, tier)] and the number of tiers."""
    step = step or (hi - lo) / 2000.0
    items = []
    for g, c, i, n in homes:
        x0, x1 = spans[c]
        items.append((g, x0, x1, x0 + (x1 - x0) * (i + 0.5) / n))
    items.sort(key=lambda t: t[3])
    tiers, placed = [], []

    def free(t, a, b):
        return all(b + gap <= l or a >= r + gap for l, r in tiers[t])
    for g, x0, x1, cx in items:
        bw = widths[g]
        a, b = (x0 + bw / 2, x1 - bw / 2) if bw <= x1 - x0 else ((x0 + x1) / 2, (x0 + x1) / 2)
        a, b = max(a, lo + bw / 2), min(b, hi - bw / 2)
        if a > b:
            a = b = min(max((x0 + x1) / 2, lo + bw / 2), hi - bw / 2)
        pref = min(max(cx, a), b)
        cands = [pref] + [x for k in range(1, 4001) for x in (pref + step * k, pref - step * k) if a <= x <= b]
        spot = None
        for k in range(max_tiers):
            if k == len(tiers):
                tiers.append([])
            spot = next(((k, x) for x in cands if free(k, x - bw / 2, x + bw / 2)), None)
            if spot:
                break
        if spot is None:
            spot = (len(tiers) - 1, pref)
        tiers[spot[0]].append((spot[1] - bw / 2, spot[1] + bw / 2))
        placed.append((g, spot[1], bw, spot[0]))
    return placed, max(1, len(tiers))


# ---------------------------------------------------------------------------
# gene-order shading (order_ref)
# ---------------------------------------------------------------------------
def order_shades(pair_tables, species_order, ref):
    """Rank (0-1) of each gene along the reference species' chromosome among the genes of its linkage group,
    carried to the other species through the ortholog pairs. pair_tables[i] is the table between
    species_order[i] and species_order[i + 1] with columns SP_gene, SP_scaf, SP_pos, gene_group, color.
    Returns {(species, gene): t}."""
    t = {}
    ri = species_order.index(ref)
    pos = collections.defaultdict(set)
    for k in (ri - 1, ri):
        if 0 <= k < len(pair_tables):
            df = pair_tables[k]
            for g, sc, p, grp, col in zip(df[ref + "_gene"], df[ref + "_scaf"], df[ref + "_pos"], df["gene_group"], df["color"]):
                if not is_unassigned(col, grp):
                    pos[(str(grp), sc)].add(float(p))
    ranks = {k: sorted(v) for k, v in pos.items()}
    for k in (ri - 1, ri):
        if 0 <= k < len(pair_tables):
            df = pair_tables[k]
            for g, sc, p, grp, col in zip(df[ref + "_gene"], df[ref + "_scaf"], df[ref + "_pos"], df["gene_group"], df["color"]):
                if not is_unassigned(col, grp):
                    v = ranks[(str(grp), sc)]
                    t[(ref, g)] = bisect.bisect_left(v, float(p)) / max(1, len(v) - 1)
    changed = True
    while changed:   # carry the shades along the chain of neighboring pairs
        changed = False
        for k, df in enumerate(pair_tables):
            a, b = species_order[k], species_order[k + 1]
            for ga, gb in zip(df[a + "_gene"], df[b + "_gene"]):
                if (a, ga) in t and (b, gb) not in t:
                    t[(b, gb)] = t[(a, ga)]
                    changed = True
                elif (b, gb) in t and (a, ga) not in t:
                    t[(a, ga)] = t[(b, gb)]
                    changed = True
    return t


def paint_halfway(panel, x0, width, chrom_len, pts, y_top, height, zorder=4):
    """Each point paints the bar halfway to its neighbors (used for gene-order shading)."""
    pts = sorted(pts)
    if not pts:
        return
    xs = [x0 + width * p / chrom_len for p, _ in pts]
    quads = []
    for k, x in enumerate(xs):
        a = x0 if k == 0 else (xs[k - 1] + x) / 2
        b = x0 + width if k == len(xs) - 1 else (x + xs[k + 1]) / 2
        quads.append([(a, y_top), (b, y_top), (b, y_top + height), (a, y_top + height)])
    cols = [c for _, c in pts]
    panel.add_collection(PolyCollection(quads, facecolors=cols, edgecolors=cols, linewidths=0.05,
                                        zorder=zorder, clip_on=False))


# ---------------------------------------------------------------------------
# follow a gene set
# ---------------------------------------------------------------------------
def parse_follow(spec):
    """spec: a list of strings or a path to a file of lines. 'SP:gene' follows one gene;
    'SP:scaffold:start-end' follows every ortholog of SP in that region.
    Returns (genes {species: set}, regions [(species, scaffold, start, end)])."""
    lines = spec
    if isinstance(spec, str):
        with open(spec) as fh:
            lines = [ln.strip() for ln in fh]
    genes, regions = collections.defaultdict(set), []
    for ln in lines:
        ln = str(ln).strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split(":")
        if len(parts) == 3 and "-" in parts[2]:
            lo, hi = (float(v.replace(",", "")) for v in parts[2].split("-"))
            regions.append((parts[0], parts[1], lo, hi))
        elif len(parts) == 2:
            genes[parts[0]].add(parts[1])
        else:
            raise ValueError("follow: cannot read %r (use SP:gene or SP:scaffold:start-end)" % ln)
    return genes, regions


def follow_mask(df, a, b, genes, regions):
    """True for rows of a pair table (species a, b) that are followed."""
    out = []
    ga, gb = genes.get(a, set()), genes.get(b, set())
    for gene_a, gene_b, sa, pa, sb, pb in zip(df[a + "_gene"], df[b + "_gene"], df[a + "_scaf"], df[a + "_pos"],
                                              df[b + "_scaf"], df[b + "_pos"]):
        hit = gene_a in ga or gene_b in gb
        if not hit:
            for sp, sc, lo, hi in regions:
                if (sp == a and sc == sa and lo <= float(pa) <= hi) or (sp == b and sc == sb and lo <= float(pb) <= hi):
                    hit = True
                    break
        out.append(hit)
    return out


# ---------------------------------------------------------------------------
# microsynteny and collinear (riparian) runs
# ---------------------------------------------------------------------------
def _ranks_along(df, sp):
    by = collections.defaultdict(list)
    for i, (sc, p) in enumerate(zip(df[sp + "_scaf"], df[sp + "_pos"])):
        by[sc].append((float(p), i))
    rank, chains = {}, []
    for v in by.values():
        v.sort()
        chains.append([i for _, i in v])
        for k, (_, i) in enumerate(v):
            rank[i] = k
    return rank, chains


def microsynteny_rows(df, a, b, window=3, min_run=3):
    """Row positions (0..n-1) in runs of >= min_run orthologs that are neighbors in a and stay within `window`
    ortholog ranks on one chromosome of b."""
    rank_b, _ = _ranks_along(df, b)
    _, chains = _ranks_along(df, a)
    scaf_b = list(df[b + "_scaf"])
    keep = set()
    for ch in chains:
        run = ch[:1]
        for u, w in zip(ch, ch[1:]):
            if scaf_b[u] == scaf_b[w] and abs(rank_b[u] - rank_b[w]) <= window:
                run.append(w)
            else:
                if len(run) >= min_run:
                    keep.update(run)
                run = [w]
        if len(run) >= min_run:
            keep.update(run)
    return keep


def collinear_runs(df, a, b, min_run=3):
    """Riparian bands only where two stretches are collinear: runs of >= min_run orthologs that are neighbors in a
    (no other ortholog of the pair between them) and stay neighbors, in the same or reversed order, on one chromosome
    of b. Returns [(row positions, reversed)]."""
    rank_b, _ = _ranks_along(df, b)
    _, chains = _ranks_along(df, a)
    scaf_b = list(df[b + "_scaf"])
    runs = []
    for ch in chains:
        run, step = ch[:1], 0
        for u, w in zip(ch, ch[1:]):
            d = rank_b[w] - rank_b[u] if scaf_b[u] == scaf_b[w] else 0
            if abs(d) == 1 and (step == 0 or d == step):
                run.append(w)
                step = d
            else:
                if len(run) >= min_run:
                    runs.append((list(run), step < 0))
                run, step = [w], 0
        if len(run) >= min_run:
            runs.append((list(run), step < 0))
    return runs


def band_patch(a0, a1, b0, b1, y0, y1, reversed_=False, **kw):
    """A riparian band from [a0, a1] at y0 to [b0, b1] at y1, twisted when reversed."""
    ym = (y0 + y1) / 2
    if reversed_:
        b0, b1 = b1, b0
    path = Path([(a0, y0), (a0, ym), (b0, ym), (b0, y1), (b1, y1), (b1, ym), (a1, ym), (a1, y0), (a0, y0)],
                [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4,
                 Path.CURVE4, Path.CLOSEPOLY])
    return PathPatch(path, **kw)


def micro_color(spec, group_color):
    """microsynteny_color: 'dark' (the linkage-group color darkened), 'gray', or any matplotlib color."""
    if spec == "dark":
        return shade(group_color if group_color and str(group_color).upper() != UNASSIGNED else "#8c8c8c", 0.85)
    if spec == "gray":
        return "#4d4d4d"
    return spec


def shuffled(n, seed=1):
    order = list(range(n))
    random.Random(seed).shuffle(order)
    return order


# ---------------------------------------------------------------------------
# chromosome names on the bars
# ---------------------------------------------------------------------------
def short_chrom_names(names):
    """Drop a prefix that all of a genome's chromosome names share when it holds no digits (PMA3 -> 3, chr10 -> 10,
    Ch_6 -> 6) and a shared version suffix (.1), but only when every shortened name is four characters or fewer.
    Accession numbers such as NC_047017.1 or CAIIXF020000008.1 stay whole, so no accession fragment passes for a
    chromosome number."""
    import re
    names = [str(n) for n in names]
    if len(names) < 2:
        return {n: n for n in names}
    pre = names[0]
    for n in names[1:]:
        while not n.startswith(pre):
            pre = pre[:-1]
    while pre and pre[-1].isdigit():
        pre = pre[:-1]
    if any(ch.isdigit() for ch in pre):
        pre = ""
    suf = names[0]
    for n in names[1:]:
        while not n.endswith(suf):
            suf = suf[1:]
    m = re.search(r"\.\d+$", suf)
    suf = m.group(0) if m else ""
    out = {}
    for n in names:
        core = n[len(pre):len(n) - len(suf)] if suf else n[len(pre):]
        if core.isdigit():
            core = core.lstrip("0") or "0"
        out[n] = core if core else n
    if len(set(out.values())) < len(out) or any(len(v) > 4 for v in out.values()):
        return {n: n for n in names}   # collisions, or accession-like names: keep the names as they are
    return out


def name_on_bar(ax, x_center, y, text, fontsize=6, zorder=8):
    """A chromosome name centered on its bar: black with a thin white edge, readable over painted colors."""
    import matplotlib.patheffects as pe
    ax.text(x_center, y, text, ha="center", va="center", fontsize=fontsize, color="black", zorder=zorder,
            clip_on=False, path_effects=[pe.withStroke(linewidth=1.2, foreground="white")])
