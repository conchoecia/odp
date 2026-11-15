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


# ---------------------------------------------------------------------------
# extra outputs: zoom, small multiples, animation, interactive page
# They use the chromosome-coordinate layout of the main plot (x from 0 to 1).
# ---------------------------------------------------------------------------
def layout_from(species_order, sp_to_chromdf, rbh_df_list):
    """{sp: {chrom: (x0, width, size)}} in the chromosome-coordinate panel, and the pair tables in order."""
    chroms = {}
    for sp in species_order:
        df = sp_to_chromdf[sp]
        chroms[sp] = {c: (float(o), float(w), float(s)) for c, o, w, s in
                      zip(df["chrom"], df["chrPlotOffset"], df["chrPlotPercent"], df["chrsize"])}
    return {"species": list(species_order), "chroms": chroms, "pairs": list(rbh_df_list)}


def _x(lay, sp, scaf, pos):
    x0, w, size = lay["chroms"][sp][scaf]
    return x0 + w * float(pos) / max(size, 1.0)


def _rounded_bar(ax, x, y, w, h, color):
    r = min(h / 2, w / 2) * 0.7
    ax.add_patch(FancyBboxPatch((x + r, y - h / 2 + r), max(w - 2 * r, 1e-9), max(h - 2 * r, 1e-9),
                                boxstyle="round,pad=%f,rounding_size=%f" % (r, r), facecolor=color, edgecolor="none",
                                zorder=5))


def zoom_figure(path, lay, top, bottom, names=None, width_in=7.087):
    """The pair top=(sp, chrom), bottom=(sp, chrom) on one Mb scale, one ribbon per ortholog between them."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import PathCollection
    (sa, ca), (sb, cb) = top, bottom
    k = next((i for i in range(len(lay["species"]) - 1) if {lay["species"][i], lay["species"][i + 1]} == {sa, sb}), None)
    if k is None:
        raise ValueError("zoom: %s and %s must be neighboring species rows" % (sa, sb))
    df = lay["pairs"][k]
    la, lb = lay["chroms"][sa][ca][2], lay["chroms"][sb][cb][2]
    mm = 1 / 25.4
    W, H = width_in / mm, 46.0
    fig = plt.figure(figsize=(width_in, H * mm))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    x0 = 32.0
    scale = (W - 4 - x0) / max(la, lb)
    ya, yb = 12.0, 36.0
    paths, cols = [], []
    for ga, sca, pa, scb, pb, grp, col in zip(df[sa + "_gene"], df[sa + "_scaf"], df[sa + "_pos"], df[sb + "_scaf"],
                                               df[sb + "_pos"], df["gene_group"], df["color"]):
        if sca == ca and scb == cb:
            xa, xb = x0 + float(pa) * scale, x0 + float(pb) * scale
            ym = (ya + yb) / 2
            paths.append(Path([(xa, ya + 0.45), (xa, ym), (xb, ym), (xb, yb - 0.45)],
                              [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]))
            cols.append(mcolors.to_rgba("#8c8c8c", 0.3) if is_unassigned(col, grp) else mcolors.to_rgba(col, 0.6))
    ax.add_collection(PathCollection(paths, facecolors="none", edgecolors=cols, linewidths=0.5, zorder=3))
    for y, sp, c, L, up in ((ya, sa, ca, la, True), (yb, sb, cb, lb, False)):
        _rounded_bar(ax, x0, y, L * scale, 0.9, "black")
        step = 10 if L > 40e6 else 5 if L > 15e6 else 2 if L > 6e6 else 1
        n = int(L / 1e6 // step)
        for j in range(n + 1):
            x = x0 + j * step * 1e6 * scale
            yy = y - 0.45 if up else y + 0.45
            ax.plot([x, x], [yy, yy + (-0.7 if up else 0.7)], color="black", linewidth=0.5)
            ax.text(x, yy + (-1.0 if up else 1.0), "%d" % (j * step) + (" Mb" if j == n else ""), ha="center",
                    va="bottom" if up else "top", fontsize=7)
        ax.text(x0 - 2, y, (names or {}).get(sp, sp), ha="right", va="center", fontsize=7, style="italic")
        ax.text(x0 - 2, y + 3.0, c, ha="right", va="center", fontsize=7)
    fig.savefig(path)
    plt.close(fig)
    return len(paths)


def small_multiples_figure(path, lay, colors, ncol=6, width_in=7.087, names=None):
    """One small panel per linkage group, with only that group's ribbons; chromosomes holding at least 5 of its
    orthologs in black, the others light gray."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import PathCollection
    groups = sorted(colors)
    sps = lay["species"]
    mm = 1 / 25.4
    W = width_in / mm
    name_w = 22.0
    pw = (W - name_w) / ncol
    gap = 8.5
    ph = 6.0 + gap * (len(sps) - 1) + 3.5
    nrow = -(-len(groups) // ncol)
    fig = plt.figure(figsize=(width_in, nrow * ph * mm))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(nrow * ph, 0)
    ax.axis("off")
    for k, g in enumerate(groups):
        cx, cy = name_w + (k % ncol) * pw, (k // ncol) * ph
        px = lambda x: cx + 0.6 + x * (pw - 1.8)
        ys = {sp: cy + 6.0 + gap * i for i, sp in enumerate(sps)}
        paths, cols = [], []
        held = collections.defaultdict(collections.Counter)
        for i, df in enumerate(lay["pairs"]):
            a, b = sps[i], sps[i + 1]
            for sca, pa, scb, pb, grp in zip(df[a + "_scaf"], df[a + "_pos"], df[b + "_scaf"], df[b + "_pos"], df["gene_group"]):
                if str(grp) != g or sca not in lay["chroms"][a] or scb not in lay["chroms"][b]:
                    continue
                held[a][sca] += 1
                held[b][scb] += 1
                xa, xb = px(_x(lay, a, sca, pa)), px(_x(lay, b, scb, pb))
                y0, y1 = ys[a] + 0.35, ys[b] - 0.35
                ym = (y0 + y1) / 2
                paths.append(Path([(xa, y0), (xa, ym), (xb, ym), (xb, y1)], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]))
                cols.append(mcolors.to_rgba(colors[g], 0.45))
        ax.add_collection(PathCollection(paths, facecolors="none", edgecolors=cols, linewidths=0.5, zorder=2))
        for sp in sps:
            for c, (x0, w, size) in lay["chroms"][sp].items():
                _rounded_bar(ax, px(x0), ys[sp], max(px(x0 + w) - px(x0), 0.2), 0.7,
                             "black" if held[sp][c] >= 5 else "#c8c8c8")
        bw = 1.8 + 1.4 * len(g)
        ax.add_patch(FancyBboxPatch((cx + 0.6 + 0.5, cy + 1.0 + 0.5), bw - 1.0, 1.6, boxstyle="round,pad=0.5,rounding_size=0.5",
                                    facecolor=colors[g], edgecolor="none", zorder=5))
        ax.text(cx + 0.6 + bw / 2, cy + 2.3, g, ha="center", va="center", fontsize=7, color=text_color_for(colors[g]), zorder=6)
        if k % ncol == 0:
            for sp in sps:
                full = (names or {}).get(sp, sp).split()
                short = "%s. %s" % (full[0][0], " ".join(full[1:])) if len(full) >= 2 else " ".join(full)
                ax.text(name_w - 2, ys[sp], short, ha="right", va="center", fontsize=7, style="italic")
    fig.savefig(path)
    plt.close(fig)
    return len(groups)


def _chains(lay):
    """Orthologs found in every species, followed through the neighboring pair tables by gene ID:
    [{sp: (scaf, pos, color)}]."""
    sps = lay["species"]
    links = []
    for i, df in enumerate(lay["pairs"]):
        a, b = sps[i], sps[i + 1]
        links.append({ga: (gb, sa, pa, sb, pb, col, grp) for ga, gb, sa, pa, sb, pb, col, grp in
                      zip(df[a + "_gene"], df[b + "_gene"], df[a + "_scaf"], df[a + "_pos"], df[b + "_scaf"],
                          df[b + "_pos"], df["color"], df["gene_group"])})
    out = []
    for g0, (g1, sa, pa, sb, pb, col, grp) in links[0].items():
        rec = {sps[0]: (sa, pa), sps[1]: (sb, pb)}
        color = "#8c8c8c" if is_unassigned(col, grp) else col
        cur = g1
        ok = True
        for i in range(1, len(links)):
            if cur not in links[i]:
                ok = False
                break
            nxt, _, _, sb2, pb2, col2, grp2 = links[i][cur]
            rec[sps[i + 1]] = (sb2, pb2)
            cur = nxt
        if ok and all(rec[sp][0] in lay["chroms"][sp] for sp in sps):
            out.append((rec, color))
    return out


def animation(path_mp4, lay, background="light", names=None, fps=30, hold_s=1.4, move_s=3.0):
    """Dots for the orthologs found in every genome travel from one genome's chromosomes to the next (1920 by 1080),
    on a light or dark background; writes an MP4 with ffmpeg and a GIF beside it."""
    import shutil
    import subprocess
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    if shutil.which("ffmpeg") is None:
        raise RuntimeError("animate needs ffmpeg on the PATH")
    sps = lay["species"]
    sel = _chains(lay)
    bg, fg = ("#111111", "white") if background == "dark" else ("white", "black")
    W, H = 406.4, 228.6
    fx = lambda x: 20.0 + x * (W - 40.0)
    yc = 150.0
    rng = np.random.default_rng(1)
    x = np.array([[fx(_x(lay, sp, *rec[sp])) for sp in sps] for rec, _ in sel])
    y = yc - 3.5 - rng.uniform(0, 1, size=(len(sel), len(sps))) * 26.0
    cols = [mcolors.to_rgb(c) for _, c in sel]
    fig = plt.figure(figsize=(16, 9), dpi=120, facecolor=bg)
    ax = fig.add_axes([0, 0, 1, 1], facecolor=bg)
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    bars, labels = [], []
    for sp in sps:
        pats = [FancyBboxPatch((fx(x0) + 1.0, yc - 0.5), max(fx(x0 + w) - fx(x0) - 2.0, 1e-6), 1.0,
                               boxstyle="round,pad=1.0,rounding_size=1.0") for x0, w, _ in lay["chroms"][sp].values()]
        pc = PatchCollection(pats, facecolor=fg, edgecolor="none", zorder=3)
        ax.add_collection(pc)
        bars.append(pc)
        labels.append(ax.text(W / 2, 190, (names or {}).get(sp, sp), ha="center", va="center", fontsize=30, style="italic", color=fg))
    sc = ax.scatter(x[:, 0], y[:, 0], s=30, c=cols, linewidths=0, zorder=4)
    ax.text(W / 2, 30, "{:,} orthologs found in every genome".format(len(sel)), ha="center", va="center", fontsize=20, color=fg)
    plan = []
    for k in range(len(sps)):
        plan.append(("hold", k, int(hold_s * fps)))
        if k < len(sps) - 1:
            plan.append(("move", k, int(move_s * fps)))
    ff = subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgba", "-s", "1920x1080",
                           "-r", str(fps), "-i", "-", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "28", "-preset", "slow", path_mp4],
                          stdin=subprocess.PIPE)
    frames = 0
    for kind, k, n in plan:
        for j in range(n):
            if kind == "hold":
                px_, py_ = x[:, k], y[:, k]
                al = [1.0 if i == k else 0.0 for i in range(len(sps))]
                nl = al
            else:
                t = j / max(n - 1, 1)
                e = t * t * (3 - 2 * t)
                dx = x[:, k + 1] - x[:, k]
                px_ = x[:, k] + dx * e
                py_ = y[:, k] + (y[:, k + 1] - y[:, k]) * e - np.minimum(np.abs(dx) * 0.22, 70.0) * math.sin(math.pi * e)
                al = [0.0] * len(sps)
                al[k], al[k + 1] = 1 - e, e
                nl = [max(0.0, 2 * v - 1) for v in al]
            sc.set_offsets(np.column_stack([px_, py_]))
            for i in range(len(sps)):
                bars[i].set_alpha(al[i])
                labels[i].set_alpha(nl[i])
            fig.canvas.draw()
            ff.stdin.write(bytes(fig.canvas.buffer_rgba()))
            frames += 1
    ff.stdin.close()
    ff.wait()
    plt.close(fig)
    gif = path_mp4.rsplit(".", 1)[0] + ".gif"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", path_mp4, "-vf",
                    "fps=10,scale=720:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=64[p];[s1][p]paletteuse=dither=none", gif], check=True)
    return {"orthologs_in_every_genome": len(sel), "frames": frames}


_PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Ribbon plot</title>
<style>
:root{--bg:#f4f4f2;--ink:#1a1a1a;--muted:#5c5c5c;--card:#fff;--line:#d6d6d2;--tip:#1a1a1a;--tipink:#fff}
@media (prefers-color-scheme:dark){:root{--bg:#151515;--ink:#ececec;--muted:#a8a8a8;--line:#333;--tip:#f2f2f2;--tipink:#111}}
body{margin:0;background:var(--bg);color:var(--ink);font-family:Helvetica,Arial,sans-serif}
main{max-width:1180px;margin:0 auto;padding:14px 16px 24px}
.help{font-size:13px;color:var(--muted)} button{font:inherit;font-size:13px;padding:3px 10px;border:1px solid var(--line);border-radius:6px;background:#fff;color:#1a1a1a;cursor:pointer}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:6px;overflow-x:auto}
svg{display:block;width:100%;min-width:640px;height:auto} svg text{font-family:Helvetica,Arial,sans-serif;font-size:2.47px;fill:#000}
#rb path{fill:none;stroke-width:0.176} #hit path{fill:none;stroke:transparent;stroke-width:0.9;pointer-events:stroke}
.bar{fill:#000;cursor:pointer} .lab{cursor:pointer}
#tip{position:fixed;pointer-events:none;background:var(--tip);color:var(--tipink);font-size:12.5px;line-height:1.35;padding:6px 9px;border-radius:6px;display:none}
</style></head><body><main>
<p><b>Ribbon plot</b> <span class="help">Point at a ribbon for its genes. Click a label or a chromosome to spotlight it; Esc or Reset clears.</span> <button id="reset">Reset</button></p>
<div class="card"><svg id="fig" viewBox="0 0 @@W@@ @@H@@"><g id="rb">@@RIB@@</g><g>@@BARS@@</g><g>@@NAMES@@</g><g>@@LABS@@</g><g id="hit">@@HITS@@</g></svg></div>
</main><div id="tip"></div>
<script>
const R=@@RECS@@;const rb=[...document.querySelectorAll('#rb path')];const tip=document.getElementById('tip');let spot=null,hot=null;
function apply(){for(const p of rb){const b=+p.dataset.o;p.setAttribute('stroke-opacity',spot?(spot(p)?Math.max(.55,b):.015):b)}}
function show(k,x,y){const r=R[k];tip.innerHTML='<b>'+(r[0]?'Linkage group '+r[0]:'No linkage group')+'</b><br>'+r[1]+'<br>'+r[2];tip.style.display='block';tip.style.left=Math.min(innerWidth-tip.offsetWidth-8,x+14)+'px';tip.style.top=Math.max(8,y-tip.offsetHeight-12)+'px'}
document.getElementById('hit').addEventListener('mousemove',e=>{const k=e.target.dataset.k;if(k===undefined)return;if(hot!==null)document.getElementById('r'+hot).style.strokeWidth='';hot=k;const p=document.getElementById('r'+k);p.style.strokeWidth=.45;p.setAttribute('stroke-opacity',1);show(+k,e.clientX,e.clientY)});
document.getElementById('hit').addEventListener('mouseleave',()=>{tip.style.display='none';if(hot!==null){document.getElementById('r'+hot).style.strokeWidth='';apply()}hot=null});
document.querySelectorAll('.lab').forEach(g=>g.addEventListener('click',()=>{const x=g.dataset.g;spot=p=>p.dataset.g===x;apply()}));
document.querySelectorAll('.bar').forEach(b=>b.addEventListener('click',()=>{const c=b.dataset.c;spot=p=>p.dataset.a===c||p.dataset.b===c;apply()}));
document.getElementById('reset').onclick=()=>{spot=null;apply()};document.addEventListener('keydown',e=>{if(e.key==='Escape'){spot=null;apply()}});
</script></body></html>"""


def interactive_page(path, lay, colors, alpha_sig=0.4, alpha_faint=0.05, names=None):
    """The chromosome-coordinate plot as a self-contained web page with gene tooltips and spotlights."""
    import html
    import json
    sps = lay["species"]
    W, gap, top, lab_col = 180.0, 24.0, 16.0, 30.0
    X = lambda x: lab_col + x * (W - lab_col - 2.0)
    ys = {sp: top + gap * i for i, sp in enumerate(sps)}
    rib, hits, recs = [], [], []
    k = 0
    for i, df in enumerate(lay["pairs"]):
        a, b = sps[i], sps[i + 1]
        for ga, gb, sa, pa, sb, pb, grp, col, fet in zip(df[a + "_gene"], df[b + "_gene"], df[a + "_scaf"], df[a + "_pos"],
                                                         df[b + "_scaf"], df[b + "_pos"], df["gene_group"], df["color"],
                                                         df["break_FET"]):
            if sa not in lay["chroms"][a] or sb not in lay["chroms"][b]:
                continue
            un = is_unassigned(col, grp)
            xa, xb = X(_x(lay, a, sa, pa)), X(_x(lay, b, sb, pb))
            y0, y1 = ys[a] + 0.45, ys[b] - 0.45
            ym = (y0 + y1) / 2
            d = "M%.2f %.2fC%.2f %.2f %.2f %.2f %.2f %.2f" % (xa, y0, xa, ym, xb, ym, xb, y1)
            o = alpha_faint if (un or not float(fet) <= 0.05) else alpha_sig
            rib.append('<path id="r%d" d="%s" stroke="%s" stroke-opacity="%.3f" data-o="%.3f" data-g="%s" data-a="%s" data-b="%s"/>'
                       % (k, d, "#8c8c8c" if un else col, o, o, "" if un else html.escape(str(grp)),
                          html.escape(a + ":" + sa), html.escape(b + ":" + sb)))
            hits.append('<path d="%s" data-k="%d"/>' % (d, k))
            recs.append(["" if un else str(grp), "%s %s, %.2f Mb (%s)" % ((names or {}).get(a, a), sa, float(pa) / 1e6, ga),
                         "%s %s, %.2f Mb (%s)" % ((names or {}).get(b, b), sb, float(pb) / 1e6, gb)])
            k += 1
    bars = ['<rect class="bar" x="%.2f" y="%.2f" width="%.2f" height="0.9" rx="0.3" data-c="%s"/>'
            % (X(x0), ys[sp] - 0.45, X(x0 + w) - X(x0), html.escape(sp + ":" + c))
            for sp in sps for c, (x0, w, _) in lay["chroms"][sp].items()]
    nms = ['<text x="%.2f" y="%.2f" text-anchor="end" font-style="italic">%s</text>'
           % (lab_col - 2, ys[sp] + 0.85, html.escape((names or {}).get(sp, sp))) for sp in sps]
    labs, x_at = [], lab_col
    for g in sorted(colors):
        bw = 1.8 + 1.5 * len(g)
        labs.append('<g class="lab" data-g="%s"><rect x="%.2f" y="2" width="%.2f" height="2.6" rx="0.5" fill="%s"/>'
                    '<text x="%.2f" y="3.85" text-anchor="middle" fill="%s">%s</text></g>'
                    % (html.escape(g), x_at, bw, colors[g], x_at + bw / 2, text_color_for(colors[g]), html.escape(g)))
        x_at += bw + 0.6
    page = _PAGE
    for key, val in (("@@W@@", "%.1f" % W), ("@@H@@", "%.1f" % (ys[sps[-1]] + 5)), ("@@RIB@@", "".join(rib)),
                     ("@@BARS@@", "".join(bars)), ("@@NAMES@@", "".join(nms)), ("@@LABS@@", "".join(labs)),
                     ("@@HITS@@", "".join(hits)), ("@@RECS@@", json.dumps(recs))):
        page = page.replace(key, val)
    with open(path, "w") as fh:
        fh.write(page)
    return k
