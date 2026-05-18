#!/usr/bin/env python3
"""
Polarity-evidence ribbon panel.

A 3-row ribbon plot showing how a hypothesized fusion / fission event
maps across three species:

    top row     unfused-side species: two distinct chromosomes carrying
                two ancestral linkage groups (e.g. ALG-A and ALG-B).
    middle row  outgroup species: the two chromosomes that preserve
                ALG-A and ALG-B as distinct linkage groups.
    bottom row  fused-side species: the single chromosome carrying
                both ALG-A and ALG-B.

Each ribbon traces one gene via 3-way RBH (top <-> middle <-> bottom).
Ribbons can be colored by:
  * `bcns`    - BCnS-Simakov2022 ALG (uses HMM-RBH calls)
  * `lg2`     - two-color split keyed to the top-row chromosomes
  * `raw`     - single neutral color

The module exposes `draw_polarity_panel(ax, ...)` so callers can compose
multiple panels into a grid (e.g. one outgroup per column, one
fused-side species per row).
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
import matplotlib.patches as mpatches
import pandas as pd


# -- defaults / constants ----------------------------------------------------

# Two-color split for the `lg2` coloring scheme.
COLOR_LG2_A = "#B7202A"   # ALG-A side (top-chrom #1)
COLOR_LG2_B = "#1F4A8F"   # ALG-B side (top-chrom #2)
FALLBACK_GREY = "#888888"

# Deep-saturated label colors picked to NOT collide with the BCnS palette
# or the lg2 two-color split. Used to tie the in-cell prefix to its
# matching grid label.
LABEL_COLOR_BOTTOM = "#8B0000"   # DarkRed: bottom-row prefix + row label
LABEL_COLOR_MIDDLE = "#00264D"   # Deep navy: middle-row prefix + col header
LABEL_COLOR_TOP    = "#000000"   # black (neutral)

ROW_Y = {"top": 0, "middle": 1, "bottom": 2}


# -- I/O helpers -------------------------------------------------------------

def load_rbh(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_bcns_map(hmm_rbh_path: str | Path, gene_col: str) -> dict:
    """Map {species_gene -> BCnS ALG name} from an odp BCnS HMM-RBH file."""
    p = Path(hmm_rbh_path)
    if not p.exists():
        return {}
    df = load_rbh(p)
    if gene_col not in df.columns or "BCnSSimakov2022_scaf" not in df.columns:
        return {}
    return dict(zip(df[gene_col], df["BCnSSimakov2022_scaf"]))


# -- string helpers ----------------------------------------------------------

def chrom_acronym(chrom: str) -> str:
    """Compact label for a chromosome ID. Strips accession version and
    common GenBank-style prefixes, keeping the distinguishing suffix."""
    s = chrom.rstrip()
    s = re.sub(r"\.\d+$", "", s)
    for pat, fmt in (
        (r"^HiC_scaffold_(\d+)$", "HiC{0}"),
        (r"^NW_\d+(\d{3})$",       "NW{0}"),
        (r"^NC_\d+(\d{2})$",       "NC{0}"),
        (r"^OZ\d+(\d{2})$",        "OZ{0}"),
        (r"^CM\d+(\d{2})$",        "CM{0}"),
        # Generic "*_chrN" / "*_chromN" — picks up haplotype-resolved
        # assemblies like eupHapAv0.3_chr1, Foo_chr14, asm_chrom5, etc.
        (r"_chr(?:om)?(\d+)$",     "chr{0}"),
    ):
        m = re.search(pat, s) if pat.startswith("_") else re.match(pat, s)
        if m:
            return fmt.format(m.group(1))
    return s


# -- layout helpers ----------------------------------------------------------

def layout_chroms(chroms: list[str], widths: dict, x_left: float, x_right: float,
                  gap_frac: float = 0.04) -> dict:
    """Return {chrom: (x_start, x_end)} sized by `widths` with gaps."""
    total = sum(widths.values()) or 1
    span = x_right - x_left
    gap_total = gap_frac * span * max(0, len(chroms) - 1)
    chrom_span = span - gap_total
    out = {}
    x = x_left
    for i, c in enumerate(chroms):
        w = chrom_span * widths[c] / total
        out[c] = (x, x + w)
        x += w + (gap_frac * span if i < len(chroms) - 1 else 0)
    return out


def gene_x(chrom: str, pos_index: int, n_on_chrom: int, layout: dict) -> float:
    x0, x1 = layout[chrom]
    if n_on_chrom <= 1:
        return (x0 + x1) / 2
    return x0 + (pos_index / (n_on_chrom - 1)) * (x1 - x0)


def bezier_ribbon(ax, x_top, y_top, x_bot, y_bot, color, alpha=0.4, lw=0.4):
    dy = y_bot - y_top
    verts = [(x_top, y_top), (x_top, y_top + dy * 0.5),
             (x_bot, y_top + dy * 0.5), (x_bot, y_bot)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    patch = mpatches.PathPatch(MplPath(verts, codes), fill=False, lw=lw,
                                edgecolor=color, alpha=alpha, capstyle="round")
    ax.add_patch(patch)


def _index_within_chrom(df: pd.DataFrame, scaf_col: str, pos_col: str,
                         idx_col: str, n_col: str) -> pd.DataFrame:
    df = df.copy()
    df[idx_col] = (df.groupby(scaf_col)[pos_col].rank(method="first")
                     .astype(int) - 1)
    df[n_col]   = df.groupby(scaf_col)[pos_col].transform("count")
    return df


# -- core rendering function -------------------------------------------------

def draw_polarity_panel(ax, *,
                         top_sp: str, top_chrA: str, top_chrB: str,
                         mid_sp: str, mid_chrA: str, mid_chrB: str,
                         bot_sp: str, bot_chr_fused: str,
                         rbh_top_mid: str | Path,
                         rbh_bot_mid: str | Path,
                         color_by: str = "bcns",
                         bcns_hmm_top: str | Path | None = None,
                         bcns_hmm_mid: str | Path | None = None,
                         bcns_hmm_bot: str | Path | None = None,
                         bcns_palette: str | Path | None = None,
                         top_abbrev: str | None = None,
                         mid_abbrev: str | None = None,
                         bot_abbrev: str | None = None,
                         show_chrom_labels: bool = True,
                         show_species_labels: bool = True,
                         show_legend: bool = True,
                         legend_outside: bool = False,
                         chrom_acronyms: bool = False,
                         chrom_label_map: dict | None = None,
                         chrom_gap_frac: float = 0.04,
                         equal_width_chroms: bool = False,
                         highlight_bottom_species: bool = False,
                         highlight_middle_species: bool = False,
                         label_fontsize: int = 7,
                         species_fontsize: int = 10,
                         lg_a_label: str = "ALG-A",
                         lg_b_label: str = "ALG-B") -> dict:
    """
    Draw one polarity-evidence panel into an existing matplotlib axes.

    Returns a dict with summary stats (n_triplets, n_unmapped).

    Required arguments
    ------------------
    top_sp / top_chrA / top_chrB
        Top-row species and its two chromosomes carrying ALG-A and ALG-B.
    mid_sp / mid_chrA / mid_chrB
        Middle (outgroup) species and the two chromosomes that preserve
        ALG-A and ALG-B as distinct linkage groups.
    bot_sp / bot_chr_fused
        Bottom species and the single chromosome carrying both ALGs.
    rbh_top_mid / rbh_bot_mid
        Paths to the pairwise reciprocal-best-hits TSV files
        (top<->mid and bot<->mid).
    """
    # -------------------------------------------------------------------
    # Build 3-way RBH triplets
    # -------------------------------------------------------------------
    TOP, MID, BOT = top_sp, mid_sp, bot_sp
    TOP_GENE, TOP_SCAF, TOP_POS = f"{TOP}_gene", f"{TOP}_scaf", f"{TOP}_pos"
    MID_GENE, MID_SCAF, MID_POS = f"{MID}_gene", f"{MID}_scaf", f"{MID}_pos"
    BOT_GENE, BOT_SCAF, BOT_POS = f"{BOT}_gene", f"{BOT}_scaf", f"{BOT}_pos"

    tm = load_rbh(rbh_top_mid)
    bm = load_rbh(rbh_bot_mid)
    triplets = tm.merge(bm, on=MID_GENE, suffixes=("_tm", "_bm"))
    triplets[MID_SCAF] = triplets[f"{MID_SCAF}_tm"]
    triplets[MID_POS]  = triplets[f"{MID_POS}_tm"]
    triplets = triplets[
        triplets[TOP_SCAF].isin([top_chrA, top_chrB]) &
        triplets[BOT_SCAF].isin([bot_chr_fused]) &
        triplets[MID_SCAF].isin([mid_chrA, mid_chrB])
    ].copy()
    n_total = len(triplets)
    if n_total == 0:
        ax.text(0.5, 0.5, f"no 3-way triplets\n{top_sp}/{mid_sp}/{bot_sp}",
                ha="center", va="center", fontsize=7, color="0.4",
                transform=ax.transAxes)
        for s in ("top", "right", "bottom", "left"):
            ax.spines[s].set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])
        return {"n_triplets": 0, "n_unmapped": 0}

    top_chroms = [top_chrA, top_chrB]
    mid_chroms = [mid_chrA, mid_chrB]
    bot_chroms = [bot_chr_fused]

    def widths_for(col, chroms):
        if equal_width_chroms:
            return {c: 1 for c in chroms}
        v = triplets[col].value_counts()
        return {c: int(v.get(c, 0)) for c in chroms}

    top_layout = layout_chroms(top_chroms, widths_for(TOP_SCAF, top_chroms),
                                0.18, 0.82, gap_frac=chrom_gap_frac)
    mid_layout = layout_chroms(mid_chroms, widths_for(MID_SCAF, mid_chroms),
                                0.18, 0.82, gap_frac=chrom_gap_frac)
    bot_layout = layout_chroms(bot_chroms, widths_for(BOT_SCAF, bot_chroms),
                                0.18, 0.82, gap_frac=chrom_gap_frac)

    triplets = triplets.sort_values([TOP_SCAF, TOP_POS]).reset_index(drop=True)
    triplets = _index_within_chrom(triplets, TOP_SCAF, TOP_POS, "_top_idx", "_top_n")
    triplets = _index_within_chrom(triplets, MID_SCAF, MID_POS, "_mid_idx", "_mid_n")
    triplets = _index_within_chrom(triplets, BOT_SCAF, BOT_POS, "_bot_idx", "_bot_n")

    # -------------------------------------------------------------------
    # Coloring
    # -------------------------------------------------------------------
    legend_elems = []
    if color_by == "lg2":
        def color_for_row(r):
            return COLOR_LG2_A if r[TOP_SCAF] == top_chrA else COLOR_LG2_B
        legend_elems = [
            mpatches.Patch(facecolor=COLOR_LG2_A, edgecolor="none", label=lg_a_label),
            mpatches.Patch(facecolor=COLOR_LG2_B, edgecolor="none", label=lg_b_label),
        ]
    elif color_by == "bcns":
        if bcns_palette is None:
            raise ValueError("bcns_palette path is required when color_by='bcns'")
        top_map = load_bcns_map(bcns_hmm_top, f"{TOP}_gene") if bcns_hmm_top else {}
        mid_map = load_bcns_map(bcns_hmm_mid, f"{MID}_gene") if bcns_hmm_mid else {}
        bot_map = load_bcns_map(bcns_hmm_bot, f"{BOT}_gene") if bcns_hmm_bot else {}
        triplets["_alg_top"] = triplets[TOP_GENE].map(top_map)
        triplets["_alg_mid"] = triplets[MID_GENE].map(mid_map)
        triplets["_alg_bot"] = triplets[BOT_GENE].map(bot_map)
        triplets["_alg"] = (triplets["_alg_top"]
                              .fillna(triplets["_alg_mid"])
                              .fillna(triplets["_alg_bot"]))
        palette_df = load_rbh(bcns_palette)
        alg_to_color = (palette_df[["gene_group", "color"]]
                          .drop_duplicates().dropna()
                          .set_index("gene_group")["color"].to_dict())

        def color_for_row(r):
            return alg_to_color.get(r["_alg"], FALLBACK_GREY)

        algs_used = triplets["_alg"].dropna().value_counts().index.tolist()
        legend_elems = [
            mpatches.Patch(facecolor=alg_to_color.get(a, FALLBACK_GREY),
                            edgecolor="none", label=f"BCnS-{a}")
            for a in algs_used
        ]
        n_unmapped = triplets["_alg"].isna().sum()
        if n_unmapped > 0:
            legend_elems.append(
                mpatches.Patch(facecolor=FALLBACK_GREY, edgecolor="none",
                                label=f"no BCnS call (n={n_unmapped})"))
    else:   # raw
        def color_for_row(r):
            return "#444444"
        legend_elems = []

    # -------------------------------------------------------------------
    # Axes setup
    # -------------------------------------------------------------------
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.4, len(ROW_Y) - 1 + 0.4)
    ax.set_ylim(ax.get_ylim()[::-1])
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])

    def draw_chroms(layout: dict, y: float, chrom_labels: dict | None):
        for c, (x0, x1) in layout.items():
            ax.plot([x0, x1], [y, y], color="black", lw=2.0,
                    solid_capstyle="round", clip_on=False)
            if chrom_labels is None or chrom_labels.get(c) is None:
                continue
            label = chrom_labels[c]
            ax.text(x0 + 0.002, y - 0.04, label,
                    ha="left", va="bottom",
                    rotation=90, rotation_mode="anchor",
                    fontsize=label_fontsize, zorder=200, clip_on=False)

    if show_chrom_labels:
        _override = chrom_label_map or {}
        def _lab(c):
            if c in _override:
                return _override[c]
            return chrom_acronym(c) if chrom_acronyms else c
        top_labels = {top_chrA: _lab(top_chrA), top_chrB: _lab(top_chrB)}
        mid_labels = {mid_chrA: _lab(mid_chrA), mid_chrB: _lab(mid_chrB)}
        bot_labels = {bot_chr_fused: _lab(bot_chr_fused)}
    else:
        top_labels = mid_labels = bot_labels = None

    draw_chroms(top_layout, ROW_Y["top"],    top_labels)
    draw_chroms(mid_layout, ROW_Y["middle"], mid_labels)
    draw_chroms(bot_layout, ROW_Y["bottom"], bot_labels)

    if show_species_labels:
        ta = top_abbrev or top_sp[:3]
        ma = mid_abbrev or mid_sp[:3]
        ba = bot_abbrev or bot_sp[:3]
        species_rows = [
            (ROW_Y["top"],    ta, "normal", LABEL_COLOR_TOP),
            (ROW_Y["middle"], ma,
                "bold" if highlight_middle_species else "normal",
                LABEL_COLOR_MIDDLE if highlight_middle_species else LABEL_COLOR_TOP),
            (ROW_Y["bottom"], ba,
                "bold" if highlight_bottom_species else "normal",
                LABEL_COLOR_BOTTOM if highlight_bottom_species else LABEL_COLOR_TOP),
        ]
        for y, name, weight, color in species_rows:
            ax.text(0.04, y, name, ha="right", va="center",
                    style="italic", weight=weight, color=color,
                    fontsize=species_fontsize,
                    clip_on=False, zorder=210)

    # -------------------------------------------------------------------
    # Ribbons
    # -------------------------------------------------------------------
    if color_by == "bcns":
        draw_order = pd.concat([triplets[triplets["_alg"].isna()],
                                 triplets[triplets["_alg"].notna()]])
        ribbon_alpha = lambda r: 0.45 if pd.isna(r.get("_alg")) else 0.65
        ribbon_lw    = lambda r: 0.55 if pd.isna(r.get("_alg")) else 0.55
    else:
        draw_order = triplets
        ribbon_alpha = lambda r: 0.5
        ribbon_lw    = lambda r: 0.5

    for _, row in draw_order.iterrows():
        color = color_for_row(row)
        tx = gene_x(row[TOP_SCAF], row["_top_idx"], row["_top_n"], top_layout)
        mx = gene_x(row[MID_SCAF], row["_mid_idx"], row["_mid_n"], mid_layout)
        bx = gene_x(row[BOT_SCAF], row["_bot_idx"], row["_bot_n"], bot_layout)
        a, l = ribbon_alpha(row), ribbon_lw(row)
        bezier_ribbon(ax, tx, ROW_Y["top"],    mx, ROW_Y["middle"], color, alpha=a, lw=l)
        bezier_ribbon(ax, mx, ROW_Y["middle"], bx, ROW_Y["bottom"], color, alpha=a, lw=l)

    if show_legend and legend_elems:
        ncol = min(len(legend_elems), 7)
        ax.legend(handles=legend_elems, loc="lower center",
                  bbox_to_anchor=(0.5, -0.10), ncol=ncol,
                  frameon=False, fontsize=label_fontsize)

    n_unmapped = (int(triplets["_alg"].isna().sum())
                   if color_by == "bcns" else 0)
    return {"n_triplets": n_total, "n_unmapped": n_unmapped}
