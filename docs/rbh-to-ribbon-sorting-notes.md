# Ribbon plot chromosome ordering — algorithm notes

Working notes from the iterative development of the chromosome-ordering and
flip-optimization pipeline for `odp_rbh_to_ribbon` / `odp_ribbon_plot.py`.
This file is a design / debugging log, not user-facing documentation.
See odp issue #127 and PR #128 for the discussion context.

## Objective

For a ribbon plot of N species stacked top-to-bottom, choose:

1. A linear ordering of each species' chromosomes within its row.
2. A per-chromosome forward/reverse (L/R) orientation.

To minimize the total number of bezier-line crossings between every
adjacent species pair. FET-significant (`whole_FET <= 0.05`) orthologs
count 1000× more than non-significant ones, so the score is a
multiplicative weighted inversion count.

## Final pipeline

The `optimize_chrom_rotation = True` path (default) runs:

1. **`optimal-lg` cluster seed** — `_order_by_lg` partitions every
   chromosome into a synteny linkage cluster using a mutual-top-1
   FET-significant graph (`_build_fet_lg_clusters`), with orphan
   rescue for chromosomes that were the runner-up partner of a
   shared hub. Within each cluster, `_seed_top_order_within` and
   `_cascade_within` apply the FET-best-partner heuristic restricted
   to that cluster.
2. **Flip cascade** — `_optimize_chrom_flips_top_down` decides each
   chromosome's L/R orientation greedily top-down, with stage-1
   independent-per-scaffold seed and stage-2 iterated flip-greedy
   refinement. Optionally seeded from `initial_flip` to preserve
   prior decisions on row 0.
3. **Brushing sweeps** — `_brushing_sweep` alternates two passes
   per generation:
   - **td (top-down)**: row 0 fixed; each row k below is re-ordered
     by FET-weighted barycenter of its partners in row k-1
     (currently pinned).
   - **bu (bottom-up)**: row N-1 fixed; each row k above is
     re-ordered by FET-weighted barycenter of its partners in row k+1.
   After each direction the flip cascade is re-run with
   `initial_flip` set to the current state. Top-row and bottom-row
   greedy moves (tm / bm) are tried after td and bu respectively,
   prioritised by per-chromosome tension. Best-seen state tracked
   across all sweeps. Patience = 4 consecutive non-improving sweeps
   before convergence.
4. **(Optional) iterated random restart** — perturb the best state
   (random shuffle of a window of chromosomes in a few rows), brush
   again, keep global best. 8 restarts × ~3 min each on a 12-species
   Lepidoptera test.

The final layout is rendered by `_render_ribbon_figure` (extracted
from `ribbon_plot` so brushing can snapshot per generation).

## Score / scoring kernel

For every adjacent species pair, every line endpoint pair (i < j)
that crosses (i.e. their top-x and bot-x orderings flip) adds
`w_i × w_j` to the score. Per-line weight = `1000 × shared_count`
if the chromosome pair's best `whole_FET <= 0.05`, else
`shared_count`. Implemented as a Fenwick-tree O(n log n) weighted
inversion count (`_fenwick_weighted_inversion`).

Two scoring resolutions:
- **Anchor coarse** (`_build_pair_anchor_data`): one record per
  (top_chrom, bot_chrom) pair at the within-chrom mean rank. ~200
  records/pair on Lep data. Used inside the brushing sweep for
  speed.
- **Per-ortholog fine** (`_build_pair_line_data`): one record per
  individual ortholog. ~2000 records/pair on Lep. Slower; used by
  the local-search refinement when it was enabled.

## What worked

| Technique | Effect |
|---|---|
| FET-graph LG clustering with within-cluster cascade | Clean LG-aligned columns across all species; correct seed for everything else. |
| Brushing td + bu alternation, FET-weighted barycenter | -66% from seed on the Lepidoptera test in the cleanest run. |
| FET 1000× weighting | Excluded trace homology from the centroid; kept LG signal dominant. |
| Best-state tracking, patience = 4 | Lets bu pass through 1-2 regressions to find a better next td. Critical. |
| Flip cascade per direction with `initial_flip` | Re-evaluated L/R in fresh context without losing previously-good decisions. |
| Top/bottom row tension-prioritised greedy moves | Marginal (~1-5M per stroke) but consistent. |
| Iterated random restart with small perturbation | -4% additional vs single descent; plateaued after 4 restarts. |
| `flush=True` on prints + per-gen PDF snapshots | Live convergence visibility on long runs. |
| Robust two-line labels (italic species + grey accession) | Fallback to id-strip if no tip_label_map. |

## What didn't work / was marginal

| Technique | Reason |
|---|---|
| Raw barycenter (no FET weighting, no LG seed) | Trace homology drowned the signal — WORSE than the legacy cascade. |
| Median heuristic; barycenter-iter (top-down only) | Same problem; cascade or brushing always beat them. |
| Pure local search (insertion + swap + flip) on crossing count | Greedy, plateaued at ~996B vs 360B for brushing. Can't propagate. |
| Adjacent-swap polish (`optimal-swap`) only | Too local. |
| Bidirectional iteration with one-anchor `_optimize_spA_based_on_rbh` cascade | Converges after 1-2 sweeps; can't escape the same local minimum the one-pass cascade hits. |
| Sticky LG-block move with looser top-k FET filter | Bridged distinct LGs via shared fusion hubs; graph collapsed to one giant component. |
| Strict convergence on first regression | Stopped at gen 2; missed deeper minima behind transient regressions. |
| Top-row feedback only (bu disabled) | 931B vs 360B for bu+td. Top-row moves only directly affect pair (0,1). |
| Flip cascade after `try_top_moves` without `initial_flip` | Reset row-0 flips, undoing accepted moves. Fixed by adding `initial_flip` to `_optimize_chrom_flips_top_down`. |
| Small-window perturbation (5 chroms × 3 rows) for restart | Plateaued at 374B after 4 successful restarts. Too small to escape the current basin. |
| First refactor of the render block via Python script | Accidentally deleted `ribbon_plot`. Restored from git; second attempt with sanity checks worked. |

## Score trajectory (Lepidoptera 12-species test, w = 1000 for FET-sig)

| Approach | Best score | Δ from `optimal-lg` seed |
|---|---|---|
| `optimal-size` legacy cascade (replaced) | ~1069B | -- |
| `optimal-lg` cluster seed | ~1069B | (seed for brushing) |
| td-only brushing | 945B | -12% |
| Pure local search on per-ortholog score | 996B | -7% |
| td + top-row feedback (no bu) | 931B | -13% |
| **td + bu brushing, patience=4** | **360B** | **-66%** ← best |
| 4-phase brushing (td + tm + bu + bm) | 389B | -64% |
| 4-phase + 8 iterated restarts | 374B | -65% |

## Key insights for the next round

1. **bu is essential.** No technique we tried can replace it. The
   alternating top-down / bottom-up anchor flip is what enables
   deep search.
2. **Track best-seen, not last-iter.** bu sweeps regress the score
   sometimes; their value is in the *next* td that follows.
3. **Patience > strict monotone.** Allow several non-improving
   sweeps before quitting; some of the best minima sit two or three
   sweeps past a regression.
4. **Small perturbations plateau.** To find truly different basins
   we likely need bigger restart perturbations (row reverse,
   multi-row shuffle) or SA-style uphill acceptance during brushing.
5. **Sugiyama needs modification.** Vanilla barycenter fails on
   data with universal trace homology. The FET-weighted seed +
   FET-cluster constraint + per-direction flip cascade is what
   makes it work here.
6. **Visibility matters.** Per-gen PDF snapshots + flushed logs
   turned algorithm debugging from blind to surgical.

## Configuration knobs

In `_brushing_sweep` (defaults shown):

- `max_iters = 50` — generations of (td + tm + bu + bm)
- `patience = 4` — consecutive non-improving sweeps before stop
- `top_moves_per_iter = 15` — how many top-row chroms to try per gen
- `top_move_offsets = (1, 2, 3, 5)` — swap distances
- `enable_bottom_up = True` — toggle bu sweep
- `bottom_moves_per_iter = 15` — symmetric to top_moves

In iterated restart wrapper (inside `ribbon_plot`):

- `n_restarts = 8`
- `perturb_rows = 3`
- `perturb_chroms_per_row = 5`
- random seed = 42 (fixed for reproducibility)

## Open problems

1. **Plateau at ~374B**. The small-perturbation restart can't escape
   it. Need larger perturbation, or SA, or a fundamentally different
   move type.
2. **bu regressions cost compute**. Every regression is a wasted sweep
   that the patience mechanism has to absorb. A predictive heuristic
   for which direction to take next would save time.
3. **bottom-up bias on this dataset**. Every bu pass on the Lep test
   regressed relative to the preceding td. Reason: Plutella (top
   row, basal) is the cleanest anchor; Danaus (bottom, ~23 chroms
   from fusions) is a poor anchor. May need row-quality weighting
   in barycenter.
4. **Local search after brushing**. We had a per-ortholog crossing-
   count local search that converged to ~1145B from the brushing
   seed, then dropped to 1110B with anchor coarsening and insertion
   moves. Currently disabled to isolate brushing. Should be re-
   enabled once the brushing+restart side stabilizes; it catches
   within-chrom detail brushing flattens.
5. **Multi-start from divergent seeds**. Try brushing from
   `optimal-size`, `optimal-lg`, `random` initial orderings;
   keep best. Cheap parallelization (each seed is independent).

## File layout

- `scripts/odp_ribbon_plot.py` — all the helpers above and
  `ribbon_plot` (the entry point).
- `scripts/odp_rbh_to_ribbon` — the snakefile; one `make_plot` rule
  that calls `ribbon_plot` with config-supplied knobs.
- This file: `docs/rbh-to-ribbon-sorting-notes.md`.
