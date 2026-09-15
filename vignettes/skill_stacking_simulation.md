---
name: haplostack-stacking-simulation
description: >
  Run haploStack recurrent breeding (stacking) simulations from a founder
  GenoMeta and marker effects — a standalone single-method run, or a comparison
  across selection schemes (TS / OHS / GBS; OCS has known issues), under the
  parents-first unit or direct pair selection — and plot/summarise the
  trajectories. Use when a user wants to project genetic gain under a
  selection strategy.
---

# Skill — Haplotype stacking simulation

Runs select → design → cross → advance over cycles and replications, scoring each
cycle on local GEBV. Full walkthrough: `vignettes/example_stacking_simulation.Rmd`.

## When to use
- You have a founder `GenoMeta` (with blocks) and marker effects.
- You want to project gain for one method, or compare several.

## Selection schemes — focus on TS, OHS, GBS
| Scheme | Meaning | Status |
|--------|---------|--------|
| `TS`  | Truncation / genomic selection (top founders by total GEBV) | ✔ |
| `OHS` | Optimal Haplotype Selection (maximise next-cross block gain) | ✔ |
| `GBS` | Genotype Building Strategy (drive block frequencies to an ideal) | ✔ |
| `OCS` | kinship-penalised truncation — **known issues**: not true optimum-contribution selection, silent truncation fall-back, replicates identical in a multi-rep call; a KCT rename + a real OPT scheme are planned | ⚠ avoid |

Every scheme runs under `selection_unit = "parents"` (default: pick a parent
subgroup, cross all pairs) or `"pairs"` (pick the families directly; needs
`pair_budget = "families"` or `"parents"`). Full combination matrix:
`vignettes/selection_grid.md`.

## Standalone run — one method
```r
meta <- load_geno_meta("founders.rds")

sim <- run_stacking_simulation(
  geno_meta     = meta,
  effect_df     = effects,          # or rely on meta$marker_effects
  scheme        = "OHS",            # single method
  ohs_method    = "Greedy",         # GA_pairwise / GA_coverage / Greedy
  replications  = 3L,
  select_cycles = 5L,
  top_founders  = 6L,
  num_offspring = 4L,
  cross_backend = "block",          # pure-R; no extra packages
  seed          = 42L,
  output_dir    = NULL, verbose = FALSE)

sim$summary           # per cycle × rep
sim$trajectories      # plotting view: Rep, Cycle, Method, MeanGEBV, ...
plot_simulation_trajectories(sim$trajectories, trait_name = "Yield", show_ribbon = TRUE)
```

## Assessment across methods
Pass a **vector** of schemes — same founders, cycles, and seeds for a fair compare.
```r
cmp <- run_stacking_simulation(
  geno_meta = meta, effect_df = effects,
  scheme = c("TS", "OHS", "GBS"),
  ohs_method = "Greedy", gbs_method = "ARCSIN",
  replications = 3L, select_cycles = 5L, top_founders = 6L,
  num_offspring = 4L, cross_backend = "block", seed = 42L,
  output_dir = NULL, verbose = FALSE)

plot_simulation_trajectories(cmp$trajectories, trait_name = "Yield", colour_by = "Method")
summ <- summarise_simulation(cmp$trajectories, group_by = c("Cycle", "Method"))
summ[Cycle == max(Cycle)][order(-MeanGEBV_mean)]   # rank methods at the final cycle
```

## Direct pair selection (0.7.2)
```r
pr <- run_stacking_simulation(
  geno_meta = meta, effect_df = effects,
  scheme = c("TS", "OHS", "GBS"),
  selection_unit = "pairs",
  pair_budget    = "families",     # same offspring population size as parents-first
                                   # ("parents": same number of distinct parents)
  ts_pair_score  = "midparent",    # or "usefulness" (+ usefulness_k * sd(F1))
  gbs_pair_level = "population",   # or "family"; pair_penalty = "loss" | "hard"
  ohs_method = "Greedy", replications = 3L, select_cycles = 5L,
  top_founders = 4L, num_offspring = 4L, output_dir = NULL, verbose = FALSE)
pr$summary[, .(cycle, sel_scheme, n_families, n_parents_used)]   # counts per cycle
pr$plans[, .(sel_scheme, Parent1, Parent2, pair_score)]           # score per family
pair_selection_budget(16, "families")                               # the numbers a budget implies
```

## Gotchas
- `top_founders` must be ≤ the number of founders; keep `replications` /
  `select_cycles` small while tuning, then scale up.
- Default `cross_backend = "block"` inherits whole parental blocks (no
  intra-block recombination — an upper bound on stacking gain). Install
  *genomicSimulation* and set `cross_backend = "genomicSimulation"` for
  marker-level recombination.
- `select_priority = "low"` minimises the trait instead of maximising.
- In one call with `replications > 1`, `OHS` GA engines and `OCS` re-seed the
  RNG every cycle, so their replicates are identical — use
  `ohs_method = "Greedy"`, or separate calls with different seeds.
- `"pairs"`: `pair_budget` is required; `OCS` and `ohs_method = "GA_coverage"`
  have no family form and are refused.
- On small datasets GBS may emit benign penalty warnings.
- Parallelise large comparisons with `n_cores > 1` (needs *future* /
  *future.apply*).
