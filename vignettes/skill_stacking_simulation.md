---
name: haplostack-stacking-simulation
description: >
  Run haploStack recurrent breeding (stacking) simulations from a founder
  GenoMeta and marker effects — a standalone single-method run, or a comparison
  across selection schemes (TS/OHS/OCS/GBS) — and plot/summarise the trajectories.
  Use when a user wants to project genetic gain under a selection strategy.
---

# Skill — Haplotype stacking simulation

Runs select → design → cross → advance over cycles and replications, scoring each
cycle on local GEBV. Full walkthrough: `vignettes/example_stacking_simulation.Rmd`.

## When to use
- You have a founder `GenoMeta` (with blocks) and marker effects.
- You want to project gain for one method, or compare several.

## Selection schemes
| Scheme | Meaning |
|--------|---------|
| `TS`  | Truncation / genomic selection (top founders by total GEBV) |
| `OHS` | Optimal Haplotype Selection (maximise next-cross block gain) |
| `OCS` | Optimal Contribution Selection (kinship-constrained) |
| `GBS` | Genotype Building Strategy (drive block frequencies to an ideal) |

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
  scheme = c("TS", "OHS", "OCS", "GBS"),
  ohs_method = "Greedy", gbs_method = "ARCSIN",
  replications = 3L, select_cycles = 5L, top_founders = 6L,
  num_offspring = 4L, cross_backend = "block", seed = 42L,
  output_dir = NULL, verbose = FALSE)

plot_simulation_trajectories(cmp$trajectories, trait_name = "Yield", colour_by = "Method")
summ <- summarise_simulation(cmp$trajectories, group_by = c("Cycle", "Method"))
summ[Cycle == max(Cycle)][order(-MeanGEBV_mean)]   # rank methods at the final cycle
```

## Gotchas
- `top_founders` must be ≤ the number of founders; keep `replications` /
  `select_cycles` small while tuning, then scale up.
- Default `cross_backend = "block"` inherits whole parental blocks (no
  intra-block recombination — an upper bound on stacking gain). Install
  *genomicSimulation* and set `cross_backend = "genomicSimulation"` for
  marker-level recombination.
- `select_priority = "low"` minimises the trait instead of maximising.
- On small datasets OCS/GBS may emit benign kinship/penalty warnings.
- Parallelise large comparisons with `n_cores > 1` (needs *future* /
  *future.apply*).
