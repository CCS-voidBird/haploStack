# Stacking simulation - choosing the settings

`run_stacking_simulation()` = select -> design crosses -> cross -> advance,
for `select_cycles` cycles x `replications`, scoring each cycle on local GEBV.

## Decide with the user first

1. **Trait** - one per run: `effect_df = meta$marker_effects[["X"]]` or
   `trait = "X"`. With several stored traits and neither, the call stops and
   names them (0.9.0; earlier versions took the first silently). For several
   traits build an index with `create_fitness_index()` and pass that.
2. **Schemes** - the default is `c("TS", "OHS", "GBS")` (0.9.0); OCS runs only
   when named. Pass `scheme =` anyway so the choice is on record.
3. **Unit** - `selection_unit = "parents"` (choose `top_founders`, cross every
   pair: C(n, 2) families) or `"pairs"` (choose families directly;
   `pair_budget` required).
4. **Direction** - `select_priority = "high"` or `"low"`. `"low"` is
   honoured by TS and by GBS with `construct_method = "argmax"`; OHS, OCS and
   GBS with the `"custom"` / `"ratio"` builder refuse it (they select the high
   tail whatever it says).
5. **Crossing backend** - `cross_backend = "block"` (whole blocks inherited;
   an upper bound on stacking) or `"genomicSimulation"` (marker-level
   recombination; needs that package).
6. **Budget** - replications x cycles; rehearse with 2 x 3 first.
7. **Seed** - fixed and recorded (`seed =`). Seeded calls leave the caller's
   random stream as it was.

## Schemes

| Scheme | Objective | Main options | Status |
|---|---|---|---|
| `TS` | truncation on total local GEBV; pairs: mid-parent or usefulness | `ts_pair_score`, `usefulness_k` | recommended |
| `OHS` | maximise the chosen set's family stacking ceilings | `ohs_method = "GA_pairwise"` (default) or `"Greedy"`; `ohs_elite_ratio` in (0, 1] with `ohs_elite_level = "pair"` or `"block"` | recommended |
| `GBS` | drive block frequencies towards an ideal genotype | `gbs_method = "ARCSIN"` or `"LOG"`, `gbs_preset = "haploStack"` or `"kemper"`, `ideal_scope`, `penalty_mode`; pairs: `gbs_pair_level`, `pair_penalty`, `pair_freq` | recommended |
| `OCS` | kinship-penalised truncation | `max_kinship` | avoid: no pair form; not in the default; its 0.8.8 numbers cannot be reproduced (a per-cycle reseed was removed in 0.9.0) |

`avoid_selfing` defaults to `TRUE`. `top_founders` above the population size
is capped.

## Pair budgets (`selection_unit = "pairs"`)

`pair_selection_budget(top_founders, pair_budget)` prints the numbers.

- `"families"` - the same offspring population as the parents-first run
  (C(top_founders, 2) families); the number of distinct parents is free.
- `"parents"` - the same number of distinct parents; the family count is free.

## Metrics and cost

`metrics` is measurement only: `"gebv"` is always on; add `"d2ideal"`,
`"gene_diversity"`, `"coancestry"`, `"ideal_hap_gebv"`, or `"all"`.
`n_cores > 1` runs replicates x schemes in parallel (needs future +
future.apply). `monitor = TRUE` writes live tables and figures (needs
`output_dir`). OHS is the slowest scheme per cycle, TS the fastest.

`output_dir` defaults to `NULL` (0.9.0): nothing is written and the result
comes back in memory. Name a directory to get `stacking_plans_<tag>.txt`,
`stacking_summary_<tag>.txt` and `Trajectories_<tag>.png`, or write later with
`write_outputs(sim, dir)`. `top_founders` above the population is an error.
`progress = function(e) ...` is called after every generation (`e$payload$summary`
is that generation's row); returning `FALSE` cancels with an `hs_error_cancelled`
whose `partial` field holds what was done.

## Step API (0.9.0)

`run_stacking_simulation()` is a wrapper over five calls that give the same
numbers: `spec <- stk_spec(...)` (every science argument, validated once),
`stk_jobs(spec)` (the replication x scheme grid), `st <- stk_init(spec, meta,
effect_df, job = k)`, `st <- stk_step(st)` one generation at a time until
`st$status != "running"`, `stk_result(states)`; `stk_close(st)` frees a
genomicSimulation session. A `"block"` state survives `saveRDS()` and resumes
identically; a genomicSimulation state does not (one gS session per R process:
`gs_load_geno_meta()` refuses a second load unless `replace = TRUE`).

## Reading results

```r
sim$summary[cycle == max(cycle),
            .(gain = mean(mean_F1mean), elite = mean(max_F1mean)), by = sel_scheme]
plot_simulation_trajectories(sim$trajectories, trait_name = "X", show_ribbon = TRUE)
plot_simulation_metrics(sim$summary, metrics = c("mean_F1mean", "mean_gene_diversity"))
```

- The result is an `hs_simulation`: `plans`, `summary`, `ideal_genotypes`,
  `final_populations`, `trajectories`, plus `hs_notes(sim)` (every silent
  fallback of the run, e.g. `fallback_truncation`, `ocs_kinship_cap`,
  `metric_failed`; one warning per code is raised at the end) and
  `hs_provenance(sim)` (the resolved spec, seed, backend, one row per job).
  Always read `hs_notes(sim)` before reporting numbers.
- `mean_F1mean` = mean offspring local GEBV; `max_F1mean` = the best single
  offspring per cycle (the top elite).
- One replicate has no ribbon and no variance; never present it as expected.
- Compare schemes at matched cycles: truncation typically leads early, fixes,
  and is overtaken later by the stacking schemes.
- `summarise_simulation(sim$trajectories)` summarises `MeanGEBV` and `SDev`;
  the all-`NA` `BlockCoverage` column of the compat view is skipped (0.9.0).
