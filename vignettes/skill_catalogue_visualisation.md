---
name: haplostack-catalogue-visualisation
description: >
  Compute per-block local GEBVs from a haploStack GenoMeta and a trait's marker
  effects, then visualise them (effect / variance Manhattan, haplotype alignment
  on subsets) or run the full catalogue orchestrator. Use when a user wants to
  score and visualise haplotype-block effects for a trait.
---

# Skill — Haplotype catalogue and visualisation

Turns a `GenoMeta` + marker effects into per-sample, per-block local GEBVs and
diagnostic plots. Full walkthrough: `vignettes/example_catalogue_visualisation.Rmd`.

## When to use
- You have a `GenoMeta` (with blocks) and effect estimates for a trait.
- You want per-block breeding values and genome-wide / alignment plots.

## Inputs
- A valid `GenoMeta` with `block_marker_map_info` (attach blocks first if absent).
- Marker effects: `Marker, Allele, Effect` — stored on `meta$marker_effects` or
  passed as `effect_df`.

## Local GEBVs
```r
meta <- load_geno_meta("catalogue_demo.rds")     # usual hand-off point
gebv <- compute_profile_localgebv(meta, level = "genotype", center_by_freq = TRUE)
# gebv: wide Block × Sample data.table (first column "Block")
```
`level = "haplotype"` keeps the two haplotypes separate; `center_by_freq = TRUE`
centres on the population mean.

## Genome-wide Manhattan plots
```r
bloc <- as.data.table(build_block_location_info(meta)$block_info)  # Block,Chrom,Pos,...
long <- data.table::melt(gebv, id.vars = "Block",
                         variable.name = "Sample", value.name = "Effect")
long <- merge(long, bloc[, .(Block, Chrom, Pos)], by = "Block")

plot_manhattan_effect(long, trait_name = "Yield")                 # needs Chrom,Pos,Effect
var_dt <- long[, .(BlockVariance = var(Effect, na.rm = TRUE)), by = .(Block, Chrom, Pos)]
plot_manhattan_variance(var_dt, trait_name = "Yield")             # needs Chrom,Pos,BlockVariance
```

## Haplotype alignment on a subset
Subset first — whole-genome alignment is dense.
```r
sub <- filter_geno_meta(meta, chromosomes = list(list(chr = "1")))
sub <- filter_geno_meta(meta, keep_samples = elite, blocks = block_ids)
plot_haplotype_alignment(sub, show_block_boundaries = TRUE, show_block_labels = TRUE)
```

## One-call catalogue (multi-trait + file export)
```r
results <- run_haplotype_catalogue(
  geno_meta   = meta,
  effect_list = list(Yield = eff_yield, Protein = eff_protein),  # named per trait
  output_dir  = "catalogue_out", save_profile = TRUE, n_cores = 1L)

results$Yield$gebv_matrix        # per-trait outputs
plot_block_counts(results$Yield$block_counts, trait_name = "Yield")

# Per-block effect summary is returned ready-made (and written to
# Block_Summary_<ld>_t<tol>_Yield.txt) — no manual variance derivation needed.
bs <- results$Yield$block_summary   # Block, BlockVariance, BlockMax, BlockMin, BlockMean, Chrom, Pos
plot_manhattan_variance(bs, trait_name = "Yield")   # variance scatter, one plot per trait
```

## Gotchas
- `plot_manhattan_effect()` needs columns `Chrom, Pos, Effect`;
  `plot_manhattan_variance()` needs `Chrom, Pos, BlockVariance`. The two share
  one internal genomic-scatter body, differing only in that value column.
- Building plots by hand (the melt + `var()` above) is only needed off the raw
  `gebv` matrix; `run_haplotype_catalogue()` already returns `block_summary`
  (with `BlockVariance/BlockMax/BlockMin/BlockMean`) per trait.
- `run_haplotype_catalogue()` errors if the `GenoMeta` has no
  `block_marker_map_info` — call `attach_block_info()` first.
- The GEBV matrix feeds directly into the stacking simulation (Example 4).
