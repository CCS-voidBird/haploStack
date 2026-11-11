---
name: haplostack
description: Usage guide for the haploStack R package (haplotype-block analysis and haplotype-stacking breeding simulation). Build a GenoMeta container from PLINK or delimited genotype files with a haplotype-block map (or LD blocks via hapS_ld_blocks), attach marker effects, compute local GEBVs and a block catalogue, count favourable blocks, and simulate TS / OHS / GBS selection over breeding cycles. Use whenever a task mentions haploStack, GenoMeta, local GEBV, haplotype blocks or catalogue, haplotype stacking, OHS / GBS / OCS selection, run_haplotype_catalogue or run_stacking_simulation.
---
<!-- haplostack-version -->

# haploStack - usage guide

haploStack turns **genotypes + a marker map + a haplotype-block map + marker
effects** into **local GEBVs per block**, a **block catalogue**, and
**stacking simulations** (selection over breeding cycles).

This guide ships inside the package, so it describes the build that is
installed. Print it any time with `haploStack::hs_agent_guide()`; the
reference pages are `hs_agent_guide("blocks")`, `"genometa"`, `"simulation"`
and `"troubleshooting"`. Every exported function has a help page: `?fn`,
`args(fn)`.

## Step 0 - orient yourself

```r
library(haploStack)            # the startup line names this guide
packageVersion("haploStack")
?haploStack                    # package page
```

Before touching real data, run the simulated route below once (seconds). It
proves the install works and shows every object's shape.

## The object: GenoMeta

A plain named list; nearly every function takes one.

```r
meta$geno_arrays            # 3-D character array [samples x markers x ploidy] - NOT a list
meta$map_info               # Chrom, Marker, Pos (bp)
meta$block_marker_map_info  # Chrom, Marker, Block  (catalogue and simulation need it)
meta$block_info             # Chrom, Block, Start, End, Pos, Num_SNPs, Block_length
meta$sample_info            # Family, Sample, Parent1, Parent2, Sex  (PLINK FAM)
meta$marker_effects         # named list: one (Marker, Allele, Effect) table per trait
meta$phenotype              # Sample + trait columns
meta$metadata               # ploidy, encoding, missing_tokens, blocks (provenance), ...
```

Inspect with `dim(meta$geno_arrays)`, `str(meta, max.level = 1)` and
`is_geno_meta(meta, verbose = TRUE)`. Never print the whole object. Details:
`hs_agent_guide("genometa")`.

## The route (tested on simulated data)

```r
library(haploStack)
options(device = function(...) grDevices::pdf(NULL))   # no stray Rplots.pdf in scripts
OUT <- "haplostack_run"; dir.create(OUT, showWarnings = FALSE)

# 0. a simulated panel, to rehearse (40 lines x 120 markers, blocks of 4)
d    <- sim_dataset(n_samples = 40L, n_markers = 120L, block_size = 4L, seed = 1L)
meta <- d$geno_meta
eff  <- sim_effect_df(meta$map_info$Marker, seed = 2L)

# 1. read real data WITH the block map (the recommended route; no LD parameter)
# meta <- read_plink_ped("panel.ped", "panel.map", block_file = "blocks.tsv")
#   or   build_geno_meta(geno_2d, map, blocks = block_df)   # rows = samples, cols = markers
#   or   read_genotype(geno_file, map_file, block_file = ...)
#   no block map at all? LD blocks are a side route: meta <- hapS_ld_blocks(meta)
#   (needs the HapSelect package; see hs_agent_guide("blocks"))

# 2. marker effects: from a file ...
# eff <- read_effect_file("effects.tsv")            # separator and header detected
meta$marker_effects <- list(Yield = eff)
stopifnot(is_geno_meta(meta))
#    ... or fitted from a phenotype (needs rrBLUP / BGLR):
# fit <- fit_genomic_prediction(meta, phenotype = ph, trait = "Yield", model = "rrblup")
# meta$marker_effects <- fit$marker_effects

# 3. phenotype (optional)
# meta <- attach_phenotype(meta, read_phenotype("pheno.tsv", sample_col = "IID"))
save_geno_meta(meta, file.path(OUT, "geno_meta.rds"))

# 4. catalogue: per-trait tables and Manhattan PNGs under output_dir
res <- run_haplotype_catalogue(meta, output_dir = file.path(OUT, "catalogue"),
                               center = "none")
res$Yield$block_summary            # Chrom, Pos, Block, BlockVariance, BlockMax, ...
lg  <- compute_profile_localgebv(meta, trait = "Yield", level = "genotype",
                                 center = "none")     # Block x samples
fav <- count_favourable_blocks(meta, trait = "Yield", prefer = "positive",
                               center = "none", plot = FALSE)
fav$table                          # per sample: GEBV, nPositive, nNegative, ...

# 5. stacking simulation: ONE trait, schemes named explicitly
sim <- run_stacking_simulation(meta, effect_df = meta$marker_effects$Yield,
                               scheme = c("TS", "OHS", "GBS"),
                               replications = 2L, select_cycles = 3L,
                               top_founders = 8L, seed = 42L,
                               metrics = c("gebv", "d2ideal", "gene_diversity"),
                               output_dir = file.path(OUT, "simulation"))
sim$summary[cycle == max(cycle),
            .(gain = mean(mean_F1mean), elite = mean(max_F1mean)), by = sel_scheme]
plot_simulation_trajectories(sim$trajectories, trait_name = "Yield")
plot_simulation_metrics(sim$summary, metrics = "mean_F1mean")
```

Rehearse a real run small (one chromosome, 2 replications x 3 cycles) before
the full budget. Simulation choices: `hs_agent_guide("simulation")`.

Genes on the blocks (0.9.0, offline, no Python):

```r
meta  <- annotate_blocks(meta, gff3 = "genes.gff3.gz",      # path or table per input
                         chr_map = c("1" = "1H"),           # container -> GFF3 spelling, if they differ
                         flank_bp = 0)                      # 0 = exact overlap; default 10 kb
paths <- export_block_annotation(meta, file.path(OUT, "annotation"))   # the five step files
plot_meta_circos(list(track_annotation(meta)), meta = meta)
```

Coordinates must be bp. A container chromosome the GFF3 never names is
warned about with the `chr_map` that would match (`strict = TRUE` stops).

## Rules that change results

1. **A block map you already have is the main route.** Give it at read time
   (`block_file =`) or with `attach_block_info()`. It must be one of two
   shapes, with these header names: a block-marker map `Chrom, Marker, Block`
   or a block-region map `Chrom, Block, Start, End`. Other headers (`Chr, SNP,
   Block`, no header) are refused with an `hs_error_block_map`. Rename columns
   first.
2. **LD blocks are a side route, and the engine decides the answer.**
   `hapS_ld_blocks()` (HapSelect, the default) and `st_ld_blocks()`
   (SelectionTools, being retired, off unless `use_st_formats(TRUE)`) give
   different partitions at the same settings. Use one engine per study and say
   which made the blocks.
3. **Centring.** `center = "none"` (default) | `"locus_mean"` (HapSelect's
   `mean_adjust`) | `"allele_freq"` (reproduces catalogues from <= 0.5.5).
   Their numbers are not comparable with each other; state which you used.
4. **The simulation is single-trait.** Pass `effect_df =` or `trait =` for
   the trait you mean; with several stored traits and neither, the call stops
   and names them (0.9.0). A multi-trait goal needs an index first
   (`create_fitness_index()`).
5. **Pass `scheme =`.** The default is `c("TS", "OHS", "GBS")` (0.9.0); OCS
   runs only when named and is not recommended (no pair form; see the
   simulation page). `output_dir` defaults to `NULL` (nothing written; use
   `write_outputs(sim, dir)` or name a directory).
6. **`metrics` never changes selection**; it only adds measurements.
7. **`ld_threshold` / `tolerance` on `run_haplotype_catalogue()` are labels**
   used in output file names (default `0.5` / `0`), whatever made the blocks.
   Where blocks came from is recorded in `meta$metadata$blocks`.
8. **Encodings.** A container stores nucleotide alleles (`A/C/G/T`, or digit
   labels such as PLINK's 1/2) or dosage. `0/1/2` dosage is `"numeric"`: only
   the user knows which allele was counted, so ask before converting.
9. **Missing calls are `NA`** inside a container (`0`, `-`, `N`, `.`, `NA`
   are normalised on read; in dosage data `0` is a genotype). Declare an odd
   token with `missing_chars = "9"` on the read.
10. **Figures take the right object:** `plot_marker_effects(meta, trait =)`,
    `plot_manhattan_variance(res[[trait]]$block_summary)`,
    `plot_simulation_metrics(sim$summary)` (not `sim`).
11. **`load_geno_meta()` is for trusted files only** (it unserialises R
    objects). Containers saved by older versions are upgraded in memory with
    one message; the file is left as it was.

More symptom-to-fix rows: `hs_agent_guide("troubleshooting")`.

## Working habits

- Never read a large genotype file whole to look at it; inspect with
  `data.table::fread(f, nrows = 5)` or `head`.
- Keep inputs read-only, write each run to its own output folder, fix and
  record the seed, and save `sessionInfo()` beside the results.
- haploStack's conditions are classed: errors are `hs_error` (input errors
  `hs_error_input`, a missing optional package `hs_error_dependency` with the
  install command in `$hint`), progress messages are `hs_message`. Silence the
  progress chatter with `hs_quiet(expr)`.
- Optional packages (Suggests) unlock routes: HapSelect (LD blocks),
  rrBLUP / BGLR (effect fitting), genomicSimulation (marker-level crossing
  backend), circlize (circos figures), future + future.apply (`n_cores`).
  A call that needs one stops and names the install command.

## Reporting a problem

If haploStack contradicts its help pages, returns wrong numbers, loses data
silently, or errors on valid input: reproduce it on simulated data
(`sim_dataset()`), note `packageVersion("haploStack")` and `sessionInfo()`,
and tell the user. File an issue at the package's `BugReports` URL
(`utils::packageDescription("haploStack")$BugReports`) only if the user asks.
