---
name: haplostack-family-tracing
description: >
  Attach an external pedigree to a haploStack GenoMeta, trace haplotype-block
  inheritance from parents to offspring, and run in-family IBD / marker-signal
  analysis with visualisation. Use when a user has a pedigree and wants to know
  which blocks/markers report inheritance.
---

# Skill — Family tracing and in-family analysis

Walks a pedigree, assigns each offspring block to a parent, and tests which
markers tag block inheritance. Full walkthrough:
`vignettes/example_family_tracing.Rmd`.

## When to use
- The population has a pedigree (a FAM-style table or any ID + parents file).
- You want per-block parental origin, or which markers reliably tag IBD.

## Load an external pedigree and attach it
```r
ped_in <- read.table("pedigree.txt", header = TRUE, stringsAsFactors = FALSE)
# columns: an individual ID + parent columns; founders use "0"/NA/"" parents.

meta <- attach_family_history(meta, ped_in,
                              parent_cols = c("Father", "Mother"),
                              id_col      = "IID")
meta$metadata$family_history          # per-generation family grouping
```
Trace without a GenoMeta: `trace_family_history(ped_in, parent_cols=..., id_col=...)`.

## Block-level inheritance
```r
ped         <- family_history_to_pedigree(meta)          # Sample, Parent1, Parent2
inheritance <- trace_haplotype_family(meta, ped)
# rows: Sample, Haplotype_idx, Block, Parent (Parent1/Parent2/ambiguous/unknown),
#       n_match, n_total
inheritance[, .N, by = Parent]                           # parental composition
```

## In-family IBD and marker-signal test
```r
ibd <- compute_block_ibd(meta, threshold = 0.5)          # parent×progeny block IBD

sig <- score_marker_signal(meta, mode = "IBD")           # "IBS" weights by similarity
sig$metrics   # per-marker TP/FP/TN/FN, precision, recall, f1, mcc, final_* (MAF-weighted)

plot_marker_signal(sig$metrics, geno_meta = meta, score = "final_mcc")
plot_signal_precision_recall(sig$metrics)                # top-right + high MAF = best tags
```

## Gotchas
- `parent_cols` / `id_col` must name real columns of your pedigree table.
- The confusion test is anchored on **parental** haplotypes vs **every** progeny
  block; `mode = "IBD"` uses the hard call, `mode = "IBS"` the similarity weight.
- `final_<metric> = metric × MAF × (1 − MAF)` down-weights near-monomorphic
  markers — prefer it over the raw metric when ranking candidate tags.
- Real data: read genotypes (`read_*()`), define/attach blocks, attach the
  pedigree, then run `score_marker_signal()`.
