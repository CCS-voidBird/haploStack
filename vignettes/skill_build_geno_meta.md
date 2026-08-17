---
name: haplostack-build-geno-meta
description: >
  Build a haploStack GenoMeta object from in-memory variables or genotype files
  (SelectionTools, genomicSimulation, PLINK, or a custom layout), attach optional
  blocks / effects / phenotypes, and write it back out in another format. Use when
  a user needs to get genotype data into haploStack or convert between formats.
---

# Skill — Build and convert a GenoMeta

`GenoMeta` is haploStack's single central object (genotypes + map + blocks +
`sample_info` + optional effects/phenotypes). Every other function consumes one.
Full walkthrough: `vignettes/example_build_geno_meta.Rmd`.

## When to use
- Genotypes are in an R matrix and need to become a `GenoMeta`.
- Genotypes are in a file (ST / GS / PLINK / custom delimited layout).
- You need to convert genotypes from one format to another.

## Build from in-memory variables
`build_geno_meta()` requires **rows = samples, columns = markers** (transpose if
yours is markers × samples).

```r
meta <- build_geno_meta(
  geno_2d       = geno_sm,          # samples × markers; cells like "1/4"
  map           = map,              # data.frame: Chrom, Marker, Pos (bp)
  fam           = fam,              # optional FAM: Family,Sample,Parent1,Parent2,Sex
  blocks        = blocks,           # optional marker→block df / coords / file
  ploidy        = 2L,
  split_char    = "/",              # "/" for "1/4"; "" for "AT"
  encoding      = "selectiontools", # or "nucleotide"
  source_format = "selectiontools")
is_geno_meta(meta)                  # must be TRUE
```

## Build from a file
```r
meta <- read_st_genotype("geno.tsv", map_file = "map.tsv", block_file = "blk.tsv")
meta <- read_plink_ped("data.ped", "data.map")
meta <- read_gs_genotype("geno_gs.csv", map_file = "map.tsv")

# Custom layout:
spec <- resolve_geno_format(list(geno_orient = "marker", split_char = "/",
                                 encoding = "selectiontools",
                                 map_order = c("Marker","Chrom","Pos")))
meta <- read_genotype("geno.tsv", map_file = "map.tsv", preset = spec)
```

## Attach optional contents
```r
meta <- attach_block_info(meta, block_source = blocks)   # or a file path
meta$marker_effects <- effects                           # Marker, Allele, Effect
meta <- attach_phenotype(meta, pheno)                    # Sample + trait columns
```

## Convert / export (format transformation)
```r
write_geno_meta(meta, "out_plink", output_format = "plink_ped")  # .ped + .map
write_geno_meta(meta, "out_gs.csv", output_format = "gs")
mat  <- geno_meta_to_matrix(meta, output_format = "st")$matrix   # inspect, no write
meta_nuc <- recode_geno_meta(meta, to_type = "nucleotide")       # 1→A, 4→T
save_geno_meta(meta, "meta.rds"); meta <- load_geno_meta("meta.rds")
```

## Gotchas
- `build_geno_meta()` orientation is **fixed** (samples × markers); it does not
  auto-transpose.
- `map` columns must be in order **Chrom, Marker, Pos**; positions are bp unless
  you pass `pos_unit = "mbp"`.
- `sample_info` is read **by position** (FAM order), and its `Sample` column must
  match the genotype sample names exactly.
- Markers are auto-intersected between genotype and map; mismatches are dropped
  with a warning.
