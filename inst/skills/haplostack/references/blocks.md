# Haplotype blocks - a block map first, an LD resolver only without one

The recommended input is a **block map you already have**; it needs no LD
parameter. Building LD blocks inside haploStack is a side route, done by
resolvers named after their backend so nobody mistakes one for the other.

| Route | Call | Status |
|---|---|---|
| 1. your block map | `read_plink_ped(block_file =)`, `read_genotype(block_file =)`, `build_geno_meta(blocks =)`, `attach_block_info()` | recommended |
| 2. LD via HapSelect | `hapS_ld_blocks(meta)` | default resolver (needs HapSelect) |
| 3. LD via SelectionTools | `st_ld_blocks(meta)` | being retired; off unless `use_st_formats(TRUE)` |

Every route ends in the same table, `Chrom, Marker, Block`. Markers absent
from the map become single-marker blocks.

## 1. A block map you already have

```r
meta <- read_plink_ped(ped, map, block_file = "blocks.tsv")   # file or data.frame
meta <- build_geno_meta(geno_2d, map, blocks = block_df)
meta <- attach_block_info(meta, "blocks.tsv")                 # later, on a container
```

**Two shapes, defined by their header names** (case does not matter):

- block-marker map: `Chrom, Marker, Block` (blocking by markers);
- block-region map: `Chrom, Block, Start, End` (blocking by position; markers
  inside each region join it, and the declared span is kept as
  `Region_Start` / `Region_End`).

Anything else stops with an `hs_error_block_map`: a headerless file, an alias
header such as `Chr, SNP, Block`, a block on two chromosomes, a map with no
marker in the panel. Rename columns before attaching. Other behaviour:

- Identifiers are read as text (`007` stays `007`).
- A one-to-one chromosome rename (`1H` in the map, `1` in the panel) is
  reconciled and reported once; otherwise pass `chrom_map =`.
- A marker in two blocks is an error unless `overlap = "first"`
  (`block_overlap = "first"` on the readers) keeps its first block.
- Fewer than half the panel's markers matched is a warning: check the marker
  ID scheme (`1_123` vs `chr1:123`).
- Re-blocking keeps the effects, phenotype and FAM. It clears what depended on
  the old blocks: the signal selection, a gene annotation, the block
  provenance.

Where the blocks came from is recorded in `meta$metadata$blocks`: `source`
(`"user_map"`, `"hapS_ld_blocks"`, `"st_ld_blocks"`, `"unknown"`), file name
and md5, `n_blocks`, `n_singletons`, `n_matched`, `n_unmatched` and the rule.

## 2. LD blocks with HapSelect - `hapS_ld_blocks()`

HapSelect is a Suggests (`remotes::install_github("wrshf7/HapSelect")`; it
compiles C++, so Windows needs Rtools).

```r
meta <- hapS_ld_blocks(meta)            # r2 >= 0.5, tolerance 3, tol_reset, flanking
blk  <- hapS_ld_blocks(meta, ld_threshold = 0.6, tolerance = 1L,
                       attach = FALSE, output_dir = OUT, prefix = "run")
```

- **The rule is passed for you.** HapSelect's own defaults are threshold 0.7,
  tolerance 1. If you call `HapSelect::def_blocks()` yourself, pass all three
  settings explicitly.
- **Unplaced markers** (no position, chromosome `0` / `UN`) are left out of the
  LD scan (HapSelect reads row order as genomic order) and become singletons.
- **Cost** grows with the square of the markers per chromosome (pure R): minutes for a few thousand
  markers, much longer for a whole genome. Rehearse on one chromosome
  (`filter_geno_meta(meta, chromosomes = ...)`). `parallelize = TRUE` uses
  HapSelect's `future` parallelism and restores your plan afterwards.
- **Many singletons are normal** at small sample sizes. Report the share.
- **Markers whose genotype never varies** (monomorphic, or every line
  heterozygous - in an inbred panel usually a collapsed paralog) give
  HapSelect's "standard deviation is zero" warning once per pair;
  `hapS_ld_blocks()` muffles exactly that warning and counts the markers.
- **Block IDs** are HapSelect's own, `"<chromosome key>:<n>"` (`"1:1"`).
- `attach = FALSE` returns the block tables and the provenance record;
  attaching that list later keeps the resolver's identity.

## 3. SelectionTools - `st_ld_blocks()` and `st_convert()`

SelectionTools is being retired from haploStack. Everything that touches its
formats is off by default: turn it on with `use_st_formats(TRUE)` for the
session, or `st_formats = TRUE` on the call.

- `st_ld_blocks(meta)` runs SelectionTools' LD blocking in a private temporary
  folder and maps the blocks back. Diploid only; `ld_threshold` in (0, 1];
  default tolerance 0. Its partition differs from HapSelect's at the same
  settings - a tolerance value does not transfer between the engines.
- `st_convert()` reads SelectionTools' own matrix + map files into a
  `GenoMeta`, turns a `st.def.hblocks()` table into a block-marker map, and
  writes a `GenoMeta` in those formats (`to = "selectiontools"`).
- SelectionTools' allele codes (`A = 1 ... T = 4`, missing `-1`) never enter a
  container; they are translated at this boundary.

## After blocking

```r
loc <- build_block_location_info(meta, output_dir = OUT, prefix = "run")
# run_block_info.csv        Chrom, Block, Start, End, Pos, Num_SNPs, Block_length
# run_block_marker_map.csv  Chrom, Marker, Block  (reads back with attach_block_info())
```

`run_haplotype_catalogue(ld_threshold =, tolerance =)` only label output file
names; pass the values the resolver used, or describe your own map's origin in
the run notes.
