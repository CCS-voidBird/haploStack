# haploStack I/O reference — v0.4.7

> **Updated 2026-07-29 for 0.5.2.** `write_block_allele_sheet()` gained a
> `fill_na` argument and **changed its output**: an unobserved allele is now
> written as `"-"` (the `haploStack.fill_na` option) in both `Signature` and
> `Alleles`, where it used to be an empty field. That old behaviour was a defect —
> an empty placeholder has zero width, so `1;1;;4` and `1;;1;4` both pasted to
> `"114"` and distinct subtypes shared one `Signature` string. Signatures are now
> fixed-width and positional. The affected rows below are corrected; nothing else
> in this document changed, and `compute_marker_stats()` is untouched. See
> [`docs/CHANGELOG.md`](CHANGELOG.md) `[0.5.2]` and thread
> [#023](project_record/threads/023-marker-effect-export-subtypes.md).
>
> The companion slide deck `haploStack_io_reference_0.4.7_slide.html` was **not**
> regenerated (no pandoc on the authoring machine) and still shows the empty-field
> behaviour on its `Alleles` / `Signature` slides.
>
> Functions added after 0.4.7 are documented in
> [`haploStack_io_reference_0.5.2.md`](haploStack_io_reference_0.5.2.md).

A user guide to the inputs and outputs of the two functions added in
**haploStack 0.4.7**:

- `compute_marker_stats()` — per-marker QC summary (MAF, heterozygosity, missing
  rate, …).
- `write_block_allele_sheet()` — per-block distinct marker-allele combinations
  (haplotype subtypes) and their frequencies.

Both read a single `GenoMeta` object and can optionally write a tab-separated
(`.txt`) file. Neither needs marker effects — they describe the genotype panel
itself.

---

## The shared input: a `GenoMeta` object

Both functions take `geno_meta`, the package's central container built by
`new_geno_meta()` or any `read_*()` reader and, for the block sheet, completed with
block assignments via `attach_block_info()`. The parts they use:

| Component | What it holds | Used by |
|---|---|---|
| `geno_arrays` | 3-D array `[n_samples × n_markers × ploidy]` of allele codes (diploid → `ploidy = 2`). | both |
| `map_info` | Marker table: `Marker`, `Chrom`, `Pos`. | both (location) |
| `block_marker_map_info` | Marker→block map: `Marker`, `Block`. `NULL` when no blocks are attached. | `compute_marker_stats` (the `Block` column); **required** by `write_block_allele_sheet` |
| `metadata$encoding` | Allele encoding (`"selectiontools"`, `"nucleotide"`, `"numeric"`); auto-detected if absent. | both (missing-token / het handling) |

> **Encoding note.** For dosage / single-copy encodings (`"numeric"`, `ploidy = 1`)
> there is no phase, so heterozygosity is reported as `NA` and each sample
> contributes one haplotype copy.

---

## 1. `compute_marker_stats()`

Per-marker quality-control summary.

### Call

```r
compute_marker_stats(geno_meta, markers = NULL, output_file = NULL)
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `geno_meta` | `GenoMeta` | — | Genotype container (see above). |
| `markers` | character vector | `NULL` | Restrict to these marker IDs. `NULL` = all markers. |
| `output_file` | path | `NULL` | When given, the table is written as a tab-separated file (default file name convention: `marker_stats.txt`). The table is **always** returned regardless. |

### Output

Returns a `data.table` (one row per marker, genomic order). When `output_file` is
set, the same table is written with `\t` separators — column names as below, one
header row.

| Column | Type | Description |
|---|---|---|
| `Marker` | character | Marker ID. |
| `Chrom` | character | Chromosome. |
| `Pos` | numeric | Genomic position. |
| `Block` | character | Assigned block ID, or `NA` if the marker is not in a block (or no blocks attached). |
| `N_obs` | integer | Number of non-missing allele observations = `n_samples × ploidy − missing cells`. |
| `Missing_rate` | numeric | Fraction of missing cells over `n_samples × ploidy` (0–1). |
| `Call_rate` | numeric | `1 − Missing_rate`. |
| `MAF` | numeric | Minor-allele frequency (`1 − major-allele frequency`). `NA` if nothing observed. |
| `Major` | character | Most frequent observed allele (as stored in the encoding). `NA` if none observed. |
| `Minor` | character | Least frequent observed allele. `NA` for a monomorphic marker (only one observed allele) or none observed. |
| `Het` | numeric | Observed heterozygosity — fraction of scored samples whose two haplotype copies differ. `NA` for single-copy (dosage) encodings. |

### Example

```r
meta  <- load_geno_meta("my_geno_meta.rds")
stats <- compute_marker_stats(meta, output_file = "marker_stats.txt")
head(stats)
#>   Marker Chrom      Pos   Block N_obs Missing_rate Call_rate    MAF Major Minor   Het
#> 1   Mkr1     1  7400000 Chr1_B1    16          0.0       1.0 0.3125     1     4 0.375
#> 2   Mkr2     1 13900000 Chr1_B1    16          0.0       1.0 0.4375     1     4 0.375

# Filter to well-called, common markers:
stats[Call_rate >= 0.95 & MAF >= 0.05]
```

---

## 2. `write_block_allele_sheet()`

For every haplotype block, enumerates the **distinct marker-allele combinations**
(subtypes) present in the panel and how common each is. A combination is the tuple
of alleles a haplotype copy carries across the block's member markers, in canonical
genomic order.

### Call

```r
write_block_allele_sheet(geno_meta, output_file = NULL, blocks = NULL,
                         count_by = c("haplotype", "genotype"),
                         fill_na = NULL)
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `geno_meta` | `GenoMeta` | — | Genotype container; **must have block assignments** (`attach_block_info()`). |
| `output_file` | path | `NULL` | When given, the table is written tab-separated (default file name convention: `block_allele_sheet.txt`). Returned either way. |
| `blocks` | character vector | `NULL` | Restrict to these block IDs. `NULL` = all blocks. |
| `count_by` | `"haplotype"` / `"genotype"` | `"haplotype"` | `"haplotype"` counts haplotype copies (Count sums to `n_samples × ploidy` per block). `"genotype"` counts distinct samples carrying the combination on ≥ 1 copy. |
| `fill_na` | single string | `NULL` | Placeholder for an unobserved allele in `Signature` / `Alleles`. `NULL` uses the `haploStack.fill_na` option, itself `"-"` (`HS_FILL_NA_DEFAULT`). **Added 0.5.2.** An empty string is refused — it would shorten `Signature` and let distinct combinations collide. |

### Output

Returns a `data.table` ordered by block then descending `Count`. **One tidy table
for all blocks** — because different blocks have different member markers, the
marker set is stored inside `Markers` / `Alleles` cells (semicolon-separated) rather
than as per-marker columns, so every row has the same shape and the sheet scales to
whole genomes.

| Column | Type | Description |
|---|---|---|
| `Block` | character | Block ID. |
| `Chrom` | character | Chromosome of the block. |
| `Signature` | character | The combination's alleles concatenated in canonical genomic order (no separator) — a compact subtype key. Always exactly as many characters as the block has member markers, with `fill_na` (`"-"`) holding the place of an unobserved allele, so two combinations differing only in *where* the gap falls stay distinct. |
| `Markers` | character | The block's member marker IDs, genomic order, joined by `;`. Identical for every row of the same block. |
| `Alleles` | character | This combination's allele **per member marker**, joined by `;`, index-aligned with `Markers`. Unobserved alleles appear as `fill_na` (`"-"`) — before 0.5.2 they were empty fields. |
| `Count` | integer | Occurrence count (per `count_by`). |
| `Freq` | numeric | `Count` ÷ the block's total count — sums to 1 within each block. |

**Reading a row.** `Markers = "Mkr1;Mkr2;Mkr3"` with `Alleles = "1;4;1"` means this
subtype carries allele `1` at Mkr1, `4` at Mkr2, `1` at Mkr3; `Signature` is `141`.
Split a cell in R with `strsplit(x, ";", fixed = TRUE)`.

A subtype with an unobserved allele reads the same way: `Alleles = "1;-;1"` →
`Signature` `"1-1"`, i.e. allele `1` at Mkr1, nothing called at Mkr2, `1` at Mkr3.
Because the placeholder occupies a position, `"1-1"` and `"11-"` are different
subtypes — under the pre-0.5.2 empty-field behaviour both collapsed to `"11"`.
Note also that the container uncalls a genotype cell when **any** copy is missing,
so on a diploid one missing copy shows as `-` on both of that sample's copies.

### Example

```r
meta  <- attach_block_info(load_geno_meta("my_geno_meta.rds"), block_map)
sheet <- write_block_allele_sheet(meta, output_file = "block_allele_sheet.txt")
head(sheet)
#>     Block Chrom Signature        Markers Alleles Count   Freq
#> 1 Chr1_B1     1       114 Mkr1;Mkr2;Mkr3   1;1;4     4 0.2500
#> 2 Chr1_B1     1       141 Mkr1;Mkr2;Mkr3   1;4;1     3 0.1875
#> 3 Chr1_B1     1       414 Mkr1;Mkr2;Mkr3   4;1;4     3 0.1875

# Most common subtype per block:
sheet[, .SD[which.max(Count)], by = Block]

# Count distinct carriers (samples) instead of haplotype copies:
write_block_allele_sheet(meta, count_by = "genotype")
```

---

## Quick reference

| Function | Default output file | Row unit | Key columns |
|---|---|---|---|
| `compute_marker_stats()` | `marker_stats.txt` | one marker | `MAF`, `Het`, `Missing_rate`, `Call_rate`, `Major`/`Minor` |
| `write_block_allele_sheet()` | `block_allele_sheet.txt` | one block × subtype | `Signature`, `Markers`, `Alleles`, `Count`, `Freq` |

Both files are tab-separated with a single header row and are read back with, e.g.,
`data.table::fread("marker_stats.txt")`.
