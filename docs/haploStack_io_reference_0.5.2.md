# haploStack I/O reference — v0.5.2

A user guide to the inputs and outputs of the seven exports added in
**haploStack 0.5.2**. They fall into three groups:

**Marker effects — get them out and look at them**

- `export_marker_effects()` — write a GenoMeta's stored marker effects to a table
  / one file per trait.
- `plot_marker_effects()` — genomic scatter of marker effects.

**Block subtypes — what haplotypes are there, and who carries them**

- `count_block_subtypes()` — *enumerate* the subtypes present in one or more
  blocks.
- `count_subtype_occurrence()` — *look up* one named subtype in one named block.
- `count_subtype_carriage()` — carriage of a *set* of named subtypes.
- `count_favourable_subtypes()` — pick the best subtype per block from a
  catalogue, then count its carriage.

**A constant**

- `HS_FILL_NA_DEFAULT` (`"-"`) — the missing-allele placeholder used when a
  haplotype signature is assembled.

Everything here reads a single `GenoMeta` and never modifies it. Only the
marker-effect functions and `count_favourable_subtypes()` need marker effects; the
subtype functions describe the genotype panel itself.

Companion sheet for the 0.4.7 additions (`compute_marker_stats()`,
`write_block_allele_sheet()`):
[`haploStack_io_reference_0.4.7.md`](haploStack_io_reference_0.4.7.md).
Container spec: [`geno_meta_structure.md`](geno_meta_structure.md).

---

## Shared vocabulary

### The subtype signature

A **subtype** is one haplotype copy's allele combination across a block's member
markers, read in canonical genomic order and pasted with no separator — its
**signature**. For a block of four markers a signature is four characters long.

An unobserved allele is written as `fill_na` (default `"-"`) and **keeps its
slot**. This matters: over markers `(M1, M2, M3)` the copies `(A, T, -)` and
`(A, -, T)` are different haplotypes, and would both paste to `"AT"` if the gap
were dropped. They are `"AT-"` and `"A-T"`.

```r
HS_FILL_NA_DEFAULT                  # "-"
options(haploStack.fill_na = "?")   # change it globally
count_block_subtypes(meta, fill_na = "?")   # or per call
```

An empty `fill_na` is refused — a zero-width placeholder is not a placeholder.

> **Missing-data note.** The container uncalls a genotype cell whenever *any*
> copy of it is missing, so on a diploid a single missing copy shows as `fill_na`
> on **both** of that sample's copies at that marker.

### Where signatures come from

Any of these produce strings you can feed back in as targets — they are all the
same string for the same haplotype:

| Source | Column |
|---|---|
| `count_block_subtypes()` | `Subtype` |
| `write_block_allele_sheet()` | `Signature` |
| `compute_block_catalogue(group_by = "signature")` | `signature` |

### `count_by`

Every counting function takes it, with the same meaning:

| Value | Counts | Sums to |
|---|---|---|
| `"haplotype"` (default) | haplotype copies — a heterozygote contributes to two subtypes | `n_samples × ploidy` per block |
| `"genotype"` | distinct samples carrying the subtype on ≥ 1 copy — a homozygote counts once | ≤ `n_samples` |

---

## 1. `export_marker_effects()`

Gets `geno_meta$marker_effects` out of the container and onto disk.

### Call

```r
export_marker_effects(geno_meta, traits = NULL, output_dir = NULL,
                      prefix = "marker_effects", output_file = NULL,
                      annotate = FALSE)
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `geno_meta` | `GenoMeta` | — | Must carry `marker_effects` (e.g. from `fit_genomic_prediction()`). |
| `traits` | character vector | `NULL` | Which traits to export, in the given order. `NULL` = all. |
| `output_dir` | path | `NULL` | Writes **one file per trait** here, created recursively if absent. Mutually exclusive with `output_file`. |
| `prefix` | string | `"marker_effects"` | File-name stem; the sanitised trait name and `.txt` are appended. |
| `output_file` | path | `NULL` | Writes a single named file. Legal only when exactly one trait is selected. |
| `annotate` | logical | `FALSE` | Adds `Chrom`/`Pos`/`Block` to the **returned** table (not to the files). |

### Output

Returns a `data.table` whether or not anything was written, with the paths
written attached as `attr(x, "files")` (`character(0)` when nothing was).

| Column | Type | Description |
|---|---|---|
| `Trait` | character | Trait the table was stored under. |
| `Chrom`, `Pos`, `Block` | character / numeric / character | Locus and block membership. **Only when `annotate = TRUE`**; `NA` for a marker absent from `map_info`. |
| `Marker` | character | Marker ID, as stored. |
| `Allele` | character | Allele state the effect applies to. |
| `Effect` | numeric | Allele effect, as stored. |

**Written files carry exactly `Marker`, `Allele`, `Effect`** — the GenoMeta's own
columns, in its own row order, so a file round-trips through `read_effect_file()`.
The trait is identified by the **file name**, not by a column. There is
deliberately no combined-file mode under `output_dir`: a combined file would need
a `Trait` column and would no longer match what the container stores.

**File naming.** `<prefix>_<trait>.txt`, with the trait sanitised for the
filesystem — a trait called `Yield (t/ha)` becomes `marker_effects_Yield_t_ha.txt`
and a message says so. The returned `Trait` column keeps the original string.

### Example

```r
meta <- load_geno_meta("my_geno_meta.rds")

eff <- export_marker_effects(meta, output_dir = "effects")
#> Marker effects (Yield) written to: effects/marker_effects_Yield.txt
#> Marker effects (Height) written to: effects/marker_effects_Height.txt
attr(eff, "files")

head(read.table("effects/marker_effects_Yield.txt", header = TRUE))
#>   Marker Allele  Effect
#> 1   Mkr1      1  0.1421
#> 2   Mkr1      4 -0.0873

# One trait to an exact path, with locus columns on the returned table
export_marker_effects(meta, traits = "Yield",
                      output_file = "yield.txt", annotate = TRUE)
```

---

## 2. `plot_marker_effects()`

The marker-level counterpart of the catalogue's block-effect scatter. It shares
the catalogue's internal plotting body, so axis, palette and theme match
`plot_manhattan_effect()` exactly.

### Call

```r
plot_marker_effects(geno_meta, trait = NULL,
                    collapse = c("max_abs", "none"),
                    point_size = 1, alpha = 0.8,
                    title = NULL, table_file = NULL)
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `geno_meta` | `GenoMeta` | — | Must carry `marker_effects`. |
| `trait` | string | `NULL` | Which stored trait. `NULL` is fine when exactly one is stored; with several it errors and lists them. |
| `collapse` | `"max_abs"` / `"none"` | `"max_abs"` | How multiple alleles per marker are reduced — see below. |
| `point_size`, `alpha` | numeric | `1`, `0.8` | Point size and transparency. |
| `title` | string | `NULL` | Replaces the default `"Marker Effect — <trait>"`. |
| `table_file` | path | `NULL` | Writes the plotted table tab-separated. |

### The `collapse` argument

Effect tables hold one row per `(Marker, Allele)`, and a genomic-prediction fit
typically stores the **reference allele at 0** and the alternate at the estimated
effect. Plotting every row therefore draws a dense, meaningless line of zeros
beneath the real signal.

- `"max_abs"` (default) keeps one row per marker — the largest effect in absolute
  value, **sign preserved**. Ties resolve on the lexicographically first allele.
- `"none"` plots every stored row.

### Output

A `ggplot` object, with the plotted table attached as `attr(p, "effect_table")`:
columns `Trait`, `Chrom`, `Pos`, `Block`, `Marker`, `Allele`, `Effect` — the same
shape as `export_marker_effects(annotate = TRUE)`.

Markers absent from `map_info` cannot be placed on the genomic axis, so they are
dropped with a warning naming the first few.

### Example

```r
p <- plot_marker_effects(meta, trait = "Yield")
print(p)
head(attr(p, "effect_table"))

ggplot2::ggsave("yield_marker_effects.png", p, width = 9, height = 4)

# Every allele row instead, and keep the backing table
plot_marker_effects(meta, trait = "Yield", collapse = "none",
                    table_file = "yield_marker_effects.txt")
```

---

## 3. `count_block_subtypes()` — enumerate

"What haplotypes exist in these blocks, how common are they, and are they good?"

### Call

```r
count_block_subtypes(geno_meta, blocks = NULL, trait = NULL,
                     count_by = c("haplotype", "genotype"),
                     center_by_freq = TRUE, output_file = NULL,
                     fill_na = NULL)
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `geno_meta` | `GenoMeta` | — | Must have block assignments. |
| `blocks` | character vector | `NULL` | Which blocks. `NULL` = all. Unknown IDs **warn and are skipped**, they do not error. |
| `trait` | string | `NULL` | Attach this trait's block effects. `NULL` = no effect computation at all. |
| `count_by` | see above | `"haplotype"` | |
| `center_by_freq` | logical | `TRUE` | Passed to `compute_block_catalogue()` when `trait` is given. |
| `output_file` | path | `NULL` | Tab-separated write; the table is returned either way. |
| `fill_na` | string | `NULL` | Missing-allele placeholder. |

### Output

A `data.table`, one row per `(Block, Subtype)`, ordered by block then descending
`Count`. **All six columns are always present**, so results of several calls
`rbindlist()` cleanly.

| Column | Type | Description |
|---|---|---|
| `Trait` | character | Trait whose effects were attached; `NA` when `trait` was `NULL`. |
| `Block` | character | Block ID. |
| `Subtype` | character | The signature. |
| `Count` | integer | Occurrences per `count_by`. |
| `Freq` | numeric | `Count` ÷ the block's total; sums to 1 within a block. Under `"genotype"` the denominator is the number of sample–subtype pairs, not the number of samples. |
| `Effect` | numeric | The subtype's block local GEBV, from its representative carrier; `NA` when `trait` was `NULL`. |

> **Cost.** `trait = NULL` does no effect computation. Naming a `trait` calls
> `compute_block_catalogue()`, which has no block filter and so runs a
> **genome-wide** local-GEBV pass even for one block. For repeated interrogation,
> subset with `filter_geno_meta()` first.

### Example

```r
count_block_subtypes(meta, blocks = c("Chr1_B3", "Chr2_B7"), trait = "Yield")
#>     Trait   Block Subtype Count   Freq  Effect
#> 1:  Yield Chr1_B3    4441     6 0.1250 -0.0581
#> 2:  Yield Chr1_B3    4414     5 0.1042 -0.3414

count_block_subtypes(meta, output_file = "block_subtypes.txt")   # whole genome
```

---

## 4. `count_subtype_occurrence()` — look up one

The elemental lookup. One block, one signature, one number.

### Call

```r
count_subtype_occurrence(geno_meta, block, subtype,
                         count_by = c("haplotype", "genotype"),
                         fill_na = NULL, .signatures = NULL)
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `geno_meta` | `GenoMeta` | — | Must have block assignments. |
| `block` | single string | — | One block ID. |
| `subtype` | single string | — | One signature. |
| `count_by` | see above | `"haplotype"` | |
| `fill_na` | string | `NULL` | Missing-allele placeholder. |
| `.signatures` | data.table | `NULL` | A pre-built signature table to reuse across calls. |

### Output

A single non-negative **integer**.

- A subtype that does not segregate returns `0` **silently** — that is a valid
  answer, and makes the function safe to sweep over candidate signatures.
- An unknown *block* returns `0` **with a warning** — that is a caller mistake.
- A signature of the wrong width will legitimately count `0`; take signatures from
  the sources listed under *Shared vocabulary* rather than typing them.

> **Cost.** Only the block's member markers are read, but every call rebuilds the
> GenoMeta's allele view, which dominates (~0.4 s on a 200 × 2000 panel). For more
> than a couple of lookups use `count_subtype_carriage()`, which answers a whole
> set in one pass, or hoist the signature table into `.signatures`.

### Example

```r
tab <- count_block_subtypes(meta, blocks = "Chr1_B3")
count_subtype_occurrence(meta, "Chr1_B3", tab$Subtype[1])
#> [1] 6

count_subtype_occurrence(meta, "Chr1_B3", "ACGT")   # not segregating here
#> [1] 0
```

---

## 5. `count_subtype_carriage()` — a set of targets

"How many of these haplotypes does this line have?" Usually pointed at a
one-sample sub-GenoMeta.

### Call

```r
count_subtype_carriage(geno_meta, subtypes,
                       count_by = c("haplotype", "genotype"),
                       fill_na = NULL, .signatures = NULL)
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `geno_meta` | `GenoMeta` | — | The panel whose carriage is counted. |
| `subtypes` | data.frame or named vector | — | The targets — see below. |
| `count_by` | see above | `"haplotype"` | Under `"genotype"` each sample's `Count` is `0`/`1`. |
| `fill_na`, `.signatures` | | `NULL` | As above. |

**Specifying targets.** Either a data.frame with a `Block` column plus one of
`Subtype` / `signature` / `Signature` — so output from `count_block_subtypes()`,
`write_block_allele_sheet()` or `compute_block_catalogue()` passes straight in —
or a named character vector, `c(Chr1_B3 = "4441", Chr2_B7 = "1114")`. Duplicates
are dropped with a warning; target blocks absent from the GenoMeta warn and
count `0`.

### Output

A named `list`:

| Element | Type | Description |
|---|---|---|
| `n_target` | integer | Number of distinct `(Block, Subtype)` targets. |
| `n_samples` | integer | Samples in the GenoMeta. |
| `n_carried` | integer | Targets carried, pooled over samples. **The headline number** — for a one-sample GenoMeta this is the answer. |
| `total` | integer | Sum of `Count` over every sample and target (matching copies under `"haplotype"`). |
| `by_sample` | data.table | `Sample`, `nTarget`, `nCarried`, `Total`, `Frac` — ordered by descending `nCarried`. `Frac` is `nCarried / nTarget`. |
| `detail` | data.table | `Block`, `Sample`, `Subtype`, `Count`, `Carried` — one row per sample × target. |
| `count_by` | string | The mode used. |

**Uncarried targets stay in `detail` as explicit zeros**, so `detail` always has
exactly `n_target × n_samples` rows and the denominator of any fraction is
visible.

### Example

```r
line <- filter_geno_meta(meta, keep_samples = "S1")
res  <- count_subtype_carriage(line, c(Chr1_B1 = "1414", Chr2_B1 = "ACGT"))

res$n_carried    # 1 of 2 targets
#> [1] 1
res$detail
#>      Block Sample Subtype Count Carried
#> 1: Chr1_B1     S1    1414     1    TRUE
#> 2: Chr2_B1     S1    ACGT     0   FALSE

# A whole panel, ranked
count_subtype_carriage(meta, targets)$by_sample
```

---

## 6. `count_favourable_subtypes()` — the pipeline

Catalogue → best subtype per block → carriage. Pointed at one line it answers
"how many favourable haplotypes does it already have?"; pointed at a panel it
ranks the lines.

### Call

```r
count_favourable_subtypes(geno_meta, catalogue = NULL, trait = NULL,
                          fitness_fn = c("max_effect", "frequent_best",
                                         "expected", "min_effect"),
                          min_count = 1, top_n = NULL,
                          count_by = c("haplotype", "genotype"),
                          center_by_freq = TRUE, fill_na = NULL,
                          output_file = NULL)
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `geno_meta` | `GenoMeta` | — | The panel whose carriage is counted. |
| `catalogue` | see below | `NULL` | The reference. `NULL` computes one from `geno_meta` + `trait`. |
| `trait` | string | `NULL` | Used only when `catalogue` is `NULL`. |
| `fitness_fn` | string or function | `"max_effect"` | Which subtype wins a block. |
| `min_count` | numeric | `1` | Minimum count for `"frequent_best"`. |
| `top_n` | integer | `NULL` | Keep only the strongest `top_n` blocks. `NULL` = all. |
| `count_by` | see above | `"haplotype"` | |
| `center_by_freq` | logical | `TRUE` | Used when the catalogue is computed here. |
| `fill_na` | string | `NULL` | Missing-allele placeholder. |
| `output_file` | path | `NULL` | Writes the **`by_sample`** table. |

**What `catalogue` accepts.** A `compute_block_catalogue(group_by = "signature")`
table; a per-trait named list as `run_haplotype_catalogue()` returns (each element
carrying `$hap_catalogue`); or an `analyze_block_catalogue()` result, whose
`selected` picks are then **honoured verbatim** — that is the escape hatch when
you need `trait_weights`, which the light path sums with equal weight.

**`fitness_fn`** is passed to the same selector `analyze_block_catalogue()` uses,
so *favourable* means one thing across the package:

| Value | Picks |
|---|---|
| `"max_effect"` (default) | highest effect |
| `"frequent_best"` | best among subtypes with `Count >= min_count` |
| `"expected"` | frequency-weighted effect (`Effect × Count / block total`) |
| `"min_effect"` | lowest effect — for counting *un*favourable haplotypes |
| a function | receives a block's subtype rows, must return exactly one |

### Output

The `count_subtype_carriage()` list, plus:

| Element | Type | Description |
|---|---|---|
| `reference` | data.table | `Block`, `Subtype`, `Effect` — the target chosen per block, strongest first (ascending under `"min_effect"`). |
| `fitness_fn` | string | The rule used (`"custom"` for a function). |

`detail` additionally carries each target's `Effect`.

### Example

```r
cat_dt <- compute_block_catalogue(meta, trait = "Yield",
                                  group_by = "signature")

res <- count_favourable_subtypes(meta, cat_dt)
res$reference          # the target chosen per block, strongest first
#>      Block Subtype    Effect
#> 1: Chr1_B2    1441 0.8670705
#> 2: Chr3_B2    1144 0.7058432

# Rank the panel on the 5 strongest blocks, and save the summary
count_favourable_subtypes(meta, cat_dt, top_n = 5,
                          output_file = "favourable_carriage.txt")$by_sample
#>    Sample nTarget nCarried Total  Frac
#> 1:     S3       5        2     3   0.4
#> 2:    S10       5        1     1   0.2

# One line on its own: n_carried is the headline number
line <- filter_geno_meta(meta, keep_samples = "S3")
count_favourable_subtypes(line, cat_dt)$n_carried
```

> **Note.** This counts *identity matches* against catalogue-chosen subtypes. It
> is a different question from counting favourable blocks by **local-GEBV sign**
> per sample (`docs/Update_notes.txt` item 1), which is planned as its own
> function.

---

## Quick reference

| Function | Returns | Default output file | Row unit |
|---|---|---|---|
| `export_marker_effects()` | data.table + `attr(,"files")` | `marker_effects_<trait>.txt` (one per trait) | marker × allele |
| `plot_marker_effects()` | ggplot + `attr(,"effect_table")` | — (`table_file` optional) | marker |
| `count_block_subtypes()` | data.table | caller-named | block × subtype |
| `count_subtype_occurrence()` | integer | — | — |
| `count_subtype_carriage()` | list | — | sample × target |
| `count_favourable_subtypes()` | list | caller-named (`by_sample`) | sample × target |

All written files are tab-separated with a single header row; read them back with
`data.table::fread()`.

> **`fread()` and numeric alleles.** A signature made only of digits (e.g. the
> SelectionTools `1`/`4` encoding) is type-guessed as an integer on re-read, which
> drops any leading placeholder distinction. Pass
> `colClasses = c(Subtype = "character")` (or `Signature`) when that matters.

### Choosing between the subtype functions

| You want | Use |
|---|---|
| everything in a block | `count_block_subtypes()` |
| one specific haplotype | `count_subtype_occurrence()` |
| a hand-picked set of haplotypes | `count_subtype_carriage()` |
| the catalogue's best per block | `count_favourable_subtypes()` |
| the marker-by-marker allele breakdown | `write_block_allele_sheet()` (0.4.7 sheet) |
