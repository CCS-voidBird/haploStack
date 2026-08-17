# haploStack I/O reference — v0.5.3

A user guide to the inputs and outputs of the two exports added in
**haploStack 0.5.3**:

- `count_favourable_blocks()` — per-sample table of how many block haplotypes
  carry a favourable local GEBV.
- `plot_favourable_blocks()` — the same judgement as a 1/0 alignment heatmap.

Both read a single `GenoMeta` and never modify it. They need marker effects,
resolved the same way as everywhere in the catalogue layer.

Companion sheets:
[`haploStack_io_reference_0.5.2.md`](haploStack_io_reference_0.5.2.md) (the
subtype-identity metric this one is deliberately separate from) and
[`haploStack_io_reference_0.4.7.md`](haploStack_io_reference_0.4.7.md).
Container spec: [`geno_meta_structure.md`](geno_meta_structure.md).

---

## What "favourable" means here

A block's **local GEBV** is the sum of its member markers' allele effects for one
haplotype copy (`compute_profile_localgebv()`). This release counts, per sample,
how many of those values fall on the side you are selecting for.

That is **not** the same question as 0.5.2's `count_favourable_subtypes()`, which
asks how many *specific catalogue-chosen haplotypes* a line carries. Same word,
two metrics:

| Function | Favourable means | Counting unit |
|---|---|---|
| `count_favourable_blocks()` (0.5.3) | the local GEBV has the preferred **sign** | block copies |
| `count_favourable_subtypes()` (0.5.2) | the copy **is** the catalogue's chosen subtype | matched targets |

They join on `IID` / `Sample` if you want both in one table.

### Two things that change the numbers

**Counting is per haplotype copy, always.** A diploid sample contributes two
copies per block, so:

```
nPositive + nNegative + nZero == nBlock * ploidy
```

`nZero` is not a column — it is that identity's remainder. Genotype-level
positive/negative block counts already exist as `run_haplotype_catalogue()`'s
`block_counts` output; there is no `count_by` argument here, deliberately, so the
two definitions cannot drift apart.

**`center_by_freq` decides what a sign means.** The default is `TRUE`, matching
`run_haplotype_catalogue()`. Leave it there unless you know why you are changing
it: uncentred, whether a block's local GEBV comes out positive mostly reflects
**which allele your effect file called the reference**, so the counts describe the
effect file's convention rather than the germplasm. Centred, a positive local GEBV
means "above the panel mean for this block".

---

## `count_favourable_blocks()`

```r
count_favourable_blocks(geno_meta,
                        effect_df      = NULL,
                        trait          = NULL,
                        prefer         = c("positive", "negative"),
                        samples        = NULL,
                        blocks         = NULL,
                        chromosomes    = NULL,
                        center_by_freq = TRUE,
                        fill_na        = 0,
                        plot           = TRUE,
                        output_file    = NULL)
```

### Inputs

| Argument | What it takes |
|---|---|
| `geno_meta` | A `GenoMeta` with block assignments (`attach_block_info()`). A family history (`attach_family_history()`) is optional for the table, required for the figure. |
| `effect_df` | `(Marker, Allele, Effect)`, read **by position**. `NULL` uses `geno_meta$marker_effects`. |
| `trait` | Which stored trait to use when `effect_df` is `NULL`. Required when several are stored. |
| `prefer` | `"positive"` (default) or `"negative"` — which sign you are selecting for. |
| `samples`, `blocks`, `chromosomes` | Optional subsets, passed to `compute_profile_localgebv()`. |
| `center_by_freq` | Default `TRUE`; see above. |
| `fill_na` | Value for a (block, copy) cell with no effect information. Default `0`. |
| `plot` | Also draw the dummy alignment. Default `TRUE`. |
| `output_file` | Optional path for the table (tab-separated). |

### Output

A list, always with these five elements:

| Element | Contents |
|---|---|
| `table` | The per-sample `data.table` (below). |
| `plots` | `plot_localgebv_alignment()`'s nested generation → family → `ggplot` list, or `NULL`. |
| `dummy` | The wide `Block` × copy matrix of `1`/`0` the figure is drawn from. |
| `prefer` | The preference used. |
| `level` | `"haplotype"` — the counting unit, always. |

`table` columns, in order:

| Column | Meaning |
|---|---|
| `GID` | Generation from the family history (`"gen_1"`, …); `NA` if unavailable. |
| `FID` | Family within that generation (`"F_1"`, …); `NA` if unavailable. |
| `IID` | Sample ID — the key, one row per sample. |
| `Label` | `"P"` parent / `"O"` offspring; `NA` if unavailable. |
| `nBlock` | Blocks evaluated for this sample. |
| `GEBV` | Total local GEBV over those blocks (all copies summed). |
| `nPositive` | Block copies with local GEBV `> 0`. |
| `nNegative` | Block copies with local GEBV `< 0`. |
| `nPreferred` | Whichever of the two `prefer` names. |

Rows are ordered by `GID`, then `FID`, then descending `nPreferred` — so each
family reads as a ranking of its members.

**`nPositive` and `nNegative` are always both reported**, whatever `prefer` is.
One call describes both directions; `prefer` only decides which one
`nPreferred` copies.

### Reading `nBlock`, and when `fill_na` matters

`nBlock` counts blocks that actually hold a value for that sample. Under the
default `fill_na = 0` that is every block in the output, so `nBlock` is the same
for every sample.

The distinction `fill_na` draws is narrower than it looks. A block appears in the
output at all only when **some** copy has effect information in it — a block whose
markers carry no effects is simply absent, not filled. So `fill_na` governs
individual (block, copy) cells: typically a sample whose alleles are missing
across that block's markers.

```r
# "no information" scored as zero effect (default): nBlock is uniform
count_favourable_blocks(meta, trait = "Yield")$table

# "no information" excluded: nBlock varies per sample; GEBV and the
# sign counts are identical either way (an empty cell contributes nothing)
count_favourable_blocks(meta, trait = "Yield", fill_na = NA_real_)$table
```

### Which pedigree row a sample gets

A sample usually appears in several generations — a parent in `gen_2` is normally
an offspring in `gen_1`. The table gives it **one** row, choosing the oldest
generation in which it is an *offspring* (its birth family); a sample that is
never an offspring takes the oldest generation in which it appears as a parent.
This is the same rule `attach_family_history()` uses to write
`sample_info$Family`. Samples absent from the history — and every sample when
there is no history at all — get `GID = FID = Label = NA`.

### Examples

```r
meta <- load_geno_meta("my_geno_meta.rds")
meta <- attach_family_history(meta)

res <- count_favourable_blocks(meta, trait = "Yield")
res$table                      # the ranking
res$plots[["gen_1"]][["F_1"]]  # the 1/0 dummy alignment

# Selecting for a low value instead (e.g. days to maturity)
count_favourable_blocks(meta, trait = "DTM", prefer = "negative")$table

# One region, written straight to file
count_favourable_blocks(meta, trait = "Yield",
                        chromosomes = list(list(chr = "1H")),
                        output_file = "favourable_1H.txt")
```

---

## `plot_favourable_blocks()`

```r
plot_favourable_blocks(x, family = NULL,
                       prefer = c("positive", "negative"), ...)
```

The viz half. A preferred-sign copy is `1`, everything else `0`, drawn **through
`plot_localgebv_alignment()`** rather than as a separate heatmap — so the
chromosome separators, per-family panels, offspring/founder boundary and genomic
block ordering are the same figure machinery as the continuous local-GEBV
alignment. The fill scale is clipped to `c(0, 1)`, putting `0` on the diverging
scale's white midpoint and `1` at full colour — the negative half of that scale
cannot be reached by a dummy, so including it would leave half the legend dead.

| Argument | What it takes |
|---|---|
| `x` | A `count_favourable_blocks()` result (its `dummy` is used) or a wide dummy table directly. |
| `family` | The `GenoMeta` the counts came from, or a `trace_family_history()` list. Required. |
| `prefer` | Legend label only; taken from `x` when `x` is a result list. |
| `...` | Passed through — e.g. `selected_blocks`, `show_block_labels`, `show_chrom_labels`. |

Returns the nested generation → family → `ggplot` list, or **`NULL` with a
warning** when the object carries no family history. That is deliberate: counting
works without a pedigree, only the per-family figure needs one, and since
`plot = TRUE` is `count_favourable_blocks()`'s default, erroring would make the
default call fail on any panel that has not had `attach_family_history()` run on
it — and lose the table with it.

```r
res   <- count_favourable_blocks(meta, trait = "Yield", plot = FALSE)
plots <- plot_favourable_blocks(res, family = meta,
                                show_block_labels = FALSE)
plots[["gen_1"]][["F_1"]]
```

> **A `0` cell means "not favourable", not "no signal".** A copy with a
> wrong-sign local GEBV and a copy with exactly zero both render as `0`. Use
> `plot_localgebv_alignment()` on the continuous matrix when you need to tell them
> apart.

---

## Related change

`plot_localgebv_alignment()` gained a **`fill_name`** argument (default
`"Local GEBV"`) so a caller drawing a derived quantity on that figure can label
the legend honestly. The default reproduces previous output exactly.
