# haploStack I/O reference — v0.5.5

A user guide to the **row-ordering contract** introduced in **haploStack 0.5.5**:
one rule that decides which sample sits on which row, in every figure that has a
sample axis.

- `use_sample_order()` — the one new export. Broadcasts an order to every figure.
- `sample_order =` — now on `plot_haplotype_alignment()` (since 0.4.10),
  `plot_localgebv_alignment()`, `plot_favourable_blocks()`,
  `plot_phenotype_bars()` and `plot_top_block_heatmap()`.
- `plot_meta_alignment(check_order =)` — tells you when composed panels disagree.

Companion sheets:
[`haploStack_io_reference_0.5.3.md`](haploStack_io_reference_0.5.3.md),
[`haploStack_io_reference_0.5.2.md`](haploStack_io_reference_0.5.2.md),
[`haploStack_io_reference_0.4.7.md`](haploStack_io_reference_0.4.7.md).
Container spec: [`geno_meta_structure.md`](geno_meta_structure.md).

---

## The problem this solves

`cluster_samples()` has produced a dendrogram leaf order since 0.4.10, and until
now exactly one figure could read it. Everything else picked its own order:

| Figure | Row order before 0.5.5 |
|---|---|
| `plot_haplotype_alignment()` | `sample_order`, else alphabetical |
| `plot_localgebv_alignment()` | family history: offspring, then founders |
| `plot_phenotype_bars()` | `samples` argument order, else alphabetical |
| `plot_top_block_heatmap()` | ggplot2's alphabetical discrete scale |

So two panels of the same population put the same line on different rows. Worse,
`plot_meta_alignment()` **hid** it: it reads the first panel's axis and imposes
that order on the others, so a mismatched pair still rendered — silently
re-ordered.

## The rule

> **A per-call `sample_order` argument wins. Then the broadcast option. Then the
> figure's own default.**

The same shape as `fill_na` (per-call argument → `haploStack.fill_na` →
`HS_FILL_NA_DEFAULT`), so there is one precedence rule in the package, not two.

## `use_sample_order(x)`

```r
cl  <- cluster_samples(meta, k = 4)
old <- use_sample_order(cl)     # broadcast; returns the previous value

plot_haplotype_alignment(meta)          # clustered rows
plot_localgebv_alignment(lg, meta)      # same rows, within each family
plot_phenotype_bars(meta, orientation = "horizontal")

use_sample_order(old)           # restore (old is NULL if nothing was set)
```

| Input | Meaning |
|---|---|
| `hs_sample_clustering` | its `$order`, **and** its `$params` for the figure label |
| character vector | the order, used verbatim |
| `NULL` | clear the broadcast |

Validation happens **on set**, not on use: a duplicated ID, an empty vector or a
non-order object fails at the `use_sample_order()` call you can see, rather than
inside whichever figure you draw next.

The return value is the previous option, invisibly, so a script can be polite:

```r
old <- use_sample_order(cl)
on.exit(use_sample_order(old), add = TRUE)
```

### Per-call overrides

| `sample_order =` | Effect |
|---|---|
| `NULL` (default) | use the broadcast if set, else the figure's own default |
| a character vector / clustering | use it for this call only |
| `TRUE` | require the broadcast — error if none is set |
| `FALSE` | ignore the broadcast for this call |

## Strict vs best-effort — the one asymmetry worth knowing

An order passed **as an argument** must be a permutation of that panel's
samples. A mismatch is an error, because you named those samples for that figure:

```r
plot_haplotype_alignment(meta, sample_order = c("L1", "L2"))
#> Error: 'sample_order' must be a permutation of the panel's samples. Missing 4: L3, L4, ...
```

A **broadcast** order is best-effort. It is intersected with whatever the panel
holds, and panel samples missing from it are appended after the ordered ones with
a message:

```r
use_sample_order(cluster_samples(full_panel))
plot_haplotype_alignment(filter_geno_meta(full_panel, keep_samples = subset))
#> (draws fine — the broadcast is intersected with the 20 samples present)
```

Without that asymmetry, one `use_sample_order()` call would turn every filtered
panel in the session into an error at a distance.

## What an order does **not** do

**It never moves a sample between groups.** In the family-faceted figures
(`plot_localgebv_alignment()`, `plot_favourable_blocks()`) an order re-arranges
rows *inside* the offspring block and *inside* the founder block. The
offspring/founder split and its dashed boundary line stay exactly where they
were — that split is the panel's main readable feature, and it is not the row
order's to move.

Haplotype-level panels keep each sample's copies adjacent and in index order
(`L1.1`, `L1.2`), so re-ordering moves whole samples, never single copies.

Blocks are untouched: the column axis stays in genomic order, which is what makes
the 0.5.0 chromosome sub-panels readable.

## Reading a figure's order back

Every figure that participates carries its provenance:

```r
p <- plot_haplotype_alignment(meta)
attr(p, "sample_order")          # the resolved order, as drawn
attr(p, "sample_order_source")   # "default" | "argument" | "option"
```

and a clustering-derived order labels itself in the subtitle:

```
6 samples · 24 markers · rows: ward.D2 / hamming, k=2
```

So a clustered panel cannot be mistaken for an alphabetical one without
re-running the code.

## Composing panels — `plot_meta_alignment(check_order =)`

```r
fam  <- plot_localgebv_alignment(lg, meta)[["gen_1"]][["F_1"]]
bars <- plot_phenotype_bars(meta, orientation = "horizontal")
plot_meta_alignment(list(fam, bars), rel_widths = c(3, 1))
#> Warning: panel(s) 2 were drawn in a different sample order from the reference
#> panel (sources: option vs default). The reference order is imposed on every
#> panel, so those panels are being re-ordered...
```

The fix is to draw both from the same order — usually by broadcasting one:

```r
use_sample_order(cluster_samples(meta, k = 4))
```

Panels are compared only on the samples they **share**, so a phenotype panel
covering a subset is not a mismatch by itself — only a different *relative* order
is. Panels with no order record are skipped rather than assumed to agree. Pass
`check_order = FALSE` to silence it.

## Notes and edge cases

- **`plot_phenotype_bars(sample_levels =)` still wins.** It is the explicit lock
  used to pin a panel to a companion figure; `sample_order` applies below it.
- **`plot_top_block_heatmap()` has no internal clustering** — its Y axis was
  always ggplot2's alphabetical scale. An order sets the factor levels, reversed
  so the first sample is the **top** row, matching the alignment figures.
- **`plot_clustered_alignment()` consumes a broadcast clustering** when
  `clustering = NULL`, instead of computing a second one that would disagree with
  every other figure in the session. A broadcast that is only a character vector
  cannot drive the dendrogram, so that case still computes one.
- **`count_favourable_blocks(plot = TRUE)`** draws through
  `plot_favourable_blocks()`, so a broadcast reaches it too.
- **The option is session state**, not container state. It is not inherited by
  `future` workers (figures are drawn in the parent), but a script that sets it
  and never clears it changes every later figure in that session.
- **Nothing about default output changed in 0.5.5.** A figure drawn with no
  `sample_order` and no broadcast produces exactly what 0.5.4 produced;
  regression tests pin the y-axis factor levels of both alignment figures.
