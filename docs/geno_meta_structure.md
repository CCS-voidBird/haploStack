# GenoMeta — structure & provenance

> **Verified against haploStack 0.5.5** (source tree, 2026-08-13). Every claim
> below was read off `R/` and, where behaviour was at stake, run against a live
> container. The 2026-07-29 pass added the [`metadata` keys](#metadata-keys)
> section; this pass added [what each transform carries
> forward](#what-each-transform-carries-forward) and covered the 0.5.2–0.5.5
> functions that read or write container fields.

`GenoMeta` is a plain named list (no S3 class) — the single source of truth for
genotype, map, block, pedigree, and effect data in haploStack. Strict validation
lives only in `new_geno_meta()` (construction) and `is_geno_meta()` (lightweight
re-check). Every legitimate change to an element flows through a named function;
direct field edits bypass validation and must be re-checked with `is_geno_meta()`.

## Elements & call-chain provenance

| Element | Structure | Created by | Updated / derived by | Validated by |
|---|---|---|---|---|
| `geno_arrays` | 3-D character array `[n_samples × n_markers × ploidy]` | `split_haplotypes()` via `build_geno_meta()` / `read_*()`; `sim_geno_meta()`; the offspring builders (`.cross_block()`, `update_geno_meta_from_matrix()`) | `filter_geno_meta()` (subset), `recode_geno_meta()` (allele recode) | `new_geno_meta()`, `is_geno_meta()` (3-D; `dim[2]==nrow(map_info)`) |
| `map_info` | data.frame `(Chrom, Marker, Pos)` | `.process_map()` via `build_geno_meta()` / `read_*()` | `filter_geno_meta()` (subset + reorder) | `new_geno_meta()`, `is_geno_meta()` |
| `block_info` | data.frame `(Chrom, Block, Start, End)` — plus the coordinate columns `(Pos, Num_SNPs, Block_length)` when built by `expand_block_info()` / `build_block_location_info()`, plus `(Signal_Marker, Signal_Score, N_Signal)` after a signal selection, see below — or `NULL` | `attach_block_info()` / `expand_block_info()`; `attach_signal_markers()` derives it from `block_marker_map_info` when absent | `filter_geno_meta()` (rows preserved; signal columns pruned), `attach_signal_markers()` (signal columns) | `new_geno_meta()` |
| `block_marker_map_info` | data.frame `(Chrom, Marker, Block)` or `NULL` | `attach_block_info()` | `filter_geno_meta()` (pruned to surviving markers) | `new_geno_meta()`, `is_geno_meta()` (markers ⊆ `map_info`) |
| `sample_info` | FAM data.frame `(Family, Sample, Parent1, Parent2, Sex)` — positional | `new_geno_meta()` (dummy from rownames) · `build_geno_meta(fam=)` · `read_plink_ped()` | `filter_geno_meta()` (row subset), `attach_family_history()` (Family/FID write-back); **reset to a dummy FAM by `attach_block_info()`** — see [carry-forward](#what-each-transform-carries-forward) | `new_geno_meta()`, `is_geno_meta()` (`Sample` col 2 aligned to genotypes) |
| `marker_effects` | named list of per-trait data.frames `(Marker, Allele, Effect)` (list name = trait), or `NULL`; a bare data.frame (optionally with a `Trait` column to split on) is accepted and normalised | `new_geno_meta(marker_effects=)`; direct assign; `fit_genomic_prediction()` (one table per fitted trait, merged in) | `filter_geno_meta()` (pruned to surviving markers); **dropped by `attach_block_info()`** | `new_geno_meta()`, `is_geno_meta()` — **optional, never required** |
| `phenotype` | data.frame `(Sample, <trait...>)` or `NULL` | `attach_phenotype()` / `read_phenotype()` · `new_geno_meta(phenotype=)` · `fit_genomic_prediction()` (GEBV trait) · `predict_gebv()` (appends `<trait>_GEBV`) | `filter_geno_meta()` (pruned to surviving samples); **dropped by `attach_block_info()` and `recode_geno_meta()`** | `new_geno_meta()`, `is_geno_meta()` (`Sample` ⊆ genotypes) |
| `metadata` | named list — five required keys plus `encoding` and two optional sub-objects; full spec in [`metadata` keys](#metadata-keys) | `build_geno_meta()` / `read_*()` / `sim_geno_meta()` | `filter_geno_meta()` (counts), `recode_geno_meta()` (`char_dict`,`encoding`), `attach_family_history()` (`family_history`), `attach_signal_markers()` (`signal_selection`) | `new_geno_meta()` (required keys present — nothing more) |

`marker_effects` also has two **read-only** consumers added in 0.5.2 —
`export_marker_effects()` (one tab-separated `(Marker, Allele, Effect)` file per
trait, trait name in the file name, so the output round-trips through
`read_effect_file()`) and `plot_marker_effects()` (marker-level scatter). Neither
writes back to the container. `predict_gebv()` reads `marker_effects` and writes
to `phenotype`; it is the only function that crosses between the two.

## What each transform carries forward

Every transform except `attach_family_history()` **rebuilds** the container
through `new_geno_meta()`, forwarding the fields it was written to forward. A
field the rebuild call does not name is not preserved — it comes back at the
constructor's default, which is `NULL` for `marker_effects` / `phenotype` and a
freshly dummy-filled FAM for `sample_info`. The result still passes
`is_geno_meta()`, because none of those fields is required.

| Transform | `sample_info` | `marker_effects` | `phenotype` | `metadata` |
|---|---|---|---|---|
| `filter_geno_meta()` | subset | pruned to surviving markers | pruned to surviving samples | counts refreshed |
| `attach_phenotype()` | kept | kept | **set** | kept |
| `attach_signal_markers()` | kept | kept | kept | `signal_selection` written |
| `attach_family_history()` | `Family` written back (edits in place, no rebuild) | kept | kept | `family_history` written |
| `recode_geno_meta()` | kept | kept | **dropped** | `char_dict` + `encoding` rewritten |
| `attach_block_info()` | **reset to dummy** | **dropped** | **dropped** | kept (incl. now-stale `family_history` / `signal_selection`) |
| offspring builders (`.cross_block()`, `update_geno_meta_from_matrix()`, the driver's advance step) | rebuilt from the cross plan, else dummy | kept | **dropped** | counts refreshed; the rest carried from the parent |

Three consequences worth planning around:

* **`attach_block_info()` is destructive beyond blocks.** `.replace_block_in_meta()`
  (`R/format_helper.R:2289`) forwards only `geno_arrays`, `map_info`, the two
  block tables and `metadata`. Re-blocking a container therefore also discards its
  effects and phenotypes and resets the whole FAM — `Family`, `Parent1`,
  `Parent2`, `Sex` all return to `"0"` — while `metadata$family_history` is
  carried forward unchanged, so the pedigree table and the tracing history end up
  **out of sync**, which is exactly what `attach_family_history()` exists to
  prevent. Only the dropped `block_info` signal columns are announced (one
  message); the rest is silent. **Attach blocks first**, then effects, phenotype
  and pedigree — or re-run `attach_phenotype()` / `attach_family_history()` and
  re-assign `marker_effects` afterwards.
  Ordering alone does not save a FAM that came **from the input file**: a PED's
  parent columns are already in `sample_info` when the reader returns, so a later
  `attach_block_info()` erases them and `attach_family_history()` then traces
  nothing. Two working routes, both verified:

  ```r
  # (a) let the reader attach the blocks — sample_info is never rebuilt
  meta <- read_plink_ped("demo.ped", "demo.map", block_file = "demo_blocks.txt")
  meta <- attach_family_history(meta)          # traces the PED pedigree

  # (b) re-block later, but carry the FAM across by hand
  si   <- meta$sample_info
  meta <- attach_block_info(meta, "demo_blocks.txt")
  meta$sample_info <- si                       # hand edit …
  stopifnot(is_geno_meta(meta))                # … so re-validate, per the rule below
  ```
* **`recode_geno_meta()` drops `phenotype`.** It forwards `sample_info` and
  `marker_effects` but not `phenotype` (`R/catalogue_profile.R:222`). Recode
  before attaching phenotypes, or re-attach after.
* **The offspring generations carry the parents' `metadata`.** A simulated
  generation refreshes `n_samples` / `n_markers` only, so `family_history` and
  `signal_selection` still describe the *founder* panel. Read `sample_info` and
  `block_info` for what the offspring container actually holds.

Only `filter_geno_meta()` is lossless across all four fields. When in doubt,
check with `names(m$marker_effects)`, `head(m$phenotype)` and `m$sample_info`
after the call rather than assuming.

## Reading `geno_arrays`

The array holds allele **tokens as character** — `"1"` / `"4"` under
`selectiontools`, `"A"` / `"T"` under `nucleotide` — one layer per copy, with
`dimnames` = (samples, markers, copy labels). It is not stored as integer codes;
`metadata$char_dict` is a label registry, not a decoder (see
[`char_dict`](#char_dict)).

Two seams read it. `geno_meta_views()` is the public one, returning the tidy
`alleles` / `marker_meta` / `block_meta` tables the catalogue layer consumes;
internally `.meta_to_allele_matrix()` rejoins the copies into a 2-D
`sample × marker` matrix of `"A/T"`-style cells. That rejoin applies the
container's **NA-any rule**: a genotype cell is uncalled when *any* copy is
missing, so on a diploid one missing copy shows as a gap on **both** copies at
that marker. This is deliberate and pre-existing; it is why block signatures
placeholder-code missing alleles (`HS_FILL_NA_DEFAULT` = `"-"`, overridable with
`options(haploStack.fill_na = )` or a per-call `fill_na=`) instead of dropping
them — dropping made distinct allele combinations collide, fixed in 0.5.2.

## `sample_info` column spec (positional FAM)

| idx | canonical name | meaning | missing / dummy |
|---|---|---|---|
| 1 | `Family`  | family / FID group label | `"0"` (unknown) |
| 2 | `Sample`  | individual ID (IID); **must equal** `dimnames(geno_arrays)[[1]]` | — (required) |
| 3 | `Parent1` | father IID | `"0"` = founder |
| 4 | `Parent2` | mother IID | `"0"` = founder |
| 5 | `Sex`     | `0` unknown · `1` male · `2` female | `"0"` |

Input fams are read **by position** (headers optional/ignored); the stored table
always carries the canonical names above, so `sample_info$Sample` is stable. A
one-column input is treated as `Sample` (back-compat), with the other columns
dummy-filled.

## `phenotype` (sample-level traits)

Optional data.frame keyed on `Sample` — the sample-side analogue of
`marker_effects` (which is keyed on `Marker`). The first column is `Sample`
(matching `dimnames(geno_arrays)[[1]]` / `sample_info$Sample`); every other
column is a numeric trait. `attach_phenotype()` filters an incoming table to the
**shared** sample set: phenotype rows for unknown samples are dropped (warning),
genotype samples lacking a phenotype row are kept with `NA` traits (left-join
semantics), and rows are reordered to the genotype sample order — all enforced by
`.standardise_phenotype()` inside `new_geno_meta()`. `filter_geno_meta()` prunes
`phenotype` in lockstep with the surviving samples; `attach_block_info()` and
`recode_geno_meta()` **drop** it (see
[carry-forward](#what-each-transform-carries-forward)). The element is optional
and absent in legacy objects, so `is_geno_meta()` only validates it when present.

Two genomic-prediction functions write into it, both via `attach_phenotype()`, so
existing trait columns survive:

* **`fit_genomic_prediction()`** returns `$geno_meta` — a copy carrying the fitted
  `marker_effects` (merged into any already stored) *and* a phenotype holding the
  modelled trait values. An externally supplied phenotype is attached first, so
  its columns persist alongside the GEBV column.
* **`predict_gebv()`** applies effects that already exist (explicit
  `marker_effects=` beats the stored ones) and appends one column named
  `paste0(trait, "_", gebv_suffix)` — `"<trait>_GEBV"` by default. Every sample in
  the container is always scored, so the column aligns one-to-one with the
  phenotype rows; `markers` / `blocks` / `chromosomes` restrict which loci
  contribute, not which samples appear.

`read_phenotype()` loads a delimited `.pheno`/CSV/TSV into the same shape
(auto-detecting the `Sample`/`IID` column and dropping `FID`). The plotting
helpers (`plot_phenotype_distribution()`, `plot_phenotype_correlation()`,
`plot_phenotype_bars()`, `phenotype_summary()`) are deliberately **portable**:
they accept either a `GenoMeta` or any data.frame carrying a `Sample` column, via
the internal `.as_phenotype_table()` seam.

`plot_phenotype_bars()` draws up to two traits as per-sample dodged bars, with a
`"vertical"` mode (samples on X, for a direct call) and a `"horizontal"` mode
(samples on Y, merge-ready). `plot_meta_alignment()` then merges any plots that
share a discrete **sample axis** — e.g. a genotype-level `plot_localgebv_alignment()`
panel (raw sample IDs on Y) plus the horizontal phenotype bars — into one figure:
it reads the sample order from the first (reference) plot, reorders the other
panels' sample axis to match (identical ticks), and pools the legends into a
shared bottom row (via `cowplot`).

Since **0.5.5** that row order is a package-wide contract rather than a
per-figure accident (`R/plot_order.R`). `use_sample_order(x)` broadcasts one order
— an `hs_sample_clustering` from `cluster_samples()`, a character vector of sample
IDs, or `NULL` to clear — through the `haploStack.sample_order` option, and
`plot_phenotype_bars()`, `plot_localgebv_alignment()`, `plot_favourable_blocks()`
and `plot_top_block_heatmap()` each take a `sample_order=` argument. Precedence is
per-call argument > option > the figure's own default; on
`plot_phenotype_bars()`, `sample_levels=` still outranks both, because it is the
explicit lock used to pin the panel to a companion figure. An order passed *to* a
figure must be a permutation of that panel's samples (mismatch = error); a
broadcast is intersected with whatever the panel holds. `plot_meta_alignment(check_order = TRUE)`
warns when composed panels were drawn in different orders instead of silently
imposing the reference panel's. None of this touches the container — the order
lives in an option and in `attr(p, "sample_order")`, never in `sample_info`.

## `metadata` keys

`metadata` is a plain named list. `new_geno_meta()` checks **only that the five
required keys are present** — not their types, not their values, and not that the
counters still describe the array. Everything else rides along unvalidated:
unknown keys survive `filter_geno_meta()`, `recode_geno_meta()`,
`attach_*()` and the `save_geno_meta()` / `load_geno_meta()` RDS round-trip
untouched, because every transform copies the list forward and edits named keys.

| Key | Required | Type | Written by | Later updated by |
|---|---|---|---|---|
| `ploidy` | ✅ | integer scalar | every reader / builder (detected by `split_haplotypes()`, or the `ploidy=` argument; `2L` for DArT/VCF) | — fixed at construction |
| `char_dict` | ✅ | named list **or** named integer vector: allele token → integer code | reader / builder (see the table below) | `recode_geno_meta()` |
| `n_samples` | ✅ | integer scalar | reader / builder | `filter_geno_meta()`; the crossing simulator's offspring builder (`stacking_cross.R`) |
| `n_markers` | ✅ | integer scalar | reader / builder | same as `n_samples` |
| `source_format` | ✅ | character scalar (free-text provenance label) | reader / builder | — |
| `encoding` | filled automatically | character scalar | the builder, else `new_geno_meta()` by `detect_genotype_encoding()` | `recode_geno_meta()` |
| `family_history` | optional | named list of `data.table(FID, IID, Label)`, one per generation | `attach_family_history()` | replaced by the next `attach_family_history()` |
| `signal_selection` | optional | list of 4 (`score`, `n`, `n_blocks`, `date`) | `attach_signal_markers()` | replaced by the next `attach_signal_markers()` |

A fully populated `metadata` — `str()` of a simulated diploid container after
`attach_family_history()` and `attach_signal_markers()`:

```
List of 8
 $ ploidy          : int 2
 $ char_dict       : Named int [1:2] 1 2
  ..- attr(*, "names")= chr [1:2] "1" "4"
 $ n_samples       : int 8
 $ n_markers       : int 12
 $ source_format   : chr "selectiontools"
 $ encoding        : chr "selectiontools"
 $ family_history  :List of 2
  ..$ gen_1:Classes 'data.table' and 'data.frame': 5 obs. of 3 variables:
  .. ..$ FID  : chr [1:5] "F_1" "F_1" "F_1" "F_1" ...
  .. ..$ IID  : chr [1:5] "S3" "S4" "S5" "S1" ...
  .. ..$ Label: chr [1:5] "O" "O" "O" "P" ...
  ..$ gen_2:Classes 'data.table' and 'data.frame': 5 obs. of 3 variables: ...
 $ signal_selection:List of 4
  ..$ score   : chr "final_mcc"
  ..$ n       : int 2
  ..$ n_blocks: int 4
  ..$ date    : chr "2026-07-29"
```

A container straight from a reader carries the first six keys only.

> **Name collision.** The `metadata` *argument* of `transform_genotype()` and
> `transform_effects()` is a different object — the `process_genotype()` result
> list (`$genotype`, `$encode`, `$map`, `$sample`, `$char_dict`) — not this
> `metadata` element, and not a `GenoMeta`. Only `char_dict` means the same thing
> in both.

### `ploidy`

Integer copy count, the third dimension of `geno_arrays`. Consumers that read it
(`stacking_cross.R`, `stacking_main.R`, `stacking_selection.R`) all default to
`2L` when it is missing, so a hand-built container without it still runs — as
diploid, silently. The package is diploid throughout; nothing writes anything
but `2L` today.

### `char_dict`

The allele **label registry**: which allele tokens this object uses and what
integer code each one stands for. It is *not* how the array is stored —
`geno_arrays` holds the tokens themselves as character (`"1"`, `"A"`).

Note what reads it: **nothing inside the package reads a `GenoMeta`'s
`char_dict`.** It is written at construction, rewritten by `recode_geno_meta()`,
and otherwise carried. The consumers of a `char_dict`
(`process_genotype()`, `transform_effects()`, `.detect_minor_alleles()`) take one
as an argument and default to `DEFAULT_ENCODING_DICT`; `transform_genotype()`
takes its own from the `process_genotype()` result. They see this field only if
you hand it to them. Treat it as a record of what the tokens mean, and
pass it explicitly when it matters.

| Producer | Content | Shape |
|---|---|---|
| `build_geno_meta()`, encoding `"selectiontools"` / `"numeric"` | `A=1, C=2, G=3, T=4, -=-1` (`DEFAULT_ENCODING_DICT`) | list |
| `build_geno_meta()`, encoding `"nucleotide"` | `A=1, T=2, C=3, G=4` (`EXTENDED_ALLELE_DICT` order) | list |
| `build_geno_meta()`, encoding not recognised | observed alleles numbered in order of appearance | list |
| `read_dart_genotype()`, `read_vcf_genotype()` | `A=1, B=2` (A/B dosage labels) | list |
| `read_plink_ped()` | observed alleles numbered, e.g. `A=1, G=2` | named integer vector |
| `sim_geno_meta()` | `"1"=1, "4"=2` (ST) or `A=1, T=2` (nucleotide) | named integer vector |
| `recode_geno_meta("nucleotide")` | the observed letters at their `DEFAULT_ENCODING_DICT` codes, e.g. `A=1, T=4` | named integer vector |
| `recode_geno_meta("selectiontools")` | the observed codes mapped to themselves, e.g. `"1"=1, "4"=4` | named integer vector |

Two things to know before you consume a stored dict:

* **The two built-in dictionaries disagree on order.**
  `DEFAULT_ENCODING_DICT` is alphabetical (`A,C,G,T,-` → `1,2,3,4,-1`);
  `EXTENDED_ALLELE_DICT` is `A,T,C,G` → `1,2,3,4`. Which one a nucleotide
  container carries depends on the path it came in through. Read the codes off
  `char_dict`; never assume a position.
* **The list form is the safe one.** `.detect_minor_alleles()` probes
  `char_dict[["na_char"]]` / `char_dict[["na_num"]]`, which returns `NULL` on a
  list but **errors** (`subscript out of bounds`) on a named integer vector. So
  feeding a `sim_geno_meta()` / `read_plink_ped()` / `recode_geno_meta()` dict
  straight into `process_genotype(char_dict=)` fails; wrap it in `as.list()`.
  No builder emits `na_char` / `na_num` today — the missing-value token lives in
  the dict as `"-" = -1L`.

### `n_samples` / `n_markers`

Derived counters, kept as a convenience for reporting (the format-helper GUI
prints them). They are **written, not derived**: `new_geno_meta()` never
recomputes them, so editing `geno_arrays` by hand leaves them stale and nothing
complains — `is_geno_meta()` checks `dim(geno_arrays)[2] == nrow(map_info)`, not
the counters. `filter_geno_meta()` and the crossing simulator refresh them.
The authoritative shape is always `dim(geno_arrays)`.

### `source_format`

Provenance label for where the object came from. No code branches on it — it is
recorded and displayed, nothing more.

| Value | Written by |
|---|---|
| `"user_supplied"` | `build_geno_meta()` default; `read_genotype(preset = "default")` |
| `"plink"` | preset `PL` |
| `"selectiontools"` | preset `ST` / `read_st_genotype()` |
| `"genomic_simulation"` | preset `GS` / `read_gs_genotype()` |
| `"plink_ped"` | `read_plink_ped()` |
| `"dart_seq"` | `read_dart_genotype()` |
| `"vcf"` | `read_vcf_genotype()` |
| `"selectiontools"` / `"nucleotide"` | `sim_geno_meta()` — which stores the *encoding* in this slot |

Any other string reaches the field if a caller passes
`build_geno_meta(source_format = ...)`.

### `encoding`

How allele tokens in `geno_arrays` should be read. Set by the builder when it
knows; otherwise `new_geno_meta()` fills it from
`detect_genotype_encoding()` — and if detection *fails*, the assignment is
`NULL`, so the key is simply **absent** and every consumer re-detects on demand.

Since **0.5.1** there is one vocabulary, `HS_ENCODINGS`
(`R/format_detect.R`): `"selectiontools"`, `"nucleotide"`, `"allele"`. Note the
edges:

* `detect_genotype_encoding()` emits `"selectiontools"`, `"nucleotide"` or
  `"numeric"` — so `"numeric"` (0/1/2 dosage) is a value that legitimately lands
  in `metadata$encoding` even though it is not in `HS_ENCODINGS`.
  `geno_meta_views()` `match.arg()`s against exactly
  `c("selectiontools", "nucleotide", "numeric")`; `"allele"` is the storage-layer
  token used by `translate_alleles()` / `.detect_minor_alleles()`, not a value
  the views layer accepts.
* `"number"` and `"char"` were **removed** in 0.5.1, not aliased. Passing either
  anywhere now errors by name — `encoding "number" was removed in 0.5.1 — use
  "selectiontools" instead` — enforced centrally by `.encode_arg()`. Objects
  saved by ≤ 0.5.0 that carry one of those strings in `metadata$encoding` will be
  rejected by `geno_meta_views()`; re-set the field (or re-read the source) with
  the canonical token.

`recode_geno_meta(to_type = )` is the supported way to change it, because it
rewrites `geno_arrays` and `char_dict` in the same call. It moves between
`"nucleotide"` and `"selectiontools"` **only** (`match.arg`), and it trusts the
stored `encoding` to say where it is starting from — it does not check that the
array actually looks that way. A container recorded as `"numeric"` therefore has
no recode path: `recode_geno_meta(m, "nucleotide")` on 0/1/2 dosages pushes them
through `REVERSED_ENCODING_DICT` and returns `1 → "A"`, `2 → "C"`, `0 → NA`
without a warning. Re-read the source instead.

### `family_history`, `signal_selection`

Two optional sub-objects, each with its own section below —
[`family_history`](#family_history-metadatafamily_history) and
[Signal markers](#signal-markers-block_infosignal_-metadatasignal_selection).
Neither is validated by `is_geno_meta()`.

## Signal markers (`block_info$Signal_*`, `metadata$signal_selection`)

`select_signal_markers()` turns the per-marker scores of `score_marker_signal()`
into a decision — the `n` marker(s) (default 1) that best tag each block's
inheritance — and `attach_signal_markers()` records that decision on the
container's block-level sheet. Three optional columns are added to `block_info`:

| Column | Type | Meaning |
|---|---|---|
| `Signal_Marker` | character | the selected marker IDs, **comma-separated in rank order** (`NA` when the block has none) |
| `Signal_Score`  | numeric   | the rank-1 marker's score — the block's signal strength |
| `N_Signal`      | integer   | how many markers were selected (`< n` for short blocks, `0` for none) |

The packing keeps `block_info` at one row per block, which `new_geno_meta()`
relies on when it re-orders the table and `rbind`s dummy singleton blocks.
**`get_signal_markers()` is the only supported parser** of the packed column —
it returns a tidy `(Block, Marker, signal_rank, Signal_Score)` table (or the
bare marker vector with `long = FALSE`).

The parameters of the selection that produced those columns are recorded
alongside them in `metadata$signal_selection`:

```
List of 4
 $ score   : chr "final_mcc"     # metric ranked on (attr(sel, "score"))
 $ n       : int 2               # markers requested per block (attr(sel, "n"))
 $ n_blocks: int 4               # blocks that got at least one marker
 $ date    : chr "2026-07-29"    # format(Sys.Date()) at attach time
```

`score` falls back to `NA_character_` and `n` to the highest observed
`signal_rank` when the selection table carries no attributes (e.g. a table
rebuilt by hand). The list is a **record of the decision, not a live index** —
`filter_geno_meta()` prunes the `block_info` columns but leaves `n_blocks`
as it was at attach time, so read the block sheet, not this list, for what is
currently selected.

Lifecycle across the container:

| Function | Effect on the selection |
|---|---|
| `filter_geno_meta()` | prunes selections whose marker left the panel, recomputes `N_Signal`, clears `Signal_Score` when the rank-1 marker is gone (messages once) |
| `attach_block_info()` | **drops** the columns with a message — new block definitions invalidate a selection made against the old ones. `metadata$signal_selection` is **not** cleared, so the record outlives the columns it describes; treat a selection record with no `Signal_*` columns as stale |
| `build_block_location_info()` | recomputes coordinates but joins the signal columns back on `Block`, and writes them to `<prefix>_block_info.csv` |
| `.read_block_file()` | reads the columns back when a coordinate block file carries them |
| `save_geno_meta()` / `load_geno_meta()` | round-trip them (RDS) |
| `plot_marker_signal()` | highlights the stored selection by default (`highlight = "auto"`) |

## `family_history` (`metadata$family_history`)

List of per-generation `data.table(FID, IID, Label)` (oldest → newest) produced by
`trace_family_history()`. `attach_family_history()` stores it here **and** writes
each sample's offspring-family FID into `sample_info$Family`, keeping the pedigree
table and the tracing history in sync.

The list is named `gen_1`, `gen_2`, … in birth order. Each generation table holds
one row per individual **per family it belongs to** — so a parent appears in the
family it founded as well as in the family it was born into:

```
    FID    IID  Label
1:  F_1     S3      O      # Label "O" = offspring of this family
2:  F_1     S4      O
3:  F_1     S5      O
4:  F_1     S1      P      # Label "P" = parent of this family
5:  F_1     S2      P
```

| Column | Content |
|---|---|
| `FID`   | family ID, `F_<i>` numbered **within the generation** — `F_1` in `gen_1` and `F_1` in `gen_2` are different families |
| `IID`   | individual ID, matching `sample_info$Sample` when the pedigree is the object's own FAM |
| `Label` | `"P"` (parent) or `"O"` (offspring) — one full parent set + its children per `FID` |

Because `FID` repeats across generations, the write-back into
`sample_info$Family` uses the **generation-qualified** form `"gen_1.F_1"`, taking
each sample's oldest birth family; samples that never appear as an `"O"` (the
founders) keep whatever `Family` they had, typically `"0"`. `attach_family_history()`
warns but still stores the history when no `IID` matches the genotype sample set.

`family_history` is supplementary: `is_geno_meta()` does not validate it, and it
rides unchanged through `filter_geno_meta()` and `recode_geno_meta()` — so after
filtering samples away it still describes the **full** pedigree it was traced
from. `family_history_to_pedigree()` flattens it back to
`Sample`/`Parent1`/`Parent2`.

It also rides through `attach_block_info()`, which **resets `sample_info`** — so
after re-blocking, the history is intact while the FAM it was written back to is
dummy again (`Family`, `Parent1`, `Parent2`, `Sex` all `"0"`). The two are then
out of sync in the one direction `attach_family_history()` exists to close; re-run
it after any `attach_block_info()` call. Simulated offspring generations carry the
parents' history for the same reason (the `metadata` list is copied forward), so
it describes the founders, not the generation in hand.

## Provenance rule

The [element table](#elements--call-chain-provenance) is the authoritative
call-chain for every element. To change an element, call the listed function
rather than editing the field directly; if a field is edited by hand, re-run
`is_geno_meta()` before passing the object to any transform.
`save_geno_meta()` / `load_geno_meta()` round-trip the whole object;
`load_geno_meta()` re-validates with `is_geno_meta()`.

Order matters as much as the call: `is_geno_meta()` cannot tell you that a
transform silently dropped an optional field, so check
[what each transform carries forward](#what-each-transform-carries-forward)
before choosing the sequence of a build script. The safe order is **read → recode
→ blocks → effects → phenotype → pedigree → signal markers**: the two lossy
steps (`recode_geno_meta()`, `attach_block_info()`) come before anything they
would discard, and `filter_geno_meta()` — the only lossless transform — can be
used at any point.
