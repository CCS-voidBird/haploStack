# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com), and the project adheres to
semantic versioning.

## [Unreleased]

_The remaining roadmap items (stacking-simulation rebuild #013, genomicSimulation
persistent session GSPS) are re-targeted to **v0.5.x**. An opening PRD for the
remaining input formats (allelic DArT, VCF) exists upstream — questions only,
nothing scoped._

### Documentation (2026-07-29, no code change)

Maintenance pass. The four findings below are documented, not fixed.

- `docs/geno_meta_structure.md` — new **`metadata` keys** section documenting
  every key the container can carry (`ploidy`, `char_dict`, `n_samples`,
  `n_markers`, `source_format`, `encoding`, `family_history`,
  `signal_selection`), each with its type, its writer, what later updates it, and
  a `str()` dump of real content. Verified against the code, not the prose: the
  eight char_dict shapes the builders actually emit (list vs named integer
  vector; `A,C,G,T,-` vs `A,T,C,G` ordering), the `"selectiontools"` /
  `"nucleotide"` / `"numeric"` values that reach `metadata$encoding` versus the
  `HS_ENCODINGS` set and the tokens removed in 0.5.1, the `source_format` labels
  each reader writes, `family_history`'s `FID`/`IID`/`Label` content and its
  `gen_1.F_1` write-back into `sample_info$Family`, and `signal_selection`'s four
  fields. Records four behaviours worth knowing at the boundary: the counters
  are written and never re-derived (so they go stale on a hand edit),
  `signal_selection$n_blocks` is not refreshed by `filter_geno_meta()`, a
  named-vector `char_dict` errors inside `.detect_minor_alleles()`
  (`char_dict[["na_char"]]` → subscript out of bounds) where the list form
  returns `NULL`, and `recode_geno_meta()` trusts the stored `encoding` rather
  than the array — on a `"numeric"` container it silently returns
  `1 → "A"`, `2 → "C"`, `0 → NA`. Documented, not fixed; no code changed.
- Same file: `block_info` row now lists the `(Pos, Num_SNPs, Block_length)`
  coordinate columns added by `expand_block_info()` /
  `build_block_location_info()`.

### Documentation (2026-08-13, no code change)

Second maintenance pass over the container manual, re-verifying
`docs/geno_meta_structure.md` against the **0.5.5** source (it was last checked at
0.5.1, before four releases of new field consumers). The three findings below are
documented, not fixed.

- `docs/geno_meta_structure.md` — new **"What each transform carries forward"**
  section: a per-transform matrix of which optional fields (`sample_info`,
  `marker_effects`, `phenotype`, `metadata`) survive each rebuild. Every
  transform except `attach_family_history()` rebuilds through `new_geno_meta()`,
  and a field the rebuild call does not name comes back at the constructor
  default — `NULL`, or a dummy FAM — while the result still passes
  `is_geno_meta()`, since none of those fields is required. Three findings,
  confirmed by running each call on a live container (`sim_geno_meta()` +
  effects + phenotype + traced pedigree + a signal selection):
  - **`attach_block_info()` is destructive beyond blocks.**
    `.replace_block_in_meta()` (`R/format_helper.R:2289`) forwards only
    `geno_arrays`, `map_info`, the two block tables and `metadata`, so
    re-blocking also **drops `marker_effects` and `phenotype` and resets the
    whole FAM** (`Family`, `Parent1`, `Parent2`, `Sex` → `"0"`) — while
    `metadata$family_history` rides on unchanged, leaving the pedigree table and
    the tracing history out of sync in the one direction
    `attach_family_history()` exists to close. Only the dropped `block_info`
    signal columns are announced.
  - **`recode_geno_meta()` drops `phenotype`.** It forwards `sample_info` and
    `marker_effects` but not `phenotype` (`R/catalogue_profile.R:222`).
  - **`metadata$signal_selection` outlives the columns it describes.**
    `attach_block_info()` drops the `block_info` `Signal_*` columns (with a
    message) but does not clear the selection record, so a container can carry a
    selection record with nothing selected.
  The section closes with the safe build order — read → recode → blocks →
  effects → phenotype → pedigree → signal markers — and notes that
  `filter_geno_meta()` is the only transform lossless across all four fields.
- Same file: coverage of the 0.5.2–0.5.5 functions that read or write container
  fields — `export_marker_effects()` / `plot_marker_effects()` as read-only
  `marker_effects` consumers, **`predict_gebv()`** as a `phenotype` writer
  (appends `paste0(trait, "_", gebv_suffix)`, default `"<trait>_GEBV"`, via
  `attach_phenotype()` so existing traits survive) alongside
  `fit_genomic_prediction()`, and the 0.5.5 row-order contract
  (`use_sample_order()`, the `haploStack.sample_order` option, `sample_order=`,
  `plot_meta_alignment(check_order=)`) with the note that none of it touches the
  container.
- Same file: new **"Reading `geno_arrays`"** section — tokens are stored as
  character, `geno_meta_views()` / `.meta_to_allele_matrix()` are the two read
  seams, and the container's **NA-any rule** (a cell is uncalled when *any* copy
  is missing) is stated where consumers meet it, with the 0.5.2 `fill_na`
  placeholder as its consequence. Also: the offspring builders
  (`.cross_block()`, `update_geno_meta_from_matrix()`, the driver's advance step)
  named as `geno_arrays` creators that carry the parents' `metadata`, and a
  version stamp at the top of the file recording what it was verified against.
- No code, tests, `NAMESPACE` or `DESCRIPTION` changed; the version stays
  **0.5.5**. The package was loaded with `pkgload::load_all()` only to exercise
  the transforms.

### Deliverable (2026-08-13) — `ready_to_use/genoMetaBuild/`, source-only container layer

A second delivery shape beside the tarball and templates, for users whose job is
only **getting data in**: the container layer of 0.5.5, `source()`d rather than
installed. No package code changed — the fork's `R/` files are verbatim copies.

- **`load_genometa.R`** sources 14 files into the session, checks `data.table` /
  `dplyr`, and sets the two `haploStack.*` options that `.onLoad()` would.
- **In scope:** build (`read_plink_ped()`, `build_geno_meta()`, `new_geno_meta()`,
  `sim_geno_meta()`), validate (`is_geno_meta()`, `geno_meta_views()`), filter,
  I/O (`save`/`load_geno_meta()`, `write_geno_meta()`, `read_effect_file()`,
  `export_marker_effects()`, `read_phenotype()`) and the container updates
  (`attach_block_info()`, `attach_phenotype()`, `attach_family_history()`,
  `attach_signal_markers()`, `recode_geno_meta()`, `build_block_location_info()`).
  **Out:** catalogue, local GEBV, plotting, genomic prediction, simulation.
- **`GenoMeta_build_guide.Rmd`** — a runnable new-user guide restricted to
  **PLINK PED/MAP and nucleotide tables**, on simulated demo data
  (`example_data/`, regenerable via `make_example_data.R`).
- **One non-verbatim file, `R/00_engine_shim.R`**: `.HS_ID_COLS` and
  `.lead_identifier_cols` lifted verbatim from `catalogue_engine.R` so the fork
  need not carry the catalogue layer, plus a **fork-specific `.hs_use_rcpp()`
  that always returns FALSE**. The package's version would return TRUE whenever
  an installed haploStack is found — and then fail, because the kernels are
  namespace-private and a *sourced* `split_haplotypes()` cannot see
  `split_genotype_cpp` (reproduced on R 4.4.3 with 0.5.5 installed). Forcing the
  pure-R fallbacks also makes the fork behave the same on every machine
  regardless of what version happens to be installed; the two paths are
  parity-tested in-package, so results are identical.
- Verified on **R 4.4.3 and R 4.2.0**, identical on both: **109** functions
  defined, `codetools::checkUsageEnv()` clean apart from the three now-unreachable
  Rcpp kernel names, all 13 copied files byte-identical (SHA-256) to their
  parents, and the guide knits end to end.

## [0.5.5] - 2026-07-31 — One row order for every figure

Before this release a sample order could reach exactly **one** figure:
`plot_haplotype_alignment(sample_order=)`. Every other sample-axis figure invented
its own row order, so two panels of the same population put the same line on
different rows — including through `plot_meta_alignment()`, whose purpose is to
stack panels that share an axis. 0.5.5 makes the order a **package-wide
contract** (spec `dev/PRD_meta_plot_clustering_order.md`, thread #027, decisions
`dev/questions_sheet.md` Q17).

**No default output changes.** Every figure drawn without a `sample_order` and
without a broadcast produces exactly what 0.5.4 produced; regression tests pin
the y-axis factor levels of both alignment figures.

### Added

- **`use_sample_order(x)`** — broadcast one order to every figure with a sample
  axis, by setting the new **`haploStack.sample_order`** option. Takes an
  `hs_sample_clustering` from `cluster_samples()`, a character vector of sample
  IDs, or `NULL` to clear; returns the previous value invisibly, so
  `on.exit(use_sample_order(old))` restores it. Validated **on set**, so a bad
  order fails at the call the user can see rather than inside the next figure.
- **`sample_order =` on four more figures**: `plot_localgebv_alignment()`,
  `plot_favourable_blocks()`, `plot_phenotype_bars()` and
  `plot_top_block_heatmap()`. Each accepts a character vector, an
  `hs_sample_clustering`, `TRUE` (require the broadcast; error when unset) or
  `FALSE` (ignore the broadcast for this call).
- **`plot_meta_alignment(check_order = TRUE)`** — warns when composed panels were
  drawn in different sample orders. The reference panel's order is imposed on
  every other panel, so without this the disagreement is resolved silently, which
  is the exact failure this contract exists to prevent. Panels carrying no order
  record are skipped rather than assumed to agree.
- **Provenance on the figure.** Plots carry `attr(p, "sample_order")` and
  `attr(p, "sample_order_source")` (`"default"` / `"argument"` / `"option"`), and
  a clustering-derived order adds `· rows: ward.D2 / hamming, k=4` to the
  subtitle. A reader must be able to tell a clustered panel from an alphabetical
  one without re-running the code.
- **`R/plot_order.R`** — the one resolver (`.hs_resolve_sample_order()`) plus
  `.hs_order_within()`, `.hs_warn_order_mismatch()` and the provenance helpers.

### Changed

- **`plot_localgebv_alignment()` orders rows *within* each family group**, never
  across it. The offspring/founder split and its boundary line are not the
  order's to move: they are the panel's main readable feature. Samples missing
  from an order follow the ordered ones, and haplotype copies stay adjacent and
  in index order.
- **`plot_clustered_alignment()` uses a broadcast clustering when
  `clustering = NULL`**, instead of computing a fresh one. Otherwise this panel
  would cluster independently of every other figure in the session and quietly
  disagree with them. A broadcast that is only a character order cannot drive the
  dendrogram, so that case still computes one.
- **`.resolve_sample_order()`** is now a thin wrapper over the shared contract.
  `plot_haplotype_alignment()`'s strict validation and alphabetical default are
  unchanged and pinned by tests.
- **An argument is strict; a broadcast is best-effort.** An order passed *to a
  figure* must be a permutation of that panel's samples — a mismatch is an error,
  because you named those samples for that figure. A *broadcast* order is
  intersected with whatever the panel holds, with unmatched panel samples
  appended. Otherwise setting one order would make every filtered panel in the
  session an error at a distance.
- **`plot_phenotype_bars(sample_levels=)` still outranks `sample_order`** — it is
  the explicit lock used to pin a panel to a companion figure.

### Notes

- Precedence, everywhere: per-call argument > `haploStack.sample_order` > the
  figure's own default. Same shape as `fill_na` (0.5.2), so there is one rule.
- The option is **unset** by default and `.onLoad()` deliberately does not set
  it: there is no sensible package-wide row order, only per-figure defaults.
- Not included (out of scope, `dev/PRD_meta_plot_clustering_order.md` §5):
  `plot_block_counts()` (its axis is a count, and it is produced inside
  `run_haplotype_catalogue()`), block/column-axis ordering (blocks stay in
  genomic order — that is what makes the 0.5.0 chromosome sub-panels readable),
  and tabular output, which keeps its current row order.
- New suite `tests/testthat/test-plot-order.R`.

---

## [0.5.4] - 2026-07-30 — Cross-trait analysis is raw by default

**Breaking (defaults).** Cross-trait analysis rescaled its inputs unless told not
to: `analyze_block_catalogue()` z-scored every trait's effects, and
`create_fitness_index()` divided the trait weights by `sum(abs(weights))`. Both
now use the values **as given**, and rescaling happens only when asked for
(thread #025).

### Changed

- **`analyze_block_catalogue(normalize = )` now defaults to `"none"`** (was
  `"zscore"`). Nothing is rescaled unless requested, so the effects reaching the
  cross-trait selection sum, the reported `S`/`N` tables and the plots are the
  local GEBVs as computed, in each trait's own units and comparable with
  `compute_block_catalogue()`'s `Effect` column. Under the default the `N` table
  is the identity transform of `S`.
  - **`value_source = "normalized"` counts as asking**, so it promotes
    `normalize = "zscore"` (the historical transform) with a message rather than
    returning raw numbers under a normalized label. Callers that were already
    passing `value_source = "normalized"` therefore see **no change**.
  - Naming **both** `value_source = "normalized"` and an explicit
    `normalize = "none"` is a contradiction — `N` would be `S` — and now errors
    naming the fix instead of silently mislabelling.
  - The result carries **`$normalize`**, the transform actually applied
    (`"none"`, a method name, or `"custom"` for a function).
- **`create_fitness_index(normalise_weights = )` now defaults to `FALSE`** (was
  `TRUE`), and likewise on `update_fitness_index()` and
  `run_fitness_pipeline()`. The composite is the plain `sum(weight * effect)`, so
  the index — and any GEBV computed from it — stays in raw units. Previously unit
  weights `c(A = 1, B = 1)` silently **halved** the composite. Rescaling is still
  one argument away, and because it is a positive constant factor it never
  changed any ranking or selection decision — only the magnitudes.
- **`plot_block_effect_genome()` no longer claims a standard score it did not
  compute.** Its subtitle reads "white = average effect (within-trait standard
  score)" only when a transform actually ran; a raw-default result re-plotted with
  `value_source = "normalized"` now reads "white = zero effect (raw)". Results
  from before 0.5.4 carry no `$normalize` field and are treated as rescaled, which
  is what they were.
- **`ready_to_use/` cross-trait templates** (`templates/` and `workshop/`
  `04_cross_trait_analysis.R`) switched from `zscore`/`normalized` to
  `none`/`raw`, with a comment naming when to switch back. They passed the old
  values explicitly, so the default change alone would not have reached them.

### Notes

- **What is *not* touched:** `center_by_freq` (default `TRUE`) still mean-centres
  allele effects by population frequency inside the catalogue. That is a
  within-trait part of the local-GEBV definition — it is what makes the sign of a
  block effect mean "above / below the panel mean" rather than "which allele the
  effect file called the reference" — and it is applied identically to every
  trait, so it is not a cross-trait rescaling. Changing it would move every
  catalogue number in the package away from the origin GRDC workflow.
- Anyone comparing **absolute** cross-trait or fitness-index magnitudes against
  output from 0.5.3 or earlier will see different numbers. Rankings, selections
  and block picks are unchanged except where a genuinely different scale changes
  the summed argmax — which is the point of the flag.

### Packaging — first distribution build since 0.4.9 (#026)

- **`.Rbuildignore` added** (package root). The package had **none**, so
  `R CMD build` shipped everything it found: 0.4.7 and 0.4.8 went out at **5.5 MB
  carrying `dev/` run outputs, `docs/`, the `.scratch/` issue tracker, and the
  compiled Windows objects `src/*.o` + `src/haploStack.dll`**. Shipping those
  makes a "source" tarball non-portable, which is the real defect; the size was
  just the symptom. The ignore list covers compiled artefacts, `dev/`, `docs/`,
  `sim/`, `.scratch/`, `CLAUDE.md`, `.agents/`, and scratch output
  (`Rplots.pdf`, `test_load.txt`, `vignettes/figure`).
- **`ready_to_use/haploStack_0.5.4.tar.gz`** — a portable **source** tarball
  (500 KB), sources only, so one file installs on Windows, Linux and macOS by
  compiling `src/*.cpp` on the target machine. Built with `vignettes = FALSE`
  because pandoc is absent here; the vignette *sources* ship, only rendered HTML
  is missing. Verified by installing it into a throwaway library — the C++
  recompiled, the Rcpp kernels loaded, and 0.5.3/0.5.4 behaviour was checked
  through the installed copy.
- **This closes INDEX discrepancy 14**, open since 2026-07-28: the distributed
  workflow had fallen **seven releases** behind source (newest tarball 0.4.9).
  `ready_to_use/README.md`'s install line now names 0.5.4.

### Verification

`devtools::document()` clean (`NAMESPACE` unchanged at 122) · `devtools::test()`
**1578 pass / 0 fail / 12 skip** (was 1556/0/12 at 0.5.3) ·
`devtools::run_examples()` **exit 0** · the 0.5.4 tarball installs from source
into a clean library and runs (`count_favourable_blocks()`,
`analyze_block_catalogue()` raw-default, `create_fitness_index()` unhalved).

---

## [0.5.3] - 2026-07-30 — Stacking evaluation: favourable blocks by local-GEBV sign

`run_haplotype_catalogue()` has always computed per-sample positive / zero /
negative block counts, but only as a by-product of a **full catalogue run** — it
needs an `output_dir`, writes seven files per trait, and returns the counts in
long `(Sample, Label, Count)` form at **genotype** level. There was no way to ask
the question on its own, at haplotype level, against a pedigree. This release is
that question asked directly (`docs/Update_notes.txt` item 1, thread #024, closing
`dev/qa_suggestions.md` **E4**).

### Added

- **`count_favourable_blocks()`** (`R/stacking_evaluate.R`) — one row per sample
  with the requested schema `GID, FID, IID, Label, nBlock, GEBV, nPositive,
  nNegative, nPreferred`, ordered by generation, family, then descending
  `nPreferred`. Standalone: it needs only a `GenoMeta` and marker effects,
  resolved exactly as elsewhere in the catalogue layer (explicit `effect_df`
  wins, else `marker_effects` picked by `trait`).
  - **Counting is always at haplotype level**, never genotype level, so each
    block contributes `ploidy` copies per sample and
    `nPositive + nNegative + nZero == nBlock * ploidy`. `nZero` is that
    identity's remainder rather than a column. Genotype-level counts remain
    `run_haplotype_catalogue()`'s `block_counts`; there is deliberately no
    `count_by` argument, so the two cannot drift.
  - **`prefer = c("positive", "negative")`** sets which sign is favourable, and
    therefore `nPreferred` — `nPositive` and `nNegative` are reported either
    way, so one call describes both directions.
  - **`center_by_freq = TRUE`** by default, matching `run_haplotype_catalogue()`
    rather than `compute_profile_localgebv()`'s bare `FALSE`. Uncentred, the sign
    of a local GEBV largely reflects which allele the effect file called the
    reference; centred, it means "above / below the panel mean for this block",
    which is the only reading under which counting positive blocks is a breeding
    statement.
  - **`nBlock` is per sample**, counted over cells that hold a value. Under the
    default `fill_na = 0` that is every block; `fill_na = NA_real_` keeps "no
    effect information" distinct from "zero effect", and `nBlock` then varies per
    sample while `GEBV` and the sign counts are unchanged.
  - Returns a list with fixed element names — `table`, `plots`, `dummy`,
    `prefer`, `level` — and writes `table` to `output_file` when given.
- **`plot_favourable_blocks()`** (`R/stacking_evaluate.R`) — the 1/0 dummy
  alignment: a preferred-sign copy is `1`, everything else `0`. It **calls
  `plot_localgebv_alignment()`** on an internally built dummy matrix rather than
  reimplementing a heatmap, so the chromosome separators, per-family panels,
  offspring/founder boundary and genomic block ordering are the same figure
  machinery as the continuous alignment. Accepts either a
  `count_favourable_blocks()` result or a bare dummy table. A `GenoMeta` with no
  family history **warns and returns `NULL`** instead of erroring — counting
  works without a pedigree, only the per-family figure needs one, and `plot =
  TRUE` is the default.

### Changed

- **`plot_localgebv_alignment()` gains `fill_name` (default `"Local GEBV"`)** —
  the one edit to an existing function, purely additive. A legend reading "Local
  GEBV" over a 1/0 dummy field is wrong; the default reproduces previous output
  and a regression test pins the default title.
- `utils::globalVariables()` gains the bare symbols the new file uses (11 column
  names + 3 update-join `i.*` symbols). "No visible binding" findings stay at
  **83**, with **0** from the new code.

### Notes

- **Pedigree columns are resolved to one row per sample.** A sample usually
  appears in several generations (a parent in `gen_2` is an offspring in
  `gen_1`), so the new internal `.sample_generation_map()` follows
  `.sample_family_map()`'s existing convention: the **oldest generation in which
  the sample is an offspring wins** (its birth family), falling back to the
  oldest generation in which it appears as a parent. Samples absent from the
  history — and every sample when there is no history — get
  `GID = FID = Label = NA`.
- **This is a different metric from 0.5.2's `count_favourable_subtypes()`**,
  which counts *identity matches* against catalogue-chosen subtypes. This one
  counts *local-GEBV signs*. The two tables join on `IID` / `Sample`.
- **The +/- classification is `run_haplotype_catalogue()`'s own `fcase()`
  idiom**, reused verbatim rather than refactored into a shared helper: the
  catalogue driver needs a long factor-ordered `(Sample, Label, Count)` table for
  `plot_block_counts()` while this needs wide per-sample counts, and the shared
  part is one three-branch expression. The duplication is deliberate and
  recorded, not accidental.

### Verification

`devtools::document()` clean (`NAMESPACE` 120 → 122) · `devtools::test()`
**1556 pass / 0 fail / 12 skip** (was 1448/0/12 at 0.5.2 — +108 assertions, no
regressions) · `devtools::run_examples()` **exit 0** · "no visible binding"
findings **83 → 83**, none from the new file.

---

## [0.5.2] - 2026-07-29 — Marker-effect export, subtype search, positional signatures

Two read-only ways to interrogate a GenoMeta that were missing. A fit attaches
`marker_effects`, but there was **no way to get them back out to a file** and no
marker-level effect figure (`plot_manhattan_effect()` hardcodes `"Block Effect"`
in its title). And while `compute_block_catalogue(group_by = "signature")`
already counts subtypes, it *requires* effects and has no block filter, so there
was no cheap way to ask "which haplotypes exist in this block, and how common are
they?". Both new files only read a GenoMeta; no existing function changed
behaviour.

### Added

- **`export_marker_effects()`** (`R/marker_effects.R`) — returns
  `geno_meta$marker_effects` as one stacked table and optionally writes **one
  tab-separated file per trait** (`output_dir` + `prefix`, trait name in the file
  name) or a single named file (`output_file`, one trait only; the two arguments
  are mutually exclusive). Written files carry exactly the stored columns
  `Marker, Allele, Effect`, so they round-trip through `read_effect_file()`; the
  trait lives in the file name, not an extra column. `annotate = TRUE` adds
  `Chrom`/`Pos`/`Block` to the *returned* table only. Paths written come back as
  `attr(x, "files")`. Trait names are sanitised for the file name (a trait called
  `"Yield (t/ha)"` no longer produces an invalid path) with a message naming the
  substitution; the returned `Trait` column keeps the original string.
- **`plot_marker_effects()`** (`R/marker_effects.R`) — the marker-level
  counterpart of the catalogue's block-effect scatter, sharing the same internal
  plotting body (`.plot_manhattan()`) so axis, palette and theme match
  `plot_manhattan_effect()` exactly. Because effect tables hold one row per
  `(Marker, Allele)` and most fits store the reference allele at `0`, the default
  `collapse = "max_abs"` keeps one point per marker — the largest effect in
  absolute value, **sign preserved** — instead of drawing a dense meaningless
  line of zeros beneath the signal; `collapse = "none"` plots every stored row.
  Markers absent from `map_info` are dropped with a counted warning. The plotted
  table comes back as `attr(p, "effect_table")`.
- **`count_block_subtypes()`** (`R/block_subtypes.R`) — per-block subtype
  (allele-combination) counts for one, several, or all blocks, with the fixed
  six-column schema `Trait, Block, Subtype, Count, Freq, Effect` (always all
  present, so results of several calls bind cleanly). `count_by = "haplotype"`
  (default) counts haplotype copies; `"genotype"` counts distinct samples
  carrying the subtype on at least one copy. `trait = NULL` (default) leaves
  `Trait`/`Effect` as `NA` and does **no** effect computation at all; naming a
  trait joins the block local GEBV from `compute_block_catalogue()`. Block IDs
  not present warn and are skipped rather than erroring — "no hits" is a valid
  answer for a search helper.

- **A three-layer subtype search** (`R/block_subtypes.R`), for asking whether a
  panel — usually a one-sample sub-GenoMeta — carries particular haplotypes:
  - **`count_subtype_occurrence()`** — the primitive: one block, one signature,
    how many times it occurs. A subtype that does not segregate returns `0`
    rather than erroring; an unknown block warns and returns `0`.
  - **`count_subtype_carriage()`** — the same question for a set of
    `(Block, Subtype)` targets, given as a table (any of `Subtype` / `signature`
    / `Signature`, so `count_block_subtypes()`, `write_block_allele_sheet()` and
    `compute_block_catalogue()` output all pass straight in) or a named vector.
    Returns a list: `n_target`, `n_samples`, `n_carried`, `total`, a per-sample
    `by_sample` table and a full `detail` grid in which uncarried targets remain
    as explicit zeros. One row per sample, so a one-sample sub-GenoMeta gives the
    headline count directly and a panel is ranked without a wrapper loop.
  - **`count_favourable_subtypes()`** — the pipeline: take a catalogue, pick the
    single most favourable subtype per block, count its carriage. The
    winner-per-block rule is `analyze_block_catalogue()`'s own
    `.select_block_subtype()`, so `fitness_fn` accepts the same values
    (`"max_effect"` default, `"frequent_best"`, `"expected"`, `"min_effect"`, or a
    custom function) and *favourable* means one thing across the package.
    `catalogue` accepts a `compute_block_catalogue()` table, a per-trait
    `run_haplotype_catalogue()`-shaped list, or an `analyze_block_catalogue()`
    result (whose `selected` picks are then honoured verbatim); `NULL` computes
    one from `geno_meta` + `trait`. `top_n` restricts to the strongest blocks.

### Changed

- **`.compute_block_signatures()` gained `blocks=`**, restricting the melt to
  those blocks' member markers. A single-block signature build on a
  200-sample x 2000-marker panel drops from 1.25 s to 0.42 s, and
  `count_block_subtypes()` uses it whenever no `trait` is given (with a `trait`
  the full table is still built, because `compute_block_catalogue()` labels every
  block). The absent-block warning now checks the GenoMeta's block map rather
  than the signature table, so it still fires when *every* requested block is
  missing and the restricted table comes back empty.

### Fixed

- **Block haplotype signatures are now positional — distinct allele combinations
  no longer collide.** An unobserved allele is written as a placeholder
  (`fill_na`, default `"-"`, the canonical missing character of
  `DEFAULT_ENCODING_DICT`) instead of being dropped, so every signature is
  exactly as wide as its block and each character keeps its marker's position.
  Previously, over markers `(M1, M2, M3)`, the copies `(A, T, -)` and
  `(A, -, T)` both pasted to `"AT"` and were merged into one subtype.
  - `.compute_block_signatures()` fills rather than drops, and **no longer drops
    the copy itself**: a haplotype copy with nothing observed in a block now
    appears as an all-`fill_na` subtype, so `sum(Count)` within a block is again
    `n_samples * ploidy`.
  - `write_block_allele_sheet()` uses the same placeholder in `Signature` and
    `Alleles`. It had the **same string collision** — `1;1;;4` also pastes to
    `"114"` — although its *counts* were always correct, because it groups by the
    per-marker columns rather than by the pasted string. Both producers now emit
    identical strings, so `Signature` and `count_block_subtypes()`'s `Subtype`
    join directly.
  - **Affects existing output:** signature/`Subtype` strings from
    `compute_block_catalogue(group_by = "signature")`,
    `analyze_block_catalogue()` and `write_block_allele_sheet()` change wherever
    a block carries missing data (clean data is unchanged), and blocks with
    missing data gain the previously-dropped copies. Note the container's NA-any
    rule: a genotype cell is uncalled when *any* copy is missing, so on a diploid
    one missing copy shows as `fill_na` on both copies at that marker.
  - Configurable globally with `options(haploStack.fill_na = )` or per call via
    the new `fill_na` argument on `write_block_allele_sheet()` and
    `count_block_subtypes()`. An empty placeholder is refused — it would
    reintroduce the collision.
  - New export: **`HS_FILL_NA_DEFAULT`** (`"-"`).
- **`R CMD check` "no visible binding" findings: 106 → 83** (measured with
  `tools:::.check_code_usage_in_package()` against the 0.5.1 tree and this one).
  Thirteen symbols added to `utils::globalVariables()` (`R/zzz.R`):
  - `"Freq"` and `"Subtype"` retire **5 pre-existing** findings —
    `write_block_allele_sheet()`, `.select_block_subtype()`,
    `analyze_block_catalogue()` (×2) and `plot_block_effect_manhattan()`.
    `"signature"` is registered alongside them for consistency; it produced no
    finding on its own (it is always `:=`-assigned before it is read).
  - `"i.Block"`, `"i.Chrom"`, `"i.Effect"`, `"i.Pos"` — the data.table
    update-join RHS symbols (`x[i, on =, col := i.col]`), supplied by the join
    and never lexically bound. Registering them retires **22** findings across
    `.hap_effects_rcpp()`, `.marker_signal_table()` and the catalogue joins.
  - `"S"` and `"Vbase"` (the block-subtype selection value columns) retire **4**
    — two pre-existing in `analyze_block_catalogue()` / `.select_block_subtype()`
    and the two `.resolve_top_subtypes()` would otherwise have added.
  - `"Carried"`, `"Copies"`, `"Frac"`, `"Target"`, `"Total"`, `"nCarried"`,
    `"nTarget"` for the new carriage tables.
  - Net effect: the new code contributes **0** findings of its own, and the
    pre-existing backlog shrinks by 23.

### Notes

- `count_block_subtypes()` derives both its counts and its effect-join key from a
  **single** `.compute_block_signatures()` call, which is then injected into
  `compute_block_catalogue(.signatures = )`, so the join key is identical by
  construction rather than by convention. This mattered more before the
  positional-signature fix above, when the package's two signature producers
  disagreed on missing data; they now agree, but the single-source construction
  is kept because it is the property that guarantees the join.
- `count_block_subtypes(trait = )` runs a **genome-wide** local-GEBV pass even for
  a single block, because `compute_block_catalogue()` has no block filter.
  Mitigations are named in the docs (`filter_geno_meta()` first, or reuse
  `.signatures=`). Left as-is rather than adding a filter to an existing
  hot-path function.
- The signature collision was first recorded here as INDEX discrepancy 17 and
  deliberately left unfixed; it was then **fixed in the same release** on a
  follow-up request (see *Fixed*, above). Discrepancy 17 is closed.
- `count_subtype_carriage()` answers for a whole target set in one keyed join
  rather than by looping `count_subtype_occurrence()`, because the irreducible
  per-call cost is `geno_meta_views()`. On a 200-sample x 2000-marker panel the
  full 400-block pipeline runs in **0.94 s** against roughly **2.7 min** for the
  loop. The primitive still defines the contract, and a test pins the two to
  agree target by target.
- Exports 113 → **120**. Verification: `document()` clean, `test()`
  **1448 pass / 0 fail / 12 skip** (from 1223/0/12 — +225 assertions, no
  regressions), `run_examples()` exit 0. `devtools::check()` cannot complete on
  this machine — pandoc is absent, so the vignette rebuild fails during the
  source-tarball step, before checking starts. The globals figures above were
  therefore taken with `tools:::.check_code_usage_in_package()`, which is the
  same pass that emits the NOTE.

## [0.5.1] - 2026-07-29 — One encoding vocabulary

The format module had **two names for one concept**. `"char"` / `"number"` was
the live vocabulary of the storage layer (`translate_alleles()`,
`.resolve_preset()`, `transform_effects()`); `"nucleotide"` / `"selectiontools"`
was the read-side and metadata vocabulary (`detect_genotype_encoding()`,
`GENO_FORMAT_PRESETS`, `metadata$encoding`). They met in `process_genotype()`,
which overwrites the preset's `encode` with the *detected* one — so an `ST`
preset produced `"selectiontools"` and handed it to `translate_alleles()`, which
only knew `"number"`. One vocabulary now.

### Removed

- **The `"number"` and `"char"` encoding tokens.** They are removed, not
  aliased: passing either now raises an error naming its replacement
  (`encoding "number" was removed in 0.5.1 — use "selectiontools" instead`).
  Silent normalisation was rejected as an option — it would have left the
  vocabulary split in fact while claiming it was unified. **Breaking:**
  `translate_alleles(arr, "char", "number")` and
  `build_geno_meta(encoding = "number")` no longer work.
- **`.resolve_geno_spec()`** — dead code (nothing called it) carrying a third
  copy of the retired vocabulary. Live counterpart: `resolve_geno_format()` /
  `GENO_FORMAT_PRESETS`.

### Added

- **`HS_ENCODINGS`** (`R/format_detect.R`) — the canonical vocabulary
  (`"selectiontools"`, `"nucleotide"`, `"allele"`), defined in one place.
- **`.encode_arg()`** — validates every encoding argument and always refuses a
  removed token by name, whatever the caller's `allowed` set. `strict = FALSE`
  serves `build_geno_meta()`, which legitimately accepts detector output it has
  no dictionary for.

### Fixed

- **`transform_genotype()`'s round-trip example, broken since ≤ 0.4.9**
  (project-record discrepancy 12). `devtools::run_examples()` now completes
  **exit 0** across the package; the ST → GS → ST round-trip is a true inverse.
  The discrepancy note had guessed this was "a real encoding-resolution defect
  behind the example" — it was exactly that.
- **`build_geno_meta()` could build a GenoMeta nothing could use**
  (discrepancy 15, F1): `"number"` was accepted for `char_dict` selection, then
  stored verbatim in `metadata$encoding`, which `geno_meta_views()` rejects. The
  object built cleanly and failed on first use.
- **`dev/alignment_check/run_alignment_check.R`** (discrepancy 15, F2) — stale
  and non-running: it passed a markers × samples matrix where `build_geno_meta()`
  fixes samples × markers, and used the retired tokens. Both corrected; its
  two-encoding rendering coverage is reachable again.

### Tests

- `tests/testthat/test-format.R` — three `translate_alleles` tests renamed to
  the canonical vocabulary; two added, covering the migration errors at every
  layer that used to accept the old tokens, and the contract that
  `build_geno_meta()` only ever stores an encoding `geno_meta_views()` accepts.

## [0.5.0] - 2026-07-29 — Chromosomes as alignment sub-panels

The haplotype alignment showed chromosomes as text and nothing else. A
whole-genome `plot_haplotype_alignment()` was one undifferentiated strip of
tiles with `Chr 1H,2H,…` in the subtitle, so nothing in the figure said where 2H
began — while blocks, the finer grouping, had both boundary lines and labels.
The two column-level groupings are now drawn at their own scales: chromosomes as
sub-panels, blocks annotated inside them.

### Changed

- **`plot_haplotype_alignment()` (`R/catalogue_viz.R`) draws one sub-panel per
  chromosome** — `facet_grid(cols = vars(chromosome), scales = "free_x",
  space = "free_x")`, so every panel is as wide as its share of the markers and
  is named by the strip above it. `x_pos` remains the global marker index, so
  `p$data` keeps its 0.4.9 contract (plus a new `chromosome` column) and the
  parity tests are untouched. New `chrom_panels` argument (default `TRUE`);
  `chrom_panels = FALSE` restores the single continuous marker axis and the
  `Chr …` subtitle, and a single-chromosome panel is never faceted regardless.
- **Block labels moved inside the panel, into reserved headroom.** They kept
  their vertical style but now anchor at `hjust = 1`, so the text grows
  *downward* from the panel top instead of upward past it, and
  `scale_y_discrete()` reserves `max(2, ceiling(0.1 * rows))` blank row-units for
  them. This both (a) keeps them clear of the chromosome strips and (b) fixes a
  standing defect: at `hjust = 0` the labels grew into a 10 pt `plot.margin` and
  were clipped by the device (`plot_marker_signal()` reserves 46 pt for the same
  construct).
- **Block boundaries and labels are layer data, not `annotate()`/bare
  `xintercept`.** `.alignment_block_boundaries()` and
  `.alignment_block_labels()` now return the owning `chromosome` (and, for
  boundaries, a `cross` flag) so ggplot routes each annotation to its own panel —
  `annotate()` replicates into *every* panel under faceting. Boundaries that
  straddle a chromosome break are dropped in panel mode, where the panel gap
  already states them, and kept in single-panel mode, which therefore keeps its
  exact pre-0.5.0 boundary set.
- **Natural chromosome order** via the new internal `.natural_chrom_levels()`:
  chromosomes whose names start with digits sort numerically first
  (`1H < 2H < 10H`), the rest alphabetically after. The previous
  `order(chromosome, position)` put `"10H"` ahead of `"2H"`, disagreeing with
  `.create_genomic_positions()` and the Manhattan plots. Identical to the old
  order for `1H`–`7H` and `1`–`9`.
- New arguments, all appended so positional calls are unaffected:
  `chrom_panels` (`TRUE`), `show_chrom_labels` (`TRUE`), `block_label_size`
  (`2.5`). The X scale expands additively by half a tile so edge tiles are no
  longer cut in half.

### Fixed

- **`plot_clustered_alignment()` dendrogram pairing.** The reserved block-label
  headroom shifts the alignment rows, which would have slid the dendrogram
  leaves out of register with the rows they index. `plot_haplotype_alignment()`
  now publishes its additive Y expansion as `attr(p, "row_expand")`,
  `plot_cluster_dendrogram()` gains a `leaf_expand` argument, and
  `plot_clustered_alignment()` mirrors one onto the other.

### Tests

- `tests/testthat/test-catalogue-viz.R` — six new tests: natural chromosome
  ordering; one panel per chromosome in natural order with every tile carrying
  its chromosome; `chrom_panels = FALSE` back to one panel + `Chr …` subtitle;
  single-chromosome never faceted; boundary/label layers routed per panel with
  the three cross-chromosome boundaries dropped; headroom reserved only when
  there are labels to place.
- `tests/testthat/test-catalogue-cluster.R` — the dendrogram leaf axis matches
  the alignment's expanded row range, and the unmirrored default does not.

## [0.4.11] - 2026-07-28 — Signal-marker selection

Scoring now ends in a decision. `score_marker_signal()` ranked every marker but
never picked one, so "which marker tags this block?" was re-derived by hand at
every use site — and the marker-signal plot invented its own top-quantile
candidate set at draw time, highlighting a different set than any downstream
analysis. The selection is now a first-class result that is recorded on the
`GenoMeta` and simply rendered by the plot.

### Added

- `select_signal_markers()` (`catalogue_similarity.R`): ranks the markers of
  each haplotype block by a signal metric (default `final_mcc`) and keeps the
  top `n` (default 1) as that block's **signal markers**. Ties break
  deterministically (`ties = "position"` by default, or `"maf"` / `"first"`);
  `min_score` leaves a block unselected rather than giving it a weak
  representative. Returns `(Block, Chrom, Marker, Pos, score_val, signal_rank)`.
- `attach_signal_markers()`: records a selection on `block_info` as
  `Signal_Marker` (comma-separated in rank order), `Signal_Score` (the rank-1
  marker's score) and `N_Signal`, plus the selection parameters in
  `metadata$signal_selection`. Derives `block_info` from
  `block_marker_map_info` when the container has none.
- `get_signal_markers()`: the single supported parser of the packed column —
  a tidy `(Block, Marker, signal_rank, Signal_Score)` table, or the bare marker
  vector with `long = FALSE`.
- `save_marker_signal()`: writes a marker-signal plot and its marker × block
  metrics table under one stem (`marker_signal.png` + `marker_signal.csv`).
- `score_marker_signal()` gained `n_signal` (default 1) / `select_score` and a
  fourth return element, `$signal`, carrying the selection.
- `plot_marker_signal()` gained `highlight`, `signal_markers`,
  `show_marker_labels`, `show_title` and `table_file`; every plot now carries
  its marker × block table as `attr(p, "signal_table")` (marker metrics joined
  to block, with `is_signal`, `block_score`, `block_n_markers`,
  `block_highlighted`).
- `tests/testthat/test-similarity.R`: 20 new cases covering selection, the
  container round-trip, the lifecycle guards below, and the plot defaults.

### Changed

- **`plot_marker_signal()` default output changed.** It now highlights the
  selected signal markers (`highlight = "auto"`: explicit `signal_markers` →
  the selection stored on `geno_meta` → the old quantile rule with a message),
  draws **no marker labels** (`show_marker_labels = FALSE`) and **no title**
  (`show_title = FALSE`) so the vertical block labels stay legible. Block
  labels, the block strip, the quantile line and the theoretical-maximum line
  are unchanged. *Migration:* `show_marker_labels = TRUE, show_title = TRUE`
  restores the 0.4.10 look; `highlight = "quantile"` restores the old
  highlighting rule. Supplying `title` still implies a title.
- `filter_geno_meta()` prunes signal markers that leave the panel, recomputing
  `N_Signal` and clearing a now-meaningless `Signal_Score` (messages once).
- `attach_block_info()` reports that it dropped a stored selection — redefining
  blocks invalidates a selection made against the old ones.
- `build_block_location_info()` joins the signal columns back on `Block` after
  recomputing coordinates, and writes them to `<prefix>_block_info.csv`;
  `.read_block_file()` reads them back when present.
- `docs/geno_meta_structure.md`: `block_info` / `metadata` rows updated and a
  new "Signal markers" section documents the columns and their lifecycle.
- `vignettes/ibd_ibs_signal.Rmd` and `vignettes/example_family_tracing.Rmd` show
  the score → select → attach → plot flow; `dev/r2_functional_test/run_functional_test.R`
  replaces its hand-rolled top-5% candidate CSV with the selection API.

### Fixed

- `combine_haplotypes()`'s `@examples` built a 1 × 2 matrix with only one column
  name (`dimnames = list("Mkr1", "S1")`), so `R CMD check` failed at
  *"length of 'dimnames' [2] not equal to array extent"*. Now
  `list("Mkr1", c("S1","S2"))`. Pre-existing (identical in the 0.4.9 tarball),
  found while checking this release. **`R CMD check` still errors on the next
  example in the chain** — `transform_genotype()`'s round-trip example hits
  *"unsupported encoding 'selectiontools' → 'char'"*, also pre-existing and
  unrelated to this release; see INDEX discrepancy 12.

Suite after the change: 22 files, 379 tests, 0 failures, 12 skips.

## [0.4.10] - 2026-07-27 — Sample clustering

> **Reconstructed from the code on 2026-07-28** — 0.4.10 was set in
> `DESCRIPTION` but shipped with no changelog entry and no tarball (newest build
> remains `haploStack_0.4.9.tar.gz`). Evidence: `R/catalogue_cluster.R` (new,
> mtime 2026-07-27) and the `catalogue_viz.R` edit of the same date; the
> `NAMESPACE` delta against the 0.4.9 tarball — `+export(cluster_samples)`,
> `+export(plot_cluster_dendrogram)`, `+export(plot_clustered_alignment)`,
> `+S3method(print, hs_sample_clustering)`; and `tests/testthat/test-catalogue-cluster.R`
> (11 tests) with a new `_snaps/` directory.

### Added

- `R/catalogue_cluster.R` — sample clustering from haplotype-level allele calls.
  Distances are computed on the haplotype copies but the clustering **leaves are
  whole samples**: the ploidy copies are collapsed to one genotype-level
  distance (`copy_merge = "average" | "pool" | "best"`). Marker-level Hamming
  distance with pairwise-missing exclusion, computed via a single BLAS
  `tcrossprod` over one-hot allele indicators rather than an R pair loop, so it
  scales to ~10K markers with no new dependencies.
  - `cluster_samples()` → `hs_sample_clustering` (with `samples`/`markers`/
    `blocks`/`chromosomes` pre-filtering, `linkage`, optional `k` cut).
  - `plot_cluster_dendrogram()` — dendrogram panel, horizontal or vertical,
    optionally coloured by cluster.
  - `plot_clustered_alignment()` — clustering → row-ordered alignment heatmap
    (+ dendrogram).
  - `print.hs_sample_clustering()`.

### Changed

- `plot_haplotype_alignment()` gained `sample_order`, accepting either a
  character permutation of the panel's samples or an `hs_sample_clustering`
  (its leaf order). Validated by `.resolve_sample_order()`, which errors on
  unknown/missing/duplicated IDs so a stale order cannot silently drop rows. In
  `"haplotype"` mode the ploidy copies stay grouped within each sample.

## [0.4.9] - 2026-07-22 — Alignment plot performance

### Changed

- `plot_haplotype_alignment()` (`catalogue_viz.R`): rebuilt the long-format
  data-prep so it scales to large panels. The `method = "genotype"` path no
  longer melts the full 3-D `geno_arrays` and then runs an N×M
  `by = .(Sample, marker)` group-by of `sort()` + `paste0()`; it now collapses
  the two haplotype copies to a sorted, separator-less genotype matrix in a
  single pass. Both methods replace `as.data.frame.table()` with a vectorised
  `rep()`/`as.character()` build. Plot output is **byte-identical** to 0.4.8.
  Benchmark (800 samples × 12 000 markers): genotype collapse **155.7 s → 2.9 s
  (~53×)**; haplotype long-format build ~2.4× faster.

### Added

- `collapse_genotype_sorted_cpp` (new Rcpp kernel in `src/format_cells.cpp`,
  Kernel A): sorts each cell's non-missing ploidy alleles and concatenates them
  with no separator — the compiled twin of `paste0(sort(allele), collapse = "")`.
  Wired through the existing `.hs_use_rcpp()` toggle with a vectorised pure-R
  fallback (`.collapse_genotype_sorted()` / `.collapse_genotype_sorted_r()` in
  `R/format_split.R`), matching the established kernel-dispatch pattern.
- `tests/testthat/test-catalogue-viz.R` (new): plot-data parity of both row
  methods against the original melt+group-by logic, under both engines.
  `test-format-cells.R` extended with sorted-collapse kernel-vs-pure-R parity
  (NA-drop, both-NA → `""`, multi-character sort order) and input-validation
  checks.
- Source tarball `haploStack_0.4.9.tar.gz` in `ready_to_use/`.

## [0.4.8] - 2026-07-21 — Catalogue output column order

> Reconstructed 2026-07-22 from the `ready_to_use/haploStack_0.4.7.tar.gz` vs
> `_0.4.8.tar.gz` source diff (build 2026-07-21 16:41). No `NAMESPACE` export
> change (106 → 106); the release was internal presentation only.

### Changed

- **Identification columns now lead every catalogue output table.** A new
  internal helper `.lead_identifier_cols()` (`catalogue_engine.R`, priority list
  `.HS_ID_COLS`) reorders the ID keys — `Trait, Chrom, Pos, Block, Sample,
  Haplotype, signature/Subtype …` — to the front of each returned/written
  `data.table`. Applied across `analyze_block_catalogue()`
  (`catalogue_analysis.R`: `subtypes`, `selected`, `S`, `N`, `top_blocks`),
  `compute_block_catalogue()` (`catalogue_gebv.R`), and the per-block effect
  summary in `run_haplotype_catalogue()` (`catalogue_main.R`). Column *values*
  are unchanged — only their order — so downstream keyed joins are unaffected.
- `tests/testthat/test-catalogue-analysis.R` / `test-catalogue.R` updated for the
  new column order; `man/compute_block_catalogue.Rd` /
  `man/run_haplotype_catalogue.Rd` regenerated.
- `DESCRIPTION`: version **0.4.7 → 0.4.8**.
- Source tarball `haploStack_0.4.8.tar.gz` in `ready_to_use/`.

## [0.4.7] - 2026-07-21 — Panel descriptors

Spec: [`dev/PRD_marker_stats_block_allele_sheet.md`](../dev/PRD_marker_stats_block_allele_sheet.md).
I/O guide: [`docs/haploStack_io_reference_0.4.7.md`](haploStack_io_reference_0.4.7.md)
(+ `_slide.html`).

### Added

- `marker_stats.R` (new): **`compute_marker_stats()`** — per-marker QC summary
  keyed on a `GenoMeta`, one tidy `data.table` row per marker in genomic order
  (`Marker, Chrom, Pos, Block, N_obs, Missing_rate, Call_rate, MAF, Major, Minor,
  Het`). All statistics come from a single integer-encoding pass over
  `geno_arrays` (promotes the internal `.encode_marker_codes` / `.marker_maf` /
  `.marker_allele_counts` helpers to a user-facing entry point). Optional TSV
  write via `data.table::fwrite` (default `marker_stats.txt`); the table is
  returned either way.
- `catalogue_gebv.R`: **`write_block_allele_sheet()`** — per haplotype block, the
  distinct marker-allele combinations (haplotype subtypes) present in the panel
  and their frequency, as one tidy `data.table`
  (`Block, Chrom, Signature, Markers, Alleles, Count, Freq`), ordered by block
  then descending `Count`. Ragged member-marker sets live inside `;`-joined
  `Markers`/`Alleles` cells (no per-marker columns), so the sheet scales to whole
  genomes. `count_by = c("haplotype","genotype")` switches copy vs distinct-carrier
  counting; missing alleles coded `""` so per-block `Count` sums exactly. Mirrors
  the internal `group_by = "signature"` catalogue path (`.compute_block_signatures`).
- `tests/testthat/test-marker-stats.R` (new): coverage for both functions.
- `docs/haploStack_io_reference_0.4.7.md` + `docs/haploStack_io_reference_0.4.7_slide.html`
  (new): user-facing I/O reference and a 16:9 analysis slide.
- Source tarball `haploStack_0.4.7.tar.gz` at the package root.

### Changed

- `DESCRIPTION`: version **0.4.6 → 0.4.7**.

## [0.4.6] - 2026-07-13

> Reconstructed from the `NAMESPACE` export delta (95 → 104) and file dates of
> `ready_to_use/haploStack_0.4.6.tar.gz` vs the 0.4.5 build.

### Added

- `simulate_dataset.R` (new): **synthetic-data generators** — `sim_dataset()`
  plus `sim_geno_meta()`, `sim_geno_matrix()`, `sim_map_data()`, `sim_block_map()`,
  and `sim_effect_df()`. Produce a self-contained `GenoMeta` (genotypes, map,
  blocks, effects) for demos, tests, and the ready-to-use templates.
- `update_geno_meta_from_matrix()` — rebuild a `GenoMeta`'s genotype arrays from a
  2-D matrix in place (round-trip helper).
- genomicSimulation session helpers **`gs_close()`** and
  **`gs_sync_group_names()`** exported (crossing-backend session management;
  groundwork for the deferred persistent-session work, `dev/PRD_gs_persistent_session.md`).

### Changed

- `DESCRIPTION`: version **0.4.5 → 0.4.6**.

## [0.4.5] - 2026-07-07

> The two blocks below are dated 2026-07-07 and were confirmed present in
> `ready_to_use/haploStack_0.4.5.tar.gz` (build 2026-07-07): `catalogue_analysis.R`
> carries `value_source`, and `fitness_functions.R` normalises by `sum(abs())`.
> This build also first bundles the `src/format_cells.cpp` + `src/marker_stats.cpp`
> Rcpp kernels and the `metrics`-gated `run_stacking_simulation()` driver (see the
> Rcpp / stacking-metrics notes at the end of this section).

### Fixed — Encoding & fitness (ready-to-use workflow review, 2026-07-07)

- **SelectionTools numeric codes are now alphabetical: `A=1, C=2, G=3, T=4,
  -=-1`** (`DEFAULT_ENCODING_DICT`). Previously `C` and `G` were swapped
  (`G=2, C=3`), so `recode_geno_meta(to_type = "selectiontools")` mis-coded C/G.
  `REVERSED_ENCODING_DICT` and the SelectionTools colour maps
  (`BASE_COLOR_MAP` / `COLOR_MAP`) were realigned to mirror the nucleotide
  palette under the new mapping. Nucleotide `A, C, G, T, -` remains the primary
  dictionary (workflow-review D2).
- **`create_fitness_index()` normalises trait weights by `sum(abs(weights))`**
  (was `sum(weights)`), so **signed** weights stay balanced and cannot blow up
  when they sum to ~0. A **negative** weight now correctly flips a
  lower-is-better trait; documented in `?create_fitness_index` (workflow-review
  D4).

### Changed — Block-catalogue analysis (BREAKING, 2026-07-07 revision)

- **`analyze_block_catalogue()` now selects ONE representative subtype per block,
  shared across all traits** (previously one subtype per *(block, trait)*). The
  block's subtype is `argmax_subtype Σ_t w_t · value(subtype, t)`, where `value`
  is the raw effect `S` or the within-trait normalized `N`, chosen by a new
  **`value_source = c("raw","normalized")`** argument. `value_source` drives both
  the summed selection and the plot fill, and is **stored on the result** so the
  plot layer honours it without the caller re-specifying (D2).
- **`normalize = "zscore"` is now a true standard score `(x − mean)/sd`** (centred
  and scaled) computed within each trait, replacing the previous divide-by-SD-only
  (zero-anchored) transform (D1). `minmax` / `rank` remain zero-anchored.
  Normalization now runs even for a single-trait analysis.
- **`trait_weights` is now signed direction/scale only** — applied to each
  subtype's value *before* the cross-trait selection sum (a negative weight flips a
  lower-is-better trait). It is no longer a weighted composite; the default is unit
  weights.
- **Removed the cross-trait adjustments:** the weighted composite `fitness`, the
  `consensus` table, the `corrected` column (`N − mean_t N`), and the
  `optimizing_potential` column / `.optimizing_potential()` helper. `fitness_fn`
  now only filters candidates / defines the per-subtype value (it no longer does
  the cross-trait combine).
- **Return shape changed (breaking):** `analyze_block_catalogue()` now returns
  `subtypes, selected` (one row per block), per-`(block, trait)` raw `S` and
  normalized `N` tables (each with a per-trait `is_top` flag), `top_blocks`,
  `gebv`, `plots`, and `value_source`. Gated writes drop `block_metric.txt` /
  `block_consensus.txt` and add `block_S.txt` / `block_N.txt`.
- **Plot layer (`catalogue_analysis_viz.R`):** `plot_block_effect_genome()` and
  `plot_block_effect_circos()` replace `metric = c("N","S","corrected")` with
  `value_source = c("raw","normalized")` (the `corrected` fill is gone);
  `plot_block_effect_manhattan()` / `plot_block_effect_density()` rename their
  `value` argument to `value_source`. Every plotter defaults `value_source` to the
  source stored on the result and accepts an explicit override for re-plotting.
  The diverging fill's white midpoint now means *average* effect under
  `normalized` (standard score) versus *zero* effect under `raw`.

### First bundled in this build (present in 0.4.5, added to the changelog 2026-07-21)

> These features carry no dated entry above but are **present in the 0.4.5
> tarball**; their exact introduction version (0.4.5 or an earlier 0.4.x) is not
> recoverable, so they are recorded against the earliest build that proves them.

- **Rcpp round 2 kernels** (`src/format_cells.cpp`, `src/marker_stats.cpp`):
  `split_genotype_cpp()` / `combine_haplotype_cpp()` (per-cell genotype split /
  combine) and `marker_allele_counts_cpp()` / `dosage_from_codes_cpp()` (per-marker
  tally / `{0,1,2}` dosage), extending the Kernel A dispatch pattern. Pure-R
  fallback retained.
- **Optional / optimised stacking metrics**: `run_stacking_simulation(metrics = )`
  (default `"gebv"`) gates the per-cycle bundle — `d2ideal`, `gene_diversity`,
  `coancestry`, and the cheap `ideal_hap_gebv` ceiling metric are computed only on
  request (measurement-only; selection and offspring are seed-identical regardless).

## [≤ 0.4.4] — earlier development history (not individually versioned)

> Everything below was accumulated across the pre-0.4.5 0.4.x line under a single
> `[Unreleased]` heading; the sources do not support splitting it into individual
> release sections. Dates on entries are their development dates.

### Performance (catalogue compute)

- `catalogue_engine.R` (new): the haplotype-catalogue compute path was
  restructured around a single **trait-invariant compute core**
  (`.catalogue_core()`). The long-format allele table, the population
  allele-frequency table, and the block/sample orderings are now built **once
  per `run_haplotype_catalogue()` call** and shared across every trait, instead
  of being rebuilt twice per trait (once at genotype level inside
  `compute_profile_localgebv()`, again at haplotype level inside
  `compute_block_catalogue()`). `compute_profile_localgebv()` and
  `compute_block_catalogue()` gained an internal `.core` argument so the run-level
  orchestrator threads the shared core into every trait (and into `future`
  workers). `compute_block_catalogue()` also gained an internal `.signatures`
  argument, and `analyze_block_catalogue()` (`.resolve_subtype_tables()`) now
  builds the core **and** the trait-invariant per-copy block-signature table once
  and shares them across all traits (previously rebuilt per trait on the
  signature + genotype-GEBV paths). Outputs are **byte-identical** to the previous
  pure-R pipeline (parity-gated by `dev/benchmark_catalogue.R`; cross-engine
  parity in `test-catalogue-engine.R`).
- `src/local_gebv.cpp` (new): **Kernel A** — `local_gebv_block_sum_cpp()`, an
  Rcpp scatter loop that sums per-block local GEBV directly from an integer-coded
  genotype matrix in one `O(n_copies × n_markers)` pass, with no long-format
  table materialised. It is **encoding-agnostic by construction**: all allele
  strings / encodings are resolved to integer codes and a numeric contribution
  table in R (`.catalogue_core()` / `.hap_effects_rcpp()`), so the kernel only
  ever indexes doubles by integer codes. Selected via the new global option
  `options(haploStack.use_rcpp = ...)` (default `TRUE`); falls back automatically
  to the pure-R `data.table` path when the option is `FALSE` **or** the kernel was
  not compiled (no compiler at install time). On `center_by_freq = TRUE` (the
  catalogue default) the per-trait block summation is ~5–6× faster than pure-R;
  results match pure-R to ~1e-13 (continuous GEBV), with all discrete structure
  (block ids, sign classification, `round(Effect, 10)` grouping, representative
  carrier) identical. Resolves dev Issue #001 (Rcpp not wired) for the catalogue
  path.
- `Rcpp` moved from `Suggests` to `Imports` + `LinkingTo`; `useDynLib` /
  `importFrom(Rcpp, sourceCpp)` wired via roxygen in `haploStack-package.R`.
- `tests/testthat/test-catalogue-engine.R` (new): core invariants, kernel input
  asserts, the `use_rcpp` toggle, and pure-R↔Rcpp parity across
  `compute_profile_localgebv()`, `compute_block_catalogue()`, and
  `run_haplotype_catalogue()`.
- `dev/benchmark_catalogue.R` (new): scalable golden-output + timing harness and
  the structure-exact / numeric-tolerant parity gate. Pre-refactor
  implementations archived verbatim in `dev/codebank/2026-06-26_catalogue_compute.R`.

### Added

- `catalogue_main.R`: `run_haplotype_catalogue()` now emits a per-block **effect
  summary** for every trait — one row per block with `BlockVariance`, `BlockMax`,
  `BlockMin` and `BlockMean` of its across-sample local GEBV, plus block location
  (`Chrom`, `Pos`). Returned as the new `block_summary` element and written to
  `Block_Summary_<ld>_t<tol>_<trait>.txt` (one file per trait). It also feeds the
  per-trait variance scatter (`plot_manhattan_variance()`), so the plotted
  variance and the on-disk table are guaranteed consistent.

### Changed

- `catalogue_viz.R`: `plot_manhattan_effect()` and `plot_manhattan_variance()`
  now share a single internal genomic-scatter body (`.plot_manhattan()`); the two
  figures differ only in the value column drawn (`Effect` vs `BlockVariance`) and
  the effect plot's dashed zero reference line. Public signatures are unchanged.

- `genomic_prediction.R`: `fit_genomic_prediction()` — bridge from a `GenoMeta`
  to the standard genomic-prediction engines. Transforms the haplo-wise
  `geno_arrays` into a numeric alt-allele dosage matrix and fits `"rrblup"`
  (`rrBLUP::mixed.solve`), `"gblup"` (`rrBLUP::A.mat` + `kin.blup`, with marker
  effects back-solved via the SNP-BLUP equivalence) or `"bayesc"`
  (`BGLR::BGLR`). Returns the estimates already shaped for re-use:
  `marker_effects` as a **named list with one `(Marker, Allele, Effect)` table
  per trait** (the effect_list shape consumed by `run_haplotype_catalogue()` /
  `compute_profile_localgebv()`), a GEBV `phenotype` table (`Sample, GEBV_<trait>`,
  drop-in for `attach_phenotype()`), the fitted `model`, and a `geno_meta` copy
  with both attached and revalidated (fitted effects merged into any already
  stored, one table per trait). `rrBLUP` / `BGLR` added to `Suggests`
  (guarded with `requireNamespace()`).
- `tests/testthat/test-genomic-prediction.R`: coverage for the dosage transform
  (`{0,1,2}` values, shape, verbatim marker names), all three model fits
  (skipped when the backend is absent), the GenoMeta round-trip through
  `compute_profile_localgebv()`, and error paths.
- `phenotype.R`: new `phenotype` element on `GenoMeta` — an optional
  `Sample` + numeric-trait data.frame, keyed on `Sample` the way `marker_effects`
  is keyed on `Marker`. Threaded through `new_geno_meta()`
  (`.standardise_phenotype()` aligns it to the genotype sample set: drop unknown
  samples, NA-fill missing rows, reorder), `is_geno_meta()` (optional validation),
  and `filter_geno_meta()` (sample-pruned). Round-trips via `save_geno_meta()` /
  `load_geno_meta()`.
- `phenotype.R`: `read_phenotype()` (delimited `.pheno`/CSV/TSV reader),
  `attach_phenotype()` (filter to shared samples + store), and `get_phenotype()`
  (accessor). Ports the data-loading half of the standalone trait-analysis script
  into the GenoMeta pipeline.
- `phenotype_viz.R`: `plot_phenotype_distribution()`, `plot_phenotype_correlation()`,
  and `phenotype_summary()` — portable (accept a `GenoMeta` or any data.frame with
  a `Sample` column) ggplot2 reimplementations of the distribution/correlation
  tasks, with no `corrplot`/`GGally`/`ggridges` dependency.
- `phenotype_viz.R`: `plot_phenotype_bars()` — per-sample dodged bar chart of up to
  two traits, with a `"vertical"` (samples on X, direct) and `"horizontal"`
  (samples on Y, merge-ready) orientation.
- `phenotype_viz.R`: `plot_meta_alignment()` — merges ggplots that share a discrete
  sample axis (e.g. a genotype-level `plot_localgebv_alignment()` panel + horizontal
  phenotype bars) into one figure: reorders every panel's sample axis to the
  reference order for identical ticks and pools the legends into a shared bottom
  row (uses `cowplot`, already in Suggests).
- `tests/testthat/test-phenotype.R`: coverage for sample intersection / NA-fill /
  warnings, validation, filter pruning, save/load round-trip, and plot outputs.

- `io/format_helper.R`: `GENO_FORMAT_PRESETS` — exported registry of genotype
  file-format *parameter sets* (`default`, `PL`, `ST`, `GS`), each describing
  on-disk orientation, header, embedded FAM columns, map column order, allele
  delimiter, and encoding. [FORMAT: PLINK/ST/GS]
- `io/format_helper.R`: `resolve_geno_format()` — dictionary resolver that maps
  a preset name (or alias, e.g. `"plink_ped"`), or a user-supplied named list,
  to a validated parameter set, applying per-field overrides. Shared by the
  reader and the writers so each format is described once. [FORMAT: PLINK/ST/GS]
- `io/format_helper.R`: `.read_genotype_files()` and `.read_map_file()` —
  internal file readers that, given a resolved parameter set, return the
  genotype (re-oriented to the fixed `sample × marker` layout), map
  (`Chrom, Marker, Pos`), and fam tables **without modifying their contents**.
- `tests/testthat/test-format.R`: tests for `resolve_geno_format()`, a PL-preset
  read with an embedded FAM block, and a preset-driven `write_geno_meta()`
  round-trip.

### Changed

- `catalogue_analysis_viz.R`: `plot_block_effect_genome()` now draws each block
  as a `geom_rect()` rectangle spanning its genomic extent
  (`Block_start`–`Block_end`, sourced from the selected-subtype table / GenoMeta
  block locations) instead of a single-point `geom_tile()`. Traits map to integer
  y-bands; point/singleton blocks get a small visible width. `geom_rect` added to
  the `ggplot2` imports.
- **Marker-effects format unified to a named list of per-trait tables.** The
  canonical shape for `geno_meta$marker_effects`, the `effect_list` of
  `run_haplotype_catalogue()`, and the `marker_effects` returned by
  `fit_genomic_prediction()` is now a **named list** (one `(Marker, Allele,
  Effect)` data.frame per trait, the list name giving the trait).
  `new_geno_meta()` / `filter_geno_meta()` / `is_geno_meta()` accept and
  normalise this list; a bare data.frame (optionally with a `Trait` column to
  split on) is still accepted for back-compatibility. Marker effects are
  **never required** for a valid `GenoMeta`.
- `catalogue_main.R`: `run_haplotype_catalogue()` — `effect_list` is now
  optional (default `NULL`) and falls back to `geno_meta$marker_effects`, so
  effects estimated by `fit_genomic_prediction()` are picked up automatically.
  `trait_names` defaults to the resolved effect traits.
- `catalogue_gebv.R`: `compute_profile_localgebv()` and
  `compute_block_catalogue()` gain a `trait` argument that selects one table
  from a stored per-trait `marker_effects` list when `effect_df` is `NULL` (a
  lone stored trait is picked automatically).
- `genomic_prediction.R`: `fit_genomic_prediction()` returns `marker_effects`
  as a one-element named list (no inline `Trait` column) and **merges** the
  fitted table into any effects already on the `GenoMeta` rather than replacing
  the whole field.
- `io/format_helper.R`: `build_geno_meta()` now enforces a **fixed input
  contract** — `geno_2d` must be `(n_samples × n_markers)` with named rows
  (samples) and columns (markers); `map` columns must be in the order
  `Chrom, Marker, Pos`; `fam` columns in the order
  `Family, Sample, Parent1, Parent2, Sex`. Orientation auto-detection and
  wide-format reshaping were removed (now handled by the reader). The function
  exposes only the changeable format parameters `split_char`, `pos_unit`, and
  `encoding`.
- `io/format_helper.R`: read-side encoding vocabulary renamed from
  `"number"`/`"char"` to `"selectiontools"`/`"nucleotide"`, matching
  `detect_genotype_encoding()` and the catalogue pipeline. The resolved encoding
  is now recorded in `metadata$encoding`. Legacy `"number"`/`"char"` remain
  accepted as silent aliases for back-compatibility.
- `io/format_helper.R`: `read_genotype()` reworked to delegate file reading to
  `.read_genotype_files()` and assembly to `build_geno_meta()`; it accepts a
  `preset`/parameter-set and the `presets` registry instead of the previous
  `row_axis`/`map_preset`/`encode` arguments.
- `io/format_helper.R`: `geno_meta_to_matrix()` and `write_geno_meta()` are now
  preset-aware (same parameter-set mechanism as the reader). They accept
  `preset`/`presets` and retain the legacy `output_format` argument
  (`"plink_ped"`/`"st"`/`"gs"`) as a deprecated alias so the Shiny GUI continues
  to work. Output orientation, allele delimiter, and map column order are taken
  from the resolved spec; the genomicSimulation write orientation is now
  marker-as-row, consistent with the GS read layout.

### Fixed

- `catalogue_profile.R`: `geno_meta_views()` no longer mangles non-syntactic
  marker IDs. The internal allele-table builder (`.bind_haplotypes()`) now calls
  `as.data.frame(..., check.names = FALSE)`, so IDs such as `"1H-123"` survive
  verbatim instead of becoming `"X1H.123"`. Previously the mangled names broke the
  `Marker` join in `compute_profile_localgebv()` (silently zeroing local GEBVs)
  and propagated mangled names into any estimated `marker_effects`.

### Migration notes

- Callers of `build_geno_meta()` must pass `geno_2d` already in
  `sample × marker` orientation. In-memory data previously relying on automatic
  orientation detection should be transposed beforehand.
- `read_genotype(..., row_axis=, map_preset=, encoding="number"/"char")` should
  migrate to `read_genotype(..., preset=)` (e.g. `preset = "ST"`), or pass a
  custom parameter-set list / `encoding = "selectiontools"|"nucleotide"`.
