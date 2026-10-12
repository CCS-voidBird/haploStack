# haploStack

**R package for haplotype block analysis and stacking simulation in diploid crop breeding programmes.**

Developed for the GRDC Barley Phase 2 project at the University of Queensland.
Generalised for use with any biparental or multiparental diploid breeding programme.

> ### 0.7.11 — a blackout is an emergency; a missing understudy is not
>
> Stage 1 of the GBS multi-subtype core upgrade (thread
> [#036](project_record/threads/036-gbs-multi-subtype.md);
> `dev/UPGRADE_PLAN_gbs_multi_subtype.md`, PRD decision D5, Q26 a–c accepted
> 2026-10-09). GBS's metric treated every target state as an independent
> coordinate, so a block whose rank-1 subtype was fixed still fired a full
> absence penalty for its rank-2 twin — picks were spent re-acquiring a +0.03
> understudy while true blackouts waited ("phantom emergencies"). The absence
> penalty now carries a within-block **rank multiplier**. Default off:
> `penalty_rank_scale = 1` is **bit-identical to 0.7.10** (T1 gate across
> method × objective × engine × penalty mode × preset, fixed seeds, plus a
> fixed-seed driver-run equality). No export change — **134 exports**.
>
> ### 0.7.10 — a block map first; LD blocking is a side route, named after its backend
>
> Owner direction of the same day: haploStack's recommended use is **a block
> map the user already has**, given to `read_plink_ped(block_file =)`,
> `build_geno_meta(blocks =)` or `attach_block_info()`, with no LD parameter
> to configure. The LD engines are built-in resolvers for users who really
> want to build blocks locally. HapSelect is the default, SelectionTools the
> alternative ("restored, but not stable"). Each gets a function name that
> cannot be mistaken for the main route or for the other engine: `hapS_` /
> `st_`. Thread [#035](project_record/threads/035-ld-block-resolvers.md). 131
> → **134 exports**.
>
> ### 0.7.9 — SelectionTools LD blocks straight from a GenoMeta; a blank line is not a record
>
> Two owner requests of the same day. First: check every bridge between a
> `GenoMeta` and SelectionTools, and give LD blocking a GenoMeta-in route that
> takes the r² threshold and tolerance and handles SelectionTools' files
> itself (thread [#034](project_record/threads/034-selectiontools-bridge.md),
> INDEX discrepancy **34**, intake 04). Second: fix the blank-line defect of
> the shared text reader (INDEX discrepancy **33**, intake 08). 130 → **131
> exports**.
>
> ### 0.7.8 — an exported effect file reads back with the defaults
>
> `export_marker_effects()` writes one **tab-separated file with a header
> row** (`Marker\tAllele\tEffect`) per trait, and its roxygen has promised
> since 0.5.2 that *"a written file round-trips through `read_effect_file()`
> unchanged"*. It never did with the reader's defaults (`sep = " ", header =
> FALSE`): the whole line landed in one column and the read stopped with
> *"transform_effects: need at least 3 columns"*. Found by the `haplostack-
> analysis` skill's contract smoke test (intake 01); closes INDEX discrepancy
> **31**. Released with the workshop blocking fix (intake 03, discrepancy
> **32**, see *Deliverables*).
>
> ### 0.7.7 — a genotype file is read as text
>
> `read_plink_ped()`, `read_effect_file()` and `read_dart_genotype()` parsed
> their file with `utils::read.table()`, whose `type.convert()` accepts `T` and
> `F` as logical literals — so a column whose every value is `"T"` came back as
> `TRUE`, and in a nucleotide genotype file that column is an allele. A PED then
> failed to read outright (*"cells contain multi-character values"*); an effect
> table keyed on `TRUE` matched no genotype and, since an absent allele scores
> 0, produced **flat GEBVs with no error at read time**; a DArT report
> mis-parsed its sample-column boundary. Every reader now goes through the new
> internal `.read_table_chr()`, which reads every column as `character` and
> leaves numeric coercion to the caller. It is the allele *column* that has to
> be uniform, not the marker — on a 1107 x 9856 barley PED, 511 of 19 712
> columns, from 511 markers, none of them monomorphic. Side effect, intended:
> identifiers keep their spelling, so a sample id written `007` no longer
> becomes `"7"`. No export change; containers on disk are unaffected. The read
> path from file to `metadata` is written up in
> `docs/haploStack_io_reference_0.7.7.md`.
>
> ### 0.7.6 — views keep the observed allele of a half-called cell
>
> `geno_meta_views()` now builds its per-copy allele table straight from the
> ploidy layers of `geno_arrays`. It used to re-split the collapsed genotype
> string, which is `NA` whenever *any* copy is missing — so a half-called cell
> `(A, NA)` came back as `(NA, NA)` and its observed allele was dropped from
> every consumer: local GEBVs at both levels, the block catalogue, block
> signatures, clustering, block IBD / signal scoring and the marker-stats
> views, and with them every simulation offspring metric from cycle 1. Founders
> read from a PED are never half-called, but derived populations carry
> 1.6–3.6 % such cells. **Numbers change only on containers that hold a
> half-called cell**; one without any is byte-identical to 0.7.5. No export
> change.
>
> ### 0.7.5 — `max_F1mean` (top elite) in the simulation summary
>
> `run_stacking_simulation()$summary` gains `max_F1mean`, the breeding value
> of the cycle's best single offspring, beside the population mean — the
> realised counterpart of a set-level "ultimate genotype" objective. The
> genomicSimulation backend now restores missing alleles to `NA` on read-back
> (they came back as a literal `"N"` allele and `is_geno_meta()` warned every
> cycle). The barley benchmark was re-run on both crossing backends at 100
> cycles with Kemper's GBS preset added; the backend does not change the
> ranking — see the changelog.
>
> ### 0.7.4 — `ohs_method = "GA_coverage"` removed
>
> The binary block-coverage engine is gone (it was the worst set on every
> metric in the 2026-09-17 review); `ohs_method` is `"GA_pairwise"` or
> `"Greedy"`, and the deprecated `select_method = "GA"` now maps to
> `GA_pairwise`. `genetic_algorithm_selection()` itself stays exported.
>
> ### 0.7.3 — OHS: elite subsetting of the families
>
> `ohs_selection()`, `ohs_pair_selection()` and `run_stacking_simulation()`
> gain `elite_ratio` / `elite_level` (`ohs_elite_ratio` / `ohs_elite_level` in
> the driver): the OHS objective may count only the best
> `ceiling(ratio · P)` of a set's `P` families — per set (`"pair"`, default;
> ratio → 0 = the best single cross) or per block (`"block"`; ratio → 0 = the
> best two haplotypes of every block from any two parents, Kemper's "ultimate
> genotype"). Ratio 1 is the previous objective, and the crossing design has
> the same size at every ratio. **Changes published OHS numbers:** the `Greedy`
> engine now starts from the best row sum of `M` (it used to start from the
> first sample column), and the OHS GA engines restore the RNG on exit, so
> replicates of a multi-rep OHS run differ again.
>
> ### 0.7.2 — direct pair selection
>
> The stacking driver can now **select the crosses themselves** instead of a
> parent subgroup: `run_stacking_simulation(selection_unit = "pairs",
> pair_budget = "families" | "parents", ...)`. Every scheme scores *families*
> (the expected F1 of one cross) and one greedy adds families until the budget
> stops it — TS by the mid-parent (or usefulness), OHS by the family's stacking
> ceiling `M[i, j]`, GBS by the ideal-genotype objective on the family's expected
> offspring, which reduces exactly to the parents-first frequency for a full
> traversal. New exports `ohs_pair_selection()`, `gbs_pair_selection()`;
> `summary` gains `n_families` / `n_parents_used`, `plans` gains `pair_score`.
> Parents-first runs are unchanged (`avoid_selfing` now defaults to `TRUE`, which
> the all-pairs design never used). Every valid combination of scheme × unit ×
> budget × option: [`vignettes/selection_grid.md`](vignettes/selection_grid.md).
>
> **Focus on `TS`, `OHS`, `GBS`.** `OCS` is kept for compatibility but has known
> issues (kinship-penalised truncation with a silent fall-back — not true
> optimum-contribution selection) and has no family form; a `KCT` rename and a
> real `OPT` scheme are planned. In a single call with `replications > 1`
> OCS re-seeds the RNG each cycle, so its replicates are identical — run one
> call per seed (the OHS GA engines restore the RNG since 0.7.3).
>
> **0.7.0 changed published GBS numbers** (`gbs_preset = "kemper"`, Kemper's
> `LOG`, the bounded risk × hope penalty scale, `ideal_scope = "founder"`), and
> **0.6.0 changed missing-allele handling and the catalogue centring default**
> — see [`docs/CHANGELOG.md`](docs/CHANGELOG.md) and
> [`docs/haploStack_io_reference_0.6.0.md`](docs/haploStack_io_reference_0.6.0.md)
> before upgrading a running analysis.

---

## Quick Start

### Install

```r
# From local source
devtools::install("path/to/Haplotype_workflow_R")
library(haploStack)
```

### Minimal Example (Simulated Data)

```r
library(haploStack)

# 1. Generate simulated dataset
geno_mat  <- sim_geno_matrix(n_samples = 10, n_markers = 20)
map_df    <- sim_map_data(rownames(geno_mat))
block_map <- sim_block_map(map_df, block_size = 4)
effect_df <- sim_effect_df(rownames(geno_mat), alleles = c("1","4"))

# 2. Assemble a GenoMeta object (genotypes + map + blocks in one container)
meta <- build_geno_meta(geno_mat, map_df, blocks = block_map)

# 3. Compute local GEBVs directly from the GenoMeta
gebv_mat <- compute_profile_localgebv(meta, effect_df)

# 4. Select parents
parents  <- greedy_selection(gebv_mat, n_parents = 6)

# 5. Design crosses
crosses  <- make_cross_design(parents, method = "all_pairs", n_offspring = 4)

# 6. Visualise
plot_haplotype_alignment(meta, show_block_labels = TRUE)
```

`GenoMeta` is the single central object: it carries the genotype array, map,
block assignments, a FAM-structured `sample_info`
(`Family, Sample, Parent1, Parent2, Sex`), and an optional `marker_effects`
copy. Derive the wide allele / block tables on demand with
`geno_meta_views(meta)`; subset with `filter_geno_meta(meta, ...)`. See
[`docs/geno_meta_structure.md`](docs/geno_meta_structure.md) for the full
element/provenance reference.

**Blocks come from a block map you supply** — `build_geno_meta(blocks =)`,
`read_plink_ped(block_file =)` or `attach_block_info()` — and need no LD
parameter. That is the recommended route. Only when no map exists, build LD
blocks locally with a resolver named after its backend:

- `hapS_ld_blocks()` uses HapSelect and is the default;
- `st_ld_blocks()` uses SelectionTools and is experimental.

Both keep the container's effects, phenotype and pedigree. They do **not**
partition alike at the same threshold and tolerance, so compare block sets only
within one backend.

Steps 3–5 above assume the marker effects are already known. To instead
**estimate** marker effects and GEBVs from phenotypes, use
`fit_genomic_prediction(meta, phenotype, model = "rrblup")` (also `"gblup"` /
`"bayesc"`, via the optional `rrBLUP` / `BGLR` packages). It returns a
`marker_effects` table and a GEBV `phenotype` table already shaped to feed back
into the GenoMeta and on into `compute_profile_localgebv()`.

---

## Pipeline Overview

```
┌────────────────────────────────────────────────────────────────────┐
│  STAGE 1: Format Standardisation                                   │
│    detect_genotype_encoding()  →  process_genotype()              │
│    transform_genotype()        →  read_effect_file()              │
├────────────────────────────────────────────────────────────────────┤
│  STAGE 2: Blocks — bring your own block map (recommended)          │
│    read_plink_ped(block_file =) / build_geno_meta(blocks =)        │
│    attach_block_info()         →  build_block_location_info()     │
│    no map? local LD resolvers (side route):                        │
│    hapS_ld_blocks() [HapSelect, default] · st_ld_blocks() [ST, exp.]│
├────────────────────────────────────────────────────────────────────┤
│  STAGE 3: Haplotype Catalogue                                      │
│    geno_meta_views()           →  filter_geno_meta()              │
│    fit_genomic_prediction()    →  compute_profile_localgebv()     │
│    run_haplotype_catalogue()   →  plot_manhattan_effect()         │
├────────────────────────────────────────────────────────────────────┤
│  STAGE 4: Stacking Simulation                                      │
│    greedy_selection()          →  genetic_algorithm_selection()   │
│    ocs_selection()             →  make_cross_design()             │
│    run_stacking_simulation()   →  plot_simulation_trajectories()  │
└────────────────────────────────────────────────────────────────────┘
```

---

## Package Structure

```
Haplotype_workflow_R/
├── DESCRIPTION                  Package metadata
├── NAMESPACE                    Exported functions
├── R/
│   ├── alphabet.R               Constants: allele dicts, color maps
│   ├── format_detect.R          detect_genotype_encoding(), orientation helpers
│   ├── format_process.R         process_genotype(), map normalisation
│   ├── format_transform.R       transform_genotype(), transform_effects(), read_effect_file()
│   ├── format_split.R           split_haplotypes(), combine_haplotypes(), translate_alleles()
│   ├── catalogue_profile.R      create/filter/save/load/recode profile
│   ├── catalogue_gebv.R         compute_profile_localgebv()
│   ├── genomic_prediction.R     fit_genomic_prediction() (rrBLUP / BGLR)
│   ├── catalogue_blockmap.R     build_block_location_info(), expand_block_info()
│   ├── blocks_ld_common.R       what the local LD resolvers share (a side route)
│   ├── blocks_hapselect.R       hapS_ld_blocks() — HapSelect, the default resolver
│   ├── blocks_selectiontools.R  st_ld_blocks(), st_ld_blocks_files() — experimental
│   ├── blocks_deprecated.R      the old define_haplotype_blocks*() names (warn, forward)
│   ├── catalogue_viz.R          plot_haplotype_alignment(), Manhattan plots, count plots
│   ├── catalogue_family.R       trace_haplotype_family()
│   ├── catalogue_main.R         run_haplotype_catalogue() orchestrator
│   ├── fitness_functions.R      create/update/write fitness index
│   ├── stacking_selection.R     greedy, GA, OCS selection; make_cross_design()
│   ├── stacking_main.R          run_stacking_simulation()
│   ├── stacking_pairs.R         direct pair selection: ohs_pair_selection(), gbs_pair_selection()
│   └── stacking_monitor.R       plot_simulation_trajectories(), summarise_simulation()
├── tests/testthat/
│   ├── helper-sim_data.R        Simulation data generators
│   ├── test-format.R            Format module QA (13 tests)
│   ├── test-catalogue.R         Catalogue module QA (22 tests)
│   ├── test-fitness.R           Fitness module QA (11 tests)
│   └── test-stacking.R          Stacking module QA (12 tests)
├── vignettes/
│   ├── quick_start.Rmd          ← START HERE: beginner guide
│   ├── simulation_catalogue.Rmd Full catalogue workflow with simulated data
│   └── simulation_stacking.Rmd  Full stacking workflow with simulated data
├── dev/
│   ├── development_sheet.md     Status, open issues, milestones
│   ├── questions_sheet.md       Design questions needing user input
│   └── qa_suggestions.md        Future development ideas
└── inst/extdata/                (Example data files — to be added)
```

---

## Key Data Structures

| Object | Description | Key columns / fields |
|--------|-------------|---------------------|
| Genotype matrix | Markers × samples character matrix | rownames = marker IDs, colnames = sample IDs, values = "1/4" or "A/T" |
| Map data frame | One row per marker | `Marker`, `Chrom`, `Pos` (bp) |
| Block map | One row per marker-block pair | `Chrom`, `Marker`, `Block` |
| Effect data frame | One row per (Marker, Allele) | col1 = Marker, col2 = Allele, col3 = Effect |
| Haplotype profile | Central object | `$alleles`, `$marker_meta`, `$block_meta`, `$samples`, `$encoding` |
| GEBV matrix | Block × Sample wide table | `Block` + one column per sample |

---

## Genotype Encoding Support

| Format | Example cell | Encoding key |
|--------|-------------|--------------|
| SelectionTools | `"1/4"` | `"selectiontools"` |
| Nucleotide (slash) | `"A/T"` | `"nucleotide"` |
| Nucleotide (concat) | `"AT"` | `"nucleotide"` |
| Dosage | `"0"`, `"1"`, `"2"` | `"numeric"` |

Auto-detection: `detect_genotype_encoding(geno_mat)` — no manual specification needed.

---

## Performance (catalogue compute)

`run_haplotype_catalogue()` builds its trait-invariant compute structures (the
long-format allele table, allele frequencies, block/sample orderings) **once per
run** and shares them across all traits. The per-block local-GEBV summation is
accelerated by an Rcpp kernel (`src/local_gebv.cpp`) that works on integer-coded
genotypes — all allele/encoding handling stays in R, so it is safe across
SelectionTools / nucleotide / numeric encodings.

The kernel is controlled by a global option (default on):

```r
options(haploStack.use_rcpp = TRUE)   # use the compiled kernel when available
options(haploStack.use_rcpp = FALSE)  # force the pure-R data.table path
```

The pure-R path is used automatically whenever the compiled kernel is not
loadable (e.g. the option is `FALSE`, or a binary build shipped without it) — the
two paths give identical results: continuous GEBV values agree to ~1e-13, and
block ids, sign classification and grouping are exact. Building from source
installs the kernel and needs the usual R toolchain (Rtools on Windows).

---

## Selection Methods

| Method | Function | Description |
|--------|----------|-------------|
| Greedy | `greedy_selection()` | Iteratively picks candidates that maximise uncovered block count |
| GA | `genetic_algorithm_selection()` | Evolutionary optimisation of block coverage |
| OCS/ARCSIN | `ocs_selection()` | Kinship-constrained GEBV maximisation |

---

## Running Tests

```r
# Run all tests
devtools::test()

# Run a specific module
testthat::test_file("tests/testthat/test-catalogue.R")
```

---

## Dependencies

**Required:** `data.table`, `dplyr`, `ggplot2`, `tidyr`, `reshape2`, `stringr`

**Optional:**
- `HapSelect` — local LD blocking when no block map exists (`hapS_ld_blocks()`,
  the default resolver); GitHub only: `remotes::install_github("wrshf7/HapSelect")`
- `SelectionTools` — the experimental LD resolver (`st_ld_blocks()`,
  `st_ld_blocks_files()`)
- `genomicSimulation` — Linkage-aware crossing simulation (`run_stacking_simulation()`)
- `rrBLUP` — RR-BLUP / G-BLUP marker-effect & GEBV estimation (`fit_genomic_prediction()`)
- `BGLR` — Bayesian (BayesC) marker-effect & GEBV estimation (`fit_genomic_prediction()`)
- `Rcpp` — C++ speed-ups for large-scale haplotype splitting and GEBV computation
- `cowplot` — Combined trajectory plots

---

## Development Status

See [`dev/development_sheet.md`](dev/development_sheet.md) for detailed status.

Open questions requiring user input: [`dev/questions_sheet.md`](dev/questions_sheet.md)

Future development suggestions: [`dev/qa_suggestions.md`](dev/qa_suggestions.md)

---

## Reference Pipeline

This package converts the workflow in:
`../Block_analysis/Haplotype_Workflow/`

Original pipeline documentation:
- `_Docs/pipeline_description.md` — Architecture overview
- `_Docs/package_conversion_guide.md` — Conversion decisions
- `_Docs/workflow_visualization.md` — Function call chains
