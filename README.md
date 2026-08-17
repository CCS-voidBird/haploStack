# haploStack

**R package for haplotype block analysis and stacking simulation in diploid crop breeding programmes.**

Developed for the GRDC Barley Phase 2 project at the University of Queensland.
Generalised for use with any biparental or multiparental diploid breeding programme.

---

## Quick Start

### Install

```r
# From GitHub
remotes::install_github("uq-bioinformatics/haploStack")

# From a local clone of this repository
devtools::install("path/to/haploStack")

# From a release source tarball
install.packages("haploStack_0.5.5.tar.gz", repos = NULL, type = "source")

library(haploStack)
```

Compiling `src/*.cpp` needs the usual R build toolchain (Rtools on Windows,
`r-base-dev` / Xcode command-line tools elsewhere).

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
│  STAGE 2: Block Definition (requires SelectionTools)               │
│    define_haplotype_blocks()   →  build_block_location_info()     │
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
haploStack/
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
│   ├── catalogue_blockmap.R     define_haplotype_blocks(), build_block_location_info()
│   ├── catalogue_viz.R          plot_haplotype_alignment(), Manhattan plots, count plots
│   ├── catalogue_family.R       trace_haplotype_family()
│   ├── catalogue_main.R         run_haplotype_catalogue() orchestrator
│   ├── fitness_functions.R      create/update/write fitness index
│   ├── stacking_selection.R     greedy, GA, OCS selection; make_cross_design()
│   ├── stacking_main.R          run_stacking_simulation()
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
├── src/                         Rcpp kernels (compiled on install)
└── docs/
    ├── CHANGELOG.md             Version history
    ├── geno_meta_structure.md   GenoMeta container reference
    └── haploStack_io_reference_0.5.5.md   I/O contract for the current release
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
- `SelectionTools` — LD-based block definition (`define_haplotype_blocks()`)
- `genomicSimulation` — Linkage-aware crossing simulation (`run_stacking_simulation()`)
- `rrBLUP` — RR-BLUP / G-BLUP marker-effect & GEBV estimation (`fit_genomic_prediction()`)
- `BGLR` — Bayesian (BayesC) marker-effect & GEBV estimation (`fit_genomic_prediction()`)
- `Rcpp` — C++ speed-ups for large-scale haplotype splitting and GEBV computation
- `cowplot` — Combined trajectory plots

---

## Documentation

- [`docs/CHANGELOG.md`](docs/CHANGELOG.md) — version history and release notes
- [`docs/geno_meta_structure.md`](docs/geno_meta_structure.md) — full `GenoMeta`
  element / provenance reference
- [`docs/haploStack_io_reference_0.5.5.md`](docs/haploStack_io_reference_0.5.5.md) —
  I/O contract for the current release
- `vignettes/` — worked examples; start with `quick_start.Rmd`

---

## Citation & Licence

Developed at the University of Queensland for the GRDC Barley Phase 2 project.
Released under the MIT licence — see [`LICENSE`](LICENSE).
