# =============================================================================
# test-stacking-driver.R
#
# Parameterised simulation driver (R/stacking_main.R, PRD #013 issue 03 /
# Phase 4): run_stacking_simulation() runs the recurrent route
# select -> design -> cross -> advance over replications x schemes x cycles,
# GenoMeta-native, and emits the origin-compatible output contract (FR5).
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

# Small founder GenoMeta + effects reused across the driver tests.
.driver_fixture <- function(seed = 5L) {
  m <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = seed)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))
  list(meta = m, eff = e)
}

PLANS_COLS <- c("CrossID", "Parent1", "Parent2", "F1count", "F1mean",
                "select_type", "cycle", "rep", "sel_scheme",
                "D2Ideal", "coancestry", "ideal_phenotype")
SUMMARY_COLS <- c("cycle", "select_type", "rep", "sel_scheme", "mean_F1mean",
                  "var_F1mean", "mean_D2Ideal", "mean_coancestry",
                  "mean_gene_diversity", "mean_ideal_phenotype", "ideal_hap_gebv",
                  "marker_size")

# ── output contract ───────────────────────────────────────────────────────────

test_that("a single-scheme run returns the origin output contract", {
  fx  <- .driver_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
    top_founders = 4L, num_offspring = 2L, scheme = "TS", metrics = "all",
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)

  expect_setequal(names(res),
                  c("plans", "summary", "ideal_genotypes",
                    "final_populations", "trajectories"))
  expect_named(res$plans, PLANS_COLS)               # exact columns + order
  expect_named(res$summary, SUMMARY_COLS)
  expect_named(res$ideal_genotypes,                 # built because d2ideal ∈ "all"
               c("Block", "Effect", "CFreq", "TargetFreq"))
  expect_true(all(vapply(res$final_populations, is_geno_meta, logical(1L))))
})

test_that("CrossID and sample names are compact and history-free across cycles", {
  fx  <- .driver_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 5L,
    top_founders = 6L, num_offspring = 3L, scheme = "TS",
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)

  # per-cross CrossID == C<cycle>F<family>; bounded, not the old PxPxP... history
  expect_true(all(grepl("^C[0-9]+F[0-9]+$", res$plans$CrossID)))
  expect_lt(max(nchar(res$plans$CrossID)), 16L)

  # offspring sample names (final populations) == C<cycle>F<family>_<child>
  fin <- res$final_populations[[1L]]
  snames <- dimnames(fin$geno_arrays)[[1L]]
  expect_true(all(grepl("^C[0-9]+F[0-9]+_[0-9]+$", snames)))
  expect_lt(max(nchar(snames)), 16L)

  # lineage is preserved in the columns, not the id: late-cycle parents are
  # themselves compact previous-generation ids.
  late <- res$plans[cycle == max(cycle)]
  expect_true(all(grepl("^C[0-9]+F[0-9]+_[0-9]+$", late$Parent1)))
})

test_that("summary carries a cycle-0 founder record", {
  fx  <- .driver_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
    top_founders = 4L, num_offspring = 2L, scheme = "TS",
    output_dir = NULL, verbose = FALSE)
  z <- res$summary[cycle == 0L]
  expect_equal(nrow(z), 1L)
  expect_equal(z$select_type, "founder")
})

# ── scheme dispatch + comparison ──────────────────────────────────────────────

test_that("each scheme runs end-to-end and produces valid plans", {
  fx <- .driver_fixture()
  for (sc in c("TS", "OHS", "OCS", "GBS")) {
    res <- run_stacking_simulation(
      fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
      top_founders = 4L, num_offspring = 2L, scheme = sc,
      ohs_method = "Greedy", select_top_blocks = 50L,
      output_dir = NULL, verbose = FALSE)
    expect_named(res$plans, PLANS_COLS, info = sc)
    expect_gt(nrow(res$plans), 0L)
    expect_true(all(res$plans$F1count >= 1L), info = sc)
    expect_identical(unique(res$summary$sel_scheme), sc, info = sc)
  }
})

test_that("a scheme vector compares schemes in one run", {
  fx  <- .driver_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 2L, select_cycles = 2L,
    top_founders = 4L, num_offspring = 2L,
    scheme = c("TS", "GBS"), select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE)

  expect_setequal(unique(res$summary$sel_scheme), c("TS", "GBS"))
  expect_setequal(names(res$final_populations),
                  c("TS_rep1", "TS_rep2", "GBS_rep1", "GBS_rep2"))
})

# ── selfing / ga_start_by route ───────────────────────────────────────────────

test_that("selfing_cycles interleaves selfing generations and ga_start_by delays the scheme", {
  fx  <- .driver_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
    selfing_cycles = 1L, ga_start_by = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "OHS", ohs_method = "Greedy",
    output_dir = NULL, verbose = FALSE)

  # whole_cycles = select_cycles * (selfing_cycles + 1) = 4, plus the founder row
  expect_equal(max(res$summary$cycle), 4L)
  st <- res$summary[order(cycle)]$select_type
  expect_equal(st, c("founder", "TS", "selfing", "OHS", "selfing"))
})

# ── back-compatibility ─────────────────────────────────────────────────────────

test_that("catalogue_results + select_method back-compat path works", {
  fx      <- .driver_fixture()
  cat_res <- list(SimTrait = list(geno_meta = fx$meta))
  res <- run_stacking_simulation(
    catalogue_results = cat_res, trait_name = "SimTrait", effect_df = fx$eff,
    replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, select_method = "GA", output_dir = NULL, verbose = FALSE)

  # GA -> OHS(GA_coverage); trajectories compat view is still emitted.
  expect_identical(unique(res$summary$sel_scheme), "OHS")
  expect_true("trajectories" %in% names(res))
  expect_named(res$trajectories,
               c("Rep", "Cycle", "Method", "MeanGEBV", "SDev", "BlockCoverage"))
  expect_gt(nrow(res$trajectories), 0L)
})

test_that("effects attached to the GenoMeta are used when effect_df is omitted", {
  fx <- .driver_fixture()
  m  <- fx$meta
  m$marker_effects <- fx$eff
  res <- run_stacking_simulation(
    m, replications = 1L, select_cycles = 1L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", output_dir = NULL, verbose = FALSE)
  expect_gt(nrow(res$plans), 0L)
})

test_that("missing effects raise a clear error", {
  fx <- .driver_fixture()
  expect_error(
    run_stacking_simulation(fx$meta, replications = 1L, select_cycles = 1L,
                            scheme = "TS", output_dir = NULL, verbose = FALSE),
    "no effects supplied")
})
