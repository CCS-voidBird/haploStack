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
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("A", "T"))
  list(meta = m, eff = e)
}

PLANS_COLS <- c("CrossID", "Parent1", "Parent2", "F1count", "F1mean",
                "select_type", "cycle", "rep", "sel_scheme",
                "D2Ideal", "coancestry", "ideal_phenotype", "pair_score")
SUMMARY_COLS <- c("cycle", "select_type", "rep", "sel_scheme", "mean_F1mean",
                  "var_F1mean", "max_F1mean", "mean_D2Ideal", "mean_coancestry",
                  "mean_gene_diversity", "mean_ideal_phenotype", "ideal_hap_gebv",
                  "marker_size", "n_families", "n_parents_used")

# ── output contract ───────────────────────────────────────────────────────────

test_that("max_F1mean is the cycle's best single offspring (0.7.5)", {
  fx  <- .driver_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, scheme = "TS", replications = 1L, select_cycles = 2L,
    top_founders = 4L, num_offspring = 3L, cross_backend = "block",
    metrics = "gebv", seed = 7L, output_dir = NULL, verbose = FALSE)
  s <- res$summary
  expect_true("max_F1mean" %in% names(s))
  expect_true(all(is.finite(s$max_F1mean)))
  expect_true(all(s$max_F1mean >= s$mean_F1mean))
  # founder row: the best founder's total local GEBV
  lg <- compute_profile_localgebv(fx$meta, fx$eff, level = "haplotype")
  expect_equal(s[s$cycle == 0L, ]$max_F1mean,
               max(haploStack:::.stk_total_gebv(lg)), tolerance = 1e-8)
  # a cycle row: the family means in `plans` cannot exceed the best offspring
  p <- res$plans
  for (cyc in unique(s$cycle[s$cycle > 0L]))
    expect_true(s[s$cycle == cyc, ]$max_F1mean >=
                  max(p[p$cycle == cyc, ]$F1mean) - 1e-8)
})

test_that("a single-scheme run returns the origin output contract", {
  fx  <- .driver_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
    top_founders = 4L, num_offspring = 2L, scheme = "TS", metrics = "all",
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)

  # 0.9.0: the five 0.8.8 elements plus notes / provenance / figures (hs_result)
  expect_setequal(names(res),
                  c("plans", "summary", "ideal_genotypes",
                    "final_populations", "trajectories",
                    "notes", "provenance", "figures"))
  expect_s3_class(res, "hs_simulation")
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
    # OCS on this toy stops at the kinship cap and falls back to truncation;
    # since 0.9.0 that is announced (hs_warning_fallback) and noted
    res <- suppressWarnings(run_stacking_simulation(
      fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
      top_founders = 4L, num_offspring = 2L, scheme = sc,
      ohs_method = "Greedy", select_top_blocks = 50L,
      output_dir = NULL, verbose = FALSE), classes = "hs_warning_fallback")
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
  # 0.9.0: catalogue_results / trait_name / select_method are on the deprecation ladder (they warn)
  res <- suppressWarnings(run_stacking_simulation(
    catalogue_results = cat_res, trait_name = "SimTrait", effect_df = fx$eff,
    replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, select_method = "GA", output_dir = NULL, verbose = FALSE),
    classes = "hs_warning_deprecated")

  # GA -> OHS(GA_pairwise) since 0.7.4 (GA_coverage removed); trajectories
  # compat view is still emitted.
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
