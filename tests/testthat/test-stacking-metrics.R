# =============================================================================
# test-stacking-metrics.R
#
# QA for the optional-metrics framework in run_stacking_simulation() (Phase M1):
#   - default metrics = "gebv": only F1mean populated; all other metric columns NA;
#     ideal_genotypes NULL (not built)
#   - metrics are measurement-only: identical selection/offspring for any `metrics`
#   - coancestry is measured on the OFFSPRING population
#   - .stk_ideal_hap_gebv: per-block top-ploidy ceiling
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

.metrics_fixture <- function(seed = 5L) {
  m <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = seed)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))
  list(meta = m, eff = e)
}

run_sim <- function(fx, metrics = "gebv", scheme = "TS", ...) {
  run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = scheme, select_top_blocks = 50L,
    metrics = metrics, output_dir = NULL, verbose = FALSE, ...)
}

# ── default is lean ───────────────────────────────────────────────────────────

test_that("default metrics = 'gebv' populates only F1mean; others NA", {
  res <- run_sim(.metrics_fixture())
  s   <- res$summary
  expect_false(anyNA(s$mean_F1mean))                    # F1mean always on
  for (col in c("mean_D2Ideal", "mean_coancestry", "mean_gene_diversity",
                "mean_ideal_phenotype", "ideal_hap_gebv"))
    expect_true(all(is.na(s[[col]])), info = col)
  expect_true(all(is.na(res$plans$D2Ideal)))
  expect_true(all(is.na(res$plans$coancestry)))
  expect_null(res$ideal_genotypes)                      # not built for lean TS
})

test_that("requested metrics are populated (offspring cycles)", {
  res <- run_sim(.metrics_fixture(), metrics = "all")
  s   <- res$summary[cycle > 0L]
  for (col in c("mean_D2Ideal", "mean_gene_diversity", "mean_coancestry",
                "ideal_hap_gebv"))
    expect_false(all(is.na(s[[col]])), info = col)
  expect_false(is.null(res$ideal_genotypes))
  # offspring coancestry is a mean block-IBS in [0, 1]
  cv <- s$mean_coancestry[!is.na(s$mean_coancestry)]
  expect_true(all(cv >= -1e-9 & cv <= 1 + 1e-9))
})

# ── metrics never change selection (the key guarantee) ────────────────────────

test_that("selection & offspring are identical for any metrics (fixed seed)", {
  for (sc in c("TS", "GBS")) {
    a <- run_sim(.metrics_fixture(), metrics = "gebv", scheme = sc)
    b <- run_sim(.metrics_fixture(), metrics = "all",  scheme = sc)
    expect_identical(a$plans[, .(CrossID, Parent1, Parent2)],
                     b$plans[, .(CrossID, Parent1, Parent2)], info = sc)
    expect_identical(dimnames(a$final_populations[[1L]]$geno_arrays)[[1L]],
                     dimnames(b$final_populations[[1L]]$geno_arrays)[[1L]], info = sc)
    # F1mean (the crossing outcome) is unaffected too
    expect_equal(a$summary$mean_F1mean, b$summary$mean_F1mean, info = sc)
  }
})

# ── ideal_hap_gebv (top-ploidy per block) ─────────────────────────────────────

test_that(".stk_ideal_hap_gebv sums the top-ploidy segments per block", {
  lg <- data.table(Block = c("B1", "B2"),
                   S1.1 = c(1, 2), S1.2 = c(3, 0),
                   S2.1 = c(5, 1), S2.2 = c(0, 4))
  # B1 top2 of {1,3,5,0} = 5+3 = 8 ; B2 top2 of {2,0,1,4} = 4+2 = 6 ; total 14
  expect_equal(haploStack:::.stk_ideal_hap_gebv(lg, ploidy = 2L), 14)
  # ploidy 1 -> per-block max: 5 + 4 = 9
  expect_equal(haploStack:::.stk_ideal_hap_gebv(lg, ploidy = 1L), 9)
})

test_that("metrics = 'all' expands and unknown tokens error", {
  fx <- .metrics_fixture()
  expect_error(run_sim(fx, metrics = "not_a_metric"), "arg")
})
