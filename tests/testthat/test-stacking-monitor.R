# =============================================================================
# test-stacking-monitor.R
#
# FR5 metrics + monitor helpers (R/stacking_monitor.R, PRD #013 issue 05 /
# Phase 5): block-IBS coancestry (compute_parent_coancestry) and the
# summary-based metrics plot (plot_simulation_metrics).
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

# ── compute_parent_coancestry() ──────────────────────────────────────────────

test_that("compute_parent_coancestry returns a symmetric matrix in [0, 1]", {
  m <- sim_geno_meta(n_samples = 8L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = 3L)
  K <- compute_parent_coancestry(m)

  expect_true(is.matrix(K))
  expect_equal(dim(K), c(8L, 8L))
  expect_equal(K, t(K))                              # symmetric
  expect_true(all(K >= 0 & K <= 1, na.rm = TRUE))
  expect_identical(rownames(K), dimnames(m$geno_arrays)[[1L]])
})

test_that("self-coancestry (diagonal) is >= mean off-diagonal", {
  m <- sim_geno_meta(n_samples = 8L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = 7L)
  K <- compute_parent_coancestry(m)
  # each individual shares its own two haplotypes (two identical pairings) so a
  # self value cannot be below the population's average relatedness.
  expect_gte(mean(diag(K)), mean(K[upper.tri(K)]))
})

test_that("compute_parent_coancestry restricts to requested samples and blocks", {
  m  <- sim_geno_meta(n_samples = 10L, n_markers = 24L, n_chrom = 2L,
                      block_size = 3L, seed = 11L)
  ss <- dimnames(m$geno_arrays)[[1L]][1:4]
  K  <- compute_parent_coancestry(m, samples = ss)
  expect_identical(rownames(K), ss)

  blk <- unique(as.character(m$block_marker_map_info$Block))[1L]
  Kb  <- compute_parent_coancestry(m, samples = ss, blocks = blk)
  expect_equal(dim(Kb), c(4L, 4L))                   # still valid, one-block IBS
})

# ── compute_gene_diversity() (Nei 1973) ──────────────────────────────────────

test_that(".nei_locus_h matches hand-computed Nei diversity", {
  expect_equal(haploStack:::.nei_locus_h(c("A", "A", "T", "T")), 0.5)
  expect_equal(haploStack:::.nei_locus_h(c("A", "A", "A")), 0)          # monomorphic
  expect_equal(haploStack:::.nei_locus_h(c("A", "T", "G")), 2 / 3)      # 3 equal alleles
  expect_equal(haploStack:::.nei_locus_h(c("A", "A", NA, "")), 0)       # NA/"" dropped
  expect_true(is.na(haploStack:::.nei_locus_h(c(NA, "", NA))))          # nothing comparable
})

test_that("compute_gene_diversity returns H in [0,1] with per-locus attribute", {
  m  <- sim_geno_meta(n_samples = 20L, n_markers = 30L, n_chrom = 2L,
                      block_size = 3L, seed = 8L)
  Hm <- compute_gene_diversity(m, level = "marker")
  expect_true(is.numeric(Hm) && length(Hm) == 1L)
  expect_gte(as.numeric(Hm), 0); expect_lte(as.numeric(Hm), 1)
  expect_length(attr(Hm, "per_locus"), 30L)                 # one per marker

  Hh <- compute_gene_diversity(m, level = "haplotype")
  expect_length(attr(Hh, "per_locus"), 10L)                 # one per block
  # block-haplotypes are multi-allelic -> at least as diverse as single markers
  expect_gte(as.numeric(Hh), as.numeric(Hm))
})

test_that("compute_gene_diversity is 0 for a monomorphic population", {
  m <- sim_geno_meta(n_samples = 10L, n_markers = 12L, block_size = 3L, seed = 2L)
  m$geno_arrays[, , ] <- "1"
  expect_equal(as.numeric(compute_gene_diversity(m, level = "marker")), 0)
  expect_equal(as.numeric(compute_gene_diversity(m, level = "haplotype")), 0)
})

# ── coancestry wired into the driver output ──────────────────────────────────

test_that("driver populates a finite block-IBS coancestry column", {
  m   <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                       block_size = 3L, seed = 5L)
  e   <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))
  res <- run_stacking_simulation(
    m, e, replications = 1L, select_cycles = 3L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
    metrics = "coancestry", output_dir = NULL, verbose = FALSE)

  coan <- res$plans$coancestry
  expect_true(any(is.finite(coan)))
  expect_true(all(coan >= 0 & coan <= 1, na.rm = TRUE))
})

test_that("driver reports Nei gene diversity per cycle, in both crossing backends", {
  m <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = 5L)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))

  run_backend <- function(backend) {
    run_stacking_simulation(
      m, e, replications = 1L, select_cycles = 4L, top_founders = 4L,
      num_offspring = 3L, scheme = "GBS", select_top_blocks = 50L,
      cross_backend = backend, metrics = "gene_diversity",
      output_dir = NULL, verbose = FALSE)$summary
  }

  for (backend in c("block", "genomicSimulation")) {
    if (backend == "genomicSimulation") skip_if_not_installed("genomicSimulation")
    s <- run_backend(backend)
    expect_true("mean_gene_diversity" %in% names(s), info = backend)
    gd <- s[cycle > 0L]$mean_gene_diversity
    expect_true(all(is.finite(gd)), info = backend)
    expect_true(all(gd >= 0 & gd <= 1), info = backend)
    # diversity is spent under selection: last cycle <= first decision cycle
    expect_lte(gd[length(gd)], gd[1L] + 1e-8)
  }
})

# ── plot_simulation_metrics() ────────────────────────────────────────────────

test_that("plot_simulation_metrics builds a ggplot from the summary contract", {
  m   <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                       block_size = 3L, seed = 5L)
  e   <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))
  res <- run_stacking_simulation(
    m, e, replications = 2L, select_cycles = 3L, top_founders = 4L,
    num_offspring = 2L, scheme = c("TS", "GBS"), select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE)

  p <- plot_simulation_metrics(res$summary)
  expect_s3_class(p, "ggplot")

  # subset of metrics + founder inclusion are honoured
  p2 <- plot_simulation_metrics(res$summary, metrics = "mean_F1mean",
                                include_founder = TRUE)
  expect_s3_class(p2, "ggplot")

  expect_error(plot_simulation_metrics(res$summary[, .(cycle)]),
               "sel_scheme")
})

# ── live monitor (issue 08) ──────────────────────────────────────────────────

test_that("monitor = TRUE writes live files, emits per-job progress, and does not change results", {
  m <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = 5L)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))
  od <- tempfile("livemon")                       # driver creates it if absent

  base <- list(m, e, replications = 2L, select_cycles = 3L, top_founders = 4L,
               num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
               metrics = "gebv")

  # monitor on (verbose so the per-job progress lines are emitted)
  msgs <- capture_messages(
    res_live <- do.call(run_stacking_simulation,
                        c(base, list(output_dir = od, monitor = TRUE,
                                     verbose = TRUE))))

  # 2 reps x 1 scheme = 2 jobs -> exactly two per-job completion lines
  expect_equal(sum(grepl("done - mean_F1mean", msgs, fixed = TRUE)), 2L)
  expect_true(file.exists(file.path(od, "live_summary.txt")))
  expect_true(file.exists(file.path(od, "live_plans.txt")))
  expect_true(file.exists(file.path(od, "live_metrics.png")))
  expect_true(file.exists(file.path(od, "live_trajectories.png")))

  # observation-only: identical return for the same seed with monitor off
  res_off <- do.call(run_stacking_simulation,
                     c(base, list(output_dir = NULL, monitor = FALSE,
                                  verbose = FALSE)))
  expect_equal(res_live$summary, res_off$summary)
  expect_equal(res_live$plans,   res_off$plans)
})

test_that("monitor = TRUE without output_dir announces and skips (no error)", {
  m <- sim_geno_meta(n_samples = 10L, n_markers = 18L, n_chrom = 2L,
                     block_size = 3L, seed = 4L)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))
  expect_message(
    res <- run_stacking_simulation(
      m, e, replications = 1L, select_cycles = 2L, top_founders = 4L,
      num_offspring = 2L, scheme = "TS", select_top_blocks = 40L,
      output_dir = NULL, monitor = TRUE, verbose = FALSE),
    "needs 'output_dir'")
  expect_s3_class(res$summary, "data.table")
})

test_that("live monitor harvests as jobs resolve under parallel execution, same result", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  m <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = 5L)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))
  od <- tempfile("livemon_par")

  base <- list(m, e, replications = 2L, select_cycles = 2L, top_founders = 4L,
               num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
               verbose = FALSE)

  res_par <- do.call(run_stacking_simulation,
                     c(base, list(output_dir = od, monitor = TRUE, n_cores = 2L)))
  expect_true(file.exists(file.path(od, "live_metrics.png")))

  # completion order may differ, but results are placed by job index -> identical
  res_seq <- do.call(run_stacking_simulation,
                     c(base, list(output_dir = NULL, monitor = FALSE, n_cores = 1L)))
  expect_equal(res_par$summary, res_seq$summary)
})
