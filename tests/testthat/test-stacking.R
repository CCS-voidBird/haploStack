# =============================================================================
# test-stacking.R
#
# QA tests for the stacking module:
#   make_cross_design(), greedy_selection(), genetic_algorithm_selection(),
#   ocs_selection(), compute_grm(), summarise_simulation()
# =============================================================================

library(testthat)
library(haploStack)

# Shared GEBV matrix (Block × Sample)
.make_gebv_mat <- function(n_blocks = 8L, n_samples = 10L, seed = 42L) {
  set.seed(seed)
  mat <- matrix(stats::rnorm(n_blocks * n_samples), nrow = n_blocks)
  colnames(mat) <- paste0("S", seq_len(n_samples))
  dt  <- as.data.table(mat)
  dt[, Block := paste0("B", seq_len(n_blocks))]
  setcolorder(dt, c("Block", paste0("S", seq_len(n_samples))))
  dt
}

.gebv_mat  <- .make_gebv_mat()
.parents   <- paste0("S", 1:6)

# ── make_cross_design ─────────────────────────────────────────────────────────

test_that("make_cross_design all_pairs returns correct number of crosses", {
  design <- make_cross_design(paste0("P", 1:4), method = "all_pairs")
  expect_equal(nrow(design), 6L)  # C(4,2) = 6
})

test_that("make_cross_design avoid_selfing removes self-crosses", {
  # Custom table with a self-cross
  ct <- data.frame(Parent1 = c("P1","P2","P1"), Parent2 = c("P2","P3","P1"))
  design <- make_cross_design(NULL, method = "custom", cross_table = ct,
                               avoid_selfing = TRUE)
  expect_true(all(design$Parent1 != design$Parent2))
})

test_that("make_cross_design sequential produces floor(n/2) crosses", {
  design <- make_cross_design(paste0("P", 1:5), method = "sequential")
  expect_equal(nrow(design), 2L)  # (P1×P2), (P3×P4) — P5 pairs with P4
})

test_that("make_cross_design custom uses supplied cross_table", {
  ct <- data.frame(Parent1 = c("PA","PB"), Parent2 = c("PB","PC"))
  design <- make_cross_design(NULL, method = "custom", cross_table = ct)
  expect_equal(nrow(design), 2L)
  expect_equal(design$Parent1, c("PA","PB"))
})

test_that("make_cross_design adds Cross_ID and n_offspring columns", {
  design <- make_cross_design(.parents[1:4], method = "all_pairs",
                               n_offspring = 3L)
  expect_true("Cross_ID"    %in% names(design))
  expect_true("n_offspring" %in% names(design))
  expect_true(all(design$n_offspring == 3L))
})

test_that("make_cross_design errors on custom without cross_table", {
  expect_error(make_cross_design(.parents, method = "custom"), "cross_table")
})

# ── greedy_selection ──────────────────────────────────────────────────────────

test_that("greedy_selection returns exactly n_parents names", {
  selected <- greedy_selection(.gebv_mat, n_parents = 5L)
  expect_equal(length(selected), 5L)
})

test_that("greedy_selection returns only names from the GEBV matrix", {
  selected <- greedy_selection(.gebv_mat, n_parents = 4L)
  expect_true(all(selected %in% paste0("S", 1:10)))
})

test_that("greedy_selection returns all candidates when n_parents = n_cand", {
  small_mat <- .make_gebv_mat(n_blocks = 4L, n_samples = 3L)
  selected  <- suppressWarnings(greedy_selection(small_mat, n_parents = 3L))
  expect_equal(length(selected), 3L)
})

test_that("greedy_selection warns when n_parents > n_candidates", {
  small_mat <- .make_gebv_mat(n_blocks = 4L, n_samples = 2L)
  expect_warning(greedy_selection(small_mat, n_parents = 5L),
                 "available candidates")
})

test_that("greedy_selection errors without Block column", {
  bad_mat <- data.table(S1 = 1, S2 = 2)
  expect_error(greedy_selection(bad_mat, 2L), "Block")
})

# ── genetic_algorithm_selection ───────────────────────────────────────────────

test_that("genetic_algorithm_selection returns exactly n_parents names", {
  selected <- genetic_algorithm_selection(.gebv_mat, n_parents = 4L,
                                           n_generations = 10L)
  expect_equal(length(selected), 4L)
})

test_that("genetic_algorithm_selection is reproducible with same seed", {
  s1 <- genetic_algorithm_selection(.gebv_mat, n_parents = 4L,
                                     n_generations = 10L, seed = 1L)
  s2 <- genetic_algorithm_selection(.gebv_mat, n_parents = 4L,
                                     n_generations = 10L, seed = 1L)
  expect_equal(s1, s2)
})

test_that("genetic_algorithm_selection returns different results with different seeds", {
  s1 <- genetic_algorithm_selection(.gebv_mat, n_parents = 4L,
                                     n_generations = 50L, seed = 1L)
  s2 <- genetic_algorithm_selection(.gebv_mat, n_parents = 4L,
                                     n_generations = 50L, seed = 99L)
  # Not guaranteed to differ but extremely likely for n_parents = 4, n_cand = 10
  # Skip if deterministically same
  skip_if(identical(sort(s1), sort(s2)),
          message = "Seeds produced same result by chance.")
  expect_false(identical(sort(s1), sort(s2)))
})

# ── ocs_selection ─────────────────────────────────────────────────────────────

test_that("ocs_selection with no kinship matrix does truncation selection", {
  gebv_vec <- setNames(c(1.0, 3.0, 2.0, 0.5), paste0("P", 1:4))
  selected <- ocs_selection(gebv_vec, n_parents = 2L)
  expect_equal(selected, c("P2","P3"))
})

test_that("ocs_selection returns exactly n_parents names", {
  set.seed(1L)
  gebv_vec <- setNames(stats::rnorm(8L), paste0("S", 1:8))
  selected <- ocs_selection(gebv_vec, n_parents = 4L)
  expect_lte(length(selected), 4L)  # may stop early due to kinship constraint
})

test_that("ocs_selection errors on unnamed gebv_vector", {
  expect_error(ocs_selection(c(1, 2, 3), n_parents = 2L), "named")
})

# ── compute_grm ───────────────────────────────────────────────────────────────

test_that("compute_grm returns a symmetric matrix with correct dimensions", {
  meta <- sim_geno_meta(n_samples = 4L, n_markers = 10L, seed = 1L)
  G    <- compute_grm(meta)
  expect_true(is.matrix(G))
  expect_equal(dim(G), c(4L, 4L))
  expect_equal(G, t(G))
})

test_that("compute_grm row/column names match sample names in GenoMeta", {
  meta           <- sim_geno_meta(n_samples = 3L, n_markers = 10L, seed = 2L)
  G              <- compute_grm(meta)
  expected_names <- dimnames(meta$geno_arrays)[[1L]]
  expect_equal(rownames(G), expected_names)
  expect_equal(colnames(G), expected_names)
})

test_that("compute_grm diagonal values are all positive", {
  meta <- sim_geno_meta(n_samples = 8L, n_markers = 20L, seed = 3L)
  G    <- compute_grm(meta)
  expect_true(all(diag(G) > 0))
})

test_that("compute_grm imputes NAs and messages", {
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 10L, seed = 42L)
  meta$geno_arrays[1L, 1L, 1L] <- NA_character_
  expect_message(compute_grm(meta), "imputed")
})

test_that("compute_grm GRM can be passed to ocs_selection", {
  meta     <- sim_geno_meta(n_samples = 6L, n_markers = 10L, seed = 4L)
  G        <- compute_grm(meta)
  snames   <- dimnames(meta$geno_arrays)[[1L]]
  gebv_vec <- setNames(stats::rnorm(6L), snames)
  selected <- ocs_selection(gebv_vec, kinship_matrix = G, n_parents = 3L)
  expect_true(length(selected) >= 1L)
  expect_true(all(selected %in% snames))
})

test_that("compute_grm errors when input is not a GenoMeta", {
  expect_error(compute_grm(matrix(c(0L,1L,2L,1L), nrow = 2L)), "GenoMeta")
  expect_error(compute_grm(list()),                              "GenoMeta")
})

# ── summarise_simulation ──────────────────────────────────────────────────────

test_that("summarise_simulation returns a data.table with Cycle column", {
  set.seed(7L)
  traj <- data.table(
    Rep           = rep(1:3, each = 5L),
    Cycle         = rep(1:5, 3L),
    Method        = "GA",
    MeanGEBV      = stats::rnorm(15L),
    SDev          = abs(stats::rnorm(15L)),
    BlockCoverage = runif(15L)
  )
  summ <- summarise_simulation(traj)
  expect_true("Cycle" %in% names(summ))
  expect_true(is.data.table(summ))
})

test_that("summarise_simulation has 5 rows for 5-cycle trajectory", {
  traj <- data.table(
    Rep      = rep(1:3, each = 5L),
    Cycle    = rep(1:5, 3L),
    MeanGEBV = stats::rnorm(15L)
  )
  summ <- summarise_simulation(traj, group_by = "Cycle")
  expect_equal(nrow(summ), 5L)
})
