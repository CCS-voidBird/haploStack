# =============================================================================
# test-stacking-gbs.R
#
# GBS scheme: ideal genotype + per-candidate selection + penalty tuning
# (R/stacking_gbs.R, PRD #013 issue 02).
#
# Fixtures are hand-built wide haplotype-level local GEBV tables (Block +
# "<Sample>.<idx>" columns, ploidy 2), as produced by
# compute_profile_localgebv(level = "haplotype").
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

# Construct-ideal fixture: 4 samples, 2 blocks, distinct effect-states.
.ideal_lg <- data.table(
  Block = c("B1", "B2"),
  A.1 = c(2, 0), A.2 = c(1, 0),
  B.1 = c(2, 5), B.2 = c(0, 5),
  C.1 = c(1, 5), C.2 = c(1, 0),
  D.1 = c(0, 0), D.2 = c(0, 0)
)

# Selection fixture: P1 carries the B1 target, P2 the B2 target, F* are fillers.
.sel_lg <- data.table(
  Block = c("B1", "B2"),
  P1.1 = c(10, 0), P1.2 = c(10, 0),
  P2.1 = c(0, 10), P2.2 = c(0, 10),
  F1.1 = c(0, 0),  F1.2 = c(0, 0),
  F2.1 = c(0, 0),  F2.2 = c(0, 0),
  F3.1 = c(0, 0),  F3.2 = c(0, 0)
)

# Rare-target fixture: P1b duplicates P1's B1 coverage; R is the lone B3 carrier.
.rare_lg <- data.table(
  Block = c("B1", "B2", "B3"),
  P1.1  = c(10, 0, 0), P1.2  = c(10, 0, 0),
  P1b.1 = c(10, 0, 0), P1b.2 = c(10, 0, 0),
  P2.1  = c(0, 10, 0), P2.2  = c(0, 10, 0),
  R.1   = c(0, 0, 7),  R.2   = c(0, 0, 0)
)

# No-op-guard fixture: X carries state 5, Y state 2, Z nothing (one block).
.ng_lg <- data.table(
  Block = c("B1"),
  X.1 = c(5), X.2 = c(5),
  Y.1 = c(2), Y.2 = c(2),
  Z.1 = c(0), Z.2 = c(0)
)

# ── ideal genotype construction ───────────────────────────────────────────────

test_that("construct_ideal_genotype returns the documented columns", {
  ig <- construct_ideal_genotype(.ideal_lg, method = "argmax")
  expect_equal(names(ig), c("Block", "Effect", "CFreq", "TargetFreq"))
})

test_that("argmax targets the top effect-state per block with TargetFreq 1", {
  ig <- construct_ideal_genotype(.ideal_lg, method = "argmax")
  expect_equal(nrow(ig), 2L)
  expect_true(all(ig$TargetFreq == 1))
  expect_equal(ig[Block == "B1", Effect], 2)   # max state at B1
  expect_equal(ig[Block == "B2", Effect], 5)   # max state at B2
  expect_equal(attr(ig, "ideal_phenotype"), 2 * 1 * 2 + 5 * 1 * 2)  # = 14
})

test_that("ratio TargetFreq sums to 1 within each block", {
  ig <- construct_ideal_genotype(.ideal_lg, method = "ratio")
  per_block <- ig[, .(s = sum(TargetFreq)), by = Block]
  expect_equal(per_block$s, rep(1, nrow(per_block)))
})

test_that("custom TargetFreq sums to 1 per block and ideal_phenotype is consistent", {
  ig <- construct_ideal_genotype(.ideal_lg, method = "custom", k_penalty = 0.5)
  per_block <- ig[, .(s = sum(TargetFreq)), by = Block]
  expect_equal(per_block$s, rep(1, nrow(per_block)))
  expect_equal(attr(ig, "ideal_phenotype"),
               sum(ig$Effect * ig$TargetFreq * 2))
})

# ── GBS selection ─────────────────────────────────────────────────────────────

test_that("gbs_selection picks the complementary target carriers over fillers", {
  sel <- gbs_selection(.sel_lg, n_parents = 2L)
  expect_setequal(sel, c("P1", "P2"))
})

test_that("gbs_selection returns n_parents names for both metrics", {
  expect_length(gbs_selection(.sel_lg, n_parents = 2L, gbs_method = "ARCSIN"), 2L)
  expect_length(gbs_selection(.sel_lg, n_parents = 2L, gbs_method = "LOG"), 2L)
})

test_that("gbs_selection returns all samples when n_parents >= n", {
  ig <- construct_ideal_genotype(.ng_lg, method = "argmax")
  expect_setequal(gbs_selection(.ng_lg, ig, n_parents = 5L), c("X", "Y", "Z"))
})

test_that("gbs_selection moves the selected set toward the ideal vs fillers", {
  ig    <- construct_ideal_genotype(.sel_lg)
  sel   <- gbs_selection(.sel_lg, ideal_genotypes = ig, n_parents = 2L)
  d_sel <- haploStack:::.distance_to_ideal(.sel_lg, ig, samples = sel)
  d_bad <- haploStack:::.distance_to_ideal(.sel_lg, ig, samples = c("F1", "F2"))
  expect_lt(d_sel, d_bad)
})

test_that("gbs_penalty protects a rare target haplotype", {
  # A large (auto) penalty recruits the lone B3 carrier R over the redundant P1b.
  sel_big <- gbs_selection(.rare_lg, n_parents = 3L, gbs_penalty = NULL)
  expect_true("R" %in% sel_big)

  # A negligible penalty over-represents the common B1 target and drops R.
  sel_tiny <- gbs_selection(.rare_lg, n_parents = 3L, gbs_penalty = -0.01)
  expect_false("R" %in% sel_tiny)
})

test_that("gbs_selection is not a no-op: selection tracks the ideal genotype", {
  ideal_hi <- data.table(Block = "B1", Effect = 5, CFreq = 2 / 6, TargetFreq = 1)
  ideal_lo <- data.table(Block = "B1", Effect = 2, CFreq = 2 / 6, TargetFreq = 1)
  expect_setequal(gbs_selection(.ng_lg, ideal_hi, n_parents = 1L), "X")
  expect_setequal(gbs_selection(.ng_lg, ideal_lo, n_parents = 1L), "Y")
})

# ── scale-relative penalty ────────────────────────────────────────────────────

test_that("auto penalty is scale-relative and a numeric override is respected", {
  auto2 <- haploStack:::.resolve_gbs_penalty(NULL, 2, "ARCSIN", c(1), 8L, 2L)
  auto4 <- haploStack:::.resolve_gbs_penalty(NULL, 4, "ARCSIN", c(1), 8L, 2L)
  expect_true(auto2 < 0)
  expect_equal(auto4, 2 * auto2)                         # linear in penalty_scale
  expect_equal(
    haploStack:::.resolve_gbs_penalty(-10, 2, "ARCSIN", c(1), 8L, 2L), -10)
})

# ── penalty tuning ────────────────────────────────────────────────────────────

test_that("tune_gbs_penalty returns a score table and a sane best penalty", {
  res <- tune_gbs_penalty(.sel_lg, n_parents = 2L, grid = c(1, 5),
                          n_rep = 2L, subsample_frac = 1)
  expect_named(res, c("best_penalty_scale", "best_penalty", "scores"))
  expect_setequal(names(res$scores), c("penalty_scale", "mean_score", "sd_score"))
  expect_true(res$best_penalty_scale %in% c(1, 5))
  expect_true(res$best_penalty < 0)
  expect_lte(min(res$scores$mean_score), max(res$scores$mean_score))
})
