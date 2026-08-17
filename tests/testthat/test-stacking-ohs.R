# =============================================================================
# test-stacking-ohs.R
#
# OHS pairwise objective + engines (R/stacking_ohs.R, PRD #013 issue 01).
#
# Fixture (2 blocks, 3 samples, ploidy 2) — complementary parents A,B vs the
# redundant A,C. Hand-computed pairwise matrix (diagonal = 0, off-diagonal only,
# matching the origin compute_pairwise.cpp):
#         A    B    C
#    A    0   20   19
#    B   20    0   19
#    C   19   19    0
# sum(M[S,S]) = 2*M[i,j]: {A,B}=40 (optimum), {A,C}=38, {B,C}=38.
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

.ohs_lg <- data.table(
  Block = c("B1", "B2"),
  A.1   = c(10, 0),  A.2 = c(0, 0),
  B.1   = c(0,  0),  B.2 = c(0, 10),
  C.1   = c(9,  0),  C.2 = c(0, 0)
)

.obj <- function(sel, M) sum(M[sel, sel])

# ── pairwise block-value matrix ───────────────────────────────────────────────

test_that(".pairwise_block_value_matrix is symmetric with the right values", {
  M <- haploStack:::.pairwise_block_value_matrix(.ohs_lg)
  expect_equal(dim(M), c(3L, 3L))
  expect_equal(rownames(M), c("A", "B", "C"))
  expect_true(isSymmetric(M))
  expect_equal(unname(diag(M)), c(0, 0, 0))       # self-pairs excluded
  expect_equal(M["A", "B"], 20)                   # complementary
  expect_equal(M["A", "C"], 19)                   # redundant at B1
  expect_equal(M["B", "C"], 19)
})

test_that(".pairwise_block_value_matrix infers ploidy from columns", {
  M <- haploStack:::.pairwise_block_value_matrix(.ohs_lg, ploidy = NULL)
  expect_equal(M["A", "B"], 20)                   # ploidy 2 inferred (top-2)
})

# ── OHS engines ───────────────────────────────────────────────────────────────

test_that("OHS Greedy picks the complementary optimum {A,B}", {
  sel <- ohs_selection(.ohs_lg, n_parents = 2L, method = "Greedy")
  expect_setequal(sel, c("A", "B"))
})

test_that("OHS GA_pairwise reaches the pairwise optimum (40)", {
  M   <- haploStack:::.pairwise_block_value_matrix(.ohs_lg)
  sel <- ohs_selection(.ohs_lg, n_parents = 2L, method = "GA_pairwise",
                       n_generations = 50L, seed = 1L)
  expect_length(sel, 2L)
  expect_equal(.obj(sel, M), 40)                  # global optimum
})

test_that("OHS GA_pairwise objective >= GA_coverage objective", {
  M    <- haploStack:::.pairwise_block_value_matrix(.ohs_lg)
  s_pw <- ohs_selection(.ohs_lg, n_parents = 2L, method = "GA_pairwise",
                        n_generations = 50L, seed = 1L)
  s_cv <- ohs_selection(.ohs_lg, n_parents = 2L, method = "GA_coverage",
                        n_generations = 50L, seed = 1L)
  expect_length(s_cv, 2L)
  expect_gte(.obj(s_pw, M), .obj(s_cv, M))
})

test_that("ohs_selection returns all samples when n_parents >= n", {
  sel <- ohs_selection(.ohs_lg, n_parents = 5L, method = "GA_pairwise")
  expect_setequal(sel, c("A", "B", "C"))
})

test_that("ohs_selection validates input", {
  expect_error(ohs_selection(data.table(x = 1), method = "Greedy"), "Block")
})

# ── haplotype -> genotype collapse ────────────────────────────────────────────

test_that(".collapse_hap_to_geno sums haplotype columns per sample", {
  g <- haploStack:::.collapse_hap_to_geno(.ohs_lg)
  expect_setequal(names(g), c("Block", "A", "B", "C"))
  expect_equal(g[Block == "B1", A], 10)    # 10 + 0
  expect_equal(g[Block == "B2", B], 10)    # 0 + 10
})
