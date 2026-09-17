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

test_that("the GA_coverage engine is gone (0.7.4)", {
  expect_error(ohs_selection(.ohs_lg, n_parents = 2L, method = "GA_coverage"),
               "arg")
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

# ── 0.7.3: elite subsetting, the greedy's first pick, the GA's RNG ───────────

# Random wide haplotype-GEBV table (ploidy 2), independent of the package.
.ohs_random_lg <- function(n_blocks, n_samples, seed = 1L) {
  set.seed(seed)
  cols <- list(Block = paste0("B", seq_len(n_blocks)))
  for (s in seq_len(n_samples)) for (h in 1:2)
    cols[[paste0("S", s, ".", h)]] <- round(rnorm(n_blocks, sd = 3), 3)
  as.data.table(cols)
}

# Oracles transcribed from the definitions, never from the package.
.oracle_pair_elite <- function(M, S, ratio) {
  m <- M[S, S]; v <- m[upper.tri(m)]
  q <- max(1L, ceiling(ratio * length(v) - 1e-9))
  sum(sort(v, decreasing = TRUE)[seq_len(q)])
}
.oracle_block_elite <- function(lg, S, ratio) {
  cols <- setdiff(names(lg), "Block"); base <- sub("[.][0-9]+$", "", cols)
  H <- as.matrix(lg[, ..cols])
  cb <- combn(S, 2L); P <- ncol(cb)
  q <- max(1L, ceiling(ratio * P - 1e-9))
  Mb <- sapply(seq_len(P), function(f)
    apply(H[, base %in% cb[, f], drop = FALSE], 1L, function(r)
      sum(sort(r, decreasing = TRUE)[1:2])))           # blocks x families
  sum(apply(matrix(Mb, ncol = P), 1L,
            function(r) sum(sort(r, decreasing = TRUE)[seq_len(q)])))
}
# HapSelect's Haploid_OHS: per block, the two best haplotype values in the set.
.oracle_haploid_ohs <- function(lg, S) {
  cols <- setdiff(names(lg), "Block"); base <- sub("[.][0-9]+$", "", cols)
  H <- as.matrix(lg[, ..cols])[, base %in% S, drop = FALSE]
  sum(apply(H, 1L, function(r) sum(sort(r, decreasing = TRUE)[1:2])))
}

test_that("elite_ratio = 1 is the full objective and the pre-0.7.3 GA picks", {
  lg <- .ohs_random_lg(20L, 25L, seed = 3L)
  for (lvl in c("pair", "block")) {
    a <- ohs_selection(lg, n_parents = 5L, method = "GA_pairwise", n_generations = 40L)
    b <- ohs_selection(lg, n_parents = 5L, method = "GA_pairwise", n_generations = 40L,
                       elite_ratio = 1, elite_level = lvl)
    expect_identical(a, b, info = lvl)
  }
  # the objective at ratio 1 is sum(M[S, S]) for both levels
  M   <- haploStack:::.pairwise_block_value_matrix(lg)
  obj <- haploStack:::.ohs_objective(M, elite_ratio = 1)
  expect_equal(obj$fit(1:5), sum(M[1:5, 1:5]))
  expect_true(obj$closed_form)
})

test_that("the pair-level elite objective is the sum of the q largest family ceilings", {
  lg  <- .ohs_random_lg(15L, 12L, seed = 4L)
  M   <- haploStack:::.pairwise_block_value_matrix(lg)
  for (r in c(0.5, 0.2, 0.001)) {
    obj <- haploStack:::.ohs_objective(M, elite_ratio = r, elite_level = "pair")
    S <- c(2L, 5L, 7L, 9L, 11L)
    expect_equal(obj$fit(S), .oracle_pair_elite(M, S, r), tolerance = 1e-5, info = r)
  }
  # q counts: 10 families -> 0.3 gives 3, 0.001 gives 1, 1 gives 10
  expect_equal(haploStack:::.ohs_elite_q(10L, 0.3), 3L)
  expect_equal(haploStack:::.ohs_elite_q(10L, 0.001), 1L)
  expect_equal(haploStack:::.ohs_elite_q(10L, 1), 10L)
})

test_that("the block-level elite objective matches its transcription, and tends to Haploid_OHS", {
  lg <- .ohs_random_lg(15L, 12L, seed = 6L)
  pm <- haploStack:::.pairwise_block_value_matrix(lg, return_array = TRUE)
  S  <- c(1L, 4L, 6L, 8L, 10L, 12L)
  for (r in c(0.5, 0.2)) {
    obj <- haploStack:::.ohs_objective(pm$M, pm$A, elite_ratio = r, elite_level = "block")
    expect_equal(obj$fit(S), .oracle_block_elite(lg, rownames(pm$M)[S], r),
                 tolerance = 1e-5, info = r)
  }
  obj1 <- haploStack:::.ohs_objective(pm$M, pm$A, elite_ratio = 0.001, elite_level = "block")
  expect_equal(obj1$fit(S), .oracle_haploid_ohs(lg, rownames(pm$M)[S]), tolerance = 1e-5)
  # a ratio just below 1 at the block level is the full family sum
  obj_full <- haploStack:::.ohs_objective(pm$M, pm$A, elite_ratio = 0.999999, elite_level = "block")
  expect_equal(obj_full$fit(S), sum(pm$M[S, S]) / 2, tolerance = 1e-5)
})

test_that("a small pair-level ratio selects a set containing the best single cross", {
  lg <- .ohs_random_lg(20L, 20L, seed = 8L)
  M  <- haploStack:::.pairwise_block_value_matrix(lg)
  best <- which(M == max(M), arr.ind = TRUE)[1L, ]
  for (m in c("GA_pairwise", "Greedy")) {
    sel <- ohs_selection(lg, n_parents = 4L, method = m, elite_ratio = 0.001)
    expect_true(all(rownames(M)[best] %in% sel), info = m)
  }
})

test_that("the greedy's first pick is the best row sum and column order does not matter", {
  lg <- .ohs_random_lg(20L, 20L, seed = 9L)
  M  <- haploStack:::.pairwise_block_value_matrix(lg)
  g  <- ohs_selection(lg, n_parents = 5L, method = "Greedy")
  expect_identical(g[1L], names(which.max(rowSums(M))))
  cols <- setdiff(names(lg), "Block")
  lg2 <- lg[, c("Block", rev(cols)), with = FALSE]
  g2  <- ohs_selection(lg2, n_parents = 5L, method = "Greedy")
  expect_setequal(g, g2)
  # never worse than the pre-0.7.3 greedy, which started from the first column
  sel <- 1L; rem <- setdiff(seq_len(nrow(M)), sel)
  for (k in 2:5) {
    mg <- 2 * colSums(M[sel, rem, drop = FALSE]); b <- rem[which.max(mg)]
    sel <- c(sel, b); rem <- setdiff(rem, b)
  }
  expect_gte(sum(M[g, g]), sum(M[sel, sel]))
})

test_that("GA_pairwise leaves the global RNG stream where it found it", {
  lg <- .ohs_random_lg(10L, 12L, seed = 2L)
  set.seed(123); before <- runif(1); set.seed(123)
  invisible(ohs_selection(lg, n_parents = 3L, method = "GA_pairwise", n_generations = 5L))
  expect_equal(runif(1), before)
})

test_that("elite arguments are validated", {
  lg <- .ohs_random_lg(5L, 6L, seed = 1L)
  expect_error(ohs_selection(lg, n_parents = 2L, elite_ratio = 0), "elite_ratio")
  expect_error(ohs_selection(lg, n_parents = 2L, elite_ratio = 1.5), "elite_ratio")
})
