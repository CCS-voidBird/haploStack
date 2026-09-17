# =============================================================================
# test-stacking-pairs.R
#
# Direct pair selection — the seam (thread #031): the budget that ties a
# direct-pairing cycle to the parents-first traversal (pair_selection_budget(),
# .stk_selection_budget()), the driver's argument validation, the engine's
# stop rules (.pairs_greedy()), and that the default parents-first path is
# untouched. The scheme scorers are tested in test-stacking-pairs-schemes.R.
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

.pairs_fixture <- function(seed = 5L) {
  m <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = seed)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))
  list(meta = m, eff = e)
}

# ── family count of a traversal ───────────────────────────────────────────────

test_that("a traversal of n parents makes n(n-1)/2 families, + n with selfing", {
  expect_equal(haploStack:::.stk_n_families(16L), 120L)
  expect_equal(haploStack:::.stk_n_families(16L, allow_selfing = TRUE), 136L)
  expect_equal(haploStack:::.stk_n_families(2L), 1L)
  expect_equal(haploStack:::.stk_n_families(4L, TRUE), 10L)
  expect_type(haploStack:::.stk_n_families(16L), "integer")
})

test_that("the all-pairs design never contains a self-cross, whatever avoid_selfing", {
  # Recorded behaviour the budget relies on (INDEX discrepancy 27 / Q24 c):
  # combn() enumerates distinct pairs only, so the flag does not add selfs.
  p <- paste0("P", 1:5)
  d_avoid <- make_cross_design(p, method = "all_pairs", avoid_selfing = TRUE)
  d_allow <- make_cross_design(p, method = "all_pairs", avoid_selfing = FALSE)
  expect_equal(nrow(d_avoid), choose(5L, 2L))
  expect_equal(nrow(d_allow), choose(5L, 2L))
  expect_false(any(d_allow$Parent1 == d_allow$Parent2))
  expect_equal(nrow(d_allow), haploStack:::.stk_n_families(5L, allow_selfing = FALSE))
})

# ── the two conflicting budgets ───────────────────────────────────────────────

test_that("pair_budget = 'families' fixes the family count and frees the parents", {
  b <- pair_selection_budget(16, "families")
  expect_equal(b$pair_budget, "families")
  expect_equal(b$n_parents_ref, 16L)
  expect_false(b$allow_selfing)
  expect_equal(b$n_families, 120L)
  expect_equal(c(b$n_families_min, b$n_families_max), c(120L, 120L))
  expect_true(is.na(b$n_parents_max))

  bs <- pair_selection_budget(16, "families", allow_selfing = TRUE)
  expect_equal(bs$n_families, 136L)
})

test_that("pair_budget = 'parents' fixes the parent count and frees the families", {
  b <- pair_selection_budget(16, "parents")
  expect_equal(b$pair_budget, "parents")
  expect_equal(b$n_parents_max, 16L)
  expect_true(is.na(b$n_families))
  expect_equal(b$n_families_min, 8L)       # every pair brings two new parents
  expect_equal(b$n_families_max, 120L)     # the full traversal
  # odd n: ceiling
  expect_equal(pair_selection_budget(5, "parents")$n_families_min, 3L)
})

test_that("pair_selection_budget validates its inputs", {
  expect_error(pair_selection_budget(1, "families"), ">= 2")
  expect_error(pair_selection_budget(16, "crosses"))
  expect_error(pair_selection_budget(c(4, 8), "families"), "single integer")
})

test_that(".stk_selection_budget resolves parents-first to the traversal", {
  b <- haploStack:::.stk_selection_budget("parents", NULL, 6L)
  expect_equal(b$selection_unit, "parents")
  expect_true(is.na(b$pair_budget))
  expect_equal(b$n_families, 15L)
  expect_equal(b$n_families_max, 15L)
  expect_equal(b$n_parents_max, 6L)

  # pairs without a budget is refused up front, naming both options
  expect_error(haploStack:::.stk_selection_budget("pairs", NULL, 6L),
               "families.*parents")
  bp <- haploStack:::.stk_selection_budget("pairs", "families", 6L)
  expect_equal(bp$selection_unit, "pairs")
  expect_equal(bp$n_families, 15L)
})

test_that(".stk_budget_text states the numbers for each unit", {
  txt <- function(unit, budget) haploStack:::.stk_budget_text(
    haploStack:::.stk_selection_budget(unit, budget, 16L))
  expect_match(txt("parents", NULL), "16 parents -> 120 families")
  expect_match(txt("pairs", "families"), "120 families per cycle")
  expect_match(txt("pairs", "parents"), "16 distinct parents.*8\\.\\.120 families")
})

# ── driver: the placeholder ───────────────────────────────────────────────────

test_that("selection_unit = 'pairs' without pair_budget fails before running", {
  fx <- .pairs_fixture()
  expect_error(
    run_stacking_simulation(
      fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
      num_offspring = 2L, scheme = "TS", selection_unit = "pairs",
      select_top_blocks = 50L, output_dir = NULL, verbose = FALSE),
    "needs 'pair_budget'")
})

test_that("selection_unit = 'pairs' runs every cycle under either budget", {
  fx <- .pairs_fixture()
  run <- function(budget) run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 3L, top_founders = 4L,
    num_offspring = 2L, scheme = "GBS", selection_unit = "pairs",
    pair_budget = budget, select_top_blocks = 50L, output_dir = NULL,
    verbose = FALSE)
  # 4 parents -> exactly 6 families under "families"; 4 distinct parents and
  # 2..6 families under "parents"
  a <- run("families")
  expect_true(all(a$summary[cycle > 0L]$n_families == 6L))
  expect_true(all(a$plans[, .N, by = cycle]$N == 6L))
  b <- run("parents")
  expect_true(all(b$summary[cycle > 0L]$n_parents_used == 4L))
  expect_true(all(b$summary[cycle > 0L]$n_families %between% c(2L, 6L)))
  expect_true(all(!is.na(b$plans$pair_score)))
})

test_that("selection_unit = 'pairs' refuses what has no family form", {
  fx <- .pairs_fixture()
  run <- function(...) run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 1L, top_founders = 4L,
    num_offspring = 2L, selection_unit = "pairs", pair_budget = "families",
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE, ...)
  expect_error(run(scheme = "OCS"), "no direct-pair form")
  expect_error(run(scheme = "OHS", ohs_method = "GA_coverage"), "arg")   # removed 0.7.4
  expect_error(run(scheme = "GBS", gbs_engine = "GA", gbs_objective = "set"),
               "family sets")
})

# ── the engine: stop rules ────────────────────────────────────────────────────

# An additive scorer over a fixed score triangle (lower = better).
.toy_scorer <- function(S, allow_selfing = FALSE) {
  n <- nrow(S)
  cand0 <- haploStack:::.pair_candidates(n, allow_selfing)
  list(scheme = "TOY", samples = rownames(S), additive = TRUE, init = NULL,
       candidates = function(state) cand0,
       score = function(state, cand) S[cbind(cand$i, cand$j)],
       update = function(state, i, j) state)
}

test_that("the engine takes the top-F families under 'families' (= order())", {
  set.seed(1)
  n <- 8L
  S <- matrix(stats::runif(n * n), n, n, dimnames = list(paste0("P", 1:n), paste0("P", 1:n)))
  S <- (S + t(S)) / 2
  sc <- .toy_scorer(S)
  b  <- pair_selection_budget(4, "families")            # 6 families
  out <- haploStack:::.pairs_greedy(sc, b)
  expect_s3_class(out, "pair_selection")
  expect_equal(nrow(out), 6L)
  expect_equal(attr(out, "n_families"), 6L)
  # equals the 6 smallest upper-triangle entries, in order
  ut <- which(upper.tri(S), arr.ind = TRUE)
  ref <- sort(S[ut])[1:6]
  expect_equal(out$score, ref)
  expect_true(all(out$Parent1 != out$Parent2))
  key <- paste(pmin(out$Parent1, out$Parent2), pmax(out$Parent1, out$Parent2))
  expect_false(anyDuplicated(key) > 0)
})

test_that("the engine stops at n distinct parents under 'parents' and honours overshoot", {
  n <- 8L
  # scores chosen so the best-first order is (1,2) (3,4) (5,6) (1,3) (2,4) ...
  S <- matrix(10, n, n, dimnames = list(paste0("P", 1:n), paste0("P", 1:n)))
  S[1, 2] <- S[2, 1] <- 1; S[3, 4] <- S[4, 3] <- 2; S[5, 6] <- S[6, 5] <- 3
  S[1, 3] <- S[3, 1] <- 4; S[2, 4] <- S[4, 2] <- 5
  sc <- .toy_scorer(S)
  b  <- pair_selection_budget(4, "parents")            # 4 names, 2..6 families
  skip <- haploStack:::.pairs_greedy(sc, b, overshoot = "skip")
  # (1,2), (3,4) fill the 4 names; (5,6) is skipped; then families among 1..4
  expect_equal(attr(skip, "n_parents_used"), 4L)
  expect_equal(nrow(skip), 6L)                          # up to the traversal cap
  expect_setequal(unique(c(skip$Parent1, skip$Parent2)), paste0("P", 1:4))
  expect_equal(skip$score[1:2], c(1, 2))
  stop_ <- haploStack:::.pairs_greedy(sc, b, overshoot = "stop")
  expect_equal(nrow(stop_), 2L)                         # ends at the first overshoot
  expect_equal(attr(stop_, "n_parents_used"), 4L)
})

test_that("the engine honours a per-parent usage cap and the selfing rule", {
  n <- 6L
  S <- matrix(10, n, n, dimnames = list(paste0("P", 1:n), paste0("P", 1:n)))
  S[1, 2] <- 1; S[1, 3] <- 2; S[1, 4] <- 3; S[2, 3] <- 4
  S <- pmin(S, t(S))
  b <- pair_selection_budget(3, "families")             # 3 families
  out <- haploStack:::.pairs_greedy(.toy_scorer(S), b, max_use = 2L)
  # P1 may appear twice only: (1,2), (1,3), then (2,3) instead of (1,4)
  expect_equal(out$score, c(1, 2, 4))
  # selfs offered only when allowed
  diag(S) <- 0
  none <- haploStack:::.pairs_greedy(.toy_scorer(S, FALSE), b)
  expect_false(any(none$Parent1 == none$Parent2))
  bs   <- pair_selection_budget(3, "families", allow_selfing = TRUE)
  some <- haploStack:::.pairs_greedy(.toy_scorer(S, TRUE), bs)
  expect_true(any(some$Parent1 == some$Parent2))
})

test_that("pair_budget is ignored (with a warning) under the parents-first unit", {
  fx <- .pairs_fixture()
  expect_warning(
    res <- run_stacking_simulation(
      fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
      num_offspring = 2L, scheme = "TS", pair_budget = "families",
      select_top_blocks = 50L, output_dir = NULL, verbose = FALSE),
    "ignored")
  expect_s3_class(res$summary, "data.table")
})

test_that("the verbose header states the parents-first budget", {
  fx <- .pairs_fixture()
  expect_message(
    run_stacking_simulation(
      fx$meta, fx$eff, replications = 1L, select_cycles = 1L, top_founders = 4L,
      num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
      output_dir = NULL, verbose = TRUE),
    "4 parents -> 6 families per cycle")
})

test_that("the default parents-first path is unchanged by the pairs unit", {
  fx <- .pairs_fixture()
  run <- function(...) run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 3L, top_founders = 4L,
    num_offspring = 2L, scheme = c("TS", "GBS"), select_top_blocks = 50L,
    output_dir = NULL, verbose = FALSE, ...)
  a <- run()
  b <- run(selection_unit = "parents")
  expect_identical(a$summary, b$summary)
  expect_identical(a$plans, b$plans)
  # one family per distinct pair of the 4 parents, every cycle
  fam <- a$plans[cycle > 0L, .N, by = .(sel_scheme, cycle)]
  expect_true(all(fam$N == choose(4L, 2L)))
  expect_false(any(a$plans$Parent1 == a$plans$Parent2))
})
