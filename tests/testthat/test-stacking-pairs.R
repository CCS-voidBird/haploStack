# =============================================================================
# test-stacking-pairs.R
#
# Direct pair selection — PLACEHOLDER (thread #031, Q23 i). What is tested is
# the seam, not a selector: the budget that ties a direct-pairing cycle to the
# parents-first traversal (pair_selection_budget(), .stk_selection_budget()),
# the driver's argument validation, the stop message of .stk_select_pairs(),
# and that the default parents-first path is untouched.
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

test_that("selection_unit = 'pairs' stops at the first cycle with the resolved budget", {
  fx <- .pairs_fixture()
  run <- function(budget) run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "GBS", selection_unit = "pairs",
    pair_budget = budget, select_top_blocks = 50L, output_dir = NULL,
    verbose = FALSE)
  # 4 parents -> 6 families under "families"; 2..6 families under "parents"
  expect_error(run("families"), "not implemented.*6 families per cycle")
  expect_error(run("parents"),  "not implemented.*4 distinct parents.*2\\.\\.6 families")
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

test_that("the default parents-first path is unchanged by the seam", {
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
