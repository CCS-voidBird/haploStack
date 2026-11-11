# =============================================================================
# test-gbs-rank-scale.R
#
# GBS 0.7.11 (dev/UPGRADE_PLAN_gbs_multi_subtype.md stage 1, thread #036):
# `penalty_rank_scale` — a within-block rank-indexed multiplier of the
# absent-state penalty beta_t only. Gates T1–T7 of the plan:
#   T1 default parity (the bit-identical gate; the 0.7.10 behaviour itself is
#      pinned by the formula oracles in test-stacking-gbs-modes.R)
#   T2 single-target blocks invariant under any vector
#   T3 c(1,1,1) == c(1) exactly
#   T4 rank assignment (ties Q26.a, recycling Q26.b)
#   T5 composition with the dynamic penalty mode (pool-extinct minor state)
#   T6 pairs parity (gbs_pair_selection, both levels)
#   T7 the phantom-emergency fix, asserted on need scores directly
# plus Q26.c (message, not error, when [1] != 1) and argument validation.
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

# ── fixtures ──────────────────────────────────────────────────────────────────

# Random wide haplotype-GEBV table (the test-stacking-gbs-modes.R generator).
.rs_random_lg <- function(n_blocks, n_samples, n_states = 4L, seed = 1L) {
  set.seed(seed)
  states <- lapply(seq_len(n_blocks), function(b) round(sort(rnorm(n_states, sd = 3)), 2))
  cols <- list(Block = paste0("B", seq_len(n_blocks)))
  for (s in seq_len(n_samples)) for (h in 1:2)
    cols[[paste0("S", s, ".", h)]] <-
      vapply(seq_len(n_blocks), function(b) sample(states[[b]], 1L), numeric(1L))
  as.data.table(cols)
}

.rs_driver_fixture <- function(seed = 5L) {
  m <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = seed)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("A", "T"))
  list(meta = m, eff = e)
}

# ── T1 — default parity across the whole option grid ─────────────────────────

test_that("T1: penalty_rank_scale default/1 leaves every selection unchanged", {
  lg <- .rs_random_lg(25L, 20L, seed = 11L)
  ig <- construct_ideal_genotype(lg, method = "custom")
  for (m in c("ARCSIN", "LOG"))
    for (o in c("candidate", "set"))
      for (e in if (o == "set") c("greedy", "GA") else "greedy")
        for (pm in c("dynamic", "fixed"))
          for (pr in c("haploStack", "kemper")) {
            info <- paste(m, o, e, pm, pr)
            base <- gbs_selection(lg, ig, n_parents = 5L, gbs_method = m,
                                  gbs_objective = o, gbs_engine = e,
                                  penalty_mode = pm, gbs_preset = pr,
                                  ga_control = list(n_generations = 10L),
                                  seed = 42L)
            on1  <- gbs_selection(lg, ig, n_parents = 5L, gbs_method = m,
                                  gbs_objective = o, gbs_engine = e,
                                  penalty_mode = pm, gbs_preset = pr,
                                  penalty_rank_scale = 1,
                                  ga_control = list(n_generations = 10L),
                                  seed = 42L)
            expect_identical(base, on1, info = info)
          }
})

test_that("T1: the driver's fixed-seed run is unchanged by the default vector", {
  fx <- .rs_driver_fixture()
  run <- function(...) run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
    top_founders = 4L, num_offspring = 2L, scheme = "GBS",
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE, ...)
  a <- run()
  b <- run(penalty_rank_scale = 1)
  expect_identical(a$plans, b$plans)
  expect_identical(a$summary, b$summary)
})

# ── T2 — single-target blocks are invariant under any vector ─────────────────

test_that("T2: with one target state per block, any vector is a no-op", {
  lg <- .rs_random_lg(20L, 15L, seed = 7L)
  ig <- construct_ideal_genotype(lg, method = "argmax")   # one state per block
  expect_true(all(as.data.table(ig)[TargetFreq > 0, .N, by = Block]$N == 1L))
  base <- gbs_selection(lg, ig, n_parents = 5L)
  for (v in list(c(1, 0.01), c(1, 0.5, 0.1), c(1, 0, 0)))
    expect_identical(gbs_selection(lg, ig, n_parents = 5L,
                                   penalty_rank_scale = v),
                     base, info = paste(v, collapse = ","))
  # and on the scores themselves: every target is rank 1, so rscale is all 1
  mm   <- haploStack:::.gbs_hapn_matrix(lg, ig)
  spec <- haploStack:::.gbs_apply_preset("haploStack")
  spec$penalty_rank_scale <- c(1, 0.01, 0.001)
  ctx  <- haploStack:::.gbs_context(mm, 5L, "ARCSIN", spec)
  expect_identical(ctx$rscale, rep(1, ncol(mm$HapN)))
})

# ── T3 — c(1,1,1) is exactly c(1) ─────────────────────────────────────────────

test_that("T3: penalty_rank_scale = c(1,1,1) equals c(1) exactly", {
  lg <- .rs_random_lg(25L, 20L, seed = 3L)
  ig <- construct_ideal_genotype(lg, method = "custom")
  for (m in c("ARCSIN", "LOG")) for (o in c("candidate", "set")) {
    s1 <- gbs_selection(lg, ig, n_parents = 5L, gbs_method = m,
                        gbs_objective = o, penalty_rank_scale = 1)
    s3 <- gbs_selection(lg, ig, n_parents = 5L, gbs_method = m,
                        gbs_objective = o, penalty_rank_scale = c(1, 1, 1))
    expect_identical(s1, s3, info = paste(m, o))
  }
  # exact equality of the scores, not just the picks
  mm   <- haploStack:::.gbs_hapn_matrix(lg, ig)
  spec <- haploStack:::.gbs_apply_preset("haploStack")
  ctx1 <- ctx3 <- NULL
  spec$penalty_rank_scale <- 1
  ctx1 <- haploStack:::.gbs_context(mm, 5L, "ARCSIN", spec)
  spec$penalty_rank_scale <- c(1, 1, 1)
  ctx3 <- haploStack:::.gbs_context(mm, 5L, "ARCSIN", spec)
  sf <- rep(0, ncol(mm$HapN))
  expect_identical(haploStack:::.gbs_need_score(ctx1, sf),
                   haploStack:::.gbs_need_score(ctx3, sf))
  expect_identical(haploStack:::.gbs_set_objective(ctx1, sf),
                   haploStack:::.gbs_set_objective(ctx3, sf))
})

# ── T4 — rank assignment: order, ties (Q26.a), recycling (Q26.b) ─────────────

test_that("T4: ranks follow TargetFreq desc, CFreq desc, Effect asc; recycling", {
  # Rows deliberately interleaved across blocks to prove the row mapping.
  targets <- data.table(
    Block      = c("B1", "B2", "B1", "B3", "B2", "B1", "B3"),
    Effect     = c(10,    5,    8,    2,    6,    7,    1),
    CFreq      = c(0.5,  0.2,  0.3,  0.4,  0.4,  0.2,  0.4),
    TargetFreq = c(0.6,  0.5,  0.3,  0.5,  0.5,  0.1,  0.5))
  # B1: plain TargetFreq order            -> ranks 1 (E10), 2 (E8), 3 (E7)
  # B2: TargetFreq tied, CFreq breaks it  -> E6 (CFreq .4) rank 1, E5 rank 2
  # B3: TargetFreq and CFreq tied, Effect ascending -> E1 rank 1, E2 rank 2
  sc <- haploStack:::.gbs_rank_scale(targets, c(1, 0.1, 0.01))
  expect_equal(sc, c(1,     # B1 E10  r1
                     0.1,   # B2 E5   r2
                     0.1,   # B1 E8   r2
                     0.1,   # B3 E2   r2 (key ascending: E1 before E2)
                     1,     # B2 E6   r1 (CFreq desc)
                     0.01,  # B1 E7   r3
                     1))    # B3 E1   r1
  # recycling: ranks beyond the vector take its last element
  expect_equal(haploStack:::.gbs_rank_scale(targets, c(1, 0.1)),
               c(1, 0.1, 0.1, 0.1, 1, 0.1, 1))
  # a length-1 vector applies to every rank (and the default is all 1)
  expect_equal(haploStack:::.gbs_rank_scale(targets, 1), rep(1, 7L))
  expect_equal(haploStack:::.gbs_rank_scale(targets, 0.5), rep(0.5, 7L))
  # a table without CFreq still ranks (tie-break falls through to Effect)
  expect_equal(haploStack:::.gbs_rank_scale(targets[, !"CFreq"], c(1, 0.1, 0.01)),
               c(1, 1, 0.1, 0.1, 0.1, 0.01, 1))
  # (B2: TargetFreq tied, no CFreq -> key ascending -> E5 r1, E6 r2;
  #  B1 and B3 unchanged)
})

# ── T5 — composition with the dynamic penalty mode ────────────────────────────

test_that("T5: a pool-extinct minor state takes -|gbs_penalty| * 1 * scale[r]", {
  # Block B1: rank-1 state (E10) carried by A/B; rank-2 state (E7) in the ideal
  # but carried by nobody -> pop_f = 0 -> risk = 1, hope = 0 -> s_t = 1 (the
  # lost-cause damping), so beta = -10, and the rank multiplier sits on top.
  lg <- data.table(Block = "B1",
                   A.1 = 10, A.2 = 10, B.1 = 10, B.2 = 0)
  ideal <- data.table(Block = c("B1", "B1"), Effect = c(10, 7),
                      CFreq = c(0.75, 0), TargetFreq = c(0.6, 0.4))
  mm   <- haploStack:::.gbs_hapn_matrix(lg, ideal)
  spec <- haploStack:::.gbs_apply_preset("haploStack")   # dynamic, -10, cap 10
  spec$penalty_rank_scale <- c(1, 0.01)
  ctx  <- haploStack:::.gbs_context(mm, 2L, "ARCSIN", spec)
  t_extinct <- which(mm$targets$Effect == 7)
  expect_equal(mm$pop_f[t_extinct], 0)
  expect_equal(ctx$beta[t_extinct], -10)                 # s_t = 1 at hope 0
  expect_equal(ctx$rscale[t_extinct], 0.01)              # rank 2
  expect_equal((ctx$beta * ctx$rscale)[t_extinct], -0.1)
  # and that is the value the need score uses while the state is absent
  sf <- rep(0, ncol(mm$HapN))
  D  <- haploStack:::.gbs_need_score(ctx, sf)
  expect_equal(D, as.vector(mm$HapN %*% (ctx$beta * ctx$rscale * ctx$w)))
})

# ── T6 — pairs parity and the scaled beta in the pair scorer ─────────────────

test_that("T6: gbs_pair_selection default parity at both levels, vector runs", {
  lg <- .rs_random_lg(20L, 14L, seed = 9L)
  ig <- construct_ideal_genotype(lg, method = "custom")
  for (lv in c("population", "family")) {
    base <- gbs_pair_selection(lg, ig, pair_budget = "families",
                               top_founders = 4L, gbs_pair_level = lv)
    on1  <- gbs_pair_selection(lg, ig, pair_budget = "families",
                               top_founders = 4L, gbs_pair_level = lv,
                               penalty_rank_scale = 1)
    expect_identical(base, on1, info = lv)
    vec  <- gbs_pair_selection(lg, ig, pair_budget = "families",
                               top_founders = 4L, gbs_pair_level = lv,
                               penalty_rank_scale = c(1, 0.01, 0.001))
    expect_s3_class(vec, "data.table")
    expect_gt(nrow(vec), 0L)
    expect_true(all(c("Parent1", "Parent2", "score", "rank") %in% names(vec)))
  }
})

# ── T7 — the phantom-emergency fix, on need scores directly ──────────────────

test_that("T7: a rank-1-held/rank-2-absent block no longer outranks a blackout", {
  # Block BP ("phantom"): rank-1 state E10 held at target, rank-2 state E9
  # absent. Block BT ("true blackout"): its only rank-1 state E8 absent.
  # Candidate X carries only BP's rank-2 state; candidate Y only BT's rank-1
  # state, one copy each. Fixed mode makes every absent beta exactly -100.
  lg <- data.table(Block = c("BP", "BT"),
                   H1.1 = c(10, 0), H1.2 = c(10, 0),   # holds BP rank-1
                   X.1  = c(9, 0),  X.2  = c(0, 0),    # carries BP rank-2
                   Y.1  = c(0, 8),  Y.2  = c(0, 0))    # carries BT rank-1
  ideal <- data.table(Block = c("BP", "BP", "BT"),
                      Effect = c(10, 9, 8),
                      CFreq = c(1/3, 1/6, 1/6),
                      TargetFreq = c(0.5, 0.5, 1))     # BP tie -> CFreq desc
  mm   <- haploStack:::.gbs_hapn_matrix(lg, ideal)
  spec <- haploStack:::.gbs_apply_preset("kemper")      # fixed -10 x 10
  n_parents <- 2L
  # SelFreq state: BP rank-1 present (H1 selected), both other targets absent.
  sf <- colSums(mm$HapN[match("H1", mm$samples), , drop = FALSE]) /
        (n_parents * mm$ploidy)
  iX <- match("X", mm$samples); iY <- match("Y", mm$samples)

  spec$penalty_rank_scale <- 1                          # current regime
  ctx0 <- haploStack:::.gbs_context(mm, n_parents, "ARCSIN", spec)
  s0 <- haploStack:::.gbs_need_score(ctx0, sf)
  expect_equal(s0[iX], s0[iY])                          # the phantom emergency:
                                                        # both absences shout -100

  spec$penalty_rank_scale <- c(1, 0.01)                 # the fix
  ctx1 <- haploStack:::.gbs_context(mm, n_parents, "ARCSIN", spec)
  s1 <- haploStack:::.gbs_need_score(ctx1, sf)
  expect_equal(s1[iY], -100)                            # blackout: full emergency
  expect_equal(s1[iX], -1)                              # rank-2 absence: demoted
  expect_lt(s1[iY], s1[iX])                             # Y is now picked first
})

# ── Q26.c message and argument validation ────────────────────────────────────

test_that("penalty_rank_scale is validated; [1] != 1 messages, never errors", {
  lg <- .rs_random_lg(8L, 8L, seed = 2L)
  ig <- construct_ideal_genotype(lg, method = "custom")
  expect_error(gbs_selection(lg, ig, n_parents = 3L,
                             penalty_rank_scale = numeric(0)), "numeric vector")
  expect_error(gbs_selection(lg, ig, n_parents = 3L,
                             penalty_rank_scale = c(1, -0.1)), ">= 0")
  expect_error(gbs_selection(lg, ig, n_parents = 3L,
                             penalty_rank_scale = c(1, NA)), "finite")
  expect_error(gbs_selection(lg, ig, n_parents = 3L,
                             penalty_rank_scale = "a"), "numeric vector")
  expect_message(gbs_selection(lg, ig, n_parents = 3L,
                               penalty_rank_scale = c(0.5, 0.1)), "rescaled")
  expect_silent(sel <- gbs_selection(lg, ig, n_parents = 3L,
                                     penalty_rank_scale = c(1, 0.1)))
  expect_length(sel, 3L)
  expect_error(gbs_pair_selection(lg, ig, pair_budget = "families",
                                  top_founders = 3L,
                                  penalty_rank_scale = c(1, -1)), ">= 0")
})

test_that("the driver validates the vector up front and threads it through", {
  fx <- .rs_driver_fixture()
  expect_error(run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 1L,
    top_founders = 4L, num_offspring = 2L, scheme = "GBS",
    penalty_rank_scale = c(1, -1),
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE), ">= 0")
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
    top_founders = 4L, num_offspring = 2L, scheme = "GBS",
    penalty_rank_scale = c(1, 0.01, 0.001),
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)
  expect_gt(nrow(res$plans), 0L)
  # pairs route
  res2 <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
    top_founders = 4L, num_offspring = 2L, scheme = "GBS",
    selection_unit = "pairs", pair_budget = "families",
    penalty_rank_scale = c(1, 0.01),
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)
  expect_gt(nrow(res2$plans), 0L)
})
