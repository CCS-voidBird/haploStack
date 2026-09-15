# =============================================================================
# test-stacking-gbs-modes.R
#
# GBS 0.7.0 (dev/UPGRADE_PLAN_opt_kemper.md §4, Phase 1 gates a–h): Kemper's
# two modes as written (`gbs_preset = "kemper"`: ARCSIN / LOG, sequential
# greedy, fixed -100 = -10 x 10) and the haploStack adjustments (bounded
# risk x hope penalty scale, set objective, GA engine, effect weight), plus the
# driver's `ideal_scope`.
#
# Oracles are independent transcriptions of the formulas — Kemper's per-candidate
# score as the origin GRDC code wrote it (Value = HapN * dev * (1 + |Effect| w),
# absent -> penalty, ranked ascending) applied sequentially — never the package
# functions under test.
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

# ── fixtures ──────────────────────────────────────────────────────────────────

# Hand toy: 3 samples x 2 blocks, ploidy 2, three target states with distinct
# TargetFreq. Step-1 scores are -beta * (target copies): A 3, B 2, C 2 -> A.
.toy_lg <- data.table(
  Block = c("B1", "B2"),
  A.1 = c(10, 8), A.2 = c(5, 0),      # A carries t1, t2, t3 (one copy each)
  B.1 = c(10, 0), B.2 = c(10, 0),     # B carries t1 twice
  C.1 = c(5, 8),  C.2 = c(0, 0)       # C carries t2, t3
)
.toy_ideal <- data.table(Block = c("B1", "B1", "B2"), Effect = c(10, 5, 8),
                         CFreq = c(0.5, 1/3, 0.5), TargetFreq = c(0.6, 0.4, 1))

# Complementary-carrier fixture (P1 -> B1 target, P2 -> B2 target, F* fillers).
.sel_lg <- data.table(
  Block = c("B1", "B2"),
  P1.1 = c(10, 0), P1.2 = c(10, 0),
  P2.1 = c(0, 10), P2.2 = c(0, 10),
  F1.1 = c(0, 0),  F1.2 = c(0, 0),
  F2.1 = c(0, 0),  F2.2 = c(0, 0),
  F3.1 = c(0, 0),  F3.2 = c(0, 0)
)

# Rare-target fixture: R is the lone carrier of the B3 target.
.rare_lg <- data.table(
  Block = c("B1", "B2", "B3"),
  P1.1  = c(10, 0, 0), P1.2  = c(10, 0, 0),
  P1b.1 = c(10, 0, 0), P1b.2 = c(10, 0, 0),
  P2.1  = c(0, 10, 0), P2.2  = c(0, 10, 0),
  R.1   = c(0, 0, 7),  R.2   = c(0, 0, 0)
)

# Random wide haplotype-GEBV table: `n_blocks` blocks x `n_samples` samples,
# ploidy 2, each haplotype drawn from `n_states` effect states per block.
.random_lg <- function(n_blocks, n_samples, n_states = 4L, seed = 1L) {
  set.seed(seed)
  states <- lapply(seq_len(n_blocks), function(b) round(sort(rnorm(n_states, sd = 3)), 2))
  cols <- list(Block = paste0("B", seq_len(n_blocks)))
  for (s in seq_len(n_samples)) for (h in 1:2)
    cols[[paste0("S", s, ".", h)]] <-
      vapply(seq_len(n_blocks), function(b) sample(states[[b]], 1L), numeric(1L))
  as.data.table(cols)
}

# ── oracle: Kemper's sequential greedy, transcribed from the formulas ─────────

# HapN[c, t], TargetFreq, Effect per target state; SelFreq projected to the
# final set size; ARCSIN asin sqrt(SelFreq) - asin sqrt(TargetFreq); LOG
# TargetFreq * log(SelFreq); absent -> penalty; Value = HapN * dev * (1 + |E| w);
# ascending rank; pick one, recompute, repeat.
.oracle_kemper <- function(localgebv, ideal, n_parents, method, penalty = -100,
                           effect_weight = 0) {
  long <- melt(as.data.table(localgebv), id.vars = "Block",
               variable.name = "SH", value.name = "Effect", variable.factor = FALSE)
  long[, Sample := sub("\\.[0-9]+$", "", SH)]
  ploidy  <- max(as.integer(sub("^.*\\.", "", long$SH)))
  samples <- sort(unique(long$Sample))
  tg      <- as.data.table(ideal)[TargetFreq > 0]
  hapn <- sapply(seq_len(nrow(tg)), function(t)
    vapply(samples, function(s)
      sum(long$Sample == s & long$Block == tg$Block[t] & long$Effect == tg$Effect[t]),
      numeric(1L)))
  hapn <- matrix(hapn, nrow = length(samples))
  w <- 1 + abs(tg$Effect) * effect_weight
  sel <- character(0L)
  for (k in seq_len(n_parents)) {
    selfreq <- if (!length(sel)) rep(0, nrow(tg))
               else colSums(hapn[match(sel, samples), , drop = FALSE]) / (n_parents * ploidy)
    dev <- vapply(seq_len(nrow(tg)), function(t) {
      if (selfreq[t] <= 0) return(penalty)
      if (method == "ARCSIN") asin(sqrt(selfreq[t])) - asin(sqrt(tg$TargetFreq[t]))
      else tg$TargetFreq[t] * log(selfreq[t])
    }, numeric(1L))
    value <- as.vector(hapn %*% (dev * w))
    value[match(sel, samples)] <- Inf
    sel <- c(sel, samples[which.min(value)])
  }
  sel
}

# ── (a) hand-computed scores on the toy ──────────────────────────────────────

test_that("(a) candidate scores match the hand formula for ARCSIN and LOG", {
  mm  <- haploStack:::.gbs_hapn_matrix(.toy_lg, .toy_ideal)
  spec <- haploStack:::.gbs_apply_preset("kemper")
  for (m in c("ARCSIN", "LOG")) {
    ctx <- haploStack:::.gbs_context(mm, 2L, m, spec)
    # step 1: every target absent -> -100 * copies carried
    s1 <- haploStack:::.gbs_need_score(ctx, rep(0, 3))
    expect_equal(unname(s1), unname(-100 * rowSums(mm$HapN)), info = m)
    # step 2 after picking A: SelFreq = HapN[A,] / (2 * 2) = 0.25 each
    p <- rep(0.25, 3); tau <- .toy_ideal$TargetFreq
    d <- if (m == "ARCSIN") asin(sqrt(p)) - asin(sqrt(tau)) else tau * log(p)
    s2 <- haploStack:::.gbs_need_score(ctx, p)
    expect_equal(unname(s2), as.vector(mm$HapN %*% d), info = m)
    expect_setequal(gbs_selection(.toy_lg, .toy_ideal, n_parents = 2L,
                                  gbs_method = m, gbs_preset = "kemper"),
                    c("A", "C"))
  }
})

test_that("(a) set objective F(S) matches the hand formula", {
  mm   <- haploStack:::.gbs_hapn_matrix(.toy_lg, .toy_ideal)
  spec <- haploStack:::.gbs_apply_preset("kemper")
  tau  <- .toy_ideal$TargetFreq
  p    <- c(0.25, 0, 0.5)                        # t2 absent -> penalty term
  ctxA <- haploStack:::.gbs_context(mm, 2L, "ARCSIN", spec)
  ctxL <- haploStack:::.gbs_context(mm, 2L, "LOG", spec)
  F_A  <- -(asin(sqrt(0.25)) - asin(sqrt(tau[1])))^2 + (-100) -
           (asin(sqrt(0.5)) - asin(sqrt(tau[3])))^2
  F_L  <- tau[1] * log(0.25) + (-100) + tau[3] * log(0.5)
  expect_equal(haploStack:::.gbs_set_objective(ctxA, p), F_A)
  expect_equal(haploStack:::.gbs_set_objective(ctxL, p), F_L)
})

# ── (b) the Kemper preset reproduces the transcribed sequential greedy ────────

test_that("(b) gbs_preset = 'kemper' reproduces Kemper's formula picks (both modes)", {
  lg <- .random_lg(30L, 25L, seed = 11L)
  ig <- construct_ideal_genotype(lg, method = "custom")
  for (m in c("ARCSIN", "LOG")) {
    got <- gbs_selection(lg, ig, n_parents = 6L, gbs_method = m, gbs_preset = "kemper")
    exp <- .oracle_kemper(lg, ig, 6L, m, penalty = -100)
    expect_identical(got, exp, info = m)
  }
  # effect weight enters exactly as in the origin transcription
  got <- gbs_selection(lg, ig, n_parents = 6L, gbs_method = "ARCSIN",
                       gbs_preset = "kemper", gbs_effect_weight = 0.5)
  expect_identical(got, .oracle_kemper(lg, ig, 6L, "ARCSIN", -100, 0.5))
})

test_that("LOG ranks ascending: a target-free candidate is never picked before a carrier", {
  sel <- gbs_selection(.sel_lg, n_parents = 2L, gbs_method = "LOG", gbs_preset = "kemper")
  expect_setequal(sel, c("P1", "P2"))
})

# ── (c) fixed penalty: base x scale ──────────────────────────────────────────

test_that("(c) fixed penalty is gbs_penalty * penalty_scale for every state", {
  pf <- c(0.5, 0.01, 0.001)
  expect_equal(haploStack:::.resolve_gbs_penalty(-10, "fixed", 1,  pf, rep(.25, 3), 16L, 2L, 500L),
               rep(-10, 3))
  expect_equal(haploStack:::.resolve_gbs_penalty(-10, "fixed", 10, pf, rep(.25, 3), 16L, 2L, 500L),
               rep(-100, 3))
  expect_equal(haploStack:::.gbs_apply_preset("kemper")$penalty_mode, "fixed")
  expect_equal(haploStack:::.gbs_apply_preset("kemper")$penalty_scale, 10)
  expect_equal(haploStack:::.gbs_apply_preset("kemper")$gbs_penalty, -10)
})

# ── (d) dynamic scale: the risk x hope curve ─────────────────────────────────

test_that("(d) dynamic scale follows the plan §4.3 table and its properties", {
  sc <- function(pf, n, N = 500L, tau = 0.25)
    haploStack:::.gbs_penalty_scale("dynamic", 10, pf, rep(tau, length(pf)), n, 2L, N)
  pf <- c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001)
  s16 <- sc(pf, 16L)
  expect_equal(s16[1], 1, tolerance = 0.01)            # common state: origin -10
  expect_gt(s16[4], 7); expect_lt(s16[4], 8)           # 10-copy state near Kemper
  expect_lt(s16[6], s16[4])                            # 1-copy (hope 1/8) below 10-copy
  expect_true(all(s16 >= 1 & s16 <= 10))               # bounded [1, cap]
  # pool vs parents: the same pop_f is safer with more parents
  s4 <- sc(pf, 4L); s64 <- sc(pf, 64L)
  expect_true(all(s4 >= s16 - 1e-12))
  expect_true(all(s16 >= s64 - 1e-12))
  # hand value at pop_f = 0.01, n = 16: 1 + 9 * 0.99^32 * min(1, 10 / 8)
  expect_equal(s16[4], 1 + 9 * 0.99^32)
})

test_that("(d) an absent target state is recruited first under the Kemper preset", {
  for (m in c("ARCSIN", "LOG")) {
    sel <- gbs_selection(.rare_lg, n_parents = 3L, gbs_method = m, gbs_preset = "kemper")
    expect_true("R" %in% sel, info = m)
    expect_setequal(sel, c("P1", "P2", "R"))
  }
})

test_that("penalty arguments are validated", {
  expect_error(gbs_selection(.sel_lg, n_parents = 2L, gbs_penalty = 5), "non-positive")
  expect_error(gbs_selection(.sel_lg, n_parents = 2L, penalty_scale = 0.5), ">= 1")
  expect_error(gbs_selection(.sel_lg, n_parents = 2L, gbs_engine = "GA"), "set")
})

test_that("an explicit argument overrides its preset value", {
  spec <- haploStack:::.gbs_apply_preset("kemper", list(penalty_mode = "dynamic",
                                                        gbs_effect_weight = 0.3))
  expect_equal(spec$penalty_mode, "dynamic")
  expect_equal(spec$gbs_effect_weight, 0.3)
  expect_equal(spec$penalty_scale, 10)                 # untouched preset entry
})

# ── (e) set-objective greedy = argmax of F(S u c) at every step ──────────────

test_that("(e) the set-objective greedy picks the brute-force best marginal at every step", {
  lg <- .random_lg(12L, 10L, seed = 3L)
  ig <- construct_ideal_genotype(lg, method = "custom")
  mm <- haploStack:::.gbs_hapn_matrix(lg, ig)
  spec <- haploStack:::.gbs_apply_preset("haploStack")
  ctx  <- haploStack:::.gbs_context(mm, 4L, "ARCSIN", spec)
  got  <- gbs_selection(lg, ig, n_parents = 4L, gbs_objective = "set")
  S <- integer(0L)
  for (k in 1:4) {
    fs <- vapply(seq_along(mm$samples), function(c) {
      if (c %in% S) return(-Inf)
      haploStack:::.gbs_set_objective(
        ctx, colSums(mm$HapN[c(S, c), , drop = FALSE]) / ctx$denom)
    }, numeric(1L))
    best <- which.max(fs)
    expect_identical(got[k], mm$samples[best], info = paste("step", k))
    S <- c(S, best)
  }
  # and the chosen set beats fillers on D2Ideal
  d_sel <- haploStack:::.distance_to_ideal(lg, ig, samples = got)
  d_all <- haploStack:::.distance_to_ideal(lg, ig)
  expect_true(is.finite(d_sel) && d_sel >= 0)
})

# ── (f) GA reaches at least the greedy set objective ─────────────────────────

test_that("(f) gbs_engine = 'GA' reaches >= the greedy F(S) (seeded)", {
  lg <- .random_lg(40L, 30L, seed = 7L)
  ig <- construct_ideal_genotype(lg, method = "custom")
  mm <- haploStack:::.gbs_hapn_matrix(lg, ig)
  spec <- haploStack:::.gbs_apply_preset("haploStack")
  ctx  <- haploStack:::.gbs_context(mm, 6L, "ARCSIN", spec)
  FS <- function(sel) haploStack:::.gbs_set_objective(
    ctx, colSums(mm$HapN[match(sel, mm$samples), , drop = FALSE]) / ctx$denom)
  g  <- gbs_selection(lg, ig, n_parents = 6L, gbs_objective = "set")
  ga <- gbs_selection(lg, ig, n_parents = 6L, gbs_objective = "set",
                      gbs_engine = "GA", ga_control = list(n_generations = 30L), seed = 9L)
  expect_length(ga, 6L)
  expect_gte(FS(ga), FS(g) - 1e-9)
})

test_that("the GA engine leaves the global RNG stream where it found it", {
  lg <- .random_lg(10L, 12L, seed = 2L)
  set.seed(123); before <- runif(1); set.seed(123)
  invisible(gbs_selection(lg, n_parents = 3L, gbs_objective = "set",
                          gbs_engine = "GA", ga_control = list(n_generations = 5L)))
  expect_equal(runif(1), before)
})

# ── (g) smoke over the whole vocabulary ──────────────────────────────────────

test_that("(g) every method x objective x penalty mode x preset runs", {
  lg <- .random_lg(15L, 12L, seed = 5L)
  ig <- construct_ideal_genotype(lg)
  for (m in c("ARCSIN", "LOG")) for (o in c("candidate", "set"))
    for (pm in c("dynamic", "fixed")) for (pr in c("haploStack", "kemper")) {
      sel <- gbs_selection(lg, ig, n_parents = 4L, gbs_method = m, gbs_objective = o,
                           penalty_mode = pm, gbs_preset = pr)
      expect_equal(length(sel), 4L, info = paste(m, o, pm, pr))
      expect_equal(length(unique(sel)), 4L, info = paste(m, o, pm, pr))
    }
})

# ── (h) D2Ideal under LOG is a divergence ────────────────────────────────────

test_that("(h) D2Ideal(LOG) is zero at target, >= 0, and worse with a lost state", {
  # A set whose pooled frequencies equal the target exactly: ideal from that set.
  lg <- .random_lg(10L, 8L, seed = 4L)
  # pooled frequency of the whole pool (8 samples x 2 haplotypes) as the target
  # -> distance 0 for both metrics
  full <- haploStack:::.localgebv_long(lg)[, .(CFreq = .N / 16), by = .(Block, Effect)]
  ideal_pool <- full[, .(Block, Effect, CFreq, TargetFreq = CFreq)]
  expect_equal(haploStack:::.distance_to_ideal(lg, ideal_pool, gbs_method = "LOG"), 0,
               tolerance = 1e-12)
  expect_equal(haploStack:::.distance_to_ideal(lg, ideal_pool, gbs_method = "ARCSIN"), 0,
               tolerance = 1e-12)
  # random subsets: non-negative at gbs_effect_weight = 0
  set.seed(8)
  for (i in 1:5) {
    s <- sample(paste0("S", 1:8), 4L)
    expect_gte(haploStack:::.distance_to_ideal(lg, ideal_pool, samples = s, gbs_method = "LOG"), -1e-12)
  }
  # a set missing one target state is worse than one carrying all of them
  ig2 <- construct_ideal_genotype(.sel_lg, method = "argmax")
  d_ok  <- haploStack:::.distance_to_ideal(.sel_lg, ig2, c("P1", "P2"), gbs_method = "LOG")
  d_bad <- haploStack:::.distance_to_ideal(.sel_lg, ig2, c("P1", "F1"), gbs_method = "LOG")
  expect_lt(d_ok, d_bad)
})

# ── tune_gbs_penalty grids the scale ─────────────────────────────────────────

test_that("tune_gbs_penalty grids penalty_scale under either mode", {
  for (pm in c("dynamic", "fixed")) {
    res <- tune_gbs_penalty(.rare_lg, n_parents = 3L, penalty_mode = pm,
                            grid = c(1, 10), n_rep = 2L, subsample_frac = 1)
    expect_true(res$best_penalty_scale %in% c(1, 10), info = pm)
    expect_equal(res$best_penalty, -10 * res$best_penalty_scale, info = pm)
  }
  expect_error(tune_gbs_penalty(.rare_lg, n_parents = 3L, grid = c(0.5, 1)), ">= 1")
})

# ── driver: ideal_scope and the Kemper preset ────────────────────────────────

.modes_fixture <- function(seed = 5L) {
  m <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = seed)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))
  list(meta = m, eff = e)
}

test_that("ideal_scope = 'founder' keeps one ideal for the run; 'cycle' rebuilds it", {
  fx <- .modes_fixture()
  run <- function(scope) run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 3L,
    top_founders = 4L, num_offspring = 3L, scheme = "GBS", metrics = "d2ideal",
    ideal_scope = scope, select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)
  a <- run("founder")
  # the founder ideal_phenotype is repeated in every cycle record
  expect_equal(length(unique(round(a$summary$mean_ideal_phenotype, 10))), 1L)
  expect_equal(a$summary$mean_ideal_phenotype[1],
               as.numeric(attr(a$ideal_genotypes, "ideal_phenotype")))
  b <- run("cycle")
  expect_named(b$summary, names(a$summary))
  expect_gt(nrow(b$plans), 0L)
})

test_that("the driver runs the Kemper preset and the set/GA options", {
  fx <- .modes_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
    top_founders = 4L, num_offspring = 2L, scheme = "GBS",
    gbs_preset = "kemper", gbs_method = "LOG",
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)
  expect_gt(nrow(res$plans), 0L)
  res2 <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L,
    top_founders = 4L, num_offspring = 2L, scheme = "GBS",
    gbs_objective = "set", gbs_engine = "GA",
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)
  expect_gt(nrow(res2$plans), 0L)
  expect_error(run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 1L,
    top_founders = 4L, num_offspring = 2L, scheme = "GBS", gbs_engine = "GA",
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE), "set")
})
