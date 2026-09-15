# =============================================================================
# test-stacking-pairs-schemes.R
#
# Direct pair selection — the scheme scorers (R/stacking_pairs.R, thread #031,
# dev/UPGRADE_PLAN_pair_selection.md Phases A-B): TS mid-parent / usefulness,
# OHS on M[i, j], GBS on the family's expected offspring — the reduction
# identity to the parents-first frequency, the separability of the candidate
# objective, the loss-weighted penalty, the two levels — and the driver under
# selection_unit = "pairs".
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

.pairs_fixture <- function(seed = 5L) {
  m <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L,
                     block_size = 3L, seed = seed)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("1", "4"))
  list(meta = m, eff = e,
       lg = compute_profile_localgebv(m, e, level = "haplotype"))
}

# 3 samples x 2 blocks, hand-computable.
.toy_lg <- function() {
  data.table(Block = c("b1", "b2"),
             S1.1 = c(5, 0), S1.2 = c(1, 0),
             S2.1 = c(4, 2), S2.2 = c(3, 2),
             S3.1 = c(0, 6), S3.2 = c(0, 1))
}

# ── OHS pairs ─────────────────────────────────────────────────────────────────

test_that("OHS pairs rank families by the hand-computed M[i, j]", {
  lg <- .toy_lg()
  # M[1,2] = top2{5,1,4,3} + top2{0,0,2,2} = 9 + 4 = 13
  # M[1,3] = top2{5,1,0,0} + top2{0,0,6,1} = 6 + 7 = 13
  # M[2,3] = top2{4,3,0,0} + top2{2,2,6,1} = 7 + 8 = 15
  M <- haploStack:::.pairwise_block_value_matrix(lg)
  expect_equal(M["S1", "S2"], 13); expect_equal(M["S1", "S3"], 13)
  expect_equal(M["S2", "S3"], 15)
  out <- ohs_pair_selection(lg, "families", top_founders = 2)
  expect_s3_class(out, "pair_selection")
  expect_equal(nrow(out), 1L)
  expect_setequal(c(out$Parent1, out$Parent2), c("S2", "S3"))
  expect_equal(out$score, -15)
  expect_error(ohs_pair_selection(lg, "families", 2, ohs_method = "GA_coverage"),
               "GA_coverage")
})

test_that("OHS pairs under 'parents' is the best set's traversal, ranked by M", {
  fx  <- .pairs_fixture()
  out <- ohs_pair_selection(fx$lg, "parents", top_founders = 4, ohs_method = "Greedy")
  expect_equal(attr(out, "n_parents_used"), 4L)
  expect_equal(nrow(out), choose(4L, 2L))
  S <- ohs_selection(fx$lg, n_parents = 4L, method = "Greedy")
  expect_setequal(unique(c(out$Parent1, out$Parent2)), S)
  expect_true(!is.unsorted(out$score))                 # -M ascending = M descending
})

test_that("a self-cross is offered under OHS only when allowed, with its doubled ceiling", {
  lg <- .toy_lg()
  none <- ohs_pair_selection(lg, "families", 3, allow_selfing = FALSE)
  expect_false(any(none$Parent1 == none$Parent2))
  # self ceiling of S1: b1 top2{5,5,1,1} = 10, b2 0 -> 10; S3: b1 0, b2 top2{6,6,1,1} = 12
  v <- haploStack:::.pair_self_ceiling(lg, c("S1", "S2", "S3"), 2L)
  expect_equal(v, c(10, 12, 12))
  some <- ohs_pair_selection(lg, "families", 3, allow_selfing = TRUE)
  expect_true(any(some$Parent1 == some$Parent2))
})

# ── TS pairs ──────────────────────────────────────────────────────────────────

test_that("TS pairs: the first family is (best, second-best) by total GEBV", {
  fx  <- .pairs_fixture()
  g   <- sort(haploStack:::.stk_total_gebv(fx$lg), decreasing = TRUE)
  sc  <- haploStack:::.pair_scorer_ts(fx$lg)
  out <- haploStack:::.pairs_greedy(sc, pair_selection_budget(4, "families"))
  expect_setequal(c(out$Parent1[1L], out$Parent2[1L]), names(g)[1:2])
  expect_equal(out$score[1L], -(g[1L] + g[2L]) / 2, ignore_attr = TRUE)
})

test_that("the F1 standard deviation of a family is exact under the block model", {
  lg <- .toy_lg()
  # per-sample sum over blocks of the population variance of its two haplotypes:
  # S1: (5-1)^2/4 + 0 = 4; S2: (4-3)^2/4 + 0 = 0.25; S3: 0 + (6-1)^2/4 = 6.25
  sdm <- haploStack:::.pair_f1_sd(lg, c("S1", "S2", "S3"))
  expect_equal(sdm["S1", "S2"], sqrt(4 + 0.25))
  expect_equal(sdm["S1", "S3"], sqrt(4 + 6.25))
  expect_equal(sdm["S1", "S1"], sqrt(8))              # a self: the parent twice
})

test_that("usefulness ranks a segregating family above a fixed one of equal mid-parent", {
  lg <- data.table(Block = "b1",
                   S1.1 = 2, S1.2 = 2, S2.1 = 2, S2.2 = 2,   # fixed, mid-parent 4
                   S3.1 = 4, S3.2 = 0, S4.1 = 4, S4.2 = 0)   # segregating, mid-parent 4
  b <- pair_selection_budget(2, "families")
  mp <- haploStack:::.pairs_greedy(haploStack:::.pair_scorer_ts(lg), b)
  uc <- haploStack:::.pairs_greedy(
    haploStack:::.pair_scorer_ts(lg, ts_pair_score = "usefulness"), b)
  expect_equal(mp$score, -4)                          # every family ties at -4; first by order
  expect_setequal(c(uc$Parent1, uc$Parent2), c("S3", "S4"))
  expect_equal(uc$score, -4 - sqrt(2 * (4 - 0)^2 / 4))
})

test_that("the mid-parent is the expected F1 mean under the block backend", {
  fx  <- .pairs_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 1L, top_founders = 2L,
    num_offspring = 200L, scheme = "TS", selection_unit = "pairs",
    pair_budget = "families", select_top_blocks = 50L, output_dir = NULL,
    verbose = FALSE)
  fam <- res$plans[cycle == 1L]
  expect_equal(nrow(fam), 1L)
  g   <- haploStack:::.stk_total_gebv(fx$lg)
  mp  <- (g[fam$Parent1] + g[fam$Parent2]) / 2
  sdm <- haploStack:::.pair_f1_sd(fx$lg, names(g))
  se  <- sdm[fam$Parent1, fam$Parent2] / sqrt(200)
  expect_lt(abs(fam$F1mean - mp), 4 * se + 1e-8)
  expect_equal(fam$pair_score, -mp, ignore_attr = TRUE)
})

# ── GBS pairs ─────────────────────────────────────────────────────────────────

.gbs_parts <- function(fx, top_founders = 4L, budget = "families") {
  ideal <- construct_ideal_genotype(fx$lg)
  mm    <- haploStack:::.gbs_hapn_matrix(fx$lg, ideal)
  spec  <- haploStack:::.gbs_apply_preset("haploStack", list())
  b     <- pair_selection_budget(top_founders, budget)
  ctx   <- haploStack:::.gbs_context(mm, b$n_parents_ref, "ARCSIN", spec)
  list(ideal = ideal, mm = mm, ctx = ctx, b = b, spec = spec)
}

test_that("reduction identity: the traversal of S has the parents-first frequency p_t(S)", {
  fx <- .pairs_fixture()
  pp <- .gbs_parts(fx, 4L)
  sc <- haploStack:::.pair_scorer_gbs(pp$mm, pp$ctx, pp$b, n_offspring = 4L,
                                      pair_pool = "all", pair_freq = "projected",
                                      pair_penalty = "hard")
  S   <- c("S1", "S4", "S7", "S10")
  idx <- match(S, pp$mm$samples)
  st  <- sc$init
  for (k in seq_along(idx)) for (l in seq_along(idx)) if (k < l)
    st <- sc$update(st, idx[k], idx[l])
  expect_equal(st$Fc, choose(4L, 2L))
  p_pairs   <- st$sumq / pp$b$n_families                       # projected to F = 6
  p_parents <- colSums(pp$ctx$HapN[idx, , drop = FALSE]) / (4L * pp$mm$ploidy)
  expect_equal(p_pairs, p_parents, tolerance = 1e-12)
})

test_that("the candidate objective over a family is separable: score(i, j) = s_i + s_j", {
  fx <- .pairs_fixture()
  pp <- .gbs_parts(fx, 4L)
  # parents-first need scores at the empty set (every state absent -> beta)
  s_par <- haploStack:::.gbs_need_score(pp$ctx, rep(0, ncol(pp$ctx$HapN)))
  names(s_par) <- pp$mm$samples
  out <- gbs_pair_selection(fx$lg, pp$ideal, "families", top_founders = 4,
                            gbs_objective = "candidate", pair_penalty = "hard",
                            pair_pool = "all", pair_freq = "projected")
  cand <- haploStack:::.pair_candidates(length(s_par))
  sums <- s_par[cand$i] + s_par[cand$j]
  expect_equal(out$score[1L], min(sums), ignore_attr = TRUE)
  expect_equal(out$score[1L], s_par[out$Parent1[1L]] + s_par[out$Parent2[1L]],
               ignore_attr = TRUE)
})

test_that("the first family of the set objective maximises F over every single family", {
  fx <- .pairs_fixture()
  pp <- .gbs_parts(fx, 4L)
  m  <- pp$mm$ploidy
  cand <- haploStack:::.pair_candidates(length(pp$mm$samples))
  Fv <- vapply(seq_len(nrow(cand)), function(r) {
    q <- (pp$ctx$HapN[cand$i[r], ] + pp$ctx$HapN[cand$j[r], ]) / (2 * m)
    haploStack:::.gbs_set_objective(pp$ctx, q / pp$b$n_families)   # projected, F = 6
  }, numeric(1L))
  best <- which.max(Fv)
  out <- gbs_pair_selection(fx$lg, pp$ideal, "families", top_founders = 4,
                            pair_penalty = "hard", pair_pool = "all",
                            pair_freq = "projected")
  expect_setequal(c(out$Parent1[1L], out$Parent2[1L]),
                  pp$mm$samples[c(cand$i[best], cand$j[best])])
  # F(empty) = sum(w * beta); the score is the negated marginal gain
  F0 <- sum(pp$ctx$w * pp$ctx$beta)
  expect_equal(out$score[1L], -(Fv[best] - F0), ignore_attr = TRUE)
})

test_that("the loss-weighted penalty tends to the hard penalty as M grows", {
  fx <- .pairs_fixture()
  pp <- .gbs_parts(fx, 4L)
  hard <- gbs_pair_selection(fx$lg, pp$ideal, "families", top_founders = 4,
                             pair_penalty = "hard", pair_pool = "all")
  loss <- gbs_pair_selection(fx$lg, pp$ideal, "families", top_founders = 4,
                             num_offspring = 1000L, pair_penalty = "loss",
                             pair_pool = "all")
  expect_equal(loss$Parent1[1L], hard$Parent1[1L])
  expect_equal(loss$Parent2[1L], hard$Parent2[1L])
  expect_equal(loss$score[1L], hard$score[1L], tolerance = 1e-6)
  # and at M = 1 a state carried by one heterozygous family is still at risk:
  # the first pick may differ, but the run is valid and complete
  one <- gbs_pair_selection(fx$lg, pp$ideal, "families", top_founders = 4,
                            num_offspring = 1L, pair_penalty = "loss",
                            pair_pool = "all")
  expect_equal(nrow(one), 6L)
})

test_that("with one family per cycle the family and population levels agree", {
  fx <- .pairs_fixture()
  pp <- .gbs_parts(fx, 2L)
  pop <- gbs_pair_selection(fx$lg, pp$ideal, "families", top_founders = 2,
                            gbs_pair_level = "population", pair_pool = "all")
  fam <- gbs_pair_selection(fx$lg, pp$ideal, "families", top_founders = 2,
                            gbs_pair_level = "family", pair_pool = "all")
  expect_equal(nrow(pop), 1L); expect_equal(nrow(fam), 1L)
  expect_setequal(c(pop$Parent1, pop$Parent2), c(fam$Parent1, fam$Parent2))
})

test_that("GBS pairs: the shortlist agrees with the full triangle on the fixture", {
  fx <- .pairs_fixture()
  pp <- .gbs_parts(fx, 4L)
  all_ <- gbs_pair_selection(fx$lg, pp$ideal, "families", top_founders = 4,
                             pair_pool = "all")
  short <- gbs_pair_selection(fx$lg, pp$ideal, "families", top_founders = 4,
                              pair_pool = "shortlist", shortlist_k = 6L)
  expect_equal(nrow(short), 6L)
  key <- function(x) paste(pmin(x$Parent1, x$Parent2), pmax(x$Parent1, x$Parent2))
  expect_gte(length(intersect(key(all_), key(short))), 4L)
})

test_that("GBS pairs run under both projections, both links, both presets, selfing", {
  fx <- .pairs_fixture()
  pp <- .gbs_parts(fx, 4L)
  for (fr in c("projected", "running")) for (lk in c("ARCSIN", "LOG"))
    for (pr in c("haploStack", "kemper")) {
      out <- gbs_pair_selection(fx$lg, pp$ideal, "parents", top_founders = 4,
                                gbs_method = lk, gbs_preset = pr, pair_freq = fr)
      expect_s3_class(out, "pair_selection")
      expect_equal(attr(out, "n_parents_used"), 4L)
      expect_true(all(is.finite(out$score)))
    }
  s <- gbs_pair_selection(fx$lg, pp$ideal, "families", top_founders = 3,
                          allow_selfing = TRUE, pair_pool = "all")
  expect_equal(nrow(s), 6L)                            # choose(3, 2) + 3 selfs
  expect_error(gbs_pair_selection(fx$lg, pp$ideal, "families", 4,
                                  gbs_objective = "candidate",
                                  gbs_pair_level = "family"), "population-level")
})

test_that("the pairs penalty scale is the one of the matched parents-first run", {
  fx <- .pairs_fixture()
  pp <- .gbs_parts(fx, 4L)
  ctx_parents <- haploStack:::.gbs_context(pp$mm, 4L, "ARCSIN", pp$spec)
  expect_equal(pp$ctx$beta, ctx_parents$beta)
})

# ── driver ────────────────────────────────────────────────────────────────────

test_that("the driver runs TS / OHS / GBS under 'pairs' with the block backend", {
  fx <- .pairs_fixture()
  for (bud in c("families", "parents")) {
    res <- run_stacking_simulation(
      fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
      num_offspring = 2L, scheme = c("TS", "OHS", "GBS"), selection_unit = "pairs",
      pair_budget = bud, ohs_method = "Greedy", metrics = "all",
      select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)
    s <- res$summary[cycle > 0L]
    expect_true(all(!is.na(s$n_families)))
    expect_true(all(!is.na(res$plans$pair_score)))
    if (bud == "families") expect_true(all(s$n_families == 6L))
    else expect_true(all(s$n_parents_used == 4L))
    expect_true(all(grepl("^C[0-9]+F[0-9]+$", res$plans$CrossID)))
    expect_false(any(res$plans$Parent1 == res$plans$Parent2))
  }
})

test_that("avoid_selfing = FALSE offers selfs under 'pairs' and counts them in the budget", {
  fx <- .pairs_fixture()
  expect_message(
    res <- run_stacking_simulation(
      fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 3L,
      num_offspring = 2L, scheme = "GBS", selection_unit = "pairs",
      pair_budget = "families", avoid_selfing = FALSE, pair_pool = "all",
      select_top_blocks = 50L, output_dir = NULL, verbose = TRUE),
    "6 families per cycle.*self-crosses counted")
  expect_true(all(res$summary[cycle > 0L]$n_families == 6L))
})

test_that("the driver runs 'pairs' with the genomicSimulation backend", {
  skip_if_not_installed("genomicSimulation")
  fx <- .pairs_fixture()
  res <- run_stacking_simulation(
    fx$meta, fx$eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", selection_unit = "pairs",
    pair_budget = "families", cross_backend = "genomicSimulation",
    select_top_blocks = 50L, output_dir = NULL, verbose = FALSE)
  expect_true(all(res$summary[cycle > 0L]$n_families == 6L))
})
