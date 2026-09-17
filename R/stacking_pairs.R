# =============================================================================
# stacking_pairs.R
#
# Direct pair (family) selection — `selection_unit = "pairs"` in
# run_stacking_simulation() (thread #031, dev/UPGRADE_PLAN_pair_selection.md).
#
# A scheme is a SCORER of families; the unit is the ENGINE. Under the
# parents-first protocol every scheme scores a parent set and the crossing
# design is the full traversal. Here every scheme scores a family — the
# expected F1 of one cross (i, j) — and one shared greedy (.pairs_greedy())
# adds families until the budget of pair_selection_budget() stops it:
#   "families"  exactly n_families families (the traversal's count);
#   "parents"   families best first until n_parents_max distinct names are in
#               use; a family that would bring more new names than the slack
#               allows is skipped (overshoot = "skip", the search then continues
#               among the names in use up to the traversal cap) or ends the
#               search (overshoot = "stop").
#
# Family scores (lower = better; the engine minimises):
#   TS   -(g_i + g_j) / 2 — the mid-parent value = the expected F1 mean;
#        optionally the usefulness (g_i + g_j) / 2 + k * sd_F1(i, j), with the
#        F1 standard deviation exact under the block model (.pair_f1_sd()).
#   OHS  -M[i, j] — the family's stacking ceiling (.pairwise_block_value_matrix()).
#        Under the "parents" budget the optimum of an additive family value
#        with a name limit is the best parent SET and every family among it,
#        i.e. sum(M[S, S]) — the parents-first objective — so that budget routes
#        through ohs_selection() and the traversal (elite_ratio / elite_level
#        passed along). With elite_level = "block" and a ratio < 1 the family
#        greedy scores the set's per-block elite instead
#        (.pair_scorer_ohs_block(), non-additive).
#   GBS  the ideal-genotype objective on the family's expected offspring,
#        q_t(i, j) = (HapN[i, t] + HapN[j, t]) / (2m). Level "population": the
#        marginal gain of the set objective F(P) at the offspring-pool frequency
#        p_t(P) = mean_f q_t(f) — which reduces to the parents-first p_t(S)
#        when P is the traversal of S. Level "family": the family's own distance
#        to the ideal. The absent-state penalty is either hard (p_t(P) == 0) or
#        weighted by the exact probability L_t = prod_f a_f(t)^M that none of the
#        M offspring of the chosen families carries state t.
#
# Functions:
#   ohs_pair_selection(), gbs_pair_selection()      exported
#   .stk_pairs_truncation()                          TS pairs (driver-only)
#   .pair_candidates(), .pairs_greedy()              the shared engine
#   .pair_scorer_ts() / _ohs() / _gbs()              scorers (one interface)
#   .pair_f1_sd(), .pair_self_ceiling()              block-model family quantities
#   .as_pair_selection(), print.pair_selection()     the output class
# =============================================================================


# ---------------------------------------------------------------------------
# Internal — candidate families
# ---------------------------------------------------------------------------

# All canonical families over the sample indices `idx` (default: all `n`):
# i < j, plus the selfs (i, i) when allowed. One place decides orientation, so
# .stk_rename_generation()'s `FamID = .GRP by (Parent1, Parent2)` never sees
# (A, B) and (B, A) as two families.
.pair_candidates <- function(n, allow_selfing = FALSE, idx = NULL) {
  idx <- if (is.null(idx)) seq_len(n) else sort(unique(as.integer(idx)))
  cand <- if (length(idx) >= 2L) {
    cb <- utils::combn(idx, 2L)
    data.table(i = cb[1L, ], j = cb[2L, ])
  } else {
    data.table(i = integer(0L), j = integer(0L))
  }
  if (isTRUE(allow_selfing))
    cand <- rbind(cand, data.table(i = idx, j = idx))
  cand
}


# ---------------------------------------------------------------------------
# Internal — the output class
# ---------------------------------------------------------------------------

.as_pair_selection <- function(dt, scheme, budget, n_parents_used) {
  dt <- as.data.table(dt)
  data.table::setattr(dt, "scheme", scheme)
  data.table::setattr(dt, "pair_budget", budget$pair_budget)
  data.table::setattr(dt, "n_families", nrow(dt))
  data.table::setattr(dt, "n_parents_used", as.integer(n_parents_used))
  data.table::setattr(dt, "class", c("pair_selection", class(dt)))
  dt
}

#' @export
print.pair_selection <- function(x, n = 10L, ...) {
  cat(sprintf("Direct pair selection (%s, pair_budget = \"%s\"): %d families, %d distinct parents\n",
              attr(x, "scheme"), attr(x, "pair_budget"),
              attr(x, "n_families"), attr(x, "n_parents_used")))
  y <- data.table::copy(x)
  data.table::setattr(y, "class", setdiff(class(y), "pair_selection"))
  print(utils::head(y, n), ...)
  if (nrow(y) > n) cat(sprintf("... %d more families\n", nrow(y) - n))
  invisible(x)
}


# ---------------------------------------------------------------------------
# Internal — the engine: greedy over families under a budget
# ---------------------------------------------------------------------------

# `scorer` is a list: scheme, samples, additive (logical), init (the state),
# candidates(state) -> data.table(i, j), score(state, cand) -> numeric (lower
# = better), update(state, i, j) -> state, and optionally stop_when_no_gain
# (logical: once the name limit is reached, stop when the best family no
# longer improves the objective — set objectives only).
#
# Additive scorers are scored once over their candidate set and taken in
# score order (exact under "families" with no usage cap; the best-first
# heuristic otherwise). Non-additive scorers are re-scored after every pick.
.pairs_greedy <- function(scorer, budget, allow_selfing = FALSE, max_use = Inf,
                          overshoot = c("skip", "stop")) {
  overshoot <- match.arg(overshoot)
  samples <- scorer$samples
  N       <- length(samples)
  usage   <- integer(N)
  used    <- logical(N)
  taken   <- matrix(FALSE, N, N)
  sel_i   <- integer(0L); sel_j <- integer(0L); sel_s <- numeric(0L)
  by_parents <- identical(budget$pair_budget, "parents")
  n_target   <- if (by_parents) budget$n_families_max else budget$n_families

  fits_usage <- function(i, j) {
    if (i == j) usage[i] + 2L <= max_use
    else usage[i] + 1L <= max_use && usage[j] + 1L <= max_use
  }
  new_names <- function(i, j) as.integer(!used[i]) + as.integer(!used[j] && i != j)
  take <- function(i, j, s) {
    sel_i <<- c(sel_i, i); sel_j <<- c(sel_j, j); sel_s <<- c(sel_s, s)
    taken[i, j] <<- TRUE; taken[j, i] <<- TRUE
    usage[i] <<- usage[i] + 1L
    usage[j] <<- usage[j] + 1L                    # a self counts twice
    used[i]  <<- TRUE; used[j] <<- TRUE
  }

  state <- scorer$init

  if (isTRUE(scorer$additive)) {
    cand <- scorer$candidates(state)
    if (nrow(cand) > 0L) {
      sc  <- scorer$score(state, cand)
      for (k in order(sc)) {
        i <- cand$i[k]; j <- cand$j[k]
        if (taken[i, j] || !fits_usage(i, j)) next
        if (by_parents) {
          slack <- budget$n_parents_max - sum(used)
          if (new_names(i, j) > slack) {
            if (overshoot == "stop") break else next
          }
        }
        take(i, j, sc[k])
        if (length(sel_i) >= n_target) break
      }
    }
  } else {
    repeat {
      cand <- scorer$candidates(state)
      if (nrow(cand) == 0L) break
      ok <- !taken[cbind(cand$i, cand$j)]
      if (is.finite(max_use))
        ok <- ok & vapply(seq_len(nrow(cand)),
                          function(r) fits_usage(cand$i[r], cand$j[r]), logical(1L))
      slack <- if (by_parents) budget$n_parents_max - sum(used) else NA_integer_
      nn <- if (by_parents) vapply(seq_len(nrow(cand)),
                                   function(r) new_names(cand$i[r], cand$j[r]),
                                   integer(1L)) else NULL
      if (by_parents && overshoot == "skip") ok <- ok & (nn <= slack)
      cand <- cand[ok]
      if (by_parents) nn <- nn[ok]
      if (nrow(cand) == 0L) break
      sc <- scorer$score(state, cand)
      b  <- which.min(sc)
      if (by_parents && overshoot == "stop" && nn[b] > slack) break
      if (by_parents && sum(used) >= budget$n_parents_max &&
          isTRUE(scorer$stop_when_no_gain) && sc[b] >= 0) break
      take(cand$i[b], cand$j[b], sc[b])
      state <- scorer$update(state, cand$i[b], cand$j[b])
      if (length(sel_i) >= n_target) break
    }
  }

  dt <- data.table(Parent1 = samples[sel_i], Parent2 = samples[sel_j],
                   score = sel_s, rank = seq_along(sel_i))
  .as_pair_selection(dt, scorer$scheme, budget, sum(used))
}


# ---------------------------------------------------------------------------
# Internal — block-model family quantities
# ---------------------------------------------------------------------------

# Standard deviation of the F1 total local GEBV of every family, exact under
# the block model: a gamete carries one of a parent's haplotypes per block, so
# a parent contributes sum_b Var(haplotype values at b) (population variance
# over its ploidy copies) and the family's variance is the sum of its two
# parents' contributions (blocks independent; a self counts its parent twice).
# Returns an n x n matrix in `samples` order.
.pair_f1_sd <- function(localgebv, samples) {
  long <- .localgebv_long(localgebv)
  per  <- long[, .(v = {
    e <- Effect[!is.na(Effect)]
    if (length(e)) mean((e - mean(e))^2) else 0
  }), by = .(Sample, Block)][, .(D = sum(v)), by = Sample]
  d <- per$D[match(samples, per$Sample)]
  d[is.na(d)] <- 0
  out <- sqrt(outer(d, d, "+"))
  dimnames(out) <- list(samples, samples)
  out
}

# The stacking ceiling of a self-cross in the OHS currency: the union of the
# parent's haplotypes with itself, i.e. every haplotype twice, top-`ploidy`
# per block (a diploid self can fix its best haplotype: 2 x max per block).
.pair_self_ceiling <- function(localgebv, samples, ploidy) {
  long <- .localgebv_long(localgebv)
  per  <- long[, .(v = {
    e <- sort(Effect[!is.na(Effect)], decreasing = TRUE)
    if (length(e)) sum(sort(rep(e, 2L), decreasing = TRUE)[seq_len(min(ploidy, 2L * length(e)))])
    else 0
  }), by = .(Sample, Block)][, .(C = sum(v)), by = Sample]
  v <- per$C[match(samples, per$Sample)]
  v[is.na(v)] <- 0
  v
}


# ---------------------------------------------------------------------------
# Internal — scorers
# ---------------------------------------------------------------------------

# TS: mid-parent value, optionally the usefulness. `select_priority = "low"`
# takes the lower tail (mid-parent - k * sd) and minimises it.
.pair_scorer_ts <- function(localgebv, allow_selfing = FALSE,
                            ts_pair_score = c("midparent", "usefulness"),
                            usefulness_k = 1,
                            select_priority = c("high", "low")) {
  ts_pair_score   <- match.arg(ts_pair_score)
  select_priority <- match.arg(select_priority)
  g       <- .stk_total_gebv(localgebv)
  samples <- names(g)
  sdm     <- if (ts_pair_score == "usefulness") .pair_f1_sd(localgebv, samples) else NULL
  sgn     <- if (select_priority == "high") -1 else 1
  cand0   <- .pair_candidates(length(samples), allow_selfing)
  list(scheme = "TS", samples = samples, additive = TRUE, init = NULL,
       candidates = function(state) cand0,
       score = function(state, cand) {
         s <- sgn * (g[cand$i] + g[cand$j]) / 2
         if (!is.null(sdm)) s <- s - usefulness_k * sdm[cbind(cand$i, cand$j)]
         unname(s)
       },
       update = function(state, i, j) state)
}

# OHS: the family's stacking ceiling M[i, j]. The diagonal of M is 0 by
# construction (no self in sum(M[S, S])); when selfs are allowed the scorer
# overlays the self ceiling without touching M.
.pair_scorer_ohs <- function(localgebv, allow_selfing = FALSE, ploidy = NULL) {
  M       <- .pairwise_block_value_matrix(localgebv, ploidy = ploidy)
  samples <- rownames(M)
  if (isTRUE(allow_selfing)) {
    if (is.null(ploidy)) {
      cols   <- setdiff(names(as.data.table(localgebv)), "Block")
      ploidy <- max(tabulate(match(sub("\\.[0-9]+$", "", cols), samples)))
    }
    diag(M) <- .pair_self_ceiling(localgebv, samples, ploidy)
  }
  cand0 <- .pair_candidates(length(samples), allow_selfing)
  list(scheme = "OHS", samples = samples, additive = TRUE, init = NULL,
       candidates = function(state) cand0,
       score = function(state, cand) -M[cbind(cand$i, cand$j)],
       update = function(state, i, j) state)
}

# OHS, block-level elite (0.7.3): the chosen family SET is scored by
# sum_b of the q largest per-block ceilings M_b[i, j] over its families,
# q = ceiling(elite_ratio * n_target) with n_target the budget's family count
# (its cap under "parents"). Non-additive: the marginal gain of a candidate
# family is sum_b max(0, v_fb - tau_b), tau_b the q-th largest chosen value at
# block b (the family's full ceiling while fewer than q families are chosen —
# so at ratio 1 the ranking is the additive one, M[i, j]). A tie-break of
# 1e-6 * M[i, j] orders families with equal gain. The per-block values are held
# for the `shortlist_families` best families by M only (default the larger of
# 4 x n_target and 100; NULL and a small pool = every family).
.pair_scorer_ohs_block <- function(localgebv, budget, elite_ratio = 1,
                                   allow_selfing = FALSE, ploidy = NULL,
                                   shortlist_families = NULL) {
  pm      <- .pairwise_block_value_matrix(localgebv, ploidy = ploidy,
                                          return_array = TRUE)
  M       <- pm$M
  A       <- pm$A
  samples <- rownames(M)
  if (isTRUE(allow_selfing))
    diag(M) <- .pair_self_ceiling(localgebv, samples, pm$ploidy)
  n_target <- if (identical(budget$pair_budget, "parents")) budget$n_families_max
              else budget$n_families
  q <- .ohs_elite_q(n_target, elite_ratio)

  cand0 <- .pair_candidates(length(samples), allow_selfing)
  if (is.null(shortlist_families)) shortlist_families <- max(4L * n_target, 100L)
  if (nrow(cand0) > shortlist_families) {
    keep  <- order(-M[cbind(cand0$i, cand0$j)])[seq_len(shortlist_families)]
    cand0 <- cand0[sort(keep)]
  }
  V   <- .ohs_pair_block_values(A, cand0$i, cand0$j)      # families x blocks
  nb  <- ncol(V)
  ckey <- paste(cand0$i, cand0$j)
  key  <- function(i, j) match(paste(i, j), ckey)
  eps  <- 1e-6

  list(scheme = "OHS", samples = samples, additive = FALSE,
       init = list(chosen = matrix(0, 0L, nb)),
       stop_when_no_gain = FALSE,
       q = q, n_candidates = nrow(cand0),
       candidates = function(state) cand0,
       score = function(state, cand) {
         Vc <- V[key(cand$i, cand$j), , drop = FALSE]
         ch <- state$chosen
         gain <- if (nrow(ch) < q) {
           rowSums(Vc)
         } else {
           tau <- apply(ch, 2L, function(x) sort(x, decreasing = TRUE)[q])
           rowSums(pmax(sweep(Vc, 2L, tau, "-"), 0))
         }
         -(gain + eps * M[cbind(cand$i, cand$j)])
       },
       update = function(state, i, j) {
         state$chosen <- rbind(state$chosen, V[key(i, j), ])
         state
       })
}

# GBS: the ideal-genotype objective on families. `mm` from .gbs_hapn_matrix(),
# `ctx` from .gbs_context() (its `beta` sized with the reference protocol's
# parent count, plan Q24 f). State: sumq = sum_f q_t(f) over the chosen
# families, Fc = their count, logL = sum_f M * log a_f(t) (log of the loss
# probability of each target state in the next pool).
.pair_scorer_gbs <- function(mm, ctx, budget, n_offspring,
                             gbs_pair_level = c("population", "family"),
                             gbs_objective  = c("set", "candidate"),
                             pair_freq      = NULL,
                             pair_penalty   = c("loss", "hard"),
                             pair_pool      = c("shortlist", "all"),
                             shortlist_k    = NULL,
                             select_priority = c("high", "low"),
                             allow_selfing  = FALSE) {
  level        <- match.arg(gbs_pair_level)
  objective    <- match.arg(gbs_objective)
  pair_penalty <- match.arg(pair_penalty)
  pair_pool    <- match.arg(pair_pool)
  select_priority <- match.arg(select_priority)
  if (is.null(pair_freq))
    pair_freq <- if (identical(budget$pair_budget, "families")) "projected" else "running"
  pair_freq <- match.arg(pair_freq, c("projected", "running"))
  if (objective == "candidate" && level == "family")
    stop("gbs_pair_selection: gbs_objective = \"candidate\" is a population-level ",
         "need score; use gbs_pair_level = \"population\" with it.")

  HapN    <- ctx$HapN
  m       <- mm$ploidy
  N       <- nrow(HapN)
  Tn      <- ncol(HapN)
  samples <- ctx$samples
  w  <- ctx$w; tf <- ctx$tf; beta <- ctx$beta; method <- ctx$method
  sgn     <- if (select_priority == "high") 1 else -1
  F_final <- if (identical(budget$pair_budget, "families")) budget$n_families
             else budget$n_families_max
  Qp    <- HapN / (2 * m)                              # one parent's share of q_t
  logvp <- n_offspring * log(pmax(1 - HapN / m, 0))    # M * log P(a gamete lacks t)
  floor_q <- if (method == "LOG") 1 / (2 * 2 * m) else 0   # half a copy in 2m haplotypes
  K <- if (is.null(shortlist_k)) 4L * budget$n_parents_ref else as.integer(shortlist_k)
  K <- max(2L, min(K, N))
  use_all <- pair_pool == "all" || K >= N

  p_cur <- function(state) {
    if (state$Fc == 0L) rep(0, Tn)
    else state$sumq / (if (pair_freq == "projected") F_final else state$Fc)
  }
  # Per-state term for frequency rows P (n x T), loss rows L, and the presence
  # of each state in the pool once the candidate is added. Population level:
  # f_t is evaluated at the pool frequency and zeroed where the state is absent
  # (the penalty takes over there; the LOG link must never see log 0).
  # Family level: f_t is the family's own distance, a state the family lacks
  # floored at half a copy under LOG so the distance stays finite.
  term_of <- function(P, L, present, family = FALSE) {
    TF <- matrix(tf,   nrow(P), Tn, byrow = TRUE)
    B  <- matrix(beta, nrow(P), Tn, byrow = TRUE)
    Pf <- P
    if (family) {
      if (method == "LOG") Pf <- pmax(Pf, floor_q)
    } else {
      Pf[!present] <- if (method == "LOG") 1 else 0
    }
    f <- .gbs_set_term(Pf, TF, method)
    if (!family) f[!present] <- 0
    if (pair_penalty == "loss") (1 - L) * f + L * B else ifelse(present, f, B)
  }
  # Objective value F(P) of the current state (population level).
  F_of_state <- function(state) {
    p  <- p_cur(state)
    L  <- exp(state$logL)
    pr <- p > 0
    pf <- p; pf[!pr] <- if (method == "LOG") 1 else 0
    f  <- .gbs_set_term(pf, tf, method); f[!pr] <- 0
    term <- if (pair_penalty == "loss") (1 - L) * f + L * beta else ifelse(pr, f, beta)
    sum(w * term)
  }
  # Score rows given each candidate's q contribution Q (n x T) and log-loss
  # increment LA (n x T); `own` = the candidate's own frequency for the family
  # level (Q for a family; the genotype frequency for a single parent).
  score_rows <- function(state, Q, LA, own) {
    n  <- nrow(Q)
    Lnew <- exp(matrix(state$logL, n, Tn, byrow = TRUE) + LA)
    Pnew <- if (pair_freq == "projected")
      matrix(p_cur(state), n, Tn, byrow = TRUE) + Q / F_final
    else
      (matrix(state$sumq, n, Tn, byrow = TRUE) + Q) / (state$Fc + 1)
    present <- Pnew > 0
    if (level == "population") {
      if (objective == "candidate") {
        p  <- p_cur(state); L <- exp(state$logL); pr <- p > 0
        pf <- p; pf[!pr] <- if (method == "LOG") 1 else 0
        d  <- .gbs_deviation(pf, tf, method); d[!pr] <- 0
        D  <- if (pair_penalty == "loss") (1 - L) * d + L * beta else ifelse(pr, d, beta)
        return(sgn * as.vector((2 * m * Q) %*% (D * w)))
      }
      F_new <- as.vector(term_of(Pnew, Lnew, present) %*% w)
      return(sgn * -(F_new - F_of_state(state)))
    }
    sgn * -as.vector(term_of(own, Lnew, present, family = TRUE) %*% w)
  }
  fam_rows <- function(cand) {
    list(Q  = Qp[cand$i, , drop = FALSE] + Qp[cand$j, , drop = FALSE],
         LA = logvp[cand$i, , drop = FALSE] + logvp[cand$j, , drop = FALSE])
  }
  cand_all <- .pair_candidates(N, allow_selfing)

  list(scheme = "GBS", samples = samples, additive = FALSE,
       stop_when_no_gain = (level == "population" && objective == "set"),
       init = list(sumq = numeric(Tn), Fc = 0L, logL = numeric(Tn)),
       candidates = function(state) {
         if (use_all) return(cand_all)
         # Shortlist: the K parents whose "half-family" (one gamete each) scores
         # best at the current state, then every family among them.
         s   <- score_rows(state, Qp, logvp, HapN / m)
         idx <- order(s)[seq_len(K)]
         .pair_candidates(N, allow_selfing, idx = idx)
       },
       score = function(state, cand) {
         r <- fam_rows(cand)
         score_rows(state, r$Q, r$LA, r$Q)
       },
       update = function(state, i, j) {
         state$sumq <- state$sumq + Qp[i, ] + Qp[j, ]
         state$Fc   <- state$Fc + 1L
         state$logL <- state$logL + logvp[i, ] + logvp[j, ]
         state
       })
}


# ---------------------------------------------------------------------------
# TS pairs (driver-only, mirrors .stk_select_truncation())
# ---------------------------------------------------------------------------

.stk_pairs_truncation <- function(localgebv, budget, opts) {
  scorer <- .pair_scorer_ts(localgebv, allow_selfing = budget$allow_selfing,
                            ts_pair_score = opts$ts_pair_score %||% "midparent",
                            usefulness_k  = opts$usefulness_k %||% 1,
                            select_priority = opts$select_priority %||% "high")
  .pairs_greedy(scorer, budget, allow_selfing = budget$allow_selfing,
                max_use = opts$max_use %||% Inf,
                overshoot = opts$overshoot %||% "skip")
}


# ---------------------------------------------------------------------------
# OHS pairs (exported)
# ---------------------------------------------------------------------------

#' Optimal Haplotype Selection of families (direct pair selection)
#'
#' Chooses the crosses themselves rather than a parent subgroup: every family
#' \code{(i, j)} is valued by its stacking ceiling \code{M[i, j]} — the sum
#' over blocks of the best \code{ploidy} haplotype values among the two parents'
#' haplotypes (\code{\link{.pairwise_block_value_matrix}}, the quantity
#' \code{\link{ohs_selection}} sums over a parent set) — and families are taken
#' under one of the two budgets of \code{\link{pair_selection_budget}}.
#'
#' @section Engines by budget:
#' Under \code{pair_budget = "families"} the objective is additive and
#' unconstrained, so the top \code{n_families} families of the triangle are the
#' exact optimum; \code{ohs_method} is irrelevant there. Under
#' \code{pair_budget = "parents"} ("at most \code{top_founders} distinct
#' parents, families free") the optimum is the best parent \emph{set} and every
#' family among it — \code{sum(M[S, S])}, the parents-first objective — so this
#' budget runs \code{\link{ohs_selection}} with the chosen \code{ohs_method}
#' (\code{"GA_pairwise"} or \code{"Greedy"}) and returns the traversal of its
#' set, ranked by \code{M}. \code{"GA_coverage"} scores a parent set by block
#' coverage and has no family reading; it is refused. A finite \code{max_use}
#' under \code{"parents"} falls back to the best-first greedy over families.
#'
#' @section Elite subsetting:
#' \code{elite_ratio} and \code{elite_level} are the knobs of
#' \code{\link{ohs_selection}}. When the set engine runs (\code{"parents"}
#' budget, no usage cap) they are passed to it. Under the family greedy:
#' \code{elite_level = "pair"} leaves the ranking unchanged — the families
#' taken best-first \emph{are} the elite subset of the pool, whatever the
#' ratio; \code{elite_level = "block"} scores the family set by
#' \code{sum_b} of the \code{q = ceiling(elite_ratio * n_families)} largest
#' per-block ceilings over its families, re-scored after every pick (the
#' family whose blocks add most to the set's "ultimate genotype"), over the
#' \code{shortlist_families} best families of the pool by \code{M}. At
#' \code{elite_ratio = 1} the block level ranks exactly as the pair level.
#' The budget, and so the number of families made, is the same at every ratio.
#'
#' @param localgebv Haplotype-level local GEBV (see \code{\link{ohs_selection}}).
#' @param pair_budget \code{"families"} or \code{"parents"} (see
#'   \code{\link{pair_selection_budget}}).
#' @param top_founders The parents-first reference the budget is derived from.
#' @param allow_selfing Logical. Offer self-crosses as candidate families (a
#'   self's ceiling is the parent's haplotypes doubled, top-\code{ploidy} per
#'   block). Default \code{FALSE}.
#' @param ohs_method Set engine used under \code{pair_budget = "parents"}.
#' @param elite_ratio,elite_level Elite subsetting (see the section above and
#'   \code{\link{ohs_selection}}). Defaults \code{1} and \code{"pair"}.
#' @param shortlist_families Block level only: candidate families held with
#'   their per-block values, the best by \code{M}. \code{NULL} (default) =
#'   the larger of \code{4 * n_families} and \code{100}.
#' @param max_use Per-parent usage cap (families a parent may appear in).
#'   Default \code{Inf}.
#' @param overshoot Under \code{"parents"}, what to do with a family that would
#'   exceed the name limit: \code{"skip"} it (default) or \code{"stop"} there.
#' @param ploidy,pop_size,n_generations,mutation_rate,seed Passed to
#'   \code{\link{ohs_selection}} / the value matrix.
#'
#' @return A \code{"pair_selection"} \code{data.table}: \code{Parent1},
#'   \code{Parent2}, \code{score} (\code{-M[i, j]}), \code{rank}; attributes
#'   \code{n_families}, \code{n_parents_used}, \code{pair_budget}, \code{scheme}.
#'
#' @seealso \code{\link{gbs_pair_selection}}, \code{\link{pair_selection_budget}},
#'   \code{\link{run_stacking_simulation}}
#' @export
ohs_pair_selection <- function(localgebv,
                               pair_budget   = c("families", "parents"),
                               top_founders  = 8L,
                               allow_selfing = FALSE,
                               ohs_method    = c("GA_pairwise", "Greedy", "GA_coverage"),
                               elite_ratio   = 1,
                               elite_level   = c("pair", "block"),
                               shortlist_families = NULL,
                               max_use       = Inf,
                               overshoot     = c("skip", "stop"),
                               ploidy        = NULL,
                               pop_size      = 100L,
                               n_generations = 200L,
                               mutation_rate = 0.1,
                               seed          = 42L) {
  ohs_method  <- match.arg(ohs_method)
  elite_level <- match.arg(elite_level)
  overshoot   <- match.arg(overshoot)
  .ohs_check_elite(elite_ratio)
  if (ohs_method == "GA_coverage")
    stop("ohs_pair_selection: 'GA_coverage' scores a parent set by block coverage ",
         "and has no family reading; use ohs_method = \"GA_pairwise\" or \"Greedy\".")
  budget <- pair_selection_budget(top_founders, pair_budget, allow_selfing)

  if (budget$pair_budget == "parents" && !is.finite(max_use)) {
    scorer  <- .pair_scorer_ohs(localgebv, allow_selfing = allow_selfing, ploidy = ploidy)
    parents <- ohs_selection(localgebv, n_parents = budget$n_parents_max,
                             method = ohs_method, ploidy = ploidy,
                             elite_ratio = elite_ratio, elite_level = elite_level,
                             pop_size = pop_size, n_generations = n_generations,
                             mutation_rate = mutation_rate, seed = seed)
    idx  <- match(parents, scorer$samples)
    cand <- .pair_candidates(length(scorer$samples), allow_selfing, idx = idx)
    sc   <- scorer$score(NULL, cand)
    o    <- order(sc)
    dt   <- data.table(Parent1 = scorer$samples[cand$i[o]],
                       Parent2 = scorer$samples[cand$j[o]],
                       score = sc[o], rank = seq_along(o))
    return(.as_pair_selection(dt, "OHS", budget, length(parents)))
  }
  scorer <- if (elite_level == "block" && elite_ratio < 1)
    .pair_scorer_ohs_block(localgebv, budget, elite_ratio = elite_ratio,
                           allow_selfing = allow_selfing, ploidy = ploidy,
                           shortlist_families = shortlist_families)
  else
    .pair_scorer_ohs(localgebv, allow_selfing = allow_selfing, ploidy = ploidy)
  .pairs_greedy(scorer, budget, allow_selfing = allow_selfing,
                max_use = max_use, overshoot = overshoot)
}


# ---------------------------------------------------------------------------
# GBS pairs (exported)
# ---------------------------------------------------------------------------

#' Genotype Building Strategy over families (direct pair selection)
#'
#' Chooses the crosses themselves, scoring every candidate family by the
#' ideal-genotype objective of \code{\link{gbs_selection}} applied to the
#' family's expected offspring: at each target state \code{t} a family
#' \code{(i, j)} transmits \code{q_t = (HapN[i, t] + HapN[j, t]) / (2 * ploidy)}
#' (exact under the block crossing model). Families are added one at a time
#' (a sequential greedy, re-scored after every pick) until the budget of
#' \code{\link{pair_selection_budget}} stops it.
#'
#' @section Two levels:
#' \describe{
#'   \item{\code{gbs_pair_level = "population"} (default)}{A family is scored by
#'     what it adds to the offspring \emph{population} being built: the marginal
#'     gain of the set objective \code{F(P) = sum_t w_t f_t(p_t(P))} at the
#'     pool frequency \code{p_t(P)}, the mean of \code{q_t} over the chosen
#'     families. When \code{P} is the full traversal of a parent set \code{S}
#'     this frequency equals the parents-first \code{p_t(S)}, so the
#'     parents-first protocol is the special case "every parent used
#'     \code{n - 1} times". \code{gbs_objective = "candidate"} instead ranks
#'     families by Kemper's need score at the current pool frequency, which is
#'     separable (\code{= s_i + s_j}: the two most-needed parents).}
#'   \item{\code{gbs_pair_level = "family"}}{A family is scored by its own
#'     distance to the ideal, \code{-sum_t w_t f_t(q_t)}: the within-family
#'     composition an offspring actually inherits. States still absent from the
#'     pool being built take the penalty, so families still complement each
#'     other on presence, not on frequency.}
#' }
#'
#' @section Frequency projection and the absent-state penalty:
#' \code{pair_freq = "projected"} divides the accumulated \code{q_t} by the
#' final family count (Kemper's projection to the final set size; the default
#' under \code{pair_budget = "families"}, where that count is known);
#' \code{"running"} uses the mean over the families chosen so far (default under
#' \code{"parents"}). With \code{pair_penalty = "loss"} (default) a state's term
#' is \code{(1 - L_t) * f_t + L_t * beta_t}, where \code{L_t = prod_f a_f(t)^M} is
#' the exact probability that none of the \code{num_offspring = M} progeny of
#' the chosen families carries \code{t} (\code{a_f(t)} = the probability that an
#' offspring of family \code{f} carries no copy) — the parents-first risk term
#' made exact for a chosen design, and the reason a state carried by a single
#' heterozygous family keeps a residual pull toward a second carrier.
#' \code{"hard"} applies \code{beta_t} only while the state is absent from the
#' pool. The penalty scale itself (\code{penalty_mode = "dynamic"}) is sized with
#' \code{top_founders} parents, so it equals the one of the parents-first run
#' the pair run is matched to.
#'
#' @section Candidate space:
#' \code{pair_pool = "shortlist"} (default) scores, at every step, the
#' \code{shortlist_k} parents (default \code{4 * top_founders}) whose single
#' gamete would help most, and considers the families among them;
#' \code{"all"} scores every family of the pool at every step (the exact greedy;
#' memory \code{~ choose(N, 2) x targets} per step).
#'
#' @inheritParams gbs_selection
#' @inheritParams ohs_pair_selection
#' @param num_offspring Offspring per family, \code{M} in the loss probability.
#'   Default \code{4}.
#' @param gbs_objective \code{"set"} (default for families, both presets) or
#'   \code{"candidate"} (population level only).
#' @param gbs_pair_level \code{"population"} or \code{"family"} (see above).
#' @param pair_freq \code{NULL} (by budget), \code{"projected"} or \code{"running"}.
#' @param pair_penalty \code{"loss"} (default) or \code{"hard"}.
#' @param pair_pool \code{"shortlist"} (default) or \code{"all"}.
#' @param shortlist_k Parents in the shortlist; \code{NULL} = \code{4 * top_founders}.
#'
#' @return A \code{"pair_selection"} \code{data.table}: \code{Parent1},
#'   \code{Parent2}, \code{score} (lower = better; at the population level the
#'   negated marginal gain of \code{F(P)}), \code{rank}; attributes
#'   \code{n_families}, \code{n_parents_used}, \code{pair_budget}, \code{scheme}.
#'
#' @references Kemper, K.E., Bowman, P.J., Pryce, J.E., Hayes, B.J., Goddard,
#'   M.E. (2012). Long-term selection strategies for complex traits using
#'   high-density genetic markers. \emph{Journal of Dairy Science}, 95(8),
#'   4646–4656.
#'
#' @seealso \code{\link{gbs_selection}}, \code{\link{ohs_pair_selection}},
#'   \code{\link{pair_selection_budget}}, \code{\link{run_stacking_simulation}}
#' @export
gbs_pair_selection <- function(localgebv,
                               ideal_genotypes   = NULL,
                               pair_budget       = c("families", "parents"),
                               top_founders      = 8L,
                               num_offspring     = 4L,
                               allow_selfing     = FALSE,
                               gbs_method        = c("ARCSIN", "LOG"),
                               gbs_preset        = c("haploStack", "kemper"),
                               gbs_objective     = c("set", "candidate"),
                               penalty_mode      = c("dynamic", "fixed"),
                               gbs_penalty       = -10,
                               penalty_scale     = 10,
                               gbs_effect_weight = 0,
                               gbs_pair_level    = c("population", "family"),
                               pair_freq         = NULL,
                               pair_penalty      = c("loss", "hard"),
                               pair_pool         = c("shortlist", "all"),
                               shortlist_k       = NULL,
                               max_use           = Inf,
                               overshoot         = c("skip", "stop"),
                               select_priority   = c("high", "low")) {
  gbs_method      <- match.arg(gbs_method)
  gbs_preset      <- match.arg(gbs_preset)
  gbs_objective   <- match.arg(gbs_objective)
  select_priority <- match.arg(select_priority)
  overshoot       <- match.arg(overshoot)

  # The preset fills the penalty knobs the caller did not pass (explicit wins);
  # the objective is the pairs default "set" unless passed.
  spec <- .gbs_apply_preset(gbs_preset, list(
    penalty_mode      = if (missing(penalty_mode))      NULL else match.arg(penalty_mode),
    gbs_penalty       = if (missing(gbs_penalty))       NULL else gbs_penalty,
    penalty_scale     = if (missing(penalty_scale))     NULL else penalty_scale,
    gbs_effect_weight = if (missing(gbs_effect_weight)) NULL else gbs_effect_weight))
  .gbs_check_penalty_args(spec$gbs_penalty, spec$penalty_mode, spec$penalty_scale)

  budget <- pair_selection_budget(top_founders, pair_budget, allow_selfing)
  if (is.null(ideal_genotypes))
    ideal_genotypes <- .construct_ideal_genotype(localgebv)
  mm <- .gbs_hapn_matrix(localgebv, ideal_genotypes)
  if (length(mm$samples) < 2L && !isTRUE(allow_selfing))
    stop("gbs_pair_selection: fewer than two candidates in 'localgebv'.")
  ctx <- .gbs_context(mm, budget$n_parents_ref, gbs_method, spec)

  scorer <- .pair_scorer_gbs(mm, ctx, budget, n_offspring = num_offspring,
                             gbs_pair_level = gbs_pair_level,
                             gbs_objective = gbs_objective,
                             pair_freq = pair_freq, pair_penalty = pair_penalty,
                             pair_pool = pair_pool, shortlist_k = shortlist_k,
                             select_priority = select_priority,
                             allow_selfing = allow_selfing)
  .pairs_greedy(scorer, budget, allow_selfing = allow_selfing,
                max_use = max_use, overshoot = overshoot)
}
