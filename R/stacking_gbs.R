# =============================================================================
# stacking_gbs.R
#
# Genotype Building Strategy (GBS) — drive a selected set's block-haplotype
# frequencies toward an "ideal genotype" (Kemper et al. 2012). Ported from the
# origin OCS_fitness_function / construct_ultimate_genotypes
# (Origin_workflow_example/_Haplotype_Stacking/_haplotype_combin_function_GA.R),
# renamed GBS to separate it from the project's kinship-constrained
# ocs_selection().
#
# 0.7.0 (dev/UPGRADE_PLAN_opt_kemper.md §4): Kemper's two modes run as written
# under `gbs_preset = "kemper"` — ARCSIN or LOG deviation, sequential greedy,
# fixed -100 (= base -10 x scale 10) for an absent target — and haploStack's
# adjustments are options layered on top: a bounded risk x hope penalty scale
# of the pool frequency (`penalty_mode = "dynamic"`, the default), effect
# weights, a set objective and a GA engine.
#
# Functions:
#   gbs_selection()            — exported: parent selection toward an ideal genotype
#   construct_ideal_genotype() — exported: ideal_genotypes(Block, Effect, CFreq, TargetFreq)
#   tune_gbs_penalty()         — exported: grid the penalty scale against D2Ideal
#   .gbs_deviation()           — internal: d_t(p) per link — the ONLY place the two
#                                formulas live (ARCSIN asin√p − asin√τ; LOG τ·log p)
#   .gbs_set_term()            — internal: f_t(p) (ARCSIN −d²; LOG d)
#   .gbs_need_score()          — internal: candidate-objective scores
#   .gbs_set_objective()       — internal: F(S)
#   .gbs_marginal_scores()     — internal: −(F(S ∪ c) − F(S)) for every candidate
#   .gbs_greedy()              — internal: sequential greedy (either objective)
#   .gbs_ga_set()              — internal: GA over subsets on F(S)
#   .gbs_penalty_scale()       — internal: s_t (fixed constant or risk × hope curve)
#   .resolve_gbs_penalty()     — internal: β_t = gbs_penalty · s_t (a vector)
#   .gbs_apply_preset()        — internal: preset → argument defaults
#   .gbs_hapn_matrix()         — internal: HapN[candidates × targets], pop_f, N
#   .construct_ideal_genotype()— internal: ideal table (custom / argmax / ratio)
#   .distance_to_ideal()       — internal: D2Ideal (ARCSIN squared distance; LOG KL)
#   .localgebv_long()          — internal: wide haplotype GEBV -> long form
#
# Input (all functions): haplotype-level local GEBV from
#   compute_profile_localgebv(geno_meta, effect_df, level = "haplotype")
# i.e. a wide data.table "Block" + one column per "<Sample>.<haplotype_idx>"
# (the same contract ohs_selection() consumes).
# =============================================================================


# ---------------------------------------------------------------------------
# Internal — wide haplotype GEBV -> long (Sample, Haplotype, Block, Effect)
# ---------------------------------------------------------------------------

# Melt the wide "Block" + "<Sample>.<idx>" table to long form. `Effect` is the
# haplotype's local-GEBV value at the block; `Haplotype` is the 1..ploidy index.
.localgebv_long <- function(localgebv) {
  dt <- as.data.table(localgebv)
  if (!"Block" %in% names(dt))
    stop(".localgebv_long: 'localgebv' must have a 'Block' column.")
  sample_cols <- setdiff(names(dt), "Block")
  if (length(sample_cols) == 0L)
    stop(".localgebv_long: 'localgebv' has no sample columns.")

  long <- data.table::melt(dt, id.vars = "Block",
                           measure.vars = sample_cols,
                           variable.name = "SampleHap",
                           value.name = "Effect",
                           variable.factor = FALSE)
  long[, Sample    := sub("\\.[0-9]+$", "", SampleHap)]
  long[, Haplotype := as.integer(sub("^.*\\.([0-9]+)$", "\\1", SampleHap))]
  long[, Effect    := as.numeric(Effect)]
  long[, SampleHap := NULL]
  long[]
}


# ---------------------------------------------------------------------------
# Internal — ideal genotype construction (origin construct_ultimate_genotypes)
# ---------------------------------------------------------------------------

# Build the ideal genotype target table from haplotype-level local GEBV.
# Returns data.table(Block, Effect, CFreq, TargetFreq); attr "ideal_phenotype"
# holds sum(Effect * TargetFreq * ploidy).
.construct_ideal_genotype <- function(localgebv,
                                      method        = c("custom", "argmax", "ratio"),
                                      select_top_n  = 10L,
                                      k_penalty     = 0.5,
                                      select_priority = c("high", "low")) {
  method          <- match.arg(method)
  select_priority <- match.arg(select_priority)

  materials <- .localgebv_long(localgebv)
  ploidy    <- max(materials$Haplotype, na.rm = TRUE)
  n_pop     <- length(unique(materials$Sample))

  # Block-wise state table: frequency of each distinct effect-state per block.
  effect.info <- materials[, .(CFreq = .N / n_pop / ploidy), by = .(Block, Effect)]

  if (method == "argmax") {
    # Target the single best (or worst, for "low") state per block.
    pick_fun <- if (select_priority == "high") max else min
    ideal <- effect.info[, .(Effect, CFreq,
                             TargetFreq = as.numeric(Effect == pick_fun(Effect))),
                         by = Block]

  } else if (method == "ratio") {
    ideal <- effect.info[, {
      shifted <- Effect - min(Effect) + 1e-10
      .(Effect = Effect, CFreq = CFreq, TargetFreq = shifted / sum(shifted))
    }, by = Block]

  } else {  # custom — Kemper diversity-penalised target
    maxN  <- select_top_n * ploidy
    parts <- vector("list", length(unique(effect.info$Block)))
    blocks <- unique(effect.info$Block)
    for (bi in seq_along(blocks)) {
      df_b <- effect.info[Block == blocks[bi]][order(-Effect)]
      eff  <- df_b$Effect
      nrb  <- length(eff)
      tf   <- numeric(nrb)
      cap  <- min(maxN, nrb)
      for (nb in seq_len(cap)) {
        xn    <- if (nb + 1L <= nrb) eff[nb + 1L] else -Inf
        xmean <- mean(eff[seq_len(nb)])
        xp    <- xmean - (-2 * k_penalty / nb)
        if (xn < xp || nb == cap) {
          t <- (eff[seq_len(nb)] - xmean) / 2 / k_penalty - 1 / nb
          t <- t - min(t) + 1e-10
          tf[seq_len(nb)] <- t / sum(t)
          break
        }
      }
      df_b[, TargetFreq := tf]
      parts[[bi]] <- df_b
    }
    ideal <- data.table::rbindlist(parts)
  }

  ideal <- ideal[!is.na(TargetFreq) & TargetFreq > 0]
  ideal <- ideal[order(Block), .(Block, Effect, CFreq, TargetFreq)]

  ideal_phenotype <- sum(ideal$Effect * ideal$TargetFreq * ploidy, na.rm = TRUE)
  data.table::setattr(ideal, "ideal_phenotype", ideal_phenotype)
  data.table::setattr(ideal, "ploidy", ploidy)
  ideal[]
}


# ---------------------------------------------------------------------------
# Internal — per-target HapN matrix (candidates x targets) for a sample set
# ---------------------------------------------------------------------------

# Returns list(HapN = matrix[samples x targets], samples, targets, ploidy,
# pop_f, n_pool). HapN[s, t] = number of haplotypes of sample s carrying target
# state t; pop_f[t] = frequency of state t in the whole candidate pool (the
# input of the dynamic penalty scale — it is computed here, from the pool
# actually supplied, never read from the ideal table's CFreq, which under a
# founder-fixed ideal is the founder frequency).
.gbs_hapn_matrix <- function(localgebv, ideal_genotypes) {
  materials <- .localgebv_long(localgebv)
  ploidy    <- max(materials$Haplotype, na.rm = TRUE)
  pop.ref   <- materials[, .(HapN = .N), by = .(Sample, Block, Effect)]

  targets <- as.data.table(ideal_genotypes)[TargetFreq > 0,
                                            .(Block, Effect, TargetFreq)]
  if (nrow(targets) == 0L)
    stop(".gbs_hapn_matrix: 'ideal_genotypes' has no positive TargetFreq rows.")
  targets[, tkey := .I]

  all_samples <- sort(unique(pop.ref$Sample))
  hits <- pop.ref[targets, on = .(Block, Effect), nomatch = NULL,
                  allow.cartesian = TRUE]

  HapN <- matrix(0, nrow = length(all_samples), ncol = nrow(targets),
                 dimnames = list(all_samples, NULL))
  if (nrow(hits) > 0L)
    HapN[cbind(match(hits$Sample, all_samples), hits$tkey)] <- hits$HapN

  n_pool <- length(all_samples)
  list(HapN = HapN, samples = all_samples, targets = targets, ploidy = ploidy,
       pop_f = colSums(HapN) / (n_pool * ploidy), n_pool = n_pool)
}


# ---------------------------------------------------------------------------
# Internal — the two deviation forms (one link function per method)
# ---------------------------------------------------------------------------

# d_t(p): per-state deviation of a selected-set frequency p from its target
# tau, on the method's link scale. This is the only place the two formulas
# live; the need score, the set objective, the penalty scale and D2Ideal all
# call it. Both are "need" scores: lower = the state is more needed.
#   ARCSIN  asin(sqrt(p)) - asin(sqrt(tau))  signed, 0 at target (Kemper 2012)
#   LOG     tau * log(p)                     <= 0, most negative for rare states
#                                            (Kemper 2012; ranked ascending —
#                                            the origin ranked it descending)
# `p` and `tau` may be vectors or conformable matrices.
.gbs_deviation <- function(p, tau, gbs_method) {
  switch(gbs_method,
         ARCSIN = asin(sqrt(p)) - asin(sqrt(tau)),
         LOG    = tau * log(p),
         stop(".gbs_deviation: unknown gbs_method '", gbs_method, "'."))
}

# f_t(p): per-state term of the set objective F(S) = sum_t w_t f_t(p_t(S)),
# to be MAXIMISED. ARCSIN's deviation is signed, so the loss is its square
# (F = -D2Ideal); Kemper's LOG deviation is already the summand of a loss
# (cross-entropy sum_t tau_t log p_t = -KL(tau || p) + const), so it is used as
# is — squaring it would invent a third objective.
.gbs_set_term <- function(p, tau, gbs_method) {
  d <- .gbs_deviation(p, tau, gbs_method)
  if (gbs_method == "ARCSIN") -d^2 else d
}


# ---------------------------------------------------------------------------
# Internal — the absent-state ("diversity") penalty  beta_t = gbs_penalty * s_t
# ---------------------------------------------------------------------------

# Scale s_t >= 1 applied to the base penalty for target state t.
#   fixed    s_t = penalty_scale for every state (10 = Kemper's -100 from a -10
#            base; 1 = the origin's flat -10).
#   dynamic  s_t = 1 + (penalty_scale - 1) * risk_t * hope_t, bounded in
#            [1, penalty_scale]:
#              risk_t = (1 - pop_f_t)^(n_parents * ploidy)
#                       P(a random parental set carries no copy of t) — the same
#                       pool frequency is dangerous with few parents, safe with many;
#              hope_t = min(1, copies_in_pool_t / (tau_t * n_parents * ploidy))
#                       can the pool supply the target? 1 for any state with at
#                       least as many copies as the parental target needs; falls
#                       for the near-hopeless (1 copy when 8 are needed -> 1/8).
#            The curve rises from the origin's -10 for common states to near
#            Kemper's -100 for "a handful of copies", then declines instead of
#            running to infinity (plan #030 §4.3; audit penalty_curves_2026-09-15).
.gbs_penalty_scale <- function(penalty_mode, penalty_scale, pop_f, target_freq,
                               n_parents, ploidy, n_pool) {
  if (penalty_mode == "fixed")
    return(rep(penalty_scale, length(pop_f)))
  nk     <- n_parents * ploidy
  risk   <- (1 - pop_f)^nk
  copies <- pop_f * n_pool * ploidy
  hope   <- pmin(1, copies / (target_freq * nk))
  1 + (penalty_scale - 1) * risk * hope
}

# beta_t: the per-state penalty vector applied where a target state is absent
# from the selected set (p_t(S) == 0). Enters the score exactly as the deviation
# does — per state, times HapN[c, t] and w_t.
.resolve_gbs_penalty <- function(gbs_penalty, penalty_mode, penalty_scale, pop_f,
                                 target_freq, n_parents, ploidy, n_pool) {
  gbs_penalty * .gbs_penalty_scale(penalty_mode, penalty_scale, pop_f,
                                   target_freq, n_parents, ploidy, n_pool)
}

# Argument validation shared by gbs_selection() and the driver.
.gbs_check_penalty_args <- function(gbs_penalty, penalty_mode, penalty_scale) {
  if (!is.numeric(gbs_penalty) || length(gbs_penalty) != 1L || gbs_penalty > 0)
    stop("'gbs_penalty' must be a single non-positive number (the base penalty, ",
         "e.g. -10).")
  if (!is.numeric(penalty_scale) || length(penalty_scale) != 1L ||
      penalty_scale < 1)
    stop("'penalty_scale' must be a single number >= 1 (the fixed scale, or the ",
         "cap of the dynamic scale).")
  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# Internal — presets
# ---------------------------------------------------------------------------

# The two presets differ in ONE knob: how the absent-state penalty is scaled.
# Everything else is Kemper's genotype building as published (sequential greedy
# on the per-candidate score, no effect weight); the haploStack adjustments
# (set objective, GA engine, effect weight) are opt-in under either preset.
.GBS_PRESETS <- list(
  haploStack = list(gbs_objective = "candidate", gbs_engine = "greedy",
                    penalty_mode = "dynamic", gbs_penalty = -10, penalty_scale = 10,
                    gbs_effect_weight = 0),
  kemper     = list(gbs_objective = "candidate", gbs_engine = "greedy",
                    penalty_mode = "fixed",   gbs_penalty = -10, penalty_scale = 10,
                    gbs_effect_weight = 0)
)

# Fill every NULL entry of `supplied` from the preset; explicit values win.
.gbs_apply_preset <- function(gbs_preset, supplied = list()) {
  base <- .GBS_PRESETS[[gbs_preset]]
  if (is.null(base))
    stop(".gbs_apply_preset: unknown gbs_preset '", gbs_preset, "'.")
  for (nm in names(supplied))
    if (!is.null(supplied[[nm]])) base[[nm]] <- supplied[[nm]]
  base
}


# ---------------------------------------------------------------------------
# Internal — scoring context and the two objectives
# ---------------------------------------------------------------------------

# Everything the greedy / GA need, computed once per call.
.gbs_context <- function(mm, n_parents, gbs_method, spec) {
  tf   <- mm$targets$TargetFreq
  beta <- .resolve_gbs_penalty(spec$gbs_penalty, spec$penalty_mode,
                               spec$penalty_scale, mm$pop_f, tf,
                               n_parents, mm$ploidy, mm$n_pool)
  list(HapN    = mm$HapN,
       samples = mm$samples,
       tf      = tf,
       w       = 1 + abs(mm$targets$Effect) * spec$gbs_effect_weight,
       beta    = beta,
       denom   = n_parents * mm$ploidy,     # SelFreq projected to the final set size
       method  = gbs_method)
}

# Candidate objective (Kemper's sequential greedy): score(c) =
# sum_t HapN[c,t] * w_t * D_t with D_t = d_t(p_t(S)) where p_t(S) > 0, else
# beta_t. Lower = more needed.
.gbs_need_score <- function(ctx, selfreq) {
  D <- ifelse(selfreq > 0,
              .gbs_deviation(pmax(selfreq, 0), ctx$tf, ctx$method),
              ctx$beta)
  as.vector(ctx$HapN %*% (D * ctx$w))
}

# Set objective F(S) for the selected-set frequency vector `p` (maximise).
.gbs_set_objective <- function(ctx, p) {
  term <- ifelse(p > 0, .gbs_set_term(pmax(p, 0), ctx$tf, ctx$method), ctx$beta)
  sum(ctx$w * term)
}

# Marginal-gain scores for every candidate: -(F(S u {c}) - F(S)), so that, as
# for the candidate objective, the minimum is the best pick. Vectorised over
# candidates: P[c, t] = p_t(S) + HapN[c, t] / denom.
.gbs_marginal_scores <- function(ctx, selfreq, F_S) {
  N  <- nrow(ctx$HapN)
  Tn <- ncol(ctx$HapN)
  P    <- ctx$HapN / ctx$denom + matrix(selfreq,  N, Tn, byrow = TRUE)
  TF   <- matrix(ctx$tf,   N, Tn, byrow = TRUE)
  BETA <- matrix(ctx$beta, N, Tn, byrow = TRUE)
  term <- ifelse(P > 0, .gbs_set_term(pmax(P, 0), TF, ctx$method), BETA)
  F_new <- as.vector(term %*% ctx$w)
  -(F_new - F_S)
}

# Sequential greedy: pick one, recompute the selected-set frequency, re-score.
# `objective = "candidate"` scores by .gbs_need_score() (Kemper); `"set"` by
# the marginal gain of F(S). `take_min` = TRUE for select_priority "high".
.gbs_greedy <- function(ctx, n_parents, objective, take_min) {
  Tn <- ncol(ctx$HapN)
  selected <- integer(0L)
  for (k in seq_len(n_parents)) {
    selfreq <- if (length(selected) == 0L) rep(0, Tn)
               else colSums(ctx$HapN[selected, , drop = FALSE]) / ctx$denom
    scores <- if (objective == "candidate")
      .gbs_need_score(ctx, selfreq)
    else
      .gbs_marginal_scores(ctx, selfreq, .gbs_set_objective(ctx, selfreq))
    scores[selected] <- if (take_min) Inf else -Inf
    best <- if (take_min) which.min(scores) else which.max(scores)
    selected <- c(selected, best)
  }
  ctx$samples[selected]
}

# GA over subsets of size n_parents on the set objective F(S), seeded with the
# set-objective greedy solution as one elite member (so it never does worse
# than the greedy). Tournament selection, union crossover, single-point
# mutation, 10 % elitism — the OHS GA recipe with a pluggable fitness. The
# global RNG stream is restored on exit so the driver's crossing step is not
# reseeded by a selection call.
.gbs_ga_set <- function(ctx, n_parents, take_min, ga_control = list(),
                        seed = 42L) {
  ctl <- utils::modifyList(list(pop_size = 100L, n_generations = 200L,
                                mutation_rate = 0.1), ga_control)
  n   <- nrow(ctx$HapN)
  sgn <- if (take_min) 1 else -1
  fit <- function(idx)
    sgn * .gbs_set_objective(ctx, colSums(ctx$HapN[idx, , drop = FALSE]) / ctx$denom)

  greedy0 <- sort(match(.gbs_greedy(ctx, n_parents, "set", take_min), ctx$samples))

  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old_seed)) rm(list = ".Random.seed", envir = globalenv())
    else assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)
  set.seed(seed)

  pop_size <- max(2L, as.integer(ctl$pop_size))
  pop <- c(list(greedy0),
           replicate(pop_size - 1L, sort(sample.int(n, n_parents)), simplify = FALSE))
  fv    <- vapply(pop, fit, numeric(1L))
  best  <- pop[[which.max(fv)]]
  bestf <- max(fv)
  elite <- max(1L, ceiling(pop_size * 0.1))

  for (g in seq_len(ctl$n_generations)) {
    ord    <- order(fv, decreasing = TRUE)
    newpop <- vector("list", pop_size)
    newpop[seq_len(elite)] <- pop[ord[seq_len(elite)]]
    for (kk in (elite + 1L):pop_size) {
      p1 <- pop[[.ohs_tournament(fv)]]
      p2 <- pop[[.ohs_tournament(fv)]]
      child <- sort(sample(union(p1, p2), n_parents))
      if (stats::runif(1L) < ctl$mutation_rate) {
        avail <- setdiff(seq_len(n), child)
        if (length(avail) > 0L) {
          child[sample.int(n_parents, 1L)] <- sample(avail, 1L)
          child <- sort(child)
        }
      }
      newpop[[kk]] <- child
    }
    pop <- newpop
    fv  <- vapply(pop, fit, numeric(1L))
    gb  <- which.max(fv)
    if (fv[gb] > bestf) { bestf <- fv[gb]; best <- pop[[gb]] }
  }
  ctx$samples[best]
}


# ---------------------------------------------------------------------------
# GBS selection (exported)
# ---------------------------------------------------------------------------

#' Genotype Building Strategy (GBS) parent selection
#'
#' Selects \code{n_parents} candidates whose pooled block-haplotype frequencies
#' are driven toward an \emph{ideal genotype}, following the genotype-building
#' strategy of Kemper et al. (2012). By default selection is Kemper's sequential
#' greedy: a running selected-set frequency \code{SelFreq} (projected to the
#' final set size, \code{sum(HapN) / (n_parents * ploidy)}) is recomputed after
#' every pick and each remaining candidate is scored by
#' \code{sum_t HapN[c, t] * w_t * D_t}, where \code{D_t} is the deviation of
#' \code{SelFreq} from \code{TargetFreq} on the chosen link scale, or the
#' absent-state penalty when the state is not yet in the set.
#'
#' @section The two deviation forms (Kemper's modes):
#' \describe{
#'   \item{\code{"ARCSIN"}}{\code{asin(sqrt(p)) - asin(sqrt(tau))} — signed, zero
#'     at target: once a state reaches its target it stops pulling, and any
#'     other state may fill a gap left by a lost one.}
#'   \item{\code{"LOG"}}{\code{tau * log(p)} — always \eqn{\le 0}, most negative
#'     for states rare in the set: keeps pulling every remaining target, harder
#'     for the ones with the larger \code{TargetFreq}. Ranked ascending (the
#'     origin GRDC transcription ranked it descending, which inverted the
#'     order; haploStack 0.6.0 replaced the formula by a signed log-ratio —
#'     dropped in 0.7.0, see the changelog).}
#' }
#' Both are \emph{need} scores: the minimum is the most needed candidate.
#'
#' @section The absent-state ("diversity") penalty:
#' A target state absent from the selected set takes
#' \code{beta_t = gbs_penalty * s_t} in place of its deviation; it enters the
#' score per state, times \code{HapN[c, t]} and \code{w_t}. With
#' \code{penalty_mode = "fixed"}, \code{s_t = penalty_scale} for every state
#' (\code{10} gives Kemper's \code{-100} from the \code{-10} base; \code{1} the
#' origin's flat \code{-10}). With \code{penalty_mode = "dynamic"} (default),
#' \code{s_t = 1 + (penalty_scale - 1) * risk_t * hope_t}, bounded in
#' \code{[1, penalty_scale]}, where \code{risk_t = (1 - pop_f_t)^(n_parents *
#' ploidy)} is the probability that a random parental set carries no copy of
#' the state (\code{pop_f_t} = its frequency in the whole candidate pool) and
#' \code{hope_t = min(1, copies_in_pool_t / (TargetFreq_t * n_parents *
#' ploidy))} is whether the pool can supply the target at all. The scale rises
#' from the origin's constant for common states to near Kemper's for states
#' with a handful of copies, then declines for the near-hopeless, and the same
#' pool frequency is penalised more with few parents than with many.
#'
#' @section Presets:
#' \code{gbs_preset = "kemper"} runs the published method: \code{"candidate"}
#' objective, \code{"greedy"} engine, \code{penalty_mode = "fixed"},
#' \code{penalty_scale = 10}, \code{gbs_penalty = -10} (i.e. \code{-100}),
#' \code{gbs_effect_weight = 0}. \code{"haploStack"} (default) is identical
#' except \code{penalty_mode = "dynamic"}. Any argument passed explicitly
#' overrides its preset value.
#'
#' @param localgebv Haplotype-level local GEBV: a wide \code{data.frame}/
#'   \code{data.table} with a \code{Block} column and one column per
#'   \code{"<Sample>.<haplotype_idx>"} (the output of
#'   \code{compute_profile_localgebv(level = "haplotype")}).
#' @param ideal_genotypes Target table with columns \code{Block}, \code{Effect},
#'   \code{CFreq}, \code{TargetFreq}, e.g. from
#'   \code{\link{construct_ideal_genotype}}. When \code{NULL} (default) it is
#'   constructed from \code{localgebv} via \code{construct_ideal_genotype()}.
#' @param n_parents Integer number of parents to select. Default \code{8}.
#' @param gbs_method Deviation form: \code{"ARCSIN"} (default) or \code{"LOG"} —
#'   Kemper's two modes, see the section above.
#' @param gbs_preset \code{"haploStack"} (default) or \code{"kemper"}; fills every
#'   knob below that is not passed explicitly.
#' @param gbs_objective \code{"candidate"} (Kemper's per-candidate score at the
#'   current \code{SelFreq}) or \code{"set"} (marginal gain of the set objective
#'   \code{F(S) = sum_t w_t f_t(p_t(S))}, with \code{f_t = -d_t^2} for ARCSIN
#'   and \code{f_t = d_t} — the cross-entropy summand — for LOG).
#' @param gbs_engine \code{"greedy"} (sequential) or \code{"GA"} (genetic
#'   algorithm over subsets maximising \code{F(S)}; requires
#'   \code{gbs_objective = "set"}; seeded with the greedy solution).
#' @param penalty_mode \code{"dynamic"} or \code{"fixed"}, see the penalty
#'   section.
#' @param gbs_penalty Base penalty \code{<= 0}. Preset default \code{-10}.
#' @param penalty_scale Scale \code{>= 1}: the fixed scale, or the cap of the
#'   dynamic scale. Preset default \code{10}.
#' @param gbs_effect_weight Weight in \code{[0, 1]} scaling each state's
#'   contribution by \code{w_t = 1 + |Effect| * gbs_effect_weight}. Preset
#'   default \code{0} (ignore effect magnitude).
#' @param select_priority \code{"high"} (default) takes the minimum score;
#'   \code{"low"} the maximum.
#' @param ga_control Optional list for \code{gbs_engine = "GA"}: \code{pop_size}
#'   (100), \code{n_generations} (200), \code{mutation_rate} (0.1).
#' @param seed Seed for the GA engine (the global RNG stream is restored on
#'   exit). Default \code{42L}.
#'
#' @return A character vector of selected sample names (length \code{n_parents},
#'   or all candidates when \code{n_parents} exceeds the pool).
#'
#' @references Kemper, K.E., Bowman, P.J., Pryce, J.E., Hayes, B.J., Goddard,
#'   M.E. (2012). Long-term selection strategies for complex traits using
#'   high-density genetic markers. \emph{Journal of Dairy Science}, 95(8),
#'   4646–4656.
#'
#' @seealso \code{\link{construct_ideal_genotype}}, \code{\link{tune_gbs_penalty}},
#'   \code{\link{ohs_selection}}, \code{\link{ocs_selection}}
#' @export
gbs_selection <- function(localgebv,
                          ideal_genotypes   = NULL,
                          n_parents         = 8L,
                          gbs_method        = c("ARCSIN", "LOG"),
                          gbs_preset        = c("haploStack", "kemper"),
                          gbs_objective     = c("candidate", "set"),
                          gbs_engine        = c("greedy", "GA"),
                          penalty_mode      = c("dynamic", "fixed"),
                          gbs_penalty       = -10,
                          penalty_scale     = 10,
                          gbs_effect_weight = 0,
                          select_priority   = c("high", "low"),
                          ga_control        = list(),
                          seed              = 42L) {
  gbs_method      <- match.arg(gbs_method)
  gbs_preset      <- match.arg(gbs_preset)
  select_priority <- match.arg(select_priority)

  # Preset fills every knob the caller did not pass; explicit arguments win.
  spec <- .gbs_apply_preset(gbs_preset, list(
    gbs_objective     = if (missing(gbs_objective))     NULL else match.arg(gbs_objective),
    gbs_engine        = if (missing(gbs_engine))        NULL else match.arg(gbs_engine),
    penalty_mode      = if (missing(penalty_mode))      NULL else match.arg(penalty_mode),
    gbs_penalty       = if (missing(gbs_penalty))       NULL else gbs_penalty,
    penalty_scale     = if (missing(penalty_scale))     NULL else penalty_scale,
    gbs_effect_weight = if (missing(gbs_effect_weight)) NULL else gbs_effect_weight))
  .gbs_check_penalty_args(spec$gbs_penalty, spec$penalty_mode, spec$penalty_scale)
  if (spec$gbs_engine == "GA" && spec$gbs_objective != "set")
    stop("gbs_selection: gbs_engine = \"GA\" optimises a set objective — pass ",
         "gbs_objective = \"set\".")

  if (is.null(ideal_genotypes))
    ideal_genotypes <- .construct_ideal_genotype(localgebv)

  mm <- .gbs_hapn_matrix(localgebv, ideal_genotypes)
  if (n_parents >= length(mm$samples)) return(mm$samples)

  ctx      <- .gbs_context(mm, n_parents, gbs_method, spec)
  take_min <- select_priority == "high"

  if (spec$gbs_engine == "GA")
    .gbs_ga_set(ctx, n_parents, take_min, ga_control, seed)
  else
    .gbs_greedy(ctx, n_parents, spec$gbs_objective, take_min)
}


#' Construct an ideal genotype target table
#'
#' Builds the block-wise ideal genotype consumed by \code{\link{gbs_selection}}
#' (and the D2Ideal metric): for each block it assigns a \code{TargetFreq} to the
#' distinct haplotype effect-states. Ported from the origin
#' \code{construct_ultimate_genotypes}.
#'
#' @param localgebv Haplotype-level local GEBV (see \code{\link{gbs_selection}}).
#' @param method \code{"custom"} (default; Kemper diversity-penalised),
#'   \code{"argmax"} (all weight on the best state per block), or \code{"ratio"}
#'   (weight proportional to shifted effect).
#' @param select_top_n Integer; caps the number of states retained per block in
#'   the \code{"custom"} method (multiplied by ploidy). Default \code{10}.
#' @param k_penalty Diversity penalty for the \code{"custom"} method (smaller =
#'   more concentrated targets). Default \code{0.5}.
#' @param select_priority \code{"high"} (default) or \code{"low"} — for
#'   \code{"argmax"}, which extreme effect-state to target.
#'
#' @return A \code{data.table} with columns \code{Block}, \code{Effect},
#'   \code{CFreq}, \code{TargetFreq}. The attribute \code{"ideal_phenotype"}
#'   holds \code{sum(Effect * TargetFreq * ploidy)}.
#'
#' @seealso \code{\link{gbs_selection}}
#' @export
construct_ideal_genotype <- function(localgebv,
                                     method          = c("custom", "argmax", "ratio"),
                                     select_top_n    = 10L,
                                     k_penalty       = 0.5,
                                     select_priority = c("high", "low")) {
  .construct_ideal_genotype(localgebv, method = match.arg(method),
                            select_top_n = select_top_n, k_penalty = k_penalty,
                            select_priority = match.arg(select_priority))
}


# ---------------------------------------------------------------------------
# Internal — distance to ideal (pooled set-level; the D2Ideal metric)
# ---------------------------------------------------------------------------

# Pooled deviation of a sample set's block-haplotype frequencies from the ideal.
# Lower = closer to the ideal genotype; zero at target.
#   ARCSIN  sum_t w_t (asin sqrt(p_t) - asin sqrt(tau_t))^2
#   LOG     sum_t w_t tau_t (log tau_t - log p_t) = KL(tau || p), the divergence
#           whose negative is the LOG set objective; a state absent from the set
#           is floored at half a copy, 1 / (2 * |S| * ploidy), so the loss is
#           finite and on the data's scale.
# `samples = NULL` pools over every sample in `localgebv`.
.distance_to_ideal <- function(localgebv, ideal_genotypes,
                               samples          = NULL,
                               gbs_method       = c("ARCSIN", "LOG"),
                               gbs_effect_weight = 0) {
  gbs_method <- match.arg(gbs_method)
  mm     <- .gbs_hapn_matrix(localgebv, ideal_genotypes)
  HapN   <- mm$HapN
  ploidy <- mm$ploidy
  tf     <- mm$targets$TargetFreq
  w      <- 1 + abs(mm$targets$Effect) * gbs_effect_weight

  if (is.null(samples)) samples <- mm$samples
  samples <- intersect(samples, mm$samples)
  if (length(samples) == 0L) return(NA_real_)

  selfreq <- colSums(HapN[samples, , drop = FALSE]) / (length(samples) * ploidy)

  term <- if (gbs_method == "ARCSIN") {
    .gbs_deviation(selfreq, tf, "ARCSIN")^2
  } else {
    eps <- 1 / (2 * length(samples) * ploidy)
    p   <- pmax(selfreq, eps)
    tf * (log(tf) - log(p))
  }
  sum(term * w)
}


# ---------------------------------------------------------------------------
# Penalty tuning (exported)
# ---------------------------------------------------------------------------

#' Tune the GBS absent-state penalty scale
#'
#' Data-optimises the \code{\link{gbs_selection}} penalty by searching over a
#' grid of \code{penalty_scale} values (the fixed scale under
#' \code{penalty_mode = "fixed"}, the cap of the dynamic scale under
#' \code{"dynamic"}) with the base \code{gbs_penalty} held fixed, scoring each
#' resulting parent set by its distance to the ideal genotype
#' (\code{\link{.distance_to_ideal}}, lower is better). Scores are averaged over
#' \code{n_rep} candidate-pool subsamples to avoid overfitting a single
#' realisation.
#'
#' @param localgebv Haplotype-level local GEBV (see \code{\link{gbs_selection}}).
#' @param ideal_genotypes Target table (see \code{\link{gbs_selection}}); built
#'   from \code{localgebv} when \code{NULL}.
#' @param n_parents Integer number of parents. Default \code{8}.
#' @param gbs_method \code{"ARCSIN"} (default) or \code{"LOG"}.
#' @param penalty_mode \code{"dynamic"} (default) or \code{"fixed"}.
#' @param gbs_penalty Base penalty held fixed across the grid. Default \code{-10}.
#' @param grid Numeric vector of \code{penalty_scale} values (\code{>= 1}) to
#'   evaluate. \code{NULL} (default) uses \code{c(1, 2, 5, 10, 20)}.
#' @param n_rep Number of subsample replicates per grid point. Default \code{5}.
#' @param subsample_frac Fraction of candidates sampled per replicate. Default
#'   \code{0.8}.
#' @param seed Random seed. Default \code{42L}.
#' @param ... Passed through to \code{\link{gbs_selection}} (e.g.
#'   \code{gbs_effect_weight}, \code{gbs_objective}, \code{select_priority}).
#'
#' @return A list with \code{best_penalty_scale}, \code{best_penalty}
#'   (\code{gbs_penalty * best_penalty_scale}: the constant under
#'   \code{"fixed"}, the most negative value the dynamic scale can reach under
#'   \code{"dynamic"}), and \code{scores}, a
#'   \code{data.table(penalty_scale, mean_score, sd_score)}.
#'
#' @seealso \code{\link{gbs_selection}}
#' @export
tune_gbs_penalty <- function(localgebv,
                             ideal_genotypes = NULL,
                             n_parents       = 8L,
                             gbs_method      = c("ARCSIN", "LOG"),
                             penalty_mode    = c("dynamic", "fixed"),
                             gbs_penalty     = -10,
                             grid            = NULL,
                             n_rep           = 5L,
                             subsample_frac  = 0.8,
                             seed            = 42L,
                             ...) {
  gbs_method   <- match.arg(gbs_method)
  penalty_mode <- match.arg(penalty_mode)
  if (is.null(grid)) grid <- c(1, 2, 5, 10, 20)
  if (any(grid < 1))
    stop("tune_gbs_penalty: every 'grid' value must be >= 1 (a penalty scale).")
  if (is.null(ideal_genotypes))
    ideal_genotypes <- .construct_ideal_genotype(localgebv)

  dt          <- as.data.table(localgebv)
  sample_cols <- setdiff(names(dt), "Block")
  all_samples <- unique(sub("\\.[0-9]+$", "", sample_cols))
  n_sub       <- max(n_parents + 1L, floor(length(all_samples) * subsample_frac))

  score_one <- function(scale_s, rep_i) {
    set.seed(seed + rep_i)
    keep_s   <- if (n_sub >= length(all_samples)) all_samples
                else sample(all_samples, n_sub)
    keep_cols <- c("Block",
                   sample_cols[sub("\\.[0-9]+$", "", sample_cols) %in% keep_s])
    sub_lg   <- dt[, keep_cols, with = FALSE]
    sel      <- gbs_selection(sub_lg, ideal_genotypes = ideal_genotypes,
                              n_parents = n_parents, gbs_method = gbs_method,
                              penalty_mode = penalty_mode,
                              gbs_penalty = gbs_penalty,
                              penalty_scale = scale_s, ...)
    .distance_to_ideal(sub_lg, ideal_genotypes, samples = sel,
                       gbs_method = gbs_method)
  }

  scores <- data.table::rbindlist(lapply(grid, function(s) {
    vals <- vapply(seq_len(n_rep), function(r) score_one(s, r), numeric(1L))
    data.table(penalty_scale = s,
               mean_score = mean(vals, na.rm = TRUE),
               sd_score   = stats::sd(vals, na.rm = TRUE))
  }))

  best_scale <- scores$penalty_scale[which.min(scores$mean_score)]

  list(best_penalty_scale = best_scale,
       best_penalty       = gbs_penalty * best_scale,
       scores             = scores)
}
