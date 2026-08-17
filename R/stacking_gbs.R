# =============================================================================
# stacking_gbs.R
#
# Genotype Building Strategy (GBS) — drive a selected set's block-haplotype
# frequencies toward an "ideal genotype" (Kemper et al. 2012). Ported from the
# origin OCS_fitness_function / construct_ultimate_genotypes
# (Origin_workflow_example/_Haplotype_Stacking/_haplotype_combin_function_GA.R),
# renamed GBS to separate it from the project's kept kinship-constrained
# ocs_selection().
#
# Functions:
#   gbs_selection()           — exported: greedy per-candidate selection toward
#                               an ideal genotype (ARCSIN / LOG).
#   tune_gbs_penalty()        — exported: data-optimise the missing-haplotype
#                               penalty against a distance-to-ideal objective.
#   .construct_ideal_genotype()— internal: ideal_genotypes(Block, Effect, CFreq,
#                               TargetFreq) (custom / argmax / ratio).
#   .distance_to_ideal()      — internal: pooled deviation of a sample set's
#                               block-haplotype frequencies from TargetFreq
#                               (reused by Phase 5 as the D2Ideal metric).
#   .localgebv_long()         — internal: wide haplotype GEBV -> long form.
#
# Input (all functions): haplotype-level local GEBV from
#   compute_profile_localgebv(geno_meta, effect_df, level = "haplotype")
# i.e. a wide data.table "Block" + one column per "<Sample>.<haplotype_idx>"
# (the same contract ohs_selection() consumes).
#
# §5.2 fix: the contribution is genuinely per-candidate (recomputed running
# SelFreq after every pick), unlike the population-scalar penalty in
# ocs_selection() which cancels out of the ranking.
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

# Returns a list(HapN = matrix[samples x targets], samples, targets, ploidy).
# HapN[s, t] = number of haplotypes of sample s carrying target state t.
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

  list(HapN = HapN, samples = all_samples, targets = targets, ploidy = ploidy)
}


# Resolve the missing-haplotype penalty. NULL -> scale-relative default
# -penalty_scale * max|dev|, the worst-case non-penalty signed deviation at the
# minimum non-zero SelFreq (one copy in the selected set).
.resolve_gbs_penalty <- function(gbs_penalty, penalty_scale, gbs_method,
                                  target_freq, n_parents, ploidy) {
  if (!is.null(gbs_penalty)) return(gbs_penalty)
  min_freq <- 1 / (n_parents * ploidy)
  if (gbs_method == "ARCSIN") {
    max_abs_dev <- max(abs(asin(sqrt(min_freq)) - asin(sqrt(target_freq))))
  } else {  # LOG
    max_abs_dev <- max(abs(log(min_freq) - log(target_freq)))
  }
  -penalty_scale * max_abs_dev
}


# ---------------------------------------------------------------------------
# GBS selection (exported)
# ---------------------------------------------------------------------------

#' Genotype Building Strategy (GBS) parent selection
#'
#' Selects \code{n_parents} candidates whose pooled block-haplotype frequencies
#' are driven toward an \emph{ideal genotype}, following the genotype-building
#' strategy of Kemper et al. (2012). Selection is greedy and genuinely
#' per-candidate: a running selected-set frequency (\code{SelFreq}) is recomputed
#' after every pick and each remaining candidate is scored by how much it pulls
#' \code{SelFreq} toward the target \code{TargetFreq}.
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
#' @param gbs_method Deviation metric comparing running \code{SelFreq} to
#'   \code{TargetFreq}: \code{"ARCSIN"} (default, variance-stabilised) or
#'   \code{"LOG"}. Both use a signed deviation (negative below target); an absent
#'   target takes \code{gbs_penalty}.
#' @param gbs_penalty Penalty applied to a target state that is still absent from
#'   the selected set (\code{SelFreq == 0}). \code{NULL} (default) resolves to a
#'   scale-relative value \code{-penalty_scale * max|dev|} (interpretable and
#'   stable across datasets); pass a numeric to override (e.g. the origin
#'   \code{-10}).
#' @param penalty_scale Multiplier for the scale-relative default penalty
#'   (worst-case deviations an uncovered target costs). Default \code{2}.
#' @param gbs_effect_weight Weight in \code{[0, 1]} scaling each state's
#'   contribution by \code{(1 + |Effect| * gbs_effect_weight)}. Default \code{0}
#'   (ignore effect magnitude).
#' @param select_priority \code{"high"} (default) drives toward high-value target
#'   states; \code{"low"} reverses the pick direction.
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
                          ideal_genotypes  = NULL,
                          n_parents        = 8L,
                          gbs_method       = c("ARCSIN", "LOG"),
                          gbs_penalty      = NULL,
                          penalty_scale    = 2,
                          gbs_effect_weight = 0,
                          select_priority  = c("high", "low")) {
  gbs_method      <- match.arg(gbs_method)
  select_priority <- match.arg(select_priority)

  if (is.null(ideal_genotypes))
    ideal_genotypes <- .construct_ideal_genotype(localgebv)

  mm      <- .gbs_hapn_matrix(localgebv, ideal_genotypes)
  HapN    <- mm$HapN
  samples <- mm$samples
  targets <- mm$targets
  ploidy  <- mm$ploidy

  if (n_parents >= length(samples)) return(samples)

  penalty <- .resolve_gbs_penalty(gbs_penalty, penalty_scale, gbs_method,
                                  targets$TargetFreq, n_parents, ploidy)
  tf <- targets$TargetFreq
  w  <- 1 + abs(targets$Effect) * gbs_effect_weight

  # Both metrics use a signed deviation that is negative when the running
  # SelFreq is below TargetFreq; an absent target (SelFreq == 0) takes the
  # (negative) penalty. "high" priority recruits the most-below targets, i.e.
  # the most negative score -> take the minimum. (The origin's LOG branch
  # ranked the opposite way to ARCSIN; this rebuild makes them consistent.)
  take_min <- select_priority == "high"
  denom    <- n_parents * ploidy

  selected <- character(0L)
  for (k in seq_len(n_parents)) {
    selfreq <- if (length(selected) == 0L)
      rep(0, length(tf))
    else
      colSums(HapN[selected, , drop = FALSE]) / denom

    if (gbs_method == "ARCSIN")
      dev <- ifelse(selfreq > 0, asin(sqrt(selfreq)) - asin(sqrt(tf)), penalty)
    else
      dev <- ifelse(selfreq > 0, log(selfreq) - log(tf), penalty)

    contrib <- dev * w
    scores  <- as.vector(HapN %*% contrib)
    names(scores) <- samples
    scores[selected] <- if (take_min) Inf else -Inf

    best <- if (take_min) samples[which.min(scores)] else samples[which.max(scores)]
    selected <- c(selected, best)
  }
  selected
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
# Internal — distance to ideal (pooled set-level; D2Ideal in Phase 5)
# ---------------------------------------------------------------------------

# Pooled deviation of a sample set's block-haplotype frequencies from the ideal.
# Lower = closer to the ideal genotype. ARCSIN uses squared arcsin-scale
# deviation; LOG uses squared log-scale deviation. `samples = NULL` pools over
# every sample in `localgebv`.
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

  if (gbs_method == "ARCSIN")
    dev2 <- (asin(sqrt(selfreq)) - asin(sqrt(tf)))^2
  else
    dev2 <- (log(selfreq + 1e-10) - log(tf))^2

  sum(dev2 * w)
}


# ---------------------------------------------------------------------------
# Penalty tuning (exported)
# ---------------------------------------------------------------------------

#' Tune the GBS missing-haplotype penalty
#'
#' Data-optimises the \code{gbs_selection} penalty by searching over a grid of
#' \code{penalty_scale} multipliers and scoring each resulting parent set by its
#' distance to the ideal genotype (\code{\link{.distance_to_ideal}}, lower is
#' better). Scores are averaged over \code{n_rep} candidate-pool subsamples to
#' avoid overfitting a single realisation.
#'
#' @param localgebv Haplotype-level local GEBV (see \code{\link{gbs_selection}}).
#' @param ideal_genotypes Target table (see \code{\link{gbs_selection}}); built
#'   from \code{localgebv} when \code{NULL}.
#' @param n_parents Integer number of parents. Default \code{8}.
#' @param gbs_method \code{"ARCSIN"} (default) or \code{"LOG"}.
#' @param grid Numeric vector of \code{penalty_scale} multipliers to evaluate.
#'   \code{NULL} (default) uses \code{c(0.5, 1, 2, 3, 5, 10)}.
#' @param n_rep Number of subsample replicates per grid point. Default \code{5}.
#' @param subsample_frac Fraction of candidates sampled per replicate. Default
#'   \code{0.8}.
#' @param seed Random seed. Default \code{42L}.
#' @param ... Passed through to \code{\link{gbs_selection}} (e.g.
#'   \code{gbs_effect_weight}, \code{select_priority}).
#'
#' @return A list with \code{best_penalty_scale}, \code{best_penalty} (the
#'   resolved numeric penalty at the full-pool ideal), and \code{scores}, a
#'   \code{data.table(penalty_scale, mean_score, sd_score)}.
#'
#' @seealso \code{\link{gbs_selection}}
#' @export
tune_gbs_penalty <- function(localgebv,
                             ideal_genotypes = NULL,
                             n_parents       = 8L,
                             gbs_method      = c("ARCSIN", "LOG"),
                             grid            = NULL,
                             n_rep           = 5L,
                             subsample_frac  = 0.8,
                             seed            = 42L,
                             ...) {
  gbs_method <- match.arg(gbs_method)
  if (is.null(grid)) grid <- c(0.5, 1, 2, 3, 5, 10)
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
                              gbs_penalty = NULL, penalty_scale = scale_s, ...)
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
  best_pen   <- .resolve_gbs_penalty(
    NULL, best_scale, gbs_method,
    as.data.table(ideal_genotypes)[TargetFreq > 0, TargetFreq],
    n_parents, attr(ideal_genotypes, "ploidy") %||% 2L)

  list(best_penalty_scale = best_scale,
       best_penalty       = best_pen,
       scores             = scores)
}
