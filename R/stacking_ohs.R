# =============================================================================
# stacking_ohs.R
#
# Optimal Haplotype Selection (OHS) — "maximise next-cross gain".
#
# OHS scores a SET of parents by the genetic value their crosses could build:
# for each block, an offspring of parents (i, j) can carry the best `ploidy`
# block haplotypes drawn from the parents' combined haplotype pool. The pairwise
# value matrix
#       M[i, j] = sum_block ( sum of top-`ploidy` block values among the
#                             haplotypes of samples i and j )
# captures that achievable gain, and the subset objective sum(M[S, S]) rewards a
# parent set whose pairwise crosses collectively maximise stacking.
#
# Functions:
#   ohs_selection()              — exported: select parents by an OHS engine
#   .pairwise_block_value_matrix()— internal: the n x n matrix M
#   .ga_pairwise() / .greedy_pairwise() — internal: GA / greedy on sum(M[S,S])
#   .collapse_hap_to_geno()      — internal: haplotype -> genotype wide GEBV
#
# Engines (PRD #013, issue 01):
#   "GA_pairwise" — origin GA maximising sum(M[S,S]).
#   "GA_coverage" — project genetic_algorithm_selection() (binary block coverage).
#   "Greedy"      — new greedy marginal-gain on sum(M[S,S]) (fixes the origin's
#                   broken 2-D/3-D greedy).
#
# Input: haplotype-level local GEBV from
#   compute_profile_localgebv(geno_meta, effect_df, level = "haplotype")
# i.e. a wide data.table: "Block" + one column per "<Sample>.<haplotype_idx>".
# =============================================================================


#' Optimal Haplotype Selection of parents
#'
#' Selects \code{n_parents} candidates that maximise expected next-cross gain,
#' using one of three engines (\code{method}). \code{GA_pairwise} and
#' \code{Greedy} optimise the pairwise objective \code{sum(M[S, S])} on the
#' block-value matrix \code{\link{.pairwise_block_value_matrix}};
#' \code{GA_coverage} delegates to \code{\link{genetic_algorithm_selection}}
#' (binary block coverage) after collapsing haplotypes to genotype level.
#'
#' @param localgebv Haplotype-level local GEBV: a wide \code{data.frame}/
#'   \code{data.table} with a \code{Block} column and one column per
#'   \code{"<Sample>.<haplotype_idx>"} (the output of
#'   \code{compute_profile_localgebv(level = "haplotype")}).
#' @param n_parents Integer number of parents to select. Default \code{8}.
#' @param method One of \code{"GA_pairwise"} (default), \code{"GA_coverage"},
#'   \code{"Greedy"}.
#' @param ploidy Haplotypes per sample. \code{NULL} (default) infers it from the
#'   column names.
#' @param pop_size,n_generations,mutation_rate GA controls (used by
#'   \code{GA_pairwise} and \code{GA_coverage}).
#' @param seed Random seed. Default \code{42L}.
#' @param verbose Logical; print GA progress. Default \code{FALSE}.
#'
#' @return A character vector of selected sample names (length \code{n_parents}).
#'
#' @seealso \code{\link{ocs_selection}}, \code{\link{genetic_algorithm_selection}},
#'   \code{\link{make_cross_design}}
#' @export
ohs_selection <- function(localgebv,
                          n_parents     = 8L,
                          method        = c("GA_pairwise", "GA_coverage", "Greedy"),
                          ploidy        = NULL,
                          pop_size      = 100L,
                          n_generations = 200L,
                          mutation_rate = 0.1,
                          seed          = 42L,
                          verbose       = FALSE) {
  method <- match.arg(method)
  dt <- as.data.table(localgebv)
  if (!"Block" %in% names(dt))
    stop("ohs_selection: 'localgebv' must have a 'Block' column.")
  if (length(setdiff(names(dt), "Block")) == 0L)
    stop("ohs_selection: 'localgebv' has no sample columns.")

  if (method == "GA_coverage") {
    geno <- .collapse_hap_to_geno(dt)
    return(genetic_algorithm_selection(
      geno, n_parents = n_parents, pop_size = pop_size,
      n_generations = n_generations, mutation_rate = mutation_rate,
      seed = seed, verbose = verbose))
  }

  M <- .pairwise_block_value_matrix(dt, ploidy = ploidy)
  if (n_parents >= nrow(M)) return(rownames(M))

  if (method == "GA_pairwise")
    .ga_pairwise(M, n_parents, pop_size = pop_size,
                 n_generations = n_generations, mutation_rate = mutation_rate,
                 seed = seed, verbose = verbose)
  else
    .greedy_pairwise(M, n_parents)
}


# ---------------------------------------------------------------------------
# Internal — pairwise block-value matrix
# ---------------------------------------------------------------------------

# Build the symmetric n_samples x n_samples matrix
#   M[i,j] = sum_block ( sum of top-`ploidy` values among the block haplotype
#                        values of samples i and j ).
# `localgebv` is the wide haplotype-level table (Block + "<sample>.<idx>" cols).
.pairwise_block_value_matrix <- function(localgebv, ploidy = NULL) {
  dt <- as.data.table(localgebv)
  if (!"Block" %in% names(dt))
    stop(".pairwise_block_value_matrix: 'localgebv' needs a 'Block' column.")
  sample_cols <- setdiff(names(dt), "Block")
  if (length(sample_cols) == 0L)
    stop(".pairwise_block_value_matrix: no sample columns.")

  base    <- sub("\\.[0-9]+$", "", sample_cols)
  samples <- unique(base)
  n       <- length(samples)
  nb      <- nrow(dt)
  if (is.null(ploidy)) ploidy <- max(tabulate(match(base, samples)))

  # A[sample, block, hap] — per-sample block values sorted descending, zero-pad.
  A <- array(0, dim = c(n, nb, ploidy))
  for (s in seq_len(n)) {
    cols <- sample_cols[base == samples[s]]
    vals <- matrix(as.numeric(as.matrix(dt[, cols, with = FALSE])),
                   nrow = nb, ncol = length(cols))
    sorted <- t(apply(vals, 1L, sort, decreasing = TRUE))
    if (nb == 1L) sorted <- matrix(sorted, nrow = 1L)
    if (ncol(sorted) < ploidy)
      sorted <- cbind(sorted, matrix(0, nb, ploidy - ncol(sorted)))
    A[s, , ] <- sorted[, seq_len(ploidy), drop = FALSE]
  }

  # Off-diagonal pairwise gains only; the diagonal stays 0 (self-pairs do not
  # contribute to sum(M[S,S]) — matches the origin compute_pairwise.cpp, which
  # fills only i < j).
  M <- matrix(0, n, n, dimnames = list(samples, samples))
  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      Ai <- matrix(A[i, , ], nb, ploidy)
      for (j in (i + 1L):n) {
        Aj <- matrix(A[j, , ], nb, ploidy)
        if (ploidy == 2L) {
          # Vectorised sum of top-2 of {a1>=a2, b1>=b2} per block.
          a1 <- Ai[, 1L]; a2 <- Ai[, 2L]; b1 <- Aj[, 1L]; b2 <- Aj[, 2L]
          top <- ifelse(a1 >= b1, a1 + pmax(a2, b1), b1 + pmax(b2, a1))
          val <- sum(top)
        } else {
          comb <- cbind(Ai, Aj)
          val  <- sum(apply(comb, 1L, function(r)
            sum(sort(r, decreasing = TRUE)[seq_len(ploidy)])))
        }
        M[i, j] <- val
        M[j, i] <- val
      }
    }
  }
  M
}


# ---------------------------------------------------------------------------
# Internal — subset optimisers on sum(M[S, S])
# ---------------------------------------------------------------------------

# Greedy marginal-gain selection. Marginal gain of adding candidate c to set S:
#   M[c,c] + 2 * sum_{s in S} M[s,c].
.greedy_pairwise <- function(M, n_select) {
  n  <- nrow(M)
  nm <- rownames(M)
  if (n_select >= n) return(nm)

  d         <- diag(M)
  selected  <- integer(0L)
  remaining <- seq_len(n)

  for (k in seq_len(n_select)) {
    if (length(selected) == 0L) {
      marg <- d[remaining]
    } else {
      marg <- d[remaining] +
        2 * colSums(M[selected, remaining, drop = FALSE])
    }
    best      <- remaining[which.max(marg)]
    selected  <- c(selected, best)
    remaining <- setdiff(remaining, best)
  }
  nm[selected]
}

# Internal: one tournament winner index (size k).
.ohs_tournament <- function(fitness, k = 3L) {
  idx <- sample.int(length(fitness), min(k, length(fitness)))
  idx[which.max(fitness[idx])]
}

# Genetic algorithm maximising sum(M[S, S]).
.ga_pairwise <- function(M, n_select, pop_size = 100L, n_generations = 200L,
                         mutation_rate = 0.1, seed = 42L, verbose = FALSE) {
  set.seed(seed)
  n  <- nrow(M)
  nm <- rownames(M)
  if (n_select >= n) return(nm)

  fit <- function(idx) sum(M[idx, idx])

  pop  <- replicate(pop_size, sort(sample.int(n, n_select)), simplify = FALSE)
  fv   <- vapply(pop, fit, numeric(1L))
  best <- pop[[which.max(fv)]]
  bestf <- max(fv)
  elite <- max(1L, ceiling(pop_size * 0.1))

  for (g in seq_len(n_generations)) {
    ord    <- order(fv, decreasing = TRUE)
    newpop <- vector("list", pop_size)
    newpop[seq_len(elite)] <- pop[ord[seq_len(elite)]]

    for (kk in (elite + 1L):pop_size) {
      p1 <- pop[[.ohs_tournament(fv)]]
      p2 <- pop[[.ohs_tournament(fv)]]
      child <- sort(sample(union(p1, p2), n_select))
      if (runif(1L) < mutation_rate) {
        avail <- setdiff(seq_len(n), child)
        if (length(avail) > 0L) {
          child[sample.int(n_select, 1L)] <- sample(avail, 1L)
          child <- sort(child)
        }
      }
      newpop[[kk]] <- child
    }

    pop <- newpop
    fv  <- vapply(pop, fit, numeric(1L))
    gb  <- which.max(fv)
    if (fv[gb] > bestf) { bestf <- fv[gb]; best <- pop[[gb]] }
    if (verbose && g %% 50L == 0L)
      message("  OHS-GA gen ", g, " — best sum(M[S,S]) = ", round(bestf, 3L))
  }
  nm[best]
}


# ---------------------------------------------------------------------------
# Internal — collapse haplotype-level GEBV to genotype level
# ---------------------------------------------------------------------------

# Sum the per-sample haplotype columns to a wide Block x Sample table (the input
# genetic_algorithm_selection() / greedy_selection() expect).
.collapse_hap_to_geno <- function(localgebv) {
  dt <- as.data.table(localgebv)
  sample_cols <- setdiff(names(dt), "Block")
  base    <- sub("\\.[0-9]+$", "", sample_cols)
  samples <- unique(base)

  out <- data.table(Block = dt$Block)
  for (s in samples) {
    cols <- sample_cols[base == s]
    out[, (s) := rowSums(matrix(as.numeric(as.matrix(dt[, cols, with = FALSE])),
                                nrow = nrow(dt)), na.rm = TRUE)]
  }
  out[]
}
