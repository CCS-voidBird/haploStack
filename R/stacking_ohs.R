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
# is the family's stacking ceiling, and the subset objective sum(M[S, S])
# rewards a parent set whose crosses collectively maximise stacking.
#
# Elite subsetting (0.7.3, `elite_ratio` / `elite_level`): instead of every
# family, the objective may count only the best `q = ceiling(ratio * P)` of the
# `P = k (k - 1) / 2` families of the set —
#   elite_level = "pair"   the q largest family ceilings M[i, j] in S
#                          (ratio -> 0: the single best cross);
#   elite_level = "block"  per block, the q largest per-block family ceilings
#                          M_b[i, j], summed over blocks (ratio -> 0: the best two
#                          haplotypes of every block from any two parents — the
#                          "ultimate genotype" of Kemper 2012 / HapSelect's
#                          Haploid_OHS).
# At ratio 1 both levels are the full-set objective (sum over every family).
#
# Functions:
#   ohs_selection()              — exported: select parents by an OHS engine
#   .pairwise_block_value_matrix()— internal: the n x n matrix M (+ the array)
#   .ohs_objective()             — internal: the set objective for a ratio/level
#   .ga_pairwise() / .greedy_pairwise() — internal: GA / greedy on the objective
#   .collapse_hap_to_geno()      — internal: haplotype -> genotype wide GEBV
#
# Engines (PRD #013, issue 01):
#   "GA_pairwise" — origin GA maximising the set objective (sum(M[S,S]) by default).
#   "GA_coverage" — project genetic_algorithm_selection() (binary block coverage).
#   "Greedy"      — greedy marginal gain on the set objective, seeded with the
#                   best row sum of M (0.7.3; the origin's first pick was the
#                   first sample column, whatever the data).
#
# Input: haplotype-level local GEBV from
#   compute_profile_localgebv(geno_meta, effect_df, level = "haplotype")
# i.e. a wide data.table: "Block" + one column per "<Sample>.<haplotype_idx>".
# =============================================================================


#' Optimal Haplotype Selection of parents
#'
#' Selects \code{n_parents} candidates that maximise expected next-cross gain,
#' using one of three engines (\code{method}). \code{GA_pairwise} and
#' \code{Greedy} optimise the set objective built on the block-value matrix
#' \code{\link{.pairwise_block_value_matrix}} — by default
#' \code{sum(M[S, S])}, every family's stacking ceiling summed;
#' \code{GA_coverage} delegates to \code{\link{genetic_algorithm_selection}}
#' (binary block coverage) after collapsing haplotypes to genotype level.
#'
#' @section Elite subsetting:
#' With \code{elite_ratio < 1} the objective counts only the best
#' \code{q = ceiling(elite_ratio * P)} of the \code{P = k (k - 1) / 2} families
#' a set of \code{k} parents implies:
#' \describe{
#'   \item{\code{elite_level = "pair"} (default)}{the \code{q} largest family
#'     ceilings \code{M[i, j]} within the set. As the ratio shrinks the objective
#'     tends to the single best cross; a set is then judged by its elite
#'     families and the other parents matter only through a tie-break on the
#'     full sum.}
#'   \item{\code{elite_level = "block"}}{per block, the \code{q} largest
#'     per-block family ceilings \code{M_b[i, j]}, summed over blocks. Elite
#'     families may differ from block to block, so as the ratio shrinks the
#'     objective tends to the best two haplotypes of every block from any two
#'     parents — the "ultimate genotype" a gene pool could eventually build
#'     (Kemper et al. 2012; HapSelect's \code{Haploid_OHS}), rather than the
#'     best single cross.}
#' }
#' At \code{elite_ratio = 1} both levels equal the full objective. The ratio
#' changes what the search values, not the crossing design: the driver still
#' crosses the selected set in its full traversal. When \code{q < k - 1} some
#' parents cannot appear in any elite family; a tie-break of \code{1e-6} times
#' the mean family ceiling keeps them ordered, but they are fillers as far as
#' the elite objective is concerned.
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
#' @param elite_ratio Fraction in \code{(0, 1]} of the set's families the
#'   objective counts (see the section above). Default \code{1} = every family.
#'   Not available with \code{method = "GA_coverage"}.
#' @param elite_level \code{"pair"} (default) or \code{"block"}: whether the
#'   elite families are chosen once per set or once per block.
#' @param pop_size,n_generations,mutation_rate GA controls (used by
#'   \code{GA_pairwise} and \code{GA_coverage}).
#' @param seed Random seed for the GA engines. \code{GA_pairwise} restores the
#'   global RNG stream on exit (0.7.3). Default \code{42L}.
#' @param verbose Logical; print GA progress. Default \code{FALSE}.
#'
#' @return A character vector of selected sample names (length \code{n_parents}).
#'
#' @seealso \code{\link{ocs_selection}}, \code{\link{genetic_algorithm_selection}},
#'   \code{\link{ohs_pair_selection}}, \code{\link{make_cross_design}}
#' @export
ohs_selection <- function(localgebv,
                          n_parents     = 8L,
                          method        = c("GA_pairwise", "GA_coverage", "Greedy"),
                          ploidy        = NULL,
                          elite_ratio   = 1,
                          elite_level   = c("pair", "block"),
                          pop_size      = 100L,
                          n_generations = 200L,
                          mutation_rate = 0.1,
                          seed          = 42L,
                          verbose       = FALSE) {
  method      <- match.arg(method)
  elite_level <- match.arg(elite_level)
  .ohs_check_elite(elite_ratio)
  dt <- as.data.table(localgebv)
  if (!"Block" %in% names(dt))
    stop("ohs_selection: 'localgebv' must have a 'Block' column.")
  if (length(setdiff(names(dt), "Block")) == 0L)
    stop("ohs_selection: 'localgebv' has no sample columns.")

  if (method == "GA_coverage") {
    if (elite_ratio < 1)
      stop("ohs_selection: 'elite_ratio' applies to the pairwise objective; ",
           "method = \"GA_coverage\" scores block coverage and has no family ",
           "reading. Use \"GA_pairwise\" or \"Greedy\".")
    geno <- .collapse_hap_to_geno(dt)
    return(genetic_algorithm_selection(
      geno, n_parents = n_parents, pop_size = pop_size,
      n_generations = n_generations, mutation_rate = mutation_rate,
      seed = seed, verbose = verbose))
  }

  need_array <- elite_level == "block" && elite_ratio < 1
  pm <- .pairwise_block_value_matrix(dt, ploidy = ploidy, return_array = need_array)
  M  <- if (need_array) pm$M else pm
  A  <- if (need_array) pm$A else NULL
  if (n_parents >= nrow(M)) return(rownames(M))

  obj <- .ohs_objective(M, A, elite_ratio = elite_ratio, elite_level = elite_level)

  if (method == "GA_pairwise")
    .ga_pairwise(M, n_parents, obj = obj, pop_size = pop_size,
                 n_generations = n_generations, mutation_rate = mutation_rate,
                 seed = seed, verbose = verbose)
  else
    .greedy_pairwise(M, n_parents, obj = obj)
}

# Argument check shared by ohs_selection(), ohs_pair_selection() and the driver.
.ohs_check_elite <- function(elite_ratio) {
  if (!is.numeric(elite_ratio) || length(elite_ratio) != 1L ||
      is.na(elite_ratio) || elite_ratio <= 0 || elite_ratio > 1)
    stop("'elite_ratio' must be a single number in (0, 1].")
  invisible(TRUE)
}

# Number of elite families out of `n_pairs` for a ratio (at least one; the
# small offset guards ceiling() against 0.3 * 10 = 3.0000000000000004).
.ohs_elite_q <- function(n_pairs, elite_ratio) {
  if (n_pairs <= 0L) return(0L)
  max(1L, as.integer(ceiling(elite_ratio * n_pairs - 1e-9)))
}


# ---------------------------------------------------------------------------
# Internal — pairwise block-value matrix
# ---------------------------------------------------------------------------

# Build the symmetric n_samples x n_samples matrix
#   M[i,j] = sum_block ( sum of top-`ploidy` values among the block haplotype
#                        values of samples i and j ).
# `localgebv` is the wide haplotype-level table (Block + "<sample>.<idx>" cols).
# With `return_array = TRUE` returns list(M, A, ploidy) where
# A[sample, block, hap] holds each sample's block values sorted descending and
# zero-padded to `ploidy` — the input of the block-level elite objective.
.pairwise_block_value_matrix <- function(localgebv, ploidy = NULL,
                                         return_array = FALSE) {
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
  if (isTRUE(return_array)) list(M = M, A = A, ploidy = ploidy) else M
}

# Per-block family ceilings M_b[i, j] for the families (i_vec[f], j_vec[f]),
# from the sorted array A: a matrix families x blocks. Vectorised for ploidy 2
# (the same top-2 rule as the matrix build); general ploidy sorts per cell.
# A self (i == j) gives 2 x the parent's best haplotype per block.
.ohs_pair_block_values <- function(A, i_vec, j_vec) {
  nb <- dim(A)[2L]; ploidy <- dim(A)[3L]
  nf <- length(i_vec)
  if (nf == 0L) return(matrix(0, 0L, nb))
  if (ploidy == 2L) {
    a1 <- matrix(A[i_vec, , 1L], nf, nb); a2 <- matrix(A[i_vec, , 2L], nf, nb)
    b1 <- matrix(A[j_vec, , 1L], nf, nb); b2 <- matrix(A[j_vec, , 2L], nf, nb)
    ge <- a1 >= b1
    out <- ifelse(ge, a1 + pmax(a2, b1), b1 + pmax(b2, a1))
    return(matrix(out, nf, nb))
  }
  out <- matrix(0, nf, nb)
  for (f in seq_len(nf)) {
    comb <- cbind(matrix(A[i_vec[f], , ], nb, ploidy),
                  matrix(A[j_vec[f], , ], nb, ploidy))
    out[f, ] <- apply(comb, 1L, function(r)
      sum(sort(r, decreasing = TRUE)[seq_len(ploidy)]))
  }
  out
}

# Sum of the q largest entries of every column of V (families x blocks), in
# one order() call over the whole matrix; q >= nrow(V) is the plain column sum.
.ohs_col_topq_sum <- function(V, q) {
  P <- nrow(V)
  if (P == 0L) return(0)
  if (q >= P) return(sum(V))
  o  <- order(col(V), -V)
  Vs <- matrix(V[o], nrow = P)
  sum(Vs[seq_len(q), , drop = FALSE])
}


# ---------------------------------------------------------------------------
# Internal — the set objective for an elite ratio and level
# ---------------------------------------------------------------------------

# Returns list(fit = function(idx) -> numeric (maximise), closed_form = logical,
# q_of = function(k) -> q). `fit` on the sample indices `idx` of a set:
#   ratio 1                 sum(M[idx, idx])                     (closed form)
#   "pair", ratio < 1       sum of the q largest M[i, j], i < j in idx,
#                           + 1e-6 * mean family ceiling (tie-break)
#   "block", ratio < 1      sum_b of the q largest M_b[i, j] over the families,
#                           + 1e-6 * mean family ceiling (tie-break)
# `A` (the sorted array) is required for the block level only.
.ohs_objective <- function(M, A = NULL, elite_ratio = 1,
                           elite_level = c("pair", "block")) {
  elite_level <- match.arg(elite_level)
  eps <- 1e-6
  if (elite_ratio >= 1)
    return(list(fit = function(idx) sum(M[idx, idx]), closed_form = TRUE,
                q_of = function(k) .ohs_elite_q(k * (k - 1L) / 2L, 1)))
  q_of <- function(k) .ohs_elite_q(as.integer(k * (k - 1L) / 2L), elite_ratio)

  if (elite_level == "pair") {
    fit <- function(idx) {
      k <- length(idx)
      if (k < 2L) return(0)
      m <- M[idx, idx]
      v <- m[upper.tri(m)]
      q <- q_of(k)
      s <- sort(v, decreasing = TRUE)
      sum(s[seq_len(q)]) + eps * mean(v)
    }
    return(list(fit = fit, closed_form = FALSE, q_of = q_of))
  }

  if (is.null(A))
    stop(".ohs_objective: the block level needs the sorted haplotype array 'A'.")
  fit <- function(idx) {
    k <- length(idx)
    if (k < 2L) return(0)
    cb <- utils::combn(idx, 2L)
    V  <- .ohs_pair_block_values(A, cb[1L, ], cb[2L, ])
    .ohs_col_topq_sum(V, q_of(k)) + eps * sum(V) / nrow(V)
  }
  list(fit = fit, closed_form = FALSE, q_of = q_of)
}


# ---------------------------------------------------------------------------
# Internal — subset optimisers on the set objective
# ---------------------------------------------------------------------------

# Greedy marginal-gain selection. The first parent is the best row sum of M
# (the best average partner); every later pick maximises the objective of the
# enlarged set. Under the full objective the marginal gain of adding c to S is
# the closed form M[c,c] + 2 * sum_{s in S} M[s,c]; otherwise the objective is
# re-evaluated for every candidate.
.greedy_pairwise <- function(M, n_select, obj = NULL) {
  n  <- nrow(M)
  nm <- rownames(M)
  if (n_select >= n) return(nm)
  if (is.null(obj)) obj <- .ohs_objective(M)

  d         <- diag(M)
  selected  <- which.max(rowSums(M))
  remaining <- setdiff(seq_len(n), selected)

  while (length(selected) < n_select) {
    if (isTRUE(obj$closed_form)) {
      marg <- d[remaining] + 2 * colSums(M[selected, remaining, drop = FALSE])
    } else {
      marg <- vapply(remaining, function(c) obj$fit(c(selected, c)), numeric(1L))
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

# Genetic algorithm maximising the set objective (sum(M[S, S]) by default).
# The global RNG stream is restored on exit (0.7.3), so a driver's crossing
# step is not reseeded by a selection call.
.ga_pairwise <- function(M, n_select, obj = NULL, pop_size = 100L,
                         n_generations = 200L, mutation_rate = 0.1, seed = 42L,
                         verbose = FALSE) {
  n  <- nrow(M)
  nm <- rownames(M)
  if (n_select >= n) return(nm)
  if (is.null(obj)) obj <- .ohs_objective(M)
  fit <- obj$fit

  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(list = ".Random.seed", envir = globalenv())
    } else assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)
  set.seed(seed)

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
      message("  OHS-GA gen ", g, " — best objective = ", round(bestf, 3L))
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
