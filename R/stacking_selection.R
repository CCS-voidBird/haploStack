# =============================================================================
# stacking_selection.R
#
# Selection algorithm functions for the haplotype stacking simulation:
#
#   greedy_selection()            — iterative greedy parent selection
#   genetic_algorithm_selection() — GA-based parent selection
#   ocs_selection()               — Optimal Contribution Selection (ARCSIN)
#   make_cross_design()           — build a crossing design table
# =============================================================================

#' Greedy parent selection for haplotype stacking
#'
#' Selects \code{n_parents} founders from a local GEBV matrix (Block × Sample)
#' using an iterative greedy strategy: at each step, the candidate that
#' maximises the number of blocks with positive effect is added to the selected
#' set.
#'
#' @param gebv_matrix  A wide \code{data.table} or matrix (Block × Sample) of
#'   local GEBV values, as produced by \code{\link{compute_profile_localgebv}}.
#'   First column must be named \code{"Block"}.
#' @param n_parents    Integer. Number of parents to select.
#' @param ideal_freq   Optional named numeric vector (Block → target frequency
#'   of the positive haplotype). When supplied, the selection criterion
#'   maximises coverage of blocks with current frequency below the target.
#'   Default \code{NULL} — maximises count of non-zero blocks per candidate.
#' @param min_effect   Numeric threshold. Only blocks with
#'   \code{|effect| > min_effect} are considered "non-zero". Default \code{0}.
#'
#' @return A character vector of selected sample names (length \code{n_parents}).
#'
#' @details
#' Selection steps:
#' \enumerate{
#'   \item Convert the GEBV matrix to long format.
#'   \item Score each candidate by the number of blocks where it carries a
#'     positive effect (effect > \code{min_effect}) not yet covered by the
#'     currently selected set.
#'   \item Add the top-scoring candidate.
#'   \item Repeat until \code{n_parents} candidates are selected.
#' }
#'
#' @examples
#' \dontrun{
#' selected <- greedy_selection(gebv_matrix, n_parents = 8)
#' }
#'
#' @seealso \code{\link{genetic_algorithm_selection}},
#'   \code{\link{make_cross_design}}
#'
#' @export
greedy_selection <- function(gebv_matrix,
                              n_parents  = 8L,
                              ideal_freq = NULL,
                              min_effect = 0) {

  mat <- as.data.table(gebv_matrix)
  if (!"Block" %in% names(mat))
    stop("greedy_selection: 'gebv_matrix' must have a 'Block' column.")

  sample_cols <- setdiff(names(mat), "Block")
  n_cand      <- length(sample_cols)

  if (n_parents > n_cand) {
    warning("greedy_selection: n_parents (", n_parents, ") > available candidates (",
            n_cand, "). Returning all candidates.")
    return(sample_cols)
  }

  # Score matrix: (n_blocks, n_candidates) — 1 if positive, else 0
  score_mat <- as.matrix(mat[, sample_cols, with = FALSE])
  score_mat <- (score_mat > min_effect) * 1L

  selected   <- character(0L)
  covered    <- rep(0L, nrow(mat))  # blocks already covered by selected set

  for (k in seq_len(n_parents)) {
    remaining <- setdiff(sample_cols, selected)
    if (length(remaining) == 0L) break

    # Marginal gain: blocks NOT yet covered that each candidate covers
    marginal_gain <- sapply(remaining, function(s) {
      col_scores <- score_mat[, s]
      sum(col_scores & !covered)
    })

    best <- remaining[which.max(marginal_gain)]
    selected <- c(selected, best)
    covered  <- covered | (score_mat[, best] > 0L)
  }

  selected
}


#' Genetic Algorithm parent selection for haplotype stacking
#'
#' Applies a simple genetic algorithm to select \code{n_parents} candidates
#' that collectively maximise a composite objective: (a) total block coverage
#' (blocks with at least one positive carrier among selected parents) and
#' optionally (b) alignment of carrier frequencies with ideal target
#' frequencies.
#'
#' @param gebv_matrix  Wide Block × Sample \code{data.table} of local GEBVs.
#'   First column must be \code{"Block"}.
#' @param n_parents    Integer. Number of parents to select.
#' @param ideal_target Optional named numeric vector (Block → target carrier
#'   proportion, 0–1). Blocks with \code{ideal_target > 0} are treated as
#'   "target blocks". Default \code{NULL}.
#' @param pop_size     GA population size (number of candidate parent subsets).
#'   Default \code{100}.
#' @param n_generations Integer. Number of GA generations. Default \code{200}.
#' @param mutation_rate Probability of flipping one parent in a solution at
#'   each mutation step. Default \code{0.1}.
#' @param crossover_rate Probability of two-point crossover between parents.
#'   Default \code{0.8}.
#' @param seed         Random seed for reproducibility. Default \code{42L}.
#' @param verbose      Logical. Print progress every 50 generations.
#'   Default \code{FALSE}.
#'
#' @return A character vector of selected sample names.
#'
#' @examples
#' \dontrun{
#' selected <- genetic_algorithm_selection(gebv_matrix, n_parents = 8,
#'                                         n_generations = 100)
#' }
#'
#' @seealso \code{\link{greedy_selection}}, \code{\link{make_cross_design}}
#'
#' @export
genetic_algorithm_selection <- function(gebv_matrix,
                                         n_parents      = 8L,
                                         ideal_target   = NULL,
                                         pop_size       = 100L,
                                         n_generations  = 200L,
                                         mutation_rate  = 0.1,
                                         crossover_rate = 0.8,
                                         seed           = 42L,
                                         verbose        = FALSE) {

  set.seed(seed)
  mat <- as.data.table(gebv_matrix)
  if (!"Block" %in% names(mat))
    stop("genetic_algorithm_selection: 'gebv_matrix' must have a 'Block' column.")

  sample_cols <- setdiff(names(mat), "Block")
  n_cand      <- length(sample_cols)
  n_blocks    <- nrow(mat)

  if (n_parents >= n_cand) {
    warning("genetic_algorithm_selection: n_parents >= n_candidates; returning greedy solution.")
    return(greedy_selection(gebv_matrix, n_parents))
  }

  # Positive-effect indicator matrix (n_blocks × n_candidates)
  effect_mat <- as.matrix(mat[, sample_cols, with = FALSE])
  pos_mat    <- (effect_mat > 0) * 1L

  # Fitness function: proportion of blocks covered (with target weighting)
  .fitness <- function(parent_idx) {
    covered <- rowSums(pos_mat[, parent_idx, drop = FALSE]) > 0L
    base_score <- mean(covered)
    if (!is.null(ideal_target)) {
      target_blocks <- names(ideal_target)[ideal_target > 0]
      t_idx  <- which(mat$Block %in% target_blocks)
      if (length(t_idx) > 0L) {
        freq_vec  <- rowMeans(pos_mat[t_idx, parent_idx, drop = FALSE])
        target_vec <- ideal_target[mat$Block[t_idx]]
        target_vec[is.na(target_vec)] <- 0
        alignment <- 1 - mean(abs(freq_vec - target_vec))
        return(0.6 * base_score + 0.4 * alignment)
      }
    }
    base_score
  }

  # Initialise population
  population <- replicate(pop_size, sample(n_cand, n_parents), simplify = FALSE)
  fitness_scores <- vapply(population, .fitness, numeric(1L))

  # GA loop
  for (gen in seq_len(n_generations)) {
    # Selection (tournament, size 2)
    new_pop <- vector("list", pop_size)
    for (i in seq_len(pop_size)) {
      t1 <- sample(pop_size, 2L)
      t2 <- sample(pop_size, 2L)
      p1 <- population[[t1[which.max(fitness_scores[t1])]]]
      p2 <- population[[t2[which.max(fitness_scores[t2])]]]

      # Crossover
      if (runif(1L) < crossover_rate) {
        pt <- sample(n_parents - 1L, 1L)
        child <- unique(c(p1[seq_len(pt)], p2[(pt + 1L):n_parents]))
        if (length(child) < n_parents)
          child <- c(child, sample(setdiff(seq_len(n_cand), child),
                                   n_parents - length(child)))
        child <- child[seq_len(n_parents)]
      } else {
        child <- p1
      }

      # Mutation
      if (runif(1L) < mutation_rate) {
        swap_out <- sample(n_parents, 1L)
        swap_in  <- sample(setdiff(seq_len(n_cand), child), 1L)
        child[swap_out] <- swap_in
      }

      new_pop[[i]] <- child
    }

    population    <- new_pop
    fitness_scores <- vapply(population, .fitness, numeric(1L))

    if (verbose && gen %% 50L == 0L)
      message("  GA gen ", gen, " — best fitness: ",
              round(max(fitness_scores), 4L))
  }

  best_idx <- population[[which.max(fitness_scores)]]
  sample_cols[best_idx]
}


#' Optimal Contribution Selection (OCS / ARCSIN method)
#'
#' Selects parents to maximise expected genetic gain while controlling
#' inbreeding, using a mean-centred GEBV and kinship constraint. Implements
#' the ARCSIN approximation described in the pipeline documentation.
#'
#' @param gebv_vector       Named numeric vector of total GEBVs per sample.
#' @param kinship_matrix    Square numeric matrix of pairwise kinship /
#'   relatedness coefficients. Row and column names must match
#'   \code{names(gebv_vector)}. When \code{NULL} (default), falls back to
#'   pure GEBV ranking (truncation selection).
#' @param n_parents         Integer. Number of parents to select.
#' @param max_kinship       Maximum allowable mean pairwise kinship among
#'   selected parents. Default \code{0.125} (equivalent to first cousins).
#' @param ideal_genotype_df Optional data.frame representing the target
#'   pan-genomic population. Required columns: \code{Chr}, \code{Block},
#'   \code{Effect}, \code{Hap} (or \code{Ploidy}), \code{TargetFreq},
#'   \code{CurrentFreq}. Optional columns: \code{Start}, \code{End}.
#'   When supplied, a frequency-deviation penalty is added to the selection
#'   objective: \eqn{score = GEBV - \lambda \cdot \Sigma (CurrentFreq -
#'   TargetFreq)^2} per block. Default \code{NULL}.
#' @param freq_penalty_weight Numeric weight \eqn{\lambda} for the
#'   frequency-deviation penalty when \code{ideal_genotype_df} is supplied.
#'   Default \code{1}.
#' @param seed              Random seed. Default \code{42L}.
#'
#' @return A character vector of selected sample names.
#'
#' @details
#' Ideal genotype file format (Q7 decision, 2026-06-08):
#' \itemize{
#'   \item Required: \code{Chr}, \code{Block}, \code{Effect},
#'     \code{Hap} (or \code{Ploidy}), \code{TargetFreq}, \code{CurrentFreq}.
#'   \item Optional: \code{Start}, \code{End} (block coordinates; ignored in
#'     selection).
#'   \item Represents a pan-genomic-like reference with target haplotype
#'     frequencies per block.
#' }
#'
#' @export
ocs_selection <- function(gebv_vector,
                           kinship_matrix      = NULL,
                           n_parents           = 8L,
                           max_kinship         = 0.125,
                           ideal_genotype_df   = NULL,
                           freq_penalty_weight = 1,
                           seed                = 42L) {

  set.seed(seed)
  if (is.null(names(gebv_vector)))
    stop("ocs_selection: 'gebv_vector' must be named.")

  samples <- names(gebv_vector)

  # Q7: validate and process ideal_genotype_df
  freq_penalty <- 0
  if (!is.null(ideal_genotype_df)) {
    freq_penalty <- .compute_freq_penalty(ideal_genotype_df, freq_penalty_weight)
  }

  # Adjust GEBVs by frequency penalty (scalar penalty is block-level, not
  # sample-level; we apply it as a single genome-wide deduction to guide
  # selection toward blocks that are furthest from their target frequency)
  adjusted_gebv <- gebv_vector - freq_penalty

  if (is.null(kinship_matrix)) {
    # Pure truncation selection on adjusted GEBVs
    selected <- names(sort(adjusted_gebv, decreasing = TRUE))[seq_len(n_parents)]
    return(selected)
  }

  # Validate kinship matrix
  if (!all(samples %in% rownames(kinship_matrix)))
    stop("ocs_selection: all samples in 'gebv_vector' must be in 'kinship_matrix'.")

  K <- kinship_matrix[samples, samples]

  # Greedy kinship-constrained selection
  selected  <- character(0L)
  remaining <- samples

  while (length(selected) < n_parents && length(remaining) > 0L) {
    if (length(selected) == 0L) {
      score <- adjusted_gebv[remaining]
    } else {
      mean_k <- rowMeans(K[remaining, selected, drop = FALSE])
      # ARCSIN transformation for better numerical behaviour
      score  <- adjusted_gebv[remaining] - asin(sqrt(pmin(pmax(mean_k, 0), 1)))
    }

    best      <- remaining[which.max(score)]
    selected  <- c(selected, best)
    remaining <- setdiff(remaining, best)

    if (length(selected) > 1L) {
      sel_K       <- K[selected, selected]
      mean_pair_k <- mean(sel_K[upper.tri(sel_K)])
      if (mean_pair_k > max_kinship) {
        message("ocs_selection: kinship constraint violated; stopping at ",
                length(selected) - 1L, " parents.")
        selected <- selected[-length(selected)]
        break
      }
    }
  }

  selected
}


# ---------------------------------------------------------------------------
# Internal: ideal genotype frequency penalty
# ---------------------------------------------------------------------------

#' Compute the frequency-deviation penalty from an ideal genotype table
#'
#' Q7 (2026-06-08): validates the ideal genotype data.frame and returns the
#' scalar penalty \eqn{\lambda \cdot \Sigma (CurrentFreq - TargetFreq)^2}
#' summed across all blocks.
#'
#' @param ideal_df          data.frame with required columns.
#' @param freq_penalty_weight Numeric weight \eqn{\lambda}.
#' @return Scalar numeric penalty.
#' @keywords internal
.compute_freq_penalty <- function(ideal_df, freq_penalty_weight = 1) {
  df <- as.data.frame(ideal_df, stringsAsFactors = FALSE)

  # Accept either "Hap" or "Ploidy" as the haplotype column name
  if ("Ploidy" %in% names(df) && !"Hap" %in% names(df))
    names(df)[names(df) == "Ploidy"] <- "Hap"

  required <- c("Chrom", "Block", "Effect", "Hap", "TargetFreq", "CurrentFreq")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0L)
    stop("ocs_selection: ideal_genotype_df is missing required column(s): ",
         paste(missing, collapse = ", "))

  df$TargetFreq  <- as.numeric(df$TargetFreq)
  df$CurrentFreq <- as.numeric(df$CurrentFreq)

  penalty <- freq_penalty_weight *
    sum((df$CurrentFreq - df$TargetFreq)^2, na.rm = TRUE)

  penalty
}


#' Build a crossing design table
#'
#' Generates a table of proposed crosses between selected parents, using one of
#' three design strategies: all pairwise combinations (\code{"all_pairs"}),
#' sequential pairs (\code{"sequential"}), or a custom user-supplied table
#' (\code{"custom"}).
#'
#' @param selected_parents Character vector of parent sample names.
#' @param method           Crossing design method:
#'   \describe{
#'     \item{\code{"all_pairs"}}{All non-redundant pairwise combinations
#'       (n*(n-1)/2 crosses).}
#'     \item{\code{"sequential"}}{Each adjacent pair: (1,2), (3,4), … Useful
#'       for half-diallel schemes.}
#'     \item{\code{"custom"}}{User-supplied \code{cross_table} argument.}
#'   }
#' @param cross_table      data.frame with columns \code{Parent1} and
#'   \code{Parent2}. Required when \code{method = "custom"}.
#' @param avoid_selfing    Logical. Remove self-crosses from the design.
#'   Default \code{TRUE}.
#' @param n_offspring      Integer. Number of offspring per cross. Stored as a
#'   column in the output for use by the simulation loop.
#'
#' @return A \code{data.table} with columns \code{Parent1}, \code{Parent2},
#'   \code{Cross_ID} (character), and \code{n_offspring}.
#'
#' @examples
#' parents <- c("P1","P2","P3","P4")
#' design  <- make_cross_design(parents, method = "all_pairs", n_offspring = 4)
#' print(design)
#'
#' @export
make_cross_design <- function(selected_parents,
                               method        = c("all_pairs","sequential","custom"),
                               cross_table   = NULL,
                               avoid_selfing = TRUE,
                               n_offspring   = 4L) {

  method <- match.arg(method)
  n      <- length(selected_parents)

  if (method == "all_pairs") {
    pairs <- as.data.table(t(combn(selected_parents, 2L)))
    setnames(pairs, c("Parent1", "Parent2"))

  } else if (method == "sequential") {
    if (n < 2L) stop("make_cross_design: need >= 2 parents for sequential design.")
    idx   <- seq(1L, n - 1L, by = 2L)
    pairs <- data.table(
      Parent1 = selected_parents[idx],
      Parent2 = selected_parents[pmin(idx + 1L, n)]
    )

  } else {  # custom
    if (is.null(cross_table))
      stop("make_cross_design: 'cross_table' required for method = 'custom'.")
    pairs <- as.data.table(cross_table)
    if (!all(c("Parent1","Parent2") %in% names(pairs)))
      setnames(pairs, names(pairs)[1:2], c("Parent1","Parent2"))
  }

  if (avoid_selfing) {
    pairs <- pairs[Parent1 != Parent2]
  }

  pairs[, Cross_ID    := paste0(Parent1, "x", Parent2)]
  pairs[, n_offspring := as.integer(n_offspring)]

  pairs
}


# ---------------------------------------------------------------------------
# Family count of a parent set, and the direct-pairing budget (placeholder seam)
# ---------------------------------------------------------------------------

# Families a full traversal of `n_parents` parents makes: n(n-1)/2 distinct
# pairs, plus the n self-crosses when the design allows them. This is the
# reference quantity for direct pair selection (pair_selection_budget()), and
# the driver sizes its compact family IDs from it.
.stk_n_families <- function(n_parents, allow_selfing = FALSE) {
  n <- as.integer(n_parents)
  as.integer(choose(n, 2L) + if (isTRUE(allow_selfing)) n else 0L)
}


#' Budget of a direct-pair selection run
#'
#' The stacking driver selects a \emph{parent subgroup} first and then crosses it
#' in a full traversal (\code{\link{make_cross_design}}, \code{"all_pairs"}).
#' Selecting \emph{pairs} directly (\code{selection_unit = "pairs"} in
#' \code{\link{run_stacking_simulation}} — a placeholder in this version) has to
#' be matched to that protocol by one of two budgets, which conflict in general
#' and must therefore be chosen explicitly:
#' \describe{
#'   \item{\code{"families"}}{The same offspring population size: as many
#'     families as the traversal of \code{top_founders} parents would make,
#'     \code{choose(top_founders, 2)} (\code{+ top_founders} when self-crosses
#'     are allowed). The number of distinct parents is then free — anywhere
#'     between roughly \code{sqrt(2 * n_families)} and \code{2 * n_families}.}
#'   \item{\code{"parents"}}{The same number of parents: pairs are added until
#'     \code{top_founders} distinct parent names are in use. The number of
#'     families is then free — between \code{ceiling(top_founders / 2)} (every
#'     pair brings two new parents) and the full traversal.}
#' }
#' This function is the "quick calculation" of what a given \code{top_founders}
#' implies under each budget.
#'
#' @param top_founders Integer. Parents per cycle of the parents-first protocol
#'   (the driver's \code{top_founders}).
#' @param pair_budget \code{"families"} or \code{"parents"} (see Details).
#' @param allow_selfing Logical. Count the \code{top_founders} self-crosses in
#'   the traversal. Default \code{FALSE} — the driver's \code{"all_pairs"} design
#'   never contains a self-cross (it enumerates \code{combn(parents, 2)}), so the
#'   driver passes \code{FALSE} whatever its \code{avoid_selfing} argument says.
#'
#' @return A list: \code{pair_budget}, \code{n_parents_ref} (the
#'   \code{top_founders} the budget is derived from), \code{allow_selfing},
#'   \code{n_families} (fixed under \code{"families"}, \code{NA} under
#'   \code{"parents"}), \code{n_families_min} / \code{n_families_max} (the range
#'   a run can produce), and \code{n_parents_max} (fixed under \code{"parents"},
#'   \code{NA} under \code{"families"}).
#'
#' @examples
#' pair_selection_budget(16, "families")                       # 120 families
#' pair_selection_budget(16, "families", allow_selfing = TRUE) # 136 families
#' pair_selection_budget(16, "parents")                        # 8..120 families
#'
#' @seealso \code{\link{run_stacking_simulation}}, \code{\link{make_cross_design}}
#' @export
pair_selection_budget <- function(top_founders,
                                  pair_budget = c("families", "parents"),
                                  allow_selfing = FALSE) {
  pair_budget <- match.arg(pair_budget)
  n <- as.integer(top_founders)
  if (length(n) != 1L || is.na(n) || n < 2L)
    stop("pair_selection_budget: 'top_founders' must be a single integer >= 2.")
  full <- .stk_n_families(n, allow_selfing)
  if (pair_budget == "families")
    list(pair_budget = "families", n_parents_ref = n,
         allow_selfing = isTRUE(allow_selfing),
         n_families = full, n_families_min = full, n_families_max = full,
         n_parents_max = NA_integer_)
  else
    list(pair_budget = "parents", n_parents_ref = n,
         allow_selfing = isTRUE(allow_selfing),
         n_families = NA_integer_,
         n_families_min = as.integer(ceiling(n / 2)), n_families_max = full,
         n_parents_max = n)
}


#' Compute a Genomic Relationship Matrix (GRM) from a GenoMeta object
#'
#' Builds the additive genomic relationship matrix via the first method of
#' VanRaden (2008). Dosage is extracted internally from the \code{GenoMeta}
#' using \code{.meta_to_dosage_matrix()} with the minor allele as the counted
#' allele for each marker. Missing values are mean-imputed before computation.
#'
#' To restrict the GRM to a marker subset, call
#' \code{\link{filter_geno_meta}} first.
#'
#' \deqn{G = \frac{Z Z'}{2 \sum_j p_j (1 - p_j)}}
#' where \eqn{Z = M - 2\mathbf{p}'} is the centred marker dosage matrix and
#' \eqn{p_j} is the minor allele frequency of marker \eqn{j}.
#'
#' @param geno_meta A valid \code{GenoMeta} list. Use
#'   \code{\link{filter_geno_meta}} to subset markers before calling if a
#'   marker-filtered GRM is needed.
#' @param min_maf   Numeric. Markers with minor allele frequency below
#'   \code{min_maf} are excluded before computing the GRM. Default
#'   \code{0.01}.
#'
#' @return A symmetric named numeric matrix (n_samples × n_samples) of
#'   genomic relationship coefficients. Diagonal values near 1 indicate
#'   expected heterozygosity; off-diagonal values reflect genetic similarity.
#'
#' @references VanRaden, P.M. (2008). Efficient methods to compute genomic
#'   predictions. \emph{Journal of Dairy Science}, 91(11), 4414–4423.
#'   \doi{10.3168/jds.2007-0980}
#'
#' @examples
#' \dontrun{
#' # Build GRM directly from a GenoMeta
#' meta <- read_st_genotype("geno.csv")
#' grm  <- compute_grm(meta)
#'
#' # Chromosome-1-only GRM
#' meta_chr1 <- filter_geno_meta(meta,
#'   keep_markers = meta$map_info$Marker[meta$map_info$Chrom == "1"])
#' grm_chr1 <- compute_grm(meta_chr1)
#'
#' selected <- ocs_selection(gebv_vec, kinship_matrix = grm, n_parents = 8)
#' }
#'
#' @seealso \code{\link{filter_geno_meta}}, \code{\link{ocs_selection}}
#'
#' @export
compute_grm <- function(geno_meta, min_maf = 0.01) {

  if (!is_geno_meta(geno_meta))
    stop("compute_grm: 'geno_meta' must be a valid GenoMeta. ",
         "Use filter_geno_meta() to subset markers before calling if needed.")

  # Extract minor-allele dosage matrix [n_samples × n_markers]
  M <- .meta_to_dosage_matrix(geno_meta)   # target_alleles = NULL → minor allele

  if (!is.numeric(M))
    stop("compute_grm: dosage matrix is not numeric; check geno_meta encoding.")
  if (nrow(M) < 2L)
    stop("compute_grm: need at least 2 samples.")
  if (ncol(M) < 1L)
    stop("compute_grm: need at least 1 marker.")

  m <- ncol(M)

  # Minor allele frequency: dosage / ploidy
  ploidy <- if (!is.null(geno_meta$metadata$ploidy)) geno_meta$metadata$ploidy else 2L
  p <- colMeans(M, na.rm = TRUE) / ploidy
  p <- pmin(pmax(p, 0), 1)

  # MAF filter (p is already the minor-allele freq)
  maf  <- pmin(p, 1 - p)
  keep <- maf >= min_maf
  if (sum(keep) == 0L)
    stop("compute_grm: no markers remain after MAF filter (min_maf = ",
         min_maf, ").")
  if (sum(keep) < m)
    message("compute_grm: ", m - sum(keep),
            " markers removed by MAF filter (", sum(keep), " retained).")

  M <- M[, keep, drop = FALSE]
  p <- p[keep]

  # Mean-impute missing values
  na_idx <- which(is.na(M), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) {
    M[na_idx] <- ploidy * p[na_idx[, 2L]]
    message("compute_grm: imputed ", nrow(na_idx),
            " missing genotype(s) to ploidy * p.")
  }

  # Centre: Z = M - ploidy * p  (broadcast across rows)
  Z <- sweep(M, 2L, ploidy * p, FUN = "-")

  # Scale factor: ploidy * sum(p_j * (1 - p_j))
  scale_factor <- ploidy * sum(p * (1 - p))
  if (scale_factor <= 0)
    stop("compute_grm: scale factor is zero; all markers may be monomorphic.")

  G <- (Z %*% t(Z)) / scale_factor

  sample_names <- dimnames(geno_meta$geno_arrays)[[1L]]
  rownames(G) <- colnames(G) <- sample_names

  G
}
