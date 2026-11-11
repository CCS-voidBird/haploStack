# =============================================================================
# stacking_monitor.R
#
# Monitoring and visualisation functions for stacking simulation results.
#
# Functions:
#   compute_parent_coancestry()    — block-IBS coancestry among a sample set
#   plot_simulation_trajectories() — GEBV and block-coverage trajectories
#   plot_simulation_metrics()      — F1mean / D2Ideal / coancestry trajectories
#   summarise_simulation()         — per-cycle summary statistics
# =============================================================================


# ---------------------------------------------------------------------------
# Block-IBS coancestry (FR5 metric) — parent-set relatedness
# ---------------------------------------------------------------------------

#' Block-IBS coancestry among a set of samples
#'
#' Computes a symmetric coancestry matrix over a sample set from block-level
#' identity-by-state (IBS). For a pair of individuals the coancestry is the mean,
#' over blocks, of the average IBS across the four haplotype pairings
#' (\code{i.1-j.1}, \code{i.1-j.2}, \code{i.2-j.1}, \code{i.2-j.2}) — i.e. the
#' probability that a random allele drawn from one individual is identical by
#' state to a random allele from the other, per block, averaged over the genome.
#' Each block's IBS reuses \code{\link{compute_segment_ibs}} on the block's
#' markers, so missing alleles are excluded pairwise.
#'
#' This is the coancestry metric reported per cycle by
#' \code{\link{run_stacking_simulation}} (FR5): the mean off-diagonal value over
#' the selected parents. Unlike \code{\link{compute_block_ibd}} it needs no
#' pedigree — it measures block-IBS directly among the genotypes.
#'
#' @param geno_meta A \code{GenoMeta} object.
#' @param samples Optional character vector of sample names to restrict to.
#'   Default \code{NULL} (all samples).
#' @param blocks Optional character vector of block ids to restrict the coancestry
#'   to (e.g. the blocks selection focuses on). Default \code{NULL} (all blocks).
#'
#' @return A symmetric numeric matrix (rows/cols named by sample) of block-IBS
#'   coancestry, diagonal included (each individual's within-genotype IBS).
#'   Returns a \code{1 x 1} matrix when only one sample is available.
#'
#' @seealso \code{\link{compute_segment_ibs}}, \code{\link{compute_block_ibd}},
#'   \code{\link{run_stacking_simulation}}
#' @export
compute_parent_coancestry <- function(geno_meta, samples = NULL, blocks = NULL) {
  .hs_check_meta(geno_meta, "compute_parent_coancestry", "geno_meta")

  arr         <- geno_meta$geno_arrays
  all_samples <- dimnames(arr)[[1L]]
  markers     <- dimnames(arr)[[2L]]

  if (is.null(samples)) samples <- all_samples
  samples <- intersect(as.character(samples), all_samples)
  P <- length(samples)
  if (P == 0L)
    .hs_abort("compute_parent_coancestry: none of 'samples' are in the GenoMeta.")

  # Marker indices per block (unblocked markers form singleton blocks so they
  # still contribute; restrict to 'blocks' when given).
  bmm <- geno_meta$block_marker_map_info
  if (is.null(bmm)) {
    mk_by_block <- list(all = seq_along(markers))
  } else {
    bmm <- as.data.table(bmm)
    if (!is.null(blocks)) bmm <- bmm[Block %in% blocks]
    idx <- match(as.character(bmm$Marker), markers)
    keep <- !is.na(idx)
    mk_by_block <- split(idx[keep], as.character(bmm$Block)[keep])
    mk_by_block <- mk_by_block[lengths(mk_by_block) > 0L]
  }
  if (length(mk_by_block) == 0L)
    .hs_abort("compute_parent_coancestry: no markers remain for the requested blocks.")

  sub <- arr[samples, , , drop = FALSE]            # [P x M x 2]

  # Compiled kernel (0.9.0): same numbers bit for bit, minus the R double loop
  # that cost ~47 s per cycle under metrics = "all" (thread 013). Alleles become
  # integer codes here (equal strings -> equal codes; NA or "" -> NA), which is
  # exactly the comparability rule of compute_segment_ibs().
  if (dim(sub)[3L] >= 2L && .hs_use_rcpp("block_ibs_coancestry_cpp")) {
    v <- as.vector(sub)
    v[!is.na(v) & !nzchar(v)] <- NA_character_
    codes <- match(v, unique(v[!is.na(v)]))
    K <- .hs_kernel("block_ibs_coancestry_cpp")(
      as.integer(codes), as.integer(dim(sub)),
      lapply(unname(mk_by_block), function(ix) as.integer(ix) - 1L))
    dimnames(K) <- list(samples, samples)
    return(K)
  }

  K <- matrix(NA_real_, P, P, dimnames = list(samples, samples))

  pair_ibs <- function(i, j) {
    bvals <- vapply(mk_by_block, function(ix) {
      mean(c(
        compute_segment_ibs(sub[i, ix, 1L], sub[j, ix, 1L]),
        compute_segment_ibs(sub[i, ix, 1L], sub[j, ix, 2L]),
        compute_segment_ibs(sub[i, ix, 2L], sub[j, ix, 1L]),
        compute_segment_ibs(sub[i, ix, 2L], sub[j, ix, 2L])), na.rm = TRUE)
    }, numeric(1L))
    mean(bvals, na.rm = TRUE)
  }

  for (i in seq_len(P)) for (j in i:P) {
    v <- pair_ibs(i, j)
    K[i, j] <- v
    K[j, i] <- v
  }
  K
}


# Mean off-diagonal (non-self) block-IBS coancestry among `parents` — the FR5
# `init_coancestry` scalar reported per cycle by run_stacking_simulation().
# Returns NA when fewer than two parents are available.
.stk_block_ibs_coancestry <- function(geno_meta, parents, blocks = NULL) {
  parents <- intersect(as.character(parents), dimnames(geno_meta$geno_arrays)[[1L]])
  if (length(parents) < 2L) return(NA_real_)
  K <- tryCatch(
    compute_parent_coancestry(geno_meta, samples = parents, blocks = blocks),
    error = function(e) NULL)
  if (is.null(K) || nrow(K) < 2L) return(NA_real_)
  mean(K[upper.tri(K)], na.rm = TRUE)
}


# ---------------------------------------------------------------------------
# Nei's (1973) gene diversity (FR5 metric) — population genetic diversity
# ---------------------------------------------------------------------------

# Nei's gene diversity for a single locus from a vector of allele copies:
# h = 1 - Σ p_i^2 over the observed allele frequencies. Missing / empty copies
# are dropped; returns NA when no copy is comparable or the locus is monomorphic
# with a single missing-free allele (h = 0).
.nei_locus_h <- function(alleles) {
  a <- alleles[!is.na(alleles) & nzchar(alleles)]
  n <- length(a)
  if (n == 0L) return(NA_real_)
  p <- tabulate(match(a, unique(a))) / n
  1 - sum(p * p)
}

#' Nei's (1973) gene diversity of a population
#'
#' Computes Nei's gene diversity — the expected heterozygosity
#' \eqn{H = 1 - \sum_i p_i^2} averaged over loci, where \eqn{p_i} are the allele
#' frequencies among all gene copies in the population (Nei 1973,
#' \emph{PNAS} 70:3321). Two resolutions are offered:
#' \describe{
#'   \item{\code{level = "marker"}}{Canonical per-marker gene diversity: each
#'     marker is a locus and its single-character alleles are the variants,
#'     pooled over both haplotype copies of every sample. \eqn{H} is the mean
#'     over markers.}
#'   \item{\code{level = "haplotype"}}{Block-haplotype gene diversity: each block
#'     is treated as a multi-allelic locus whose "alleles" are the distinct
#'     haplotypes carried across the block's markers (a gene copy contributing a
#'     haplotype with any missing marker is dropped). \eqn{H} is the mean over
#'     blocks — the diversity resolution most relevant to haplotype stacking.}
#' }
#'
#' This is the \code{gene_diversity} metric reported per cycle by
#' \code{\link{run_stacking_simulation}} (FR5): as the population converges under
#' selection it falls toward 0 (loss of diversity), complementing the
#' \code{\link{compute_parent_coancestry}} coancestry (which rises toward 1).
#'
#' @param geno_meta A \code{GenoMeta} object.
#' @param level \code{"marker"} (default, canonical Nei 1973) or
#'   \code{"haplotype"} (per block-haplotype).
#' @param samples Optional character vector of samples to restrict the population
#'   to. Default \code{NULL} (all samples).
#' @param blocks Optional character vector of block ids (used when
#'   \code{level = "haplotype"}, or to restrict markers to those blocks when
#'   \code{level = "marker"}). Default \code{NULL} (all).
#'
#' @return A single numeric: mean gene diversity over loci, with the per-locus
#'   diversities attached as \code{attr(., "per_locus")} (named by marker or
#'   block).
#'
#' @references Nei, M. (1973) Analysis of gene diversity in subdivided
#'   populations. \emph{PNAS} 70(12):3321-3323.
#'
#' @seealso \code{\link{compute_parent_coancestry}},
#'   \code{\link{run_stacking_simulation}}
#' @export
compute_gene_diversity <- function(geno_meta, level = c("marker", "haplotype"),
                                   samples = NULL, blocks = NULL) {
  level <- match.arg(level)
  .hs_check_meta(geno_meta, "compute_gene_diversity", "geno_meta")

  arr     <- geno_meta$geno_arrays
  all_s   <- dimnames(arr)[[1L]]
  if (!is.null(samples)) {
    samples <- intersect(as.character(samples), all_s)
    if (length(samples) == 0L)
      .hs_abort("compute_gene_diversity: none of 'samples' are in the GenoMeta.")
    arr <- arr[samples, , , drop = FALSE]
  }
  markers <- dimnames(arr)[[2L]]

  # Block -> marker indices (shared by both levels; scopes markers for "marker"
  # and defines the loci for "haplotype").
  bmm <- geno_meta$block_marker_map_info
  mk_by_block <- NULL
  if (!is.null(bmm)) {
    bmm <- as.data.table(bmm)
    if (!is.null(blocks)) bmm <- bmm[Block %in% blocks]
    ix   <- match(as.character(bmm$Marker), markers)
    keep <- !is.na(ix)
    mk_by_block <- split(ix[keep], as.character(bmm$Block)[keep])
    mk_by_block <- mk_by_block[lengths(mk_by_block) > 0L]
  }

  if (level == "marker") {
    keep_mk <- if (!is.null(mk_by_block)) .hs_sort(unique(unlist(mk_by_block)))
               else seq_along(markers)
    h <- vapply(keep_mk, function(m)
      .nei_locus_h(c(arr[, m, 1L], arr[, m, 2L])), numeric(1L))
    names(h) <- markers[keep_mk]
  } else {
    if (is.null(mk_by_block))
      .hs_abort("compute_gene_diversity: level = 'haplotype' needs block assignments ",
           "(attach_block_info() / block_marker_map_info).")
    hap_string <- function(layer_mat) {
      # one string per sample; NA if any marker in the block is missing
      apply(layer_mat, 1L, function(z)
        if (anyNA(z) || any(!nzchar(z))) NA_character_
        else paste(z, collapse = "_"))
    }
    h <- vapply(mk_by_block, function(ixb) {
      m1 <- matrix(arr[, ixb, 1L], nrow = dim(arr)[1L])
      m2 <- matrix(arr[, ixb, 2L], nrow = dim(arr)[1L])
      .nei_locus_h(c(hap_string(m1), hap_string(m2)))
    }, numeric(1L))
    names(h) <- names(mk_by_block)
  }

  out <- mean(h, na.rm = TRUE)
  attr(out, "per_locus") <- h
  out
}


#' Plot simulation GEBV trajectories across cycles
#'
#' Creates a line plot of mean GEBV (and optionally ± 1 SD) across breeding
#' cycles, faceted or coloured by selection method and replication.
#'
#' @param trajectories data.table from \code{\link{run_stacking_simulation}}
#'   with columns \code{Rep}, \code{Cycle}, \code{Method}, \code{MeanGEBV},
#'   \code{SDev}, \code{BlockCoverage}.
#' @param trait_name   Character label for plot title.
#' @param show_ribbon  Logical. Draw ± 1 SD ribbon. Default \code{TRUE}.
#' @param show_coverage Logical. Add a second panel for block coverage.
#'   Default \code{FALSE}.
#' @param colour_by    \code{"Method"} (default) or \code{"Rep"}.
#'
#' @return A \code{ggplot} object (or a combined \code{cowplot} figure when
#'   \code{show_coverage = TRUE}).
#'
#' @examples
#' \dontrun{
#' p <- plot_simulation_trajectories(sim_results$trajectories,
#'                                   trait_name = "Yield")
#' ggplot2::ggsave("trajectories.png", p, width = 10, height = 5)
#' }
#'
#' @export
plot_simulation_trajectories <- function(trajectories,
                                          trait_name    = "Trait",
                                          show_ribbon   = TRUE,
                                          show_coverage = FALSE,
                                          colour_by     = "Method") {

  traj <- as.data.table(trajectories)
  req_cols <- c("Rep","Cycle","MeanGEBV")
  if (!all(req_cols %in% names(traj)))
    .hs_abort("plot_simulation_trajectories: trajectories must have columns ",
         paste(req_cols, collapse = ", "))

  # Aggregate across reps for ribbon
  agg_col <- if ("Method" %in% names(traj) && colour_by == "Method")
    "Method" else "Rep"

  agg <- traj[, .(
    mu     = mean(MeanGEBV, na.rm = TRUE),
    lo     = mean(MeanGEBV, na.rm = TRUE) - stats::sd(MeanGEBV, na.rm = TRUE),
    hi     = mean(MeanGEBV, na.rm = TRUE) + stats::sd(MeanGEBV, na.rm = TRUE)
  ), by = c("Cycle", agg_col)]

  p_gebv <- ggplot(agg, aes(x = Cycle, y = mu,
                             colour = .data[[agg_col]],
                             fill   = .data[[agg_col]])) +
    geom_line(linewidth = 1) +
    labs(
      title  = paste0("GEBV Trajectories \u2014 ", trait_name),
      x      = "Breeding Cycle",
      y      = "Mean Local GEBV",
      colour = agg_col, fill = agg_col
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  if (show_ribbon) {
    p_gebv <- p_gebv +
      ggplot2::geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2,
                           colour = NA)
  }

  if (!show_coverage || !"BlockCoverage" %in% names(traj)) return(p_gebv)

  # Block coverage panel
  agg_cov <- traj[, .(
    cov_mu = mean(BlockCoverage, na.rm = TRUE)
  ), by = c("Cycle", agg_col)]

  p_cov <- ggplot(agg_cov, aes(x = Cycle, y = cov_mu,
                                colour = .data[[agg_col]])) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title  = "Block Coverage",
      x      = "Breeding Cycle",
      y      = "Proportion of Blocks with \u22651 Positive Carrier",
      colour = agg_col
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  if (requireNamespace("cowplot", quietly = TRUE)) {
    cowplot::plot_grid(p_gebv, p_cov, nrow = 2L, align = "v")
  } else {
    .hs_inform("Install 'cowplot' for a combined two-panel plot. Returning GEBV plot only.")
    p_gebv
  }
}


#' Summarise simulation results across replications
#'
#' Computes per-cycle summary statistics (mean, SD, min, max) for GEBV and
#' block coverage across all replications, optionally grouped by selection
#' method.
#'
#' @param trajectories data.table from \code{\link{run_stacking_simulation}}.
#' @param group_by     Character vector of grouping columns. Default
#'   \code{c("Cycle","Method")} when a \code{Method} column is present,
#'   otherwise \code{"Cycle"}.
#'
#' @return A \code{data.table} with one row per (Cycle × Method) combination
#'   and columns for mean/SD/min/max of MeanGEBV and BlockCoverage.
#'
#' @export
summarise_simulation <- function(trajectories,
                                  group_by = NULL) {

  traj <- as.data.table(trajectories)

  if (is.null(group_by)) {
    group_by <- if ("Method" %in% names(traj)) c("Cycle","Method")
                else "Cycle"
  }

  metrics <- intersect(c("MeanGEBV","SDev","BlockCoverage"), names(traj))
  # A metric with no value at all (BlockCoverage is all NA in the compat view)
  # is skipped (0.9.0, sc/F23 / I02): up to 0.8.8 it gave +-Inf columns and
  # two "no non-missing arguments" warnings per group.
  metrics <- metrics[vapply(metrics, function(m) any(!is.na(traj[[m]])), logical(1L))]
  if (!length(metrics))
    .hs_abort("summarise_simulation: no metric column with a value ",
              "(MeanGEBV, SDev, BlockCoverage are all missing or all NA).")

  result_list <- lapply(metrics, function(m) {
    traj[, .(
      mean_val = mean(get(m), na.rm = TRUE),
      sd_val   = stats::sd(get(m), na.rm = TRUE),
      min_val  = min(get(m), na.rm = TRUE),
      max_val  = max(get(m), na.rm = TRUE)
    ), by = group_by] |>
      setnames(c(group_by, paste0(m, c("_mean","_sd","_min","_max"))))
  })

  Reduce(function(a, b) merge(a, b, by = group_by), result_list)
}


# ---------------------------------------------------------------------------
# Metrics trajectories on the origin-compatible summary (FR5)
# ---------------------------------------------------------------------------

#' Plot per-cycle simulation metrics (F1mean / D2Ideal / coancestry)
#'
#' Faceted trajectories of the FR5 cycle metrics directly from the
#' \code{summary} table returned by \code{\link{run_stacking_simulation}} (the
#' origin-compatible output), averaged across replications and coloured by
#' selection scheme. This complements \code{\link{plot_simulation_trajectories}}
#' (which plots only mean GEBV from the back-compat \code{trajectories} view) by
#' surfacing the ideal-genotype distance and coancestry alongside genetic gain.
#'
#' @param summary The \code{summary} data.table from
#'   \code{\link{run_stacking_simulation}} (columns \code{cycle, sel_scheme,
#'   mean_F1mean, mean_D2Ideal, mean_coancestry, ...}).
#' @param metrics Character vector of summary columns to panel. Default
#'   \code{c("mean_F1mean", "mean_D2Ideal", "mean_coancestry",
#'   "mean_gene_diversity")}; any subset of the numeric \code{mean_*} columns is
#'   accepted (missing columns are skipped).
#' @param include_founder Logical. Keep the \code{cycle = 0} founder row.
#'   Default \code{FALSE} (start panels at the first decision cycle).
#'
#' @return A \code{ggplot} object (metrics faceted, free y-scales).
#'
#' @seealso \code{\link{plot_simulation_trajectories}},
#'   \code{\link{run_stacking_simulation}}
#' @export
plot_simulation_metrics <- function(
    summary,
    metrics = c("mean_F1mean", "mean_D2Ideal", "mean_coancestry",
                "mean_gene_diversity"),
    include_founder = FALSE) {

  s <- as.data.table(summary)
  if (!all(c("cycle", "sel_scheme") %in% names(s)))
    .hs_abort("plot_simulation_metrics: 'summary' needs columns 'cycle' and ",
         "'sel_scheme' (pass run_stacking_simulation()$summary).")

  metrics <- intersect(metrics, names(s))
  if (length(metrics) == 0L)
    .hs_abort("plot_simulation_metrics: none of the requested 'metrics' are in ",
         "'summary'.")
  if (!include_founder) s <- s[cycle > 0L]

  long <- data.table::melt(
    s, id.vars = c("cycle", "sel_scheme"), measure.vars = metrics,
    variable.name = "metric", value.name = "value")
  agg <- long[, .(mu = mean(value, na.rm = TRUE)),
              by = .(cycle, sel_scheme, metric)]

  ggplot(agg, aes(x = cycle, y = mu, colour = sel_scheme)) +
    geom_line(linewidth = 1) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    labs(x = "Breeding cycle", y = "Mean (across replications)",
         colour = "Scheme") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
}


# ---------------------------------------------------------------------------
# Internal — live-monitor render (issue 08)
# ---------------------------------------------------------------------------

# Re-render the accumulating simulation record + figures after each completed
# (rep, scheme) job so a long run shows progress mid-flight. Called by
# run_stacking_simulation(monitor = TRUE) with the records collected so far;
# strictly observation-only (it reads accumulated outputs, never touches
# selection or the return value). Writes into `output_dir`:
#   live_summary.txt / live_plans.txt  — the partial FR5 tables
#   live_metrics.png                    — plot_simulation_metrics()
#   live_trajectories.png               — plot_simulation_trajectories()
# Every write/render is guarded so a transient plotting/IO hiccup degrades to a
# skipped refresh rather than aborting the whole simulation — and, since
# 0.9.0, says so: each failure is a note (hs_warning_fallback, code
# monitor_write_failed) that the driver's result carries. PNGs go through
# hs_render_figure() (the one device opener) at the 0.8.8 sizes.
.stk_render_live <- function(summary_acc, plans_acc, output_dir, tag = "sim") {
  if (is.null(output_dir)) return(invisible(NULL))
  summ <- as.data.table(summary_acc)
  plns <- as.data.table(plans_acc)
  guard <- function(what, expr) tryCatch(expr, error = function(e) {
    .hs_note("monitor", "monitor_write_failed",
             paste0("live monitor: ", what, " not written (", conditionMessage(e), ")."))
    NULL
  })

  guard("live_summary.txt",
        data.table::fwrite(summ, file.path(output_dir, "live_summary.txt"), sep = "\t"))
  guard("live_plans.txt",
        data.table::fwrite(plns, file.path(output_dir, "live_plans.txt"), sep = "\t"))

  # Metrics panels straight off the origin summary contract. Only panel the
  # metrics that actually carry data this run (unrequested ones are all-NA), so
  # the live figure has no empty facets.
  guard("live_metrics.png", {
    cand <- c("mean_F1mean", "mean_D2Ideal", "mean_coancestry",
              "mean_gene_diversity")
    have <- cand[vapply(cand, function(cn)
      cn %in% names(summ) && any(is.finite(summ[[cn]])), logical(1L))]
    if (length(have)) {
      pm <- plot_simulation_metrics(summ, metrics = have, include_founder = FALSE)
      suppressWarnings(suppressMessages(hs_render_figure(
        pm, file.path(output_dir, "live_metrics.png"),
        width_mm = 254, height_mm = 152.4, dpi = 150, overwrite = TRUE)))
    }
  })

  # GEBV trajectories via the back-compat view (the one transform the final
  # combine uses too, .stk_trajectories_view()).
  guard("live_trajectories.png", {
    traj <- .stk_trajectories_view(summ)
    if (nrow(traj) > 0L) {
      pt <- plot_simulation_trajectories(traj, trait_name = tag)
      suppressWarnings(suppressMessages(hs_render_figure(
        pt, file.path(output_dir, "live_trajectories.png"),
        width_mm = 254, height_mm = 127, dpi = 150, overwrite = TRUE)))
    }
  })

  invisible(NULL)
}
