# =============================================================================
# catalogue_cluster.R
#
# Sample (genotype-level) clustering from haplotype-level allele calls, plus the
# portal that feeds a clustering into the haplotype-alignment plot.
#
# The distance is computed on the HAPLOTYPE-level allele calls (one row per
# haplotype copy) but the clustering LEAVES are whole samples: the ploidy copies
# of a sample are collapsed into a single genotype-level distance. So a user
# clusters "by the haplotypes a sample carries" yet reads the result — dendrogram
# leaves, cluster assignments, alignment row order — at the sample level.
#
# Distance (user-selected): marker-level Hamming (matching) distance, missing
# excluded pairwise, computed with a single BLAS `tcrossprod` over one-hot allele
# indicators (agreements) and over the non-missing mask (comparable markers), so
# the hot path is C/BLAS rather than an R pair loop — this is the ~10K-marker
# scale lever. No new package dependencies (stats + base only for the numerics).
#
# Public:
#   cluster_samples()            — hierarchical clustering -> hs_sample_clustering
#   plot_cluster_dendrogram()    — ggplot dendrogram panel (leaves = samples)
#   plot_clustered_alignment()   — cluster -> row-ordered alignment (+ dendrogram)
#   print.hs_sample_clustering()
# =============================================================================


# ---------------------------------------------------------------------------
# Internal — integer-coded haplotype-copy matrix
# ---------------------------------------------------------------------------

# Build the [n_copies x n_markers] integer allele-code matrix used for distance.
# Rows are haplotype copies (Sample x Haplotype_idx); NA marks a missing/absent
# allele (encoding-aware, via .missing_tokens). Mirrors the integer-encoding step
# of .catalogue_core() but over ALL requested markers, so clustering works with or
# without haplotype blocks attached.
# Returns list(code = matrix, copy_sample = chr, copy_hap = int, markers = chr,
#              samples = chr, encoding = chr, ploidy = int).
# @keywords internal
.cluster_copy_codes <- function(geno_meta) {
  prof    <- geno_meta_views(geno_meta)
  alleles <- prof$alleles                                # Sample, Haplotype_idx, <markers...>
  markers <- setdiff(names(alleles), c("Sample", "Haplotype_idx"))
  if (nrow(alleles) == 0L || length(markers) == 0L)
    .hs_abort("cluster_samples: no genotype data to cluster after filtering.")

  miss  <- .missing_tokens(prof$encoding, .meta_extra_missing(geno_meta))
  n_cp  <- nrow(alleles)
  n_mk  <- length(markers)
  code  <- matrix(NA_integer_, n_cp, n_mk, dimnames = list(NULL, markers))
  for (j in seq_len(n_mk)) {
    v <- as.character(alleles[[markers[j]]])
    v[is.na(v) | !nzchar(trimws(v)) | v %in% miss] <- NA_character_
    code[, j] <- as.integer(factor(v))                  # per-marker allele codes
  }
  list(code        = code,
       copy_sample = as.character(alleles$Sample),
       copy_hap    = as.integer(alleles$Haplotype_idx),
       markers     = markers,
       samples     = prof$samples,
       encoding    = prof$encoding,
       ploidy      = max(alleles$Haplotype_idx, na.rm = TRUE))
}


# Matching (Hamming) distance between the ROWS of an integer code matrix, missing
# (NA) excluded pairwise: d(i,j) = 1 - agreements / comparable over markers scored
# in both rows. Agreements and comparable counts are each a single `tcrossprod`
# (BLAS), giving an O(n^2) memory / BLAS-bound hot path rather than an R loop over
# pairs. Cells with no comparable marker are NA. Returns a full symmetric matrix
# with a zero diagonal.
# @keywords internal
.matching_distance <- function(code) {
  n <- nrow(code); m <- ncol(code)
  nm <- !is.na(code)                                     # [n x m] non-missing mask
  comparable <- tcrossprod(nm * 1)                       # both scored, per pair

  # One-hot indicator over (marker, allele-code): a 1 at column off[j]+code for
  # every scored (copy, marker). tcrossprod counts, per pair, the markers where
  # BOTH carry the same code == agreements.
  kj  <- apply(code, 2L, function(c) { mx <- suppressWarnings(max(c, na.rm = TRUE))
                                       if (is.finite(mx)) mx else 0L })
  off <- c(0L, cumsum(kj)[-m])
  P   <- sum(kj)
  agree <- matrix(0, n, n)
  if (P > 0L) {
    idx <- which(nm, arr.ind = TRUE)                     # (row=copy, col=marker)
    col_g <- off[idx[, 2L]] + code[nm]                   # global one-hot column
    onehot <- matrix(0, n, P)
    onehot[cbind(idx[, 1L], col_g)] <- 1
    agree <- tcrossprod(onehot)
  }

  d <- 1 - agree / comparable
  d[comparable == 0] <- NA_real_
  diag(d) <- 0
  d
}


# Collapse a copy-level distance matrix to a sample-level one. "average" is the
# masked mean over the ploidy x ploidy copy pairings between two samples (NA copy
# pairs skipped); "best" is the minimum copy-pair distance. Returns a symmetric
# sample x sample matrix (zero diagonal) with dimnames = samples.
# @keywords internal
.merge_copy_distance <- function(dcopy, copy_sample, method = c("average", "best")) {
  method  <- match.arg(method)
  samples <- unique(copy_sample)
  G <- outer(copy_sample, samples, `==`) * 1            # [n_copies x n_samples] 0/1

  if (method == "average") {
    dz   <- dcopy; ok <- is.finite(dz); dz[!ok] <- 0
    num  <- crossprod(G, dz %*% G)                       # sum of comparable copy dists
    den  <- crossprod(G, (ok * 1) %*% G)                 # # comparable copy pairs
    ds   <- num / den
    ds[den == 0] <- NA_real_
  } else {                                               # "best" = min copy-pair dist
    ns <- length(samples)
    ds <- matrix(NA_real_, ns, ns)
    grp <- lapply(samples, function(s) which(copy_sample == s))
    for (a in seq_len(ns)) for (b in a:ns) {
      block <- dcopy[grp[[a]], grp[[b]], drop = FALSE]
      v <- suppressWarnings(min(block[is.finite(block)]))
      ds[a, b] <- ds[b, a] <- if (is.finite(v)) v else NA_real_
    }
  }
  diag(ds) <- 0
  dimnames(ds) <- list(samples, samples)
  ds
}


# "pool" merge: collapse each sample's copies to one per-marker genotype profile
# first (sorted-allele string via .collapse_genotype_sorted, diploid), then a
# single sample-level matching distance. Falls back to the copy path for non-
# diploid data (handled by the caller). Returns a sample x sample matrix.
# @keywords internal
.pool_sample_distance <- function(geno_meta, samples) {
  arr <- geno_meta$geno_arrays
  gmat <- .collapse_genotype_sorted(arr)                 # [n_samples x n_markers] strings
  gmat <- gmat[samples, , drop = FALSE]
  code <- matrix(NA_integer_, nrow(gmat), ncol(gmat), dimnames = list(rownames(gmat), NULL))
  for (j in seq_len(ncol(gmat))) {
    v <- as.character(gmat[, j]); v[is.na(v) | !nzchar(trimws(v))] <- NA_character_
    code[, j] <- as.integer(factor(v))
  }
  d <- .matching_distance(code)
  dimnames(d) <- list(rownames(gmat), rownames(gmat))
  d
}


# ---------------------------------------------------------------------------
# cluster_samples()
# ---------------------------------------------------------------------------

#' Hierarchically cluster samples from their haplotype-level alleles
#'
#' Clusters the samples of a \code{GenoMeta} using a distance computed on the
#' \strong{haplotype-level} allele calls (one row per haplotype copy), collapsed to
#' the \strong{genotype (sample) level} so that every dendrogram leaf is one sample.
#'
#' The distance is a marker-level Hamming (matching) distance — the fraction of
#' jointly-scored markers whose alleles differ, missing calls excluded pairwise —
#' computed with a single BLAS \code{tcrossprod} over one-hot allele indicators, so
#' it scales to panels of many thousands of markers without a per-pair R loop. The
#' ploidy copies of a sample are then merged into one per-sample distance
#' (\code{copy_merge}). Hierarchical clustering uses \code{stats::hclust}.
#'
#' @param geno_meta A valid \code{GenoMeta} list.
#' @param samples,markers Optional character vectors restricting the samples /
#'   markers used (passed to \code{\link{filter_geno_meta}} as \code{keep_samples}
#'   / \code{keep_markers}).
#' @param blocks Optional block IDs restricting the markers used (block vocabulary
#'   of \code{\link{filter_geno_meta}}).
#' @param chromosomes Optional chromosome-region list (see
#'   \code{\link{filter_geno_meta}}).
#' @param distance Distance metric. Currently \code{"hamming"} (marker-level
#'   matching distance); reserved for future metrics.
#' @param copy_merge How the ploidy haplotype copies of a sample collapse into one
#'   sample distance: \code{"average"} (default; masked mean of the
#'   copy-vs-copy distances), \code{"pool"} (collapse copies to one genotype
#'   profile first; diploid only, else falls back to \code{"average"}), or
#'   \code{"best"} (minimum copy-pair distance).
#' @param linkage Agglomeration method passed to \code{stats::hclust}. Default
#'   \code{"ward.D2"}.
#' @param k Optional number of clusters to cut the tree into (\code{stats::cutree});
#'   when supplied, the result carries a per-sample \code{clusters} assignment.
#' @param verbose Logical; print a progress line via \code{message()}.
#'
#' @return An \code{hs_sample_clustering} object: a list with \code{order} (sample
#'   IDs in dendrogram leaf order), \code{hclust}, \code{dist} (sample-level),
#'   \code{clusters} (named sample -> cluster integer vector, or \code{NULL}),
#'   \code{params}, \code{n_samples}, \code{n_markers}, and \code{call}.
#'
#' @examples
#' \dontrun{
#' meta <- sim_geno_meta(n_samples = 30, n_markers = 200, block_size = 5)
#' cl   <- cluster_samples(meta, linkage = "ward.D2", k = 4)
#' plot_haplotype_alignment(meta, method = "genotype", sample_order = cl)
#' plot_clustered_alignment(meta, cl)
#' }
#'
#' @seealso \code{\link{plot_clustered_alignment}},
#'   \code{\link{plot_cluster_dendrogram}}, \code{\link{plot_haplotype_alignment}}
#' @export
cluster_samples <- function(geno_meta, samples = NULL, markers = NULL,
                            blocks = NULL, chromosomes = NULL,
                            distance   = "hamming",
                            copy_merge = c("average", "pool", "best"),
                            linkage    = "ward.D2",
                            k          = NULL,
                            verbose    = FALSE) {
  .hs_check_meta(geno_meta, "cluster_samples", "geno_meta")
  distance   <- match.arg(distance, c("hamming"))
  copy_merge <- match.arg(copy_merge)

  gm <- geno_meta
  if (!is.null(samples) || !is.null(markers) || !is.null(blocks) ||
      !is.null(chromosomes))
    gm <- filter_geno_meta(gm, keep_markers = markers, keep_samples = samples,
                           blocks = blocks, chromosomes = chromosomes)

  cc <- .cluster_copy_codes(gm)
  if (length(cc$samples) < 2L)
    .hs_abort("cluster_samples: need at least 2 samples to cluster (got ",
         length(cc$samples), ").")

  if (copy_merge == "pool" && cc$ploidy == 2L) {
    dsamp <- .pool_sample_distance(gm, cc$samples)
    used_merge <- "pool"
  } else {
    if (copy_merge == "pool") {
      .hs_warn("cluster_samples: copy_merge = 'pool' needs diploid data; ",
              "falling back to 'average'.")
      used_merge <- "average"
    } else used_merge <- copy_merge
    dcopy <- .matching_distance(cc$code)
    dsamp <- .merge_copy_distance(dcopy, cc$copy_sample,
                                  method = if (used_merge == "best") "best" else "average")
  }

  if (anyNA(dsamp)) {
    n_na <- sum(is.na(dsamp[upper.tri(dsamp)]))
    .hs_warn("cluster_samples: ", n_na, " sample pair(s) share no comparable ",
            "marker; their distance was set to the maximum observed.")
    dsamp[is.na(dsamp)] <- max(dsamp[is.finite(dsamp)], 0)
    diag(dsamp) <- 0
  }

  d  <- stats::as.dist(dsamp)
  hc <- stats::hclust(d, method = linkage)
  hc$labels <- cc$samples
  ordn <- hc$labels[hc$order]

  clusters <- NULL
  if (!is.null(k)) {
    k <- as.integer(k)
    if (k < 1L || k > length(cc$samples))
      .hs_abort("cluster_samples: 'k' must be between 1 and n_samples (",
           length(cc$samples), ").")
    clusters <- stats::cutree(hc, k = k)
    names(clusters) <- cc$samples
  }

  if (isTRUE(verbose))
    .hs_inform(sprintf("[cluster_samples] %d samples x %d markers | distance=%s | ",
                    length(cc$samples), length(cc$markers), distance),
            sprintf("copy_merge=%s | linkage=%s%s", used_merge, linkage,
                    if (!is.null(k)) sprintf(" | k=%d", k) else ""))

  structure(list(
    order     = ordn,
    hclust    = hc,
    dist      = d,
    clusters  = clusters,
    params    = list(distance = distance, copy_merge = used_merge,
                     linkage = linkage, k = k, n_markers_used = length(cc$markers)),
    n_samples = length(cc$samples),
    n_markers = length(cc$markers),
    call      = match.call()),
    class = "hs_sample_clustering")
}


#' @export
print.hs_sample_clustering <- function(x, ...) {
  cat(sprintf("<hs_sample_clustering> %d samples, %d markers\n",
              x$n_samples, x$n_markers))
  cat(sprintf("  distance=%s  copy_merge=%s  linkage=%s%s\n",
              x$params$distance, x$params$copy_merge, x$params$linkage,
              if (!is.null(x$params$k)) sprintf("  k=%d", x$params$k) else ""))
  cat(sprintf("  leaf order: %s%s\n",
              paste(utils::head(x$order, 6L), collapse = ", "),
              if (length(x$order) > 6L) ", ..." else ""))
  invisible(x)
}


# ---------------------------------------------------------------------------
# Dendrogram segments (base R; no ggdendro dependency)
# ---------------------------------------------------------------------------

# Turn an hclust into leaf positions + link segments. Leaves are placed at
# integer positions 1..n in dendrogram order; each internal node sits at the mean
# position of its children and at height = merge height. Returns
# list(segments = data.frame(x, y, xend, yend), leaves = data.frame(label, pos)).
# @keywords internal
.hclust_segments <- function(hc) {
  n     <- length(hc$order)
  pos   <- numeric(n + nrow(hc$merge))     # x-position of each node
  hgt   <- numeric(n + nrow(hc$merge))     # height (y) of each node
  leaf_pos <- match(seq_len(n), hc$order)  # leaf i -> its slot in the ordering
  node_x <- function(id) if (id < 0L) leaf_pos[-id] else pos[n + id]
  node_h <- function(id) if (id < 0L) 0 else hgt[n + id]

  segs <- list()
  for (m in seq_len(nrow(hc$merge))) {
    a <- hc$merge[m, 1L]; b <- hc$merge[m, 2L]
    xa <- node_x(a); xb <- node_x(b); ha <- node_h(a); hb <- node_h(b)
    h  <- hc$height[m]
    pos[n + m] <- (xa + xb) / 2
    hgt[n + m] <- h
    # two verticals up to the join height + the horizontal connector
    segs[[length(segs) + 1L]] <- data.frame(x = xa, y = ha, xend = xa, yend = h)
    segs[[length(segs) + 1L]] <- data.frame(x = xb, y = hb, xend = xb, yend = h)
    segs[[length(segs) + 1L]] <- data.frame(x = xa, y = h,  xend = xb, yend = h)
  }
  list(segments = do.call(rbind, segs),
       leaves   = data.frame(label = hc$labels[hc$order],
                             pos   = seq_len(n), stringsAsFactors = FALSE))
}


#' Plot a sample-clustering dendrogram
#'
#' Renders the dendrogram of a \code{\link{cluster_samples}} result as a
#' \code{ggplot}, with leaves ordered to line up with a haplotype alignment's
#' sample rows. Intended as the side panel of \code{\link{plot_clustered_alignment}}
#' but usable on its own.
#'
#' @param x An \code{hs_sample_clustering} (from \code{\link{cluster_samples}}).
#' @param orientation \code{"horizontal"} (default; leaves down the y-axis, tree
#'   growing left, to sit left of a row-based alignment) or \code{"vertical"}
#'   (leaves along the x-axis).
#' @param k Optional number of clusters for colouring the leaves/branches; defaults
#'   to the \code{k} stored on \code{x} (if any).
#' @param color_clusters Logical; colour by cluster when \code{k} is available.
#' @param line_width Dendrogram line width. Default \code{0.4}.
#' @param label_leaves Logical; draw leaf (sample) tick labels. Default
#'   \code{FALSE} (the paired alignment panel carries the labels).
#' @param leaf_expand Numeric length-2 additive expansion \code{c(low, high)} of
#'   the leaf axis, in leaf units. Default \code{c(0.5, 0.5)}. Raise the upper
#'   value to match a companion panel that reserves headroom — this is how
#'   \code{\link{plot_clustered_alignment}} keeps the leaves level with the
#'   alignment rows once the alignment reserves room for its block labels.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{cluster_samples}}, \code{\link{plot_clustered_alignment}}
#' @export
plot_cluster_dendrogram <- function(x, orientation = c("horizontal", "vertical"),
                                    k = NULL, color_clusters = TRUE,
                                    line_width = 0.4, label_leaves = FALSE,
                                    leaf_expand = c(0.5, 0.5)) {
  if (!inherits(x, "hs_sample_clustering"))
    .hs_abort("plot_cluster_dendrogram: 'x' must be an hs_sample_clustering ",
         "(from cluster_samples()).")
  orientation <- match.arg(orientation)
  leaf_expand <- as.numeric(leaf_expand)
  if (length(leaf_expand) == 1L) leaf_expand <- rep(leaf_expand, 2L)
  if (length(leaf_expand) != 2L || anyNA(leaf_expand))
    .hs_abort("plot_cluster_dendrogram: 'leaf_expand' must be one or two finite ",
         "numbers.")
  leaf_exp <- ggplot2::expansion(mult = 0, add = leaf_expand)
  seg <- .hclust_segments(x$hclust)
  segments <- seg$segments; leaves <- seg$leaves

  k <- k %||% x$params$k
  leaf_col <- NULL
  if (isTRUE(color_clusters) && !is.null(k)) {
    cl <- stats::cutree(x$hclust, k = k)
    leaves$cluster <- factor(cl[leaves$label])
    leaf_col <- leaves
  }

  p <- ggplot2::ggplot()
  if (orientation == "horizontal") {
    # leaf position -> y; height -> x (tree grows toward -x, leaves at x = 0)
    p <- p + ggplot2::geom_segment(
      data = segments,
      ggplot2::aes(x = -y, y = x, xend = -yend, yend = xend),
      linewidth = line_width, colour = "grey25")
    if (!is.null(leaf_col))
      p <- p + ggplot2::geom_point(
        data = leaf_col, ggplot2::aes(x = 0, y = pos, colour = cluster), size = 1.6)
    if (isTRUE(label_leaves))
      p <- p + ggplot2::scale_y_continuous(
        breaks = leaves$pos, labels = leaves$label, expand = leaf_exp)
    else
      p <- p + ggplot2::scale_y_continuous(expand = leaf_exp) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
    p <- p + ggplot2::labs(x = "Height", y = NULL)
  } else {
    p <- p + ggplot2::geom_segment(
      data = segments,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = line_width, colour = "grey25")
    if (!is.null(leaf_col))
      p <- p + ggplot2::geom_point(
        data = leaf_col, ggplot2::aes(x = pos, y = 0, colour = cluster), size = 1.6)
    if (isTRUE(label_leaves))
      p <- p + ggplot2::scale_x_continuous(
        breaks = leaves$pos, labels = leaves$label, expand = leaf_exp)
    else
      p <- p + ggplot2::scale_x_continuous(expand = leaf_exp)
    p <- p + ggplot2::labs(x = NULL, y = "Height")
  }

  p + ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = "Cluster",
                                 na.translate = FALSE)
}


#' Cluster samples and plot the row-ordered haplotype alignment (with dendrogram)
#'
#' The one-call portal from \code{\link{cluster_samples}} to
#' \code{\link{plot_haplotype_alignment}}: it orders the alignment's sample rows by
#' a clustering and, optionally, places the clustering dendrogram as a side panel
#' whose leaves line up with the alignment rows.
#'
#' @param geno_meta A valid \code{GenoMeta} list.
#' @param clustering Optional \code{hs_sample_clustering}; computed with
#'   \code{\link{cluster_samples}} (defaults) when \code{NULL}.
#' @param method Alignment method: \code{"genotype"} (default, one row per sample —
#'   required for the dendrogram panel) or \code{"haplotype"} (rows grouped by the
#'   clustered sample order; no dendrogram panel).
#' @param show_dendrogram Logical; compose the dendrogram side panel (needs the
#'   \pkg{cowplot} package and \code{method = "genotype"}).
#' @param rel_widths Relative widths of \code{c(dendrogram, alignment)}. Default
#'   \code{c(1, 4)}.
#' @param k Optional number of clusters (passed through to
#'   \code{\link{cluster_samples}} / colouring) when \code{clustering} is computed
#'   here.
#' @param ... Further arguments forwarded to \code{\link{plot_haplotype_alignment}}.
#'
#' @return A \code{ggplot} (no dendrogram) or a composed \pkg{cowplot} object
#'   (with dendrogram), saveable with \code{ggplot2::ggsave}.
#'
#' @seealso \code{\link{cluster_samples}}, \code{\link{plot_cluster_dendrogram}},
#'   \code{\link{plot_haplotype_alignment}}, \code{\link{plot_meta_alignment}}
#' @export
plot_clustered_alignment <- function(geno_meta, clustering = NULL,
                                     method = c("genotype", "haplotype"),
                                     show_dendrogram = TRUE, rel_widths = c(1, 4),
                                     k = NULL, ...) {
  .hs_check_meta(geno_meta, "plot_clustered_alignment", "geno_meta")
  method <- match.arg(method)
  if (is.null(clustering)) {
    # A broadcast clustering is used in preference to computing a fresh one:
    # otherwise this panel would cluster independently of every other figure in
    # the session and quietly disagree with them. A broadcast that is only a
    # character order cannot drive the dendrogram, so that case still computes.
    bcast <- .hs_sample_order()
    if (inherits(bcast, "hs_sample_clustering")) {
      clustering <- bcast
      .hs_inform("plot_clustered_alignment: using the clustering broadcast by ",
              "use_sample_order(); pass 'clustering' to override.")
    } else {
      clustering <- cluster_samples(geno_meta, k = k)
    }
  }
  if (!inherits(clustering, "hs_sample_clustering"))
    .hs_abort("plot_clustered_alignment: 'clustering' must be an hs_sample_clustering.")

  aln <- plot_haplotype_alignment(geno_meta, method = method,
                                  sample_order = clustering, ...)

  want_dendro <- isTRUE(show_dendrogram) && method == "genotype"
  # viz/F16 rule (0.9.0): a missing optional package is an error that names the
  # install, never a silently different return class.
  if (want_dendro && !requireNamespace("cowplot", quietly = TRUE))
    .hs_abort("plot_clustered_alignment: the 'cowplot' package is required for ",
              "show_dendrogram = TRUE. Install it, or pass show_dendrogram = FALSE.",
              class = "hs_error_dependency",
              data = list(hint = "install.packages(\"cowplot\")"))
  if (!want_dendro) {
    if (isTRUE(show_dendrogram) && method == "haplotype")
      .hs_inform("plot_clustered_alignment: dendrogram panel is genotype-level; ",
              "returning the haplotype alignment ordered by the clustering.")
    return(aln)
  }

  # Mirror the alignment's own Y expansion onto the dendrogram's leaf axis: the
  # alignment reserves headroom above the top row for its block labels, and
  # without the same reservation the leaves drift upward against the rows.
  den <- plot_cluster_dendrogram(clustering, orientation = "horizontal", k = k,
                                 leaf_expand = attr(aln, "row_expand") %||%
                                               c(0.5, 0.5))
  # Alignment rows run bottom->top over sample_order == clustering$order, i.e. the
  # dendrogram leaf order, so leaf position i already matches alignment row i.
  plot_meta_alignment(list(den, aln), rel_widths = rel_widths,
                      sync_axis = "none", legend = "bottom", labels_on = 2L)
}
