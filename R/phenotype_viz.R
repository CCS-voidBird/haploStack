# =============================================================================
# phenotype_viz.R
#
# Portable visualisation of sample-level phenotype data.
#
# The plotting functions accept EITHER a GenoMeta (their `$phenotype` element is
# used) OR any data.frame / object carrying a `Sample` column plus numeric trait
# columns — so they work on a stored phenotype, on a bare phenotype table from
# read_phenotype(), or on any other alignment-style table that contains sample
# data. This portability is provided by the internal .as_phenotype_table() seam.
#
# Functions:
#   plot_phenotype_distribution()  — per-trait histogram / density / boxplot
#   plot_phenotype_correlation()   — trait correlation heatmap
#   plot_phenotype_bars()          — per-sample dodged bar chart (<=2 traits)
#   phenotype_summary()            — per-trait numeric summary table
#   plot_meta_alignment()          — merge plots that share a discrete sample axis
#   .as_phenotype_table()          — internal: coerce input to long Sample/Trait/Value
#   .discrete_axis_levels()        — internal: read a ggplot's discrete axis order
# =============================================================================


# Internal: coerce a phenotype-bearing input to a long data.frame with columns
# Sample, Trait, Value (+ a wide attribute for correlation).
#
#   x       a GenoMeta list (uses $phenotype), or any data.frame/object with a
#           Sample column + numeric trait columns.
#   traits  optional character vector selecting / ordering trait columns.
#
# Returns a long data.frame (Sample, Trait, Value); the wide table is attached as
# attr(., "wide") for callers that need a sample × trait matrix (e.g.
# correlation). Implemented in base R (no data.table `:=`) so the plotting helpers
# work regardless of whether the caller is in a data.table-aware namespace.
.as_phenotype_table <- function(x, traits = NULL) {
  # Resolve a wide Sample + traits data.frame from the input.
  if (is.list(x) && !is.data.frame(x) && !is.null(x[["phenotype"]])) {
    wide <- x[["phenotype"]]
  } else if (is.list(x) && !is.data.frame(x) &&
             ("geno_arrays" %in% names(x))) {
    stop(".as_phenotype_table: GenoMeta has no phenotype attached. ",
         "Use attach_phenotype() first.")
  } else if (is.data.frame(x)) {
    wide <- x
  } else {
    stop(".as_phenotype_table: 'x' must be a GenoMeta or a data.frame with a ",
         "'Sample' column.")
  }

  wide <- as.data.frame(wide, stringsAsFactors = FALSE)
  sample_col <- if ("Sample" %in% names(wide)) "Sample" else names(wide)[1L]

  trait_cols <- setdiff(names(wide), sample_col)
  if (!is.null(traits)) {
    missing_t <- setdiff(traits, trait_cols)
    if (length(missing_t) > 0L)
      stop(".as_phenotype_table: requested trait(s) not found: ",
           paste(missing_t, collapse = ", "), ".")
    trait_cols <- traits
  }
  # Keep only numeric trait columns.
  is_num <- vapply(trait_cols, function(tc) is.numeric(wide[[tc]]) ||
                     !anyNA(suppressWarnings(as.numeric(wide[[tc]]))),
                   logical(1L))
  trait_cols <- trait_cols[is_num]
  if (length(trait_cols) == 0L)
    stop(".as_phenotype_table: no numeric trait columns to plot.")

  w <- data.frame(Sample = as.character(wide[[sample_col]]),
                  stringsAsFactors = FALSE)
  for (tc in trait_cols) w[[tc]] <- suppressWarnings(as.numeric(wide[[tc]]))

  # Stack to long form in base R (column-major: trait blocks of nrow(w) rows).
  long <- data.frame(
    Sample = rep(w$Sample, times = length(trait_cols)),
    Trait  = factor(rep(trait_cols, each = nrow(w)), levels = trait_cols),
    Value  = unlist(w[trait_cols], use.names = FALSE),
    stringsAsFactors = FALSE
  )
  attr(long, "wide") <- w
  long
}


#' Plot per-trait phenotype distributions
#'
#' Faceted distribution plot of every (selected) trait. Portable: \code{x} may be
#' a \code{GenoMeta} (its \code{$phenotype} is used) or any data.frame with a
#' \code{Sample} column plus numeric trait columns.
#'
#' @param x      A \code{GenoMeta} with an attached phenotype, or a phenotype
#'   data.frame (\code{Sample} + traits).
#' @param traits Optional character vector selecting / ordering traits.
#'   \code{NULL} (default) uses all numeric trait columns.
#' @param type   \code{"histogram"} (histogram + density overlay, the default),
#'   \code{"density"}, or \code{"box"} (boxplot with jittered points).
#' @param bins   Number of histogram bins (\code{type = "histogram"}).
#'   Default \code{30}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' p <- plot_phenotype_distribution(meta, type = "histogram")
#' ggplot2::ggsave("trait_distributions.png", p, width = 10, height = 6)
#' }
#'
#' @seealso \code{\link{plot_phenotype_correlation}}, \code{\link{phenotype_summary}}
#'
#' @export
plot_phenotype_distribution <- function(x, traits = NULL,
                                        type = c("histogram", "density", "box"),
                                        bins = 30L) {
  type <- match.arg(type)
  long <- .as_phenotype_table(x, traits = traits)
  df   <- as.data.frame(long)

  if (type == "box") {
    p <- ggplot(df, aes(x = Trait, y = Value, fill = Trait)) +
      ggplot2::geom_boxplot(outlier.shape = 21, outlier.fill = "white",
                            outlier.colour = "grey40", width = 0.5) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.15, size = 0.8) +
      labs(x = NULL, y = "Value") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold"))
    return(p)
  }

  p <- ggplot(df, aes(x = Value))
  if (type == "histogram") {
    p <- p +
      ggplot2::geom_histogram(aes(y = ggplot2::after_stat(density)),
                              bins = bins, fill = "#4C72B0",
                              colour = "white", alpha = 0.8) +
      ggplot2::geom_density(colour = "#C44E52", linewidth = 1.1)
  } else {
    p <- p +
      ggplot2::geom_density(fill = "#4C72B0", colour = "#C44E52",
                            alpha = 0.5, linewidth = 1.0)
  }
  p +
    facet_wrap(~Trait, scales = "free") +
    labs(x = "Value", y = "Density") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}


#' Plot a phenotype trait-correlation heatmap
#'
#' Lower-triangular correlation heatmap with coefficient labels. Portable in the
#' same way as \code{\link{plot_phenotype_distribution}}.
#'
#' @param x      A \code{GenoMeta} with an attached phenotype, or a phenotype
#'   data.frame (\code{Sample} + traits).
#' @param traits Optional character vector selecting / ordering traits.
#' @param method Correlation method passed to \code{stats::cor}:
#'   \code{"pearson"} (default), \code{"spearman"}, or \code{"kendall"}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' plot_phenotype_correlation(meta, method = "spearman")
#' }
#'
#' @seealso \code{\link{plot_phenotype_distribution}}
#'
#' @export
plot_phenotype_correlation <- function(x, traits = NULL, method = "pearson") {
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  long <- .as_phenotype_table(x, traits = traits)
  wide <- attr(long, "wide")

  trait_cols <- setdiff(names(wide), "Sample")
  if (length(trait_cols) < 2L)
    stop("plot_phenotype_correlation: need at least two traits.")

  cmat <- stats::cor(as.matrix(wide[, trait_cols, drop = FALSE]),
                     method = method, use = "pairwise.complete.obs")

  lv <- trait_cols
  dt <- as.data.frame.table(cmat, stringsAsFactors = FALSE)
  names(dt) <- c("Var1", "Var2", "Corr")
  # Keep lower triangle (incl. diagonal) for a clean half-matrix.
  dt <- dt[match(dt$Var1, lv) >= match(dt$Var2, lv), , drop = FALSE]
  dt$Var1 <- factor(dt$Var1, levels = lv)
  dt$Var2 <- factor(dt$Var2, levels = rev(lv))

  ggplot(dt, aes(x = Var1, y = Var2, fill = Corr)) +
    geom_tile(colour = "white") +
    ggplot2::geom_text(aes(label = sprintf("%.2f", Corr)), size = 4) +
    scale_fill_gradient2(low = "#C44E52", mid = "white", high = "#4C72B0",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    ggplot2::coord_fixed() +
    labs(x = NULL, y = NULL,
         title = paste0(tools::toTitleCase(method), " correlation")) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
}


#' Per-trait numeric summary of a phenotype table
#'
#' Returns mean / sd / median / skew / min / max for each (selected) trait.
#' Portable in the same way as \code{\link{plot_phenotype_distribution}}.
#'
#' @param x      A \code{GenoMeta} with an attached phenotype, or a phenotype
#'   data.frame (\code{Sample} + traits).
#' @param traits Optional character vector selecting / ordering traits.
#' @param digits Rounding digits. Default \code{3}.
#'
#' @return A data.frame with one row per trait and columns \code{Trait},
#'   \code{Mean}, \code{SD}, \code{Median}, \code{Skew}, \code{Min}, \code{Max},
#'   \code{N}.
#'
#' @export
phenotype_summary <- function(x, traits = NULL, digits = 3L) {
  long <- .as_phenotype_table(x, traits = traits)
  wide <- attr(long, "wide")
  trait_cols <- setdiff(names(wide), "Sample")

  rows <- lapply(trait_cols, function(tc) {
    v  <- wide[[tc]]
    sd <- stats::sd(v, na.rm = TRUE)
    data.frame(
      Trait  = tc,
      Mean   = round(mean(v, na.rm = TRUE), digits),
      SD     = round(sd, digits),
      Median = round(stats::median(v, na.rm = TRUE), digits),
      Skew   = round((mean(v, na.rm = TRUE) - stats::median(v, na.rm = TRUE)) /
                       sd, digits),
      Min    = round(min(v, na.rm = TRUE), digits),
      Max    = round(max(v, na.rm = TRUE), digits),
      N      = sum(!is.na(v)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


#' Per-sample phenotype bar chart (up to two traits, dodged)
#'
#' Draws one (or two, side-by-side / dodged) bar(s) per sample. Portable: \code{x}
#' may be a \code{GenoMeta} (its \code{$phenotype} is used) or any data.frame with a
#' \code{Sample} column plus numeric trait columns.
#'
#' The function has two modes. In a \strong{direct} call use
#' \code{orientation = "vertical"} (the default): samples run along the X-axis and
#' trait values up the Y-axis. For a \strong{merged} call — placing the bars beside
#' an alignment heatmap that carries samples on its Y-axis (see
#' \code{\link{plot_meta_alignment}}) — use \code{orientation = "horizontal"}:
#' samples run down the Y-axis and values along the X-axis, so the sample rows line
#' up with the alignment rows.
#'
#' @param x        A \code{GenoMeta} with an attached phenotype, or a phenotype
#'   data.frame (\code{Sample} + traits).
#' @param traits   Optional character vector selecting / ordering traits. At most
#'   two are drawn; if more are available the first two are kept (with a warning).
#' @param samples  Optional character vector selecting / ordering the samples to
#'   show. \code{NULL} (default) keeps all samples.
#' @param orientation \code{"vertical"} (default, samples on X) or
#'   \code{"horizontal"} (samples on Y, merge-ready).
#' @param sample_levels Optional explicit factor levels for the sample axis, used
#'   to lock the order to a companion plot. \code{NULL} (default) uses
#'   \code{samples}, else \code{sample_order}, else the sorted unique sample IDs.
#'   Supplied levels win over \code{sample_order}: this is the explicit lock used
#'   when pinning this panel beside another.
#' @param sample_order Optional row order: a character vector, an
#'   \code{hs_sample_clustering} from \code{\link{cluster_samples}}, \code{TRUE}
#'   to require the broadcast order (\code{\link{use_sample_order}}), or
#'   \code{FALSE} to ignore it. \code{NULL} (default) uses the broadcast order
#'   when one is set. Applied to whichever samples remain after \code{samples}
#'   filtering; unmatched ones follow the ordered ones.
#' @param bar_width Width passed to \code{geom_col}. Default \code{0.7}.
#' @param palette  Fill colours for the (up to two) traits. Default a blue/red pair.
#' @param legend_position Passed to \code{theme(legend.position=)}. Default
#'   \code{"right"}; \code{plot_meta_alignment()} overrides this when collecting
#'   legends.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' plot_phenotype_bars(meta, traits = c("Pheno1", "Pheno2"))            # vertical
#' plot_phenotype_bars(meta, orientation = "horizontal")               # merge-ready
#' }
#'
#' @seealso \code{\link{plot_meta_alignment}}, \code{\link{plot_phenotype_distribution}}
#'
#' @export
plot_phenotype_bars <- function(x, traits = NULL, samples = NULL,
                                orientation = c("vertical", "horizontal"),
                                sample_levels = NULL, sample_order = NULL,
                                bar_width = 0.7,
                                palette = c("#4C72B0", "#C44E52"),
                                legend_position = "right") {
  orientation <- match.arg(orientation)
  long <- .as_phenotype_table(x, traits = traits)
  trait_cols <- levels(long$Trait)

  # "max col 2": keep at most the first two traits.
  if (length(trait_cols) > 2L) {
    warning("plot_phenotype_bars: more than two traits supplied; keeping the ",
            "first two (", paste(trait_cols[1:2], collapse = ", "), ").")
    trait_cols <- trait_cols[1:2]
    long <- long[long$Trait %in% trait_cols, , drop = FALSE]
    long$Trait <- factor(as.character(long$Trait), levels = trait_cols)
  }

  # Optional sample subset.
  if (!is.null(samples)) {
    samples <- as.character(samples)
    missing_s <- setdiff(samples, unique(long$Sample))
    if (length(missing_s) > 0L)
      warning("plot_phenotype_bars: ", length(missing_s),
              " requested sample(s) not present: ",
              paste(utils::head(missing_s, 5L), collapse = ", "),
              if (length(missing_s) > 5L) " ..." else "")
    long <- long[long$Sample %in% samples, , drop = FALSE]
  }

  # Resolve the sample-axis factor order. `sample_levels` is the explicit lock
  # (it exists to pin this panel to a companion figure) and therefore outranks
  # both `sample_order` and the broadcast; below it the shared contract applies.
  order_res <- NULL
  if (is.null(sample_levels)) {
    base_levels <- if (!is.null(samples)) samples else sort(unique(long$Sample))
    order_res <- .hs_resolve_sample_order(
      sample_order, base_levels, caller = "plot_phenotype_bars",
      subset = "intersect", default = "keep")
    sample_levels <- order_res$order %||% base_levels
  }
  sample_levels <- as.character(sample_levels)
  long$Sample <- factor(long$Sample, levels = sample_levels)

  dodge <- ggplot2::position_dodge2(preserve = "single", padding = 0.1)

  if (orientation == "horizontal") {
    p <- ggplot(long, aes(x = Value, y = Sample, fill = Trait)) +
      ggplot2::geom_col(position = dodge, width = bar_width) +
      labs(x = "Phenotype value", y = "Sample", fill = "Trait")
  } else {
    p <- ggplot(long, aes(x = Sample, y = Value, fill = Trait)) +
      ggplot2::geom_col(position = dodge, width = bar_width) +
      labs(x = "Sample", y = "Phenotype value", fill = "Trait")
  }

  p <- p +
    scale_fill_manual(values = palette, drop = FALSE) +
    theme_minimal(base_size = 12) +
    theme(legend.position = legend_position,
          panel.grid.major.y = if (orientation == "horizontal")
            element_blank() else element_line(),
          axis.text.x = if (orientation == "vertical")
            element_text(angle = 45, hjust = 1) else element_text())

  .hs_tag_sample_order(p, order_res)
}


# Internal: read a ggplot's discrete axis labels in plot order (bottom->top for y,
# left->right for x), via ggplot_build introspection. Returns character(0) when the
# axis is not discrete or cannot be read (callers then skip reordering).
.discrete_axis_levels <- function(p, which = c("y", "x")) {
  which <- match.arg(which)
  lv <- tryCatch({
    pp   <- ggplot2::ggplot_build(p)$layout$panel_params[[1L]]
    labs <- pp[[which]]$get_labels()
    as.character(labs[!is.na(labs)])
  }, error = function(e) character(0L))
  if (is.null(lv)) character(0L) else lv
}


#' Merge plots that share a discrete sample axis into one figure
#'
#' Combines several \code{ggplot} panels that share a discrete sample axis — for
#' example a genotype-level local-GEBV alignment (\code{\link{plot_localgebv_alignment}})
#' or a genotype alignment (\code{\link{plot_haplotype_alignment}}) on the left and a
#' phenotype bar chart (\code{\link{plot_phenotype_bars}}, \code{orientation =
#' "horizontal"}) on the right — into a single figure with \strong{identical sample
#' ticks} across panels and the panel legends collected into one shared row.
#'
#' Tick sharing is achieved by reading the sample order from the first
#' (\emph{reference}) plot and reordering every other panel's sample axis to match.
#' This requires the panels to use the \strong{same raw sample IDs} on that axis,
#' i.e. one row per sample (genotype-level alignments). Haplotype-level alignments
#' (two rows per sample, labelled \code{"<Sample>_hapN"}) are not supported here.
#'
#' Requires the \pkg{cowplot} package.
#'
#' @param plots     A list of \code{ggplot} objects ordered left-to-right;
#'   \code{plots[[1]]} is the reference whose sample axis defines the shared order.
#' @param rel_widths Optional numeric vector of relative panel widths (length
#'   \code{length(plots)}). \code{NULL} (default) gives equal widths.
#' @param sync_axis Discrete axis to synchronise across panels: \code{"y"}
#'   (default), \code{"x"}, or \code{"none"} to skip reordering.
#' @param legend    Where to place the collected legends: \code{"bottom"}
#'   (default, one shared row beneath the panels), \code{"right"} (each panel keeps
#'   its own legend), or \code{"none"}.
#' @param align     Passed to \code{cowplot::plot_grid} (\code{"h"} by default, so
#'   panel rows line up).
#' @param labels_on Index of the single panel that keeps its sample-axis tick
#'   labels/title (default \code{1}, the leftmost); the rest are blanked to avoid
#'   duplicate ticks.
#' @param check_order Logical; warn when the panels were drawn with
#'   \emph{different} sample orders (default \code{TRUE}). The reference panel's
#'   order silently wins the sync, so a clustered panel composed with an
#'   alphabetical one would be re-ordered without saying so — the failure the
#'   \code{sample_order} contract exists to prevent. Panels that carry no order
#'   record (drawn by a function outside the contract, or with the default order)
#'   are skipped.
#'
#' @return A combined plot object from \code{cowplot} (saveable with
#'   \code{ggplot2::ggsave}).
#'
#' @examples
#' \dontrun{
#' lg    <- compute_profile_localgebv(meta, eff, level = "genotype")
#' fam   <- plot_localgebv_alignment(lg, meta)[["gen_1"]][["F_1"]]
#' fams  <- .discrete_axis_levels(fam, "y")               # family sample order
#' bars  <- plot_phenotype_bars(meta, samples = fams, orientation = "horizontal")
#' plot_meta_alignment(list(fam, bars), rel_widths = c(3, 1))
#' }
#'
#' @seealso \code{\link{plot_phenotype_bars}}, \code{\link{plot_localgebv_alignment}},
#'   \code{\link{plot_haplotype_alignment}}
#'
#' @export
plot_meta_alignment <- function(plots, rel_widths = NULL, sync_axis = "y",
                                legend = c("bottom", "right", "none"),
                                align = "h", labels_on = 1L,
                                check_order = TRUE) {
  legend    <- match.arg(legend)
  sync_axis <- match.arg(sync_axis, c("y", "x", "none"))

  if (!requireNamespace("cowplot", quietly = TRUE))
    stop("plot_meta_alignment: the 'cowplot' package is required. ",
         "Install it with install.packages('cowplot').")
  if (!is.list(plots) || length(plots) < 2L)
    stop("plot_meta_alignment: 'plots' must be a list of at least two ggplots.")
  is_gg <- vapply(plots, function(p) inherits(p, "ggplot"), logical(1L))
  if (!all(is_gg))
    stop("plot_meta_alignment: every element of 'plots' must be a ggplot.")
  if (!is.null(rel_widths) && length(rel_widths) != length(plots))
    stop("plot_meta_alignment: 'rel_widths' must have one value per plot.")

  # 0. Do the panels agree on the row order? The sync below imposes the
  # reference's order on everyone, so a disagreement is silently resolved unless
  # it is reported here.
  if (isTRUE(check_order) && sync_axis != "none")
    .hs_warn_order_mismatch(plots)

  # 1. Synchronise the shared discrete axis to the reference's order.
  if (sync_axis != "none") {
    ref_levels <- .discrete_axis_levels(plots[[1L]], sync_axis)
    if (length(ref_levels) > 0L) {
      scale_fn <- if (sync_axis == "y") ggplot2::scale_y_discrete
                  else ggplot2::scale_x_discrete
      for (i in seq_along(plots)[-1L])
        plots[[i]] <- plots[[i]] + scale_fn(limits = ref_levels)
    }
  }

  # 2. Blank the shared-axis labels on every panel except `labels_on`.
  blank_axis <- if (sync_axis == "y")
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  else
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  for (i in seq_along(plots))
    if (i != labels_on && sync_axis != "none")
      plots[[i]] <- plots[[i]] + blank_axis

  # 3. Collect legends (when pooling at the bottom) and strip them from panels.
  legends <- NULL
  if (legend == "bottom") {
    legends <- lapply(plots, function(p)
      tryCatch(cowplot::get_legend(
        p + theme(legend.position = "bottom", legend.box = "horizontal")),
        error = function(e) NULL))
    legends <- Filter(Negate(is.null), legends)
    plots   <- lapply(plots, function(p) p + theme(legend.position = "none"))
  }

  # 4. Arrange the panels in one row, aligned so the sample rows match.
  panel_row <- cowplot::plot_grid(plotlist = plots, nrow = 1L, align = align,
                                  axis = "tb", rel_widths = rel_widths)

  if (legend != "bottom" || length(legends) == 0L)
    return(panel_row)

  legend_row <- cowplot::plot_grid(plotlist = legends, nrow = 1L)
  cowplot::plot_grid(panel_row, legend_row, ncol = 1L,
                     rel_heights = c(1, 0.12))
}
