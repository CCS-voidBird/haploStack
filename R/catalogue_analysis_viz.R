# =============================================================================
# catalogue_analysis_viz.R
#
# Plot layer for the cross-trait block-catalogue analysis. Each plotter consumes
# the list returned by analyze_block_catalogue() and draws one figure from the one
# subtype selected per block (shared across traits), filling by that subtype's
# per-trait value — raw effect S or within-trait standard score N, chosen by
# `value_source` (defaulting to the source stored on the analysis result, D2):
#
#   plot_block_effect_genome()    — rectangular genome (block extent x trait,
#                                   geom_rect filled by value_source)
#   plot_block_effect_manhattan() — multi-trait Manhattan, top blocks labelled
#   plot_block_effect_density()   — density of selected-subtype values per trait
#   plot_top_block_heatmap()      — sample x top-block heatmap, faceted by trait
#   plot_block_effect_circos()    — circular genome, one ring per trait (falls
#                                   back to the linear genome plot without circlize)
# =============================================================================


# Internal: validate the analysis list carries the named tables a plotter needs.
.check_analysis <- function(analysis, need, fn) {
  if (!is.list(analysis) || !all(need %in% names(analysis)))
    stop(fn, ": 'analysis' must be the list returned by ",
         "analyze_block_catalogue() (missing: ",
         paste(setdiff(need, names(analysis)), collapse = ", "), ").")
  invisible(TRUE)
}

# Internal: resolve the plot value source. A NULL default falls back to the
# `value_source` stored on the analysis result (D2), so the built plots honour the
# analyze-chosen source without the caller re-specifying it; an explicit
# "raw"/"normalized" overrides it for a re-plot.
#
# `rescaled` reports whether N is actually a transform of S: since 0.5.4 the
# analysis default is normalize = "none", so a result's N column can be the raw
# values. A caller overriding value_source = "normalized" on such a result must
# not get a subtitle claiming a standard score. Older results carry no
# `normalize` field; they predate the raw default, so their N is a transform.
#
# Returns list(source, col, table, rescaled).
.resolve_value_source <- function(analysis, value_source) {
  if (is.null(value_source))
    value_source <- if (!is.null(analysis$value_source)) analysis$value_source
                    else "raw"
  value_source <- match.arg(value_source, c("raw", "normalized"))
  norm <- analysis$normalize
  list(source   = value_source,
       col      = if (value_source == "normalized") "N" else "S",
       table    = if (value_source == "normalized") analysis$N else analysis$S,
       rescaled = is.null(norm) || !identical(norm, "none"))
}

# Internal: short, readable subtype label for plot annotation.
.subtype_label <- function(x, n = 6L) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1L, n), "…"), x)
}


#' Rectangular genome plot of per-block cross-trait effects
#'
#' Draws each block as a rectangle spanning its \strong{genomic extent}
#' (\code{Block_start} to \code{Block_end}, faceted by chromosome) coloured by its
#' \strong{selected subtype}'s per-trait value, with one band per trait. Because
#' selection is per block, a block's subtype is identical across its trait bands.
#' The block extents are the origin \code{Block_start}/\code{Block_end} carried on
#' the analysis value tables — derived from the GenoMeta block locations — so each
#' rectangle covers the block's real span rather than a single point.
#'
#' The diverging fill is anchored at \code{0} (white). Under
#' \code{value_source = "raw"} white means a \strong{zero-effect} block; under
#' \code{"normalized"} (a within-trait standard score) white means an
#' \strong{average-effect} block for that trait.
#'
#' @param analysis The list returned by \code{\link{analyze_block_catalogue}}.
#' @param value_source Which per-(block, trait) value to colour by: \code{"raw"}
#'   (the selected subtype's raw effect \code{S}) or \code{"normalized"} (its
#'   within-trait standard score \code{N}). \code{NULL} (default) uses the
#'   \code{value_source} stored on \code{analysis}, so the plot matches the
#'   analyze-time choice; pass one explicitly to re-plot the other source without
#'   recomputing.
#' @param title Optional plot title.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{analyze_block_catalogue}}
#' @export
plot_block_effect_genome <- function(analysis, value_source = NULL,
                                     title = NULL) {
  .check_analysis(analysis, c("S", "N"), "plot_block_effect_genome")
  vs <- .resolve_value_source(analysis, value_source)
  dt <- as.data.table(vs$table)
  dt[, PlotVal := get(vs$col)]

  # Block extents (Block_start, Block_end) are carried on the value table — the
  # origin block locations sourced from the GenoMeta. Fall back to the block
  # midpoint Pos when an extent is missing (e.g. singleton / dummy blocks).
  if (!all(c("Block_start", "Block_end") %in% names(dt)))
    dt[, `:=`(Block_start = NA_real_, Block_end = NA_real_)]
  dt[is.na(Block_start), Block_start := Pos]
  dt[is.na(Block_end),   Block_end   := Pos]

  # Give zero-width (point) blocks a small visible width so they still render.
  span <- suppressWarnings(
    diff(range(c(dt$Block_start, dt$Block_end), na.rm = TRUE)))
  if (!is.finite(span) || span == 0) span <- 1
  half <- span * 0.002
  dt[Block_end <= Block_start,
     `:=`(Block_start = Block_start - half, Block_end = Block_end + half)]

  # Map traits to integer y-bands for geom_rect (rectangles need ymin/ymax).
  dt[, Trait := factor(Trait)]
  trait_levels <- levels(dt$Trait)
  dt[, yc := as.integer(Trait)]

  # One chromosome per row (vertical panel stack); strip label in the right
  # margin so each chromosome reads as its own track.
  n_chr <- data.table::uniqueN(dt$Chrom)

  # Diverging fill anchored at zero (white): the scale spans the real data range
  # across zero, (-a, 0, b), so negative blocks read blue and positive red while
  # zero stays white regardless of any asymmetry between a and b.
  rng <- suppressWarnings(range(dt$PlotVal, na.rm = TRUE))
  if (!all(is.finite(rng))) rng <- c(-1, 1)
  fill_lim <- c(min(rng[1], 0), max(rng[2], 0))

  ggplot(dt, aes(xmin = Block_start, xmax = Block_end,
                 ymin = yc - 0.45, ymax = yc + 0.45, fill = PlotVal)) +
    geom_rect(colour = NA) +
    facet_wrap(~ Chrom, nrow = n_chr,) +
    scale_y_continuous(breaks = seq_along(trait_levels), labels = trait_levels,
                       expand = c(0, 0.6)) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = fill_lim,
                         breaks = c(fill_lim[1], 0, fill_lim[2])) +
    labs(title = title %||% paste0("Block effect across traits (", vs$source, ")"),
         subtitle = if (vs$source == "normalized" && vs$rescaled)
           "white = average effect (within-trait standard score)"
         else "white = zero effect (raw)",
         x = "Genomic position", y = "Trait", fill = vs$col) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 6))
}


#' Multi-trait Manhattan of selected-subtype effects
#'
#' One Manhattan panel per trait (shared cumulative genomic axis), each block
#' drawn from its selected subtype. Per-trait top blocks are enlarged and
#' labelled with the subtype id and its population \code{Count}.
#'
#' @param analysis The list returned by \code{\link{analyze_block_catalogue}}.
#' @param value_source Which value to plot on y: \code{"raw"} (the selected
#'   subtype's raw effect \code{S}) or \code{"normalized"} (its within-trait
#'   standard score \code{N}). \code{NULL} (default) uses the \code{value_source}
#'   stored on \code{analysis}.
#' @param label_top Logical; label each trait's top block(s). Default \code{TRUE}.
#' @param point_size Base point size. Default \code{1.4}.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{analyze_block_catalogue}},
#'   \code{\link{plot_manhattan_effect}}
#' @export
plot_block_effect_manhattan <- function(analysis, value_source = NULL,
                                        label_top = TRUE, point_size = 1.4) {
  .check_analysis(analysis, c("S", "N"), "plot_block_effect_manhattan")
  vs <- .resolve_value_source(analysis, value_source)
  dt <- as.data.table(vs$table)
  dt[, PlotVal := get(vs$col)]

  pos_info <- .create_genomic_positions(dt, "Chrom", "Pos")
  pdt      <- pos_info$data
  chr_info <- pos_info$chr_info

  p <- ggplot(pdt, aes(x = GenomicPos, y = PlotVal, colour = factor(Chrom))) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(size = point_size, alpha = 0.85) +
    facet_wrap(~ Trait, ncol = 1L, scales = "free_y") +
    scale_x_continuous(breaks = chr_info$mid, labels = chr_info$Chrom,
                       expand = c(0.01, 0)) +
    labs(title = "Cross-trait block Manhattan", x = "Chromosome",
         y = vs$col, colour = "Chr") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  if (isTRUE(label_top) && "is_top" %in% names(pdt)) {
    tops <- pdt[is_top == TRUE]
    if (nrow(tops) > 0L) {
      tops[, lab := paste0(.subtype_label(Subtype), " (", Count, ")")]
      p <- p +
        geom_point(data = tops, size = point_size + 1.6, colour = "black") +
        geom_text(data = tops, aes(label = lab), vjust = -0.8, size = 2.6,
                  colour = "black")
    }
  }
  p
}


#' Density of selected-subtype effects across traits
#'
#' @param analysis The list returned by \code{\link{analyze_block_catalogue}}.
#' @param value_source Which value to show: \code{"raw"} (the selected subtype's
#'   raw effect \code{S}) or \code{"normalized"} (its within-trait standard score
#'   \code{N}). \code{NULL} (default) uses the \code{value_source} stored on
#'   \code{analysis}.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{analyze_block_catalogue}}
#' @export
plot_block_effect_density <- function(analysis, value_source = NULL) {
  .check_analysis(analysis, c("S", "N"), "plot_block_effect_density")
  vs <- .resolve_value_source(analysis, value_source)
  dt <- as.data.table(vs$table)
  dt[, PlotVal := get(vs$col)]
  dt[, Trait := factor(Trait)]

  ggplot(dt, aes(x = PlotVal, colour = Trait, fill = Trait)) +
    geom_density(alpha = 0.25, na.rm = TRUE) +
    labs(title = "Distribution of selected-subtype effects",
         x = vs$col, y = "Density") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


#' Sample x top-block heatmap, faceted by trait
#'
#' For each trait, takes its \code{top_n} top blocks (by selected effect) and
#' draws a sample x block heatmap of the genotype-level local GEBVs, optionally
#' z-scored within each block.
#'
#' @param analysis The list returned by \code{\link{analyze_block_catalogue}}
#'   (must carry \code{gebv}).
#' @param top_n Number of top blocks per trait. Default \code{1}.
#' @param scale Logical; z-score each block's GEBVs across samples. Default
#'   \code{TRUE}.
#' @param sample_order Optional order for the sample (Y) axis: a character
#'   vector, an \code{hs_sample_clustering} from \code{\link{cluster_samples}},
#'   \code{TRUE} to require the broadcast order (\code{\link{use_sample_order}}),
#'   or \code{FALSE} to ignore it. \code{NULL} (default) uses the broadcast order
#'   when one is set, else the alphabetical order \pkg{ggplot2} gives a discrete
#'   axis. Rows read bottom-to-top, so the first sample in the order is the
#'   \emph{top} row — matching the alignment figures.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{analyze_block_catalogue}}, \code{\link{use_sample_order}}
#' @export
plot_top_block_heatmap <- function(analysis, top_n = 1L, scale = TRUE,
                                   sample_order = NULL) {
  .check_analysis(analysis, c("S", "gebv"), "plot_top_block_heatmap")
  bm   <- as.data.table(analysis$S)
  gebv <- analysis$gebv

  parts <- lapply(unique(bm$Trait), function(tr) {
    gm <- gebv[[tr]]
    if (is.null(gm)) {
      message("plot_top_block_heatmap: no GEBV matrix for trait '", tr,
              "'; skipped.")
      return(NULL)
    }
    gm   <- as.data.table(gm)
    blks <- bm[Trait == tr][order(-S)][seq_len(min(top_n, .N)), Block]
    sub  <- gm[Block %in% blks]
    if (nrow(sub) == 0L) return(NULL)
    scols <- setdiff(names(sub), "Block")
    long  <- data.table::melt(sub, id.vars = "Block", measure.vars = scols,
                              variable.name = "Sample", value.name = "GEBV",
                              variable.factor = FALSE)
    if (isTRUE(scale))
      long[, GEBV := {
        s <- stats::sd(GEBV, na.rm = TRUE)
        if (is.na(s) || s == 0) GEBV * 0 else (GEBV - mean(GEBV, na.rm = TRUE)) / s
      }, by = Block]
    long[, Trait := tr]
    long
  })
  dt <- data.table::rbindlist(Filter(Negate(is.null), parts))
  if (nrow(dt) == 0L)
    stop("plot_top_block_heatmap: no GEBV matrices available to plot.")

  # Row order. Without one the Y axis is ggplot2's alphabetical discrete scale,
  # which is what this figure always had; `rev()` because a discrete Y axis draws
  # its first level at the BOTTOM and the alignment figures put the first sample
  # at the top.
  order_res <- .hs_resolve_sample_order(
    sample_order, unique(dt$Sample), caller = "plot_top_block_heatmap",
    subset = "intersect", default = "keep")
  if (!is.null(order_res$order))
    dt[, Sample := factor(Sample, levels = rev(order_res$order))]

  p <- ggplot(dt, aes(x = Block, y = Sample, fill = GEBV)) +
    geom_tile() +
    facet_wrap(~ Trait, scales = "free_x") +
    scale_fill_gradient2(low = "#2166ac", mid = "grey92", high = "#b2182b",
                         midpoint = 0) +
    labs(title = paste0("Top-", top_n, " block GEBVs per trait",
                        if (scale) " (z-scored)" else ""),
         x = "Block", y = "Sample",
         fill = if (scale) "z(GEBV)" else "GEBV") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 7),
          panel.grid = element_blank())

  .hs_tag_sample_order(p, order_res)
}


#' Circular (circos) plot of per-block cross-trait effects
#'
#' Draws the genome as chromosome arcs with \strong{one ring per trait}; each
#' block is a segment spanning its genomic range, coloured by its selected
#' subtype's per-trait value (\code{value_source}). Per-trait top blocks are labelled,
#' and (optionally) the traits' top blocks are chord-linked to highlight shared /
#' cross-trait signal.
#'
#' Requires the \pkg{circlize} package (in \code{Suggests}). When it is not
#' installed the function emits a \code{message()} and returns the linear
#' equivalent (\code{\link{plot_block_effect_genome}}) instead — there is no
#' in-package polar approximation.
#'
#' @param analysis The list returned by \code{\link{analyze_block_catalogue}}.
#' @param value_source Which per-(block, trait) value to colour by: \code{"raw"}
#'   (the selected subtype's raw effect \code{S}) or \code{"normalized"} (its
#'   within-trait standard score \code{N}). \code{NULL} (default) uses the
#'   \code{value_source} stored on \code{analysis}. The diverging scale is centred
#'   at \code{0} (grey): under \code{"raw"} that means zero effect, under
#'   \code{"normalized"} it means average effect for the trait.
#' @param link_top Logical; chord-link the traits' top blocks. Default
#'   \code{FALSE}.
#' @param palette Length-3 colour vector (low, mid, high) for the diverging
#'   value scale, centred at 0. Default blue–grey–red.
#'
#' @return When \pkg{circlize} is available, invisibly a recorded base plot
#'   (\code{grDevices::recordPlot()}); otherwise the \code{ggplot} fallback.
#' @seealso \code{\link{analyze_block_catalogue}},
#'   \code{\link{plot_block_effect_genome}}
#' @export
plot_block_effect_circos <- function(analysis, value_source = NULL,
                                     link_top = FALSE,
                                     palette = c("#2166ac", "grey92", "#b2182b")) {
  .check_analysis(analysis, c("S", "N"), "plot_block_effect_circos")
  vs <- .resolve_value_source(analysis, value_source)

  if (!requireNamespace("circlize", quietly = TRUE)) {
    message("plot_block_effect_circos: package 'circlize' is not installed; ",
            "returning the linear equivalent (plot_block_effect_genome()).")
    return(plot_block_effect_genome(analysis, value_source = vs$source))
  }

  dt <- as.data.table(vs$table)
  dt[, MetricVal := get(vs$col)]
  if (!all(c("Block_start", "Block_end") %in% names(dt)))
    dt[, `:=`(Block_start = NA_real_, Block_end = NA_real_)]
  dt[is.na(Block_start), Block_start := Pos]
  dt[is.na(Block_end),   Block_end   := Pos]

  # Chromosome sectors (natural order) and their genomic spans.
  chr_levels <- {
    u <- unique(as.character(dt$Chrom))
    n <- suppressWarnings(as.integer(u))
    u[order(!is.na(n), n, u, na.last = TRUE)]
  }
  xlim <- do.call(rbind, lapply(chr_levels, function(ch) {
    z <- dt[Chrom == ch]
    range(c(z$Block_start, z$Block_end), na.rm = TRUE)
  }))

  amax    <- max(abs(dt$MetricVal), na.rm = TRUE)
  if (!is.finite(amax) || amax == 0) amax <- 1
  col_fun <- circlize::colorRamp2(c(-amax, 0, amax), palette)
  traits  <- unique(as.character(dt$Trait))

  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  circlize::circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 3,
                       start.degree = 90)
  circlize::circos.initialize(sectors = chr_levels, xlim = xlim)

  # Chromosome name axis (outermost, thin).
  circlize::circos.track(ylim = c(0, 1), track.height = 0.06, bg.border = NA,
    panel.fun = function(x, y) {
      circlize::circos.text(circlize::CELL_META$xcenter, 0.5,
                            circlize::CELL_META$sector.index,
                            cex = 0.7, facing = "inside", niceFacing = TRUE)
    })

  # One ring per trait: each block a coloured segment; top blocks labelled.
  for (tr in traits) local({
    tdat <- dt[Trait == tr]
    circlize::circos.track(ylim = c(0, 1), track.height = 0.12, bg.border = NA,
      panel.fun = function(x, y) {
        sd <- tdat[Chrom == circlize::CELL_META$sector.index]
        if (nrow(sd) > 0L) {
          circlize::circos.rect(sd$Block_start, rep(0, nrow(sd)),
                                sd$Block_end, rep(1, nrow(sd)),
                                col = col_fun(sd$MetricVal), border = NA)
          tp <- sd[is_top == TRUE]
          if (nrow(tp) > 0L)
            circlize::circos.text((tp$Block_start + tp$Block_end) / 2,
                                  0.5, .subtype_label(tp$Subtype),
                                  cex = 0.5, col = "black")
        }
      })
  })

  # Optional chord links between the traits' top blocks.
  if (isTRUE(link_top) && length(traits) >= 2L) {
    tops <- dt[is_top == TRUE, .SD[1L], by = Trait]
    if (nrow(tops) >= 2L)
      for (i in seq_len(nrow(tops) - 1L))
        circlize::circos.link(
          tops$Chrom[i],     (tops$Block_start[i]     + tops$Block_end[i])     / 2,
          tops$Chrom[i + 1L], (tops$Block_start[i + 1L] + tops$Block_end[i + 1L]) / 2,
          col = "#00000040", lwd = 1.5)
  }

  invisible(grDevices::recordPlot())
}


# ---------------------------------------------------------------------------
# Internal: build the requested plot set for analyze_block_catalogue()
# ---------------------------------------------------------------------------

# Returns a named list of plot objects for the requested `plots`. Each plotter is
# wrapped so a single failure degrades to a warning rather than aborting the whole
# analysis. Every plotter is passed the analysis-chosen `value_source` (read off
# the result, D2) so the built figures fill by the same value the selection used.
# "circos" builds via plot_block_effect_circos() (falling back to the linear
# genome plot with a message when 'circlize' is not installed).
.build_analysis_plots <- function(result, plots, multi) {
  if (length(plots) == 0L) return(list())
  plots <- unique(plots)
  out   <- list()
  vsrc  <- result$value_source

  safe <- function(nm, expr) {
    tryCatch(out[[nm]] <<- expr,
             error = function(e)
               warning("analyze_block_catalogue: plot '", nm, "' skipped: ",
                       conditionMessage(e)))
  }

  if ("rectangular" %in% plots) safe("rectangular", plot_block_effect_genome(result, value_source = vsrc))
  if ("manhattan"   %in% plots) safe("manhattan",   plot_block_effect_manhattan(result, value_source = vsrc))
  if ("density"     %in% plots) safe("density",     plot_block_effect_density(result, value_source = vsrc))
  if ("heatmap"     %in% plots) safe("heatmap",     plot_top_block_heatmap(result))
  if ("circos"      %in% plots) safe("circos",      plot_block_effect_circos(result, value_source = vsrc))

  out
}
