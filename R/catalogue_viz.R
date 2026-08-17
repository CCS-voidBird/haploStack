# =============================================================================
# catalogue_viz.R
#
# Visualisation functions for the haplotype catalogue module.
#
# Functions:
#   plot_haplotype_alignment()  — alignment heatmap of a haplotype profile
#   plot_localgebv_alignment()  — per-family local-GEBV heatmaps (nested list)
#   plot_manhattan_effect()     — Manhattan plot of block/marker effects
#   plot_manhattan_variance()   — Manhattan plot of block effect variance
#   plot_block_counts()         — stacked bar chart of positive/negative blocks
#   plot_marker_signal()        — Manhattan of marker signal + block strip
#   save_marker_signal()        — write that plot + its marker x block table
#   .create_genomic_positions() — internal: cumulative genomic X positions
#   .natural_chrom_levels()     — internal: numeric-then-alphabetic chrom order
#   .resolve_signal_markers()   — internal: selection argument -> marker IDs
#   .marker_signal_table()      — internal: marker x block table behind the plot
# =============================================================================

#' Plot a haplotype alignment-style heatmap from a GenoMeta object
#'
#' Generates a rectangular heatmap where each cell represents one marker allele
#' of one sample-haplotype, coloured by allele value. Markers run along the
#' X-axis (grouped within blocks); sample-haplotype pairs run along the Y-axis.
#'
#' @section Chromosomes and blocks:
#' The figure carries two column-level groupings. **Chromosomes** are drawn as
#' sub-panels — one facet column per chromosome, each as wide as its marker
#' count, named by the strip above it (\code{chrom_panels = TRUE}, the default;
#' a single-chromosome panel is never faceted and keeps the \code{Chr …}
#' subtitle instead). **Blocks** sit inside each panel: dot-dash boundary lines
#' between adjacent blocks, and vertical block IDs hanging from the top of the
#' panel into headroom reserved on the Y scale, so they clear both the tiles and
#' the chromosome strip. Set \code{chrom_panels = FALSE} for the pre-0.5.0
#' single continuous marker axis.
#'
#' @param geno_meta A valid \code{GenoMeta} list from
#'   \code{\link{new_geno_meta}} or any \code{read_*} function.
#' @param color_map Optional named character vector mapping allele values to hex
#'   colours. When \code{NULL} (default), a palette from \code{COLOR_MAP} is
#'   chosen based on the genotype encoding.
#' @param method Row resolution. \code{"haplotype"} (default) shows each
#'   haplotype copy on its own row, labelled \code{"<Sample>_hap<idx>"};
#'   \code{"genotype"} collapses the two haplotypes of a diploid sample into one
#'   sorted two-allele genotype per marker, on a single row labelled
#'   \code{"<Sample>_g"} (diploid only).
#' @param sample_order Optional ordering of the sample rows: either a character
#'   vector of sample IDs (a permutation of the panel's samples) or an
#'   \code{hs_sample_clustering} object from \code{\link{cluster_samples}} (its
#'   leaf \code{order} is used). \code{NULL} (default) keeps the alphabetical
#'   order. In \code{"haplotype"} mode the haplotype copies stay grouped within
#'   each sample, in this order.
#' @param show_block_boundaries Logical; draw vertical lines at block
#'   boundaries. Default \code{TRUE}.
#' @param show_block_labels Logical; label block IDs at the top of the plot.
#'   Default \code{TRUE}.
#' @param show_individual_borders Logical; draw horizontal lines between
#'   individual samples. Default \code{FALSE}.
#' @param show_individual_names Logical; label sample names on the Y-axis.
#'   Default \code{FALSE}.
#' @param show_marker_labels Logical; label marker IDs on the X-axis.
#'   Default \code{FALSE}.
#' @param chrom_panels Logical; split the alignment into one sub-panel (facet
#'   column) per chromosome, each panel as wide as its share of the markers.
#'   Default \code{TRUE}. Automatically a no-op when the data hold a single
#'   chromosome, which is then drawn as one panel with a \code{Chr …} subtitle.
#' @param show_chrom_labels Logical; draw the chromosome name in the strip above
#'   each sub-panel. Default \code{TRUE}. Ignored when \code{chrom_panels} is
#'   \code{FALSE} or only one chromosome is present.
#' @param block_label_size Font size of the vertical block-ID labels. Default
#'   \code{2.5}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' # Whole genome: one sub-panel per chromosome, blocks labelled within each.
#' p <- plot_haplotype_alignment(meta, show_block_labels = TRUE)
#' ggplot2::ggsave("alignment.png", p, width = 14, height = 8)
#'
#' # Pre-0.5.0 look: one continuous marker axis, chromosomes only in the subtitle.
#' plot_haplotype_alignment(meta, chrom_panels = FALSE)
#' }
#'
#' @seealso \code{\link{plot_clustered_alignment}},
#'   \code{\link{plot_localgebv_alignment}}
#'
# Resolve the `sample_order` argument of plot_haplotype_alignment() to a validated
# permutation of the panel's samples. Since 0.5.5 this is a thin wrapper over the
# shared contract in R/plot_order.R (which added the broadcast option, the
# provenance record and the intersect mode); the strict/alphabetical behaviour it
# had before is preserved exactly by `subset = "strict"`, `default = "sorted"`.
# @keywords internal
.resolve_sample_order <- function(sample_order, samples,
                                  caller = "plot_haplotype_alignment") {
  .hs_resolve_sample_order(sample_order, samples, caller = caller,
                           subset = "strict", default = "sorted")
}


# Chromosome names in natural order: those starting with an integer sort
# numerically first ("1H" < "2H" < "10H"), the rest alphabetically after. Plain
# alphabetical order — what order(chromosome, position) used before 0.5.0 —
# puts "10H" ahead of "2H", which then reorders both the marker axis and the
# chromosome sub-panels. Same intent as the chromosome ordering in
# .create_genomic_positions(), but keyed on the leading digits so barley-style
# "1H".."7H" labels sort as well as bare "1".."7".
# @keywords internal
.natural_chrom_levels <- function(x) {
  u <- unique(as.character(x))
  u <- u[!is.na(u)]
  if (length(u) == 0L) return(character(0L))
  num <- suppressWarnings(as.numeric(
    ifelse(grepl("^[0-9]+", u), sub("^([0-9]+).*$", "\\1", u), NA)))
  u[order(is.na(num), num, u)]
}


#' @export
plot_haplotype_alignment <- function(geno_meta,
                                      color_map               = NULL,
                                      method                  = "haplotype",
                                      sample_order            = NULL,
                                      show_block_boundaries   = TRUE,
                                      show_block_labels       = TRUE,
                                      show_individual_borders = FALSE,
                                      show_individual_names   = FALSE,
                                      show_marker_labels      = FALSE,
                                      chrom_panels            = TRUE,
                                      show_chrom_labels       = TRUE,
                                      block_label_size        = 2.5) {

  if (!is_geno_meta(geno_meta))
    stop("plot_haplotype_alignment: 'geno_meta' is not a valid GenoMeta.")
  method <- match.arg(method, c("haplotype", "genotype"))

  # Read the genotype array directly — no profile materialisation.
  arr <- geno_meta$geno_arrays
  if (dim(arr)[1L] == 0L || dim(arr)[2L] == 0L)
    stop("plot_haplotype_alignment: GenoMeta contains no data to plot.")
  samples <- dimnames(arr)[[1L]]
  markers <- dimnames(arr)[[2L]]
  ploidy  <- dim(arr)[3L]

  if (method == "genotype" && ploidy != 2L)
    stop("plot_haplotype_alignment: method = 'genotype' supports diploid data ",
         "only (ploidy = ", ploidy, ").")

  encoding <- geno_meta$metadata$encoding %||%
    detect_genotype_encoding(.meta_to_allele_matrix(geno_meta))

  # Marker order + per-marker block_id, straight from map_info / block map
  marker_meta <- as.data.table(geno_meta$map_info)[Marker %in% markers,
                   .(marker     = as.character(Marker),
                     chromosome = as.character(Chrom),
                     position   = as.numeric(Pos))]
  bmm <- geno_meta$block_marker_map_info
  if (!is.null(bmm)) {
    blk_lookup  <- as.data.table(bmm)[Marker %in% markers,
                                      .(marker = as.character(Marker),
                                        block_id = as.character(Block))]
    marker_meta <- merge(marker_meta, blk_lookup, by = "marker", all.x = TRUE)
  } else {
    marker_meta[, block_id := NA_character_]
  }
  # Chromosome is now a plotted dimension (facet columns), not just subtitle
  # text, so it is ordered naturally and carried through to the long table.
  marker_meta[is.na(chromosome), chromosome := "unplaced"]
  chrom_levels    <- .natural_chrom_levels(marker_meta$chromosome)
  marker_meta[, chrom_f := factor(chromosome, levels = chrom_levels)]
  marker_meta     <- marker_meta[order(chrom_f, position)]
  ordered_markers <- marker_meta$marker
  marker_x        <- data.table(marker = ordered_markers,
                                x_pos  = seq_along(ordered_markers))
  chr <- paste(chrom_levels, collapse = ",")

  # Sub-panels only earn their keep with more than one chromosome; a filtered
  # single-chromosome GenoMeta keeps the historic one-panel figure + subtitle.
  facet_chrom <- isTRUE(chrom_panels) && length(chrom_levels) > 1L

  # Long format built directly from the 3-D array with vectorised rep()/
  # as.character() — avoids as.data.frame.table()'s slow per-element factor
  # expansion. as.character(arr) yields the column-major flat vector (sample
  # fastest, then marker, then haplotype copy), which the rep() patterns mirror.
  n_samp     <- length(samples)
  order_res    <- .resolve_sample_order(sample_order, samples)
  sample_order <- order_res$order
  if (method == "genotype") {
    # Diploid: collapse the two haplotype alleles into one sorted genotype
    # string per (sample, marker) — e.g. "A","T" -> "AT"; "1","4" -> "14".
    # Done straight on the array via the compiled kernel / vectorised fallback
    # (.collapse_genotype_sorted), so there is no melt of the full 3-D array and
    # no N x M R group-by. One row per sample.
    geno_mat <- .collapse_genotype_sorted(arr)[, ordered_markers, drop = FALSE]
    long_dt  <- data.table(
      Sample = rep(samples, times = length(ordered_markers)),
      marker = rep(ordered_markers, each = n_samp),
      allele = as.character(geno_mat))
    long_dt[, y_label := paste0(Sample, "_g")]
    y_levels <- paste0(sample_order, "_g")
  } else {
    n_mark_all <- length(markers)
    long_dt <- data.table(
      Sample        = rep(samples, times = n_mark_all * ploidy),
      marker        = rep(rep(markers, each = n_samp), times = ploidy),
      Haplotype_idx = rep(seq_len(ploidy), each = n_samp * n_mark_all),
      allele        = as.character(arr))
    long_dt <- long_dt[marker %in% ordered_markers]
    long_dt[, y_label := paste0(Sample, "_hap", Haplotype_idx)]
    y_levels <- as.character(unlist(lapply(sample_order,
                  function(s) paste0(s, "_hap", seq_len(ploidy)))))
  }
  long_dt[, y_label := factor(y_label, levels = y_levels)]
  long_dt[marker_x, on = "marker", x_pos := i.x_pos]
  long_dt[marker_meta, on = "marker", chromosome := i.chromosome]
  long_dt[, chromosome := factor(chromosome, levels = chrom_levels)]

  # Colour palette: single-allele (haplotype) or two-allele genotype, by length.
  if (is.null(color_map))
    color_map <- make_allele_color_map(long_dt$allele, encoding)

  name_theme <- if (show_individual_names)
    theme(axis.text.y = element_text(size = 8),
          axis.ticks.y = element_line(linewidth = 0.5))
  else
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  marker_theme <- if (show_marker_labels)
    theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5),
          axis.ticks.x = element_line(linewidth = 0.5))
  else
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # Block labels hang down from the top of the panel (see below); the Y scale
  # reserves that many blank row-units so no label ever covers a tile, and so
  # the labels stay clear of the chromosome strip above the panel. 10% of the
  # rows, floored at 2, keeps the reservation proportional on tall panels
  # without eating a small one.
  lab <- if (show_block_labels)
    .alignment_block_labels(marker_meta, marker_x) else NULL
  if (!is.null(lab) && nrow(lab) > 0L)
    lab[, chromosome := factor(chromosome, levels = chrom_levels)]
  head_rows <- if (!is.null(lab) && nrow(lab) > 0L)
    max(2L, ceiling(0.10 * length(y_levels))) else 0L

  p <- ggplot(long_dt, aes(x = x_pos, y = y_label, fill = allele)) +
    geom_tile() +
    scale_fill_manual(values = color_map, na.value = "grey80", name = "Allele") +
    scale_x_continuous(
      breaks = if (show_marker_labels) marker_x$x_pos else NULL,
      labels = if (show_marker_labels) marker_x$marker else NULL,
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(add = c(0.6, head_rows + 0.6))) +
    labs(subtitle = paste0(
           if (facet_chrom)
             sprintf("%d chromosomes · %d markers · %d samples",
                     length(chrom_levels), length(ordered_markers), n_samp)
           else sprintf("Chr %s", chr),
           # Provenance: a clustered panel must not be mistaken for an
           # alphabetical one (empty string when the order is the default).
           .hs_order_suffix(order_res$desc)),
         x = "Markers (genomic position)", y = "Sample / Haplotype") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank()) + name_theme + marker_theme

  # Chromosomes as sub-panels: one facet column each, width proportional to the
  # chromosome's marker count (space = "free_x"). x_pos stays the global marker
  # index — scales = "free_x" clips each panel to its own range, so nothing
  # downstream of `p$data` has to change.
  if (facet_chrom)
    p <- p + ggplot2::facet_grid(cols = ggplot2::vars(chromosome),
                                 scales = "free_x", space = "free_x")

  # Block boundaries (between markers where the block_id changes). Carried as a
  # data layer rather than a bare xintercept vector so each line lands in its
  # own chromosome panel.
  if (show_block_boundaries) {
    bnd <- .alignment_block_boundaries(marker_meta, marker_x)
    # A chromosome change is also a block change, but under faceting that line
    # falls on a panel edge and merely restates the panel gap.
    if (facet_chrom) bnd <- bnd[cross == FALSE]
    if (nrow(bnd) > 0L)
      bnd[, chromosome := factor(chromosome, levels = chrom_levels)]
    if (nrow(bnd) > 0L)
      p <- p + geom_vline(data = bnd, aes(xintercept = xintercept),
                          inherit.aes = FALSE, colour = "grey10",
                          linewidth = 0.8, linetype = "dotdash")
  }

  # Block labels, vertical, hanging from the panel top into the headroom
  # reserved above. hjust = 1 with angle = 90 anchors the END of the string at
  # y = Inf so the text grows downward, inside the panel — before 0.5.0 it grew
  # upward past the panel into a 10 pt margin and was clipped by the device.
  if (!is.null(lab) && nrow(lab) > 0L)
    p <- p + geom_text(data = lab,
                       aes(x = x_center, y = Inf, label = block_id),
                       inherit.aes = FALSE, vjust = 0.5, hjust = 1,
                       size = block_label_size, angle = 90, colour = "grey15")

  # Horizontal separators between samples
  if (show_individual_borders) {
    pl  <- if (method == "genotype") 1L else ploidy
    hbnd <- (seq_along(sample_order) - 1L) * pl + 0.5
    p <- p + geom_hline(yintercept = hbnd, colour = "black",
                        linewidth = 0.3, linetype = "dashed")
  }

  strip_theme <- if (facet_chrom)
    theme(strip.background = ggplot2::element_rect(fill = "grey92",
                                                   colour = NA),
          strip.text = if (show_chrom_labels)
            element_text(size = 9, face = "bold", colour = "grey15")
            else element_blank())
  else theme()

  p <- p +
    ggplot2::coord_cartesian(clip = "off") +
    strip_theme +
    theme(
      plot.margin         = ggplot2::margin(10, 10, 10, 10, "pt"),
      plot.title.position = "plot",
      plot.subtitle       = element_text(size = 12, colour = "grey45", hjust = 0)
    )

  # The additive Y expansion is published so a plot pinned beside this one can
  # match it and keep its rows aligned — the block-label headroom otherwise
  # slides the alignment rows relative to, say, the dendrogram in
  # plot_clustered_alignment(). Set last: `p + <layer>` would drop the attribute.
  attr(p, "row_expand") <- c(0.6, head_rows + 0.6)
  .hs_tag_sample_order(p, order_res)
}


# Internal: X positions of block boundaries (midpoint between adjacent markers
# whose block_id differs), in genomic order. Returns one row per boundary with
# the chromosome it belongs to (the left-hand marker's), so the caller can route
# each line to its facet panel, and `cross` flagging boundaries that straddle a
# chromosome break.
.alignment_block_boundaries <- function(marker_meta, marker_x) {
  empty <- data.table(chromosome = character(0L), xintercept = numeric(0L),
                      cross = logical(0L))
  mm <- merge(marker_meta[, .(marker, block_id, chromosome)], marker_x,
              by = "marker")
  mm <- mm[!is.na(block_id)][order(x_pos)]
  if (nrow(mm) < 2L) return(empty)
  chg <- which(mm$block_id[-nrow(mm)] != mm$block_id[-1L])
  if (length(chg) == 0L) return(empty)
  data.table(chromosome = mm$chromosome[chg],
             xintercept = (mm$x_pos[chg] + mm$x_pos[chg + 1L]) / 2,
             cross      = mm$chromosome[chg] != mm$chromosome[chg + 1L])
}

# Internal: per-block label X centre (mean marker position within the block),
# with the block's chromosome for facet routing.
.alignment_block_labels <- function(marker_meta, marker_x) {
  mm <- merge(marker_meta[, .(marker, block_id, chromosome)], marker_x,
              by = "marker")
  mm <- mm[!is.na(block_id)]
  if (nrow(mm) == 0L)
    return(data.table(block_id = character(0L), x_center = numeric(0L),
                      chromosome = character(0L)))
  mm[, .(x_center = mean(x_pos), chromosome = chromosome[1L]), by = block_id]
}

# Internal: per-block genomic span + highlight flag for plot_marker_signal().
# Mirrors .alignment_block_labels(): collapse markers to their block, here using
# the cumulative GenomicPos already on `pdt`. `pdt` must carry GenomicPos and
# `candidate`; `bmm` is a (Marker, Block) map. A block is "highlighted" when it
# holds at least one candidate marker.
.signal_block_bands <- function(pdt, bmm) {
  empty <- data.table(Block = character(0L), xmin = numeric(0L),
                       xmax = numeric(0L), xcenter = numeric(0L),
                       highlighted = logical(0L))
  bmm <- as.data.table(bmm)
  if (!all(c("Marker", "Block") %in% names(bmm))) return(empty)
  bmm <- bmm[, .(Marker = as.character(Marker), Block = as.character(Block))]
  d <- merge(pdt[, .(Marker = as.character(Marker), GenomicPos, candidate)],
             bmm, by = "Marker")
  d <- d[!is.na(Block)]
  if (nrow(d) == 0L) return(empty)
  d[, .(xmin = min(GenomicPos), xmax = max(GenomicPos),
        xcenter = mean(GenomicPos), highlighted = any(candidate)),
    by = Block]
}


#' Per-family local-GEBV alignment heatmaps
#'
#' Visualises the output of \code{\link{compute_profile_localgebv}} as a nested
#' list of heatmaps, one per family per generation, mirroring the family-history
#' structure produced by \code{\link{trace_family_history}}. Each heatmap shows
#' per-block local GEBVs (continuous, signed) with blocks on the X-axis (one
#' block per cell, in genomic order) and samples on the Y-axis. Offspring rows
#' are drawn above founder rows, separated by a boundary line.
#'
#' @param localgebv A wide \code{data.frame}/\code{data.table} as returned by
#'   \code{\link{compute_profile_localgebv}}: first column \code{Block}, then one
#'   column per sample (\code{level = "genotype"}) or per \code{"<Sample>.<idx>"}
#'   (\code{level = "haplotype"}).
#' @param family Either a full \code{GenoMeta} (supplies both the family grouping
#'   via \code{metadata$family_history} and \code{block_info} for chromosome
#'   labels) or a bare family-history list as returned by
#'   \code{\link{trace_family_history}} (grouping only; chromosome labels are
#'   omitted with a warning). Auto-detected.
#' @param selected_blocks Optional character vector of block IDs to keep, using
#'   the same vocabulary as \code{filter_geno_meta(blocks = )}. Genomic order is
#'   preserved.
#' @param method Row resolution. \code{"auto"} (default) infers from the
#'   \code{localgebv} column names; \code{"genotype"} = one row per sample,
#'   \code{"haplotype"} = one row per \code{"<Sample>.<idx>"}.
#' @param sample_order Optional row order \strong{within} each family panel: a
#'   character vector of sample IDs, an \code{hs_sample_clustering} from
#'   \code{\link{cluster_samples}}, \code{TRUE} to require the broadcast order
#'   (see \code{\link{use_sample_order}}), or \code{FALSE} to ignore it.
#'   \code{NULL} (default) uses the broadcast order when one is set, else the
#'   family-history order. The offspring / founder split is \strong{never}
#'   overridden — an order re-arranges rows inside each of the two groups, and
#'   samples missing from it follow the ordered ones. In \code{"haplotype"} mode
#'   the copies of a sample stay adjacent and in index order.
#' @param fill_limits Optional numeric length-2 vector clipping the diverging
#'   fill scale. \code{NULL} (default) uses a symmetric range derived from all
#'   values, so every family plot shares one comparable scale.
#' @param fill_name Legend title for the fill scale. Default
#'   \code{"Local GEBV"}. Exists so callers drawing a derived quantity on this
#'   figure can label it honestly — see \code{\link{plot_favourable_blocks}},
#'   which passes a 1/0 dummy matrix.
#' @param show_block_labels Logical; label block IDs along the bottom X-axis.
#'   Default \code{TRUE}.
#' @param show_chrom_labels Logical; label chromosomes along the top of the
#'   panel (requires \code{block_info}, i.e. a \code{GenoMeta}). Default
#'   \code{TRUE}.
#' @param show_founder_boundary Logical; draw a horizontal line between
#'   offspring and founder rows. Default \code{TRUE}.
#'
#' @return A named list keyed by generation (\code{gen_1}, \code{gen_2}, ...);
#'   each element is a named list keyed by family ID (\code{F_1}, ...) whose
#'   values are \code{ggplot} objects. Generations/families with no matching
#'   samples are dropped.
#'
#' @examples
#' \dontrun{
#' lg    <- compute_profile_localgebv(meta, eff, level = "genotype")
#' plots <- plot_localgebv_alignment(lg, meta, selected_blocks = my_blocks)
#' plots[["gen_1"]][["F_1"]]
#' }
#'
#' @seealso \code{\link{compute_profile_localgebv}},
#'   \code{\link{plot_haplotype_alignment}}, \code{\link{trace_family_history}}
#'
#' @export
plot_localgebv_alignment <- function(localgebv,
                                     family,
                                     selected_blocks       = NULL,
                                     method                = c("auto", "genotype",
                                                               "haplotype"),
                                     sample_order          = NULL,
                                     fill_limits           = NULL,
                                     fill_name             = "Local GEBV",
                                     show_block_labels     = TRUE,
                                     show_chrom_labels     = TRUE,
                                     show_founder_boundary = TRUE) {

  method <- match.arg(method)

  # 1. Validate localgebv -----------------------------------------------------
  if (!is.data.frame(localgebv) && !inherits(localgebv, "data.table"))
    stop("plot_localgebv_alignment: 'localgebv' must be a data.frame/data.table ",
         "(the output of compute_profile_localgebv()).")
  lg <- as.data.table(localgebv)
  if (ncol(lg) < 2L || names(lg)[1L] != "Block")
    stop("plot_localgebv_alignment: 'localgebv' must have 'Block' as its first ",
         "column followed by sample columns.")
  sample_cols <- setdiff(names(lg), "Block")
  if (length(sample_cols) == 0L)
    stop("plot_localgebv_alignment: 'localgebv' has no sample columns.")

  # 2. Resolve method + column -> base-sample map -----------------------------
  is_hap_cols <- all(grepl("^.*\\.\\d+$", sample_cols))
  if (method == "auto")
    method <- if (is_hap_cols) "haplotype" else "genotype"
  if (method == "haplotype" && !is_hap_cols)
    stop("plot_localgebv_alignment: method = 'haplotype' but columns are not in ",
         "'<Sample>.<idx>' form; pass a haplotype-level localgebv.")

  base_sample <- if (method == "haplotype")
    sub("\\.\\d+$", "", sample_cols) else sample_cols
  col_map <- data.table(col = sample_cols, sample = base_sample)

  # Row order, resolved ONCE for the whole nested figure: a resolution per family
  # would emit the "not in the supplied order" message once per panel. `default =
  # "keep"` means "no order supplied and none broadcast" returns NULL, and every
  # panel keeps the family-history order it had before 0.5.5, byte for byte.
  order_res <- .hs_resolve_sample_order(
    sample_order, unique(col_map$sample),
    caller = "plot_localgebv_alignment",
    subset = "intersect", default = "keep")

  # 3. Resolve family (auto-detect GenoMeta vs family-history list) -----------
  if (is_geno_meta(family)) {
    history    <- family$metadata$family_history
    block_info <- family$block_info
    if (is.null(history))
      stop("plot_localgebv_alignment: GenoMeta carries no family history. ",
           "Run attach_family_history() first.")
  } else if (is.list(family) && length(family) > 0L &&
             all(vapply(family, function(g)
               (is.data.frame(g) || inherits(g, "data.table")) &&
               all(c("FID", "IID", "Label") %in% names(g)), logical(1L)))) {
    history    <- family
    block_info <- NULL
    warning("plot_localgebv_alignment: a bare family-history list was supplied; ",
            "chromosome labels need a GenoMeta (block_info) and will be omitted.")
  } else {
    stop("plot_localgebv_alignment: 'family' must be a GenoMeta or a ",
         "trace_family_history() list (named gen_*, each with FID/IID/Label).")
  }
  if (is.null(names(history)))
    names(history) <- paste0("gen_", seq_along(history))

  # 4. Block selection (preserve the genomic order already in lg) -------------
  if (!is.null(selected_blocks)) {
    selected_blocks <- as.character(selected_blocks)
    missing_blk <- setdiff(selected_blocks, lg$Block)
    if (length(missing_blk) > 0L)
      warning("plot_localgebv_alignment: requested block(s) not present: ",
              paste(missing_blk, collapse = ", "))
    lg <- lg[Block %in% selected_blocks]
  }
  if (nrow(lg) == 0L)
    stop("plot_localgebv_alignment: no blocks to plot after selection.")

  # 5. Block -> chromosome map + X positions ----------------------------------
  ordered_blocks <- lg$Block
  blk_x <- data.table(block = ordered_blocks, x_pos = seq_along(ordered_blocks))
  if (!is.null(block_info)) {
    bi <- as.data.table(block_info)[, .(block      = as.character(Block),
                                        chromosome = as.character(Chrom))]
    blk_x[bi, on = "block", chromosome := i.chromosome]
  } else {
    blk_x[, chromosome := NA_character_]
  }
  chrom_bnd <- .localgebv_chrom_boundaries(blk_x)
  chrom_lab <- .localgebv_chrom_labels(blk_x)

  # 6. Melt to long once: Block x col x GEBV ----------------------------------
  long_all <- data.table::melt(lg, id.vars = "Block",
                               measure.vars   = sample_cols,
                               variable.name  = "col",
                               value.name     = "GEBV",
                               variable.factor = FALSE)
  long_all[, GEBV := as.numeric(GEBV)]
  long_all[blk_x,   on = c("Block" = "block"), x_pos  := i.x_pos]
  long_all[col_map, on = "col",                sample := i.sample]

  if (is.null(fill_limits)) {
    rng <- suppressWarnings(max(abs(long_all$GEBV), na.rm = TRUE))
    if (!is.finite(rng) || rng == 0) rng <- 1
    fill_limits <- c(-rng, rng)
  }

  # 7. Iterate generations -> families ----------------------------------------
  warned_empty <- FALSE
  out <- list()
  for (gen_name in names(history)) {
    gen       <- as.data.table(history[[gen_name]])
    gen_plots <- list()

    for (fid in unique(gen$FID)) {
      grp     <- gen[FID == fid]
      off_ids <- unique(grp$IID[grp$Label == "O"])
      fnd_ids <- unique(grp$IID[grp$Label == "P"])

      # A resolved order re-arranges rows INSIDE each group; the offspring /
      # founder split and its boundary line are not the order's to move.
      off_ids <- .hs_order_within(off_ids, order_res$order)
      fnd_ids <- .hs_order_within(fnd_ids, order_res$order)

      off_cols <- col_map[sample %in% off_ids][order(match(sample, off_ids), col)]
      fnd_cols <- col_map[sample %in% fnd_ids][order(match(sample, fnd_ids), col)]
      fam_cols <- c(off_cols$col, fnd_cols$col)

      if (length(fam_cols) == 0L) {
        if (!warned_empty) {
          warning("plot_localgebv_alignment: one or more families have no ",
                  "samples in 'localgebv'; those plots are skipped.")
          warned_empty <- TRUE
        }
        next
      }

      n_off <- nrow(off_cols)
      n_fnd <- nrow(fnd_cols)
      # ggplot draws the first factor level at the bottom: founders at the
      # bottom, offspring on top => levels run founders-up-then-offspring.
      y_levels <- rev(c(off_cols$col, fnd_cols$col))

      pdat <- long_all[col %in% fam_cols]
      pdat[, y_label := factor(col, levels = y_levels)]

      p <- ggplot(pdat, aes(x = x_pos, y = y_label, fill = GEBV)) +
        geom_tile() +
        scale_fill_gradient2(
          low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
          midpoint = 0, limits = fill_limits, na.value = "grey50",
          name = fill_name) +
        scale_x_continuous(
          breaks = if (show_block_labels) blk_x$x_pos else NULL,
          labels = if (show_block_labels) blk_x$block else NULL,
          expand = c(0, 0)) +
        labs(
          title    = sprintf("%s — %s", gen_name, fid),
          subtitle = paste0(sprintf("%d offspring · %d founder(s)", n_off, n_fnd),
                            .hs_order_suffix(order_res$desc)),
          x = "Block (genomic position)", y = "Sample") +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid   = element_blank(),
          axis.text.y  = element_text(size = 8),
          axis.text.x  = if (show_block_labels)
            element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1)
            else element_blank(),
          axis.ticks.x = if (show_block_labels)
            element_line(linewidth = 0.4) else element_blank())

      # Chromosome separators + top labels (only when block_info is available)
      if (length(chrom_bnd) > 0L)
        p <- p + geom_vline(xintercept = chrom_bnd, colour = "grey10",
                            linewidth = 0.6, linetype = "dotdash")
      if (show_chrom_labels && !is.null(block_info) && nrow(chrom_lab) > 0L)
        p <- p + annotate("text", x = chrom_lab$x_center, y = Inf,
                          label = chrom_lab$chromosome, vjust = -0.4,
                          hjust = 0.5, size = 3, fontface = "bold")

      # Offspring / founder boundary
      if (show_founder_boundary && n_off > 0L && n_fnd > 0L)
        p <- p + geom_hline(yintercept = n_fnd + 0.5, colour = "black",
                            linewidth = 0.5, linetype = "dashed")

      p <- p +
        ggplot2::coord_cartesian(clip = "off") +
        theme(
          plot.margin         = ggplot2::margin(16, 10, 10, 10, "pt"),
          plot.title.position = "plot",
          plot.subtitle       = element_text(size = 10, colour = "grey45"))

      # Per-panel provenance: the order actually drawn here is this family's
      # rows, not the whole figure's, so plot_meta_alignment() compares like
      # with like when a single family panel is composed with another figure.
      gen_plots[[fid]] <- .hs_tag_sample_order(
        p, list(order = unique(c(off_cols$sample, fnd_cols$sample)),
                source = order_res$source))
    }

    if (length(gen_plots) > 0L) out[[gen_name]] <- gen_plots
  }

  if (length(out) == 0L)
    warning("plot_localgebv_alignment: no family matched the localgebv sample ",
            "columns; returning an empty list.")
  out
}


# Internal: X positions of chromosome boundaries (midpoint between adjacent
# blocks whose chromosome differs), in genomic order. `blk_dt` has columns
# block, x_pos, chromosome.
.localgebv_chrom_boundaries <- function(blk_dt) {
  bd <- blk_dt[!is.na(chromosome)][order(x_pos)]
  if (nrow(bd) < 2L) return(numeric(0))
  chg <- which(bd$chromosome[-nrow(bd)] != bd$chromosome[-1L])
  if (length(chg) == 0L) return(numeric(0))
  (bd$x_pos[chg] + bd$x_pos[chg + 1L]) / 2
}

# Internal: per-chromosome label X centre (mean block position within chrom).
.localgebv_chrom_labels <- function(blk_dt) {
  bd <- blk_dt[!is.na(chromosome)]
  if (nrow(bd) == 0L)
    return(data.table(chromosome = character(0L), x_center = numeric(0L)))
  bd[, .(x_center = mean(x_pos)), by = chromosome]
}


# Internal: shared body for the per-block/marker genomic scatter ("Manhattan").
# `dt` carries Chrom, Pos and the value column `value_col`; points run along the
# cumulative genomic axis, coloured by chromosome. `zero_line` adds a dashed
# y = 0 guide (used for signed effects, not for the non-negative variance).
# plot_manhattan_effect()/plot_manhattan_variance() are thin wrappers over this —
# the only real difference between the two figures is which value column is drawn.
.plot_manhattan <- function(dt, value_col, title, y_lab,
                            point_size = 1, alpha = 0.8, zero_line = FALSE) {
  dt <- as.data.table(dt)

  pos_info <- .create_genomic_positions(dt, "Chrom", "Pos")
  plot_dt  <- pos_info$data
  chr_info <- pos_info$chr_info
  plot_dt[, PlotVal := get(value_col)]

  p <- ggplot(plot_dt, aes(x = GenomicPos, y = PlotVal,
                           colour = factor(Chrom))) +
    geom_point(size = point_size, alpha = alpha) +
    scale_x_continuous(
      breaks = chr_info$mid,
      labels = chr_info$Chrom,
      expand = c(0.01, 0)
    ) +
    labs(title = title, x = "Chromosome", y = y_lab, colour = "Chr") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank())

  if (zero_line)
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50")
  p
}


#' Manhattan plot of per-block or per-marker effects
#'
#' @param effect_data  data.frame / data.table with columns \code{Chrom},
#'   \code{Pos}, \code{Effect}. Additional columns are ignored.
#' @param trait_name   Character label used in the plot title.
#' @param point_size   Point size (default \code{1}).
#' @param alpha        Point transparency (default \code{0.8}).
#' @param title_suffix Optional suffix appended to the plot title.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{plot_manhattan_variance}}
#' @export
plot_manhattan_effect <- function(effect_data,
                                   trait_name   = "Trait",
                                   point_size   = 1,
                                   alpha        = 0.8,
                                   title_suffix = "") {

  dt <- as.data.table(effect_data)
  if (!all(c("Chrom","Pos","Effect") %in% names(dt)))
    stop("plot_manhattan_effect: need columns Chrom, Pos, Effect.")

  .plot_manhattan(dt, value_col = "Effect",
                  title = paste0("Block Effect — ", trait_name, title_suffix),
                  y_lab = "Effect", point_size = point_size, alpha = alpha,
                  zero_line = TRUE)
}


#' Genomic scatter of per-block effect variance
#'
#' Draws one point per block at its genomic position, y = \code{BlockVariance}
#' (the across-sample variance of the block's local GEBV). Same genomic-scatter
#' geometry as \code{\link{plot_manhattan_effect}} — it shares the internal
#' plotting body, differing only in the value column.
#'
#' @param variance_data data.frame / data.table with columns \code{Chrom},
#'   \code{Pos}, \code{BlockVariance} (e.g. the block-summary table from
#'   \code{\link{run_haplotype_catalogue}}).
#' @param trait_name   Character label for the plot title.
#' @param point_size   Point size (default \code{1}).
#' @param alpha        Point transparency (default \code{0.8}).
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{plot_manhattan_effect}}
#' @export
plot_manhattan_variance <- function(variance_data,
                                     trait_name = "Trait",
                                     point_size = 1,
                                     alpha      = 0.8) {

  dt <- as.data.table(variance_data)
  if (!all(c("Chrom","Pos","BlockVariance") %in% names(dt)))
    stop("plot_manhattan_variance: need columns Chrom, Pos, BlockVariance.")

  .plot_manhattan(dt, value_col = "BlockVariance",
                  title = paste0("Block Variance — ", trait_name),
                  y_lab = "Block Effect Variance",
                  point_size = point_size, alpha = alpha)
}


#' Stacked bar chart of positive / zero / negative block counts per sample
#'
#' @param haplotype_counts data.frame with columns \code{Sample},
#'   \code{Label} (values: \code{"Positive"}, \code{"Zero"},
#'   \code{"Negative"}), \code{Count}, and optionally \code{GEBV}.
#' @param trait_name Character label for the plot title.
#' @param x_text_size Font size for X-axis sample labels (default \code{6}).
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_block_counts <- function(haplotype_counts,
                               trait_name  = "Trait",
                               x_text_size = 6) {

  haplotype_counts$Label <- factor(
    haplotype_counts$Label,
    levels = c("Positive", "Zero", "Negative")
  )

  ggplot(haplotype_counts,
         aes(x = Sample, y = Count, fill = Label)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = c(Positive = "#2ca25f", Zero = "#d9d9d9", Negative = "#de2d26")
    ) +
    labs(
      title = paste0("Block Haplotype Counts — ", trait_name),
      x     = "Sample",
      y     = "Number of Blocks",
      fill  = "Block Effect"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,
                                  vjust = 0.5, size = x_text_size)
    )
}


#' Manhattan-style plot of marker signal scores
#'
#' Visualises a per-marker signal metric (from \code{\link{score_marker_signal}})
#' along the genome and highlights the \strong{signal markers} — the markers the
#' selection step (\code{\link{select_signal_markers}}) chose to represent their
#' blocks. Points run by genomic position. When block membership is available, a
#' strip along the top groups markers into their haplotype blocks: every block is
#' drawn as a shaded rectangle spanning its genomic range, and blocks holding a
#' highlighted marker are emphasised in colour and labelled — making it clear
#' which block(s) the signal markers belong to. A blue dashed line marks the
#' theoretical maximum of the score — \code{1} for raw metrics and \code{0.25}
#' (\eqn{1 \times 0.5 \times 0.5}) for MAF-weighted \code{final_*} scores — so
#' achieved values can be read against the best attainable.
#'
#' Marker labels and the plot title are \strong{off by default} so the vertical
#' block labels stay readable; set \code{show_marker_labels = TRUE} /
#' \code{show_title = TRUE} for the pre-0.4.11 look.
#'
#' @param metrics A \code{data.frame}/\code{data.table} with \code{Chrom},
#'   \code{Marker} and the score column (e.g. the \code{metrics} element of
#'   \code{\link{score_marker_signal}}). If it lacks a \code{Pos} column,
#'   positions are looked up from \code{geno_meta} or \code{map_info}.
#' @param geno_meta Optional \code{GenoMeta}; its \code{map_info} supplies
#'   per-marker \code{Pos} (and \code{Chrom}) when \code{metrics} has no
#'   \code{Pos}, and its \code{block_info} supplies the stored signal-marker
#'   selection (see \code{\link{attach_signal_markers}}).
#' @param map_info Optional map \code{data.frame} (\code{Chrom, Marker, Pos});
#'   overrides \code{geno_meta} for position lookup.
#' @param score Name of the column to plot. Default \code{"final_mcc"}
#'   (MAF-weighted MCC); any metric / \code{final_*} column is accepted.
#' @param quantile Upper-tail quantile (0-1) drawn as a dashed reference line,
#'   and — under \code{highlight = "quantile"} — the rule defining which markers
#'   are highlighted. \code{NULL} draws no line. Default \code{0.95}.
#' @param highlight Which markers to colour: \code{"auto"} (default) uses
#'   \code{signal_markers} if given, else the selection stored on
#'   \code{geno_meta}, else falls back to the \code{quantile} rule with a
#'   message; \code{"signal"} requires a selection (errors without one);
#'   \code{"quantile"} forces the top-\code{quantile} rule; \code{"none"}
#'   highlights nothing.
#' @param signal_markers Optional explicit selection: a character vector of
#'   marker IDs, or a selection table (the \code{signal} element of
#'   \code{\link{score_marker_signal}} / the output of
#'   \code{\link{select_signal_markers}}).
#' @param top_n Number of highest-scoring markers to label — only used when
#'   \code{show_marker_labels = TRUE}. Default \code{10}.
#' @param show_marker_labels Logical; print marker IDs next to the labelled
#'   points. Default \code{FALSE} (they crowd the block labels).
#' @param show_title Logical; draw a title. Default \code{FALSE}. Supplying
#'   \code{title} implies \code{TRUE}.
#' @param point_size Point size. Default \code{1.2}.
#' @param show_blocks Logical; draw the haplotype-block strip along the top.
#'   Requires a block map (see \code{block_marker_map}). When \code{TRUE} but no
#'   map is available, the strip is skipped with a message. Default \code{TRUE}.
#' @param block_marker_map Optional \code{data.frame}/\code{data.table} with
#'   \code{Marker} and \code{Block} columns mapping markers to haplotype blocks;
#'   overrides \code{geno_meta$block_marker_map_info} (e.g. from
#'   \code{\link{attach_block_info}}).
#' @param table_file Optional path; when given, the marker × block metrics table
#'   backing the plot is written there with \code{fwrite()}. The table is always
#'   available as \code{attr(p, "signal_table")}, and
#'   \code{\link{save_marker_signal}} writes it beside the image under a matching
#'   name.
#' @param title Optional plot title (implies \code{show_title = TRUE}).
#'
#' @return A \code{ggplot} object carrying the marker × block metrics table as
#'   the \code{"signal_table"} attribute.
#'
#' @examples
#' \dontrun{
#' sig  <- score_marker_signal(meta, mode = "IBD", n_signal = 1)
#' meta <- attach_signal_markers(meta, sig$signal)
#' p    <- plot_marker_signal(sig$metrics, geno_meta = meta, score = "final_mcc")
#' save_marker_signal(p, "marker_signal_IBD.png")   # + marker_signal_IBD.csv
#' }
#'
#' @seealso \code{\link{score_marker_signal}}, \code{\link{compute_signal_metrics}},
#'   \code{\link{select_signal_markers}}, \code{\link{save_marker_signal}}
#' @export
plot_marker_signal <- function(metrics,
                               geno_meta          = NULL,
                               map_info           = NULL,
                               score              = "final_mcc",
                               quantile           = 0.95,
                               highlight          = c("auto", "signal",
                                                      "quantile", "none"),
                               signal_markers     = NULL,
                               top_n              = 10,
                               show_marker_labels = FALSE,
                               show_title         = FALSE,
                               point_size         = 1.2,
                               show_blocks        = TRUE,
                               block_marker_map   = NULL,
                               table_file         = NULL,
                               title              = NULL) {

  highlight <- match.arg(highlight)

  dt <- as.data.table(data.table::copy(metrics))
  if (!"Marker" %in% names(dt))
    stop("plot_marker_signal: 'metrics' needs a 'Marker' column.")
  if (!score %in% names(dt))
    stop("plot_marker_signal: score column '", score, "' not found in 'metrics'.")

  # Resolve genomic position / chromosome when absent.
  if (!all(c("Chrom", "Pos") %in% names(dt))) {
    mi <- if (!is.null(map_info)) as.data.table(map_info)
          else if (!is.null(geno_meta)) as.data.table(geno_meta$map_info)
          else stop("plot_marker_signal: 'metrics' lacks Chrom/Pos; supply ",
                    "'geno_meta' or 'map_info'.")
    mi <- mi[, .(Marker = as.character(Marker),
                 Chrom  = as.character(Chrom),
                 Pos    = as.numeric(Pos))]
    dt[, Marker := as.character(Marker)]
    if (!"Chrom" %in% names(dt)) dt[mi, on = "Marker", Chrom := i.Chrom]
    dt[mi, on = "Marker", Pos := i.Pos]
  }

  dt[, score_val := as.numeric(get(score))]
  dt <- dt[is.finite(score_val) & !is.na(Pos)]
  if (nrow(dt) == 0L)
    stop("plot_marker_signal: no markers with a finite '", score, "' and position.")

  pos_info <- .create_genomic_positions(dt, "Chrom", "Pos")
  pdt      <- pos_info$data
  chr_info <- pos_info$chr_info

  cut <- if (!is.null(quantile))
    as.numeric(stats::quantile(pdt$score_val, probs = quantile, na.rm = TRUE))
    else NA_real_

  # ── Which markers are highlighted ──────────────────────────────────────────
  # The selection made by select_signal_markers() is the default source: the
  # plot renders a decision rather than re-deriving one at draw time. The
  # quantile rule remains available (and is the fallback) for metrics tables
  # that were never put through a selection.
  sel_markers <- .resolve_signal_markers(signal_markers, geno_meta)
  use_signal  <- switch(
    highlight,
    signal   = if (length(sel_markers) == 0L)
                 stop("plot_marker_signal: highlight = 'signal' but no signal ",
                      "markers were supplied or stored on 'geno_meta' — run ",
                      "select_signal_markers() / attach_signal_markers() first.")
               else TRUE,
    quantile = FALSE,
    none     = FALSE,
    auto     = length(sel_markers) > 0L)

  if (highlight == "auto" && !use_signal)
    message("plot_marker_signal: no signal-marker selection found; ",
            "highlighting the top ", if (is.null(quantile)) "(none)"
            else paste0(100 * (1 - quantile), "%"),
            " by ", score, " instead. See select_signal_markers().")

  pdt[, candidate := if (use_signal) as.character(Marker) %in% sel_markers
                     else if (highlight == "none") FALSE
                     else !is.na(cut) & score_val >= cut]
  pdt[, is_signal := as.character(Marker) %in% sel_markers]

  # Resolve marker -> block map for the top block strip (arg overrides geno_meta).
  bands <- NULL
  bmm   <- if (!is.null(block_marker_map)) block_marker_map
           else if (!is.null(geno_meta)) geno_meta$block_marker_map_info
           else NULL
  if (isTRUE(show_blocks)) {
    if (is.null(bmm)) {
      message("plot_marker_signal: no block map (block_marker_map / ",
              "geno_meta$block_marker_map_info); skipping block strip.")
    } else {
      bands <- .signal_block_bands(pdt, bmm)
      if (nrow(bands) == 0L) bands <- NULL
    }
  }

  lab <- if (isTRUE(show_marker_labels) && top_n > 0L)
    pdt[order(-score_val)][seq_len(min(top_n, nrow(pdt)))] else pdt[0L]

  top_pts <- max(pdt$score_val, na.rm = TRUE)

  # Theoretical maximum of the score: raw metrics (mcc, f1, ...) cap at 1, while
  # MAF-weighted final_* scores are multiplied by MAF * (1 - MAF), which peaks at
  # 0.25 (MAF = 0.5). Drawn as a dashed reference line so the achieved scores can
  # be read against the best attainable value.
  score_max <- if (grepl("^final_", score)) 0.25 else 1

  p <- ggplot(pdt, aes(x = GenomicPos, y = score_val, colour = candidate))

  # Chromosome barriers: vertical separators between chromosomes (drawn first,
  # behind the points). Reuses the boundary helper from the local-GEBV plot.
  chrom_bnd <- .localgebv_chrom_boundaries(
    data.table(x_pos = pdt$GenomicPos, chromosome = as.character(pdt$Chrom)))
  if (length(chrom_bnd) > 0L)
    p <- p + geom_vline(xintercept = chrom_bnd, colour = "grey70",
                        linewidth = 0.4, linetype = "dashed")

  # Faint full-height shading behind highlighted blocks (drawn first, so points
  # sit on top) — visually ties each highlighted marker to its top-strip block.
  if (!is.null(bands)) {
    hb <- bands[highlighted == TRUE]
    if (nrow(hb) > 0L)
      p <- p + ggplot2::geom_rect(data = hb,
                         aes(xmin = xmin, xmax = xmax, ymin = -Inf,
                             ymax = top_pts),
                         inherit.aes = FALSE, fill = "#d7301f", alpha = 0.07)
  }

  # Title is opt-in (an explicit `title` implies it): a title eats the headroom
  # the vertical block labels grow into.
  draw_title <- isTRUE(show_title) || !is.null(title)
  plot_title <- if (draw_title) title %||% paste0("Marker signal — ", score)
                else NULL

  p <- p +
    geom_point(size = point_size, alpha = 0.8) +
    scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#d7301f"),
                       guide = "none") +
    scale_x_continuous(breaks = chr_info$mid, labels = chr_info$Chrom,
                       expand = c(0.01, 0)) +
    labs(title = plot_title, x = "Chromosome", y = score) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  if (!is.na(cut))
    p <- p + geom_hline(yintercept = cut, linetype = "dashed",
                        colour = "grey40")

  # Reference line at the theoretical maximum attainable score.
  p <- p +
    geom_hline(yintercept = score_max, linetype = "dashed",
               colour = "#08519c", linewidth = 0.4) +
    annotate("text", x = -Inf, y = score_max,
             label = paste0("max = ", score_max),
             hjust = -0.05, vjust = -0.4, size = 2.6, colour = "#08519c")

  if (nrow(lab) > 0L)
    p <- p + annotate("text", x = lab$GenomicPos, y = lab$score_val,
                      label = lab$Marker, size = 2.6, vjust = -0.6,
                      hjust = 0.5, colour = "grey15")

  # Block strip across the top: every block as a shaded rectangle over its
  # genomic span; blocks holding a highlighted marker are coloured + labelled.
  if (!is.null(bands)) {
    yr   <- range(pdt$score_val, na.rm = TRUE)
    span <- diff(yr)
    if (!is.finite(span) || span == 0)
      span <- if (is.finite(yr[2]) && yr[2] != 0) abs(yr[2]) else 1
    strip_ymin <- yr[2] + 0.08 * span
    strip_ymax <- yr[2] + 0.16 * span

    p <- p +
      ggplot2::geom_rect(data = bands,
                aes(xmin = xmin, xmax = xmax,
                    ymin = strip_ymin, ymax = strip_ymax, fill = highlighted),
                inherit.aes = FALSE, colour = "grey60", linewidth = 0.2) +
      scale_fill_manual(values = c(`FALSE` = "grey85", `TRUE` = "#d7301f"),
                        guide = "none")

    # Vertical block labels above the strip (same style as
    # plot_haplotype_alignment): far more legible than text inside the band.
    hl <- bands[highlighted == TRUE]
    if (nrow(hl) > 0L)
      p <- p + ggplot2::geom_text(data = hl,
                         aes(x = xcenter, y = strip_ymax, label = Block),
                         inherit.aes = FALSE, size = 2.5, angle = 90,
                         hjust = 0, vjust = 0.5, colour = "grey15")

    # Headroom for the vertical block labels; a title needs its own band on top.
    p <- p +
      ggplot2::coord_cartesian(clip = "off") +
      theme(plot.margin = ggplot2::margin(if (draw_title) 58 else 46,
                                          10, 10, 10, "pt"))
  }

  # Companion marker x block table (req 4): the numbers behind the figure, so
  # the plot and its table can be written side by side under matching names.
  tab <- .marker_signal_table(pdt, bands, bmm, score)
  attr(p, "signal_table") <- tab
  if (!is.null(table_file)) {
    data.table::fwrite(tab, table_file)
    message("Marker-signal table written to: ", table_file)
  }

  p
}


# Internal: normalise the `signal_markers` argument / the selection stored on a
# GenoMeta to a plain character vector of marker IDs. Accepts a character
# vector, a select_signal_markers() table (or any data.frame with a Marker
# column), or NULL -> falls back to get_signal_markers(geno_meta).
.resolve_signal_markers <- function(signal_markers, geno_meta) {
  if (is.null(signal_markers)) {
    if (is.null(geno_meta)) return(character(0L))
    return(as.character(get_signal_markers(geno_meta, long = FALSE)))
  }
  if (is.data.frame(signal_markers)) {
    if (!"Marker" %in% names(signal_markers))
      stop("plot_marker_signal: 'signal_markers' table needs a 'Marker' column.")
    return(unique(as.character(signal_markers$Marker)))
  }
  unique(as.character(signal_markers))
}


# Internal: the marker x block table backing plot_marker_signal(). One row per
# plotted marker: its metrics, its block, whether it was selected as a signal
# marker, and the block-level aggregates drawn as the top strip.
.marker_signal_table <- function(pdt, bands, bmm, score) {
  tab <- data.table::copy(pdt)
  tab[, Marker := as.character(Marker)]

  if (!"Block" %in% names(tab)) {
    blk <- if (!is.null(bmm)) as.data.table(bmm) else NULL
    if (!is.null(blk) && all(c("Marker", "Block") %in% names(blk))) {
      blk <- unique(blk[, .(Marker = as.character(Marker),
                            Block  = as.character(Block))])
      tab[blk, on = "Marker", Block := i.Block]
    } else {
      tab[, Block := NA_character_]
    }
  }

  if (!is.null(bands)) {
    b <- as.data.table(bands)[, .(Block, block_start = xmin, block_end = xmax,
                                  block_highlighted = highlighted)]
    tab[b, on = "Block", `:=`(block_start = i.block_start,
                              block_end   = i.block_end,
                              block_highlighted = i.block_highlighted)]
  }
  tab[, block_score     := max(score_val, na.rm = TRUE), by = Block]
  tab[, block_n_markers := .N, by = Block]
  tab[, score_column    := score]

  lead <- intersect(c("Chrom", "Pos", "GenomicPos", "Marker", "Block",
                      "is_signal", "candidate", "score_column", "score_val",
                      "block_score", "block_n_markers", "block_highlighted"),
                    names(tab))
  data.table::setcolorder(tab, lead)
  data.table::setorderv(tab, intersect(c("Chrom", "Pos"), names(tab)),
                        na.last = TRUE)
  tab[]
}


#' Save a marker-signal plot together with its metrics table
#'
#' Writes the plot to \code{file} and the marker × block table behind it
#' (\code{attr(plot, "signal_table")}, attached by
#' \code{\link{plot_marker_signal}}) to a CSV of the \strong{same stem}, so the
#' figure and its numbers stay together — \code{"marker_signal.png"} is
#' accompanied by \code{"marker_signal.csv"}.
#'
#' @param plot A \code{ggplot} from \code{\link{plot_marker_signal}}.
#' @param file Image path (any extension \code{ggplot2::ggsave} understands).
#' @param width,height,dpi Passed to \code{ggplot2::ggsave}. Defaults suit the
#'   wide Manhattan layout.
#' @param table_ext Extension used for the companion table. Default
#'   \code{".csv"}.
#' @param ... Further arguments for \code{ggplot2::ggsave}.
#'
#' @return Invisibly, a character vector of the two paths written (or just the
#'   image when the plot carries no table).
#'
#' @seealso \code{\link{plot_marker_signal}}
#' @export
save_marker_signal <- function(plot, file, width = 11, height = 4.5, dpi = 130,
                               table_ext = ".csv", ...) {
  if (!inherits(plot, "ggplot"))
    stop("save_marker_signal: 'plot' must be a ggplot from plot_marker_signal().")

  ggplot2::ggsave(file, plot, width = width, height = height, dpi = dpi, ...)

  tab <- attr(plot, "signal_table")
  if (is.null(tab)) {
    warning("save_marker_signal: the plot carries no 'signal_table' attribute; ",
            "only the image was written.")
    return(invisible(file))
  }

  stem       <- sub("\\.[^./\\\\]+$", "", file)
  table_path <- paste0(stem, table_ext)
  data.table::fwrite(tab, table_path)
  message("Written: ", file, " and ", table_path)
  invisible(c(file, table_path))
}


#' Precision-recall scatter of marker signal, coloured by MAF
#'
#' A non-genomic companion to \code{\link{plot_marker_signal}}: each marker is
#' placed by its recall (sensitivity) and precision, coloured by minor allele
#' frequency. Markers toward the top-right are both sensitive and specific tags
#' of block inheritance; the MAF colour separates informative (high-MAF) from
#' near-monomorphic markers. The top \code{top_n} markers by \code{score} are
#' labelled.
#'
#' @param metrics A \code{data.frame}/\code{data.table} with \code{Marker},
#'   \code{precision}, \code{recall} and \code{MAF} (the \code{metrics} element
#'   of \code{\link{score_marker_signal}}).
#' @param score Column used to pick the markers to label. Default
#'   \code{"final_mcc"}; ignored if absent.
#' @param top_n Number of highest-scoring markers to label. Default \code{10}.
#' @param point_size Point size. Default \code{2.2}.
#' @param title Optional plot title.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{plot_marker_signal}}, \code{\link{score_marker_signal}}
#' @export
plot_signal_precision_recall <- function(metrics,
                                         score      = "final_mcc",
                                         top_n      = 10,
                                         point_size = 2.2,
                                         title      = NULL) {

  dt <- as.data.table(data.table::copy(metrics))
  need <- c("Marker", "precision", "recall", "MAF")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0L)
    stop("plot_signal_precision_recall: 'metrics' missing column(s): ",
         paste(miss, collapse = ", "))

  dt <- dt[is.finite(precision) & is.finite(recall)]
  if (nrow(dt) == 0L)
    stop("plot_signal_precision_recall: no markers with finite precision & recall.")

  lab <- if (top_n > 0L && score %in% names(dt))
    dt[is.finite(get(score))][order(-get(score))][seq_len(min(top_n, .N))]
    else dt[0L]

  p <- ggplot(dt, aes(x = recall, y = precision, colour = MAF)) +
    geom_point(size = point_size, alpha = 0.85) +
    scale_color_gradient(low = "grey80", high = "#08519c", name = "MAF") +
    labs(title = title %||% "Marker signal — precision vs recall",
         x = "Recall (sensitivity)", y = "Precision") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank()) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

  if (nrow(lab) > 0L)
    p <- p + annotate("text", x = lab$recall, y = lab$precision,
                      label = lab$Marker, size = 2.6, vjust = -0.6,
                      hjust = 0.5, colour = "grey15")

  p
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Compute cumulative genomic positions for Manhattan plots
#'
#' Converts (Chrom, Pos) pairs into a single cumulative X coordinate and
#' returns per-chromosome midpoints for X-axis labelling.
#'
#' @param dt        data.table with at least \code{chrom_col} and
#'   \code{pos_col} columns.
#' @param chrom_col Column name for chromosome (character).
#' @param pos_col   Column name for position (character).
#'
#' @return Named list with elements \code{data} (original \code{dt} plus a
#'   \code{GenomicPos} column) and \code{chr_info} (data.table with
#'   \code{Chrom} and \code{mid} for X-axis breaks).
#' @keywords internal
.create_genomic_positions <- function(dt, chrom_col = "Chrom",
                                      pos_col = "Pos") {
  dt <- as.data.table(dt)
  dt[, Chrom_ := as.character(get(chrom_col))]
  dt[, Pos_   := as.numeric(get(pos_col))]

  # Sort chromosomes by natural order (numeric then alphanumeric)
  chr_order <- unique(dt$Chrom_)
  chr_nums  <- suppressWarnings(as.integer(chr_order))
  chr_order <- chr_order[order(!is.na(chr_nums), chr_nums,
                                chr_order, na.last = FALSE)]

  # Cumulative offset per chromosome
  # NB: precompute the reorder index in this frame. Doing the match() inside
  # chr_max[...] would let data.table's NSE bind `chr_max` to the same-named
  # column (an atomic vector), breaking `chr_max$Chrom_`.
  chr_max <- dt[, .(chr_max = max(Pos_, na.rm = TRUE)), by = Chrom_]
  chr_idx <- match(chr_order, chr_max$Chrom_)
  chr_max <- chr_max[chr_idx]
  offsets <- cumsum(c(0, chr_max$chr_max[-nrow(chr_max)]))
  names(offsets) <- chr_order

  dt[, GenomicPos := Pos_ + offsets[Chrom_]]

  chr_info <- dt[, .(mid = mean(GenomicPos, na.rm = TRUE)), by = .(Chrom = Chrom_)]
  info_idx <- match(chr_order, chr_info$Chrom)
  chr_info <- chr_info[info_idx]

  dt[, c("Chrom_","Pos_") := NULL]
  list(data = dt, chr_info = chr_info)
}
