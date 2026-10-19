# =============================================================================
# catalogue_circos_meta.R
#
# BAN (thread #033) stage 6 — the catalogue-side circos meta plot (Q25.n) and
# the journal figure standard (dev/figure_standards_nature.md):
#
#   circos_track()            — one ring's spec (data + mark type + palette)
#   track_annotation()        — gene / QTL counts from an annotation df
#   track_block_effect()      — one trait ring from analyze_block_catalogue()
#   track_block_variance()    — per-block variance (line or rect)
#   track_marker_stat()       — marker-level points
#   track_marker_labels()     — named loci with leader lines + deconfliction
#   plot_meta_circos()        — chr–pos bar outline + the ranked track list
#   theme_hs_journal()        — the journal ggplot theme (post hoc)
#   save_journal_figure()     — mm-width cairo_pdf/PNG saver; re-renders a
#                               circos at the declared width, never rescales
#
# The layer boundary (Q25 owner clarification): these functions CONSUME
# annotation/catalogue tables in the standard formats — from
# compute_block_annotation(), a GenoMeta, re-read step files, or a hand-built
# data.frame — and never call back into the annotation module.
#
# Sizing contract (dev/figure_standards_nature.md): figures are authored at
# final print width. Text is specified in POINTS AT THE DECLARED WIDTH
# (.pt_cex), ring heights in millimetres (.mm_frac), strokes in points
# (.lwd_pt) — so the 5 pt floor holds when the same object is re-rendered at
# another column width.
# =============================================================================


# --- journal sizing helpers --------------------------------------------------

# cex that renders `pt` points on the CURRENT device (device pointsize-aware).
.pt_cex <- function(pt) pt / par("ps")

# lwd that renders `pt` points (R lwd 1 = 1/96 in = 0.75 pt).
.lwd_pt <- function(pt) pt / 0.75

# circlize track.height fraction for a physical ring height in mm at a
# declared figure width. The unit circle (radius 1) spans width_mm / (2 *
# canvas_half) mm, canvas_half being the canvas.xlim half-width.
.mm_frac <- function(mm, width_mm, canvas_half = 1.06) {
  mm * 2 * canvas_half / width_mm
}

# Okabe–Ito, the colourblind-safe categorical order (fixed, never cycled) —
# the chromosome-bar palette.
.hs_okabe_ito <- function() c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                              "#56B4E9", "#D55E00", "#F0E442", "#999999")

# The package's diverging pair (blue / neutral grey / red, anchored at 0 —
# matches plot_block_effect_circos()) and a single-hue sequential ramp.
.hs_diverging <- function() c("#2166AC", "grey92", "#B2182B")
.hs_sequential <- function() c("grey96", "#08519C")

# Natural chromosome order: numeric where possible, then string.
.hs_chr_order <- function(x) {
  u <- unique(as.character(x))
  n <- suppressWarnings(as.integer(u))
  u[order(!is.na(n), n, u, na.last = TRUE)]
}


# --- track specs -------------------------------------------------------------

#' One ring of a circos meta plot
#'
#' The track spec consumed by \code{\link{plot_meta_circos}} — data plus a mark
#' type, in the catalogue layer's standard coordinate columns: \code{Chrom}
#' with either \code{Start}/\code{End} (interval marks) or \code{Pos} (point
#' marks). The five \code{track_*()} constructors build these from the common
#' sources; call this directly for anything else.
#'
#' @param data data.frame/data.table with \code{Chrom} and \code{Start}/
#'   \code{End} (for \code{"rect"}) or \code{Pos} (for \code{"points"},
#'   \code{"line"}, \code{"labels"}; \code{"rect"} data carrying only
#'   \code{Pos} is widened to a visible sliver). \code{"labels"} also needs
#'   \code{Name}.
#' @param type Mark: \code{"rect"} (filled block segments), \code{"points"},
#'   \code{"line"}, or \code{"labels"} (leader line + text).
#' @param value_col Column mapped to the mark's value: fill for \code{"rect"},
#'   y for \code{"points"}/\code{"line"}, the ranking for \code{"labels"}
#'   (optional there).
#' @param name   Ring name for the centre legend.
#' @param palette For \code{"rect"}: length-2 (sequential low/high) or
#'   length-3 (diverging low/mid/high, anchored at 0) colour vector; defaults
#'   to diverging when the values cross 0, sequential otherwise. For point /
#'   line / label marks: a single colour (defaults drawn from the Okabe–Ito
#'   order).
#' @param height_mm Ring height in mm at the declared figure width (defaults
#'   per type: rect 3.5, line 4.5, points 5.5, labels 8).
#' @param label_top_n \code{"labels"} only: at most this many labels, highest
#'   \code{value_col} first, after angular deconfliction. Default 10.
#' @param min_gap_rad \code{"labels"} only: minimum angular gap between kept
#'   labels (default 0.14 rad, the Python circos' deconfliction gap).
#'
#' @return An \code{hs_circos_track} spec.
#' @seealso \code{\link{plot_meta_circos}}
#' @export
circos_track <- function(data, type = c("rect", "points", "line", "labels"),
                         value_col = NULL, name = NULL, palette = NULL,
                         height_mm = NULL, label_top_n = 10L,
                         min_gap_rad = 0.14) {
  type <- match.arg(type)
  dt <- as.data.table(data)
  if (!"Chrom" %in% names(dt))
    stop("circos_track: 'data' needs a Chrom column.")
  dt[, Chrom := as.character(Chrom)]

  has_int <- all(c("Start", "End") %in% names(dt))
  has_pos <- "Pos" %in% names(dt)
  if (type == "rect") {
    if (!has_int && !has_pos)
      stop("circos_track: a 'rect' track needs Start/End (or Pos) columns.")
    if (!has_int) dt[, `:=`(Start = as.numeric(Pos), End = as.numeric(Pos))]
    dt[, `:=`(Start = as.numeric(Start), End = as.numeric(End))]
  } else {
    if (!has_pos) {
      if (!has_int)
        stop("circos_track: a '", type, "' track needs a Pos column ",
             "(or Start/End to take midpoints from).")
      dt[, Pos := (as.numeric(Start) + as.numeric(End)) / 2]
    }
    dt[, Pos := as.numeric(Pos)]
  }
  if (type == "labels" && !"Name" %in% names(dt))
    stop("circos_track: a 'labels' track needs a Name column.")
  if (!is.null(value_col)) {
    if (!value_col %in% names(dt))
      stop("circos_track: value_col '", value_col, "' not in 'data'.")
    dt[, .val := as.numeric(get(value_col))]
  } else if (type %in% c("points", "line")) {
    stop("circos_track: '", type, "' tracks need a value_col.")
  } else {
    dt[, .val := NA_real_]
  }

  if (is.null(height_mm))
    height_mm <- c(rect = 3.5, points = 5.5, line = 4.5, labels = 8)[[type]]

  structure(list(data = dt[], type = type, value_col = value_col,
                 name = name %||% value_col %||% type, palette = palette,
                 height_mm = as.numeric(height_mm),
                 label_top_n = as.integer(label_top_n),
                 min_gap_rad = as.numeric(min_gap_rad)),
            class = "hs_circos_track")
}


#' Annotation count ring (genes / QTLs per block)
#'
#' Builds a \code{"rect"} track from a block-annotation table in the standard
#' format — the \code{blocks} table of \code{\link{compute_block_annotation}},
#' a \code{GenoMeta} carrying \code{$annotation}, or any data.frame with
#' \code{Chrom, Block, Start, End} plus the count column (e.g. a re-read
#' \code{step3_local_annotation.tsv}, renamed to the container convention).
#' Per the thread-#033 layer boundary, the source never matters.
#'
#' @param x GenoMeta with \code{$annotation}, a
#'   \code{compute_block_annotation()} result, or a blocks-format data.frame.
#' @param value Count column to fill by: \code{"gene_count"} (default) or
#'   \code{"qtl_count"}.
#' @param name,palette,height_mm Passed to \code{\link{circos_track}}
#'   (sequential palette by default).
#' @return An \code{hs_circos_track}.
#' @export
track_annotation <- function(x, value = "gene_count", name = value,
                             palette = NULL, height_mm = NULL) {
  bl <- if (is.data.frame(x)) x
        else if (.hs_is_annotation(x)) x$blocks
        else if (is.list(x) && !is.null(x$annotation)) x$annotation$blocks
        else stop("track_annotation: 'x' must be a blocks-format data.frame, ",
                  "a compute_block_annotation() result, or a GenoMeta ",
                  "carrying $annotation.")
  if (!value %in% names(bl))
    stop("track_annotation: no '", value, "' column on the blocks table.")
  circos_track(bl, type = "rect", value_col = value, name = name,
               palette = palette %||% .hs_sequential(), height_mm = height_mm)
}


#' Block-effect ring for one trait
#'
#' Builds a \code{"rect"} track from an \code{\link{analyze_block_catalogue}}
#' result: each block's selected-subtype value for one trait, on the diverging
#' scale anchored at 0 (the same encoding as
#' \code{\link{plot_block_effect_circos}}). One call per ring — a multi-trait
#' figure stacks one of these per trait in the \code{tracks} list.
#'
#' @param analysis The \code{analyze_block_catalogue()} list.
#' @param trait    Trait to draw (must exist in the analysis tables).
#' @param value_source \code{"raw"} / \code{"normalized"} / \code{NULL}
#'   (the analysis' own choice), as in the other catalogue plotters.
#' @param name,palette,height_mm Passed to \code{\link{circos_track}}.
#' @return An \code{hs_circos_track}.
#' @export
track_block_effect <- function(analysis, trait, value_source = NULL,
                               name = trait, palette = NULL,
                               height_mm = NULL) {
  .check_analysis(analysis, c("S", "N"), "track_block_effect")
  vs <- .resolve_value_source(analysis, value_source)
  dt <- as.data.table(vs$table)
  if (!trait %in% unique(as.character(dt$Trait)))
    stop("track_block_effect: trait '", trait, "' not in the analysis (has: ",
         paste(unique(dt$Trait), collapse = ", "), ").")
  dt <- dt[Trait == trait]
  dt[, .plotval := get(vs$col)]
  if (!all(c("Block_start", "Block_end") %in% names(dt)))
    dt[, `:=`(Block_start = NA_real_, Block_end = NA_real_)]
  dt[is.na(Block_start), Block_start := Pos]
  dt[is.na(Block_end),   Block_end   := Pos]
  circos_track(dt[, .(Chrom, Start = Block_start, End = Block_end,
                      value = .plotval)],
               type = "rect", value_col = "value", name = name,
               palette = palette %||% .hs_diverging(), height_mm = height_mm)
}


#' Block-variance ring
#'
#' Builds a \code{"line"} (default) or \code{"rect"} track from any per-block
#' variance table carrying \code{Chrom} + \code{Start}/\code{End} (the
#' \code{Start_Pos}/\code{End_Pos} spellings of HapSelect-style sheets are
#' accepted) and a variance column.
#'
#' @param x data.frame with the block coordinates and \code{value_col}.
#' @param value_col Variance column (default \code{"Block_Var"}).
#' @param trait Optional: filter \code{x} to \code{Trait == trait} first.
#' @param type \code{"line"} (midpoint profile; default) or \code{"rect"}.
#' @param name,palette,height_mm Passed to \code{\link{circos_track}}.
#' @return An \code{hs_circos_track}.
#' @export
track_block_variance <- function(x, value_col = "Block_Var", trait = NULL,
                                 type = c("line", "rect"),
                                 name = "block variance", palette = NULL,
                                 height_mm = NULL) {
  type <- match.arg(type)
  dt <- as.data.table(x)
  if (!is.null(trait)) {
    if (!"Trait" %in% names(dt))
      stop("track_block_variance: 'trait' given but 'x' has no Trait column.")
    dt <- dt[as.character(Trait) == trait]
  }
  if ("Start_Pos" %in% names(dt) && !"Start" %in% names(dt))
    setnames(dt, c("Start_Pos", "End_Pos"), c("Start", "End"))
  if (!value_col %in% names(dt))
    stop("track_block_variance: no '", value_col, "' column.")
  circos_track(dt, type = type, value_col = value_col, name = name,
               palette = palette %||%
                 (if (type == "rect") .hs_sequential() else "#009E73"),
               height_mm = height_mm)
}


#' Marker-statistic ring
#'
#' Builds a \code{"points"} track from marker-level values: either a
#' data.frame with \code{Chrom, Pos} + \code{stat}, or a \code{GenoMeta} plus
#' the name of a stored effect trait (then \code{stat} values are
#' \code{|Effect|} per marker, joined to \code{map_info} positions; zero
#' effects are dropped).
#'
#' @param x data.frame or GenoMeta.
#' @param stat Column name in \code{x} (data.frame route) or trait name in
#'   \code{x$marker_effects} (GenoMeta route).
#' @param name,palette,height_mm Passed to \code{\link{circos_track}}.
#' @return An \code{hs_circos_track}.
#' @export
track_marker_stat <- function(x, stat, name = stat, palette = NULL,
                              height_mm = NULL) {
  if (is.data.frame(x)) {
    dt <- as.data.table(x)
    if (!all(c("Chrom", "Pos", stat) %in% names(dt)))
      stop("track_marker_stat: 'x' needs Chrom, Pos and '", stat,
           "' columns.")
    dt <- dt[, .(Chrom, Pos, value = as.numeric(get(stat)))]
  } else if (is.list(x) && !is.null(x$map_info)) {
    eff <- .coerce_single_effect_df(x$marker_effects, trait = stat)
    if (is.null(eff))
      stop("track_marker_stat: the GenoMeta carries no marker effects.")
    map <- as.data.table(x$map_info)
    dt <- as.data.table(eff)[, .(Marker, value = abs(as.numeric(Effect)))]
    dt <- dt[, .(value = max(value)), by = Marker]
    dt <- merge(dt, map[, .(Marker, Chrom, Pos = as.numeric(Pos))],
                by = "Marker")[value > 0]
    if (nrow(dt) == 0L)
      stop("track_marker_stat: every |effect| for '", stat, "' is zero.")
  } else {
    stop("track_marker_stat: 'x' must be a data.frame or a GenoMeta.")
  }
  circos_track(dt, type = "points", value_col = "value",
               name = name %||% stat, palette = palette %||% "#E69F00",
               height_mm = height_mm)
}


#' Named-locus label ring
#'
#' Builds a \code{"labels"} track (leader line + text, angularly deconflicted)
#' from a marker-level annotation data.frame — the Q25.o contract: \code{Name}
#' plus either \code{Chrom}+\code{Pos} or \code{Marker} (resolved through
#' \code{\link{standardize_marker_annotations}} against \code{meta}).
#'
#' @param x data.frame with \code{Name} and \code{Chrom}+\code{Pos} or
#'   \code{Marker}; or an already-standardized table.
#' @param meta GenoMeta for the \code{Marker} route (and for \code{match_bp}
#'   matching); optional when \code{x} carries positions.
#' @param value_col Optional ranking column — the top \code{label_top_n}
#'   survive deconfliction first.
#' @param match_bp Passed to the standardizer (default the shared 10 kb).
#' @param label_top_n,min_gap_rad,name,palette,height_mm Passed to
#'   \code{\link{circos_track}}.
#' @return An \code{hs_circos_track}.
#' @export
track_marker_labels <- function(x, meta = NULL, value_col = NULL,
                                match_bp = .HS_ANNOT_DEFAULT_BP,
                                label_top_n = 10L, min_gap_rad = 0.14,
                                name = "loci", palette = NULL,
                                height_mm = NULL) {
  dt <- as.data.table(x)
  if (!all(c("Chrom", "Pos") %in% names(dt)) || !"Name" %in% names(dt)) {
    std <- standardize_marker_annotations(x, meta = meta, match_bp = match_bp)
    dt <- as.data.table(std)
  }
  keep <- c("Chrom", "Pos", "Name", intersect(value_col, names(dt)))
  circos_track(dt[, keep, with = FALSE], type = "labels",
               value_col = if (!is.null(value_col) &&
                               value_col %in% names(dt)) value_col else NULL,
               name = name, palette = palette %||% "grey20",
               height_mm = height_mm, label_top_n = label_top_n,
               min_gap_rad = min_gap_rad)
}


# --- the meta plot -----------------------------------------------------------

# Greedy angular deconfliction (the Python circos' rule): candidates sorted by
# rank value descending; keep one when it is >= min_gap radians from everything
# already kept; at most n labels.
.hs_deconflict <- function(theta, value, n, min_gap) {
  ord <- order(-ifelse(is.na(value), -Inf, value))
  kept <- integer(0L)
  for (i in ord) {
    if (length(kept) >= n) break
    if (all(abs(theta[kept] - theta[i]) >= min_gap)) kept <- c(kept, i)
  }
  sort(kept)
}

# Auto tick-label step (bp): the smallest "nice" step giving >= min_mm of arc
# between labels at the declared width. The WHOLE genome wraps the circle, so
# the arc per bp is set by the total size, not the largest chromosome.
.hs_tick_step <- function(total_bp, width_mm, min_mm = 13) {
  circum_mm <- pi * width_mm * 0.82         # ~chr-bar circumference
  mm_per_bp <- circum_mm / total_bp
  for (step in c(10, 20, 50, 100, 200, 250, 500, 1000, 2000) * 1e6)
    if (step * mm_per_bp >= min_mm) return(step)
  5000e6
}

#' Circos meta plot — a ranked track list inside a chr–pos bar outline
#'
#' The polar analogue of \code{\link{plot_meta_alignment}} (Q25.n): the caller
#' supplies an ordered list of \code{\link{circos_track}} specs and this
#' function owns the shared frame. The outermost ring is always the
#' \strong{chr–pos bar} — one sector per chromosome in natural order,
#' Okabe–Ito colours, Mb ticks — and \code{tracks[[1]]} is the outermost data
#' ring, the last element the innermost. A centre legend names each ring
#' (with a mini colourbar for \code{"rect"} rings).
#'
#' Sizing follows \code{dev/figure_standards_nature.md}: the figure is
#' authored at \code{width_mm} (89 mm single column by default) — ring heights
#' are physical mm and text sizes are points \emph{at that width}, so nothing
#' falls below the 5 pt floor. \code{\link{save_journal_figure}} re-renders
#' the returned object at the saved width rather than rescaling it. The
#' function warns when a ring lands below ~2.5 mm and suggests double column.
#'
#' Requires \pkg{circlize} (Suggests); there is no polar approximation — the
#' function stops with an install hint.
#'
#' @param tracks A \code{circos_track()} spec or a list of them, outermost
#'   first.
#' @param meta  Optional GenoMeta: supplies chromosome extents (max
#'   \code{map_info$Pos} per chromosome) when \code{chr_sizes} is absent.
#' @param chr_sizes Optional named numeric vector of chromosome lengths (bp);
#'   beats \code{meta}; the fallback is the union of the track data.
#' @param rel_heights Optional numeric multiplier per track (like
#'   \code{plot_meta_alignment()}'s \code{rel_widths}).
#' @param width_mm Declared print width in mm (default 89, single column).
#' @param gap_degree,start_degree \code{circlize::circos.par()} frame
#'   parameters.
#' @param chr_colors Optional named colours per chromosome (default
#'   Okabe–Ito in natural order).
#' @param legend Draw the centre legend (default \code{TRUE}).
#' @param draw  Render now on the current device (default \code{TRUE}).
#'
#' @return Invisibly, an \code{hs_meta_circos} object: \code{$render(width_mm)}
#'   re-draws the figure at any declared width (what
#'   \code{\link{save_journal_figure}} calls), plus the frame parameters.
#'
#' @seealso \code{\link{circos_track}}, \code{\link{track_annotation}},
#'   \code{\link{track_block_effect}}, \code{\link{track_block_variance}},
#'   \code{\link{track_marker_stat}}, \code{\link{track_marker_labels}},
#'   \code{\link{save_journal_figure}}, \code{\link{plot_block_effect_circos}}
#' @export
plot_meta_circos <- function(tracks, meta = NULL, chr_sizes = NULL,
                             rel_heights = NULL, width_mm = 89,
                             gap_degree = 2, start_degree = 90,
                             chr_colors = NULL, legend = TRUE, draw = TRUE) {
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("plot_meta_circos: the 'circlize' package is required. ",
         "Install it with install.packages('circlize').")
  if (inherits(tracks, "hs_circos_track")) tracks <- list(tracks)
  if (!is.list(tracks) || length(tracks) == 0L ||
      !all(vapply(tracks, inherits, logical(1L), "hs_circos_track")))
    stop("plot_meta_circos: 'tracks' must be circos_track() specs ",
         "(see the track_*() constructors).")
  if (!is.null(rel_heights) && length(rel_heights) != length(tracks))
    stop("plot_meta_circos: 'rel_heights' must have one value per track.")

  # --- chromosome frame: sizes + order --------------------------------------
  union_max <- function() {
    dd <- data.table::rbindlist(lapply(tracks, function(tr) {
      d <- tr$data
      hi <- if (all(c("Start", "End") %in% names(d))) pmax(d$Start, d$End)
            else d$Pos
      data.table::data.table(Chrom = d$Chrom, hi = as.numeric(hi))
    }))
    v <- dd[, max(hi, na.rm = TRUE), by = Chrom]
    stats::setNames(v$V1, v$Chrom)
  }
  sizes <- if (!is.null(chr_sizes)) {
    if (is.null(names(chr_sizes)))
      stop("plot_meta_circos: 'chr_sizes' must be named by chromosome.")
    chr_sizes
  } else if (!is.null(meta) && !is.null(meta$map_info)) {
    mp <- as.data.table(meta$map_info)
    v <- mp[, max(as.numeric(Pos), na.rm = TRUE), by = Chrom]
    stats::setNames(v$V1, as.character(v$Chrom))
  } else union_max()
  chr_levels <- .hs_chr_order(names(sizes))
  sizes <- as.numeric(sizes[chr_levels])
  names(sizes) <- chr_levels

  cols <- chr_colors %||%
    stats::setNames(rep_len(.hs_okabe_ito(), length(chr_levels)), chr_levels)

  heights_mm <- vapply(tracks, `[[`, numeric(1L), "height_mm")
  if (!is.null(rel_heights)) heights_mm <- heights_mm * rel_heights
  track_names <- vapply(tracks, `[[`, character(1L), "name")

  # --- one ring's palette ----------------------------------------------------
  ring_fill <- function(tr) {
    v <- tr$data$.val
    pal <- tr$palette
    if (tr$type != "rect") return(list(fun = NULL, col = pal[1L]))
    rng <- range(v, na.rm = TRUE)
    if (is.null(pal)) pal <- if (rng[1L] < 0 && rng[2L] > 0) .hs_diverging()
                             else .hs_sequential()
    if (length(pal) >= 3L) {
      a <- max(abs(rng)); if (!is.finite(a) || a == 0) a <- 1
      list(fun = circlize::colorRamp2(c(-a, 0, a), pal[c(1L, 2L, 3L)]),
           lab = c(-a, 0, a))
    } else {
      hi <- max(rng[2L], 1e-12)
      list(fun = circlize::colorRamp2(c(0, hi), pal[c(1L, 2L)]),
           lab = c(0, hi))
    }
  }

  render <- function(width_mm_now = width_mm) {
    small <- heights_mm[heights_mm < 2.5]
    if (length(small) > 0L)
      message("plot_meta_circos: ", length(small), " ring(s) under 2.5 mm at ",
              width_mm_now, " mm — consider double column (183 mm) or fewer ",
              "rings for journal legibility.")

    # Radial budget: ring heights are physical mm, so the centre hole shrinks
    # as rings stack. A cramped hole makes the innermost ring (and the centre
    # legend) collide — say so, and drop the legend rather than overprint it.
    radius_mm <- width_mm_now / (2 * 1.06)
    used_mm   <- 2.2 + 2.5 + sum(heights_mm) + 0.8 * (length(tracks) + 1L)
    hole_mm   <- radius_mm - used_mm
    if (hole_mm < 0.35 * radius_mm)
      message("plot_meta_circos: ", length(tracks), " rings (",
              round(sum(heights_mm), 1), " mm) leave a ",
              max(0, round(hole_mm)), " mm centre at ", width_mm_now,
              " mm width — use double column (183 mm) or fewer rings.")
    do_legend <- isTRUE(legend) && hole_mm >= 18
    if (isTRUE(legend) && !do_legend)
      message("plot_meta_circos: the centre is too small for the legend at ",
              width_mm_now, " mm — legend skipped (it fits at double ",
              "column).")

    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)
    circlize::circos.par(cell.padding = c(0, 0, 0, 0),
                         gap.degree = gap_degree,
                         start.degree = start_degree,
                         track.margin = c(.mm_frac(0.4, width_mm_now),
                                          .mm_frac(0.4, width_mm_now)),
                         canvas.xlim = c(-1.06, 1.06),
                         canvas.ylim = c(-1.06, 1.06),
                         points.overflow.warning = FALSE)
    circlize::circos.initialize(sectors = factor(chr_levels,
                                                 levels = chr_levels),
                                xlim = cbind(0, sizes))

    # Chr–pos bar outline: coloured bar, sparse labelled Mb ticks (unlabelled
    # minors between), chromosome name inside the bar.
    tick <- .hs_tick_step(sum(sizes), width_mm_now)
    circlize::circos.track(
      ylim = c(0, 1), track.height = .mm_frac(2.2, width_mm_now),
      bg.border = NA,
      panel.fun = function(x, y) {
        ch <- circlize::CELL_META$sector.index
        xl <- circlize::CELL_META$xlim
        circlize::circos.rect(xl[1L], 0, xl[2L], 1, col = cols[[ch]],
                              border = NA)
        at_minor <- seq(0, xl[2L], by = tick / 2)
        circlize::circos.axis(h = "top", major.at = at_minor, labels = FALSE,
                              col = "grey30",
                              major.tick.length = circlize::mm_y(0.5),
                              lwd = .lwd_pt(0.4))
        at <- seq(tick, xl[2L], by = tick)
        if (length(at) > 0L)
          circlize::circos.axis(h = "top", major.at = at,
                                labels = at / 1e6,
                                labels.cex = .pt_cex(5),
                                labels.col = "black",
                                # "downward" = screen-upright numerals at any
                                # angle — a tangential label at exactly 3/9
                                # o'clock sits on circlize's niceFacing flip
                                # boundary and can render bottom-to-top
                                # ("600" -> "009"); upright text cannot.
                                col = "grey30", labels.facing = "downward",
                                major.tick.length = circlize::mm_y(0.9),
                                lwd = .lwd_pt(0.5))
        # Upright ("downward") so a 6 or 9 can never flip into the other at
        # the bending-text boundaries (same trap as the axis numerals).
        circlize::circos.text(mean(xl), 0.5, ch, cex = .pt_cex(6),
                              col = "white", font = 2, facing = "downward")
      })

    # Data rings, outermost first.
    for (k in seq_along(tracks)) {
      tr <- tracks[[k]]
      th <- .mm_frac(heights_mm[k], width_mm_now)
      d  <- tr$data
      if (tr$type == "rect") {
        rf <- ring_fill(tr)
        circlize::circos.track(
          ylim = c(0, 1), track.height = th, bg.border = "grey85",
          bg.lwd = .lwd_pt(0.25),
          panel.fun = function(x, y) {
            ch <- circlize::CELL_META$sector.index
            sd <- d[d$Chrom == ch, ]
            if (nrow(sd) == 0L) return()
            xl <- circlize::CELL_META$xlim
            min_bp <- (xl[2L] - xl[1L]) * 0.004
            s <- sd$Start; e <- sd$End
            narrow <- (e - s) < min_bp
            mid <- (s + e) / 2
            s[narrow] <- pmax(xl[1L], mid[narrow] - min_bp / 2)
            e[narrow] <- pmin(xl[2L], mid[narrow] + min_bp / 2)
            circlize::circos.rect(s, rep(0, nrow(sd)), e, rep(1, nrow(sd)),
                                  col = rf$fun(sd$.val), border = NA)
          })
      } else if (tr$type %in% c("points", "line")) {
        rng <- range(d$.val, na.rm = TRUE)
        if (!all(is.finite(rng)) || diff(rng) == 0)
          rng <- rng[1L] + c(-0.5, 0.5)
        pad <- diff(rng) * 0.06
        yl <- c(rng[1L] - pad, rng[2L] + pad)
        col1 <- (tr$palette %||% "#0072B2")[1L]
        circlize::circos.track(
          ylim = yl, track.height = th, bg.border = "grey85",
          bg.lwd = .lwd_pt(0.25),
          panel.fun = function(x, y) {
            ch <- circlize::CELL_META$sector.index
            sd <- d[d$Chrom == ch, ]
            if (nrow(sd) == 0L) return()
            if (tr$type == "points") {
              circlize::circos.points(sd$Pos, sd$.val, pch = 16,
                                      cex = .pt_cex(2.2),
                                      col = grDevices::adjustcolor(col1, 0.55))
            } else {
              o <- order(sd$Pos)
              circlize::circos.lines(sd$Pos[o], sd$.val[o], col = col1,
                                     lwd = .lwd_pt(0.5))
            }
          })
      } else {                                            # labels
        # theta per candidate, deconflicted globally, then drawn per sector.
        gaps  <- gap_degree * length(chr_levels)
        avail <- 360 - gaps
        off   <- c(0, cumsum(sizes))[seq_along(sizes)]
        frac  <- (off[match(d$Chrom, chr_levels)] + d$Pos) / sum(sizes)
        theta <- (frac * avail +
                    (match(d$Chrom, chr_levels) - 1L) * gap_degree) * pi / 180
        kept  <- .hs_deconflict(theta, d$.val, tr$label_top_n, tr$min_gap_rad)
        dk    <- d[kept, ]
        col1  <- (tr$palette %||% "grey20")[1L]
        circlize::circos.track(
          ylim = c(0, 1), track.height = th, bg.border = NA,
          panel.fun = function(x, y) {
            ch <- circlize::CELL_META$sector.index
            sd <- dk[dk$Chrom == ch, ]
            if (nrow(sd) == 0L) return()
            circlize::circos.segments(sd$Pos, 1, sd$Pos, 0.82, col = col1,
                                      lwd = .lwd_pt(0.5))
            circlize::circos.text(sd$Pos, 0.76, sd$Name, cex = .pt_cex(5),
                                  col = "black", facing = "clockwise",
                                  niceFacing = TRUE, adj = c(1, 0.5))
          })
      }
    }

    # Centre legend: one row per ring, in ink (5 pt), spaced so nothing
    # collides — a mini colourbar with min/max beneath for rect rings, a
    # chip + name for the rest.
    if (do_legend) {
      row_h <- vapply(tracks, function(tr)
        if (tr$type == "rect") 6 else 3.4, numeric(1L))
      row_h <- .mm_frac(row_h, width_mm_now)
      bar_w <- .mm_frac(10, width_mm_now)
      bar_h <- .mm_frac(1.6, width_mm_now)
      tops <- cumsum(c(0, row_h[-length(row_h)]))
      total_h <- sum(row_h)
      for (k in seq_along(tracks)) {
        tr <- tracks[[k]]
        y_top <- total_h / 2 - tops[k]
        if (tr$type == "rect") {
          rf <- ring_fill(tr)
          y_bar <- y_top - .mm_frac(2.2, width_mm_now)
          graphics::text(-bar_w / 2, y_top - .mm_frac(0.7, width_mm_now),
                         tr$name, adj = c(0, 0.5), cex = .pt_cex(5),
                         col = "black")
          brks <- seq(rf$lab[1L], rf$lab[length(rf$lab)], length.out = 24L)
          xs <- seq(-bar_w / 2, bar_w / 2, length.out = 25L)
          graphics::rect(xs[-25L], y_bar - bar_h / 2, xs[-1L],
                         y_bar + bar_h / 2, col = rf$fun(brks), border = NA)
          graphics::text(-bar_w / 2 - .mm_frac(0.8, width_mm_now), y_bar,
                         signif(rf$lab[1L], 2), adj = c(1, 0.5),
                         cex = .pt_cex(5), col = "black")
          graphics::text(bar_w / 2 + .mm_frac(0.8, width_mm_now), y_bar,
                         signif(rf$lab[length(rf$lab)], 2), adj = c(0, 0.5),
                         cex = .pt_cex(5), col = "black")
        } else {
          col1 <- (tr$palette %||% "#0072B2")[1L]
          y <- y_top - .mm_frac(1.2, width_mm_now)
          if (tr$type == "points")
            graphics::points(-bar_w / 2, y, pch = 16, cex = .pt_cex(4),
                             col = col1)
          else
            graphics::segments(-bar_w / 2 - .mm_frac(1, width_mm_now), y,
                               -bar_w / 2 + .mm_frac(1, width_mm_now), y,
                               col = col1, lwd = .lwd_pt(0.75))
          graphics::text(-bar_w / 2 + .mm_frac(2, width_mm_now), y, tr$name,
                         adj = c(0, 0.5), cex = .pt_cex(5), col = "black")
        }
      }
    }
    invisible(NULL)
  }

  out <- structure(list(render = render, width_mm = width_mm,
                        tracks = track_names, chr_levels = chr_levels),
                   class = "hs_meta_circos")
  if (isTRUE(draw)) render(width_mm)
  invisible(out)
}


# --- the journal standard ----------------------------------------------------

#' The journal (Nature-quality) ggplot theme
#'
#' The \code{dev/figure_standards_nature.md} contract as a post-hoc theme:
#' 5–7 pt sans lettering, black axis text, 0.5 pt axis lines with outward
#' ticks, no grid, content-cropped 2 pt margins, borderless bottom legend,
#' and \strong{no in-figure titles} (captions carry them — a surviving title
#' is nulled). Apply with \code{p + theme_hs_journal()}; existing plotters
#' keep their interactive defaults.
#'
#' @param base_size Base text size in pt (default 7; axis text renders at
#'   ~6 pt, never below the 5 pt floor).
#' @param base_family Font family (default \code{"sans"} — Arial/Helvetica;
#'   save through \code{\link{save_journal_figure}} so cairo embeds it).
#'
#' @return A \code{ggplot2} theme object.
#' @seealso \code{\link{save_journal_figure}}
#' @export
theme_hs_journal <- function(base_size = 7, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title      = element_blank(),
      plot.subtitle   = element_blank(),
      plot.margin     = ggplot2::margin(2, 2, 2, 2, "pt"),
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      # ggplot linewidth 1 = 0.75 mm ≈ 2.13 pt, so 0.5 pt = 0.5 / 2.13
      axis.line       = element_line(colour = "black", linewidth = 0.5 / 2.13),
      axis.ticks      = element_line(colour = "black", linewidth = 0.5 / 2.13),
      axis.ticks.length = grid::unit(1.5, "mm"),
      axis.text       = element_text(colour = "black",
                                     size = max(5, base_size - 1)),
      axis.title      = element_text(colour = "black", size = base_size),
      legend.position = "bottom",
      legend.background = element_blank(),
      legend.key      = element_blank(),
      legend.key.height = grid::unit(2.5, "mm"),
      legend.key.width  = grid::unit(5, "mm"),
      legend.text     = element_text(size = max(5, base_size - 2)),
      legend.title    = element_text(size = max(5, base_size - 1)),
      legend.margin   = ggplot2::margin(1, 1, 1, 1, "pt"),
      strip.background = element_blank(),
      strip.text      = element_text(size = base_size, colour = "black",
                                     margin = ggplot2::margin(1, 1, 1, 1,
                                                              "pt"))
    )
}


#' Save a figure at journal size
#'
#' The \code{dev/figure_standards_nature.md} saver: figures land at a named
#' print width — \code{"single"} (89 mm), \code{"onehalf"} (120 mm),
#' \code{"double"} (183 mm), or a numeric mm — as vector PDF through
#' \code{grDevices::cairo_pdf} (fonts embedded) or, for \code{.png}, raster at
#' \code{dpi} through \pkg{ragg} when installed (else \code{grDevices::png}).
#'
#' \itemize{
#'   \item A \code{ggplot} with a surviving title/subtitle triggers a warning
#'     (the caption owns the title at a journal).
#'   \item An \code{hs_meta_circos} is \strong{re-rendered} at the saved width
#'     — text re-sized to its pt spec — never rescaled.
#'   \item A plain \code{recordedplot} is refused: rescaling recorded output
#'     breaks the pt contract (re-render at the right width instead).
#' }
#'
#' @param plot A \code{ggplot}, an \code{hs_meta_circos} from
#'   \code{\link{plot_meta_circos}}, or any object with a \code{plot()}
#'   method drawn at author size.
#' @param file Output path; the extension picks the device (\code{.pdf} or
#'   \code{.png}).
#' @param width \code{"single"}, \code{"onehalf"}, \code{"double"}, or mm.
#' @param height Height in mm; \code{NULL} = the width (square — the circos
#'   default) for \code{hs_meta_circos}, else 0.75 × width.
#' @param dpi Raster resolution for \code{.png} (default 600).
#'
#' @return Invisibly, \code{file}.
#' @seealso \code{\link{theme_hs_journal}}, \code{\link{plot_meta_circos}}
#' @export
save_journal_figure <- function(plot, file,
                                width = c("single", "onehalf", "double"),
                                height = NULL, dpi = 600) {
  if (inherits(plot, "recordedplot"))
    stop("save_journal_figure: refusing to rescale a recorded plot — text ",
         "sizes would drift from their pt spec. Re-render at the target ",
         "width instead (plot_meta_circos() objects do this automatically).")
  w_mm <- if (is.numeric(width)) width
          else c(single = 89, onehalf = 120, double = 183)[[match.arg(width)]]
  is_circos <- inherits(plot, "hs_meta_circos")
  h_mm <- height %||% if (is_circos) w_mm else w_mm * 0.75

  if (inherits(plot, "ggplot")) {
    tt <- plot$labels[c("title", "subtitle")]
    if (any(vapply(tt, function(x) !is.null(x) && nzchar(paste(x)[1L]),
                   logical(1L))))
      warning("save_journal_figure: the plot carries a title/subtitle — a ",
              "journal figure's caption does that job. theme_hs_journal() ",
              "removes them.")
  }

  ext <- tolower(tools::file_ext(file))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  in_w <- w_mm / 25.4; in_h <- h_mm / 25.4
  if (ext == "pdf") {
    grDevices::cairo_pdf(file, width = in_w, height = in_h, pointsize = 12)
  } else if (ext == "png") {
    if (requireNamespace("ragg", quietly = TRUE))
      ragg::agg_png(file, width = in_w, height = in_h, units = "in",
                    res = dpi)
    else
      grDevices::png(file, width = in_w, height = in_h, units = "in",
                     res = dpi, type = "cairo")
  } else {
    stop("save_journal_figure: unsupported extension '.", ext,
         "' (use .pdf or .png).")
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  if (is_circos) {
    graphics::par(mar = c(0, 0, 0, 0))
    plot$render(w_mm)
  } else if (inherits(plot, "ggplot")) {
    print(plot)
  } else {
    plot(plot)
  }
  invisible(file)
}
