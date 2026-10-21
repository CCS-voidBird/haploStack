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
#   track_marker_labels()     — named loci as regions (+ optional search
#                               window band) with deconflicted text
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
#'   \code{"line"}; \code{"rect"} data carrying only \code{Pos} is widened to
#'   a visible sliver). \code{"labels"} needs \code{Name} and a
#'   \strong{region} per locus: \code{Start}/\code{End} (a gene model, a QTL
#'   interval); a locus reported at a single coordinate \code{a} (\code{Pos}
#'   alone, only one of \code{Start}/\code{End}, or \code{Start == End})
#'   becomes \code{Start = a}, \code{End = a + 1}, and a reversed span is
#'   swapped. Optional \code{Window_start}/\code{Window_end} give each
#'   locus's search window (see \code{window_bp}).
#' @param type Mark: \code{"rect"} (filled block segments), \code{"points"},
#'   \code{"line"}, or \code{"labels"} (each locus's region as a band at the
#'   ring's outer edge, its text just inside, and — when present — its
#'   search window as a lighter, taller band behind the region, with no
#'   text).
#' @param value_col Column mapped to the mark's value: fill for \code{"rect"},
#'   y for \code{"points"}/\code{"line"}, the ranking for \code{"labels"}
#'   (optional there). For \code{"rect"}, a character or factor column makes
#'   the ring \strong{categorical}: one colour per category (see
#'   \code{palette}), drawn in level order so later levels sit on top, and
#'   one legend chip per category in place of a colourbar — e.g. every block
#'   in grey with a selected subset highlighted.
#' @param name   Ring name for the legend.
#' @param palette For \code{"rect"}: length-2 (sequential low/high) or
#'   length-3 (diverging low/mid/high, anchored at 0) colour vector; defaults
#'   to diverging when the values cross 0, sequential otherwise. For a
#'   categorical \code{"rect"} ring: a colour vector named by category (it
#'   must name every category; its order sets the level order unless the
#'   column is a factor), or an unnamed vector recycled in level order
#'   (default Okabe–Ito). For point / line marks: a single colour (defaults
#'   drawn from the Okabe–Ito order). For \code{"labels"}: the region colour,
#'   optionally followed by the search-window colour (defaults
#'   \code{"grey20"} and Okabe–Ito sky blue \code{"#56B4E9"}).
#' @param height_mm Ring height in mm at the declared figure width (defaults
#'   per type: rect 3.5, line 4.5, points 5.5, labels 8).
#' @param label_top_n \code{"labels"} only: at most this many labels, highest
#'   \code{value_col} first, after angular deconfliction. Default 10.
#' @param min_gap_rad \code{"labels"} only: minimum angular gap between kept
#'   labels (default 0.14 rad, the Python circos' deconfliction gap).
#' @param window_bp \code{"labels"} only: a search window of the region
#'   \eqn{\pm} \code{window_bp} bp (floored at 0) for every locus without its
#'   own \code{Window_start}/\code{Window_end}. \code{NULL} (default): only
#'   those columns give windows; with neither, no window is drawn.
#' @param window_name \code{"labels"} only: the search window's legend name
#'   (default \code{"search window"}). The window has no text of its own.
#'
#' @return An \code{hs_circos_track} spec.
#' @seealso \code{\link{plot_meta_circos}}
#' @export
circos_track <- function(data, type = c("rect", "points", "line", "labels"),
                         value_col = NULL, name = NULL, palette = NULL,
                         height_mm = NULL, label_top_n = 10L,
                         min_gap_rad = 0.14, window_bp = NULL,
                         window_name = "search window") {
  type <- match.arg(type)
  dt <- as.data.table(data)
  if (!"Chrom" %in% names(dt))
    stop("circos_track: 'data' needs a Chrom column.")
  dt[, Chrom := as.character(Chrom)]
  if (!is.null(window_bp)) {
    if (type != "labels")
      stop("circos_track: 'window_bp' applies to 'labels' tracks only.")
    if (!is.numeric(window_bp) || length(window_bp) != 1L ||
        !is.finite(window_bp) || window_bp < 0)
      stop("circos_track: 'window_bp' must be a single number >= 0.")
  }

  has_int <- all(c("Start", "End") %in% names(dt))
  has_pos <- "Pos" %in% names(dt)
  if (type == "labels") {
    if (!any(c("Start", "End", "Pos") %in% names(dt)))
      stop("circos_track: a 'labels' track needs Start/End (a region) or ",
           "Pos (a single position) columns.")
    dt <- .hs_label_regions(dt, window_bp)
  } else if (type == "rect") {
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
  categories <- NULL
  if (!is.null(value_col)) {
    if (!value_col %in% names(dt))
      stop("circos_track: value_col '", value_col, "' not in 'data'.")
    v <- dt[[value_col]]
    if (type == "rect" && (is.character(v) || is.factor(v))) {
      ct <- .hs_categories(v, palette)
      categories <- ct$levels
      palette <- ct$palette
      dt[, `:=`(.cat = as.character(v), .val = NA_real_)]
    } else {
      dt[, .val := as.numeric(get(value_col))]
    }
  } else if (type %in% c("points", "line")) {
    stop("circos_track: '", type, "' tracks need a value_col.")
  } else {
    dt[, .val := NA_real_]
  }

  if (is.null(height_mm))
    height_mm <- c(rect = 3.5, points = 5.5, line = 4.5, labels = 8)[[type]]

  structure(list(data = dt[], type = type, value_col = value_col,
                 name = name %||% value_col %||% type, palette = palette,
                 categories = categories,
                 height_mm = as.numeric(height_mm),
                 label_top_n = as.integer(label_top_n),
                 min_gap_rad = as.numeric(min_gap_rad),
                 window_name = as.character(window_name)),
            class = "hs_circos_track")
}

# A label ring's loci as regions. Start/End where given; a locus reported at a
# single coordinate a (Pos alone, one of Start/End, or Start == End) becomes
# [a, a + 1]; a reversed span is swapped. Pos becomes the region midpoint (the
# text anchor and deconfliction angle). The window: Window_start/Window_end
# where given (a missing side takes the region's bound), else the region +-
# window_bp floored at 0, else none (NA).
.hs_label_regions <- function(dt, window_bp = NULL) {
  num <- function(col) if (col %in% names(dt)) as.numeric(dt[[col]])
                       else rep(NA_real_, nrow(dt))
  s <- num("Start"); e <- num("End"); p <- num("Pos")
  one <- is.na(s) | is.na(e)
  a <- ifelse(is.na(s), ifelse(is.na(e), p, e), s)
  s[one] <- a[one]; e[one] <- a[one] + 1
  lo <- pmin(s, e); hi <- pmax(s, e)
  pt <- !is.na(lo) & lo == hi
  hi[pt] <- lo[pt] + 1

  ws <- num("Window_start"); we <- num("Window_end")
  half <- xor(is.na(ws), is.na(we))
  ws[half & is.na(ws)] <- lo[half & is.na(ws)]
  we[half & is.na(we)] <- hi[half & is.na(we)]
  if (!is.null(window_bp)) {
    fill <- is.na(ws) & is.na(we)
    ws[fill] <- pmax(0, lo[fill] - window_bp)
    we[fill] <- hi[fill] + window_bp
  }
  dt[, c("Start", "End", "Pos", "Window_start", "Window_end") :=
       list(lo, hi, (lo + hi) / 2, pmin(ws, we), pmax(ws, we))]
  dt
}

# A label ring's colours: c(region, window).
.hs_label_cols <- function(tr) {
  pal <- tr$palette %||% "grey20"
  c(pal[1L], if (length(pal) >= 2L) pal[2L] else "#56B4E9")
}

# A label ring's radial layout (ylim 0–1, 1 = the outer edge): the window band
# stands proud of the region band on both sides, so the two read apart even
# when both are at their minimum arc; text starts at text_y and runs inward.
# Arcs narrower than reg_mm are widened to it; a window is drawn at least
# halo_mm wider than its region on each side.
.hs_lab_geom <- list(win = c(0.82, 1), reg = c(0.86, 0.96), text_y = 0.76,
                     reg_mm = 0.35, halo_mm = 0.3)

# Widen each [s, e] to at least w (bp, recycled) about its midpoint, then
# clamp to the sector's xlim.
.hs_widen <- function(s, e, w, xl) {
  w <- rep_len(w, length(s))
  narrow <- (e - s) < w
  mid <- (s + e) / 2
  s[narrow] <- mid[narrow] - w[narrow] / 2
  e[narrow] <- mid[narrow] + w[narrow] / 2
  list(s = pmax(xl[1L], s), e = pmin(xl[2L], e))
}

# A categorical rect ring's levels and colours. Levels: factor order, else the
# named palette's order (categories it does not name follow, sorted), else
# sorted. Colours: the named palette (it must name every level), else an
# unnamed palette — or Okabe–Ito — in level order.
.hs_categories <- function(v, palette) {
  present <- unique(as.character(v[!is.na(v)]))
  lev <- if (is.factor(v)) intersect(levels(v), present)
         else if (!is.null(names(palette)))
           c(intersect(names(palette), present),
             sort(setdiff(present, names(palette))))
         else sort(present)
  if (!is.null(names(palette))) {
    miss <- setdiff(lev, names(palette))
    if (length(miss) > 0L)
      stop("circos_track: the named 'palette' has no colour for ",
           paste0("'", miss, "'", collapse = ", "), ".")
    pal <- palette[lev]
  } else {
    pal <- stats::setNames(rep_len(palette %||% .hs_okabe_ito(), length(lev)),
                           lev)
  }
  list(levels = lev, palette = pal)
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
#' Builds a \code{"labels"} track from a locus-level annotation data.frame:
#' each locus is drawn as its \strong{region} (a band at the ring's outer
#' edge) with its name just inside, the names angularly deconflicted. The
#' Q25.o contract plus regions: \code{Name} with \code{Chrom} +
#' \code{Start}/\code{End} (a gene model, a QTL interval), or \code{Chrom} +
#' \code{Pos}, or \code{Marker} (resolved through
#' \code{\link{standardize_marker_annotations}} against \code{meta}). A locus
#' reported at a single position \code{a} is drawn as \code{[a, a + 1]}.
#'
#' A \strong{search window} — the stretch searched around each locus, e.g.
#' for markers or blocks — is drawn behind its region as a lighter, taller
#' band with no text; the legend names it (\code{window_name}). Give it per
#' locus as \code{Window_start}/\code{Window_end} columns, or for every locus
#' as \code{window_bp} either side of the region.
#'
#' @param x data.frame with \code{Name} and \code{Chrom} + \code{Start}/
#'   \code{End}, \code{Chrom} + \code{Pos}, or \code{Marker}; optional
#'   \code{Window_start}/\code{Window_end}; or an already-standardized table.
#' @param meta GenoMeta for the \code{Marker} route (and for \code{match_bp}
#'   matching); optional when \code{x} carries positions.
#' @param value_col Optional ranking column — the top \code{label_top_n}
#'   survive deconfliction first.
#' @param match_bp Passed to the standardizer (default the shared 10 kb).
#' @param window_bp Optional search-window flank (bp) either side of each
#'   region, for loci without \code{Window_start}/\code{Window_end}.
#'   \code{NULL} (default) draws windows only from those columns.
#' @param label_top_n,min_gap_rad,name,palette,height_mm,window_name Passed
#'   to \code{\link{circos_track}} (\code{palette}: the region colour, then
#'   optionally the window colour).
#' @return An \code{hs_circos_track}.
#' @export
track_marker_labels <- function(x, meta = NULL, value_col = NULL,
                                match_bp = .HS_ANNOT_DEFAULT_BP,
                                label_top_n = 10L, min_gap_rad = 0.14,
                                name = "loci", palette = NULL,
                                height_mm = NULL, window_bp = NULL,
                                window_name = "search window") {
  dt <- as.data.table(x)
  located <- "Chrom" %in% names(dt) &&
    any(c("Start", "End", "Pos") %in% names(dt))
  if (!located || !"Name" %in% names(dt)) {
    std <- as.data.table(standardize_marker_annotations(x, meta = meta,
                                                        match_bp = match_bp))
    # one row per input row, in order: carry region, window and rank over
    carry <- setdiff(intersect(c("Start", "End", "Window_start", "Window_end",
                                 value_col), names(dt)), names(std))
    for (cc in carry) std[, (cc) := dt[[cc]]]
    dt <- std
  }
  keep <- intersect(c("Chrom", "Start", "End", "Pos", "Window_start",
                      "Window_end", "Name", value_col), names(dt))
  circos_track(dt[, keep, with = FALSE], type = "labels",
               value_col = if (!is.null(value_col) &&
                               value_col %in% names(dt)) value_col else NULL,
               name = name, palette = palette %||% "grey20",
               height_mm = height_mm, label_top_n = label_top_n,
               min_gap_rad = min_gap_rad, window_bp = window_bp,
               window_name = window_name)
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


# --- the legend: measured in mm, in the hole or outside the circle -------------

# Helvetica advance widths (1/1000 em) for ASCII 32–126, from R's own
# grDevices/afm/Helvetica.afm.gz; Arial (the cairo / ragg "sans") is
# metric-compatible. They size the legend box in mm before any device exists.
.hs_helv_afm <- c(
  278L, 278L, 355L, 556L, 556L, 889L, 667L, 222L, 333L, 333L, 389L, 584L,
  278L, 333L, 278L, 278L, 556L, 556L, 556L, 556L, 556L, 556L, 556L, 556L,
  556L, 556L, 278L, 278L, 584L, 584L, 584L, 556L, 1015L, 667L, 667L, 722L,
  722L, 667L, 611L, 778L, 722L, 278L, 500L, 667L, 556L, 833L, 722L, 778L,
  667L, 778L, 722L, 667L, 611L, 722L, 667L, 944L, 667L, 667L, 611L, 278L,
  278L, 278L, 469L, 556L, 222L, 556L, 556L, 500L, 556L, 556L, 278L, 556L,
  556L, 222L, 222L, 500L, 222L, 833L, 556L, 556L, 556L, 556L, 333L, 500L,
  278L, 556L, 500L, 722L, 500L, 500L, 500L, 334L, 260L, 334L, 584L)

# Printed width (mm) of each string at `pt` points; a character outside
# ASCII counts as a digit (556).
.hs_text_mm <- function(x, pt = 5) {
  vapply(as.character(x), function(s) {
    if (is.na(s)) return(0)
    cp <- utf8ToInt(s)
    w <- rep(556L, length(cp))
    ok <- cp >= 32L & cp <= 126L
    w[ok] <- .hs_helv_afm[cp[ok] - 31L]
    sum(w) / 1000 * pt * 25.4 / 72
  }, numeric(1L), USE.NAMES = FALSE)
}

# Legend geometry, mm at print size (text is 5 pt). A colourbar entry: title
# centred 0.9 mm under its top, a 10 x 1.6 mm bar from 2.3 mm, min / max
# beside the bar — 3.9 mm tall. A chip entry (point, line, or a square: a
# categorical ring's level, a label ring's region or search window): chip and
# name centred 1.3 mm under its top — 2.6 mm tall. `reach` is how far
# the tick labels stand outside the chr bar, `air` the clearance kept from
# them, `margin` the gap to the page edge.
.hs_lg <- list(bar_w = 10, bar_h = 1.6, bar_top = 2.3, title_y = 0.9,
               bar_item_h = 3.9, chip_y = 1.3, chip_item_h = 2.6,
               lab_gap = 0.8, chip_name = 2, row_gap = 1.4, col_gap = 3.5,
               margin = 1.5, reach = 5, air = 1.5)

# Lay the entries out (mm). "stack" = one column (centre and corners), every
# bar and chip on one shared gutter; "flow" = rows across the page (bottom),
# wrapping at avail_mm, each row centred, a chip in a row with a colourbar
# sitting on the bar's centre line. Returns the box size and each entry's
# top-left offset and gutter, in mm from the box's top-left.
.hs_legend_layout <- function(items, mode = c("stack", "flow"),
                              avail_mm = Inf) {
  mode <- match.arg(mode)
  g <- .hs_lg
  n <- length(items)
  is_bar <- vapply(items, function(it) it$kind == "bar", logical(1L))
  own <- vapply(items, function(it)
    if (it$kind == "bar") it$wl + g$lab_gap else 1, numeric(1L))
  gutter <- if (mode == "stack") rep(max(own), n) else own
  width <- vapply(seq_len(n), function(i) {
    it <- items[[i]]
    if (it$kind == "bar") gutter[i] + max(g$bar_w + g$lab_gap + it$wr, it$wt)
    else gutter[i] + g$chip_name + it$wt
  }, numeric(1L))
  height <- ifelse(is_bar, g$bar_item_h, g$chip_item_h)

  if (mode == "stack")
    return(list(w = max(width), h = sum(height) + g$row_gap * (n - 1L),
                x = rep(0, n), y = cumsum(c(0, height[-n] + g$row_gap)),
                gutter = gutter))

  row <- integer(n); r <- 1L; run <- 0
  for (i in seq_len(n)) {
    if (run > 0 && run + g$col_gap + width[i] > avail_mm) {
      r <- r + 1L; run <- 0
    }
    run <- run + (if (run > 0) g$col_gap else 0) + width[i]
    row[i] <- r
  }
  has_bar <- vapply(seq_len(r), function(k) any(is_bar[row == k]),
                    logical(1L))
  off <- ifelse(!is_bar & has_bar[row],
                g$bar_top + g$bar_h / 2 - g$chip_y, 0)
  row_w <- vapply(seq_len(r), function(k)
    sum(width[row == k]) + g$col_gap * (sum(row == k) - 1L), numeric(1L))
  row_h <- vapply(seq_len(r), function(k) max((off + height)[row == k]),
                  numeric(1L))
  box_w <- max(row_w)
  x <- numeric(n)
  for (k in seq_len(r)) {
    idx <- which(row == k)
    x[idx] <- (box_w - row_w[k]) / 2 +
      cumsum(c(0, width[idx] + g$col_gap))[seq_along(idx)]
  }
  list(w = box_w, h = sum(row_h) + g$row_gap * (r - 1L), x = x,
       y = cumsum(c(0, row_h[-r] + g$row_gap))[row] + off, gutter = gutter)
}

# The canvas for a legend position at a declared width. circlize plots with
# asp = 1 and R's 4 % axis padding per side, so the page shows a half-width of
# 1.06 * 1.08 canvas units across width_mm. An outside legend extends the
# canvas down (bottom*) or up (top*) — never sideways, so the declared width
# and the circle's size hold — only as far as the box needs to clear the
# circle (chr bar + tick labels + air); a corner box uses the free corner
# first. Returns the canvas ylim and the page height (mm) that keeps 1:1.
.hs_circos_frame <- function(position, box_w, box_h, width_mm) {
  g <- .hs_lg
  V  <- 1.06 * 1.08
  k  <- 2 * V / width_mm                         # canvas units per mm
  rc <- 1 + (g$reach + g$air) * k
  w  <- box_w * k; h <- box_h * k; m <- g$margin * k
  # |y| the box's inner edge must reach, given its inner |x|
  clear_y <- function(x_in)
    if (x_in >= rc) -Inf else if (x_in <= 0) rc else sqrt(rc^2 - x_in^2)
  lo <- -V; hi <- V
  if (position == "bottom") lo <- min(lo, -rc - h - m)
  if (position %in% c("bottomright", "bottomleft"))
    lo <- min(lo, -clear_y(V - m - w) - h - m)
  if (position %in% c("topright", "topleft"))
    hi <- max(hi, clear_y(V - m - w) + h + m)
  span <- (hi - lo) / 1.08
  list(ylim = c(lo + 0.04 * span, hi - 0.04 * span),
       height_mm = width_mm * (hi - lo) / (2 * V))
}

# Top-left corner (canvas units) of the legend box on the drawn plot: `k`
# canvas units per mm on the live device, `usr` its visible extent.
.hs_legend_anchor <- function(position, box_w, box_h, k, usr) {
  g <- .hs_lg
  w <- box_w * k; h <- box_h * k; m <- g$margin * k
  switch(position,
         center      = c(-w / 2, h / 2),
         bottom      = c(-w / 2, -(1 + (g$reach + g$air) * k)),
         bottomright = c(usr[2L] - m - w, usr[3L] + m + h),
         bottomleft  = c(usr[1L] + m, usr[3L] + m + h),
         topright    = c(usr[2L] - m - w, usr[4L] - m),
         topleft     = c(usr[1L] + m, usr[4L] - m))
}

# Draw the laid-out entries with the box's top-left at (left, top).
.hs_draw_legend <- function(items, lay, left, top, k) {
  g <- .hs_lg
  cex <- .pt_cex(5)
  for (i in seq_along(items)) {
    it <- items[[i]]
    y0 <- top - lay$y[i] * k
    gx <- left + (lay$x[i] + lay$gutter[i]) * k
    if (it$kind == "bar") {
      br <- gx + g$bar_w * k
      yb <- y0 - (g$bar_top + g$bar_h / 2) * k
      hb <- g$bar_h / 2 * k
      graphics::text(gx, y0 - g$title_y * k, it$name, adj = c(0, 0.5),
                     cex = cex, col = "black")
      xs <- seq(gx, br, length.out = 25L)
      graphics::rect(xs[-25L], yb - hb, xs[-1L], yb + hb,
                     col = it$fun(seq(it$lab[1L], it$lab[2L],
                                      length.out = 24L)),
                     border = NA)
      graphics::text(gx - g$lab_gap * k, yb, it$lo, adj = c(1, 0.5),
                     cex = cex, col = "black")
      graphics::text(br + g$lab_gap * k, yb, it$hi, adj = c(0, 0.5),
                     cex = cex, col = "black")
    } else {
      yc <- y0 - g$chip_y * k
      if (it$mark == "points")
        graphics::points(gx, yc, pch = 16, cex = .pt_cex(4), col = it$col)
      else if (it$mark == "square")       # a category, a region, a window
        graphics::rect(gx - k, yc - 0.8 * k, gx + k, yc + 0.8 * k,
                       col = it$col, border = NA)
      else
        graphics::segments(gx - k, yc, gx + k, yc, col = it$col,
                           lwd = .lwd_pt(0.75))
      graphics::text(gx + g$chip_name * k, yc, it$name, adj = c(0, 0.5),
                     cex = cex, col = "black")
    }
  }
}

#' Circos meta plot — a ranked track list inside a chr–pos bar outline
#'
#' The polar analogue of \code{\link{plot_meta_alignment}} (Q25.n): the caller
#' supplies an ordered list of \code{\link{circos_track}} specs and this
#' function owns the shared frame. The outermost ring is always the
#' \strong{chr–pos bar} — one sector per chromosome in natural order,
#' Okabe–Ito colours, Mb ticks — and \code{tracks[[1]]} is the outermost data
#' ring, the last element the innermost. A legend names each ring (with a
#' mini colourbar for \code{"rect"} rings, one chip per category for a
#' categorical one, and for a \code{"labels"} ring one chip for the regions
#' plus one for the search windows when any drawn locus has one) — in the
#' centre by default, or outside the circle via \code{legend_position}.
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
#' @param legend Draw the legend (default \code{TRUE}); \code{FALSE} is
#'   \code{legend_position = "none"}.
#' @param draw  Render now on the current device (default \code{TRUE}).
#' @param legend_position Where the legend goes: \code{"center"} (default;
#'   \code{"centre"} also accepted) inside the hole, skipped with a message
#'   when the hole is too small for it; \code{"bottomright"},
#'   \code{"bottomleft"}, \code{"topright"} or \code{"topleft"} — one column
#'   tucked into that corner, outside the circle; \code{"bottom"} — the
#'   entries in rows under the circle; or \code{"none"}. An outside legend
#'   keeps the circle at its declared size and makes the figure taller by
#'   only what the legend needs — often nothing for a corner at double
#'   column — never wider; \code{\link{save_journal_figure}} takes the page
#'   height from \code{$height_mm()}.
#'
#' @return Invisibly, an \code{hs_meta_circos} object: \code{$render(width_mm)}
#'   re-draws the figure at any declared width (what
#'   \code{\link{save_journal_figure}} calls), \code{$height_mm(width_mm)}
#'   gives the figure height at that width (the width itself unless an
#'   outside legend needs room), \code{$legend} the legend entry names in
#'   order, plus the frame parameters.
#'
#' @seealso \code{\link{circos_track}}, \code{\link{track_annotation}},
#'   \code{\link{track_block_effect}}, \code{\link{track_block_variance}},
#'   \code{\link{track_marker_stat}}, \code{\link{track_marker_labels}},
#'   \code{\link{save_journal_figure}}, \code{\link{plot_block_effect_circos}}
#' @export
plot_meta_circos <- function(tracks, meta = NULL, chr_sizes = NULL,
                             rel_heights = NULL, width_mm = 89,
                             gap_degree = 2, start_degree = 90,
                             chr_colors = NULL, legend = TRUE, draw = TRUE,
                             legend_position = c("center", "bottomright",
                                                 "bottomleft", "topright",
                                                 "topleft", "bottom",
                                                 "none")) {
  if (identical(legend_position, "centre")) legend_position <- "center"
  legend_position <- match.arg(legend_position)
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
      if ("Window_end" %in% names(d))
        hi <- pmax(hi, d$Window_end, na.rm = TRUE)
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

  # --- label rings: which loci are drawn (width-independent) -----------------
  # theta per locus (region midpoint), deconflicted globally; loci with no
  # position, or on no drawn chromosome, are dropped first.
  label_kept <- lapply(tracks, function(tr) {
    if (tr$type != "labels") return(NULL)
    d <- tr$data
    ci <- match(d$Chrom, chr_levels)
    off <- c(0, cumsum(sizes))[seq_along(sizes)]
    theta <- ((off[ci] + d$Pos) / sum(sizes) * (360 - gap_degree *
                                                  length(chr_levels)) +
                (ci - 1L) * gap_degree) * pi / 180
    ok <- which(is.finite(theta))
    ok[.hs_deconflict(theta[ok], d$.val[ok], tr$label_top_n,
                      tr$min_gap_rad)]
  })

  # --- one ring's palette ----------------------------------------------------
  ring_fill <- function(tr) {
    v <- tr$data$.val
    pal <- tr$palette
    if (tr$type != "rect" || !is.null(tr$categories))
      return(list(fun = NULL, col = pal[1L]))
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

  # --- the legend: one entry per ring, measured in mm (width-independent) ---
  use_legend <- isTRUE(legend) && legend_position != "none"
  legend_items <- unlist(lapply(seq_along(tracks), function(i) {
    tr <- tracks[[i]]
    if (!is.null(tr$categories))                  # one chip per category
      return(lapply(tr$categories, function(l)
        list(kind = "chip", name = l, mark = "square",
             col = unname(tr$palette[[l]]), wt = .hs_text_mm(l))))
    if (tr$type == "labels") {                    # region (+ search window)
      lc <- .hs_label_cols(tr)
      dk <- tr$data[label_kept[[i]]]
      win <- any(is.finite(dk$Window_start) & is.finite(dk$Window_end))
      return(c(list(list(kind = "chip", name = tr$name, mark = "square",
                         col = lc[1L], wt = .hs_text_mm(tr$name))),
               if (win) list(list(kind = "chip", name = tr$window_name,
                                  mark = "square", col = lc[2L],
                                  wt = .hs_text_mm(tr$window_name)))))
    }
    list(if (tr$type == "rect") {
      rf  <- ring_fill(tr)
      lab <- c(rf$lab[1L], rf$lab[length(rf$lab)])
      lo  <- as.character(signif(lab[1L], 2))
      hi  <- as.character(signif(lab[2L], 2))
      list(kind = "bar", name = tr$name, fun = rf$fun, lab = lab,
           lo = lo, hi = hi, wt = .hs_text_mm(tr$name),
           wl = .hs_text_mm(lo), wr = .hs_text_mm(hi))
    } else {
      list(kind = "chip", name = tr$name,
           mark = if (tr$type == "points") "points" else "line",
           col = (tr$palette %||% "#0072B2")[1L],
           wt = .hs_text_mm(tr$name))
    })
  }), recursive = FALSE)
  layout_at <- function(w_mm)
    .hs_legend_layout(legend_items,
                      if (legend_position == "bottom") "flow" else "stack",
                      avail_mm = w_mm - 2 * .hs_lg$margin)
  outside <- use_legend && legend_position != "center"
  frame_at <- function(w_mm) {
    if (!outside) return(list(ylim = c(-1.06, 1.06), height_mm = w_mm))
    lay <- layout_at(w_mm)
    .hs_circos_frame(legend_position, lay$w, lay$h, w_mm)
  }

  render <- function(width_mm_now = width_mm) {
    small <- heights_mm[heights_mm < 2.5]
    if (length(small) > 0L)
      message("plot_meta_circos: ", length(small), " ring(s) under 2.5 mm at ",
              width_mm_now, " mm — consider double column (183 mm) or fewer ",
              "rings for journal legibility.")

    # Radial budget: ring heights are physical mm, so the centre hole shrinks
    # as rings stack. A cramped hole makes the innermost ring (and a centre
    # legend) collide — say so, and drop a centre legend rather than overprint
    # it, pointing at the outside positions.
    radius_mm <- width_mm_now / (2 * 1.06)
    used_mm   <- 2.2 + 2.5 + sum(heights_mm) + 0.8 * (length(tracks) + 1L)
    hole_mm   <- radius_mm - used_mm
    if (hole_mm < 0.35 * radius_mm)
      message("plot_meta_circos: ", length(tracks), " rings (",
              round(sum(heights_mm), 1), " mm) leave a ",
              max(0, round(hole_mm)), " mm centre at ", width_mm_now,
              " mm width — use double column (183 mm) or fewer rings.")
    lay <- if (use_legend) layout_at(width_mm_now)
    do_legend <- use_legend
    if (use_legend && !outside) {
      # The box's half-diagonal (+ air) must fit inside the hole, less what
      # an innermost label ring's text reaches into it.
      k0 <- 2 * 1.06 * 1.08 / width_mm_now
      hole_u <- 1 - .mm_frac(used_mm, width_mm_now)
      last <- tracks[[length(tracks)]]
      if (last$type == "labels")
        hole_u <- min(hole_u, hole_u +
                        .hs_lab_geom$text_y *
                          .mm_frac(heights_mm[length(tracks)],
                                   width_mm_now) -
                        max(.hs_text_mm(last$data$Name), 0) * k0)
      do_legend <- (sqrt(lay$w^2 + lay$h^2) / 2 + .hs_lg$air) * k0 <= hole_u
      if (!do_legend)
        message("plot_meta_circos: the centre is too small for the legend ",
                "at ", width_mm_now, " mm — legend skipped; ",
                "legend_position = \"bottomright\" (or another outside ",
                "position) draws it outside the circle.")
    }
    fr <- frame_at(width_mm_now)

    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)
    circlize::circos.par(cell.padding = c(0, 0, 0, 0),
                         gap.degree = gap_degree,
                         start.degree = start_degree,
                         track.margin = c(.mm_frac(0.4, width_mm_now),
                                          .mm_frac(0.4, width_mm_now)),
                         canvas.xlim = c(-1.06, 1.06),
                         canvas.ylim = fr$ylim,
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
            if (is.null(tr$categories)) {
              fill <- rf$fun(sd$.val)
              o <- seq_len(nrow(sd))
            } else {                     # later levels draw on top
              fill <- unname(tr$palette[sd$.cat])
              o <- order(match(sd$.cat, tr$categories))
            }
            circlize::circos.rect(s[o], rep(0, nrow(sd)), e[o],
                                  rep(1, nrow(sd)), col = fill[o],
                                  border = NA)
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
        # The deconflicted loci, drawn per sector: search windows first, the
        # regions over them, then the text (windows carry none).
        dk <- d[label_kept[[k]], ]
        lc <- .hs_label_cols(tr)
        g  <- .hs_lab_geom
        circlize::circos.track(
          ylim = c(0, 1), track.height = th, bg.border = NA,
          panel.fun = function(x, y) {
            ch <- circlize::CELL_META$sector.index
            sd <- dk[dk$Chrom == ch, ]
            if (nrow(sd) == 0L) return()
            xl <- circlize::CELL_META$xlim
            n  <- nrow(sd)
            reg <- .hs_widen(sd$Start, sd$End,
                             circlize::mm_x(g$reg_mm, h = mean(g$reg)), xl)
            w <- which(is.finite(sd$Window_start) & is.finite(sd$Window_end))
            if (length(w) > 0L) {
              halo <- circlize::mm_x(g$halo_mm, h = mean(g$win))
              win <- .hs_widen(sd$Window_start[w], sd$Window_end[w],
                               reg$e[w] - reg$s[w] + 2 * halo, xl)
              circlize::circos.rect(pmin(win$s, reg$s[w]),
                                    rep(g$win[1L], length(w)),
                                    pmax(win$e, reg$e[w]),
                                    rep(g$win[2L], length(w)),
                                    col = lc[2L], border = NA)
            }
            circlize::circos.rect(reg$s, rep(g$reg[1L], n), reg$e,
                                  rep(g$reg[2L], n), col = lc[1L],
                                  border = NA)
            circlize::circos.text(sd$Pos, g$text_y, sd$Name,
                                  cex = .pt_cex(5), col = "black",
                                  facing = "clockwise", niceFacing = TRUE,
                                  adj = c(1, 0.5))
          })
      }
    }

    # Legend: one entry per ring, in ink (5 pt) — a mini colourbar with
    # min / max beside it for rect rings, a chip + name for the rest. Drawn
    # in true mm off the live device (R's axis padding makes .mm_frac's mm
    # ~8 % short of the page, which put the old centre legend's titles onto
    # their bars), anchored by legend_position.
    if (do_legend) {
      usr <- graphics::par("usr")
      k <- diff(usr[1:2]) / (graphics::par("pin")[1L] * 25.4)
      a <- .hs_legend_anchor(legend_position, lay$w, lay$h, k, usr)
      .hs_draw_legend(legend_items, lay, a[1L], a[2L], k)
    }
    invisible(NULL)
  }

  out <- structure(list(render = render, width_mm = width_mm,
                        height_mm = function(width_mm)
                          frame_at(width_mm)$height_mm,
                        legend_position = if (use_legend) legend_position
                                          else "none",
                        tracks = track_names,
                        legend = vapply(legend_items, `[[`, character(1L),
                                        "name"),
                        chr_levels = chr_levels),
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
#' @param height Height in mm; \code{NULL} = the object's own height at that
#'   width for \code{hs_meta_circos} (\code{$height_mm()}: square unless an
#'   outside legend adds room), else 0.75 × width.
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
  h_mm <- height %||% if (!is_circos) w_mm * 0.75
                      else if (is.function(plot$height_mm)) plot$height_mm(w_mm)
                      else w_mm                     # an object from <= 0.8.1

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
