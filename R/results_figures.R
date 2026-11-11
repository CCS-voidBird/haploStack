# =============================================================================
# results_figures.R — figure recipes and the one renderer (0.9.0, design D3.2 /
# framework F5)
#
# Compute functions never touch a graphics device. A figure inside a result is
# an `hs_figure`: a data-only recipe (plotter name + arguments + size + file
# stem) that survives saveRDS() and a JSON conversion, and is built on demand
# by hs_build_figure(). An argument may be an `hs_ref` — a pointer to one of
# the owning result's tables — so a catalogue result does not keep a second
# copy of the blocks × samples data its Manhattan plot is drawn from.
#
# hs_render_figure() is the only place a file device is opened for a figure:
# it accepts every kind the package produces (ggplot, hs_figure, hs_meta_circos,
# a nested hs_figure_set, a recordedplot), picks the device from the file
# extension with capability fallbacks, and restores the caller's device, par()
# and circlize state afterwards. save_journal_figure() and write_outputs() call
# it; nothing else opens a device.
# =============================================================================

# ---- recipes ----------------------------------------------------------------

# A figure recipe. `plotter` is the name of an exported haploStack function;
# `args` its arguments, where any element may be an `hs_ref` (see below).
# `size` is the legacy render size write_outputs() reproduces; `file` the
# stem (with extension) the writer uses inside the output directory.
.hs_figure <- function(name, plotter, args = list(), width_mm = 183,
                       height_mm = NULL, dpi = 150, file = NULL) {
  if (!is.character(plotter) || length(plotter) != 1L || !nzchar(plotter))
    .hs_abort(".hs_figure: 'plotter' must be the name of one exported function.",
              class = "hs_error_internal")
  if (!is.list(args))
    .hs_abort(".hs_figure: 'args' must be a list.", class = "hs_error_internal")
  structure(list(name    = as.character(name)[1L],
                 plotter = plotter,
                 args    = args,
                 size    = list(width_mm  = as.numeric(width_mm),
                                height_mm = as.numeric(height_mm %||% width_mm * 0.75),
                                dpi       = as.numeric(dpi)),
                 file    = file %||% paste0(.hs_file_stem(name), ".png")),
            class = "hs_figure")
}

# A pointer into the owning result: `what` names an element, `trait` the trait
# for a trait-keyed result (hs_catalogue). `derive` names a registered
# derivation instead (a plot table that is computed from the result's tables
# rather than stored twice). Resolved by hs_build_figure().
.hs_ref <- function(what, trait = NULL, derive = NULL) {
  if (!is.null(derive) && is.null(.HS_REF_DERIVE[[derive]]))
    .hs_abort(".hs_ref: unknown derivation '", derive, "'.", class = "hs_error_internal")
  structure(list(what = if (is.null(what)) NA_character_ else as.character(what)[1L],
                 trait = trait, derive = derive),
            class = "hs_ref")
}

# Derivations: each takes the result element the ref points at (a trait
# element for a catalogue, the result itself otherwise) and returns the table
# the plotter wants. Data only; reproduces exactly what 0.8.8 built inline.
.HS_REF_DERIVE <- list(
  # The element itself (an analysis result is what its plotters take).
  self = function(el) el,
  # plot_manhattan_effect(): the long Block x Sample effects joined with the
  # block positions (Chrom, Pos) that block_summary carries.
  manhattan_effect_data = function(el) {
    long <- .melt_gebv(el$gebv_matrix)
    pos  <- data.table::as.data.table(el$block_summary)[, c("Block", "Chrom", "Pos"), with = FALSE]
    merge(long, unique(pos), by = "Block", all.x = TRUE)
  },
  # plot_meta_circos(): one block-effect ring per trait of an analysis result.
  block_effect_tracks = function(el) {
    traits <- unique(as.character(el$S$Trait))
    lapply(traits, function(tr) track_block_effect(el, trait = tr))
  }
)

.hs_deref <- function(ref, result) {
  what  <- if (!is.null(ref$derive)) paste0("derivation '", ref$derive, "'")
           else sprintf("'%s'", ref$what)
  where <- if (is.null(ref$trait)) what else sprintf("%s of trait '%s'", what, ref$trait)
  if (is.null(result))
    .hs_abort("hs_build_figure: the recipe refers to ", where, " of its result; ",
              "pass result = <the hs_result>, or build through ",
              "hs_figures(res, build = TRUE).", call. = FALSE)
  el <- if (is.null(ref$trait)) result
        else if (is.list(result[[ref$trait]])) result[[ref$trait]] else NULL
  if (is.null(el))
    .hs_abort("hs_build_figure: the result has no trait '", ref$trait, "'.", call. = FALSE)
  if (!is.null(ref$derive)) return(.HS_REF_DERIVE[[ref$derive]](el))
  x <- el[[ref$what]]
  if (is.null(x))
    .hs_abort("hs_build_figure: the result has no element ", where, ".", call. = FALSE)
  x
}

.hs_resolve_args <- function(args, result) {
  lapply(args, function(a) if (inherits(a, "hs_ref")) .hs_deref(a, result) else a)
}

#' Build a figure recipe into a plot object
#'
#' A result's figures (\code{\link{hs_figures}}) are stored as recipes of
#' class \code{"hs_figure"}: the name of the plotter, its arguments and the
#' size the figure is written at. Nothing is drawn until a recipe is built.
#' \code{hs_build_figure()} calls the plotter with the stored arguments; where
#' an argument points into the owning result (so the result does not carry a
#' second copy of a table), \code{result} supplies it. Printing a recipe builds
#' it and draws it on the current device.
#'
#' @param fig An \code{hs_figure} recipe.
#' @param result The \code{\link[=hs_result]{hs_result}} the recipe belongs
#'   to; needed only when the recipe refers to one of its tables.
#' @return The built figure: a ggplot, or an \code{hs_meta_circos}.
#' @seealso \code{\link{hs_render_figure}}, \code{\link{hs_figures}}
#' @export
hs_build_figure <- function(fig, result = NULL) {
  if (!inherits(fig, "hs_figure"))
    .hs_abort("hs_build_figure: 'fig' must be an hs_figure recipe (see hs_figures()).")
  fn <- tryCatch(getExportedValue("haploStack", fig$plotter), error = function(e) NULL)
  if (!is.function(fn))
    .hs_abort("hs_build_figure: unknown plotter '", fig$plotter,
              "' (a recipe names an exported haploStack function).")
  do.call(fn, .hs_resolve_args(fig$args, result))
}

#' @export
print.hs_figure <- function(x, ...) {
  .hs_draw_now(hs_build_figure(x))
  invisible(x)
}

# Draw a built figure on the current device (console use: print methods).
.hs_draw_now <- function(obj) {
  if (inherits(obj, "ggplot") || inherits(obj, "hs_meta_circos")) print(obj)
  else if (inherits(obj, "recordedplot")) grDevices::replayPlot(obj)
  else graphics::plot(obj)
  invisible(obj)
}

# ---- nested figure sets -----------------------------------------------------

# plot_localgebv_alignment() / plot_favourable_blocks() return a nested list
# (generation -> family -> ggplot). The shape is kept (viz/F08); the class lets
# the renderer write one file per leaf and a writer dispatch on it.
.hs_figure_set <- function(x) {
  if (inherits(x, "hs_figure_set")) return(x)
  structure(x, class = c("hs_figure_set", "list"))
}

.hs_is_figure_like <- function(x) {
  inherits(x, "ggplot") || inherits(x, "hs_meta_circos") ||
    inherits(x, "hs_figure") || inherits(x, "recordedplot")
}

# Flatten a set to its leaves, named by the path of list names joined by "_"
# (unnamed levels contribute their index).
.hs_set_leaves <- function(x, prefix = character()) {
  if (.hs_is_figure_like(x)) {
    out <- list(x); names(out) <- paste(prefix, collapse = "_"); return(out)
  }
  if (!is.list(x)) return(list())
  nm <- names(x) %||% rep("", length(x))
  nm[!nzchar(nm)] <- as.character(seq_along(x))[!nzchar(nm)]
  out <- list()
  for (i in seq_along(x))
    out <- c(out, .hs_set_leaves(x[[i]], c(prefix, nm[i])))
  out
}

.hs_leaf_files <- function(file, leaf_names) {
  ext  <- tools::file_ext(file)
  stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", file) else file
  paste0(stem, "_", gsub("[^A-Za-z0-9._-]+", "_", leaf_names),
         if (nzchar(ext)) paste0(".", ext) else "")
}

# ---- the renderer -----------------------------------------------------------

# Open the file device chosen by the extension. ragg when installed, else the
# cairo bitmap devices; pdf through cairo_pdf (fonts embedded) when cairo is
# available, else grDevices::pdf() with a note; svg needs cairo.
.hs_open_device <- function(file, width_mm, height_mm, dpi) {
  ext   <- tolower(tools::file_ext(file))
  in_w  <- width_mm / 25.4
  in_h  <- height_mm / 25.4
  cairo <- isTRUE(unname(capabilities("cairo")))
  ragg  <- requireNamespace("ragg", quietly = TRUE)
  # Without cairo the bitmap devices take R's own default type (no option read:
  # options are read only through .hs_setting(), hygiene H1).
  bitmap <- function(dev) {
    if (cairo) dev(file, width = in_w, height = in_h, units = "in", res = dpi, type = "cairo")
    else       dev(file, width = in_w, height = in_h, units = "in", res = dpi)
  }
  switch(ext,
    png = if (ragg) ragg::agg_png(file, width = in_w, height = in_h, units = "in", res = dpi)
          else bitmap(grDevices::png),
    tiff = if (ragg) ragg::agg_tiff(file, width = in_w, height = in_h, units = "in", res = dpi)
           else bitmap(grDevices::tiff),
    pdf = if (cairo) grDevices::cairo_pdf(file, width = in_w, height = in_h, pointsize = 12)
          else {
            .hs_note("hs_render_figure", "fonts_not_embedded",
                     "cairo is not available: the PDF was written with grDevices::pdf(), fonts not embedded")
            grDevices::pdf(file, width = in_w, height = in_h, pointsize = 12)
          },
    svg = if (cairo) grDevices::svg(file, width = in_w, height = in_h, pointsize = 12)
          else .hs_abort("hs_render_figure: writing .svg needs cairo (capabilities(\"cairo\") is FALSE).",
                         class = "hs_error_dependency",
                         data = list(hint = "install an R build with cairo, or write .png / .pdf")),
    .hs_abort("hs_render_figure: unsupported file extension '.", ext,
              "' (use .png, .pdf, .svg or .tiff)."))
  invisible(grDevices::dev.cur())
}

#' Render a figure to a file
#'
#' The one place haploStack opens a file device for a figure. It takes every
#' kind of figure the package produces — a ggplot, an \code{hs_figure} recipe
#' (built first), an \code{hs_meta_circos} from \code{\link{plot_meta_circos}},
#' a nested \code{hs_figure_set} (one file per leaf, the leaf path appended to
#' the file stem) or a \code{recordedplot} (replayed as recorded) — picks the
#' device from the file extension, draws, and closes the device. The caller's
#' current device, \code{par()} and circlize state are left as they were.
#'
#' Devices: \code{.png} / \code{.tiff} use \pkg{ragg} when installed, else the
#' cairo bitmap devices; \code{.pdf} uses \code{cairo_pdf()} (fonts embedded)
#' when cairo is available, else \code{pdf()} with a \code{fonts_not_embedded}
#' note; \code{.svg} needs cairo. Any other extension is an error.
#'
#' @param x A ggplot, \code{hs_figure}, \code{hs_meta_circos},
#'   \code{hs_figure_set} or \code{recordedplot}.
#' @param file Output path; the extension picks the device. For a set, the
#'   leaf path is inserted before the extension (\code{stem_gen_family.png}).
#' @param width_mm,height_mm Page size in mm. Defaults: an
#'   \code{hs_meta_circos} uses its declared width and its own height at that
#'   width; a recipe uses its stored size; anything else 183 mm by 0.75 of the
#'   width.
#' @param dpi Raster resolution (default 300; a recipe's stored dpi when not
#'   given).
#' @param create_dir Create a missing output directory (default \code{FALSE}:
#'   a missing directory is an error that names it).
#' @param overwrite Replace an existing file (default \code{TRUE}).
#' @return Invisibly, the normalised path(s) written.
#' @seealso \code{\link{hs_build_figure}}, \code{\link{save_journal_figure}},
#'   \code{\link{write_outputs}}
#' @export
hs_render_figure <- function(x, file, width_mm = NULL, height_mm = NULL, dpi = 300,
                             create_dir = FALSE, overwrite = TRUE) {
  if (!is.character(file) || length(file) != 1L || !nzchar(file))
    .hs_abort("hs_render_figure: 'file' must be a single path.")

  if (inherits(x, "hs_figure_set") ||
      (is.list(x) && !.hs_is_figure_like(x) && !is.data.frame(x))) {
    leaves <- .hs_set_leaves(x)
    if (!length(leaves))
      .hs_abort("hs_render_figure: the figure set holds no figure.")
    files <- .hs_leaf_files(file, names(leaves))
    out <- character(length(leaves))
    for (i in seq_along(leaves))
      out[i] <- hs_render_figure(leaves[[i]], files[i], width_mm = width_mm,
                                 height_mm = height_mm, dpi = dpi,
                                 create_dir = create_dir, overwrite = overwrite)
    return(invisible(out))
  }

  if (inherits(x, "hs_figure")) {
    if (missing(dpi) && !is.null(x$size$dpi)) dpi <- x$size$dpi
    width_mm  <- width_mm  %||% x$size$width_mm
    height_mm <- height_mm %||% x$size$height_mm
    x <- hs_build_figure(x)
  }
  if (inherits(x, "hs_meta_circos")) {
    width_mm  <- width_mm %||% x$width_mm
    height_mm <- height_mm %||%
      if (is.function(x$height_mm)) x$height_mm(width_mm) else width_mm
  } else if (inherits(x, "ggplot") || inherits(x, "recordedplot")) {
    width_mm  <- width_mm  %||% 183
    height_mm <- height_mm %||% width_mm * 0.75
  } else {
    .hs_abort("hs_render_figure: cannot render an object of class ",
              paste(class(x), collapse = "/"),
              " (expected a ggplot, hs_figure, hs_meta_circos, hs_figure_set or recordedplot).")
  }
  if (!is.numeric(width_mm) || !is.numeric(height_mm) ||
      !is.finite(width_mm) || !is.finite(height_mm) || width_mm <= 0 || height_mm <= 0)
    .hs_abort("hs_render_figure: 'width_mm' and 'height_mm' must be positive numbers.")

  .hs_dir_ready(dirname(file), create_dir = create_dir)
  .hs_check_overwrite(file, overwrite)

  prev <- grDevices::dev.cur()
  .hs_open_device(file, width_mm, height_mm, dpi)
  this <- grDevices::dev.cur()
  on.exit({
    if (this %in% grDevices::dev.list()) grDevices::dev.off(this)
    if (prev > 1L && prev %in% grDevices::dev.list()) grDevices::dev.set(prev)
  }, add = TRUE)

  if (inherits(x, "hs_meta_circos")) {
    if (requireNamespace("circlize", quietly = TRUE))
      on.exit(circlize::circos.clear(), add = TRUE, after = FALSE)
    graphics::par(mar = c(0, 0, 0, 0), cex = 1)
    x$render(width_mm)
  } else if (inherits(x, "ggplot")) {
    print(x)
  } else {
    grDevices::replayPlot(x)
  }
  invisible(normalizePath(file, winslash = "/", mustWork = FALSE))
}
