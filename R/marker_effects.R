# =============================================================================
# marker_effects.R
#
# .hs_safe_file_part()     — internal: sanitise a string for use in a file name
# .marker_locus_table()    — internal: Marker -> Chrom/Pos/Block lookup
# .effect_trait_list()     — internal: resolve geno_meta$marker_effects to traits
# export_marker_effects()  — exported: marker effects as a table / per-trait files
# plot_marker_effects()    — exported: genomic scatter of marker effects
#
# These functions only READ a GenoMeta's stored marker effects; none of them
# modify a GenoMeta or alter the behaviour of the catalogue / GEBV layer.
# =============================================================================

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Sanitise a string so it is safe to embed in a file name on every platform.
# Runs of characters outside [A-Za-z0-9._-] collapse to a single underscore.
# A message is emitted when the result differs from the input so a user whose
# trait is called e.g. "Yield (t/ha)" can see which file was actually written.
# @keywords internal
.hs_safe_file_part <- function(x) {
  orig <- as.character(x)[1L]
  safe <- gsub("[^A-Za-z0-9._-]+", "_", orig)
  safe <- gsub("_+", "_", safe)
  safe <- gsub("^[_.]+|_+$", "", safe)
  if (!nzchar(safe))
    stop("cannot build a file name from '", orig,
         "': no usable characters remain after sanitising.")
  if (!identical(safe, orig))
    message("Trait '", orig, "' sanitised to '", safe, "' for the file name.")
  safe
}


# Marker -> genomic locus + block lookup, built from the GenoMeta's map tables.
#
# Deliberately NOT geno_meta_views(): that helper rebuilds the entire wide allele
# table (an O(n_samples x n_markers) melt over the genotype array) to expose the
# same three fields. Marker effects need only the map, so reading map_info and
# block_marker_map_info directly keeps this lookup linear in the marker count.
# @keywords internal
.marker_locus_table <- function(geno_meta) {
  mi <- data.table::as.data.table(geno_meta$map_info)
  loc <- mi[, .(Marker = as.character(Marker),
                Chrom  = as.character(Chrom),
                Pos    = as.numeric(Pos))]

  bmm <- geno_meta$block_marker_map_info
  if (is.null(bmm)) {
    loc[, Block := NA_character_]
  } else {
    bt <- data.table::as.data.table(bmm)
    bt <- unique(bt[, .(Marker = as.character(Marker),
                        Block  = as.character(Block))], by = "Marker")
    loc[, Block := NA_character_]
    loc[bt, on = "Marker", Block := i.Block]
  }

  data.table::setkeyv(loc, "Marker")
  loc[]
}


# Resolve a GenoMeta's marker_effects to a named list of per-trait tables.
#
# Normalisation is delegated to .as_effect_list() (R/format_helper.R), so the
# canonical named list, a bare data.frame, and a data.frame carrying a `Trait`
# column are all accepted exactly as new_geno_meta() accepts them.
#
#   traits  character vector selecting traits, or NULL for all (order preserved)
#   single  when TRUE, require exactly one trait and return
#           list(name = <trait>, table = <data.frame>); the NULL-trait
#           disambiguation matches .resolve_effect_df() (R/catalogue_gebv.R).
# @keywords internal
.effect_trait_list <- function(geno_meta, traits = NULL, caller,
                               single = FALSE) {
  me <- geno_meta$marker_effects
  if (is.null(me))
    stop(caller, ": no marker effects on this GenoMeta — attach ",
         "'marker_effects' (e.g. via fit_genomic_prediction()) before exporting.")

  lst <- .as_effect_list(me)
  if (length(lst) == 0L)
    stop(caller, ": 'marker_effects' is present but empty.")

  if (!is.null(traits)) {
    want <- as.character(traits)
    bad  <- setdiff(want, names(lst))
    if (length(bad) > 0L)
      stop(caller, ": trait(s) not in geno_meta$marker_effects: ",
           paste(bad, collapse = ", "), ". Available: ",
           paste(names(lst), collapse = ", "), ".")
    lst <- lst[want]
  }

  if (!single) return(lst)

  if (length(lst) > 1L)
    stop(caller, ": geno_meta$marker_effects holds multiple traits (",
         paste(names(lst), collapse = ", "), "); set 'trait' to choose one.")
  list(name = names(lst)[1L], table = lst[[1L]])
}


# ---------------------------------------------------------------------------
# Exported: export_marker_effects
# ---------------------------------------------------------------------------

#' Export a GenoMeta's marker effects as a table, optionally to file
#'
#' Reads \code{geno_meta$marker_effects} — the per-trait allele-effect tables a
#' genomic-prediction fit attaches (see \code{\link{fit_genomic_prediction}}) —
#' and returns them as one tidy table, optionally writing one tab-separated file
#' per trait. The written files carry \strong{exactly the columns stored on the
#' GenoMeta} (\code{Marker}, \code{Allele}, \code{Effect}); the trait is
#' identified by the file name rather than by an extra column, so a written file
#' round-trips through \code{\link{read_effect_file}} unchanged.
#'
#' This function only reads the GenoMeta — it never modifies it.
#'
#' @param geno_meta A valid \code{GenoMeta} carrying \code{marker_effects}.
#' @param traits Optional character vector of trait names selecting which
#'   per-trait tables to export, in the given order. Default \code{NULL} (all
#'   traits). Unlike \code{\link{plot_marker_effects}}, several traits are not
#'   ambiguous here — each becomes its own file.
#' @param output_dir Optional directory. When non-\code{NULL}, \strong{one file
#'   per trait} is written there (created recursively when absent). Mutually
#'   exclusive with \code{output_file}. Default \code{NULL}.
#' @param prefix File-name stem used with \code{output_dir}; the sanitised trait
#'   name and \code{".txt"} are appended. Default \code{"marker_effects"}.
#' @param output_file Optional single path, legal only when exactly one trait is
#'   selected. Mutually exclusive with \code{output_dir}. Default \code{NULL}.
#' @param annotate Logical; when \code{TRUE} the \emph{returned} table also
#'   carries \code{Chrom}, \code{Pos} and \code{Block} joined from
#'   \code{map_info} / \code{block_marker_map_info}, sorted into genomic order.
#'   Written files are unaffected — they always keep the GenoMeta's own columns.
#'   Default \code{FALSE}.
#'
#' @return A \code{data.table}, returned whether or not anything was written,
#'   with the paths written attached as \code{attr(x, "files")}
#'   (\code{character(0)} when nothing was written). Columns:
#' \describe{
#'   \item{\code{Trait}}{Trait the effect table was stored under.}
#'   \item{\code{Chrom}, \code{Pos}, \code{Block}}{Locus and block membership;
#'     present only when \code{annotate = TRUE}, \code{NA} for a marker absent
#'     from \code{map_info}.}
#'   \item{\code{Marker}}{Marker ID, as stored.}
#'   \item{\code{Allele}}{Allele state the effect applies to, as stored.}
#'   \item{\code{Effect}}{Allele effect, as stored.}
#' }
#'
#' @examples
#' \dontrun{
#' meta <- load_geno_meta("my_geno_meta.rds")
#'
#' # Return only
#' eff <- export_marker_effects(meta)
#'
#' # One file per trait: marker_effects_Yield.txt, marker_effects_Height.txt
#' eff <- export_marker_effects(meta, output_dir = "effects")
#' attr(eff, "files")
#'
#' # A single named trait to an exact path, with locus annotation returned
#' export_marker_effects(meta, traits = "Yield",
#'                       output_file = "yield_effects.txt", annotate = TRUE)
#' }
#'
#' @seealso \code{\link{plot_marker_effects}},
#'   \code{\link{fit_genomic_prediction}}, \code{\link{compute_marker_stats}},
#'   \code{\link{read_effect_file}}
#' @export
export_marker_effects <- function(geno_meta,
                                  traits      = NULL,
                                  output_dir  = NULL,
                                  prefix      = "marker_effects",
                                  output_file = NULL,
                                  annotate    = FALSE) {

  if (!is_geno_meta(geno_meta))
    stop("export_marker_effects: 'geno_meta' is not a valid GenoMeta.")
  if (!is.null(output_dir) && !is.null(output_file))
    stop("export_marker_effects: 'output_dir' and 'output_file' are mutually ",
         "exclusive — 'output_dir' writes one file per trait, 'output_file' ",
         "writes a single named file.")

  eff_list <- .effect_trait_list(geno_meta, traits, "export_marker_effects")

  out <- data.table::rbindlist(lapply(eff_list, data.table::as.data.table),
                               idcol = "Trait", use.names = TRUE)

  if (annotate) {
    loc <- .marker_locus_table(geno_meta)
    out[loc, on = "Marker",
        `:=`(Chrom = i.Chrom, Pos = i.Pos, Block = i.Block)]
    .lead_identifier_cols(out)
    data.table::setorder(out, Trait, Chrom, Pos, Marker, Allele, na.last = TRUE)
  }

  # Written files keep the GenoMeta's own columns; the trait lives in the name.
  stored_cols <- c("Marker", "Allele", "Effect")
  written     <- character(0L)

  if (!is.null(output_file)) {
    if (length(eff_list) > 1L)
      stop("export_marker_effects: 'output_file' writes a single file but ",
           length(eff_list), " traits were selected (",
           paste(names(eff_list), collapse = ", "),
           "). Use 'output_dir' (one file per trait) or restrict 'traits'.")
    data.table::fwrite(out[, stored_cols, with = FALSE], output_file, sep = "\t")
    message("Marker effects (", names(eff_list)[1L], ") written to: ",
            output_file)
    written <- output_file

  } else if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    for (tr in names(eff_list)) {
      path <- file.path(output_dir,
                        paste0(prefix, "_", .hs_safe_file_part(tr), ".txt"))
      data.table::fwrite(out[Trait == tr, stored_cols, with = FALSE], path,
                         sep = "\t")
      message("Marker effects (", tr, ") written to: ", path)
      written <- c(written, path)
    }
  }

  out <- out[]
  data.table::setattr(out, "files", written)
  out
}


# ---------------------------------------------------------------------------
# Exported: plot_marker_effects
# ---------------------------------------------------------------------------

#' Genomic scatter of marker effects
#'
#' Draws one point per marker at its genomic position, y = the allele effect,
#' coloured by chromosome with a dashed zero line — the marker-level counterpart
#' of the haplotype catalogue's block-effect scatter. It shares the catalogue's
#' internal plotting body, so the axis, palette and theme match
#' \code{\link{plot_manhattan_effect}} and \code{\link{plot_manhattan_variance}}
#' exactly.
#'
#' @section Multiple alleles per marker:
#' Effect tables hold one row per \code{(Marker, Allele)}; a
#' \code{\link{fit_genomic_prediction}} fit typically stores the reference allele
#' at \code{0} and the alternate at the estimated effect. Plotting every row
#' therefore draws a dense, meaningless line of zeros beneath the real signal, so
#' \code{collapse = "max_abs"} (the default) keeps one row per marker — the
#' largest effect in absolute value, \strong{sign preserved} — and
#' \code{collapse = "none"} draws every stored row. Ties resolve on the
#' lexicographically first allele.
#'
#' @param geno_meta A valid \code{GenoMeta} carrying \code{marker_effects}.
#' @param trait Optional trait name selecting one stored effect table. When
#'   \code{NULL} (default) and exactly one trait is stored, that one is used;
#'   with several stored traits an error asks which.
#' @param collapse How multiple alleles per marker are reduced.
#'   \code{"max_abs"} (default) keeps the largest \code{|Effect|} per marker;
#'   \code{"none"} plots every \code{(Marker, Allele)} row.
#' @param point_size Point size. Default \code{1}.
#' @param alpha Point transparency. Default \code{0.8}.
#' @param title Optional plot title, replacing the default
#'   \code{"Marker Effect — <trait>"}.
#' @param table_file Optional path; when non-\code{NULL} the plotted table is
#'   written there as a tab-separated file. Default \code{NULL}.
#'
#' @return A \code{ggplot} object. The plotted table is attached as
#'   \code{attr(p, "effect_table")}, a \code{data.table} with columns
#'   \code{Trait}, \code{Chrom}, \code{Pos}, \code{Block}, \code{Marker},
#'   \code{Allele}, \code{Effect} — the same shape as
#'   \code{export_marker_effects(annotate = TRUE)}.
#'
#' @examples
#' \dontrun{
#' meta <- load_geno_meta("my_geno_meta.rds")
#' p <- plot_marker_effects(meta, trait = "Yield")
#' print(p)
#' head(attr(p, "effect_table"))
#'
#' # Every allele row, and save the backing table
#' plot_marker_effects(meta, trait = "Yield", collapse = "none",
#'                     table_file = "yield_marker_effects.txt")
#' }
#'
#' @seealso \code{\link{export_marker_effects}},
#'   \code{\link{plot_manhattan_effect}}, \code{\link{plot_marker_signal}}
#' @export
plot_marker_effects <- function(geno_meta,
                                trait      = NULL,
                                collapse   = c("max_abs", "none"),
                                point_size = 1,
                                alpha      = 0.8,
                                title      = NULL,
                                table_file = NULL) {

  collapse <- match.arg(collapse)
  if (!is_geno_meta(geno_meta))
    stop("plot_marker_effects: 'geno_meta' is not a valid GenoMeta.")

  res <- .effect_trait_list(geno_meta, trait, "plot_marker_effects",
                            single = TRUE)
  dt  <- data.table::as.data.table(res$table)
  dt[, Trait := res$name]

  loc <- .marker_locus_table(geno_meta)
  dt[loc, on = "Marker", `:=`(Chrom = i.Chrom, Pos = i.Pos, Block = i.Block)]

  off <- dt[is.na(Chrom) | is.na(Pos)]
  if (nrow(off) > 0L) {
    nm <- unique(off$Marker)
    warning("plot_marker_effects: ", nrow(off), " effect row(s) for ",
            length(nm), " marker(s) absent from map_info were dropped: ",
            paste(utils::head(nm, 5L), collapse = ", "),
            if (length(nm) > 5L) ", ..." else "")
    dt <- dt[!is.na(Chrom) & !is.na(Pos)]
  }
  if (nrow(dt) == 0L)
    stop("plot_marker_effects: no marker effects have a genomic position — ",
         "check that the effect table's Marker IDs match geno_meta$map_info$Marker.")

  data.table::setorder(dt, Chrom, Pos, Marker, Allele)
  if (collapse == "max_abs") {
    idx <- dt[, .I[which.max(abs(Effect))], by = Marker][[2L]]
    dt  <- dt[idx]
  }
  .lead_identifier_cols(dt)

  p <- .plot_manhattan(
    dt, value_col = "Effect",
    title      = title %||% paste0("Marker Effect — ", res$name),
    y_lab      = "Marker effect",
    point_size = point_size, alpha = alpha, zero_line = TRUE)

  if (!is.null(table_file)) {
    data.table::fwrite(dt, table_file, sep = "\t")
    message("Marker effect table written to: ", table_file)
  }

  attr(p, "effect_table") <- dt[]
  p
}
