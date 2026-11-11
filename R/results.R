# =============================================================================
# results.R — the result contract (0.9.0, design D3 §3.1 / framework F4, A1)
#
# Compute functions return data; they write files only when the caller names a
# destination. Every result carries three things beside its 0.8.8 elements:
#
#   notes       data.frame(stage, code, message, n): every haploStack warning
#               and silent fallback raised while the result was computed
#   provenance  list(package, result_schema, created, fun, ...)
#   figures     unrendered figure recipes (hs_figure, D3.2) or ggplots
#
# reached through the accessors hs_notes(), hs_provenance(), hs_figures(),
# which is the contract. Where they live depends on the result's shape:
# `hs_catalogue` is a list keyed by trait, so it keeps notes and provenance as
# attributes and per-trait figures inside each trait element; every other
# result (`hs_block_analysis`, `hs_simulation`) keeps them as elements.
# `.hs_result()` is the one constructor; write_outputs() is the one writer
# (an S3 generic, one method per result class, landing with each class).
#
# Results must survive saveRDS()/readRDS() and a JSON conversion of their
# tables: no closures, environments, recorded plots, external pointers.
# =============================================================================

.HS_RESULT_SCHEMA <- 1L

# ---- notes ------------------------------------------------------------------

.hs_notes_df <- function(stage = character(), code = character(),
                         message = character(), n = integer()) {
  data.frame(stage = as.character(stage), code = as.character(code),
             message = as.character(message), n = as.integer(n),
             stringsAsFactors = FALSE)
}

# Raise a silent fallback as a classed warning so that it (a) reaches the
# caller like any other warning and (b) lands as a row in the result's notes
# when the entry point runs under .hs_with_notes(). `...` are extra fields
# (D4 adds job / first_cycle).
.hs_note <- function(stage, code, message, ...) {
  .hs_warn(message, call. = FALSE, class = "hs_warning_fallback",
           data = c(list(stage = stage, code = code), list(...)))
}

# One notes row from a condition: stage = the condition's own `stage` field,
# else the function that raised it; code = its own `code`, else its first
# class.
.hs_note_row <- function(w) {
  stage <- w$stage
  if (is.null(stage)) {
    cl <- conditionCall(w)
    stage <- if (is.null(cl)) NA_character_ else as.character(cl[[1L]])[1L]
  }
  code <- w$code
  if (is.null(code)) code <- class(w)[1L]
  extra <- w[setdiff(names(w), c("message", "call", "stage", "code"))]
  extra <- extra[vapply(extra, function(v) is.atomic(v) && length(v) == 1L, TRUE)]
  row <- .hs_notes_df(stage, code, conditionMessage(w), 1L)
  for (nm in names(extra)) row[[nm]] <- extra[[nm]]
  row
}

# Collapse identical rows (same stage, code, message and extra columns) into
# one with their count in `n`.
.hs_collapse_notes <- function(notes) {
  if (!NROW(notes)) return(.hs_notes_df())
  key_cols <- setdiff(names(notes), "n")
  key <- do.call(paste, c(lapply(notes[key_cols], as.character), sep = "\r"))
  first <- !duplicated(key)
  out <- notes[first, , drop = FALSE]
  out$n <- as.integer(tapply(notes$n, factor(key, levels = unique(key)), sum))
  rownames(out) <- NULL
  out
}

# Run `expr` and collect every haploStack warning raised inside it as a notes
# row; the warnings themselves propagate unchanged (muffle = FALSE) or are
# swallowed (muffle = TRUE, for drivers that re-raise once per code later).
# Returns list(value, notes).
.hs_with_notes <- function(expr, muffle = FALSE) {
  rows <- list()
  value <- withCallingHandlers(expr, hs_warning = function(w) {
    rows[[length(rows) + 1L]] <<- .hs_note_row(w)
    if (isTRUE(muffle)) invokeRestart("muffleWarning")
  })
  notes <- if (length(rows)) {
    cols <- unique(unlist(lapply(rows, names)))
    rows <- lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] })
    .hs_collapse_notes(do.call(rbind, rows))
  } else .hs_notes_df()
  list(value = value, notes = notes)
}

# ---- constructor and accessors ---------------------------------------------

.hs_provenance <- function(fun, ...) {
  c(list(package = as.character(utils::packageVersion("haploStack")),
         result_schema = .HS_RESULT_SCHEMA,
         created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
         fun = fun),
    list(...))
}

# The one constructor. `class` is the result's own class ("hs_catalogue",
# "hs_block_analysis", "hs_simulation"); "hs_result" and "list" are appended.
# An `hs_catalogue` keeps notes and provenance as attributes (its top level is
# keyed by trait); every other result keeps them as elements.
.hs_result <- function(x, class, notes = NULL, provenance = NULL) {
  if (!is.list(x))
    .hs_abort(class = "hs_error_internal", ".hs_result: 'x' must be a list.")
  if (is.null(notes)) notes <- .hs_notes_df()
  if (is.null(provenance)) provenance <- .hs_provenance(fun = NA_character_)
  if ("hs_catalogue" %in% class) {
    attr(x, "hs_notes") <- notes
    attr(x, "hs_provenance") <- provenance
  } else {
    x$notes <- notes
    x$provenance <- provenance
    if (is.null(x$figures)) x$figures <- list()
  }
  class(x) <- unique(c(class, "hs_result", "list"))
  x
}

.hs_is_catalogue <- function(x) inherits(x, "hs_catalogue")

#' Notes, provenance and figures of a haploStack result
#'
#' Since 0.9.0 the analysis entry points return a result object (class
#' \code{"hs_result"} plus its own class: \code{"hs_catalogue"} from
#' \code{\link{run_haplotype_catalogue}}, \code{"hs_block_analysis"} from
#' \code{\link{analyze_block_catalogue}}, \code{"hs_simulation"} from
#' \code{\link{run_stacking_simulation}}). The 0.8.8 elements are unchanged
#' (\code{res$summary}, \code{res[["Yield"]]$gebv_matrix}, ... keep working);
#' beside them every result carries notes, provenance and unrendered figures,
#' reached through these accessors rather than by position, because where they
#' are stored depends on the result's shape.
#'
#' \describe{
#'   \item{\code{hs_notes(x)}}{A data.frame \code{stage, code, message, n}
#'     with one row per distinct haploStack warning or silent fallback raised
#'     while the result was computed (\code{n} counts repeats); zero rows when
#'     the run was clean. The same warnings were also signalled at the time.}
#'   \item{\code{hs_provenance(x)}}{A list: \code{package} version,
#'     \code{result_schema}, \code{created} (UTC), \code{fun}, and what the
#'     entry point records (centring, traits, file stems, label, block
#'     provenance, engine, resolved settings).}
#'   \item{\code{hs_figures(x, build = FALSE)}}{A named list of unrendered
#'     figures: figure recipes (\code{hs_figure}) or ggplots. For an
#'     \code{hs_catalogue} the list is keyed by trait. \code{build = TRUE}
#'     builds every recipe against the result (\code{\link{hs_build_figure}})
#'     and returns plot objects. Rendering to a file is
#'     \code{\link{write_outputs}}'s or \code{\link{hs_render_figure}}'s job.}
#' }
#'
#' @param x An \code{hs_result}.
#' @param build Build the recipes into plot objects (default \code{FALSE}:
#'   return the recipes).
#' @return See the description.
#' @seealso \code{\link{write_outputs}}
#' @examples
#' \dontrun{
#' res <- run_haplotype_catalogue(meta, effects)   # in memory, nothing written
#' hs_notes(res)
#' hs_provenance(res)$label
#' write_outputs(res, "catalogue_output")
#' }
#' @name hs_result
#' @export
hs_notes <- function(x) {
  .hs_check_result(x, "hs_notes")
  n <- if (.hs_is_catalogue(x)) attr(x, "hs_notes", exact = TRUE) else x$notes
  if (is.null(n)) .hs_notes_df() else n
}

#' @rdname hs_result
#' @export
hs_provenance <- function(x) {
  .hs_check_result(x, "hs_provenance")
  if (.hs_is_catalogue(x)) attr(x, "hs_provenance", exact = TRUE) else x$provenance
}

#' @rdname hs_result
#' @export
hs_figures <- function(x, build = FALSE) {
  .hs_check_result(x, "hs_figures")
  bld <- function(fig) if (inherits(fig, "hs_figure")) hs_build_figure(fig, result = x) else fig
  if (.hs_is_catalogue(x)) {
    f <- lapply(unclass(x), function(tr) if (is.list(tr)) tr$figures)
    f <- f[!vapply(f, is.null, TRUE)]
    if (isTRUE(build)) f <- lapply(f, function(figs) lapply(figs, bld))
    f
  } else {
    f <- if (is.null(x$figures)) list() else x$figures
    if (isTRUE(build)) f <- lapply(f, bld)
    f
  }
}

.hs_check_result <- function(x, fn) {
  if (!inherits(x, "hs_result"))
    .hs_abort(fn, ": 'x' is not a haploStack result (class \"hs_result\").",
              call. = FALSE)
  invisible(TRUE)
}

#' @export
print.hs_result <- function(x, ...) {
  own <- setdiff(class(x), c("hs_result", "list"))
  prov <- hs_provenance(x)
  notes <- hs_notes(x)
  cat("<haploStack result: ", paste(own, collapse = "/"), ">\n", sep = "")
  if (.hs_is_catalogue(x)) {
    cat("  traits : ", paste(names(x), collapse = ", "), "\n", sep = "")
  } else {
    cat("  fields : ", paste(setdiff(names(x), c("notes", "provenance", "figures")),
                           collapse = ", "), "\n", sep = "")
  }
  if (!is.null(prov$label)) cat("  label  : ", prov$label, "\n", sep = "")
  if (!is.null(prov$fun) && !is.na(prov$fun))
    cat("  from   : ", prov$fun, "() in haploStack ", prov$package, "\n", sep = "")
  nf <- length(hs_figures(x))
  cat("  notes  : ", NROW(notes), "   figures: ", nf, "\n", sep = "")
  if (NROW(notes)) {
    show <- utils::head(notes, 3L)
    for (i in seq_len(nrow(show)))
      cat("    [", show$code[i], "] ", show$message[i],
          if (show$n[i] > 1L) sprintf(" (x%d)", show$n[i]) else "", "\n", sep = "")
    if (nrow(notes) > 3L) cat("    ... and ", nrow(notes) - 3L, " more (hs_notes())\n", sep = "")
  }
  cat("  Use write_outputs(x, dir) to write tables and figures.\n")
  invisible(x)
}

#' @export
summary.hs_result <- function(object, ...) {
  tbl <- function(lst, prefix) {
    rows <- lapply(names(lst), function(nm) {
      v <- lst[[nm]]
      if (is.data.frame(v))
        data.frame(object = paste0(prefix, nm), class = class(v)[1L],
                   nrow = nrow(v), ncol = ncol(v), stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }
  out <- if (.hs_is_catalogue(object)) {
    do.call(rbind, lapply(names(object), function(tr)
      tbl(object[[tr]], paste0(tr, "/"))))
  } else {
    tbl(object[setdiff(names(object), c("notes", "provenance", "figures"))], "")
  }
  if (is.null(out)) out <- data.frame(object = character(), class = character(),
                                      nrow = integer(), ncol = integer(),
                                      stringsAsFactors = FALSE)
  structure(list(tables = out, notes = hs_notes(object),
                 provenance = hs_provenance(object)),
            class = "summary.hs_result")
}

#' @export
print.summary.hs_result <- function(x, ...) {
  cat("Tables:\n"); print(x$tables, row.names = FALSE)
  cat("Notes: ", NROW(x$notes), "\n", sep = "")
  invisible(x)
}

# ---- file stems and the table writer ---------------------------------------

# Sanitise a trait name or label for use inside a file name. A changed stem
# is reported as a note (hs_warning_fallback, code "file_stem_renamed") so it
# reaches the result's notes; nothing is printed.
.hs_file_stem <- function(x, what = "trait") {
  orig <- as.character(x)[1L]
  safe <- gsub("[^A-Za-z0-9._-]+", "_", orig)
  safe <- gsub("_+", "_", safe)
  safe <- gsub("^[_.]+|_+$", "", safe)
  if (!nzchar(safe))
    .hs_abort("cannot build a file name from ", what, " '", orig,
              "': no usable characters remain after sanitising.", call. = FALSE)
  if (!identical(safe, orig))
    .hs_note("file_stem", "file_stem_renamed",
             paste0(what, " '", orig, "' sanitised to '", safe, "' for the file name."))
  safe
}

# Stems for a vector of names; two names that collide after sanitising are an
# input error (their files would overwrite each other).
.hs_file_stems <- function(x, what = "trait") {
  x <- as.character(x)
  stems <- vapply(x, .hs_file_stem, character(1L), what = what, USE.NAMES = FALSE)
  if (anyDuplicated(stems)) {
    dup <- stems[duplicated(stems)]
    .hs_abort("these ", what, " names give the same file name after sanitising: ",
              paste(sprintf("'%s'", x[stems %in% dup]), collapse = ", "),
              " -> '", paste(unique(dup), collapse = "', '"), "'. Rename one of them.",
              call. = FALSE, data = list(bad = x[stems %in% dup]))
  }
  stats::setNames(stems, x)
}

# The 0.8.8 name, kept for export_marker_effects() (which still reports the
# rename as a message); new code uses .hs_file_stem().
.hs_safe_file_part <- function(x) {
  orig <- as.character(x)[1L]
  safe <- gsub("[^A-Za-z0-9._-]+", "_", orig)
  safe <- gsub("_+", "_", safe)
  safe <- gsub("^[_.]+|_+$", "", safe)
  if (!nzchar(safe))
    .hs_abort("cannot build a file name from '", orig,
         "': no usable characters remain after sanitising.")
  if (!identical(safe, orig))
    .hs_inform("Trait '", orig, "' sanitised to '", safe, "' for the file name.")
  safe
}

.hs_manifest <- function(path = character(), kind = character(),
                         object = character(), bytes = numeric()) {
  data.frame(path = as.character(path), kind = as.character(kind),
             object = as.character(object), bytes = as.numeric(bytes),
             stringsAsFactors = FALSE)
}

# Refuse to overwrite existing targets unless asked; names every file.
.hs_check_overwrite <- function(paths, overwrite) {
  if (isTRUE(overwrite)) return(invisible(TRUE))
  there <- paths[file.exists(paths)]
  if (length(there))
    .hs_abort("these output files exist already (pass overwrite = TRUE to replace them): ",
              paste(there, collapse = ", "), call. = FALSE,
              data = list(bad = there))
  invisible(TRUE)
}

# Make sure `dir` can receive files: create it when asked, else it must exist.
.hs_dir_ready <- function(dir, create_dir = TRUE) {
  if (!is.character(dir) || length(dir) != 1L || !nzchar(dir))
    .hs_abort("'dir' must be a single directory path.", call. = FALSE)
  if (!dir.exists(dir)) {
    if (!isTRUE(create_dir))
      .hs_abort("output directory does not exist: ", dir,
                " (pass create_dir = TRUE to create it).", call. = FALSE)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(normalizePath(dir, winslash = "/", mustWork = TRUE))
}

# Write one table the way 0.8.8 wrote them (fwrite, tab-separated) and return
# its manifest row. `object` names what was written ("Yield/gebv_matrix").
.hs_write_table <- function(dt, path, overwrite = TRUE, object = NA_character_,
                            sep = "\t", ...) {
  .hs_check_overwrite(path, overwrite)
  data.table::fwrite(dt, file = path, sep = sep, ...)
  .hs_manifest(normalizePath(path, winslash = "/", mustWork = TRUE), "table",
               object, file.info(path)$size)
}

# Write several tables into one directory (the exporters' writer, D5 §3.6 h):
# `tables` is a named list, the names are the file names (`prefix_` goes in
# front when given). Returns the absolute paths (named as `tables`) and the
# manifest. `...` goes to fwrite, so each exporter keeps its 0.8.8 bytes.
.hs_write_tables <- function(tables, dir, prefix = NULL, sep = "\t",
                             overwrite = TRUE, ...) {
  if (!is.list(tables) || is.null(names(tables)) || any(!nzchar(names(tables))))
    .hs_abort(".hs_write_tables: 'tables' must be a named list.", call. = FALSE)
  dir <- .hs_dir_ready(dir, create_dir = TRUE)
  pfx <- if (is.null(prefix) || !nzchar(prefix)) "" else paste0(prefix, "_")
  paths <- file.path(dir, paste0(pfx, names(tables)))
  .hs_check_overwrite(paths, overwrite)
  for (i in seq_along(tables))
    data.table::fwrite(tables[[i]], file = paths[i], sep = sep, ...)
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  names(paths) <- names(tables)
  list(paths = paths,
       manifest = .hs_manifest(unname(paths), "table", names(tables),
                               file.info(paths)$size))
}

# An ISO-8601 UTC timestamp for provenance (replaces Sys.Date() stamps).
.hs_utc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

# ---- write_outputs() ---------------------------------------------------------

#' Write a haploStack result to disk
#'
#' Writes the tables, figures (and, for a catalogue, the container) of an
#' \code{\link[=hs_result]{hs_result}} into \code{dir}, with the same file names
#' and contents the entry point wrote into \code{output_dir} up to 0.8.8. So a
#' result kept in memory (\code{output_dir = NULL}) can be written later, and
#' written again elsewhere, and gives the same files. One method per result
#' class; they land with each class (\code{hs_catalogue},
#' \code{hs_block_analysis}, \code{hs_simulation}).
#'
#' @param x An \code{hs_result}.
#' @param dir Output directory.
#' @param ... Method arguments: \code{what} (which kinds to write),
#'   \code{overwrite} (default \code{FALSE} when called directly: existing
#'   files are an error that names them; the entry points pass \code{TRUE}, as
#'   0.8.8 did), \code{create_dir} (default \code{TRUE}).
#' @return Invisibly, a manifest \code{data.frame(path, kind, object, bytes)}
#'   with one row per file written (\code{bytes} is \code{NA} for a figure that
#'   failed to render; the failure is also a note).
#' @seealso \code{\link{hs_result}}
#' @export
write_outputs <- function(x, dir, ...) UseMethod("write_outputs")

#' @export
write_outputs.default <- function(x, dir, ...) {
  .hs_abort("write_outputs: no method for an object of class ",
            paste(class(x), collapse = "/"),
            "; expected a haploStack result (hs_catalogue, hs_block_analysis, ",
            "hs_simulation).", call. = FALSE)
}
