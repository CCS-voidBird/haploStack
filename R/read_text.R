# =============================================================================
# read_text.R — the one text-table reader (thread #032)
#
# .read_table_chr() reads every haploStack text input as character columns —
# an identifier like "007" keeps its spelling, an allele column never turns
# logical — and refuses a ragged file with the offending line named; blank
# lines are skipped. .read_map_file() reads a map through it and reorders the
# columns to a preset's map_order. Every reader descends from these two.
# Split out of format_helper.R in 0.9.0 (D2-9).
# =============================================================================

# .read_table_chr(file, header = FALSE, sep = "auto", strict = FALSE)
#
# Read a whitespace/delimited text table with EVERY column as character, so
# that R's type inference can never reinterpret a cell's contents.
#
# This exists because of one specific, silent failure. utils::read.table() runs
# type.convert() on each column, and type.convert() maps a column whose every
# value is "T" (or "F") to LOGICAL. In a nucleotide PED that is an ordinary
# allele column: a T/A or T/G marker that happens to be monomorphic-T across
# the panel. The column comes back TRUE, `mode(x) <- "character"` then yields
# "TRUE", and .reshape_wide_alleles() refuses it —
#
#   .reshape_wide_alleles: cells contain multi-character values.
#   Wide format requires exactly one character per cell.
#
# — so the whole file fails to read. It is not rare: on the GRDC barley AGG1
# panel 511 of 19,712 allele columns are affected (and the count depends on
# which lines are in the panel, so a subset can make it larger).
#
# The same trap reaches two more readers, for the same reason and with worse
# consequences, so they share this primitive:
#   * read_effect_file()   — an all-"T" Allele column becomes TRUE, and an
#     allele spelled "TRUE" matches nothing, so every effect silently
#     contributes 0. No error is raised anywhere.
#   * read_dart_genotype() — a logical column breaks the metadata/sample column
#     heuristic (as.numeric("TRUE") is NA), shifting where the sample block is
#     judged to start.
#
# The rule this encodes: a genotype file is read as TEXT. Every column arrives
# as character and the caller coerces, explicitly, the columns it knows to be
# numeric (positions, effects, dosages). Nothing about a genotype's meaning is
# left to a heuristic that cannot tell an allele from a boolean.
#
# data.table::fread is used rather than read.table(colClasses = "character")
# because it is the package's existing reader everywhere else, it does not
# type-convert when colClasses is given, and it is far faster on a wide PED
# (a 43 MB / 19,718-column file reads in seconds rather than minutes).
#
# RAGGED INPUT — read this before changing `fill`.
#
# `fread(fill = FALSE)` does NOT reliably error on a ragged file. Given a PED
# whose first line is longer than the rest, it settles on the narrower width and
# **silently drops the wider row**, returning a table with one fewer sample and
# no warning. That is strictly worse than the defect this function exists to
# fix, so the fread call always uses `fill = TRUE` (short rows padded with NA,
# no row ever discarded) and raggedness is detected here instead:
#
#   * every data line in the file must appear as a row (`n_lines` guard), and
#   * with `strict = TRUE` the table must be rectangular — no padded cell in the
#     last column. A PED/MAP line is fixed-width by construction, so a short one
#     is a truncated file, not missing data, and must fail loudly rather than
#     become a row of missing calls.
#
# `strict = FALSE` keeps the lenient read.table(fill = TRUE) behaviour the
# effect-file and DArT readers always had.
#
# BLANK LINES (0.7.9). A blank line — empty, or nothing but spaces and tabs —
# is not a record, so the guard counts it out, and fread must not return it
# either. With fread's default blank.lines.skip = FALSE a blank line between
# records, or a second trailing one, came back as a row, the counts disagreed,
# and a well-formed file was refused as "ragged" (INDEX discrepancy 33).
# blank.lines.skip = TRUE handles an EMPTY line. It does not handle a line of
# spaces or tabs: with sep = "\t" a "\t\t" line is two empty fields to fread,
# so it comes back as a row. A file holding such a line is therefore parsed
# from its lines with the blank ones taken out (fread(text = )); every other
# file is parsed from disk exactly as before.
#
# Deliberately NOT blank: a line of separators such as ",,". It holds
# non-whitespace characters, and up to 0.7.8 it read as a row of missing
# values; it still does. (Dropping all-empty rows instead would refuse such a
# file as ragged, or shift the rows a caller indexes — the first 0.7.9 build did
# the former.)
#
# @keywords internal
.read_table_chr <- function(file, header = FALSE, sep = "auto",
                            strict = FALSE) {
  file <- as.character(file)
  if (!file.exists(file))
    .hs_abort(".read_table_chr: file not found: ", file)
  if (file.size(file) == 0L)
    .hs_abort(".read_table_chr: file is empty: ", file)

  # Guard input: the records are counted from the file rather than trusted from
  # the parser, because a parser that drops rows is exactly what we are
  # defending against. One pass — a 43 MB PED is expensive to re-read.
  lines      <- readLines(file, warn = FALSE)
  blank      <- !nzchar(trimws(lines))
  if (all(blank))
    .hs_abort(".read_table_chr: no data read from ", file)
  data_lines <- which(!blank)                          # file line numbers
  if (isTRUE(header)) data_lines <- data_lines[-1L]
  expected   <- length(data_lines)
  # A line of spaces or tabs is parsed from the lines, without it (see BLANK
  # LINES above); an empty line fread skips itself.
  text <- if (any(blank & nzchar(lines))) lines[!blank] else NULL
  rm(lines)

  # fill = TRUE always: see the note above. Dropping a row silently is the one
  # failure mode this function must never have.
  rd <- function(...) data.table::fread(..., header = header, sep = sep,
                                        fill = TRUE, blank.lines.skip = TRUE,
                                        colClasses = "character",
                                        data.table = FALSE, check.names = FALSE,
                                        showProgress = FALSE)
  out <- if (is.null(text)) rd(file) else rd(text = text)
  rm(text)
  if (nrow(out) == 0L || ncol(out) == 0L)
    .hs_abort(".read_table_chr: no data read from ", file)

  # fread returns "" for an empty field; every caller treats an absent value as
  # NA, and no allele vocabulary contains the empty string.
  for (j in seq_len(ncol(out))) {
    v <- out[[j]]
    if (is.character(v)) out[[j]] <- ifelse(nzchar(v), v, NA_character_)
  }

  # Guard: no line may go missing.
  if (nrow(out) != expected)
    .hs_abort(".read_table_chr: read ", nrow(out), " row(s) from ", file,
         " but the file holds ", expected, " data line(s). The file is ",
         "ragged \u2014 lines do not all carry the same number of fields.")

  if (isTRUE(strict) && ncol(out) > 0L) {
    short <- which(is.na(out[[ncol(out)]]))
    if (length(short) > 0L)
      .hs_abort(".read_table_chr: ", length(short), " line(s) in ", file,
           " are shorter than the widest line (first: line ",
           data_lines[short[1L]],
           "). A fixed-width genotype file must be rectangular; padding a ",
           "truncated line into missing calls would hide the damage.")
  }
  out
}


#' Read a map file and reorder its columns to (Chrom, Marker, Pos) (internal)
#'
#' Reads a delimited map file and returns a three-column data.frame in the
#' canonical \code{Chrom, Marker, Pos} order. Column \emph{roles} are taken
#' positionally from \code{map_order} (the file column order). Contents are not
#' modified — position-unit conversion is left to \code{\link{build_geno_meta}}.
#'
#' @param map_file  Path to the map file (assumed to have a header).
#' @param map_order Length-3 permutation of \code{c("Chrom","Marker","Pos")}
#'   giving the role of each of the first three file columns.
#' @return data.frame with columns \code{Chrom}, \code{Marker}, \code{Pos}.
#' @keywords internal
#' @noRd
.read_map_file <- function(map_file, map_order) {
  map_file <- as.character(map_file)
  if (!file.exists(map_file))
    .hs_abort(".read_map_file: map file not found: ", map_file)
  pos_of <- function(role) which(map_order == role)[1L]
  # Chrom and Marker are read as text ("007" stays "007"); Pos keeps fread's type.
  id_cols <- c(pos_of("Chrom"), pos_of("Marker"))
  raw_map <- data.table::fread(map_file, data.table = FALSE,
                               colClasses = list(character = id_cols[!is.na(id_cols)]))
  if (ncol(raw_map) < 3L)
    .hs_abort(".read_map_file: map file must have at least three columns (roles: ",
         paste(map_order, collapse = ", "), ").")
  data.frame(
    Chrom  = as.character(raw_map[[pos_of("Chrom")]]),
    Marker = as.character(raw_map[[pos_of("Marker")]]),
    Pos    = raw_map[[pos_of("Pos")]],
    stringsAsFactors = FALSE
  )
}
