# =============================================================================
# format_presets.R — the file-format presets and their resolver
#
# GENO_FORMAT_PRESETS is the registry of parameter sets (orientation, allele
# separator, map column order, FAM columns, encoding) that drive both the
# readers (read_genotype() and its wrappers) and the writers
# (geno_meta_to_matrix(), write_geno_meta()); resolve_geno_format() turns a
# preset name or a custom list into one validated spec. Split out of
# format_helper.R in 0.9.0 (D2-9); the history of the preset model is in
# docs/CHANGELOG.md (2026-06-12 onwards).
# =============================================================================

#' Genotype file-format presets (parameter-set registry)
#'
#' A named list of \emph{parameter sets}, one per supported genotype file
#' layout. Each parameter set is a named list describing how a delimited
#' genotype/map file pair is laid out on disk, and is consumed by both the
#' reader (\code{\link{read_genotype}}) and the writers
#' (\code{\link{geno_meta_to_matrix}}, \code{\link{write_geno_meta}}) so that a
#' format is described in exactly one place.
#'
#' @section Parameter-set fields:
#' \describe{
#'   \item{\code{geno_orient}}{\code{"sample"} (rows = samples, the fixed
#'     internal orientation) or \code{"marker"} (rows = markers; transposed on
#'     read and on write).}
#'   \item{\code{geno_header}}{Logical; whether the genotype file has a header
#'     row of column identifiers.}
#'   \item{\code{fam_cols}}{Integer count of leading FAM/pedigree columns
#'     embedded in the genotype file (PLINK-style). \code{0} means none; the
#'     first column is then the row identifier. When \code{> 0}, columns
#'     \code{1:fam_cols} are the FAM block and genotype values start at column
#'     \code{fam_cols + 1}.}
#'   \item{\code{map_order}}{Character length-3 permutation of
#'     \code{c("Chrom", "Marker", "Pos")} giving the \emph{file} column order of
#'     the map.}
#'   \item{\code{split_char}}{Allele delimiter within a genotype cell:
#'     \code{"/"} (e.g. \code{"1/4"}) or \code{""} (one character per allele,
#'     e.g. \code{"AT"}).}
#'   \item{\code{encoding}}{Allele vocabulary of the cells:
#'     \code{"nucleotide"} (A/T/C/G, or digits used as allele labels) or
#'     \code{"selectiontools"} (SelectionTools' codes 1/2/3/4 — the read is
#'     then behind \code{\link{use_st_formats}} and the codes are translated to
#'     nucleotides; see \code{\link{read_genotype}}).}
#'   \item{\code{source_format}}{Label stored in
#'     \code{metadata$source_format}.}
#' }
#'
#' @section Built-in presets:
#' \describe{
#'   \item{\code{default}}{Sample-as-row table with a header; first column =
#'     sample IDs, remaining column names = marker IDs. Map columns ordered
#'     \code{Chrom, Marker, Pos}.}
#'   \item{\code{PL}}{PLINK-style: sample-as-row, header, the first six columns
#'     are FAM fields (FID, IID, PAT, MAT, SEX, PHENO) and genotype values begin
#'     at column seven. Map columns ordered \code{Chrom, Marker, Pos}.}
#'   \item{\code{ST}}{SelectionTools: marker-as-row, header; first column =
#'     marker IDs, remaining column names = sample IDs; slash-separated
#'     SelectionTools codes (\code{1/4}), translated to nucleotides on the read
#'     behind \code{\link{use_st_formats}}. Map columns ordered \code{Marker,
#'     Chrom, Pos}.}
#'   \item{\code{GS}}{genomicSimulation: marker-as-row, header; two-character
#'     nucleotide cells. Map columns ordered \code{Marker, Chrom, Pos}.}
#' }
#'
#' Users may define their own presets by adding entries to a copy of this list
#' (or by passing a named list directly to the \code{preset} argument); see
#' \code{\link{resolve_geno_format}}.
#'
#' @export
GENO_FORMAT_PRESETS <- list(
  default = list(
    geno_orient   = "sample",
    geno_header   = TRUE,
    fam_cols      = 0L,
    map_order     = c("Chrom", "Marker", "Pos"),
    split_char    = "",
    encoding      = "nucleotide",
    source_format = "user_supplied"
  ),
  PL = list(
    geno_orient   = "sample",
    geno_header   = TRUE,
    fam_cols      = 6L,
    map_order     = c("Chrom", "Marker", "Pos"),
    split_char    = "",
    encoding      = "nucleotide",
    source_format = "plink"
  ),
  ST = list(
    geno_orient   = "marker",
    geno_header   = TRUE,
    fam_cols      = 0L,
    map_order     = c("Marker", "Chrom", "Pos"),
    split_char    = "/",
    encoding      = "selectiontools",
    source_format = "selectiontools"
  ),
  GS = list(
    geno_orient   = "marker",
    geno_header   = TRUE,
    fam_cols      = 0L,
    map_order     = c("Marker", "Chrom", "Pos"),
    split_char    = "",
    encoding      = "nucleotide",
    source_format = "genomic_simulation"
  )
)


# Canonical parameter-set field names (used to screen custom presets).
.GENO_FORMAT_FIELDS <- c("geno_orient", "geno_header", "fam_cols",
                         "map_order", "split_char", "encoding",
                         "source_format")


# Internal: normalise a preset name (case-insensitive) and resolve common
# aliases to a registry key.
.normalise_preset_key <- function(x) {
  switch(tolower(x),
    "default"            = "default",
    "pl"                 = ,
    "plink"              = ,
    "plink_ped"          = "PL",
    "st"                 = ,
    "selectiontools"     = "ST",
    "gs"                 = ,
    "genomicsimulation"  = ,
    "genomic_simulation" = "GS",
    toupper(x))
}


# Internal: validate a fully resolved format parameter set.
.validate_geno_format <- function(spec) {
  if (!is.list(spec))
    .hs_abort("resolve_geno_format: parameter set must be a named list.")
  if (!identical(spec$geno_orient, "sample") &&
      !identical(spec$geno_orient, "marker"))
    .hs_abort("resolve_geno_format: 'geno_orient' must be \"sample\" or \"marker\".")
  if (!is.logical(spec$geno_header) || length(spec$geno_header) != 1L)
    .hs_abort("resolve_geno_format: 'geno_header' must be a single logical.")
  if (!is.numeric(spec$fam_cols) || length(spec$fam_cols) != 1L ||
      spec$fam_cols < 0L)
    .hs_abort("resolve_geno_format: 'fam_cols' must be a single non-negative integer.")
  if (length(spec$map_order) != 3L ||
      !setequal(spec$map_order, c("Chrom", "Marker", "Pos")))
    .hs_abort("resolve_geno_format: 'map_order' must be a permutation of ",
         "c(\"Chrom\", \"Marker\", \"Pos\").")
  if (!spec$encoding %in% c("selectiontools", "nucleotide"))
    .hs_abort("resolve_geno_format: 'encoding' must be \"selectiontools\" or ",
         "\"nucleotide\" (got '", spec$encoding, "').")
  invisible(spec)
}


#' Resolve a genotype file-format parameter set
#'
#' Looks up a preset by name in a registry (a \emph{dictionary resolver}),
#' merges it over the \code{default} parameter set, and then applies any
#' explicit per-field overrides. A custom parameter set can be supplied
#' directly as a named list in \code{preset}, so callers are never restricted to
#' the built-in presets.
#'
#' @param preset  One of: \code{NULL} (use \code{default}); a preset name
#'   (\code{"PL"}, \code{"ST"}, \code{"GS"}, \code{"default"}, or an alias such
#'   as \code{"plink_ped"}/\code{"selectiontools"}); or a named list giving a
#'   custom parameter set (merged over \code{default}).
#' @param presets Registry of presets to resolve against. Default
#'   \code{\link{GENO_FORMAT_PRESETS}}; supply an extended copy to register
#'   bespoke formats.
#' @param split_char,encoding,geno_orient,geno_header,fam_cols,map_order,source_format
#'   Optional per-field overrides. Any non-\code{NULL} value replaces the
#'   corresponding field after the preset is applied.
#'
#' @return A validated named list with the fields documented in
#'   \code{\link{GENO_FORMAT_PRESETS}}.
#'
#' @seealso \code{\link{GENO_FORMAT_PRESETS}}, \code{\link{read_genotype}},
#'   \code{\link{write_geno_meta}}
#'
#' @examples
#' resolve_geno_format("ST")
#' # Custom parameter set: SelectionTools layout but nucleotide encoding
#' resolve_geno_format(list(geno_orient = "marker", split_char = "/",
#'                          encoding = "nucleotide"))
#'
#' @export
resolve_geno_format <- function(preset        = NULL,
                                presets       = GENO_FORMAT_PRESETS,
                                split_char    = NULL,
                                encoding      = NULL,
                                geno_orient   = NULL,
                                geno_header   = NULL,
                                fam_cols      = NULL,
                                map_order     = NULL,
                                source_format = NULL) {

  spec <- presets[["default"]]
  if (is.null(spec))
    .hs_abort("resolve_geno_format: 'presets' must contain a 'default' entry.")

  if (!is.null(preset)) {
    if (is.list(preset)) {
      p <- preset
    } else if (is.character(preset) && length(preset) == 1L) {
      key <- .normalise_preset_key(preset)
      p   <- presets[[key]]
      if (is.null(p))
        .hs_abort("resolve_geno_format: unknown preset '", preset, "'. Known ",
             "presets: ", paste(setdiff(names(presets), "default"),
                                 collapse = ", "),
             ". Pass a named list to define a custom parameter set.")
    } else {
      .hs_abort("resolve_geno_format: 'preset' must be NULL, a preset name, or a ",
           "named list parameter set.")
    }
    unknown <- setdiff(names(p), .GENO_FORMAT_FIELDS)
    if (length(unknown) > 0L)
      .hs_warn("resolve_geno_format: ignoring unknown parameter-set field(s): ",
              paste(unknown, collapse = ", "), ".")
    spec <- modifyList(spec, p[intersect(names(p), .GENO_FORMAT_FIELDS)])
  }

  # Explicit overrides (any non-NULL value wins over the preset/default)
  if (!is.null(split_char))    spec$split_char    <- split_char
  if (!is.null(encoding))      spec$encoding      <- encoding
  if (!is.null(geno_orient))   spec$geno_orient   <- geno_orient
  if (!is.null(geno_header))   spec$geno_header   <- as.logical(geno_header)
  if (!is.null(fam_cols))      spec$fam_cols      <- as.integer(fam_cols)
  if (!is.null(map_order))     spec$map_order     <- map_order
  if (!is.null(source_format)) spec$source_format <- source_format

  spec$fam_cols <- as.integer(spec$fam_cols)
  .validate_geno_format(spec)
  spec
}


# ---------------------------------------------------------------------------
# Missing-allele normalisation (MAN, 0.6.0)
# ---------------------------------------------------------------------------
