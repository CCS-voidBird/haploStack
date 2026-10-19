# =============================================================================
# phenotype.R
#
# Sample-level phenotype I/O for the GenoMeta container.
#
# A phenotype table is a data.frame with a `Sample` column plus one or more
# numeric trait columns. It is stored as the optional `phenotype` element of a
# GenoMeta (see docs/geno_meta_structure.md) and is keyed on `Sample` the way
# `marker_effects` is keyed on `Marker`.
#
#   read_phenotype()     — read a delimited .pheno/CSV/TSV into a phenotype table
#   attach_phenotype()   — filter a phenotype to shared samples + store on a meta
#   get_phenotype()      — accessor for the stored phenotype (optional subset)
# =============================================================================


#' Read a delimited phenotype file into a phenotype table
#'
#' Reads a \code{.pheno}/CSV/TSV table of sample-level trait values and returns a
#' tidy data.frame with a \code{Sample} column followed by numeric trait columns.
#' Leading pedigree/identifier columns commonly shipped alongside phenotypes
#' (\code{FID}) are dropped; the remaining non-sample columns are coerced to
#' numeric.
#'
#' @param file       Path to the phenotype file.
#' @param sample_col Column to use as the sample identifier. \code{NULL}
#'   (default) auto-detects: the first of \code{Sample}, \code{IID}, \code{ID}
#'   present, else the first column. May be a column name or a 1-based index.
#' @param drop_cols  Character vector of additional id columns to drop before
#'   treating the rest as traits. Default \code{c("FID")}.
#' @param traits     Optional character vector of trait columns to keep (in this
#'   order). \code{NULL} (default) keeps every non-id column.
#' @param sep        Field separator passed to \code{data.table::fread}.
#'   \code{NULL} (default) lets \code{fread} auto-detect.
#'
#' @return A data.frame with column \code{Sample} (character) followed by numeric
#'   trait columns — the shape of the GenoMeta \code{phenotype} element.
#'
#' @seealso \code{\link{attach_phenotype}}, \code{\link{get_phenotype}}
#'
#' @examples
#' \dontrun{
#' ph <- read_phenotype("traits.pheno")
#' ph <- read_phenotype("traits.csv", traits = c("NFNB", "SFNB", "Scald"))
#' }
#'
#' @export
read_phenotype <- function(file,
                           sample_col = NULL,
                           drop_cols  = c("FID"),
                           traits     = NULL,
                           sep        = NULL) {
  file <- as.character(file)
  if (!file.exists(file))
    stop("read_phenotype: file not found: ", file)

  raw <- data.table::fread(file, header = TRUE, data.table = FALSE,
                           sep = if (is.null(sep)) "auto" else sep)
  if (ncol(raw) < 2L)
    stop("read_phenotype: file must have an identifier column plus at least ",
         "one trait column.")

  # Resolve the sample-id column.
  if (is.null(sample_col)) {
    sample_col <- intersect(c("Sample", "IID", "ID"), names(raw))[1L]
    if (is.na(sample_col)) sample_col <- names(raw)[1L]
  } else if (is.numeric(sample_col)) {
    sample_col <- names(raw)[as.integer(sample_col)]
  }
  if (!sample_col %in% names(raw))
    stop("read_phenotype: sample column '", sample_col, "' not found in file.")

  sample_ids <- as.character(raw[[sample_col]])

  # Trait columns = everything except the sample column and any drop_cols.
  trait_cols <- setdiff(names(raw), c(sample_col, drop_cols))
  if (!is.null(traits)) {
    missing_t <- setdiff(traits, trait_cols)
    if (length(missing_t) > 0L)
      stop("read_phenotype: requested trait(s) not found: ",
           paste(missing_t, collapse = ", "), ".")
    trait_cols <- traits
  }
  if (length(trait_cols) == 0L)
    stop("read_phenotype: no trait columns remain after dropping id columns.")

  ph <- data.frame(Sample = sample_ids, stringsAsFactors = FALSE)
  for (tc in trait_cols)
    ph[[tc]] <- suppressWarnings(as.numeric(raw[[tc]]))
  ph
}


#' Attach a phenotype table to a GenoMeta, filtered to shared samples
#'
#' Stores a phenotype table on a GenoMeta's \code{phenotype} element after
#' restricting it to the samples shared with the genotype data
#' (\code{intersect(phenotype$Sample, sample_info$Sample)}). Genotype samples
#' lacking a phenotype row are retained with \code{NA} trait values; phenotype
#' rows for unknown samples are dropped with a warning. The result is rebuilt
#' through \code{\link{new_geno_meta}} so the whole object is revalidated.
#'
#' @param geno_meta A valid \code{GenoMeta} list.
#' @param phenotype A phenotype data.frame (\code{Sample} + numeric trait
#'   columns), e.g. from \code{\link{read_phenotype}}, or a path to a file to be
#'   read with \code{\link{read_phenotype}}.
#' @param traits    Optional character vector of trait columns to keep.
#'   \code{NULL} (default) keeps all.
#'
#' @return The \code{GenoMeta} with its \code{phenotype} element populated.
#'
#' @seealso \code{\link{read_phenotype}}, \code{\link{get_phenotype}},
#'   \code{\link{plot_phenotype_distribution}}
#'
#' @examples
#' \dontrun{
#' meta <- attach_phenotype(meta, read_phenotype("traits.pheno"))
#' meta <- attach_phenotype(meta, "traits.pheno", traits = c("NFNB", "SFNB"))
#' }
#'
#' @export
attach_phenotype <- function(geno_meta, phenotype, traits = NULL) {
  if (!is_geno_meta(geno_meta))
    stop("attach_phenotype: 'geno_meta' is not a valid GenoMeta.")

  if (is.character(phenotype) && length(phenotype) == 1L)
    phenotype <- read_phenotype(phenotype, traits = traits)
  if (!is.data.frame(phenotype))
    stop("attach_phenotype: 'phenotype' must be a data.frame or a file path.")
  phenotype <- as.data.frame(phenotype, stringsAsFactors = FALSE)

  # Resolve / normalise the Sample column to be named "Sample".
  sample_col <- if ("Sample" %in% names(phenotype)) "Sample"
                else names(phenotype)[1L]
  if (!identical(sample_col, "Sample"))
    names(phenotype)[match(sample_col, names(phenotype))] <- "Sample"

  # Optional trait subset (when not already applied during read).
  if (!is.null(traits)) {
    missing_t <- setdiff(traits, setdiff(names(phenotype), "Sample"))
    if (length(missing_t) > 0L)
      stop("attach_phenotype: requested trait(s) not found: ",
           paste(missing_t, collapse = ", "), ".")
    phenotype <- phenotype[, c("Sample", traits), drop = FALSE]
  }

  geno_samples <- dimnames(geno_meta$geno_arrays)[[1L]]
  shared <- intersect(as.character(phenotype$Sample), geno_samples)
  if (length(shared) == 0L)
    stop("attach_phenotype: no samples shared between the phenotype and the ",
         "genotype data.")

  # .standardise_phenotype() (via new_geno_meta) performs the drop/NA-fill/reorder
  # and warns; rebuild so the entire object is revalidated.
  out <- new_geno_meta(
    geno_arrays           = geno_meta$geno_arrays,
    map_info              = geno_meta$map_info,
    block_info            = geno_meta$block_info,
    block_marker_map_info = geno_meta$block_marker_map_info,
    metadata              = geno_meta$metadata,
    sample_info           = geno_meta$sample_info,
    marker_effects        = geno_meta$marker_effects,
    phenotype             = phenotype
  )
  # Blocks are unchanged, so a stored block annotation stays valid (#033).
  out$annotation <- geno_meta$annotation
  out
}


#' Retrieve the phenotype table stored on a GenoMeta
#'
#' @param geno_meta A valid \code{GenoMeta} list.
#' @param traits    Optional character vector of trait columns to return (with
#'   \code{Sample}). \code{NULL} (default) returns all columns.
#'
#' @return A data.frame with \code{Sample} + trait columns.
#'
#' @seealso \code{\link{attach_phenotype}}
#'
#' @export
get_phenotype <- function(geno_meta, traits = NULL) {
  if (!is_geno_meta(geno_meta))
    stop("get_phenotype: 'geno_meta' is not a valid GenoMeta.")
  ph <- geno_meta$phenotype
  if (is.null(ph))
    stop("get_phenotype: no phenotype attached. Use attach_phenotype() first.")
  if (!is.null(traits)) {
    missing_t <- setdiff(traits, setdiff(names(ph), "Sample"))
    if (length(missing_t) > 0L)
      stop("get_phenotype: requested trait(s) not found: ",
           paste(missing_t, collapse = ", "), ".")
    ph <- ph[, c("Sample", traits), drop = FALSE]
  }
  ph
}
