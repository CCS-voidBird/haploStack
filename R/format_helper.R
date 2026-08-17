# =============================================================================
# format_helper.R
#
# GenoMeta — unified internal metadata container for genotype data.
#
# A GenoMeta object is a plain named list (no formal class) that all format
# readers return and all transform functions validate on entry.
#
#   new_geno_meta()          — constructor + dim validator
#   is_geno_meta()           — lightweight dim-alignment check
#   GENO_FORMAT_PRESETS      — preset registry (parameter sets per file format)
#   resolve_geno_format()    — preset/parameter-set resolver (read + write)
#   .read_genotype_files()   — internal: read geno/map/fam files unmodified
#   .read_map_file()         — internal: read + reorder a map file
#   read_genotype()          — reader: delimited geno table via preset
#   read_plink_ped()         — reader: Plink 1.9 PED + MAP
#   read_st_genotype()       — reader: SelectionTools CSV/tab
#   read_gs_genotype()       — reader: genomicSimulation CSV
#   build_geno_meta()        — assemble GenoMeta from fixed-shape inputs
#   .read_block_file()       — internal: load block_info + block_marker_map
#   filter_geno_meta()       — intersect-filter markers and/or samples
#   geno_meta_to_matrix()    — export 3D → 2D for a given output format
#   write_geno_meta()        — write converted files to disk
#
# Q5 (2026-06-08): Mbp detection warning delegated to .process_map().
# Q8 (2026-06-08): Plain named list; strict validation only in constructors
#                  and transform functions.
# 2026-06-12: Format handling unified on a preset/parameter-set model. A single
#             GENO_FORMAT_PRESETS registry + resolve_geno_format() drives both
#             reading and writing. build_geno_meta() now takes FIXED-orientation
#             inputs (geno [sample × marker], map [Chrom, Marker, Pos], fam
#             [Family, Sample, Parent1, Parent2, Sex]) and only the changeable
#             format parameters split_char / pos_unit / encoding. The encoding
#             vocabulary is "selectiontools" / "nucleotide", matching
#             detect_genotype_encoding() and the catalogue pipeline.
# 0.5.1     — the former "number" / "char" tokens are REMOVED, not just
#             superseded. They had lingered as accepted aliases here and as the
#             live vocabulary of translate_alleles()/.resolve_preset(), so one
#             concept had two names that met in process_genotype() and did not
#             match. See .encode_arg() in format_detect.R.
# =============================================================================


# ---------------------------------------------------------------------------
# Format presets (parameter sets) and resolver
# ---------------------------------------------------------------------------

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
#'   \item{\code{encoding}}{Allele encoding: \code{"selectiontools"} (numeric
#'     codes 1/2/3/4) or \code{"nucleotide"} (A/T/C/G).}
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
#'     marker IDs, remaining column names = sample IDs; slash-separated numeric
#'     cells. Map columns ordered \code{Marker, Chrom, Pos}.}
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
    stop("resolve_geno_format: parameter set must be a named list.")
  if (!identical(spec$geno_orient, "sample") &&
      !identical(spec$geno_orient, "marker"))
    stop("resolve_geno_format: 'geno_orient' must be \"sample\" or \"marker\".")
  if (!is.logical(spec$geno_header) || length(spec$geno_header) != 1L)
    stop("resolve_geno_format: 'geno_header' must be a single logical.")
  if (!is.numeric(spec$fam_cols) || length(spec$fam_cols) != 1L ||
      spec$fam_cols < 0L)
    stop("resolve_geno_format: 'fam_cols' must be a single non-negative integer.")
  if (length(spec$map_order) != 3L ||
      !setequal(spec$map_order, c("Chrom", "Marker", "Pos")))
    stop("resolve_geno_format: 'map_order' must be a permutation of ",
         "c(\"Chrom\", \"Marker\", \"Pos\").")
  if (!spec$encoding %in% c("selectiontools", "nucleotide"))
    stop("resolve_geno_format: 'encoding' must be \"selectiontools\" or ",
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
    stop("resolve_geno_format: 'presets' must contain a 'default' entry.")

  if (!is.null(preset)) {
    if (is.list(preset)) {
      p <- preset
    } else if (is.character(preset) && length(preset) == 1L) {
      key <- .normalise_preset_key(preset)
      p   <- presets[[key]]
      if (is.null(p))
        stop("resolve_geno_format: unknown preset '", preset, "'. Known ",
             "presets: ", paste(setdiff(names(presets), "default"),
                                 collapse = ", "),
             ". Pass a named list to define a custom parameter set.")
    } else {
      stop("resolve_geno_format: 'preset' must be NULL, a preset name, or a ",
           "named list parameter set.")
    }
    unknown <- setdiff(names(p), .GENO_FORMAT_FIELDS)
    if (length(unknown) > 0L)
      warning("resolve_geno_format: ignoring unknown parameter-set field(s): ",
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
# GenoMeta structure
# ---------------------------------------------------------------------------

#' Create and validate a GenoMeta internal metadata object
#'
#' A \strong{GenoMeta} is a plain named list that consolidates all genotype
#' and map data into a single structure consumed by the haploStack pipeline.
#' No formal class system is used — any field can be read or modified directly.
#' Strict dimension validation is applied only here and in transform functions.
#'
#' @section Structure:
#' \describe{
#'   \item{\code{geno_arrays}}{3-D integer/character array
#'     \code{[n_samples × n_markers × ploidy]}. Values are encoded per
#'     \code{metadata$char_dict}.}
#'   \item{\code{map_info}}{data.frame \code{(Chrom, Marker, Pos)}.
#'     \code{nrow(map_info) == dim(geno_arrays)[2]}.}
#'   \item{\code{block_info}}{data.frame \code{(Chr, Block, Start, End)} or
#'     \code{NULL}. Block rows are preserved even when no markers remain.}
#'   \item{\code{block_marker_map_info}}{data.frame \code{(Chr, Marker, Block)}
#'     or \code{NULL}. \code{Marker} column must be a subset of
#'     \code{map_info$Marker}.}
#'   \item{\code{sample_info}}{Required FAM-structured data.frame with positional
#'     columns \code{(Family, Sample, Parent1, Parent2, Sex)}. Column \emph{index}
#'     is the contract — an input fam is read by position and stored under these
#'     canonical names (so \code{$Sample} is stable). \code{Sample} (position 2)
#'     must match \code{dimnames(geno_arrays)[[1]]}. When \code{NULL}, a dummy
#'     table is built from the genotype sample names (other columns \code{"0"}).
#'     See \code{docs/geno_meta_structure.md}.}
#'   \item{\code{marker_effects}}{Optional allele effects, or \code{NULL}.
#'     Stored canonically as a \strong{named list of per-trait data.frames},
#'     each \code{(Marker, Allele, Effect)} with the list name giving the trait
#'     — the same shape as the \code{effect_list} consumed by
#'     \code{\link{run_haplotype_catalogue}} and produced by
#'     \code{\link{fit_genomic_prediction}}. A bare data.frame (optionally with a
#'     \code{Trait} column to split on) is accepted and normalised to this list.
#'     Used as the default effect source by
#'     \code{\link{compute_profile_localgebv}} and
#'     \code{\link{run_haplotype_catalogue}}; an explicit effect argument always
#'     overrides it. Never required for a valid GenoMeta.}
#'   \item{\code{phenotype}}{Optional data.frame of sample-level trait values, or
#'     \code{NULL}. First column \code{Sample} (matching
#'     \code{dimnames(geno_arrays)[[1]]}); remaining columns are numeric trait
#'     measurements. Keyed on \code{Sample} the way \code{marker_effects} is keyed
#'     on \code{Marker}: standardised on construction, pruned to the surviving
#'     samples by \code{\link{filter_geno_meta}}. Attach via
#'     \code{\link{attach_phenotype}} / \code{\link{read_phenotype}}.}
#'   \item{\code{metadata}}{Named list:
#'     \code{ploidy} (int), \code{char_dict} (named char vector),
#'     \code{n_samples} (int), \code{n_markers} (int),
#'     \code{source_format} (character), and \code{encoding} (auto-detected
#'     when absent).}
#' }
#'
#' @param geno_arrays         3-D array \code{[n_samples × n_markers × ploidy]}.
#' @param map_info            data.frame with \code{Chr}, \code{Marker},
#'   \code{Pos}.
#' @param block_info          data.frame with \code{Chr}, \code{Block},
#'   \code{Start}, \code{End}, or \code{NULL}.
#' @param block_marker_map_info data.frame with \code{Chr}, \code{Marker},
#'   \code{Block}, or \code{NULL}.
#' @param metadata            Named list with at minimum \code{ploidy},
#'   \code{char_dict}, \code{n_samples}, \code{n_markers},
#'   \code{source_format}.
#' @param sample_info         FAM-structured data.frame read by position
#'   (\code{Family}, \code{Sample}, \code{Parent1}, \code{Parent2}, \code{Sex}),
#'   or \code{NULL} (default) to auto-build a dummy table from the genotype
#'   sample names. A one-column input is treated as the \code{Sample} list.
#' @param marker_effects      Optional allele effects, or \code{NULL} (default).
#'   Either a named list of per-trait data.frames (\code{Marker}, \code{Allele},
#'   \code{Effect}) or a single data.frame (optionally with a \code{Trait}
#'   column); normalised on construction to the named-list form.
#' @param phenotype          Optional data.frame of sample-level trait values
#'   (\code{Sample} + numeric trait columns), or \code{NULL} (default). Aligned to
#'   the genotype sample set on construction (see \code{\link{attach_phenotype}}).
#'
#' @return A plain named list of class \code{"list"}. Access fields directly
#'   with \code{$}.
#'
#' @export
new_geno_meta <- function(geno_arrays,
                           map_info,
                           block_info              = NULL,
                           block_marker_map_info   = NULL,
                           metadata,
                           sample_info             = NULL,
                           marker_effects          = NULL,
                           phenotype               = NULL) {

  # ---- Validate required fields ----------------------------------------
  if (!is.array(geno_arrays) || length(dim(geno_arrays)) != 3L)
    stop("new_geno_meta: 'geno_arrays' must be a 3-D array.")
  if (!is.data.frame(map_info))
    stop("new_geno_meta: 'map_info' must be a data.frame.")
  required_map <- c("Chrom", "Marker", "Pos")
  missing_map  <- setdiff(required_map, names(map_info))
  if (length(missing_map) > 0L)
    stop("new_geno_meta: 'map_info' missing columns: ",
         paste(missing_map, collapse = ", "))
  if (dim(geno_arrays)[2L] != nrow(map_info))
    stop("new_geno_meta: dim(geno_arrays)[2] (", dim(geno_arrays)[2L],
         ") must equal nrow(map_info) (", nrow(map_info), ").")

  if (!is.null(block_marker_map_info)) {
    if (!is.data.frame(block_marker_map_info))
      stop("new_geno_meta: 'block_marker_map_info' must be a data.frame or NULL.")

    # Filter: keep only markers present in map_info
    bad_markers <- setdiff(block_marker_map_info$Marker, map_info$Marker)
    if (length(bad_markers) > 0L) {
      warning("new_geno_meta: ", length(bad_markers),
              " marker(s) in block_marker_map_info absent from map_info — dropped: ",
              paste(head(bad_markers, 5L), collapse = ", "),
              if (length(bad_markers) > 5L) " ..." else "")
      block_marker_map_info <- block_marker_map_info[
        block_marker_map_info$Marker %in% map_info$Marker, , drop = FALSE
      ]
    }

    # Filter block_info to only blocks still represented after marker pruning
    if (!is.null(block_info)) {
      retained_blocks <- unique(block_marker_map_info$Block)
      block_info <- block_info[
        block_info$Block %in% retained_blocks, , drop = FALSE
      ]
    }

    # Coerce identifier columns now so rbind() below stays type-consistent
    block_marker_map_info$Chrom  <- as.character(block_marker_map_info$Chrom)
    block_marker_map_info$Marker <- as.character(block_marker_map_info$Marker)
    block_marker_map_info$Block  <- as.character(block_marker_map_info$Block)

    # Dummy blocks: every map_info marker absent from block_marker_map_info
    # becomes its own singleton block (one marker per dummy block).
    absent_markers <- setdiff(map_info$Marker, block_marker_map_info$Marker)
    if (length(absent_markers) > 0L) {
      message("new_geno_meta: ", length(absent_markers),
              " map marker(s) absent from block_marker_map_info — ",
              "created singleton dummy block(s).")

      map_idx     <- match(absent_markers, map_info$Marker)
      dummy_chrom <- as.character(map_info$Chrom[map_idx])
      dummy_pos   <- as.numeric(map_info$Pos[map_idx])
      dummy_block <- paste0("dummy_", absent_markers)

      block_marker_map_info <- rbind(
        block_marker_map_info[, c("Chrom", "Marker", "Block"), drop = FALSE],
        data.frame(Chrom  = dummy_chrom,
                   Marker = as.character(absent_markers),
                   Block  = dummy_block,
                   stringsAsFactors = FALSE)
      )

      # Update block_info with one row per dummy block (Start = End = Pos)
      if (!is.null(block_info)) {
        block_info$Chrom <- as.character(block_info$Chrom)
        block_info$Block <- as.character(block_info$Block)
        dummy_bi <- data.frame(Chrom = dummy_chrom,
                               Block = dummy_block,
                               Start = dummy_pos,
                               End   = dummy_pos,
                               stringsAsFactors = FALSE)
        # Fill any extra block_info columns (e.g. Num_SNPs) so rbind aligns
        extra_cols <- setdiff(names(block_info), names(dummy_bi))
        for (cc in extra_cols) dummy_bi[[cc]] <- NA
        dummy_bi <- dummy_bi[, names(block_info), drop = FALSE]
        block_info <- rbind(block_info, dummy_bi)
      }
    }

    # Reorder block_marker_map_info to follow map_info marker order
    block_marker_map_info <- block_marker_map_info[
      order(match(block_marker_map_info$Marker, map_info$Marker)), , drop = FALSE
    ]
    rownames(block_marker_map_info) <- NULL

    # Reorder block_info by genomic position (Chrom, then Start)
    if (!is.null(block_info)) {
      block_info <- block_info[
        order(block_info$Chrom, block_info$Start), , drop = FALSE
      ]
      rownames(block_info) <- NULL
    }
  }

  required_meta <- c("ploidy", "char_dict", "n_samples", "n_markers", "source_format")
  if (!is.list(metadata) || !all(required_meta %in% names(metadata)))
    stop("new_geno_meta: 'metadata' must be a named list with fields: ",
         paste(required_meta, collapse = ", "))

  # ---- sample_info: required FAM-structured field, aligned to genotypes ------
  geno_samples <- dimnames(geno_arrays)[[1L]]
  sample_info  <- .standardise_sample_info(sample_info, geno_samples,
                                           caller = "new_geno_meta")

  # ---- marker_effects: optional; normalise to the canonical named list of
  #      per-trait effect tables (Marker, Allele, Effect). A bare data.frame
  #      (optionally carrying a Trait column) is accepted for convenience. ------
  if (!is.null(marker_effects)) {
    marker_effects <- .as_effect_list(marker_effects)
  }

  # ---- phenotype: optional, align to the genotype sample set -----------------
  if (!is.null(phenotype)) {
    phenotype <- .standardise_phenotype(phenotype, geno_samples,
                                        caller = "new_geno_meta")
  }

  # ---- encoding: store in metadata when absent (best-effort detection) -------
  if (is.null(metadata$encoding)) {
    metadata$encoding <- tryCatch(
      detect_genotype_encoding(.meta_to_allele_matrix(list(geno_arrays = geno_arrays))),
      error = function(e) NULL
    )
  }

  # ---- Coerce Chrom / Marker / Block columns to character --------------------
  map_info$Chrom  <- as.character(map_info$Chrom)
  map_info$Marker <- as.character(map_info$Marker)

  if (!is.null(block_marker_map_info)) {
    block_marker_map_info$Chrom  <- as.character(block_marker_map_info$Chrom)
    block_marker_map_info$Marker <- as.character(block_marker_map_info$Marker)
    block_marker_map_info$Block  <- as.character(block_marker_map_info$Block)
  }

  if (!is.null(block_info)) {
    block_info$Chrom <- as.character(block_info$Chrom)
    block_info$Block <- as.character(block_info$Block)
  }

  # Return plain named list — no class attribute
  list(
    geno_arrays           = geno_arrays,
    map_info              = map_info,
    block_info            = block_info,
    block_marker_map_info = block_marker_map_info,
    sample_info           = sample_info,
    marker_effects        = marker_effects,
    phenotype             = phenotype,
    metadata              = metadata
  )
}


# Internal: normalise any marker-effects object to the canonical named list of
# per-trait effect tables. Each element is a data.frame (Marker, Allele, Effect)
# and the element name is the trait. Accepts:
#   * NULL                                    -> NULL
#   * a named list of data.frames             -> standardised element-wise
#   * a single data.frame with a Trait column -> split into one table per trait
#   * a single data.frame without Trait       -> wrapped under `default_name`
.as_effect_list <- function(x, default_name = "Trait") {
  if (is.null(x)) return(NULL)

  if (is.data.frame(x)) {
    if ("Trait" %in% names(x)) {
      keep <- setdiff(names(x), "Trait")
      sp   <- split(x[, keep, drop = FALSE], as.character(x$Trait))
      return(lapply(sp, .standardise_effect_df))
    }
    return(stats::setNames(list(.standardise_effect_df(x)), default_name))
  }

  if (is.list(x)) {
    if (length(x) == 0L) return(list())
    if (is.null(names(x)) || any(!nzchar(names(x))))
      stop("marker effects given as a list must be fully named ",
           "(one element name per trait).")
    return(lapply(x, .standardise_effect_df))
  }

  stop("marker effects must be a data.frame or a named list of data.frames.")
}

# Internal: standardise a single effect table to (Marker, Allele, Effect),
# reading the first three columns by position.
.standardise_effect_df <- function(d) {
  if (!is.data.frame(d) && !inherits(d, "data.table"))
    stop("each marker-effect table must be a data.frame.")
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  if (ncol(d) < 3L)
    stop("each marker-effect table needs >= 3 columns (Marker, Allele, Effect).")
  data.frame(
    Marker = as.character(d[[1L]]),
    Allele = as.character(d[[2L]]),
    Effect = as.numeric(d[[3L]]),
    stringsAsFactors = FALSE
  )
}

# Internal: collapse a marker-effects object (named list of per-trait tables or
# a bare data.frame) to a single effect table, or NULL when empty. Used where a
# single-trait table is required (e.g. founder-pool breeding values). `trait`
# selects a named element; otherwise the first table is used.
.coerce_single_effect_df <- function(x, trait = NULL) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    if (length(x) == 0L) return(NULL)
    if (!is.null(trait) && trait %in% names(x)) return(x[[trait]])
    return(x[[1L]])
  }
  NULL
}

# Internal: standardise a phenotype table to (Sample, <trait...>) aligned to the
# genotype sample set.
#
#   ph           data.frame with a Sample column (or a first column treated as
#                Sample) plus one or more numeric trait columns
#   geno_samples character vector dimnames(geno_arrays)[[1]]
#
# The Sample column is matched against geno_samples: rows for unknown samples are
# dropped with a warning, genotype samples lacking a phenotype row are added with
# NA trait values (left-join semantics), and the result is reordered to the
# genotype sample order. Trait columns are coerced to numeric.
.standardise_phenotype <- function(ph, geno_samples, caller = "new_geno_meta") {
  if (!is.data.frame(ph))
    stop(caller, ": 'phenotype' must be a data.frame or NULL.")
  ph <- as.data.frame(ph, stringsAsFactors = FALSE)
  if (ncol(ph) < 2L)
    stop(caller, ": 'phenotype' must have a Sample column plus at least one ",
         "trait column.")

  # Resolve the Sample column: prefer an explicit "Sample" name, else position 1.
  sample_col <- if ("Sample" %in% names(ph)) "Sample" else names(ph)[1L]
  sample_ids <- as.character(ph[[sample_col]])
  trait_cols <- setdiff(names(ph), sample_col)

  std <- data.frame(Sample = sample_ids, stringsAsFactors = FALSE)
  for (tc in trait_cols)
    std[[tc]] <- suppressWarnings(as.numeric(ph[[tc]]))

  geno_samples <- as.character(geno_samples)

  # Drop phenotype rows for samples absent from the genotype set.
  extra <- setdiff(std$Sample, geno_samples)
  if (length(extra) > 0L) {
    warning(caller, ": ", length(extra),
            " phenotype sample(s) absent from the genotype set — dropped: ",
            paste(utils::head(extra, 5L), collapse = ", "),
            if (length(extra) > 5L) " ..." else "")
    std <- std[std$Sample %in% geno_samples, , drop = FALSE]
  }
  if (anyDuplicated(std$Sample))
    stop(caller, ": 'phenotype' has duplicate Sample IDs.")

  # Left-join onto the genotype sample order: genotype samples without a
  # phenotype row get NA traits.
  missing_pheno <- setdiff(geno_samples, std$Sample)
  if (length(missing_pheno) > 0L) {
    warning(caller, ": ", length(missing_pheno),
            " genotype sample(s) lack a phenotype row — filled with NA: ",
            paste(utils::head(missing_pheno, 5L), collapse = ", "),
            if (length(missing_pheno) > 5L) " ..." else "")
    pad <- std[seq_len(length(missing_pheno)), , drop = FALSE]
    pad$Sample <- missing_pheno
    for (tc in setdiff(names(std), "Sample")) pad[[tc]] <- NA_real_
    std <- rbind(std, pad)
  }

  std <- std[match(geno_samples, std$Sample), , drop = FALSE]
  rownames(std) <- NULL
  std
}

# Canonical FAM column names for sample_info (positional contract; see
# docs/geno_meta_structure.md).
.SAMPLE_INFO_COLS <- c("Family", "Sample", "Parent1", "Parent2", "Sex")

# Internal: coerce any fam-like input to the canonical 5-column FAM sample_info,
# aligned (by the Sample column = position 2) to the genotype sample set.
#
#   x            data.frame fam, or NULL (build a dummy table from geno_samples)
#   geno_samples character vector dimnames(geno_arrays)[[1]]
#
# Columns are read BY POSITION (1=Family, 2=Sample, 3=Parent1, 4=Parent2,
# 5=Sex); a one-column input is treated as Sample for back-compatibility.
# Missing/short columns are dummy-filled with "0". The result always carries the
# canonical names so sample_info$Sample stays stable.
.standardise_sample_info <- function(x, geno_samples, caller = "new_geno_meta") {
  n <- length(geno_samples)

  if (is.null(x)) {
    si <- data.frame(
      Family  = rep("0", n),
      Sample  = as.character(geno_samples),
      Parent1 = rep("0", n),
      Parent2 = rep("0", n),
      Sex     = rep("0", n),
      stringsAsFactors = FALSE
    )
    return(si)
  }

  if (!is.data.frame(x))
    stop(caller, ": 'sample_info' must be a data.frame or NULL.")
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  if (ncol(x) == 1L) {
    # Back-compat: a single column is the Sample/IID list.
    si <- data.frame(
      Family  = rep("0", nrow(x)),
      Sample  = as.character(x[[1L]]),
      Parent1 = rep("0", nrow(x)),
      Parent2 = rep("0", nrow(x)),
      Sex     = rep("0", nrow(x)),
      stringsAsFactors = FALSE
    )
  } else {
    # Read up to five columns by POSITION; pad any missing with "0".
    col <- function(i) if (ncol(x) >= i) as.character(x[[i]]) else rep("0", nrow(x))
    si <- data.frame(
      Family  = col(1L),
      Sample  = col(2L),
      Parent1 = col(3L),
      Parent2 = col(4L),
      Sex     = col(5L),
      stringsAsFactors = FALSE
    )
  }

  if (!setequal(si$Sample, geno_samples))
    stop(caller, ": 'sample_info' Sample column (position 2) must match the ",
         "genotype sample set exactly (dimnames(geno_arrays)[[1]]).")

  # Reorder rows to genotype sample order
  si <- si[match(geno_samples, si$Sample), , drop = FALSE]
  rownames(si) <- NULL
  si
}


#' Check whether an object is a valid GenoMeta
#'
#' Performs a lightweight dim-alignment check: verifies that \code{geno_arrays}
#' dim 1 equals \code{nrow(map_info)}, and that \code{block_marker_map_info}
#' markers are a subset of \code{map_info} markers. Does not revalidate
#' the full structure.
#'
#' When the check fails and \code{verbose = TRUE}, a \code{message()} describing
#' the first failing condition is emitted before \code{FALSE} is returned.
#'
#' @param x       Any R object.
#' @param verbose Logical. When \code{TRUE} (default), emit a diagnostic
#'   \code{message()} explaining why \code{x} is not a valid GenoMeta.
#' @return \code{TRUE} if \code{x} appears to be a valid GenoMeta;
#'   \code{FALSE} otherwise.
#'
#' @export
is_geno_meta <- function(x, verbose = TRUE) {
  fail <- function(msg) {
    if (verbose) message("is_geno_meta: ", msg)
    FALSE
  }

  if (!is.list(x))
    return(fail("'x' is not a list."))

  required <- c("geno_arrays","map_info","block_info",
                "block_marker_map_info","sample_info","metadata")
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0L)
    return(fail(paste0("missing required field(s): ",
                       paste(missing_fields, collapse = ", "), ".")))

  if (!is.array(x$geno_arrays) || length(dim(x$geno_arrays)) != 3L)
    return(fail("'geno_arrays' must be a 3-D array."))

  if (!is.data.frame(x$map_info))
    return(fail("'map_info' must be a data.frame."))

  if (dim(x$geno_arrays)[2L] != nrow(x$map_info))
    return(fail(paste0("dim(geno_arrays)[2] (", dim(x$geno_arrays)[2L],
                       ") does not equal nrow(map_info) (",
                       nrow(x$map_info), ").")))

  if (!is.null(x$block_marker_map_info)) {
    bad <- setdiff(x$block_marker_map_info$Marker, x$map_info$Marker)
    if (length(bad) > 0L)
      return(fail(paste0(length(bad),
                         " marker(s) in block_marker_map_info absent from map_info: ",
                         paste(head(bad, 5L), collapse = ", "),
                         if (length(bad) > 5L) " ..." else "", ".")))
  }

  # sample_info: required, must align with the genotype sample set
  if (!is.data.frame(x$sample_info))
    return(fail("'sample_info' must be a data.frame."))
  if (!"Sample" %in% names(x$sample_info))
    return(fail("'sample_info' must contain a 'Sample' column."))
  if (!setequal(as.character(x$sample_info$Sample),
                dimnames(x$geno_arrays)[[1L]]))
    return(fail(paste0("'sample_info$Sample' does not match the genotype ",
                       "sample set (dimnames(geno_arrays)[[1]]).")))

  # marker_effects: always optional (never required for a valid GenoMeta). The
  # canonical form is a named list of per-trait effect tables; a bare data.frame
  # is also accepted. When present, every table must carry Marker/Allele/Effect.
  if (!is.null(x$marker_effects)) {
    me_list <- if (is.data.frame(x$marker_effects)) list(x$marker_effects)
               else x$marker_effects
    if (!is.list(me_list))
      return(fail(paste0("'marker_effects' must be a data.frame, a named list ",
                         "of data.frames, or NULL.")))
    need <- c("Marker", "Allele", "Effect")
    for (me in me_list) {
      if (!is.data.frame(me))
        return(fail("each 'marker_effects' table must be a data.frame."))
      miss <- setdiff(need, names(me))
      if (length(miss) > 0L)
        return(fail(paste0("'marker_effects' table missing column(s): ",
                           paste(miss, collapse = ", "), ".")))
    }
  }

  # phenotype: optional (absent in legacy objects); when present must be a
  # data.frame with a Sample column ⊆ the genotype sample set.
  if ("phenotype" %in% names(x) && !is.null(x$phenotype)) {
    if (!is.data.frame(x$phenotype))
      return(fail("'phenotype' must be a data.frame or NULL."))
    if (!"Sample" %in% names(x$phenotype))
      return(fail("'phenotype' must contain a 'Sample' column."))
    bad_p <- setdiff(as.character(x$phenotype$Sample),
                     dimnames(x$geno_arrays)[[1L]])
    if (length(bad_p) > 0L)
      return(fail(paste0(length(bad_p),
                         " phenotype sample(s) absent from the genotype set: ",
                         paste(utils::head(bad_p, 5L), collapse = ", "),
                         if (length(bad_p) > 5L) " ..." else "", ".")))
  }

  TRUE
}


# ---------------------------------------------------------------------------
# Format readers
# ---------------------------------------------------------------------------

#' Read a Plink 1.9 PED + MAP file pair into a GenoMeta object
#'
#' Parses a \code{.ped} file (6-column header + allele pairs per marker) and
#' its companion \code{.map} file into the unified \code{GenoMeta} structure.
#'
#' @param ped_file  Path to the \code{.ped} file.
#' @param map_file  Path to the \code{.map} file (4-column PLINK format:
#'   chr, snp_id, cM, bp).
#' @param block_file Optional path to a block coordinate file
#'   (\code{Chr, Block, Start, End}) or block-marker map
#'   (\code{Chr, Marker, Block}). When supplied, \code{block_info} and
#'   \code{block_marker_map_info} are populated.
#' @param pos_unit  Position unit in \code{map_file}: \code{"bp"} (default)
#'   or \code{"mbp"}.
#'
#' @return A valid \code{GenoMeta} list.
#'
#' @export
read_plink_ped <- function(ped_file, map_file, block_file = NULL,
                            pos_unit = "bp") {
  ped_file <- as.character(ped_file)
  map_file <- as.character(map_file)

  # Parse MAP file (chr, snp_id, cM, bp) — use process_genotype map path
  raw_map <- utils::read.table(map_file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(raw_map) >= 4L) {
    map_df <- data.frame(
      Chrom  = as.character(raw_map[[1L]]),
      Marker = as.character(raw_map[[2L]]),
      Pos    = as.numeric(raw_map[[4L]]),
      stringsAsFactors = FALSE
    )
  } else if (ncol(raw_map) == 3L) {
    map_df <- data.frame(
      Chrom  = as.character(raw_map[[1L]]),
      Marker = as.character(raw_map[[2L]]),
      Pos    = as.numeric(raw_map[[3L]]),
      stringsAsFactors = FALSE
    )
  } else {
    stop("read_plink_ped: MAP file must have 3 or 4 columns.")
  }

  # Q5: delegate Mbp warning, the preset of map is "PL"
  map_df <- .process_map(map_df, marker_names = map_df$Marker,
                         map_preset = "PL", pos_unit = pos_unit)

  n_markers <- nrow(map_df)

  # Parse PED file
  ped <- utils::read.table(ped_file, header = FALSE, stringsAsFactors = FALSE)
  # Columns 1-6: FID, IID, father, mother, sex, phenotype; then ploidy*n_markers allele cols
  n_allele_cols <- ncol(ped) - 6L
  if (n_allele_cols %% n_markers != 0L)
    stop("read_plink_ped: n_allele_cols (", n_allele_cols,
         ") is not divisible by n_markers (", n_markers,
         "). Cannot determine ploidy.")
  ploidy_detected <- as.integer(n_allele_cols %/% n_markers)

  sample_names <- as.character(ped[[2L]])
  n_samples    <- length(sample_names)

  allele_mat <- as.matrix(ped[, 7L:ncol(ped)])
  mode(allele_mat) <- "character"
  rownames(allele_mat) <- sample_names

  # Collapse ploidy consecutive allele cols → (n_samples × n_markers)
  geno_concat <- .reshape_wide_alleles(allele_mat, n_markers, ploidy_detected)
  colnames(geno_concat) <- map_df$Marker

  split_res <- split_haplotypes(geno_concat,
                                 split_char   = "",
                                 orientation  = 1L,
                                 ploidy       = ploidy_detected,
                                 sample_names = sample_names,
                                 marker_names = map_df$Marker)
  geno_arr <- split_res$array   # (n_samples × n_markers × ploidy_detected)

  # Simple char_dict from observed alleles
  obs <- unique(c(geno_arr))
  obs <- obs[!is.na(obs) & obs != "0"]
  char_dict <- setNames(seq_along(obs), obs)

  blocks <- if (!is.null(block_file)) .read_block_file(block_file, map_df) else list()

  # PED columns 1-5 are the FAM fields (FID, IID, father, mother, sex) — feed
  # them straight into sample_info (positional FAM order).
  fam_si <- data.frame(
    Family  = as.character(ped[[1L]]),
    Sample  = as.character(ped[[2L]]),
    Parent1 = as.character(ped[[3L]]),
    Parent2 = as.character(ped[[4L]]),
    Sex     = as.character(ped[[5L]]),
    stringsAsFactors = FALSE
  )

  new_geno_meta(
    geno_arrays           = geno_arr,
    map_info              = map_df,
    block_info            = blocks$block_info,
    block_marker_map_info = blocks$block_marker_map_info,
    metadata = list(
      ploidy        = ploidy_detected,
      char_dict     = char_dict,
      n_samples     = n_samples,
      n_markers     = n_markers,
      source_format = "plink_ped"
    ),
    sample_info = fam_si
  )
}


#' Read genotype, map, and fam files into structured variables (internal)
#'
#' The single place where genotype files are \emph{physically read}. Given a
#' resolved format parameter set (see \code{\link{resolve_geno_format}}), it
#' loads the genotype table (and optional map), applies the file layout encoded
#' in the parameter set, and returns three plain variables — \code{geno},
#' \code{map}, \code{fam} — \strong{without modifying their contents} (no
#' encoding translation, no position-unit conversion, no marker filtering). All
#' downstream normalisation happens in \code{\link{build_geno_meta}}.
#'
#' The returned \code{geno} is always re-oriented to the fixed internal layout
#' \code{(n_samples × n_markers)} with sample IDs as row names and marker IDs as
#' column names, regardless of the on-disk orientation.
#'
#' @param geno_file Path to the genotype file.
#' @param map_file  Optional path to a map file.
#' @param spec      A resolved format parameter set from
#'   \code{\link{resolve_geno_format}}.
#'
#' @return Named list: \code{geno} (character matrix, samples × markers),
#'   \code{map} (data.frame \code{Chrom, Marker, Pos}, or \code{NULL}), and
#'   \code{fam} (data.frame of the embedded FAM columns, or \code{NULL}).
#' @keywords internal
.read_genotype_files <- function(geno_file, map_file = NULL, spec) {
  geno_file <- as.character(geno_file)
  if (!file.exists(geno_file))
    stop(".read_genotype_files: genotype file not found: ", geno_file)

  raw <- data.table::fread(geno_file, header = spec$geno_header,
                           data.table = FALSE)
  if (ncol(raw) < 2L)
    stop(".read_genotype_files: genotype file must have an identifier column ",
         "plus at least one value column.")

  fam <- NULL
  if (spec$fam_cols > 0L) {
    # PLINK-style embedded FAM block: columns 1:fam_cols are pedigree fields;
    # genotype values begin at column fam_cols + 1. Sample IDs are the IID
    # column (position 2 of the FAM block, or position 1 when fam_cols == 1).
    if (ncol(raw) <= spec$fam_cols)
      stop(".read_genotype_files: file has ", ncol(raw), " column(s) but the ",
           "preset reserves the first ", spec$fam_cols, " as FAM fields — no ",
           "genotype columns remain.")
    fam        <- raw[, seq_len(spec$fam_cols), drop = FALSE]
    iid_col    <- if (spec$fam_cols >= 2L) 2L else 1L
    sample_ids <- as.character(raw[[iid_col]])
    geno_2d    <- as.matrix(raw[, (spec$fam_cols + 1L):ncol(raw), drop = FALSE])
    mode(geno_2d)     <- "character"
    rownames(geno_2d) <- sample_ids
    # Column names (marker IDs) come straight from the header.
  } else {
    row_ids   <- as.character(raw[[1L]])
    geno_2d   <- as.matrix(raw[, -1L, drop = FALSE])
    mode(geno_2d)     <- "character"
    rownames(geno_2d) <- row_ids
  }

  # Re-orient to the FIXED internal layout (rows = samples, cols = markers).
  if (spec$geno_orient == "marker") geno_2d <- t(geno_2d)

  map <- if (!is.null(map_file)) {
    .read_map_file(map_file, spec$map_order)
  } else {
    NULL
  }

  list(geno = geno_2d, map = map, fam = fam)
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
.read_map_file <- function(map_file, map_order) {
  map_file <- as.character(map_file)
  if (!file.exists(map_file))
    stop(".read_map_file: map file not found: ", map_file)
  raw_map <- data.table::fread(map_file, data.table = FALSE)
  if (ncol(raw_map) < 3L)
    stop(".read_map_file: map file must have at least three columns (roles: ",
         paste(map_order, collapse = ", "), ").")
  pos_of <- function(role) which(map_order == role)[1L]
  data.frame(
    Chrom  = as.character(raw_map[[pos_of("Chrom")]]),
    Marker = as.character(raw_map[[pos_of("Marker")]]),
    Pos    = raw_map[[pos_of("Pos")]],
    stringsAsFactors = FALSE
  )
}


#' Read a delimited genotype table into a GenoMeta object
#'
#' Unified reader for delimited genotype/map file pairs. The on-disk file layout
#' is described by a \code{preset} parameter set (see
#' \code{\link{GENO_FORMAT_PRESETS}}); any field can be overridden, or a fully
#' custom parameter set can be supplied. File reading is delegated to
#' \code{\link{.read_genotype_files}} (which returns the genotype, map, and fam
#' unchanged), and assembly to \code{\link{build_geno_meta}} (which builds the
#' \code{(n_samples, n_markers, ploidy)} array, intersects markers with the map,
#' resolves the encoding dictionary, attaches blocks, and validates the result).
#'
#' PLINK PED (\code{\link{read_plink_ped}}), DArTseq
#' (\code{\link{read_dart_genotype}}), and VCF (\code{\link{read_vcf_genotype}})
#' have their own dedicated readers.
#'
#' @param geno_file  Path to the genotype file.
#' @param map_file   Optional path to a map file. The file column order is taken
#'   from the preset's \code{map_order} field.
#' @param preset     File-format preset name (\code{"ST"}, \code{"GS"},
#'   \code{"PL"}, \code{"default"}), a named list giving a custom parameter set,
#'   or \code{NULL} for the \code{default} layout. See
#'   \code{\link{GENO_FORMAT_PRESETS}}.
#' @param block_file Optional path to a block coordinate or block-marker file.
#' @param presets    Preset registry to resolve \code{preset} against. Default
#'   \code{\link{GENO_FORMAT_PRESETS}}.
#' @param split_char Optional override for the allele delimiter (\code{"/"} or
#'   \code{""}).
#' @param encoding   Optional override for the allele encoding
#'   (\code{"selectiontools"} or \code{"nucleotide"}).
#' @param ploidy     Integer ploidy level (default \code{2}).
#' @param pos_unit   Position unit in \code{map_file}: \code{"bp"} (default) or
#'   \code{"mbp"}.
#'
#' @return A valid \code{GenoMeta} list.
#'
#' @seealso \code{\link{resolve_geno_format}}, \code{\link{build_geno_meta}},
#'   \code{\link{read_st_genotype}}, \code{\link{read_gs_genotype}},
#'   \code{\link{read_plink_ped}}
#'
#' @export
read_genotype <- function(geno_file,
                          map_file   = NULL,
                          preset     = NULL,
                          block_file = NULL,
                          presets    = GENO_FORMAT_PRESETS,
                          split_char = NULL,
                          encoding   = NULL,
                          ploidy     = 2L,
                          pos_unit   = "bp") {

  spec  <- resolve_geno_format(preset, presets = presets,
                               split_char = split_char, encoding = encoding)
  parts <- .read_genotype_files(geno_file, map_file, spec)

  # Synthesise a placeholder map when none is supplied (positions unknown).
  map <- parts$map
  if (is.null(map))
    map <- data.frame(Chrom = NA_character_, Marker = colnames(parts$geno),
                      Pos = NA_real_, stringsAsFactors = FALSE)

  build_geno_meta(
    geno_2d       = parts$geno,
    map           = map,
    fam           = parts$fam,
    blocks        = block_file,
    ploidy        = as.integer(ploidy),
    split_char    = spec$split_char,
    encoding      = spec$encoding,
    source_format = spec$source_format,
    pos_unit      = pos_unit
  )
}


#' Read a SelectionTools genotype file into a GenoMeta object
#'
#' Thin wrapper over \code{\link{read_genotype}} with \code{preset = "ST"}
#' (markers as rows, slash-separated numeric cells, \code{name | chr | pos} map).
#'
#' @param file       Path to the SelectionTools genotype file.
#' @param map_file   Optional path to a map file (\code{name | chr | pos}).
#' @param block_file Optional path to a block coordinate or block-marker file.
#' @param pos_unit   Position unit: \code{"bp"} (default) or \code{"mbp"}.
#'
#' @return A valid \code{GenoMeta} list.
#'
#' @export
read_st_genotype <- function(file, map_file = NULL, block_file = NULL,
                              pos_unit = "bp") {
  read_genotype(file, map_file = map_file, preset = "ST",
                block_file = block_file, pos_unit = pos_unit)
}


#' Read a genomicSimulation genotype file into a GenoMeta object
#'
#' Thin wrapper over \code{\link{read_genotype}} with \code{preset = "GS"}
#' (markers as rows, two-character allele cells, e.g. \code{"AA"}).
#'
#' @param file       Path to the genomicSimulation genotype file.
#' @param map_file   Optional path to a map file (\code{name | chr | pos}).
#' @param block_file Optional path to a block coordinate or block-marker file.
#' @param pos_unit   Position unit: \code{"bp"} (default) or \code{"mbp"}.
#'
#' @return A valid \code{GenoMeta} list.
#'
#' @export
read_gs_genotype <- function(file, map_file = NULL, block_file = NULL,
                              pos_unit = "bp") {
  read_genotype(file, map_file = map_file, preset = "GS",
                block_file = block_file, pos_unit = pos_unit)
}


#' Read a DArTseq SNP report into a GenoMeta object
#'
#' Parses a DArTseq 1-row SNP report (wide format: one row per marker, one
#' column per sample, dosage values 0/1/2) into the unified \code{GenoMeta}
#' structure. Marker metadata columns are detected automatically, and the
#' dosage is split into two haplotype slices for consistency with the rest of
#' the haploStack pipeline (2 → ["A","A"], 1 → ["A","B"], 0 → ["B","B"]).
#'
#' @section DArT report format:
#' The file must be a CSV or tab-delimited text file. Multiple leading metadata
#' rows (rows before the marker data) are handled via \code{n_meta_rows} or
#' auto-detection. Sample columns are identified as the first contiguous block
#' of columns containing only numeric values (0, 1, 2, or NA).
#'
#' @param file               Path to the DArTseq SNP report file.
#' @param n_meta_rows        Integer number of leading metadata rows to skip,
#'   or \code{NULL} (default) to auto-detect by scanning for the row where
#'   non-metadata columns first become all-numeric.
#' @param marker_col         Column name used as the marker ID. Searched in
#'   order: \code{"AlleleID"}, \code{"TargetID"}, \code{"CloneID"},
#'   \code{"MarkerName"}, \code{"Marker"}. Or supply a column name directly.
#'   Default \code{NULL} = auto.
#' @param chrom_col          Column name for chromosome. Searched in order:
#'   first column matching \code{"[Cc]hrom"} but not \code{"[Pp]os"}.
#'   Default \code{NULL} = auto.
#' @param pos_col            Column name for SNP position. Searched in order:
#'   first column matching \code{"[Pp]os"}, \code{"ChromPos"}, \code{"SnpPos"}.
#'   Default \code{NULL} = auto.
#' @param sample_start_col   Integer index of the first sample column. When
#'   \code{NULL} (default), auto-detected as the first column after all
#'   character-type marker-metadata columns.
#' @param block_file         Optional path to a block coordinate or
#'   block-marker file.
#' @param pos_unit           Position unit: \code{"bp"} (default) or
#'   \code{"mbp"}.
#' @param call_rate_threshold Numeric (0–1). Markers with a call rate below
#'   this threshold are removed. Default \code{0.8}.
#' @param maf_threshold      Numeric (0–1). Markers with minor allele
#'   frequency below this threshold are removed. Default \code{0.01}.
#'
#' @return A valid \code{GenoMeta} list with:
#' \describe{
#'   \item{\code{geno_arrays}}{3-D character array \code{[n_samples × n_markers × 2]}.
#'     Allele encoding: \code{"A"} = reference, \code{"B"} = alternative.
#'     \code{NA} = missing call.}
#'   \item{\code{map_info}}{data.frame \code{(Chr, Marker, Pos)}.}
#'   \item{\code{metadata$source_format}}{\code{"dart_seq"}.}
#' }
#'
#' @seealso \code{\link{read_st_genotype}}, \code{\link{read_plink_ped}},
#'   \code{\link{filter_geno_meta}}, \code{\link{compute_grm}}
#'
#' @export
read_dart_genotype <- function(file,
                                n_meta_rows         = NULL,
                                marker_col          = NULL,
                                chrom_col           = NULL,
                                pos_col             = NULL,
                                sample_start_col    = NULL,
                                block_file          = NULL,
                                pos_unit            = "bp",
                                call_rate_threshold = 0.8,
                                maf_threshold       = 0.01) {

  file <- as.character(file)
  if (!file.exists(file))
    stop("read_dart_genotype: file not found: ", file)

  # ---- Read raw file -------------------------------------------------------
  sep <- if (grepl("[.]csv$", file, ignore.case = TRUE)) "," else "\t"
  raw <- utils::read.table(file, header = TRUE, sep = sep,
                           stringsAsFactors = FALSE, fill = TRUE,
                           check.names = FALSE, comment.char = "")

  if (nrow(raw) == 0L || ncol(raw) < 2L)
    stop("read_dart_genotype: file appears empty or has fewer than 2 columns.")

  # ---- Handle leading metadata rows ----------------------------------------
  # DArT reports sometimes have extra rows at the top (callrate rows, etc.)
  # before the actual per-marker data. Auto-detect by finding the first row
  # where the majority of sample-candidate columns are numeric 0/1/2.
  if (!is.null(n_meta_rows) && n_meta_rows > 0L) {
    colnames(raw) <- as.character(unlist(raw[n_meta_rows, ]))
    raw <- raw[(n_meta_rows + 1L):nrow(raw), , drop = FALSE]
    rownames(raw) <- NULL
  }

  # ---- Identify marker metadata vs sample columns --------------------------
  # Heuristic: marker metadata columns are character; sample columns are
  # coercible to 0/1/2/NA. Walk from right until we find a non-numeric column.
  is_numeric_col <- vapply(seq_len(ncol(raw)), function(j) {
    v <- suppressWarnings(as.numeric(as.character(raw[[j]])))
    non_na <- v[!is.na(v)]
    length(non_na) > 0L && all(non_na %in% c(0, 1, 2))
  }, logical(1L))

  if (!is.null(sample_start_col)) {
    first_sample <- as.integer(sample_start_col)
  } else {
    # Find the first column that looks like a sample column; require that all
    # subsequent columns also pass the numeric test.
    first_sample <- NULL
    for (j in seq_len(ncol(raw))) {
      if (all(is_numeric_col[j:ncol(raw)])) {
        first_sample <- j
        break
      }
    }
    if (is.null(first_sample))
      stop("read_dart_genotype: could not auto-detect sample columns. ",
           "Please supply 'sample_start_col'.")
  }

  meta_cols_df <- raw[, seq_len(first_sample - 1L), drop = FALSE]
  geno_df      <- raw[, first_sample:ncol(raw), drop = FALSE]
  sample_names <- colnames(geno_df)

  # ---- Identify marker ID column -------------------------------------------
  candidate_marker_cols <- c("AlleleID","TargetID","CloneID","MarkerName","Marker")
  if (!is.null(marker_col)) {
    mk_col <- marker_col
  } else {
    mk_col <- NULL
    for (cn in candidate_marker_cols) {
      if (cn %in% names(meta_cols_df)) { mk_col <- cn; break }
    }
    if (is.null(mk_col)) mk_col <- names(meta_cols_df)[1L]
  }
  marker_names <- as.character(meta_cols_df[[mk_col]])

  # Deduplicate marker names (rare but possible)
  if (anyDuplicated(marker_names)) {
    marker_names <- make.unique(marker_names, sep = "_")
    warning("read_dart_genotype: duplicate marker IDs found; made unique.")
  }

  # ---- Identify chromosome and position columns ----------------------------
  all_meta_cols <- names(meta_cols_df)

  if (!is.null(chrom_col)) {
    chr_col_name <- chrom_col
  } else {
    chr_col_name <- NULL
    for (cn in all_meta_cols) {
      if (grepl("[Cc]hrom", cn) && !grepl("[Pp]os", cn)) {
        chr_col_name <- cn; break
      }
    }
  }

  if (!is.null(pos_col)) {
    pos_col_name <- pos_col
  } else {
    pos_col_name <- NULL
    pos_patterns <- c("ChromPosSnp", "SnpPosition", "ChromPos",
                      "[Pp]os", "[Pp]osition")
    for (pat in pos_patterns) {
      for (cn in all_meta_cols) {
        if (grepl(pat, cn)) { pos_col_name <- cn; break }
      }
      if (!is.null(pos_col_name)) break
    }
  }

  chr_vec <- if (!is.null(chr_col_name) && chr_col_name %in% all_meta_cols) {
    # Normalise H-suffix (1H → 1) common in barley
    gsub("H$", "", as.character(meta_cols_df[[chr_col_name]]))
  } else {
    rep(NA_character_, nrow(meta_cols_df))
  }

  pos_vec <- if (!is.null(pos_col_name) && pos_col_name %in% all_meta_cols) {
    suppressWarnings(as.numeric(meta_cols_df[[pos_col_name]]))
  } else {
    rep(NA_real_, nrow(meta_cols_df))
  }

  # ---- Build numeric dosage matrix -----------------------------------------
  n_all_markers <- length(marker_names)
  n_samples     <- length(sample_names)

  dosage <- matrix(NA_real_, nrow = n_all_markers, ncol = n_samples,
                   dimnames = list(marker_names, sample_names))
  for (j in seq_len(n_samples)) {
    v <- suppressWarnings(as.numeric(as.character(geno_df[[j]])))
    # DArT missing: -1, blank strings, or literal NA
    v[!is.na(v) & v < 0L] <- NA_real_
    dosage[, j] <- v
  }

  # ---- Filter by call rate --------------------------------------------------
  call_rate <- rowMeans(!is.na(dosage))
  pass_cr   <- call_rate >= call_rate_threshold
  if (sum(pass_cr) < n_all_markers)
    message("read_dart_genotype: ", sum(!pass_cr),
            " marker(s) removed (call rate < ", call_rate_threshold, ").")

  dosage       <- dosage[pass_cr, , drop = FALSE]
  chr_vec      <- chr_vec[pass_cr]
  pos_vec      <- pos_vec[pass_cr]
  marker_names <- rownames(dosage)
  n_markers    <- nrow(dosage)

  # ---- Filter by MAF -------------------------------------------------------
  freq <- rowMeans(dosage, na.rm = TRUE) / 2   # allele freq of alt allele
  maf  <- pmin(freq, 1 - freq)
  pass_maf <- !is.na(maf) & maf >= maf_threshold
  if (sum(pass_maf) < n_markers)
    message("read_dart_genotype: ", sum(!pass_maf),
            " marker(s) removed (MAF < ", maf_threshold, ").")

  dosage       <- dosage[pass_maf, , drop = FALSE]
  chr_vec      <- chr_vec[pass_maf]
  pos_vec      <- pos_vec[pass_maf]
  marker_names <- rownames(dosage)
  n_markers    <- nrow(dosage)

  if (n_markers == 0L)
    stop("read_dart_genotype: no markers remain after filtering.")

  # ---- Split dosage to two haplotype slices --------------------------------
  # 2 (AA) → hap1 = "A", hap2 = "A"
  # 1 (AB) → hap1 = "A", hap2 = "B"
  # 0 (BB) → hap1 = "B", hap2 = "B"
  # NA     → hap1 = NA,  hap2 = NA
  hap1 <- matrix(ifelse(dosage == 0L, "B", "A"),
                 nrow = n_markers, ncol = n_samples)
  hap2 <- matrix(ifelse(dosage >= 2L, "A", "B"),
                 nrow = n_markers, ncol = n_samples)
  hap1[is.na(dosage)] <- NA_character_
  hap2[is.na(dosage)] <- NA_character_

  # [n_samples × n_markers × ploidy=2]; hap1/hap2 are [n_markers × n_samples]
  geno_arr <- array(NA_character_,
                    dim      = c(n_samples, n_markers, 2L),
                    dimnames = list(sample_names, marker_names, c("1","2")))
  geno_arr[, , 1L] <- t(hap1)
  geno_arr[, , 2L] <- t(hap2)

  # ---- Build map_info -------------------------------------------------------
  map_df <- data.frame(
    Chrom  = chr_vec,
    Marker = marker_names,
    Pos    = pos_vec,
    stringsAsFactors = FALSE
  )
  map_df <- .process_map(map_df, marker_names = marker_names,
                         map_preset = NULL, pos_unit = pos_unit)

  # ---- Block files ----------------------------------------------------------
  blocks <- if (!is.null(block_file)) .read_block_file(block_file, map_df) else list()

  # ---- Assemble GenoMeta ---------------------------------------------------
  char_dict <- list("A" = 1L, "B" = 2L)

  new_geno_meta(
    geno_arrays           = geno_arr,
    map_info              = map_df,
    block_info            = blocks$block_info,
    block_marker_map_info = blocks$block_marker_map_info,
    metadata = list(
      ploidy        = 2L,
      char_dict     = char_dict,
      n_samples     = n_samples,
      n_markers     = n_markers,
      source_format = "dart_seq"
    )
  )
}


#' Read a VCF file into a GenoMeta object
#'
#' Parses a standard VCF (Variant Call Format) file into the unified
#' \code{GenoMeta} structure using the \pkg{vcfR} package. Genotype calls in
#' the \code{GT} FORMAT field (phased or unphased) are converted to dosage
#' (0/1/2) and split into two haplotype slices
#' (2 -> [\code{"A"},\code{"A"}], 1 -> [\code{"A"},\code{"B"}],
#' 0 -> [\code{"B"},\code{"B"}]).
#'
#' @section VCF requirements:
#' \itemize{
#'   \item VCF version 4.0+ with a \code{GT} FORMAT field.
#'   \item Diploid samples only (exactly two allele calls per GT).
#'   \item Multi-allelic sites are filtered to biallelic only when
#'     \code{biallelic_only = TRUE} (default).
#' }
#'
#' @param file               Path to the VCF file (plain or bgzipped).
#' @param block_file         Optional path to a block coordinate or
#'   block-marker file.
#' @param pos_unit           Position unit: \code{"bp"} (default) or
#'   \code{"mbp"}.
#' @param call_rate_threshold Numeric (0-1). Markers below this call rate are
#'   removed. Default \code{0.8}.
#' @param maf_threshold      Numeric (0-1). Markers below this minor allele
#'   frequency are removed. Default \code{0.01}.
#' @param biallelic_only     Logical. Drop multi-allelic sites. Default
#'   \code{TRUE}.
#'
#' @return A valid \code{GenoMeta} list with \code{metadata$source_format =
#'   "vcf"} and A/B haplotype encoding.
#'
#' @seealso \code{\link{read_dart_genotype}}, \code{\link{read_plink_ped}},
#'   \code{\link{filter_geno_meta}}, \code{\link{compute_grm}}
#'
#' @export
read_vcf_genotype <- function(file,
                               block_file          = NULL,
                               pos_unit            = "bp",
                               call_rate_threshold = 0.8,
                               maf_threshold       = 0.01,
                               biallelic_only      = TRUE) {

  if (!file.exists(file))
    stop("read_vcf_genotype: file not found: ", file)
  if (!requireNamespace("vcfR", quietly = TRUE))
    stop("read_vcf_genotype requires the 'vcfR' package. ",
         "Install with: install.packages('vcfR')")

  vcf <- vcfR::read.vcfR(file, verbose = FALSE)

  # ---- Extract fixed fields (CHROM, POS, ID, ALT) --------------------------
  fix_df <- as.data.frame(vcfR::getFIX(vcf), stringsAsFactors = FALSE)
  colnames(fix_df) <- toupper(colnames(fix_df))

  if (!all(c("CHROM","POS") %in% names(fix_df)))
    stop("read_vcf_genotype: VCF is missing CHROM or POS columns.")

  # ---- Biallelic filter (before GT extraction) -----------------------------
  if (biallelic_only) {
    alt_col      <- fix_df[["ALT"]]
    is_biallelic <- !grepl(",", alt_col, fixed = TRUE) &
                    !is.na(alt_col) & alt_col != "."
    n_multi <- sum(!is_biallelic)
    if (n_multi > 0L)
      message("read_vcf_genotype: removed ", n_multi, " multi-allelic site(s).")
    vcf    <- vcf[is_biallelic, ]
    fix_df <- fix_df[is_biallelic, , drop = FALSE]
  }

  if (nrow(fix_df) == 0L)
    stop("read_vcf_genotype: no sites remaining after biallelic filter.")

  # ---- Build marker IDs ----------------------------------------------------
  raw_ids    <- fix_df[["ID"]]
  missing_id <- is.na(raw_ids) | raw_ids == "." | raw_ids == ""
  raw_ids[missing_id] <- paste0(fix_df[["CHROM"]][missing_id], ":",
                                 fix_df[["POS"]][missing_id])
  marker_ids <- make.unique(raw_ids, sep = "_dup")

  # ---- Clean chromosome names ----------------------------------------------
  chrom_raw   <- fix_df[["CHROM"]]
  chrom_clean <- gsub("^[Cc]hr", "", chrom_raw)   # strip chr/Chr prefix
  chrom_clean <- gsub("H$",      "", chrom_clean)  # strip barley H-suffix

  pos_num <- suppressWarnings(as.numeric(fix_df[["POS"]]))
  if (tolower(pos_unit) == "mbp") pos_num <- pos_num * 1e6

  # ---- Extract GT matrix ---------------------------------------------------
  gt_mat <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)
  if (is.null(gt_mat) || nrow(gt_mat) == 0L)
    stop("read_vcf_genotype: no GT data found in VCF.")

  rownames(gt_mat) <- marker_ids

  # ---- Parse GT strings to dosage (count of alt alleles) ------------------
  .gt_to_dosage <- function(gt_vec) {
    gt_clean <- gsub("|", "/", gt_vec, fixed = TRUE)  # phased to unphased
    parts    <- strsplit(gt_clean, "/", fixed = TRUE)
    vapply(parts, function(x) {
      if (length(x) != 2L || any(x == ".")) return(NA_integer_)
      a1 <- suppressWarnings(as.integer(x[1L]))
      a2 <- suppressWarnings(as.integer(x[2L]))
      if (is.na(a1) || is.na(a2)) return(NA_integer_)
      pmin(a1 + a2, 2L)
    }, integer(1L))
  }

  dos_mat <- matrix(
    .gt_to_dosage(as.vector(gt_mat)),
    nrow     = nrow(gt_mat),
    ncol     = ncol(gt_mat),
    dimnames = list(marker_ids, colnames(gt_mat))
  )

  n_mk  <- nrow(dos_mat)
  n_smp <- ncol(dos_mat)
  if (n_mk  == 0L) stop("read_vcf_genotype: no markers after GT parsing.")
  if (n_smp == 0L) stop("read_vcf_genotype: no samples found in VCF.")

  sample_names <- colnames(dos_mat)

  # ---- Call rate filter ----------------------------------------------------
  call_rate <- rowMeans(!is.na(dos_mat))
  keep_cr   <- call_rate >= call_rate_threshold
  n_fail_cr <- sum(!keep_cr)
  if (n_fail_cr > 0L)
    message("read_vcf_genotype: removed ", n_fail_cr,
            " marker(s) below call rate ", call_rate_threshold, ".")
  dos_mat     <- dos_mat[keep_cr, , drop = FALSE]
  chrom_clean <- chrom_clean[keep_cr]
  pos_num     <- pos_num[keep_cr]

  # ---- MAF filter ----------------------------------------------------------
  p          <- rowMeans(dos_mat, na.rm = TRUE) / 2
  maf        <- pmin(p, 1 - p)
  keep_maf   <- !is.na(maf) & maf >= maf_threshold
  n_fail_maf <- sum(!keep_maf)
  if (n_fail_maf > 0L)
    message("read_vcf_genotype: removed ", n_fail_maf,
            " marker(s) below MAF ", maf_threshold, ".")
  dos_mat     <- dos_mat[keep_maf, , drop = FALSE]
  chrom_clean <- chrom_clean[keep_maf]
  pos_num     <- pos_num[keep_maf]

  n_mk_final <- nrow(dos_mat)
  if (n_mk_final == 0L)
    stop("read_vcf_genotype: no markers remaining after call rate / MAF filtering.")

  marker_names <- rownames(dos_mat)

  # ---- Mean-impute NA genotypes --------------------------------------------
  na_rows <- which(rowSums(is.na(dos_mat)) > 0L)
  if (length(na_rows) > 0L) {
    message("read_vcf_genotype: mean-imputing NA genotypes in ",
            length(na_rows), " marker(s).")
    row_means <- rowMeans(dos_mat, na.rm = TRUE)
    for (r in na_rows) {
      na_cols <- which(is.na(dos_mat[r, ]))
      dos_mat[r, na_cols] <- as.integer(round(row_means[r]))
    }
  }

  # ---- Split dosage to A/B haplotype slices --------------------------------
  arr <- array(
    NA_character_,
    dim      = c(n_smp, n_mk_final, 2L),
    dimnames = list(sample_names, marker_names, c("H1", "H2"))
  )
  arr[, , 1L] <- t(ifelse(dos_mat == 0L, "B", "A"))  # hap1: "B" only if dosage=0
  arr[, , 2L] <- t(ifelse(dos_mat >= 2L, "A", "B"))  # hap2: "A" only if dosage=2

  # ---- Map info & block files ----------------------------------------------
  map_df <- data.frame(
    Chrom  = chrom_clean,
    Marker = marker_names,
    Pos    = pos_num,
    stringsAsFactors = FALSE
  )

  blocks <- if (!is.null(block_file)) .read_block_file(block_file, map_df) else list()

  # ---- Assemble GenoMeta ---------------------------------------------------
  new_geno_meta(
    geno_arrays           = arr,
    map_info              = map_df,
    block_info            = blocks$block_info,
    block_marker_map_info = blocks$block_marker_map_info,
    metadata = list(
      ploidy        = 2L,
      char_dict     = list("A" = 1L, "B" = 2L),
      n_samples     = n_smp,
      n_markers     = n_mk_final,
      source_format = "vcf"
    )
  )
}


#' Build a GenoMeta object from fixed-shape in-memory variables
#'
#' Assembles a valid \code{GenoMeta} from three primary elements with a
#' \strong{fixed input contract}: a genotype matrix, a map, and an optional fam
#' table. This is the single assembly entry point shared by
#' \code{\link{read_genotype}} and by callers holding in-memory data. The file
#' layout — orientation, headers, embedded FAM columns, map column order — is
#' resolved \emph{before} this point (in the reader); \code{build_geno_meta()}
#' itself takes only the changeable \emph{format} parameters \code{split_char},
#' \code{pos_unit}, and \code{encoding}.
#'
#' @section Fixed input contract:
#' \describe{
#'   \item{\code{geno_2d}}{Genotype matrix in \code{(n_samples × n_markers)}
#'     orientation only — rows are samples, columns are markers. Row names must
#'     be sample IDs and column names must be marker IDs. Orientation is
#'     \strong{not} auto-detected.}
#'   \item{\code{map}}{At least three columns in the order
#'     \code{Chrom, Marker, Pos}.}
#'   \item{\code{fam}}{Columns in the order \code{Family, Sample, Parent1,
#'     Parent2, Sex} (the \code{Sample} column, position 2, identifies samples).}
#' }
#'
#' Marker and sample name intersections are resolved internally with informative
#' warnings when names are dropped.
#'
#' @param geno_2d      Character matrix or data.frame of raw genotype calls in
#'   \code{(n_samples × n_markers)} orientation with named rows (samples) and
#'   columns (markers). Cells may be slash-separated haplotypes (\code{"1/4"}),
#'   two-character concatenations (\code{"AA"}), or single codes.
#' @param map          data.frame whose first three columns are, in order,
#'   \code{Chrom}, \code{Marker}, \code{Pos}.
#' @param fam          Optional fam data.frame with columns in the order
#'   \code{Family}, \code{Sample}, \code{Parent1}, \code{Parent2}, \code{Sex}.
#'   The \code{Sample} column (position 2) supplies the sample IDs; a
#'   single-column fam is treated as the sample list. Samples absent from
#'   \code{fam} are dropped with a warning. \code{NULL} (default) keeps all
#'   samples and builds a dummy fam.
#' @param blocks       Optional block source (data.frame or file path) passed to
#'   \code{\link{attach_block_info}}. \code{NULL} (default) skips block
#'   attachment.
#' @param ploidy       Integer ploidy. Default \code{2L}.
#' @param split_char   Allele separator used in genotype cells: \code{"/"}
#'   (default) or \code{""} (one character per allele).
#' @param pos_unit     Position unit of \code{map$Pos}: \code{"bp"} (default) or
#'   \code{"mbp"}. Internally all positions are stored in bp.
#' @param encoding     Allele encoding: \code{"selectiontools"} or
#'   \code{"nucleotide"}. \code{NULL} (default) auto-detects.
#' @param source_format Label stored in \code{metadata$source_format}.
#'   Default \code{"user_supplied"}.
#'
#' @return A valid \code{GenoMeta} list (see \code{\link{new_geno_meta}}).
#'
#' @details
#' \strong{Automatic marker filtering.} \code{build_geno_meta()} intersects the
#' genotype markers with \code{map$Marker} and keeps only the shared set — no
#' manual pre-filtering is needed. Genotype markers absent from the map are
#' dropped (with a warning), map markers absent from the genotypes are excluded,
#' and a complete mismatch (no shared markers) raises an error. Every reader
#' (\code{\link{read_genotype}} and its wrappers) delegates here, so this is the
#' single place marker intersection happens. An optional \code{fam} table applies
#' the analogous automatic intersection on the sample axis.
#'
#' @examples
#' \dontrun{
#' # geno_2d: rows = samples, cols = markers (fixed orientation)
#' geno <- matrix(c("1/1", "1/4", "4/4",
#'                  "4/4", "1/4", "1/1"),
#'                nrow = 2L, byrow = TRUE,
#'                dimnames = list(c("S1", "S2"), c("M1", "M2", "M3")))
#' # map: columns in order Chrom, Marker, Pos
#' map  <- data.frame(Chrom = c("1", "1", "2"),
#'                    Marker = c("M1", "M2", "M3"),
#'                    Pos = c(1e6, 2e6, 3e6), stringsAsFactors = FALSE)
#' meta <- build_geno_meta(geno, map)
#' is_geno_meta(meta)  # TRUE
#' }
#'
#' @seealso \code{\link{new_geno_meta}}, \code{\link{attach_block_info}},
#'   \code{\link{read_genotype}}, \code{\link{resolve_geno_format}}
#'
#' @export
build_geno_meta <- function(geno_2d,
                             map,
                             fam           = NULL,
                             blocks        = NULL,
                             ploidy        = 2L,
                             split_char    = "/",
                             pos_unit      = "bp",
                             encoding      = NULL,
                             source_format = "user_supplied") {

  # ---- Validate genotype matrix (FIXED sample × marker orientation) ----------
  if (!is.matrix(geno_2d) && !is.data.frame(geno_2d))
    stop("build_geno_meta: 'geno_2d' must be a matrix or data.frame in ",
         "(n_samples × n_markers) orientation.")
  geno_2d <- as.matrix(geno_2d)
  mode(geno_2d) <- "character"
  if (is.null(rownames(geno_2d)) || is.null(colnames(geno_2d)))
    stop("build_geno_meta: 'geno_2d' must have row names (sample IDs) and ",
         "column names (marker IDs). Orientation is fixed: rows = samples, ",
         "columns = markers. Re-orient the input before calling.")
  sample_names <- rownames(geno_2d)
  marker_names <- colnames(geno_2d)

  # ---- Validate / normalise map (FIXED column order: Chrom, Marker, Pos) -----
  if (missing(map) || is.null(map))
    stop("build_geno_meta: 'map' is required (columns in order ",
         "Chrom, Marker, Pos).")
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  if (ncol(map) < 3L)
    stop("build_geno_meta: 'map' must have at least three columns in the ",
         "order Chrom, Marker, Pos.")
  # Positional contract: the first three columns ARE Chrom, Marker, Pos. Honour
  # an already-canonical naming, otherwise rename the first three by position.
  if (all(c("Chrom", "Marker", "Pos") %in% names(map))) {
    map <- map[, c("Chrom", "Marker", "Pos"), drop = FALSE]
  } else {
    map <- map[, seq_len(3L), drop = FALSE]
    names(map) <- c("Chrom", "Marker", "Pos")
  }
  # Convert position units (→ bp) and emit the Mbp-scale warning (Q5).
  map <- .process_map(map, marker_names = map$Marker, map_preset = NULL,
                      pos_unit = pos_unit)
  map$Chrom  <- as.character(map$Chrom)
  map$Marker <- as.character(map$Marker)

  # ---- Marker intersection (canonical auto-filter: geno ∩ map) ---------------
  # Single source of truth for marker filtering — every reader delegates here,
  # so callers never need to pre-intersect geno_2d and map by hand.
  map_markers <- map$Marker
  extra_geno  <- setdiff(marker_names, map_markers)
  extra_map   <- setdiff(map_markers,  marker_names)
  if (length(extra_geno) > 0L)
    warning("build_geno_meta: ", length(extra_geno),
            " marker(s) in geno_2d not found in map and will be dropped: ",
            paste(head(extra_geno, 5L), collapse = ", "),
            if (length(extra_geno) > 5L) " ..." else "")
  if (length(extra_map) > 0L)
    warning("build_geno_meta: ", length(extra_map),
            " marker(s) in map not found in geno_2d (absent from array): ",
            paste(head(extra_map, 5L), collapse = ", "),
            if (length(extra_map) > 5L) " ..." else "")
  keep_mk <- intersect(marker_names, map_markers)
  if (length(keep_mk) == 0L)
    stop("build_geno_meta: no markers in common between geno_2d and map.")

  # ---- Sample intersection (optional fam) → becomes sample_info --------------
  # The fam is a positional FAM table (Family, Sample, Parent1, Parent2, Sex);
  # the Sample column is read by position (col 2, or col 1 if single-column).
  keep_smp <- sample_names
  fam_si   <- NULL
  if (!is.null(fam)) {
    fam <- as.data.frame(fam, stringsAsFactors = FALSE)
    if (ncol(fam) < 1L)
      stop("build_geno_meta: 'fam' must have at least one column (Sample).")
    sample_col <- if (ncol(fam) == 1L) 1L else 2L
    fam_ids    <- as.character(fam[[sample_col]])
    extra_smp  <- setdiff(sample_names, fam_ids)
    if (length(extra_smp) > 0L)
      warning("build_geno_meta: ", length(extra_smp),
              " sample(s) in geno_2d not found in fam and will be dropped: ",
              paste(head(extra_smp, 5L), collapse = ", "),
              if (length(extra_smp) > 5L) " ..." else "")
    keep_smp <- intersect(sample_names, fam_ids)
    if (length(keep_smp) == 0L)
      stop("build_geno_meta: no samples in common between geno_2d and fam.")
    # Retain fam rows for the kept samples; new_geno_meta() standardises and
    # aligns them to the genotype array order.
    fam_si <- fam[fam_ids %in% keep_smp, , drop = FALSE]
  }

  # ---- Subset and split to 3D array ------------------------------------------
  # Orientation is fixed (rows = samples): always orientation = 1L.
  geno_sub  <- geno_2d[keep_smp, keep_mk, drop = FALSE]
  split_res <- split_haplotypes(geno_sub, split_char = split_char,
                                orientation = 1L, ploidy = as.integer(ploidy),
                                sample_names = keep_smp, marker_names = keep_mk)
  # split_res$array is (n_samples, n_markers, ploidy)

  # ---- Resolve encoding and char_dict ----------------------------------------
  enc <- encoding %||% split_res$encode
  if (is.null(enc)) enc <- "unknown"
  message("build_geno_meta: detected encoding '", enc, "'")
  # Canonical encodings: "selectiontools" (numeric allele codes) and
  # "nucleotide" (A/T/C/G); "numeric" is what detect_genotype_encoding() emits
  # for 0/1/2-looking data and is kept as an alias here. The legacy "number" /
  # "char" tokens were removed in 0.5.1 — they used to be accepted here and then
  # stored verbatim in metadata$encoding, producing a GenoMeta that
  # geno_meta_views() rejected.
  enc <- .encode_arg(enc, "build_geno_meta",
                     allowed = c("selectiontools", "nucleotide", "numeric"),
                     strict  = FALSE)
  char_dict <- switch(enc,
    selectiontools = ,
    numeric        = as.list(DEFAULT_ENCODING_DICT),
    nucleotide     = setNames(as.list(seq_along(EXTENDED_ALLELE_DICT)),
                               EXTENDED_ALLELE_DICT),
    {
      # Derive from observed alleles
      obs <- as.vector(split_res$array)
      obs <- unique(obs[!is.na(obs)])
      setNames(as.list(seq_along(obs)), obs)
    }
  )

  # ---- Subset and reorder map to match array marker order --------------------
  map_sub <- map[map$Marker %in% keep_mk, , drop = FALSE]
  map_sub <- map_sub[match(split_res$markers, map_sub$Marker), , drop = FALSE]
  rownames(map_sub) <- NULL

  # ---- Build base GenoMeta ---------------------------------------------------
  meta <- new_geno_meta(
    geno_arrays = split_res$array,
    map_info    = map_sub,
    metadata    = list(
      ploidy        = as.integer(ploidy),
      char_dict     = char_dict,
      n_samples     = length(split_res$samples),
      n_markers     = length(split_res$markers),
      source_format = source_format,
      encoding      = enc
    ),
    sample_info = fam_si   # NULL → dummy FAM auto-built in new_geno_meta()
  )

  # ---- Attach blocks if provided ---------------------------------------------
  if (!is.null(blocks)) {
    meta <- attach_block_info(meta, blocks)
  }

  meta
}


#' Load block_info and block_marker_map_info from a block file
#'
#' Accepts either a \code{(Chr, Block, Start, End)} file (block coordinates)
#' or a \code{(Chr, Marker, Block)} file (marker-level). Detects format by
#' column count and content. Returns a list with \code{block_info} and
#' \code{block_marker_map_info} (either may be \code{NULL} if not derivable).
#'
#' @param block_file Path to the block file.
#' @param map_df     \code{map_info} data.frame used to expand coordinate
#'   blocks into marker-level entries via \code{\link{expand_block_info}}.
#' @return Named list: \code{block_info}, \code{block_marker_map_info}.
#' @keywords internal
.read_block_file <- function(block_file, map_df) {
  raw <- utils::read.table(block_file, header = TRUE, stringsAsFactors = FALSE,
                           sep = "", fill = TRUE)
  if (ncol(raw) < 3L)
    raw <- utils::read.csv(block_file, stringsAsFactors = FALSE)

  cn <- tolower(names(raw))

  # Detect format: coordinate file (4 cols with Start/End) vs marker-level (3 cols)
  has_start <- any(cn %in% c("start","start_pos","begin"))
  has_end   <- any(cn %in% c("end","end_pos"))

  if (has_start && has_end && ncol(raw) >= 4L) {
    # Coordinate format — expand to marker-level
    # Try case-insensitive name matching first, then positional fallback
    cn_lower  <- tolower(names(raw))
    idx_chr   <- which(cn_lower %in% c("chrom","chr","chromosome"))[1L]
    idx_blk   <- which(cn_lower %in% c("block","block_id","hap_block","hapblock"))[1L]
    idx_start <- which(cn_lower %in% c("start","start_pos","begin"))[1L]
    idx_end   <- which(cn_lower %in% c("end","end_pos"))[1L]
    if (!is.na(idx_chr) && !is.na(idx_blk) && !is.na(idx_start) && !is.na(idx_end)) {
      names(raw)[c(idx_chr, idx_blk, idx_start, idx_end)] <- c("Chrom","Block","Start","End")
    } else {
      names(raw)[1:4] <- c("Chrom","Block","Start","End")
    }
    # Carry a signal-marker selection through when the file has one (written by
    # a previous build_block_location_info(); see attach_signal_markers()).
    keep_cols  <- c("Chrom", "Block", "Start", "End",
                    intersect(.signal_marker_cols(), names(raw)))
    block_info <- raw[, keep_cols, drop = FALSE]
    block_marker_map_info <- expand_block_info(map_df, block_info)
  } else {
    # Assume marker-level format: (Chrom, Marker, Block)
    # Try case-insensitive name matching first, then positional fallback
    cn_lower   <- tolower(names(raw))
    idx_chrom  <- which(cn_lower %in% c("chrom","chr","chromosome"))[1L]
    idx_marker <- which(cn_lower %in% c("marker","snp","snpid","snp_id","mkr"))[1L]
    idx_block  <- which(cn_lower %in% c("block","block_id","hap_block","hapblock"))[1L]

    if (!is.na(idx_chrom) && !is.na(idx_marker) && !is.na(idx_block)) {
      names(raw)[c(idx_chrom, idx_marker, idx_block)] <- c("Chrom","Marker","Block")
    } else if (ncol(raw) >= 3L) {
      names(raw)[1:3] <- c("Chrom","Marker","Block")
    } else {
      stop(".read_block_file: could not identify Chrom/Marker/Block columns in '",
           block_file, "'. File has ", ncol(raw), " column(s) after parsing. ",
           "Expected a marker-level block file with columns: Chrom, Marker, Block.")
    }
    block_marker_map_info <- raw[, c("Chrom","Marker","Block"), drop = FALSE]
    # Derive block_info from marker positions
    bi_list    <- .compute_block_location(block_marker_map_info, map_df)
    block_info <- as.data.frame(bi_list$block_info[, c("Chrom","Block","Start","End"),
                                                   with = FALSE])
  }

  # Coerce identifier columns to character
  block_info$Chrom <- as.character(block_info$Chrom)
  block_info$Block <- as.character(block_info$Block)
  block_marker_map_info$Chrom  <- as.character(block_marker_map_info$Chrom)
  block_marker_map_info$Marker <- as.character(block_marker_map_info$Marker)
  block_marker_map_info$Block  <- as.character(block_marker_map_info$Block)

  list(block_info = block_info, block_marker_map_info = block_marker_map_info)
}


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

#' Filter a GenoMeta object by samples, markers, blocks, or chromosomes
#'
#' Intersects \code{geno_arrays} dim 1 with \code{map_info$Marker} and
#' optionally restricts to user-supplied samples, markers, haplotype blocks, or
#' chromosome regions. Multiple marker-level criteria (\code{keep_markers},
#' \code{blocks}, \code{chromosomes}) are combined with \strong{AND}
#' (intersection) logic. \code{block_marker_map_info} is pruned to the retained
#' markers; \code{block_info} rows are always preserved (blocks may span beyond
#' markers). \code{sample_info} rows and \code{marker_effects} rows are subset to
#' the survivors. \code{metadata} is updated with the new counts.
#'
#' @param meta         A \code{GenoMeta} list (from \code{\link{new_geno_meta}}
#'   or any reader function).
#' @param keep_markers Optional character vector of marker names to retain.
#'   Further intersected with the existing marker set.
#' @param keep_samples Optional character vector of sample names to retain.
#' @param blocks       Optional character vector of block IDs to retain
#'   (resolved to markers via \code{block_marker_map_info}).
#' @param chromosomes  Optional list of chromosome region specs, each a list
#'   with a required \code{chr} field and optional \code{start}/\code{end}
#'   fields, e.g. \code{list(list(chr = "1", start = 1e6, end = 5e6))}.
#'
#' @return A new \code{GenoMeta} list (input is not modified).
#'
#' @export
filter_geno_meta <- function(meta, keep_markers = NULL, keep_samples = NULL,
                             blocks = NULL, chromosomes = NULL) {
  if (!is_geno_meta(meta))
    stop("filter_geno_meta: 'meta' is not a valid GenoMeta.")

  current_samples <- dimnames(meta$geno_arrays)[[1L]]
  current_markers <- dimnames(meta$geno_arrays)[[2L]]

  # Intersect genotype dim2 with map_info
  valid_markers <- intersect(current_markers, meta$map_info$Marker)

  # --- Marker-level criteria (AND logic across keep_markers/blocks/chromosomes)
  marker_sets <- list()

  if (!is.null(keep_markers)) {
    not_found <- setdiff(keep_markers, valid_markers)
    if (length(not_found) > 0L)
      warning("filter_geno_meta: ", length(not_found),
              " requested marker(s) not present: ",
              paste(head(not_found, 5L), collapse = ", "),
              if (length(not_found) > 5L) " ..." else "")
    marker_sets$markers <- intersect(valid_markers, keep_markers)
  }

  if (!is.null(blocks)) {
    bmm <- meta$block_marker_map_info
    avail_blocks <- if (!is.null(bmm)) unique(bmm$Block) else character(0L)
    miss <- setdiff(blocks, avail_blocks)
    if (length(miss) > 0L)
      warning("filter_geno_meta: ", length(miss),
              " requested block(s) not present: ",
              paste(head(miss, 5L), collapse = ", "),
              if (length(miss) > 5L) " ..." else "")
    blk_markers <- if (!is.null(bmm)) bmm$Marker[bmm$Block %in% blocks]
                   else character(0L)
    marker_sets$blocks <- intersect(valid_markers, blk_markers)
  }

  if (!is.null(chromosomes)) {
    map <- meta$map_info
    chr_markers <- character(0L)
    for (region in chromosomes) {
      if (is.null(region$chr)) {
        warning("filter_geno_meta: region missing 'chr' field, skipping.")
        next
      }
      chr_val <- as.character(region$chr)
      hit     <- map$Chrom == chr_val
      if (!is.null(region$start)) hit <- hit & as.numeric(map$Pos) >= region$start
      if (!is.null(region$end))   hit <- hit & as.numeric(map$Pos) <= region$end
      chr_markers <- union(chr_markers, map$Marker[hit])
    }
    marker_sets$chromosomes <- intersect(valid_markers, chr_markers)
  }

  if (length(marker_sets) > 0L)
    valid_markers <- Reduce(intersect, marker_sets)

  valid_samples <- current_samples
  if (!is.null(keep_samples)) {
    not_found_s <- setdiff(keep_samples, current_samples)
    if (length(not_found_s) > 0L)
      warning("filter_geno_meta: ", length(not_found_s),
              " requested sample(s) not present: ",
              paste(head(not_found_s, 5L), collapse = ", "),
              if (length(not_found_s) > 5L) " ..." else "")
    valid_samples <- intersect(current_samples, keep_samples)
  }

  # Subset array
  new_arr <- meta$geno_arrays[valid_samples, valid_markers, , drop = FALSE]

  # Subset map_info
  new_map <- meta$map_info[meta$map_info$Marker %in% valid_markers, , drop = FALSE]

  # Prune block_marker_map_info — retain only rows whose markers survived
  new_bmm <- if (!is.null(meta$block_marker_map_info)) {
    meta$block_marker_map_info[
      meta$block_marker_map_info$Marker %in% valid_markers, , drop = FALSE
    ]
  } else NULL

  # block_info: rows always preserved (blocks may extend beyond remaining
  # markers), but a stored signal-marker selection must not point at a marker
  # that just left the panel — .prune_signal_markers() drops those and
  # recomputes N_Signal (no-op when no selection is stored).
  new_bi <- .prune_signal_markers(meta$block_info, valid_markers,
                                  caller = "filter_geno_meta")

  # Subset sample_info to surviving samples
  new_si <- meta$sample_info[
    meta$sample_info$Sample %in% valid_samples, , drop = FALSE
  ]

  # Prune marker_effects to surviving markers (when present). Handles both the
  # canonical named list of per-trait tables and a bare data.frame.
  new_me <- if (is.null(meta$marker_effects)) {
    NULL
  } else if (is.data.frame(meta$marker_effects)) {
    meta$marker_effects[meta$marker_effects$Marker %in% valid_markers, ,
                        drop = FALSE]
  } else {
    lapply(meta$marker_effects, function(df)
      df[df$Marker %in% valid_markers, , drop = FALSE])
  }

  # Prune phenotype to surviving samples (when present)
  new_pheno <- if (!is.null(meta$phenotype)) {
    meta$phenotype[meta$phenotype$Sample %in% valid_samples, , drop = FALSE]
  } else NULL

  # Update metadata
  new_meta              <- meta$metadata
  new_meta$n_markers    <- length(valid_markers)
  new_meta$n_samples    <- length(valid_samples)

  new_geno_meta(
    geno_arrays           = new_arr,
    map_info              = new_map,
    block_info            = new_bi,
    block_marker_map_info = new_bmm,
    metadata              = new_meta,
    sample_info           = new_si,
    marker_effects        = new_me,
    phenotype             = new_pheno
  )
}


# ---------------------------------------------------------------------------
# Internal: 3-D array → 2-D conversions
# ---------------------------------------------------------------------------

# .meta_to_allele_matrix(geno_meta, sep = "/")
#
# Reconstruct a 2-D character genotype matrix from the 3-D split array.
# Each cell is the allele calls from all ploidy layers joined by `sep`
# (e.g. "A/T" for diploid). NA in any ploidy layer → NA in the output cell.
# No subsetting is performed; call filter_geno_meta() first if needed.
#
# Returns: character matrix [n_samples × n_markers].
.meta_to_allele_matrix <- function(geno_meta, sep = "/") {
  arr          <- geno_meta$geno_arrays
  n_samples    <- dim(arr)[1L]
  n_markers    <- dim(arr)[2L]
  n_ploidy     <- dim(arr)[3L]
  sample_names <- dimnames(arr)[[1L]]
  marker_names <- dimnames(arr)[[2L]]

  result  <- arr[, , 1L, drop = TRUE]          # [n_samples × n_markers] character
  na_mask <- is.na(result)

  for (k in seq(2L, n_ploidy)) {
    layer_k  <- arr[, , k, drop = TRUE]
    na_mask  <- na_mask | is.na(layer_k)
    result   <- matrix(paste(result, layer_k, sep = sep),
                       nrow = n_samples, ncol = n_markers)
  }

  result[na_mask] <- NA_character_
  dimnames(result) <- list(sample_names, marker_names)
  result   # [n_samples × n_markers]
}


# .meta_to_dosage_matrix(geno_meta, target_alleles = NULL)
#
# Count how many ploidy layers carry a specified allele per marker × sample.
#
# target_alleles: character vector, length == n_markers (positionally matched
#   to dim 1 of geno_arrays). No names required.
#   NULL → auto minor allele per marker (least frequent across all layers).
#
# Strict: if length(target_alleles) != n_markers an error is thrown.
# Subsetting must be done by the caller via filter_geno_meta() first.
# NA in any ploidy layer → NA for that marker × sample cell.
#
# Returns: integer matrix [n_samples × n_markers].
.meta_to_dosage_matrix <- function(geno_meta, target_alleles = NULL) {
  arr          <- geno_meta$geno_arrays
  n_samples    <- dim(arr)[1L]
  n_markers    <- dim(arr)[2L]
  n_ploidy     <- dim(arr)[3L]
  sample_names <- dimnames(arr)[[1L]]
  marker_names <- dimnames(arr)[[2L]]

  if (!is.null(target_alleles)) {
    if (length(target_alleles) != n_markers)
      stop(".meta_to_dosage_matrix: length(target_alleles) (", length(target_alleles),
           ") must equal n_markers (", n_markers, "). ",
           "Call filter_geno_meta() to subset the GenoMeta to the desired markers first.")
    target_alleles <- as.character(target_alleles)
  } else {
    # Auto: minor allele = least frequent OBSERVED allele per marker (ties broken
    # alphabetically, matching names(sort(table(obs)))[1L]). Detection is in R;
    # counting is the marker_allele_counts kernel.
    counts <- .marker_allele_counts(arr)
    levs   <- colnames(counts)
    target_alleles <- vapply(seq_len(n_markers), function(m) {
      cj <- counts[m, ]
      if (all(cj == 0L)) return(NA_character_)
      cj[cj == 0L] <- .Machine$integer.max         # ignore alleles unobserved here
      levs[which.min(cj)]                          # least frequent; ties -> alpha
    }, character(1L))
  }

  # Dosage = count of target allele per cell; NA if ANY ploidy layer is missing
  # (na_if_any = TRUE preserves the original has_na semantics).
  .dosage_matrix(arr, target_alleles, na_if_any = TRUE)
}


# ---------------------------------------------------------------------------
# Block attachment
# ---------------------------------------------------------------------------

#' Attach or update block information in a GenoMeta object
#'
#' Updates an existing \code{GenoMeta} with block assignment data and returns
#' a new \code{GenoMeta} with both \code{block_info} and
#' \code{block_marker_map_info} populated. The original object is not modified.
#' Any previously stored block data is replaced.
#'
#' @section Source types:
#' \describe{
#'   \item{\code{"block_marker_map"}}{A data.frame with columns
#'     \code{Chrom}, \code{Marker}, \code{Block} — one row per marker per
#'     block. All markers must be present in \code{geno_meta$map_info}.}
#'   \item{\code{"block_info"}}{A data.frame with columns
#'     \code{Chrom}, \code{Block}, \code{Start}, \code{End} — block coordinate
#'     table. \code{\link{expand_block_info}} is called internally to derive
#'     the marker-level assignment.}
#'   \item{\code{"file"}}{Path to a block file. Format is detected
#'     automatically: four-column files are treated as \code{block_info};
#'     three-column files as \code{block_marker_map}.}
#'   \item{\code{"auto"}}{Detect from column names of \code{block_source}:
#'     presence of \code{"Start"} or \code{"End"} selects \code{block_info};
#'     otherwise \code{block_marker_map}.}
#' }
#'
#' @param geno_meta    A valid \code{GenoMeta} list.
#' @param block_source A data.frame supplying block assignments, or a file
#'   path. See \emph{Source types} above.
#' @param source_type  Character. One of \code{"auto"} (default),
#'   \code{"block_marker_map"}, \code{"block_info"}, or \code{"file"}.
#'
#' @return A new \code{GenoMeta} with \code{block_info} and
#'   \code{block_marker_map_info} populated.
#'
#' @examples
#' \dontrun{
#' # Typical workflow: read data, then attach blocks from a marker-level table
#' meta  <- read_st_genotype("geno.csv")
#' bmap  <- read.csv("block_marker_map.csv")           # Chrom, Marker, Block
#' meta  <- attach_block_info(meta, bmap)
#'
#' # From block coordinates (Start / End); markers are matched by position
#' binfo <- data.frame(Chrom = "1", Block = "B1", Start = 1e5, End = 5e5)
#' meta  <- attach_block_info(meta, binfo, source_type = "block_info")
#'
#' # After LD-based block definition
#' ld_out <- define_haplotype_blocks("geno.st", "map.txt", output_dir = "blocks/")
#' meta   <- attach_block_info(meta, ld_out$block_marker_map)
#' }
#'
#' @seealso \code{\link{expand_block_info}}, \code{\link{define_haplotype_blocks}},
#'   \code{\link{filter_geno_meta}}, \code{\link{new_geno_meta}}
#'
#' @export
attach_block_info <- function(geno_meta,
                               block_source,
                               source_type = c("auto", "block_marker_map",
                                               "block_info", "file")) {
  if (!is_geno_meta(geno_meta))
    stop("attach_block_info: 'geno_meta' is not a valid GenoMeta.")
  source_type <- match.arg(source_type)

  # ── File path ──────────────────────────────────────────────────────────────
  if (source_type == "file" ||
      (source_type == "auto" && is.character(block_source) &&
       length(block_source) == 1L)) {
    if (!file.exists(block_source))
      stop("attach_block_info: file not found: ", block_source)
    block_list <- .read_block_file(block_source, geno_meta$map_info)
    return(.replace_block_in_meta(geno_meta,
                                  block_list$block_marker_map_info,
                                  block_list$block_info))
  }

  # ── Data.frame source ──────────────────────────────────────────────────────
  if (!is.data.frame(block_source) && !inherits(block_source, "data.table"))
    stop("attach_block_info: 'block_source' must be a data.frame or a file path.")

  resolved_type <- if (source_type == "auto") {
    if (any(c("Start", "End") %in% names(block_source))) "block_info"
    else "block_marker_map"
  } else source_type

  if (resolved_type == "block_marker_map") {
    bmm <- as.data.frame(block_source, stringsAsFactors = FALSE)
    if (!all(c("Chrom", "Marker", "Block") %in% names(bmm)))
      names(bmm)[seq_len(min(3L, ncol(bmm)))] <-
        c("Chrom", "Marker", "Block")[seq_len(min(3L, ncol(bmm)))]

    bad <- setdiff(bmm$Marker, geno_meta$map_info$Marker)
    if (length(bad) > 0L)
      warning("attach_block_info: ", length(bad),
              " marker(s) in block_source absent from GenoMeta map_info: ",
              paste(head(bad, 5L), collapse = ", "),
              if (length(bad) > 5L) " ..." else "")

    loc <- .compute_block_location(bmm, geno_meta$map_info)
    return(.replace_block_in_meta(geno_meta, bmm, loc$block_info))
  }

  # resolved_type == "block_info" — expand coordinates to marker-level first
  bi_df <- as.data.frame(block_source, stringsAsFactors = FALSE)
  bmm   <- expand_block_info(geno_meta$map_info, bi_df)
  loc   <- .compute_block_location(bmm, geno_meta$map_info)
  .replace_block_in_meta(geno_meta, bmm, loc$block_info)
}


# Internal: replace block fields and return updated GenoMeta copy
# Rebuilds via new_geno_meta() so the same block filtering, singleton dummy-block
# creation, and block ordering are applied as in any other GenoMeta construction.
#
# Note on signal markers: the incoming block_info is recomputed from marker
# positions, so any Signal_Marker/Signal_Score/N_Signal columns previously
# written by attach_signal_markers() are lost here. That is correct —
# redefining the blocks invalidates a selection made against the old ones — but
# it must not be silent, hence the message.
.replace_block_in_meta <- function(geno_meta, block_marker_map, block_info) {
  old_sig <- intersect(.signal_marker_cols(), names(geno_meta$block_info %||%
                                                      data.frame()))
  new_sig <- intersect(.signal_marker_cols(), names(block_info))
  if (length(old_sig) > 0L && length(new_sig) == 0L)
    message("attach_block_info: block definitions replaced — the stored ",
            "signal-marker selection was dropped. Re-run ",
            "score_marker_signal()/select_signal_markers() against the new ",
            "blocks.")

  new_geno_meta(
    geno_arrays           = geno_meta$geno_arrays,
    map_info              = geno_meta$map_info,
    block_info            = as.data.frame(block_info,
                                          stringsAsFactors = FALSE),
    block_marker_map_info = as.data.frame(block_marker_map,
                                          stringsAsFactors = FALSE),
    metadata              = geno_meta$metadata
  )
}


# ---------------------------------------------------------------------------
# Export (3D → 2D)
# ---------------------------------------------------------------------------

#' Convert a GenoMeta 3-D array to a 2-D genotype matrix
#'
#' Collapses the haplotype dimension of \code{geno_arrays} into genotype calls
#' laid out per a format \code{preset} — the \emph{same} preset/parameter-set
#' mechanism used by \code{\link{read_genotype}}, so reading and writing share
#' one description of each format. The preset's \code{geno_orient} sets the
#' output orientation and \code{split_char} the in-cell allele separator; the
#' map is returned with its columns reordered to the preset's \code{map_order}.
#'
#' @param meta          A valid \code{GenoMeta} list.
#' @param preset        Format preset name, a custom parameter-set list, or
#'   \code{NULL}. See \code{\link{GENO_FORMAT_PRESETS}}.
#' @param output_format Deprecated alias for \code{preset} accepting the legacy
#'   values \code{"plink_ped"}, \code{"st"}, \code{"gs"}. Used only when
#'   \code{preset} is \code{NULL}.
#' @param presets       Preset registry. Default \code{\link{GENO_FORMAT_PRESETS}}.
#' @param split_char    Optional override for the in-cell allele separator.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{matrix}}{2-D character matrix in the preset's orientation.}
#'   \item{\code{map}}{data.frame with columns ordered per the preset's
#'     \code{map_order}.}
#'   \item{\code{fam}}{The \code{sample_info} FAM table.}
#'   \item{\code{output_format}}{The resolved preset name/spec source.}
#'   \item{\code{spec}}{The resolved format parameter set.}
#' }
#'
#' @seealso \code{\link{resolve_geno_format}}, \code{\link{write_geno_meta}}
#'
#' @export
geno_meta_to_matrix <- function(meta,
                                 preset        = NULL,
                                 output_format = NULL,
                                 presets       = GENO_FORMAT_PRESETS,
                                 split_char    = NULL) {
  if (!is_geno_meta(meta))
    stop("geno_meta_to_matrix: 'meta' is not a valid GenoMeta.")

  # Back-compat: fall back to the legacy output_format argument (default "PL"
  # mirrors the former "plink_ped" default).
  if (is.null(preset)) preset <- output_format %||% "PL"
  spec <- resolve_geno_format(preset, presets = presets,
                              split_char = split_char)

  # Join ploidy layers with the preset separator → [n_samples × n_markers].
  am <- .meta_to_allele_matrix(meta, sep = spec$split_char)

  # Orient: "sample" keeps samples-as-rows (the native layout); "marker"
  # transposes to markers-as-rows.
  mat <- if (spec$geno_orient == "marker") t(am) else am

  # Map columns reordered to the preset's file order.
  base_map <- meta$map_info[, c("Chrom", "Marker", "Pos"), drop = FALSE]
  map_out  <- base_map[, spec$map_order, drop = FALSE]

  list(matrix        = mat,
       map           = map_out,
       fam           = meta$sample_info,
       output_format = if (is.character(preset)) preset else "custom",
       spec          = spec)
}


#' Write a GenoMeta object to disk in a specified format
#'
#' Resolves a format \code{preset} (the same parameter-set mechanism as
#' \code{\link{read_genotype}}) and writes the genotype to \code{path}. When the
#' preset embeds FAM columns (\code{fam_cols > 0}, e.g. \code{"PL"}), a PLINK
#' \code{.ped}/\code{.map} pair is written: the FAM block is taken from
#' \code{meta$sample_info} and alleles are space-separated. Otherwise a single
#' delimited file is written in the preset's orientation, with the row
#' identifier (marker or sample IDs) in the first column.
#'
#' @param meta          A valid \code{GenoMeta} list.
#' @param path          Output file path (without extension for PLINK output;
#'   full path including extension otherwise).
#' @param preset        Format preset name, a custom parameter-set list, or
#'   \code{NULL}. See \code{\link{GENO_FORMAT_PRESETS}}.
#' @param output_format Deprecated alias for \code{preset} accepting the legacy
#'   values \code{"plink_ped"}, \code{"st"}, \code{"gs"}. Used only when
#'   \code{preset} is \code{NULL}.
#' @param presets       Preset registry. Default \code{\link{GENO_FORMAT_PRESETS}}.
#' @param split_char    Optional override for the in-cell allele separator.
#' @param pos_unit      Unit for written map positions: \code{"bp"} (default) or
#'   \code{"mbp"}. Positions are stored internally in bp and converted on write.
#'
#' @return Invisibly, the result of \code{\link{geno_meta_to_matrix}}.
#'
#' @seealso \code{\link{resolve_geno_format}}, \code{\link{geno_meta_to_matrix}}
#'
#' @export
write_geno_meta <- function(meta, path,
                            preset        = NULL,
                            output_format = NULL,
                            presets       = GENO_FORMAT_PRESETS,
                            split_char    = NULL,
                            pos_unit      = "bp") {
  if (is.null(preset)) preset <- output_format %||% "PL"
  result <- geno_meta_to_matrix(meta, preset = preset, presets = presets,
                                split_char = split_char)
  spec   <- result$spec

  if (spec$fam_cols > 0L) {
    # ---- PLINK-style output: .ped (FAM block + space-separated alleles) ------
    arr     <- meta$geno_arrays
    sn      <- dimnames(arr)[[1L]]
    mn      <- dimnames(arr)[[2L]]
    ploidy  <- dim(arr)[3L]

    allele_block <- matrix(NA_character_, nrow = length(sn), ncol = length(mn),
                           dimnames = list(sn, mn))
    for (j in seq_along(mn)) {
      layers <- vapply(seq_len(ploidy),
                       function(h) as.character(arr[, j, h]),
                       character(length(sn)))
      if (is.null(dim(layers))) layers <- matrix(layers, nrow = length(sn))
      allele_block[, j] <- apply(layers, 1L, function(r) paste(r, collapse = " "))
    }

    si <- meta$sample_info
    ped_df <- cbind(
      si$Family, si$Sample, si$Parent1, si$Parent2, si$Sex, -9L,
      allele_block
    )
    ped_path <- paste0(path, ".ped")
    utils::write.table(ped_df, ped_path, sep = " ",
                       row.names = FALSE, col.names = FALSE, quote = FALSE)

    # Companion .map: chr, marker, cM (0), position (in requested unit).
    map_out <- data.frame(
      Chrom  = meta$map_info$Chrom,
      Marker = meta$map_info$Marker,
      cM     = 0,
      Pos    = .convert_pos_unit(suppressWarnings(as.numeric(meta$map_info$Pos)),
                                 from_unit = "bp", to_unit = pos_unit),
      stringsAsFactors = FALSE
    )
    map_path <- paste0(path, ".map")
    utils::write.table(map_out, map_path, sep = "\t",
                       row.names = FALSE, col.names = FALSE, quote = FALSE)
    message("Written: ", ped_path, " and ", map_path)

  } else {
    # ---- Single delimited file in the preset's orientation -------------------
    mat      <- result$matrix
    id_label <- if (spec$geno_orient == "marker") "Marker" else "Sample"
    out_df   <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
    out_df   <- cbind(setNames(
                        data.frame(rownames(mat), stringsAsFactors = FALSE),
                        id_label),
                      out_df)
    utils::write.csv(out_df, path, row.names = FALSE, quote = FALSE)
    message("Written: ", path)
  }

  invisible(result)
}
