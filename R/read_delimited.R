# =============================================================================
# read_delimited.R — delimited genotype tables through a preset
#
# read_genotype() reads any delimited genotype table whose layout a preset
# describes (GENO_FORMAT_PRESETS, format_presets.R); read_st_genotype() and
# read_gs_genotype() are its SelectionTools and genomicSimulation wrappers
# (the SelectionTools one behind use_st_formats()). .read_genotype_files()
# reads the geno / map / FAM files unmodified. Split out of format_helper.R
# in 0.9.0 (D2-9).
# =============================================================================

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
#' @noRd
.read_genotype_files <- function(geno_file, map_file = NULL, spec) {
  geno_file <- as.character(geno_file)
  if (!file.exists(geno_file))
    .hs_abort(".read_genotype_files: genotype file not found: ", geno_file)

  # Identifier columns (the row ID, or the whole FAM block) are read as text so
  # a sample or marker ID such as "007" or "1e5" keeps its spelling, as
  # read_plink_ped() does; the value columns keep fread's type, as before.
  raw <- data.table::fread(geno_file, header = spec$geno_header,
                           colClasses = list(character = seq_len(max(1L, spec$fam_cols))),
                           data.table = FALSE)
  if (ncol(raw) < 2L)
    .hs_abort(".read_genotype_files: genotype file must have an identifier column ",
         "plus at least one value column.")

  fam <- NULL
  if (spec$fam_cols > 0L) {
    # PLINK-style embedded FAM block: columns 1:fam_cols are pedigree fields;
    # genotype values begin at column fam_cols + 1. Sample IDs are the IID
    # column (position 2 of the FAM block, or position 1 when fam_cols == 1).
    if (ncol(raw) <= spec$fam_cols)
      .hs_abort(".read_genotype_files: file has ", ncol(raw), " column(s) but the ",
           "preset reserves the first ", spec$fam_cols, " as FAM fields \u2014 no ",
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


#' Read a delimited genotype table into a GenoMeta object
#'
#' Unified reader for delimited genotype/map file pairs. The on-disk file layout
#' is described by a \code{preset} parameter set (see
#' \code{\link{GENO_FORMAT_PRESETS}}); any field can be overridden, or a fully
#' custom parameter set can be supplied. File reading is delegated to
#' \code{.read_genotype_files()} (which returns the genotype, map, and fam
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
#' @param block_file Optional block map: a path or a data.frame in one of the
#'   two shapes \code{\link{attach_block_info}} accepts (block-marker map
#'   \code{Chrom, Marker, Block} or block-region map \code{Chrom, Block, Start,
#'   End}; header required). Attached after the read.
#' @param block_overlap Overlap policy for the block map, as
#'   \code{attach_block_info(overlap =)}: \code{"error"} (default) or
#'   \code{"first"}.
#' @param presets    Preset registry to resolve \code{preset} against. Default
#'   \code{\link{GENO_FORMAT_PRESETS}}.
#' @param split_char Optional override for the allele delimiter (\code{"/"} or
#'   \code{""}).
#' @param encoding   Optional override for the allele encoding
#'   (\code{"nucleotide"}, or \code{"selectiontools"} for a file in
#'   SelectionTools' 1/2/3/4 codes, which are translated on the read; see
#'   \code{\link{build_geno_meta}}).
#' @param ploidy     Integer ploidy level (default \code{2}).
#' @param pos_unit   Position unit in \code{map_file}: \code{"bp"} (default) or
#'   \code{"mbp"}.
#' @param extra_na   Optional character vector of additional missing sentinels
#'   this file uses, passed through to the normalisation portal (see
#'   \code{\link{build_geno_meta}}).
#' @param missing_chars Extra missing-allele characters for this read:
#'   \code{NULL} (default) takes the session registry
#'   (\code{\link{use_missing_chars}}); a character vector names them directly
#'   (\code{character(0)} ignores the registry). See
#'   \code{\link{new_geno_meta}}.
#' @param st_formats For the \code{"ST"} preset (or any parameter set whose
#'   \code{encoding} is \code{"selectiontools"}): \code{TRUE} enables
#'   SelectionTools' allele codes for this call; \code{NULL} (default) takes
#'   the session switch (\code{\link{use_st_formats}}). The codes are
#'   translated to nucleotides on the read.
#'
#' @section SelectionTools' allele codes:
#' A parameter set whose \code{encoding} is \code{"selectiontools"} (the
#' \code{"ST"} preset) describes cells in SelectionTools' codes — \code{A = 1,
#' C = 2, G = 3, T = 4}, missing \code{-1}, written \code{1/-1}. Since 0.9.0 a
#' \code{GenoMeta} stores nucleotide alleles only, so such a read is behind
#' \code{\link{use_st_formats}} (SelectionTools is being retired) and the
#' cells are translated before the container is built; the container is
#' stamped \code{encoding = "nucleotide"} with \code{source_format =
#' "selectiontools"} as its provenance. A file whose digits are allele labels
#' of your own (PLINK's 1/2 coding) is read with \code{encoding =
#' "nucleotide"}, which keeps the digits as labels.
#'
#' @return A valid \code{GenoMeta} list.
#'
#' @seealso \code{\link{resolve_geno_format}}, \code{\link{build_geno_meta}},
#'   \code{\link{read_st_genotype}}, \code{\link{read_gs_genotype}},
#'   \code{\link{read_plink_ped}}, \code{\link{st_convert}}
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
                          pos_unit   = "bp",
                          extra_na   = NULL,
                          missing_chars = NULL,
                          block_overlap = c("error", "first"),
                          st_formats = NULL) {

  spec  <- resolve_geno_format(preset, presets = presets,
                               split_char = split_char, encoding = encoding)
  # SelectionTools' allele codes (0.9.0, C-5): gated, translated on the read.
  st_codes <- identical(spec$encoding, "selectiontools")
  if (st_codes) .st_gate("read_genotype", st_formats)
  parts <- .read_genotype_files(geno_file, map_file, spec)
  if (st_codes) {
    parts$geno    <- .st_decode_cells(parts$geno, spec$split_char, missing_chars)
    spec$encoding <- "nucleotide"
  }

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
    block_overlap = block_overlap,
    ploidy        = as.integer(ploidy),
    split_char    = spec$split_char,
    encoding      = spec$encoding,
    source_format = spec$source_format,
    pos_unit      = pos_unit,
    extra_na      = extra_na,
    missing_chars = missing_chars
  )
}


#' Read a SelectionTools genotype file into a GenoMeta object
#'
#' Thin wrapper over \code{\link{read_genotype}} with \code{preset = "ST"}
#' (markers as rows, slash-separated numeric cells, \code{name | chr | pos}
#' map). The cells are SelectionTools' allele codes, so the read is behind
#' \code{\link{use_st_formats}} and the codes are translated to nucleotides
#' (see \code{\link{read_genotype}}, section \emph{SelectionTools' allele
#' codes}); \code{\link{st_convert}} reads SelectionTools' own
#' tab-separated matrix format.
#'
#' @param file       Path to the SelectionTools genotype file.
#' @param map_file   Optional path to a map file (\code{name | chr | pos}).
#' @param block_file Optional block map (a path or a data.frame; see
#'   \code{\link{attach_block_info}} for the two shapes), attached after the
#'   read.
#' @param block_overlap Overlap policy for the block map (\code{"error"},
#'   the default, or \code{"first"}); see \code{\link{attach_block_info}}.
#' @param pos_unit   Position unit: \code{"bp"} (default) or \code{"mbp"}.
#' @param extra_na   Optional character vector of additional missing sentinels
#'   this file uses (see \code{\link{build_geno_meta}}).
#' @param missing_chars Extra missing-allele characters for this read (see
#'   \code{\link{read_genotype}}).
#' @param st_formats \code{TRUE} enables SelectionTools' allele codes for this
#'   call; \code{NULL} (default) takes the session switch
#'   (\code{\link{use_st_formats}}).
#'
#' @return A valid \code{GenoMeta} list.
#'
#' @export
read_st_genotype <- function(file, map_file = NULL, block_file = NULL,
                              pos_unit = "bp", extra_na = NULL,
                              missing_chars = NULL,
                              block_overlap = c("error", "first"),
                              st_formats = NULL) {
  read_genotype(file, map_file = map_file, preset = "ST",
                block_file = block_file, pos_unit = pos_unit,
                extra_na = extra_na, missing_chars = missing_chars,
                block_overlap = block_overlap, st_formats = st_formats)
}


#' Read a genomicSimulation genotype file into a GenoMeta object
#'
#' Thin wrapper over \code{\link{read_genotype}} with \code{preset = "GS"}
#' (markers as rows, two-character allele cells, e.g. \code{"AA"}).
#'
#' @param file       Path to the genomicSimulation genotype file.
#' @param map_file   Optional path to a map file (\code{name | chr | pos}).
#' @param block_file Optional block map (a path or a data.frame; see
#'   \code{\link{attach_block_info}} for the two shapes), attached after the
#'   read.
#' @param block_overlap Overlap policy for the block map (\code{"error"},
#'   the default, or \code{"first"}); see \code{\link{attach_block_info}}.
#' @param pos_unit   Position unit: \code{"bp"} (default) or \code{"mbp"}.
#' @param extra_na   Optional character vector of additional missing sentinels
#'   this file uses (see \code{\link{build_geno_meta}}).
#' @param missing_chars Extra missing-allele characters for this read (see
#'   \code{\link{read_genotype}}).
#'
#' @return A valid \code{GenoMeta} list.
#'
#' @export
read_gs_genotype <- function(file, map_file = NULL, block_file = NULL,
                              pos_unit = "bp", extra_na = NULL,
                              missing_chars = NULL,
                              block_overlap = c("error", "first")) {
  read_genotype(file, map_file = map_file, preset = "GS",
                block_file = block_file, pos_unit = pos_unit,
                extra_na = extra_na, missing_chars = missing_chars,
                block_overlap = block_overlap)
}
