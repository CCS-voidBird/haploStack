# =============================================================================
# format_process.R
#
# process_genotype() — convert a raw genotype data frame into a standardised
# metadata list (3D array + map + sample names + encoding info).
# =============================================================================

#' Process a raw genotype data frame into standardised metadata
#'
#' Reads any common genotype format, splits diploid calls into a 3D array,
#' normalises the marker map, detects minor alleles, and returns a single
#' metadata list consumed by \code{\link{transform_genotype}} and the catalogue
#' pipeline functions.
#'
#' @param genotype    A matrix or data frame of raw genotype values. Cells may
#'   be strings like \code{"1/4"} (SelectionTools), \code{"A/T"} (nucleotide),
#'   or integers \code{0}/\code{1}/\code{2} (dosage).
#' @param map_info    Optional map data frame. Accepted formats are
#'   \itemize{
#'     \item \strong{SelectionTools} (\code{map_preset = "ST"}): columns
#'       \code{name | chr | pos}
#'     \item \strong{PLINK} (\code{map_preset = "PL"}): columns
#'       \code{chr | name | cM | bp} (4-col) or \code{chr | name | bp} (3-col)
#'     \item \strong{generic}: auto-detection via column-name matching
#'       (\code{marker/name/id}, \code{chr/chrom}, \code{pos/bp/position}).
#'   }
#' @param preset      Short-cut preset for common tools:
#'   \describe{
#'     \item{\code{"ST"}}{SelectionTools — row-marker, slash-numeric encoding}
#'     \item{\code{"GS"}}{genomicSimulation — row-marker, nucleotide chars}
#'     \item{\code{"MA"}}{MAGIC / generic — row-sample, slash-separated chars}
#'     \item{\code{"PL"}}{PLINK — row-sample, dosage (0/1/2)}
#'     \item{\code{NULL}}{manual specification via the remaining arguments}
#'   }
#' @param orientation Integer. \code{1} = rows are samples, cols are markers;
#'   \code{2} = rows are markers, cols are samples. Auto-detected when
#'   \code{NULL} (default).
#' @param split_char  Allele separator within each cell (e.g. \code{"/"}). The
#'   preset overrides this argument.
#' @param char_dict   Encoding dictionary (default: \code{DEFAULT_ENCODING_DICT}).
#' @param encode      Encoding of stored values: \code{"nucleotide"},
#'   \code{"selectiontools"}, or \code{"allele"}. Overridden by \code{preset}.
#' @param ploidy      Integer ploidy level (default \code{2}). Validated against
#'   the data.
#' @param map_preset  Format preset for \code{map_info}: \code{"ST"},
#'   \code{"PL"}, or \code{NULL} (auto-detect).
#' @param pos_unit    Position unit in \code{map_info}: \code{"bp"} (default)
#'   or \code{"mbp"}. Internally all positions are stored in bp.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{genotype}}{3-D character array, \code{(n_samples, n_markers, ploidy)}.}
#'   \item{\code{map}}{data.frame with columns \code{Chr}, \code{Marker},
#'     \code{Pos} (bp), \code{Minor}.}
#'   \item{\code{sample}}{Character vector of sample names.}
#'   \item{\code{split_char}}{Resolved allele separator.}
#'   \item{\code{orientation}}{Resolved orientation (1 or 2).}
#'   \item{\code{char_dict}}{Encoding dictionary used.}
#'   \item{\code{encode}}{Resolved encoding.}
#' }
#'
#' @seealso \code{\link{transform_genotype}}, \code{\link{detect_genotype_encoding}}
#'
#' @examples
#' # Minimal example — SelectionTools format
#' geno <- matrix(c("1/1","1/4","4/4","1/1"), nrow = 2,
#'                dimnames = list(c("Mkr1","Mkr2"), c("S1","S2")))
#' meta <- process_genotype(geno, preset = "ST")
#' str(meta)
#'
#' @export
process_genotype <- function(genotype,
                              map_info    = NULL,
                              preset      = NULL,
                              orientation = NULL,
                              split_char  = NULL,
                              char_dict   = DEFAULT_ENCODING_DICT,
                              encode      = NULL,
                              ploidy      = 2L,
                              map_preset  = NULL,
                              pos_unit    = "bp") {

  params      <- .resolve_preset(preset, orientation, split_char, encode)
  orientation <- params$orientation
  split_char  <- params$split_char
  encode      <- params$encode

  if (is.na(orientation) || is.null(orientation)) {
    map_tmp     <- .process_map(map_info, NULL, map_preset, pos_unit)
    orientation <- .detect_orientation(genotype, map_tmp)
  }

  split_result <- split_haplotypes(
    geno_data   = genotype,
    split_char  = split_char,
    orientation = orientation,
    ploidy      = as.integer(ploidy)
  )

  geno_array   <- split_result$array
  marker_names <- split_result$markers
  sample_names <- split_result$samples
  if (!is.null(split_result$encode)) encode <- split_result$encode

  map_df <- .process_map(map_info, marker_names, map_preset, pos_unit)
  if (!is.null(map_df)) .validate_geno_map(geno_array, map_df)

  minor_alleles <- .detect_minor_alleles(geno_array, encode, char_dict)

  if (!is.null(map_df)) {
    map_df$Minor <- minor_alleles[match(map_df$Marker, names(minor_alleles))]
  } else {
    map_df <- data.frame(
      Chr    = NA_character_,
      Marker = marker_names,
      Pos    = NA_real_,
      Minor  = minor_alleles[marker_names],
      stringsAsFactors = FALSE
    )
  }

  list(
    genotype    = geno_array,
    map         = map_df,
    sample      = sample_names,
    split_char  = split_char,
    orientation = orientation,
    char_dict   = char_dict,
    encode      = encode
  )
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Standardise a map data frame to columns Chr, Marker, Pos
#'
#' All positions are stored internally in bp. Use \code{.convert_pos_unit()}
#' when writing output in mbp.
#'
#' @param map_info    Raw map data frame, or \code{NULL}.
#' @param marker_names Character vector (for ordering / cross-checking).
#' @param map_preset  \code{"ST"}, \code{"PL"}, or \code{NULL}.
#' @param pos_unit    Input position unit: \code{"bp"} or \code{"mbp"}.
#' @return Standardised data.frame with columns Chr, Marker, Pos, or NULL.
#' @keywords internal
.process_map <- function(map_info,
                          marker_names = NULL,
                          map_preset   = NULL,
                          pos_unit     = "bp") {
  if (is.null(map_info)) return(NULL)
  df <- as.data.frame(map_info, stringsAsFactors = FALSE)

  if (!is.null(map_preset) && map_preset == "ST") {
    if (ncol(df) >= 3L) colnames(df)[seq_len(3L)] <- c("Marker", "Chrom", "Pos")

  } else if (!is.null(map_preset) && map_preset == "PL") {
    if (ncol(df) >= 4L) {
      df <- df[, c(2L, 1L, 4L)]
    } else if (ncol(df) >= 3L) {
      df <- df[, c(2L, 1L, 3L)]
    }
    colnames(df)[seq_len(3L)] <- c("Marker", "Chrom", "Pos")

  } else {
    cn         <- tolower(colnames(df))
    col_marker <- which(cn %in% c("marker","name","snp","id","rs","locus"))[1L]
    col_chr    <- which(cn %in% c("chr","chrom","chromosome"))[1L]
    col_pos    <- which(cn %in% c("pos","bp","position","mbp","loc"))[1L]
    if (!is.na(col_marker)) colnames(df)[col_marker] <- "Marker"
    if (!is.na(col_chr))    colnames(df)[col_chr]    <- "Chrom"
    if (!is.na(col_pos))    colnames(df)[col_pos]    <- "Pos"
  }

  for (col in c("Chrom", "Marker", "Pos")) {
    if (!col %in% colnames(df)) {
      warning("Map column '", col, "' not found; setting to NA.")
      df[[col]] <- NA
    }
  }

  df$Pos <- .convert_pos_unit(suppressWarnings(as.numeric(df$Pos)),
                              from_unit = pos_unit, to_unit = "bp")

  # Q5: warn when bp positions look like Mbp (barley genome is ~5.1 Gbp)
  if (pos_unit == "bp" && !all(is.na(df$Pos))) {
    pos_max <- max(df$Pos, na.rm = TRUE)
    if (is.finite(pos_max) && pos_max < 10000) {
      warning(
        "Maximum marker position is ", pos_max, " (pos_unit = \"bp\"). ",
        "Values this small look like Mbp scale. ",
        "If positions are in megabase-pairs, set pos_unit = \"mbp\"."
      )
    }
  }

  df[, c("Chrom", "Marker", "Pos"), drop = FALSE]
}


#' Convert genomic position values between units
#'
#' All map positions are stored internally in bp. Use this when reading or
#' writing maps that use megabase-pair coordinates.
#'
#' @param pos       Numeric vector.
#' @param from_unit \code{"bp"} or \code{"mbp"}.
#' @param to_unit   \code{"bp"} or \code{"mbp"}.
#' @return Numeric vector.
#' @keywords internal
.convert_pos_unit <- function(pos, from_unit = "bp", to_unit = "bp") {
  if (identical(from_unit, to_unit)) return(pos)
  pos <- suppressWarnings(as.numeric(pos))
  if (from_unit == "mbp" && to_unit == "bp")  return(pos * 1e6)
  if (from_unit == "bp"  && to_unit == "mbp") return(pos / 1e6)
  stop("Unsupported unit conversion: '", from_unit, "' -> '", to_unit, "'")
}


#' Warn on marker-set mismatches between genotype array and map
#' @keywords internal
.validate_geno_map <- function(geno_array, map_df) {
  geno_markers <- dimnames(geno_array)[[2L]]
  map_markers  <- map_df$Marker
  miss_map  <- setdiff(geno_markers, map_markers)
  miss_geno <- setdiff(map_markers, geno_markers)
  if (length(miss_map)  > 0L)
    warning(length(miss_map), " genotype marker(s) absent from map.")
  if (length(miss_geno) > 0L)
    warning(length(miss_geno), " map marker(s) absent from genotype.")
  invisible(NULL)
}
