# =============================================================================
# format_transform.R
#
# transform_genotype() — convert processed metadata to any target format.
# read_effect_file()   — read and standardise marker effect tables.
# transform_effects()  — re-encode allele states in an effect data frame.
# =============================================================================

#' Convert genotype metadata to a target output format
#'
#' Takes the metadata list produced by \code{\link{process_genotype}} and
#' reformats the genotype data to a different encoding, orientation, or
#' output preset (e.g. SelectionTools → PLINK dosage).
#'
#' @param metadata           Named list from \code{\link{process_genotype}}.
#' @param output_preset      Short-cut preset: \code{"ST"}, \code{"GS"},
#'   \code{"MA"}, \code{"PL"}, or \code{NULL} (manual specification).
#' @param output_orientation \code{1} (rows = samples) or \code{2} (rows =
#'   markers); ignored when \code{output_preset} is set.
#' @param output_split_char  Allele separator for the output; ignored when
#'   preset is set or when \code{output_encode = "allele"}.
#' @param output_char_dict   Encoding dictionary for output (default: same as
#'   input).
#' @param output_encode      \code{"nucleotide"}, \code{"selectiontools"},
#'   or \code{"allele"}; ignored when preset is set.
#' @param name_as_col        \code{TRUE} — marker/sample name as first column;
#'   \code{FALSE} — stored as row/column names; \code{NULL} — preset default.
#' @param output_pos_unit    Position unit for the output map: \code{"bp"}
#'   (default) or \code{"mbp"}.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{genotype}}{data.frame ready to write.}
#'   \item{\code{map}}{data.frame formatted for the output preset/unit.}
#' }
#'
#' @seealso \code{\link{process_genotype}}
#'
#' @examples
#' # Round-trip ST → nucleotide → ST
#' geno <- matrix(c("1/1","1/4","4/4","1/1"), nrow = 2,
#'                dimnames = list(c("Mkr1","Mkr2"), c("S1","S2")))
#' meta <- process_genotype(geno, preset = "ST")
#' out  <- transform_genotype(meta, output_preset = "GS")
#' head(out$genotype)
#'
#' @export
transform_genotype <- function(metadata,
                                output_preset      = NULL,
                                output_orientation = NULL,
                                output_split_char  = NULL,
                                output_char_dict   = NULL,
                                output_encode      = NULL,
                                name_as_col        = NULL,
                                output_pos_unit    = "bp") {

  if (is.null(output_char_dict)) output_char_dict <- metadata$char_dict

  out_p              <- .resolve_preset(output_preset, output_orientation,
                                        output_split_char, output_encode)
  output_orientation <- out_p$orientation
  output_split_char  <- out_p$split_char
  output_encode      <- out_p$encode

  geno_array   <- metadata$genotype
  input_encode <- metadata$encode
  map_df       <- metadata$map
  sample_names <- metadata$sample

  # Translate encoding if needed
  needs_translate <- !is.null(input_encode) && !is.null(output_encode) &&
    !identical(input_encode, output_encode)

  if (needs_translate) {
    minor_alleles <- if (!is.null(map_df) && "Minor" %in% colnames(map_df)) {
      setNames(map_df$Minor, map_df$Marker)
    } else NULL

    geno_array <- .genotype_translate(
      geno_array    = geno_array,
      from_encode   = input_encode,
      to_encode     = output_encode,
      char_dict     = output_char_dict,
      minor_alleles = minor_alleles
    )
  }

  out_df <- .combine_genotype(
    geno_array         = geno_array,
    map_df             = map_df,
    sample_names       = sample_names,
    output_preset      = output_preset,
    output_orientation = output_orientation,
    output_split_char  = output_split_char,
    output_encode      = output_encode,
    output_char_dict   = output_char_dict,
    name_as_col        = name_as_col
  )

  out_map <- .format_map_for_preset(map_df, output_preset, output_pos_unit)
  list(genotype = out_df, map = out_map)
}


#' Read a marker effect file into a standardised data frame
#'
#' Reads a whitespace- or tab-separated effect file and returns a data frame
#' with columns \code{Marker}, \code{Allele}, and \code{Effect}. Additional
#' trait columns are preserved as \code{Effect2}, \code{Effect3}, etc.
#'
#' Allele encoding is auto-detected (SelectionTools numeric or nucleotide
#' characters) and can optionally be converted to a target encoding.
#'
#' @param path          Path to the effect file. Should have at least three
#'   columns: marker ID, allele state, effect size.
#' @param sep           Column separator. Default \code{" "} (space). Use
#'   \code{"\t"} for tab-separated files.
#' @param header        Logical; whether the file has a header row. Default
#'   \code{FALSE} (effect files typically have no header).
#' @param output_encode Optional target encoding: \code{"nucleotide"} or
#'   \code{"selectiontools"}. When \code{NULL} (default) the original encoding
#'   is preserved.
#' @param char_dict     Encoding dictionary (default:
#'   \code{DEFAULT_ENCODING_DICT}).
#' @param extend        Logical. When \code{TRUE} (default) fills in zero-effect
#'   rows for any (Marker, Allele) combination missing from the file.
#'
#' @return data.frame with columns \code{Marker}, \code{Allele},
#'   \code{Effect} [, \code{Effect2}, ...].
#'
#' @seealso \code{\link{transform_effects}}
#'
#' @examples
#' \dontrun{
#' eff <- read_effect_file("effects/trait1.GS.Eff")
#' head(eff)
#' }
#'
#' @export
read_effect_file <- function(path,
                              sep           = " ",
                              header        = FALSE,
                              output_encode = NULL,
                              char_dict     = DEFAULT_ENCODING_DICT,
                              extend        = TRUE) {
  if (!file.exists(path)) stop("read_effect_file: File not found: ", path)

  # Read as text, then coerce the effect columns back. The allele column (2)
  # must never be type-inferred: a table whose alleles are all "T" came back as
  # logical TRUE, and an allele spelled "TRUE" matches no genotype, so every
  # effect silently contributed 0 — no error, just a flat GEBV. See
  # .read_table_chr(). Column 1 (Marker) stays character so an id like "0100"
  # keeps its spelling; columns 3+ are the numeric effects.
  raw <- .read_table_chr(path, header = header, sep = sep)
  if (ncol(raw) >= 3L)
    for (j in 3:ncol(raw))
      raw[[j]] <- suppressWarnings(as.numeric(raw[[j]]))

  transform_effects(
    effects_df    = raw,
    metadata      = NULL,
    output_encode = output_encode,
    char_dict     = char_dict,
    extend        = extend
  )
}


#' Re-encode allele states in a (Marker, Allele, Effect) data frame
#'
#' Converts the allele state column from one encoding to another and optionally
#' fills in zero-effect rows for any (Marker, Allele) combination absent from
#' the input.
#'
#' @param effects_df    data.frame with columns: Marker (col 1), Allele state
#'   (col 2), Effect (col 3). Extra numeric columns are treated as additional
#'   traits.
#' @param metadata      Optional list from \code{\link{process_genotype}} —
#'   used for map-order sorting and dictionary fallback.
#' @param input_encode  Encoding of the State column: \code{"nucleotide"} or
#'   \code{"selectiontools"}. Auto-detected when \code{NULL}.
#' @param output_encode Target encoding: \code{"nucleotide"} or
#'   \code{"selectiontools"}. When \code{NULL} the original encoding is
#'   preserved.
#' @param char_dict     Encoding dictionary.
#' @param extend        When \code{TRUE} (default) adds zero-effect rows for
#'   missing (Marker, Allele) combinations.
#'
#' @return data.frame with columns \code{Marker}, \code{Allele},
#'   \code{Effect} [, extra effect columns].
#'
#' @export
transform_effects <- function(effects_df,
                               metadata      = NULL,
                               input_encode  = NULL,
                               output_encode = NULL,
                               char_dict     = NULL,
                               extend        = TRUE) {

  if (is.null(char_dict)) {
    char_dict <- if (!is.null(metadata$char_dict)) metadata$char_dict
                 else DEFAULT_ENCODING_DICT
  }

  df <- as.data.frame(effects_df, stringsAsFactors = FALSE)
  if (ncol(df) < 3L) stop("transform_effects: need at least 3 columns.")

  colnames(df)[1:3] <- c("Marker", "Allele", "Effect")
  extra_cols <- if (ncol(df) > 3L) seq(4L, ncol(df)) else integer(0L)

  df$Marker <- as.character(df$Marker)
  df$Allele <- as.character(df$Allele)
  df$Effect <- as.numeric(df$Effect)

  # Translate allele encoding if requested
  if (!is.null(input_encode) && !is.null(output_encode) &&
      !identical(input_encode, output_encode)) {

    state_map <- .eff_build_state_map(input_encode, output_encode, char_dict)
    translated <- state_map[df$Allele]
    bad <- is.na(translated) & !is.na(df$Allele)
    if (any(bad)) translated[bad] <- df$Allele[bad]
    df$Allele <- as.character(translated)
  }

  # Extend with zero rows for missing allele × marker combinations
  if (extend && !is.null(char_dict[["valid_chars"]])) {
    valid_states <- if (!is.null(output_encode) &&
                        output_encode == "selectiontools")
      as.character(unlist(char_dict[vapply(names(char_dict), function(k)
        k %in% c("valid_nums","A","G","C","T"), logical(1L))]))
    else char_dict[["valid_chars"]]

    if (!is.null(valid_states) && length(valid_states) > 0L) {
      all_markers <- unique(df$Marker)
      fill_rows <- lapply(all_markers, function(m) {
        present <- df$Allele[df$Marker == m]
        missing  <- setdiff(valid_states, present)
        if (length(missing) == 0L) return(NULL)
        data.frame(Marker = m, Allele = missing, Effect = 0,
                   stringsAsFactors = FALSE)
      })
      fill_df <- do.call(rbind, Filter(Negate(is.null), fill_rows))
      if (!is.null(fill_df) && nrow(fill_df) > 0L) {
        df <- rbind(df[, c("Marker","Allele","Effect")], fill_df)
      }
    }
  }

  rownames(df) <- NULL
  df
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Build allele state translation map
#' @keywords internal
.eff_build_state_map <- function(from_encode, to_encode, char_dict) {
  from_encode <- .encode_arg(from_encode, "transform_effects")
  to_encode   <- .encode_arg(to_encode,   "transform_effects")
  if (from_encode == "nucleotide" && to_encode == "selectiontools") {
    m <- setNames(as.character(unlist(DEFAULT_ENCODING_DICT)),
                  names(DEFAULT_ENCODING_DICT))
    return(m)
  }
  if (from_encode == "selectiontools" && to_encode == "nucleotide") {
    return(REVERSED_ENCODING_DICT)
  }
  stop("transform_effects: unsupported encoding '", from_encode,
       "' → '", to_encode, "'.")
}


#' Format map for target preset / position unit
#' @keywords internal
.format_map_for_preset <- function(map_df, output_preset = NULL,
                                   output_pos_unit = "bp") {
  if (is.null(map_df)) return(NULL)
  map_out <- map_df
  map_out$Pos <- .convert_pos_unit(map_df$Pos, from_unit = "bp",
                                   to_unit = output_pos_unit)

  if (!is.null(output_preset) && output_preset == "ST") {
    df <- map_out[, intersect(c("Marker","Chrom","Pos"), colnames(map_out))]
    colnames(df) <- c("name","chr","pos")[seq_along(colnames(df))]
    return(df)
  }
  if (!is.null(output_preset) && output_preset == "PL") {
    return(data.frame(chr  = map_out$Chrom,
                      name = map_out$Marker,
                      cM   = 0,
                      pos  = map_out$Pos,
                      stringsAsFactors = FALSE))
  }
  map_out[, intersect(c("Chrom","Marker","Pos","Minor"), colnames(map_out)),
          drop = FALSE]
}


#' Placeholder for .genotype_translate (implemented in format_split.R)
#' @keywords internal
.genotype_translate <- function(geno_array, from_encode, to_encode,
                                char_dict, minor_alleles = NULL) {
  # Implemented in format_split.R via translate_alleles()
  translate_alleles(geno_array, from_encode, to_encode, char_dict,
                    minor_alleles)
}


#' Placeholder for .combine_genotype (implemented in format_split.R)
#' @keywords internal
.combine_genotype <- function(geno_array, map_df, sample_names,
                               output_preset, output_orientation,
                               output_split_char, output_encode,
                               output_char_dict, name_as_col) {
  combine_haplotypes(geno_array, sample_names, output_split_char,
                     output_orientation)
}
