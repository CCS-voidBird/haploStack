# =============================================================================
# format_detect.R
#
# Functions for auto-detecting genotype matrix properties: encoding format,
# orientation (rows = samples vs rows = markers), and minor alleles per marker.
# =============================================================================

#' Detect genotype encoding format
#'
#' Samples up to 100 non-NA values from a genotype matrix and classifies the
#' encoding as \code{"nucleotide"} (allele labels, e.g. \code{"A/T"} or
#' \code{"AT"}) or \code{"numeric"} (dosage counts 0/1/2). Slash-separated
#' integer codes (\code{"1/4"}, \code{"1/-1"}) are SelectionTools' allele
#' encoding, which a \code{GenoMeta} no longer stores (0.9.0): they are
#' refused with an \code{hs_error_format} that names
#' \code{\link{st_convert}} and the two declarations that resolve the
#' ambiguity (\code{encoding = "selectiontools"} to translate the codes,
#' \code{encoding = "nucleotide"} to keep digits as labels).
#'
#' @param geno_data A matrix, data.frame, or data.table. Values should be
#'   genotype call strings. Either orientation (markers × samples or samples ×
#'   markers) is accepted.
#'
#' @return A character scalar: \code{"nucleotide"} or \code{"numeric"}.
#'
#' @details
#' Detection priority:
#' 1. If >50 % of sampled values contain \code{"/"} and the allele parts are
#'    all numeric → SelectionTools codes: refused (see above).
#' 2. If >50 % contain \code{"/"} and allele parts are alphabetic → \code{"nucleotide"}.
#' 3. If no slash but >50 % contain A/G/C/T characters → \code{"nucleotide"}.
#' 4. If values are simple integers 0/1/2 → \code{"numeric"}.
#'
#' @examples
#' # Nucleotide format
#' detect_genotype_encoding(matrix(c("A/A", "A/T", "T/T"), 1))
#' # => "nucleotide"
#'
#' # Dosage numeric
#' detect_genotype_encoding(matrix(c(0, 1, 2, 1), 2))
#' # => "numeric"
#'
#' # SelectionTools codes are refused, naming st_convert()
#' try(detect_genotype_encoding(matrix(c("1/1", "1/4", "4/4", "-1/-1"), 2)))
#'
#' @export
detect_genotype_encoding <- function(geno_data) {
  enc <- .detect_encoding(geno_data)
  if (identical(enc, "selectiontools")) .st_refuse_codes("detect_genotype_encoding")
  enc
}

# The detector proper. Still returns "selectiontools" for slash-separated
# integer codes, so the portal (new_geno_meta()), build_geno_meta() and the
# legacy process_genotype() layer can tell codes from labels; the exported
# wrapper above refuses that verdict.
# @keywords internal
.detect_encoding <- function(geno_data) {
  if (is.data.frame(geno_data) || inherits(geno_data, "data.table")) {
    geno_mat <- as.matrix(geno_data)
  } else {
    geno_mat <- geno_data
  }

  all_vals <- as.character(as.vector(geno_mat))
  all_vals <- all_vals[!is.na(all_vals) & nzchar(all_vals) & all_vals != "NA"]

  if (length(all_vals) == 0L) {
    .hs_abort("detect_genotype_encoding: No non-NA genotype values found.")
  }

  # Classify on the calls, not the absences. A cell whose every allele is a
  # missing token ("00", "--", "0/0", "./.") says nothing about which
  # vocabulary the file is written in, and on a panel with a high missing rate
  # enough of them will dilute the ratios below to misclassify it — a
  # heavily-missing PED reads as "numeric" without this filter. The union of
  # every mode's missing set is used because the mode is what we are trying to
  # determine; the filter only ever *removes* candidates, and if it removes
  # everything the unfiltered values are used instead.
  all_vals <- .drop_all_missing_cells(all_vals)

  # A fixed local seed keeps the verdict reproducible; the caller's own
  # random stream is restored afterwards (R/state.R).
  sampled <- .with_seed(42L, sample(all_vals, min(100L, length(all_vals))))

  has_slash  <- grepl("/", sampled, fixed = TRUE)

  if (mean(has_slash) > 0.5) {
    slash_vals <- sampled[has_slash]
    parts <- unlist(strsplit(slash_vals, "/", fixed = TRUE))
    if (mean(grepl("[AGCTagct]", parts)) > 0.5) {
      return("nucleotide")
    } else {
      return("selectiontools")
    }
  }

  if (mean(grepl("[AGCTagct]", sampled)) > 0.5) return("nucleotide")
  if (mean(grepl("^[012]$", trimws(sampled))) > 0.5) return("numeric")

  nums <- suppressWarnings(as.numeric(trimws(sampled)))
  if (mean(!is.na(nums)) > 0.5) return("numeric")

  .hs_abort(
    "detect_genotype_encoding: Unable to determine encoding. ",
    "Sample values: ", paste(head(sampled, 5L), collapse = ", ")
  )
}


# ── Internal helpers ──────────────────────────────────────────────────────────

# Internal: drop genotype cells whose every allele is a missing token.
#
# Works on the unique values rather than the full vector — an allele vocabulary
# is tiny next to a 20-million-cell panel, so this is one small loop plus one
# %in%. A cell splits on "/" when it has one and into characters otherwise,
# which is the same pair of writings split_haplotypes() accepts.
#
# Returns the input unchanged when the filter would empty it (an all-missing
# marker set, or a dosage panel that happens to be all zeros).
# @keywords internal
.drop_all_missing_cells <- function(vals) {
  if (length(vals) == 0L) return(vals)

  miss <- unique(c(.missing_tokens("nucleotide"),
                   .missing_tokens("selectiontools"),
                   .missing_tokens("numeric")))

  u <- unique(vals)
  all_missing <- vapply(u, function(v) {
    if (v %in% miss) return(TRUE)
    parts <- if (grepl("/", v, fixed = TRUE))
               strsplit(v, "/", fixed = TRUE)[[1L]]
             else
               strsplit(v, "", fixed = TRUE)[[1L]]
    length(parts) > 0L && all(parts %in% miss)
  }, logical(1L), USE.NAMES = FALSE)

  if (!any(all_missing)) return(vals)
  keep <- vals[!vals %in% u[all_missing]]
  if (length(keep) == 0L) vals else keep
}


#' Detect whether rows are samples or markers
#'
#' Uses map row count for an exact check; falls back to aspect-ratio heuristic.
#'
#' @param df        Data frame or matrix of genotype data.
#' @param map_info  Optional map data frame (one row per marker).
#' @return Integer: \code{1L} (rows = samples) or \code{2L} (rows = markers).
#' @keywords internal
#' @noRd
.detect_orientation <- function(df, map_info = NULL) {
  if (!is.null(map_info)) {
    n_markers <- nrow(map_info)
    if (nrow(df) %in% c(n_markers, n_markers + 1L)) return(2L)
    if (ncol(df) %in% c(n_markers, n_markers + 1L)) return(1L)
  }
  if (nrow(df) > ncol(df)) return(2L)
  return(1L)
}


#' Detect the minor allele for every marker in a 3D genotype array
#'
#' The minor allele is the non-NA allele with the lowest count across all
#' samples and haplotype copies. For dosage (\code{"allele"}) encoding the
#' minor allele is always \code{"1"}.
#'
#' @param geno_array 3D array \code{(n_samples, n_markers, ploidy)}.
#' @param encode     Storage encoding: \code{"nucleotide"},
#'   \code{"selectiontools"}, or \code{"allele"}.
#' @param char_dict  Encoding dictionary; used to identify missing value codes.
#' @return Named character vector: marker name → minor allele (as stored).
#' @keywords internal
#' @noRd
.detect_minor_alleles <- function(geno_array,
                                  encode    = "nucleotide",
                                  char_dict = DEFAULT_ENCODING_DICT) {
  n_markers    <- dim(geno_array)[2L]
  marker_names <- dimnames(geno_array)[[2L]]
  if (is.null(marker_names)) marker_names <- paste0("M", seq_len(n_markers))

  na_vals <- unique(c(NA_character_,
                      as.character(char_dict[["na_char"]]),
                      as.character(char_dict[["na_num"]])))

  if (encode == "allele") {
    return(setNames(rep("1", n_markers), marker_names))
  }

  # Minor allele = least frequent OBSERVED allele per marker (ties broken
  # alphabetically, matching names(sort(table(alleles)))[1L]). Detection of the
  # missing tokens + dictionary is in R; the kernel only counts codes.
  counts <- .marker_allele_counts(geno_array, na_vals = na_vals)
  levs   <- colnames(counts)
  na_char <- as.character(char_dict[["na_char"]])
  result <- vapply(seq_len(n_markers), function(j) {
    cj <- counts[j, ]
    if (all(cj == 0L)) return(na_char)
    cj[cj == 0L] <- .Machine$integer.max          # ignore alleles unobserved here
    levs[which.min(cj)]                            # least frequent; ties -> alpha
  }, character(1L))

  setNames(result, marker_names)
}


# Canonical allele-encoding vocabulary, shared by every layer of the format
# module: the storage encodings that translate_alleles() converts between, and
# the same tokens detect_genotype_encoding() emits and metadata$encoding stores.
#
# Before 0.5.1 there were two vocabularies for one concept: the low-level
# "char" / "number" used by .resolve_preset() / translate_alleles() /
# transform_effects(), and the read-side "nucleotide" / "selectiontools" used by
# detect_genotype_encoding(), GENO_FORMAT_PRESETS and geno_meta$metadata. They
# met in process_genotype(), which overwrites the preset's encode with
# split_haplotypes()'s detected one — so an ST preset produced "selectiontools"
# and then handed it to translate_alleles(), which only knew "number". That is
# what broke the transform_genotype() round-trip example. One vocabulary now.
#
# 0.9.0 (C-5): "selectiontools" is no longer a vocabulary the container
# stores — a GenoMeta holds nucleotide alleles (and dosage). The token stays
# accepted here for the legacy process_genotype() / translate_alleles() layer
# (deprecated, warns) and for the declared-encoding deprecation rung of the
# builders; it joins .HS_ENCODING_REMOVED in 0.10.
HS_ENCODINGS <- c("selectiontools", "nucleotide", "allele")

# Tokens removed in 0.5.1, mapped to what replaced them, so a caller still on
# the old vocabulary gets told what to write instead of a downstream
# "unsupported encoding" error from whichever helper happens to see it first.
.HS_ENCODING_REMOVED <- c(number = "selectiontools", char = "nucleotide")

# Internal: validate an encoding argument. `what` names the calling function for
# the error message. NULL passes through — callers use it for "not specified".
#
# A removed token is ALWAYS an error, whatever `allowed` and `strict` say: that
# is what makes the removal loud instead of surfacing later as an "unsupported
# encoding" from whichever helper happens to see it first.
#
# `strict = FALSE` lets an unrecognised value through untouched, for the one
# caller (build_geno_meta) that legitimately accepts detector output it does not
# have a dictionary for and falls back to deriving one from the observed alleles.
# @keywords internal
.encode_arg <- function(x, what, allowed = HS_ENCODINGS, strict = TRUE) {
  if (is.null(x)) return(NULL)
  x <- as.character(x)
  if (length(x) != 1L || is.na(x))
    .hs_abort(what, ": encoding must be a single string, one of ",
         paste0('"', allowed, '"', collapse = ", "), ".")
  if (x %in% names(.HS_ENCODING_REMOVED))
    .hs_abort(what, ": encoding \"", x, "\" was removed in 0.5.1 \u2014 use \"",
         .HS_ENCODING_REMOVED[[x]], "\" instead.")
  if (x %in% allowed || !strict) return(x)
  .hs_abort(what, ": unknown encoding \"", x, "\". Use one of ",
       paste0('"', allowed, '"', collapse = ", "), ".")
}


#' Resolve format preset parameters
#'
#' Maps a short preset name to the canonical orientation, split character, and
#' encoding values used internally by \code{process_genotype} /
#' \code{transform_genotype}.
#'
#' @param preset      One of \code{"ST"}, \code{"GS"}, \code{"MA"},
#'   \code{"PL"}, or \code{NULL}.
#' @param orientation Integer 1 or 2, or \code{NULL} to auto-detect.
#' @param split_char  Allele separator string or \code{NULL}.
#' @param encode      Encoding name (\code{"selectiontools"},
#'   \code{"nucleotide"}, \code{"allele"}) or \code{NULL}.
#' @return Named list with elements \code{orientation}, \code{split_char},
#'   \code{encode}.
#' @keywords internal
#' @noRd
.resolve_preset <- function(preset, orientation, split_char, encode) {
  encode <- .encode_arg(encode, ".resolve_preset")
  if (is.null(preset)) {
    return(list(orientation = orientation,
                split_char  = split_char %||% "/",
                encode      = encode %||% "nucleotide"))
  }

  preset_map <- list(
    ST  = list(orientation = 2L, split_char = "/",  encode = "selectiontools"),
    GS  = list(orientation = 2L, split_char = NULL, encode = "nucleotide"),
    MA  = list(orientation = 1L, split_char = "/",  encode = "nucleotide"),
    PL  = list(orientation = 1L, split_char = NULL, encode = "allele")
  )

  p <- preset_map[[toupper(preset)]]
  if (is.null(p)) {
    .hs_abort(".resolve_preset: Unknown preset '", preset,
         "'. Choose from: ST, GS, MA, PL.")
  }
  list(
    orientation = p$orientation,
    split_char  = p$split_char %||% split_char %||% "/",
    encode      = p$encode %||% encode %||% "nucleotide"
  )
}


