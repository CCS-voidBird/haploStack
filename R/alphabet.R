# =============================================================================
# alphabet.R
#
# Package-level constants: allele dictionaries, encoding maps, and color
# palettes used throughout haploStack. These objects are exported so users
# can inspect or extend them when integrating custom allele encodings.
# =============================================================================

#' Default two-allele dictionary
#'
#' Character vector of the two canonical allelic states used in standard
#' biallelic SNP panels (A = reference, T = alternate).
#'
#' @export
DEFAULT_ALLELE_DICT <- c("A", "T")

#' Extended four-allele dictionary
#'
#' All four nucleotide bases. Used when working with multi-allelic markers or
#' re-sequencing panels.
#'
#' @export
EXTENDED_ALLELE_DICT <- c("A", "T", "C", "G")

#' Default numeric encoding dictionary
#'
#' Named list mapping nucleotide characters to their SelectionTools integer
#' codes. Missing alleles ("-") map to -1L.
#'
#' @details
#' The canonical nucleotide dictionary (\code{A, C, G, T, -}) maps to the
#' SelectionTools integer codes in \strong{alphabetical} order:
#'
#' | Character | Code |
#' |-----------|------|
#' | A         | 1    |
#' | C         | 2    |
#' | G         | 3    |
#' | T         | 4    |
#' | -         | -1   |
#'
#' \code{recode_geno_meta(to_type = "selectiontools")} therefore encodes
#' \code{A -> 1, C -> 2, G -> 3, T -> 4} (missing \code{-> -1}); nucleotide
#' \code{A, C, G, T, -} is the primary dictionary and SelectionTools numeric codes
#' are the alternate encoding (used for tools that ingest numeric genotypes, e.g.
#' LD estimation).
#'
#' @export
DEFAULT_ENCODING_DICT <- list(
  "A" = 1L,
  "C" = 2L,
  "G" = 3L,
  "T" = 4L,
  "-" = -1L
)

#' Default placeholder for a missing allele in a haplotype signature
#'
#' The character written in place of an unobserved allele when a block's
#' haplotype signature is assembled (see \code{\link{write_block_allele_sheet}},
#' \code{\link{count_block_subtypes}}). \code{"-"} is the canonical missing-allele
#' character of \code{\link{DEFAULT_ENCODING_DICT}}, where it maps to \code{-1L}.
#'
#' @details
#' A signature is the alleles of a haplotype copy pasted together in genomic
#' order, so it is only a faithful identifier if every member marker occupies a
#' position. Dropping a missing allele (or writing an empty string for it) makes
#' the string \strong{shorter} and loses the position, so distinct combinations
#' collapse onto one another: over markers \code{(M1, M2, M3)}, the copies
#' \code{(A, T, -)} and \code{(A, -, T)} both paste to \code{"AT"}. Filling with a
#' placeholder keeps every signature the width of its block and keeps distinct
#' combinations distinct.
#'
#' Override globally with \code{options(haploStack.fill_na = "?")}, or per call
#' via the \code{fill_na} argument.
#'
#' @seealso \code{\link{DEFAULT_ENCODING_DICT}},
#'   \code{\link{write_block_allele_sheet}}, \code{\link{count_block_subtypes}}
#' @export
HS_FILL_NA_DEFAULT <- "-"

# Resolve the active missing-allele placeholder: the `haploStack.fill_na` option
# when set (see zzz.R), else HS_FILL_NA_DEFAULT. Must be a single non-NA string;
# an empty string is refused because it reintroduces the collision the
# placeholder exists to prevent.
# @keywords internal
.hs_fill_na <- function(fill_na = NULL) {
  x <- fill_na %||% getOption("haploStack.fill_na", HS_FILL_NA_DEFAULT)
  if (!is.character(x) || length(x) != 1L || is.na(x))
    stop("fill_na must be a single non-NA character string; got ",
         paste(utils::capture.output(utils::str(x)), collapse = " "), ".")
  if (!nzchar(x))
    stop('fill_na must not be the empty string — an empty placeholder shortens ',
         'the signature and lets distinct allele combinations collide. Use "-" ',
         "(the default) or another single character.")
  x
}

#' Reversed encoding dictionary
#'
#' Inverse of \code{DEFAULT_ENCODING_DICT}: maps integer codes back to
#' nucleotide character strings. Useful for converting SelectionTools numeric
#' genotypes back to nucleotide format.
#'
#' @export
REVERSED_ENCODING_DICT <- setNames(
  names(DEFAULT_ENCODING_DICT),
  as.character(unname(unlist(DEFAULT_ENCODING_DICT)))
)

#' Single-allele base color map
#'
#' Named list with two encoding keys (\code{"nucleotide"} and
#' \code{"selectiontools"}). Each element is a named character vector mapping
#' individual allele values to hex color codes.
#'
#' @export
BASE_COLOR_MAP <- list(
  "nucleotide" = c(
    "A" = "#6a3d9a",  # purple
    "G" = "#009E73",  # green
    "C" = "#1f78b4",  # blue
    "T" = "#ff7f00",  # orange
    "-" = "#D3D3D3"   # light grey (missing)
  ),
  "selectiontools" = c(          # codes A=1, C=2, G=3, T=4 (mirror the nucleotide colours)
    "1"  = "#6a3d9a",            # A purple
    "2"  = "#1f78b4",            # C blue
    "3"  = "#009E73",            # G green
    "4"  = "#ff7f00",            # T orange
    "0"  = "#D3D3D3",
    "-1" = "#D3D3D3",
    "-"  = "#D3D3D3"
  )
)

#' Diploid genotype color map
#'
#' Named list with two encoding keys (\code{"nucleotide"} and
#' \code{"selectiontools"}). Each element maps two-allele genotype strings
#' (e.g. \code{"11"}, \code{"14"}, \code{"AA"}, \code{"AT"}) to hex colors.
#' Homozygotes are saturated; heterozygotes are intermediate hues; missing
#' calls are light grey.
#'
#' @export
COLOR_MAP <- list(
  "selectiontools" = c(          # codes A=1, C=2, G=3, T=4 (mirror the nucleotide genotypes)
    "00"  = "#D3D3D3",
    "01"  = "#c9aede", "11" = "#6a3d9a",   # -A / AA
    "02"  = "#80d4b8", "22" = "#009E73",   # -C / CC
    "03"  = "#a6c8e8", "33" = "#1f78b4",   # -G / GG
    "04"  = "#ffd199", "44" = "#ff7f00",   # -T / TT
    "12"  = "#3d7a6a",                     # AC
    "13"  = "#3d5a8a",                     # AG
    "14"  = "#b35e00",                     # AT
    "23"  = "#0a6678",                     # CG
    "24"  = "#a05a1f",                     # CT
    "34"  = "#7a8a00"                      # GT
  ),
  "nucleotide" = c(
    "--"  = "#D3D3D3",
    "-A"  = "#c9aede", "AA" = "#6a3d9a",
    "-G"  = "#a6c8e8", "GG" = "#1f78b4",
    "-C"  = "#80d4b8", "CC" = "#009E73",
    "-T"  = "#ffd199", "TT" = "#ff7f00",
    "AG"  = "#3d5a8a",
    "AC"  = "#3d7a6a",
    "AT"  = "#b35e00",
    "CG"  = "#0a6678",
    "CT"  = "#a05a1f",
    "GT"  = "#7a8a00"
  )
)


#' Build a discrete allele / genotype colour map for observed values
#'
#' Chooses the single-allele palette (\code{\link{BASE_COLOR_MAP}}) when every
#' observed value is one character (haplotype view) or the two-allele genotype
#' palette (\code{\link{COLOR_MAP}}) when every value is two characters (diploid
#' genotype view), keyed by \code{encoding}. Observed values absent from the
#' chosen palette receive distinct fallback colours; a mix of one- and
#' two-character values is an error.
#'
#' @param allele_values Character vector of observed allele / genotype values.
#' @param encoding      \code{"nucleotide"}, \code{"selectiontools"} or
#'   \code{"numeric"}.
#'
#' @return A named character vector mapping each observed value to a hex colour,
#'   suitable for \code{ggplot2::scale_fill_manual(values = ...)}.
#'
#' @export
make_allele_color_map <- function(allele_values, encoding) {
  u <- sort(unique(allele_values[!is.na(allele_values) & nzchar(allele_values)]))
  if (length(u) == 0L) return(character(0))

  if (identical(encoding, "numeric")) {
    cm <- c("0" = "#f0f0f0", "1" = "#6baed6", "2" = "#08519c")
    cm <- cm[names(cm) %in% u]
  } else {
    lens <- nchar(u)
    base <- if (all(lens == 1L)) {
      BASE_COLOR_MAP
    } else if (all(lens == 2L)) {
      COLOR_MAP
    } else {
      stop("make_allele_color_map: mixed allele-string lengths (",
           paste(sort(unique(lens)), collapse = ", "),
           "); all values must be single alleles or all two-allele genotypes.")
    }
    palette <- if (encoding %in% names(base)) base[[encoding]]
               else base[["nucleotide"]]
    cm <- palette[names(palette) %in% u]
  }

  # Distinct fallback colours for any observed value not in the palette
  missing <- setdiff(u, names(cm))
  if (length(missing) > 0L) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Dark 3")
    names(extra) <- missing
    cm <- c(cm, extra)
  }
  cm
}
