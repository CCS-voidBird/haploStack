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

# The SelectionTools code table (DEFAULT_ENCODING_DICT, REVERSED_ENCODING_DICT)
# moved to R/st_formats.R in 0.9.0 (C-5): a GenoMeta stores nucleotide alleles
# only, and the table is applied at the SelectionTools boundary.

#' Default placeholder for a missing allele in a haplotype signature
#'
#' The character written in place of an unobserved allele when a block's
#' haplotype signature is assembled (see \code{\link{write_block_allele_sheet}},
#' \code{\link{count_block_subtypes}}). \code{"-"} is the canonical missing-allele
#' character of the nucleotide vocabulary.
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
#' @seealso \code{\link{write_block_allele_sheet}},
#'   \code{\link{count_block_subtypes}}
#' @export
HS_FILL_NA_DEFAULT <- "-"

# Resolve the active missing-allele placeholder: the call's `fill_na` argument,
# else the `haploStack.fill_na` option (see zzz.R), else HS_FILL_NA_DEFAULT
# (.hs_setting(), R/state.R). Must be a single non-NA string; an empty string
# is refused because it reintroduces the collision the placeholder exists to
# prevent.
# @keywords internal
.hs_fill_na <- function(fill_na = NULL) {
  .hs_setting("fill_na", fill_na)
}

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
  )
)

#' Diploid genotype color map
#'
#' Named list with one encoding key, \code{"nucleotide"} (the
#' \code{"selectiontools"} palette went with the encoding in 0.9.0). It maps
#' two-allele genotype strings (e.g. \code{"AA"}, \code{"AT"}) to hex colors.
#' Homozygotes are saturated; heterozygotes are intermediate hues; missing
#' calls are light grey.
#'
#' @export
COLOR_MAP <- list(
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
  u <- .hs_sort(unique(allele_values[!is.na(allele_values) & nzchar(allele_values)]))
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
      .hs_abort("make_allele_color_map: mixed allele-string lengths (",
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


# ---------------------------------------------------------------------------
# Missing-allele vocabulary (MAN, 0.6.0)
# ---------------------------------------------------------------------------

# The three source-allele vocabularies a missing character can be registered
# against. These are source *formats* (what a file's alleles mean), which is the
# same vocabulary normalise_missing_alleles() takes — not the container's
# provenance record. "allele" is accepted as an alias of "numeric" because
# .detect_minor_alleles() uses it for dosage (the split is DVF Q3, still open).
HS_MISSING_MODES <- c("nucleotide", "selectiontools", "numeric")

# Characters that are REAL alleles in a mode and can therefore never be
# registered as missing there. The numeric row is the load-bearing one: 0/1/2
# are genotypes, not absences. Digits 5-9 are barred there too — a dosage panel
# that genuinely writes 9 for missing should declare it per read
# (`extra_na = "9"`), where the caller owns the risk, rather than package-wide.
.hs_reserved_alleles <- function(mode) {
  switch(mode,
    nucleotide     = c("A", "C", "G", "T", "a", "c", "g", "t"),
    selectiontools = as.character(unlist(DEFAULT_ENCODING_DICT)),
    numeric        = as.character(0:9),
    character(0L))
}

.hs_missing_mode <- function(mode, caller) {
  if (missing(mode) || is.null(mode))
    .hs_abort(caller, ": 'mode' is required \u2014 name the source format the character ",
         "belongs to: ", paste0('"', HS_MISSING_MODES, '"', collapse = ", "),
         ". A missing character is only ever missing for the format that ",
         "writes it; \"0\" is an absence in a PED and a genotype in a dosage ",
         "matrix.", call. = FALSE)
  if (!is.character(mode) || length(mode) != 1L || is.na(mode))
    .hs_abort(caller, ": 'mode' must be a single non-NA string.", call. = FALSE)
  if (identical(mode, "allele")) mode <- "numeric"     # dosage alias
  if (!mode %in% HS_MISSING_MODES)
    .hs_abort(caller, ': unknown mode "', mode, '". Use one of ',
         paste0('"', HS_MISSING_MODES, '"', collapse = ", "), ".", call. = FALSE)
  mode
}

# The built-in missing-allele vocabulary for a source format, plus `extra`
# (characters the caller resolved: a container's stamp, a registry read at the
# portal). "0" is missing for the two allele vocabularies and never for a
# dosage one, where it is a real genotype — that exclusion is load-bearing, see
# normalise_missing_alleles(). Since 0.9.0 this function reads no option: the
# session registry (use_missing_chars()) reaches a container only where it is
# built or repaired, through .extra_missing_chars(); afterwards the container's
# metadata$missing_tokens stamp is the record every consumer uses
# (.meta_extra_missing()).
# @keywords internal
.missing_tokens <- function(encoding, extra = character(0L)) {
  base <- c("", "NA", "N", "n", ".", "-")
  if (identical(encoding, "nucleotide") ||
      identical(encoding, "selectiontools")) base <- c(base, "0")
  if (length(extra)) {
    extra <- as.character(extra)
    base  <- c(base, extra[!is.na(extra)])
  }
  unique(base)
}

# The characters a container was normalised against beyond the built-in set:
# its stamp minus the vocabulary of its recorded encoding. character(0) for
# an unstamped (pre-0.6.0 or hand-built) container.
# @keywords internal
.meta_extra_missing <- function(meta) {
  stamp <- meta$metadata$missing_tokens
  if (!is.character(stamp)) return(character(0L))
  setdiff(stamp[!is.na(stamp)], .missing_tokens(meta$metadata$encoding))
}

# The extra missing characters in force for one mode at a build or repair.
# `missing_chars` is the call's argument: NULL takes the session registry
# (use_missing_chars()), a character vector names the characters for this
# mode, a per-mode list is read like the registry. An argument is validated
# as use_missing_chars() validates (single characters, never a real allele of
# the mode); an unusable option value is ignored, because options() set by
# hand must not be able to change what is missing. Returns character(0) when
# nothing applies or the mode is NULL/unknown.
# @keywords internal
.extra_missing_chars <- function(mode, missing_chars = NULL) {
  if (is.null(mode) || !is.character(mode) || length(mode) != 1L || is.na(mode))
    return(character(0L))
  if (identical(mode, "allele")) mode <- "numeric"
  reg <- .hs_setting("missing_chars", missing_chars)
  if (is.null(reg)) return(character(0L))
  x <- if (is.character(reg)) reg else reg[[mode]]
  if (is.null(x) || !is.character(x)) return(character(0L))
  if (!is.null(missing_chars)) {
    if (length(x) == 0L) return(character(0L))
    return(.hs_check_missing_chars(x, mode, "missing_chars"))
  }
  unique(x[!is.na(x) & nchar(x) == 1L])
}

# The validation use_missing_chars() and the `missing_chars` argument share:
# a non-empty character vector of single, non-blank, non-NA characters, none
# of them a real allele of `mode`. Returns the unique characters.
.hs_check_missing_chars <- function(chars, mode, caller) {
  if (!is.character(chars) || length(chars) == 0L)
    .hs_abort(caller, ": 'chars' must be a non-empty character vector.",
         call. = FALSE)
  if (anyNA(chars))
    .hs_abort(caller, ": 'chars' must not contain NA — NA is already the ",
         "one representation of missing inside a GenoMeta.", call. = FALSE)

  bad_len <- chars[nchar(chars) != 1L]
  if (length(bad_len) > 0L)
    .hs_abort(caller, ": each character must be exactly one character; ",
         "got ", paste0('"', bad_len, '"', collapse = ", "),
         ". Missing is decided per allele, so a multi-character token would ",
         "never match a single allele.", call. = FALSE)

  blank <- chars[!nzchar(trimws(chars))]
  if (length(blank) > 0L)
    .hs_abort(caller, ": blank and whitespace-only cells are already ",
         "missing under every mode — nothing to register.", call. = FALSE)

  reserved <- intersect(chars, .hs_reserved_alleles(mode))
  if (length(reserved) > 0L)
    .hs_abort(caller, ": ", paste0('"', reserved, '"', collapse = ", "),
         if (length(reserved) > 1L) " are real alleles of " else " is a real allele of ",
         '"', mode, '" and cannot mean missing there.',
         if (identical(mode, "numeric"))
           " For a dosage panel that genuinely writes a digit for an absent call, declare it on the read instead (extra_na=), where the caller owns the risk."
         else "", call. = FALSE)
  unique(chars)
}


#' Register an extra missing-allele character for a source format
#'
#' Attaches one or more characters to the package's missing-allele vocabulary
#' for the duration of the session, so a panel that marks absent calls with a
#' character haploStack does not know about — \code{"#"}, \code{"?"}, \code{"*"}
#' — is read as missing rather than as a fifth allele.
#'
#' @section The default vocabulary:
#' Missing is decided \strong{per allele}, and the defaults are:
#'
#' \tabular{lll}{
#'   \strong{mode} \tab \strong{missing characters} \tab \strong{notes} \cr
#'   \code{"nucleotide"} \tab \code{"" "NA" "N" "n" "." "-" "0"} \tab
#'     \strong{both} \code{"0"} (PLINK's code) and \code{"-"} \cr
#'   \code{"selectiontools"} \tab \code{"" "NA" "N" "n" "." "-" "0"} \tab
#'     \code{"0"} since 0.6.0 — a 1/2/3/4 file writes it for an absent call;
#'     \strong{deprecated as a registration mode since 0.9.0}: it reaches only
#'     \code{\link{st_convert}} and the \code{"ST"} preset, which translate the
#'     codes on the way in \cr
#'   \code{"numeric"} \tab \code{"" "NA" "N" "n" "." "-"} \tab
#'     \strong{never} \code{"0"} — that is a real genotype
#' }
#'
#' Because the rule is per allele, a diploid genotype cell written \code{"00"}
#' or \code{"--"} is a missing \emph{genotype}: both copies normalise to
#' \code{NA}. A half-called cell (\code{"A0"}, \code{"A-"}) keeps its observed
#' allele and marks only the absent copy — the observed \code{A} still counts
#' toward frequencies, and the NA-any policy at consumption
#' (\code{\link{count_block_subtypes}}) is what decides whether the cell is
#' usable.
#'
#' @section What is refused:
#' A registered character must be a \strong{single} character, and must not be a
#' real allele of that mode: A/C/G/T for \code{"nucleotide"}, the codes 1-4 for
#' \code{"selectiontools"}, and \strong{any digit} for \code{"numeric"}. A
#' dosage panel that really does write \code{9} for missing should declare it on
#' the read that parses it (\code{extra_na = "9"}), where the caller owns the
#' risk, rather than package-wide.
#'
#' @section Scope and lifetime:
#' The registration is a session option (\code{haploStack.missing_chars}) and is
#' gone on restart — deliberately, because it is a statement about a file being
#' read, not about the package.
#'
#' Since 0.9.0 it is read in exactly one place: where a container is
#' \strong{built or repaired} — the readers, \code{\link{build_geno_meta}},
#' \code{\link{new_geno_meta}} and \code{\link{repair_geno_meta}}. Each of
#' those also takes the same characters as its own \code{missing_chars}
#' argument, which wins over the registry (\code{character(0)} ignores it for
#' that call). The vocabulary a container was normalised against is recorded in
#' \code{metadata$missing_tokens}, and from then on that stamp is what every
#' consumer — \code{\link{is_geno_meta}}, the frequency, IBS and clustering
#' paths — goes by. A character registered \emph{after} a container was built
#' does not change how that container is read; repair it with
#' \code{repair_geno_meta(meta, missing_chars = "#")} if it should. This is
#' what lets two jobs in one process read different panels without one
#' registration leaking into the other's numbers.
#'
#' @param chars Character vector of single characters to register.
#' @param mode  Source format the characters belong to: \code{"nucleotide"},
#'   \code{"selectiontools"} or \code{"numeric"}. \strong{Required} — there is
#'   no package-wide default, because a character that means "absent" in one
#'   vocabulary can be a genotype in another.
#'
#' @return Invisibly, the full effective missing set for \code{mode}.
#'
#' @examples
#' old <- getOption("haploStack.missing_chars")
#' use_missing_chars("#", mode = "nucleotide")
#' missing_chars("nucleotide")
#' missing_chars("numeric")          # unaffected
#' reset_missing_chars()
#' options(haploStack.missing_chars = old)
#'
#' @seealso \code{\link{missing_chars}}, \code{\link{reset_missing_chars}},
#'   \code{\link{repair_geno_meta}}, \code{\link{HS_FILL_NA_DEFAULT}}
#' @export
use_missing_chars <- function(chars, mode) {
  mode  <- .hs_missing_mode(mode, "use_missing_chars")
  chars <- .hs_check_missing_chars(chars, mode, "use_missing_chars")
  if (identical(mode, "selectiontools")) .st_gate("use_missing_chars", per_call = FALSE)  # C-6
  if (identical(mode, "selectiontools"))
    .st_encoding_deprecated("use_missing_chars", what = "mode = \"selectiontools\"",
                            then = paste0("the registration reaches only st_convert() ",
                                          "and the \"ST\" preset, which translate the ",
                                          "codes on the way in"))

  reg <- .hs_option_raw("missing_chars")
  if (!is.list(reg)) reg <- list()
  reg[[mode]] <- unique(c(.extra_missing_chars(mode), chars))
  options(haploStack.missing_chars = reg)

  eff <- .missing_tokens(mode, .extra_missing_chars(mode))
  .hs_inform("use_missing_chars: ", paste0('"', chars, '"', collapse = ", "),
          " now read as missing for '", mode, "' (session only). Effective set: ",
          paste0('"', eff, '"', collapse = " "), ".")
  invisible(eff)
}


#' The effective missing-allele characters for a source format
#'
#' @param chars_only Logical. \code{FALSE} (default) returns the full effective
#'   set — the built-in vocabulary plus anything registered with
#'   \code{\link{use_missing_chars}}. \code{TRUE} returns only what was
#'   registered this session, which is what \code{\link{reset_missing_chars}}
#'   clears.
#' @inheritParams use_missing_chars
#'
#' @return Character vector.
#'
#' @examples
#' missing_chars("nucleotide")
#' missing_chars("numeric")
#'
#' @seealso \code{\link{use_missing_chars}}
#' @export
missing_chars <- function(mode, chars_only = FALSE) {
  mode  <- .hs_missing_mode(mode, "missing_chars")
  extra <- .extra_missing_chars(mode)
  if (isTRUE(chars_only)) extra else .missing_tokens(mode, extra)
}


#' Clear registered missing-allele characters
#'
#' Removes what \code{\link{use_missing_chars}} registered. The built-in
#' vocabulary is not affected and cannot be cleared.
#'
#' @param mode Source format to clear, or \code{NULL} (default) to clear every
#'   mode.
#'
#' @return Invisibly, the registry as it stands after clearing.
#'
#' @examples
#' use_missing_chars("#", mode = "nucleotide")
#' reset_missing_chars("nucleotide")
#' missing_chars("nucleotide", chars_only = TRUE)
#'
#' @seealso \code{\link{use_missing_chars}}
#' @export
reset_missing_chars <- function(mode = NULL) {
  if (is.null(mode)) {
    options(haploStack.missing_chars = NULL)
    return(invisible(list()))
  }
  mode <- .hs_missing_mode(mode, "reset_missing_chars")
  reg  <- .hs_option_raw("missing_chars")
  if (!is.list(reg)) reg <- list()
  reg[[mode]] <- NULL
  options(haploStack.missing_chars = reg)
  invisible(reg)
}
