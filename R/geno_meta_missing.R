# =============================================================================
# geno_meta_missing.R — the missing-allele portal (MAN, 0.6.0; settings 0.9.0)
#
# Inside a GenoMeta a missing allele is NA and nothing else. Every route into a
# container passes .normalise_at_portal(): the source format's own missing
# tokens (and the extra characters a call declares, or the session registry
# at build time) become NA, and the container is stamped with the vocabulary
# it was built under (metadata$missing_tokens). repair_geno_meta() runs the
# portal over a container that still carries tokens; the token report and
# its notice live here too. Split out of format_helper.R in 0.9.0 (D2-9).
# =============================================================================

#' Normalise a source format's missing tokens to NA (internal)
#'
#' Inside a \code{GenoMeta} a missing allele is \code{NA} and nothing else. This
#' is the primitive that enforces the invariant: it maps every missing token of
#' the format being \emph{parsed} to \code{NA_character_} in an allele array.
#'
#' @section Why the argument is a source format:
#' The token set is keyed on the format the caller is parsing, supplied by the
#' caller — \strong{not} read back out of a container's \code{metadata$encoding}.
#' \code{metadata$encoding} is a record of provenance, not an instruction
#' consumers follow to decide what an allele means. Normalisation therefore
#' happens once, at the boundary where the raw tokens and the format that
#' explains them are both in hand; after that no layer needs to be
#' encoding-aware to know what is missing.
#'
#' @section The dosage rule:
#' \code{"0"} is a missing code for \code{"nucleotide"} and
#' \code{"selectiontools"} sources, and is \strong{never} missing for a
#' dosage/\code{"numeric"} source, where \code{0} is a real genotype
#' (homozygous reference). The rule is enforced here even against an explicit
#' \code{extra_na}, because corrupting a dosage panel would be a worse failure
#' than the bug this function exists to fix.
#'
#' @section Per allele, and what that means for a genotype:
#' Normalisation acts on the \strong{allele} array, after the reader has split
#' each genotype cell into its copies. For a nucleotide source both \code{"0"}
#' and \code{"-"} are missing, so a diploid cell written \code{"00"} or
#' \code{"--"} arrives here as two missing alleles and becomes a missing
#' genotype by consequence. A half-called cell (\code{"A0"}, \code{"A-"}) keeps
#' its observed allele: only the absent copy becomes \code{NA}, the \code{A}
#' still counts toward allele frequencies, and whether the cell is usable is
#' decided later by the NA-any policy at consumption
#' (\code{\link{count_block_subtypes}}).
#'
#' Characters registered with \code{\link{use_missing_chars}} join the set for
#' their mode and are honoured here like any built-in one, unless the call
#' names its own \code{missing_chars}.
#'
#' @param arr           Allele array (any dimension) or matrix. Non-character
#'   input is returned unchanged — a numeric dosage array carries no tokens.
#' @param source_format Allele vocabulary of the source being parsed:
#'   \code{"nucleotide"}, \code{"selectiontools"}, or \code{"numeric"}. Any
#'   other value falls back to the encoding-independent token set
#'   (\code{""}, \code{"NA"}, \code{"N"}, \code{"n"}, \code{"."}, \code{"-"}).
#' @param extra_na      Optional character vector of additional sentinels the
#'   caller knows its own format uses (e.g. \code{read_plink_ped()} declares
#'   \code{"0"}).
#' @param missing_chars Extra missing-allele characters for this call:
#'   \code{NULL} (default) takes the session registry
#'   (\code{\link{use_missing_chars}}), a character vector names the characters
#'   for \code{source_format} (validated as the registry is; \code{character(0)}
#'   ignores the registry), a per-mode list is read like the registry.
#'
#' @return Named list: \code{array} (the normalised array, dimensions and
#'   dimnames preserved) and \code{n_changed} (integer count of cells set to
#'   \code{NA}). The count is returned rather than reported so callers can
#'   state what they changed instead of mutating silently.
#'
#' @seealso \code{\link{repair_geno_meta}}, \code{\link{is_geno_meta}}
#' @keywords internal
normalise_missing_alleles <- function(arr, source_format, extra_na = NULL,
                                      missing_chars = NULL) {
  n_zero <- list(array = arr, n_changed = 0L)
  if (is.null(arr) || length(arr) == 0L) return(n_zero)
  # A numeric/integer array holds no tokens to map — nothing to do, and nothing
  # to risk.
  if (!is.character(arr)) return(n_zero)

  is_dosage <- identical(source_format, "numeric")
  tokens    <- .missing_tokens(source_format,
                               .extra_missing_chars(source_format, missing_chars))

  if (!is.null(extra_na)) {
    extra_na <- as.character(extra_na)
    extra_na <- extra_na[!is.na(extra_na)]
    if (is_dosage && "0" %in% extra_na)
      .hs_warn("normalise_missing_alleles: extra_na = \"0\" ignored for a ",
              "dosage ('numeric') source \u2014 0 is a real genotype there.")
    tokens <- unique(c(tokens, extra_na))
  }
  # Load-bearing: the token set never carries "0" for a dosage source.
  if (is_dosage) tokens <- setdiff(tokens, "0")

  hit <- !is.na(arr) & trimws(arr) %in% tokens
  n   <- sum(hit)
  if (n > 0L) arr[hit] <- NA_character_

  list(array = arr, n_changed = as.integer(n))
}


# Internal: the single normalisation portal, called from new_geno_meta() and
# from nowhere else.
#
# Why one portal and why this one: read_plink_ped(), read_dart_genotype() and
# read_vcf_genotype() bypass build_geno_meta() and construct directly, so
# build_geno_meta() is not a chokepoint — new_geno_meta() is the only function
# every route to a container passes through, including the internal transforms
# (filter_geno_meta(), .replace_block_in_meta(), recode_geno_meta(), the
# simulation and crossing paths). "100 % of source-file -> GenoMeta passes a
# normalisation" is a property of the call graph, not of a convention.
#
# The vocabulary is keyed on `source_format` when the caller declares one (a
# reader knows the format it parsed) and otherwise on the encoding the
# constructor just resolved. Both are the caller's declaration at construction
# time, which is what Q21.d asks for; what Q21.d forbids is a *consumer* reading
# metadata$encoding later to decide what an allele means.
#
# The stamp makes the portal cheap. A container whose metadata$missing_tokens
# already covers the vocabulary in force has been normalised against it, so a
# derived container (a filter, a block swap, a crossing generation) skips the
# scan instead of re-walking a 20-million-cell array to change nothing.
.normalise_at_portal <- function(arr, metadata, source_format = NULL,
                                 extra_na = NULL, missing_chars = NULL,
                                 caller = "new_geno_meta") {
  src <- source_format %||% metadata$encoding

  extra_na <- if (is.null(extra_na)) character(0L)
              else as.character(extra_na)[!is.na(as.character(extra_na))]

  # The session registry (use_missing_chars()) is read here and in
  # repair_geno_meta() only (0.9.0, form row 2): the call's own `missing_chars`
  # wins, and the stamp written below is what every later consumer goes by.
  registered <- .extra_missing_chars(src, missing_chars)
  effective  <- unique(c(.missing_tokens(src, registered), extra_na))
  # The dosage rule, restated at the portal: 0 is a genotype, never an absence.
  if (identical(src, "numeric") || identical(src, "allele"))
    effective <- setdiff(effective, "0")

  stamp <- metadata$missing_tokens
  stamp <- if (is.character(stamp)) unique(stamp[!is.na(stamp)]) else NULL

  # Already normalised against at least this vocabulary — nothing to do.
  if (!is.null(stamp) && length(setdiff(effective, stamp)) == 0L) {
    metadata$missing_tokens <- stamp
    return(list(array = arr, metadata = metadata, n_changed = 0L))
  }

  res <- normalise_missing_alleles(arr, src, extra_na = c(extra_na, registered),
                                   missing_chars = character(0L))
  if (res$n_changed > 0L)
    .hs_inform(caller, ": normalised ", res$n_changed,
            " missing allele cell(s) to NA (source format '",
            src %||% "unknown", "').")

  metadata$missing_tokens <- unique(c(stamp, effective))
  list(array = res$array, metadata = metadata, n_changed = res$n_changed)
}


# Internal: report the missing tokens a container still carries as alleles.
#
# Legacy compatibility, not the mechanism (Q21.d): a container built by 0.6.0 or
# later is normalised at the portal, so this returns NULL. It exists so a
# hand-built container, or one restored from a pre-0.6.0 .rds, is detected
# rather than silently mis-counted.
#
# The token set is keyed on the recorded provenance (metadata$encoding), which
# for this one purpose is the best available statement of what the source file
# meant. A dosage container never has "0" in its set — see .missing_tokens().
#
# It is the UNION of that set and metadata$missing_tokens, the vocabulary the
# container records having been normalised against (written wherever
# normalisation runs). A detector should be more sensitive, not less: the union
# catches both a token this session has registered and one a past session did.
# The dosage exclusion survives the union.
#
# Returns NULL when clean, else list(tokens, n_cells, n_markers, markers).
.missing_token_report <- function(x, source_format = NULL) {
  arr <- x$geno_arrays
  if (is.null(arr) || !is.character(arr) || length(arr) == 0L) return(NULL)

  enc    <- source_format %||% x$metadata$encoding
  tokens <- .missing_tokens(enc)
  stamp  <- x$metadata$missing_tokens
  if (is.character(stamp)) tokens <- unique(c(tokens, stamp[!is.na(stamp)]))
  if (identical(enc, "numeric") || identical(enc, "allele"))
    tokens <- setdiff(tokens, "0")

  # Cheap first pass: the allele vocabulary is tiny next to the array.
  present <- unique(as.vector(arr))
  hit_tok <- intersect(present[!is.na(present)], tokens)
  if (length(hit_tok) == 0L) return(NULL)

  hit      <- !is.na(arr) & arr %in% hit_tok
  dim(hit) <- dim(arr)
  by_mk    <- apply(hit, 2L, any)
  markers  <- dimnames(arr)[[2L]]
  if (is.null(markers)) markers <- as.character(seq_along(by_mk))

  list(tokens    = hit_tok,
       n_cells   = as.integer(sum(hit)),
       n_markers = as.integer(sum(by_mk)),
       markers   = markers[by_mk])
}


# Internal: the FR4 notice. Warns by default, errors when `strict_missing` is
# TRUE (the call's argument, else options(haploStack.strict_missing), resolved
# by .hs_setting()). The warning carries the class "haploStack_missing_token"
# so a caller that reports the condition its own way (load_geno_meta) can
# muffle it without muffling anything else.
.notify_missing_tokens <- function(rep, caller, strict_missing = NULL) {
  if (is.null(rep)) return(invisible(FALSE))
  msg <- paste0(
    caller, ": ", rep$n_cells, " allele cell(s) across ", rep$n_markers,
    " marker(s) hold a missing token (",
    paste0("\"", rep$tokens, "\"", collapse = ", "),
    ") instead of NA. Inside a GenoMeta a missing allele is NA and nothing ",
    "else; a container that predates 0.6.0 (or was built by hand) can still ",
    "carry the source token, which shifts marker effects and frequency-centred ",
    "GEBVs. Repair it with repair_geno_meta().")
  if (isTRUE(.hs_setting("strict_missing", strict_missing)))
    .hs_abort(msg, call. = FALSE)
  .hs_warn(structure(
    class   = c("haploStack_missing_token", "warning", "condition"),
    list(message = paste0(msg, "\n"), call = NULL)))
  invisible(TRUE)
}


#' Normalise the missing alleles of an existing GenoMeta
#'
#' The migration path for containers built before 0.6.0, when a missing allele
#' could be stored as the source file's token (PLINK's \code{"0"}, \code{"N"},
#' \code{"."}, \code{"-"}) rather than as \code{NA}. Maps every such token to
#' \code{NA_character_} and returns the repaired container.
#'
#' Containers built by 0.6.0 or later are normalised at construction, so
#' repairing one is a no-op and says so.
#'
#' @section What it does not do:
#' It does not touch the file a container was loaded from —
#' \code{\link{load_geno_meta}} never rewrites what you saved. Persist the
#' repair yourself with \code{\link{save_geno_meta}}. It is also not imputation:
#' a missing allele is marked missing, not filled.
#'
#' @section The vocabulary it repairs against:
#' The union of the built-in set for the source format (see
#' \code{\link{use_missing_chars}} for the table), anything registered in this
#' session, and \code{metadata$missing_tokens} if the container carries one — so
#' a character registered in the session that \emph{built} a container is still
#' honoured when it is repaired somewhere else. The dosage exclusion survives
#' the union: \code{"0"} is never removed from a \code{"numeric"} container.
#'
#' On return, \code{metadata$missing_tokens} records the vocabulary the
#' container has been cleaned against, which is what keeps it self-describing in
#' a session where nothing is registered.
#'
#' @param meta          A valid \code{GenoMeta} list.
#' @param verbose       Logical. When \code{TRUE} (default), report how many
#'   cells changed.
#' @param source_format Allele vocabulary of the container's source
#'   (\code{"nucleotide"}, \code{"selectiontools"}, \code{"numeric"}). Defaults
#'   to the recorded provenance \code{metadata$encoding}. Supply it explicitly
#'   when that record is absent or wrong — \code{"0"} is a missing code for
#'   nucleotide and SelectionTools sources and a \strong{real genotype} for a
#'   dosage one.
#' @param missing_chars Extra missing-allele characters for this repair:
#'   \code{NULL} (default) takes the session registry
#'   (\code{\link{use_missing_chars}}); a character vector names them for
#'   \code{source_format} directly (\code{character(0)} ignores the registry).
#'   This and the build portal are the only two places the registry is read.
#'
#' @return The \code{GenoMeta} with \code{geno_arrays} normalised.
#'
#' @examples
#' \dontrun{
#' meta <- load_geno_meta("panel.rds")   # warns if it predates 0.6.0
#' meta <- repair_geno_meta(meta)
#' save_geno_meta(meta, "panel.rds")     # your call, not the loader's
#' }
#'
#' @seealso \code{\link{is_geno_meta}}, \code{\link{load_geno_meta}}
#' @export
repair_geno_meta <- function(meta, verbose = TRUE, source_format = NULL,
                             missing_chars = NULL) {
  if (!is.list(meta) || is.null(meta$geno_arrays))
    .hs_abort("repair_geno_meta: 'meta' is not a GenoMeta (no 'geno_arrays').")

  enc <- source_format %||% meta$metadata$encoding
  if (is.null(enc) && isTRUE(verbose))
    .hs_inform("repair_geno_meta: no recorded encoding \u2014 using the ",
            "encoding-independent token set only (\"0\" is left alone).")

  # Repair against the union of this call's vocabulary (the registry, or the
  # `missing_chars` argument) and whatever the container records having been
  # normalised against, so a rule registered in the session that built it is
  # not lost when it is repaired in another.
  stamp <- meta$metadata$missing_tokens
  stamp <- if (is.character(stamp)) unique(stamp[!is.na(stamp)]) else character(0L)
  if (identical(enc, "numeric") || identical(enc, "allele"))
    stamp <- setdiff(stamp, "0")
  registered <- .extra_missing_chars(enc, missing_chars)

  res <- normalise_missing_alleles(meta$geno_arrays, enc,
                                   extra_na = c(stamp, registered),
                                   missing_chars = character(0L))
  meta$geno_arrays <- res$array

  # Stamp the vocabulary this container has now been cleaned against. Written
  # here because this is where normalisation actually runs; the constructor
  # stamps it too when it normalises.
  meta$metadata$missing_tokens <- unique(c(.missing_tokens(enc, registered), stamp))

  if (isTRUE(verbose)) {
    if (res$n_changed == 0L)
      .hs_inform("repair_geno_meta: already normalised \u2014 no cells changed.")
    else
      .hs_inform("repair_geno_meta: normalised ", res$n_changed,
              " missing allele cell(s) to NA.")
  }
  meta
}



# ---------------------------------------------------------------------------
# GenoMeta structure
# ---------------------------------------------------------------------------
