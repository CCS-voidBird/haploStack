# =============================================================================
# plot_order.R
#
# The sample row-ordering contract shared by every figure with a sample axis.
#
#   use_sample_order()             — broadcast an order to every figure
#   .hs_sample_order()             — read the haploStack.sample_order option
#   .hs_resolve_sample_order()     — the single resolver (argument > option >
#                                    the caller's own default)
#   .hs_order_within()             — apply a resolved order inside one group
#   .hs_order_suffix()             — provenance text for a figure subtitle
#   .hs_tag_sample_order()         — stamp the provenance onto a ggplot
#
# Before 0.5.5 an order could reach exactly one figure: plot_haplotype_alignment()
# took `sample_order` and .resolve_sample_order() (this file's ancestor) validated
# it. Every other sample-axis figure invented its own row order, so two panels of
# the same population put the same line on different rows. This file is that one
# rule, written once.
#
# Precedence, everywhere: an explicit `sample_order` argument wins, then the
# option, then the function's own default — the same shape as `fill_na` (per-call
# argument > haploStack.fill_na > HS_FILL_NA_DEFAULT, R/alphabet.R).
# =============================================================================


#' Broadcast a sample order to every haploStack figure
#'
#' Sets the \code{haploStack.sample_order} option, which every figure with a
#' sample axis reads when its own \code{sample_order} argument is \code{NULL}.
#' The usual argument is a clustering from \code{\link{cluster_samples}}, so one
#' call makes the alignment, the local-GEBV panels, the phenotype bars and the
#' top-block heatmap agree on which row is which.
#'
#' @section Precedence:
#' An explicit \code{sample_order} argument on a plot function wins; then this
#' option; then the function's own default (alphabetical for
#' \code{\link{plot_haplotype_alignment}}, family-history order for
#' \code{\link{plot_localgebv_alignment}}). Pass \code{sample_order = FALSE} to a
#' figure to make it ignore the broadcast for that one call.
#'
#' @section Strict vs best-effort:
#' An order given \emph{as an argument} must be a permutation of the panel's
#' samples — a mismatch is an error, because you named those samples for that
#' figure. A \emph{broadcast} order is best-effort: it is intersected with
#' whatever the panel holds and unmatched panel samples are appended after the
#' ordered ones. A broadcast is a preference for many figures, not an assertion
#' about one, and a filtered panel must not become an error at a distance.
#'
#' @param x An \code{hs_sample_clustering} from \code{\link{cluster_samples}}, a
#'   character vector of sample IDs, or \code{NULL} to clear the broadcast.
#'
#' @return The previous option value, invisibly (so
#'   \code{on.exit(use_sample_order(old))} restores it).
#'
#' @examples
#' \dontrun{
#' cl <- cluster_samples(meta, k = 4)
#' old <- use_sample_order(cl)
#' plot_haplotype_alignment(meta)                 # clustered rows
#' plot_localgebv_alignment(lg, meta)             # same rows, within each family
#' use_sample_order(old)                          # restore
#' }
#'
#' @seealso \code{\link{cluster_samples}}, \code{\link{plot_clustered_alignment}},
#'   \code{\link{plot_haplotype_alignment}}, \code{\link{plot_localgebv_alignment}}
#' @export
use_sample_order <- function(x) {
  old <- .hs_option_raw("sample_order")

  if (is.null(x)) {
    options(haploStack.sample_order = NULL)
    .hs_inform("use_sample_order: sample-order broadcast cleared; figures fall ",
            "back to their own defaults.")
    return(invisible(old))
  }

  # Validate on set, not on use: a bad order must fail at the call the user can
  # see, not three figures later.
  if (inherits(x, "hs_sample_clustering")) {
    ord <- as.character(x$order)
  } else if (is.character(x) || is.factor(x)) {
    ord <- as.character(x)
  } else {
    .hs_abort("use_sample_order: 'x' must be an hs_sample_clustering (from ",
         "cluster_samples()), a character vector of sample IDs, or NULL.")
  }
  if (length(ord) == 0L)
    .hs_abort("use_sample_order: the order is empty.")
  if (anyNA(ord))
    .hs_abort("use_sample_order: the order contains NA sample IDs.")
  if (anyDuplicated(ord))
    .hs_abort("use_sample_order: the order has duplicated sample IDs: ",
         paste(utils::head(unique(ord[duplicated(ord)]), 5L), collapse = ", "),
         ".")

  options(haploStack.sample_order = x)
  desc <- .hs_order_describe(x)
  .hs_inform("use_sample_order: broadcasting a row order over ", length(ord),
          " sample(s)",
          if (is.null(desc)) "." else paste0(" (", desc, ")."))
  invisible(old)
}


# Read the broadcast order. Returns NULL when unset. Validated on read as well
# as in use_sample_order(), because options() can be set directly: an unusable
# value is ignored with a warning (.check_sample_order(), R/state.R).
# @keywords internal
.hs_sample_order <- function() {
  .hs_setting("sample_order")
}


# Human-readable description of where an order came from, for a subtitle or a
# message. NULL for anything unlabelled.
# @keywords internal
.hs_order_describe <- function(x) {
  if (inherits(x, "hs_sample_clustering")) {
    p <- x$params
    return(sprintf("%s / %s%s", p$linkage %||% "hclust",
                   p$distance %||% "distance",
                   if (!is.null(p$k)) sprintf(", k=%d", as.integer(p$k)) else ""))
  }
  NULL
}


# Subtitle/message fragment for a described order; "" when there is nothing to say.
# @keywords internal
.hs_order_suffix <- function(desc, prefix = " \u00B7 rows: ") {
  if (is.null(desc) || !nzchar(desc)) "" else paste0(prefix, desc)
}


#' Resolve a sample row order (internal)
#'
#' The single seam behind every figure's \code{sample_order} argument.
#'
#' @param order_spec \code{NULL} (use the option, else the caller's default),
#'   \code{FALSE} (the caller's default, ignoring the option), \code{TRUE} (the
#'   option; an error when unset), a character vector, or an
#'   \code{hs_sample_clustering}.
#' @param samples Character vector of the sample IDs present in the panel.
#' @param caller Name used in messages, so an error names the function the user
#'   called.
#' @param subset \code{"strict"} — an explicitly supplied order must be a
#'   permutation of \code{samples} (the pre-0.5.5 contract); \code{"intersect"} —
#'   keep what matches and append the rest. A broadcast order is always resolved
#'   as \code{"intersect"} regardless of this argument (see
#'   \code{\link{use_sample_order}}).
#' @param default \code{"sorted"} — with nothing supplied, return
#'   \code{sort(unique(samples))}; \code{"keep"} — return \code{NULL}, meaning
#'   "caller keeps its own order untouched".
#'
#' @return A list: \code{order} (a full permutation of \code{samples}, or
#'   \code{NULL} under \code{default = "keep"} with nothing supplied),
#'   \code{source} (\code{"default"} / \code{"argument"} / \code{"option"}),
#'   \code{clustering} (the object when one was supplied, else \code{NULL}) and
#'   \code{desc} (provenance text, or \code{NULL}).
#' @keywords internal
#' @noRd
.hs_resolve_sample_order <- function(order_spec = NULL, samples,
                                     caller  = "plot",
                                     subset  = c("strict", "intersect"),
                                     default = c("sorted", "keep")) {
  subset  <- match.arg(subset)
  default <- match.arg(default)
  samples <- unique(as.character(samples))

  fallback <- function(src) {
    list(order      = if (default == "sorted") .hs_sort(samples) else NULL,
         source     = src,
         clustering = NULL,
         desc       = NULL)
  }

  # ---- Where does the order come from? --------------------------------------
  source <- "argument"
  if (is.null(order_spec)) {
    order_spec <- .hs_sample_order()
    source     <- if (is.null(order_spec)) "default" else "option"
  } else if (is.logical(order_spec) && length(order_spec) == 1L) {
    if (isTRUE(order_spec)) {
      order_spec <- .hs_sample_order()
      if (is.null(order_spec))
        .hs_abort(caller, ": sample_order = TRUE, but no order is broadcast. ",
             "Call use_sample_order(cluster_samples(meta)) first, or pass an ",
             "order directly.")
      source <- "option"
    } else {
      return(fallback("default"))                 # FALSE = ignore the broadcast
    }
  }
  if (source == "default") return(fallback("default"))

  clustering <- NULL
  if (inherits(order_spec, "hs_sample_clustering")) {
    clustering <- order_spec
    order_spec <- order_spec$order
  }
  desc <- .hs_order_describe(clustering) %||%
    (if (source == "option") "broadcast order" else "custom order")

  ord <- as.character(order_spec)
  if (anyDuplicated(ord))
    .hs_abort(caller, ": 'sample_order' has duplicated sample IDs.")

  # A broadcast is best-effort; an explicit argument is checked (see the
  # function docs for why the two differ).
  mode <- if (source == "option") "intersect" else subset

  missing_ord <- setdiff(samples, ord)     # in the panel, not in the order
  extra_ord   <- setdiff(ord, samples)     # in the order, not in the panel

  if (mode == "strict" && (length(missing_ord) || length(extra_ord)))
    .hs_abort(caller, ": 'sample_order' must be a permutation of the panel's ",
         "samples. ",
         if (length(missing_ord))
           sprintf("Missing %d: %s. ", length(missing_ord),
                   paste(utils::head(missing_ord, 5L), collapse = ", ")) else "",
         if (length(extra_ord))
           sprintf("Unknown %d: %s.", length(extra_ord),
                   paste(utils::head(extra_ord, 5L), collapse = ", ")) else "")

  if (length(extra_ord) && length(extra_ord) == length(ord))
    .hs_warn(caller, ": the supplied sample order matches none of this panel's ",
            "samples; the panel keeps its own order.")

  # One warning per call, never per panel: a nested figure would otherwise emit
  # the same message once per family.
  if (length(missing_ord))
    .hs_inform(caller, ": ", length(missing_ord), " panel sample(s) are not in the ",
            "supplied order and are placed after the ordered ones: ",
            paste(utils::head(missing_ord, 5L), collapse = ", "),
            if (length(missing_ord) > 5L) " ..." else "")

  kept <- ord[ord %in% samples]
  tail <- if (default == "sorted") .hs_sort(missing_ord) else missing_ord
  list(order      = c(kept, tail),
       source     = source,
       clustering = clustering,
       desc       = desc)
}


# Apply a resolved order INSIDE one group, preserving the group's membership.
# `ids` keep their relative order when they are absent from `ord`, and unordered
# ids always follow ordered ones — the rule that lets a clustering re-order the
# rows of a family panel without dissolving the family.
# @keywords internal
.hs_order_within <- function(ids, ord) {
  ids <- as.character(ids)
  if (is.null(ord) || length(ids) == 0L) return(ids)
  hit <- match(ids, ord)
  ids[order(is.na(hit), hit, seq_along(ids))]
}


# Warn when composed panels were drawn with different sample orders. Compared on
# the samples the panels SHARE, so a phenotype panel covering a subset of the
# alignment's samples is not a mismatch by itself — only a different relative
# order is. Panels with no order record (drawn outside the contract, or with the
# default order) are skipped rather than assumed to agree.
# @keywords internal
.hs_warn_order_mismatch <- function(plots) {
  ords <- lapply(plots, function(p) attr(p, "sample_order"))
  keep <- which(vapply(ords, function(o) length(o) > 0L, logical(1L)))
  if (length(keep) < 2L) return(invisible(FALSE))

  ref  <- ords[[keep[1L]]]
  bad  <- integer(0L)
  for (i in keep[-1L]) {
    shared <- intersect(ref, ords[[i]])
    if (length(shared) < 2L) next
    if (!identical(ref[ref %in% shared], ords[[i]][ords[[i]] %in% shared]))
      bad <- c(bad, i)
  }
  if (length(bad) == 0L) return(invisible(FALSE))

  src <- vapply(plots, function(p)
    attr(p, "sample_order_source") %||% "unknown", character(1L))
  .hs_warn("plot_meta_alignment: panel(s) ", paste(bad, collapse = ", "),
          " were drawn in a different sample order from the reference panel ",
          "(sources: ", src[keep[1L]], " vs ",
          paste(unique(src[bad]), collapse = "/"),
          "). The reference order is imposed on every panel, so those panels ",
          "are being re-ordered. Draw them with the same 'sample_order' \u2014 or ",
          "broadcast one with use_sample_order() \u2014 to make the figure honest. ",
          "Pass check_order = FALSE to silence this.")
  invisible(TRUE)
}


# Stamp the resolved order onto a ggplot so a reader (and plot_meta_alignment())
# can tell a clustered panel from an alphabetical one. Set LAST: `p + <layer>`
# drops attributes.
# @keywords internal
# ---- chromosome order (0.9.0, design D3 §3.6, decision-form row 29 a) --------

#' One natural chromosome order
#'
#' The single rule every haploStack figure axis and table row order follows
#' since 0.9.0 (decision-form row 29 a). Up to 0.8.8 four rules disagreed:
#' the circos and Manhattan plots put \code{"Un"} first and \code{"10H"} before
#' \code{"1H"}, the alignment panels used a natural order, and tables sorted
#' chromosomes as plain strings (\code{"10"} before \code{"2"}).
#'
#' The rule: strip a leading \code{chr}, \code{chrom} or \code{chromosome}
#' (any case, with an optional \code{_}, \code{.} or \code{-}) to get the key;
#' names whose key starts with an integer come first, ordered by that integer
#' and then by the rest of the key (\code{1A} before \code{1B}, \code{2A} before
#' \code{10A}); every other name follows, alphabetically by key; ties keep the
#' original spelling order. The comparison is locale-free.
#'
#' Who sees a change from 0.8.8: panels with ten or more numbered chromosomes,
#' bare numbers plus a scaffold (\code{"Un"}), or \code{chr}-prefixed names.
#' Barley \code{1H}–\code{7H} (with or without \code{Un}) is identical under
#' every rule. Values never change, only the order of rows and axes;
#' HapSelect block ids embed the key, so they change on exactly those panels.
#'
#' @param x Chromosome names (any vector; coerced to character, \code{NA}
#'   dropped).
#' @return The distinct names in natural order.
#' @examples
#' hs_chrom_levels(c("10H", "2H", "1H", "Un"))   # 1H 2H 10H Un
#' hs_chrom_levels(c("chr10", "chr2", "chr1"))   # chr1 chr2 chr10
#' hs_chrom_levels(c("X", "2", "1", "Y", "10", "MT"))  # 1 2 10 MT X Y
#' @export
hs_chrom_levels <- function(x) {
  u <- unique(as.character(x))
  u <- u[!is.na(u)]
  if (!length(u)) return(character(0L))
  key  <- sub("^(chromosome|chrom|chr)[_.-]?", "", u, ignore.case = TRUE)
  lead <- suppressWarnings(as.numeric(
    ifelse(grepl("^[0-9]+", key), sub("^([0-9]+).*$", "\\1", key), NA)))
  rest <- ifelse(is.na(lead), key, sub("^[0-9]+", "", key))
  u[order(is.na(lead), lead, rest, u, method = "radix")]
}

# Integer rank of each name under hs_chrom_levels(); for data.table sorts.
.hs_chrom_rank <- function(x) {
  x <- as.character(x)
  match(x, hs_chrom_levels(x))
}

.hs_tag_sample_order <- function(p, res) {
  if (is.null(res)) return(p)
  attr(p, "sample_order")        <- res$order
  attr(p, "sample_order_source") <- res$source
  p
}
