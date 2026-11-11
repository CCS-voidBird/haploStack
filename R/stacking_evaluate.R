# =============================================================================
# stacking_evaluate.R
#
# Stacking evaluation — how many of a line's block haplotypes carry a favourable
# local GEBV:
#   count_favourable_blocks()  — exported: the per-sample evaluation table
#   plot_favourable_blocks()   — exported: the 1/0 dummy alignment heatmap
#   .sample_generation_map()   — internal: family history -> (IID, GID, FID, Label)
#   .favourable_dummy_matrix() — internal: long counts -> wide 1/0 matrix
#
# The metric here is the SIGN of a block's local GEBV, counted per haplotype
# copy. It is deliberately distinct from count_favourable_subtypes()
# (R/block_subtypes.R), which counts identity matches against catalogue-chosen
# subtypes. Both answer "how favourable is this line", by different definitions;
# they join on Sample / IID.
#
# Nothing here reimplements the local-GEBV computation: it is
# compute_profile_localgebv(level = "haplotype"), and the sign classification is
# run_haplotype_catalogue()'s own fcase() idiom (R/catalogue_main.R).
# =============================================================================


# ---------------------------------------------------------------------------
# Internal: family history -> one pedigree row per sample
# ---------------------------------------------------------------------------

# Resolve each sample's (GID, FID, Label) from a trace_family_history() result.
#
# A sample can appear in several generations — a parent in gen_2 is usually an
# offspring in gen_1 — so one membership has to be chosen. Follows
# .sample_family_map()'s convention (R/catalogue_family.R): the OLDEST
# generation in which the sample is an offspring wins (its birth family). A
# sample that is never an offspring takes the oldest generation in which it
# appears as a parent, so founders still carry a family.
#
# Returns data.table(IID, GID, FID, Label), one row per sample seen in the
# history; empty (typed) table for a NULL/empty history.
.sample_generation_map <- function(history) {
  empty <- data.table(IID   = character(0L), GID = character(0L),
                      FID   = character(0L), Label = character(0L))
  if (is.null(history) || length(history) == 0L) return(empty)

  gen_names <- names(history)
  if (is.null(gen_names)) gen_names <- paste0("gen_", seq_along(history))

  parts <- lapply(seq_along(history), function(g) {
    gen <- as.data.table(history[[g]])
    if (nrow(gen) == 0L ||
        !all(c("FID", "IID", "Label") %in% names(gen))) return(NULL)
    data.table(IID   = as.character(gen$IID),
               GID   = gen_names[g],
               FID   = as.character(gen$FID),
               Label = as.character(gen$Label),
               gen_i = g)
  })
  long <- data.table::rbindlist(parts)
  if (nrow(long) == 0L) return(empty)

  # Offspring memberships first, then generation order: the first row per IID is
  # its oldest birth family, falling back to its oldest parent appearance.
  long[, is_off := Label == "O"]
  data.table::setorder(long, IID, -is_off, gen_i)
  out <- unique(long, by = "IID")
  out[, .(IID, GID, FID, Label)]
}


# ---------------------------------------------------------------------------
# Internal: long per-copy counts -> wide 1/0 dummy matrix
# ---------------------------------------------------------------------------

# `long` carries Block, SampleHap and Preferred (logical). `block_order` is the
# genomic block order of the source local-GEBV matrix; `col_order` its column
# order. The result has the same shape as the local-GEBV matrix it came from, so
# plot_localgebv_alignment() accepts it unchanged.
.favourable_dummy_matrix <- function(long, block_order, col_order) {
  d <- data.table::dcast(long, Block ~ SampleHap, value.var = "Preferred",
                         fun.aggregate = function(z) as.integer(any(z)),
                         fill = 0L)
  # Restore the source ordering (dcast sorts both margins).
  d <- d[match(block_order, Block)]
  data.table::setcolorder(d, c("Block", intersect(col_order, names(d))))
  d[]
}


# ---------------------------------------------------------------------------
# Exported: the evaluation table
# ---------------------------------------------------------------------------

#' Count favourable blocks per sample by local-GEBV sign
#'
#' Scores every line in a \code{GenoMeta} by how many of its block haplotypes
#' carry a favourable local GEBV. One row per sample, with the pedigree
#' coordinates (\code{GID}/\code{FID}/\code{Label}) taken from the object's
#' family history when it has one.
#'
#' Counting is always at \strong{haplotype level}: each block contributes
#' \code{ploidy} copies per sample, so
#' \code{nPositive + nNegative + nZero == nBlock * ploidy}. \code{nZero} is not a
#' column — it is that identity's remainder. Genotype-level positive/negative
#' block counts already exist as \code{\link{run_haplotype_catalogue}}'s
#' \code{block_counts}.
#'
#' \strong{Which sign is favourable} is set by \code{prefer}; \code{nPositive}
#' and \code{nNegative} are reported either way, so one call describes both
#' directions and \code{nPreferred} simply names the one being selected for.
#'
#' \strong{Centering matters to the sign.} The default since 0.9.0 is
#' \code{center = "locus_mean"} (decision-form row 15 c): each marker's allele
#' effects are shifted by their frequency-weighted mean, so a positive local GEBV
#' means \dQuote{above the panel mean for this block} while every haplotype of a
#' block moves by the same amount and the block's subtype ranking is untouched.
#' Without centering (\code{"none"}) the sign mostly reflects which allele the
#' effect file called the reference, and counting positive blocks stops being a
#' breeding statement. Up to 0.8.8 the default was \code{"allele_freq"}
#' (\code{center_by_freq = TRUE}); pass \code{center = "allele_freq"} to
#' reproduce those counts. A value within \code{tol} of zero is neither sign.
#'
#' @param geno_meta A valid \code{GenoMeta} with block assignments (see
#'   \code{\link{attach_block_info}}).
#' @param effect_df Marker effects as \code{(Marker, Allele, Effect)}, read by
#'   position. \code{NULL} (default) uses \code{geno_meta$marker_effects} (see
#'   \code{trait}); an explicit table always overrides it.
#' @param trait Optional trait name selecting one table from a stored per-trait
#'   \code{geno_meta$marker_effects} list when \code{effect_df} is \code{NULL}.
#' @param prefer Which sign counts as favourable: \code{"positive"} (default) or
#'   \code{"negative"}. Sets \code{nPreferred} only.
#' @param samples,blocks,chromosomes Optional subsets, passed to
#'   \code{\link{compute_profile_localgebv}}.
#' @param center        Centring mode applied to the allele effects:
#'   \code{"none"}, \code{"locus_mean"} or \code{"allele_freq"}. See
#'   \code{\link{compute_profile_localgebv}} for what each preserves.
#' @param center_by_freq \strong{Deprecated} logical form of \code{center}
#'   (\code{TRUE} = \code{"allele_freq"}, \code{FALSE} = \code{"none"}); warns
#'   on every call. The default centring here is \code{"locus_mean"} since 0.9.0
#'   (see above).
#' @param tol Non-negative number. A local GEBV with \code{|value| <= tol} counts
#'   as \code{Zero} rather than \code{Positive} / \code{Negative}. Default
#'   \code{0}: the strict sign rule of 0.8.8.
#' @param fill_na Value for (block, copy) cells with no effect information.
#'   Default \code{0}, matching \code{\link{compute_profile_localgebv}} — every
#'   block then holds a value and \code{nBlock} is the same for every sample.
#'   Pass \code{NA_real_} to keep \dQuote{no information} distinct from
#'   \dQuote{zero effect}: \code{nBlock} then varies per sample and \code{GEBV}
#'   sums the observed cells only.
#' @param plot Logical; also draw the 1/0 dummy alignment via
#'   \code{\link{plot_favourable_blocks}}. Default \code{TRUE}. Needs a family
#'   history; without one it warns and \code{plots} is \code{NULL}.
#' @param output_file Optional path; writes \code{table} as a tab-separated
#'   file. The list is returned either way.
#'
#' @return A named list, always with the same elements:
#' \describe{
#'   \item{\code{table}}{\code{data.table}, one row per sample:
#'     \code{GID} (generation, \code{NA} without a family history), \code{FID}
#'     (family, \code{NA} likewise), \code{IID} (sample), \code{Label}
#'     (\code{"P"}/\code{"O"}, \code{NA} likewise), \code{nBlock} (blocks
#'     evaluated), \code{GEBV} (total local GEBV over them), \code{nPositive},
#'     \code{nNegative}, \code{nPreferred}. Ordered by \code{GID}, \code{FID},
#'     then descending \code{nPreferred}.}
#'   \item{\code{plots}}{The \code{\link{plot_localgebv_alignment}} nested list
#'     (generation -> family -> \code{ggplot}) for the dummy matrix, or
#'     \code{NULL}.}
#'   \item{\code{dummy}}{\code{data.table} — the wide \code{Block} x copy 1/0
#'     matrix the figure is drawn from, \code{1} marking a preferred-sign copy.}
#'   \item{\code{prefer}}{The preference used.}
#'   \item{\code{level}}{\code{"haplotype"} — the counting unit, always.}
#' }
#'
#' @examples
#' \dontrun{
#' meta <- load_geno_meta("my_geno_meta.rds")
#' meta <- attach_family_history(meta)
#'
#' res <- count_favourable_blocks(meta, trait = "Yield")
#' res$table                      # the ranking table
#' res$plots[["gen_1"]][["F_1"]]  # the 1/0 dummy alignment
#'
#' # Selecting for a low value instead (e.g. days to maturity)
#' count_favourable_blocks(meta, trait = "DTM", prefer = "negative")$table
#' }
#'
#' @seealso \code{\link{compute_profile_localgebv}} (the local GEBVs counted
#'   here), \code{\link{plot_favourable_blocks}},
#'   \code{\link{count_favourable_subtypes}} (the same question by subtype
#'   identity rather than sign; joins on \code{IID} / \code{Sample}),
#'   \code{\link{run_haplotype_catalogue}} (genotype-level \code{block_counts})
#'
#' @export
count_favourable_blocks <- function(geno_meta,
                                    effect_df      = NULL,
                                    trait          = NULL,
                                    prefer         = c("positive", "negative"),
                                    samples        = NULL,
                                    blocks         = NULL,
                                    chromosomes    = NULL,
                                    center         = NULL,
                                    center_by_freq = NULL,
                                    tol            = 0,
                                    fill_na        = 0,
                                    plot           = TRUE,
                                    output_file    = NULL) {

  center <- .resolve_center(center, center_by_freq, "count_favourable_blocks")
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0)
    .hs_abort("count_favourable_blocks: 'tol' must be a single non-negative number.")
  prefer <- match.arg(prefer)
  .hs_check_meta(geno_meta, "count_favourable_blocks", "geno_meta")

  # 1. Haplotype-level local GEBVs (request 1e: never genotype level) ----------
  lg <- compute_profile_localgebv(
    geno_meta, effect_df = effect_df, trait = trait,
    samples = samples, blocks = blocks, chromosomes = chromosomes,
    level = "haplotype", center = center, fill_na = fill_na)
  lg <- as.data.table(lg)

  hap_cols <- setdiff(names(lg), "Block")
  if (nrow(lg) == 0L || length(hap_cols) == 0L)
    .hs_abort("count_favourable_blocks: no block x haplotype local GEBVs to count ",
         "(empty genotype data after filtering, or no markers with block ",
         "assignments).")

  # 2. Long form: one row per (Block, Sample, Haplotype) copy ------------------
  #    Same split as compute_block_catalogue() (R/catalogue_gebv.R).
  long <- data.table::melt(lg, id.vars = "Block", measure.vars = hap_cols,
                           variable.name = "SampleHap", value.name = "GEBV",
                           variable.factor = FALSE)
  hc <- .hap_cols(hap_cols)                            # one parser (sc/F24)
  long[, Sample := hc$base[match(SampleHap, hc$cols)]]
  long[, GEBV   := as.numeric(GEBV)]

  # 3. Sign classification — run_haplotype_catalogue()'s idiom (request 1d) ----
  # C-prefer (0.9.0): a value within `tol` of zero is Zero (default tol = 0 is
  # the 0.8.8 rule: strict sign).
  long[, Sign := data.table::fcase(
    GEBV >  tol, "Positive",
    GEBV < -tol, "Negative",
    default  = "Zero")]                     # NA (no information) falls here too
  long[is.na(GEBV), Sign := NA_character_]
  long[, Preferred := Sign == if (prefer == "positive") "Positive" else "Negative"]
  long[is.na(Preferred), Preferred := FALSE]

  # 4. Per-sample counts ------------------------------------------------------
  #    nBlock counts blocks that actually hold a value, so fill_na = NA_real_
  #    makes it per-sample meaningful; with the default 0 it is every block.
  tab <- long[, .(nBlock    = data.table::uniqueN(Block[!is.na(GEBV)]),
                  GEBV      = sum(GEBV, na.rm = TRUE),
                  nPositive = sum(Sign == "Positive", na.rm = TRUE),
                  nNegative = sum(Sign == "Negative", na.rm = TRUE),
                  nPreferred = sum(Preferred)),
              by = .(IID = Sample)]

  # 5. Pedigree coordinates ---------------------------------------------------
  ped <- .sample_generation_map(geno_meta$metadata$family_history)
  if (nrow(ped) > 0L) {
    tab[ped, on = "IID", `:=`(GID = i.GID, FID = i.FID, Label = i.Label)]
  } else {
    tab[, `:=`(GID = NA_character_, FID = NA_character_,
               Label = NA_character_)]
  }

  data.table::setcolorder(tab, c("GID", "FID", "IID", "Label", "nBlock",
                                 "GEBV", "nPositive", "nNegative",
                                 "nPreferred"))
  data.table::setorder(tab, GID, FID, -nPreferred, IID, na.last = TRUE)

  # 6. Dummy 1/0 matrix + figure ----------------------------------------------
  dummy <- .favourable_dummy_matrix(long, lg$Block, hap_cols)

  plots <- NULL
  if (isTRUE(plot)) {
    plots <- plot_favourable_blocks(dummy, family = geno_meta, prefer = prefer)
  }

  if (!is.null(output_file)) {
    data.table::fwrite(tab, output_file, sep = "\t")
    .hs_inform("Favourable-block counts written to: ", output_file)
  }

  structure(hs_center = center,                     # D3.5: the centre travels with the result
  list(table  = tab[],
       plots  = plots,
       dummy  = dummy,
       prefer = prefer,
       level  = "haplotype"))
}


# ---------------------------------------------------------------------------
# Exported: the 1/0 dummy alignment
# ---------------------------------------------------------------------------

#' Dummy alignment heatmap of favourable block haplotypes
#'
#' The viz half of \code{\link{count_favourable_blocks}}: the same per-copy
#' judgement drawn as an alignment, with a preferred-sign local GEBV marked
#' \code{1} and everything else \code{0}.
#'
#' This is \code{\link{plot_localgebv_alignment}} called on the 1/0 dummy matrix
#' rather than a second heatmap implementation, so the chromosome separators,
#' per-family panels, offspring/founder boundary and genomic block ordering are
#' the same figure machinery as the continuous local-GEBV alignment. The fill
#' scale is clipped to \code{c(0, 1)}, which puts \code{0} on the diverging
#' scale's white midpoint and \code{1} at full colour; the scale's negative half
#' is unreachable for a dummy, so including it would leave half the legend dead.
#'
#' @param x Either a \code{\link{count_favourable_blocks}} result (its
#'   \code{dummy} matrix is used) or a wide dummy table directly: first column
#'   \code{Block}, then one 0/1 column per haplotype copy.
#' @param family A \code{GenoMeta} carrying \code{metadata$family_history}, or a
#'   \code{\link{trace_family_history}} list. Required — counting works without a
#'   pedigree but the per-family figure does not; a missing history warns and
#'   returns \code{NULL} rather than erroring, so \code{count_favourable_blocks()}
#'   still returns its table.
#' @param prefer The preference the dummy encodes, used in the legend title
#'   only. Taken from \code{x} when \code{x} is a result list.
#' @param sample_order Row order within each family panel, forwarded to
#'   \code{\link{plot_localgebv_alignment}}: a character vector, an
#'   \code{hs_sample_clustering}, \code{TRUE}/\code{FALSE} for the broadcast
#'   order (\code{\link{use_sample_order}}), or \code{NULL} (default). Named
#'   explicitly rather than left to \code{...} so it is documented where it is
#'   used.
#' @param ... Passed to \code{\link{plot_localgebv_alignment}} (e.g.
#'   \code{selected_blocks}, \code{show_block_labels}).
#'
#' @return The \code{\link{plot_localgebv_alignment}} nested list (generation ->
#'   family -> \code{ggplot}), or \code{NULL} when no family history is
#'   available.
#'
#' @examples
#' \dontrun{
#' res <- count_favourable_blocks(meta, trait = "Yield", plot = FALSE)
#' plots <- plot_favourable_blocks(res, family = meta)
#' plots[["gen_1"]][["F_1"]]
#' }
#'
#' @seealso \code{\link{count_favourable_blocks}},
#'   \code{\link{plot_localgebv_alignment}}
#'
#' @export
plot_favourable_blocks <- function(x, family = NULL,
                                   prefer = c("positive", "negative"),
                                   sample_order = NULL, ...) {

  # Result list vs bare dummy table.
  if (is.list(x) && !is.data.frame(x) && !is.null(x$dummy)) {
    if (missing(prefer) && !is.null(x$prefer)) prefer <- x$prefer
    dummy <- x$dummy
  } else {
    dummy <- x
  }
  prefer <- match.arg(prefer)

  if (!is.data.frame(dummy) && !inherits(dummy, "data.table"))
    .hs_abort("plot_favourable_blocks: 'x' must be a count_favourable_blocks() ",
         "result or a wide dummy table with 'Block' as its first column.")

  if (is.null(family))
    .hs_abort("plot_favourable_blocks: 'family' is required \u2014 pass the GenoMeta the ",
         "counts were computed from, or a trace_family_history() list.")

  has_history <- if (is_geno_meta(family))
    !is.null(family$metadata$family_history) else TRUE
  if (!has_history) {
    .hs_warn("plot_favourable_blocks: the GenoMeta carries no family history, ",
            "so the per-family dummy alignment cannot be drawn. Run ",
            "attach_family_history() first; returning NULL.")
    return(NULL)
  }

  plot_localgebv_alignment(
    dummy, family,
    sample_order = sample_order,
    fill_limits  = c(0, 1),
    fill_name    = sprintf("Favourable\n(%s)", prefer),
    ...)
}
