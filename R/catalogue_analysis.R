# =============================================================================
# catalogue_analysis.R
#
# Cross-trait block-catalogue analysis (PRD: dev/PRD_block_catalogue_analysis.md,
# revised 2026-07-07).
#
#   analyze_block_catalogue()  — driver: assemble per-trait block SUBTYPES,
#                                optionally rescale each trait's effects, select
#                                ONE representative subtype per block (shared
#                                across traits) by the highest summed value
#                                across traits, and emit each trait's own raw (S)
#                                / normalized (N) values plus per-trait top
#                                blocks.
#
# RAW BY DEFAULT (0.5.4). `normalize` defaults to "none" and `value_source` to
# "raw", so the effects that reach the cross-trait sum, the reported tables and
# the plots are the local GEBVs as computed — nothing is rescaled unless the
# caller asks. Asking for value_source = "normalized" counts as asking and
# promotes the "zscore" transform; combining it with an explicit
# normalize = "none" is a contradiction and errors.
#
# Sequencing is load-bearing (PRD R1): rescale WITHIN each trait (when asked) ->
# sum the per-trait value across (signed) traits to select one subtype per block
# -> plot each trait's own value. The cross-trait sum is a SELECTION criterion
# only; no composite fitness, corrected effect, or optimizing-potential is
# reported (those cross-trait adjustments were removed in the 2026-07-07
# revision).
# =============================================================================


# ---------------------------------------------------------------------------
# Internal: one-subtype-per-block selection (Step 1)
# ---------------------------------------------------------------------------

# Choose a single representative subtype per block, shared across every trait, as
#   selected(b) = argmax_subtype  Σ_t  w_t · value(subtype, t)
# where `value` is the per-(subtype, trait) selection value carried on the long
# `subtypes` table in column `Vbase` (raw effect S or within-trait standard score
# N, per value_source), and `w` is the signed per-trait weight vector.
#
# `fitness_fn` only filters the candidate subtypes / redefines the per-subtype
# value BEFORE the sum (it no longer does any cross-trait combine):
#   max_effect (default) — candidates = all; value = Vbase.
#   frequent_best        — candidates restricted to Count >= min_count (per block,
#                          falling back to all subtypes when none qualify).
#   expected             — value = Vbase * Count / block-total Count (the
#                          frequency-weighted contribution).
#   min_effect           — as max_effect but the block picks the LOWEST summed
#                          value (argmin).
#   a function(subtypes_dt) -> one row — full generality; receives a block's long
#                          subtype rows and must return exactly one, whose
#                          `Subtype` id represents the block.
#
# Returns a data.table(Block, Subtype) — one row per block.
.select_block_subtype <- function(subtypes, fitness_fn, min_count, w) {
  if (is.function(fitness_fn)) {
    return(subtypes[, {
      row <- as.data.table(fitness_fn(.SD))
      if (nrow(row) != 1L)
        stop("analyze_block_catalogue: a custom 'fitness_fn' must return exactly ",
             "one subtype row (got ", nrow(row), ").")
      .(Subtype = as.character(row$Subtype[1L]))
    }, by = Block])
  }

  st <- data.table::copy(subtypes)
  st[, Vsel := Vbase]

  if (fitness_fn == "expected") {
    tot <- unique(st[, .(Block, Subtype, Count)])[
      , .(BlockTot = sum(Count, na.rm = TRUE)), by = Block]
    st  <- merge(st, tot, by = "Block")
    st[, Vsel := Vbase * Count / BlockTot]
  } else if (fitness_fn == "frequent_best") {
    st[, keep := Count >= min_count]
    st[, keep := keep | !any(keep), by = Block]      # fallback: all if none qualify
    st  <- st[keep == TRUE]
  } else if (!fitness_fn %in% c("max_effect", "min_effect")) {
    stop("analyze_block_catalogue: unknown fitness_fn '", fitness_fn, "'.")
  }

  # Cross-trait summed value per (Block, Subtype), then one subtype per block.
  score <- st[, .(score = sum(w[Trait] * Vsel, na.rm = TRUE)),
              by = .(Block, Subtype)]
  pick  <- if (fitness_fn == "min_effect")
    score[, .SD[which.min(score)], by = Block]
  else
    score[, .SD[which.max(score)], by = Block]
  pick[, .(Block, Subtype)]
}


# ---------------------------------------------------------------------------
# Internal: per-trait normaliser (Step 2 — normalize)
# ---------------------------------------------------------------------------

# Fit a per-trait normaliser on the selected-effect vector `x` (one value per
# block). Returns list(values = normalised x, apply = function reusable on other
# values from the same trait or NULL when the transform cannot be reapplied,
# label = method name). `apply` lets a reuse step normalise other values from the
# same trait on the identical scale.
#
# `zscore` is a TRUE standard score `(x - mean)/sd` (D1, 2026-07-07): it CENTRES
# and scales within the trait, so a normalised value of 0 means an AVERAGE-effect
# subtype for that trait (the diverging genome/circos fill's white midpoint then
# reads "average" under a normalized value, versus "zero effect" under a raw one).
# `minmax` and `rank` remain ZERO-ANCHORED (origin and sign preserved), keeping
# their spread comparable to the white-at-zero fills.
.fit_normalizer <- function(x, method) {
  if (is.function(method))
    return(list(values = as.numeric(method(x)), apply = NULL, label = "custom"))

  switch(method,
    none = list(values = x, apply = function(v) v, label = "none"),
    zscore = {
      # True standard score: centre by the mean AND scale to unit SD, so the
      # mean maps to 0 (average effect) and the spread to unit SD.
      m <- mean(x, na.rm = TRUE)
      s <- stats::sd(x, na.rm = TRUE)
      f <- if (is.na(s) || s == 0) function(v) v * 0
           else function(v) (v - m) / s
      list(values = f(x), apply = f, label = "zscore")
    },
    minmax = {
      # Symmetric scale by the largest absolute effect: maps to [-1, 1] with
      # 0 fixed at 0, rather than re-basing on the (non-zero) minimum.
      mx <- max(abs(x), na.rm = TRUE)
      f  <- if (!is.finite(mx) || mx == 0) function(v) v * 0
            else function(v) v / mx
      list(values = f(x), apply = f, label = "minmax")
    },
    rank = {
      # Sign-preserving magnitude rank: rank by |x| within each sign, scaled to
      # (0, 1] and re-signed, so zero-effect blocks sit at 0 and the two tails
      # spread either side of it.
      n <- sum(!is.na(x))
      r <- data.table::frank(abs(x), ties.method = "average", na.last = "keep")
      v <- if (n > 1L) sign(x) * (r / n) else x * 0
      list(values = as.numeric(v), apply = NULL, label = "rank")
    },
    stop("analyze_block_catalogue: unknown normalize '", method, "'.")
  )
}


# ---------------------------------------------------------------------------
# Internal: resolve per-trait subtype tables (FR1 data-source resolution)
# ---------------------------------------------------------------------------

# Returns list(tabs, gebv): `tabs` is a named (by trait) list of subtype
# data.tables (compute_block_catalogue output); `gebv` is a named list of wide
# Block x Sample genotype-level local-GEBV matrices (for the top-block heatmap).
# Accepts a pre-computed `catalogue` (run_haplotype_catalogue result) or computes
# from geno_meta + effect_list with signature-based subtypes.
.resolve_subtype_tables <- function(geno_meta, effect_list, catalogue,
                                    trait_names, center) {
  if (!is.null(catalogue)) {
    if (is.null(trait_names)) trait_names <- names(catalogue)
    tabs <- lapply(trait_names, function(tr) {
      el <- catalogue[[tr]]
      if (is.null(el) || is.null(el$hap_catalogue))
        stop("analyze_block_catalogue: catalogue['", tr,
             "'] has no 'hap_catalogue'.")
      as.data.table(el$hap_catalogue)
    })
    gebv <- lapply(trait_names, function(tr) {
      gm <- catalogue[[tr]]$gebv_matrix
      if (is.null(gm)) NULL else as.data.table(gm)
    })
    names(tabs) <- names(gebv) <- trait_names
    return(list(tabs = tabs, gebv = gebv))
  }

  if (is.null(geno_meta) || is.null(effect_list))
    stop("analyze_block_catalogue: supply either 'catalogue', or 'geno_meta' + ",
         "'effect_list'.")
  if (is.data.frame(effect_list)) effect_list <- list(Trait1 = effect_list)
  if (is.null(names(effect_list)) || any(names(effect_list) == ""))
    stop("analyze_block_catalogue: 'effect_list' must be a named list ",
         "(one element per trait).")
  if (is.null(trait_names)) trait_names <- names(effect_list)

  bloc <- build_block_location_info(geno_meta)$block_info

  # The long-format allele table and the per-copy block signatures are
  # trait-invariant, so build them ONCE and share across every trait instead of
  # rebuilding them inside each compute_block_catalogue() / compute_profile_localgebv()
  # call (the catalogue-perf core; see R/catalogue_engine.R).
  core <- .catalogue_core(geno_meta)
  sigs <- .compute_block_signatures(geno_meta)

  tabs <- lapply(trait_names, function(tr) {
    eff <- effect_list[[tr]]
    if (is.null(eff))
      stop("analyze_block_catalogue: no effects for trait '", tr, "'.")
    compute_block_catalogue(geno_meta, effect_df = eff, block_info = bloc,
                            center   = center,
                            group_by = "signature",
                            .core = core, .signatures = sigs)
  })
  gebv <- lapply(trait_names, function(tr) {
    compute_profile_localgebv(geno_meta, effect_df = effect_list[[tr]],
                              level = "genotype", center = center,
                              .core = core)
  })
  names(tabs) <- names(gebv) <- trait_names
  list(tabs = tabs, gebv = gebv)
}


# ---------------------------------------------------------------------------
# Exported driver
# ---------------------------------------------------------------------------

#' Cross-trait block-catalogue analysis
#'
#' Stacks the per-trait haplotype-block catalogues and represents each block by
#' the \strong{single subtype whose value, summed across traits, is highest} — one
#' allele combination shared across every trait. Each trait's own value of that
#' selected subtype is then reported so cross-trait trade-offs stay visible. It
#' generalises the origin \code{Block_catalogue_analysis.R} script into a tested,
#' reusable object.
#'
#' @section Raw effects by default:
#' The cross-trait sum, the reported tables and the plots all use the
#' \strong{local GEBVs as computed} — \code{normalize = "none"} and
#' \code{value_source = "raw"}. Nothing is rescaled unless you ask, so a block
#' effect in the output is in the trait's own units and comparable with
#' \code{\link{compute_block_catalogue}}'s \code{Effect} column.
#'
#' Ask for rescaling when the traits sit on scales different enough that the
#' summed selection criterion would otherwise be decided by the largest-scale
#' trait alone: name a \code{normalize} transform, and/or set
#' \code{value_source = "normalized"}. Setting \code{value_source =
#' "normalized"} counts as asking — it promotes \code{normalize = "zscore"} (with
#' a message) rather than returning raw numbers under a normalized label.
#' Combining it with an explicit \code{normalize = "none"} is a contradiction and
#' errors.
#'
#' Under the default the returned \code{N} table is the identity transform of
#' \code{S}, i.e. the same raw values; \code{$normalize} on the result records
#' which transform was actually applied.
#'
#' @section The unit of analysis:
#' A \emph{subtype} is a distinct allele-combination state of a block (one row of
#' \code{\link{compute_block_catalogue}(group_by = "signature")}), carrying its
#' local-GEBV \code{Effect}, a representative carrier, and the population
#' \code{Count}. A subtype's identity is trait-invariant, so the same subtype
#' carries a different \code{Effect} in each trait. Each block is represented by a
#' single selected subtype — a block's value is the value of that one allele
#' combination, not a population average.
#'
#' @section Method (load-bearing sequence):
#' \enumerate{
#'   \item \strong{Rescale within each trait}, if asked: transform every subtype's
#'     effect \code{S(subtype,t)} with \code{normalize} to \code{N(subtype,t)}
#'     (\code{"zscore"} = a true standard score). Under the default
#'     \co\code{"none"} this is the identity, so \code{N == S} and the analysis runs
#'     entirely on raw effects.
#'   \item \strong{Select one subtype per block} shared across traits:
#'     \code{selected(b) = argmax_subtype Σ_t w_t · value(subtype, t)}, where
#'     \code{value} is \code{S} (\code{value_source = "raw"}) or \code{N}
#'     (\code{"normalized"}), \code{w_t} is the signed \code{trait_weights}, and
#'     \code{fitness_fn} filters candidates / defines the per-subtype value.
#'   \item \strong{Report each trait's own value}: the selected subtype's per-trait
#'     \code{S} and \code{N}, and per-trait top blocks (\code{argmax_b value}). No
#'     cross-trait composite, corrected effect, or optimizing potential is formed.
#' }
#'
#' @param geno_meta A valid \code{GenoMeta} (with block assignments). Required
#'   unless \code{catalogue} is given.
#' @param effect_list Named list of per-trait marker-effect data.frames
#'   (\code{Marker, Allele, Effect}). Required unless \code{catalogue} is given.
#' @param catalogue Optional result list from \code{\link{run_haplotype_catalogue}}
#'   (one element per trait, each carrying \code{hap_catalogue}). When supplied it
#'   is used directly; otherwise the catalogues are computed from \code{geno_meta}
#'   + \code{effect_list}.
#' @param trait_names Optional character vector restricting / ordering the traits.
#'   Defaults to the names of \code{catalogue} or \code{effect_list}.
#' @param selected_blocks Optional character vector of block IDs to keep.
#' @param fitness_fn Candidate filter / per-subtype value definition used before
#'   the cross-trait selection sum (it no longer does the cross-trait combine):
#'   \code{"max_effect"} (default; value = the subtype's effect),
#'   \code{"frequent_best"} (restrict candidates to \code{Count >= min_count}),
#'   \code{"expected"} (value = frequency-weighted effect \code{Effect*Count/N}),
#'   \code{"min_effect"} (pick the lowest summed value), or a function
#'   \code{function(subtypes_dt) -> one row} receiving a block's long subtype rows.
#' @param min_count Minimum \code{Count} for \code{fitness_fn = "frequent_best"}.
#'   Default \code{1}.
#' @param value_source The value used \strong{both} to sum-select the
#'   representative subtype and to fill the rectangular / circos plots:
#'   \code{"raw"} (default; the raw effect \code{S}) or \code{"normalized"} (the
#'   rescaled value \code{N}). Use \code{"normalized"} when traits sit on very
#'   different scales, so the selection sum (and plot colour) is not dominated by
#'   the largest-scale trait; it promotes \code{normalize = "zscore"} when
#'   \code{normalize} was left at its default. Stored on the result as
#'   \code{value_source} so the plot layer fills by the same value without the
#'   caller re-specifying it.
#' @param normalize Per-trait scale transform applied \strong{within each trait}:
#'   \co\code{"none"} (default — no rescaling; \code{N == S} in raw units),
#'   \code{"zscore"}, \code{"minmax"}, \code{"rank"}, or a function
#'   \code{function(numeric) -> numeric}. \code{"zscore"} is a \strong{true
#'   standard score} \code{(x - mean)/sd} — it centres and scales, so a normalised
#'   \code{0} means an \emph{average}-effect subtype for that trait (the diverging
#'   plot fill's white midpoint then reads "average" under a normalized value,
#'   versus "zero effect" under a raw one). \code{"minmax"} and \code{"rank"} stay
#'   \strong{zero-anchored} (origin and sign preserved): \code{"minmax"} divides by
#'   the largest absolute effect (range \code{[-1, 1]}), and \code{"rank"} is a
#'   sign-preserving magnitude rank.
#' @param trait_weights Optional named numeric vector giving each trait's
#'   \strong{signed direction / scale}, applied to \code{value(subtype, t)} before
#'   the cross-trait selection sum — \strong{not} a composite metric. \code{NULL}
#'   (default) uses equal unit weights (\code{1}); a negative weight flips a
#'   lower-is-better trait so its subtypes are summed in the intended direction.
#'   When supplied it must be named and cover every trait.
#' @param center        Centring mode applied to the allele effects:
#'   \code{"none"}, \code{"locus_mean"} or \code{"allele_freq"}. See
#'   \code{\link{compute_profile_localgebv}} for what each preserves.
#' @param center_by_freq Passed to \code{\link{compute_block_catalogue}} on the
#'   compute path. Default \code{TRUE}.
#' @param output_dir Directory for output \code{*.txt} tables. \code{NULL}
#'   (default) writes nothing and emits a \code{warning()}; the full result is
#'   still returned.
#' @param write_outputs Logical. When \code{FALSE}, suppress disk writes even if
#'   \code{output_dir} is set (no warning). Default \code{TRUE}.
#' @param plots Character vector naming the figures to build into \code{$plots}:
#'   any of \code{"rectangular"} (\code{\link{plot_block_effect_genome}}),
#'   \code{"manhattan"} (\code{\link{plot_block_effect_manhattan}}),
#'   \code{"density"} (\code{\link{plot_block_effect_density}}), and
#'   \code{"heatmap"} (\code{\link{plot_top_block_heatmap}}). \code{"circos"} is
#'   reserved for a later phase and currently falls back with a message. Pass
#'   \code{character(0)} to skip plotting.
#' @param ... Reserved for forward compatibility.
#'
#' @return Invisibly, a named list:
#' \describe{
#'   \item{\code{subtypes}}{All block subtypes, long, with a \code{Trait} column
#'     and per-(subtype, trait) raw \code{S} and normalized \code{N}.}
#'   \item{\code{selected}}{The one representative subtype chosen per block (one
#'     row per block), with its carrier and block location.}
#'   \item{\code{S}}{Per-(block, trait) \strong{raw} effect of the selected
#'     subtype, with the per-trait top-block flag \code{is_top}.}
#'   \item{\code{N}}{Per-(block, trait) \strong{rescaled} value of the selected
#'     subtype, with \code{is_top}. Under the default \code{normalize = "none"}
#'     this is \code{S} unchanged.}
#'   \item{\code{top_blocks}}{Per-trait top block(s) (\code{argmax_b} the selected
#'     subtype's value).}
#'   \item{\code{gebv}}{Named list of per-trait wide Block x Sample genotype-level
#'     local-GEBV matrices (used by the top-block heatmap).}
#'   \item{\code{plots}}{Named list of the \code{ggplot} figures requested via
#'     \code{plots}.}
#'   \item{\code{value_source}}{The \code{value_source} used, so the plot layer
#'     fills by the same value.}
#'   \item{\code{normalize}}{The transform actually applied
#'     (\co\code{"none"} by default, or \code{"custom"} for a function), so the plot
#'     layer can label \code{N} honestly rather than assuming a standard score.}
#' }
#'
#' @seealso \code{\link{compute_block_catalogue}},
#'   \code{\link{run_haplotype_catalogue}}, \code{\link{build_block_location_info}}
#'
#' @export
analyze_block_catalogue <- function(geno_meta       = NULL,
                                    effect_list     = NULL,
                                    catalogue       = NULL,
                                    trait_names     = NULL,
                                    selected_blocks = NULL,
                                    fitness_fn      = c("max_effect",
                                                        "frequent_best",
                                                        "expected", "min_effect"),
                                    min_count       = 1,
                                    normalize       = c("none", "zscore",
                                                        "minmax", "rank"),
                                    value_source    = c("raw", "normalized"),
                                    trait_weights   = NULL,
                                    center          = NULL,
                                    center_by_freq  = NULL,
                                    output_dir      = NULL,
                                    write_outputs   = TRUE,
                                    plots           = c("rectangular", "manhattan",
                                                        "density", "heatmap",
                                                        "circos"),
                                    ...) {

  if (!is.function(fitness_fn)) fitness_fn <- match.arg(fitness_fn)
  normalize_given <- !missing(normalize)
  if (!is.function(normalize))  normalize  <- match.arg(normalize)
  value_source <- match.arg(value_source)

  # Raw by default: nothing is rescaled unless the caller asks for it. Asking for
  # the normalized value_source IS asking, so it promotes the "zscore" transform
  # rather than handing back raw numbers under a "normalized" label. Naming
  # normalize = "none" alongside it is a contradiction, so say so.
  if (identical(value_source, "normalized") && !is.function(normalize) &&
      identical(normalize, "none")) {
    if (normalize_given)
      stop("analyze_block_catalogue: value_source = \"normalized\" with ",
           "normalize = \"none\" has no meaning — the normalized value N would ",
           "be the raw value S. Name a transform (\"zscore\", \"minmax\", ",
           "\"rank\", or a function), or use value_source = \"raw\".")
    normalize <- "zscore"
    message("analyze_block_catalogue: value_source = \"normalized\" requested, ",
            "so normalize = \"zscore\" is used (raw is the default otherwise). ",
            "Pass 'normalize' to choose another transform.")
  }
  norm_label <- if (is.function(normalize)) "custom" else normalize

  # ── FR1: resolve per-trait subtype tables (+ GEBV matrices for heatmap) ────
  center      <- .resolve_center(center, center_by_freq,
                                 "analyze_block_catalogue",
                                 default = "allele_freq")
  resolved    <- .resolve_subtype_tables(geno_meta, effect_list, catalogue,
                                         trait_names, center)
  tabs        <- resolved$tabs
  gebv_mats   <- resolved$gebv
  trait_names <- names(tabs)
  n_traits    <- length(trait_names)
  if (n_traits == 0L)
    stop("analyze_block_catalogue: no traits to analyse.")
  multi <- n_traits > 1L

  # ── FR2: stack into one long subtype table with a Trait column ─────────────
  loc_cols <- c("Chrom", "Pos", "Block_start", "Block_end",
                "Num_SNPs_in_Block", "Block_length")
  subtypes <- data.table::rbindlist(lapply(trait_names, function(tr) {
    dt <- data.table::copy(tabs[[tr]])
    dt[, Trait := tr]
    # A stable subtype id: content signature when available, else carrier label.
    if ("signature" %in% names(dt)) dt[, Subtype := signature]
    else dt[, Subtype := paste0("h", Haplotype, "_e", formatC(Effect, format = "g"))]
    dt
  }), fill = TRUE)

  if (!is.null(selected_blocks))
    subtypes <- subtypes[Block %in% selected_blocks]
  if (nrow(subtypes) == 0L)
    stop("analyze_block_catalogue: no subtypes after filtering.")

  data.table::setorder(subtypes, Block, Trait, -Count)

  # Raw per-(subtype, trait) effect S, and its within-trait rescaled value N.
  # The transform runs per trait whatever it is — under the default
  # normalize = "none" it is the identity, so N == S and the reported N table is
  # the raw one in raw units.
  subtypes[, S := as.numeric(Effect)]
  subtypes[, N := .fit_normalizer(S, normalize)$values, by = Trait]

  # value(subtype, t) used for both selection and plot fill (raw S or normalized N).
  subtypes[, Vbase := if (value_source == "normalized") N else S]

  # Trait weights: signed direction/scale only (default equal unit weights). A
  # constant scaling of every weight leaves the argmax selection unchanged.
  if (is.null(trait_weights)) {
    w <- stats::setNames(rep(1, n_traits), trait_names)
  } else {
    if (is.null(names(trait_weights)))
      stop("analyze_block_catalogue: 'trait_weights' must be named.")
    miss <- setdiff(trait_names, names(trait_weights))
    if (length(miss) > 0L)
      stop("analyze_block_catalogue: 'trait_weights' missing trait(s): ",
           paste(miss, collapse = ", "), ".")
    w <- trait_weights[trait_names]
  }

  # ── Step 1: one representative subtype per block, shared across traits ──────
  chosen <- .select_block_subtype(subtypes, fitness_fn, min_count, w)

  bl <- unique(subtypes[, c("Block", loc_cols), with = FALSE], by = "Block")

  # `selected`: one row per block — the chosen subtype, a representative carrier,
  # and the block location.
  sub_rep  <- unique(subtypes[, .(Block, Subtype, Count, Sample, Haplotype)],
                     by = c("Block", "Subtype"))
  selected <- merge(chosen, sub_rep, by = c("Block", "Subtype"))
  selected <- merge(selected, bl, by = "Block", all.x = TRUE)
  data.table::setorder(selected, Block)

  # ── Step 2/3: expand the selected subtype to its per-trait S and N ─────────
  bt <- merge(chosen,
              subtypes[, .(Block, Subtype, Trait, S, N, Count, Chrom, Pos,
                           Block_start, Block_end)],
              by = c("Block", "Subtype"))
  # Per-trait top block(s): argmax over blocks of the selected subtype's value
  # (raw and normalized agree within a trait, so flagged on S).
  bt[, is_top := S == max(S, na.rm = TRUE), by = Trait]
  data.table::setorder(bt, Block, Trait)

  S_tab <- bt[, .(Block, Trait, Subtype, Count, S, Chrom, Pos,
                  Block_start, Block_end, is_top)]
  N_tab <- bt[, .(Block, Trait, Subtype, Count, N, Chrom, Pos,
                  Block_start, Block_end, is_top)]

  top_blocks <- data.table::copy(bt[is_top == TRUE])
  data.table::setorder(top_blocks, Trait, -S)

  # Identification columns lead every returned/written analysis table
  # (Trait, Chrom, Pos, Block, Sample, Haplotype, signature/Subtype …).
  for (tb in list(subtypes, selected, S_tab, N_tab, top_blocks))
    .lead_identifier_cols(tb)

  result <- list(subtypes     = subtypes[],
                 selected     = selected[],
                 S            = S_tab[],
                 N            = N_tab[],
                 top_blocks   = top_blocks[],
                 gebv         = gebv_mats,
                 plots        = list(),
                 value_source = value_source,
                 normalize    = norm_label)

  # ── FR3/FR4: build the requested plots over the analysis result ────────────
  result$plots <- .build_analysis_plots(result, plots, multi)

  # ── FR5: gated disk writes ─────────────────────────────────────────────────
  if (is.null(output_dir)) {
    warning("analyze_block_catalogue: no 'output_dir' given; results returned ",
            "in-session only, nothing written.")
  } else if (isTRUE(write_outputs)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    data.table::fwrite(subtypes,   file.path(output_dir, "block_subtypes.txt"), sep = "\t")
    data.table::fwrite(selected,   file.path(output_dir, "block_selected.txt"), sep = "\t")
    data.table::fwrite(S_tab,      file.path(output_dir, "block_S.txt"),        sep = "\t")
    data.table::fwrite(N_tab,      file.path(output_dir, "block_N.txt"),        sep = "\t")
    data.table::fwrite(top_blocks, file.path(output_dir, "block_top.txt"),      sep = "\t")
    for (nm in names(result$plots)) {
      p   <- result$plots[[nm]]
      png <- file.path(output_dir, paste0("block_", nm, ".png"))
      if (inherits(p, "ggplot")) {
        suppressMessages(ggplot2::ggsave(png, p, width = 11, height = 6,
                                         dpi = 150))
      } else if (inherits(p, "recordedplot")) {
        grDevices::png(png, width = 1100, height = 1100, res = 150)
        grDevices::replayPlot(p)
        grDevices::dev.off()
      }
    }
  }

  invisible(result)
}
