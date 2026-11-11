# =============================================================================
# catalogue_engine.R
#
# Shared, trait-invariant compute core for the haplotype-catalogue pipeline, plus
# the local-GEBV engine that turns one trait's marker effects into a wide
# Block x {Sample | SampleHap} matrix.
#
# Why this file exists (catalogue-perf upgrade, 2026-06-26):
#   The pre-refactor pipeline rebuilt the long-format allele table 2*T times per
#   run (once per trait at genotype level inside compute_profile_localgebv(),
#   again at haplotype level inside compute_block_catalogue()). Every rebuild
#   re-ran geno_meta_views() and re-melted the whole genotype matrix, even though
#   none of that work depends on the marker effects. .catalogue_core() does it
#   ONCE; per-trait work then collapses to an effect join + a per-block sum.
#
# Engine dispatch:
#   .local_gebv_engine() computes the per-block haplotype effects via either a
#   compiled Rcpp kernel (local_gebv_block_sum_cpp, src/local_gebv.cpp) or a pure
#   data.table fallback, selected by .hs_use_rcpp(). Both paths are required to
#   produce identical results (parity-tested in test-catalogue-engine.R). The
#   pre-refactor implementations are archived verbatim in
#   dev/codebank/2026-06-26_catalogue_compute.R.
# =============================================================================


# ---------------------------------------------------------------------------
# Shared column-ordering convention for block/catalogue output tables
# ---------------------------------------------------------------------------

# Identification columns, in the order they should lead an output table. Any of
# these present in a table are moved to the front in this order; every other
# column keeps its existing relative order behind them. Genomic coordinates lead
# (matching a Chrom/Pos/Block/Marker/Sample reading order); Trait, when present,
# is the outermost key and comes first.
# @keywords internal

# Identification columns, in the order they should lead an output table. Any of
# these present in a table are moved to the front in this order; every other
# column keeps its existing relative order behind them. Genomic coordinates lead
# (matching a Chrom/Pos/Block/Marker/Sample reading order); Trait, when present,
# is the outermost key and comes first.
# @keywords internal
.HS_ID_COLS <- c("Trait", "Chrom", "Pos", "Block", "Marker", "Sample",
                 "Haplotype", "signature", "Subtype")


# ---------------------------------------------------------------------------
# Rcpp availability / global toggle
# ---------------------------------------------------------------------------

# Resolve ANY Rcpp-exported kernel from the package namespace by symbol name, or
# NULL when it was never compiled in (the package was installed without a
# compiler / without the src/ sources). One lookup strategy for every kernel.
# @keywords internal
.hs_kernel <- function(name) {
  fn <- tryCatch(
    get0(name, envir = asNamespace("haploStack"), inherits = FALSE,
         ifnotfound = NULL),
    error = function(e) NULL)
  if (is.function(fn)) fn else NULL
}


# Is the named Rcpp kernel usable right now? TRUE only when (a) the user has not
# disabled compiled kernels via options(haploStack.use_rcpp = FALSE) and (b) the
# symbol is actually loadable. The single option governs ALL kernels; any
# consumer MUST retain a pure-R fallback for when this returns FALSE. The option
# is read through .hs_setting() (R/state.R); it stays a process option because
# it picks a code path, never a result.
# @keywords internal
.hs_use_rcpp <- function(name) {
  isTRUE(.hs_setting("use_rcpp")) && is.function(.hs_kernel(name))
}


# ---------------------------------------------------------------------------
# Effect-table standardisation (shared by both public wrappers)
# ---------------------------------------------------------------------------


# Standardise a raw effect_df to a data.table(Marker, Allele[char], Effect[num])
# read BY POSITION from its first three columns, unique by (Marker, Allele).
.standardise_effect_dt <- function(effect_df, caller) {
  if (!is.data.frame(effect_df) && !inherits(effect_df, "data.table"))
    .hs_abort(caller, ": 'effect_df' must be a data.frame or data.table.")
  if (ncol(effect_df) < 3L)
    .hs_abort(caller, ": 'effect_df' must have >= 3 columns.")
  eff <- as.data.table(effect_df[, 1:3, drop = FALSE])
  setnames(eff, c("Marker", "Allele", "Effect"))
  eff[, Allele := as.character(Allele)]
  eff[, Effect := as.numeric(Effect)]
  unique(eff, by = c("Marker", "Allele"))
}


# ---------------------------------------------------------------------------
# Trait-invariant compute core
# ---------------------------------------------------------------------------


#' Build the trait-invariant catalogue compute core (once per run)
#'
#' Derives everything the per-block GEBV summation needs that does NOT depend on
#' marker effects: the long-format allele table (Sample, Haplotype_idx, Marker,
#' Allele, Block) with NA/blank alleles dropped and unassigned markers removed,
#' the population allele-frequency table, the block/sample output orderings, and
#' — when \code{use_rcpp} is \code{TRUE} — the integer-coded matrix + per-marker
#' allele dictionaries + marker→block index consumed by the Rcpp kernel.
#'
#' All allele-encoding handling lives here, in R: the genotype cells and (later)
#' the effect alleles are mapped into the SAME per-marker integer code space, so
#' the compiled kernel only ever sees integers and doubles and is therefore
#' encoding-agnostic. An effect allele absent from a marker's observed states
#' simply has no code and contributes zero — exactly the silent-drop behaviour of
#' the pre-refactor update join.
#'
#' @param geno_meta A valid \code{GenoMeta} (already filtered if the caller wants
#'   a subset; the core is built over whatever \code{geno_meta} contains).
#' @param use_rcpp  Logical; build the integer-coded structures for the Rcpp
#'   kernel. Default: the compiled kernel when enabled and loadable.
#' @return A list with elements: \code{empty} (logical), \code{encoding},
#'   \code{long_core}, \code{freq_core}, \code{block_order}, \code{sample_order},
#'   \code{present_copies}, and (when \code{use_rcpp}) \code{code_mat},
#'   \code{contrib_levels}, \code{block_idx}, \code{block_levels},
#'   \code{copy_sample}, \code{copy_hap}, \code{n_code}.
#' @keywords internal
#' @noRd
.catalogue_core <- function(geno_meta,
                            use_rcpp = .hs_use_rcpp("local_gebv_block_sum_cpp"),
                            views = NULL) {
  # `views`: a geno_meta_views(geno_meta) the caller already holds (cc/F09:
  # run_haplotype_catalogue() builds the views once for the core and the
  # signatures). NULL builds them here.
  prof      <- if (is.null(views)) geno_meta_views(geno_meta) else views
  meta_cols <- c("Sample", "Haplotype_idx")
  all_mks   <- setdiff(names(prof$alleles), meta_cols)
  block_mks <- prof$marker_meta[!is.na(block_id), marker]
  markers   <- intersect(all_mks, block_mks)

  empty_core <- list(empty = TRUE, encoding = prof$encoding)
  if (nrow(prof$alleles) == 0L || nrow(prof$marker_meta) == 0L ||
      length(markers) == 0L)
    return(empty_core)

  # 0.9.0 (row 29 a): blocks by the one natural chromosome order, then position
  block_order  <- prof$block_meta[order(.hs_chrom_rank(chromosome), start_pos), block_id]
  sample_order <- prof$samples

  # ── Long-format allele table (the single expensive melt, done ONCE) ────────
  long_core <- data.table::melt(
    prof$alleles[, c(meta_cols, markers), with = FALSE],
    id.vars         = meta_cols,
    measure.vars    = markers,
    variable.name   = "Marker",
    value.name      = "Allele",
    variable.factor = FALSE)
  setDT(long_core)
  long_core[, Allele := as.character(Allele)]
  long_core[, Marker := as.character(Marker)]
  # The third of the three missing-detection idioms this release removes. It
  # now uses the shared token set rather than its own `nzchar(trimws())` rule
  # (FR7). On a container that came through the portal the two agree exactly —
  # the portal already normalises whitespace-only cells — so this is a no-op.
  # On a legacy container it stops a retained token being counted in
  # `MarkerTotal` below, which is what deflated every real allele's frequency
  # and shifted local GEBVs under `center_by_freq = TRUE`.
  .core_miss <- .missing_tokens(prof$encoding, .meta_extra_missing(geno_meta))
  long_core  <- long_core[!is.na(Allele) & !trimws(Allele) %chin% .core_miss]

  blk_map <- prof$marker_meta[, .(Marker = marker, Block = block_id)]
  long_core[, Block := NA_character_]
  long_core[blk_map, on = "Marker", Block := i.Block]
  long_core <- long_core[!is.na(Block)]

  if (nrow(long_core) == 0L) return(empty_core)

  # ── Allele-frequency table (trait-invariant; reused across traits) ─────────
  freq_core <- long_core[, .(AlleleCount = .N), by = .(Marker, Allele)]
  tot_core  <- long_core[, .(MarkerTotal = .N), by = Marker]
  freq_core <- freq_core[tot_core, on = "Marker"][
    , AlleleFreq := AlleleCount / MarkerTotal][
    , c("AlleleCount", "MarkerTotal") := NULL]

  present_copies <- unique(long_core[, paste0(Sample, ".", Haplotype_idx)])

  core <- list(
    empty          = FALSE,
    encoding       = prof$encoding,
    long_core      = long_core,
    freq_core      = freq_core,
    block_order    = block_order,
    sample_order   = sample_order,
    present_copies = present_copies)

  if (!isTRUE(use_rcpp)) return(core)

  # ── Integer-coded structures for the Rcpp kernel (built only when needed) ──
  # code_mat[i, j] = 1-based allele code at copy i, marker j (NA where absent).
  # block_idx[j]   = 1-based index of marker j's block into block_levels (0=skip,
  #                  i.e. the block has no present allele anywhere).
  # contrib_levels = once-built (j, Marker, code, Allele) table; per trait it is
  #                  joined to the effect table to form the kernel's contrib_mat.
  alleles_dt <- prof$alleles
  n_copies   <- nrow(alleles_dt)
  n_markers  <- length(markers)
  code_mat   <- matrix(NA_integer_, n_copies, n_markers)
  lev_list   <- vector("list", n_markers)
  for (j in seq_len(n_markers)) {
    vals <- as.character(alleles_dt[[markers[j]]])
    vals[is.na(vals) | nchar(trimws(vals)) == 0L] <- NA_character_
    f             <- factor(vals)
    code_mat[, j] <- as.integer(f)
    lev_list[[j]] <- levels(f)
  }

  block_levels <- unique(long_core$Block)
  block_of_mk  <- prof$marker_meta[match(markers, marker), block_id]
  block_idx    <- match(block_of_mk, block_levels)
  block_idx[is.na(block_idx)] <- 0L

  contrib_levels <- data.table::rbindlist(lapply(seq_len(n_markers), function(j) {
    lv <- lev_list[[j]]
    if (!length(lv)) return(NULL)
    data.table(j = j, Marker = markers[j], code = seq_along(lv), Allele = lv)
  }))

  core$code_mat       <- code_mat
  core$contrib_levels <- contrib_levels
  core$block_idx      <- as.integer(block_idx)
  core$block_levels   <- block_levels
  core$copy_sample    <- as.character(alleles_dt$Sample)
  core$copy_hap       <- as.integer(alleles_dt$Haplotype_idx)
  core$n_code         <- max(vapply(lev_list, length, integer(1L)), 0L)
  core
}


# ---------------------------------------------------------------------------
# Per-block haplotype effects (the dispatched hot loop)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Centring vocabulary (CEN, 0.6.0)
# ---------------------------------------------------------------------------


# The three centrings a local GEBV can be computed under. They are different
# operations, not settings of one:
#
#   "none"        contrib(a) = Effect(a)
#   "locus_mean"  contrib(a) = Effect(a) - mean_i,  mean_i = sum_a p_a * e_a
#                 ONE value per marker, subtracted from every allele there. A
#                 pure location shift: the A-vs-B contrast is preserved and
#                 block variance is unchanged. Reproduces HapSelect's
#                 mean_adjust = TRUE value for value.
#   "allele_freq" contrib(a) = Effect(a) * (1 - p_a)
#                 Each allele centred by ITS OWN frequency, so the contrast
#                 moves from e_A - e_B to q*e_A - p*e_B. A reweighting: block
#                 variances change (median factor ~3 on the barley panel).
HS_CENTER_MODES <- c("none", "locus_mean", "allele_freq")


# One-time-per-session deprecation notice for the retired logical.
# One table of centring defaults, looked up by the calling entry point (design
# D3 §3.7, decision-form row 15 c, 0.9.0). Every catalogue / selection function
# defaults to "none" (Q22.h); count_favourable_blocks() to "locus_mean": its
# output is a SIGN count, and an uncentred sign mostly reflects which allele the
# effect file called the reference, while locus_mean is a per-block location
# shift that leaves every subtype argmax unchanged. Pin the 0.8.8 numbers with
# center = "allele_freq". An unknown caller gets "none".
.HS_CENTER_DEFAULT <- c(
  run_haplotype_catalogue   = "none",
  compute_profile_localgebv = "none",
  predict_gebv              = "none",
  compute_block_catalogue   = "none",
  analyze_block_catalogue   = "none",
  count_block_subtypes      = "none",
  count_favourable_subtypes = "none",
  count_favourable_blocks   = "locus_mean"
)


# Resolve the centring from the new argument and the deprecated logical.
# `default` is what applies when the caller supplies neither — it differs by
# entry point, which is the disagreement CEN removes.
.resolve_center <- function(center = NULL, center_by_freq = NULL,
                            caller = "compute_profile_localgebv") {
  # Defunct since 0.9.0 (row 26 a; deprecated since 0.6.0). The argument stays
  # in the eight signatures for this release so the caller is told the
  # replacement; the message names the mode the logical used to mean.
  if (!is.null(center_by_freq)) {
    new <- if (isTRUE(center_by_freq)) "center = \"allele_freq\""
           else if (isFALSE(center_by_freq)) "center = \"none\""
           else "center"
    .hs_defunct("center_by_freq", new, "0.9.0", dep_since = "0.6.0",
                what = "argument", fn = caller)
  }

  if (is.null(center)) {
    default <- .HS_CENTER_DEFAULT[caller]
    return(if (is.na(default)) "none" else unname(default))
  }
  if (is.logical(center))
    .hs_abort(caller, ": 'center' takes a mode name, not a logical. Use ",
         "center = \"allele_freq\" or center = \"none\".", call. = FALSE)
  match.arg(as.character(center)[1L], HS_CENTER_MODES)
}


# Per-marker locus mean:  mean_i = sum_a AlleleFreq(a) * Effect(a).
# Built from the same frequency table "allele_freq" uses, aggregated one level
# up (per marker rather than per allele). Alleles with no supplied effect count
# as effect 0, matching the effect join everywhere else.
.locus_means <- function(freq_dt, eff) {
  lm_dt <- data.table::copy(as.data.table(freq_dt))
  lm_dt[, Effect := NA_real_]
  lm_dt[eff, on = .(Marker, Allele), Effect := i.Effect]
  lm_dt[is.na(Effect), Effect := 0]
  lm_dt[, .(LocusMean = sum(AlleleFreq * Effect, na.rm = TRUE)), by = Marker]
}


# Pure-R path: reuse the archived update-join + grouped-sum on a COPY of the
# cached long table (so the shared core is never mutated). Identical arithmetic
# and grouping/summation order to the pre-refactor .compute_effect_size().
.hap_effects_r <- function(core, eff, center) {
  .compute_effect_size(data.table::copy(core$long_core), eff,
                       center = center, freq_dt = core$freq_core)
}


# Rcpp path: build a per-(marker, code) contribution matrix in R (this is where
# encoding is resolved), scatter-sum it per block in C++, then reshape to the
# same (Sample, Haplotype_idx, Block, HapEffect) long form the pure-R path emits.
.hap_effects_rcpp <- function(core, eff, center) {
  ld <- data.table::copy(core$contrib_levels)
  ld[eff, on = .(Marker, Allele), Effect := i.Effect]
  ld[is.na(Effect), Effect := 0]
  if (identical(center, "allele_freq")) {
    ld[core$freq_core, on = .(Marker, Allele), AlleleFreq := i.AlleleFreq]
    ld[is.na(AlleleFreq), AlleleFreq := 0]
    ld[, contrib := Effect * (1 - AlleleFreq)]          # == Eff - freq*Eff
  } else if (identical(center, "locus_mean")) {
    ld[.locus_means(core$freq_core, eff), on = "Marker",
       LocusMean := i.LocusMean]
    ld[is.na(LocusMean), LocusMean := 0]
    ld[, contrib := Effect - LocusMean]
  } else {
    ld[, contrib := Effect]
  }

  n_markers   <- ncol(core$code_mat)
  contrib_mat <- matrix(0, nrow = n_markers, ncol = max(core$n_code, 1L))
  contrib_mat[cbind(ld$j, ld$code)] <- ld$contrib

  res <- .hs_kernel("local_gebv_block_sum_cpp")(core$code_mat, contrib_mat, core$block_idx,
                        length(core$block_levels))
  colnames(res) <- core$block_levels
  dt <- as.data.table(res)
  dt[, Sample := core$copy_sample]
  dt[, Haplotype_idx := core$copy_hap]
  hap_dt <- data.table::melt(dt, id.vars = c("Sample", "Haplotype_idx"),
                             variable.name = "Block", value.name = "HapEffect",
                             variable.factor = FALSE)
  # Drop entirely-absent haplotype copies so haplotype-level columns match the
  # pure-R output (genotype level is unaffected; absent copies sum to 0).
  hap_dt[paste0(Sample, ".", Haplotype_idx) %in% core$present_copies]
}


.hap_effects <- function(core, eff, center, use_rcpp) {
  if (isTRUE(use_rcpp)) .hap_effects_rcpp(core, eff, center)
  else                  .hap_effects_r(core, eff, center)
}


# ---------------------------------------------------------------------------
# Wide-matrix assembly (identical for both engines and both levels)
# ---------------------------------------------------------------------------


# Turn (Sample, Haplotype_idx, Block, HapEffect) into the wide Block x {Sample |
# SampleHap} data.table, preserving the pre-refactor ordering and fill semantics.
# keep_blocks (optional) restricts the output blocks (used when a shared full
# core serves a block-subsetted compute_profile_localgebv call).
.assemble_localgebv_wide <- function(hap_dt, core, level, fill_na,
                                     keep_blocks = NULL) {
  if (!is.null(keep_blocks)) hap_dt <- hap_dt[Block %in% keep_blocks]
  block_order  <- core$block_order
  sample_order <- core$sample_order

  if (level == "haplotype") {
    n_hap <- max(hap_dt$Haplotype_idx, na.rm = TRUE)
    hap_dt[, SampleHap := paste0(Sample, ".", Haplotype_idx)]
    gebv_wide <- data.table::dcast(hap_dt, Block ~ SampleHap,
                                   value.var = "HapEffect", fill = fill_na)
    all_hap_keys <- as.character(unlist(lapply(sample_order,
      function(s) paste0(s, ".", seq_len(n_hap)))))
    valid_bo  <- block_order[block_order %in% gebv_wide$Block]
    valid_hk  <- all_hap_keys[all_hap_keys %in% names(gebv_wide)]
    gebv_wide <- gebv_wide[match(valid_bo, gebv_wide[["Block"]])]
    setcolorder(gebv_wide, c("Block", valid_hk))
  } else {
    geno_dt <- hap_dt[, .(GenoEffect = sum(HapEffect, na.rm = TRUE)),
                      by = .(Sample, Block)]
    gebv_wide <- data.table::dcast(geno_dt, Block ~ Sample,
                                   value.var = "GenoEffect", fill = fill_na)
    valid_bo  <- block_order[block_order %in% gebv_wide$Block]
    valid_so  <- sample_order[sample_order %in% names(gebv_wide)]
    gebv_wide <- gebv_wide[match(valid_bo, gebv_wide[["Block"]])]
    setcolorder(gebv_wide, c("Block", valid_so))
  }
  gebv_wide
}


#' Local-GEBV engine: one trait's effects -> wide Block x Sample (or SampleHap) table
#'
#' The single seam both \code{compute_profile_localgebv()} and
#' \code{compute_block_catalogue()} call. Dispatches the per-block summation to
#' the Rcpp kernel or the pure-R path, then assembles the wide matrix.
#'
#' @param core     A \code{.catalogue_core()} result.
#' @param eff      Standardised effect data.table (Marker, Allele, Effect).
#' @param center   Centring mode: \code{"none"}, \code{"locus_mean"} or
#'   \code{"allele_freq"} (see \code{HS_CENTER_MODES}). Already resolved by the
#'   caller — the engine does not read the deprecated logical.
#' @param level    \code{"genotype"} or \code{"haplotype"}.
#' @param fill_na  Fill value for absent (Block, Sample) cells.
#' @param keep_blocks Optional character vector restricting output blocks.
#' @keywords internal
#' @noRd
.local_gebv_engine <- function(core, eff, center = "none",
                               level = c("genotype", "haplotype"),
                               fill_na = 0, keep_blocks = NULL) {
  level  <- match.arg(level)
  center <- match.arg(as.character(center)[1L], HS_CENTER_MODES)
  if (isTRUE(core$empty)) return(data.table(Block = character(0L)))
  # The Rcpp path emits a dense [copy x block] result, so absent (copy, block)
  # cells read as 0; that matches the pure-R fill ONLY when fill_na == 0. For any
  # other fill value (a non-default direct call), use the pure-R path which omits
  # absent cells and lets dcast fill them. Requires the kernel + coded core.
  use_rcpp <- .hs_use_rcpp("local_gebv_block_sum_cpp") && fill_na == 0 &&
    !is.null(core$code_mat)
  hap_dt   <- .hap_effects(core, eff, center, use_rcpp)
  if (nrow(hap_dt) == 0L) return(data.table(Block = character(0L)))
  .assemble_localgebv_wide(hap_dt, core, level, fill_na, keep_blocks)
}
