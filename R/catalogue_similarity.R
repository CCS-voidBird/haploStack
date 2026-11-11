# =============================================================================
# catalogue_similarity.R
#
# Haplotype-block similarity analysis (R2): IBS of haplotype segments,
# family-history-driven block IBD, and a marker-level signal test.
#
# Functions:
#   compute_segment_ibs()   — IBS ratio of two haplotype segments (R2.1)
#   compute_block_ibd()     — per-block parent→progeny IBD, all pairwise (R2.2)
#   score_marker_signal()   — marker-level TP/FP/TN/FN signal test (R2.3)
#   compute_signal_metrics()— classification metrics + MAF final score (R2.4)
#   select_signal_markers() — top-N signal marker(s) per block (R2.5)
#   attach_signal_markers() — record the selection on geno_meta$block_info
#   get_signal_markers()    — read the selection back off a GenoMeta
#
# Internal:
#   .resolve_families()        — family-history / pedigree → family list
#   .compute_similarity_long() — shared per-comparison marker-level table
#   .marker_maf()              — per-marker genome-wide minor allele frequency
#   .signal_marker_cols()      — the block_info column names written by attach
#   .prune_signal_markers()    — drop selections whose marker no longer exists
#
# Comparison model (Q10): anchored on parental haplotypes — for each of a
# family's n_parents × ploidy parental haplotypes (e.g. a double-cross family
# with 4 parents → up to 8), compare against every progeny haplotype block and
# keep every pairwise (ParentHap × ProgenyHap) comparison (no best-match
# collapsing). IBD is inclusive: Similarity >= threshold.
# =============================================================================


# ---------------------------------------------------------------------------
# R2.1 — segment IBS
# ---------------------------------------------------------------------------

#' Identity-by-state of two haplotype segments
#'
#' Compares two equal-length haplotype segments (allele vectors) and returns the
#' fraction of comparable markers whose alleles are identical. Positions where
#' either segment is \code{NA}/empty are excluded from both numerator and
#' denominator.
#'
#' @param seg1,seg2 Character vectors of alleles of equal length (one haplotype
#'   copy across the markers of a block/segment).
#' @param block_length Optional numeric physical (bp) length of the block, used
#'   only when \code{weight_by_length = TRUE}.
#' @param weight_by_length Logical. When \code{TRUE}, the match/total ratio is
#'   scaled by \code{block_length}. Default \code{FALSE}.
#'
#' @return A single numeric: the IBS ratio in \code{[0, 1]} (or
#'   \code{ratio * block_length} when weighting); \code{NA_real_} when no
#'   markers are comparable.
#'
#' @examples
#' \dontrun{
#' compute_segment_ibs(c("A","T","G"), c("A","C","G"))            # 0.667
#' compute_segment_ibs(c("A","T"), c("A","T"), block_length = 1e5,
#'                     weight_by_length = TRUE)                    # 1e5
#' }
#'
#' @seealso \code{\link{compute_block_ibd}}, \code{\link{score_marker_signal}}
#' @export
compute_segment_ibs <- function(seg1, seg2, block_length = NULL,
                                weight_by_length = FALSE) {
  if (length(seg1) != length(seg2))
    .hs_abort("compute_segment_ibs: 'seg1' and 'seg2' must have equal length.")

  s1 <- as.character(seg1)
  s2 <- as.character(seg2)
  comparable <- !is.na(s1) & !is.na(s2) & nzchar(s1) & nzchar(s2)
  n <- sum(comparable)
  if (n == 0L) return(NA_real_)

  ratio <- sum(s1[comparable] == s2[comparable]) / n
  if (weight_by_length) {
    if (is.null(block_length))
      .hs_abort("compute_segment_ibs: 'block_length' is required when ",
           "weight_by_length = TRUE.")
    return(ratio * as.numeric(block_length))
  }
  ratio
}


# ---------------------------------------------------------------------------
# Internal — missing-allele tokens (encoding-aware)
# ---------------------------------------------------------------------------

# The missing-allele vocabulary (.missing_tokens(), .meta_extra_missing()) lives
# in R/alphabet.R since 0.9.0. The rule is per ALLELE, so a diploid cell written
# "00" or "--" is a missing genotype by consequence — both copies land in the
# set — while a half-called cell ("A0") keeps its observed allele and marks only
# the absent copy. Since 0.6.0 a container built by the package is already
# normalised, so the set is legacy compatibility for hand-built or pre-0.6.0
# containers, not the mechanism that decides what is missing.


# Emit a progress message (to stderr) only when verbose is TRUE.
.sim_msg <- function(verbose, ...) {
  if (isTRUE(verbose)) .hs_inform(sprintf(...))
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Internal — family resolution
# ---------------------------------------------------------------------------

# Build a uniform list of families: each element is
# list(gen, fid, parents, offspring), filtered to samples present in the
# genotype set. Source is either a family-history list (named gen_*, each a
# data.table(FID, IID, Label)) or a flat pedigree (Sample, Parent1, Parent2).
.resolve_families <- function(history, pedigree, geno_samples) {
  fams <- list()

  if (!is.null(history)) {
    if (is.null(names(history)))
      names(history) <- paste0("gen_", seq_along(history))
    for (gen_name in names(history)) {
      g <- as.data.table(history[[gen_name]])
      for (fid in unique(g$FID)) {
        grp <- g[FID == fid]
        fams[[length(fams) + 1L]] <- list(
          gen       = gen_name,
          fid       = as.character(fid),
          parents   = unique(grp$IID[grp$Label == "P"]),
          offspring = unique(grp$IID[grp$Label == "O"]))
      }
    }
  } else {
    ped <- as.data.table(pedigree)
    if (!all(c("Sample", "Parent1", "Parent2") %in% names(ped)))
      .hs_abort(".resolve_families: pedigree needs columns Sample, Parent1, Parent2.")
    pset <- lapply(seq_len(nrow(ped)), function(i) {
      p <- c(as.character(ped$Parent1[i]), as.character(ped$Parent2[i]))
      .hs_sort(unique(p[!is.na(p) & nzchar(p) & p != "0"]))
    })
    keys <- vapply(pset, paste, collapse = "|", FUN.VALUE = character(1L))
    uk <- unique(keys[nzchar(keys)])
    for (i in seq_along(uk)) {
      idx <- which(keys == uk[i])
      fams[[length(fams) + 1L]] <- list(
        gen       = "ped",
        fid       = paste0("F_", i),
        parents   = pset[[idx[1L]]],
        offspring = unique(as.character(ped$Sample[idx])))
    }
  }

  # Keep only genotyped parents/offspring; drop families that lose either side.
  dropped <- FALSE
  kept <- list()
  for (f in fams) {
    p <- intersect(f$parents,   geno_samples)
    o <- intersect(f$offspring, geno_samples)
    if (length(p) == 0L || length(o) == 0L) { dropped <- TRUE; next }
    f$parents <- p; f$offspring <- o
    kept[[length(kept) + 1L]] <- f
  }
  attr(kept, "dropped") <- dropped
  kept
}


# ---------------------------------------------------------------------------
# Internal — shared per-comparison marker-level similarity table
# ---------------------------------------------------------------------------

# Returns a marker-level data.table with one row per
# (Gen, FID, Parent, ParentHap, Progeny, ProgenyHap, Block, Marker) carrying:
#   Chrom, Pos, marker_same (logical), comparable (logical),
#   Similarity (block-level IBS for that pairwise comparison), IBD (logical).
# Block similarity is the mean per-marker identity over comparable markers in
# the block; IBD = Similarity >= threshold (inclusive).
.compute_similarity_long <- function(geno_meta, family_history, pedigree,
                                     blocks, threshold, verbose = FALSE) {
  gm <- if (!is.null(blocks)) filter_geno_meta(geno_meta, blocks = blocks)
        else geno_meta
  prof <- geno_meta_views(gm)

  empty <- data.table(
    Gen = character(0L), FID = character(0L), Parent = character(0L),
    ParentHap = integer(0L), Progeny = character(0L), ProgenyHap = integer(0L),
    Block = character(0L), Chrom = character(0L), Pos = numeric(0L),
    Marker = character(0L), marker_same = logical(0L),
    comparable = logical(0L), Similarity = numeric(0L), IBD = logical(0L))

  mk_block <- prof$marker_meta[!is.na(block_id)]
  if (nrow(mk_block) == 0L || nrow(prof$alleles) == 0L) {
    .sim_msg(verbose, "[similarity] no block-assigned markers / no genotypes; ",
             "nothing to compare.")
    return(empty)
  }

  ploidy   <- max(prof$alleles$Haplotype_idx, na.rm = TRUE)
  n_blocks <- length(unique(mk_block$block_id))
  .sim_msg(verbose,
           "[similarity] encoding=%s | %d markers in %d blocks | ploidy=%d (num_hap/sample)",
           prof$encoding, nrow(mk_block), n_blocks, ploidy)

  miss <- .missing_tokens(prof$encoding, .meta_extra_missing(gm))
  marker_names <- mk_block$marker
  long <- data.table::melt(
    prof$alleles[, c("Sample", "Haplotype_idx", marker_names), with = FALSE],
    id.vars       = c("Sample", "Haplotype_idx"),
    measure.vars  = marker_names,
    variable.name = "Marker",
    value.name    = "Allele",
    variable.factor = FALSE)
  data.table::setDT(long)
  long[, Allele := as.character(Allele)]
  bm <- mk_block[, .(Marker = marker, Block = block_id,
                     Chrom = chromosome, Pos = position)]
  long[bm, on = "Marker", `:=`(Block = i.Block, Chrom = i.Chrom, Pos = i.Pos)]

  fams <- .resolve_families(family_history, pedigree, prof$samples)
  if (isTRUE(attr(fams, "dropped")))
    .hs_warn("compute similarity: one or more families had no genotyped ",
            "parent/offspring and were skipped.")
  if (length(fams) == 0L) {
    .sim_msg(verbose, "[similarity] no usable families after matching to ",
             "genotypes.")
    return(empty)
  }
  .sim_msg(verbose, "[similarity] %d family(ies) to process; IBD threshold = %.3f (inclusive).",
           length(fams), threshold)

  n_fam <- length(fams)
  per_family <- lapply(seq_len(n_fam), function(fi) {
    f <- fams[[fi]]
    n_par <- length(f$parents); n_off <- length(f$offspring)
    .sim_msg(verbose,
             "  (%d/%d) [%s/%s] %d parent(s) -> %d parental haplotypes x %d offspring (%d progeny haplotypes)",
             fi, n_fam, f$gen, f$fid, n_par, n_par * ploidy, n_off,
             n_off * ploidy)

    par <- long[Sample %in% f$parents,
                .(Parent = Sample, ParentHap = Haplotype_idx,
                  Marker, Block, Chrom, Pos, PAllele = Allele)]
    off <- long[Sample %in% f$offspring,
                .(Progeny = Sample, ProgenyHap = Haplotype_idx,
                  Marker, PrAllele = Allele)]
    if (nrow(par) == 0L || nrow(off) == 0L) {
      .sim_msg(verbose, "      skipped (no genotyped alleles on one side).")
      return(NULL)
    }

    m <- merge(par, off, by = "Marker", allow.cartesian = TRUE)
    if (nrow(m) == 0L) {
      .sim_msg(verbose, "      skipped (no shared markers).")
      return(NULL)
    }

    m[, comparable := !is.na(PAllele) & !is.na(PrAllele) &
                      !(PAllele %in% miss) & !(PrAllele %in% miss)]
    m[, marker_same := comparable & (PAllele == PrAllele)]

    # Block-level similarity per pairwise comparison (over comparable markers)
    m[, Similarity := NA_real_]
    blk <- m[comparable == TRUE,
             .(sim = mean(marker_same)),
             by = .(Parent, ParentHap, Progeny, ProgenyHap, Block)]
    m[blk, on = .(Parent, ParentHap, Progeny, ProgenyHap, Block),
      Similarity := i.sim]
    m[, IBD := !is.na(Similarity) & Similarity >= threshold]
    m[, `:=`(Gen = f$gen, FID = f$fid)]

    .sim_msg(verbose,
             "      -> %d block comparisons across %d block(s); %.1f%% IBD (mean similarity %.3f)",
             nrow(blk), length(unique(blk$Block)),
             100 * mean(blk$sim >= threshold, na.rm = TRUE),
             mean(blk$sim, na.rm = TRUE))

    m[, .(Gen, FID, Parent, ParentHap, Progeny, ProgenyHap, Block, Chrom, Pos,
          Marker, marker_same, comparable, Similarity, IBD)]
  })

  out <- data.table::rbindlist(Filter(Negate(is.null), per_family))
  if (nrow(out) == 0L) return(empty)
  .sim_msg(verbose, "[similarity] done: %d marker-level rows from %d family(ies).",
           nrow(out), n_fam)
  out
}


# ---------------------------------------------------------------------------
# R2.2 — block IBD
# ---------------------------------------------------------------------------

#' Per-block identity-by-descent between parents and progenies
#'
#' Computes, for every family in a family history (or pedigree), the pairwise
#' block similarity between each parental haplotype and every progeny haplotype,
#' and flags identity-by-descent when the similarity meets a threshold.
#'
#' The comparison is anchored on the parental haplotypes: for each of a family's
#' \code{n_parents x ploidy} parental haplotypes (a family may have more than two
#' parents — e.g. a double-cross family with 4 parents has up to 8), the block
#' segment is compared against every progeny haplotype block. All pairwise
#' comparisons are kept (no best-match collapsing).
#'
#' @param geno_meta A valid \code{GenoMeta} list.
#' @param family_history Optional family-history list (named \code{gen_*}, each a
#'   \code{data.table(FID, IID, Label)}). Defaults to
#'   \code{geno_meta$metadata$family_history}.
#' @param pedigree Optional flat pedigree \code{data.frame(Sample, Parent1,
#'   Parent2)} used when no family history is available (groups offspring by
#'   their parent set; supports up to two parents).
#' @param blocks Optional character vector of block IDs to restrict the analysis
#'   (same vocabulary as \code{filter_geno_meta(blocks = )}).
#' @param threshold IBD similarity threshold, inclusive. Default \code{0.5}
#'   (\code{IBD = Similarity >= threshold}).
#' @param verbose Logical; print progress (families, haplotype counts, per-family
#'   IBD summaries) via \code{message()}. Default \code{FALSE}.
#'
#' @return A \code{data.table} with one row per pairwise comparison:
#'   \code{Gen}, \code{FID}, \code{Parent}, \code{ParentHap}, \code{Progeny},
#'   \code{ProgenyHap}, \code{Block}, \code{Chrom}, \code{Similarity},
#'   \code{IBD}.
#'
#' @seealso \code{\link{compute_segment_ibs}}, \code{\link{score_marker_signal}},
#'   \code{\link{trace_family_history}}
#' @export
compute_block_ibd <- function(geno_meta, family_history = NULL,
                              pedigree = NULL, blocks = NULL, threshold = 0.5,
                              verbose = FALSE) {
  .hs_check_meta(geno_meta, "compute_block_ibd", "geno_meta")
  if (is.null(family_history)) family_history <- geno_meta$metadata$family_history
  if (is.null(family_history) && is.null(pedigree))
    .hs_abort("compute_block_ibd: supply 'family_history' (or attach one via ",
         "attach_family_history()) or a 'pedigree'.")

  .sim_msg(verbose, "[compute_block_ibd] starting...")
  sim <- .compute_similarity_long(geno_meta, family_history, pedigree,
                                  blocks, threshold, verbose = verbose)

  block_tab <- unique(sim[, .(Gen, FID, Parent, ParentHap, Progeny, ProgenyHap,
                              Block, Chrom, Similarity, IBD)])
  data.table::setorder(block_tab, Gen, FID, Parent, ParentHap, Progeny,
                       ProgenyHap, Block)
  .sim_msg(verbose,
           "[compute_block_ibd] %d pairwise block comparisons; %.1f%% IBD.",
           nrow(block_tab), 100 * mean(block_tab$IBD, na.rm = TRUE))
  block_tab[]
}


# ---------------------------------------------------------------------------
# R2.3 — marker-level signal test
# ---------------------------------------------------------------------------

#' Marker-level inheritance-signal test
#'
#' Scores how reliably each marker reports the IBD status of its block. Block
#' IBD (per pairwise parent-haplotype × progeny-haplotype comparison) is the
#' ground truth; per-marker allele identity (parent vs progeny on the haplotypes
#' of that comparison) is the prediction. Builds a TP/FP/TN/FN confusion matrix
#' (standard orientation) summed over all comparisons.
#'
#' \describe{
#'   \item{\code{mode = "IBD"}}{Boolean / conditional. Each comparison×marker
#'     contributes integer 1 to one cell, gated by the boolean \code{IBD} flag:
#'     \code{same & IBD = TP}, \code{same & !IBD = FP}, \code{diff & IBD = FN},
#'     \code{diff & !IBD = TN}.}
#'   \item{\code{mode = "IBS"}}{Numeric / similarity-weighted. The same cells
#'     receive fractional weights: \code{TP += same × Similarity},
#'     \code{FP += same × (1 - Similarity)}, \code{FN += diff × Similarity},
#'     \code{TN += diff × (1 - Similarity)}.}
#' }
#'
#' @inheritParams compute_block_ibd
#' @param mode Measurement mode: \code{"IBD"} (boolean counts, default) or
#'   \code{"IBS"} (similarity-weighted).
#' @param n_signal Number of signal marker(s) to select per block once scoring
#'   finishes (see \code{\link{select_signal_markers}}). Default \code{1}; set
#'   \code{0} to skip the selection.
#' @param select_score Metric used to rank markers within a block when selecting.
#'   Default \code{"final_mcc"}.
#'
#' @return A named \code{list}:
#' \describe{
#'   \item{\code{detail}}{per (Gen, FID, IID_parent, ParentHap, Block, Marker)
#'     TP/FP/TN/FN (summed over that parental haplotype's progeny comparisons).}
#'   \item{\code{marker}}{aggregated per marker: \code{Chrom, Marker, TP, FP,
#'     TN, FN}.}
#'   \item{\code{metrics}}{the \code{marker} table passed through
#'     \code{\link{compute_signal_metrics}} (metrics + MAF-weighted final
#'     scores).}
#'   \item{\code{signal}}{the selected signal marker(s) per block — the
#'     \code{\link{select_signal_markers}} table. Attach it to the container with
#'     \code{\link{attach_signal_markers}}.}
#' }
#'
#' @seealso \code{\link{compute_block_ibd}}, \code{\link{compute_signal_metrics}},
#'   \code{\link{select_signal_markers}}, \code{\link{attach_signal_markers}}
#' @export
score_marker_signal <- function(geno_meta, family_history = NULL,
                                pedigree = NULL, blocks = NULL,
                                threshold = 0.5, mode = c("IBD", "IBS"),
                                n_signal = 1L, select_score = "final_mcc",
                                verbose = FALSE) {
  mode <- match.arg(mode)
  .hs_check_meta(geno_meta, "score_marker_signal", "geno_meta")
  if (is.null(family_history)) family_history <- geno_meta$metadata$family_history
  if (is.null(family_history) && is.null(pedigree))
    .hs_abort("score_marker_signal: supply 'family_history' (or attach one via ",
         "attach_family_history()) or a 'pedigree'.")

  .sim_msg(verbose, "[score_marker_signal] mode = %s; building comparisons...",
           mode)
  sim <- .compute_similarity_long(geno_meta, family_history, pedigree,
                                  blocks, threshold, verbose = verbose)

  empty_marker <- data.table(Chrom = character(0L), Marker = character(0L),
                             TP = numeric(0L), FP = numeric(0L),
                             TN = numeric(0L), FN = numeric(0L))
  if (nrow(sim) == 0L) {
    empty_metrics <- compute_signal_metrics(empty_marker, geno_meta = geno_meta)
    return(list(
      detail  = data.table(Gen = character(0L), FID = character(0L),
                           IID_parent = character(0L), ParentHap = integer(0L),
                           Block = character(0L), Chrom = character(0L),
                           Marker = character(0L), TP = numeric(0L),
                           FP = numeric(0L), TN = numeric(0L), FN = numeric(0L)),
      marker  = empty_marker,
      metrics = empty_metrics,
      signal  = data.table(Block = character(0L), Chrom = character(0L),
                           Marker = character(0L), Pos = numeric(0L),
                           score_val = numeric(0L),
                           signal_rank = integer(0L))))
  }

  # Only comparable marker rows contribute to the confusion matrix.
  d <- sim[comparable == TRUE]
  .sim_msg(verbose,
           "[score_marker_signal] %d comparable comparison-marker rows; tallying %s confusion matrix...",
           nrow(d), if (mode == "IBD") "integer" else "similarity-weighted")

  if (mode == "IBD") {
    d[, `:=`(
      TP = as.numeric(marker_same & IBD),
      FP = as.numeric(marker_same & !IBD),
      FN = as.numeric(!marker_same & IBD),
      TN = as.numeric(!marker_same & !IBD))]
  } else {
    d[, `:=`(
      TP = ifelse(marker_same,  Similarity,       0),
      FP = ifelse(marker_same,  1 - Similarity,   0),
      FN = ifelse(!marker_same, Similarity,       0),
      TN = ifelse(!marker_same, 1 - Similarity,   0))]
  }

  detail <- d[, .(TP = sum(TP), FP = sum(FP), TN = sum(TN), FN = sum(FN)),
              by = .(Gen, FID, IID_parent = Parent, ParentHap, Block, Chrom,
                     Marker, Pos)]
  data.table::setorder(detail, Chrom, Pos, IID_parent, ParentHap)
  detail[, Pos := NULL]

  marker <- d[, .(TP = sum(TP), FP = sum(FP), TN = sum(TN), FN = sum(FN)),
              by = .(Chrom, Marker, Pos)]
  data.table::setorder(marker, Chrom, Pos)
  marker[, Pos := NULL]

  metrics <- compute_signal_metrics(marker, geno_meta = geno_meta,
                                    verbose = verbose)

  # Selection: the block -> representative marker(s) decision. Uses the block
  # membership already carried by `detail`, so no extra map lookup is needed.
  n_signal <- suppressWarnings(as.integer(n_signal))
  if (is.na(n_signal) || n_signal < 0L)
    .hs_abort("score_marker_signal: 'n_signal' must be a non-negative integer.")
  signal <- if (n_signal > 0L) {
    select_signal_markers(
      metrics, geno_meta = geno_meta,
      block_marker_map = unique(detail[, .(Marker, Block)]),
      n = n_signal, score = select_score)
  } else {
    select_signal_markers(metrics[0L], geno_meta = geno_meta,
                          block_marker_map = unique(detail[, .(Marker, Block)]),
                          n = 1L, score = select_score)
  }
  .sim_msg(verbose,
           "[score_marker_signal] done: %d markers scored (detail rows: %d); %d signal marker(s) selected across %d block(s).",
           nrow(marker), nrow(detail), nrow(signal),
           length(unique(signal$Block)))

  list(detail = detail[], marker = marker[], metrics = metrics,
       signal = signal)
}


# ---------------------------------------------------------------------------
# R2.4 — classification metrics + MAF final score
# ---------------------------------------------------------------------------

#' Classification metrics and MAF-weighted final scores
#'
#' Turns a per-marker TP/FP/TN/FN confusion-matrix table into accuracy,
#' precision, recall (sensitivity), specificity, F1 and MCC, then multiplies
#' every metric by \code{MAF * (1 - MAF)} to produce \code{final_*} columns.
#' Undefined metrics (zero denominator, never \code{Inf}/\code{NaN}) are filled
#' with \code{zero_division} (default \code{0}); markers with no computable MAF
#' get \code{MAF = 0} (so their \code{final_*} scores become 0).
#'
#' @param marker_table A \code{data.frame}/\code{data.table} with \code{TP},
#'   \code{FP}, \code{TN}, \code{FN} (and usually \code{Chrom}, \code{Marker}).
#' @param geno_meta Optional \code{GenoMeta} used to compute per-marker
#'   genome-wide minor allele frequency. Required unless \code{maf} is supplied.
#' @param maf Optional named numeric vector of per-marker MAF (names = marker
#'   IDs). Overrides \code{geno_meta}.
#' @param zero_division Value used for metrics whose denominator is zero (e.g.
#'   \code{precision} when a marker is never called "same"). Default \code{0}
#'   (scikit-learn convention); set \code{NA} to restore the old NA behaviour.
#' @param verbose Logical; print a one-line progress summary via
#'   \code{message()}. Default \code{FALSE}.
#'
#' @return The input table (as a \code{data.table}) with added columns
#'   \code{accuracy}, \code{precision}, \code{recall}, \code{specificity},
#'   \code{f1}, \code{mcc}, \code{MAF}, and \code{final_*} for each metric.
#'
#' @seealso \code{\link{score_marker_signal}}
#' @export
compute_signal_metrics <- function(marker_table, geno_meta = NULL, maf = NULL,
                                   zero_division = 0, verbose = FALSE) {
  dt <- as.data.table(data.table::copy(marker_table))
  if (!all(c("TP", "FP", "TN", "FN") %in% names(dt)))
    .hs_abort("compute_signal_metrics: 'marker_table' needs TP, FP, TN, FN columns.")
  if (is.null(maf) && is.null(geno_meta))
    .hs_abort("compute_signal_metrics: supply 'geno_meta' or a named 'maf' vector.")

  zero_division <- as.numeric(zero_division)
  fill <- function(x) { x[!is.finite(x)] <- zero_division; x }
  metric_cols <- c("accuracy", "precision", "recall", "specificity", "f1", "mcc")

  if (nrow(dt) == 0L) {
    for (m in c(metric_cols, "MAF", paste0("final_", metric_cols)))
      dt[, (m) := numeric(0L)]
    return(dt[])
  }

  TP <- as.numeric(dt$TP); FP <- as.numeric(dt$FP)
  TN <- as.numeric(dt$TN); FN <- as.numeric(dt$FN)

  dt[, accuracy    := fill((TP + TN) / (TP + TN + FP + FN))]
  dt[, precision   := fill(TP / (TP + FP))]
  dt[, recall      := fill(TP / (TP + FN))]
  dt[, specificity := fill(TN / (TN + FP))]
  # f1 from the already-filled precision/recall; a residual 0/0 -> zero_division.
  dt[, f1          := fill(2 * precision * recall / (precision + recall))]
  mcc_den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  dt[, mcc         := fill((TP * TN - FP * FN) / mcc_den)]

  if (is.null(maf)) maf <- .marker_maf(geno_meta)
  if (is.null(names(maf)))
    .hs_abort("compute_signal_metrics: 'maf' must be a named vector (marker IDs).")
  if (!"Marker" %in% names(dt))
    .hs_abort("compute_signal_metrics: 'marker_table' needs a 'Marker' column to ",
         "join MAF.")
  dt[, MAF := unname(maf[as.character(Marker)])]
  n_no_maf <- sum(is.na(dt$MAF))
  if (n_no_maf > 0L)
    .hs_warn("compute_signal_metrics: ", n_no_maf,
            " marker(s) have no MAF; set to 0 (their final_* scores become 0).")
  dt[is.na(MAF), MAF := 0]

  w <- dt$MAF * (1 - dt$MAF)
  for (m in metric_cols)
    dt[, (paste0("final_", m)) := fill(get(m) * w)]

  .sim_msg(verbose,
           "[compute_signal_metrics] %d markers; MAF source: %s; %d marker(s) without MAF.",
           nrow(dt), if (is.null(geno_meta)) "supplied vector" else "geno_meta",
           n_no_maf)
  dt[]
}


# ---------------------------------------------------------------------------
# R2.5 — signal-marker selection (block -> representative marker[s])
# ---------------------------------------------------------------------------

# Internal: the block_info columns written by attach_signal_markers(). Kept in
# one place because filter_geno_meta(), attach_block_info() and
# build_block_location_info() all have to recognise them.
.signal_marker_cols <- function() c("Signal_Marker", "Signal_Score", "N_Signal")

# Internal: separator used to pack >1 selected marker into one block_info cell.
.signal_sep <- ","


#' Select the signal marker(s) of each haplotype block
#'
#' Turns the per-marker metrics of \code{\link{score_marker_signal}} into a
#' decision: the \code{n} marker(s) that best tag each block's inheritance. These
#' are the markers worth genotyping (or plotting) to follow the block through a
#' pedigree.
#'
#' Within each block, markers are ranked by \code{score} (descending) and the top
#' \code{n} are kept. Blocks with fewer scored markers yield fewer rows. Ties are
#' broken deterministically — see \code{ties} — so repeated runs select the same
#' markers.
#'
#' @param metrics A \code{data.frame}/\code{data.table} of per-marker metrics —
#'   normally the \code{metrics} element of \code{\link{score_marker_signal}}.
#'   Needs \code{Marker} and the \code{score} column; \code{Chrom}/\code{Pos} are
#'   filled from \code{geno_meta} when absent.
#' @param geno_meta Optional \code{GenoMeta} supplying the marker→block map
#'   (\code{block_marker_map_info}) and \code{map_info} positions.
#' @param block_marker_map Optional \code{data.frame}/\code{data.table} with
#'   \code{Marker} and \code{Block}; overrides \code{geno_meta}.
#' @param n Number of signal markers to select per block. Default \code{1}.
#' @param score Metric used to rank markers within a block. Default
#'   \code{"final_mcc"} (MAF-weighted MCC).
#' @param min_score Optional lower bound: markers scoring below it are never
#'   selected, so a block with no informative marker is left unselected rather
#'   than being given a weak representative. Default \code{NULL} (no floor).
#' @param ties Tie-break rule for markers with an identical \code{score}:
#'   \code{"position"} (default, lowest \code{Pos}), \code{"maf"} (highest
#'   \code{MAF}; falls back to position when \code{metrics} has no \code{MAF}),
#'   or \code{"first"} (row order of \code{metrics}).
#'
#' @return A \code{data.table} with one row per selected marker:
#'   \code{Block}, \code{Chrom}, \code{Marker}, \code{Pos}, \code{score_val},
#'   \code{signal_rank} (1 = best in its block). The \code{score} name and
#'   \code{n} are kept as attributes (\code{attr(x, "score")}, \code{attr(x, "n")}).
#'
#' @examples
#' \dontrun{
#' sig <- score_marker_signal(meta, mode = "IBD")
#' sel <- select_signal_markers(sig$metrics, geno_meta = meta, n = 1)
#' meta <- attach_signal_markers(meta, sel)
#' }
#'
#' @seealso \code{\link{score_marker_signal}}, \code{\link{attach_signal_markers}},
#'   \code{\link{get_signal_markers}}, \code{\link{plot_marker_signal}}
#' @export
select_signal_markers <- function(metrics, geno_meta = NULL,
                                  block_marker_map = NULL, n = 1L,
                                  score = "final_mcc", min_score = NULL,
                                  ties = c("position", "maf", "first")) {
  ties <- match.arg(ties)
  n    <- suppressWarnings(as.integer(n))
  if (is.na(n) || n < 1L)
    .hs_abort("select_signal_markers: 'n' must be a positive integer.")

  dt <- as.data.table(data.table::copy(metrics))
  if (!"Marker" %in% names(dt))
    .hs_abort("select_signal_markers: 'metrics' needs a 'Marker' column.")
  if (!score %in% names(dt))
    .hs_abort("select_signal_markers: score column '", score,
         "' not found in 'metrics'.")
  dt[, Marker := as.character(Marker)]
  dt[, score_val := as.numeric(get(score))]

  empty <- data.table(Block = character(0L), Chrom = character(0L),
                      Marker = character(0L), Pos = numeric(0L),
                      score_val = numeric(0L), signal_rank = integer(0L))
  data.table::setattr(empty, "score", score)
  data.table::setattr(empty, "n", n)

  # ── marker -> block map (argument wins over the container) ──────────────────
  bmm <- block_marker_map %||%
    (if (!is.null(geno_meta)) geno_meta$block_marker_map_info else NULL)
  if (is.null(bmm))
    .hs_abort("select_signal_markers: no marker->block map; supply ",
         "'block_marker_map' or a 'geno_meta' carrying block_marker_map_info ",
         "(see attach_block_info()).")
  bmm <- as.data.table(bmm)
  if (!all(c("Marker", "Block") %in% names(bmm)))
    .hs_abort("select_signal_markers: 'block_marker_map' needs Marker and Block ",
         "columns.")
  bmm <- bmm[, .(Marker = as.character(Marker), Block = as.character(Block))]
  if ("Block" %in% names(dt)) dt[, Block := NULL]
  dt[bmm, on = "Marker", Block := i.Block]

  # ── Chrom / Pos: from `metrics`, else from the map ──────────────────────────
  if (!all(c("Chrom", "Pos") %in% names(dt)) && !is.null(geno_meta)) {
    mi <- as.data.table(geno_meta$map_info)[
      , .(Marker = as.character(Marker), Chrom = as.character(Chrom),
          Pos = as.numeric(Pos))]
    if (!"Chrom" %in% names(dt)) dt[mi, on = "Marker", Chrom := i.Chrom]
    if (!"Pos"   %in% names(dt)) dt[mi, on = "Marker", Pos   := i.Pos]
  }
  if (!"Chrom" %in% names(dt)) dt[, Chrom := NA_character_]
  if (!"Pos"   %in% names(dt)) dt[, Pos   := NA_real_]

  sel <- dt[!is.na(Block) & is.finite(score_val)]
  if (!is.null(min_score)) sel <- sel[score_val >= as.numeric(min_score)]
  if (nrow(sel) == 0L) return(empty)

  # ── deterministic tie-break ────────────────────────────────────────────────
  sel[, tie_key := switch(
    ties,
    position = as.numeric(Pos),
    maf      = if ("MAF" %in% names(sel)) -as.numeric(MAF) else as.numeric(Pos),
    first    = seq_len(.N))]
  if (ties == "maf" && !"MAF" %in% names(sel))
    .hs_warn("select_signal_markers: ties = 'maf' but 'metrics' has no MAF ",
            "column; falling back to position.")
  data.table::setorderv(sel, c("Block", "score_val", "tie_key", "Marker"),
                        c(1L, -1L, 1L, 1L), na.last = TRUE)

  out <- sel[, utils::head(.SD, n), by = Block,
             .SDcols = c("Chrom", "Marker", "Pos", "score_val")]
  out[, signal_rank := seq_len(.N), by = Block]
  out <- out[, .(Block, Chrom = as.character(Chrom), Marker, Pos = as.numeric(Pos),
                 score_val, signal_rank = as.integer(signal_rank))]
  data.table::setorderv(out, c("Chrom", "Pos", "signal_rank"), na.last = TRUE)
  data.table::setattr(out, "score", score)
  data.table::setattr(out, "n", n)
  out[]
}


#' Record selected signal markers on a GenoMeta
#'
#' Writes a \code{\link{select_signal_markers}} selection into the
#' \code{block_info} table of a \code{GenoMeta} — the container's block-level
#' sheet — so the decision travels with the data instead of living in a
#' side file.
#'
#' Three columns are added (or overwritten): \code{Signal_Marker} (the selected
#' marker IDs, comma-separated in rank order), \code{Signal_Score} (the rank-1
#' marker's score — the block's signal strength) and \code{N_Signal} (how many
#' markers were selected; may be smaller than the requested \code{n} for short
#' blocks). Blocks with no selection get \code{NA} / \code{0}. The selection
#' parameters are recorded in \code{metadata$signal_selection}.
#'
#' \strong{Lifecycle.} \code{\link{filter_geno_meta}} prunes selections whose
#' marker no longer survives the filter; \code{\link{attach_block_info}} drops
#' the columns (with a message) because redefining blocks invalidates any
#' selection made against the old ones. \code{\link{save_geno_meta}} /
#' \code{\link{load_geno_meta}} round-trip them; a CSV block file does not.
#'
#' @param geno_meta A valid \code{GenoMeta} carrying block information. When it
#'   has \code{block_marker_map_info} but no \code{block_info}, the coordinate
#'   table is derived on the fly (same computation as
#'   \code{\link{build_block_location_info}}).
#' @param signal A selection table from \code{\link{select_signal_markers}} (or
#'   the \code{signal} element of \code{\link{score_marker_signal}}): needs
#'   \code{Block}, \code{Marker} and, for \code{Signal_Score}, \code{score_val}
#'   and \code{signal_rank}.
#' @param overwrite Logical; replace an existing selection. When \code{FALSE} and
#'   the columns are already present, the call errors. Default \code{TRUE}.
#'
#' @return The \code{GenoMeta} with the signal columns on \code{block_info}.
#'
#' @seealso \code{\link{select_signal_markers}}, \code{\link{get_signal_markers}}
#' @export
attach_signal_markers <- function(geno_meta, signal, overwrite = TRUE) {
  .hs_check_meta(geno_meta, "attach_signal_markers", "geno_meta")

  # block_info is optional on a GenoMeta (a block map alone is enough for most
  # of the package), so derive the coordinate table when it is missing rather
  # than making the caller run attach_block_info() first.
  if (is.null(geno_meta$block_info)) {
    if (is.null(geno_meta$block_marker_map_info))
      .hs_abort("attach_signal_markers: 'geno_meta' has no block information; call ",
           "attach_block_info() first.")
    .hs_inform("attach_signal_markers: no block_info on the GenoMeta \u2014 derived it ",
            "from block_marker_map_info.")
    geno_meta$block_info <- as.data.frame(
      .compute_block_location(geno_meta$block_marker_map_info,
                              geno_meta$map_info)$block_info,
      stringsAsFactors = FALSE)
  }

  sig <- as.data.table(data.table::copy(signal))
  if (!all(c("Block", "Marker") %in% names(sig)))
    .hs_abort("attach_signal_markers: 'signal' needs Block and Marker columns ",
         "(see select_signal_markers()).")
  if (!"signal_rank" %in% names(sig)) sig[, signal_rank := seq_len(.N), by = Block]
  if (!"score_val"   %in% names(sig)) sig[, score_val := NA_real_]

  bi   <- as.data.frame(geno_meta$block_info, stringsAsFactors = FALSE)
  cols <- .signal_marker_cols()
  if (!isTRUE(overwrite) && any(cols %in% names(bi)))
    .hs_abort("attach_signal_markers: 'geno_meta' already carries a signal-marker ",
         "selection and overwrite = FALSE.")

  unknown <- setdiff(unique(sig$Block), as.character(bi$Block))
  if (length(unknown) > 0L) {
    .hs_warn("attach_signal_markers: ", length(unknown),
            " block(s) in 'signal' absent from block_info \u2014 dropped: ",
            paste(utils::head(unknown, 5L), collapse = ", "),
            if (length(unknown) > 5L) " ..." else "")
    sig <- sig[Block %in% as.character(bi$Block)]
  }

  data.table::setorder(sig, Block, signal_rank)
  packed <- sig[, .(Signal_Marker = paste(as.character(Marker),
                                          collapse = .signal_sep),
                    Signal_Score  = as.numeric(score_val[1L]),
                    N_Signal      = as.integer(.N)),
                by = Block]

  idx <- match(as.character(bi$Block), packed$Block)
  bi$Signal_Marker <- packed$Signal_Marker[idx]
  bi$Signal_Score  <- packed$Signal_Score[idx]
  bi$N_Signal      <- ifelse(is.na(idx), 0L, packed$N_Signal[idx])

  # The block definitions are unchanged (only Signal_* columns were added), so
  # everything else is forwarded (one rebuild path, 0.9.0). No timestamp in
  # the record: containers built from the same inputs stay identical().
  .rebuild_geno_meta(
    geno_meta,
    block_info = bi,
    metadata   = list(signal_selection = list(
      score    = attr(signal, "score") %||% NA_character_,
      n        = attr(signal, "n") %||% max(sig$signal_rank, 0L),
      n_blocks = nrow(packed))),
    .caller = "attach_signal_markers")
}


#' Read the signal-marker selection stored on a GenoMeta
#'
#' The single parser of the packed \code{block_info$Signal_Marker} column — every
#' consumer (plots, exports, downstream filters) should go through here rather
#' than splitting the string itself.
#'
#' @param geno_meta A \code{GenoMeta}, or any \code{block_info}-shaped
#'   \code{data.frame} carrying \code{Block} and \code{Signal_Marker}.
#' @param long Logical. \code{TRUE} (default) returns a tidy
#'   \code{data.table(Block, Marker, signal_rank, Signal_Score)}; \code{FALSE}
#'   returns just the character vector of selected marker IDs.
#'
#' @return A \code{data.table} (\code{long = TRUE}) or a character vector.
#'   Empty / \code{character(0)} when no selection is stored.
#'
#' @seealso \code{\link{select_signal_markers}}, \code{\link{attach_signal_markers}}
#' @export
get_signal_markers <- function(geno_meta, long = TRUE) {
  bi <- if (is.data.frame(geno_meta)) geno_meta else geno_meta$block_info
  empty <- data.table(Block = character(0L), Marker = character(0L),
                      signal_rank = integer(0L), Signal_Score = numeric(0L))

  if (is.null(bi) || !all(c("Block", "Signal_Marker") %in% names(bi)))
    return(if (long) empty else character(0L))

  bi <- as.data.table(bi)[!is.na(Signal_Marker) & nzchar(Signal_Marker)]
  if (nrow(bi) == 0L) return(if (long) empty else character(0L))

  parts <- strsplit(as.character(bi$Signal_Marker), .signal_sep, fixed = TRUE)
  out <- data.table(
    Block        = rep(as.character(bi$Block), lengths(parts)),
    Marker       = trimws(unlist(parts, use.names = FALSE)),
    Signal_Score = rep(if ("Signal_Score" %in% names(bi))
                         as.numeric(bi$Signal_Score) else NA_real_,
                       lengths(parts)))
  out[, signal_rank := seq_len(.N), by = Block]
  if (!long) return(out$Marker)
  out[, .(Block, Marker, signal_rank = as.integer(signal_rank), Signal_Score)][]
}


# Internal: drop stored selections whose marker no longer exists in `markers`
# (e.g. after filter_geno_meta() prunes the panel), recomputing the packed
# column and N_Signal. Returns the block_info data.frame unchanged when it
# carries no selection. `caller` only shapes the message.
.prune_signal_markers <- function(block_info, markers,
                                  caller = "filter_geno_meta") {
  if (is.null(block_info) || !"Signal_Marker" %in% names(block_info))
    return(block_info)

  sel <- get_signal_markers(block_info, long = TRUE)
  if (nrow(sel) == 0L) return(block_info)

  gone <- setdiff(sel$Marker, as.character(markers))
  if (length(gone) == 0L) return(block_info)

  kept <- sel[!Marker %in% gone]
  .hs_inform(caller, ": ", length(gone),
          " signal marker(s) no longer in the panel \u2014 dropped from block_info.")

  if (nrow(kept) == 0L) {
    block_info$Signal_Marker <- NA_character_
    block_info$Signal_Score  <- NA_real_
    block_info$N_Signal      <- 0L
    return(block_info)
  }

  data.table::setorder(kept, Block, signal_rank)
  packed <- kept[, .(Signal_Marker = paste(Marker, collapse = .signal_sep),
                     N_Signal      = as.integer(.N)), by = Block]
  idx <- match(as.character(block_info$Block), packed$Block)
  block_info$Signal_Marker <- packed$Signal_Marker[idx]
  block_info$N_Signal      <- ifelse(is.na(idx), 0L, packed$N_Signal[idx])
  # The rank-1 marker may itself have gone; the stored score then no longer
  # describes the surviving top marker, so it is cleared rather than left stale.
  if ("Signal_Score" %in% names(block_info)) {
    lost_top <- unique(sel[signal_rank == 1L & Marker %in% gone]$Block)
    block_info$Signal_Score[as.character(block_info$Block) %in% lost_top] <- NA_real_
    block_info$Signal_Score[is.na(idx)] <- NA_real_
  }
  block_info
}


# Internal: per-marker genome-wide minor allele frequency (biallelic:
# 1 - max(allele freq)). Returns a named numeric vector over all markers.
.marker_maf <- function(geno_meta) {
  arr     <- geno_meta$geno_arrays
  markers <- dimnames(arr)[[2L]]
  enc <- geno_meta$metadata$encoding %||%
    detect_genotype_encoding(.meta_to_allele_matrix(geno_meta))
  miss <- .missing_tokens(enc, .meta_extra_missing(geno_meta))

  # Per-marker allele counts (missing tokens resolved in R); MAF = 1 - major
  # allele frequency. Markers with no observed allele stay NA.
  counts <- .marker_allele_counts(arr, na_vals = miss)
  tot    <- rowSums(counts)
  mx     <- if (ncol(counts) > 0L) apply(counts, 1L, max) else integer(nrow(counts))
  maf    <- ifelse(tot > 0L, 1 - mx / tot, NA_real_)
  stats::setNames(maf[match(markers, rownames(counts))], markers)
}
