# =============================================================================
# marker_stats.R
#
# R-side allele/DICT detection + dispatch for the per-marker allele-count and
# per-cell dosage kernels (src/marker_stats.cpp, Kernel B). This layer owns ALL
# allele semantics; the compiled kernels only count integer codes and match a
# target code (see the header of src/marker_stats.cpp).
#
# Design contract (mirrors src/local_gebv.cpp): every decision about which string
# values are missing/blank, the GLOBAL allele dictionary (factor levels), and
# each marker's target allele (minor/major/alt) is resolved HERE, in R. The
# genotype array is reduced to 1-based integer codes (NA = missing) before any
# kernel is called, so the kernels stay encoding-agnostic and identical results
# are produced by the Rcpp and pure-R paths (parity-tested in test-marker-stats.R).
# =============================================================================


# ---------------------------------------------------------------------------
# Global integer coding of a 3-D character genotype array (ALL detection here)
# ---------------------------------------------------------------------------

# Encode a 3-D character genotype array [n_samples x n_markers x ploidy] to
# 1-based integer codes against a GLOBAL, alphabetically-sorted allele dictionary
# built here in R. Every missingness decision happens here:
#   * NA cells                          -> NA
#   * cells whose value is in `na_vals` -> NA
#   * (optional) blank / whitespace-only cells, when strip_blank = TRUE -> NA
# Returns list(codes = <integer array, same dims>, levels = <character dict>).
.encode_marker_codes <- function(arr, na_vals = character(0L),
                                  strip_blank = FALSE) {
  flat   <- as.character(arr)
  na_hit <- is.na(flat) | flat %in% na_vals
  if (isTRUE(strip_blank)) na_hit <- na_hit | !nzchar(trimws(flat))
  flat[na_hit] <- NA_character_
  f <- factor(flat)                       # levels() are sorted -> the dictionary
  list(codes  = array(as.integer(f), dim = dim(arr)),
       levels = levels(f))
}


# ---------------------------------------------------------------------------
# Per-marker allele counts  [n_markers x n_levels]
# ---------------------------------------------------------------------------

# Pure-R fallback: tabulate integer codes per marker across all ploidy layers.
.marker_allele_counts_r <- function(codes, n_mark, n_lev) {
  out <- matrix(0L, n_mark, n_lev)
  for (j in seq_len(n_mark)) {
    v <- codes[, j, ]
    v <- v[!is.na(v)]
    if (length(v)) out[j, ] <- tabulate(v, nbins = n_lev)
  }
  out
}

# Count each dictionary allele per marker. Returns an [n_markers x n_levels]
# integer matrix with rownames = markers and colnames = the allele dictionary.
# Dispatches to marker_allele_counts_cpp or the pure-R fallback; identical result.
.marker_allele_counts <- function(arr, na_vals = character(0L),
                                   strip_blank = FALSE,
                                   use_rcpp = .hs_use_rcpp("marker_allele_counts_cpp"),
                                   enc = NULL) {
  # `enc`: a pre-built .encode_marker_codes() result, so a caller that already
  # encoded the array (compute_marker_stats(), cc/F19) does not encode twice.
  if (is.null(enc)) enc <- .encode_marker_codes(arr, na_vals, strip_blank)
  markers <- dimnames(arr)[[2L]]
  n_mark  <- dim(arr)[2L]
  n_lev   <- length(enc$levels)
  if (n_lev == 0L)
    return(matrix(0L, n_mark, 0L, dimnames = list(markers, NULL)))

  counts <- if (isTRUE(use_rcpp)) {
    .hs_kernel("marker_allele_counts_cpp")(as.integer(enc$codes),
                                           as.integer(dim(arr)), n_lev)
  } else {
    .marker_allele_counts_r(enc$codes, n_mark, n_lev)
  }
  dimnames(counts) <- list(markers, enc$levels)
  counts
}

# Major / minor allele per marker from an [n_markers x n_levels] count matrix
# (levels in dictionary = alphabetical order). The one tie rule (0.9.0,
# decision-form row 16 a, D22):
#   Minor = the least frequent observed allele; ties -> first in allele order
#           (the same allele the dosage counts, so MAF and dosage agree)
#   Major = the most frequent among the REMAINING observed alleles; ties ->
#           first in allele order
# so on an exact two-way tie Major != Minor (0.8.8 reported the same allele for
# both). A monomorphic marker has Major only; no observed allele gives NA/NA.
.hs_major_minor <- function(counts) {
  n_mark <- nrow(counts)
  levs   <- colnames(counts)
  major  <- rep(NA_character_, n_mark)
  minor  <- rep(NA_character_, n_mark)
  if (length(levs) == 0L || n_mark == 0L)
    return(list(Major = major, Minor = minor))
  for (j in seq_len(n_mark)) {
    cj  <- counts[j, ]
    obs <- which(cj > 0L)
    if (length(obs) == 0L) next
    if (length(obs) == 1L) { major[j] <- levs[obs]; next }
    i_min <- obs[which.min(cj[obs])]            # first in allele order on a tie
    rest  <- setdiff(obs, i_min)
    i_max <- rest[which.max(cj[rest])]          # first in allele order on a tie
    minor[j] <- levs[i_min]
    major[j] <- levs[i_max]
  }
  list(Major = major, Minor = minor)
}


# ---------------------------------------------------------------------------
# Per-cell dosage  [n_samples x n_markers]
# ---------------------------------------------------------------------------

# Pure-R fallback for the dosage kernel (see kernel header for the NA policy).
.dosage_matrix_r <- function(codes, tcode, na_if_any) {
  d      <- dim(codes)
  n_samp <- d[1L]; n_mark <- d[2L]; ploidy <- d[3L]
  out    <- matrix(0L, n_samp, n_mark)
  for (j in seq_len(n_mark)) {
    slab   <- matrix(codes[, j, ], nrow = n_samp, ncol = ploidy)
    n_miss <- rowSums(is.na(slab))
    tj     <- tcode[j]
    res    <- if (is.na(tj)) integer(n_samp)
              else as.integer(rowSums(slab == tj, na.rm = TRUE))
    if (isTRUE(na_if_any)) res[n_miss > 0L] <- NA_integer_
    res[n_miss == ploidy]  <- NA_integer_          # fully missing
    out[, j] <- res
  }
  out
}

# Dosage = number of ploidy layers equal to each marker's target allele. All
# allele/missing resolution is done here; `target_alleles` is a length-n_markers
# character vector of target allele strings (NA = no target -> 0 where observed).
# `na_if_any` selects the missing-cell policy (see src/marker_stats.cpp header).
# Returns an [n_samples x n_markers] integer matrix (dimnames from the array).
.dosage_matrix <- function(arr, target_alleles, na_if_any = TRUE,
                           na_vals = character(0L), strip_blank = FALSE,
                           use_rcpp = .hs_use_rcpp("dosage_from_codes_cpp")) {
  enc     <- .encode_marker_codes(arr, na_vals, strip_blank)
  samples <- dimnames(arr)[[1L]]
  markers <- dimnames(arr)[[2L]]
  tcode   <- as.integer(match(as.character(target_alleles), enc$levels))

  dose <- if (isTRUE(use_rcpp)) {
    .hs_kernel("dosage_from_codes_cpp")(as.integer(enc$codes),
                                        as.integer(dim(arr)), tcode,
                                        isTRUE(na_if_any))
  } else {
    .dosage_matrix_r(enc$codes, tcode, isTRUE(na_if_any))
  }
  dimnames(dose) <- list(samples, markers)
  dose
}


# ---------------------------------------------------------------------------
# Exported: per-marker QC summary
# ---------------------------------------------------------------------------

#' Per-marker quality-control summary (MAF, heterozygosity, missing rate)
#'
#' Computes a tidy per-marker statistics table from a \code{GenoMeta} object:
#' minor-allele frequency, observed heterozygosity, missing / call rate, the
#' observed count, and the major / minor allele. All allele semantics (missing
#' tokens, the allele dictionary) are resolved from the GenoMeta encoding, so the
#' result is independent of how alleles happen to be stored.
#'
#' Statistics are derived from a single integer-encoding pass over
#' \code{geno_meta$geno_arrays} (\code{[n_samples x n_markers x ploidy]}) and reuse
#' the internal counting kernels (\code{.marker_allele_counts()},
#' \code{.marker_maf()}). Marker location and block membership are taken from
#' the GenoMeta map and block assignments.
#'
#' @param geno_meta A valid \code{GenoMeta} list (see \code{\link{new_geno_meta}}).
#' @param markers   Optional character vector of marker IDs to restrict the
#'   summary to. Default \code{NULL} (all markers).
#' @param output_file Optional path. When non-\code{NULL}, the table is written as
#'   a tab-separated file via \code{data.table::fwrite()}; the table is returned
#'   either way. Default \code{NULL} (return only, nothing written).
#'
#' @return A \code{data.table}, one row per marker in genomic order, with columns:
#' \describe{
#'   \item{\code{Marker}}{Marker ID.}
#'   \item{\code{Chrom}}{Chromosome.}
#'   \item{\code{Pos}}{Genomic position.}
#'   \item{\code{Block}}{Assigned block ID, or \code{NA} when the marker is not in
#'     a block (or no block assignments are attached).}
#'   \item{\code{N_obs}}{Number of non-missing allele observations (sample x
#'     ploidy cells) at the marker.}
#'   \item{\code{Missing_rate}}{Fraction of missing cells (over
#'     \code{n_samples x ploidy}).}
#'   \item{\code{Call_rate}}{\code{1 - Missing_rate}.}
#'   \item{\code{MAF}}{Minor-allele frequency (\code{1 - major-allele frequency});
#'     \code{NA} for a marker with no observed allele.}
#'   \item{\code{Major}}{Most frequent observed allele (as stored); \code{NA} when
#'     nothing is observed.}
#'   \item{\code{Minor}}{Least frequent observed allele; \code{NA} for a
#'     monomorphic marker (only one observed allele) or none observed.}
#'   \item{\code{Het}}{Observed heterozygosity — the fraction of scored samples
#'     whose two haplotype copies differ. \code{NA} for single-copy (dosage /
#'     \code{"numeric"}) encodings, which carry no phase.}
#' }
#'
#' @examples
#' \dontrun{
#' meta  <- load_geno_meta("my_geno_meta.rds")
#' stats <- compute_marker_stats(meta, output_file = "marker_stats.txt")
#' head(stats)
#' }
#'
#' @seealso \code{\link{geno_meta_views}}, \code{\link{write_block_allele_sheet}}
#' @section Major, minor and ploidy (0.9.0):
#' One rule names the alleles (decision-form row 16 a): \code{Minor} is the
#' least frequent observed allele, ties broken by allele order (so it is the
#' allele the dosage counts and \code{MAF} agrees with it); \code{Major} is the
#' most frequent of the remaining alleles, ties broken the same way. On an exact
#' two-way tie the two therefore differ (0.8.8 reported the same allele for
#' both). \code{MAF} is \code{1 - } the major-allele frequency from the same
#' counts. \code{Het} is defined for diploid data only and is \code{NA} at any
#' other ploidy; a dosage-coded container (\code{metadata$encoding ==
#' "numeric"}) is read as one layer, so its \code{Missing_rate} is the true one.
#'
#' @export
compute_marker_stats <- function(geno_meta, markers = NULL, output_file = NULL) {
  .hs_check_meta(geno_meta, "compute_marker_stats", "geno_meta")

  arr <- geno_meta$geno_arrays
  if (!is.null(markers)) {
    keep <- dimnames(arr)[[2L]] %in% as.character(markers)
    if (!any(keep))
      .hs_abort("compute_marker_stats: none of 'markers' are present in geno_meta.")
    arr <- arr[, keep, , drop = FALSE]
  }
  marker_names <- dimnames(arr)[[2L]]
  d            <- dim(arr)
  n_samp <- d[1L]; n_mark <- d[2L]; ploidy <- d[3L]

  # Encoding + missing tokens (all missingness decided in R, as elsewhere).
  enc_str <- geno_meta$metadata$encoding %||%
    detect_genotype_encoding(.meta_to_allele_matrix(geno_meta))
  miss <- .missing_tokens(enc_str, .meta_extra_missing(geno_meta))

  # cc/F21: a dosage container ("numeric", the phantom second layer that P1
  # retires) carries one real layer, so its cell count is n_samp x 1 and the
  # missing rate no longer reads 0.5 on every marker.
  eff_ploidy <- if (identical(enc_str, "numeric")) 1L else ploidy
  if (eff_ploidy < ploidy) arr <- arr[, , seq_len(eff_ploidy), drop = FALSE]

  # ONE integer-encoding pass (cc/F19); every statistic below reads it.
  enc    <- .encode_marker_codes(arr, na_vals = miss)
  codes  <- enc$codes                                   # [n_samp x n_mark x eff_ploidy]
  is_na  <- is.na(codes)

  # Missing / call rate and observed count per marker.
  n_cells      <- n_samp * eff_ploidy
  miss_cells   <- apply(is_na, 2L, sum)
  missing_rate <- miss_cells / n_cells
  n_obs        <- n_cells - miss_cells

  # Observed heterozygosity: diploid only (cc/F21) — scored samples whose two
  # copies differ. NA for any other ploidy rather than a two-layer guess.
  if (eff_ploidy == 2L) {
    l1 <- codes[, , 1L]; l2 <- codes[, , 2L]      # [n_samp x n_mark]
    if (n_mark == 1L) { l1 <- matrix(l1, n_samp, 1L); l2 <- matrix(l2, n_samp, 1L) }
    both_obs <- !is.na(l1) & !is.na(l2)
    het_ij   <- both_obs & (l1 != l2)
    scored   <- colSums(both_obs)
    het      <- ifelse(scored > 0L, colSums(het_ij) / scored, NA_real_)
  } else {
    het <- rep(NA_real_, n_mark)
  }

  # Per-marker allele counts from the same encoding; MAF = 1 - major frequency
  # (the .marker_maf() rule, computed here from these counts, cc/F19); Major /
  # Minor by the one tie rule (.hs_major_minor(), row 16 a).
  counts <- .marker_allele_counts(arr, na_vals = miss, enc = enc)   # [n_mark x n_lev]
  tot    <- rowSums(counts)
  mx     <- if (ncol(counts) > 0L) apply(counts, 1L, max) else integer(nrow(counts))
  maf    <- ifelse(tot > 0L, 1 - mx / tot, NA_real_)
  mm_al  <- .hs_major_minor(counts)

  stats <- data.table::data.table(
    Marker       = marker_names,
    N_obs        = as.integer(n_obs),
    Missing_rate = as.numeric(missing_rate),
    Call_rate    = as.numeric(1 - missing_rate),
    MAF          = as.numeric(maf),
    Major        = mm_al$Major,
    Minor        = mm_al$Minor,
    Het          = as.numeric(het))

  # Location + block membership from the map and the block map directly
  # (cc/F19): geno_meta_views() would melt the whole allele table for the same
  # four fields.
  mm <- .marker_locus_table(geno_meta)[, .(Marker, Chrom, Pos, Block)]
  stats <- merge(mm, stats, by = "Marker", all.y = TRUE, sort = FALSE)
  # 0.9.0 (row 29 a): the one natural chromosome order, then position
  stats[, .chr_key := .hs_chrom_rank(Chrom)]
  data.table::setorder(stats, .chr_key, Pos)
  stats[, .chr_key := NULL]
  data.table::setcolorder(stats,
    c("Marker", "Chrom", "Pos", "Block", "N_obs", "Missing_rate",
      "Call_rate", "MAF", "Major", "Minor", "Het"))

  if (!is.null(output_file))
    data.table::fwrite(stats, output_file, sep = "\t")
  stats[]
}
