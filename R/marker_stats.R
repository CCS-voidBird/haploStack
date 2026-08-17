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
                                   use_rcpp = .hs_use_rcpp("marker_allele_counts_cpp")) {
  enc     <- .encode_marker_codes(arr, na_vals, strip_blank)
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
#' the internal counting kernels (\code{\link{.marker_allele_counts}},
#' \code{\link{.marker_maf}}). Marker location and block membership are taken from
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
#' @export
compute_marker_stats <- function(geno_meta, markers = NULL, output_file = NULL) {
  if (!is_geno_meta(geno_meta))
    stop("compute_marker_stats: 'geno_meta' is not a valid GenoMeta.")

  arr <- geno_meta$geno_arrays
  if (!is.null(markers)) {
    keep <- dimnames(arr)[[2L]] %in% as.character(markers)
    if (!any(keep))
      stop("compute_marker_stats: none of 'markers' are present in geno_meta.")
    arr <- arr[, keep, , drop = FALSE]
  }
  marker_names <- dimnames(arr)[[2L]]
  d            <- dim(arr)
  n_samp <- d[1L]; n_mark <- d[2L]; ploidy <- d[3L]

  # Encoding + missing tokens (all missingness decided in R, as elsewhere).
  enc_str <- geno_meta$metadata$encoding %||%
    detect_genotype_encoding(.meta_to_allele_matrix(geno_meta))
  miss <- .missing_tokens(enc_str)

  # Single integer-encoding pass; NA = missing.
  codes  <- .encode_marker_codes(arr, na_vals = miss)$codes   # [n_samp x n_mark x ploidy]
  is_na  <- is.na(codes)

  # Missing / call rate and observed count per marker.
  n_cells      <- n_samp * ploidy
  miss_cells   <- apply(is_na, 2L, sum)
  missing_rate <- miss_cells / n_cells
  n_obs        <- n_cells - miss_cells

  # Observed heterozygosity (diploid only): scored samples whose two copies differ.
  if (ploidy >= 2L) {
    l1 <- codes[, , 1L]; l2 <- codes[, , 2L]      # [n_samp x n_mark]
    if (n_mark == 1L) { l1 <- matrix(l1, n_samp, 1L); l2 <- matrix(l2, n_samp, 1L) }
    both_obs <- !is.na(l1) & !is.na(l2)
    het_ij   <- both_obs & (l1 != l2)
    scored   <- colSums(both_obs)
    het      <- ifelse(scored > 0L, colSums(het_ij) / scored, NA_real_)
  } else {
    het <- rep(NA_real_, n_mark)
  }

  # MAF (reuse the shared kernel-backed helper).
  maf_all <- .marker_maf(geno_meta)
  maf     <- unname(maf_all[marker_names])

  # Major / minor allele from per-marker observed counts.
  counts <- .marker_allele_counts(arr, na_vals = miss)   # [n_mark x n_lev]
  levs   <- colnames(counts)
  major  <- rep(NA_character_, n_mark)
  minor  <- rep(NA_character_, n_mark)
  if (length(levs) > 0L) {
    for (j in seq_len(n_mark)) {
      cj  <- counts[j, ]
      obs <- which(cj > 0L)
      if (length(obs) == 0L) next
      major[j] <- levs[obs[which.max(cj[obs])]]
      if (length(obs) >= 2L) minor[j] <- levs[obs[which.min(cj[obs])]]
    }
  }

  stats <- data.table::data.table(
    Marker       = marker_names,
    N_obs        = as.integer(n_obs),
    Missing_rate = as.numeric(missing_rate),
    Call_rate    = as.numeric(1 - missing_rate),
    MAF          = as.numeric(maf),
    Major        = major,
    Minor        = minor,
    Het          = as.numeric(het))

  # Location + block membership, in genomic order, via the shared views.
  mm <- geno_meta_views(geno_meta)$marker_meta[
    , .(Marker = marker, Chrom = chromosome, Pos = position, Block = block_id)]
  stats <- merge(mm, stats, by = "Marker", all.y = TRUE, sort = FALSE)
  data.table::setorder(stats, Chrom, Pos)
  data.table::setcolorder(stats,
    c("Marker", "Chrom", "Pos", "Block", "N_obs", "Missing_rate",
      "Call_rate", "MAF", "Major", "Minor", "Het"))

  if (!is.null(output_file))
    data.table::fwrite(stats, output_file, sep = "\t")
  stats[]
}
