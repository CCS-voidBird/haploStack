// =============================================================================
// marker_stats.cpp  —  Kernel B: per-marker allele counts + per-cell dosage
//
// Two compiled hot loops that replace four duplicated per-marker R loops
// (table()/sort() over an N x ploidy slab): the allele-frequency / minor-allele /
// dosage machinery in .marker_allele_counts() and .dosage_matrix() (R/marker_stats.R),
// which back .meta_to_dosage_matrix(), .gp_allele_dosage(), .detect_minor_alleles(),
// and .marker_maf().
//
// ENCODING SAFETY (same philosophy as src/local_gebv.cpp): these kernels are
// deliberately encoding-AGNOSTIC and do NO allele/DICT detection. ALL of that —
// building the global allele dictionary (the factor levels), deciding which
// string values are missing / blank, and choosing each marker's target allele
// (minor / major / alt) — is resolved IN R, inside R/marker_stats.R, BEFORE these
// kernels are called. The genotype array arrives here already reduced to 1-based
// INTEGER CODES (NA where R decided the cell is missing).
//
// The ONLY "detection" performed in C++, stated explicitly:
//   1. NA-code  == "missing" for this cell     -> skipped / counted as missing.
//   2. code == target_code                     -> a dosage match (integer ==).
//   3. range checks that a code is a valid dictionary index (fail loudly).
// No string is ever parsed here and no allele semantics are inferred.
//
// ARRAY LAYOUT: the 3-D genotype array [n_samples x n_markers x ploidy] is passed
// as its flat, column-major backing vector `codes` plus `dims = c(n_samples,
// n_markers, ploidy)`. Element [i, j, k] (all 0-based) lives at
//   idx = i + n_samples * (j + n_markers * k).
// =============================================================================

#include <Rcpp.h>
using namespace Rcpp;

// Validate and unpack the shared (codes, dims) array contract.
static inline void check_dims(const IntegerVector& codes,
                              const IntegerVector& dims) {
  if (dims.size() != 3)
    stop("marker_stats: 'dims' must have length 3 (n_samples, n_markers, "
         "ploidy); got %d.", dims.size());
  const double expected = (double)dims[0] * (double)dims[1] * (double)dims[2];
  if ((double)codes.size() != expected)
    stop("marker_stats: length(codes) = %d but prod(dims) = %.0f; the flat code "
         "vector must match the array dimensions.", codes.size(), expected);
  if (dims[0] < 0 || dims[1] < 0 || dims[2] < 0)
    stop("marker_stats: all dims must be >= 0.");
}

// ---------------------------------------------------------------------------
// marker_allele_counts_cpp
//
// Count, per marker, how many (sample, ploidy-layer) cells carry each dictionary
// allele code. NA cells are skipped (detection rule 1). Returns an
// [n_markers x n_codes] IntegerMatrix; result(j, c-1) = number of cells at marker
// j whose code == c. R turns these counts into minor / major / MAF / target.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerMatrix marker_allele_counts_cpp(IntegerVector codes,
                                       IntegerVector dims,
                                       int n_codes) {
  check_dims(codes, dims);
  if (n_codes < 0)
    stop("marker_allele_counts_cpp: n_codes must be >= 0 (got %d).", n_codes);

  const int n_samples = dims[0];
  const int n_markers = dims[1];
  const int ploidy    = dims[2];

  IntegerMatrix out(n_markers, n_codes);

  for (int k = 0; k < ploidy; ++k) {
    const R_xlen_t base_k = (R_xlen_t)n_samples * (R_xlen_t)n_markers * k;
    for (int j = 0; j < n_markers; ++j) {
      if ((j & 255) == 0) Rcpp::checkUserInterrupt();   // 0.9.0 D1-5
      const R_xlen_t base_jk = base_k + (R_xlen_t)n_samples * j;
      for (int i = 0; i < n_samples; ++i) {
        const int code = codes[base_jk + i];
        if (IntegerVector::is_na(code)) continue;          // rule 1: missing
        if (code < 1 || code > n_codes)
          stop("marker_allele_counts_cpp: code %d at [%d, %d, %d] is out of "
               "range [1, %d].", code, i + 1, j + 1, k + 1, n_codes);
        out(j, code - 1) += 1;
      }
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// dosage_from_codes_cpp
//
// Per-cell dosage = number of a marker's ploidy layers whose code equals that
// marker's target allele code (detection rule 2). target_code[j] is a 1-based
// dictionary code, or NA meaning "no target" (dosage 0 wherever observed).
//
// Missing-cell (NA) policy is decided in R and passed as `na_if_any`:
//   na_if_any == true  -> the (sample, marker) cell is NA if ANY layer is missing
//                         (matches .meta_to_dosage_matrix).
//   na_if_any == false -> the cell is NA only if ALL layers are missing; a
//                         partially-missing cell still counts its observed layers
//                         (matches .gp_allele_dosage).
// A fully-missing cell is always NA regardless of na_if_any.
//
// Returns an [n_samples x n_markers] IntegerMatrix.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerMatrix dosage_from_codes_cpp(IntegerVector codes,
                                    IntegerVector dims,
                                    IntegerVector target_code,
                                    bool na_if_any) {
  check_dims(codes, dims);
  const int n_samples = dims[0];
  const int n_markers = dims[1];
  const int ploidy    = dims[2];

  if (target_code.size() != n_markers)
    stop("dosage_from_codes_cpp: length(target_code) = %d but n_markers = %d; "
         "they must match.", target_code.size(), n_markers);

  IntegerMatrix out(n_samples, n_markers);

  for (int j = 0; j < n_markers; ++j) {
    if ((j & 255) == 0) Rcpp::checkUserInterrupt();   // 0.9.0 D1-5
    const int tj    = target_code[j];
    const bool tj_na = IntegerVector::is_na(tj);
    if (!tj_na && tj < 1)
      stop("dosage_from_codes_cpp: target_code[%d] = %d must be >= 1 or NA.",
           j + 1, tj);

    for (int i = 0; i < n_samples; ++i) {
      int n_missing = 0, n_match = 0, n_obs = 0;
      for (int k = 0; k < ploidy; ++k) {
        const R_xlen_t idx =
          (R_xlen_t)i +
          (R_xlen_t)n_samples * ((R_xlen_t)j + (R_xlen_t)n_markers * k);
        const int code = codes[idx];
        if (IntegerVector::is_na(code)) {                  // rule 1: missing
          ++n_missing;
        } else {
          ++n_obs;
          if (!tj_na && code == tj) ++n_match;             // rule 2: dosage match
        }
      }
      if (n_obs == 0 || (na_if_any && n_missing > 0))
        out(i, j) = NA_INTEGER;
      else
        out(i, j) = n_match;
    }
  }
  return out;
}
