// =============================================================================
// coancestry.cpp  —  Kernel D: block-IBS coancestry matrix (0.9.0)
//
// Compiled twin of the pure-R double loop in compute_parent_coancestry()
// (R/stacking_monitor.R), which the record measured at 47 of ~49 s per cycle
// under metrics = "all" (thread 013). The R loop stays as the fallback and as
// the parity oracle (tests/testthat/test-maintenance-0.9.0.R).
//
// ENCODING SAFETY (same philosophy as local_gebv.cpp / marker_stats.cpp): no
// allele semantics here. R reduces the genotype slab to INTEGER CODES first
// (equal strings -> equal codes; a missing or empty allele -> NA). The only
// tests in C++ are "code is NA" and integer equality.
//
// BIT-FOR-BIT PARITY with the R loop:
//   * per block and haplotype pairing, IBS = (#equal) / (#comparable) as a
//     double division, NA when nothing is comparable (compute_segment_ibs());
//   * the four pairings and then the blocks are averaged with r_mean(), which
//     reproduces R's mean(): NA/NaN dropped (na.rm = TRUE), long-double sum,
//     division, then the second-pass correction; NaN when nothing remains;
//   * blocks are visited in the order of the list R passes in.
//
// ARRAY LAYOUT: codes is the column-major backing vector of the
// [n_samples x n_markers x ploidy] slab; [i, j, k] (0-based) is at
// i + n_samples * (j + n_markers * k). Only haplotype slots 0 and 1 are used,
// as in the R loop.
// =============================================================================

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

static double r_mean(const std::vector<double>& x) {
  const R_xlen_t n = static_cast<R_xlen_t>(x.size());
  long double s = 0.0;
  for (double v : x) s += v;
  s /= n;
  if (R_FINITE(static_cast<double>(s))) {
    long double t = 0.0;
    for (double v : x) t += (v - s);
    s += t / n;
  }
  return static_cast<double>(s);
}

// [[Rcpp::export]]
NumericMatrix block_ibs_coancestry_cpp(IntegerVector codes, IntegerVector dims,
                                       List blocks) {
  if (dims.size() != 3) stop("block_ibs_coancestry_cpp: dims must have length 3");
  const int P = dims[0], M = dims[1], H = dims[2];
  if (H < 2) stop("block_ibs_coancestry_cpp: need at least two haplotype slots");
  if (static_cast<double>(codes.size()) !=
      static_cast<double>(P) * M * H)
    stop("block_ibs_coancestry_cpp: length(codes) != prod(dims)");

  const int B = blocks.size();
  std::vector<IntegerVector> blk(B);
  for (int b = 0; b < B; ++b) {
    blk[b] = as<IntegerVector>(blocks[b]);
    for (int m : blk[b])
      if (m < 0 || m >= M) stop("block_ibs_coancestry_cpp: marker index out of range");
  }

  const int* cd = codes.begin();
  auto at = [&](int i, int j, int k) -> int {
    return cd[i + static_cast<R_xlen_t>(P) * (j + static_cast<R_xlen_t>(M) * k)];
  };

  NumericMatrix K(P, P);
  std::vector<double> four, bvals;
  four.reserve(4);
  bvals.reserve(B);

  for (int i = 0; i < P; ++i) {
    checkUserInterrupt();
    for (int j = i; j < P; ++j) {
      bvals.clear();
      for (int b = 0; b < B; ++b) {
        four.clear();
        for (int hi = 0; hi < 2; ++hi) for (int hj = 0; hj < 2; ++hj) {
          int n = 0, eq = 0;
          for (int m : blk[b]) {
            const int a1 = at(i, m, hi), a2 = at(j, m, hj);
            if (a1 == NA_INTEGER || a2 == NA_INTEGER) continue;
            ++n;
            if (a1 == a2) ++eq;
          }
          if (n > 0) four.push_back(static_cast<double>(eq) / static_cast<double>(n));
        }
        const double bv = r_mean(four);
        if (!ISNAN(bv)) bvals.push_back(bv);
      }
      const double v = r_mean(bvals);
      K(i, j) = v;
      K(j, i) = v;
    }
  }
  return K;
}
