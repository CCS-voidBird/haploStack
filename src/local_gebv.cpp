// =============================================================================
// local_gebv.cpp  —  Kernel A: block-summed local GEBV
//
// One compiled hot loop for the haplotype-catalogue compute path. Given a
// genotype matrix already reduced to integer ALLELE CODES, a per-(marker, code)
// CONTRIBUTION table, and a marker->block map, it accumulates each haplotype
// copy's per-block local GEBV in a single O(n_copies x n_markers) scatter pass.
// No long-format table is ever materialised.
//
// USAGE (R side): called only by haploStack:::.hap_effects_rcpp() in
//   R/catalogue_engine.R, which is reached via .local_gebv_engine() ->
//   .hap_effects() when .use_rcpp() is TRUE (option haploStack.use_rcpp, default
//   TRUE) and the kernel symbol is loadable. When unavailable, the pure-R
//   data.table path (.hap_effects_r) is used instead — results are identical and
//   parity-tested in tests/testthat/test-catalogue-engine.R.
//
// ENCODING SAFETY: this kernel is deliberately encoding-AGNOSTIC. All allele
// strings and the genotype encoding ("selectiontools" / "nucleotide" / numeric)
// are resolved to integer codes and a numeric contribution table IN R, inside
// .catalogue_core() / .hap_effects_rcpp(). The kernel never parses an allele or
// knows the encoding; it only indexes doubles by integer codes. This is what
// makes it safe to feed differently-encoded GenoMeta objects through the same
// compiled routine.
//
// INPUT CONTRACT (all 1-based indices come from R; validated below):
//   code_mat    IntegerMatrix [n_copies x n_markers]
//                 code_mat(i, j) = 1-based allele code of copy i at marker j,
//                 or NA_integer_ when that cell has no observed allele
//                 (missing / blank). NA cells contribute nothing.
//   contrib_mat NumericMatrix [n_markers x n_code]
//                 contrib_mat(j, c - 1) = the contribution of allele code c at
//                 marker j (the marker effect, or Effect*(1 - allele_freq) when
//                 centred). Codes with no matching effect are pre-filled 0 in R.
//   block_idx   IntegerVector [n_markers]
//                 block_idx[j] = 1-based block index of marker j into the block
//                 axis of the output, or 0 to skip marker j (its block has no
//                 observed allele anywhere).
//   n_blocks    int  number of output block columns (>= 0).
//
// OUTPUT:
//   NumericMatrix [n_copies x n_blocks]  result(i, b) = sum over markers j with
//   block_idx[j] == b+1 of contrib_mat(j, code_mat(i, j) - 1), skipping NA cells.
// =============================================================================

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix local_gebv_block_sum_cpp(IntegerMatrix code_mat,
                                       NumericMatrix contrib_mat,
                                       IntegerVector block_idx,
                                       int n_blocks) {
  const int n_copies  = code_mat.nrow();
  const int n_markers = code_mat.ncol();
  const int n_code    = contrib_mat.ncol();

  // ── Input asserts (fail loudly rather than corrupt a GEBV) ────────────────
  if (contrib_mat.nrow() != n_markers)
    stop("local_gebv_block_sum_cpp: contrib_mat has %d rows but code_mat has %d "
         "markers (columns); they must match.", contrib_mat.nrow(), n_markers);
  if (block_idx.size() != n_markers)
    stop("local_gebv_block_sum_cpp: block_idx has length %d but code_mat has %d "
         "markers (columns); they must match.", block_idx.size(), n_markers);
  if (n_blocks < 0)
    stop("local_gebv_block_sum_cpp: n_blocks must be >= 0 (got %d).", n_blocks);

  NumericMatrix result(n_copies, n_blocks);

  for (int j = 0; j < n_markers; ++j) {
    const int b = block_idx[j];              // 1-based block, 0 = skip
    if (b == 0) continue;
    if (b < 0 || b > n_blocks)
      stop("local_gebv_block_sum_cpp: block_idx[%d] = %d is out of range "
           "[0, %d].", j + 1, b, n_blocks);
    const int bcol = b - 1;

    for (int i = 0; i < n_copies; ++i) {
      const int code = code_mat(i, j);       // 1-based code, NA = no allele
      if (IntegerMatrix::is_na(code)) continue;
      if (code < 1 || code > n_code)
        stop("local_gebv_block_sum_cpp: code_mat(%d, %d) = %d is out of range "
             "[1, %d] (contrib_mat columns).", i + 1, j + 1, code, n_code);
      result(i, bcol) += contrib_mat(j, code - 1);
    }
  }

  return result;
}
