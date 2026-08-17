// =============================================================================
// format_cells.cpp  —  Kernel A: haplotype cell split / combine
//
// Two compiled hot loops that replace the per-cell string transforms in
// split_haplotypes() (a per-marker strsplit/substr loop) and combine_haplotypes()
// (a double N x M per-cell paste loop) in R/format_split.R. Both scale with
// N samples x M markers x ploidy.
//
// ENCODING SAFETY (same philosophy as src/local_gebv.cpp and src/marker_stats.cpp):
// these kernels do NO allele/DICT/encoding detection. They are purely STRUCTURAL
// string transforms — split each genotype cell on a caller-supplied separator, or
// join ploidy allele cells with it. The separator `sep`, the ploidy, and the
// missing token `na_token` ("-1") are all decided in R (R/format_split.R). No
// allele is ever interpreted and no encoding is inferred here.
//
// The ONLY per-cell "detection" performed in C++, stated explicitly:
//   1. NA cell            -> propagated (split: all layers NA; combine: na_token).
//   2. sep found / absent -> where a genotype string is cut (fixed strsplit).
//   3. position >= length -> a fixed-width layer is missing (-> NA).
//
// R-SEMANTICS REPRODUCED (verified against R 4.4.3 before implementation):
//   * Delimiter mode (sep != ""): matches strsplit(x, sep, fixed = TRUE) followed
//     by p[[h]] (else NA). R drops a TRAILING empty field but keeps leading /
//     internal empties: "1/"->["1"], "1/4/"->["1","4"], "/4"->["","4"],
//     "1//4"->["1","","4"], ""->[] (all layers NA). An empty token ("") is kept
//     as "" (NOT turned into NA) in delimiter mode.
//   * Fixed-width mode (sep == ""): matches substr(x, h, h) with "" -> NA, i.e.
//     layer h is the h-th character, or NA when the string is shorter.
//
// ARRAY LAYOUT: the 3-D haplotype array [n_samples x n_markers x ploidy] is the
// column-major flat vector; element [i, j, h] (0-based) is at
//   idx = i + n_samples * (j + n_markers * h) = (i + n_samples*j) + n_cell*h,
// where n_cell = n_samples * n_markers.
// =============================================================================

#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <vector>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// split_genotype_cpp
//
// Split an [n_samples x n_markers] character genotype matrix into the flat
// [n_samples x n_markers x ploidy] haplotype array (column-major, with the dim
// attribute set so R receives a ready 3-D array; R attaches dimnames).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
CharacterVector split_genotype_cpp(CharacterMatrix geno,
                                   std::string sep,
                                   int ploidy) {
  if (ploidy < 1)
    stop("split_genotype_cpp: ploidy must be >= 1 (got %d).", ploidy);

  const int n_samp = geno.nrow();
  const int n_mark = geno.ncol();
  const R_xlen_t n_cell = (R_xlen_t)n_samp * (R_xlen_t)n_mark;
  const bool fixed_width = sep.empty();

  CharacterVector out(n_cell * (R_xlen_t)ploidy);

  for (int j = 0; j < n_mark; ++j) {
    for (int i = 0; i < n_samp; ++i) {
      const R_xlen_t cell = (R_xlen_t)i + (R_xlen_t)n_samp * j;

      // Rule 1: a missing genotype cell -> every layer NA.
      if (CharacterVector::is_na(geno[cell])) {
        for (int h = 0; h < ploidy; ++h) out[cell + n_cell * h] = NA_STRING;
        continue;
      }

      const std::string v = as<std::string>(geno[cell]);

      if (fixed_width) {
        const int len = (int)v.size();
        for (int h = 0; h < ploidy; ++h)                       // rule 3
          out[cell + n_cell * h] =
            (h < len) ? String(std::string(1, v[h])) : String(NA_STRING);
      } else {
        // Rule 2: R fixed strsplit — emit the token before each separator; emit
        // the final remainder ONLY if non-empty (drops a trailing empty field);
        // an empty string yields zero tokens.
        std::vector<std::string> parts;
        if (!v.empty()) {
          std::size_t start = 0, pos;
          while ((pos = v.find(sep, start)) != std::string::npos) {
            parts.push_back(v.substr(start, pos - start));
            start = pos + sep.size();
          }
          const std::string rem = v.substr(start);
          if (!rem.empty()) parts.push_back(rem);
        }
        const int np = (int)parts.size();
        for (int h = 0; h < ploidy; ++h)
          out[cell + n_cell * h] =
            (h < np) ? String(parts[h]) : String(NA_STRING);
      }
    }
  }

  out.attr("dim") = IntegerVector::create(n_samp, n_mark, ploidy);
  return out;
}

// ---------------------------------------------------------------------------
// combine_haplotype_cpp
//
// Collapse the flat [n_samples x n_markers x ploidy] haplotype array back to an
// [n_samples x n_markers] character genotype matrix by joining each cell's ploidy
// alleles with `sep`; missing alleles are written as `na_token` (matches the R
// path's paste() with NA -> "-1").
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
CharacterMatrix combine_haplotype_cpp(CharacterVector hap,
                                      IntegerVector dims,
                                      std::string sep,
                                      std::string na_token) {
  if (dims.size() != 3)
    stop("combine_haplotype_cpp: 'dims' must have length 3 (n_samples, "
         "n_markers, ploidy); got %d.", dims.size());
  const int n_samp = dims[0], n_mark = dims[1], ploidy = dims[2];
  const R_xlen_t n_cell = (R_xlen_t)n_samp * (R_xlen_t)n_mark;
  if ((R_xlen_t)hap.size() != n_cell * (R_xlen_t)ploidy)
    stop("combine_haplotype_cpp: length(hap) = %d but prod(dims) = %.0f.",
         hap.size(), (double)n_cell * ploidy);

  CharacterMatrix out(n_samp, n_mark);

  for (int j = 0; j < n_mark; ++j) {
    for (int i = 0; i < n_samp; ++i) {
      const R_xlen_t cell = (R_xlen_t)i + (R_xlen_t)n_samp * j;
      std::string s;
      for (int h = 0; h < ploidy; ++h) {
        if (h) s += sep;
        const R_xlen_t idx = cell + n_cell * h;
        if (CharacterVector::is_na(hap[idx])) s += na_token;    // rule 1
        else                                  s += as<std::string>(hap[idx]);
      }
      out[cell] = s;                            // column-major -> out(i, j)
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// collapse_genotype_sorted_cpp
//
// Collapse the flat [n_samples x n_markers x ploidy] haplotype array to an
// [n_samples x n_markers] character genotype matrix by SORTING each cell's
// non-missing ploidy alleles and concatenating them with NO separator. This is
// the compiled twin of the per-group R expression used by the "genotype" row
// mode in plot_haplotype_alignment(): paste0(sort(allele), collapse = "").
//
// Same structural-only philosophy as combine_haplotype_cpp above: no allele /
// encoding is interpreted. The ONLY per-cell rules, stated explicitly:
//   1. NA layer -> dropped (matches R's sort(), which discards NA by default).
//   2. all layers NA -> "" (empty string; sort() -> character(0), paste0 -> "").
//   3. surviving alleles sorted ascending, then concatenated.
//
// SORT-ORDER PARITY: sorting is byte order via std::sort on std::string. This
// matches R's locale-aware sort() only because allele tokens here are single
// characters (nucleotides A/C/G/T or numeric codes), for which C-locale and
// R-locale collation agree — the same assumption the rest of this file relies
// on. Multi-character tokens are NOT guaranteed to match a non-C R locale.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
CharacterMatrix collapse_genotype_sorted_cpp(CharacterVector hap,
                                             IntegerVector dims) {
  if (dims.size() != 3)
    stop("collapse_genotype_sorted_cpp: 'dims' must have length 3 (n_samples, "
         "n_markers, ploidy); got %d.", dims.size());
  const int n_samp = dims[0], n_mark = dims[1], ploidy = dims[2];
  const R_xlen_t n_cell = (R_xlen_t)n_samp * (R_xlen_t)n_mark;
  if ((R_xlen_t)hap.size() != n_cell * (R_xlen_t)ploidy)
    stop("collapse_genotype_sorted_cpp: length(hap) = %d but prod(dims) = %.0f.",
         hap.size(), (double)n_cell * ploidy);

  CharacterMatrix out(n_samp, n_mark);
  std::vector<std::string> parts;
  parts.reserve(ploidy);

  for (int j = 0; j < n_mark; ++j) {
    for (int i = 0; i < n_samp; ++i) {
      const R_xlen_t cell = (R_xlen_t)i + (R_xlen_t)n_samp * j;
      parts.clear();
      for (int h = 0; h < ploidy; ++h) {
        const R_xlen_t idx = cell + n_cell * h;
        if (!CharacterVector::is_na(hap[idx]))        // rule 1: drop NA layers
          parts.push_back(as<std::string>(hap[idx]));
      }
      std::sort(parts.begin(), parts.end());          // rule 3
      std::string s;                                  // rule 2: all-NA -> ""
      for (std::size_t k = 0; k < parts.size(); ++k) s += parts[k];
      out[cell] = s;                            // column-major -> out(i, j)
    }
  }
  return out;
}
