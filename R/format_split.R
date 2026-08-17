# =============================================================================
# format_split.R
#
# split_haplotypes()   — separate diploid genotype strings into a 3D array.
# combine_haplotypes() — collapse a 3D haplotype array back to a 2D data frame.
# translate_alleles()  — re-encode a 3D genotype array between formats.
# =============================================================================

#' Split diploid genotype strings into a 3D haplotype array
#'
#' Converts a 2D genotype matrix (samples × markers OR markers × samples) into
#' a 3D array \code{(n_samples, n_markers, ploidy)} by splitting each cell on
#' the specified delimiter.
#'
#' @param geno_data A matrix or data.frame of genotype strings.
#' @param split_char Allele delimiter (default \code{"/"}). Use \code{""} for
#'   no delimiter — each character is treated as one allele (2-char cells for
#'   diploids, e.g. \code{"AT"}).
#' @param orientation Integer: \code{1} = rows are samples, \code{2} = rows are
#'   markers. Auto-detected when \code{NULL}.
#' @param ploidy Integer ploidy level (default \code{2}).
#' @param sample_names Optional character vector of sample names. Inferred from
#'   the matrix dimension names when \code{NULL}.
#' @param marker_names Optional character vector of marker names.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{array}}{3D character array \code{(n_samples, n_markers, ploidy)}.}
#'   \item{\code{samples}}{Character vector of sample names.}
#'   \item{\code{markers}}{Character vector of marker names.}
#'   \item{\code{encode}}{Detected encoding (\code{"selectiontools"},
#'     \code{"nucleotide"}, or \code{"numeric"}).}
#' }
#'
#' @examples
#' geno <- matrix(c("1/1","1/4","4/4","1/1"), nrow = 2,
#'                dimnames = list(c("Mkr1","Mkr2"), c("S1","S2")))
#' res <- split_haplotypes(geno, orientation = 2)
#' dim(res$array)  # 2 x 2 x 2
#'
#' @export
split_haplotypes <- function(geno_data,
                              split_char   = "/",
                              orientation  = NULL,
                              ploidy       = 2L,
                              sample_names = NULL,
                              marker_names = NULL) {

  if (is.data.frame(geno_data) || inherits(geno_data, "data.table")) {
    geno_mat <- as.matrix(geno_data)
  } else {
    geno_mat <- as.matrix(geno_data)
  }

  # Detect orientation
  if (is.null(orientation)) {
    orientation <- .detect_orientation(geno_mat)
  }

  # Normalise to (rows = samples, cols = markers)
  if (orientation == 2L) geno_mat <- t(geno_mat)

  n_samp <- nrow(geno_mat)
  n_mark <- ncol(geno_mat)

  if (is.null(sample_names)) {
    sample_names <- rownames(geno_mat) %||% paste0("Sample", seq_len(n_samp))
  }
  if (is.null(marker_names)) {
    marker_names <- colnames(geno_mat) %||% paste0("Mkr", seq_len(n_mark))
  }

  # Work on a character matrix so both the kernel and the fallback see allele
  # strings (the caller may pass a numeric matrix).
  if (!is.character(geno_mat)) storage.mode(geno_mat) <- "character"

  ploidy       <- as.integer(ploidy)
  hap_dimnames <- list(sample_names, marker_names, paste0("H", seq_len(ploidy)))

  # Split each cell into its ploidy alleles: compiled kernel when available, the
  # pure-R loop (identical results) otherwise — see .split_genotype_r().
  if (.hs_use_rcpp("split_genotype_cpp")) {
    hap_array <- split_genotype_cpp(geno_mat, split_char, ploidy)
    dimnames(hap_array) <- hap_dimnames
  } else {
    hap_array <- .split_genotype_r(geno_mat, split_char, ploidy, hap_dimnames)
  }

  encoding <- tryCatch(detect_genotype_encoding(geno_mat), error = function(e) "unknown")

  list(array   = hap_array,
       samples = sample_names,
       markers = marker_names,
       encode  = encoding)
}


# Pure-R fallback for split_genotype_cpp: the original per-marker strsplit /
# substr loop. Kept as the guaranteed path when the compiled kernel is
# unavailable (package built without a compiler) or disabled via
# options(haploStack.use_rcpp = FALSE). Results are identical to the kernel
# (parity-tested in test-format-cells.R).
.split_genotype_r <- function(geno_mat, split_char, ploidy, hap_dimnames) {
  n_samp <- nrow(geno_mat)
  n_mark <- ncol(geno_mat)
  hap_array <- array(NA_character_, dim = c(n_samp, n_mark, ploidy),
                     dimnames = hap_dimnames)
  for (j in seq_len(n_mark)) {
    vals <- as.character(geno_mat[, j])
    if (nzchar(split_char)) {
      parts <- strsplit(vals, split_char, fixed = TRUE)
      for (h in seq_len(ploidy)) {
        hap_array[, j, h] <- vapply(parts, function(p)
          if (length(p) >= h) p[[h]] else NA_character_, character(1L))
      }
    } else {
      # No delimiter: split by character position
      for (h in seq_len(ploidy)) {
        hap_array[, j, h] <- substr(vals, h, h)
        hap_array[hap_array[, j, h] == "", j, h] <- NA_character_
      }
    }
  }
  hap_array
}


# .reshape_wide_alleles(mat, n_markers, ploidy)
#
# Converts a wide allele matrix (n_samples × ploidy*n_markers), where each
# cell is a single character, into a (n_samples × n_markers) matrix where
# each cell is ploidy characters concatenated (e.g. "A"+"T" → "AT").
# Suitable for split_haplotypes(split_char = "").
#
# Errors if any non-NA cell has nchar > 1 (ambiguous multi-char wide input).
.reshape_wide_alleles <- function(mat, n_markers, ploidy) {
  ploidy    <- as.integer(ploidy)
  n_samples <- nrow(mat)

  probe <- as.vector(mat[seq_len(min(3L, n_samples)), , drop = FALSE])
  probe <- probe[!is.na(probe) & nzchar(probe)]
  if (length(probe) && any(nchar(probe) > 1L))
    stop(".reshape_wide_alleles: cells contain multi-character values. ",
         "Wide format requires exactly one character per cell.")

  out <- matrix(NA_character_, nrow = n_samples, ncol = n_markers)
  for (k in seq_len(n_markers)) {
    cols     <- seq.int((k - 1L) * ploidy + 1L, k * ploidy)
    out[, k] <- apply(mat[, cols, drop = FALSE], 1L,
                      function(r) paste0(r, collapse = ""))
  }
  rownames(out) <- rownames(mat)
  out
}


# Alias used internally by process_genotype
.split_genotype <- split_haplotypes


# Pure-R fallback for combine_haplotype_cpp: the original per-cell paste double
# loop. Guaranteed path when the compiled kernel is unavailable or disabled;
# identical results to the kernel (parity-tested in test-format-cells.R).
.combine_genotype_r <- function(hap_array, sample_names, marker_names,
                                split_char, n_samp, n_mark, ploidy) {
  geno_mat <- matrix(NA_character_, nrow = n_samp, ncol = n_mark,
                     dimnames = list(sample_names, marker_names))
  for (i in seq_len(n_samp)) {
    for (j in seq_len(n_mark)) {
      alleles <- as.character(hap_array[i, j, ])
      alleles[is.na(alleles)] <- "-1"
      geno_mat[i, j] <- paste(alleles, collapse = split_char)
    }
  }
  geno_mat
}


# Collapse the 3-D haplotype array [n_samples x n_markers x ploidy] into an
# [n_samples x n_markers] matrix of SORTED, separator-less genotype strings —
# the structural transform behind the "genotype" row mode of
# plot_haplotype_alignment() (paste0(sort(allele), collapse = "") per cell).
# Compiled kernel when available, the vectorised pure-R path otherwise; results
# are identical (parity-tested in test-format-cells.R). NA alleles are dropped
# (matching sort()), an all-NA cell yields "".
# @keywords internal
.collapse_genotype_sorted <- function(arr,
                                      use_rcpp =
                                        .hs_use_rcpp("collapse_genotype_sorted_cpp")) {
  if (use_rcpp) {
    g <- collapse_genotype_sorted_cpp(as.character(arr), dim(arr))
    dimnames(g) <- dimnames(arr)[1:2]
    g
  } else {
    .collapse_genotype_sorted_r(arr)
  }
}


# Pure-R fallback for collapse_genotype_sorted_cpp. Diploid fast path uses
# vectorised pmin()/pmax() (lexicographic on characters, matching sort()) with
# explicit NA handling; a general per-cell loop covers ploidy != 2 (which the
# genotype row mode never reaches — it is diploid-only). Identical results to
# the kernel.
.collapse_genotype_sorted_r <- function(arr) {
  d      <- dim(arr)
  n_samp <- d[1L]; n_mark <- d[2L]; ploidy <- d[3L]

  if (ploidy == 2L) {
    h1  <- matrix(arr[, , 1L], n_samp, n_mark)
    h2  <- matrix(arr[, , 2L], n_samp, n_mark)
    na1 <- is.na(h1); na2 <- is.na(h2)
    g   <- paste0(pmin(h1, h2), pmax(h1, h2))   # NA wherever either allele is NA
    one <- xor(na1, na2)                        # exactly one present -> keep it
    g[one] <- ifelse(na1[one], h2[one], h1[one])
    g[na1 & na2] <- ""                          # both NA -> ""
    dim(g)      <- c(n_samp, n_mark)
    dimnames(g) <- dimnames(arr)[1:2]
    return(g)
  }

  g <- matrix(NA_character_, n_samp, n_mark, dimnames = dimnames(arr)[1:2])
  for (j in seq_len(n_mark)) {
    for (i in seq_len(n_samp)) {
      a <- arr[i, j, ]
      g[i, j] <- paste0(sort(a[!is.na(a)]), collapse = "")
    }
  }
  g
}


#' Collapse a 3D haplotype array back to a 2D genotype data frame
#'
#' Reverses \code{\link{split_haplotypes}} by joining haplotype layers with
#' a delimiter. Optionally adds sample or marker names as the first column.
#'
#' @param hap_array 3D character array \code{(n_samples, n_markers, ploidy)}.
#' @param sample_names Character vector of sample names.
#' @param split_char   Allele separator to insert between haplotype alleles
#'   (default \code{"/"}). Use \code{""} for concatenation.
#' @param output_orientation \code{1} = rows are samples (default), \code{2} =
#'   rows are markers.
#' @param name_as_col  Logical. When \code{TRUE} the row identifier (sample or
#'   marker name) is placed in the first column. Default \code{FALSE}.
#'
#' @return A data.frame ready for writing.
#'
#' @examples
#' res  <- split_haplotypes(matrix(c("1/1","4/4"), 1,
#'                 dimnames = list("Mkr1", c("S1","S2"))), orientation = 2)
#' df   <- combine_haplotypes(res$array, res$samples)
#' print(df)
#'
#' @export
combine_haplotypes <- function(hap_array,
                                sample_names       = NULL,
                                split_char         = "/",
                                output_orientation = 1L,
                                name_as_col        = FALSE) {

  n_samp <- dim(hap_array)[1L]
  n_mark <- dim(hap_array)[2L]
  ploidy <- dim(hap_array)[3L]

  if (is.null(sample_names))
    sample_names <- dimnames(hap_array)[[1L]] %||% paste0("S", seq_len(n_samp))
  marker_names <- dimnames(hap_array)[[2L]] %||% paste0("M", seq_len(n_mark))

  # Collapse haplotype layers to genotype strings: compiled kernel when available,
  # the pure-R double loop (identical results) otherwise — see .combine_genotype_r().
  ploidy <- as.integer(ploidy)
  if (.hs_use_rcpp("combine_haplotype_cpp")) {
    geno_mat <- combine_haplotype_cpp(as.character(hap_array),
                                      as.integer(c(n_samp, n_mark, ploidy)),
                                      split_char, "-1")
    dimnames(geno_mat) <- list(sample_names, marker_names)
  } else {
    geno_mat <- .combine_genotype_r(hap_array, sample_names, marker_names,
                                    split_char, n_samp, n_mark, ploidy)
  }

  if (output_orientation == 2L) geno_mat <- t(geno_mat)

  df <- as.data.frame(geno_mat, stringsAsFactors = FALSE)

  if (name_as_col) {
    row_ids <- rownames(df)
    df <- cbind(data.frame(Name = row_ids, stringsAsFactors = FALSE), df)
    rownames(df) <- NULL
  }

  df
}


#' Translate a 3D genotype array between allele encodings
#'
#' Converts all allele values in a \code{(n_samples, n_markers, ploidy)} array
#' from one encoding to another using \code{DEFAULT_ENCODING_DICT} or a custom
#' \code{char_dict}.
#'
#' @param geno_array  3D character array.
#' @param from_encode Source encoding: \code{"nucleotide"} (A/C/G/T),
#'   \code{"selectiontools"} (numeric allele codes), or \code{"allele"} (dosage;
#'   requires \code{minor_alleles}). The legacy names \code{"char"} and
#'   \code{"number"} were removed in 0.5.1.
#' @param to_encode   Target encoding: same options.
#' @param char_dict   Encoding dictionary.
#' @param minor_alleles Optional named character vector (marker → minor allele),
#'   required when converting \emph{to} or \emph{from} \code{"allele"} encoding.
#'
#' @return 3D character array with translated allele values.
#'
#' @export
translate_alleles <- function(geno_array,
                               from_encode,
                               to_encode,
                               char_dict     = DEFAULT_ENCODING_DICT,
                               minor_alleles = NULL) {

  from_encode <- .encode_arg(from_encode, "translate_alleles")
  to_encode   <- .encode_arg(to_encode,   "translate_alleles")
  if (identical(from_encode, to_encode)) return(geno_array)

  # Build lookup vector
  if (from_encode == "nucleotide" && to_encode == "selectiontools") {
    lut <- setNames(as.character(unlist(DEFAULT_ENCODING_DICT)),
                    names(DEFAULT_ENCODING_DICT))
  } else if (from_encode == "selectiontools" && to_encode == "nucleotide") {
    lut <- REVERSED_ENCODING_DICT
  } else if (to_encode == "allele") {
    if (is.null(minor_alleles)) {
      stop("translate_alleles: minor_alleles required for dosage encoding.")
    }
    # Convert to dosage (count of minor allele) per marker
    n_samp <- dim(geno_array)[1L]
    n_mark <- dim(geno_array)[2L]
    ploidy <- dim(geno_array)[3L]
    out <- array(NA_character_, dim = dim(geno_array),
                 dimnames = dimnames(geno_array))
    marker_names <- dimnames(geno_array)[[2L]]
    for (j in seq_len(n_mark)) {
      mn <- minor_alleles[marker_names[j]]
      for (i in seq_len(n_samp)) {
        alleles <- geno_array[i, j, ]
        dosage  <- sum(alleles == mn, na.rm = TRUE)
        out[i, j, ] <- as.character(c(dosage, rep(NA, ploidy - 1L)))
      }
    }
    return(out)
  } else {
    stop("translate_alleles: unsupported encoding '", from_encode,
         "' → '", to_encode, "'.")
  }

  # Vectorised LUT translation
  orig_dim  <- dim(geno_array)
  orig_dims <- dimnames(geno_array)
  flat      <- as.character(geno_array)
  translated <- lut[flat]
  translated[is.na(translated)] <- flat[is.na(translated)]
  array(translated, dim = orig_dim, dimnames = orig_dims)
}
