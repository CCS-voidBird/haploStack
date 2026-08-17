# =============================================================================
# test-marker-stats.R
#
# QA for Kernel B (R/marker_stats.R + src/marker_stats.cpp):
#   - marker_allele_counts_cpp / dosage_from_codes_cpp input asserts
#   - pure-R vs Rcpp parity for .marker_allele_counts() and .dosage_matrix()
#   - end-to-end parity across the four refactored call sites
#     (.meta_to_dosage_matrix, .gp_allele_dosage, .detect_minor_alleles,
#      .marker_maf), toggled via options(haploStack.use_rcpp)
#   - edge cases: NA-any vs NA-all policy, monomorphic markers, count ties,
#     all-missing markers.
# =============================================================================

library(testthat)
library(haploStack)

counts_kernel <- function() is.function(haploStack:::.hs_kernel("marker_allele_counts_cpp"))
dosage_kernel <- function() is.function(haploStack:::.hs_kernel("dosage_from_codes_cpp"))

# A small diploid array with deliberate edge cases:
#   Mkr1: clean biallelic          Mkr2: monomorphic ("1" only)
#   Mkr3: partial-missing cells     Mkr4: fully-missing (all NA)
#   Mkr5: count tie between "1"/"4"
.mk_arr <- function() {
  s <- paste0("S", 1:6); m <- paste0("Mkr", 1:5)
  arr <- array(NA_character_, dim = c(6L, 5L, 2L),
               dimnames = list(s, m, c("1", "2")))
  arr[, 1, ] <- c("1","1","4","4","1","4",  "4","1","4","1","1","4")   # biallelic
  arr[, 2, ] <- "1"                                                     # monomorphic
  arr[, 3, 1] <- c("1","4",NA,"1","4","1"); arr[, 3, 2] <- c(NA,"4","4","1",NA,"4")
  arr[, 4, ] <- NA_character_                                           # all missing
  arr[, 5, ] <- c("1","1","1","4","4","4",  "1","1","1","4","4","4")   # 6x"1" 6x"4" tie
  arr
}

# ── kernel input asserts ──────────────────────────────────────────────────────

test_that("marker_allele_counts_cpp validates its inputs", {
  skip_if_not(counts_kernel(), "kernel not compiled")
  k <- haploStack:::.hs_kernel("marker_allele_counts_cpp")
  codes <- c(1L, 2L, NA, 1L)                     # 2 samples x 2 markers x 1
  expect_silent(k(codes, c(2L, 2L, 1L), 2L))
  expect_error(k(codes, c(2L, 2L), 2L), "length 3")
  expect_error(k(codes, c(2L, 2L, 2L), 2L), "prod\\(dims\\)")
  expect_error(k(c(1L, 9L, 1L, 1L), c(2L, 2L, 1L), 2L), "out of range")
})

test_that("dosage_from_codes_cpp validates target length", {
  skip_if_not(dosage_kernel(), "kernel not compiled")
  k <- haploStack:::.hs_kernel("dosage_from_codes_cpp")
  codes <- c(1L, 2L, 1L, 2L)                     # 2 samples x 2 markers x 1
  expect_silent(k(codes, c(2L, 2L, 1L), c(1L, 2L), TRUE))
  expect_error(k(codes, c(2L, 2L, 1L), c(1L), TRUE), "must match")
})

# ── low-level parity ──────────────────────────────────────────────────────────

test_that(".marker_allele_counts: pure-R and Rcpp agree", {
  skip_if_not(counts_kernel(), "kernel not compiled")
  arr <- .mk_arr()
  a <- haploStack:::.marker_allele_counts(arr, use_rcpp = FALSE)
  b <- haploStack:::.marker_allele_counts(arr, use_rcpp = TRUE)
  expect_identical(a, b)
  # spot-check: Mkr1 has six "1" and six "4"; Mkr4 all zero
  expect_identical(unname(a["Mkr1", ]), c(6L, 6L))
  expect_identical(unname(a["Mkr4", ]), c(0L, 0L))
})

test_that(".dosage_matrix: pure-R and Rcpp agree under both NA policies", {
  skip_if_not(dosage_kernel(), "kernel not compiled")
  arr <- .mk_arr()
  tgt <- rep("4", 5L)                            # count copies of "4"
  for (pol in c(TRUE, FALSE)) {
    a <- haploStack:::.dosage_matrix(arr, tgt, na_if_any = pol, use_rcpp = FALSE)
    b <- haploStack:::.dosage_matrix(arr, tgt, na_if_any = pol, use_rcpp = TRUE)
    expect_identical(a, b)
  }
})

test_that(".dosage_matrix NA policy differs on partial-missing cells", {
  arr <- .mk_arr()
  tgt <- rep("4", 5L)
  any_na <- haploStack:::.dosage_matrix(arr, tgt, na_if_any = TRUE,  use_rcpp = FALSE)
  all_na <- haploStack:::.dosage_matrix(arr, tgt, na_if_any = FALSE, use_rcpp = FALSE)
  # Mkr3 sample S1: layers ("1", NA) -> any-policy NA, all-policy counts 0
  expect_true(is.na(any_na["S1", "Mkr3"]))
  expect_identical(all_na["S1", "Mkr3"], 0L)
  # fully-missing marker is NA under both policies
  expect_true(all(is.na(any_na[, "Mkr4"])))
  expect_true(all(is.na(all_na[, "Mkr4"])))
})

# ── end-to-end parity across the toggle ───────────────────────────────────────

# A GenoMeta with supported edge cases: a partial-missing cell and a monomorphic
# marker. (A *fully*-missing marker is exercised at the low-level .dosage_matrix /
# .marker_allele_counts tests above; .detect_minor_alleles does not support it —
# neither the original nor the refactor, since this sim char_dict has no na_char.)
.mk_meta_with_na <- function() {
  meta <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 3L,
                        block_size = 4L, seed = 7L)
  meta$geno_arrays[2, 3, 1] <- NA_character_     # partial-missing cell
  meta$geno_arrays[, 5, ]   <- "1"               # monomorphic marker
  meta
}

engine_toggle <- function(rcpp, expr) {
  withr::with_options(list(haploStack.use_rcpp = rcpp), expr)
}

test_that(".meta_to_dosage_matrix: pure-R and Rcpp agree (auto + supplied target)", {
  skip_if_not(dosage_kernel(), "kernel not compiled")
  meta <- .mk_meta_with_na()
  a <- engine_toggle(FALSE, haploStack:::.meta_to_dosage_matrix(meta))
  b <- engine_toggle(TRUE,  haploStack:::.meta_to_dosage_matrix(meta))
  expect_identical(a, b)
  tgt <- rep("4", dim(meta$geno_arrays)[2L])
  a2 <- engine_toggle(FALSE, haploStack:::.meta_to_dosage_matrix(meta, tgt))
  b2 <- engine_toggle(TRUE,  haploStack:::.meta_to_dosage_matrix(meta, tgt))
  expect_identical(a2, b2)
})

test_that(".gp_allele_dosage: pure-R and Rcpp agree (M + allele_map)", {
  skip_if_not(dosage_kernel(), "kernel not compiled")
  meta <- .mk_meta_with_na()
  a <- engine_toggle(FALSE, haploStack:::.gp_allele_dosage(meta, impute = "none"))
  b <- engine_toggle(TRUE,  haploStack:::.gp_allele_dosage(meta, impute = "none"))
  expect_equal(a$M, b$M)
  expect_identical(a$allele_map, b$allele_map)
})

test_that(".detect_minor_alleles: pure-R and Rcpp agree", {
  skip_if_not(counts_kernel(), "kernel not compiled")
  meta <- .mk_meta_with_na()
  a <- engine_toggle(FALSE, haploStack:::.detect_minor_alleles(meta$geno_arrays))
  b <- engine_toggle(TRUE,  haploStack:::.detect_minor_alleles(meta$geno_arrays))
  expect_identical(a, b)
})

test_that(".marker_maf: pure-R and Rcpp agree", {
  skip_if_not(counts_kernel(), "kernel not compiled")
  meta <- .mk_meta_with_na()
  a <- engine_toggle(FALSE, haploStack:::.marker_maf(meta))
  b <- engine_toggle(TRUE,  haploStack:::.marker_maf(meta))
  expect_equal(a, b)
})

# ── compute_marker_stats (exported) ───────────────────────────────────────────

test_that("compute_marker_stats: shape, columns, and MAF agree with .marker_maf", {
  meta  <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 3L,
                         block_size = 4L, seed = 7L)
  stats <- compute_marker_stats(meta)

  expect_s3_class(stats, "data.table")
  expect_identical(names(stats),
    c("Marker", "Chrom", "Pos", "Block", "N_obs", "Missing_rate",
      "Call_rate", "MAF", "Major", "Minor", "Het"))
  expect_equal(nrow(stats), dim(meta$geno_arrays)[2L])

  # MAF matches the shared helper, marker-for-marker.
  maf <- haploStack:::.marker_maf(meta)
  expect_equal(stats$MAF, unname(maf[stats$Marker]))

  # Block membership populated when block assignments exist.
  expect_true(all(!is.na(stats$Block)))
  expect_equal(stats$Call_rate, 1 - stats$Missing_rate)
})

test_that("compute_marker_stats: missing rate matches an injected NA cell", {
  meta <- sim_geno_meta(n_samples = 10L, n_markers = 6L, n_chrom = 2L,
                        block_size = 3L, seed = 11L)
  d      <- dim(meta$geno_arrays)                 # [10 x 6 x 2]
  n_cell <- d[1L] * d[3L]
  meta$geno_arrays[2, 3, 1] <- NA_character_      # one missing cell in marker 3
  stats  <- compute_marker_stats(meta)
  m3     <- stats[Marker == dimnames(meta$geno_arrays)[[2L]][3L]]
  expect_equal(m3$Missing_rate, 1 / n_cell)
  expect_equal(m3$N_obs, n_cell - 1L)
})

test_that("compute_marker_stats: Het is 1 for all-het, 0 for monomorphic", {
  meta <- sim_geno_meta(n_samples = 8L, n_markers = 4L, n_chrom = 1L,
                        block_size = 2L, seed = 3L)
  meta$geno_arrays[, 1, 1] <- "1"; meta$geno_arrays[, 1, 2] <- "4"  # all het
  meta$geno_arrays[, 2, ]  <- "1"                                   # monomorphic
  stats <- compute_marker_stats(meta)
  mk    <- dimnames(meta$geno_arrays)[[2L]]
  expect_equal(stats[Marker == mk[1L]]$Het, 1)
  expect_equal(stats[Marker == mk[2L]]$Het, 0)
})

test_that("compute_marker_stats: output_file round-trips", {
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 8L, n_chrom = 2L,
                        block_size = 4L, seed = 5L)
  f     <- tempfile(fileext = ".txt")
  stats <- compute_marker_stats(meta, output_file = f)
  expect_true(file.exists(f))
  back <- data.table::fread(f)
  expect_equal(nrow(back), nrow(stats))
  expect_equal(back$Marker, stats$Marker)
})

# ── write_block_allele_sheet (exported) ───────────────────────────────────────

test_that("write_block_allele_sheet: uniform columns, counts and freq per block", {
  meta  <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 3L,
                         block_size = 4L, seed = 7L)
  sheet <- write_block_allele_sheet(meta)

  expect_s3_class(sheet, "data.table")
  expect_identical(names(sheet),
    c("Block", "Chrom", "Signature", "Markers", "Alleles", "Count", "Freq"))

  ploidy <- dim(meta$geno_arrays)[3L]
  n_samp <- dim(meta$geno_arrays)[1L]
  # Count sums to n_samples * ploidy within every block; Freq sums to 1.
  by_blk <- sheet[, .(tot = sum(Count), fsum = sum(Freq)), by = Block]
  expect_true(all(by_blk$tot == n_samp * ploidy))
  expect_true(all(abs(by_blk$fsum - 1) < 1e-9))

  # Alleles field has one entry per member marker (semicolon-separated).
  n_mk <- lengths(strsplit(sheet$Markers, ";", fixed = TRUE))
  n_al <- lengths(strsplit(sheet$Alleles, ";", fixed = TRUE))
  expect_equal(n_mk, n_al)
})

test_that("write_block_allele_sheet: genotype counting and file write", {
  meta <- sim_geno_meta(n_samples = 10L, n_markers = 12L, n_chrom = 2L,
                        block_size = 3L, seed = 9L)
  n_samp <- dim(meta$geno_arrays)[1L]
  geno  <- write_block_allele_sheet(meta, count_by = "genotype")
  by_blk <- geno[, .(tot = sum(Count)), by = Block]
  expect_true(all(by_blk$tot <= n_samp * dim(meta$geno_arrays)[3L]))
  expect_true(all(geno[, .(fsum = sum(Freq)), by = Block]$fsum - 1 < 1e-9))

  f <- tempfile(fileext = ".txt")
  out <- write_block_allele_sheet(meta, output_file = f)
  expect_true(file.exists(f) && file.info(f)$size > 0)
  expect_equal(nrow(data.table::fread(f)), nrow(out))
})
