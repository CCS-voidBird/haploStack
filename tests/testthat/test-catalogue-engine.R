# =============================================================================
# test-catalogue-engine.R
#
# QA for the catalogue compute engine (R/catalogue_engine.R + src/local_gebv.cpp,
# catalogue-perf upgrade):
#   - .catalogue_core() invariants
#   - Rcpp kernel local_gebv_block_sum_cpp input asserts
#   - the global haploStack.use_rcpp toggle + pure-R fallback
#   - pure-R vs Rcpp parity (structure exact, continuous GEBV within tolerance)
#     across compute_profile_localgebv(), compute_block_catalogue(), and
#     run_haplotype_catalogue().
# =============================================================================

library(testthat)
library(haploStack)

kernel_available <- function() is.function(haploStack:::.hs_kernel("local_gebv_block_sum_cpp"))

# Compare two wide/long catalogue tables: structure (dims, names, non-numeric
# columns) must be IDENTICAL; numeric columns may differ only by `tol` (the Rcpp
# kernel sums per block in a different order than the pure-R per-copy grouping).
expect_engine_parity <- function(a, b, tol = 1e-8) {
  a <- as.data.table(a); b <- as.data.table(b)
  expect_identical(dim(a), dim(b))
  expect_identical(names(a), names(b))
  numc <- names(a)[vapply(a, is.numeric, logical(1L))]
  for (cc in setdiff(names(a), numc)) expect_identical(a[[cc]], b[[cc]])
  for (cc in numc)
    expect_lt(max(abs(a[[cc]] - b[[cc]]), na.rm = TRUE), tol)
}

.eng_meta <- sim_geno_meta(n_samples = 8L, n_markers = 18L,
                           n_chrom = 2L, block_size = 3L, seed = 21L)
.eng_eff  <- sim_effect_df(dimnames(.eng_meta$geno_arrays)[[2L]],
                           alleles = c("A", "T"), seed = 21L)

# ── .catalogue_core ───────────────────────────────────────────────────────────

test_that(".catalogue_core builds the trait-invariant long table + frequencies", {
  core <- haploStack:::.catalogue_core(.eng_meta, use_rcpp = FALSE)
  expect_false(core$empty)
  expect_true(all(c("Sample", "Haplotype_idx", "Marker", "Allele", "Block")
                  %in% names(core$long_core)))
  expect_true(all(c("Marker", "Allele", "AlleleFreq") %in% names(core$freq_core)))
  # frequencies sum to 1 within each marker
  s <- core$freq_core[, .(tot = sum(AlleleFreq)), by = Marker]
  expect_true(all(abs(s$tot - 1) < 1e-9))
})

test_that(".catalogue_core builds integer-coded structures when use_rcpp=TRUE", {
  core <- haploStack:::.catalogue_core(.eng_meta, use_rcpp = TRUE)
  expect_true(is.matrix(core$code_mat))
  expect_equal(ncol(core$code_mat), length(core$block_idx))
  expect_true(all(core$block_idx >= 0L))
  expect_true(all(core$block_idx <= length(core$block_levels)))
  # codes are 1-based within each marker's level count (or NA)
  expect_true(all(core$code_mat >= 1L, na.rm = TRUE))
})

test_that(".catalogue_core flags empty input", {
  m0 <- sim_geno_meta(n_samples = 4L, n_markers = 6L, block_size = NULL,
                      seed = 5L)                       # no blocks attached
  core <- haploStack:::.catalogue_core(m0, use_rcpp = FALSE)
  expect_true(isTRUE(core$empty))
})

# ── Rcpp kernel input asserts ─────────────────────────────────────────────────

test_that("local_gebv_block_sum_cpp validates its inputs", {
  skip_if_not(kernel_available(), "Rcpp kernel not compiled")
  k <- haploStack:::.hs_kernel("local_gebv_block_sum_cpp")
  code    <- matrix(1L, 2L, 3L)
  contrib <- matrix(0,  3L, 2L)
  blk     <- c(1L, 1L, 2L)
  expect_silent(k(code, contrib, blk, 2L))
  # mismatched marker count
  expect_error(k(code, matrix(0, 4L, 2L), blk, 2L), "markers")
  expect_error(k(code, contrib, c(1L, 1L), 2L), "length")
  # out-of-range block index / allele code
  expect_error(k(code, contrib, c(1L, 1L, 9L), 2L), "out of range")
  expect_error(k(matrix(5L, 2L, 3L), contrib, blk, 2L), "out of range")
})

test_that("kernel skips NA codes and 0 block indices", {
  skip_if_not(kernel_available(), "Rcpp kernel not compiled")
  k <- haploStack:::.hs_kernel("local_gebv_block_sum_cpp")
  code    <- matrix(c(1L, NA_integer_, 2L, 1L), nrow = 2L)  # 2 copies x 2 markers
  contrib <- matrix(c(10, 20, 30, 40), nrow = 2L)           # 2 markers x 2 codes
  res <- k(code, contrib, c(1L, 0L), 1L)                    # marker 2 skipped
  # copy1: marker1 code1 -> contrib(1,1)=10 ; copy2: marker1 NA -> 0
  expect_equal(as.numeric(res), c(10, 0))
})

# ── use_rcpp toggle ───────────────────────────────────────────────────────────

test_that(".hs_use_rcpp honours the option and kernel availability", {
  withr::local_options(haploStack.use_rcpp = FALSE)
  expect_false(haploStack:::.hs_use_rcpp("local_gebv_block_sum_cpp"))
  withr::local_options(haploStack.use_rcpp = TRUE)
  expect_equal(haploStack:::.hs_use_rcpp("local_gebv_block_sum_cpp"), kernel_available())
})

# ── pure-R vs Rcpp parity ─────────────────────────────────────────────────────

test_that("compute_profile_localgebv: pure-R and Rcpp agree (genotype)", {
  skip_if_not(kernel_available(), "Rcpp kernel not compiled")
  a <- withr::with_options(list(haploStack.use_rcpp = FALSE),
         compute_profile_localgebv(.eng_meta, .eng_eff, center = "allele_freq"))
  b <- withr::with_options(list(haploStack.use_rcpp = TRUE),
         compute_profile_localgebv(.eng_meta, .eng_eff, center = "allele_freq"))
  expect_engine_parity(a, b)
})

test_that("compute_profile_localgebv: pure-R and Rcpp agree (haplotype, no centre)", {
  skip_if_not(kernel_available(), "Rcpp kernel not compiled")
  a <- withr::with_options(list(haploStack.use_rcpp = FALSE),
         compute_profile_localgebv(.eng_meta, .eng_eff, level = "haplotype",
                                   center = "none"))
  b <- withr::with_options(list(haploStack.use_rcpp = TRUE),
         compute_profile_localgebv(.eng_meta, .eng_eff, level = "haplotype",
                                   center = "none"))
  expect_engine_parity(a, b)
})

test_that("compute_block_catalogue: pure-R and Rcpp agree", {
  skip_if_not(kernel_available(), "Rcpp kernel not compiled")
  a <- withr::with_options(list(haploStack.use_rcpp = FALSE),
         compute_block_catalogue(.eng_meta, .eng_eff))
  b <- withr::with_options(list(haploStack.use_rcpp = TRUE),
         compute_block_catalogue(.eng_meta, .eng_eff))
  expect_engine_parity(a, b)
})

test_that("run_haplotype_catalogue: pure-R and Rcpp agree on every table", {
  skip_if_not(kernel_available(), "Rcpp kernel not compiled")
  out_r <- withr::local_tempdir(); out_c <- withr::local_tempdir()
  withr::with_options(list(haploStack.use_rcpp = FALSE),
    run_haplotype_catalogue(.eng_meta, list(T1 = .eng_eff), output_dir = out_r,
                            save_profile = FALSE))
  withr::with_options(list(haploStack.use_rcpp = TRUE),
    run_haplotype_catalogue(.eng_meta, list(T1 = .eng_eff), output_dir = out_c,
                            save_profile = FALSE))
  for (f in list.files(out_r, pattern = "[.]txt$"))
    expect_engine_parity(data.table::fread(file.path(out_r, f)),
                         data.table::fread(file.path(out_c, f)))
})

# ── shared-core reuse equals the standalone path ──────────────────────────────

test_that("passing a shared .core matches building one internally", {
  withr::local_options(haploStack.use_rcpp = FALSE)
  core <- haploStack:::.catalogue_core(.eng_meta, use_rcpp = FALSE)
  a <- compute_profile_localgebv(.eng_meta, .eng_eff, center = "allele_freq")
  b <- compute_profile_localgebv(.eng_meta, .eng_eff, center = "allele_freq",
                                 .core = core)
  expect_identical(a, b)
  ca <- compute_block_catalogue(.eng_meta, .eng_eff)
  cb <- compute_block_catalogue(.eng_meta, .eng_eff, .core = core)
  expect_identical(ca, cb)
})

test_that("shared .core + .signatures match the standalone signature path", {
  withr::local_options(haploStack.use_rcpp = FALSE)
  core <- haploStack:::.catalogue_core(.eng_meta, use_rcpp = FALSE)
  sigs <- haploStack:::.compute_block_signatures(.eng_meta)
  a <- compute_block_catalogue(.eng_meta, .eng_eff, group_by = "signature")
  b <- compute_block_catalogue(.eng_meta, .eng_eff, group_by = "signature",
                               .core = core, .signatures = sigs)
  expect_identical(a, b)
})

test_that("analyze_block_catalogue agrees across engines (shared core path)", {
  skip_if_not(kernel_available(), "Rcpp kernel not compiled")
  eff2 <- list(T1 = .eng_eff,
               T2 = sim_effect_df(dimnames(.eng_meta$geno_arrays)[[2L]],
                                  alleles = c("A", "T"), seed = 99L))
  ar <- suppressWarnings(withr::with_options(list(haploStack.use_rcpp = FALSE),
          analyze_block_catalogue(.eng_meta, eff2, plots = character(0))))
  ac <- suppressWarnings(withr::with_options(list(haploStack.use_rcpp = TRUE),
          analyze_block_catalogue(.eng_meta, eff2, plots = character(0))))
  expect_engine_parity(ar$selected, ac$selected)
  expect_engine_parity(ar$S, ac$S)
  expect_engine_parity(ar$N, ac$N)
})

test_that(".core cannot be combined with sample/marker/chromosome filters", {
  core <- haploStack:::.catalogue_core(.eng_meta, use_rcpp = FALSE)
  expect_error(
    compute_profile_localgebv(.eng_meta, .eng_eff, samples = "S1", .core = core),
    "cannot be combined")
})
