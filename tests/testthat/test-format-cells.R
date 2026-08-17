# =============================================================================
# test-format-cells.R
#
# QA for Kernel A (R/format_split.R + src/format_cells.cpp):
#   - split_genotype_cpp / combine_haplotype_cpp input asserts
#   - pure-R vs Rcpp parity for split_haplotypes() and combine_haplotypes(),
#     toggled via options(haploStack.use_rcpp), including the tricky R strsplit
#     edge cases (trailing/leading/internal empties), fixed-width mode, and NA.
# =============================================================================

library(testthat)
library(haploStack)

split_kernel   <- function() is.function(haploStack:::.hs_kernel("split_genotype_cpp"))
combine_kernel <- function() is.function(haploStack:::.hs_kernel("combine_haplotype_cpp"))

toggle <- function(rcpp, expr) withr::with_options(list(haploStack.use_rcpp = rcpp), expr)

# Genotype matrix exercising every split edge case R reproduces:
#   "1/4" clean; "1/" trailing-empty drop; "/4" leading empty; "1//4" internal
#   empty; "1" single token; "" empty; NA.
.edge_geno <- function() {
  matrix(c("1/4", "/4",
           "1/",  "1//4",
           "1",   "4",
           "",    "A/T",
           NA,    "9/9"),
         nrow = 2L, byrow = FALSE,
         dimnames = list(c("S1", "S2"), paste0("M", 1:5)))
}

# ── kernel input asserts ──────────────────────────────────────────────────────

test_that("split_genotype_cpp validates ploidy", {
  skip_if_not(split_kernel(), "kernel not compiled")
  k <- haploStack:::.hs_kernel("split_genotype_cpp")
  g <- matrix(c("1/4", "4/4"), 1L, 2L)
  expect_silent(k(g, "/", 2L))
  expect_error(k(g, "/", 0L), "ploidy")
})

test_that("combine_haplotype_cpp validates dims / length", {
  skip_if_not(combine_kernel(), "kernel not compiled")
  k <- haploStack:::.hs_kernel("combine_haplotype_cpp")
  hap <- c("1", "4", "4", "1")                    # 1 x 2 x 2
  expect_silent(k(hap, c(1L, 2L, 2L), "/", "-1"))
  expect_error(k(hap, c(1L, 2L), "/", "-1"), "length 3")
  expect_error(k(hap, c(1L, 2L, 3L), "/", "-1"), "prod\\(dims\\)")
})

# ── split parity (delimiter + fixed-width) ────────────────────────────────────

test_that("split_haplotypes: pure-R and Rcpp agree (delimiter, all edges)", {
  skip_if_not(split_kernel(), "kernel not compiled")
  g <- .edge_geno()
  a <- toggle(FALSE, split_haplotypes(g, split_char = "/", orientation = 1L))
  b <- toggle(TRUE,  split_haplotypes(g, split_char = "/", orientation = 1L))
  expect_identical(a$array, b$array)
  # spot-check the tricky ones against the documented R semantics
  expect_identical(unname(b$array["S2", "M1", ]), c("", "4"))            # "/4"
  expect_identical(unname(b$array["S1", "M2", ]), c("1", NA_character_)) # "1/" trailing drop
  expect_identical(unname(b$array["S2", "M2", ]), c("1", ""))            # "1//4" -> h2 ""
  expect_identical(unname(b$array["S1", "M4", ]), c(NA_character_, NA_character_)) # "" -> all NA
  expect_identical(unname(b$array["S1", "M5", ]), c(NA_character_, NA_character_)) # NA -> all NA
})

test_that("split_haplotypes: pure-R and Rcpp agree (fixed-width, sep = '')", {
  skip_if_not(split_kernel(), "kernel not compiled")
  g <- matrix(c("AT", "A", "", "AG", "TT", NA),
              nrow = 2L, dimnames = list(c("S1", "S2"), paste0("M", 1:3)))
  a <- toggle(FALSE, split_haplotypes(g, split_char = "", orientation = 1L))
  b <- toggle(TRUE,  split_haplotypes(g, split_char = "", orientation = 1L))
  expect_identical(a$array, b$array)
  expect_identical(unname(b$array["S2", "M1", ]), c("A", NA_character_)) # "A" -> h2 NA
  expect_identical(unname(b$array["S1", "M2", ]), c(NA_character_, NA_character_)) # "" -> all NA
})

# ── combine parity ────────────────────────────────────────────────────────────

test_that("combine_haplotypes: pure-R and Rcpp agree (delimiter + concat, NA->-1)", {
  skip_if_not(combine_kernel(), "kernel not compiled")
  arr <- toggle(FALSE, split_haplotypes(.edge_geno(), split_char = "/",
                                        orientation = 1L))$array
  for (sc in c("/", "")) {
    a <- toggle(FALSE, combine_haplotypes(arr, split_char = sc))
    b <- toggle(TRUE,  combine_haplotypes(arr, split_char = sc))
    expect_identical(a, b)
  }
})

test_that("split_haplotypes agrees across engines on realistic simulated data", {
  skip_if_not(split_kernel(), "kernel not compiled")
  g <- sim_geno_matrix(n_samples = 15L, n_markers = 40L, encoding = "selectiontools")
  # sim matrix is markers x samples -> orientation = 2
  a <- toggle(FALSE, split_haplotypes(g, split_char = "/", orientation = 2L))
  b <- toggle(TRUE,  split_haplotypes(g, split_char = "/", orientation = 2L))
  expect_identical(a$array, b$array)
})

# ── sorted genotype collapse parity (.collapse_genotype_sorted / kernel) ───────

collapse_kernel <- function()
  is.function(haploStack:::.hs_kernel("collapse_genotype_sorted_cpp"))

# 3-D [n_samp x n_mark x 2] array covering: sorted swap, equal alleles, one NA
# on each side, both NA, and a multi-character token pair.
.collapse_arr <- function() {
  a1 <- c("T", "A", NA,  "A", NA,  "12")
  a2 <- c("A", "A", "G", NA,  NA,  "3")
  array(c(a1, a2), dim = c(3L, 2L, 2L),
        dimnames = list(paste0("S", 1:3), paste0("M", 1:2), c("1", "2")))
}

test_that(".collapse_genotype_sorted: pure-R matches sort()+paste0 (all NA cases)", {
  arr <- .collapse_arr()
  ref <- matrix(NA_character_, 3L, 2L, dimnames = dimnames(arr)[1:2])
  for (i in 1:3) for (j in 1:2) {
    v <- arr[i, j, ]; ref[i, j] <- paste0(sort(v[!is.na(v)]), collapse = "")
  }
  expect_identical(toggle(FALSE, haploStack:::.collapse_genotype_sorted(arr)), ref)
  expect_identical(unname(ref[1, 1]), "AT")    # sorted swap "T","A" -> "AT"
  expect_identical(unname(ref[2, 1]), "AA")    # equal alleles
  expect_identical(unname(ref[3, 1]), "G")     # h1 NA -> keep present allele
  expect_identical(unname(ref[1, 2]), "A")     # h2 NA -> keep present allele
  expect_identical(unname(ref[2, 2]), "")      # both NA -> ""
  expect_identical(unname(ref[3, 2]), "123")   # multi-char sorted: "12" < "3"
})

test_that(".collapse_genotype_sorted: Rcpp and pure-R agree", {
  skip_if_not(collapse_kernel(), "kernel not compiled")
  arr <- .collapse_arr()
  a <- toggle(FALSE, haploStack:::.collapse_genotype_sorted(arr))
  b <- toggle(TRUE,  haploStack:::.collapse_genotype_sorted(arr))
  expect_identical(a, b)
  # realistic-scale agreement
  meta <- sim_geno_meta(n_samples = 20L, n_markers = 50L)
  ra <- toggle(FALSE, haploStack:::.collapse_genotype_sorted(meta$geno_arrays))
  rb <- toggle(TRUE,  haploStack:::.collapse_genotype_sorted(meta$geno_arrays))
  expect_identical(ra, rb)
})

test_that("collapse_genotype_sorted_cpp validates dims / length", {
  skip_if_not(collapse_kernel(), "kernel not compiled")
  k   <- haploStack:::.hs_kernel("collapse_genotype_sorted_cpp")
  hap <- c("1", "4", "4", "1")                     # 1 x 2 x 2
  expect_silent(k(hap, c(1L, 2L, 2L)))
  expect_error(k(hap, c(1L, 2L)),     "length 3")
  expect_error(k(hap, c(1L, 2L, 3L)), "prod\\(dims\\)")
})
