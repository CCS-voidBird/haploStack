# Extracted from test-090-simulation.R:49

# prequel ----------------------------------------------------------------------
library(testthat)
library(haploStack)
library(data.table)
.e_fixture <- function(seed = 5L, n_samples = 12L, n_markers = 24L) {
  m <- sim_geno_meta(n_samples = n_samples, n_markers = n_markers, n_chrom = 2L,
                     block_size = 3L, seed = seed)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("A", "T"))
  list(meta = m, eff = e,
       lg = compute_profile_localgebv(m, e, level = "haplotype"))
}

# test -------------------------------------------------------------------------
reg <- haploStack:::.STK_SCHEMES
expect_named(reg, c("TS", "OHS", "OCS", "GBS"))
for (s in names(reg)) {
    expect_true(is.function(reg[[s]]$parents), info = s)
    expect_true(is.null(reg[[s]]$pairs) || is.function(reg[[s]]$pairs), info = s)
    expect_true(is.logical(reg[[s]]$needs_ideal), info = s)
    expect_true(reg[[s]]$direction %in% c("both", "high", "argmax"), info = s)
  }
expect_null(reg$OCS$pairs)
expect_true(reg$GBS$needs_ideal)
expect_true(reg$GBS$pin_founder_blocks)
expect_identical(haploStack:::.STK_BASELINE, "TS")
expect_identical(haploStack:::.stk_scheme_flag(c("TS", "GBS"), "needs_ideal"),
                   c(FALSE, TRUE))
expect_error(haploStack:::.stk_scheme("XYZ"), class = "hs_error_internal")
src <- testthat::test_path("..", "..", "R", "stacking_main.R")
skip_if_not(file.exists(src), "R/ not reachable")
code <- readLines(src, warn = FALSE)
code <- code[!grepl("^\\s*#", code)]
expect_false(any(grepl('scheme == "(TS|OHS|OCS|GBS)"', code)))
