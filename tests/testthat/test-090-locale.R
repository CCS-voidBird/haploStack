# =============================================================================
# test-090-locale.R — locale-independent ordering (design D1 §3.11, step
# D1-11, decision-form row 33 a)
#
#   O1  .hs_sort() / .hs_order() give C-locale byte order whatever the
#       session collation is (and differ from sort() under a non-C locale
#       on a mixed-case set)
#   O2  the identifier orders of the package follow it: a mixed-case sample
#       set orders the same under "C" and under the session locale
#   O3  no locale-dependent sort() on identifiers is left in R/ (every
#       remaining sort()/order() call is numeric or already radix)
# =============================================================================

test_that("O1: the sorters are C-locale under any collation", {
  x <- c("b", "A", "a", "B", "S10", "S2", "s1", "Chr1", "chr1")
  want <- c("A", "B", "Chr1", "S10", "S2", "a", "b", "chr1", "s1")
  expect_identical(haploStack:::.hs_sort(x), want)
  expect_identical(x[haploStack:::.hs_order(x)], want)
  expect_identical(haploStack:::.hs_sort(x, decreasing = TRUE), rev(want))
  expect_identical(haploStack:::.hs_sort(factor(x)), want)
  for (loc in c("C", "en_US.UTF-8", "Chinese (Simplified)_China.utf8")) {
    ok <- tryCatch({ withr::local_collate(loc); TRUE }, warning = function(w) FALSE,
                   error = function(e) FALSE)
    if (!ok) next
    expect_identical(haploStack:::.hs_sort(x), want, label = loc)
  }
  # under a non-C collation base sort() orders case-insensitively: that is
  # the difference row 33 a removes (skipped where the locale is C already)
  cur <- Sys.getlocale("LC_COLLATE")
  if (!cur %in% c("C", "POSIX")) {
    expect_false(identical(sort(x), want))
  }
})

test_that("O2: identifier orders in the package are the same under C and the session locale", {
  m <- sim_geno_meta(n_samples = 6L, n_markers = 8L, n_chrom = 1L, seed = 4L,
                     encoding = "nucleotide")
  ids <- c("beta", "Alpha", "alpha", "Beta", "s10", "S2")
  old <- dimnames(m$geno_arrays)[[1L]]
  dimnames(m$geno_arrays)[[1L]] <- ids
  for (cc in names(m$sample_info))          # the sample sheet names them too
    if (identical(as.character(m$sample_info[[cc]]), old)) m$sample_info[[cc]] <- ids
  m$phenotype <- NULL
  expect_true(is_geno_meta(m))
  run <- function() {
    list(views   = geno_meta_views(m)$samples,
         palette = names(make_allele_color_map(c("t", "A", "g", "C", "a"), "nucleotide")),
         levels  = haploStack:::.hs_sort(unique(ids)))
  }
  a <- run()
  withr::with_collate("C", b <- run())
  expect_identical(a, b)
  expect_identical(a$levels, c("Alpha", "Beta", "S2", "alpha", "beta", "s10"))
})

test_that("O3: no locale-dependent identifier sort is left in R/", {
  rdir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(rdir), "package R/ sources not reachable from the test dir")
  files <- list.files(rdir, pattern = "[.]R$", full.names = TRUE)
  # sort() calls that are not .hs_sort(): each must be on a numeric quantity
  # or an integer index (reviewed list, D1-11)
  allowed <- c(
    "alphabet.R:sort(unique(lens))",                      # integer lengths
    "catalogue_circos_meta.R:sort(kept)",                 # integer indices
    "simulate_dataset.R:sort(rep(seq_len",                # integers
    "simulate_dataset.R:sort(sample(seq(",                # integers
    "stacking_gbs.R:sort(match(",                         # integers
    "stacking_gbs.R:sort(sample.int(",                    # integers
    "stacking_gbs.R:sort(sample(union(",                  # integers
    "stacking_gbs.R:sort(child)",                         # integers
    "stacking_main.R:sort(r, decreasing",                 # numeric
    "stacking_main.R:sort(v, decreasing",                 # numeric
    "stacking_main.R:sort(step_time",                     # numeric
    "stacking_ohs.R:sort(r, decreasing",                  # numeric
    "stacking_ohs.R:sort(v, decreasing",                  # numeric
    "stacking_ohs.R:sort(sample.int(",                    # integers
    "stacking_ohs.R:sort(sample(union(",                  # integers
    "stacking_ohs.R:sort(child)",                         # integers
    "stacking_pairs.R:sort(unique(as.integer(",           # integers
    "stacking_pairs.R:sort(Effect[",                      # numeric
    "stacking_pairs.R:sort(rep(e,",                       # numeric
    "stacking_pairs.R:sort(keep)",                        # integers
    "stacking_pairs.R:sort(x, decreasing",                # numeric
    "stacking_selection.R:sort(adjusted_gebv",            # numeric
    "state.R:sort(x, decreasing = decreasi"               # the helper itself
  )
  hits <- character(0L)
  for (f in files) {
    src <- readLines(f, warn = FALSE)
    src <- sub("#.*$", "", src)                          # strip comments
    i <- grep("(^|[^._a-zA-Z])sort\\(", src)
    for (k in i) {
      frag <- regmatches(src[k], regexpr("sort\\([^\n]{0,24}", src[k]))
      key  <- paste0(basename(f), ":", frag)
      if (!any(startsWith(key, allowed))) hits <- c(hits, paste0(key, " (line ", k, ")"))
    }
  }
  expect_identical(hits, character(0L))
})
