# =============================================================================
# test-maintenance-0.9.0.R
#
# Regression tests for the 0.9.0 maintenance fixes that need no owner decision
# (.scratch/maintenance-0.9.0/DECISION_FORM.md, "Done without asking"):
#   M1  write_geno_meta() -> read_plink_ped() keeps missing calls missing
#   M2  read_genotype() / read_phenotype() keep identifier spelling ("007")
#   M3  recode_geno_meta() forwards the phenotype (discrepancy 20)
#   M4  export_block_annotation() on a custom-only container names the fix
#   M5  .left_join_dt() has dplyr::left_join() semantics (dplyr dropped)
# =============================================================================

# Nucleotide container with a fully missing call, a half call and an ID that
# a numeric guess would rewrite ("007", marker "0012").
.m090_container <- function() {
  s <- c("S1", "S2", "007"); m <- c("M1", "M2", "M3", "0012")
  h1 <- matrix(c("A", "A", NA,  "T",
                 "T", NA,  "A", "A",
                 "A", "T", "T", "A"), nrow = 3, byrow = TRUE,
               dimnames = list(s, m))
  h2 <- matrix(c("A", "T", NA,  "T",
                 "T", NA,  "T", "A",
                 "T", "T", NA,  "A"), nrow = 3, byrow = TRUE,
               dimnames = list(s, m))
  map <- data.frame(Chrom = c("1", "1", "2", "2"), Marker = m,
                    Pos = c(1e6, 2e6, 3e6, 4e6), stringsAsFactors = FALSE)
  arr <- array(c(h1, h2), dim = c(dim(h1), 2L),
               dimnames = c(dimnames(h1), list(c("1", "2"))))
  suppressMessages(new_geno_meta(
    arr, map,
    metadata = list(ploidy = 2L,
                    char_dict = stats::setNames(as.list(c(1L, 4L)), c("A", "T")),
                    n_samples = 3L, n_markers = 4L,
                    source_format = "hand_built", encoding = "nucleotide")))
}

test_that("M1: a PED written by write_geno_meta() reads back with the same missing cells", {
  gm <- .m090_container()
  stem <- tempfile("m090_ped_")
  on.exit(unlink(paste0(stem, c(".ped", ".map"))), add = TRUE)
  suppressMessages(write_geno_meta(gm, stem, preset = "PL"))

  # Missing alleles are written as PLINK's 0, never as the text NA.
  ped <- readLines(paste0(stem, ".ped"))
  expect_false(any(grepl("\\bNA\\b", ped)))

  back <- suppressMessages(suppressWarnings(
    read_plink_ped(paste0(stem, ".ped"), paste0(stem, ".map"))))
  a0 <- gm$geno_arrays; a1 <- back$geno_arrays[rownames(a0), colnames(a0), ]
  dimnames(a1)[[3L]] <- dimnames(a0)[[3L]]
  # Fully missing (S1, M3) stays fully missing: up to 0.8.8 it came back (NA, "A").
  expect_true(all(is.na(a1["S1", "M3", ])))
  expect_identical(is.na(a1), is.na(a0))
  expect_identical(a1[!is.na(a1)], a0[!is.na(a0)])
})

test_that("M2: read_genotype() and read_phenotype() keep ID spelling", {
  gf <- tempfile(fileext = ".txt"); mf <- tempfile(fileext = ".txt")
  pf <- tempfile(fileext = ".txt")
  on.exit(unlink(c(gf, mf, pf)), add = TRUE)
  writeLines(c("ID\tS1\t007",
               "0012\t1/1\t1/4",
               "M2\t4/4\t1/1"), gf)
  writeLines(c("Marker\tChrom\tPos",
               "0012\t01\t1000000",
               "M2\t01\t2000000"), mf)
  gm <- suppressMessages(suppressWarnings(
    read_genotype(gf, map_file = mf, preset = "ST", st_formats = TRUE)))
  expect_setequal(dimnames(gm$geno_arrays)[[1L]], c("S1", "007"))
  expect_setequal(dimnames(gm$geno_arrays)[[2L]], c("0012", "M2"))
  expect_setequal(gm$map_info$Marker, c("0012", "M2"))
  expect_true(all(gm$map_info$Chrom == "01"))
  expect_true(is.numeric(gm$map_info$Pos))

  writeLines(c("Sample\tYield\tHeight", "007\t1.5\t90", "S1\tNA\t85"), pf)
  ph <- read_phenotype(pf)
  expect_identical(ph$Sample, c("007", "S1"))
  expect_identical(ph$Yield, c(1.5, NA))
  expect_identical(ph$Height, c(90, 85))
})

test_that("M3: recode_geno_meta() keeps the phenotype", {
  gm <- .m090_container()
  ph <- data.frame(Sample = dimnames(gm$geno_arrays)[[1L]], Yield = c(1, 2, 3))
  gp <- suppressMessages(attach_phenotype(gm, ph))
  rc <- suppressMessages(recode_geno_meta(gp, to_type = "nucleotide"))
  expect_identical(rc$phenotype, gp$phenotype)
  expect_true(is_geno_meta(rc))
})

test_that("M4: export_block_annotation() refuses a custom-only container with the fix", {
  gm <- .m090_container()
  df <- data.frame(Marker = c("M1", "0012"), Name = c("L1", "L2"))
  m1 <- suppressMessages(attach_marker_annotations(gm, df, source = "unit test"))
  expect_identical(names(m1$annotation), "custom")
  od <- tempfile("m090_ann_")
  on.exit(unlink(od, recursive = TRUE), add = TRUE)
  expect_error(export_block_annotation(m1, od), "attach_block_annotation")
})

test_that("M5: .left_join_dt() keeps x order, repeats multi-matches and suffixes shared names", {
  x <- data.frame(k = c("b", "a", "z", "a"), v = 1:4, s = "x")
  y <- data.frame(k = c("a", "b", "a"), w = c(10, 20, 30), s = "y")
  out <- haploStack:::.left_join_dt(x, y, by = "k")
  expect_identical(names(out), c("k", "v", "s.x", "w", "s.y"))
  expect_identical(out$k, c("b", "a", "a", "z", "a", "a"))
  expect_identical(out$v, c(1L, 2L, 2L, 3L, 4L, 4L))
  expect_identical(out$w, c(20, 10, 30, NA, 10, 30))
  # inputs are not modified by reference
  xd <- data.table::as.data.table(x)
  invisible(haploStack:::.left_join_dt(xd, y, by = "k"))
  expect_identical(names(xd), c("k", "v", "s"))
})

test_that("M6: package conditions are classed; text and call are unchanged", {
  # error: class, message and the calling function as the call
  e <- tryCatch(read_phenotype(tempfile("no_such_")), error = function(e) e)
  expect_s3_class(e, "hs_error")
  expect_match(conditionMessage(e), "^read_phenotype: file not found: ")
  expect_identical(as.character(conditionCall(e)[[1L]]), "read_phenotype")
  # message: class, and hs_quiet() muffles it while other messages pass
  f <- function() { haploStack:::.hs_inform("from haploStack"); message("from R") }
  expect_message(f(), "from haploStack", class = "hs_message")
  msgs <- character()
  withCallingHandlers(hs_quiet(f()), message = function(m) {
    msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
  })
  expect_identical(msgs, "from R\n")
  # warning: class
  expect_warning(haploStack:::.hs_warn("w"), "w", class = "hs_warning")
  # a condition object passes through untouched
  cnd <- errorCondition("custom", class = "my_error")
  expect_error(haploStack:::.hs_abort(cnd), class = "my_error")
})

test_that("M7: compiled coancestry kernel equals the R loop bit for bit", {
  skip_if_not(haploStack:::.hs_use_rcpp("block_ibs_coancestry_cpp"))
  gm <- sim_geno_meta(n_samples = 10, n_markers = 36, block_size = 4, seed = 5)
  set.seed(9)
  a <- gm$geno_arrays
  a[runif(length(a)) < 0.1] <- NA          # missing cells, incl. whole blocks
  a[1, 1:4, ] <- NA
  gm$geno_arrays <- a
  kc <- compute_parent_coancestry(gm)
  withr::local_options(haploStack.use_rcpp = FALSE)
  kr <- compute_parent_coancestry(gm)
  expect_identical(kc, kr)
  blk <- unique(gm$block_marker_map_info$Block)[c(1, 3)]
  withr::local_options(haploStack.use_rcpp = TRUE)
  kc2 <- compute_parent_coancestry(gm, samples = rownames(kc)[c(4, 2, 9)], blocks = blk)
  withr::local_options(haploStack.use_rcpp = FALSE)
  kr2 <- compute_parent_coancestry(gm, samples = rownames(kc)[c(4, 2, 9)], blocks = blk)
  expect_identical(kc2, kr2)
})

test_that("M8: OHS pair matrix equals the per-pair top-2 rule; NA is an absent copy", {
  set.seed(3)
  cols <- as.vector(t(outer(sprintf("S%02d", 1:7), 1:2, paste, sep = ".")))
  v <- matrix(round(rnorm(25 * length(cols)), 1), 25)   # rounding makes ties
  dt <- data.table::data.table(Block = sprintf("b%02d", 1:25))
  dt[, (cols) := as.data.frame(v)]
  r <- haploStack:::.pairwise_block_value_matrix(dt, return_array = TRUE)
  A <- r$A
  ref <- matrix(0, 7, 7)
  for (i in 1:6) for (j in (i + 1):7) {
    ref[i, j] <- ref[j, i] <- sum(apply(cbind(A[i, , ], A[j, , ]), 1L,
                                        function(x) sum(sort(x, decreasing = TRUE)[1:2])))
  }
  expect_identical(unname(r$M), ref)
  # an NA no longer crashes: (NA, 3) -> (3, 0); (NA, NA) -> (0, 0)
  d <- data.table::data.table(Block = c("b1", "b2"), S1.1 = c(1, NA), S1.2 = c(2, 3),
                              S2.1 = c(1, NA), S2.2 = c(0, NA))
  r2 <- haploStack:::.pairwise_block_value_matrix(d, return_array = TRUE)
  expect_identical(r2$A[, 2L, ], matrix(c(3, 0, 0, 0), 2))
  expect_identical(r2$M[1, 2], 2 + 1 + 3 + 0)
})
