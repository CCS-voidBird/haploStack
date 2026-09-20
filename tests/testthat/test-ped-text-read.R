# =============================================================================
# test-ped-text-read.R
#
# 0.7.7 — a genotype file is read as TEXT.
#
# utils::read.table() runs type.convert() per column, and type.convert() maps a
# column whose every value is "T" (or "F") to LOGICAL. In a nucleotide PED that
# is an ordinary allele column — a T/A or T/G marker that is monomorphic-T in
# the panel — and after mode(x) <- "character" the cells read "TRUE", which
# .reshape_wide_alleles() rejects. The whole file then fails to read.
#
# read_plink_ped() now goes through .read_table_chr(), which reads every column
# as character and leaves numeric coercion to the caller. These tests pin:
#   T1  the monomorphic-T column that used to abort the read
#   T2  the F/T pair, the other type.convert() trap
#   T3  identifiers are preserved verbatim (leading zeros, "TRUE" as a name)
#   T4  a clean panel is unchanged by the new reader
#   T5  the primitive itself
#
# Spec: INDEX discrepancy 30.
# =============================================================================

library(testthat)
library(haploStack)

# -- helper -------------------------------------------------------------------

# Write a PED/MAP pair from an explicit (n_samples x 2*n_markers) allele matrix,
# so a test can place the exact column pattern it wants to exercise.
write_ped_from <- function(alleles, markers = NULL, samples = NULL,
                           fam = NULL) {
  n_samples <- nrow(alleles)
  n_markers <- ncol(alleles) / 2L
  if (is.null(markers)) markers <- paste0("M", seq_len(n_markers))
  if (is.null(samples)) samples <- paste0("S", seq_len(n_samples))
  if (is.null(fam))     fam     <- paste0("F", seq_len(n_samples))

  ped <- cbind(
    data.frame(FID = fam, IID = samples, PAT = "0", MAT = "0",
               SEX = "0", PHENO = "-9", stringsAsFactors = FALSE),
    as.data.frame(alleles, stringsAsFactors = FALSE)
  )
  ped_file <- tempfile(fileext = ".ped")
  map_file <- tempfile(fileext = ".map")
  utils::write.table(ped, ped_file, quote = FALSE, sep = " ",
                     row.names = FALSE, col.names = FALSE)
  utils::write.table(
    data.frame(Chrom = "1", Marker = markers, cM = 0,
               Pos = seq_len(n_markers) * 1000000L, stringsAsFactors = FALSE),
    map_file, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
  list(ped = ped_file, map = map_file, markers = markers, samples = samples)
}

# -- T1: the column that used to abort the read -------------------------------

test_that("T1: a monomorphic-T marker does not break the PED read", {
  n <- 12L
  # marker 1 is T/T for every line (both allele columns are all "T"), marker 2
  # segregates, marker 3 is missing for half the panel.
  al <- cbind(
    rep("T", n), rep("T", n),                       # M1 — the trap
    rep(c("A", "T"), length.out = n), rep("A", n),  # M2
    c(rep("0", n / 2L), rep("G", n / 2L)),          # M3 allele 1
    c(rep("0", n / 2L), rep("T", n / 2L))           # M3 allele 2
  )
  fx   <- write_ped_from(al)
  meta <- suppressMessages(read_plink_ped(fx$ped, fx$map))

  expect_true(is_geno_meta(meta))
  expect_equal(dim(meta$geno_arrays), c(n, 3L, 2L))

  # The whole point: M1 is the allele T on both copies, for every sample —
  # not "TRUE", not NA.
  expect_true(all(meta$geno_arrays[, "M1", ] == "T"))
  expect_equal(meta$metadata$encoding, "nucleotide")
  expect_setequal(names(meta$metadata$char_dict), c("A", "G", "T"))

  # PLINK's "0" is still a missing call, not a fourth allele.
  expect_true(all(is.na(meta$geno_arrays[seq_len(n / 2L), "M3", ])))
  expect_false("0" %in% names(meta$metadata$char_dict))
})

# -- T2: F is the other half of the trap --------------------------------------

test_that("T2: an all-F column and a mixed T/F column read as alleles", {
  # "F" is not a base, but type.convert() treats T/F as logical either way, and
  # a 1/2/3/4 or A/B/C-coded panel can legitimately contain either letter. The
  # reader must not care what the letters mean.
  n  <- 8L
  al <- cbind(rep("F", n), rep("F", n),
              rep(c("T", "F"), length.out = n),
              rep(c("F", "T"), length.out = n))
  fx   <- write_ped_from(al, markers = c("M1", "M2"))
  meta <- suppressMessages(read_plink_ped(fx$ped, fx$map))

  expect_true(all(meta$geno_arrays[, "M1", ] == "F"))
  expect_setequal(unique(as.vector(meta$geno_arrays[, "M2", ])), c("T", "F"))
  expect_false(any(as.vector(meta$geno_arrays) %in% c("TRUE", "FALSE")))
})

# -- T3: identifiers survive verbatim ----------------------------------------

test_that("T3: sample and marker identifiers are not type-converted", {
  n  <- 4L
  al <- cbind(rep("A", n), rep("T", n))
  fx <- write_ped_from(
    al,
    markers = "M1",
    # a numeric-looking ID with a leading zero, and one that reads as a boolean
    samples = c("007", "0012", "TRUE", "S4"),
    fam     = c("01", "01", "02", "02")
  )
  meta <- suppressMessages(read_plink_ped(fx$ped, fx$map))

  expect_equal(dimnames(meta$geno_arrays)[[1L]],
               c("007", "0012", "TRUE", "S4"))
  expect_equal(meta$sample_info$Sample, c("007", "0012", "TRUE", "S4"))
  expect_equal(meta$sample_info$Family, c("01", "01", "02", "02"))
})

test_that("T3b: a numeric-looking marker name keeps its exact spelling", {
  n  <- 4L
  al <- cbind(rep("A", n), rep("T", n), rep("C", n), rep("C", n))
  fx <- write_ped_from(al, markers = c("0100", "11_10694"))
  meta <- suppressMessages(read_plink_ped(fx$ped, fx$map))

  expect_equal(meta$map_info$Marker, c("0100", "11_10694"))
  expect_equal(dimnames(meta$geno_arrays)[[2L]], c("0100", "11_10694"))
  # Pos is still numeric — the caller coerces it explicitly.
  expect_type(meta$map_info$Pos, "double")
  expect_equal(meta$map_info$Pos, c(1e6, 2e6))
})

# -- T4: a clean panel is unaffected -----------------------------------------

test_that("T4: a panel with no trap column reads exactly as before", {
  set.seed(7L)
  n <- 20L; m <- 5L
  a1 <- matrix(sample(c("A", "C"), n * m, replace = TRUE), nrow = n)
  a2 <- matrix(sample(c("A", "C"), n * m, replace = TRUE), nrow = n)
  al <- matrix(NA_character_, nrow = n, ncol = 2L * m)
  al[, seq(1L, 2L * m, by = 2L)] <- a1
  al[, seq(2L, 2L * m, by = 2L)] <- a2

  fx   <- write_ped_from(al)
  meta <- suppressMessages(read_plink_ped(fx$ped, fx$map))

  expect_equal(dim(meta$geno_arrays), c(n, m, 2L))
  expect_equal(unname(meta$geno_arrays[, , 1L]), unname(a1))
  expect_equal(unname(meta$geno_arrays[, , 2L]), unname(a2))
  expect_equal(meta$metadata$ploidy, 2L)
  expect_equal(meta$metadata$source_format, "plink_ped")
})

test_that("T4b: a 1/2/3/4-coded PED is still declared selectiontools", {
  n  <- 6L
  al <- cbind(rep("1", n), rep("1", n), rep(c("1", "4"), length.out = n),
              rep("4", n))
  fx   <- write_ped_from(al, markers = c("M1", "M2"))
  meta <- suppressMessages(read_plink_ped(fx$ped, fx$map))

  expect_equal(meta$metadata$encoding, "selectiontools")
  expect_true(all(meta$geno_arrays[, "M1", ] == "1"))
})

# -- T6: read_effect_file() — the silent one, and its first tests ------------

# This reader had no test coverage at all before 0.7.7, which is how an allele
# column could become TRUE and zero out every effect without anything noticing.

test_that("T6: an all-T Allele column stays the allele T, not TRUE", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("M1 T 0.4", "M2 T -0.2", "M3 T 0.9"), f)
  eff <- read_effect_file(f)

  expect_equal(as.character(eff$Allele), rep("T", 3L))
  expect_false(any(as.character(eff$Allele) %in% c("TRUE", "FALSE")))
  expect_type(eff$Effect, "double")
  expect_equal(eff$Effect, c(0.4, -0.2, 0.9))
  expect_equal(as.character(eff$Marker), c("M1", "M2", "M3"))
})

test_that("T6b: the zeroed-GEBV consequence is gone", {
  # The failure this guards: an effect table keyed on TRUE matches no genotype,
  # and .compute_effect_size() scores an absent allele as 0 — so every local
  # GEBV came out 0 with nothing raised at read time.
  f <- tempfile(fileext = ".txt")
  writeLines(c("M1 T 0.4", "M2 T -0.2", "M3 T 0.9"), f)

  geno <- matrix(c("T/T", "T/A", "A/A",
                   "T/T", "A/A", "T/A",
                   "T/A", "T/T", "A/A"), nrow = 3L, byrow = TRUE,
                 dimnames = list(c("S1", "S2", "S3"), c("M1", "M2", "M3")))
  map  <- data.frame(Chrom = "1", Marker = c("M1", "M2", "M3"),
                     Pos = c(1e6, 2e6, 3e6), stringsAsFactors = FALSE)
  meta <- suppressMessages(build_geno_meta(
    geno_2d = geno, map = map, ploidy = 2L, split_char = "/",
    encoding = "nucleotide", source_format = "nucleotide"))
  meta <- attach_block_info(meta, data.frame(
    Chrom = "1", Marker = c("M1", "M2", "M3"), Block = "b1",
    stringsAsFactors = FALSE))

  eff <- read_effect_file(f)
  lg  <- suppressMessages(compute_profile_localgebv(
    meta, effect_df = as.data.frame(eff[, c("Marker", "Allele", "Effect")]),
    level = "genotype"))
  vals <- unlist(as.data.frame(lg)[, c("S1", "S2", "S3")], use.names = FALSE)

  # The alleles now match, so the block is not uniformly zero.
  expect_false(all(vals == 0))
  expect_true(any(vals != 0))
})

test_that("T6c: a mixed-allele effect file is read the same as it always was", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("M1 A 0.1", "M1 T -0.1", "M2 C 0.3", "M2 G -0.3"), f)
  eff <- read_effect_file(f)

  expect_setequal(as.character(eff$Allele), c("A", "T", "C", "G"))
  expect_equal(eff$Effect[order(eff$Marker, eff$Allele)],
               c(0.1, -0.1, 0.3, -0.3)[order(c("M1", "M1", "M2", "M2"),
                                             c("A", "T", "C", "G"))])
})

# -- T5: the primitive --------------------------------------------------------

test_that("T5: .read_table_chr returns every column as character", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("1 T 007 3.5 TRUE", "2 T 008 4.5 FALSE"), f)
  out <- haploStack:::.read_table_chr(f, header = FALSE)

  expect_equal(ncol(out), 5L)
  expect_true(all(vapply(out, is.character, logical(1L))))
  expect_equal(out[[2L]], c("T", "T"))
  expect_equal(out[[3L]], c("007", "008"))
  expect_equal(out[[5L]], c("TRUE", "FALSE"))
})

test_that("T5b: .read_table_chr errors clearly on a missing or empty file", {
  expect_error(haploStack:::.read_table_chr(tempfile()), "file not found")

  empty <- tempfile(fileext = ".txt")
  file.create(empty)
  expect_error(haploStack:::.read_table_chr(empty))
})

# -- T7: ragged input must never lose a line ---------------------------------

# fread(fill = FALSE) does NOT reliably error on a ragged file: given a first
# line wider than the rest it settles on the narrower width and SILENTLY DROPS
# the wider row. That is worse than the defect this whole change fixes, so the
# reader pads instead and checks the row count itself.

test_that("T7: a ragged file never silently loses a row", {
  f <- tempfile(fileext = ".txt")
  # line 1 is wider than line 2 — the shape that made fread drop a row
  writeLines(c("F1 S1 0 0 0 -9 A T C G",
               "F2 S2 0 0 0 -9 A T"), f)

  expect_error(haploStack:::.read_table_chr(f, strict = TRUE),
               "shorter than the widest line")

  # non-strict: both rows survive, the short one padded — nothing is discarded
  out <- haploStack:::.read_table_chr(f, strict = FALSE)
  expect_equal(nrow(out), 2L)
  expect_equal(out[[2L]], c("S1", "S2"))
  expect_true(is.na(out[[ncol(out)]][2L]))
})

test_that("T7b: a truncated PED line is refused, not padded into missing calls", {
  ped <- tempfile(fileext = ".ped")
  map <- tempfile(fileext = ".map")
  writeLines(c("F1 S1 0 0 0 -9 A T C G",
               "F2 S2 0 0 0 -9 A T"), ped)          # one whole marker short
  writeLines(c("1 M1 0 1000000", "1 M2 0 2000000"), map)

  # The danger: the short line is a clean multiple of ploidy, so the ploidy
  # check cannot catch it. Before this guard the sample was dropped outright
  # and read_plink_ped() happily returned a 1-sample, ploidy-1 container.
  expect_error(suppressWarnings(read_plink_ped(ped, map)),
               "rectangular|shorter than the widest line")
})

test_that("T7c: a well-formed PED with a blank trailing line still reads", {
  ped <- tempfile(fileext = ".ped")
  map <- tempfile(fileext = ".map")
  writeLines(c("F1 S1 0 0 0 -9 A T C G",
               "F2 S2 0 0 0 -9 A A C C",
               ""), ped)
  writeLines(c("1 M1 0 1000000", "1 M2 0 2000000"), map)

  meta <- suppressMessages(suppressWarnings(read_plink_ped(ped, map)))
  expect_equal(dim(meta$geno_arrays), c(2L, 2L, 2L))
  expect_equal(dimnames(meta$geno_arrays)[[1L]], c("S1", "S2"))
})
