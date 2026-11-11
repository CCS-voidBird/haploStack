# =============================================================================
# test-090-marker-stats.R — compute_marker_stats() route numbers (design D3.7;
# §3.7 row 16 a / D22, cc/F19, cc/F21)
#
#   R4  exact two-way tie: Major != Minor; Minor = first in allele order, the
#       allele the dosage counts; MAF unchanged and equal to .marker_maf()
#   R6  Het is NA for ploidy != 2; a "numeric" (dosage) container is read as one
#       layer so Missing_rate is the true one, not 0.5
#   F19 one encoding pass: .marker_allele_counts(enc =) equals a fresh call;
#       the locus table equals the views' marker_meta on the same fields
# =============================================================================

.ms_tie_meta <- function() {
  s <- paste0("S", 1:6); m <- paste0("Mkr", 1:5)
  arr <- array(NA_character_, dim = c(6L, 5L, 2L), dimnames = list(s, m, c("H1", "H2")))
  arr[, 1, ] <- c("A","A","T","T","A","T",  "A","A","T","A","A","T")   # biallelic: 7 A / 5 T
  arr[, 2, ] <- "A"                                                     # monomorphic
  arr[, 3, 1] <- c("A","T",NA,"A","T","A"); arr[, 3, 2] <- c(NA,"T","T","A",NA,"T")
  arr[, 4, ] <- NA_character_                                           # all missing
  arr[, 5, ] <- c("A","A","A","T","T","T",  "A","A","A","T","T","T")   # 6 A / 6 T: exact tie
  map <- data.frame(Chrom = "1", Marker = m, Pos = seq_along(m) * 100, stringsAsFactors = FALSE)
  suppressMessages(new_geno_meta(arr, map_info = map,
                                 metadata = list(source_format = "user_supplied",
                                                 encoding = "nucleotide")))
}

test_that("R4: an exact two-way tie gives Major != Minor, Minor first in allele order", {
  m  <- .ms_tie_meta()
  st <- compute_marker_stats(m)
  tie <- st[Marker == "Mkr5"]
  expect_identical(tie$Minor, "A")                       # first in allele order
  expect_identical(tie$Major, "T")                       # the other allele, not "A" again
  expect_false(identical(tie$Major, tie$Minor))
  expect_equal(tie$MAF, 0.5)
  # the non-tied markers are as before
  expect_identical(st[Marker == "Mkr1"]$Major, "A"); expect_identical(st[Marker == "Mkr1"]$Minor, "T")
  expect_identical(st[Marker == "Mkr2"]$Major, "A"); expect_true(is.na(st[Marker == "Mkr2"]$Minor))
  expect_true(is.na(st[Marker == "Mkr4"]$Major)); expect_true(is.na(st[Marker == "Mkr4"]$Minor))
  # MAF still equals the shared helper, on every marker
  expect_equal(st$MAF, unname(haploStack:::.marker_maf(m)[st$Marker]))
  # the Minor on the tie is the allele the dosage counts
  mn <- haploStack:::.detect_minor_alleles(m$geno_arrays[, "Mkr5", , drop = FALSE],
                                           encode = "nucleotide")
  expect_identical(unname(mn["Mkr5"]), tie$Minor)
})

test_that("R4: .hs_major_minor() on hand-built counts", {
  cnt <- matrix(c(6L, 6L, 0L,   3L, 9L, 0L,   0L, 0L, 0L,   4L, 0L, 0L,   2L, 2L, 2L),
                ncol = 3L, byrow = TRUE, dimnames = list(paste0("m", 1:5), c("A", "C", "T")))
  mm <- haploStack:::.hs_major_minor(cnt)
  expect_identical(mm$Minor, c("A", "A", NA, NA, "A"))
  expect_identical(mm$Major, c("C", "C", NA, "A", "C"))   # tie: the other allele; three-way: first of the rest
  e <- haploStack:::.hs_major_minor(matrix(integer(0), 0L, 0L))
  expect_identical(e, list(Major = character(0), Minor = character(0)))
})

test_that("R6: Het is NA at ploidy 4; a dosage container is read as one layer", {
  ds <- sim_dataset(n_samples = 8L, n_markers = 12L, block_size = 4L, seed = 9L)
  m2 <- ds$geno_meta
  st2 <- compute_marker_stats(m2)
  expect_true(all(!is.na(st2$Het)))
  # tetraploid: duplicate the two layers
  arr4 <- array(m2$geno_arrays[, , c(1L, 2L, 1L, 2L)], dim = c(dim(m2$geno_arrays)[1:2], 4L),
                dimnames = c(dimnames(m2$geno_arrays)[1:2], list(NULL)))
  m4 <- suppressMessages(new_geno_meta(arr4, map_info = m2$map_info,
                                       metadata = list(source_format = "user_supplied",
                                                       encoding = "nucleotide")))
  expect_identical(m4$metadata$ploidy, 4L)
  st4 <- compute_marker_stats(m4)
  expect_true(all(is.na(st4$Het)))
  expect_equal(st4$MAF, st2$MAF)                          # same allele proportions
  expect_identical(st4$Major, st2$Major)
  # dosage ("numeric") container: layer 2 is the phantom NA layer
  dos <- matrix(sample(0:2, 8L * 6L, replace = TRUE), 8L, 6L,
                dimnames = list(paste0("S", 1:8), paste0("D", 1:6)))
  dos[1, 1] <- NA
  map <- data.frame(Chrom = "1", Marker = colnames(dos), Pos = 1:6 * 10)
  md <- suppressMessages(suppressWarnings(build_geno_meta(dos, map, encoding = "numeric")))
  expect_true(all(is.na(md$geno_arrays[, , 2L])))
  sd <- compute_marker_stats(md)
  expect_equal(sd[Marker == "D1"]$Missing_rate, 1 / 8)   # not 0.5 + 1/16
  expect_equal(sd[Marker == "D2"]$Missing_rate, 0)
  expect_equal(sd[Marker == "D2"]$Call_rate, 1)
  expect_true(all(is.na(sd$Het)))
})

test_that("F19: one encoding pass and the locus table give the same statistics", {
  ds <- sim_dataset(n_samples = 10L, n_markers = 16L, block_size = 4L, seed = 2L)
  m  <- ds$geno_meta
  miss <- haploStack:::.missing_tokens("nucleotide", haploStack:::.meta_extra_missing(m))
  enc  <- haploStack:::.encode_marker_codes(m$geno_arrays, na_vals = miss)
  c1 <- haploStack:::.marker_allele_counts(m$geno_arrays, na_vals = miss)
  c2 <- haploStack:::.marker_allele_counts(m$geno_arrays, na_vals = miss, enc = enc)
  expect_identical(c1, c2)
  loc <- haploStack:::.marker_locus_table(m)
  mm  <- geno_meta_views(m)$marker_meta
  expect_setequal(loc$Marker, mm$marker)
  expect_identical(loc[match(mm$marker, Marker), Chrom], as.character(mm$chromosome))
  expect_equal(loc[match(mm$marker, Marker), Pos], as.numeric(mm$position))
  expect_identical(loc[match(mm$marker, Marker), Block], as.character(mm$block_id))
  st <- compute_marker_stats(m)
  expect_identical(unname(as.character(st$Chrom)), as.character(mm$chromosome))   # same natural order
  expect_identical(st$Marker, mm$marker)
})
