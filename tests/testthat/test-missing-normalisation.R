# =============================================================================
# test-missing-normalisation.R
#
# MAN (0.6.0) — one representation of a missing allele inside a GenoMeta.
# Spec: dev/PRD_missing_allele_normalisation.md; plan:
# dev/UPGRADE_PLAN_missing_allele_normalisation.md.
#
# Stage 1 covers the primitive and the token map (FR1, FR2) with the two gates
# that decide whether the rest of the work is safe:
#   V2  a dosage source keeps its "0"        (the regression that matters most)
#   V8  a selectiontools source loses it     (Q21.a)
# plus the first tests read_plink_ped() has ever had.
# =============================================================================

library(testthat)
library(haploStack)

# -- helpers ------------------------------------------------------------------

# Write a minimal PED/MAP pair. `miss_frac` of marker 1's genotypes are written
# as the PLINK missing code "0 0".
write_ped_fixture <- function(n_samples = 60L, n_markers = 8L,
                              miss_frac = 0.4, seed = 11L) {
  set.seed(seed)
  markers <- paste0("M", seq_len(n_markers))
  samples <- paste0("S", seq_len(n_samples))

  a1 <- matrix(sample(c("A", "T"), n_samples * n_markers, replace = TRUE),
               nrow = n_samples)
  a2 <- matrix(sample(c("A", "T"), n_samples * n_markers, replace = TRUE),
               nrow = n_samples)
  n_miss <- as.integer(round(miss_frac * n_samples))
  if (n_miss > 0L) {
    miss_rows <- seq_len(n_miss)
    a1[miss_rows, 1L] <- "0"
    a2[miss_rows, 1L] <- "0"
  }

  # Interleave the two strands: marker 1 allele 1, marker 1 allele 2, ...
  alleles <- matrix(NA_character_, nrow = n_samples, ncol = 2L * n_markers)
  alleles[, seq(1L, 2L * n_markers, by = 2L)] <- a1
  alleles[, seq(2L, 2L * n_markers, by = 2L)] <- a2

  ped <- cbind(
    data.frame(FID = paste0("F", seq_len(n_samples)), IID = samples,
               PAT = "0", MAT = "0", SEX = "0", PHENO = "-9",
               stringsAsFactors = FALSE),
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

  list(ped = ped_file, map = map_file, markers = markers, samples = samples,
       n_miss = n_miss)
}

# A (samples x markers x ploidy) character array from a per-strand matrix pair.
as_allele_array <- function(m1, m2) {
  arr <- array(NA_character_,
               dim = c(nrow(m1), ncol(m1), 2L),
               dimnames = list(rownames(m1), colnames(m1), c("H1", "H2")))
  arr[, , 1L] <- m1
  arr[, , 2L] <- m2
  arr
}

# -- FR2: the token map -------------------------------------------------------

test_that("FR2: '0' is missing for nucleotide and selectiontools, never for numeric", {
  expect_true("0" %in% haploStack:::.missing_tokens("nucleotide"))
  expect_true("0" %in% haploStack:::.missing_tokens("selectiontools"))
  expect_false("0" %in% haploStack:::.missing_tokens("numeric"))
})

test_that("FR2: the encoding-independent token set is unchanged", {
  expect_setequal(haploStack:::.missing_tokens("numeric"),
                  c("", "NA", "N", "n", ".", "-"))
})

# -- FR1: the primitive -------------------------------------------------------

test_that("FR1: nucleotide tokens map to NA and the count is reported", {
  arr <- as_allele_array(
    matrix(c("A", "0", "T", "-"), nrow = 2,
           dimnames = list(c("S1", "S2"), c("M1", "M2"))),
    matrix(c("A", "N", "T", "T"), nrow = 2,
           dimnames = list(c("S1", "S2"), c("M1", "M2"))))

  res <- haploStack:::normalise_missing_alleles(arr, "nucleotide")

  expect_equal(res$n_changed, 3L)
  expect_true(is.na(res$array["S2", "M1", "H1"]))
  expect_true(is.na(res$array["S2", "M1", "H2"]))
  expect_true(is.na(res$array["S2", "M2", "H1"]))
  expect_equal(unname(res$array["S1", "M1", "H1"]), "A")
  expect_equal(unname(res$array["S2", "M2", "H2"]), "T")
})

test_that("FR1: dimensions and dimnames survive normalisation", {
  arr <- as_allele_array(
    matrix(c("A", "0"), nrow = 1, dimnames = list("S1", c("M1", "M2"))),
    matrix(c("A", "0"), nrow = 1, dimnames = list("S1", c("M1", "M2"))))
  res <- haploStack:::normalise_missing_alleles(arr, "nucleotide")
  expect_equal(dim(res$array), dim(arr))
  expect_equal(dimnames(res$array), dimnames(arr))
})

test_that("FR1: a clean array is returned untouched with n_changed = 0", {
  arr <- as_allele_array(
    matrix(c("A", "T"), nrow = 1, dimnames = list("S1", c("M1", "M2"))),
    matrix(c("T", "A"), nrow = 1, dimnames = list("S1", c("M1", "M2"))))
  res <- haploStack:::normalise_missing_alleles(arr, "nucleotide")
  expect_equal(res$n_changed, 0L)
  expect_identical(res$array, arr)
})

test_that("FR1: existing NA cells are left alone and not counted", {
  arr <- as_allele_array(
    matrix(c("A", NA), nrow = 1, dimnames = list("S1", c("M1", "M2"))),
    matrix(c("A", NA), nrow = 1, dimnames = list("S1", c("M1", "M2"))))
  res <- haploStack:::normalise_missing_alleles(arr, "nucleotide")
  expect_equal(res$n_changed, 0L)
  expect_true(all(is.na(res$array[, "M2", ])))
})

test_that("FR1: whitespace-only cells are missing under every source format", {
  arr <- as_allele_array(
    matrix(c("A", "  "), nrow = 1, dimnames = list("S1", c("M1", "M2"))),
    matrix(c("A", ""),   nrow = 1, dimnames = list("S1", c("M1", "M2"))))
  res <- haploStack:::normalise_missing_alleles(arr, "numeric")
  expect_equal(res$n_changed, 2L)
  expect_true(all(is.na(res$array[, "M2", ])))
})

test_that("FR1: extra_na adds a caller-declared sentinel", {
  arr <- as_allele_array(
    matrix(c("A", "?"), nrow = 1, dimnames = list("S1", c("M1", "M2"))),
    matrix(c("A", "?"), nrow = 1, dimnames = list("S1", c("M1", "M2"))))
  expect_equal(
    haploStack:::normalise_missing_alleles(arr, "nucleotide")$n_changed, 0L)
  res <- haploStack:::normalise_missing_alleles(arr, "nucleotide",
                                                extra_na = "?")
  expect_equal(res$n_changed, 2L)
})

test_that("FR1: a non-character array is returned unchanged", {
  arr <- array(c(0L, 1L, 2L, 0L), dim = c(2L, 1L, 2L))
  res <- haploStack:::normalise_missing_alleles(arr, "numeric")
  expect_identical(res$array, arr)
  expect_equal(res$n_changed, 0L)
})

# -- V2: the regression that matters most, a dosage source keeps its 0 --------

test_that("V2: '0' survives normalisation for a dosage ('numeric') source", {
  arr <- as_allele_array(
    matrix(c("0", "1", "2", "0"), nrow = 2,
           dimnames = list(c("S1", "S2"), c("M1", "M2"))),
    matrix(c("0", "1", "2", "0"), nrow = 2,
           dimnames = list(c("S1", "S2"), c("M1", "M2"))))

  res <- haploStack:::normalise_missing_alleles(arr, "numeric")

  expect_equal(res$n_changed, 0L)
  expect_identical(res$array, arr)
  expect_false(any(is.na(res$array)))
})

test_that("V2: extra_na cannot smuggle '0' into a dosage source", {
  arr <- as_allele_array(
    matrix(c("0", "2"), nrow = 1, dimnames = list("S1", c("M1", "M2"))),
    matrix(c("0", "2"), nrow = 1, dimnames = list("S1", c("M1", "M2"))))
  expect_warning(
    haploStack:::normalise_missing_alleles(arr, "numeric", extra_na = "0"),
    "dosage")
  res <- suppressWarnings(
    haploStack:::normalise_missing_alleles(arr, "numeric", extra_na = "0"))
  expect_equal(res$n_changed, 0L)
  expect_equal(unname(res$array["S1", "M1", "H1"]), "0")
})

# -- V8: a selectiontools source loses it (Q21.a) -----------------------------

test_that("V8: '0' is normalised away for a selectiontools source", {
  arr <- as_allele_array(
    matrix(c("1", "0", "4", "2"), nrow = 2,
           dimnames = list(c("S1", "S2"), c("M1", "M2"))),
    matrix(c("1", "0", "4", "3"), nrow = 2,
           dimnames = list(c("S1", "S2"), c("M1", "M2"))))

  res <- haploStack:::normalise_missing_alleles(arr, "selectiontools")

  expect_equal(res$n_changed, 2L)
  expect_true(all(is.na(res$array["S2", "M1", ])))
  expect_setequal(setdiff(unique(as.vector(res$array)), NA_character_),
                  c("1", "4", "2", "3"))
})

test_that("V8: a 1/2/3/4-coded panel read as ST has its '0' removed by the primitive", {
  geno <- matrix(c("1/1", "0/0", "4/4", "1/4"), nrow = 2,
                 dimnames = list(c("M1", "M2"), c("S1", "S2")))
  split_res <- split_haplotypes(geno, orientation = 2L)
  res <- haploStack:::normalise_missing_alleles(split_res$array,
                                                "selectiontools")
  expect_equal(res$n_changed, 2L)
  expect_false(any(as.vector(res$array) %in% "0"))
})

# -- read_plink_ped(): its first tests (stage 1, before stage 3 touches it) ---

test_that("read_plink_ped parses dimensions, ploidy and identifiers", {
  fx   <- write_ped_fixture(n_samples = 20L, n_markers = 6L, miss_frac = 0)
  meta <- read_plink_ped(fx$ped, fx$map)

  expect_true(is_geno_meta(meta))
  expect_equal(dim(meta$geno_arrays), c(20L, 6L, 2L))
  expect_equal(meta$metadata$ploidy, 2L)
  expect_equal(meta$metadata$n_samples, 20L)
  expect_equal(meta$metadata$n_markers, 6L)
  expect_equal(meta$metadata$source_format, "plink_ped")
  expect_equal(dimnames(meta$geno_arrays)[[1L]], paste0("S", 1:20))
  expect_equal(dimnames(meta$geno_arrays)[[2L]], paste0("M", 1:6))
})

test_that("read_plink_ped carries the PED's FAM columns into sample_info", {
  fx   <- write_ped_fixture(n_samples = 8L, n_markers = 4L, miss_frac = 0)
  meta <- read_plink_ped(fx$ped, fx$map)
  expect_equal(names(meta$sample_info)[1:5],
               c("Family", "Sample", "Parent1", "Parent2", "Sex"))
  expect_equal(meta$sample_info$Family, paste0("F", 1:8))
  expect_equal(meta$sample_info$Sample, paste0("S", 1:8))
})

test_that("read_plink_ped reads the MAP file's fourth column as position", {
  fx   <- write_ped_fixture(n_samples = 6L, n_markers = 5L, miss_frac = 0)
  meta <- read_plink_ped(fx$ped, fx$map)
  expect_equal(meta$map_info$Marker, paste0("M", 1:5))
  expect_equal(meta$map_info$Pos, seq_len(5L) * 1000000)
  expect_equal(unique(meta$map_info$Chrom), "1")
})

test_that("read_plink_ped keeps the missing code out of char_dict", {
  fx   <- write_ped_fixture(n_samples = 60L, n_markers = 8L, miss_frac = 0.4)
  meta <- read_plink_ped(fx$ped, fx$map)
  expect_false("0" %in% names(meta$metadata$char_dict))
  expect_setequal(names(meta$metadata$char_dict), c("A", "T"))
})

test_that("read_plink_ped errors on a MAP file with too few columns", {
  fx  <- write_ped_fixture(n_samples = 4L, n_markers = 3L, miss_frac = 0)
  bad <- tempfile(fileext = ".map")
  utils::write.table(data.frame(Chrom = "1", Marker = c("M1", "M2", "M3")),
                     bad, quote = FALSE, sep = " ",
                     row.names = FALSE, col.names = FALSE)
  expect_error(read_plink_ped(fx$ped, bad), "3 or 4 columns")
})

# -- stage 2: detection before correction (FR4, FR5, FR6) --------------------

# A container whose marker M1 is written with a source token instead of NA —
# i.e. every pre-0.6.0 container on disk.
#
# Since stage 3 the constructor IS the normalisation portal, so this cannot be
# built through the API any more: the token would be normalised away on the way
# in. It is therefore assembled clean and the token written into geno_arrays
# afterwards, which is the only way such a container can still exist — restored
# from an old .rds, or assembled by hand. That population is exactly what FR4,
# FR5 and FR6 are for.
token_container <- function(encoding = "nucleotide", token = "0") {
  clean <- if (identical(encoding, "numeric")) c("1", "2") else c("A", "T")
  arr <- as_allele_array(
    matrix(c(clean[1L], clean[1L], clean[2L], clean[1L]), nrow = 2,
           dimnames = list(c("S1", "S2"), c("M1", "M2"))),
    matrix(c(clean[1L], clean[1L], clean[2L], clean[2L]), nrow = 2,
           dimnames = list(c("S1", "S2"), c("M1", "M2"))))
  map <- data.frame(Chrom = "1", Marker = c("M1", "M2"),
                    Pos = c(1e6, 2e6), stringsAsFactors = FALSE)
  m <- suppressMessages(new_geno_meta(
    arr, map,
    metadata = list(ploidy = 2L,
                    char_dict = stats::setNames(as.list(c(1L, 4L)), clean),
                    n_samples = 2L, n_markers = 2L,
                    source_format = "hand_built", encoding = encoding)))

  m$geno_arrays["S2", "M1", ] <- token
  m$metadata$missing_tokens   <- NULL   # a legacy container carries no stamp
  m
}

test_that("V5: is_geno_meta warns on a token container and repair clears it", {
  m <- token_container()

  expect_warning(is_geno_meta(m), "missing token")
  expect_warning(is_geno_meta(m), "repair_geno_meta")
  expect_warning(is_geno_meta(m), "1 marker")
  expect_true(suppressWarnings(is_geno_meta(m)))

  expect_message(m2 <- repair_geno_meta(m), "normalised 2")
  expect_true(all(is.na(m2$geno_arrays["S2", "M1", ])))
  expect_silent(is_geno_meta(m2))
  expect_true(is_geno_meta(m2))
})

test_that("V5: a dosage container with real 0 genotypes never warns", {
  m <- token_container(encoding = "numeric", token = "0")
  expect_silent(is_geno_meta(m))
  expect_true(is_geno_meta(m))
  expect_message(m2 <- repair_geno_meta(m), "already normalised")
  expect_identical(m2$geno_arrays, m$geno_arrays)
})

test_that("FR5: repair_geno_meta reports and is idempotent", {
  m  <- token_container()
  m2 <- suppressMessages(repair_geno_meta(m))
  expect_message(m3 <- repair_geno_meta(m2), "already normalised")
  expect_identical(m2$geno_arrays, m3$geno_arrays)
  expect_silent(repair_geno_meta(m2, verbose = FALSE))
})

test_that("FR5: repair_geno_meta follows an explicit source_format over the record", {
  # A dosage panel mislabelled 'nucleotide' would lose its 0 genotypes to the
  # record; the caller's declaration governs.
  m <- token_container(encoding = "nucleotide", token = "0")
  m2 <- suppressMessages(repair_geno_meta(m, source_format = "numeric"))
  expect_equal(unname(m2$geno_arrays["S2", "M1", "H1"]), "0")
})

test_that("V9: options(haploStack.strict_missing) turns the warning into an error", {
  m <- token_container()
  withr::with_options(list(haploStack.strict_missing = TRUE), {
    expect_error(is_geno_meta(m), "missing token")
  })
  # Default stays a warning — a hard stop by default would break every
  # container already on disk.
  expect_warning(is_geno_meta(m), "missing token")
})

test_that("V10: load_geno_meta mentions the repair and leaves the file untouched", {
  m    <- token_container()
  path <- tempfile(fileext = ".rds")
  saveRDS(m, path)
  before <- tools::md5sum(path)

  expect_message(load_geno_meta(path), "repair_geno_meta")
  expect_message(load_geno_meta(path), "NOT been modified")

  loaded <- suppressMessages(load_geno_meta(path))
  expect_equal(unname(loaded$geno_arrays["S2", "M1", "H1"]), "0")
  expect_equal(unname(tools::md5sum(path)), unname(before))
})

test_that("V10: loading a normalised container is silent", {
  m    <- suppressMessages(repair_geno_meta(token_container()))
  path <- tempfile(fileext = ".rds")
  saveRDS(m, path)
  expect_silent(load_geno_meta(path))
})

# -- the nucleotide default: 0 and - are both missing, per allele -------------

test_that("nucleotide treats both '0' and '-' as missing alleles", {
  toks <- haploStack:::.missing_tokens("nucleotide")
  expect_true(all(c("0", "-") %in% toks))

  arr <- as_allele_array(
    matrix(c("A", "0", "T", "-"), nrow = 2,
           dimnames = list(c("S1", "S2"), c("M1", "M2"))),
    matrix(c("A", "0", "T", "-"), nrow = 2,
           dimnames = list(c("S1", "S2"), c("M1", "M2"))))
  res <- haploStack:::normalise_missing_alleles(arr, "nucleotide")
  expect_equal(res$n_changed, 4L)
  expect_true(all(is.na(res$array["S2", , ])))
})

test_that("a '00' or '--' genotype cell becomes a missing genotype", {
  # Both the concatenated and the slash-separated writing of the same cell.
  for (g in list(
    matrix(c("AA", "00", "AT", "--"), nrow = 2,
           dimnames = list(c("M1", "M2"), c("S1", "S2"))),
    matrix(c("A/A", "0/0", "A/T", "-/-"), nrow = 2,
           dimnames = list(c("M1", "M2"), c("S1", "S2"))))) {

    sc  <- if (any(grepl("/", g))) "/" else ""
    spl <- split_haplotypes(g, split_char = sc, orientation = 2L)
    res <- haploStack:::normalise_missing_alleles(spl$array, "nucleotide")

    expect_equal(res$n_changed, 4L)
    expect_true(all(is.na(res$array["S1", "M2", ])))   # "00"
    expect_true(all(is.na(res$array["S2", "M2", ])))   # "--"
    expect_false(any(is.na(res$array[, "M1", ])))      # real genotypes intact
  }
})

test_that("a half-called genotype keeps its observed allele", {
  g   <- matrix(c("A0", "A-"), nrow = 1,
                dimnames = list("M1", c("S1", "S2")))
  spl <- split_haplotypes(g, split_char = "", orientation = 2L)
  res <- haploStack:::normalise_missing_alleles(spl$array, "nucleotide")

  expect_equal(res$n_changed, 2L)
  expect_equal(unname(res$array["S1", "M1", 1L]), "A")
  expect_true(is.na(res$array["S1", "M1", 2L]))
  expect_equal(unname(res$array["S2", "M1", 1L]), "A")
  expect_true(is.na(res$array["S2", "M1", 2L]))
})

# -- use_missing_chars(): the session registry -------------------------------

test_that("use_missing_chars extends one mode and leaves the others alone", {
  withr::local_options(list(haploStack.missing_chars = NULL))

  expect_message(use_missing_chars("#", mode = "nucleotide"),
                 "now read as missing")

  expect_true("#" %in% missing_chars("nucleotide"))
  expect_false("#" %in% missing_chars("selectiontools"))
  expect_false("#" %in% missing_chars("numeric"))
  expect_equal(missing_chars("nucleotide", chars_only = TRUE), "#")
  expect_equal(missing_chars("numeric", chars_only = TRUE), character(0))
})

test_that("a registered character is honoured by the primitive", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  arr <- as_allele_array(
    matrix(c("A", "#"), nrow = 1, dimnames = list("S1", c("M1", "M2"))),
    matrix(c("A", "#"), nrow = 1, dimnames = list("S1", c("M1", "M2"))))

  expect_equal(
    haploStack:::normalise_missing_alleles(arr, "nucleotide")$n_changed, 0L)

  suppressMessages(use_missing_chars("#", mode = "nucleotide"))
  res <- haploStack:::normalise_missing_alleles(arr, "nucleotide")
  expect_equal(res$n_changed, 2L)
  expect_true(all(is.na(res$array[, "M2", ])))
})

test_that("a registered character reaches a container through repair, not through the check", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  m <- token_container(encoding = "nucleotide", token = "#")

  expect_silent(is_geno_meta(m))                 # not missing until declared

  # 0.9.0 (form row 2): the registry is read where a container is built or
  # repaired, never by a consumer. is_geno_meta() goes by the built-in set and
  # the container's own stamp, so a later registration does not change what
  # it says; repair_geno_meta() does read the registry.
  suppressMessages(use_missing_chars("#", mode = "nucleotide"))
  expect_silent(is_geno_meta(m))
  expect_message(m2 <- repair_geno_meta(m), "normalised 2")
  expect_true(all(is.na(m2$geno_arrays["S2", "M1", ])))
  expect_true("#" %in% m2$metadata$missing_tokens)
  # and once stamped, a full scan sees "#" in a session with nothing registered
  reset_missing_chars()
  m3 <- m2; m3$geno_arrays["S1", "M2", ] <- "#"
  expect_warning(is_geno_meta(m3, strict = TRUE), "missing token")
})

test_that("multiple characters can be registered, and registration is unique", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  suppressMessages(use_missing_chars(c("#", "?"), mode = "nucleotide"))
  suppressMessages(use_missing_chars("#", mode = "nucleotide"))
  expect_setequal(missing_chars("nucleotide", chars_only = TRUE), c("#", "?"))
})

test_that("reset_missing_chars clears one mode or all of them", {
  withr::local_options(list(haploStack.missing_chars = NULL,
                            haploStack.st_formats = TRUE))      # the ST mode is gated (0.9.0, C-6)
  suppressMessages(use_missing_chars("#", mode = "nucleotide"))
  suppressWarnings(suppressMessages(use_missing_chars("#", mode = "selectiontools")))

  reset_missing_chars("nucleotide")
  expect_equal(missing_chars("nucleotide", chars_only = TRUE), character(0))
  expect_equal(missing_chars("selectiontools", chars_only = TRUE), "#")

  reset_missing_chars()
  expect_equal(missing_chars("selectiontools", chars_only = TRUE), character(0))
  # The built-in vocabulary is not clearable.
  expect_setequal(missing_chars("nucleotide"),
                  c("", "NA", "N", "n", ".", "-", "0"))
})

test_that("use_missing_chars refuses anything that is not one character", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  expect_error(use_missing_chars("##", mode = "nucleotide"),
               "exactly one character")
  expect_error(use_missing_chars("-9", mode = "nucleotide"),
               "exactly one character")
  expect_error(use_missing_chars(NA_character_, mode = "nucleotide"), "NA")
  expect_error(use_missing_chars(" ", mode = "nucleotide"), "already")
  expect_error(use_missing_chars(9, mode = "nucleotide"), "character vector")
})

test_that("use_missing_chars refuses a real allele of the mode", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  expect_error(use_missing_chars("A", mode = "nucleotide"), "real allele")
  expect_error(use_missing_chars("t", mode = "nucleotide"), "real allele")
  expect_error(use_missing_chars("1", mode = "selectiontools"), "real allele")
  expect_error(use_missing_chars("0", mode = "numeric"), "real allele")
  expect_error(use_missing_chars("9", mode = "numeric"), "extra_na")
  # "0" is not an allele of nucleotide — it is already missing there.
  expect_message(use_missing_chars("0", mode = "nucleotide"),
                 "now read as missing")
})

test_that("mode is required and validated", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  expect_error(use_missing_chars("#"), "'mode' is required")
  expect_error(use_missing_chars("#", mode = "dosage"), "unknown mode")
  # "allele" is accepted as the dosage alias of "numeric".
  suppressMessages(use_missing_chars("#", mode = "allele"))
  expect_equal(missing_chars("numeric", chars_only = TRUE), "#")
})

test_that("a malformed option cannot change what is missing", {
  # missing_chars() reads the registry the way the build portal does; the
  # built-in set itself (.missing_tokens()) reads no option at all (0.9.0).
  withr::local_options(list(haploStack.missing_chars = "nonsense"))
  expect_setequal(missing_chars("nucleotide"),
                  c("", "NA", "N", "n", ".", "-", "0"))
  withr::local_options(
    list(haploStack.missing_chars = list(nucleotide = c("##", NA))))
  expect_setequal(missing_chars("nucleotide"),
                  c("", "NA", "N", "n", ".", "-", "0"))
  expect_setequal(haploStack:::.missing_tokens("nucleotide"),
                  c("", "NA", "N", "n", ".", "-", "0"))
})

# -- metadata$missing_tokens: the container stays self-describing -------------

test_that("repair stamps the vocabulary it cleaned against", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  suppressMessages(use_missing_chars("#", mode = "nucleotide"))
  m2 <- suppressMessages(repair_geno_meta(token_container(token = "#")))

  expect_true("#" %in% m2$metadata$missing_tokens)
  expect_true(all(c("0", "-") %in% m2$metadata$missing_tokens))
})

test_that("the stamp keeps detection working with nothing registered", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  m <- token_container(token = "#")
  m$metadata$missing_tokens <- c("", "NA", "N", "n", ".", "-", "0", "#")

  # "#" is not registered here; the container's own record supplies it. The
  # stamp is trusted by default (0.9.0), so a hand-edited array needs the
  # full scan to be seen.
  expect_silent(is_geno_meta(m))
  expect_warning(is_geno_meta(m, strict = TRUE), "missing token")
  expect_message(m2 <- repair_geno_meta(m), "normalised 2")
  expect_true(all(is.na(m2$geno_arrays["S2", "M1", ])))
})

test_that("a dosage container's 0 survives both the stamp and the registry", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  m <- token_container(encoding = "numeric", token = "0")
  m$metadata$missing_tokens <- c("", "NA", "N", "n", ".", "-", "0")

  expect_silent(is_geno_meta(m))
  m2 <- suppressMessages(repair_geno_meta(m))
  expect_equal(unname(m2$geno_arrays["S2", "M1", "H1"]), "0")
  expect_false("0" %in% m2$metadata$missing_tokens)
})

# =============================================================================
# Stage 3 — normalise at the portal (FR3), and stage 4 (FR7)
# =============================================================================

# -- the portal is the constructor, and every route reaches it ---------------

test_that("V1: read_plink_ped leaves no missing token in geno_arrays", {
  fx   <- write_ped_fixture(n_samples = 60L, n_markers = 8L, miss_frac = 0.4)
  meta <- suppressMessages(read_plink_ped(fx$ped, fx$map))

  expect_false(any(meta$geno_arrays %in% "0"))
  expect_equal(sum(is.na(meta$geno_arrays["S1", "M1", ])), 2L)
  expect_setequal(setdiff(unique(as.vector(meta$geno_arrays)), NA_character_),
                  c("A", "T"))
  # 24 samples x 1 marker x 2 copies were written "0 0" in the PED.
  expect_equal(sum(is.na(meta$geno_arrays)), 2L * fx$n_miss)
})

test_that("V1: a normalised PED fits marker effects without the >2 alleles warning", {
  skip_if_not_installed("rrBLUP")
  fx   <- write_ped_fixture(n_samples = 60L, n_markers = 8L, miss_frac = 0.4)
  meta <- suppressMessages(read_plink_ped(fx$ped, fx$map))

  set.seed(3L)
  ph <- data.frame(Sample = dimnames(meta$geno_arrays)[[1L]],
                   Y = rnorm(dim(meta$geno_arrays)[1L]),
                   stringsAsFactors = FALSE)

  expect_no_warning(
    fit <- suppressMessages(
      fit_genomic_prediction(meta, phenotype = ph, model = "rrblup",
                             trait = "Y")))

  eff <- fit$marker_effects
  if (is.list(eff) && !is.data.frame(eff)) eff <- eff[[1L]]
  expect_false(any(as.character(eff$Allele) == "0"))
})

test_that("read_plink_ped declares its format, so a 1/2-coded PED is not dosage", {
  # PLINK's numeric allele coding. After splitting, the alleles are single
  # digits — detection would call that "numeric", and the dosage rule would
  # then keep "0" as a real genotype.
  set.seed(5L)
  n <- 40L; m <- 4L
  a1 <- matrix(sample(c("1", "2"), n * m, replace = TRUE), nrow = n)
  a2 <- matrix(sample(c("1", "2"), n * m, replace = TRUE), nrow = n)
  a1[1:10, 1L] <- "0"; a2[1:10, 1L] <- "0"
  al <- matrix(NA_character_, nrow = n, ncol = 2L * m)
  al[, seq(1L, 2L * m, by = 2L)] <- a1
  al[, seq(2L, 2L * m, by = 2L)] <- a2
  ped <- cbind(data.frame(FID = paste0("F", 1:n), IID = paste0("S", 1:n),
                          PAT = "0", MAT = "0", SEX = "0", PHENO = "-9",
                          stringsAsFactors = FALSE),
               as.data.frame(al, stringsAsFactors = FALSE))
  pf <- tempfile(fileext = ".ped"); mf <- tempfile(fileext = ".map")
  utils::write.table(ped, pf, quote = FALSE, sep = " ",
                     row.names = FALSE, col.names = FALSE)
  utils::write.table(
    data.frame(C = "1", M = paste0("M", 1:m), cM = 0, P = (1:m) * 1000000L),
    mf, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

  meta <- suppressMessages(read_plink_ped(pf, mf))
  expect_equal(meta$metadata$encoding, "nucleotide")   # labels, since 0.9.0 (C-5)
  expect_false(any(meta$geno_arrays %in% "0"))
  expect_equal(sum(is.na(meta$geno_arrays)), 20L)
})

test_that("build_geno_meta normalises, and records the vocabulary it used", {
  geno <- matrix(c("A/A", "0/0", "A/T", "T/T"), nrow = 2,
                 dimnames = list(c("S1", "S2"), c("M1", "M2")))
  map  <- data.frame(Chrom = "1", Marker = c("M1", "M2"),
                     Pos = c(1e6, 2e6), stringsAsFactors = FALSE)

  expect_message(
    meta <- build_geno_meta(geno, map, split_char = "/",
                            encoding = "nucleotide"),
    "normalised 2 missing allele cell")

  expect_true(all(is.na(meta$geno_arrays["S2", "M1", ])))
  expect_true(all(c("0", "-") %in% meta$metadata$missing_tokens))
  expect_silent(is_geno_meta(meta))
})

test_that("the dosage rule holds at the portal: a 0/1/2 container keeps its 0", {
  d <- matrix(c("0", "1", "2", "0"), nrow = 2,
              dimnames = list(c("S1", "S2"), c("M1", "M2")))
  map <- data.frame(Chrom = "1", Marker = c("M1", "M2"),
                    Pos = c(1e6, 2e6), stringsAsFactors = FALSE)
  meta <- suppressMessages(build_geno_meta(d, map, split_char = "", ploidy = 1L,
                                           encoding = "numeric"))
  expect_true("0" %in% meta$geno_arrays)
  expect_false("0" %in% meta$metadata$missing_tokens)
})

test_that("extra_na is the per-read escape hatch for an unknown sentinel", {
  geno <- matrix(c("A/A", "9/9", "A/T", "T/T"), nrow = 2,
                 dimnames = list(c("S1", "S2"), c("M1", "M2")))
  map  <- data.frame(Chrom = "1", Marker = c("M1", "M2"),
                     Pos = c(1e6, 2e6), stringsAsFactors = FALSE)

  plain <- suppressMessages(build_geno_meta(geno, map, split_char = "/",
                                            encoding = "nucleotide"))
  expect_true("9" %in% plain$geno_arrays)

  declared <- suppressMessages(build_geno_meta(geno, map, split_char = "/",
                                               encoding = "nucleotide",
                                               extra_na = "9"))
  expect_true(all(is.na(declared$geno_arrays["S2", "M1", ])))
  expect_true("9" %in% declared$metadata$missing_tokens)
})

test_that("a derived container carries the stamp and is not re-scanned", {
  geno <- matrix(c("A/A", "0/0", "A/T", "T/T"), nrow = 2,
                 dimnames = list(c("S1", "S2"), c("M1", "M2")))
  map  <- data.frame(Chrom = "1", Marker = c("M1", "M2"),
                     Pos = c(1e6, 2e6), stringsAsFactors = FALSE)
  meta <- suppressMessages(build_geno_meta(geno, map, split_char = "/",
                                           encoding = "nucleotide"))
  # A second pass through the portal has nothing to say — no message, and the
  # stamp survives.
  sub <- suppressWarnings(filter_geno_meta(meta, keep_markers = "M1"))
  expect_true(all(c("0", "-") %in% sub$metadata$missing_tokens))
  expect_silent(is_geno_meta(sub))
})

# -- V3 / V4: what the normalisation fixes downstream ------------------------

test_that("V3: a token source and an NA source now build the same container", {
  map <- data.frame(Chrom = "1", Marker = c("M1", "M2", "M3", "M4"),
                    Pos = (1:4) * 1e6, stringsAsFactors = FALSE)
  set.seed(8L)
  cells <- matrix(sample(c("A/A", "A/T", "T/T"), 40L, replace = TRUE), nrow = 10,
                  dimnames = list(paste0("S", 1:10), map$Marker))
  miss  <- matrix(runif(40L) < 0.25, nrow = 10)

  tok <- cells; tok[miss] <- "0/0"
  nas <- cells; nas[miss] <- NA_character_

  m_tok <- suppressMessages(build_geno_meta(tok, map, split_char = "/",
                                            encoding = "nucleotide"))
  m_na  <- suppressMessages(build_geno_meta(nas, map, split_char = "/",
                                            encoding = "nucleotide"))

  expect_identical(m_tok$geno_arrays, m_na$geno_arrays)

  eff <- data.frame(Marker = rep(map$Marker, each = 2L),
                    Allele = rep(c("A", "T"), 4L),
                    Effect = c(1, -1, 0.5, -0.5, 2, -2, 0.25, -0.25),
                    stringsAsFactors = FALSE)
  g_tok <- suppressMessages(suppressWarnings(compute_profile_localgebv(
    m_tok, effect_df = eff, center = "allele_freq")))
  g_na <- suppressMessages(suppressWarnings(compute_profile_localgebv(
    m_na, effect_df = eff, center = "allele_freq")))
  expect_equal(g_tok, g_na)
})

test_that("V4: signatures use fill_na and respond to the option", {
  map <- data.frame(Chrom = "1", Marker = c("M1", "M2", "M3"),
                    Pos = (1:3) * 1e6, stringsAsFactors = FALSE)
  bmap <- data.frame(Chrom = "1", Marker = map$Marker, Block = "B1",
                     stringsAsFactors = FALSE)
  cells <- matrix(c("A/A", "A/T", "T/T",
                    "A/A", "0/0", "T/T"), nrow = 2, byrow = TRUE,
                  dimnames = list(c("S1", "S2"), map$Marker))
  meta <- suppressMessages(build_geno_meta(cells, map, blocks = bmap,
                                           split_char = "/",
                                           encoding = "nucleotide"))
  sig <- function() {
    st <- suppressMessages(suppressWarnings(count_block_subtypes(meta)))
    unique(as.character(as.data.frame(st)$Subtype))
  }
  expect_true(any(grepl("-", sig(), fixed = TRUE)))
  expect_false(any(grepl("0", sig(), fixed = TRUE)))

  withr::with_options(list(haploStack.fill_na = "?"), {
    expect_true(any(grepl("?", sig(), fixed = TRUE)))
  })
})

# -- V11: recode declares both ends, and renders missing one way -------------

test_that("V11: recode_geno_meta refuses a dosage source for an allele target", {
  d <- matrix(c("0", "1", "2", "0"), nrow = 2,
              dimnames = list(c("S1", "S2"), c("M1", "M2")))
  map <- data.frame(Chrom = "1", Marker = c("M1", "M2"),
                    Pos = c(1e6, 2e6), stringsAsFactors = FALSE)
  meta <- suppressMessages(build_geno_meta(d, map, split_char = "", ploidy = 1L,
                                           encoding = "numeric"))
  expect_error(recode_geno_meta(meta, to_type = "nucleotide"), "dosage scale")
})

test_that("V11: a declared from_type is validated against the data", {
  geno <- matrix(c("A/A", "A/T"), nrow = 1,
                 dimnames = list("S1", c("M1", "M2")))
  map  <- data.frame(Chrom = "1", Marker = c("M1", "M2"),
                     Pos = c(1e6, 2e6), stringsAsFactors = FALSE)
  meta <- suppressMessages(build_geno_meta(geno, map, split_char = "/",
                                           encoding = "nucleotide"))
  expect_error(
    recode_geno_meta(meta, to_type = "nucleotide",
                     from_type = "selectiontools"),
    "alphabetic")
})

test_that("V11: recode follows from_type over a wrong recorded encoding", {
  geno <- matrix(c("A/A", "A/T"), nrow = 1,
                 dimnames = list("S1", c("M1", "M2")))
  map  <- data.frame(Chrom = "1", Marker = c("M1", "M2"),
                     Pos = c(1e6, 2e6), stringsAsFactors = FALSE)
  meta <- suppressMessages(build_geno_meta(geno, map, split_char = "/",
                                           encoding = "nucleotide"))
  # a container still holding SelectionTools codes, with a stale, wrong record
  meta$geno_arrays[] <- chartr("AT", "14", meta$geno_arrays)
  meta$metadata$encoding <- "numeric"

  out <- suppressMessages(recode_geno_meta(meta, to_type = "nucleotide",
                                           from_type = "selectiontools"))
  expect_equal(unname(out$geno_arrays["S1", "M1", 1L]), "A")
  expect_equal(out$metadata$encoding, "nucleotide")
})

test_that("na_as renders missing back out, one way only", {
  geno <- matrix(c("A/A", "0/0"), nrow = 1,
                 dimnames = list("S1", c("M1", "M2")))
  map  <- data.frame(Chrom = "1", Marker = c("M1", "M2"),
                     Pos = c(1e6, 2e6), stringsAsFactors = FALSE)
  meta <- suppressMessages(build_geno_meta(geno, map, split_char = "/",
                                           encoding = "nucleotide"))
  expect_true(all(is.na(meta$geno_arrays["S1", "M2", ])))

  expect_message(
    out <- recode_geno_meta(meta, to_type = "nucleotide", na_as = "-"),
    "rendered 2 missing allele")
  expect_equal(unname(out$geno_arrays["S1", "M2", 1L]), "-")

  # The rendering is deliberate, so the result is no longer normalised — and
  # repairing it returns NA, not the "0" the source file wrote.
  expect_warning(is_geno_meta(out), "missing token")
  back <- suppressMessages(repair_geno_meta(out))
  expect_true(all(is.na(back$geno_arrays["S1", "M2", ])))
  expect_false(any(back$geno_arrays %in% "0"))
})

test_that("detection classifies on the calls, not the absences", {
  # A heavily-missing nucleotide panel: without the filter the "00" cells
  # outvote the real calls and the panel reads as dosage.
  cells <- c(rep("00", 70L), rep("AT", 20L), rep("AA", 10L))
  expect_equal(detect_genotype_encoding(matrix(cells, nrow = 10)), "nucleotide")

  # A dosage panel is still dosage after its 0s are set aside.
  expect_equal(
    detect_genotype_encoding(matrix(c(rep("0", 50L), rep("1", 30L),
                                      rep("2", 20L)), nrow = 10)),
    "numeric")

  # An all-missing matrix falls back rather than erroring differently.
  expect_error(detect_genotype_encoding(matrix(rep(NA_character_, 4L), 2)),
               "No non-NA")
})
