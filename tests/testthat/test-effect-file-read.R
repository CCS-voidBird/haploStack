# =============================================================================
# test-effect-file-read.R  (0.7.8)
#
# read_effect_file() with sep = "auto", header = "auto":
#   - E1  export_marker_effects() output round-trips with the defaults — the
#         contract its roxygen always promised (intake 01)
#   - E2  the long-standing layout (space-separated, no header) is read exactly
#         as the pre-0.7.8 defaults (sep = " ", header = FALSE) read it
#   - E3  tab / comma, with / without a header, all give the same table
#   - E4  a first data row with an NA effect is data, not a header
#   - E5  the separator is decided on the first line, so "|" / ":" inside a
#         DArT-style marker ID cannot split a space-separated file
#   - E6  explicit arguments keep their old meaning
#   - E7  a header-only file and a bad 'header' argument fail loudly
# =============================================================================

library(testthat)
library(haploStack)

eff_lines <- function(file, lines) { writeLines(lines, file); file }

test_that("E1: export_marker_effects() output round-trips with the defaults", {
  m  <- sim_geno_meta(n_samples = 10L, n_markers = 20L, n_chrom = 2L, block_size = 4L)
  mk <- dimnames(m$geno_arrays)[[2L]]
  m$marker_effects <- list(Yield = sim_effect_df(mk, alleles = c("1", "4"), seed = 1L))
  out <- tempfile("fx_"); dir.create(out)
  invisible(export_marker_effects(m, traits = "Yield", output_dir = out))
  f <- file.path(out, "marker_effects_Yield.txt")
  expect_identical(readLines(f, n = 1L), "Marker\tAllele\tEffect")   # tab + header

  back <- read_effect_file(f)
  stored <- as.data.frame(m$marker_effects$Yield)
  expect_identical(nrow(back), nrow(stored))
  key_b <- paste(back$Marker, back$Allele)
  key_s <- paste(stored$Marker, stored$Allele)
  expect_setequal(key_b, key_s)
  expect_equal(back$Effect[match(key_s, key_b)], stored$Effect)
})

test_that("E2: a space-separated file without a header reads as the old defaults did", {
  f <- eff_lines(tempfile(fileext = ".txt"),
                 c("M1 A 0.1", "M1 T -0.1", "M2 C 0.3", "M2 G -0.3"))
  expect_identical(read_effect_file(f), read_effect_file(f, sep = " ", header = FALSE))

  f4 <- eff_lines(tempfile(fileext = ".txt"), c("M1 A 0.1 1.5", "M1 T -0.1 -1.5"))
  expect_identical(read_effect_file(f4), read_effect_file(f4, sep = " ", header = FALSE))
  expect_identical(ncol(read_effect_file(f4)), 4L)                  # extra trait kept
})

test_that("E3: tab / comma, with or without a header, give the same table", {
  ref <- read_effect_file(eff_lines(tempfile(), c("M1 A 0.1", "M1 T -0.1", "M2 C 0.3")))
  variants <- list(
    tab_noheader   = c("M1\tA\t0.1", "M1\tT\t-0.1", "M2\tC\t0.3"),
    tab_header     = c("Marker\tAllele\tEffect", "M1\tA\t0.1", "M1\tT\t-0.1", "M2\tC\t0.3"),
    comma_header   = c("Marker,Allele,Effect", "M1,A,0.1", "M1,T,-0.1", "M2,C,0.3"),
    comma_noheader = c("M1,A,0.1", "M1,T,-0.1", "M2,C,0.3"),
    space_header   = c("Marker Allele Effect", "M1 A 0.1", "M1 T -0.1", "M2 C 0.3"),
    renamed_header = c("snp allele beta", "M1 A 0.1", "M1 T -0.1", "M2 C 0.3"))
  for (v in names(variants)) {
    got <- read_effect_file(eff_lines(tempfile(), variants[[v]]))
    expect_identical(got, ref, info = v)
  }
})

test_that("E4: a first data row with an NA effect is data, not a header", {
  f <- eff_lines(tempfile(), c("M1 A NA", "M1 T -0.1", "M2 C 0.3"))
  got <- read_effect_file(f)
  expect_true("M1" %in% got$Marker)
  expect_true(is.na(got$Effect[got$Marker == "M1" & got$Allele == "A"]))
  expect_identical(got, read_effect_file(f, sep = " ", header = FALSE))

  fh <- eff_lines(tempfile(), c("Marker\tAllele\tEffect", "M1\tA\tNA", "M1\tT\t-0.1"))
  goth <- read_effect_file(fh)
  expect_identical(nrow(goth), 2L)
  expect_true(is.na(goth$Effect[goth$Allele == "A"]))
})

test_that("E5: '|' or ':' inside a marker ID does not split a space-separated file", {
  f <- eff_lines(tempfile(), c("100013|F|0-12:A>G A 0.5", "100021|F|0-9:C>T T -0.2"))
  got <- read_effect_file(f)
  expect_identical(ncol(got), 3L)
  expect_setequal(unique(got$Marker), c("100013|F|0-12:A>G", "100021|F|0-9:C>T"))
  expect_equal(got$Effect[got$Marker == "100013|F|0-12:A>G" & got$Allele == "A"], 0.5)
})

test_that("E6: explicit sep / header keep their old meaning", {
  f <- eff_lines(tempfile(), c("Marker Allele Effect", "M1 A 0.1", "M1 T -0.1"))
  asked_false <- read_effect_file(f, sep = " ", header = FALSE)   # header line is data, as always
  expect_true("Marker" %in% asked_false$Marker)
  expect_true(is.na(asked_false$Effect[asked_false$Marker == "Marker"][1L]))
  expect_false("Marker" %in% read_effect_file(f)$Marker)          # auto drops it

  ft <- eff_lines(tempfile(), c("Marker\tAllele\tEffect", "M1\tA\t0.1"))
  expect_identical(read_effect_file(ft, sep = "\t", header = TRUE), read_effect_file(ft))
})

test_that("E7: a header-only file and a bad 'header' argument fail loudly", {
  f <- eff_lines(tempfile(), "Marker\tAllele\tEffect")
  expect_error(read_effect_file(f), "header line and no effects")
  g <- eff_lines(tempfile(), c("M1 A 0.1"))
  expect_error(read_effect_file(g, header = "yes"), "must be TRUE, FALSE or \"auto\"")
})
