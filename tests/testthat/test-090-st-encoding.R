# =============================================================================
# test-090-st-encoding.R — C-5: the SelectionTools allele encoding is retired
# from the container (0.9.0). A GenoMeta stores nucleotide alleles only; the
# code table A = 1, C = 2, G = 3, T = 4, missing -1 is applied at the
# SelectionTools boundary (R/st_formats.R).
#
#   S1  the translation primitives
#   S2  the public detector refuses codes; the internal one still names them
#   S3  build_geno_meta(): undeclared codes refused, declared ones deprecated
#   S4  new_geno_meta(): the portal does the same
#   S5  the "ST" preset is gated and translates on the read
#   S6  st_convert(): files in are nucleotide; dict = "fixed" out
#   S7  upgrade_geno_meta() / load_geno_meta() recode a saved container
#   S8  the sim_* generators default to nucleotide; "selectiontools" warns
#   S9  recode_geno_meta(to_type = "selectiontools") is defunct; the recode
#       to nucleotide stays as the explicit route
#   S10 the legacy layer and the registry warn
#   S11 a PLINK 1/2-coded PED keeps its labels and its effects join
#   S12 C-6: the "selectiontools" token is behind use_st_formats() like the
#       files (builders, generators, legacy layer, registry); the one exception
#       is the migration of a saved container (upgrade / load / recode)
#   S13 C-6: the "ST" writer preset is gated; every recode is the adapter's table
# =============================================================================

library(testthat)
library(haploStack)

.s090_map <- function(markers, chrom = "1") {
  data.frame(Chrom = chrom, Marker = markers,
             Pos = seq_along(markers) * 1e6, stringsAsFactors = FALSE)
}

.s090_codes <- function() {
  matrix(c("1/1", "1/4", "4/4",
           "1/-1", "0/0", "2/3"),
         nrow = 2L, byrow = TRUE,
         dimnames = list(c("S1", "S2"), c("M1", "M2", "M3")))
}

# -- S1 -----------------------------------------------------------------------

test_that("S1: the code table translates tokens and cells", {
  f <- haploStack:::.st_to_nucleotide
  expect_identical(f(c("1", "2", "3", "4")), c("A", "C", "G", "T"))
  # a change of vocabulary, not a repair: only -1 (SelectionTools' missing
  # code) is rewritten, to the nucleotide missing character
  expect_identical(f(c("-1", "0", "-", NA, "")), c("-", "0", "-", NA, ""))
  expect_identical(f("5"), "5")                       # unknown: kept verbatim
  arr <- array(c("1", "4", "-1", "2"), dim = c(1, 2, 2),
               dimnames = list("S1", c("M1", "M2"), c("H1", "H2")))
  out <- f(arr)
  expect_identical(dim(out), dim(arr))
  expect_identical(dimnames(out), dimnames(arr))
  expect_identical(as.vector(out), c("A", "T", "-", "C"))

  g <- haploStack:::.st_decode_cells
  cells <- .s090_codes()
  out <- g(cells)
  expect_identical(dimnames(out), dimnames(cells))
  expect_identical(as.vector(out),
                   c("A/A", "A/-", "A/T", "-/-", "T/T", "C/G"))
  expect_identical(as.vector(g(matrix(c("14", "41", "00", "23"), 2L), split_char = "")),
                   c("AT", "TA", "--", "CG"))
  # a digit inside a longer token is not a code
  expect_identical(as.vector(g(matrix("12/1", 1L))), "12/A")
  # a character registered for the SelectionTools mode is missing too
  expect_identical(as.vector(g(matrix("#/1", 1L), missing_chars = "#")), "-/A")

  expect_true(haploStack:::.st_is_coded(c("1", "4", NA, "0", "-1")))
  expect_true(haploStack:::.st_is_coded(c("1", "2")))      # PLINK 1/2 looks the same
  expect_false(haploStack:::.st_is_coded(c("A", "T")))
  expect_false(haploStack:::.st_is_coded(c("1", "5")))
  expect_false(haploStack:::.st_is_coded(c(NA, "0")))
})

# -- S2 -----------------------------------------------------------------------

test_that("S2: detect_genotype_encoding() refuses codes, the internal detector names them", {
  codes <- matrix(c("1/1", "1/4", "4/4", "-1/-1"), nrow = 2)
  e <- tryCatch(detect_genotype_encoding(codes), error = function(e) e)
  expect_s3_class(e, "hs_error_format")
  expect_s3_class(e, "hs_error_input")
  expect_identical(e$hint, "st_convert()")
  expect_match(conditionMessage(e), "st_convert\\(\\)")
  expect_match(conditionMessage(e), "encoding = \"nucleotide\"")
  expect_identical(haploStack:::.detect_encoding(codes), "selectiontools")
  # the other verdicts are untouched
  expect_identical(detect_genotype_encoding(matrix(c("A/A", "A/T", "T/T", "G/C"), 2)),
                   "nucleotide")
  expect_identical(detect_genotype_encoding(matrix(c("AA", "AT", "TT", "AA"), 2)),
                   "nucleotide")
  expect_identical(detect_genotype_encoding(matrix(c("0", "1", "2", "0", "1", "2"), 3)),
                   "numeric")
  # split_haplotypes() reports the vocabulary of the cells, codes included
  sp <- split_haplotypes(matrix(c("1/1", "1/4"), 1L, dimnames = list("S1", c("M1", "M2"))))
  expect_identical(sp$encode, "selectiontools")
})

# -- S3 -----------------------------------------------------------------------

test_that("S3: build_geno_meta() refuses undeclared codes and translates declared ones", {
  cells <- .s090_codes()
  map   <- .s090_map(colnames(cells))

  e <- tryCatch(build_geno_meta(cells, map), error = function(e) e)
  expect_s3_class(e, "hs_error_format")
  expect_match(conditionMessage(e), "^build_geno_meta: ")

  # declared: deprecated, translated, stamped nucleotide (behind the switch, C-6)
  withr::local_options(list(haploStack.st_formats = TRUE))
  w <- NULL
  m <- withCallingHandlers(
    suppressMessages(build_geno_meta(cells, map, encoding = "selectiontools")),
    warning = function(cnd) {
      if (inherits(cnd, "hs_warning_deprecated")) {
        w <<- cnd
        invokeRestart("muffleWarning")
      }
    })
  expect_s3_class(w, "hs_warning_deprecated")
  expect_s3_class(w, "deprecatedWarning")
  expect_identical(w$new, "st_convert()")
  expect_match(conditionMessage(w), "new_geno_meta\\(encoding = \"selectiontools\"\\)")
  expect_identical(m$metadata$encoding, "nucleotide")
  expect_identical(as.vector(m$geno_arrays["S1", , ]), c("A", "A", "T", "A", "T", "T"))
  expect_true(all(is.na(m$geno_arrays["S2", "M2", ])))            # 0/0
  expect_identical(as.vector(m$geno_arrays["S2", "M1", ]), c("A", NA))  # 1/-1
  expect_setequal(names(m$metadata$char_dict), c("A", "C", "G", "T"))
  expect_true(is_geno_meta(m, verbose = FALSE))
  expect_true(all(c("0", "-") %in% m$metadata$missing_tokens))

  # digits declared as labels are kept as labels
  lab <- matrix(c("1/1", "1/2", "2/2", "0/0"), nrow = 2L, byrow = TRUE,
                dimnames = list(c("S1", "S2"), c("M1", "M2")))
  m2 <- suppressMessages(build_geno_meta(lab, .s090_map(c("M1", "M2")),
                                         encoding = "nucleotide"))
  expect_identical(m2$metadata$encoding, "nucleotide")
  expect_setequal(stats::na.omit(unique(as.vector(m2$geno_arrays))), c("1", "2"))
  expect_true(all(is.na(m2$geno_arrays["S2", "M2", ])))
  expect_no_error(geno_meta_views(m2))
})

# -- S4 -----------------------------------------------------------------------

test_that("S4: new_geno_meta() is the portal: same ladder, same refusal", {
  arr <- array(c("1", "4", "4", "1", "-1", "2", "0", "3"), dim = c(2, 2, 2),
               dimnames = list(c("S1", "S2"), c("M1", "M2"), c("H1", "H2")))
  map <- .s090_map(c("M1", "M2"))

  e <- tryCatch(new_geno_meta(arr, map, metadata = list(source_format = "x")),
                error = function(e) e)
  expect_s3_class(e, "hs_error_format")
  expect_match(conditionMessage(e), "^new_geno_meta: ")

  withr::local_options(list(haploStack.st_formats = TRUE))          # C-6
  expect_warning(
    m <- suppressMessages(new_geno_meta(arr, map,
                                        metadata = list(source_format = "x",
                                                        encoding = "selectiontools",
                                                        char_dict = c("1" = 1L, "4" = 2L)))),
    class = "hs_warning_deprecated")
  expect_identical(m$metadata$encoding, "nucleotide")
  expect_identical(as.vector(m$geno_arrays), c("A", "T", "T", "A", NA, "C", NA, "G"))
  expect_setequal(names(m$metadata$char_dict), c("A", "C", "G", "T"))

  # labels declared stay labels; a nucleotide array never touches the ladder
  expect_no_warning(suppressMessages(new_geno_meta(
    array(c("1", "2", "2", "1"), dim = c(1, 2, 2),
          dimnames = list("S1", c("M1", "M2"), c("H1", "H2"))),
    map, metadata = list(source_format = "x", encoding = "nucleotide"))))
})

# -- S5 -----------------------------------------------------------------------

test_that("S5: the ST preset is gated and translates on the read", {
  d <- tempfile("st_"); dir.create(d)
  gf <- file.path(d, "geno.csv"); mf <- file.path(d, "map.csv")
  writeLines(c("marker,S1,S2,S3",
               "M1,1/1,1/4,4/4",
               "M2,1/-1,0/0,2/3",
               "M3,4/4,4/1,1/1"), gf)
  writeLines(c("name,chrom,pos",
               "M1,1,1000000", "M2,1,2000000", "M3,2,3000000"), mf)

  withr::local_options(list(haploStack.st_formats = NULL))
  e <- tryCatch(read_genotype(gf, mf, preset = "ST"), error = function(e) e)
  expect_s3_class(e, "hs_error_input")
  expect_identical(e$hint, "use_st_formats(TRUE)")
  expect_match(conditionMessage(e), "^read_genotype: SelectionTools formats are off")
  expect_error(read_st_genotype(gf, mf), class = "hs_error_input")

  m <- suppressMessages(read_genotype(gf, mf, preset = "ST", st_formats = TRUE))
  expect_identical(m$metadata$encoding, "nucleotide")
  expect_identical(m$metadata$source_format, "selectiontools")
  expect_identical(as.vector(m$geno_arrays["S3", , ]), c("T", "C", "A", "T", "G", "A"))
  expect_identical(as.vector(m$geno_arrays["S1", "M2", ]), c("A", NA))
  expect_true(all(is.na(m$geno_arrays["S2", "M2", ])))
  expect_true(is_geno_meta(m, verbose = FALSE))

  withr::local_options(list(haploStack.st_formats = TRUE))
  m2 <- suppressMessages(read_st_genotype(gf, mf))
  expect_identical(m2$geno_arrays, m$geno_arrays)

  # a nucleotide file through the same preset layout is not gated
  withr::local_options(list(haploStack.st_formats = NULL))
  nf <- file.path(d, "nuc.csv")
  writeLines(c("marker,S1,S2", "M1,A/A,A/T", "M2,T/T,A/A", "M3,A/T,T/T"), nf)
  n <- suppressMessages(read_genotype(nf, mf, preset = "ST", encoding = "nucleotide"))
  expect_identical(n$metadata$encoding, "nucleotide")
})

# -- S6 -----------------------------------------------------------------------

test_that("S6: st_convert() gives a nucleotide container; dict = \"fixed\" writes the table", {
  withr::local_options(list(haploStack.st_formats = TRUE))
  d <- tempfile("stc_"); dir.create(d)
  gf <- file.path(d, "geno.txt"); mf <- file.path(d, "map.txt")
  writeLines(c("m/i\tS1\tS2",
               "M1\t1/1\t1/4",
               "M2\t-1/1\t2/3",
               "M3\t4/4\t-/-"), gf)
  writeLines(c("name\tchrom\tpos", "M1\t1\t1000000", "M2\t1\t2000000", "M3\t2\t3000000"), mf)

  m <- suppressMessages(st_convert(gf, map = mf))
  expect_identical(m$metadata$encoding, "nucleotide")
  expect_identical(m$metadata$source_format, "selectiontools")
  expect_identical(as.vector(m$geno_arrays["S1", , ]), c("A", NA, "T", "A", "A", "T"))
  expect_identical(as.vector(m$geno_arrays["S2", , ]), c("A", "C", NA, "T", "G", NA))
  expect_false(any(stats::na.omit(as.vector(m$geno_arrays)) %in% as.character(1:4)))

  # out, with the fixed table: A = 1, C = 2, G = 3, T = 4
  out <- suppressMessages(st_convert(m, to = "selectiontools", output_dir = d,
                                     prefix = "fx", dict = "fixed"))
  g <- utils::read.table(out$geno_file, header = TRUE, sep = "\t",
                         check.names = FALSE, stringsAsFactors = FALSE)
  expect_identical(g[[1L]], c("M1", "M2", "M3"))
  expect_identical(g$S1, c("1/1", "-/1", "4/4"))
  expect_identical(g$S2, c("1/4", "2/3", "-/-"))
  # the per-marker default stays what the resolver needs
  out2 <- suppressMessages(st_convert(m, to = "selectiontools", output_dir = d, prefix = "pm"))
  g2 <- utils::read.table(out2$geno_file, header = TRUE, sep = "\t",
                          check.names = FALSE, stringsAsFactors = FALSE)
  expect_identical(g2$S2, c("1/2", "2/3", "-/-"))        # M2 holds A, C, G
  # fixed refuses labels outside the table
  lab <- suppressMessages(build_geno_meta(
    matrix(c("1/1", "1/2"), 1L, dimnames = list("S1", c("M1", "M2"))),
    .s090_map(c("M1", "M2")), encoding = "nucleotide"))
  expect_error(st_convert(lab, to = "selectiontools", output_dir = d, dict = "fixed"),
               "dict = \"fixed\"")
})

# -- S7 -----------------------------------------------------------------------

test_that("S7: upgrade_geno_meta() and load_geno_meta() recode a saved container", {
  m <- sim_geno_meta(4, 6, n_chrom = 2L, block_size = 3L, seed = 7L)
  eff <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("A", "T"), seed = 7L)
  m$marker_effects <- list(Y = eff)
  # as 0.8.8 (or a 0.9.0 development build) saved it: codes, the old token
  old <- m
  old$geno_arrays[] <- chartr("AT", "14", m$geno_arrays)
  old$metadata$encoding  <- "selectiontools"
  old$metadata$char_dict <- list("1" = 1L, "4" = 2L)
  old$marker_effects$Y$Allele <- chartr("AT", "14", old$marker_effects$Y$Allele)
  old$metadata$schema_version <- NULL                  # schema 1 ...

  expect_message(up <- upgrade_geno_meta(old), "recoded to nucleotides")
  expect_identical(up$metadata$encoding, "nucleotide")
  expect_identical(up$geno_arrays, m$geno_arrays)
  expect_identical(up$marker_effects$Y$Allele, m$marker_effects$Y$Allele)
  expect_setequal(names(up$metadata$char_dict), c("A", "T"))
  expect_true(is_geno_meta(up, verbose = FALSE, strict = TRUE))

  old2 <- old; old2$metadata$schema_version <- 2L       # ... and schema 2
  up2 <- suppressMessages(upgrade_geno_meta(old2))
  expect_identical(up2$geno_arrays, m$geno_arrays)
  expect_identical(up2$metadata$encoding, "nucleotide")

  f <- tempfile(fileext = ".rds")
  saveRDS(old, f)
  ld <- suppressMessages(suppressWarnings(load_geno_meta(f)))
  expect_identical(ld$geno_arrays, m$geno_arrays)
  expect_identical(ld$metadata$encoding, "nucleotide")
  expect_no_error(geno_meta_views(ld))
})

# -- S8 -----------------------------------------------------------------------

test_that("S8: the simulator emits nucleotides; the old encoding warns and gives the same draws", {
  m <- sim_geno_meta(5, 8, n_chrom = 2L, block_size = 4L, seed = 12L)
  expect_identical(m$metadata$encoding, "nucleotide")
  expect_setequal(stats::na.omit(unique(as.vector(m$geno_arrays))), c("A", "T"))
  expect_identical(as.vector(m$geno_arrays[1, , ]),
    c("T", "T", "A", "A", "A", "T", "T", "A", "T", "A", "T", "T",
      "T", "A", "A", "T"))                             # the 0.8.8 draw, 1 -> A, 4 -> T

  withr::local_options(list(haploStack.st_formats = TRUE))          # C-6
  expect_warning(ms <- sim_geno_meta(5, 8, "selectiontools", n_chrom = 2L,
                                     block_size = 4L, seed = 12L),
                 class = "hs_warning_deprecated")
  expect_identical(ms$geno_arrays, m$geno_arrays)
  expect_identical(ms$metadata$encoding, "nucleotide")

  g <- sim_geno_matrix(3, 4, seed = 2L)
  expect_true(all(grepl("^[AT]/[AT]$", g)))
  expect_warning(gs <- sim_geno_matrix(3, 4, "selectiontools", seed = 2L),
                 class = "hs_warning_deprecated")
  expect_identical(chartr("14", "AT", gs), g)

  ds <- sim_dataset(4, 6, block_size = 3L, seed = 14L)
  expect_identical(ds$geno_meta$metadata$encoding, "nucleotide")
  expect_identical(as.vector(ds$geno_meta$geno_arrays[, , 1L]),
                   as.vector(t(chartr("14", "AT", sub("/.*$", "", ds$geno_st)))))
  expect_true(all(grepl("^[14]/[14]$", ds$geno_st)))    # the legacy element stays
})

# -- S9 -----------------------------------------------------------------------

test_that("S9: recode to the old encoding is defunct; recode from it stays", {
  m <- sim_geno_meta(3, 4, seed = 3L)
  e <- tryCatch(recode_geno_meta(m, to_type = "selectiontools"), error = function(e) e)
  expect_s3_class(e, "hs_error_defunct")
  expect_s3_class(e, "defunctError")
  expect_identical(e$new, "st_convert(to = \"selectiontools\")")
  expect_match(conditionMessage(e), "haploStack-defunct")

  old <- m
  old$geno_arrays[] <- chartr("AT", "14", m$geno_arrays)
  old$metadata$encoding <- "selectiontools"
  rc <- suppressMessages(recode_geno_meta(old, to_type = "nucleotide"))
  expect_identical(rc$geno_arrays, m$geno_arrays)
  expect_identical(rc$metadata$encoding, "nucleotide")
  # the explicit source works too, and the data check still guards it
  expect_identical(suppressMessages(recode_geno_meta(old, "nucleotide",
                                                     from_type = "selectiontools"))$geno_arrays,
                   m$geno_arrays)
  expect_error(recode_geno_meta(m, "nucleotide", from_type = "selectiontools"),
               "alphabetic")
})

# -- S10 ----------------------------------------------------------------------

test_that("S10: the legacy layer and the registry warn on the old token", {
  withr::local_options(list(haploStack.st_formats = TRUE))          # C-6
  arr <- array(c("A", "T"), dim = c(1, 1, 2))
  expect_warning(out <- translate_alleles(arr, "nucleotide", "selectiontools"),
                 class = "hs_warning_deprecated")
  expect_identical(as.vector(out), c("1", "4"))
  expect_no_warning(translate_alleles(arr, "nucleotide", "nucleotide"))

  eff <- data.frame(Marker = "M1", Allele = c("A", "T"), Effect = c(1, -1))
  expect_warning(te <- transform_effects(eff, input_encode = "nucleotide",
                                         output_encode = "selectiontools"),
                 class = "hs_warning_deprecated")
  expect_identical(te$Allele, c("1", "4"))

  old <- getOption("haploStack.missing_chars")
  on.exit(options(haploStack.missing_chars = old), add = TRUE)
  expect_warning(suppressMessages(use_missing_chars("#", mode = "selectiontools")),
                 class = "hs_warning_deprecated")
  expect_no_warning(suppressMessages(use_missing_chars("#", mode = "nucleotide")))
  reset_missing_chars()
})

# -- S11 ----------------------------------------------------------------------

test_that("S11: a PLINK 1/2-coded PED keeps its labels, and effects on them join", {
  d <- tempfile("ped12_"); dir.create(d)
  pf <- file.path(d, "x.ped"); mf <- file.path(d, "x.map")
  writeLines(c("F1 S1 0 0 1 -9 1 1 1 2 2 2 1 1",
               "F1 S2 0 0 2 -9 1 2 2 2 0 0 1 2",
               "F2 S3 0 0 1 -9 2 2 1 1 1 2 2 2"), pf)
  writeLines(c("1 M1 0 1000000", "1 M2 0 2000000", "1 M3 0 3000000", "2 M4 0 4000000"), mf)
  bmap <- data.frame(Chrom = c("1", "1", "1", "2"), Marker = paste0("M", 1:4),
                     Block = c("B1", "B1", "B2", "B3"), stringsAsFactors = FALSE)
  m <- suppressMessages(read_plink_ped(pf, mf, block_file = bmap))
  expect_identical(m$metadata$encoding, "nucleotide")
  expect_setequal(stats::na.omit(unique(as.vector(m$geno_arrays))), c("1", "2"))
  expect_true(all(is.na(m$geno_arrays["S2", "M3", ])))
  eff <- data.frame(Marker = rep(paste0("M", 1:4), each = 2L),
                    Allele = rep(c("1", "2"), 4L),
                    Effect = c(1, -1, 0.5, -0.5, 2, -2, 0.1, -0.1))
  expect_no_warning(g <- suppressMessages(compute_profile_localgebv(m, eff, center = "none")))
  expect_true(is.data.frame(g))
  vals <- as.data.frame(g)[, setdiff(names(g), "Block"), drop = FALSE]
  expect_true(any(abs(as.matrix(vals)) > 0, na.rm = TRUE))
})

# -- S12 ----------------------------------------------------------------------

test_that("S12: the token is behind use_st_formats() like the files (C-6)", {
  withr::local_options(list(haploStack.st_formats = NULL))
  cells <- .s090_codes()
  map   <- .s090_map(colnames(cells))
  off <- function(expr, caller) {
    e <- tryCatch(expr, error = function(e) e)
    expect_s3_class(e, "hs_error_input")
    expect_identical(e$hint, "use_st_formats(TRUE)")
    expect_match(conditionMessage(e),
                 paste0("^", caller, ": SelectionTools formats are off"))
    # the token takes no per-call argument: the session switch is the only way in
    expect_false(grepl("st_formats = TRUE", conditionMessage(e), fixed = TRUE))
  }
  off(build_geno_meta(cells, map, encoding = "selectiontools"), "build_geno_meta")
  arr <- array(c("1", "4", "4", "1"), dim = c(1, 2, 2),
               dimnames = list("S1", c("M1", "M2"), c("H1", "H2")))
  off(new_geno_meta(arr, .s090_map(c("M1", "M2")),
                    metadata = list(source_format = "x", encoding = "selectiontools")),
      "new_geno_meta")
  off(sim_geno_meta(3, 4, "selectiontools", seed = 1L), "sim_geno_meta")
  off(sim_geno_matrix(3, 4, "selectiontools", seed = 1L), "sim_geno_matrix")
  nuc <- array(c("A", "T"), dim = c(1, 1, 2))
  off(translate_alleles(nuc, "nucleotide", "selectiontools"), "translate_alleles")
  eff <- data.frame(Marker = "M1", Allele = c("A", "T"), Effect = c(1, -1))
  off(transform_effects(eff, input_encode = "nucleotide",
                        output_encode = "selectiontools"),
      "transform_effects")
  old_reg <- getOption("haploStack.missing_chars")
  on.exit(options(haploStack.missing_chars = old_reg), add = TRUE)
  off(use_missing_chars("#", mode = "selectiontools"), "use_missing_chars")

  # nothing else is touched: label and nucleotide routes run without the switch
  expect_no_error(suppressMessages(build_geno_meta(
    matrix(c("1/1", "1/2"), 1L, dimnames = list("S1", c("M1", "M2"))),
    .s090_map(c("M1", "M2")), encoding = "nucleotide")))
  expect_no_error(translate_alleles(nuc, "nucleotide", "nucleotide"))
  expect_no_error(suppressMessages(use_missing_chars("#", mode = "nucleotide")))
  reset_missing_chars()

  # the exception: a container saved with the old vocabulary still opens and recodes
  m <- sim_geno_meta(3, 4, seed = 3L)
  old <- m
  old$geno_arrays[] <- chartr("AT", "14", m$geno_arrays)
  old$metadata$encoding  <- "selectiontools"
  old$metadata$char_dict <- list("1" = 1L, "4" = 2L)
  expect_identical(suppressMessages(upgrade_geno_meta(old))$geno_arrays, m$geno_arrays)
  expect_identical(suppressMessages(recode_geno_meta(old, to_type = "nucleotide"))$geno_arrays,
                   m$geno_arrays)
  f <- tempfile(fileext = ".rds")
  saveRDS(old, f)
  expect_identical(suppressMessages(suppressWarnings(load_geno_meta(f)))$geno_arrays,
                   m$geno_arrays)

  # with the switch on, the deprecated behaviour is what C-5 built
  withr::local_options(list(haploStack.st_formats = TRUE))
  expect_warning(m2 <- suppressMessages(build_geno_meta(cells, map, encoding = "selectiontools")),
                 class = "hs_warning_deprecated")
  expect_identical(m2$metadata$encoding, "nucleotide")
  expect_warning(g <- sim_geno_matrix(3, 4, "selectiontools", seed = 2L),
                 class = "hs_warning_deprecated")
  expect_true(all(grepl("^[14]/[14]$", g)))
})

# -- S14 (C-6 item 5): the dead "selectiontools" branches are gone -----------

test_that("S14: no code path past the portal understands the retired token", {
  expect_false("selectiontools" %in% names(BASE_COLOR_MAP))
  expect_false("selectiontools" %in% names(COLOR_MAP))
  # the colour-map builder falls back to the nucleotide palette for it
  expect_identical(make_allele_color_map(c("A", "T"), "selectiontools"),
                   make_allele_color_map(c("A", "T"), "nucleotide"))
  m <- sim_geno_meta(3, 4, n_chrom = 1L, seed = 5L)
  expect_error(geno_meta_views(m, encoding = "selectiontools"), "arg")
  bad <- m; bad$metadata$encoding <- "selectiontools"
  expect_error(haploStack:::.gp_allele_dosage(bad), "'nucleotide'")
  # the exported constants keep the ladder to 0.10
  expect_true("selectiontools" %in% HS_ENCODINGS)
  expect_true("selectiontools" %in% HS_MISSING_MODES)
})

# -- S13 ----------------------------------------------------------------------

test_that("S13: the ST writer preset is gated; every recode is the adapter's table", {
  withr::local_options(list(haploStack.st_formats = NULL))
  m <- sim_geno_meta(3, 4, n_chrom = 1L, seed = 5L)
  e <- tryCatch(geno_meta_to_matrix(m, preset = "ST"), error = function(e) e)
  expect_s3_class(e, "hs_error_input")
  expect_identical(e$hint, "use_st_formats(TRUE)")
  expect_match(conditionMessage(e), "^geno_meta_to_matrix: .*st_formats = TRUE")
  f <- tempfile(fileext = ".csv")
  expect_error(write_geno_meta(m, f, preset = "ST"), class = "hs_error_input")
  expect_false(file.exists(f))
  # the other presets are not gated
  expect_no_error(geno_meta_to_matrix(m, preset = "GS"))

  out <- geno_meta_to_matrix(m, preset = "ST", st_formats = TRUE)
  expect_identical(dim(out$matrix), c(4L, 3L))
  expect_true(all(grepl("^[AT]/[AT]$", out$matrix)))      # the layout, nucleotides
  expect_no_error(suppressMessages(write_geno_meta(m, f, preset = "ST", st_formats = TRUE)))
  expect_true(file.exists(f))
  withr::local_options(list(haploStack.st_formats = TRUE))
  expect_identical(geno_meta_to_matrix(m, preset = "ST")$matrix, out$matrix)

  # one table: the legacy layer, recode_geno_meta() and the adapter agree
  codes <- array(c("1", "2", "3", "4", "-1", "0"), dim = c(1, 3, 2),
                 dimnames = list("S1", c("M1", "M2", "M3"), c("H1", "H2")))
  via_legacy <- suppressWarnings(translate_alleles(codes, "selectiontools", "nucleotide"))
  expect_identical(as.vector(via_legacy), as.vector(haploStack:::.st_to_nucleotide(codes)))
  expect_identical(as.vector(via_legacy), c("A", "C", "G", "T", "-", "0"))
  back <- suppressWarnings(translate_alleles(via_legacy, "nucleotide", "selectiontools"))
  expect_identical(as.vector(back), c("1", "2", "3", "4", "-1", "0"))
  expect_identical(
    suppressWarnings(haploStack:::.eff_build_state_map("selectiontools", "nucleotide", NULL)),
    haploStack:::.ST_TO_NUCLEOTIDE)
  expect_identical(
    suppressWarnings(haploStack:::.eff_build_state_map("nucleotide", "selectiontools", NULL)),
    haploStack:::.ST_FROM_NUCLEOTIDE)
})
