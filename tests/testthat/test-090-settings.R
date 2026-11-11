# =============================================================================
# test-090-settings.R — settings as arguments, option as default (design D1-3)
#
#   S1  .hs_setting(): argument > option > default; an argument is validated
#       like the setter; an unusable option is ignored, never obeyed
#   S2  the missing-character registry reaches a container only when it is
#       built: a registration in "job A" does not change "job B" read with
#       missing_chars = character(0); missing_chars = "#" works without it
#   S3  a container built before a registration keeps its stamp: the later
#       registration changes neither validation nor the numbers; repair is
#       the route, and honours both the registry and its own argument
#   S4  strict_missing as an argument: errors without touching the option, and
#       the argument overrides an option set the other way
#   S5  missing_chars is refused where use_missing_chars() refuses it
#   S6  fill_na and sample_order resolve through the same rule
# =============================================================================

.s090_cells <- function(token = "#") {
  geno <- matrix(c("A/A", paste0(token, "/", token), "A/T", "T/T",
                   "T/T", "A/A", "A/T", "A/A"),
                 nrow = 2, byrow = TRUE,
                 dimnames = list(c("S1", "S2"), c("M1", "M2", "M3", "M4")))
  map <- data.frame(Chrom = "1", Marker = colnames(geno),
                    Pos = c(1e6, 2e6, 3e6, 4e6), stringsAsFactors = FALSE)
  list(geno = geno, map = map)
}
.s090_build <- function(...) {
  d <- .s090_cells()
  suppressMessages(build_geno_meta(d$geno, d$map, split_char = "/",
                                   encoding = "nucleotide", ...))
}

test_that("S1: argument > option > default, validated like the setters", {
  withr::local_options(list(haploStack.strict_missing = NULL,
                            haploStack.fill_na = NULL))
  st <- haploStack:::.hs_setting
  expect_false(st("strict_missing"))                       # built-in default
  withr::local_options(list(haploStack.strict_missing = TRUE))
  expect_true(st("strict_missing"))                        # option
  expect_false(st("strict_missing", FALSE))                # argument wins
  expect_error(st("strict_missing", "yes"), class = "hs_error_input")
  expect_identical(st("fill_na"), HS_FILL_NA_DEFAULT)
  expect_identical(st("fill_na", "?"), "?")
  expect_error(st("fill_na", ""), "empty string", class = "hs_error_input")
  expect_error(st("fill_na", c("a", "b")), class = "hs_error_input")
  expect_error(st("no_such_setting"), class = "hs_error_internal")
  # an unusable option is ignored, not obeyed: the registry must be a list
  withr::local_options(list(haploStack.missing_chars = "nonsense"))
  expect_null(st("missing_chars"))
  expect_setequal(missing_chars("nucleotide"), c("", "NA", "N", "n", ".", "-", "0"))
  withr::local_options(list(haploStack.missing_chars = list(nucleotide = c("##", NA))))
  expect_setequal(missing_chars("nucleotide"), c("", "NA", "N", "n", ".", "-", "0"))
  # a per-mode list argument reads like the registry
  expect_identical(haploStack:::.extra_missing_chars("nucleotide", list(nucleotide = "#")), "#")
  expect_identical(haploStack:::.extra_missing_chars("numeric", list(nucleotide = "#")),
                   character(0))
})

test_that("S2: a registration in job A does not reach job B", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  # job A registers "#"; its build sees it
  suppressMessages(use_missing_chars("#", mode = "nucleotide"))
  a <- .s090_build()
  expect_true(all(is.na(a$geno_arrays["S1", "M2", ])))
  expect_true("#" %in% a$metadata$missing_tokens)
  # job B reads with an explicit empty set: "#" stays an allele, no stamp of it
  b <- .s090_build(missing_chars = character(0))
  expect_identical(unname(b$geno_arrays["S1", "M2", ]), c("#", "#"))
  expect_false("#" %in% b$metadata$missing_tokens)
  # and with nothing registered, the argument alone does the job
  reset_missing_chars()
  c1 <- .s090_build(missing_chars = "#")
  expect_true(all(is.na(c1$geno_arrays["S1", "M2", ])))
  expect_true("#" %in% c1$metadata$missing_tokens)
  # the same argument on new_geno_meta() and repair_geno_meta()
  c2 <- suppressMessages(new_geno_meta(b$geno_arrays, b$map_info,
                                       metadata = list(source_format = "t", encoding = "nucleotide"),
                                       missing_chars = "#"))
  expect_true(all(is.na(c2$geno_arrays["S1", "M2", ])))
  c3 <- suppressMessages(repair_geno_meta(b, missing_chars = "#"))
  expect_true(all(is.na(c3$geno_arrays["S1", "M2", ])))
  expect_true("#" %in% c3$metadata$missing_tokens)
})

test_that("S3: a container keeps the vocabulary it was built with", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  m <- .s090_build()                                   # "#" is an allele here
  expect_identical(unname(m$geno_arrays["S1", "M2", ]), c("#", "#"))
  before <- compute_marker_stats(m)
  suppressMessages(use_missing_chars("#", mode = "nucleotide"))
  # the later registration changes neither validation nor the numbers
  expect_silent(is_geno_meta(m))
  expect_identical(compute_marker_stats(m), before)
  expect_true(is_geno_meta(m))
  expect_true("#" %in% before$Alleles[before$Marker == "M2"] ||
                any(grepl("#", before[before$Marker == "M2", ], fixed = TRUE)))
  # repair is the route, and it reads the registry
  r <- suppressMessages(repair_geno_meta(m))
  expect_true(all(is.na(r$geno_arrays["S1", "M2", ])))
  # a derived container (filter) keeps the original stamp, not the registry
  f <- suppressMessages(suppressWarnings(filter_geno_meta(m, keep_markers = c("M1", "M2"))))
  expect_identical(unname(f$geno_arrays["S1", "M2", ]), c("#", "#"))
  expect_false("#" %in% f$metadata$missing_tokens)
  reset_missing_chars()
})

test_that("S4: strict_missing is an argument, the option only its default", {
  withr::local_options(list(haploStack.strict_missing = NULL))
  m <- .s090_build()
  m$geno_arrays["S2", "M1", ] <- "0"
  m$metadata$missing_tokens <- NULL                    # legacy: unstamped
  expect_warning(is_geno_meta(m), "missing token")
  expect_error(is_geno_meta(m, strict_missing = TRUE), "missing token",
               class = "hs_error_input")
  expect_null(getOption("haploStack.strict_missing"))  # untouched
  withr::local_options(list(haploStack.strict_missing = TRUE))
  expect_error(is_geno_meta(m), "missing token")
  expect_warning(is_geno_meta(m, strict_missing = FALSE), "missing token")
  expect_error(is_geno_meta(m, strict_missing = "no"), class = "hs_error_input")
})

test_that("S5: missing_chars is validated as use_missing_chars() validates", {
  withr::local_options(list(haploStack.missing_chars = NULL))
  expect_error(.s090_build(missing_chars = "A"), "real allele", class = "hs_error_input")
  expect_error(.s090_build(missing_chars = "##"), "exactly one character",
               class = "hs_error_input")
  expect_error(.s090_build(missing_chars = NA_character_), "NA", class = "hs_error_input")
  expect_error(.s090_build(missing_chars = 9), class = "hs_error_input")
  d <- .s090_cells()
  expect_error(suppressMessages(build_geno_meta(d$geno, d$map, split_char = "/",
                                                encoding = "numeric", missing_chars = "0")),
               "real allele", class = "hs_error_input")
  # the registry and the setter agree on the effective set
  suppressMessages(use_missing_chars(c("#", "?"), mode = "nucleotide"))
  expect_setequal(missing_chars("nucleotide", chars_only = TRUE), c("#", "?"))
  expect_setequal(haploStack:::.extra_missing_chars("nucleotide"), c("#", "?"))
  reset_missing_chars()
})

test_that("S6: fill_na and sample_order resolve by the same rule", {
  withr::local_options(list(haploStack.fill_na = NULL, haploStack.sample_order = NULL))
  expect_identical(haploStack:::.hs_fill_na(), HS_FILL_NA_DEFAULT)
  withr::local_options(list(haploStack.fill_na = "?"))
  expect_identical(haploStack:::.hs_fill_na(), "?")
  expect_identical(haploStack:::.hs_fill_na("*"), "*")
  expect_error(haploStack:::.hs_fill_na(""), "empty string")
  expect_null(haploStack:::.hs_sample_order())
  withr::local_options(list(haploStack.sample_order = c("b", "a")))
  expect_identical(haploStack:::.hs_sample_order(), c("b", "a"))
  withr::local_options(list(haploStack.sample_order = 42))
  expect_warning(x <- haploStack:::.hs_sample_order(), "unusable")
  expect_null(x)
  expect_true(haploStack:::.hs_setting("use_rcpp"))
  withr::local_options(list(haploStack.use_rcpp = FALSE))
  expect_false(haploStack:::.hs_use_rcpp("local_gebv_block_sum_cpp"))
})
