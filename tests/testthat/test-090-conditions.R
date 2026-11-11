# =============================================================================
# test-090-conditions.R — 0.9.0 condition vocabulary (design D1-1)
#
#   C1  .hs_abort() defaults to hs_error_input; a named kind replaces it;
#       input subkinds stack
#   C2  `data` fields reach the handler; reserved and unnamed fields refused
#   C3  every dependency guard in R/ carries an install hint
#   C4  text and call unchanged for three sampled errors (0.8.8 wording)
#   C5  assertion sites raise hs_error_internal, not an input error
#   C6  a condition object passes through untouched
# =============================================================================

test_that("C1: .hs_abort() defaults to hs_error_input; other kinds replace it", {
  e <- tryCatch(haploStack:::.hs_abort("x"), error = identity)
  expect_identical(class(e), c("hs_error_input", "hs_error", "error", "condition"))

  for (k in c("hs_error_dependency", "hs_error_state", "hs_error_cancelled",
              "hs_error_internal")) {
    e <- tryCatch(haploStack:::.hs_abort(class = k, "x"), error = identity)
    expect_identical(class(e)[1:2], c(k, "hs_error"))
    expect_false(inherits(e, "hs_error_input"))
  }

  e <- tryCatch(haploStack:::.hs_abort(class = c("hs_error_block_map", "hs_error_input"), "x"),
                error = identity)
  expect_true(inherits(e, "hs_error_block_map"))
  expect_true(inherits(e, "hs_error_input"))
  expect_true(inherits(e, "hs_error"))

  # warnings keep an empty default class
  w <- tryCatch(haploStack:::.hs_warn("w"), warning = identity)
  expect_identical(class(w), c("hs_warning", "warning", "condition"))
})

test_that("C2: data fields reach the handler; reserved and unnamed fields refused", {
  e <- tryCatch(haploStack:::.hs_abort("x", data = list(hint = "h", partial = list(1L))),
                error = identity)
  expect_identical(e$hint, "h")
  expect_identical(e$partial, list(1L))
  expect_identical(conditionMessage(e), "x")

  w <- tryCatch(haploStack:::.hs_warn("w", class = "hs_warning_fallback",
                                      data = list(stage = "s", code = "c")),
                warning = identity)
  expect_identical(class(w)[1:2], c("hs_warning_fallback", "hs_warning"))
  expect_identical(w$stage, "s")
  expect_identical(w$code, "c")
  expect_identical(conditionMessage(w), "w")

  expect_error(haploStack:::.hs_abort("x", data = list(message = "no")),
               class = "hs_error_internal")
  expect_error(haploStack:::.hs_abort("x", data = list(call = "no")),
               class = "hs_error_internal")
  expect_error(haploStack:::.hs_abort("x", data = list("unnamed")),
               class = "hs_error_internal")
  # a data problem that does not stop the call
  expect_warning(haploStack:::.hs_warn("w", class = "hs_warning_input"),
                 class = "hs_warning_input")
})

test_that("C3: every dependency guard carries an install hint", {
  rdir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(rdir), "package R/ sources not reachable from the test dir")
  files <- list.files(rdir, pattern = "[.]R$", full.names = TRUE)
  n_sites <- 0L
  for (f in files) {
    src <- readLines(f, warn = FALSE)
    hits <- grep('class = "hs_error_dependency"', src, fixed = TRUE)
    for (i in hits) {
      n_sites <- n_sites + 1L
      window <- paste(src[i:min(length(src), i + 2L)], collapse = "\n")
      expect_true(grepl("data = list(hint = ", window, fixed = TRUE),
                  label = sprintf("%s:%d has an install hint", basename(f), i))
    }
  }
  expect_gte(n_sites, 9L)
})

test_that("C4: text and call unchanged for sampled errors (0.8.8 wording)", {
  e <- tryCatch(fit_genomic_prediction(list()), error = identity)
  expect_s3_class(e, "hs_error_input")
  # phase B (D2-3): the guard keeps the 0.8.8 prefix and appends the reason
  expect_match(conditionMessage(e),
               "^fit_genomic_prediction: 'geno_meta' is not a valid GenoMeta: ")
  expect_identical(as.character(conditionCall(e)[[1L]]), "fit_genomic_prediction")

  e <- tryCatch(ocs_selection(c(1, 2)), error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_identical(conditionMessage(e), "ocs_selection: 'gebv_vector' must be named.")
  expect_identical(as.character(conditionCall(e)[[1L]]), "ocs_selection")

  e <- tryCatch(read_phenotype(tempfile("no_such_")), error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "^read_phenotype: file not found: ")
  expect_identical(as.character(conditionCall(e)[[1L]]), "read_phenotype")
})

test_that("C5: an assertion failure is hs_error_internal, not an input error", {
  e <- tryCatch(haploStack:::.st_parse_blocks(data.frame(a = 1), NULL), error = identity)
  expect_s3_class(e, "hs_error_internal")
  expect_false(inherits(e, "hs_error_input"))
  expect_match(conditionMessage(e), "^st_ld_blocks: unexpected SelectionTools block table")
})

test_that("C6: a condition object passes through .hs_abort()/.hs_warn() untouched", {
  cnd <- errorCondition("custom", class = "my_error")
  e <- tryCatch(haploStack:::.hs_abort(cnd), error = identity)
  expect_identical(class(e), c("my_error", "error", "condition"))
  wc <- warningCondition("custom w", class = "my_warning")
  w <- tryCatch(haploStack:::.hs_warn(wc), warning = identity)
  expect_identical(class(w), c("my_warning", "warning", "condition"))
})
