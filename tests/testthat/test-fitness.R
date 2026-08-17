# =============================================================================
# test-fitness.R
#
# QA tests for the fitness module:
#   create_fitness_index(), update_fitness_index(), write_fitness_file(),
#   run_fitness_pipeline()
# =============================================================================

library(testthat)
library(haploStack)

# Shared test data
.eff1 <- data.frame(
  Marker = c("M1","M1","M2","M2"),
  Allele = c("A","T","A","T"),
  Effect = c(0.5, -0.5, 0.2, -0.2),
  stringsAsFactors = FALSE
)
.eff2 <- data.frame(
  Marker = c("M1","M1","M2","M2"),
  Allele = c("A","T","A","T"),
  Effect = c(0.1, -0.1, 0.8, -0.8),
  stringsAsFactors = FALSE
)
.weights <- c(Yield = 0.7, Protein = 0.3)
.eff_list <- list(Yield = .eff1, Protein = .eff2)

# ── create_fitness_index ──────────────────────────────────────────────────────

test_that("create_fitness_index returns a data.table with required columns", {
  fi <- create_fitness_index(.eff_list, .weights)
  expect_true(is.data.table(fi))
  expect_true(all(c("Marker","Allele","Effect") %in% names(fi)))
})

test_that("create_fitness_index has correct number of rows", {
  fi <- create_fitness_index(.eff_list, .weights)
  # M1: 2 alleles (A,T), M2: 2 alleles → 4 rows total
  expect_equal(nrow(fi), 4L)
})

test_that("create_fitness_index computes weighted effect correctly", {
  fi <- create_fitness_index(.eff_list, .weights,
                              normalise_weights = FALSE)
  # M1/A: 0.7 * 0.5 + 0.3 * 0.1 = 0.35 + 0.03 = 0.38
  m1a <- fi[Marker == "M1" & Allele == "A", Effect]
  expect_equal(m1a, 0.7 * 0.5 + 0.3 * 0.1, tolerance = 1e-9)
})

test_that("create_fitness_index uses the weights AS GIVEN by default (0.5.4)", {
  # Raw by default: the composite is the plain sum(weight * effect), so doubling
  # every weight doubles the composite instead of cancelling out.
  fi_1 <- create_fitness_index(.eff_list, c(Yield = 1, Protein = 1))
  fi_2 <- create_fitness_index(.eff_list, c(Yield = 2, Protein = 2))
  expect_equal(fi_2$Effect, 2 * fi_1$Effect, tolerance = 1e-9)

  # M1/A: 1 * 0.5 + 1 * 0.1 = 0.6, not 0.3
  expect_equal(fi_1[Marker == "M1" & Allele == "A", Effect], 0.6,
               tolerance = 1e-9)
})

test_that("normalise_weights = TRUE rescales by sum(abs) without changing ranking", {
  fi_norm <- create_fitness_index(.eff_list, c(Yield = 2, Protein = 2),
                                  normalise_weights = TRUE)
  fi_base <- create_fitness_index(.eff_list, c(Yield = 1, Protein = 1),
                                  normalise_weights = TRUE)
  # After normalisation (50/50), both weight vectors give the same result.
  expect_equal(fi_norm$Effect, fi_base$Effect, tolerance = 1e-9)
  # It is a positive rescale of the raw composite, so the order is untouched.
  raw <- create_fitness_index(.eff_list, c(Yield = 2, Protein = 2))
  expect_identical(order(fi_norm$Effect), order(raw$Effect))
})

test_that("create_fitness_index normalises SIGNED weights by sum(abs) and flips direction", {
  # A negative weight makes a 'lower is better' trait subtract from the composite.
  fi <- create_fitness_index(.eff_list, c(Yield = 1, Protein = -1),
                             normalise_weights = TRUE)          # sum(abs) = 2
  # M1/A: (1 * 0.5 + (-1) * 0.1) / 2 = 0.4 / 2 = 0.20
  m1a <- fi[Marker == "M1" & Allele == "A", Effect]
  expect_equal(m1a, (1 * 0.5 + (-1) * 0.1) / 2, tolerance = 1e-9)

  # The sign flip does not depend on normalisation: raw keeps the subtraction.
  expect_equal(
    create_fitness_index(.eff_list,
                         c(Yield = 1, Protein = -1))[Marker == "M1" &
                                                       Allele == "A", Effect],
    1 * 0.5 + (-1) * 0.1, tolerance = 1e-9)

  # Weights that sum to ~0 must NOT blow up (guarded by sum(abs) not sum()).
  expect_true(all(is.finite(
    create_fitness_index(.eff_list, c(Yield = 1, Protein = -1),
                         normalise_weights = TRUE)$Effect)))
  expect_error(create_fitness_index(.eff_list, c(Yield = 0, Protein = 0),
                                    normalise_weights = TRUE),
               "all zero")
})

test_that("create_fitness_index errors when effect_list is empty", {
  expect_error(create_fitness_index(list(), c()), "at least one trait")
})

test_that("create_fitness_index errors when weights do not cover all traits", {
  expect_error(
    create_fitness_index(.eff_list, c(Yield = 1)),
    "must have names"
  )
})

test_that("create_fitness_index fills missing Marker×Allele with fill_missing", {
  # eff2 has no missing, but use different allele_dict to trigger fill
  fi <- create_fitness_index(list(Yield = .eff1), c(Yield = 1),
                              allele_dict = c("A","T","C"),
                              normalise_weights = FALSE)
  # Should include C allele rows with Effect = fill_missing (0 by default)
  c_rows <- fi[Allele == "C"]
  expect_true(nrow(c_rows) > 0L)
  expect_true(all(c_rows$Effect == 0))
})

# ── update_fitness_index ──────────────────────────────────────────────────────

test_that("update_fitness_index returns same structure as original", {
  fi      <- create_fitness_index(.eff_list, .weights)
  fi_new  <- update_fitness_index(fi, .eff_list, .weights)
  expect_equal(names(fi_new), names(fi))
  expect_equal(nrow(fi_new), nrow(fi))
})

test_that("update_fitness_index errors on non-data.table input", {
  expect_error(update_fitness_index(data.frame(), .eff_list, .weights),
               "data.table")
})

# ── write_fitness_file ────────────────────────────────────────────────────────

test_that("write_fitness_file writes a readable file with correct columns", {
  fi  <- create_fitness_index(.eff_list, .weights)
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp))
  write_fitness_file(fi, tmp)
  expect_true(file.exists(tmp))
  read_back <- utils::read.table(tmp, header = TRUE, stringsAsFactors = FALSE)
  expect_true(all(c("Marker","Allele","Effect") %in% names(read_back)))
})

test_that("write_fitness_file number of rows matches fitness index", {
  fi  <- create_fitness_index(.eff_list, .weights)
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp))
  write_fitness_file(fi, tmp)
  read_back <- utils::read.table(tmp, header = TRUE)
  expect_equal(nrow(read_back), nrow(fi))
})

# ── run_fitness_pipeline ──────────────────────────────────────────────────────

test_that("run_fitness_pipeline returns a fitness data.table invisibly", {
  fi <- run_fitness_pipeline(.eff_list, .weights)
  expect_true(is.data.table(fi))
  expect_true("Effect" %in% names(fi))
})

test_that("run_fitness_pipeline writes output_file when specified", {
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp))
  run_fitness_pipeline(.eff_list, .weights, output_file = tmp)
  expect_true(file.exists(tmp))
})
