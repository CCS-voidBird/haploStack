# =============================================================================
# test-090-plan.R — the future plan has one owner (design D1 §3.7, step D1-4)
#
#   P1  .hs_parallel_ok(): n_cores <= 1 is FALSE; > 1 needs future + future.apply
#   P2  .hs_with_plan(): n_cores <= 1 leaves the plan alone; under a sequential
#       plan it runs multisession with n_cores workers and puts sequential back
#   P3  a plan the caller set is used as it is (one message) and left in place
#   P4  .hs_keep_plan() puts the caller's plan back whatever the code did
#   P5  the entry points go through it: no future::plan() call outside state.R
#       (the hygiene test) and n_cores = 1 never touches the plan
# =============================================================================

library(testthat)
library(haploStack)

# the specific strategy in a plan's class vector (it starts with "FutureStrategy")
.plan_kind <- function() {
  k <- class(future::plan())
  intersect(c("multisession", "multicore", "cluster", "sequential"), k)[1L]
}

test_that("P1: .hs_parallel_ok() is the one rule for 'can we go parallel'", {
  expect_false(haploStack:::.hs_parallel_ok(1L))
  expect_false(haploStack:::.hs_parallel_ok(0L))
  expect_false(haploStack:::.hs_parallel_ok(NA))
  skip_if_not_installed("future"); skip_if_not_installed("future.apply")
  expect_true(haploStack:::.hs_parallel_ok(2L))
  expect_true(haploStack:::.hs_parallel_ok("3"))
})

test_that("P2: .hs_with_plan() sets multisession for the call and restores sequential", {
  skip_if_not_installed("future"); skip_if_not_installed("future.apply")
  withr::local_options(list(future.rng.onMisuse = "ignore"))
  future::plan(future::sequential)
  withr::defer(future::plan(future::sequential))
  # n_cores = 1: untouched, the code runs inline
  seen <- NULL
  v <- haploStack:::.hs_with_plan(1L, { seen <- .plan_kind(); 42L })
  expect_identical(v, 42L)
  expect_identical(seen, "sequential")
  # n_cores = 2: multisession with 2 workers inside, sequential again after
  inside <- haploStack:::.hs_with_plan(2L, {
    list(cls = .plan_kind(), n = future::nbrOfWorkers(),
         res = unlist(future.apply::future_lapply(1:3, function(i) i * 2L)))
  })
  expect_identical(inside$cls, "multisession")
  expect_identical(as.integer(inside$n), 2L)
  expect_identical(inside$res, c(2L, 4L, 6L))
  expect_s3_class(future::plan(), "sequential")
  # restored on error too
  expect_error(haploStack:::.hs_with_plan(2L, stop("inside")), "inside")
  expect_s3_class(future::plan(), "sequential")
})

test_that("P3: a plan the caller set is used as it is, with one message", {
  skip_if_not_installed("future"); skip_if_not_installed("future.apply")
  withr::local_options(list(future.rng.onMisuse = "ignore"))
  future::plan(future::multisession, workers = 2L)
  withr::defer(future::plan(future::sequential))
  expect_message(
    inside <- haploStack:::.hs_with_plan(3L, list(cls = .plan_kind(),
                                                  n = future::nbrOfWorkers())),
    "using the caller's future plan")
  expect_identical(inside$cls, "multisession")
  expect_identical(as.integer(inside$n), 2L)          # not 3: the caller's count
  expect_s3_class(future::plan(), "multisession")       # left in place
})

test_that("P4: .hs_keep_plan() puts the caller's plan back", {
  skip_if_not_installed("future")
  future::plan(future::sequential)
  withr::defer(future::plan(future::sequential))
  v <- haploStack:::.hs_keep_plan({ future::plan(future::multisession, workers = 2L); "x" })
  expect_identical(v, "x")
  expect_s3_class(future::plan(), "sequential")
})

test_that("P5: n_cores = 1 runs never touch the plan; the entry points share the helper", {
  skip_if_not_installed("future")
  future::plan(future::sequential)
  withr::defer(future::plan(future::sequential))
  m <- sim_geno_meta(n_samples = 10L, n_markers = 18L, n_chrom = 2L,
                     block_size = 3L, seed = 4L)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("A", "T"))
  before <- future::plan()
  invisible(run_haplotype_catalogue(m, list(A = e), figures = FALSE, n_cores = 1L))
  invisible(run_stacking_simulation(m, e, replications = 1L, select_cycles = 1L,
                                    top_founders = 4L, num_offspring = 2L, scheme = "TS",
                                    select_top_blocks = 40L, output_dir = NULL,
                                    verbose = FALSE, n_cores = 1L))
  expect_identical(class(future::plan()), class(before))
  src <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(src), "R/ not reachable")
  # the catalogue dispatches in place; the simulation's dispatcher is in
  # stacking_step.R (D4-3), the wrapper decides parallel or not
  for (f in c("catalogue_main.R", "stacking_step.R")) {
    code <- readLines(file.path(src, f), warn = FALSE)
    expect_true(any(grepl(".hs_parallel_ok(", code, fixed = TRUE)), info = f)
    expect_true(any(grepl(".hs_with_plan(", code, fixed = TRUE)), info = f)
  }
  code <- readLines(file.path(src, "stacking_main.R"), warn = FALSE)
  expect_true(any(grepl(".hs_parallel_ok(", code, fixed = TRUE)))
  expect_false(any(grepl("future::plan(", code, fixed = TRUE)))
})
