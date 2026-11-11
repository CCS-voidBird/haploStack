# =============================================================================
# test-090-gs-session.R — one genomicSimulation session per R process
# (design D1 §3.6, form row 4 a, build step D1-6)
#
#   S1  a second load while a session is open is an hs_error_state;
#       replace = TRUE closes the open one and invalidates its handles
#   S2  work_dir = NULL is a private directory, removed by gs_close();
#       a caller's work_dir is kept
#   S3  gs_close() acts only for the live session: a stale handle is a no-op,
#       a second close is a no-op, NULL / non-handles are no-ops
#   S4  the four bridge functions refuse a closed or replaced handle
#   S5  the driver leaves no session open and no private directory behind
#
# Everything below needs genomicSimulation (Suggests) and skips without it.
# =============================================================================

library(testthat)
library(haploStack)

.gss_meta <- function(seed = 7L) {
  m <- sim_geno_meta(n_samples = 8L, n_markers = 12L, n_chrom = 2L,
                     block_size = 2L, seed = seed)
  m$marker_effects <- sim_effect_df(dimnames(m$geno_arrays)[[2L]],
                                    alleles = c("A", "T"))
  m
}
.gss_clean <- function() withr::defer(haploStack:::.gs_session_clear(), envir = parent.frame())
.gss_dirs <- function() list.files(tempdir(), pattern = "^hs_gs_", full.names = TRUE)

test_that("S1: one session per process; replace = TRUE replaces and invalidates", {
  skip_if_not_installed("genomicSimulation")
  .gss_clean()
  m  <- .gss_meta()
  h1 <- gs_load_geno_meta(m, map_unit = "bp")
  expect_true(haploStack:::.gs_session_open())
  expect_type(h1$token, "integer")
  expect_true(h1$own_dir)
  expect_error(gs_load_geno_meta(m, map_unit = "bp"), class = "hs_error_state")
  expect_error(gs_load_geno_meta(m, map_unit = "bp"), "one session per R process")
  d1 <- h1$work_dir
  expect_true(dir.exists(d1))
  h2 <- gs_load_geno_meta(m, replace = TRUE, map_unit = "bp")
  expect_true(haploStack:::.gs_session_open())
  expect_identical(h2$token, h1$token + 1L)
  expect_false(dir.exists(d1))                      # old private dir gone
  plan <- make_cross_design(c("S1", "S2", "S3"), method = "all_pairs", n_offspring = 1L)
  expect_error(gs_cross_plan(h1, plan), class = "hs_error_state")   # stale
  cp <- gs_cross_plan(h2, plan)                                       # live
  expect_s3_class(cp$handle, "gs_handle")
  expect_identical(cp$handle$token, h2$token)
  gs_close(cp$handle)
  expect_false(haploStack:::.gs_session_open())
})

test_that("S2: private scratch is removed on close; a caller's work_dir is kept", {
  skip_if_not_installed("genomicSimulation")
  .gss_clean()
  m <- .gss_meta()
  h <- gs_load_geno_meta(m, map_unit = "bp")
  expect_match(basename(h$work_dir), "^hs_gs_")
  expect_true(file.exists(h$files[["geno"]]))
  gs_close(h)
  expect_false(dir.exists(h$work_dir))
  own <- withr::local_tempdir()
  h <- gs_load_geno_meta(m, work_dir = own, map_unit = "bp")
  expect_false(h$own_dir)
  gs_close(h)
  expect_true(dir.exists(own))
  expect_true(file.exists(file.path(own, "gs_geno.txt")))
})

test_that("S3: gs_close() is idempotent and ignores stale handles", {
  skip_if_not_installed("genomicSimulation")
  .gss_clean()
  m  <- .gss_meta()
  h1 <- gs_load_geno_meta(m, map_unit = "bp")
  h2 <- gs_load_geno_meta(m, replace = TRUE, map_unit = "bp")
  expect_null(gs_close(h1))                          # stale: no-op
  expect_true(haploStack:::.gs_session_open())       # h2 untouched
  plan <- make_cross_design(c("S1", "S2"), method = "all_pairs", n_offspring = 1L)
  expect_s3_class(gs_cross_plan(h2, plan, advance = FALSE)$handle, "gs_handle")
  expect_null(gs_close(h2))
  expect_false(haploStack:::.gs_session_open())
  expect_null(gs_close(h2))                          # twice: no-op
  expect_null(gs_close(NULL))
  expect_null(gs_close(list()))
})

test_that("S4: the bridge functions refuse a closed handle", {
  skip_if_not_installed("genomicSimulation")
  .gss_clean()
  m <- .gss_meta()
  h <- gs_load_geno_meta(m, map_unit = "bp")
  gs_close(h)
  plan <- make_cross_design(c("S1", "S2"), method = "all_pairs", n_offspring = 1L)
  expect_error(gs_cross_plan(h, plan), class = "hs_error_state")
  expect_error(gs_group_to_geno_meta(h), class = "hs_error_state")
  expect_error(gs_group_gebv(h), class = "hs_error_state")
  expect_error(gs_sync_group_names(h, "a", "b"), class = "hs_error_state")
  expect_error(gs_cross_plan(list(), plan), "must be a gs_handle")
})

test_that("S5: a driver run leaves no session open and no private directory", {
  skip_if_not_installed("genomicSimulation")
  .gss_clean()
  m <- .gss_meta(seed = 11L)
  before <- .gss_dirs()
  res <- run_stacking_simulation(
    m, replications = 2L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", select_top_blocks = 20L,
    cross_backend = "genomicSimulation", output_dir = NULL, verbose = FALSE,
    backend_opts = list(map_unit = "bp"))
  expect_true(all(is.finite(res$summary[cycle > 0L]$mean_F1mean)))
  expect_false(haploStack:::.gs_session_open())
  expect_setequal(.gss_dirs(), before)
  # and a new session can open right after
  h <- gs_load_geno_meta(m, map_unit = "bp")
  expect_true(haploStack:::.gs_session_open())
  gs_close(h)
})
