# =============================================================================
# test-stacking-gs-session.R
#
# Verification for the persistent genomicSimulation session backend (PRD
# dev/PRD_gs_persistent_session.md, issues 02/03): the founder trio is written
# ONCE, output shape matches the block backend, and the new path is deterministic
# within itself (same seed -> same result). All gS-dependent, so skipped unless
# genomicSimulation is installed.
# =============================================================================

library(testthat)
library(haploStack)

.gss_inputs <- function(n_samples = 12L, n_markers = 16L, seed = 7L) {
  m <- sim_geno_meta(n_samples = n_samples, n_markers = n_markers,
                     n_chrom = 2L, block_size = 2L, seed = seed)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("A", "T"))
  list(m = m, e = e, n_samples = n_samples, n_markers = n_markers)
}

test_that("persistent gS session writes the founder trio once (no per-cycle write)", {
  skip_if_not_installed("genomicSimulation")
  io <- .gss_inputs()

  # 0.9.0: each job's scratch is a private hs_gs_* directory under tempdir(),
  # removed when the job closes unless backend_opts$keep_files keeps it.
  before <- list.files(tempdir(), pattern = "^hs_gs_", full.names = TRUE)
  run_stacking_simulation(
    io$m, io$e, replications = 1L, select_cycles = 4L, top_founders = 4L,
    num_offspring = 3L, scheme = "TS", select_top_blocks = 50L,
    cross_backend = "genomicSimulation", output_dir = NULL, verbose = FALSE,
    backend_opts = list(keep_files = TRUE, map_unit = "bp"))
  kept <- setdiff(list.files(tempdir(), pattern = "^hs_gs_", full.names = TRUE), before)
  withr::defer(unlink(kept, recursive = TRUE))
  expect_length(kept, 1L)
  geno <- file.path(kept, "gs_geno.txt")
  expect_true(file.exists(geno))

  # The genotype file must still hold the FOUNDER matrix (written once at load):
  # markers x founders. A per-cycle rewrite would instead leave the last
  # generation (choose(4,2)*3 = 18 offspring != 12 founders).
  tab <- utils::read.delim(geno, row.names = 1L, check.names = FALSE)
  expect_equal(ncol(tab), io$n_samples)     # founders only -> written once
  expect_equal(nrow(tab), io$n_markers)

  # map / eff are constant across cycles too -> present and written once.
  expect_true(file.exists(file.path(kept, "gs_map.txt")))
  expect_true(file.exists(file.path(kept, "gs_eff.txt")))

  # without keep_files nothing is left behind
  before <- list.files(tempdir(), pattern = "^hs_gs_", full.names = TRUE)
  run_stacking_simulation(
    io$m, io$e, replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", select_top_blocks = 50L,
    cross_backend = "genomicSimulation", output_dir = NULL, verbose = FALSE,
    backend_opts = list(map_unit = "bp"))
  expect_setequal(list.files(tempdir(), pattern = "^hs_gs_", full.names = TRUE), before)
})

test_that("persistent gS backend matches the block backend's output shape", {
  skip_if_not_installed("genomicSimulation")
  io <- .gss_inputs(seed = 11L)
  run <- function(backend) run_stacking_simulation(
    io$m, io$e, replications = 1L, select_cycles = 3L, top_founders = 4L,
    num_offspring = 2L, scheme = "GBS", select_top_blocks = 50L,
    metrics = c("d2ideal", "coancestry", "gene_diversity"),
    cross_backend = backend, output_dir = NULL, verbose = FALSE, seed = 3L,
    backend_opts = list(map_unit = "bp"))

  rb <- run("block")
  rg <- run("genomicSimulation")

  expect_equal(names(rg$summary), names(rb$summary))
  expect_equal(names(rg$plans),   names(rb$plans))
  expect_equal(nrow(rg$summary),  nrow(rb$summary))
  expect_equal(nrow(rg$plans),    nrow(rb$plans))
  # both produce finite F1 means; values themselves are a valid-but-different stream
  expect_true(all(is.finite(rg$summary[cycle > 0L]$mean_F1mean)))
})

test_that("persistent gS path is deterministic within itself (selfing + interlude)", {
  skip_if_not_installed("genomicSimulation")
  io <- .gss_inputs(seed = 13L)
  run <- function() run_stacking_simulation(
    io$m, io$e, replications = 1L, select_cycles = 3L, top_founders = 4L,
    num_offspring = 3L, scheme = "GBS", selfing_cycles = 1L, ga_start_by = 2L,
    select_top_blocks = 50L, metrics = "gene_diversity",
    cross_backend = "genomicSimulation", output_dir = NULL, verbose = FALSE,
    seed = 42L, backend_opts = list(map_unit = "bp"))

  r1 <- run(); r2 <- run()
  # same seed -> identical summary AND plans across full runs (decision cycles,
  # the ga_start_by TS->GBS interlude, and the selfing sub-generations).
  expect_equal(r1$summary, r2$summary)
  expect_equal(r1$plans,   r2$plans)
  expect_true(any(r1$plans$select_type == "selfing"))   # selfing branch exercised
})
