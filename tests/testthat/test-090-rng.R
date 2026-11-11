# =============================================================================
# test-090-rng.R — local seeds (design D1-2, framework F2)
#
# Every seeded call gives the output a bare set.seed(seed) gave up to 0.8.8
# (values pinned below from the 0.8.8 snapshot) AND leaves the caller's random
# stream exactly as it was. "Untouched" is checked two ways: a stream that has
# drawn already continues identically, and a session that has not drawn yet
# stays that way (.Random.seed absent).
#
#   N1  .with_seed(): value, restore, NULL passthrough, restore on error
#   N2  .local_seed() and .with_rng_state() round trips
#   N3  sim_geno_matrix / sim_effect_df / sim_geno_meta pinned to 0.8.8
#   N4  sim_dataset: pinned, and $map equals the container's map (changed)
#   N5  the encoding detector restores the caller's stream
#   N6  tune_gbs_penalty() and the three GA selectors restore the stream
#   N7  fit_genomic_prediction(seed =) restores the stream (BGLR, skip-gated)
# =============================================================================

.rng_absent <- function() !exists(".Random.seed", envir = globalenv(), inherits = FALSE)
.rng_clear  <- function() {
  if (!.rng_absent()) rm(list = ".Random.seed", envir = globalenv())
}
# Next draws of the untouched stream, without disturbing it for the caller.
.rng_peek <- function(n = 3L) {
  s <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", s, envir = globalenv()))
  stats::runif(n)
}
# Run `expr` and check the caller's stream continues as if `expr` had not run.
expect_stream_untouched <- function(expr) {
  set.seed(2026)
  want <- .rng_peek()
  force(expr)
  expect_identical(stats::runif(3L), want)
  .rng_clear()
  force(expr)
  expect_true(.rng_absent())
}

test_that("N1: .with_seed() reproduces set.seed() and restores the caller's stream", {
  set.seed(5L); want <- stats::runif(4L)
  expect_identical(haploStack:::.with_seed(5L, stats::runif(4L)), want)
  expect_stream_untouched(haploStack:::.with_seed(5L, stats::runif(4L)))

  # NULL evaluates on the caller's stream
  set.seed(9L); want <- stats::runif(2L)
  set.seed(9L)
  expect_identical(haploStack:::.with_seed(NULL, stats::runif(2L)), want)

  # restored even when the code throws
  set.seed(2026); want <- .rng_peek()
  expect_error(haploStack:::.with_seed(1L, stop("boom")), "boom")
  expect_identical(stats::runif(3L), want)
})

test_that("N2: .local_seed() and .with_rng_state() round trips", {
  f <- function(seed) { haploStack:::.local_seed(seed); stats::runif(3L) }
  set.seed(7L); want <- stats::runif(3L)
  expect_identical(f(7L), want)
  expect_stream_untouched(f(7L))
  set.seed(3L); want <- stats::runif(3L); set.seed(3L)
  expect_identical(f(NULL), want)          # NULL: caller's stream, no restore

  # a job resumed from its captured state continues the same stream
  set.seed(1L); r1 <- stats::runif(2L); r2 <- stats::runif(2L)
  s1 <- haploStack:::.with_rng_state(1L, stats::runif(2L))
  s2 <- haploStack:::.with_rng_state(s1$state, stats::runif(2L))
  expect_identical(s1$value, r1)
  expect_identical(s2$value, r2)
  expect_true(is.integer(s1$state) && length(s1$state) > 1L)
  expect_stream_untouched(haploStack:::.with_rng_state(1L, stats::runif(2L)))
})

test_that("N3: sim_geno_matrix / sim_effect_df / sim_geno_meta pinned to 0.8.8", {
  g <- sim_geno_matrix(4, 6, "nucleotide", seed = 11L)
  expect_identical(as.vector(g),
    c("T/T", "T/T", "A/A", "A/T", "A/T", "T/A", "T/A", "A/T", "A/A",
      "A/A", "T/A", "T/A", "T/A", "T/T", "A/A", "T/T", "T/A", "T/T",
      "A/T", "A/T", "T/A", "A/A", "T/T", "A/T"))
  expect_identical(dimnames(g), list(paste0("Mkr", 1:6), paste0("S", 1:4)))

  e <- sim_effect_df(paste0("Mkr", 1:3), seed = 13L)
  expect_identical(round(e$Effect, 6),
                   c(0.277163, 0.09366, -0.140136, 0.571263, 0.887582, 0.207763))

  # the same draws as the 0.8.8 "selectiontools" default, 1 -> A, 4 -> T
  # (the container stores nucleotides only since 0.9.0, C-5)
  m <- sim_geno_meta(5, 8, "nucleotide", n_chrom = 2L, block_size = 4L, seed = 12L)
  expect_identical(m$map_info$Pos,
    c(19300000L, 56000000L, 70100000L, 95800000L, 19600000L, 26200000L,
      56100000L, 80500000L))
  expect_identical(as.vector(m$geno_arrays[1, , ]),
    c("T", "T", "A", "A", "A", "T", "T", "A", "T", "A", "T", "T",
      "T", "A", "A", "T"))
  expect_true(is_geno_meta(m))

  expect_stream_untouched(sim_geno_matrix(4, 6, "nucleotide", seed = 11L))
  expect_stream_untouched(sim_effect_df(paste0("Mkr", 1:3), seed = 13L))
  expect_stream_untouched(sim_geno_meta(5, 8, block_size = 4L, seed = 12L))
  expect_stream_untouched(sim_map_data(paste0("Mkr", 1:6), seed = 4L))

  # a seeded map is reproducible on its own; the unseeded map draws from the
  # caller's stream, as before
  expect_identical(sim_map_data(paste0("Mkr", 1:6), seed = 4L),
                   sim_map_data(paste0("Mkr", 1:6), seed = 4L))
  set.seed(4L); want <- sim_map_data(paste0("Mkr", 1:6))
  expect_identical(sim_map_data(paste0("Mkr", 1:6), seed = 4L), want)

  expect_error(sim_geno_matrix(2, 2, encoding = "dosage"), "'arg' should be one of")
})

test_that("N4: sim_dataset() pinned to 0.8.8, and $map is the container's map", {
  ds <- sim_dataset(4, 6, block_size = 3L, seed = 14L)
  expect_identical(as.vector(ds$geno_st),
    c("1/1", "4/4", "1/1", "1/4", "1/4", "1/4", "1/4", "4/1", "1/4",
      "1/4", "4/4", "4/4", "1/4", "1/1", "1/1", "4/1", "4/4", "4/4",
      "4/4", "1/4", "1/4", "4/4", "4/1", "1/4"))
  expect_identical(as.vector(ds$geno_nuc),
    c("A/A", "A/A", "T/A", "A/A", "A/A", "A/A", "T/T", "T/A", "A/A",
      "A/A", "T/A", "T/T", "A/T", "T/T", "T/A", "T/T", "T/A", "T/A",
      "A/A", "A/T", "T/T", "T/T", "A/A", "T/T"))
  expect_identical(round(ds$eff_st$Effect, 6),
    c(-0.330925, -0.03244, 0.859477, 0.534497, 1.060833, -0.188483,
      0.748577, 0.521592, -0.01807, -0.191411, 0.615973, 0.149711))
  # the container is unchanged from 0.8.8 ...
  expect_identical(ds$geno_meta$map_info$Pos,
    c(3100000L, 71800000L, 20000000L, 60200000L, 25000000L, 57800000L))
  # ... and $map now carries the same positions (0.8.8 drew it separately:
  # 1900000, 72200000, 69400000, 73000000, 19800000, 29800000)
  expect_identical(ds$map$Pos, ds$geno_meta$map_info$Pos)
  expect_identical(ds$map$Marker, ds$geno_meta$map_info$Marker)
  expect_identical(names(ds$map), c("Marker", "Chrom", "Pos"))
  expect_identical(ds$block_map, ds$geno_meta$block_marker_map_info[, c("Chrom", "Marker", "Block")])
  expect_stream_untouched(sim_dataset(4, 6, block_size = 3L, seed = 14L))
})

test_that("N5: the encoding detector keeps its verdicts and restores the stream", {
  v <- c(rep("A/T", 60), rep("G/G", 50), rep("0/0", 10))
  expect_identical(detect_genotype_encoding(v), "nucleotide")
  expect_identical(haploStack:::.detect_encoding(c(rep("1/4", 30), rep("4/4", 40))),
                   "selectiontools")
  expect_stream_untouched(detect_genotype_encoding(v))
  # a read/build goes through the detector; up to 0.8.8 it left seed 42 behind
  arr <- sim_geno_matrix(6, 8, "nucleotide", seed = 2L)
  expect_stream_untouched(suppressMessages(split_haplotypes(arr)))
})

test_that("N6: tune_gbs_penalty() and the GA selectors restore the caller's stream", {
  # a wide local-GEBV table (Block + <sample>.<haplotype> columns), as the
  # selectors take it; rounding makes ties so the GAs have real choices
  cols <- as.vector(t(outer(sprintf("S%02d", 1:8), 1:2, paste, sep = ".")))
  v    <- haploStack:::.with_seed(21L, matrix(round(stats::rnorm(12 * length(cols)), 1), 12))
  lg   <- data.table::data.table(Block = sprintf("b%02d", 1:12))
  lg[, (cols) := as.data.frame(v)]
  ga   <- list(pop_size = 10L, n_generations = 5L)

  tg <- function() suppressWarnings(suppressMessages(
    tune_gbs_penalty(lg, n_parents = 4L, grid = c(1, 2), n_rep = 2L,
                     subsample_frac = 0.8, seed = 3L)))
  r1 <- tg(); r2 <- tg()
  expect_identical(r1, r2)
  expect_stream_untouched(tg())

  g1 <- function() suppressWarnings(suppressMessages(
    suppressWarnings(classes = "hs_warning_deprecated",      # row 27 a
    genetic_algorithm_selection(lg, n_parents = 4L, pop_size = 10L,
                                n_generations = 5L, seed = 4L))))
  expect_identical(g1(), g1())
  expect_stream_untouched(g1())

  o1 <- function() suppressWarnings(suppressMessages(
    ohs_selection(lg, n_parents = 4L, pop_size = 10L, n_generations = 5L, seed = 5L)))
  expect_identical(o1(), o1())
  expect_stream_untouched(o1())

  b1 <- function() suppressWarnings(suppressMessages(
    gbs_selection(lg, n_parents = 4L, gbs_objective = "set", gbs_engine = "GA",
                  ga_control = ga, seed = 6L)))
  expect_identical(b1(), b1())
  expect_stream_untouched(b1())
})

test_that("N7: fit_genomic_prediction(seed =) is reproducible and restores the stream", {
  skip_if_not_installed("rrBLUP")
  meta <- sim_geno_meta(20, 30, "nucleotide", block_size = 5L, seed = 31L)
  ph <- data.frame(Sample = dimnames(meta$geno_arrays)[[1L]],
                   Yield = haploStack:::.with_seed(1L, stats::rnorm(20)))
  f <- function(...) suppressWarnings(suppressMessages(
    fit_genomic_prediction(meta, ph, model = "rrblup", ...)))
  expect_identical(f()$gebv, f(seed = 2L)$gebv)      # deterministic model
  expect_stream_untouched(f(seed = 2L))

  skip_if_not_installed("BGLR")
  td <- withr::local_tempdir()
  b <- function(seed) suppressWarnings(suppressMessages(
    fit_genomic_prediction(meta, ph, model = "bayesc", nIter = 300L, burnIn = 100L,
                           seed = seed, saveAt = file.path(td, "b_"))))
  expect_identical(b(8L)$gebv, b(8L)$gebv)
  expect_false(identical(b(8L)$gebv, b(9L)$gebv))
  expect_stream_untouched(b(8L))
})
