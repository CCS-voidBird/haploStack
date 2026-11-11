# =============================================================================
# test-090-schema.R — container schema, upgrade and the RDS pre-read check
# (designs D2-4 and D1-9)
#
#   C6a new containers are schema 2 and stamped with the package version
#   C6b a 0.8.8 fixture (schema 1) upgrades: counters, char_dict list, block
#       provenance "unknown"; load_geno_meta() upgrades in memory with one
#       message and leaves the file unchanged
#   C6c schema 3 is refused (hs_error_schema) by upgrade and load; is_geno_meta
#       returns FALSE on it
#   C6d a pre-0.6.0 fixture keeps its "0" tokens after load (no silent repair)
#       and strict_missing = TRUE makes the load an error
#   C6e upgrade refuses a container two blocks share a marker in
#       (hs_error_block_map) and a stale annotation (hs_error_schema)
#   R1  load_geno_meta(): a text file, a size over max_bytes and a file that
#       is not a container are refused with the right classes; the pre-read
#       check accepts every RDS compression
#   G1  .hs_check_meta(): reason and hint in the message, hint as a field,
#       the caller's call on the condition; the token notice only on an
#       unstamped container
# =============================================================================

.fx090 <- function(name) {
  p <- testthat::test_path("fixtures", name)
  skip_if_not(file.exists(p), paste("fixture", name, "missing"))
  p
}

test_that("C6a: new containers carry schema 2 and the building version", {
  m <- sim_geno_meta(6, 12, block_size = 4, seed = 1)
  expect_identical(m$metadata$schema_version, 2L)
  expect_identical(m$metadata$built_with, as.character(utils::packageVersion("haploStack")))
  expect_identical(haploStack:::.hs_schema_of(m), 2L)
  expect_identical(haploStack:::.HS_SCHEMA_VERSION, 2L)
  # a transform keeps the stamp; an upgrade of a current container is a no-op
  f <- suppressWarnings(filter_geno_meta(m, keep_samples = c("S1", "S2")))
  expect_identical(f$metadata$schema_version, 2L)
  expect_identical(upgrade_geno_meta(m), m)
})

test_that("C6b: the 0.8.8 fixture is schema 1 and still a valid container", {
  p <- .fx090("genometa_v1_0.8.8.rds")
  v1 <- readRDS(p)
  expect_null(v1$metadata$schema_version)
  expect_identical(haploStack:::.hs_schema_of(v1), 1L)
  expect_false(is.list(v1$metadata$char_dict))         # 0.8.8 wrote a vector
  expect_true(is_geno_meta(v1, verbose = FALSE))        # schema 1 is readable
})

test_that("C6b: upgrade_geno_meta() on the 0.8.8 fixture", {
  p <- .fx090("genometa_v1_0.8.8.rds")
  v1 <- readRDS(p)
  v1$metadata$n_markers <- 999L                       # a stale counter, as a hand edit leaves
  expect_message(upgrade_geno_meta(v1), "schema 1 -> 2", class = "hs_message")
  expect_message(upgrade_geno_meta(v1), "n_markers re-derived", class = "hs_message")
  up <- suppressMessages(upgrade_geno_meta(v1))
  expect_identical(up$metadata$schema_version, 2L)
  expect_identical(up$metadata$n_markers, 12L)
  expect_true(is.list(up$metadata$char_dict))
  expect_identical(up$metadata$blocks$source, "unknown")
  expect_identical(up$metadata$blocks$n_blocks,
                   length(unique(v1$block_marker_map_info$Block)))
  # the 0.8.8 fixture holds SelectionTools codes: the one change to the array
  # is the recode to nucleotides (0.9.0, C-5), and the effects follow
  expect_message(upgrade_geno_meta(v1), "recoded to nucleotides", class = "hs_message")
  expect_identical(up$metadata$encoding, "nucleotide")
  expect_identical(up$geno_arrays, haploStack:::.st_to_nucleotide(v1$geno_arrays))
  expect_identical(up$phenotype, v1$phenotype)
  want_eff <- lapply(v1$marker_effects, function(d) {
    d$Allele <- haploStack:::.st_to_nucleotide(d$Allele); d })
  expect_identical(up$marker_effects, want_eff)
  expect_true(is_geno_meta(up, strict = TRUE))
  # the loader does the same, says so once, and leaves the file as it was
  md5 <- tools::md5sum(p)
  msgs <- character()
  withCallingHandlers(l <- load_geno_meta(p),
                      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_length(grep("saved at container schema 1", msgs), 1L)
  expect_identical(l$metadata$schema_version, 2L)
  expect_identical(tools::md5sum(p), md5)
  expect_identical(l$geno_arrays, haploStack:::.st_to_nucleotide(v1$geno_arrays))
})

test_that("C6c: a newer schema is refused", {
  m <- sim_geno_meta(4, 8, block_size = 2, seed = 3)
  m$metadata$schema_version <- 3L
  m$metadata$built_with <- "9.9.9"
  expect_false(is_geno_meta(m, verbose = FALSE))
  expect_message(is_geno_meta(m), "newer than this haploStack", class = "hs_message")
  e <- tryCatch(upgrade_geno_meta(m), error = identity)
  expect_s3_class(e, "hs_error_schema")
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "saved by haploStack 9.9.9")
  tf <- tempfile(fileext = ".rds"); saveRDS(m, tf)
  expect_error(load_geno_meta(tf), class = "hs_error_schema")
  # not a container at all
  saveRDS(list(a = 1), tf)
  expect_error(load_geno_meta(tf), "not a GenoMeta container", class = "hs_error_schema")
})

test_that("C6d: a pre-0.6.0 container keeps its tokens after load; strict_missing errors", {
  withr::local_options(list(haploStack.strict_missing = NULL))
  p <- .fx090("genometa_pre060.rds")
  msgs <- character()
  withCallingHandlers(l <- load_geno_meta(p),
                      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(l$geno_arrays == "0", na.rm = TRUE))          # no silent repair
  expect_null(l$metadata$missing_tokens)                         # still unstamped
  expect_length(grep("predates the 0.6.0 invariant", msgs), 1L)
  expect_identical(l$metadata$schema_version, 2L)
  expect_error(load_geno_meta(p, strict_missing = TRUE), "missing token",
               class = "hs_error_input")
  r <- suppressMessages(repair_geno_meta(l))
  expect_false(any(r$geno_arrays == "0", na.rm = TRUE))
  expect_true(is_geno_meta(r, strict = TRUE))
})

test_that("C6e: upgrade refuses an overlap and a stale annotation", {
  m <- sim_geno_meta(4, 8, block_size = 2, seed = 4)
  m$metadata$schema_version <- NULL
  dup <- m; dup$block_marker_map_info <- rbind(dup$block_marker_map_info,
                                               dup$block_marker_map_info[1, ])
  e <- tryCatch(upgrade_geno_meta(dup, verbose = FALSE), error = identity)
  expect_s3_class(e, "hs_error_block_map")
  expect_identical(e$bad, as.character(m$block_marker_map_info$Marker[1]))
  expect_false(is_geno_meta(dup, verbose = FALSE))
  expect_message(is_geno_meta(dup), "more than one block", class = "hs_message")
  bad <- m
  bad$annotation <- list(blocks = data.frame(Block = "NOPE"),
                         features = data.frame(Block = "NOPE"),
                         marker_gene = data.frame(Marker = character(0), Block = character(0)))
  e <- tryCatch(upgrade_geno_meta(bad, verbose = FALSE), error = identity)
  expect_s3_class(e, "hs_error_schema")
  expect_match(conditionMessage(e), "compute_block_annotation")
})

test_that("R1: load_geno_meta() checks the file before unserialising it", {
  tf <- tempfile(fileext = ".rds")
  writeLines(c("Sample\tT1", "S1\t1"), tf)
  e <- tryCatch(load_geno_meta(tf), error = identity)
  expect_s3_class(e, "hs_error_format")
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "does not start like an R serialisation")
  m <- sim_geno_meta(4, 8, block_size = 2, seed = 5)
  for (comp in list(TRUE, "bzip2", "xz", FALSE)) {
    saveRDS(m, tf, compress = comp)
    expect_true(haploStack:::.hs_check_rds_file(tf))
    expect_identical(suppressMessages(load_geno_meta(tf))$geno_arrays, m$geno_arrays)
  }
  saveRDS(m, tf, ascii = TRUE)
  expect_true(haploStack:::.hs_check_rds_file(tf))
  saveRDS(m, tf)
  e <- tryCatch(load_geno_meta(tf, max_bytes = 100), error = identity)
  expect_s3_class(e, "hs_error_format")
  expect_match(conditionMessage(e), "over max_bytes")
  expect_error(load_geno_meta(tempfile()), "File not found")
  # the help page says so
  rd <- testthat::test_path("..", "..", "man", "load_geno_meta.Rd")
  skip_if_not(file.exists(rd), "man/ not reachable")
  expect_true(any(grepl("Trusted input only", readLines(rd, warn = FALSE), fixed = TRUE)))
})

test_that("G1: the entry-point guard names the caller, the reason and the hint", {
  e <- tryCatch(compute_grm(list(a = 1)), error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e),
               "^compute_grm: 'geno_meta' is not a valid GenoMeta: missing required field\\(s\\)")
  expect_match(conditionMessage(e), "filter_geno_meta\\(\\) to subset markers")
  expect_identical(e$hint, "Use filter_geno_meta() to subset markers before calling if needed.")
  expect_identical(as.character(conditionCall(e)[[1L]]), "compute_grm")
  e <- tryCatch(compute_marker_stats(42), error = identity)
  expect_match(conditionMessage(e), "^compute_marker_stats: 'geno_meta' is not a valid GenoMeta: 'x' is not a list\\.$")
  expect_null(e$hint)
  # the same core answers is_geno_meta()
  m <- sim_geno_meta(4, 8, block_size = 2, seed = 6)
  m$metadata$n_samples <- 3L
  expect_identical(haploStack:::.hs_meta_reason(m),
                   "metadata$n_samples (3) does not equal dim(geno_arrays)[1] (4).")
  e <- tryCatch(compute_marker_stats(m), error = identity)
  expect_match(conditionMessage(e), "metadata\\$n_samples \\(3\\)")
  # the token notice: raised on an unstamped container only
  u <- sim_geno_meta(4, 8, block_size = 2, seed = 6)
  u$geno_arrays["S1", 1L, ] <- "0"
  expect_silent(haploStack:::.hs_check_meta(u, "f"))            # stamped: trusted
  u$metadata$missing_tokens <- NULL
  expect_warning(haploStack:::.hs_check_meta(u, "f"), "missing token")
})
