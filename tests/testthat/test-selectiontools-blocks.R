# =============================================================================
# test-selectiontools-blocks.R  (0.7.9; names from 0.7.10)
#
# The SelectionTools resolver: the GenoMeta route st_ld_blocks() (0.7.9's
# define_haplotype_blocks_meta()) and the file route st_ld_blocks_files()
# (define_haplotype_blocks() up to 0.7.9). Toy panels: helper-ld-toy.R.
#   B1  planted LD groups come back as exactly those blocks
#   B2  attach = TRUE keeps effects / phenotype / FAM, and chromosome NAMES
#   B3  tolerance absorbs an LD-breaking marker inside a block
#   B4  unplaced markers are left out of blocking and become singletons
#   B5  the blocking rule is validated before SelectionTools is touched
#   B6  no file lands in the working directory; the private dir is removed
#       (or kept on request); output_dir writes the two canonical files
#   B7  the file route on a named-chromosome map returns names and positions
#       (it returned "1" and NA Start/End up to 0.7.8) and agrees with B1
#   B8  SelectionTools' global st.input.dir / st.output.dir: gone (after
#       rm(list = ls())) -> supplied for the call and removed again; set to a
#       directory -> refused with the remedy
#   B9  the pre-0.7.10 names still work: they warn (deprecated) and return
#       exactly what st_ld_blocks() / st_ld_blocks_files() return
# =============================================================================

library(testthat)
library(haploStack)
withr::local_options(list(haploStack.st_formats = TRUE))   # the gate (0.9.0)

test_that("B1: planted LD groups come back as exactly those blocks", {
  skip_if_not_installed("SelectionTools")
  toy <- st_toy()
  r <- st_ld_blocks(toy, ld_threshold = 0.5, tolerance = 0,
                                    attach = FALSE, verbose = FALSE)
  expect_named(r, c("block_marker_map", "block_info", "provenance"))
  bm <- r$block_marker_map
  expect_setequal(bm$Marker, dimnames(toy$geno_arrays)[[2L]])
  expect_identical(unname(partition(bm)),
                   unname(partition(data.frame(Marker = bm$Marker,
                                               Block = planted_group(bm$Marker)))))
  expect_true(all(grepl("^b[0-9]{6}$", bm$Block)))        # SelectionTools' IDs
  expect_true(all(c("Chrom", "Block", "Start", "End", "Num_SNPs") %in% names(r$block_info)))
  expect_false(anyNA(r$block_info$Start))
})

test_that("B2: attach = TRUE keeps effects, phenotype, FAM and chromosome names", {
  skip_if_not_installed("SelectionTools")
  toy <- st_toy()
  ids <- dimnames(toy$geno_arrays)[[1L]]
  mk  <- dimnames(toy$geno_arrays)[[2L]]
  toy$marker_effects <- list(Yield = sim_effect_df(mk, alleles = c("A", "T")))
  toy <- suppressMessages(suppressWarnings(attach_phenotype(
    toy, data.frame(Sample = ids, Yield = seq_along(ids)))))
  toy$sample_info$Parent1[3:4] <- ids[1L]
  toy$sample_info$Parent2[3:4] <- ids[2L]

  out <- suppressMessages(st_ld_blocks(toy, 0.5, 0, verbose = FALSE))
  expect_true(is_geno_meta(out, verbose = FALSE))
  expect_identical(out$marker_effects, toy$marker_effects)
  expect_identical(out$phenotype, toy$phenotype)
  expect_identical(out$sample_info, toy$sample_info)
  bmi <- out$block_marker_map_info
  expect_setequal(unique(bmi$Chrom), c("1H", "2H", "3H"))
  expect_identical(bmi$Chrom[match(mk, bmi$Marker)], toy$map_info$Chrom)
  expect_identical(length(unique(bmi$Block)), 11L)
})

test_that("B3: tolerance absorbs an LD-breaking marker inside a block", {
  skip_if_not_installed("SelectionTools")
  noisy <- st_toy(noise_in = 5L)            # group 5: 5 markers, the 3rd is noise
  g5 <- sprintf("g05_m%d", 1:5)
  blocks_of <- function(tol) {
    bm <- st_ld_blocks(noisy, 0.5, tol, attach = FALSE,
                                       verbose = FALSE)$block_marker_map
    bm$Block[match(g5, bm$Marker)]
  }
  expect_length(unique(blocks_of(0L)), 3L)   # m1+m2 | m3 | m4+m5
  expect_length(unique(blocks_of(1L)), 1L)   # one block of five
})

test_that("B4: unplaced markers are left out of blocking and become singletons", {
  skip_if_not_installed("SelectionTools")
  toy <- st_toy()
  toy$map_info$Pos[toy$map_info$Marker == "g01_m2"]   <- NA
  toy$map_info$Chrom[toy$map_info$Marker == "g05_m3"] <- "UN"
  toy <- suppressMessages(new_geno_meta(
    geno_arrays = toy$geno_arrays, map_info = toy$map_info,
    metadata = toy$metadata, sample_info = toy$sample_info))

  r <- st_ld_blocks(toy, 0.5, 0, attach = FALSE, verbose = FALSE)
  expect_false(any(c("g01_m2", "g05_m3") %in% r$block_marker_map$Marker))
  expect_identical(nrow(r$block_marker_map), length(dimnames(toy$geno_arrays)[[2L]]) - 2L)

  out <- suppressMessages(st_ld_blocks(toy, 0.5, 0, verbose = FALSE))
  bmi <- out$block_marker_map_info
  expect_identical(bmi$Block[bmi$Marker == "g01_m2"], "dummy_g01_m2")
  expect_identical(bmi$Block[bmi$Marker == "g05_m3"], "dummy_g05_m3")
})

test_that("B5: the blocking rule is validated before SelectionTools is touched", {
  toy <- st_toy()
  for (thr in list(0, -0.1, 1.2, NA_real_, c(0.3, 0.5), "0.5"))
    expect_error(st_ld_blocks(toy, ld_threshold = thr),
                 "ld_threshold", info = format(thr))
  for (tol in list(-1, 1.5, NA_real_, c(1, 2)))
    expect_error(st_ld_blocks(toy, 0.5, tolerance = tol),
                 "tolerance", info = format(tol))
  expect_error(st_ld_blocks(list(a = 1)), "valid GenoMeta")
  f <- tempfile()
  expect_error(suppressWarnings(st_ld_blocks_files(f, f, ld_threshold = 0),
                                classes = "hs_warning_deprecated"), "ld_threshold")
})

test_that("B6: no file lands in the working directory; temp dir handling; output_dir", {
  skip_if_not_installed("SelectionTools")
  toy <- st_toy()
  wd <- tempfile("wd_"); dir.create(wd)
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  before_tmp <- list.files(tempdir(), pattern = "^hs_st_")

  invisible(st_ld_blocks(toy, 0.5, 0, attach = FALSE, verbose = FALSE))
  expect_length(list.files(wd, all.files = TRUE, no.. = TRUE), 0L)
  expect_setequal(list.files(tempdir(), pattern = "^hs_st_"), before_tmp)
  expect_identical(normalizePath(getwd()), normalizePath(wd))  # working dir restored

  kept <- NULL
  withCallingHandlers(
    invisible(st_ld_blocks(toy, 0.5, 0, attach = FALSE,
                                           keep_files = TRUE, verbose = FALSE)),
    message = function(m) {
      txt <- trimws(conditionMessage(m))
      if (grepl("kept in ", txt, fixed = TRUE)) kept <<- sub(".*kept in ", "", txt)
      invokeRestart("muffleMessage")
    })
  expect_true(dir.exists(kept))
  expect_true(all(c("hs_geno.txt", "hs_map.txt") %in% list.files(kept)))
  expect_identical(readLines(file.path(kept, "hs_map.txt"), n = 1L), "name\tchrom\tpos")
  unlink(kept, recursive = TRUE)

  od <- file.path(wd, "blocks")
  out <- suppressMessages(st_ld_blocks(toy, 0.5, 0, output_dir = od,
                                                       prefix = "toy", verbose = FALSE))
  expect_true(all(c("toy_block_info.csv", "toy_block_marker_map.csv") %in% list.files(od)))
})

test_that("B7: the file route keeps chromosome names and positions, and agrees with B1", {
  skip_if_not_installed("SelectionTools")
  toy <- st_toy()
  d <- tempfile("st_files_"); dir.create(d)
  arr <- toy$geno_arrays
  cells <- matrix(paste0(arr[, , 1], "/", arr[, , 2]), nrow = dim(arr)[1L],
                  dimnames = dimnames(arr)[1:2])
  gf <- file.path(d, "geno.txt"); mf <- file.path(d, "map.txt")
  writeLines(c(paste(c("m/i", rownames(cells)), collapse = "\t"),
               paste(colnames(cells), apply(t(cells), 1L, paste, collapse = "\t"),
                     sep = "\t")), gf)
  writeLines(c("name\tchrom\tpos",
               paste(toy$map_info$Marker, toy$map_info$Chrom, toy$map_info$Pos,
                     sep = "\t")), mf)

  fr <- suppressWarnings(classes = "hs_warning_deprecated",
        suppressMessages(st_ld_blocks_files(gf, mf, 0.5, 0,
                                                 output_dir = file.path(d, "out"))))
  expect_setequal(unique(fr$block_marker_map$Chrom), c("1H", "2H", "3H"))
  expect_false(anyNA(fr$block_info$Start))
  expect_false(anyNA(fr$block_info$End))
  mr <- st_ld_blocks(toy, 0.5, 0, attach = FALSE, verbose = FALSE)
  expect_identical(unname(partition(fr$block_marker_map)),
                   unname(partition(mr$block_marker_map)))
})

test_that("B8: SelectionTools' global dirs — supplied when gone, restored when set; cwd untouched", {
  skip_if_not_installed("SelectionTools")
  requireNamespace("SelectionTools", quietly = TRUE)
  ge   <- globalenv()
  vars <- c("st.input.dir", "st.output.dir")
  kept <- mget(vars, envir = ge, ifnotfound = list(NULL))
  on.exit({
    suppressWarnings(rm(list = vars, envir = ge))
    for (v in vars) if (!is.null(kept[[v]])) assign(v, kept[[v]], envir = ge)
  }, add = TRUE)
  toy <- st_toy()

  # gone, as after rm(list = ls()): SelectionTools alone would stop with
  # "object 'st.input.dir' not found"
  suppressWarnings(rm(list = vars, envir = ge))
  r <- st_ld_blocks(toy, 0.5, 0, attach = FALSE, verbose = FALSE)
  expect_identical(length(unique(r$block_marker_map$Block)), 11L)
  expect_false(any(vapply(vars, exists, logical(1L), envir = ge, inherits = FALSE)))

  # set to a directory, as SelectionTools' manual shows: used to be refused;
  # since 0.9.0 (D1-8) the bridge points both at its private directory for the
  # call and restores the user's values afterwards
  assign("st.input.dir", "input", envir = ge)
  assign("st.output.dir", "", envir = ge)
  r2 <- st_ld_blocks(toy, 0.5, 0, attach = FALSE, verbose = FALSE)
  expect_identical(r2$block_marker_map, r$block_marker_map)
  expect_identical(get("st.input.dir", envir = ge), "input")   # restored
  expect_identical(get("st.output.dir", envir = ge), "")

  # the working directory is neither changed nor written to
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  r3 <- st_ld_blocks(toy, 0.5, 0, attach = FALSE, verbose = FALSE)
  expect_identical(normalizePath(getwd()), normalizePath(wd))
  expect_identical(list.files(wd, all.files = TRUE, no.. = TRUE), character(0L))
  expect_identical(r3$block_marker_map, r$block_marker_map)
})

test_that("B9: the pre-0.7.10 names are defunct (0.9.0); the file route still matches st_ld_blocks()", {
  skip_if_not_installed("SelectionTools")
  toy <- st_toy()
  expect_error(define_haplotype_blocks_meta(toy, 0.5, 0, attach = FALSE, verbose = FALSE),
               "st_ld_blocks()", class = "hs_error_defunct", fixed = TRUE)

  d <- tempfile("st_dep_"); dir.create(d)
  arr <- toy$geno_arrays
  cells <- matrix(paste0(arr[, , 1], "/", arr[, , 2]), nrow = dim(arr)[1L],
                  dimnames = dimnames(arr)[1:2])
  gf <- file.path(d, "geno.txt"); mf <- file.path(d, "map.txt")
  writeLines(c(paste(c("m/i", rownames(cells)), collapse = "\t"),
               paste(colnames(cells), apply(t(cells), 1L, paste, collapse = "\t"),
                     sep = "\t")), gf)
  writeLines(c("name\tchrom\tpos",
               paste(toy$map_info$Marker, toy$map_info$Chrom, toy$map_info$Pos,
                     sep = "\t")), mf)
  expect_warning(new_f <- suppressMessages(st_ld_blocks_files(gf, mf, 0.5, 0,
                                               output_dir = file.path(d, "new"))),
                 class = "hs_warning_deprecated")      # row 27 a
  expect_error(define_haplotype_blocks(gf, mf, 0.5, 0, output_dir = file.path(d, "old")),
               class = "hs_error_defunct")
  # the file route is st_convert() + st_ld_blocks(): same tables as the container route
  new_r <- st_ld_blocks(toy, 0.5, 0, attach = FALSE, verbose = FALSE)
  expect_identical(new_f$block_marker_map, new_r$block_marker_map)
})
