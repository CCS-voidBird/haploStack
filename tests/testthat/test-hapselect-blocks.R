# =============================================================================
# test-hapselect-blocks.R  (0.7.10)
#
# hapS_ld_blocks(), the default local LD resolver (HapSelect), promoted from
# dev/hapselect_bridge/. Toy panels: helper-ld-toy.R.
#   H1  planted LD groups come back as exactly those blocks, at tolerance 0 and
#       at the default 3; HapSelect's "<key>:<n>" block IDs
#   H2  attach = TRUE keeps effects / phenotype / FAM, and chromosome NAMES
#   H3  tolerance absorbs an LD-breaking marker inside a block
#   H4  unplaced markers are left out of blocking and become singletons
#   H5  the rule is validated before HapSelect is touched
#   H6  markers without genotype variation: no warning escapes, every marker
#       is still blocked, and the summary names them
#   H7  side effects contained: progressr's handlers option, the working
#       directory, silence with verbose = FALSE; output_dir writes two files
#   H8  on planted structure the two resolvers agree block for block
# =============================================================================

library(testthat)
library(haploStack)
withr::local_options(list(haploStack.st_formats = TRUE))   # H8 compares with SelectionTools (gate, 0.9.0)

test_that("H1: planted LD groups come back as exactly those blocks", {
  skip_if_not_installed("HapSelect")
  toy <- st_toy()
  planted <- function(bm) unname(partition(data.frame(
    Marker = bm$Marker, Block = planted_group(bm$Marker))))
  for (tol in c(0L, 3L)) {
    r <- hapS_ld_blocks(toy, tolerance = tol, attach = FALSE, verbose = FALSE)
    expect_named(r, c("block_marker_map", "block_info", "provenance"))
    bm <- r$block_marker_map
    expect_setequal(bm$Marker, dimnames(toy$geno_arrays)[[2L]])
    expect_identical(unname(partition(bm)), planted(bm), info = tol)
    expect_true(all(grepl("^[0-9]+:[0-9]+$", bm$Block)))     # HapSelect's IDs
    expect_false(anyNA(r$block_info$Start))
  }
})

test_that("H2: attach = TRUE keeps effects, phenotype, FAM and chromosome names", {
  skip_if_not_installed("HapSelect")
  toy <- st_toy()
  ids <- dimnames(toy$geno_arrays)[[1L]]
  mk  <- dimnames(toy$geno_arrays)[[2L]]
  toy$marker_effects <- list(Yield = sim_effect_df(mk, alleles = c("A", "T")))
  toy <- suppressMessages(suppressWarnings(attach_phenotype(
    toy, data.frame(Sample = ids, Yield = seq_along(ids)))))
  toy$sample_info$Parent1[3:4] <- ids[1L]
  toy$sample_info$Parent2[3:4] <- ids[2L]

  out <- hapS_ld_blocks(toy, verbose = FALSE)
  expect_true(is_geno_meta(out, verbose = FALSE))
  expect_identical(out$marker_effects, toy$marker_effects)
  expect_identical(out$phenotype, toy$phenotype)
  expect_identical(out$sample_info, toy$sample_info)
  bmi <- out$block_marker_map_info
  expect_setequal(unique(bmi$Chrom), c("1H", "2H", "3H"))
  expect_identical(bmi$Chrom[match(mk, bmi$Marker)], toy$map_info$Chrom)
  expect_identical(length(unique(bmi$Block)), 11L)
})

test_that("H3: tolerance absorbs an LD-breaking marker inside a block", {
  skip_if_not_installed("HapSelect")
  noisy <- st_toy(noise_in = 5L)            # group 5: 5 markers, the 3rd is noise
  g5 <- sprintf("g05_m%d", 1:5)
  blocks_of <- function(tol) {
    bm <- hapS_ld_blocks(noisy, tolerance = tol, attach = FALSE,
                         verbose = FALSE)$block_marker_map
    bm$Block[match(g5, bm$Marker)]
  }
  expect_length(unique(blocks_of(0L)), 3L)   # m1+m2 | m3 | m4+m5
  expect_length(unique(blocks_of(1L)), 1L)   # one block of five
})

test_that("H4: unplaced markers are left out of blocking and become singletons", {
  skip_if_not_installed("HapSelect")
  toy <- st_toy()
  toy$map_info$Pos[toy$map_info$Marker == "g01_m2"]   <- NA
  toy$map_info$Chrom[toy$map_info$Marker == "g05_m3"] <- "UN"
  toy <- suppressMessages(new_geno_meta(
    geno_arrays = toy$geno_arrays, map_info = toy$map_info,
    metadata = toy$metadata, sample_info = toy$sample_info))

  r <- hapS_ld_blocks(toy, attach = FALSE, verbose = FALSE)
  expect_false(any(c("g01_m2", "g05_m3") %in% r$block_marker_map$Marker))
  expect_identical(nrow(r$block_marker_map),
                   length(dimnames(toy$geno_arrays)[[2L]]) - 2L)

  bmi <- hapS_ld_blocks(toy, verbose = FALSE)$block_marker_map_info
  expect_identical(bmi$Block[bmi$Marker == "g01_m2"], "dummy_g01_m2")
  expect_identical(bmi$Block[bmi$Marker == "g05_m3"], "dummy_g05_m3")
})

test_that("H5: the rule is validated before HapSelect is touched", {
  toy <- st_toy()
  for (thr in list(0, -0.1, 1.2, NA_real_, c(0.3, 0.5), "0.5"))
    expect_error(hapS_ld_blocks(toy, ld_threshold = thr), "ld_threshold",
                 info = format(thr))
  for (tol in list(-1, 1.5, NA_real_, c(1, 2)))
    expect_error(hapS_ld_blocks(toy, tolerance = tol), "tolerance",
                 info = format(tol))
  expect_error(hapS_ld_blocks(toy, tol_reset = NA), "tol_reset")
  expect_error(hapS_ld_blocks(toy, method = "median"), "should be one of")
  expect_error(hapS_ld_blocks(list(a = 1)), "valid GenoMeta")
})

test_that("H6: markers without genotype variation are blocked quietly and named", {
  skip_if_not_installed("HapSelect")
  set.seed(2)
  base  <- sample(c("A", "T"), 30L, replace = TRUE)
  cells <- cbind(paste0(base, "/", base), paste0(base, "/", base),
                 rep("A/A", 30L), rep("A/T", 30L), paste0(base, "/", base))
  dimnames(cells) <- list(sprintf("L%02d", 1:30),
                          c("a1", "a2", "mono", "allhet", "a3"))
  cm <- suppressMessages(build_geno_meta(
    cells, data.frame(Chrom = "1", Marker = colnames(cells), Pos = (1:5) * 1e6),
    split_char = "/", encoding = "nucleotide"))

  msgs <- character(0)
  r <- NULL
  expect_no_warning(r <- withCallingHandlers(
    hapS_ld_blocks(cm, tolerance = 0L, attach = FALSE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    }))
  expect_setequal(r$block_marker_map$Marker, colnames(cells))
  expect_true(any(grepl("2 marker(s) without genotype variation", msgs,
                        fixed = TRUE)))
  # not singletons: at tolerance 0 the run of two becomes a block of its own
  blk <- r$block_marker_map$Block[match(c("mono", "allhet"),
                                        r$block_marker_map$Marker)]
  expect_identical(blk[1L], blk[2L])
})

test_that("H7: HapSelect's side effects are contained; output_dir writes two files", {
  skip_if_not_installed("HapSelect")
  toy <- st_toy()
  wd <- tempfile("wd_"); dir.create(wd)
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  handlers_before <- getOption("progressr.handlers")

  expect_silent(r <- hapS_ld_blocks(toy, attach = FALSE, verbose = FALSE))
  expect_identical(getOption("progressr.handlers"), handlers_before)
  expect_length(list.files(wd, all.files = TRUE, no.. = TRUE), 0L)

  od <- file.path(wd, "blocks")
  invisible(hapS_ld_blocks(toy, attach = FALSE, output_dir = od,
                           prefix = "toy", verbose = FALSE))
  expect_setequal(list.files(od),
                  c("toy_block_info.csv", "toy_block_marker_map.csv"))
})

test_that("H8: on planted structure the two resolvers agree block for block", {
  skip_if_not_installed("HapSelect")
  skip_if_not_installed("SelectionTools")
  toy <- st_toy()
  h <- hapS_ld_blocks(toy, tolerance = 0L, attach = FALSE, verbose = FALSE)
  s <- st_ld_blocks(toy, tolerance = 0L, attach = FALSE, verbose = FALSE)
  expect_identical(unname(partition(h$block_marker_map)),
                   unname(partition(s$block_marker_map)))
})
