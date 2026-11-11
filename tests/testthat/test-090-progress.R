# =============================================================================
# test-090-progress.R — the progress and cancel protocol on the LD resolvers
# and the compiled kernels (design D1 §3.5, step D1-5)
#
#   Q1  hapS_ld_blocks(progress =): one event per chromosome, in the container's
#       chromosome names, with the block count; the same blocks and IDs as the
#       call without a callback (per-chromosome == whole-panel gate)
#   Q2  FALSE cancels with hs_error_cancelled; `partial` holds the finished
#       chromosomes' block marker map, `at` says where
#   Q3  st_ld_blocks(progress =): start / done events; FALSE at start cancels
#   Q4  the kernels carry an interrupt check (static), and the containment of
#       progressr's options still holds through the per-chromosome loop
#   Q5  a non-function `progress` is refused before any work
# =============================================================================

library(testthat)
library(haploStack)

.pg_panel <- function(seed = 3L, n_chrom = 3L, n_markers = 120L) {
  sim_geno_meta(n_samples = 30L, n_markers = n_markers, n_chrom = n_chrom,
                block_size = 4L, encoding = "nucleotide", seed = seed)
}

test_that("Q1: hapS_ld_blocks() ticks per chromosome and gives the whole-panel blocks", {
  skip_if_not_installed("HapSelect")
  m  <- .pg_panel()
  ev <- list()
  plain <- hapS_ld_blocks(m, attach = FALSE, verbose = FALSE)
  withp <- hapS_ld_blocks(m, attach = FALSE, verbose = FALSE,
                          progress = function(e) { ev[[length(ev) + 1L]] <<- e; TRUE })
  expect_identical(withp$block_marker_map, plain$block_marker_map)
  expect_identical(withp$block_info, plain$block_info)
  chroms <- unique(as.character(m$map_info$Chrom))
  expect_length(ev, length(chroms))
  expect_identical(vapply(ev, `[[`, character(1L), "label"), chroms)
  expect_identical(vapply(ev, `[[`, integer(1L), "i"), seq_along(chroms))
  expect_true(all(vapply(ev, function(e) e$n == length(chroms), logical(1L))))
  expect_identical(unique(vapply(ev, `[[`, character(1L), "stage")), "ld_blocks")
  expect_identical(unique(vapply(ev, `[[`, character(1L), "unit")), "chromosome")
  n_blk <- vapply(ev, function(e) e$payload$n_blocks, integer(1L))
  per_chrom <- table(unique(plain$block_marker_map[, c("Chrom", "Block")])$Chrom)
  expect_identical(unname(n_blk), as.integer(per_chrom[chroms]))
  expect_true(all(diff(vapply(ev, `[[`, numeric(1L), "elapsed")) >= 0))
})

test_that("Q2: FALSE cancels hapS_ld_blocks() with the finished chromosomes as partial", {
  skip_if_not_installed("HapSelect")
  m <- .pg_panel()
  plain <- hapS_ld_blocks(m, attach = FALSE, verbose = FALSE)
  err <- tryCatch(hapS_ld_blocks(m, attach = FALSE, verbose = FALSE,
                                 progress = function(e) e$i < 2L),
                  error = identity)
  expect_s3_class(err, "hs_error_cancelled")
  expect_identical(err$at$unit, "chromosome")
  expect_identical(err$at$i, 2L)
  part <- err$partial
  expect_s3_class(part, "data.frame")
  expect_named(part, c("Chrom", "Marker", "Block"))
  done <- unique(as.character(m$map_info$Chrom))[1:2]
  expect_setequal(unique(part$Chrom), done)
  ref <- plain$block_marker_map[plain$block_marker_map$Chrom %in% done, ]
  expect_identical(part[order(part$Marker), ], ref[order(ref$Marker), ],
                   ignore_attr = TRUE)
  # the container is untouched and the options were put back
  expect_true(is_geno_meta(m))
})

test_that("Q3: st_ld_blocks() emits start / done; FALSE at start cancels", {
  skip_if_not_installed("SelectionTools")
  m <- .pg_panel(seed = 5L, n_chrom = 2L, n_markers = 60L)
  ev <- list()
  res <- st_ld_blocks(m, attach = FALSE, verbose = FALSE, st_formats = TRUE,
                      progress = function(e) { ev[[length(ev) + 1L]] <<- e; TRUE })
  expect_length(ev, 2L)
  expect_identical(vapply(ev, `[[`, character(1L), "label"), c("start", "done"))
  expect_identical(vapply(ev, `[[`, character(1L), "unit"), c("run", "run"))
  expect_identical(ev[[2L]]$payload$n_blocks,
                   length(unique(res$block_marker_map$Block)))
  err <- tryCatch(st_ld_blocks(m, attach = FALSE, verbose = FALSE, st_formats = TRUE,
                               progress = function(e) FALSE), error = identity)
  expect_s3_class(err, "hs_error_cancelled")
  expect_null(err$partial)
})

test_that("Q4: the kernels check for interrupts; progressr options are put back", {
  src <- testthat::test_path("..", "..", "src")
  skip_if_not(dir.exists(src), "src/ not reachable")
  for (f in c("local_gebv.cpp", "marker_stats.cpp")) {
    code <- readLines(file.path(src, f), warn = FALSE)
    expect_gte(sum(grepl("checkUserInterrupt", code, fixed = TRUE)),
               if (f == "marker_stats.cpp") 2L else 1L, label = f)
  }
  skip_if_not_installed("HapSelect")
  m <- .pg_panel()
  withr::local_options(list(progressr.handlers = NULL, progressr.enable = NULL))
  invisible(hapS_ld_blocks(m, attach = FALSE, verbose = FALSE,
                           progress = function(e) TRUE))
  expect_null(getOption("progressr.handlers"))
  expect_null(getOption("progressr.enable"))
})

test_that("Q5: a non-function progress is refused before any work", {
  skip_if_not_installed("HapSelect")
  m <- .pg_panel()
  expect_error(hapS_ld_blocks(m, attach = FALSE, verbose = FALSE, progress = 1),
               "'progress' must be NULL or a function", class = "hs_error_input")
})
