# =============================================================================
# test-marker-effects.R
#
# QA for R/marker_effects.R:
#   - export_marker_effects(): trait resolution, column contract, per-trait
#     file writing, file-name sanitising, argument exclusivity, annotate=
#   - plot_marker_effects(): multi-allele collapse, off-map markers, title
#   - a non-disturbance guard proving plot_manhattan_effect() is unchanged
# =============================================================================

library(testthat)
library(haploStack)

# Two-trait fixture: 20 markers x 3 chromosomes, blocks of 4.
eff_meta <- function(n_markers = 20L) {
  m  <- sim_geno_meta(n_samples = 10L, n_markers = n_markers,
                      n_chrom = 3L, block_size = 4L)
  mk <- dimnames(m$geno_arrays)[[2L]]
  m$marker_effects <- list(
    Yield  = sim_effect_df(mk, alleles = c("1", "4"), seed = 1L),
    Height = sim_effect_df(mk, alleles = c("1", "4"), seed = 2L))
  m
}


# ---------------------------------------------------------------------------
# export_marker_effects — table contract
# ---------------------------------------------------------------------------

test_that("a single-trait GenoMeta exports the stored columns", {
  m <- eff_meta()
  m$marker_effects <- m$marker_effects["Yield"]

  out <- export_marker_effects(m)
  expect_true(data.table::is.data.table(out))
  expect_identical(names(out), c("Trait", "Marker", "Allele", "Effect"))
  expect_identical(nrow(out), nrow(m$marker_effects$Yield))
  expect_identical(unique(out$Trait), "Yield")
  expect_identical(out$Marker, m$marker_effects$Yield$Marker)
  expect_equal(out$Effect, m$marker_effects$Yield$Effect)
  expect_identical(attr(out, "files"), character(0L))
})

test_that("multi-trait export stacks every trait and 'traits' subsets it", {
  m <- eff_meta()

  all_out <- export_marker_effects(m)
  expect_identical(nrow(all_out), 80L)
  expect_setequal(unique(all_out$Trait), c("Yield", "Height"))

  one <- export_marker_effects(m, traits = "Height")
  expect_identical(unique(one$Trait), "Height")
  expect_identical(nrow(one), 40L)

  # The caller's order is preserved, not the stored order.
  ordered <- export_marker_effects(m, traits = c("Height", "Yield"))
  expect_identical(unique(ordered$Trait), c("Height", "Yield"))
})

test_that("an unknown trait errors and lists what is available", {
  m <- eff_meta()
  expect_error(export_marker_effects(m, traits = "Protein"),
               "not in geno_meta\\$marker_effects")
  expect_error(export_marker_effects(m, traits = "Protein"), "Yield")
})

test_that("a missing marker_effects slot errors and names the fitting function", {
  m <- sim_geno_meta(block_size = 4L)
  expect_error(export_marker_effects(m), "fit_genomic_prediction")
})

test_that("bare-data.frame and Trait-column effect tables are accepted", {
  m  <- eff_meta()
  mk <- dimnames(m$geno_arrays)[[2L]]

  m$marker_effects <- sim_effect_df(mk, alleles = c("1", "4"))
  expect_identical(unique(export_marker_effects(m)$Trait), "Trait")

  stacked <- rbind(
    cbind(sim_effect_df(mk, alleles = c("1", "4"), seed = 1L), Trait = "A"),
    cbind(sim_effect_df(mk, alleles = c("1", "4"), seed = 2L), Trait = "B"))
  m$marker_effects <- stacked
  expect_setequal(unique(export_marker_effects(m)$Trait), c("A", "B"))
})


# ---------------------------------------------------------------------------
# export_marker_effects — file writing
# ---------------------------------------------------------------------------

test_that("output_dir writes one file per trait with the stored columns", {
  m   <- eff_meta()
  dir <- withr::local_tempdir()

  out <- suppressMessages(export_marker_effects(m, output_dir = dir))
  files <- attr(out, "files")
  expect_length(files, 2L)
  expect_true(all(file.exists(files)))
  expect_setequal(basename(files),
                  c("marker_effects_Yield.txt", "marker_effects_Height.txt"))

  back <- data.table::fread(file.path(dir, "marker_effects_Yield.txt"))
  expect_identical(names(back), c("Marker", "Allele", "Effect"))
  expect_identical(back$Marker, m$marker_effects$Yield$Marker)
  expect_equal(back$Effect, m$marker_effects$Yield$Effect)
})

test_that("a trait name with unsafe characters is sanitised for the file name", {
  m <- eff_meta()
  names(m$marker_effects)[1L] <- "Yield (t/ha)"
  dir <- withr::local_tempdir()

  expect_message(export_marker_effects(m, traits = "Yield (t/ha)",
                                       output_dir = dir),
                 "sanitised")
  out <- suppressMessages(
    export_marker_effects(m, traits = "Yield (t/ha)", output_dir = dir))

  expect_identical(basename(attr(out, "files")), "marker_effects_Yield_t_ha.txt")
  # The returned table keeps the original trait string.
  expect_identical(unique(out$Trait), "Yield (t/ha)")
})

test_that("output_file takes one trait only and output_dir is exclusive to it", {
  m   <- eff_meta()
  dir <- withr::local_tempdir()
  f   <- file.path(dir, "one.txt")

  expect_error(export_marker_effects(m, output_file = f), "output_dir")
  expect_error(export_marker_effects(m, output_dir = dir, output_file = f),
               "mutually exclusive")

  out <- suppressMessages(
    export_marker_effects(m, traits = "Yield", output_file = f))
  expect_identical(attr(out, "files"), f)
  expect_identical(names(data.table::fread(f)), c("Marker", "Allele", "Effect"))
})

test_that("a nested output_dir is created recursively", {
  m   <- eff_meta()
  dir <- file.path(withr::local_tempdir(), "a", "b")
  expect_false(dir.exists(dir))
  suppressMessages(export_marker_effects(m, traits = "Yield",
                                         output_dir = dir))
  expect_true(file.exists(file.path(dir, "marker_effects_Yield.txt")))
})

test_that("annotate = TRUE adds the locus columns and NA for off-map markers", {
  m <- eff_meta()
  m$marker_effects$Yield <- rbind(
    m$marker_effects$Yield,
    data.frame(Marker = "Ghost", Allele = "1", Effect = 0.5,
               stringsAsFactors = FALSE))

  out <- export_marker_effects(m, traits = "Yield", annotate = TRUE)
  expect_identical(names(out),
                   c("Trait", "Chrom", "Pos", "Block", "Marker", "Allele",
                     "Effect"))
  expect_true(all(is.na(out[Marker == "Ghost", Chrom])))
  expect_false(any(is.na(out[Marker != "Ghost", Pos])))
  # Positioned rows are in genomic order.
  pos <- out[Marker != "Ghost"]
  expect_false(is.unsorted(pos$Chrom))
})


# ---------------------------------------------------------------------------
# plot_marker_effects
# ---------------------------------------------------------------------------

test_that("the default collapse draws one point per marker", {
  m <- eff_meta()
  p <- plot_marker_effects(m, trait = "Yield")

  expect_s3_class(p, "ggplot")
  tbl <- attr(p, "effect_table")
  expect_identical(nrow(tbl), 20L)
  expect_identical(data.table::uniqueN(tbl$Marker), 20L)
  expect_true(all(c("GenomicPos", "PlotVal") %in% names(p$data)))
  expect_identical(nrow(p$data), 20L)
})

test_that("collapse = 'none' keeps every (Marker, Allele) row", {
  m <- eff_meta()
  p <- plot_marker_effects(m, trait = "Yield", collapse = "none")
  expect_identical(nrow(attr(p, "effect_table")), 40L)
})

test_that("max_abs picks the largest |Effect| per marker and keeps its sign", {
  m <- sim_geno_meta(n_samples = 6L, n_markers = 2L, n_chrom = 1L,
                     block_size = 2L)
  m$marker_effects <- list(Yield = data.frame(
    Marker = c("Mkr1", "Mkr1", "Mkr2", "Mkr2"),
    Allele = c("1", "4", "1", "4"),
    Effect = c(0.2, -0.9, 0.5, 0.1),
    stringsAsFactors = FALSE))

  tbl <- attr(plot_marker_effects(m), "effect_table")
  expect_identical(nrow(tbl), 2L)
  expect_equal(tbl[Marker == "Mkr1", Effect], -0.9)   # sign preserved
  expect_identical(tbl[Marker == "Mkr1", Allele], "4")
  expect_equal(tbl[Marker == "Mkr2", Effect], 0.5)
})

test_that("markers absent from map_info are dropped with a warning", {
  m <- eff_meta()
  m$marker_effects$Yield <- rbind(
    m$marker_effects$Yield,
    data.frame(Marker = "Ghost", Allele = "1", Effect = 9,
               stringsAsFactors = FALSE))

  expect_warning(p <- plot_marker_effects(m, trait = "Yield"), "map_info")
  expect_false("Ghost" %in% attr(p, "effect_table")$Marker)
})

test_that("an effect table with no mappable marker errors", {
  m <- eff_meta()
  m$marker_effects <- list(Yield = data.frame(
    Marker = "Ghost", Allele = "1", Effect = 1, stringsAsFactors = FALSE))
  expect_error(suppressWarnings(plot_marker_effects(m)), "map_info\\$Marker")
})

test_that("several stored traits require an explicit trait", {
  m <- eff_meta()
  expect_error(plot_marker_effects(m), "multiple traits")
  expect_error(plot_marker_effects(m), "Height")
})

test_that("the title says Marker Effect, not Block Effect", {
  m <- eff_meta()
  expect_match(plot_marker_effects(m, trait = "Yield")$labels$title,
               "^Marker Effect")
  expect_identical(plot_marker_effects(m, trait = "Yield", title = "Custom")$
                     labels$title, "Custom")
})

test_that("plot_manhattan_effect is left untouched by the new marker plot", {
  dt <- data.table::data.table(Chrom = c("1", "1", "2"),
                               Pos = c(1e6, 2e6, 1e6),
                               Effect = c(0.1, -0.2, 0.3))
  expect_match(plot_manhattan_effect(dt, trait_name = "Yield")$labels$title,
               "^Block Effect")
})

test_that("table_file still writes the plotted table but is deprecated (0.9.0)", {
  m <- eff_meta()
  f <- file.path(withr::local_tempdir(), "eff.txt")
  expect_warning(p <- suppressMessages(plot_marker_effects(m, trait = "Yield",
                                                           table_file = f)),
                 class = "hs_warning_deprecated")
  back <- data.table::fread(f)
  expect_identical(nrow(back), nrow(attr(p, "effect_table")))
  expect_equal(back$Effect, attr(p, "effect_table")$Effect)
  # the accessor reads the effect table too
  expect_identical(marker_signal_table(p), attr(p, "effect_table"))
})
