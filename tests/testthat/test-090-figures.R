# =============================================================================
# test-090-figures.R — figure recipes and the one renderer (design D3.2; F5)
#
#   G1  no plotter opens a device with its defaults; plot_meta_circos() returns
#       visibly (draw = FALSE) and invisibly with draw = TRUE
#   G2  hs_render_figure() restores the previous device, par() and circlize
#   G3  png / pdf / svg written; unknown extension, missing directory and
#       overwrite = FALSE are classed errors; create_dir = TRUE creates
#   G4  an hs_figure recipe rebuilds an identical ggplot data layer; hs_ref is
#       resolved against the owning result; a recipe prints and renders;
#       hs_figures(build = TRUE)
#   G4b a nested hs_figure_set renders one file per leaf, named by its path
#   G4c save_journal_figure() goes through the renderer and restores the device
# (G5–G8 land with D3.6 / D3.4 / D3.9.)
# =============================================================================

.fig_fixture <- function() {
  ds <- sim_dataset(n_samples = 12L, n_markers = 24L, block_size = 4L, seed = 11L)
  meta <- ds$geno_meta
  meta$marker_effects <- list(Trait = ds$eff_nuc)
  analysis <- suppressWarnings(suppressMessages(analyze_block_catalogue(
    geno_meta = meta, effect_list = list(T1 = ds$eff_nuc, T2 = ds$eff_nuc))))
  bi <- as.data.frame(build_block_location_info(meta)$block_info)
  genes <- data.table::data.table(
    chrom = bi$Chrom[1L], start = bi$Start[1L], end = bi$End[1L],
    strand = "+", gene_id = "gene:G1", gene_name = "GeneOne",
    biotype = "protein_coding")
  ann <- suppressWarnings(compute_block_annotation(meta, genes, flank_bp = 0))
  list(meta = meta, analysis = analysis, ann = ann)
}

.err_class <- function(expr) tryCatch(expr, error = identity)

test_that("G1: plotters open no device; plot_meta_circos() returns visibly by default", {
  fx <- .fig_fixture()
  before <- grDevices::dev.cur(); devs <- grDevices::dev.list()

  p <- plot_block_effect_genome(fx$analysis)
  expect_s3_class(p, "ggplot")
  expect_identical(grDevices::dev.cur(), before)

  skip_if_not_installed("circlize")
  v <- withVisible(suppressMessages(
    plot_meta_circos(track_annotation(fx$ann), meta = fx$meta)))
  expect_true(v$visible)
  expect_s3_class(v$value, "hs_meta_circos")
  expect_identical(grDevices::dev.cur(), before)
  expect_identical(grDevices::dev.list(), devs)

  # draw = TRUE still draws at once (on the test null device) and returns invisibly
  v2 <- withVisible(suppressMessages(
    plot_meta_circos(track_annotation(fx$ann), meta = fx$meta, draw = TRUE)))
  expect_false(v2$visible)
  expect_s3_class(v2$value, "hs_meta_circos")
  expect_identical(grDevices::dev.cur(), before)
})

test_that("G2: hs_render_figure() restores device, par and circlize state", {
  skip_if_not_installed("circlize")
  fx  <- .fig_fixture()
  cir <- suppressMessages(plot_meta_circos(track_annotation(fx$ann), meta = fx$meta))
  before <- grDevices::dev.cur(); mar <- graphics::par("mar")
  f <- tempfile(fileext = ".png")
  out <- suppressMessages(hs_render_figure(cir, f))
  expect_true(file.exists(f) && file.size(f) > 0)
  expect_identical(normalizePath(out, winslash = "/"), normalizePath(f, winslash = "/"))
  expect_identical(grDevices::dev.cur(), before)
  expect_identical(graphics::par("mar"), mar)
  expect_equal(circlize::circos.par("start.degree"), 0)   # cleared, back to default
  expect_false(file.exists(file.path(getwd(), "Rplots.pdf")))
})

test_that("G3: devices by extension; classed errors for extension, directory, overwrite", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  d <- withr::local_tempdir()
  for (ext in c("png", "pdf")) {
    f <- file.path(d, paste0("fig.", ext))
    hs_render_figure(p, f, width_mm = 89, height_mm = 60, dpi = 72)
    expect_true(file.size(f) > 0, info = ext)
  }
  if (isTRUE(unname(capabilities("cairo")))) {
    f <- file.path(d, "fig.svg")
    hs_render_figure(p, f, width_mm = 89, height_mm = 60)
    expect_true(file.size(f) > 0)
  }
  if (requireNamespace("pdftools", quietly = TRUE)) {
    sz <- pdftools::pdf_pagesize(file.path(d, "fig.pdf"))
    expect_equal(sz$width,  89 / 25.4 * 72, tolerance = 0.01)
    expect_equal(sz$height, 60 / 25.4 * 72, tolerance = 0.01)
  }

  e <- .err_class(hs_render_figure(p, file.path(d, "fig.bmp")))
  expect_s3_class(e, "hs_error_input"); expect_match(conditionMessage(e), "extension")

  e <- .err_class(hs_render_figure(p, file.path(d, "nodir", "fig.png")))
  expect_s3_class(e, "hs_error_input"); expect_match(conditionMessage(e), "create_dir")
  hs_render_figure(p, file.path(d, "nodir", "fig.png"), create_dir = TRUE, dpi = 72)
  expect_true(file.exists(file.path(d, "nodir", "fig.png")))

  e <- .err_class(hs_render_figure(p, file.path(d, "fig.png"), overwrite = FALSE))
  expect_s3_class(e, "hs_error_input"); expect_match(conditionMessage(e), "fig.png", fixed = TRUE)

  e <- .err_class(hs_render_figure(list(a = 1), file.path(d, "x.png")))
  expect_s3_class(e, "hs_error_input")
  e <- .err_class(hs_render_figure(data.frame(a = 1), file.path(d, "x.png")))
  expect_s3_class(e, "hs_error_input")
})

test_that("G4: a recipe rebuilds the same ggplot; hs_ref resolves against the result", {
  fx <- .fig_fixture()
  fig <- haploStack:::.hs_figure("genome", "plot_block_effect_genome",
                                 list(analysis = haploStack:::.hs_ref("analysis")),
                                 width_mm = 280, height_mm = 150, dpi = 150)
  expect_s3_class(fig, "hs_figure")
  expect_identical(fig$file, "genome.png")
  res <- haploStack:::.hs_result(list(analysis = fx$analysis, figures = list(g = fig)),
                                 "hs_block_analysis")
  expect_s3_class(hs_figures(res)$g, "hs_figure")
  built  <- hs_figures(res, build = TRUE)$g
  direct <- plot_block_effect_genome(fx$analysis)
  expect_s3_class(built, "ggplot")
  expect_identical(ggplot2::layer_data(built), ggplot2::layer_data(direct))
  expect_identical(hs_build_figure(fig, result = res)$labels, direct$labels)

  e <- .err_class(hs_build_figure(fig))                 # the ref needs its result
  expect_s3_class(e, "hs_error_input"); expect_match(conditionMessage(e), "result")
  e <- .err_class(hs_build_figure(haploStack:::.hs_figure("x", "no_such_plotter_xyz")))
  expect_s3_class(e, "hs_error_input"); expect_match(conditionMessage(e), "plotter")
  e <- .err_class(hs_build_figure(fig, result = haploStack:::.hs_result(list(a = 1), "hs_block_analysis")))
  expect_s3_class(e, "hs_error_input"); expect_match(conditionMessage(e), "no element")

  # a self-contained recipe renders at its stored size and prints on the current device
  fig2 <- haploStack:::.hs_figure("genome", "plot_block_effect_genome",
                                  list(analysis = fx$analysis),
                                  width_mm = 120, height_mm = 80, dpi = 72)
  f <- tempfile(fileext = ".pdf")
  hs_render_figure(fig2, f)
  expect_true(file.size(f) > 0)
  if (requireNamespace("pdftools", quietly = TRUE)) {
    sz <- pdftools::pdf_pagesize(f)
    expect_equal(sz$width,  120 / 25.4 * 72, tolerance = 0.01)
    expect_equal(sz$height,  80 / 25.4 * 72, tolerance = 0.01)
  }
  before <- grDevices::dev.cur()
  expect_invisible(print(fig2))
  expect_identical(grDevices::dev.cur(), before)

  # survives a round trip as data
  rt <- readRDS({ tf <- tempfile(fileext = ".rds"); saveRDS(res, tf); tf })
  expect_identical(hs_figures(rt)$g, fig)
})

test_that("G4b: an hs_figure_set renders one file per leaf, named by its path", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  set <- haploStack:::.hs_figure_set(list(gen_1 = list(F01 = p, F02 = p), gen_2 = list(F01 = p)))
  expect_s3_class(set, "hs_figure_set")
  expect_true(is.list(set)); expect_length(set, 2L)
  d <- withr::local_tempdir()
  files <- hs_render_figure(set, file.path(d, "aln.png"), width_mm = 80, height_mm = 50, dpi = 72)
  expect_length(files, 3L)
  expect_true(all(file.exists(files)))
  expect_identical(basename(files), c("aln_gen_1_F01.png", "aln_gen_1_F02.png", "aln_gen_2_F01.png"))
  e <- .err_class(hs_render_figure(haploStack:::.hs_figure_set(list(gen_1 = list())),
                                   file.path(d, "none.png")))
  expect_s3_class(e, "hs_error_input")
})

test_that("G4c: save_journal_figure() renders through hs_render_figure() and restores the device", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + theme_hs_journal()
  before <- grDevices::dev.cur()
  d <- withr::local_tempdir()
  f <- file.path(d, "sub", "fig.pdf")                 # create_dir = TRUE kept
  expect_identical(save_journal_figure(p, f, width = "single"), f)
  expect_true(file.size(f) > 0)
  expect_identical(grDevices::dev.cur(), before)
  e <- .err_class(save_journal_figure(p, file.path(d, "fig.svg")))
  expect_s3_class(e, "hs_error_input"); expect_match(conditionMessage(e), "unsupported")
})

test_that("G6: hs_chrom_levels() follows the D3 §3.6 table", {
  f <- hs_chrom_levels
  expect_identical(f(c("2H", "1H", "7H", "Un")),            c("1H", "2H", "7H", "Un"))
  expect_identical(f(c("10H", "2H", "1H", "Un")),           c("1H", "2H", "10H", "Un"))
  expect_identical(f(c("10", "2", "1", "Un")),              c("1", "2", "10", "Un"))
  expect_identical(f(c("chr10", "chr2", "chr1")),           c("chr1", "chr2", "chr10"))
  expect_identical(f(c("chr2", "chr1", "chrUn", "chr10")),  c("chr1", "chr2", "chr10", "chrUn"))
  expect_identical(f(c("10A", "1B", "1A", "2A")),           c("1A", "1B", "2A", "10A"))
  expect_identical(f(c("X", "2", "1", "Y", "10", "MT")),    c("1", "2", "10", "MT", "X", "Y"))
  expect_identical(f(c("Chr1", "Chr10", "Chr2", "ChrSy", "ChrUn")),
                   c("Chr1", "Chr2", "Chr10", "ChrSy", "ChrUn"))
  expect_identical(f(c("1", NA, "1", "3")), c("1", "3"))
  expect_identical(f(character(0)), character(0))
  expect_identical(f(c(3L, 10L, 1L)), c("1", "3", "10"))
  expect_identical(haploStack:::.hs_chrom_rank(c("10", "Un", "2")), c(2L, 3L, 1L))
  # the old helpers are aliases now
  expect_identical(haploStack:::.natural_chrom_levels(c("10H", "2H")), c("2H", "10H"))
  expect_identical(haploStack:::.hs_chr_order(c("10", "2", "Un")), c("2", "10", "Un"))
})

test_that("G7: Manhattan, circos, block key and tables agree on c('10','1','2','Un')", {
  chroms <- c("10", "1", "2", "Un")
  want   <- c("1", "2", "10", "Un")
  d <- data.table::data.table(Chrom = chroms, Pos = c(5, 7, 3, 9))
  expect_identical(as.character(haploStack:::.create_genomic_positions(d)$chr_info$Chrom), want)
  expect_identical(haploStack:::.hapS_chrom_key(chroms)$Chrom, want)
  # an 11-chromosome panel: every table in natural chromosome order
  m <- sim_geno_meta(n_samples = 8L, n_markers = 44L, n_chrom = 11L, block_size = 4L, seed = 5L)
  lev <- as.character(1:11)
  st <- compute_marker_stats(m)
  expect_identical(unique(as.character(st$Chrom)), lev)
  mk <- geno_meta_views(m)$marker_meta
  expect_identical(unique(as.character(mk$chromosome)), lev)
  bi <- build_block_location_info(m)$block_info
  expect_identical(unique(as.character(bi$Chrom)), lev)
  core <- haploStack:::.catalogue_core(m)
  blk_chr <- as.character(bi$Chrom[match(core$block_order, bi$Block)])
  expect_identical(unique(blk_chr), lev)
})

# ---- D3.9 small items (design D3 §3.8) -----------------------------------------

test_that("G5: an all-negative rect ring gets a non-degenerate ramp; limits fix the scale", {
  skip_if_not_installed("circlize")
  rf <- haploStack:::.hs_ring_fill
  neg <- circos_track(data.frame(Chrom = "1", Start = c(1, 10, 20), End = c(5, 15, 25),
                                 v = c(-3, -1, -0.2)), type = "rect", value_col = "v")
  r <- rf(neg)
  expect_equal(r$lab, c(-3, -0.2))
  expect_false(identical(r$fun(-3), r$fun(-0.2)))              # viz/F06: not one flat grey
  pos <- circos_track(data.frame(Chrom = "1", Start = 1, End = 5, v = 2), type = "rect", value_col = "v")
  expect_equal(rf(pos)$lab, c(0, 2))
  mix <- circos_track(data.frame(Chrom = "1", Start = c(1, 10), End = c(5, 15), v = c(-2, 4)),
                      type = "rect", value_col = "v")
  expect_equal(rf(mix)$lab, c(-4, 0, 4))                      # symmetric about zero
  # intake 12: a shared scale
  lim <- circos_track(data.frame(Chrom = "1", Start = c(1, 10), End = c(5, 15), v = c(-2, 4)),
                      type = "rect", value_col = "v", limits = c(-10, 10))
  expect_equal(rf(lim)$lab, c(-10, 0, 10))
  expect_identical(lim$limits, c(-10, 10))
  expect_identical(rf(lim)$fun(4), rf(lim)$fun(4))
  expect_false(identical(rf(lim)$fun(4), rf(mix)$fun(4)))     # 4 is no longer the end colour
  e <- .err_class(circos_track(data.frame(Chrom = "1", Pos = 1, v = 1), type = "points",
                               value_col = "v", limits = c(0, 1)))
  expect_s3_class(e, "hs_error_input")
  e <- .err_class(circos_track(data.frame(Chrom = "1", Start = 1, End = 2, v = 1), type = "rect",
                               value_col = "v", limits = c(1, 0)))
  expect_s3_class(e, "hs_error_input")
})

test_that("viz/F21: circos_track() validates its geometry arguments", {
  d <- data.frame(Chrom = "1", Start = 1, End = 2, v = 1)
  for (bad in list(list(height_mm = 0), list(height_mm = -1), list(height_mm = Inf),
                   list(min_gap_rad = Inf), list(min_gap_rad = -0.1),
                   list(label_top_n = -1L), list(label_top_n = 2.5)))
    expect_s3_class(.err_class(do.call(circos_track, c(list(data = d, type = "rect",
                                                             value_col = "v"), bad))),
                    "hs_error_input")
  expect_s3_class(circos_track(d, type = "rect", value_col = "v", height_mm = 4), "hs_circos_track")
})

test_that("viz/F18: .pt_cex divides by the device cex as well as the point size", {
  op <- graphics::par(ps = 12, cex = 1); on.exit(graphics::par(op), add = TRUE)
  expect_equal(haploStack:::.pt_cex(6), 0.5)
  graphics::par(cex = 2)
  expect_equal(haploStack:::.pt_cex(6), 0.25)
})

test_that("viz/F23: plot_phenotype_bars(sample_levels =) tags the order as an argument", {
  df <- data.frame(Sample = c("A", "B", "C"), T1 = 1:3, stringsAsFactors = FALSE)
  p <- plot_phenotype_bars(df, sample_levels = c("C", "B", "A"))
  expect_identical(attr(p, "sample_order"), c("C", "B", "A"))
  expect_identical(attr(p, "sample_order_source"), "argument")
})

test_that("io/F26: .discrete_axis_levels() reads the tag first and notes a failed sync", {
  df <- data.frame(Sample = c("A", "B", "C"), T1 = 1:3, stringsAsFactors = FALSE)
  p <- plot_phenotype_bars(df, sample_levels = c("C", "B", "A"), orientation = "horizontal")
  expect_identical(haploStack:::.discrete_axis_levels(p, "y"), c("C", "B", "A"))
  # a continuous axis still has limits to read (no note); an empty plot has nothing
  q1 <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_length(haploStack:::.discrete_axis_levels(q1, "y"), 2L)
  # (an empty ggplot still builds default panels; a plot that cannot be built has none)
  q <- ggplot2::ggplot(data.frame(x = 1), ggplot2::aes(nope, x)) + ggplot2::geom_point()
  w <- tryCatch(haploStack:::.discrete_axis_levels(q, "y"), warning = identity)
  expect_s3_class(w, "hs_warning_fallback")
  expect_identical(w$code, "axis_sync_skipped")
  expect_identical(suppressWarnings(haploStack:::.discrete_axis_levels(q, "y")), character(0))
})

test_that("C-run-fitness-print: run_fitness_pipeline(plot_landscape = TRUE) draws nothing in a script", {
  mk  <- paste0("M", 1:4)
  eff <- list(T1 = data.frame(Marker = rep(mk, each = 2L), Allele = rep(c("A", "T"), 4L),
                              Effect = c(1, -1, 0.5, -0.5, 2, -2, 0.25, -0.25)))
  before <- grDevices::dev.cur(); devs <- grDevices::dev.list()
  fi <- suppressMessages(run_fitness_pipeline(eff, c(T1 = 1), plot_landscape = TRUE))
  expect_s3_class(attr(fi, "landscape"), "ggplot")
  expect_identical(grDevices::dev.cur(), before)
  expect_identical(grDevices::dev.list(), devs)
  fi0 <- suppressMessages(run_fitness_pipeline(eff, c(T1 = 1)))
  expect_null(attr(fi0, "landscape"))
})

test_that("marker_signal_table() returns NULL for a plot without a table", {
  q <- ggplot2::ggplot(data.frame(x = 1), ggplot2::aes(x, x)) + ggplot2::geom_point()
  expect_null(marker_signal_table(q))
})
