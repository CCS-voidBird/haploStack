# =============================================================================
# test-circos-meta.R — BAN (#033) stage 6: the circos meta plot family, the
# marker-annotation overlay normaliser, and the journal figure standard
# (dev/figure_standards_nature.md). The visual pass on real data is the SIM
# test drive (data/SIM/.../haplostack_annotation/circos/).
# =============================================================================

# A blocked container + analysis + annotation shared across the tests.
.circos_fixture <- function() {
  ds <- sim_dataset(n_samples = 12L, n_markers = 24L, block_size = 4L,
                    seed = 11L)
  meta <- ds$geno_meta
  meta$marker_effects <- list(Trait = ds$eff_st)   # for the GenoMeta stat route
  analysis <- suppressWarnings(suppressMessages(analyze_block_catalogue(
    geno_meta = meta, effect_list = list(T1 = ds$eff_st, T2 = ds$eff_st))))
  bi <- as.data.frame(build_block_location_info(meta)$block_info)
  genes <- data.table::data.table(
    chrom = bi$Chrom[1L], start = bi$Start[1L], end = bi$End[1L],
    strand = "+", gene_id = "gene:G1", gene_name = "GeneOne",
    biotype = "protein_coding")
  ann <- suppressWarnings(
    compute_block_annotation(meta, genes, flank_bp = 0))
  list(meta = meta, analysis = analysis, ann = ann, bi = bi)
}


# --- track specs ---------------------------------------------------------------

test_that("circos_track validates its data contract", {
  expect_error(circos_track(data.frame(x = 1), type = "rect"), "Chrom")
  expect_error(circos_track(data.frame(Chrom = "1"), type = "rect"),
               "Start/End")
  expect_error(circos_track(data.frame(Chrom = "1", Pos = 5),
                            type = "points"), "value_col")
  expect_error(circos_track(data.frame(Chrom = "1", Pos = 5, v = 1),
                            type = "points", value_col = "nope"), "not in")
  expect_error(circos_track(data.frame(Chrom = "1", Pos = 5, v = 1),
                            type = "labels", value_col = "v"), "Name")

  tr <- circos_track(data.frame(Chrom = "1", Start = 10, End = 20, v = 2),
                     type = "rect", value_col = "v", name = "ring")
  expect_s3_class(tr, "hs_circos_track")
  expect_identical(tr$name, "ring")
  expect_equal(tr$data$.val, 2)

  # a rect track from point data widens later at draw time; Pos accepted
  tr2 <- circos_track(data.frame(Chrom = "1", Pos = 10, v = 1),
                      type = "rect", value_col = "v")
  expect_equal(tr2$data$Start, tr2$data$End)
  # a point track from interval data takes midpoints
  tr3 <- circos_track(data.frame(Chrom = "1", Start = 10, End = 30, v = 1),
                      type = "line", value_col = "v")
  expect_equal(tr3$data$Pos, 20)
})

test_that("the five constructors build consuming the standard formats", {
  fx <- .circos_fixture()

  t_ann <- track_annotation(fx$ann)                    # annotation list
  expect_s3_class(t_ann, "hs_circos_track")
  expect_identical(t_ann$type, "rect")
  expect_identical(t_ann$name, "gene_count")
  t_ann2 <- track_annotation(as.data.frame(fx$ann$blocks))   # bare df source
  expect_equal(nrow(t_ann2$data), nrow(t_ann$data))
  m2 <- attach_block_annotation(fx$meta, fx$ann)       # container source
  t_ann3 <- track_annotation(m2)
  expect_equal(nrow(t_ann3$data), nrow(t_ann$data))
  expect_error(track_annotation(fx$ann, value = "nope"), "no 'nope' column")

  t_eff <- track_block_effect(fx$analysis, trait = "T1")
  expect_identical(t_eff$type, "rect")
  expect_error(track_block_effect(fx$analysis, trait = "zzz"), "not in the")

  vdf <- data.frame(Chrom = fx$bi$Chrom, Start_Pos = fx$bi$Start,
                    End_Pos = fx$bi$End, Block_Var = seq_len(nrow(fx$bi)),
                    Trait = "T1")
  t_var <- track_block_variance(vdf, trait = "T1")
  expect_identical(t_var$type, "line")
  expect_equal(nrow(t_var$data), nrow(fx$bi))

  stat_df <- data.frame(Chrom = "1", Pos = c(1e6, 2e6), s = c(0.1, 0.4))
  t_stat <- track_marker_stat(stat_df, stat = "s")
  expect_identical(t_stat$type, "points")
  t_stat2 <- track_marker_stat(fx$meta, stat = "Trait")  # GenoMeta route
  expect_true(all(t_stat2$data$.val > 0))

  lab_df <- data.frame(Chrom = "1", Pos = 5e6, Name = "QTL_A", w = 1)
  t_lab <- track_marker_labels(lab_df, value_col = "w")
  expect_identical(t_lab$type, "labels")
})


# --- the meta plot -------------------------------------------------------------

test_that("plot_meta_circos renders, orders tracks, and re-renders by width", {
  skip_if_not_installed("circlize")
  fx <- .circos_fixture()
  tracks <- list(track_annotation(fx$ann),
                 track_block_effect(fx$analysis, trait = "T1"),
                 track_marker_labels(
                   data.frame(Chrom = "1", Pos = 5e6, Name = "L1", w = 1),
                   value_col = "w"))

  expect_error(plot_meta_circos(list(1)), "circos_track")
  expect_error(plot_meta_circos(tracks, rel_heights = 1), "one value per")

  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 89 / 25.4, height = 89 / 25.4, units = "in",
                 res = 150, type = "cairo")
  obj <- plot_meta_circos(tracks, meta = fx$meta)
  grDevices::dev.off()
  expect_s3_class(obj, "hs_meta_circos")
  expect_identical(obj$tracks, c("gene_count", "T1", "loci"))
  expect_identical(obj$chr_levels, c("1", "2", "3"))
  expect_true(file.exists(f) && file.size(f) > 0)

  # render() re-draws at another declared width (the saver's contract)
  f2 <- tempfile(fileext = ".png")
  grDevices::png(f2, width = 183 / 25.4, height = 183 / 25.4, units = "in",
                 res = 150, type = "cairo")
  suppressMessages(obj$render(183))
  grDevices::dev.off()
  expect_true(file.size(f2) > 0)
})

test_that("the angular deconfliction keeps top-ranked, gap-separated labels", {
  theta <- c(0, 0.05, 1, 2)
  value <- c(5, 10, 1, 3)
  kept <- .hs_deconflict(theta, value, n = 10L, min_gap = 0.14)
  expect_true(2L %in% kept)        # highest value wins its neighbourhood
  expect_false(1L %in% kept)       # 0.05 rad from a kept label
  expect_setequal(kept, c(2L, 3L, 4L))
  expect_length(.hs_deconflict(theta, value, n = 2L, min_gap = 0.14), 2L)
})


# --- the journal standard ------------------------------------------------------

test_that("theme_hs_journal applies the figure-standards contract", {
  th <- theme_hs_journal()
  expect_s3_class(th, "theme")
  expect_s3_class(th$plot.title, "element_blank")
  expect_s3_class(th$panel.grid, "element_blank")
  expect_identical(th$axis.text$colour, "black")
  expect_gte(th$axis.text$size, 5)
  p <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point() +
    theme_hs_journal()
  expect_s3_class(p, "ggplot")
})

test_that("save_journal_figure sizes, warns on titles, refuses recorded plots", {
  p <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()

  f <- tempfile(fileext = ".pdf")
  save_journal_figure(p + theme_hs_journal(), f, width = "single")
  expect_true(file.exists(f) && file.size(f) > 0)
  skip_on_cran()
  if (requireNamespace("pdftools", quietly = TRUE)) {
    sz <- pdftools::pdf_pagesize(f)       # points: 89 mm = 252.3 pt
    expect_equal(sz$width, 89 / 25.4 * 72, tolerance = 0.01)
  }

  expect_warning(save_journal_figure(p + ggplot2::ggtitle("NO"),
                                     tempfile(fileext = ".pdf")),
                 "title")

  grDevices::png(tempfile(fileext = ".png"))
  plot(1:3)
  rp <- grDevices::recordPlot()
  grDevices::dev.off()
  expect_error(save_journal_figure(rp, tempfile(fileext = ".pdf")),
               "recorded")

  skip_if_not_installed("circlize")
  fx <- .circos_fixture()
  cir <- suppressMessages(
    plot_meta_circos(track_annotation(fx$ann), meta = fx$meta, draw = FALSE))
  f3 <- tempfile(fileext = ".pdf")
  suppressMessages(save_journal_figure(cir, f3, width = "double"))
  expect_true(file.size(f3) > 0)
  if (requireNamespace("pdftools", quietly = TRUE)) {
    sz <- pdftools::pdf_pagesize(f3)
    expect_equal(sz$width, 183 / 25.4 * 72, tolerance = 0.01)
    expect_equal(sz$height, 183 / 25.4 * 72, tolerance = 0.01)
  }
})


# --- the Q25.o overlay normaliser ---------------------------------------------

test_that("standardize_marker_annotations resolves both routes and warns unplaced", {
  fx <- .circos_fixture()
  map <- fx$meta$map_info

  # Marker route
  s1 <- standardize_marker_annotations(
    data.frame(Marker = map$Marker[2L], Name = "A"), meta = fx$meta)
  expect_identical(s1$Chrom, map$Chrom[2L])
  expect_equal(s1$Pos, as.numeric(map$Pos[2L]))
  expect_equal(s1$matched_dist_bp, 0)
  expect_false(is.na(s1$Block))

  # positional route: within threshold -> nearest marker; beyond -> unplaced
  near <- data.frame(Chrom = map$Chrom[3L], Pos = as.numeric(map$Pos[3L]) + 500,
                     Name = "B", Type = "QTL")
  far  <- data.frame(Chrom = "1", Pos = 1, Name = "C")
  expect_warning(
    s2 <- standardize_marker_annotations(rbind(near, transform(far, Type = NA)),
                                         meta = fx$meta, match_bp = 1000),
    "match no marker")
  expect_identical(s2$Marker[1L], map$Marker[3L])
  expect_equal(s2$matched_dist_bp[1L], 500)
  expect_identical(s2$Type[1L], "QTL")
  expect_true(is.na(s2$Marker[2L]) && s2$Name[2L] == "C")

  # no meta: positional table normalised without matching
  s3 <- standardize_marker_annotations(far)
  expect_true(is.na(s3$Marker) && s3$Pos == 1)
  expect_error(standardize_marker_annotations(data.frame(Marker = "x",
                                                         Name = "n")),
               "needs 'meta'")
  expect_error(standardize_marker_annotations(data.frame(Pos = 1)), "Name")
})

test_that("attach_marker_annotations stores custom, coexists and carries forward", {
  fx <- .circos_fixture()
  map <- fx$meta$map_info
  df <- data.frame(Marker = map$Marker[c(1L, 10L)], Name = c("L1", "L2"))

  # custom alone on a container with no gene annotation
  m1 <- attach_marker_annotations(fx$meta, df, source = "unit test")
  expect_true(is_geno_meta(m1))
  expect_identical(names(m1$annotation), "custom")
  expect_equal(m1$metadata$annotation$custom$n_matched, 2L)

  # gene tables + custom coexist, in either attach order
  m2 <- attach_block_annotation(m1, fx$ann)
  expect_setequal(names(m2$annotation),
                  c("blocks", "features", "marker_gene", "custom"))
  expect_identical(m2$annotation$custom, m1$annotation$custom)
  m2b <- attach_marker_annotations(attach_block_annotation(fx$meta, fx$ann), df)
  expect_setequal(names(m2b$annotation),
                  c("blocks", "features", "marker_gene", "custom"))

  # filter prunes custom rows whose marker left; re-blocking drops everything
  keep <- map$Marker[1:6]
  f <- suppressWarnings(filter_geno_meta(m2, keep_markers = keep))
  expect_identical(f$annotation$custom$Name, "L1")
  bmap <- sim_block_map(m2$map_info, block_size = 6L)
  expect_message(m3 <- attach_block_info(m2, bmap), "annotation was dropped")
  expect_null(m3$annotation)

  # a corrupted custom marker fails validation
  m4 <- m1
  m4$annotation$custom$Marker[1L] <- "NOT_A_MARKER"
  expect_false(is_geno_meta(m4, verbose = FALSE))
})
