# =============================================================================
# test-090-analysis.R — analyze_block_catalogue() as an hs_block_analysis
# (design D3.4; §3.3, §3.4, row 12 b, intake-12 hooks)
#
#   A1  class, fields, figure recipes; nothing drawn (device untouched)
#   A2  output_dir writes the 0.8.8 set through write_outputs(); writing later
#       gives the same table bytes; overwrite = FALSE refuses; what =
#   A3  "circos" is built on the engine (hs_meta_circos in $plots, recipe in
#       $figures) and rendered to block_circos.png
#   A4  centre check with catalogue = (R2): explicit mismatch errors; the default
#       adopts the catalogue's with a note; a pre-0.9.0 catalogue -> center_unknown
#   A5  gebv_matrices = FALSE: no matrices, heatmap skipped with a note
#   A6  candidates = "complete": the candidate filter and its fallback
#   A7  unused ... warns hs_warning_input; write_outputs= deprecated
#   A8  (G8) plot_block_effect_circos() warns deprecated and, with no device
#       open, leaves no Rplots.pdf
# =============================================================================

.an_fixture <- function(seed = 5L) {
  meta <- sim_geno_meta(n_samples = 12L, n_markers = 18L, n_chrom = 3L,
                        block_size = 3L, seed = seed)
  mk <- dimnames(meta$geno_arrays)[[2L]]
  set.seed(seed)
  eff1 <- data.frame(Marker = rep(mk, each = 2L), Allele = rep(c("A", "T"), length(mk)),
                     Effect = round(rnorm(2L * length(mk), 0, 0.5), 3), stringsAsFactors = FALSE)
  eff2 <- data.frame(Marker = rep(mk, each = 2L), Allele = rep(c("A", "T"), length(mk)),
                     Effect = round(rnorm(2L * length(mk), 0, 3.0), 3), stringsAsFactors = FALSE)
  list(meta = meta, eff1 = eff1, eff2 = eff2)
}
.an_quiet <- function(expr) suppressMessages(suppressWarnings(expr))
.an_err   <- function(expr) tryCatch(expr, error = identity)

test_that("A1: an hs_block_analysis with fields, recipes and notes; no device touched", {
  fx <- .an_fixture()
  before <- grDevices::dev.cur(); devs <- grDevices::dev.list()
  res <- .an_quiet(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1, Protein = fx$eff2),
                                           plots = c("rectangular", "manhattan", "density", "heatmap")))
  expect_s3_class(res, "hs_block_analysis"); expect_s3_class(res, "hs_result")
  expect_true(all(c("subtypes", "selected", "S", "N", "top_blocks", "gebv", "plots",
                    "value_source", "normalize", "figures", "notes", "provenance") %in% names(res)))
  expect_named(res$plots, c("rectangular", "manhattan", "density", "heatmap"))
  expect_true(all(vapply(res$plots, inherits, TRUE, "ggplot")))
  expect_named(hs_figures(res), names(res$plots))
  expect_true(all(vapply(hs_figures(res), inherits, TRUE, "hs_figure")))
  expect_identical(hs_figures(res)$heatmap$file, "block_heatmap.png")
  pv <- hs_provenance(res)
  expect_identical(pv$fun, "analyze_block_catalogue")
  expect_identical(pv$traits, c("Yield", "Protein"))
  expect_identical(pv$center, "none")               # row 15 c (D3.5): "none" everywhere
  expect_identical(pv$candidates, "all"); expect_true(pv$gebv_matrices)
  expect_null(pv$catalogue_center)
  expect_identical(grDevices::dev.cur(), before); expect_identical(grDevices::dev.list(), devs)
  built <- hs_figures(res, build = TRUE)
  expect_identical(ggplot2::layer_data(built$rectangular), ggplot2::layer_data(res$plots$rectangular))
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f)); saveRDS(res, f)
  expect_identical(readRDS(f), res)
})

test_that("A2: output_dir writes the 0.8.8 set via write_outputs(); later writes match", {
  fx <- .an_fixture()
  d <- withr::local_tempdir()
  res <- .an_quiet(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1), output_dir = d,
                                           plots = c("density", "rectangular")))
  expect_setequal(list.files(d), c("block_subtypes.txt", "block_selected.txt", "block_S.txt",
                                   "block_N.txt", "block_top.txt", "block_density.png",
                                   "block_rectangular.png"))
  d2 <- withr::local_tempdir()
  m <- write_outputs(res, d2)
  expect_identical(names(m), c("path", "kind", "object", "bytes"))
  expect_setequal(basename(m$path), list.files(d))
  expect_identical(sort(m$object[m$kind == "table"]),
                   sort(c("subtypes", "selected", "S", "N", "top_blocks")))
  for (f in grep("[.]txt$", list.files(d), value = TRUE))
    expect_identical(unname(tools::md5sum(file.path(d2, f))),
                     unname(tools::md5sum(file.path(d, f))), info = f)
  expect_true(all(m$bytes[m$kind == "figure"] > 0))
  e <- .an_err(write_outputs(res, d2))
  expect_s3_class(e, "hs_error_input")
  m3 <- write_outputs(res, withr::local_tempdir(), what = "tables")
  expect_identical(nrow(m3), 5L); expect_true(all(m3$kind == "table"))
  # the table written equals the result's table
  st <- as.data.frame(data.table::fread(file.path(d, "block_top.txt")))
  expect_identical(st$Block, as.character(res$top_blocks$Block))
})

test_that("A3: 'circos' is built on the engine and rendered by write_outputs()", {
  skip_if_not_installed("circlize")
  fx <- .an_fixture()
  before <- grDevices::dev.cur()
  res <- .an_quiet(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1, Protein = fx$eff2),
                                           plots = "circos"))
  expect_s3_class(res$plots$circos, "hs_meta_circos")
  expect_setequal(res$plots$circos$tracks, c("Yield", "Protein"))   # outermost first
  expect_identical(grDevices::dev.cur(), before)
  fig <- hs_figures(res)$circos
  expect_s3_class(fig, "hs_figure"); expect_identical(fig$plotter, "plot_meta_circos")
  obj <- hs_build_figure(fig, result = res)
  expect_s3_class(obj, "hs_meta_circos")
  d <- withr::local_tempdir()
  m <- .an_quiet(write_outputs(res, d))
  expect_true("block_circos.png" %in% basename(m$path))
  expect_gt(file.size(file.path(d, "block_circos.png")), 0)
  expect_identical(grDevices::dev.cur(), before)
  expect_identical(nrow(hs_notes(res)[hs_notes(res)$code == "plot_skipped", ]), 0L)
})

test_that("A4: the centre is checked against a passed catalogue (R2)", {
  fx <- .an_fixture()
  # a catalogue computed under the 0.8.8 centring; the analysis default is "none" (row 15 c)
  cat_ <- .an_quiet(run_haplotype_catalogue(fx$meta, effect_list = list(Yield = fx$eff1),
                                            center = "allele_freq"))
  expect_identical(hs_provenance(cat_)$center, "allele_freq")
  # explicit, different -> error
  e <- .an_err(.an_quiet(analyze_block_catalogue(catalogue = cat_, center = "none",
                                                 plots = character(0))))
  expect_s3_class(e, "hs_error_input"); expect_match(conditionMessage(e), "center")
  # explicit, same -> fine, no note
  r <- .an_quiet(analyze_block_catalogue(catalogue = cat_, center = "allele_freq",
                                         plots = character(0)))
  expect_identical(hs_provenance(r)$center, "allele_freq")
  expect_false("center_adopted" %in% hs_notes(r)$code)
  # default -> adopts the catalogue's with a note
  r <- .an_quiet(analyze_block_catalogue(catalogue = cat_, plots = character(0)))
  expect_identical(hs_provenance(r)$center, "allele_freq")
  expect_identical(hs_provenance(r)$catalogue_center, "allele_freq")
  expect_true("center_adopted" %in% hs_notes(r)$code)
  # same default on both sides -> no note
  cat_n <- .an_quiet(run_haplotype_catalogue(fx$meta, effect_list = list(Yield = fx$eff1)))
  r <- .an_quiet(analyze_block_catalogue(catalogue = cat_n, plots = character(0)))
  expect_identical(hs_provenance(r)$center, "none")
  expect_false("center_adopted" %in% hs_notes(r)$code)
  # a pre-0.9.0 catalogue (plain list) -> center_unknown
  old <- lapply(unclass(cat_), function(tr) tr[c("hap_catalogue", "gebv_matrix")])
  r <- .an_quiet(analyze_block_catalogue(catalogue = old, plots = character(0)))
  expect_true("center_unknown" %in% hs_notes(r)$code)
  expect_null(hs_provenance(r)$catalogue_center)
})

test_that("A5: gebv_matrices = FALSE skips the matrices and the heatmap, with a note", {
  fx <- .an_fixture()
  res <- .an_quiet(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                           plots = c("density", "heatmap"), gebv_matrices = FALSE))
  expect_true(all(vapply(res$gebv, is.null, TRUE)))
  expect_named(res$plots, "density")
  nt <- hs_notes(res)
  expect_true(any(nt$code == "plot_skipped" & grepl("heatmap", nt$message)))
  expect_false(hs_provenance(res)$gebv_matrices)
  # on the catalogue path too
  cat_ <- .an_quiet(run_haplotype_catalogue(fx$meta, effect_list = list(Yield = fx$eff1)))
  r2 <- .an_quiet(analyze_block_catalogue(catalogue = cat_, plots = character(0), gebv_matrices = FALSE))
  expect_true(all(vapply(r2$gebv, is.null, TRUE)))
  # the selection itself is unchanged by the opt-out
  r1 <- .an_quiet(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1), plots = character(0)))
  expect_identical(r1$selected, res$selected)
})

test_that("A6: candidates = 'complete' drops gapped signatures per block, with fallback", {
  st <- data.table::data.table(
    Block   = c("b1", "b1", "b1", "b2", "b2"),
    Subtype = c("AT", "A-", "TT", "A-", "-T"),
    Trait   = "Y", Vbase = c(1, 5, 2, 3, 4), Count = c(3L, 1L, 2L, 1L, 1L))
  out <- haploStack:::.filter_complete_candidates(st, "-")
  expect_identical(out[Block == "b1", Subtype], c("AT", "TT"))   # gapped one dropped
  expect_identical(out[Block == "b2", Subtype], c("A-", "-T"))   # none complete -> all kept
  expect_false("ok" %in% names(out))
  fx <- .an_fixture()
  r_all <- .an_quiet(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1), plots = character(0)))
  r_cmp <- .an_quiet(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1), plots = character(0),
                                             candidates = "complete"))
  expect_identical(hs_provenance(r_cmp)$candidates, "complete")
  expect_identical(nrow(r_cmp$subtypes), nrow(r_all$subtypes))   # the reported table is not filtered
  fill <- haploStack:::.hs_fill_na(NULL)
  expect_false(any(grepl(fill, r_cmp$selected$Subtype, fixed = TRUE) &
                     r_cmp$selected$Block %in%
                       r_all$subtypes[!grepl(fill, Subtype, fixed = TRUE), unique(Block)]))
  expect_error(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1), candidates = "nope"))
})

test_that("A7: unused ... warns; write_outputs= is deprecated", {
  fx <- .an_fixture()
  expect_warning(suppressMessages(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                          plots = character(0), bogus = 1)),
                 class = "hs_warning_input")
  expect_warning(suppressMessages(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                          plots = character(0), write_outputs = TRUE)),
                 class = "hs_warning_deprecated")
  e <- .an_err(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1), output_dir = 1))
  expect_s3_class(e, "hs_error_input")
})

test_that("A8 (G8): plot_block_effect_circos() warns deprecated and leaves no Rplots.pdf", {
  skip_if_not_installed("circlize")
  fx <- .an_fixture()
  res <- .an_quiet(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1, Protein = fx$eff2),
                                           plots = character(0)))
  d <- withr::local_tempdir(); withr::local_dir(d)
  # close the run-wide null device so no device is open, then restore one for teardown
  grDevices::graphics.off()
  on.exit(grDevices::pdf(NULL), add = TRUE)
  expect_identical(unname(grDevices::dev.cur()), 1L)
  expect_warning(p <- plot_block_effect_circos(res, link_top = TRUE),
                 class = "hs_warning_deprecated")
  expect_s3_class(p, "recordedplot")
  expect_identical(unname(grDevices::dev.cur()), 1L)                  # its null device closed again
  expect_false(file.exists(file.path(d, "Rplots.pdf")))
})
