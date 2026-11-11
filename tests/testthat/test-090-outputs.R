# =============================================================================
# test-090-outputs.R — the result contract base (design D3.1; F4 / A1)
#
# The entry-point tests (O1–O3, O6, O7, O9) land with D3.3; this file covers
# the base on stub results:
#   O4  overwrite = FALSE errors with hs_error_input naming the files
#   O5  file stems: "Yield (t/ha)" -> "Yield_t_ha" with a note; collisions error
#   O8  notes: a forced fallback appears once, with n; muffle swallows it
#   O10 constructor + accessors store/read either way (attributes vs elements)
#   O11 print()/summary() open no device; write_outputs() default errors;
#       a result survives saveRDS()/readRDS()
# =============================================================================

.stub_catalogue <- function() {
  tr <- list(gebv_matrix = data.frame(Block = c("b1", "b2"), S1 = c(1, 2)),
             block_summary = data.frame(Block = "b1", Effect = 0.5),
             figures = list(manhattan = structure(list(name = "m"), class = "hs_figure")))
  haploStack:::.hs_result(list(Yield = tr, Height = tr), "hs_catalogue",
                          notes = haploStack:::.hs_notes_df("x", "c", "m", 2L),
                          provenance = haploStack:::.hs_provenance("run_haplotype_catalogue",
                                                                   label = "0.5_t0"))
}
.stub_analysis <- function() {
  haploStack:::.hs_result(list(subtypes = data.frame(Block = "b1"), top = data.frame(a = 1:3),
                               figures = list(circos = structure(list(name = "c"), class = "hs_figure"))),
                          "hs_block_analysis",
                          provenance = haploStack:::.hs_provenance("analyze_block_catalogue"))
}

test_that("O4: .hs_write_table() refuses to overwrite unless asked, naming the files", {
  d <- withr::local_tempdir()
  p <- file.path(d, "t.txt")
  m <- haploStack:::.hs_write_table(data.frame(a = 1:2, b = c("x", "y")), p, object = "T/a")
  expect_true(file.exists(p))
  expect_identical(names(m), c("path", "kind", "object", "bytes"))
  expect_identical(m$kind, "table"); expect_identical(m$object, "T/a")
  expect_identical(m$bytes, as.numeric(file.info(p)$size))
  expect_identical(readLines(p), c("a\tb", "1\tx", "2\ty"))       # fwrite, tab, as 0.8.8
  e <- tryCatch(haploStack:::.hs_write_table(data.frame(a = 1), p, overwrite = FALSE),
                error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "t.txt", fixed = TRUE)
  expect_identical(e$bad, p)
  expect_identical(readLines(p), c("a\tb", "1\tx", "2\ty"))       # untouched
  # directory readiness
  expect_error(haploStack:::.hs_dir_ready(file.path(d, "no"), create_dir = FALSE),
               class = "hs_error_input")
  expect_true(dir.exists(haploStack:::.hs_dir_ready(file.path(d, "yes"))))
})

test_that("O5: file stems sanitise with a note; collisions are an input error", {
  r <- haploStack:::.hs_with_notes(haploStack:::.hs_file_stem("Yield (t/ha)"), muffle = TRUE)
  expect_identical(r$value, "Yield_t_ha")
  expect_identical(r$notes$code, "file_stem_renamed")
  expect_identical(r$notes$stage, "file_stem")
  expect_false(grepl("/", r$value, fixed = TRUE))
  # a clean name gives no note
  r <- haploStack:::.hs_with_notes(haploStack:::.hs_file_stem("Yield"), muffle = TRUE)
  expect_identical(r$value, "Yield"); expect_identical(nrow(r$notes), 0L)
  # the note is a real warning outside .hs_with_notes()
  expect_warning(haploStack:::.hs_file_stem("a b"), class = "hs_warning_fallback")
  expect_error(haploStack:::.hs_file_stem("(/)"), class = "hs_error_input")
  e <- tryCatch(suppressWarnings(haploStack:::.hs_file_stems(c("Yield (t/ha)", "Yield_t_ha", "H"))),
                error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_identical(e$bad, c("Yield (t/ha)", "Yield_t_ha"))
  expect_identical(unname(suppressWarnings(haploStack:::.hs_file_stems(c("A B", "C")))),
                   c("A_B", "C"))
  # the 0.8.8 alias still reports as a message
  expect_message(haploStack:::.hs_safe_file_part("a b"), "sanitised", class = "hs_message")
})

test_that("O8: .hs_with_notes() records every haploStack warning once, with n", {
  f <- function() {
    haploStack:::.hs_note("stage1", "fallback_a", "fell back")
    haploStack:::.hs_note("stage1", "fallback_a", "fell back")
    haploStack:::.hs_warn("plain warning")
    haploStack:::.hs_note("stage2", "fallback_b", "other", extra = "e")
    42L
  }
  # not muffled: warnings propagate unchanged
  w <- character()
  r <- withCallingHandlers(haploStack:::.hs_with_notes(f()),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_identical(r$value, 42L)
  expect_identical(w, c("fell back", "fell back", "plain warning", "other"))
  n <- r$notes
  expect_identical(n$code, c("fallback_a", "hs_warning", "fallback_b"))
  expect_identical(n$n, c(2L, 1L, 1L))
  expect_identical(n$stage[1L], "stage1")
  expect_identical(n$stage[2L], "f")                 # calling function for a plain warning
  expect_identical(n$extra, c(NA, NA, "e"))          # extra columns allowed
  # muffled: nothing escapes
  expect_silent(r2 <- haploStack:::.hs_with_notes(f(), muffle = TRUE))
  expect_identical(r2$notes, n)
  # a base warning is not a haploStack note
  r3 <- suppressWarnings(haploStack:::.hs_with_notes({ warning("base"); 1 }))
  expect_identical(nrow(r3$notes), 0L)
})

test_that("O10: the constructor stores and the accessors read either shape", {
  cat_ <- .stub_catalogue()
  expect_identical(class(cat_), c("hs_catalogue", "hs_result", "list"))
  expect_identical(names(cat_), c("Yield", "Height"))        # top level stays keyed by trait
  expect_identical(hs_notes(cat_)$n, 2L)
  expect_identical(hs_provenance(cat_)$label, "0.5_t0")
  expect_identical(hs_provenance(cat_)$result_schema, 1L)
  expect_identical(names(hs_figures(cat_)), c("Yield", "Height"))
  expect_s3_class(hs_figures(cat_)$Yield$manhattan, "hs_figure")
  expect_identical(cat_$Yield$gebv_matrix$S1, c(1, 2))        # 0.8.8 access unchanged

  an <- .stub_analysis()
  expect_identical(class(an), c("hs_block_analysis", "hs_result", "list"))
  expect_identical(nrow(hs_notes(an)), 0L)
  expect_identical(hs_provenance(an)$fun, "analyze_block_catalogue")
  expect_identical(names(hs_figures(an)), "circos")
  expect_true(all(c("notes", "provenance", "figures") %in% names(an)))

  expect_error(hs_notes(list()), class = "hs_error_input")
  expect_error(haploStack:::.hs_result(1, "hs_simulation"), class = "hs_error_internal")
})

test_that("O11: print/summary open no device; default write_outputs errors; serialisable", {
  cat_ <- .stub_catalogue(); an <- .stub_analysis()
  before <- grDevices::dev.list()
  out <- capture.output(print(cat_))
  expect_match(out[1L], "hs_catalogue")
  expect_true(any(grepl("traits : Yield, Height", out, fixed = TRUE)))
  expect_true(any(grepl("notes  : 1", out, fixed = TRUE)))
  s <- summary(cat_)
  expect_s3_class(s, "summary.hs_result")
  expect_identical(s$tables$object, c("Yield/gebv_matrix", "Yield/block_summary",
                                      "Height/gebv_matrix", "Height/block_summary"))
  expect_identical(s$tables$nrow, c(2L, 1L, 2L, 1L))
  out <- capture.output(print(an))
  expect_true(any(grepl("fields : subtypes, top", out, fixed = TRUE)))
  expect_identical(summary(an)$tables$object, c("subtypes", "top"))
  expect_identical(grDevices::dev.list(), before)

  # D3.3: the catalogue method writes what the stub holds; its bogus recipe
  # becomes a figure_failed note (bytes NA), never an error.
  d_stub <- withr::local_tempdir()
  m <- suppressWarnings(write_outputs(cat_, d_stub))
  expect_identical(sum(m$kind == "table"), 4L)
  expect_identical(sum(m$kind == "figure"), 2L)
  expect_true(all(is.na(m$bytes[m$kind == "figure"])))
  expect_true(file.exists(file.path(d_stub, "Haplotype_GA_Matrix_0.5_t0_Yield.txt")))
  expect_error(write_outputs(list(), tempdir()), class = "hs_error_input")

  f <- tempfile(fileext = ".rds"); on.exit(unlink(f))
  saveRDS(cat_, f)
  expect_identical(readRDS(f), cat_)
  saveRDS(an, f)
  expect_identical(readRDS(f), an)
})

# =============================================================================
# D3.3 — the entry point: run_haplotype_catalogue() + write_outputs.hs_catalogue()
#   O1  output_dir = NULL writes nothing, returns hs_catalogue visibly
#   O2  output_dir writes exactly the 0.8.8 file set; tables identical to res
#   O3  write_outputs() later gives the same files and table bytes
#   O5' trait stems on a real run; collisions error before any file
#   O6  label: LD provenance, user map file, table, none (+ note), overrides
#   O7  empty core -> hs_error_input before any file
#   O9  survives saveRDS; recipes build against the result; derivation = inline
#   P   progress events per trait; FALSE cancels with the partial result
# =============================================================================

.cat_fixture <- function() {
  ds <- sim_dataset(n_samples = 12L, n_markers = 24L, block_size = 4L, seed = 11L)
  list(meta = ds$geno_meta, eff = ds$eff_nuc)
}
.quiet <- function(expr) suppressMessages(suppressWarnings(expr))
.as_df <- function(x) {
  x <- as.data.frame(x)
  x[] <- lapply(x, function(col) if (is.factor(col)) as.character(col) else col)
  x
}

test_that("O1: output_dir = NULL writes nothing and returns an hs_catalogue visibly", {
  fx <- .cat_fixture()
  d <- withr::local_tempdir(); withr::local_dir(d)
  v <- withVisible(.quiet(run_haplotype_catalogue(fx$meta, effect_list = list(T1 = fx$eff))))
  expect_true(v$visible)
  res <- v$value
  expect_s3_class(res, "hs_catalogue"); expect_s3_class(res, "hs_result")
  expect_identical(list.files(d, recursive = TRUE), character(0))
  expect_identical(names(res), "T1")
  expect_setequal(names(res$T1), c("gebv_matrix", "gebv_summary", "block_counts",
                                   "hap_catalogue", "block_summary", "figures"))
  expect_null(res$T1$geno_meta)
  expect_type(res$T1$block_counts$Sample, "character")
  expect_named(hs_figures(res)$T1, c("block_counts", "manhattan_variance", "manhattan_effect"))
  expect_s3_class(hs_figures(res)$T1$block_counts, "hs_figure")
  pv <- hs_provenance(res)
  expect_identical(pv$fun, "run_haplotype_catalogue")
  expect_identical(pv$center, "none")
  expect_identical(pv$traits, "T1")
  expect_identical(pv$file_stems, c(T1 = "T1"))
  expect_true(is.character(pv$label) && nzchar(pv$label))
  expect_false("file_stem_renamed" %in% hs_notes(res)$code)

  r2 <- .quiet(run_haplotype_catalogue(fx$meta, effect_list = list(T1 = fx$eff), figures = FALSE))
  expect_length(hs_figures(r2), 0L)
  expect_null(r2$T1$figures)

  blocks <- unique(as.character(fx$meta$block_marker_map_info$Block))
  r3 <- .quiet(run_haplotype_catalogue(fx$meta, effect_list = list(T1 = fx$eff),
                                       selected_blocks = data.frame(Block = blocks[1:2])))
  expect_setequal(r3$T1$gebv_matrix$Block, blocks[1:2])

  e <- tryCatch(.quiet(run_haplotype_catalogue(fx$meta, effect_list = list(T1 = fx$eff),
                                               trait_names = "nope")), error = identity)
  expect_s3_class(e, "hs_error_input"); expect_match(conditionMessage(e), "nope")
})

test_that("O2/O3: output_dir writes the 0.8.8 file set; write_outputs() later gives the same bytes", {
  fx <- .cat_fixture()
  d <- withr::local_tempdir()
  v <- withVisible(.quiet(run_haplotype_catalogue(fx$meta, effect_list = list(T1 = fx$eff),
                                                  output_dir = d)))
  expect_false(v$visible)
  res <- v$value
  L <- hs_provenance(res)$label
  expected <- c(sprintf("Haplotype_GA_Matrix_%s_T1.txt", L), sprintf("GEBV_%s_T1.txt", L),
                sprintf("Block_Counts_%s_T1.txt", L), "Block_Counts_T1.png",
                sprintf("Haplotype_Catalogue_%s_T1.txt", L), sprintf("Block_Summary_%s_T1.txt", L),
                sprintf("Manhattan_Variance_%s_T1.png", L), sprintf("Manhattan_Effect_%s_T1.png", L),
                sprintf("geno_meta_%s.rds", L))
  expect_setequal(list.files(d), expected)
  rt <- function(f) {                       # fread types "1","2" as integer; Chrom is text
    x <- as.data.frame(data.table::fread(file.path(d, f)))
    if ("Chrom" %in% names(x)) x$Chrom <- as.character(x$Chrom)
    x
  }
  expect_equal(rt(sprintf("Haplotype_GA_Matrix_%s_T1.txt", L)), .as_df(res$T1$gebv_matrix),
               ignore_attr = TRUE)
  expect_equal(rt(sprintf("Block_Summary_%s_T1.txt", L)), .as_df(res$T1$block_summary),
               ignore_attr = TRUE)
  expect_equal(rt(sprintf("GEBV_%s_T1.txt", L)), .as_df(res$T1$gebv_summary), ignore_attr = TRUE)
  bc <- rt(sprintf("Block_Counts_%s_T1.txt", L))
  expect_identical(bc$Sample, res$T1$block_counts$Sample)
  expect_identical(bc$Label, as.character(res$T1$block_counts$Label))
  expect_equal(bc$Count, res$T1$block_counts$Count)
  expect_true(is_geno_meta(readRDS(file.path(d, sprintf("geno_meta_%s.rds", L)))))
  for (f in grep("[.]png$", expected, value = TRUE))
    expect_gt(file.size(file.path(d, f)), 0)

  # O3: written later, elsewhere, from the in-memory result
  d2 <- withr::local_tempdir()
  m <- write_outputs(res, d2, geno_meta = fx$meta)
  expect_s3_class(m, "data.frame")
  expect_identical(names(m), c("path", "kind", "object", "bytes"))
  expect_setequal(basename(m$path), expected)
  expect_identical(sort(m$object[m$kind == "table"]),
                   sort(paste0("T1/", c("gebv_matrix", "gebv_summary", "block_counts",
                                        "hap_catalogue", "block_summary"))))
  for (f in grep("[.]txt$", expected, value = TRUE))
    expect_identical(unname(tools::md5sum(file.path(d2, f))),
                     unname(tools::md5sum(file.path(d, f))), info = f)
  expect_true(all(m$bytes[m$kind == "figure"] > 0))
  e <- tryCatch(write_outputs(res, d2), error = identity)         # overwrite = FALSE
  expect_s3_class(e, "hs_error_input")
  d3 <- withr::local_tempdir()
  m3 <- write_outputs(res, d3, what = "tables")
  expect_true(all(m3$kind == "table")); expect_identical(nrow(m3), 5L)
  m4 <- write_outputs(res, withr::local_tempdir(), what = "container")   # no geno_meta given
  expect_identical(nrow(m4), 0L)
  d5 <- withr::local_tempdir()
  .quiet(run_haplotype_catalogue(fx$meta, effect_list = list(T1 = fx$eff), output_dir = d5,
                                 save_profile = FALSE))
  expect_false(any(grepl("[.]rds$", list.files(d5))))
})

test_that("O5': trait names become file stems on a real run; collisions error before any file", {
  fx <- .cat_fixture()
  d <- withr::local_tempdir()
  r <- .quiet(run_haplotype_catalogue(fx$meta, effect_list = list(`Yield (t/ha)` = fx$eff),
                                      output_dir = d, save_profile = FALSE))
  expect_true(any(grepl("_Yield_t_ha[.]txt$", list.files(d))))
  expect_length(list.dirs(d, recursive = FALSE), 0L)
  expect_identical(unname(hs_provenance(r)$file_stems["Yield (t/ha)"]), "Yield_t_ha")
  expect_true("file_stem_renamed" %in% hs_notes(r)$code)
  d2 <- withr::local_tempdir()
  e <- tryCatch(.quiet(run_haplotype_catalogue(fx$meta, effect_list = list(`A B` = fx$eff, A_B = fx$eff),
                                               output_dir = d2)), error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_identical(list.files(d2), character(0))
})

test_that("O6: the output label follows the block provenance, with overrides", {
  fx <- .cat_fixture()
  run <- function(m, ...) .quiet(run_haplotype_catalogue(m, effect_list = list(T1 = fx$eff), ...))
  nb <- length(unique(fx$meta$block_marker_map_info$Block))
  m_ld <- fx$meta
  m_ld$metadata$blocks <- list(source = "hapS_ld_blocks", label = "hapS_ld_blocks_r0.5_t3",
                               params = list(ld_threshold = 0.5, tolerance = 3L), n_blocks = nb)
  expect_identical(hs_provenance(run(m_ld))$label, "0.5_t3")
  m_um <- fx$meta
  m_um$metadata$blocks <- list(source = "user_map", label = "my blocks.tsv", params = list(),
                               file = "D:/somewhere/my blocks.tsv", n_blocks = nb)
  r <- run(m_um)
  expect_identical(hs_provenance(r)$label, "my_blocks")
  expect_true("file_stem_renamed" %in% hs_notes(r)$code)
  m_df <- fx$meta
  m_df$metadata$blocks <- list(source = "user_map", label = "data.frame", params = list(),
                               file = NA_character_, n_blocks = nb)
  expect_identical(hs_provenance(run(m_df))$label, "usermap")
  m_no <- fx$meta
  m_no$metadata$blocks <- NULL
  r <- run(m_no)
  expect_identical(hs_provenance(r)$label, "0.5_t0")
  expect_true("label_legacy_default" %in% hs_notes(r)$code)
  expect_identical(hs_provenance(run(m_no, label = "mine"))$label, "mine")
  expect_identical(hs_provenance(run(m_no, ld_threshold = 0.5, tolerance = 3L))$label, "0.5_t3")
  expect_identical(hs_provenance(run(m_no, tolerance = 2))$label, "0.5_t2")
  expect_error(run(m_no, label = ""), class = "hs_error_input")
  d <- withr::local_tempdir()
  run(m_ld, output_dir = d, save_profile = FALSE)
  expect_true(file.exists(file.path(d, "Haplotype_Catalogue_0.5_t3_T1.txt")))
})

test_that("O7: an empty catalogue core is an input error before any file is written", {
  fx <- .cat_fixture()
  m <- fx$meta
  m$block_marker_map_info$Marker <- paste0("nope_", seq_len(nrow(m$block_marker_map_info)))
  d <- withr::local_tempdir()
  e <- tryCatch(.quiet(run_haplotype_catalogue(m, effect_list = list(T1 = fx$eff), output_dir = d)),
                error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "marker")
  expect_identical(list.files(d), character(0))
})

test_that("O9: the result survives saveRDS and its recipes build against it", {
  fx <- .cat_fixture()
  res <- .quiet(run_haplotype_catalogue(fx$meta, effect_list = list(T1 = fx$eff)))
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f))
  saveRDS(res, f)
  rt <- readRDS(f)
  expect_identical(rt, res)
  built <- hs_figures(rt, build = TRUE)$T1
  expect_s3_class(built$block_counts, "ggplot")
  expect_s3_class(built$manhattan_variance, "ggplot")
  expect_s3_class(built$manhattan_effect, "ggplot")
  expect_identical(attr(built$block_counts, "sample_order_source"), "positive_count")
  long <- haploStack:::.melt_gebv(res$T1$gebv_matrix)
  pos  <- unique(as.data.frame(res$T1$block_summary)[, c("Block", "Chrom", "Pos")])
  direct <- plot_manhattan_effect(merge(long, pos, by = "Block", all.x = TRUE), trait_name = "T1")
  expect_equal(ggplot2::layer_data(built$manhattan_effect), ggplot2::layer_data(direct))
  for (nm in c("gebv_matrix", "gebv_summary", "block_counts", "hap_catalogue", "block_summary"))
    expect_true(is.data.frame(rt$T1[[nm]]), info = nm)
  expect_s3_class(summary(rt), "summary.hs_result")
  expect_identical(nrow(hs_notes(rt)), nrow(hs_notes(res)))
})

test_that("progress: one event per trait; FALSE cancels with the partial result", {
  fx <- .cat_fixture()
  ev <- list()
  res <- .quiet(run_haplotype_catalogue(
    fx$meta, effect_list = list(A = fx$eff, B = fx$eff),
    progress = function(e) { ev[[length(ev) + 1L]] <<- e; TRUE }))
  expect_length(ev, 2L)
  expect_identical(vapply(ev, `[[`, "", "label"), c("A", "B"))
  expect_identical(ev[[2L]]$i, 2L); expect_identical(ev[[2L]]$n, 2L)
  expect_identical(ev[[1L]]$stage, "catalogue"); expect_identical(ev[[1L]]$unit, "trait")
  expect_identical(ev[[1L]]$payload, list(trait = "A"))
  expect_true(is.numeric(ev[[1L]]$elapsed))
  e <- tryCatch(.quiet(run_haplotype_catalogue(
    fx$meta, effect_list = list(A = fx$eff, B = fx$eff),
    progress = function(e) e$i < 1L)), error = identity)
  expect_s3_class(e, "hs_error_cancelled")
  expect_identical(names(e$partial), "A")
  expect_identical(e$at, list(unit = "trait", i = 1L, n = 2L))
  expect_error(run_haplotype_catalogue(fx$meta, effect_list = list(A = fx$eff), progress = 1),
               class = "hs_error_input")
})
