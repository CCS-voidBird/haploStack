# =============================================================================
# test-090-converter.R — the Shiny format converter as an app directory
# (design D5-4, decision-form row 28 a)
#
#   K1  the app directory ships under inst/shiny and the launcher finds it
#   K2  shiny::testServer: a PED + MAP read completes when the reader warns
#       (the warning goes to the log, cont/F23) and the output uses preset
#       names (no output_format)
# =============================================================================

.k_fixture <- function() {
  m <- sim_geno_meta(4, 6, n_chrom = 1L, seed = 11L)
  d <- withr::local_tempdir(.local_envir = parent.frame())
  suppressMessages(write_geno_meta(m, file.path(d, "g"), preset = "PL"))
  # a block file whose markers are mostly unknown: the resolver warns (< 50 % matched)
  blk <- data.frame(Chrom = m$map_info$Chrom[1L],
                    Marker = c(m$map_info$Marker[1L], "not_a_marker_1", "not_a_marker_2"),
                    Block = "B1")
  utils::write.csv(blk, file.path(d, "blocks.csv"), row.names = FALSE)
  list(ped = file.path(d, "g.ped"), map = file.path(d, "g.map"),
       blk = file.path(d, "blocks.csv"), meta = m)
}

test_that("K1: the app directory is in inst/shiny and the launcher resolves it", {
  d <- haploStack:::.hs_converter_dir()
  expect_true(file.exists(file.path(d, "app.R")))
  expect_match(d, "shiny[/\\\\]format_converter$")
  # the old in-namespace UI / server closures are gone
  expect_false(exists(".fca_ui", envir = asNamespace("haploStack"), inherits = FALSE))
  expect_false(exists(".fca_server", envir = asNamespace("haploStack"), inherits = FALSE))
})

test_that("K2: a read that warns completes, the warning is logged, output uses presets", {
  skip_if_not_installed("shiny")
  fx <- .k_fixture()
  app <- shiny::shinyAppDir(haploStack:::.hs_converter_dir())
  shiny::testServer(app, {
    session$setInputs(input_format = "plink_ped", pos_unit = "bp",
                      ped_file = data.frame(datapath = fx$ped, stringsAsFactors = FALSE),
                      map_file = data.frame(datapath = fx$map, stringsAsFactors = FALSE),
                      block_file = data.frame(datapath = fx$blk, stringsAsFactors = FALSE),
                      output_format = "GS", out_prefix = "x")
    session$setInputs(btn_run = 1L)
    res <- geno_meta_r()
    expect_false(is.null(res$meta))                    # the read completed
    expect_true(any(grepl("^\\[WARNING\\]", res$log)))   # ... and the warning was logged
    expect_true(any(grepl("marker", res$log[grepl("WARNING", res$log)])))
    expect_identical(res$meta$metadata$n_markers, 6L)
    cv <- converted_r()
    expect_identical(cv$output_format, "GS")
    expect_match(output$summary_text, "Output format:  GS", fixed = TRUE)
  })
})
