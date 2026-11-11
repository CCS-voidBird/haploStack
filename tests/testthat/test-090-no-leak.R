# =============================================================================
# test-090-no-leak.R — the server contract (design D1 §3.1 / §7, step D1-10)
#
# One job per worker: a route through the package must leave the R process as
# it found it. A snapshot of the process state is taken, "job A" runs the block
# route end to end (text read, blocks, catalogue in memory, analysis, a
# rendered figure, a two-cycle simulation, an rrBLUP fit), and the snapshot is
# taken again: identical. Gated variants add genomicSimulation, HapSelect,
# SelectionTools and BayesC where installed. Then hs_reset_session() puts a
# dirtied process back.
#
#   NL1  the block route leaks nothing
#   NL2  genomicSimulation (gated)
#   NL3  HapSelect (gated)
#   NL4  SelectionTools (gated)
#   NL5  BayesC (gated)
#   NL6  hs_reset_session()
# =============================================================================

.nl_snapshot <- function() {
  o <- options()
  keep <- grep("^(haploStack|progressr|future)\\.", names(o), value = TRUE)
  list(
    options = o[sort(keep)],
    seed    = if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
                get(".Random.seed", envir = globalenv()) else NULL,
    plan    = if (requireNamespace("future", quietly = TRUE)) class(future::plan())[1L] else NA,
    wd      = normalizePath(getwd()),
    globals = sort(ls(globalenv(), all.names = TRUE)),
    devices = names(grDevices::dev.list()),
    circos  = if ("circlize" %in% loadedNamespaces())
                tryCatch(circlize::circos.par("gap.degree"), error = function(e) NA) else NA,
    scratch = sort(list.files(tempdir(), pattern = "^hs_")),
    gs_open = haploStack:::.hs_gs$open,
    sinks   = sink.number(),
    cons    = nrow(showConnections(all = FALSE))
  )
}

# Loading a namespace is not a leak: future sets its options and circlize its
# parameters on first use, so both are warmed before the first snapshot.
.nl_warm <- function() {
  if (requireNamespace("future", quietly = TRUE)) invisible(future::plan())
  if (requireNamespace("circlize", quietly = TRUE)) invisible(circlize::circos.par("gap.degree"))
  if (requireNamespace("rrBLUP", quietly = TRUE)) loadNamespace("rrBLUP")
  if (requireNamespace("BGLR", quietly = TRUE)) loadNamespace("BGLR")
  if (requireNamespace("ggplot2", quietly = TRUE)) loadNamespace("ggplot2")
  invisible(.nl_snapshot())
}

.nl_job <- function() {
  m0  <- sim_geno_meta(n_samples = 12L, n_markers = 24L, n_chrom = 2L, block_size = 4L,
                       encoding = "nucleotide", seed = 31L)
  eff <- sim_effect_df(m0$map_info$Marker, alleles = c("A", "T"), seed = 31L)
  d   <- withr::local_tempdir(.local_envir = parent.frame())
  # a text read (the detector runs under a local seed)
  suppressMessages(write_geno_meta(m0, file.path(d, "g"), preset = "PL"))
  m   <- suppressMessages(read_plink_ped(file.path(d, "g.ped"), file.path(d, "g.map")))
  m   <- suppressMessages(attach_block_info(m, sim_block_map(m$map_info, block_size = 4L)))
  # the catalogue and its analysis, in memory
  cat_res <- suppressMessages(suppressWarnings(run_haplotype_catalogue(
    m, effect_list = list(T1 = eff), output_dir = NULL, save_profile = FALSE)))
  ana <- suppressMessages(suppressWarnings(analyze_block_catalogue(catalogue = cat_res)))
  # a rendered figure (circlize + a device)
  if (requireNamespace("circlize", quietly = TRUE)) {
    fig <- plot_meta_circos(list(track_block_effect(ana, trait = "T1")), meta = m, draw = FALSE)
    suppressMessages(hs_render_figure(fig, file.path(d, "ring.png"), width_mm = 89))
  }
  # a two-cycle simulation
  sim <- suppressMessages(suppressWarnings(run_stacking_simulation(
    m, effect_df = eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", output_dir = NULL, verbose = FALSE)))
  # an rrBLUP fit
  if (requireNamespace("rrBLUP", quietly = TRUE)) {
    ph <- data.frame(Sample = dimnames(m$geno_arrays)[[1L]],
                     Y = haploStack:::.with_seed(2L, stats::rnorm(12L)))
    suppressMessages(fit_genomic_prediction(m, phenotype = ph, model = "rrblup", trait = "Y"))
  }
  invisible(list(meta = m, eff = eff, sim = sim))
}

test_that("NL1: the block route leaves the process as it found it", {
  skip_if_not_installed("rrBLUP")
  set.seed(7L); runif(1L)                      # a drawn stream, as a worker would have
  .nl_warm()
  before <- .nl_snapshot()
  .nl_job()
  after <- .nl_snapshot()
  for (k in names(before))
    expect_identical(after[[k]], before[[k]], label = paste("process state:", k))
})

test_that("NL2: a genomicSimulation job leaves no session and no scratch", {
  skip_if_not_installed("genomicSimulation")
  m0  <- sim_geno_meta(n_samples = 8L, n_markers = 12L, n_chrom = 1L, block_size = 3L,
                       encoding = "nucleotide", seed = 33L)
  eff <- sim_effect_df(m0$map_info$Marker, alleles = c("A", "T"), seed = 33L)
  .nl_warm()
  before <- .nl_snapshot()
  suppressMessages(suppressWarnings(run_stacking_simulation(
    m0, effect_df = eff, replications = 1L, select_cycles = 2L, top_founders = 4L,
    num_offspring = 2L, scheme = "TS", cross_backend = "genomicSimulation",
    backend_opts = list(map_unit = "bp"), output_dir = NULL, verbose = FALSE)))
  after <- .nl_snapshot()
  for (k in names(before))
    expect_identical(after[[k]], before[[k]], label = paste("process state:", k))
  expect_false(haploStack:::.hs_gs$open)
})

test_that("NL3: HapSelect blocking restores the plan and the progressr options", {
  skip_if_not_installed("HapSelect")
  m0 <- sim_geno_meta(n_samples = 10L, n_markers = 16L, n_chrom = 2L, block_size = 4L,
                      encoding = "nucleotide", seed = 35L)
  .nl_warm()
  before <- .nl_snapshot()
  suppressMessages(hapS_ld_blocks(m0, ld_threshold = 0.5, tolerance = 1L, verbose = FALSE))
  after <- .nl_snapshot()
  for (k in names(before))
    expect_identical(after[[k]], before[[k]], label = paste("process state:", k))
})

test_that("NL4: SelectionTools blocking restores its globals and the working directory", {
  skip_if_not_installed("SelectionTools")
  withr::local_options(list(haploStack.st_formats = TRUE))
  m0 <- sim_geno_meta(n_samples = 10L, n_markers = 16L, n_chrom = 2L, block_size = 4L,
                      encoding = "nucleotide", seed = 37L)
  .nl_warm()
  before <- .nl_snapshot()
  suppressMessages(st_ld_blocks(m0, ld_threshold = 0.5, tolerance = 0L, attach = FALSE,
                                verbose = FALSE))
  after <- .nl_snapshot()
  for (k in names(before))
    expect_identical(after[[k]], before[[k]], label = paste("process state:", k))
})

test_that("NL5: a BayesC fit leaves no file behind", {
  skip_if_not_installed("BGLR")
  m0 <- sim_geno_meta(n_samples = 12L, n_markers = 20L, n_chrom = 1L, block_size = 4L,
                      encoding = "nucleotide", seed = 39L)
  ph <- data.frame(Sample = dimnames(m0$geno_arrays)[[1L]],
                   Y = haploStack:::.with_seed(3L, stats::rnorm(12L)))
  .nl_warm()
  before <- .nl_snapshot()
  suppressMessages(fit_genomic_prediction(m0, phenotype = ph, model = "bayesc", trait = "Y",
                                          nIter = 200L, burnIn = 50L))
  after <- .nl_snapshot()
  for (k in names(before))
    expect_identical(after[[k]], before[[k]], label = paste("process state:", k))
})

test_that("NL6: hs_reset_session() puts a dirtied process back and says what it did", {
  withr::local_options(list(haploStack.sample_order = NULL, haploStack.missing_chars = NULL,
                            haploStack.strict_missing = NULL, haploStack.use_rcpp = NULL,
                            haploStack.fill_na = NULL, haploStack.st_formats = NULL))
  options(haploStack.strict_missing = TRUE, haploStack.use_rcpp = FALSE,
          haploStack.st_formats = TRUE)
  junk <- tempfile("hs_junk_"); dir.create(junk); writeLines("x", file.path(junk, "f"))
  expect_message(done <- hs_reset_session(quiet = FALSE), "hs_reset_session", class = "hs_message")
  expect_true(all(c("haploStack.strict_missing", "haploStack.use_rcpp",
                    "haploStack.st_formats") %in% done$options))   # (+ fill_na, nulled above)
  expect_null(getOption("haploStack.strict_missing"))
  expect_true(getOption("haploStack.use_rcpp"))
  expect_false(getOption("haploStack.st_formats"))
  expect_true(basename(junk) %in% done$scratch)
  expect_false(dir.exists(junk))
  expect_false(done$gs_session)
  # nothing to do: an empty report, no error
  done2 <- hs_reset_session()
  expect_identical(done2$options, character(0L))
  expect_identical(done2$scratch, character(0L))
  # a genomicSimulation session is closed
  if (requireNamespace("genomicSimulation", quietly = TRUE)) {
    m0 <- sim_geno_meta(n_samples = 4L, n_markers = 6L, n_chrom = 1L, seed = 41L,
                        encoding = "nucleotide")
    h <- suppressMessages(suppressWarnings(gs_load_geno_meta(m0, map_unit = "bp")))
    expect_true(haploStack:::.hs_gs$open)
    done3 <- hs_reset_session()
    expect_true(done3$gs_session)
    expect_false(haploStack:::.hs_gs$open)
  }
})
