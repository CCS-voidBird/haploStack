# Extracted from test-090-no-leak.R:132

# prequel ----------------------------------------------------------------------
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
  ana <- suppressMessages(suppressWarnings(analyze_block_catalogue(cat_res)))
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

# test -------------------------------------------------------------------------
skip_if_not_installed("BGLR")
m0 <- sim_geno_meta(n_samples = 12L, n_markers = 20L, n_chrom = 1L, block_size = 4L,
                      encoding = "nucleotide", seed = 39L)
ph <- data.frame(Sample = dimnames(m0$geno_arrays)[[1L]],
                   Y = haploStack:::.with_seed(3L, stats::rnorm(12L)))
before <- .nl_snapshot()
suppressMessages(fit_genomic_prediction(m0, phenotype = ph, model = "bayesc", trait = "Y",
                                          nIter = 200L, burnIn = 50L))
after <- .nl_snapshot()
for (k in names(before))
    expect_identical(after[[k]], before[[k]], label = paste("process state:", k))
