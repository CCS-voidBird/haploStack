# =============================================================================
# test-090-smoke-exports.R — exports that had no test of their own
# (design D5 §3.7, deps/F16; step D5-12): each runs on simulated data and
# returns the documented shape.
#
#   SM1  expand_block_info()
#   SM2  geno_meta_to_matrix() for the three presets
#   SM3  plot_block_counts(), plot_simulation_trajectories()
#   SM4  read_gs_genotype() and read_st_genotype() round trips
# =============================================================================

.sm_meta <- function() {
  sim_geno_meta(n_samples = 6L, n_markers = 12L, n_chrom = 2L, block_size = 3L,
                encoding = "nucleotide", seed = 51L)
}

test_that("SM1: expand_block_info() turns block spans into a marker-level map", {
  m  <- .sm_meta()
  bi <- as.data.frame(build_block_location_info(m)$block_info)
  bm <- expand_block_info(m$map_info, bi[, c("Chrom", "Block", "Start", "End")])
  expect_true(is.data.frame(bm))
  expect_true(all(c("Chrom", "Marker", "Block") %in% names(bm)))
  expect_setequal(bm$Marker, m$map_info$Marker)
  expect_identical(sort(unique(bm$Block)), sort(unique(bi$Block)))
  # a span covering nothing contributes no marker
  bi0 <- data.frame(Chrom = "1", Block = "none", Start = 1, End = 2)
  w <- testthat::capture_warnings(e0 <- expand_block_info(m$map_info, bi0))
  expect_true(any(grepl("contains no markers", w)))
  expect_identical(nrow(e0), 0L)
})

test_that("SM2: geno_meta_to_matrix() lays the container out in every preset", {
  m <- .sm_meta()
  pl <- geno_meta_to_matrix(m, preset = "PL")
  expect_true(is.matrix(pl$matrix))
  expect_identical(pl$output_format, "PL")
  expect_identical(dim(pl$matrix), c(6L, 12L))                     # samples x markers
  gs <- geno_meta_to_matrix(m, preset = "GS")
  expect_identical(dim(gs$matrix), c(12L, 6L))                     # markers x samples
  st <- geno_meta_to_matrix(m, preset = "ST", st_formats = TRUE)
  expect_identical(dim(st$matrix), c(12L, 6L))
  expect_true(all(grepl("^[ACGT]/[ACGT]$", st$matrix)))
  expect_error(geno_meta_to_matrix(m, preset = "ST"), class = "hs_error_input")  # gated
})

test_that("SM3: the two count / trajectory plotters return ggplots", {
  m   <- .sm_meta()
  eff <- sim_effect_df(m$map_info$Marker, alleles = c("A", "T"), seed = 51L)
  cat_res <- suppressMessages(suppressWarnings(run_haplotype_catalogue(
    m, effect_list = list(T1 = eff), output_dir = NULL, save_profile = FALSE)))
  p <- plot_block_counts(cat_res$T1$block_counts, trait_name = "T1")
  expect_s3_class(p, "ggplot")
  p2 <- plot_block_counts(cat_res$T1$block_counts, trait_name = "T1",
                          sample_order = rev(dimnames(m$geno_arrays)[[1L]]))
  expect_s3_class(p2, "ggplot")
  sim <- suppressMessages(suppressWarnings(run_stacking_simulation(
    m, effect_df = eff, replications = 2L, select_cycles = 2L, top_founders = 3L,
    num_offspring = 2L, scheme = "TS", output_dir = NULL, verbose = FALSE)))
  tr <- plot_simulation_trajectories(sim$trajectories, trait_name = "T1")
  expect_s3_class(tr, "ggplot")
  tr2 <- suppressWarnings(plot_simulation_trajectories(sim$trajectories, trait_name = "T1",
                                                       show_ribbon = FALSE, show_coverage = TRUE))
  expect_s3_class(tr2, "ggplot")
})

test_that("SM4: the genomicSimulation and SelectionTools readers read what the writers wrote", {
  m <- .sm_meta()
  d <- withr::local_tempdir()
  # each reader takes its preset's map layout (resolve_geno_format()$map_order)
  write_map <- function(preset, path) {
    ord <- resolve_geno_format(preset)$map_order
    data.table::fwrite(m$map_info[, ord, drop = FALSE], path)
    path
  }
  gs_map <- write_map("GS", file.path(d, "map_gs.csv"))
  st_map <- write_map("ST", file.path(d, "map_st.csv"))
  gs_file <- file.path(d, "g_gs.csv")
  suppressMessages(write_geno_meta(m, gs_file, preset = "GS"))
  g <- suppressMessages(read_gs_genotype(gs_file, map_file = gs_map))
  expect_true(is_geno_meta(g))
  expect_identical(dim(g$geno_arrays), dim(m$geno_arrays))
  expect_identical(dimnames(g$geno_arrays)[[2L]], dimnames(m$geno_arrays)[[2L]])
  expect_identical(unname(g$geno_arrays[, , 1L]), unname(m$geno_arrays[, , 1L]))
  st_file <- file.path(d, "g_st.csv")
  suppressMessages(write_geno_meta(m, st_file, preset = "ST", st_formats = TRUE))
  expect_error(read_st_genotype(st_file, map_file = st_map),
               class = "hs_error_input")                            # the gate
  s <- suppressMessages(suppressWarnings(read_st_genotype(
    st_file, map_file = st_map, st_formats = TRUE)))
  expect_true(is_geno_meta(s))
  expect_identical(dim(s$geno_arrays), dim(m$geno_arrays))
  expect_identical(unname(s$geno_arrays[, , 1L]), unname(m$geno_arrays[, , 1L]))
})
