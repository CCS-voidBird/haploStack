# =============================================================================
# test-genomic-prediction.R
#
# QA tests for the genomic-prediction bridge:
#   .gp_allele_dosage()  — haplo-wise array -> {0,1,2} dosage matrix
#   fit_genomic_prediction() — rrblup / gblup / bayesc fits + GenoMeta round-trip
#
# Backend-dependent blocks skip when rrBLUP / BGLR are not installed.
# =============================================================================

library(testthat)
library(haploStack)

# A fixture whose trait is a linear function of marker dosages, so the estimated
# GEBVs should correlate with the (known) true genetic values.
.gp_fixture <- local({
  meta <- sim_geno_meta(n_samples = 40L, n_markers = 60L, n_chrom = 3L,
                        block_size = 5L, seed = 1L)
  d    <- haploStack:::.gp_allele_dosage(meta)
  set.seed(7L)
  b_true <- stats::rnorm(ncol(d$M), sd = 1)
  g_true <- as.numeric(scale(d$M, center = TRUE, scale = FALSE) %*% b_true)
  y      <- g_true + stats::rnorm(nrow(d$M), sd = 0.3 * stats::sd(g_true))
  ph     <- data.frame(Sample = rownames(d$M), Yield = y,
                       stringsAsFactors = FALSE)
  list(meta = meta, M = d$M, g_true = g_true, ph = ph)
})


test_that(".gp_allele_dosage returns a {0,1,2} matrix with correct shape", {
  d <- haploStack:::.gp_allele_dosage(.gp_fixture$meta, impute = "none")
  expect_true(is.matrix(d$M))
  expect_equal(dim(d$M), c(40L, 60L))
  expect_identical(rownames(d$M), dimnames(.gp_fixture$meta$geno_arrays)[[1L]])
  expect_true(all(d$M[!is.na(d$M)] %in% c(0, 1, 2)))
  expect_identical(names(d), c("M", "allele_map", "samples"))
  expect_true(all(c("Marker", "Ref", "Alt") %in% names(d$allele_map)))
})


test_that("marker names are preserved verbatim (non-syntactic IDs)", {
  meta <- sim_geno_meta(n_samples = 8L, n_markers = 6L, seed = 3L)
  ids  <- paste0(seq_len(6L), "H-", 100L + seq_len(6L))   # e.g. "1H-101"
  dimnames(meta$geno_arrays)[[2L]] <- ids
  meta$map_info$Marker <- ids
  expect_true(is_geno_meta(meta))

  d <- haploStack:::.gp_allele_dosage(meta)
  expect_identical(colnames(d$M), ids)
  expect_identical(d$allele_map$Marker, ids)

  skip_if_not_installed("rrBLUP")
  ph  <- data.frame(Sample = dimnames(meta$geno_arrays)[[1L]],
                    Y = stats::rnorm(8L), stringsAsFactors = FALSE)
  res <- fit_genomic_prediction(meta, ph, model = "rrblup")
  eff <- res$marker_effects[[1L]]
  expect_true(all(eff$Marker %in% meta$map_info$Marker))
  expect_setequal(unique(eff$Marker), ids)
})


test_that(".gp_allele_dosage rejects non-diploid encodings", {
  meta <- .gp_fixture$meta
  meta$metadata$encoding <- "numeric"
  expect_error(haploStack:::.gp_allele_dosage(meta), "diploid")
})


test_that("rrblup fit returns attachable marker effects and GEBVs", {
  skip_if_not_installed("rrBLUP")
  res <- fit_genomic_prediction(.gp_fixture$meta, .gp_fixture$ph,
                                model = "rrblup")

  # marker_effects contract: a named list with one table per trait, each table
  # carrying the (Marker, Allele, Effect) effect_list columns.
  expect_type(res$marker_effects, "list")
  expect_named(res$marker_effects, "Yield")
  eff <- res$marker_effects$Yield
  expect_s3_class(eff, "data.frame")
  expect_identical(names(eff), c("Marker", "Allele", "Effect"))
  expect_true(is.numeric(eff$Effect))

  # GEBV phenotype table ready for attach_phenotype()
  expect_identical(names(res$gebv), c("Sample", "GEBV_Yield"))
  expect_true(is.numeric(res$gebv$GEBV_Yield))

  # attached + revalidated GenoMeta
  expect_true(is_geno_meta(res$geno_meta))
  expect_false(is.null(res$geno_meta$marker_effects))
  expect_false(is.null(res$geno_meta$phenotype))
  expect_true("GEBV_Yield" %in% names(res$geno_meta$phenotype))

  # predictive: GEBV correlates with the true genetic value
  expect_gt(stats::cor(res$gebv$GEBV_Yield, .gp_fixture$g_true), 0.5)
})


test_that("external phenotype is attached (shared individuals) beside the GEBV", {
  skip_if_not_installed("rrBLUP")
  res <- fit_genomic_prediction(.gp_fixture$meta, .gp_fixture$ph,
                                model = "rrblup")

  ph <- res$geno_meta$phenotype
  # The modelled trait AND its GEBV both persist on the GenoMeta.
  expect_true(all(c("Yield", "GEBV_Yield") %in% names(ph)))
  # Restricted to the genotype sample set (shared individuals).
  expect_setequal(ph$Sample, dimnames(.gp_fixture$meta$geno_arrays)[[1L]])
  # The attached trait values match the supplied external phenotype.
  got <- ph$Yield[match(.gp_fixture$ph$Sample, ph$Sample)]
  expect_equal(got, .gp_fixture$ph$Yield, tolerance = 1e-8)
})


test_that("external-phenotype attach keeps only individuals shared with geno", {
  skip_if_not_installed("rrBLUP")
  ph_extra <- rbind(
    .gp_fixture$ph,
    data.frame(Sample = c("ghost1", "ghost2"), Yield = c(99, -99),
               stringsAsFactors = FALSE))
  res <- suppressWarnings(
    fit_genomic_prediction(.gp_fixture$meta, ph_extra, model = "rrblup"))

  expect_false(any(c("ghost1", "ghost2") %in% res$geno_meta$phenotype$Sample))
  expect_setequal(res$geno_meta$phenotype$Sample,
                  dimnames(.gp_fixture$meta$geno_arrays)[[1L]])
})


test_that("fitted effects round-trip through compute_profile_localgebv", {
  skip_if_not_installed("rrBLUP")
  res <- fit_genomic_prediction(.gp_fixture$meta, .gp_fixture$ph,
                                model = "rrblup")
  lg <- compute_profile_localgebv(res$geno_meta)
  expect_s3_class(lg, "data.table")
  expect_true("Block" %in% names(lg))
  expect_gt(nrow(lg), 0L)
})


test_that("gblup back-solves marker effects and predicts GEBVs", {
  skip_if_not_installed("rrBLUP")
  res <- fit_genomic_prediction(.gp_fixture$meta, .gp_fixture$ph,
                                model = "gblup")
  expect_named(res$marker_effects, "Yield")
  expect_gt(nrow(res$marker_effects$Yield), 0L)
  expect_identical(names(res$gebv), c("Sample", "GEBV_Yield"))
  expect_true(is_geno_meta(res$geno_meta))
  expect_gt(stats::cor(res$gebv$GEBV_Yield, .gp_fixture$g_true), 0.5)
})


test_that("bayesc fit returns marker effects and GEBVs", {
  skip_if_not_installed("BGLR")
  res <- fit_genomic_prediction(
    .gp_fixture$meta, .gp_fixture$ph, model = "bayesc",
    nIter = 1500L, burnIn = 500L,
    saveAt = file.path(tempdir(), "bglr_test_"))
  expect_named(res$marker_effects, "Yield")
  expect_identical(names(res$marker_effects$Yield),
                   c("Marker", "Allele", "Effect"))
  expect_identical(names(res$gebv), c("Sample", "GEBV_Yield"))
  expect_true(is_geno_meta(res$geno_meta))
  expect_gt(stats::cor(res$gebv$GEBV_Yield, .gp_fixture$g_true), 0)
})


test_that("D1-7: BayesC chain files go to a private directory, never getwd()", {
  skip_if_not_installed("BGLR")
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  before <- list.files(tempdir(), pattern = "^hs_bglr_")
  res <- fit_genomic_prediction(.gp_fixture$meta, .gp_fixture$ph, model = "bayesc",
                                nIter = 300L, burnIn = 100L)
  expect_identical(list.files(wd), character(0L))                      # nothing in getwd()
  expect_identical(list.files(tempdir(), pattern = "^hs_bglr_"), before)  # private dir removed
  expect_true(is.na(res$model$ETA[[1L]]$NamefileOut))
  expect_identical(res$summary$nIter, 300L)
  expect_identical(res$summary$burnIn, 100L)
  # a caller's saveAt is passed through untouched
  keep <- file.path(withr::local_tempdir(), "chain_")
  res2 <- fit_genomic_prediction(.gp_fixture$meta, .gp_fixture$ph, model = "bayesc",
                                 nIter = 300L, burnIn = 100L, saveAt = keep)
  expect_true(length(list.files(dirname(keep))) > 0L)
  expect_false(is.na(res2$model$ETA[[1L]]$NamefileOut))
})

test_that("missing-backend and missing-phenotype paths error clearly", {
  expect_error(
    fit_genomic_prediction(.gp_fixture$meta, model = "rrblup"),
    "no phenotype")
})


# ---------------------------------------------------------------------------
# predict_gebv(): apply existing marker effects -> whole-genome GEBV in pheno
# ---------------------------------------------------------------------------

# Encoding-compatible effects: put a known effect on each marker's Alt allele
# (read from the dosage allele map), so the whole-genome GEBV must equal the
# dosage matrix times that effect vector (Ref effect = 0).
.gp_eff_fixture <- local({
  meta <- attach_phenotype(.gp_fixture$meta, .gp_fixture$ph)
  d    <- haploStack:::.gp_allele_dosage(meta)
  set.seed(11L)
  eff  <- data.frame(Marker = d$allele_map$Marker,
                     Allele = d$allele_map$Alt,
                     Effect = stats::rnorm(nrow(d$allele_map)),
                     stringsAsFactors = FALSE)
  eff  <- eff[!is.na(eff$Allele), , drop = FALSE]
  b    <- stats::setNames(numeric(ncol(d$M)), colnames(d$M))
  b[eff$Marker] <- eff$Effect
  expected      <- stats::setNames(as.numeric(d$M %*% b), rownames(d$M))
  list(meta = meta, eff_list = list(Yield = eff), expected = expected)
})


test_that("predict_gebv appends a [Trait]_GEBV column matching the dosage sum", {
  f  <- .gp_eff_fixture
  expect_message(
    gm <- predict_gebv(f$meta, marker_effects = f$eff_list, trait = "Yield"),
    "externally supplied")

  expect_true(is_geno_meta(gm))
  expect_true("Yield_GEBV" %in% names(gm$phenotype))
  expect_true("Yield" %in% names(gm$phenotype))       # original trait preserved

  got <- stats::setNames(gm$phenotype$Yield_GEBV, gm$phenotype$Sample)
  expect_equal(got[names(f$expected)], f$expected, tolerance = 1e-8)
})


test_that("predict_gebv honours a custom GEBV suffix and scores all samples", {
  f  <- .gp_eff_fixture
  gm <- predict_gebv(f$meta, marker_effects = f$eff_list, trait = "Yield",
                     gebv_suffix = "GEBV_51b")
  expect_true("Yield_GEBV_51b" %in% names(gm$phenotype))
  expect_setequal(gm$phenotype$Sample, dimnames(f$meta$geno_arrays)[[1L]])
  expect_false(anyNA(gm$phenotype$Yield_GEBV_51b))
})


test_that("predict_gebv falls back to stored effects and reports the source", {
  f  <- .gp_eff_fixture
  gm <- f$meta
  gm$marker_effects <- haploStack:::.as_effect_list(f$eff_list)
  expect_true(is_geno_meta(gm))

  expect_message(out <- predict_gebv(gm, trait = "Yield"),
                 "stored on the GenoMeta")
  expect_true("Yield_GEBV" %in% names(out$phenotype))
})


test_that("predict_gebv requires a trait when several effect tables exist", {
  f        <- .gp_eff_fixture
  two      <- f$eff_list
  two$Other <- f$eff_list$Yield
  expect_error(predict_gebv(f$meta, marker_effects = two), "multiple traits")
})


test_that("predict_gebv errors when no effects are available", {
  expect_error(predict_gebv(.gp_fixture$meta), "no marker effects")
})
