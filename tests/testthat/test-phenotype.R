# =============================================================================
# test-phenotype.R
#
# QA tests for the phenotype module:
#   .standardise_phenotype() via new_geno_meta(), attach_phenotype(),
#   read_phenotype(), get_phenotype(), is_geno_meta() validation,
#   filter_geno_meta() pruning, save/load round-trip, and the portable plots.
# =============================================================================

library(testthat)
library(haploStack)

# Shared GenoMeta (6 samples S1..S6) + a phenotype with one extra sample (SX,
# absent from genotypes) and one missing sample (S6 has no phenotype row).
.pheno_data <- local({
  # No block_size: phenotype handling is independent of block_info, and this
  # keeps the fixture off the data.table `:=` path in attach_block_info().
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 12L,
                        n_chrom = 2L, seed = 42L)
  samples <- dimnames(meta$geno_arrays)[[1L]]          # S1..S6
  ph <- data.frame(
    Sample = c(samples[1:5], "SX"),                    # drop S6, add unknown SX
    Pheno1 = c(10, 12, 14, 16, 18, 99),
    Pheno2 = c(1.0, 1.5, 2.0, 2.5, 3.0, 9.9),
    stringsAsFactors = FALSE
  )
  list(meta = meta, samples = samples, ph = ph)
})

# ── attach_phenotype / .standardise_phenotype ──────────────────────────────────

test_that("attach_phenotype intersects samples, drops unknowns, NA-fills missing", {
  expect_warning(
    expect_warning(
      m <- attach_phenotype(.pheno_data$meta, .pheno_data$ph),
      "absent from the genotype set"          # SX dropped
    ),
    "lack a phenotype row"                      # S6 NA-filled
  )
  expect_false(is.null(m$phenotype))
  # Reordered to the genotype sample order, unknown sample gone.
  expect_identical(m$phenotype$Sample, .pheno_data$samples)
  # S6 (position 6) NA-filled.
  expect_true(is.na(m$phenotype$Pheno1[6]))
  expect_true(is.na(m$phenotype$Pheno2[6]))
  # S1 value preserved.
  expect_equal(m$phenotype$Pheno1[1], 10)
})

test_that("attach_phenotype errors with no shared samples", {
  bad <- data.frame(Sample = c("ZZ1", "ZZ2"), Pheno1 = c(1, 2),
                    stringsAsFactors = FALSE)
  expect_error(suppressWarnings(attach_phenotype(.pheno_data$meta, bad)),
               "no samples shared")
})

test_that("attach_phenotype trait subset keeps only requested traits", {
  m <- suppressWarnings(
    attach_phenotype(.pheno_data$meta, .pheno_data$ph, traits = "Pheno1"))
  expect_identical(setdiff(names(m$phenotype), "Sample"), "Pheno1")
})

# ── is_geno_meta validation ────────────────────────────────────────────────────

test_that("is_geno_meta accepts attached phenotype and legacy objects", {
  m <- suppressWarnings(attach_phenotype(.pheno_data$meta, .pheno_data$ph))
  expect_true(is_geno_meta(m))
  # Legacy object without the phenotype field still validates.
  legacy <- .pheno_data$meta
  legacy$phenotype <- NULL
  expect_true(is_geno_meta(legacy))
})

test_that("is_geno_meta rejects a phenotype with an unknown sample", {
  m <- suppressWarnings(attach_phenotype(.pheno_data$meta, .pheno_data$ph))
  m$phenotype$Sample[1] <- "NOT_A_SAMPLE"
  expect_false(is_geno_meta(m, verbose = FALSE))
})

# ── filter_geno_meta pruning ───────────────────────────────────────────────────

test_that("filter_geno_meta prunes phenotype in lockstep with samples", {
  m    <- suppressWarnings(attach_phenotype(.pheno_data$meta, .pheno_data$ph))
  keep <- .pheno_data$samples[1:3]
  f    <- filter_geno_meta(m, keep_samples = keep)
  expect_identical(sort(f$phenotype$Sample), sort(keep))
  expect_true(is_geno_meta(f))
})

# ── save / load round-trip ─────────────────────────────────────────────────────

test_that("save/load round-trips the phenotype element", {
  m   <- suppressWarnings(attach_phenotype(.pheno_data$meta, .pheno_data$ph))
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  save_geno_meta(m, tmp)
  m2 <- load_geno_meta(tmp)
  expect_equal(m2$phenotype, m$phenotype)
})

# ── read_phenotype ─────────────────────────────────────────────────────────────

test_that("read_phenotype reads a delimited file and drops FID", {
  tmp <- tempfile(fileext = ".pheno")
  on.exit(unlink(tmp), add = TRUE)
  write.table(
    data.frame(FID = c("F1", "F2"), IID = c("S1", "S2"),
               NFNB = c(3.1, 4.2), Scald = c(0.5, 0.7)),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t")
  ph <- read_phenotype(tmp)
  expect_identical(names(ph), c("Sample", "NFNB", "Scald"))
  expect_identical(ph$Sample, c("S1", "S2"))
  expect_true(is.numeric(ph$NFNB))
})

# ── portable plots ─────────────────────────────────────────────────────────────

test_that("phenotype plots return ggplot objects from a GenoMeta", {
  m <- suppressWarnings(attach_phenotype(.pheno_data$meta, .pheno_data$ph))
  expect_s3_class(plot_phenotype_distribution(m), "ggplot")
  expect_s3_class(plot_phenotype_distribution(m, type = "box"), "ggplot")
  expect_s3_class(plot_phenotype_correlation(m), "ggplot")
})

test_that("phenotype plots are portable to a bare data.frame with Sample", {
  df <- data.frame(Sample = c("A", "B", "C", "D"),
                   T1 = c(1, 2, 3, 4), T2 = c(4, 3, 2, 1),
                   stringsAsFactors = FALSE)
  expect_s3_class(plot_phenotype_distribution(df), "ggplot")
  expect_s3_class(plot_phenotype_correlation(df), "ggplot")
  s <- phenotype_summary(df)
  expect_identical(s$Trait, c("T1", "T2"))
  expect_equal(s$Mean, c(2.5, 2.5))
})

# ── plot_phenotype_bars ────────────────────────────────────────────────────────

test_that("plot_phenotype_bars returns a ggplot in both orientations", {
  m <- suppressWarnings(attach_phenotype(.pheno_data$meta, .pheno_data$ph))
  expect_s3_class(plot_phenotype_bars(m), "ggplot")
  expect_s3_class(plot_phenotype_bars(m, orientation = "horizontal"), "ggplot")
})

test_that("plot_phenotype_bars keeps at most two traits (with warning)", {
  df <- data.frame(Sample = c("A", "B", "C"),
                   T1 = 1:3, T2 = 3:1, T3 = c(2, 2, 2),
                   stringsAsFactors = FALSE)
  expect_warning(p <- plot_phenotype_bars(df), "keeping the first two")
  expect_setequal(levels(p$data$Trait), c("T1", "T2"))
})

test_that("plot_phenotype_bars locks the sample axis order", {
  df <- data.frame(Sample = c("A", "B", "C"), T1 = 1:3,
                   stringsAsFactors = FALSE)
  p <- plot_phenotype_bars(df, sample_levels = c("C", "B", "A"))
  expect_identical(levels(p$data$Sample), c("C", "B", "A"))
})

# ── plot_meta_alignment ────────────────────────────────────────────────────────

# A minimal alignment-like reference plot: raw sample IDs on the discrete y-axis
# (mirrors a genotype-level plot_localgebv_alignment() panel) — no heavy pipeline.
.ref_align_plot <- function(samples) {
  d <- expand.grid(Sample = samples, Block = paste0("B", 1:3),
                   stringsAsFactors = FALSE)
  d$GEBV <- stats::rnorm(nrow(d))
  ggplot2::ggplot(d, ggplot2::aes(x = Block, y = Sample, fill = GEBV)) +
    ggplot2::geom_tile()
}

test_that("plot_meta_alignment merges and syncs the sample axis", {
  skip_if_not_installed("cowplot")
  samples <- c("S1", "S2", "S3", "S4")
  ref  <- .ref_align_plot(samples)
  bars <- plot_phenotype_bars(
    data.frame(Sample = rev(samples), P1 = 1:4, P2 = 4:1,
               stringsAsFactors = FALSE),
    orientation = "horizontal")
  combined <- plot_meta_alignment(list(ref, bars), rel_widths = c(3, 1))
  expect_true(inherits(combined, "ggplot") || inherits(combined, "gtable") ||
                inherits(combined, "gg"))
  # Reference order should be readable and the bars reordered to match it.
  ref_levels <- haploStack:::.discrete_axis_levels(ref, "y")
  expect_setequal(ref_levels, samples)
})

test_that("plot_meta_alignment validates its inputs", {
  skip_if_not_installed("cowplot")
  expect_error(plot_meta_alignment(list(ggplot2::ggplot())),
               "at least two")
  expect_error(plot_meta_alignment(list(ggplot2::ggplot(), 1)),
               "must be a ggplot")
})
