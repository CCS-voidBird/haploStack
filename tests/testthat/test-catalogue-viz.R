# =============================================================================
# test-catalogue-viz.R
#
# QA for plot_haplotype_alignment() (R/catalogue_viz.R), in two parts:
#
#  1. Parity (0.4.9) — the vectorised / Rcpp-backed long-format build must
#     reproduce, tile for tile, the original as.data.frame.table() +
#     per-(Sample,marker) group-by logic, under both the compiled kernel and the
#     pure-R fallback (options(haploStack.use_rcpp)).
#  2. Chromosome sub-panels (0.5.0) — chromosomes are facet columns in natural
#     order, blocks are annotated per panel, and the single-panel path is still
#     reachable via chrom_panels = FALSE.
# =============================================================================

library(testthat)
library(haploStack)

toggle <- function(rcpp, expr)
  withr::with_options(list(haploStack.use_rcpp = rcpp), expr)

# Reference: the ORIGINAL long-format logic, returning (x_pos, allele) sorted so
# the per-marker allele multiset can be compared against the plot's data.
.ref_align_data <- function(meta, method) {
  arr     <- meta$geno_arrays
  markers <- dimnames(arr)[[2L]]
  mm <- as.data.table(meta$map_info)[Marker %in% markers,
          .(marker = as.character(Marker), chromosome = as.character(Chrom),
            position = as.numeric(Pos))][order(chromosome, position)]
  mx <- data.table(marker = mm$marker, x_pos = seq_along(mm$marker))

  ld <- as.data.table(as.data.frame.table(arr, responseName = "allele",
                                          stringsAsFactors = FALSE))
  data.table::setnames(ld, c("Sample", "marker", "Haplotype_idx", "allele"))
  ld[, allele := as.character(allele)]
  ld <- ld[marker %in% mm$marker]
  if (method == "genotype")
    ld <- ld[, .(allele = paste0(sort(allele), collapse = "")),
             by = .(Sample, marker)]
  ld[mx, on = "marker", x_pos := i.x_pos]
  ld[order(x_pos, allele), .(x_pos, allele)]
}

.plot_align_data <- function(meta, method) {
  p <- plot_haplotype_alignment(meta, method = method)
  as.data.table(p$data)[order(x_pos, allele),
                        .(x_pos, allele = as.character(allele))]
}

test_that("plot_haplotype_alignment data matches the original logic", {
  meta <- sim_geno_meta(n_samples = 12L, n_markers = 40L,
                        n_chrom = 3L, block_size = 4L)
  for (m in c("genotype", "haplotype")) {
    ref <- .ref_align_data(meta, m)
    for (u in c(TRUE, FALSE)) {
      got <- toggle(u, .plot_align_data(meta, m))
      expect_identical(got$x_pos,  ref$x_pos,  info = paste(m, "use_rcpp =", u))
      expect_identical(got$allele, ref$allele, info = paste(m, "use_rcpp =", u))
    }
  }
})

test_that("genotype and haplotype engines agree with each other", {
  meta <- sim_geno_meta(n_samples = 15L, n_markers = 30L, block_size = 5L)
  for (m in c("genotype", "haplotype")) {
    a <- toggle(FALSE, .plot_align_data(meta, m))
    b <- toggle(TRUE,  .plot_align_data(meta, m))
    expect_identical(a, b, info = m)
  }
})


# =============================================================================
# 0.5.0 — chromosome sub-panels
# =============================================================================

# Panel table of a built plot: one row per facet panel, carrying the facet
# variable. Unfaceted plots build to a single panel.
.panels <- function(p) ggplot2::ggplot_build(p)$layout$layout

# Data of the first layer drawn by `geom`, e.g. the block-boundary vlines.
.layer_data_for <- function(p, geom) {
  hit <- vapply(p$layers, function(l) inherits(l$geom, geom), logical(1L))
  if (!any(hit)) return(NULL)
  p$layers[[which(hit)[1L]]]$data
}

test_that(".natural_chrom_levels sorts numeric-led names before the rest", {
  nat <- haploStack:::.natural_chrom_levels
  expect_identical(nat(c("10H", "2H", "1H", "7H")),
                   c("1H", "2H", "7H", "10H"))
  expect_identical(nat(c("2", "10", "1")), c("1", "2", "10"))
  # Non-numeric names sort alphabetically, after every numeric-led one.
  expect_identical(nat(c("Un", "2H", "chrX", "1H")),
                   c("1H", "2H", "Un", "chrX"))
  expect_identical(nat(c(NA, "1H", "1H")), "1H")
  expect_identical(nat(character(0L)), character(0L))
})

test_that("chromosomes become one sub-panel each, in natural order", {
  meta <- sim_geno_meta(n_samples = 8L, n_markers = 48L, n_chrom = 4L,
                        block_size = 4L)
  p   <- plot_haplotype_alignment(meta, method = "genotype")
  lay <- .panels(p)

  expect_true("chromosome" %in% names(lay))
  expect_identical(as.character(lay$chromosome), c("1", "2", "3", "4"))
  expect_identical(nrow(lay), 4L)          # one column of panels per chromosome
  expect_true(all(lay$ROW == 1L))
  expect_identical(lay$COL, 1:4)

  # Every plotted tile carries the chromosome its marker sits on.
  d <- as.data.table(p$data)
  mp <- as.data.table(meta$map_info)[, .(marker = as.character(Marker),
                                         chrom  = as.character(Chrom))]
  d[mp, on = "marker", chrom_ref := i.chrom]
  expect_identical(as.character(d$chromosome), d$chrom_ref)
})

test_that("chrom_panels = FALSE keeps the pre-0.5.0 single-panel figure", {
  meta <- sim_geno_meta(n_samples = 8L, n_markers = 48L, n_chrom = 4L,
                        block_size = 4L)
  p <- plot_haplotype_alignment(meta, method = "genotype",
                                chrom_panels = FALSE)
  expect_identical(nrow(.panels(p)), 1L)
  expect_identical(p$labels$subtitle, "Chr 1,2,3,4")

  # ...and the faceted figure names the chromosomes in its strips instead.
  pf <- plot_haplotype_alignment(meta, method = "genotype")
  expect_false(identical(pf$labels$subtitle, "Chr 1,2,3,4"))
})

test_that("a single-chromosome GenoMeta is never faceted", {
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 12L, n_chrom = 1L,
                        block_size = 3L)
  p <- plot_haplotype_alignment(meta, method = "genotype",
                                chrom_panels = TRUE)
  expect_identical(nrow(.panels(p)), 1L)
  expect_identical(p$labels$subtitle, "Chr 1")
})

test_that("block boundaries and labels are routed to their own panel", {
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 48L, n_chrom = 4L,
                        block_size = 4L)

  vl_f <- .layer_data_for(plot_haplotype_alignment(meta, method = "genotype"),
                          "GeomVline")
  vl_1 <- .layer_data_for(
    plot_haplotype_alignment(meta, method = "genotype", chrom_panels = FALSE),
    "GeomVline")

  # Both carry the facet variable so ggplot can place each line.
  expect_true(all(c("chromosome", "xintercept") %in% names(vl_f)))
  # The single-panel figure keeps every block change, including the three that
  # straddle a chromosome break; the faceted one drops those (the panel gap
  # already says it).
  expect_identical(nrow(vl_1) - nrow(vl_f), 3L)
  expect_false(any(vl_f$cross))
  expect_identical(sum(vl_1$cross), 3L)

  lab <- .layer_data_for(plot_haplotype_alignment(meta, method = "genotype"),
                         "GeomText")
  expect_true(all(c("chromosome", "x_center", "block_id") %in% names(lab)))
  # Each block label sits on the chromosome its markers belong to.
  bmm <- as.data.table(meta$block_marker_map_info)
  mp  <- as.data.table(meta$map_info)[, .(Marker = as.character(Marker),
                                          Chrom  = as.character(Chrom))]
  ref <- unique(merge(bmm[, .(Marker = as.character(Marker),
                              Block  = as.character(Block))],
                      mp, by = "Marker")[, .(Block, Chrom)])
  chk <- merge(as.data.table(lab)[, .(Block = as.character(block_id),
                                      got   = as.character(chromosome))],
               ref, by = "Block")
  expect_identical(chk$got, chk$Chrom)
})

test_that("block labels reserve headroom, and companions can mirror it", {
  meta <- sim_geno_meta(n_samples = 20L, n_markers = 24L, n_chrom = 2L,
                        block_size = 4L)

  p  <- plot_haplotype_alignment(meta, method = "genotype")
  re <- attr(p, "row_expand")
  expect_identical(length(re), 2L)
  expect_gt(re[2L], re[1L])                     # headroom only at the top
  expect_equal(re[2L], ceiling(0.10 * 20) + 0.6)

  # No labels requested -> no headroom, so an unannotated alignment keeps the
  # plain discrete expansion on both ends.
  p0 <- plot_haplotype_alignment(meta, method = "genotype",
                                 show_block_labels = FALSE)
  expect_equal(attr(p0, "row_expand"), c(0.6, 0.6))

  # ...and a GenoMeta with no block map has nothing to label either.
  no_blk <- sim_geno_meta(n_samples = 20L, n_markers = 24L, n_chrom = 2L)
  expect_equal(attr(plot_haplotype_alignment(no_blk, method = "genotype"),
                    "row_expand"), c(0.6, 0.6))
})
