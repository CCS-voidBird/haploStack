# Tests for sample clustering + the clustering -> alignment portal
# (R/catalogue_cluster.R and the plot_haplotype_alignment sample_order hook).

test_that(".matching_distance is a correct pairwise matching distance", {
  # 3 units x 3 markers; row 3 is missing at marker 1.
  code <- matrix(c(1L, 1L, 2L,
                   1L, 2L, 2L,
                   NA, 1L, 2L), nrow = 3L, byrow = TRUE)
  d <- haploStack:::.matching_distance(code)

  expect_equal(dim(d), c(3L, 3L))
  expect_equal(diag(d), c(0, 0, 0))
  expect_true(isSymmetric(d))
  # d(1,2): all 3 comparable, mismatch only at marker 2 -> 1/3
  expect_equal(d[1, 2], 1 / 3)
  # d(1,3): marker 1 excluded (NA); markers 2,3 identical -> 0
  expect_equal(d[1, 3], 0)
  # d(2,3): marker 1 excluded; marker 2 differs, marker 3 same -> 1/2
  expect_equal(d[2, 3], 0.5)
})

test_that(".matching_distance excludes pairs with no comparable marker", {
  code <- matrix(c(1L, NA,
                   NA, 1L), nrow = 2L, byrow = TRUE)
  d <- haploStack:::.matching_distance(code)
  expect_true(is.na(d[1, 2]))       # no jointly-scored marker
  expect_equal(diag(d), c(0, 0))
})

test_that("copy_merge='average' is the mean of the copy-pair distances", {
  # 2 samples x 2 copies each (copies 1,2 = sample A; 3,4 = sample B).
  copy_sample <- c("A", "A", "B", "B")
  dcopy <- matrix(0, 4, 4)
  vals <- c(`13` = 0.2, `14` = 0.4, `23` = 0.6, `24` = 0.8)
  dcopy[1, 3] <- dcopy[3, 1] <- vals[["13"]]
  dcopy[1, 4] <- dcopy[4, 1] <- vals[["14"]]
  dcopy[2, 3] <- dcopy[3, 2] <- vals[["23"]]
  dcopy[2, 4] <- dcopy[4, 2] <- vals[["24"]]

  ds <- haploStack:::.merge_copy_distance(dcopy, copy_sample, method = "average")
  expect_equal(dimnames(ds)[[1]], c("A", "B"))
  expect_equal(ds["A", "B"], mean(vals))        # (0.2+0.4+0.6+0.8)/4 = 0.5
  expect_equal(unname(diag(ds)), c(0, 0))

  ds_best <- haploStack:::.merge_copy_distance(dcopy, copy_sample, method = "best")
  expect_equal(ds_best["A", "B"], min(vals))    # 0.2
})

test_that("cluster_samples returns a valid hs_sample_clustering", {
  meta <- sim_geno_meta(n_samples = 12, n_markers = 40, block_size = 5, seed = 11)
  cl   <- cluster_samples(meta, linkage = "ward.D2")

  expect_s3_class(cl, "hs_sample_clustering")
  expect_setequal(cl$order, dimnames(meta$geno_arrays)[[1L]])   # a permutation
  expect_false(anyDuplicated(cl$order) > 0)
  expect_s3_class(cl$hclust, "hclust")
  expect_equal(cl$n_samples, 12L)
  expect_null(cl$clusters)

  # Deterministic: same inputs -> same order.
  cl2 <- cluster_samples(meta, linkage = "ward.D2")
  expect_identical(cl$order, cl2$order)
})

test_that("cluster_samples k cuts into k clusters", {
  meta <- sim_geno_meta(n_samples = 15, n_markers = 40, block_size = 5, seed = 7)
  cl   <- cluster_samples(meta, k = 4)
  expect_equal(length(unique(cl$clusters)), 4L)
  expect_setequal(names(cl$clusters), dimnames(meta$geno_arrays)[[1L]])
})

test_that("cluster_samples honours copy_merge variants", {
  meta <- sim_geno_meta(n_samples = 10, n_markers = 30, block_size = 5, seed = 3)
  expect_s3_class(cluster_samples(meta, copy_merge = "average"), "hs_sample_clustering")
  expect_s3_class(cluster_samples(meta, copy_merge = "pool"),    "hs_sample_clustering")
  expect_s3_class(cluster_samples(meta, copy_merge = "best"),    "hs_sample_clustering")
})

test_that("plot_haplotype_alignment sample_order reorders rows (genotype)", {
  meta  <- sim_geno_meta(n_samples = 8, n_markers = 30, block_size = 5, seed = 5)
  cl    <- cluster_samples(meta)
  p     <- plot_haplotype_alignment(meta, method = "genotype", sample_order = cl)
  expect_equal(levels(p$data$y_label), paste0(cl$order, "_g"))

  # A bare character vector works too.
  ord   <- rev(sort(dimnames(meta$geno_arrays)[[1L]]))
  p2    <- plot_haplotype_alignment(meta, method = "genotype", sample_order = ord)
  expect_equal(levels(p2$data$y_label), paste0(ord, "_g"))
})

test_that("plot_haplotype_alignment sample_order reorders rows (haplotype)", {
  meta   <- sim_geno_meta(n_samples = 6, n_markers = 24, block_size = 4, seed = 9)
  cl     <- cluster_samples(meta)
  ploidy <- dim(meta$geno_arrays)[3L]
  p      <- plot_haplotype_alignment(meta, method = "haplotype", sample_order = cl)
  expect_equal(
    levels(p$data$y_label),
    as.character(unlist(lapply(cl$order,
                               function(s) paste0(s, "_hap", seq_len(ploidy))))))
})

test_that("default order is unchanged (alphabetical) and bad orders error", {
  meta <- sim_geno_meta(n_samples = 8, n_markers = 20, block_size = 5, seed = 1)
  p    <- plot_haplotype_alignment(meta, method = "genotype")
  expect_equal(levels(p$data$y_label),
               paste0(sort(dimnames(meta$geno_arrays)[[1L]]), "_g"))

  expect_error(
    plot_haplotype_alignment(meta, method = "genotype",
                             sample_order = c("nope1", "nope2")),
    "permutation")
})

test_that("dendrogram leaf order matches the clustering order", {
  meta <- sim_geno_meta(n_samples = 10, n_markers = 30, block_size = 5, seed = 4)
  cl   <- cluster_samples(meta, k = 3)
  seg  <- haploStack:::.hclust_segments(cl$hclust)
  expect_equal(seg$leaves$label, cl$order)

  p <- plot_cluster_dendrogram(cl)
  expect_s3_class(p, "ggplot")
})

test_that("the dendrogram's leaf axis mirrors the alignment's row expansion", {
  # 0.5.0: plot_haplotype_alignment() reserves headroom above the top row for
  # its block labels. Unless the dendrogram reserves the same, the leaves drift
  # upward against the alignment rows they are supposed to index.
  meta <- sim_geno_meta(n_samples = 10, n_markers = 30, n_chrom = 2,
                        block_size = 5, seed = 4)
  cl   <- cluster_samples(meta, k = 3)

  aln <- plot_haplotype_alignment(meta, method = "genotype",
                                  sample_order = cl)
  re  <- attr(aln, "row_expand")
  expect_gt(re[2L], 0.6)                       # there IS headroom to mirror

  leaf_rng <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1L]]$y.range
  aln_rng  <- ggplot2::ggplot_build(aln)$layout$panel_params[[1L]]$y.range
  expect_equal(leaf_rng(plot_cluster_dendrogram(cl, leaf_expand = re)),
               aln_rng)

  # The default (unmirrored) dendrogram is the one that would misalign.
  expect_false(isTRUE(all.equal(leaf_rng(plot_cluster_dendrogram(cl)), aln_rng)))

  expect_error(plot_cluster_dendrogram(cl, leaf_expand = c(1, 2, 3)),
               "leaf_expand")
})

test_that("plot_clustered_alignment composes without error", {
  skip_if_not_installed("cowplot")
  meta <- sim_geno_meta(n_samples = 8, n_markers = 24, block_size = 4, seed = 6)
  cl   <- cluster_samples(meta, k = 2)
  fig  <- plot_clustered_alignment(meta, cl, method = "genotype",
                                   show_dendrogram = TRUE)
  expect_true(inherits(fig, "gg") || inherits(fig, "ggplot") ||
              inherits(fig, "gtable") || inherits(fig, "cowplot"))

  # haplotype method returns the plain (ordered) alignment
  fig2 <- plot_clustered_alignment(meta, cl, method = "haplotype")
  expect_s3_class(fig2, "ggplot")
})
