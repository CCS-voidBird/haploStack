# Extracted from test-090-locale.R:48

# test -------------------------------------------------------------------------
m <- sim_geno_meta(n_samples = 6L, n_markers = 8L, n_chrom = 1L, seed = 4L,
                     encoding = "nucleotide")
ids <- c("beta", "Alpha", "alpha", "Beta", "s10", "S2")
old <- dimnames(m$geno_arrays)[[1L]]
dimnames(m$geno_arrays)[[1L]] <- ids
for (cc in names(m$sample_info))          # the sample sheet names them too
    if (identical(as.character(m$sample_info[[cc]]), old)) m$sample_info[[cc]] <- ids
m$phenotype <- NULL
expect_true(is_geno_meta(m))
run <- function() {
    list(views   = geno_meta_views(m)$samples,
         palette = names(make_allele_color_map(ids, "nucleotide")),
         levels  = haploStack:::.hs_sort(unique(ids)))
  }
