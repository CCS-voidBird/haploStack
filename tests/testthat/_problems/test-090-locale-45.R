# Extracted from test-090-locale.R:45

# test -------------------------------------------------------------------------
m <- sim_geno_meta(n_samples = 6L, n_markers = 8L, n_chrom = 1L, seed = 4L,
                     encoding = "nucleotide")
ids <- c("beta", "Alpha", "alpha", "Beta", "s10", "S2")
dimnames(m$geno_arrays)[[1L]] <- ids
m$sample_info <- NULL
m$phenotype   <- NULL
run <- function() {
    list(views   = geno_meta_views(m)$samples,
         palette = names(make_allele_color_map(ids, "nucleotide")),
         levels  = haploStack:::.hs_sort(unique(ids)))
  }
