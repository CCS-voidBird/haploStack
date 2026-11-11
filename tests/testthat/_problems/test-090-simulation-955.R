# Extracted from test-090-simulation.R:955

# prequel ----------------------------------------------------------------------
library(testthat)
library(haploStack)
library(data.table)
.e_fixture <- function(seed = 5L, n_samples = 12L, n_markers = 24L) {
  m <- sim_geno_meta(n_samples = n_samples, n_markers = n_markers, n_chrom = 2L,
                     block_size = 3L, seed = seed)
  e <- sim_effect_df(dimnames(m$geno_arrays)[[2L]], alleles = c("A", "T"))
  list(meta = m, eff = e,
       lg = compute_profile_localgebv(m, e, level = "haplotype"))
}
.e_stream_untouched <- function(expr) {
  set.seed(314159)
  pre <- .Random.seed
  force(expr)
  identical(pre, .Random.seed)
}

# test -------------------------------------------------------------------------
skip_if_not_installed("rrBLUP")
m  <- sim_geno_meta(n_samples = 20L, n_markers = 30L, n_chrom = 2L, block_size = 3L, seed = 9L)
ph <- data.frame(Sample = dimnames(m$geno_arrays)[[1L]],
                   Yield = haploStack:::.with_seed(1L, stats::rnorm(20L)))
fit <- suppressMessages(fit_genomic_prediction(m, phenotype = ph, model = "rrblup",
                                                 trait = "Yield"))
s <- fit$summary
expect_named(s, c("model", "n_train", "var_u", "var_e", "h2", "seed", "backend",
                    "backend_version"))
