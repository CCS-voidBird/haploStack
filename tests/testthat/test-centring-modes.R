# =============================================================================
# test-centring-modes.R
#
# CEN (0.6.0) — one centring vocabulary: none / locus_mean / allele_freq.
# Spec: dev/PRD_missing_allele_normalisation.md, *Second workstream*.
# Gates V12-V16.
#
# The three modes are different operations, not settings of one:
#   none         contrib(a) = e_a
#   locus_mean   contrib(a) = e_a - mean_i     (location shift; variance kept)
#   allele_freq  contrib(a) = e_a * (1 - p_a)  (reweighting; variance moves)
# =============================================================================

library(testthat)
library(haploStack)

# -- fixture ------------------------------------------------------------------

# One block of `n_markers` biallelic markers, no missing data, so every mode is
# computable on complete cells and the algebra below has no NA to hide behind.
cen_fixture <- function(n_samples = 40L, n_markers = 6L, seed = 4L) {
  set.seed(seed)
  markers <- sprintf("M%d", seq_len(n_markers))
  samples <- sprintf("S%02d", seq_len(n_samples))

  a1 <- matrix(sample(c("A", "T"), n_samples * n_markers, replace = TRUE),
               nrow = n_samples, dimnames = list(samples, markers))
  a2 <- matrix(sample(c("A", "T"), n_samples * n_markers, replace = TRUE),
               nrow = n_samples, dimnames = list(samples, markers))
  cells <- matrix(paste(a1, a2, sep = "/"), nrow = n_samples,
                  dimnames = list(samples, markers))

  map  <- data.frame(Chrom = "1", Marker = markers,
                     Pos = seq_len(n_markers) * 1e6, stringsAsFactors = FALSE)
  bmap <- data.frame(Chrom = "1", Marker = markers, Block = "B1",
                     stringsAsFactors = FALSE)

  meta <- suppressMessages(build_geno_meta(
    cells, map, blocks = bmap, ploidy = 2L, split_char = "/",
    encoding = "nucleotide", source_format = "cen_fixture"))

  eff <- data.frame(
    Marker = rep(markers, each = 2L),
    Allele = rep(c("A", "T"), times = n_markers),
    Effect = round(rnorm(2L * n_markers), 6),
    stringsAsFactors = FALSE)

  list(meta = meta, eff = eff, markers = markers, samples = samples,
       a1 = a1, a2 = a2)
}

gebv <- function(fx, center, level = "genotype") {
  suppressMessages(suppressWarnings(compute_profile_localgebv(
    fx$meta, effect_df = fx$eff, level = level, center = center)))
}

as_vec <- function(dt) {
  d <- as.data.frame(dt)
  unlist(d[, setdiff(names(d), "Block"), drop = FALSE], use.names = FALSE)
}

# -- V12: locus_mean is the mode that agrees with HapSelect, value for value --

test_that("V12: locus_mean equals the substitution-effect form HapSelect computes", {
  fx <- cen_fixture()

  # HapSelect's mean_adjust = TRUE computes (d - mean_d) * beta per marker, with
  # d the counted-allele dosage and beta the substitution effect e_A - e_B. The
  # PRD's identity says locus-mean centring produces exactly that:
  #
  #   d*e_A + (ploidy - d)*e_B - ploidy*mean_i  ==  (d - ploidy*p_A)*(e_A - e_B)
  #
  # Computing the right-hand side here from first principles gives a check that
  # needs no second package installed.
  d  <- (fx$a1 == "A") + (fx$a2 == "A")           # copies of A per (sample, marker)
  pA <- colMeans(d) / 2                            # population frequency of A
  eA <- setNames(fx$eff$Effect[fx$eff$Allele == "A"], fx$markers)
  eT <- setNames(fx$eff$Effect[fx$eff$Allele == "T"], fx$markers)
  beta <- eA - eT

  expected <- as.vector(
    (d - matrix(2 * pA, nrow(d), ncol(d), byrow = TRUE)) %*% beta)
  names(expected) <- rownames(d)

  got <- gebv(fx, "locus_mean")
  got_v <- as_vec(got)
  names(got_v) <- setdiff(names(as.data.frame(got)), "Block")

  expect_equal(max(abs(got_v[names(expected)] - expected)), 0, tolerance = 1e-12)
})

# -- V13: locus_mean is a pure location shift -------------------------------

test_that("V13: locus_mean leaves block variance identical to uncentred", {
  fx <- cen_fixture()
  v_none <- var(as_vec(gebv(fx, "none")))
  v_locus <- var(as_vec(gebv(fx, "locus_mean")))
  expect_equal(v_locus, v_none, tolerance = 1e-12)
})

test_that("V13: locus_mean differs from uncentred only by a per-block constant", {
  fx <- cen_fixture()
  d <- as_vec(gebv(fx, "none")) - as_vec(gebv(fx, "locus_mean"))
  expect_equal(max(d) - min(d), 0, tolerance = 1e-12)
})

test_that("allele_freq DOES move block variance — the difference that matters", {
  fx <- cen_fixture()
  v_none <- var(as_vec(gebv(fx, "none")))
  v_freq <- var(as_vec(gebv(fx, "allele_freq")))
  expect_false(isTRUE(all.equal(v_freq, v_none, tolerance = 1e-6)))
})

# -- V16: parameterisation invariance ---------------------------------------

test_that("V16: adding a constant to both alleles leaves locus_mean unchanged", {
  fx  <- cen_fixture()
  fx2 <- fx
  # Shift BOTH alleles of marker M1 by the same constant: the allele contrast is
  # untouched, so any mode that claims parameterisation-invariance must not move.
  fx2$eff$Effect[fx2$eff$Marker == "M1"] <-
    fx2$eff$Effect[fx2$eff$Marker == "M1"] + 5

  expect_equal(as_vec(gebv(fx2, "locus_mean")), as_vec(gebv(fx, "locus_mean")),
               tolerance = 1e-12)

  # allele_freq — the mode 0.5.3 chose FOR that invariance — does move.
  expect_false(isTRUE(all.equal(as_vec(gebv(fx2, "allele_freq")),
                                as_vec(gebv(fx, "allele_freq")),
                                tolerance = 1e-6)))
})

# -- V15: the deprecation shim ----------------------------------------------

test_that("V15: center_by_freq maps to the new vocabulary", {
  expect_equal(
    suppressMessages(haploStack:::.resolve_center(NULL, TRUE, "f")),
    "allele_freq")
  expect_equal(
    suppressMessages(haploStack:::.resolve_center(NULL, FALSE, "f")),
    "none")
  expect_equal(haploStack:::.resolve_center(NULL, NULL, "f", default = "none"),
               "none")
  expect_equal(
    haploStack:::.resolve_center(NULL, NULL, "f", default = "allele_freq"),
    "allele_freq")
})

test_that("V15: supplying both center and center_by_freq is an error", {
  fx <- cen_fixture()
  expect_error(
    compute_profile_localgebv(fx$meta, effect_df = fx$eff,
                              center = "none", center_by_freq = TRUE),
    "not both")
})

test_that("V15: the deprecated logical still produces the old numbers", {
  fx <- cen_fixture()
  old_true <- suppressMessages(suppressWarnings(compute_profile_localgebv(
    fx$meta, effect_df = fx$eff, center_by_freq = TRUE)))
  expect_equal(as_vec(old_true), as_vec(gebv(fx, "allele_freq")))

  old_false <- suppressMessages(suppressWarnings(compute_profile_localgebv(
    fx$meta, effect_df = fx$eff, center_by_freq = FALSE)))
  expect_equal(as_vec(old_false), as_vec(gebv(fx, "none")))
})

test_that("V15: a logical passed as `center` is refused with the fix named", {
  fx <- cen_fixture()
  expect_error(
    compute_profile_localgebv(fx$meta, effect_df = fx$eff, center = TRUE),
    "mode name")
  expect_error(
    compute_profile_localgebv(fx$meta, effect_df = fx$eff, center = "middle"),
    "should be one of")
})

# -- the two engines must agree on all three modes ---------------------------

test_that("the Rcpp and pure-R paths agree on every centring mode", {
  fx <- cen_fixture()
  for (mode in c("none", "locus_mean", "allele_freq")) {
    a <- withr::with_options(list(haploStack.use_rcpp = FALSE),
                             as_vec(gebv(fx, mode)))
    b <- withr::with_options(list(haploStack.use_rcpp = TRUE),
                             as_vec(gebv(fx, mode)))
    expect_equal(a, b, tolerance = 1e-12,
                 info = paste("centring mode:", mode))
  }
})

test_that("haplotype level agrees with genotype level under every mode", {
  fx <- cen_fixture()
  for (mode in c("none", "locus_mean", "allele_freq")) {
    g <- as.data.frame(gebv(fx, mode, level = "genotype"))
    h <- as.data.frame(gebv(fx, mode, level = "haplotype"))
    # genotype = sum over the ploidy copies (verified in the 2026-08-21
    # investigation as max |genotype - (hap1 + hap2)| = 0).
    hs <- rowSums(h[, setdiff(names(h), "Block"), drop = FALSE])
    gs <- rowSums(g[, setdiff(names(g), "Block"), drop = FALSE])
    expect_equal(unname(hs), unname(gs), tolerance = 1e-12,
                 info = paste("centring mode:", mode))
  }
})

# -- C2: the default flips, and the old numbers stay reachable ---------------

test_that("C2: run_haplotype_catalogue defaults to no centring", {
  expect_equal(formals(run_haplotype_catalogue)$center_by_freq, NULL)
  fx  <- cen_fixture()
  out <- withr::local_tempdir()

  res_default <- suppressMessages(suppressWarnings(run_haplotype_catalogue(
    fx$meta, effect_list = list(T1 = fx$eff), output_dir = out,
    save_profile = FALSE)))
  res_none <- suppressMessages(suppressWarnings(run_haplotype_catalogue(
    fx$meta, effect_list = list(T1 = fx$eff), output_dir = out,
    save_profile = FALSE, center = "none")))

  expect_equal(res_default, res_none)
})

test_that("C2: center = 'allele_freq' reproduces the 0.5.5 catalogue numbers", {
  fx  <- cen_fixture()
  out <- withr::local_tempdir()

  # 0.5.5's default was center_by_freq = TRUE; the shim must land on the same
  # arithmetic as the named mode.
  old <- suppressMessages(suppressWarnings(run_haplotype_catalogue(
    fx$meta, effect_list = list(T1 = fx$eff), output_dir = out,
    save_profile = FALSE, center_by_freq = TRUE)))
  named <- suppressMessages(suppressWarnings(run_haplotype_catalogue(
    fx$meta, effect_list = list(T1 = fx$eff), output_dir = out,
    save_profile = FALSE, center = "allele_freq")))

  expect_equal(old, named)

  # And it is genuinely a different answer from the new default.
  new_default <- suppressMessages(suppressWarnings(run_haplotype_catalogue(
    fx$meta, effect_list = list(T1 = fx$eff), output_dir = out,
    save_profile = FALSE)))
  expect_false(isTRUE(all.equal(old, new_default)))
})

test_that("C2: compute_profile_localgebv's default is unchanged", {
  # It already defaulted to no centring; the flip makes the two entry points
  # agree rather than introducing a new answer anywhere.
  fx <- cen_fixture()
  expect_equal(as_vec(suppressMessages(suppressWarnings(
    compute_profile_localgebv(fx$meta, effect_df = fx$eff)))),
    as_vec(gebv(fx, "none")))
})
