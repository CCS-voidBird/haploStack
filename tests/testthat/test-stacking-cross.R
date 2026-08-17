# =============================================================================
# test-stacking-cross.R
#
# Crossing backends (R/stacking_cross.R, PRD #013 issue 04):
#   .expand_cross_plan / .cross_block / .cross_backend_dispatch
#
# The block backend is pure-R and always tested; the genomicSimulation path is
# skipped when the (Suggests) package is absent. Backend parity asserts both
# return the same-shaped GenoMeta + plans.
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

# Hand-built diploid GenoMeta with one 3-marker block where each parent's two
# haplotypes are fully distinguishable ("1.." vs "4.."), so intra-block
# recombination (if any) would be directly visible.
.make_block_parents <- function() {
  samples <- c("P1", "P2")
  markers <- c("M1", "M2", "M3")
  arr <- array(NA_character_, dim = c(2L, 3L, 2L),
               dimnames = list(samples, markers, c("1", "2")))
  for (s in samples) {
    arr[s, , 1L] <- c("1", "1", "1")   # haplotype 1: all "1"
    arr[s, , 2L] <- c("4", "4", "4")   # haplotype 2: all "4"
  }
  new_geno_meta(
    geno_arrays           = arr,
    map_info              = data.frame(Marker = markers, Chrom = "1",
                                       Pos = c(1, 2, 3), stringsAsFactors = FALSE),
    block_marker_map_info = data.frame(Chrom = "1", Marker = markers,
                                       Block = "B1", stringsAsFactors = FALSE),
    metadata = list(ploidy = 2L, char_dict = c("1" = 1L, "4" = 2L),
                    n_samples = 2L, n_markers = 3L,
                    source_format = "selectiontools"))
}

# ── .expand_cross_plan ────────────────────────────────────────────────────────

test_that(".expand_cross_plan expands and names like gs_cross_plan", {
  cd  <- make_cross_design(c("A", "B", "C"), method = "all_pairs",
                           n_offspring = 2L)
  exp <- haploStack:::.expand_cross_plan(cd)
  expect_equal(nrow(exp), 6L)                       # 3 pairs x 2 offspring
  expect_setequal(names(exp), c("Parent1", "Parent2", "CrossID"))
  expect_true(all(grepl("n0[12]$", exp$CrossID)))   # "<id>n01"/"<id>n02"
  expect_equal(exp[Parent1 == "A" & Parent2 == "B", CrossID],
               c("AxBn01", "AxBn02"))
})

# ── .cross_block ──────────────────────────────────────────────────────────────

test_that(".cross_block returns a valid offspring GenoMeta of the right size", {
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 12L, n_chrom = 2L,
                        block_size = 3L, seed = 11L)
  cd   <- make_cross_design(c("S1", "S2", "S3", "S4"), method = "all_pairs",
                            n_offspring = 2L)
  set.seed(1L)
  res  <- haploStack:::.cross_block(meta, cd)

  expect_true(is_geno_meta(res$geno_meta))
  expect_equal(dim(res$geno_meta$geno_arrays)[1L], nrow(cd) * 2L)  # offspring
  expect_equal(dimnames(res$geno_meta$geno_arrays)[[2L]],
               dimnames(meta$geno_arrays)[[2L]])                   # same markers
  expect_setequal(names(res$plans),
                  c("Parent1", "Parent2", "CrossID", "F1mean"))
})

test_that(".cross_block: each offspring haplotype's alleles come from its parent", {
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 12L, n_chrom = 2L,
                        block_size = 3L, seed = 7L)
  cd   <- make_cross_design(c("S1", "S2", "S3"), method = "all_pairs",
                            n_offspring = 1L)
  set.seed(2L)
  res  <- haploStack:::.cross_block(meta, cd)

  par_arr <- meta$geno_arrays
  off_arr <- res$geno_meta$geno_arrays
  for (i in seq_len(nrow(res$plans))) {
    p1 <- res$plans$Parent1[i]; p2 <- res$plans$Parent2[i]
    for (m in dimnames(off_arr)[[2L]]) {
      expect_true(off_arr[i, m, 1L] %in% par_arr[p1, m, ])  # hap1 from Parent1
      expect_true(off_arr[i, m, 2L] %in% par_arr[p2, m, ])  # hap2 from Parent2
    }
  }
})

test_that(".cross_block computes F1mean from local GEBV when effects are given", {
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 12L, n_chrom = 2L,
                        block_size = 3L, seed = 11L)
  eff  <- sim_effect_df(dimnames(meta$geno_arrays)[[2L]], alleles = c("1", "4"))
  cd   <- make_cross_design(c("S1", "S2", "S3"), method = "all_pairs",
                            n_offspring = 2L)
  set.seed(1L)
  with_eff <- haploStack:::.cross_block(meta, cd, effect_df = eff)
  no_eff   <- haploStack:::.cross_block(meta, cd)            # no effects anywhere

  expect_false(anyNA(with_eff$plans$F1mean))                 # F1mean populated
  expect_type(with_eff$plans$F1mean, "double")
  expect_length(with_eff$plans$F1mean, nrow(with_eff$plans))
  expect_true(all(is.na(no_eff$plans$F1mean)))               # NA without effects
})

test_that(".cross_block does not recombine within a block", {
  meta <- .make_block_parents()
  cd   <- make_cross_design(c("P1", "P2"), method = "all_pairs",
                            n_offspring = 8L, avoid_selfing = TRUE)
  set.seed(3L)
  res  <- haploStack:::.cross_block(meta, cd)
  off  <- res$geno_meta$geno_arrays

  # Within the single block, each offspring haplotype must be homogeneous
  # (entirely "1" or entirely "4") — a crossover inside the block would mix them.
  for (i in seq_len(dim(off)[1L])) {
    expect_length(unique(off[i, , 1L]), 1L)
    expect_length(unique(off[i, , 2L]), 1L)
  }
})

test_that(".cross_block rejects unknown parents and non-diploid ploidy", {
  meta <- sim_geno_meta(n_samples = 4L, n_markers = 8L, block_size = 2L, seed = 5L)
  bad  <- data.table(Parent1 = "S1", Parent2 = "GHOST", n_offspring = 1L)
  expect_error(haploStack:::.cross_block(meta, bad), "not in the pool")

  meta3 <- meta; meta3$metadata$ploidy <- 3L
  cd    <- make_cross_design(c("S1", "S2"), n_offspring = 1L)
  expect_error(haploStack:::.cross_block(meta3, cd), "diploid")
})

# ── .cross_backend_dispatch ───────────────────────────────────────────────────

test_that("dispatch 'block' returns geno_meta + plans", {
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 12L, block_size = 3L, seed = 9L)
  cd   <- make_cross_design(c("S1", "S2", "S3"), method = "all_pairs",
                            n_offspring = 2L)
  set.seed(4L)
  res  <- haploStack:::.cross_backend_dispatch("block", meta, cd)
  expect_setequal(names(res), c("geno_meta", "plans", "session"))
  expect_null(res$session)                       # block backend keeps no session
  expect_true(is_geno_meta(res$geno_meta))
  expect_setequal(names(res$plans),
                  c("Parent1", "Parent2", "CrossID", "F1mean"))
})

test_that("dispatch 'genomicSimulation' matches the block backend's shape", {
  skip_if_not_installed("genomicSimulation")
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 12L, block_size = 3L, seed = 13L)
  eff  <- sim_effect_df(dimnames(meta$geno_arrays)[[2L]], alleles = c("1", "4"))
  cd   <- make_cross_design(c("S1", "S2", "S3"), method = "all_pairs",
                            n_offspring = 2L)

  set.seed(6L)
  rb <- haploStack:::.cross_backend_dispatch("block", meta, cd, effect_df = eff)
  rg <- haploStack:::.cross_backend_dispatch("genomicSimulation", meta, cd,
                                             effect_df = eff)

  # Same-shaped output: both are GenoMeta + identical plan columns over the same
  # marker set and offspring count.
  expect_true(is_geno_meta(rg$geno_meta))
  expect_equal(names(rg$plans), names(rb$plans))
  expect_equal(dim(rg$geno_meta$geno_arrays)[1L],
               dim(rb$geno_meta$geno_arrays)[1L])
  expect_setequal(dimnames(rg$geno_meta$geno_arrays)[[2L]],
                  dimnames(rb$geno_meta$geno_arrays)[[2L]])
  expect_false(anyNA(rg$plans$F1mean))   # gS supplies native breeding values
  expect_false(anyNA(rb$plans$F1mean))   # block computes F1mean from local GEBV
})

# ── .stk_rename_generation — offspring FAM parity (flag 1) ─────────────────────

test_that(".stk_rename_generation writes offspring FAM with real parents", {
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 12L, block_size = 3L, seed = 21L)
  cd   <- make_cross_design(c("S1", "S2", "S3"), method = "all_pairs",
                            n_offspring = 2L)
  set.seed(3L)
  cross_res <- haploStack:::.cross_block(meta, cd, "F")
  opts <- list(id_fmt_sample = "C%02dF%02d_%01d", id_fmt_cross = "C%02dF%02d")

  gen <- haploStack:::.stk_rename_generation(cross_res, cycle = 1L, opts)
  si  <- gen$off$sample_info

  # progeny FAM: real founder parents (not dummy "0"), a non-trivial Family id,
  # and rows aligned to the genotype sample order.
  expect_true(all(si$Parent1 %in% c("S1", "S2", "S3")))
  expect_true(all(si$Parent2 %in% c("S1", "S2", "S3")))
  expect_false(any(si$Family == "0"))
  expect_equal(si$Sample, dimnames(gen$off$geno_arrays)[[1L]])
  # full-sibs (same parent pair) share a Family; distinct pairs differ
  fam_by_pair <- unique(si[, c("Family", "Parent1", "Parent2")])
  expect_equal(nrow(fam_by_pair), length(unique(paste(si$Parent1, si$Parent2))))
})
