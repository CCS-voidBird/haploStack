# =============================================================================
# test-half-called-views.R
#
# 0.7.6 — a half-called cell (one copy observed, one missing) keeps its
# observed allele in geno_meta_views() and in everything built on the views.
#
# Background (docs/haploStack_TAG_manuscript/06_bench_review_2026-09-17.md §3,
# INDEX discrepancy 29): the views used to rebuild hap1 / hap2 by re-splitting
# .meta_to_allele_matrix()'s collapsed genotype string, which is NA whenever
# ANY ploidy layer is NA. So (A, NA) became (NA, NA): the observed allele was
# dropped from local GEBVs, signatures, the catalogue and the stacking
# ceiling. Founders read from PED are never half-called (missing = "0 0"),
# but every block- or genomicSimulation-backend offspring population is.
#
# The rule pinned here: a half-called cell contributes its observed allele
# and nothing for the absent copy — at the views, at both local-GEBV levels,
# on both compute engines, in the block signatures, and in the simulation's
# stacking ceiling (which can never rise above the founders' under the block
# model).
# =============================================================================

library(testthat)
library(haploStack)

kernel_available <- function() is.function(haploStack:::.hs_kernel("local_gebv_block_sum_cpp"))

# Fixture: a small blocked panel in the requested encoding, with (i) a few
# both-missing cells (what a PED reader produces) and (ii) a set of
# half-called cells (what a crossing backend produces). Returns the intact
# container, the punctured one, and the index table of the punctured cells.
.hc_fixture <- function(encoding = "nucleotide", seed = 7L) {
  meta <- sim_geno_meta(n_samples = 8L, n_markers = 18L, encoding = encoding,
                        n_chrom = 2L, block_size = 3L, seed = seed)
  arr  <- meta$geno_arrays
  s    <- dimnames(arr)[[1L]]; m <- dimnames(arr)[[2L]]
  set.seed(seed)
  half <- data.frame(Sample = sample(s, 12L, replace = TRUE),
                     Marker = sample(m, 12L, replace = TRUE),
                     Copy   = sample(1:2, 12L, replace = TRUE),
                     stringsAsFactors = FALSE)
  half <- half[!duplicated(half[, c("Sample", "Marker")]), ]
  both <- data.frame(Sample = s[1:2], Marker = m[c(3L, 9L)],
                     stringsAsFactors = FALSE)
  punct <- meta
  for (i in seq_len(nrow(half)))
    punct$geno_arrays[half$Sample[i], half$Marker[i], half$Copy[i]] <- NA_character_
  for (i in seq_len(nrow(both)))
    punct$geno_arrays[both$Sample[i], both$Marker[i], ] <- NA_character_
  list(intact = meta, punct = punct, half = half, both = both,
       alleles = if (encoding == "nucleotide") c("A", "T") else c("1", "4"))
}

# Hand computation of the per-copy block sums straight from the array:
# effect(marker, allele) summed over a block's markers, NA contributing 0.
.hc_hand_hap <- function(gm, eff) {
  bmm <- as.data.frame(gm$block_marker_map_info)
  key <- setNames(eff$Effect, paste(eff$Marker, eff$Allele))
  arr <- gm$geno_arrays
  blocks <- unique(bmm$Block)
  out <- list()
  for (b in blocks) {
    mk <- bmm$Marker[bmm$Block == b]
    for (k in 1:2) {
      a <- matrix(arr[, mk, k], nrow = dim(arr)[1L])
      v <- matrix(key[paste(rep(mk, each = nrow(a)), a)], nrow(a))
      v[is.na(a)] <- 0
      v[is.na(v)] <- 0                     # allele without an effect row
      out[[paste(b, k)]] <- data.frame(Block = b, Copy = k,
                                       Sample = dimnames(arr)[[1L]],
                                       hand = rowSums(v),
                                       stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

.hc_engine_hap <- function(gm, eff, use_rcpp) {
  old <- getOption("haploStack.use_rcpp"); on.exit(options(haploStack.use_rcpp = old))
  options(haploStack.use_rcpp = use_rcpp)
  lg <- as.data.table(compute_profile_localgebv(gm, eff, level = "haplotype"))
  long <- data.table::melt(lg, id.vars = "Block", variable.name = "col",
                           value.name = "engine", variable.factor = FALSE)
  long[, Sample := sub("[.][0-9]+$", "", col)]
  long[, Copy := as.integer(sub("^.*[.]", "", col))]
  long[, .(Block, Copy, Sample, engine)]
}

for (enc in "nucleotide") {     # the one categorical vocabulary since 0.9.0 (C-5)

  test_that(paste0("views keep the observed copy of a half-called cell (", enc, ")"), {
    fx <- .hc_fixture(enc)
    v  <- geno_meta_views(fx$punct)$alleles
    arr <- fx$punct$geno_arrays
    for (i in seq_len(nrow(fx$half))) {
      s <- fx$half$Sample[i]; m <- fx$half$Marker[i]; k <- fx$half$Copy[i]
      other <- 3L - k
      expect_true(is.na(v[Sample == s & Haplotype_idx == k][[m]]))
      expect_identical(v[Sample == s & Haplotype_idx == other][[m]],
                       unname(arr[s, m, other]))
      expect_false(is.na(arr[s, m, other]))      # the fixture really is half-called
    }
    # a both-missing cell is NA on both copies
    for (i in seq_len(nrow(fx$both))) {
      s <- fx$both$Sample[i]; m <- fx$both$Marker[i]
      expect_true(all(is.na(v[Sample == s][[m]])))
    }
    # every other cell is exactly the array
    for (k in 1:2) {
      vk <- as.matrix(v[Haplotype_idx == k, dimnames(arr)[[2L]], with = FALSE])
      rownames(vk) <- v[Haplotype_idx == k, Sample]
      expect_identical(unname(vk[dimnames(arr)[[1L]], ]),
                       unname(matrix(arr[, , k], nrow = dim(arr)[1L])))
    }
  })

  test_that(paste0("views are byte-identical to the pre-0.7.6 builder when no cell is half-called (", enc, ")"), {
    fx  <- .hc_fixture(enc)
    gm  <- fx$intact
    # inject only both-missing cells: the old and new routes must agree
    for (i in seq_len(nrow(fx$both)))
      gm$geno_arrays[fx$both$Sample[i], fx$both$Marker[i], ] <- NA_character_
    new <- geno_meta_views(gm)$alleles
    # reference = the 0.7.5 route: collapse to "a/b" (NA if any copy NA), re-split
    gt  <- haploStack:::.meta_to_allele_matrix(gm)
    sp  <- dimnames(gm$geno_arrays); n <- length(sp[[1L]]); p <- length(sp[[2L]])
    h1 <- matrix(NA_character_, n, p); h2 <- h1
    for (j in seq_len(p)) {
      parts <- strsplit(as.character(gt[, j]), "/", fixed = TRUE)
      h1[, j] <- vapply(parts, function(q) if (length(q) >= 1L) q[1L] else NA_character_, "")
      h2[, j] <- vapply(parts, function(q) if (length(q) >= 2L) q[2L] else NA_character_, "")
    }
    colnames(h1) <- sp[[2L]]; colnames(h2) <- sp[[2L]]
    d1 <- as.data.table(as.data.frame(h1, check.names = FALSE)); d1[, Sample := sp[[1L]]]; d1[, Haplotype_idx := 1L]
    d2 <- as.data.table(as.data.frame(h2, check.names = FALSE)); d2[, Sample := sp[[1L]]]; d2[, Haplotype_idx := 2L]
    old <- rbindlist(list(d1, d2)); setcolorder(old, names(new))
    expect_identical(as.data.frame(new), as.data.frame(old))
  })

  test_that(paste0("local GEBVs honour the observed copy at both levels, both engines (", enc, ")"), {
    fx  <- .hc_fixture(enc)
    eff <- sim_effect_df(dimnames(fx$punct$geno_arrays)[[2L]], alleles = fx$alleles, seed = 3L)
    hand <- as.data.table(.hc_hand_hap(fx$punct, eff))
    for (rc in unique(c(FALSE, kernel_available()))) {
      eng <- .hc_engine_hap(fx$punct, eff, use_rcpp = rc)
      cmp <- merge(hand, eng, by = c("Block", "Copy", "Sample"))
      expect_equal(nrow(cmp), nrow(hand))
      expect_lt(max(abs(cmp$hand - cmp$engine)), 1e-9)
      # genotype level = the two copies summed
      old <- getOption("haploStack.use_rcpp"); options(haploStack.use_rcpp = rc)
      lg <- as.data.table(compute_profile_localgebv(fx$punct, eff, level = "genotype"))
      options(haploStack.use_rcpp = old)
      gh <- hand[, .(hand = sum(hand)), by = .(Block, Sample)]
      for (i in seq_len(nrow(gh)))
        expect_equal(lg[Block == gh$Block[i]][[gh$Sample[i]]], gh$hand[i], tolerance = 1e-9)
    }
    # and the half-called copies really carry a non-zero contribution somewhere
    # (otherwise this test would pass on the broken code by accident)
    hc <- fx$half
    contrib <- vapply(seq_len(nrow(hc)), function(i) {
      k <- 3L - hc$Copy[i]; a <- fx$punct$geno_arrays[hc$Sample[i], hc$Marker[i], k]
      eff$Effect[eff$Marker == hc$Marker[i] & eff$Allele == a]
    }, numeric(1L))
    expect_true(any(abs(contrib) > 1e-6))
  })

  test_that(paste0("block signatures keep the observed copy and fill only the absent one (", enc, ")"), {
    fx   <- .hc_fixture(enc)
    fill <- haploStack:::.hs_fill_na(NULL)
    sig0 <- haploStack:::.compute_block_signatures(fx$intact)
    sig1 <- haploStack:::.compute_block_signatures(fx$punct)
    # signature slot of each marker: block members in genomic order (as the
    # signature builder orders them)
    mm   <- geno_meta_views(fx$punct)$marker_meta[!is.na(block_id)]
    slot <- function(m) { b <- mm[marker == m, block_id]; which(mm[block_id == b, marker] == m) }
    for (i in seq_len(nrow(fx$half))) {
      s <- fx$half$Sample[i]; m <- fx$half$Marker[i]; k <- fx$half$Copy[i]
      b <- mm[marker == m, block_id]; j <- slot(m)
      # (several punctures may hit one sample x block, so compare the slot,
      # not the whole signature)
      a_keep <- sig0[Block == b & Sample == s & Haplotype_idx == 3L - k, signature]
      z_keep <- sig1[Block == b & Sample == s & Haplotype_idx == 3L - k, signature]
      a_gone <- sig0[Block == b & Sample == s & Haplotype_idx == k, signature]
      z_gone <- sig1[Block == b & Sample == s & Haplotype_idx == k, signature]
      expect_equal(nchar(z_keep), nchar(a_keep)); expect_equal(nchar(z_gone), nchar(a_gone))
      # the copy that kept its allele: that slot unchanged (and not the fill token)
      expect_identical(substr(z_keep, j, j), substr(a_keep, j, j))
      expect_false(identical(substr(z_keep, j, j), fill))
      # the absent copy: the fill token at that slot, an allele before
      expect_identical(substr(z_gone, j, j), fill)
      expect_false(identical(substr(a_gone, j, j), fill))
    }
  })
}

test_that("the stacking ceiling never rises above the founders' under the block backend", {
  # Founders with both-missing cells (as a PED gives); the block backend then
  # produces half-called offspring. Before 0.7.6 those offspring scored
  # phantom haplotype values and ideal_hap_gebv rose from cycle 1.
  meta <- sim_geno_meta(n_samples = 12L, n_markers = 24L, encoding = "nucleotide",
                        n_chrom = 2L, block_size = 4L, seed = 11L)
  set.seed(11L)
  cells <- cbind(sample(12L, 30L, TRUE), sample(24L, 30L, TRUE))
  for (i in seq_len(nrow(cells))) meta$geno_arrays[cells[i, 1L], cells[i, 2L], ] <- NA_character_
  eff <- sim_effect_df(dimnames(meta$geno_arrays)[[2L]], alleles = c("A", "T"), seed = 11L)
  res <- run_stacking_simulation(geno_meta = meta, effect_df = eff, replications = 1L,
                                 select_cycles = 3L, top_founders = 4L, num_offspring = 3L,
                                 select_top_blocks = 6L, scheme = "TS",
                                 metrics = c("gebv", "ideal_hap_gebv"),
                                 cross_backend = "block", seed = 11L,
                                 output_dir = NULL, verbose = FALSE)
  s <- as.data.table(res$summary)
  ceiling0 <- s[cycle == 0L, ideal_hap_gebv][1L]
  expect_true(all(s$ideal_hap_gebv <= ceiling0 + 1e-9))
  # and the offspring container really is half-called somewhere
  oa <- res$final_populations[[1L]]$geno_arrays
  expect_gt(sum(xor(is.na(oa[, , 1L]), is.na(oa[, , 2L]))), 0L)
})
