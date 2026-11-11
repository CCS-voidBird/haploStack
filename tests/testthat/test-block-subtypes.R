# =============================================================================
# test-block-subtypes.R
#
# QA for R/block_subtypes.R (count_block_subtypes):
#   - the fixed six-column schema and Block / -Count ordering
#   - count_by = "haplotype" vs "genotype", checked against an independent
#     recomputation from geno_meta_views()
#   - Effect join agreement with compute_block_catalogue(group_by="signature")
#   - agreement with write_block_allele_sheet() on clean AND missing data
#   - positional signatures: complementary gaps stay distinct (discrepancy 17),
#     fixed signature width, configurable fill_na
#   - block filtering, absent blocks, empty results, file output
# =============================================================================

library(testthat)
library(haploStack)

SCHEMA <- c("Trait", "Block", "Subtype", "Count", "Freq", "Effect")

sub_meta <- function() {
  m  <- sim_geno_meta(n_samples = 10L, n_markers = 20L,
                      n_chrom = 3L, block_size = 4L)
  mk <- dimnames(m$geno_arrays)[[2L]]
  m$marker_effects <- list(
    Yield = sim_effect_df(mk, alleles = c("A", "T"), seed = 1L))
  m
}

# Member markers of a block, in the genomic order the signature uses.
block_markers <- function(m, blk) {
  mm <- geno_meta_views(m)$marker_meta
  as.character(mm[block_id == blk, marker])
}


# ---------------------------------------------------------------------------
# Schema and counting
# ---------------------------------------------------------------------------

test_that("haplotype counting yields the fixed schema and covers every copy", {
  m   <- sub_meta()
  blk <- "Chr1_B1"
  out <- suppressMessages(count_block_subtypes(m, blocks = blk))

  expect_true(data.table::is.data.table(out))
  expect_identical(names(out), SCHEMA)
  expect_identical(unique(out$Block), blk)
  # No missing data in the fixture: every copy of every sample is counted.
  expect_identical(sum(out$Count), 20L)
  expect_equal(sum(out$Freq), 1)
  expect_true(all(is.na(out$Trait)))
  expect_true(all(is.na(out$Effect)))
})

test_that("genotype counting matches an independent per-sample recomputation", {
  m   <- sub_meta()
  blk <- "Chr1_B1"
  mk  <- block_markers(m, blk)

  al  <- geno_meta_views(m)$alleles
  sub <- al[, c("Sample", mk), with = FALSE]
  sig <- do.call(paste0, sub[, mk, with = FALSE])
  ref <- data.table::data.table(Sample = sub$Sample, Subtype = sig)[
    , .(Count = data.table::uniqueN(Sample)), by = Subtype]
  data.table::setorder(ref, Subtype)

  out <- suppressMessages(
    count_block_subtypes(m, blocks = blk, count_by = "genotype"))
  data.table::setorder(out, Subtype)

  expect_identical(out$Subtype, ref$Subtype)
  expect_identical(out$Count, as.integer(ref$Count))
  expect_equal(sum(out$Freq), 1)
})

test_that("genotype counts never exceed haplotype counts on the same subtype", {
  m <- sub_meta()
  h <- suppressMessages(count_block_subtypes(m, count_by = "haplotype"))
  g <- suppressMessages(count_block_subtypes(m, count_by = "genotype"))
  j <- merge(h[, .(Block, Subtype, H = Count)],
             g[, .(Block, Subtype, G = Count)], by = c("Block", "Subtype"))
  expect_identical(nrow(j), nrow(h))
  expect_true(all(j$G <= j$H))
})

test_that("results are ordered by block then descending count", {
  m   <- sub_meta()
  out <- suppressMessages(
    count_block_subtypes(m, blocks = c("Chr1_B1", "Chr2_B1")))
  expect_setequal(unique(out$Block), c("Chr1_B1", "Chr2_B1"))
  expect_false(is.unsorted(out$Block))
  expect_true(all(out[, all(diff(Count) <= 0L), by = Block]$V1))
})

test_that("blocks = NULL covers every block and Freq sums to 1 in each", {
  m    <- sub_meta()
  out  <- suppressMessages(count_block_subtypes(m))
  sigs <- haploStack:::.compute_block_signatures(m)

  expect_identical(nrow(out),
                   nrow(unique(sigs[, .(Block, signature)])))
  expect_setequal(unique(out$Block), unique(sigs$Block))
  expect_true(all(abs(out[, sum(Freq), by = Block]$V1 - 1) < 1e-12))
})


# ---------------------------------------------------------------------------
# Block selection
# ---------------------------------------------------------------------------

test_that("absent blocks warn, present ones still return", {
  m <- sub_meta()
  expect_warning(
    out <- count_block_subtypes(m, blocks = c("Chr1_B1", "Nope")),
    "not found")
  expect_identical(unique(out$Block), "Chr1_B1")
})

test_that("an all-absent block set returns the empty schema, not an error", {
  m <- sub_meta()
  expect_warning(out <- count_block_subtypes(m, blocks = "Nope"), "not found")
  expect_identical(names(out), SCHEMA)
  expect_identical(nrow(out), 0L)
})

test_that("a GenoMeta with no block assignments returns the empty schema", {
  m <- sim_geno_meta(block_size = NULL)
  out <- count_block_subtypes(m)
  expect_identical(names(out), SCHEMA)
  expect_identical(nrow(out), 0L)
})

test_that("factor block IDs are coerced", {
  m <- sub_meta()
  a <- suppressMessages(count_block_subtypes(m, blocks = "Chr1_B1"))
  b <- suppressMessages(count_block_subtypes(m, blocks = factor("Chr1_B1")))
  expect_equal(a, b)
})


# ---------------------------------------------------------------------------
# Effects
# ---------------------------------------------------------------------------

test_that("trait = NULL attaches no effects but points at what is available", {
  m <- sub_meta()
  expect_message(count_block_subtypes(m, blocks = "Chr1_B1"), "Yield")

  m$marker_effects <- NULL
  expect_silent(out <- count_block_subtypes(m, blocks = "Chr1_B1"))
  expect_true(all(is.na(out$Effect)))
})

test_that("Effect matches compute_block_catalogue on the same signature key", {
  m   <- sub_meta()
  out <- count_block_subtypes(m, blocks = c("Chr1_B1", "Chr3_B2"),
                              trait = "Yield")

  expect_identical(unique(out$Trait), "Yield")
  expect_false(any(is.na(out$Effect)))

  cat_dt <- compute_block_catalogue(m, trait = "Yield",
                                    group_by = "signature")
  ref <- cat_dt[Block %in% c("Chr1_B1", "Chr3_B2"),
                .(Block, Subtype = signature, Effect, Count)]

  j <- merge(out[, .(Block, Subtype, Effect, Count)], ref,
             by = c("Block", "Subtype"), suffixes = c("", ".ref"))
  expect_identical(nrow(j), nrow(out))
  expect_equal(j$Effect, j$Effect.ref)
  # Counting is the catalogue's own, not a parallel implementation.
  expect_identical(j$Count, j$Count.ref)
})

test_that("trait errors propagate from the effect resolver", {
  m <- sub_meta()
  expect_error(count_block_subtypes(m, blocks = "Chr1_B1", trait = "Protein"),
               "not in geno_meta\\$marker_effects")

  m$marker_effects <- NULL
  expect_error(count_block_subtypes(m, blocks = "Chr1_B1", trait = "Yield"),
               "fit_genomic_prediction")
})


# ---------------------------------------------------------------------------
# Agreement with, and documented divergence from, the allele sheet
# ---------------------------------------------------------------------------

test_that("counts agree with write_block_allele_sheet when nothing is missing", {
  m   <- sub_meta()
  blk <- "Chr1_B1"

  for (mode in c("haplotype", "genotype")) {
    ours <- suppressMessages(
      count_block_subtypes(m, blocks = blk, count_by = mode))
    sheet <- write_block_allele_sheet(m, blocks = blk, count_by = mode)

    j <- merge(ours[, .(Subtype, Count)],
               sheet[, .(Subtype = Signature, Sheet = Count)], by = "Subtype")
    expect_identical(nrow(j), nrow(ours), info = mode)
    expect_identical(j$Count, as.integer(j$Sheet), info = mode)
  }
})

test_that("counts agree with write_block_allele_sheet on MISSING data too", {
  m   <- sub_meta()
  blk <- "Chr1_B1"
  mk  <- block_markers(m, blk)
  m$geno_arrays["S1", mk, 1L] <- NA_character_

  ours  <- suppressMessages(count_block_subtypes(m, blocks = blk))
  sheet <- write_block_allele_sheet(m, blocks = blk, count_by = "haplotype")

  # Both producers fill the gap with the same placeholder, so no copy is
  # dropped and the two Signature/Subtype strings are identical.
  expect_identical(sum(sheet$Count), 20L)
  expect_identical(sum(ours$Count), 20L)
  expect_setequal(ours$Subtype, sheet$Signature)

  j <- merge(ours[, .(Subtype, Count)],
             sheet[, .(Subtype = Signature, Sheet = Count)], by = "Subtype")
  expect_identical(nrow(j), nrow(ours))
  expect_identical(j$Count, as.integer(j$Sheet))

  # S1 is still represented — as an all-fill subtype, not as a dropped copy.
  fill <- haploStack:::.hs_fill_na()
  expect_true(strrep(fill, length(mk)) %in% ours$Subtype)
})


# ---------------------------------------------------------------------------
# Positional signatures (discrepancy 17)
# ---------------------------------------------------------------------------

test_that("complementary gaps stay distinct instead of collapsing", {
  m   <- sub_meta()
  blk <- "Chr1_B1"
  mk  <- block_markers(m, blk)

  # Force two samples onto the same observed alleles but with the gap in
  # different positions: (a, -, c, d) vs (a, b, -, d) with b == c.
  for (s in c("S1", "S2")) m$geno_arrays[s, mk, ] <- "1"
  m$geno_arrays["S1", mk[2], ] <- NA_character_
  m$geno_arrays["S2", mk[3], ] <- NA_character_

  sigs <- haploStack:::.compute_block_signatures(m)[Block == blk]
  s1 <- unique(sigs[Sample == "S1", signature])
  s2 <- unique(sigs[Sample == "S2", signature])

  expect_length(s1, 1L); expect_length(s2, 1L)
  expect_false(identical(s1, s2))          # the collision this fixes
  expect_identical(s1, "1-11")
  expect_identical(s2, "11-1")

  # Both are full block width, and every copy is still counted.
  out <- suppressMessages(count_block_subtypes(m, blocks = blk))
  expect_true(all(nchar(out$Subtype) == length(mk)))
  expect_identical(sum(out$Count), 20L)
})

test_that("every signature is exactly the block's marker count wide", {
  m <- sub_meta()
  mk_n <- vapply(unique(m$block_marker_map_info$Block),
                 function(b) length(block_markers(m, b)), integer(1L))
  # Punch holes at random so most blocks carry a gap somewhere.
  set.seed(3)
  am <- dimnames(m$geno_arrays)[[2L]]
  for (k in sample(am, 8L)) m$geno_arrays[sample(dimnames(m$geno_arrays)[[1L]], 3L), k, ] <- NA

  out <- suppressMessages(count_block_subtypes(m))
  w   <- out[, .(w = unique(nchar(Subtype))), by = Block]
  expect_identical(nrow(w), data.table::uniqueN(out$Block))   # one width per block
  expect_identical(w$w, unname(mk_n[w$Block]))
  expect_true(all(out[, sum(Count), by = Block]$V1 == 20L))
})

test_that("fill_na is configurable per call and by option, and cannot be empty", {
  m   <- sub_meta()
  blk <- "Chr1_B1"
  mk  <- block_markers(m, blk)
  m$geno_arrays["S1", mk[1], ] <- NA_character_

  per_call <- suppressMessages(
    count_block_subtypes(m, blocks = blk, fill_na = "?"))
  expect_true(any(grepl("?", per_call$Subtype, fixed = TRUE)))
  expect_false(any(grepl("-", per_call$Subtype, fixed = TRUE)))

  withr::with_options(list(haploStack.fill_na = "N"), {
    opt <- suppressMessages(count_block_subtypes(m, blocks = blk))
    expect_true(any(grepl("N", opt$Subtype, fixed = TRUE)))
  })

  expect_error(count_block_subtypes(m, blocks = blk, fill_na = ""), "empty string")
  expect_error(count_block_subtypes(m, blocks = blk, fill_na = c("-", "?")),
               "single non-NA character")
  expect_error(write_block_allele_sheet(m, blocks = blk, fill_na = ""),
               "empty string")
})


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

test_that("a pre-built signature table gives the same answer", {
  m    <- sub_meta()
  sigs <- haploStack:::.compute_block_signatures(m)
  a <- suppressMessages(count_block_subtypes(m, blocks = "Chr1_B1"))
  b <- suppressMessages(
    count_block_subtypes(m, blocks = "Chr1_B1", .signatures = sigs))
  expect_equal(a, b)
})

test_that("output_file round-trips the returned table", {
  m <- sub_meta()
  f <- file.path(withr::local_tempdir(), "subtypes.txt")
  out  <- suppressMessages(count_block_subtypes(m, blocks = "Chr1_B1",
                                                output_file = f))
  back <- data.table::fread(f)
  expect_identical(names(back), SCHEMA)
  # fread type-guesses an all-digit signature as integer, as it does for
  # write_block_allele_sheet()'s Signature column; compare on the text.
  expect_identical(as.character(back$Subtype), out$Subtype)
  expect_identical(back$Count, out$Count)
})
