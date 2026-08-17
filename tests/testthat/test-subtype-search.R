# =============================================================================
# test-subtype-search.R
#
# QA for the subtype-search layer in R/block_subtypes.R:
#   count_subtype_occurrence()  — elemental lookup, 0 for an absent subtype
#   count_subtype_carriage()    — target-set carriage, per sample + pooled;
#                                 pinned to AGREE with looping the primitive
#   count_favourable_subtypes() — catalogue -> top subtype per block -> carriage,
#                                 incl. agreement with analyze_block_catalogue()
#   .compute_block_signatures(blocks=) — the block-restricted melt
# =============================================================================

library(testthat)
library(haploStack)

srch_meta <- function(n_samples = 12L) {
  m  <- sim_geno_meta(n_samples = n_samples, n_markers = 24L,
                      n_chrom = 3L, block_size = 4L)
  mk <- dimnames(m$geno_arrays)[[2L]]
  m$marker_effects <- list(
    Yield = sim_effect_df(mk, alleles = c("1", "4"), seed = 1L))
  m
}

# One sample's sub-GenoMeta, the intended target of functions 2 and 3.
one_sample <- function(m, s = dimnames(m$geno_arrays)[[1L]][1L]) {
  filter_geno_meta(m, keep_samples = s)
}


# ---------------------------------------------------------------------------
# .compute_block_signatures(blocks =)
# ---------------------------------------------------------------------------

test_that("the block filter returns exactly the full-genome rows for that block", {
  m    <- srch_meta()
  full <- haploStack:::.compute_block_signatures(m)
  one  <- haploStack:::.compute_block_signatures(m, blocks = "Chr1_B1")

  expect_identical(unique(one$Block), "Chr1_B1")
  ref <- full[Block == "Chr1_B1"]
  data.table::setorder(ref, Sample, Haplotype_idx)
  data.table::setorder(one, Sample, Haplotype_idx)
  expect_equal(one, ref, ignore_attr = TRUE)
})

test_that("the block filter accepts several blocks and ignores unknown ones", {
  m   <- srch_meta()
  two <- haploStack:::.compute_block_signatures(
    m, blocks = c("Chr1_B1", "Chr2_B1", "Nope"))
  expect_setequal(unique(two$Block), c("Chr1_B1", "Chr2_B1"))

  none <- haploStack:::.compute_block_signatures(m, blocks = "Nope")
  expect_identical(nrow(none), 0L)
  expect_identical(names(none),
                   c("Block", "Sample", "Haplotype_idx", "signature"))
})


# ---------------------------------------------------------------------------
# count_subtype_occurrence — the primitive
# ---------------------------------------------------------------------------

test_that("occurrence matches the enumeration it is derived from", {
  m   <- srch_meta()
  tab <- suppressMessages(count_block_subtypes(m, blocks = "Chr1_B1"))

  for (i in seq_len(nrow(tab)))
    expect_identical(
      count_subtype_occurrence(m, "Chr1_B1", tab$Subtype[i]),
      tab$Count[i],
      info = tab$Subtype[i])

  # And the copies sum to the panel's copy count.
  total <- sum(vapply(tab$Subtype,
                      function(s) count_subtype_occurrence(m, "Chr1_B1", s),
                      integer(1L)))
  expect_identical(total, 24L)                       # 12 samples x 2 copies
})

test_that("genotype mode counts distinct samples, never more than haplotype", {
  m   <- srch_meta()
  tab <- suppressMessages(
    count_block_subtypes(m, blocks = "Chr1_B1", count_by = "genotype"))
  for (i in seq_len(nrow(tab)))
    expect_identical(
      count_subtype_occurrence(m, "Chr1_B1", tab$Subtype[i],
                               count_by = "genotype"),
      tab$Count[i])

  s <- tab$Subtype[1L]
  expect_lte(count_subtype_occurrence(m, "Chr1_B1", s, count_by = "genotype"),
             count_subtype_occurrence(m, "Chr1_B1", s))
})

test_that("an absent subtype returns 0 without complaint", {
  m <- srch_meta()
  expect_identical(count_subtype_occurrence(m, "Chr1_B1", "1414141414"), 0L)
  expect_silent(count_subtype_occurrence(m, "Chr1_B1", "ACGT"))
  # Right width, not observed: still a clean 0.
  expect_identical(count_subtype_occurrence(m, "Chr1_B1", "ACGT"), 0L)
})

test_that("an absent block warns and returns 0", {
  m <- srch_meta()
  expect_warning(v <- count_subtype_occurrence(m, "Nope", "1111"), "not found")
  expect_identical(v, 0L)
})

test_that("block and subtype must be single values", {
  m <- srch_meta()
  expect_error(count_subtype_occurrence(m, c("Chr1_B1", "Chr1_B2"), "1111"),
               "single non-empty block ID")
  expect_error(count_subtype_occurrence(m, "Chr1_B1", c("1111", "4444")),
               "single non-NA signature")
})

test_that("a supplied .signatures table gives the same answer", {
  m    <- srch_meta()
  sigs <- haploStack:::.compute_block_signatures(m)
  tab  <- suppressMessages(count_block_subtypes(m, blocks = "Chr2_B1"))
  expect_identical(
    count_subtype_occurrence(m, "Chr2_B1", tab$Subtype[1L], .signatures = sigs),
    count_subtype_occurrence(m, "Chr2_B1", tab$Subtype[1L]))
})


# ---------------------------------------------------------------------------
# count_subtype_carriage — the set
# ---------------------------------------------------------------------------

test_that("carriage agrees with looping the primitive (the layering contract)", {
  m    <- srch_meta()
  tgt  <- suppressMessages(count_block_subtypes(m))[, .SD[1L], by = Block][
    , .(Block, Subtype)]

  res <- count_subtype_carriage(m, tgt)
  expect_identical(res$n_target, nrow(tgt))
  expect_identical(nrow(res$detail), nrow(tgt) * res$n_samples)

  # Pooled copies per target must equal the primitive's answer.
  pooled <- res$detail[, .(Count = sum(Count)), by = .(Block, Subtype)]
  for (i in seq_len(nrow(pooled)))
    expect_identical(
      count_subtype_occurrence(m, pooled$Block[i], pooled$Subtype[i]),
      as.integer(pooled$Count[i]),
      info = paste(pooled$Block[i], pooled$Subtype[i]))
})

test_that("uncarried targets stay in the table as explicit zeros", {
  m   <- srch_meta()
  res <- count_subtype_carriage(m, c(Chr1_B1 = "ACGT"))      # never observed
  expect_identical(res$n_target, 1L)
  expect_identical(res$n_carried, 0L)
  expect_identical(res$total, 0L)
  expect_identical(nrow(res$detail), res$n_samples)
  expect_true(all(!res$detail$Carried))
  expect_true(all(res$by_sample$nCarried == 0L))
})

test_that("a one-sample sub-GenoMeta gives the headline count directly", {
  m   <- srch_meta()
  one <- one_sample(m)
  tgt <- suppressMessages(count_block_subtypes(m))[, .SD[1L], by = Block][
    , .(Block, Subtype)]

  res <- count_subtype_carriage(one, tgt)
  expect_identical(res$n_samples, 1L)
  expect_identical(nrow(res$by_sample), 1L)
  expect_identical(res$by_sample$nCarried, res$n_carried)
  expect_equal(res$by_sample$Frac, res$n_carried / res$n_target)
  # A diploid carries 0, 1 or 2 copies of a given target.
  expect_true(all(res$detail$Count %in% 0:2))
})

test_that("genotype mode collapses each sample to carried / not carried", {
  m   <- srch_meta()
  tgt <- suppressMessages(count_block_subtypes(m))[, .SD[1L], by = Block][
    , .(Block, Subtype)]

  h <- count_subtype_carriage(m, tgt, count_by = "haplotype")
  g <- count_subtype_carriage(m, tgt, count_by = "genotype")

  expect_true(all(g$detail$Count %in% 0:1))
  expect_identical(g$n_carried, h$n_carried)      # same carriers either way
  expect_lte(g$total, h$total)
  expect_identical(g$total, g$n_carried)
  expect_identical(g$count_by, "genotype")
})

test_that("targets accept a table or a named vector, and aliases for Subtype", {
  m   <- srch_meta()
  tab <- suppressMessages(count_block_subtypes(m, blocks = "Chr1_B1"))
  s   <- tab$Subtype[1L]

  from_vec <- count_subtype_carriage(m, stats::setNames(s, "Chr1_B1"))
  from_dt  <- count_subtype_carriage(
    m, data.table::data.table(Block = "Chr1_B1", Subtype = s))
  expect_identical(from_vec$n_carried, from_dt$n_carried)

  # `signature` (catalogue) and `Signature` (allele sheet) are accepted too.
  sig_dt <- data.table::data.table(Block = "Chr1_B1", signature = s)
  cap_dt <- data.table::data.table(Block = "Chr1_B1", Signature = s)
  expect_identical(count_subtype_carriage(m, sig_dt)$n_carried,
                   from_dt$n_carried)
  expect_identical(count_subtype_carriage(m, cap_dt)$n_carried,
                   from_dt$n_carried)

  # A whole allele sheet row-set passes straight in.
  sheet <- write_block_allele_sheet(m, blocks = "Chr1_B1")
  expect_identical(count_subtype_carriage(m, sheet)$n_target, nrow(sheet))
})

test_that("bad target specifications are rejected or warned", {
  m <- srch_meta()
  expect_error(count_subtype_carriage(m, data.frame(Subtype = "1111")),
               "needs a 'Block' column")
  expect_error(count_subtype_carriage(m, data.frame(Block = "Chr1_B1")),
               "needs one of 'Subtype'")
  expect_error(count_subtype_carriage(m, c("1111", "4444")), "fully named")
  expect_error(count_subtype_carriage(m, 1:3), "must be a data.frame")

  dup <- data.table::data.table(Block = c("Chr1_B1", "Chr1_B1"),
                                Subtype = c("1111", "1111"))
  expect_warning(res <- count_subtype_carriage(m, dup), "duplicate")
  expect_identical(res$n_target, 1L)
})

test_that("an absent target block warns and counts zero", {
  m <- srch_meta()
  expect_warning(
    res <- count_subtype_carriage(m, c(Chr1_B1 = "1111", Nope = "1111")),
    "not found")
  expect_identical(res$n_target, 2L)
  expect_true(all(res$detail[Block == "Nope", Count] == 0L))
})


# ---------------------------------------------------------------------------
# count_favourable_subtypes — the pipeline
# ---------------------------------------------------------------------------

test_that("the reference is one top-effect subtype per block", {
  m   <- srch_meta()
  cat_dt <- compute_block_catalogue(m, trait = "Yield", group_by = "signature")
  res <- count_favourable_subtypes(m, cat_dt)

  ref <- res$reference
  expect_identical(data.table::uniqueN(ref$Block), nrow(ref))
  expect_setequal(ref$Block, unique(cat_dt$Block))
  expect_identical(res$fitness_fn, "max_effect")

  # Each pick really is its block's maximum effect.
  best <- cat_dt[, .(Best = max(Effect)), by = Block]
  chk  <- merge(ref, best, by = "Block")
  expect_equal(chk$Effect, chk$Best)

  # reference is sorted strongest first
  expect_false(is.unsorted(rev(ref$Effect)))
})

test_that("min_effect selects the least favourable subtype instead", {
  m   <- srch_meta()
  cat_dt <- compute_block_catalogue(m, trait = "Yield", group_by = "signature")
  res <- count_favourable_subtypes(m, cat_dt, fitness_fn = "min_effect")

  worst <- cat_dt[, .(Worst = min(Effect)), by = Block]
  chk   <- merge(res$reference, worst, by = "Block")
  expect_equal(chk$Effect, chk$Worst)
  expect_identical(res$fitness_fn, "min_effect")
  expect_false(is.unsorted(res$reference$Effect))   # ascending under min_effect
})

test_that("top_n keeps the strongest blocks only", {
  m   <- srch_meta()
  cat_dt <- compute_block_catalogue(m, trait = "Yield", group_by = "signature")
  all_res <- count_favourable_subtypes(m, cat_dt)
  cut     <- count_favourable_subtypes(m, cat_dt, top_n = 3L)

  expect_identical(nrow(cut$reference), 3L)
  expect_identical(cut$n_target, 3L)
  expect_identical(cut$reference$Block, utils::head(all_res$reference$Block, 3L))
  # top_n above the block count is harmless
  big <- count_favourable_subtypes(m, cat_dt, top_n = 10000L)
  expect_identical(big$n_target, all_res$n_target)
  expect_error(count_favourable_subtypes(m, cat_dt, top_n = 0L), "positive integer")
})

test_that("carriage matches count_subtype_carriage on the same reference", {
  m   <- srch_meta()
  cat_dt <- compute_block_catalogue(m, trait = "Yield", group_by = "signature")
  res <- count_favourable_subtypes(m, cat_dt)
  ref <- count_subtype_carriage(m, res$reference[, .(Block, Subtype)])

  expect_identical(res$n_carried, ref$n_carried)
  expect_identical(res$total, ref$total)
  expect_equal(res$by_sample, ref$by_sample)
  # detail additionally carries the target's effect
  expect_true("Effect" %in% names(res$detail))
  expect_false(anyNA(res$detail$Effect))
})

test_that("a catalogue is computed from the GenoMeta when none is given", {
  m <- srch_meta()
  from_meta <- count_favourable_subtypes(m, trait = "Yield")
  cat_dt <- compute_block_catalogue(m, trait = "Yield", group_by = "signature")
  from_cat  <- count_favourable_subtypes(m, cat_dt)

  expect_equal(from_meta$reference, from_cat$reference)
  expect_identical(from_meta$n_carried, from_cat$n_carried)
})

test_that("an analyze_block_catalogue result has its picks honoured verbatim", {
  m <- srch_meta()
  an <- suppressWarnings(analyze_block_catalogue(
    geno_meta = m, effect_list = list(Yield = m$marker_effects$Yield),
    plots = character(0), write_outputs = FALSE))

  res <- count_favourable_subtypes(m, an)
  sel <- unique(data.table::as.data.table(an$selected)[, .(Block, Subtype)],
                by = "Block")
  expect_setequal(paste(res$reference$Block, res$reference$Subtype),
                  paste(sel$Block, sel$Subtype))
})

test_that("a per-trait catalogue list is accepted", {
  m <- srch_meta()
  cat_dt <- compute_block_catalogue(m, trait = "Yield", group_by = "signature")
  as_list <- list(Yield = list(hap_catalogue = cat_dt))

  expect_equal(count_favourable_subtypes(m, as_list)$reference,
               count_favourable_subtypes(m, cat_dt)$reference)
  expect_error(count_favourable_subtypes(m, list(cat_dt)), "must be named")
})

test_that("a catalogue without an Effect column is rejected", {
  m <- srch_meta()
  bad <- compute_block_catalogue(m, trait = "Yield",
                                 group_by = "signature")[, .(Block, signature)]
  expect_error(count_favourable_subtypes(m, bad), "needs an 'Effect' column")
})

test_that("one line's favourable count is reproducible and writable", {
  m   <- srch_meta()
  cat_dt <- compute_block_catalogue(m, trait = "Yield", group_by = "signature")
  line <- one_sample(m)

  f   <- file.path(withr::local_tempdir(), "carriage.txt")
  res <- suppressMessages(
    count_favourable_subtypes(line, cat_dt, output_file = f))

  expect_identical(res$n_samples, 1L)
  expect_lte(res$n_carried, res$n_target)
  back <- data.table::fread(f)
  expect_identical(names(back),
                   c("Sample", "nTarget", "nCarried", "Total", "Frac"))
  expect_identical(back$nCarried, res$by_sample$nCarried)

  # The panel total is the sum over its lines.
  panel <- count_favourable_subtypes(m, cat_dt)
  per   <- vapply(dimnames(m$geno_arrays)[[1L]], function(s)
    count_favourable_subtypes(one_sample(m, s), cat_dt)$n_carried, integer(1L))
  expect_identical(sum(per), panel$n_carried)
})
