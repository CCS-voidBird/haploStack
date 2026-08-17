# =============================================================================
# test-favourable-blocks.R
#
# QA for R/stacking_evaluate.R (count_favourable_blocks / plot_favourable_blocks):
#   - the fixed nine-column schema and its ordering
#   - the counting identity nPositive + nNegative + nZero == nBlock * ploidy
#   - counts against an independent recomputation from the wide matrix
#   - GEBV equals the genotype-level local GEBV column sum (copies sum to it)
#   - prefer = "negative" moves nPreferred only
#   - pedigree columns: NA without a family history, matched with one, one row
#     per sample even for a parent appearing in two generations
#   - the dummy matrix (shape, 1/0 only) and the reused figure
#   - fill_na = NA_real_ makes nBlock per-sample
#   - plot_localgebv_alignment()'s default legend title is unchanged (FR4)
# =============================================================================

library(testthat)
library(haploStack)

SCHEMA <- c("GID", "FID", "IID", "Label", "nBlock", "GEBV",
            "nPositive", "nNegative", "nPreferred")

# A panel with block assignments and one stored trait, no pedigree.
fb_meta <- function(n_samples = 8L, seed = 7L) {
  m  <- sim_geno_meta(n_samples = n_samples, n_markers = 24L,
                      n_chrom = 3L, block_size = 4L, seed = seed)
  mk <- dimnames(m$geno_arrays)[[2L]]
  m$marker_effects <- list(
    Yield = sim_effect_df(mk, alleles = c("1", "4"), seed = 1L))
  m
}

# The same panel with a two-generation pedigree: s1,s2 -> s3,s4 -> s5,s6 ...
# so s3/s4 are offspring in gen_1 and parents in gen_2.
fb_meta_ped <- function(...) {
  m    <- fb_meta(...)
  samp <- dimnames(m$geno_arrays)[[1L]]
  n    <- length(samp)
  fam  <- data.frame(
    IID    = samp,
    Father = c(NA, NA, rep(samp[1], 2L), rep(samp[3], n - 4L)),
    Mother = c(NA, NA, rep(samp[2], 2L), rep(samp[4], n - 4L)),
    stringsAsFactors = FALSE)
  attach_family_history(m, fam, parent_cols = c("Father", "Mother"))
}


# ---------------------------------------------------------------------------
# Schema, list shape and the counting identity
# ---------------------------------------------------------------------------

test_that("the result is the documented list and the table the fixed schema", {
  m   <- fb_meta()
  res <- suppressWarnings(count_favourable_blocks(m, trait = "Yield"))

  expect_identical(names(res),
                   c("table", "plots", "dummy", "prefer", "level"))
  expect_identical(res$prefer, "positive")
  expect_identical(res$level, "haplotype")
  expect_true(data.table::is.data.table(res$table))
  expect_identical(names(res$table), SCHEMA)

  # One row per sample, keyed on IID.
  samp <- dimnames(m$geno_arrays)[[1L]]
  expect_setequal(res$table$IID, samp)
  expect_identical(nrow(res$table), length(samp))
})

test_that("nPositive + nNegative + nZero == nBlock * ploidy", {
  m   <- fb_meta()
  res <- suppressWarnings(count_favourable_blocks(m, trait = "Yield"))
  tab <- res$table
  ploidy <- dim(m$geno_arrays)[3L]

  # With the default fill_na = 0 every block holds a value for every sample.
  n_blocks <- data.table::uniqueN(res$dummy$Block)
  expect_true(all(tab$nBlock == n_blocks))

  n_zero <- tab$nBlock * ploidy - tab$nPositive - tab$nNegative
  expect_true(all(n_zero >= 0L))
  expect_true(all(tab$nPositive + tab$nNegative + n_zero ==
                    tab$nBlock * ploidy))
})

test_that("counts match an independent recomputation from the wide matrix", {
  m  <- fb_meta()
  lg <- compute_profile_localgebv(m, trait = "Yield", level = "haplotype",
                                  center_by_freq = TRUE)
  lg <- data.table::as.data.table(lg)

  cols <- setdiff(names(lg), "Block")
  ref  <- data.table::rbindlist(lapply(unique(sub("\\.[0-9]+$", "", cols)),
    function(s) {
      v <- unlist(lg[, cols[sub("\\.[0-9]+$", "", cols) == s], with = FALSE],
                  use.names = FALSE)
      data.table::data.table(IID = s, GEBV = sum(v),
                             nPositive = sum(v > 0), nNegative = sum(v < 0))
    }))

  got <- suppressWarnings(count_favourable_blocks(m, trait = "Yield"))$table
  j   <- merge(got, ref, by = "IID", suffixes = c("", ".ref"))

  expect_equal(j$GEBV, j$GEBV.ref)
  expect_identical(j$nPositive, as.integer(j$nPositive.ref))
  expect_identical(j$nNegative, as.integer(j$nNegative.ref))
})

test_that("GEBV equals the genotype-level local GEBV column sum", {
  # Haplotype copies of a block sum to that block's genotype-level value, so the
  # per-sample total must agree across the two levels.
  m  <- fb_meta()
  gt <- data.table::as.data.table(
    compute_profile_localgebv(m, trait = "Yield", level = "genotype",
                              center_by_freq = TRUE))
  ref <- vapply(setdiff(names(gt), "Block"),
                function(s) sum(gt[[s]]), numeric(1L))

  got <- suppressWarnings(count_favourable_blocks(m, trait = "Yield"))$table
  expect_equal(got$GEBV, unname(ref[got$IID]))
})


# ---------------------------------------------------------------------------
# prefer
# ---------------------------------------------------------------------------

test_that("prefer selects nPreferred and changes nothing else", {
  m   <- fb_meta()
  pos <- suppressWarnings(
    count_favourable_blocks(m, trait = "Yield", prefer = "positive"))$table
  neg <- suppressWarnings(
    count_favourable_blocks(m, trait = "Yield", prefer = "negative"))$table

  data.table::setorder(pos, IID)
  data.table::setorder(neg, IID)

  expect_identical(pos$nPreferred, pos$nPositive)
  expect_identical(neg$nPreferred, neg$nNegative)
  # Both directions are always reported, and only nPreferred differs.
  expect_identical(pos$nPositive, neg$nPositive)
  expect_identical(pos$nNegative, neg$nNegative)
  expect_equal(pos$GEBV, neg$GEBV)
})

test_that("the dummy matrix follows prefer", {
  m  <- fb_meta()
  lg <- data.table::as.data.table(
    compute_profile_localgebv(m, trait = "Yield", level = "haplotype",
                              center_by_freq = TRUE))
  cols <- setdiff(names(lg), "Block")

  for (p in c("positive", "negative")) {
    d <- suppressWarnings(
      count_favourable_blocks(m, trait = "Yield", prefer = p))$dummy
    for (cc in cols) {
      want <- if (p == "positive") as.integer(lg[[cc]] > 0)
              else                 as.integer(lg[[cc]] < 0)
      expect_identical(d[[cc]], want, info = paste(p, cc))
    }
  }
})


# ---------------------------------------------------------------------------
# Pedigree columns
# ---------------------------------------------------------------------------

test_that("pedigree columns are NA without a family history", {
  m   <- fb_meta()
  tab <- suppressWarnings(count_favourable_blocks(m, trait = "Yield"))$table

  expect_true(all(is.na(tab$GID)))
  expect_true(all(is.na(tab$FID)))
  expect_true(all(is.na(tab$Label)))
})

test_that("pedigree columns match the family history, one row per sample", {
  m    <- fb_meta_ped()
  samp <- dimnames(m$geno_arrays)[[1L]]
  res  <- count_favourable_blocks(m, trait = "Yield", plot = FALSE)
  tab  <- res$table

  # Every sample once, even the parents that appear in two generations.
  expect_identical(nrow(tab), length(samp))
  expect_identical(anyDuplicated(tab$IID), 0L)

  hist <- m$metadata$family_history
  expect_true(length(hist) >= 2L)

  # Founders are gen_1 parents; the middle generation is offspring in gen_1
  # even though it also parents gen_2 (oldest birth family wins).
  expect_identical(tab[IID == samp[1], GID], "gen_1")
  expect_identical(tab[IID == samp[1], Label], "P")
  expect_identical(tab[IID == samp[3], GID], "gen_1")
  expect_identical(tab[IID == samp[3], Label], "O")
  expect_identical(tab[IID == samp[5], GID], "gen_2")
  expect_identical(tab[IID == samp[5], Label], "O")

  # Every (GID, FID, IID, Label) row exists in the history it came from.
  for (i in seq_len(nrow(tab))) {
    g <- data.table::as.data.table(hist[[tab$GID[i]]])
    expect_true(nrow(g[IID == tab$IID[i] & FID == tab$FID[i] &
                         Label == tab$Label[i]]) == 1L)
  }
})

test_that("a sample outside the history keeps NA pedigree columns", {
  m    <- fb_meta_ped()
  samp <- dimnames(m$geno_arrays)[[1L]]
  # Drop the last sample from the history, keep it in the genotypes.
  h <- lapply(m$metadata$family_history,
              function(g) data.table::as.data.table(g)[IID != samp[length(samp)]])
  m$metadata$family_history <- h

  tab <- count_favourable_blocks(m, trait = "Yield", plot = FALSE)$table
  expect_true(is.na(tab[IID == samp[length(samp)], GID]))
  expect_false(any(is.na(tab[IID == samp[1], GID])))
})

test_that("rows are ordered by GID, FID, then descending nPreferred", {
  m   <- fb_meta_ped(n_samples = 10L)
  tab <- count_favourable_blocks(m, trait = "Yield", plot = FALSE)$table

  # Each (GID, FID) group is one contiguous run of rows ...
  tab[, grp := paste(GID, FID, sep = "\r")]
  expect_identical(data.table::uniqueN(tab$grp),
                   length(unique(data.table::rleid(tab$grp))))
  # ... and within a run, nPreferred descends.
  for (g in unique(tab$grp)) {
    v <- tab[grp == g, nPreferred]
    expect_identical(v, sort(v, decreasing = TRUE), info = g)
  }
})


# ---------------------------------------------------------------------------
# Dummy matrix and the reused figure
# ---------------------------------------------------------------------------

test_that("the dummy matrix has the local-GEBV matrix's shape and is 1/0 only", {
  m  <- fb_meta()
  lg <- data.table::as.data.table(
    compute_profile_localgebv(m, trait = "Yield", level = "haplotype",
                              center_by_freq = TRUE))
  d  <- suppressWarnings(count_favourable_blocks(m, trait = "Yield"))$dummy

  expect_identical(names(d), names(lg))
  expect_identical(d$Block, lg$Block)          # genomic block order preserved
  vals <- unlist(d[, setdiff(names(d), "Block"), with = FALSE],
                 use.names = FALSE)
  expect_true(all(vals %in% c(0L, 1L)))
})

test_that("nPreferred equals the dummy matrix's per-sample column sum", {
  m   <- fb_meta()
  res <- suppressWarnings(count_favourable_blocks(m, trait = "Yield"))
  d   <- res$dummy
  cols <- setdiff(names(d), "Block")
  base <- sub("\\.[0-9]+$", "", cols)

  ref <- vapply(unique(base), function(s) {
    sum(unlist(d[, cols[base == s], with = FALSE], use.names = FALSE))
  }, numeric(1L))

  tab <- res$table
  expect_identical(tab$nPreferred, as.integer(unname(ref[tab$IID])))
})

test_that("plot = TRUE returns the nested plot list on a pedigreed object", {
  m   <- fb_meta_ped()
  res <- count_favourable_blocks(m, trait = "Yield", plot = TRUE)

  expect_type(res$plots, "list")
  expect_true(all(grepl("^gen_", names(res$plots))))
  expect_s3_class(res$plots[["gen_1"]][[1L]], "ggplot")
  # Drawn from the dummy: every fill value is 0 or 1.
  expect_true(all(res$plots[["gen_1"]][[1L]]$data$GEBV %in% c(0, 1)))
})

test_that("plot without a family history warns and yields NULL, table intact", {
  m   <- fb_meta()
  res <- NULL
  expect_warning(res <- count_favourable_blocks(m, trait = "Yield",
                                                plot = TRUE),
                 "no family history")
  expect_null(res$plots)
  expect_identical(names(res$table), SCHEMA)
})

test_that("plot_favourable_blocks accepts a result list or a bare dummy", {
  m   <- fb_meta_ped()
  res <- count_favourable_blocks(m, trait = "Yield", plot = FALSE)
  expect_null(res$plots)

  from_res <- plot_favourable_blocks(res, family = m)
  from_dt  <- plot_favourable_blocks(res$dummy, family = m)
  expect_s3_class(from_res[["gen_1"]][[1L]], "ggplot")
  expect_s3_class(from_dt[["gen_1"]][[1L]], "ggplot")

  # The legend is labelled for what is actually drawn, not "Local GEBV".
  lab <- from_res[["gen_1"]][[1L]]$scales$scales[[1L]]$name
  expect_true(grepl("Favourable", lab))
  expect_true(grepl("positive", lab))
})

test_that("plot_favourable_blocks requires a family", {
  m   <- fb_meta()
  res <- suppressWarnings(count_favourable_blocks(m, trait = "Yield",
                                                 plot = FALSE))
  expect_error(plot_favourable_blocks(res), "'family' is required")
})


# ---------------------------------------------------------------------------
# fill_na, subsets, file output
# ---------------------------------------------------------------------------

test_that("fill_na = NA_real_ makes nBlock the per-sample informative count", {
  # A block row only exists at all when SOME copy has effect information in it,
  # so the difference fill_na makes is per (block, copy): punch out one sample's
  # alleles across one block's markers and that sample alone loses the block.
  m    <- fb_meta()
  samp <- dimnames(m$geno_arrays)[[1L]][1L]
  mm   <- geno_meta_views(m)$marker_meta
  hole <- as.character(mm[block_id == "Chr1_B1", marker])
  m$geno_arrays[samp, hole, ] <- NA_character_

  full <- suppressWarnings(
    count_favourable_blocks(m, trait = "Yield", fill_na = 0))$table
  part <- suppressWarnings(
    count_favourable_blocks(m, trait = "Yield", fill_na = NA_real_))$table
  data.table::setorder(full, IID)
  data.table::setorder(part, IID)

  # The holed sample loses exactly that block; everyone else is unchanged.
  expect_identical(part[IID == samp, nBlock], full[IID == samp, nBlock] - 1L)
  expect_identical(part[IID != samp, nBlock], full[IID != samp, nBlock])
  # An empty cell contributes nothing either way, so the totals must agree.
  expect_equal(full$GEBV, part$GEBV)
  expect_identical(full$nPositive, part$nPositive)
  expect_identical(full$nNegative, part$nNegative)
})

test_that("blocks / samples subsets are honoured", {
  m    <- fb_meta()
  blks <- c("Chr1_B1", "Chr1_B2")
  samp <- dimnames(m$geno_arrays)[[1L]][1:3]

  res <- suppressWarnings(count_favourable_blocks(
    m, trait = "Yield", blocks = blks, samples = samp))

  expect_setequal(res$dummy$Block, blks)
  expect_setequal(res$table$IID, samp)
  expect_true(all(res$table$nBlock == length(blks)))
})

test_that("output_file writes the table and the list is still returned", {
  m   <- fb_meta()
  f   <- withr::local_tempfile(fileext = ".txt")
  res <- suppressWarnings(suppressMessages(
    count_favourable_blocks(m, trait = "Yield", output_file = f)))

  expect_true(file.exists(f))
  back <- data.table::fread(f)
  expect_identical(names(back), SCHEMA)
  expect_identical(nrow(back), nrow(res$table))
  expect_identical(back$nPreferred, res$table$nPreferred)
})

test_that("an invalid prefer or GenoMeta errors", {
  m <- fb_meta()
  expect_error(count_favourable_blocks(list(), trait = "Yield"),
               "not a valid GenoMeta")
  expect_error(count_favourable_blocks(m, trait = "Yield", prefer = "high"),
               "should be one of")
})


# ---------------------------------------------------------------------------
# FR4 — the additive fill_name must not change existing output
# ---------------------------------------------------------------------------

test_that("plot_localgebv_alignment keeps its default legend title", {
  m  <- fb_meta_ped()
  lg <- compute_profile_localgebv(m, trait = "Yield")
  p  <- plot_localgebv_alignment(lg, m)[["gen_1"]][[1L]]

  expect_identical(p$scales$scales[[1L]]$name, "Local GEBV")
})
