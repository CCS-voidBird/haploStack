# =============================================================================
# test-plot-order.R
#
# QA for the shared sample row-ordering contract (R/plot_order.R, 0.5.5):
#   - .hs_resolve_sample_order(): precedence, the two defaults, strict vs
#     intersect, TRUE/FALSE, provenance
#   - use_sample_order(): validation, clearing, returning the previous value
#   - .hs_order_within(): ordered ids first, unmatched ones after, order kept
#   - REGRESSION: every figure's default output is byte-identical to 0.5.4
#   - plot_localgebv_alignment(sample_order=): rows re-ordered INSIDE a family,
#     founder boundary and offspring/founder split untouched
#   - plot_meta_alignment(check_order=): warns on a genuine disagreement only
# =============================================================================

library(testthat)
library(haploStack)

# Always leave the session's broadcast as we found it — an option leaking out of
# a test would silently re-order every later figure in the suite.
local_no_broadcast <- function(env = parent.frame()) {
  old <- getOption("haploStack.sample_order")
  withr::defer(options(haploStack.sample_order = old), envir = env)
  options(haploStack.sample_order = NULL)
  invisible(old)
}

po_meta <- function(n_samples = 6L, seed = 11L) {
  m  <- sim_geno_meta(n_samples = n_samples, n_markers = 24L,
                      n_chrom = 2L, block_size = 4L, seed = seed)
  mk <- dimnames(m$geno_arrays)[[2L]]
  m$marker_effects <- list(
    Yield = sim_effect_df(mk, alleles = c("A", "T"), seed = 3L))
  m
}

# Two-generation pedigree over the panel: s1,s2 are founders, the rest offspring.
po_meta_ped <- function(...) {
  m    <- po_meta(...)
  samp <- dimnames(m$geno_arrays)[[1L]]
  fam  <- data.frame(
    IID    = samp,
    Father = c(NA, NA, rep(samp[1], length(samp) - 2L)),
    Mother = c(NA, NA, rep(samp[2], length(samp) - 2L)),
    stringsAsFactors = FALSE)
  attach_family_history(m, fam, parent_cols = c("Father", "Mother"))
}

# The y-axis factor levels a ggplot was built with, bottom -> top.
y_levels_of <- function(p) levels(p$data$y_label)


# ── .hs_resolve_sample_order ─────────────────────────────────────────────────

test_that("nothing supplied gives the caller's own default", {
  local_no_broadcast()
  s <- c("b", "a", "c")

  r1 <- haploStack:::.hs_resolve_sample_order(NULL, s, default = "sorted")
  expect_equal(r1$order, c("a", "b", "c"))
  expect_equal(r1$source, "default")
  expect_null(r1$desc)

  # "keep" means the caller keeps its own order untouched — NULL, not a sort.
  r2 <- haploStack:::.hs_resolve_sample_order(NULL, s, default = "keep")
  expect_null(r2$order)
  expect_equal(r2$source, "default")
})

test_that("an explicit character order is used verbatim and is checked", {
  local_no_broadcast()
  s <- c("a", "b", "c")

  r <- haploStack:::.hs_resolve_sample_order(c("c", "a", "b"), s)
  expect_equal(r$order, c("c", "a", "b"))
  expect_equal(r$source, "argument")

  # strict (the pre-0.5.5 contract): not a permutation -> error
  expect_error(
    haploStack:::.hs_resolve_sample_order(c("a", "b"), s, caller = "f",
                                          subset = "strict"),
    "permutation")
  expect_error(
    haploStack:::.hs_resolve_sample_order(c("a", "b", "c", "z"), s,
                                          caller = "f", subset = "strict"),
    "Unknown")
  expect_error(
    haploStack:::.hs_resolve_sample_order(c("a", "a", "b"), s, caller = "f"),
    "duplicated")
})

test_that("intersect mode keeps what matches and appends the rest", {
  local_no_broadcast()
  s <- c("a", "b", "c", "d")

  expect_message(
    r <- haploStack:::.hs_resolve_sample_order(c("z", "c", "a"), s,
                                               caller = "f",
                                               subset = "intersect"),
    "not in the supplied order")
  expect_equal(r$order, c("c", "a", "b", "d"))   # ordered first, then the rest
  expect_setequal(r$order, s)                    # still a full permutation
})

test_that("the broadcast is used when no argument is given, and is best-effort", {
  local_no_broadcast()
  s <- c("a", "b", "c")
  use_sample_order(c("c", "b", "a"))

  r <- haploStack:::.hs_resolve_sample_order(NULL, s)
  expect_equal(r$order, c("c", "b", "a"))
  expect_equal(r$source, "option")

  # A broadcast covering samples this panel does not hold is NOT an error, even
  # in a strict caller: it is a preference for many figures, not an assertion
  # about one. This is the filtered-panel case.
  r2 <- haploStack:::.hs_resolve_sample_order(NULL, c("a", "c"),
                                              subset = "strict")
  expect_equal(r2$order, c("c", "a"))
  expect_equal(r2$source, "option")
})

test_that("FALSE ignores the broadcast, TRUE requires one", {
  local_no_broadcast()
  s <- c("a", "b", "c")

  use_sample_order(c("c", "b", "a"))
  rf <- haploStack:::.hs_resolve_sample_order(FALSE, s, default = "sorted")
  expect_equal(rf$order, c("a", "b", "c"))
  expect_equal(rf$source, "default")

  rt <- haploStack:::.hs_resolve_sample_order(TRUE, s)
  expect_equal(rt$order, c("c", "b", "a"))

  use_sample_order(NULL)
  expect_error(haploStack:::.hs_resolve_sample_order(TRUE, s, caller = "f"),
               "no order is broadcast")
})

test_that("a clustering contributes its leaf order and its provenance text", {
  local_no_broadcast()
  meta <- po_meta()
  cl   <- cluster_samples(meta, k = 2L)

  r <- haploStack:::.hs_resolve_sample_order(cl, cl$order)
  expect_equal(r$order, cl$order)
  expect_s3_class(r$clustering, "hs_sample_clustering")
  expect_match(r$desc, "ward.D2")
  expect_match(r$desc, "hamming")
  expect_match(r$desc, "k=2")
})


# ── use_sample_order ─────────────────────────────────────────────────────────

test_that("use_sample_order validates, clears, and returns the previous value", {
  local_no_broadcast()

  expect_error(use_sample_order(42), "must be an hs_sample_clustering")
  expect_error(use_sample_order(c("a", "a")), "duplicated")
  expect_error(use_sample_order(character(0L)), "empty")

  expect_message(old <- use_sample_order(c("a", "b")), "broadcasting")
  expect_null(old)                                     # nothing was set before
  expect_equal(getOption("haploStack.sample_order"), c("a", "b"))

  expect_message(prev <- use_sample_order(NULL), "cleared")
  expect_equal(prev, c("a", "b"))                      # the value it replaced
  expect_null(getOption("haploStack.sample_order"))
})

test_that("an unusable option value is ignored with a warning, not obeyed", {
  local_no_broadcast()
  options(haploStack.sample_order = 42)
  expect_warning(x <- haploStack:::.hs_sample_order(), "unusable")
  expect_null(x)
})


# ── .hs_order_within ─────────────────────────────────────────────────────────

test_that(".hs_order_within orders members and appends the unmatched", {
  expect_equal(haploStack:::.hs_order_within(c("a", "b", "c"),
                                             c("c", "x", "a")),
               c("c", "a", "b"))
  # No order -> untouched (the default path every family panel took before).
  expect_equal(haploStack:::.hs_order_within(c("b", "a"), NULL), c("b", "a"))
  # Unmatched ids keep their relative order.
  expect_equal(haploStack:::.hs_order_within(c("q", "p", "a"), c("a")),
               c("a", "q", "p"))
})


# ── Regression: the defaults must not move ───────────────────────────────────

test_that("plot_haplotype_alignment's default rows are still alphabetical", {
  local_no_broadcast()
  meta <- po_meta()
  p    <- plot_haplotype_alignment(meta, method = "genotype")
  expect_equal(y_levels_of(p), paste0(sort(dimnames(meta$geno_arrays)[[1L]]), "_g"))
  expect_equal(attr(p, "sample_order_source"), "default")
})

test_that("plot_localgebv_alignment's default rows are the family-history order", {
  local_no_broadcast()
  meta <- po_meta_ped()
  lg   <- compute_profile_localgebv(meta, level = "genotype")
  pl   <- plot_localgebv_alignment(lg, meta)
  p    <- pl[[1L]][[1L]]

  hist <- meta$metadata$family_history[[1L]]
  off  <- unique(hist$IID[hist$Label == "O"])
  fnd  <- unique(hist$IID[hist$Label == "P"])
  # Levels run founders-up-then-offspring (ggplot draws level 1 at the bottom).
  expect_equal(y_levels_of(p), rev(c(off, fnd)))
  expect_equal(attr(p, "sample_order_source"), "default")
})

test_that("a broadcast does not leak into a figure asked to ignore it", {
  local_no_broadcast()
  meta <- po_meta()
  samp <- dimnames(meta$geno_arrays)[[1L]]
  use_sample_order(rev(samp))

  p <- plot_haplotype_alignment(meta, method = "genotype", sample_order = FALSE)
  expect_equal(y_levels_of(p), paste0(sort(samp), "_g"))
})


# ── The broadcast reaching the figures ───────────────────────────────────────

test_that("a broadcast order reaches plot_haplotype_alignment", {
  local_no_broadcast()
  meta <- po_meta()
  samp <- dimnames(meta$geno_arrays)[[1L]]
  use_sample_order(rev(samp))

  p <- plot_haplotype_alignment(meta, method = "genotype")
  expect_equal(y_levels_of(p), paste0(rev(samp), "_g"))
  expect_equal(attr(p, "sample_order_source"), "option")
  expect_equal(attr(p, "sample_order"), rev(samp))
})

test_that("an order re-arranges rows inside a family, not across the split", {
  local_no_broadcast()
  meta <- po_meta_ped()
  lg   <- compute_profile_localgebv(meta, level = "genotype")

  base <- plot_localgebv_alignment(lg, meta)[[1L]][[1L]]
  samp <- dimnames(meta$geno_arrays)[[1L]]
  rev_p <- plot_localgebv_alignment(lg, meta, sample_order = rev(samp))[[1L]][[1L]]

  hist <- meta$metadata$family_history[[1L]]
  off  <- unique(hist$IID[hist$Label == "O"])
  fnd  <- unique(hist$IID[hist$Label == "P"])

  # Same membership, same offspring/founder partition, different order within.
  expect_setequal(y_levels_of(rev_p), y_levels_of(base))
  expect_equal(y_levels_of(rev_p), rev(c(rev(off), rev(fnd))))
  expect_false(identical(y_levels_of(rev_p), y_levels_of(base)))
  expect_equal(attr(rev_p, "sample_order_source"), "argument")

  # The founder boundary sits at the same row count as before: the split moved
  # no one between the two groups.
  bnd <- function(p) {
    hl <- Filter(function(l) inherits(l$geom, "GeomHline"), p$layers)
    if (length(hl) == 0L) NA_real_ else hl[[1L]]$data$yintercept[1L]
  }
  expect_equal(bnd(rev_p), bnd(base))
})

test_that("haplotype copies stay adjacent and in index order after re-ordering", {
  local_no_broadcast()
  meta <- po_meta_ped()
  lg   <- compute_profile_localgebv(meta, level = "haplotype")
  samp <- dimnames(meta$geno_arrays)[[1L]]

  p  <- plot_localgebv_alignment(lg, meta, sample_order = rev(samp))[[1L]][[1L]]
  lv <- rev(y_levels_of(p))                       # top -> bottom, as listed
  base_of <- sub("\\.\\d+$", "", lv)
  idx_of  <- as.integer(sub("^.*\\.(\\d+)$", "\\1", lv))

  # Each sample's copies form one contiguous run, numbered 1..ploidy in order.
  expect_equal(rle(base_of)$lengths, rep(2L, length(unique(base_of))))
  expect_equal(idx_of, rep(c(1L, 2L), length(unique(base_of))))
})

test_that("a clustering broadcast is labelled on the figure", {
  local_no_broadcast()
  meta <- po_meta()
  cl   <- cluster_samples(meta, k = 2L)
  use_sample_order(cl)

  p <- plot_haplotype_alignment(meta, method = "genotype")
  expect_match(p$labels$subtitle, "rows: ward.D2 / hamming, k=2", fixed = TRUE)
})

test_that("plot_favourable_blocks forwards sample_order", {
  local_no_broadcast()
  meta <- po_meta_ped()
  res  <- count_favourable_blocks(meta, trait = "Yield", plot = FALSE)
  samp <- dimnames(meta$geno_arrays)[[1L]]

  pl <- plot_favourable_blocks(res, family = meta, sample_order = rev(samp))
  expect_equal(attr(pl[[1L]][[1L]], "sample_order_source"), "argument")
})

test_that("plot_top_block_heatmap orders its sample axis first-at-top", {
  local_no_broadcast()
  meta <- po_meta()
  ana  <- suppressWarnings(
    analyze_block_catalogue(meta, meta$marker_effects))
  samp <- dimnames(meta$geno_arrays)[[1L]]

  p <- plot_top_block_heatmap(ana, sample_order = rev(samp))
  # A discrete Y axis draws level 1 at the bottom, so the levels are reversed
  # relative to the requested top-down order.
  expect_equal(levels(p$data$Sample), samp)
})

test_that("plot_phenotype_bars: sample_levels still outranks sample_order", {
  local_no_broadcast()
  ph <- data.frame(Sample = c("a", "b", "c"), Trait1 = c(1, 2, 3),
                   stringsAsFactors = FALSE)

  p <- plot_phenotype_bars(ph, sample_order = c("c", "b", "a"))
  expect_equal(levels(p$data$Sample), c("c", "b", "a"))

  # The explicit lock wins — it exists to pin this panel to a companion figure.
  p2 <- plot_phenotype_bars(ph, sample_order = c("c", "b", "a"),
                            sample_levels = c("a", "b", "c"))
  expect_equal(levels(p2$data$Sample), c("a", "b", "c"))
})


# ── plot_meta_alignment(check_order=) ────────────────────────────────────────

test_that("composing panels with different orders warns; agreeing ones do not", {
  skip_if_not_installed("cowplot")
  local_no_broadcast()
  ph <- data.frame(Sample = c("a", "b", "c"), Trait1 = c(1, 2, 3),
                   stringsAsFactors = FALSE)
  p1 <- plot_phenotype_bars(ph, sample_order = c("a", "b", "c"),
                            orientation = "horizontal")
  p2 <- plot_phenotype_bars(ph, sample_order = c("c", "b", "a"),
                            orientation = "horizontal")
  p3 <- plot_phenotype_bars(ph, sample_order = c("a", "b", "c"),
                            orientation = "horizontal")

  # legend = "none" throughout: collecting legends from a 3-sample fixture makes
  # cowplot warn about empty panels, which is noise unrelated to the order check.
  expect_warning(plot_meta_alignment(list(p1, p2), legend = "none"),
                 "different sample order")

  # Matched on the message, not on silence: cowplot's own empty-panel warnings
  # ("no non-missing arguments to min") fire on fixtures this small and would
  # make expect_silent() fail for an unrelated reason.
  mismatch_warnings <- function(...) {
    w <- testthat::capture_warnings(plot_meta_alignment(...))
    sum(grepl("different sample order", w))
  }
  expect_equal(mismatch_warnings(list(p1, p3), legend = "none"), 0L)
  expect_equal(mismatch_warnings(list(p1, p2), legend = "none",
                                 check_order = FALSE), 0L)
})
