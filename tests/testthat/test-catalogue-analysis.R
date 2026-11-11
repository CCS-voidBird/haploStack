# =============================================================================
# test-catalogue-analysis.R
#
# QA for the cross-trait block-catalogue analysis (2026-07-07 revision):
#   analyze_block_catalogue() — within-trait standard score, ONE representative
#   subtype per block chosen by the highest summed value across traits
#   (value_source = raw S or normalized N), signed trait_weights, per-trait top
#   blocks, data-source resolution, gated writes, and the plot layer that fills by
#   the analysis-chosen value_source.
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

# Shared fixture: a small GenoMeta with blocks + two trait effect tables on
# deliberately different scales (so within-trait normalisation matters).
.analysis_fixture <- function(seed = 5L) {
  meta <- sim_geno_meta(n_samples = 12L, n_markers = 18L, n_chrom = 3L,
                        block_size = 3L, seed = seed)
  mk   <- dimnames(meta$geno_arrays)[[2L]]
  set.seed(seed)
  eff1 <- data.frame(Marker = rep(mk, each = 2L),
                     Allele = rep(c("A", "T"), length(mk)),
                     Effect = round(rnorm(2L * length(mk), 0, 0.5), 3),
                     stringsAsFactors = FALSE)
  eff2 <- data.frame(Marker = rep(mk, each = 2L),
                     Allele = rep(c("A", "T"), length(mk)),
                     Effect = round(rnorm(2L * length(mk), 0, 3.0), 3),
                     stringsAsFactors = FALSE)
  list(meta = meta, eff1 = eff1, eff2 = eff2)
}

# Expected one-subtype-per-block pick: argmax over subtypes of the summed value.
.expected_pick <- function(res, w = NULL, col = "S") {
  sub <- data.table::copy(res$subtypes)
  trs <- unique(sub$Trait)
  if (is.null(w)) w <- stats::setNames(rep(1, length(trs)), trs)
  sub[, V := get(col)]
  sc  <- sub[, .(sc = sum(w[Trait] * V)), by = .(Block, Subtype)]
  sc[, .SD[which.max(sc)], by = Block][, .(Block, wantSub = Subtype)]
}

# ── return structure ──────────────────────────────────────────────────────────

test_that("analyze_block_catalogue returns the documented components", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                  plots = character(0)))
  expect_setequal(names(res),
    c("subtypes", "selected", "S", "N", "top_blocks", "gebv", "plots",
      "value_source", "normalize",
      "figures", "notes", "provenance"))          # 0.9.0: hs_block_analysis (D3.4)
  expect_s3_class(res, "hs_block_analysis")
  expect_s3_class(res$selected, "data.table")
  expect_true(all(c("Block", "Trait", "S", "is_top") %in% names(res$S)))
  expect_true(all(c("Block", "Trait", "N", "is_top") %in% names(res$N)))
  # identification columns lead the analysis tables (Trait/Chrom/Pos/Block first)
  expect_identical(names(res$S)[1:2], c("Trait", "Chrom"))
  expect_true(which(names(res$S) == "Block") < which(names(res$S) == "S"))
  expect_type(res$plots, "list")
  # removed with the 2026-07-07 revision:
  expect_false(any(c("consensus", "block_metric") %in% names(res)))
})

# ── one subtype per block (shared across traits) ──────────────────────────────

test_that("each block has exactly one selected subtype, identical across traits", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2), plots = character(0)))
  # one row per block in `selected`, one distinct subtype per block in S/N
  expect_equal(nrow(res$selected), data.table::uniqueN(res$selected$Block))
  perblk <- res$S[, .(n = data.table::uniqueN(Subtype)), by = Block]
  expect_true(all(perblk$n == 1L))
})

# ── selection by summed value + fitness_fn ────────────────────────────────────

test_that("max_effect selects the highest summed-value subtype per block", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2),
            fitness_fn = "max_effect", plots = character(0)))
  got <- unique(res$selected[, .(Block, Subtype)])
  chk <- merge(got, .expected_pick(res), by = "Block")
  expect_equal(chk$Subtype, chk$wantSub)
})

test_that("single-trait max_effect picks the highest-Effect subtype per block", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                  plots = character(0)))
  chk <- merge(res$S[, .(Block, S)],
               res$subtypes[, .(maxE = max(Effect)), by = Block], by = "Block")
  expect_equal(chk$S, chk$maxE)
})

test_that("frequent_best restricts the selected subtype to Count >= min_count", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
            fitness_fn = "frequent_best", min_count = 3, plots = character(0)))
  sm  <- merge(res$selected[, .(Block, Subtype, Count)], res$S[, .(Block, S)],
               by = "Block")
  for (i in seq_len(nrow(sm))) {
    pool <- res$subtypes[Block == sm$Block[i] & Count >= 3]
    if (nrow(pool) > 0L) {
      expect_gte(sm$Count[i], 3)
      expect_equal(sm$S[i], max(pool$Effect))
    }
  }
})

test_that("expected selects the highest frequency-weighted contribution per block", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                  fitness_fn = "expected",
                                                  plots = character(0)))
  sub <- data.table::copy(res$subtypes)
  sub[, contr := Effect * Count]                       # blockTot cancels in argmax
  want <- sub[, .SD[which.max(contr)], by = Block][, .(Block, wantSub = Subtype)]
  got  <- unique(res$selected[, .(Block, Subtype)])
  expect_equal(merge(got, want, by = "Block")[, Subtype],
               merge(got, want, by = "Block")[, wantSub])
})

test_that("min_effect selects the lowest summed-value subtype per block", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                  fitness_fn = "min_effect",
                                                  plots = character(0)))
  chk <- merge(res$S[, .(Block, S)],
               res$subtypes[, .(minE = min(Effect)), by = Block], by = "Block")
  expect_equal(chk$S, chk$minE)
})

test_that("a custom fitness_fn (one row per block) is honoured", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
            fitness_fn = function(dt) dt[which.min(dt$Effect)],
            plots = character(0)))
  chk <- merge(res$S[, .(Block, S)],
               res$subtypes[, .(minE = min(Effect)), by = Block], by = "Block")
  expect_equal(chk$S, chk$minE)
})

test_that("a custom fitness_fn returning != 1 row errors", {
  fx <- .analysis_fixture()
  expect_error(
    suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                        fitness_fn = function(dt) dt, plots = character(0))),
    "exactly")
})

# ── value_source: raw vs normalized selection ─────────────────────────────────

test_that("value_source is stored on the result", {
  fx  <- .analysis_fixture()
  r1  <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                  plots = character(0)))
  r2  <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                            value_source = "normalized", plots = character(0)))
  expect_equal(r1$value_source, "raw")
  expect_equal(r2$value_source, "normalized")
})

test_that("raw vs normalized can select a different subtype under mixed scales", {
  fx  <- .analysis_fixture()
  rr  <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2),
            value_source = "raw", plots = character(0)))
  rn  <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2),
            value_source = "normalized", plots = character(0)))
  # each honours its own summed criterion,
  expect_equal(merge(unique(rr$selected[, .(Block, Subtype)]),
                     .expected_pick(rr, col = "S"), by = "Block")[, Subtype],
               merge(unique(rr$selected[, .(Block, Subtype)]),
                     .expected_pick(rr, col = "S"), by = "Block")[, wantSub])
  expect_equal(merge(unique(rn$selected[, .(Block, Subtype)]),
                     .expected_pick(rn, col = "N"), by = "Block")[, Subtype],
               merge(unique(rn$selected[, .(Block, Subtype)]),
                     .expected_pick(rn, col = "N"), by = "Block")[, wantSub])
  # and the raw-scale sum (dominated by the large-scale trait) picks differently
  # from the normalized sum for at least one block.
  d <- merge(unique(rr$selected[, .(Block, a = Subtype)]),
             unique(rn$selected[, .(Block, b = Subtype)]), by = "Block")
  expect_gt(sum(d$a != d$b), 0L)
})

# ── signed trait weights (direction only) ─────────────────────────────────────

test_that("a signed trait_weight flips a lower-is-better trait in the sum", {
  fx  <- .analysis_fixture()
  w   <- c(Yield = 1, Protein = -1)
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2),
            trait_weights = w, plots = character(0)))
  got <- unique(res$selected[, .(Block, Subtype)])
  chk <- merge(got, .expected_pick(res, w = w, col = "S"), by = "Block")
  expect_equal(chk$Subtype, chk$wantSub)                # negative weight subtracts
})

test_that("trait_weights must cover every trait", {
  fx <- .analysis_fixture()
  expect_error(
    suppressWarnings(analyze_block_catalogue(fx$meta,
        list(Yield = fx$eff1, Protein = fx$eff2),
        trait_weights = c(Yield = 1), plots = character(0))),
    "missing trait")
})

# ── within-trait normalisation math (over all subtypes) ───────────────────────

test_that("zscore normalisation is a true standard score (mean 0, SD 1 per trait)", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2), normalize = "zscore",
            plots = character(0)))
  agg <- res$subtypes[, .(m = mean(N), s = sd(N)), by = Trait]
  expect_true(all(abs(agg$m) < 1e-8))                  # centred: mean 0 per trait
  expect_true(all(abs(agg$s - 1) < 1e-8))              # scaled: unit SD per trait
})

test_that("minmax normalisation scales by largest |effect| onto [-1, 1]", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2), normalize = "minmax",
            plots = character(0)))
  n <- res$subtypes
  expect_true(all(n$N >= -1 - 1e-9 & n$N <= 1 + 1e-9))
  agg <- n[, .(mx = max(abs(N))), by = Trait]
  expect_true(all(abs(agg$mx - 1) < 1e-8))
  expect_true(all(sign(n$N) == sign(n$S)))             # zero-anchored
})

test_that("rank normalisation is sign-preserving on [-1, 1]", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2), normalize = "rank",
            plots = character(0)))
  n <- res$subtypes
  expect_true(all(n$N >= -1 - 1e-9 & n$N <= 1 + 1e-9))
  expect_true(all(sign(n$N) == sign(n$S)))
})

test_that("normalize = 'none' keeps N == S (raw units)", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2), normalize = "none",
            plots = character(0)))
  expect_equal(res$subtypes$N, res$subtypes$S)
  expect_equal(res$N$N, res$S$S)
})

test_that("single-trait analysis still selects, and stays in raw units", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                  plots = character(0)))
  expect_equal(nrow(res$selected), data.table::uniqueN(res$selected$Block))
  # Raw by default: no transform runs, so N is S untouched.
  expect_identical(res$normalize, "none")
  expect_equal(res$subtypes$N, res$subtypes$S)
})

test_that("an explicitly requested transform still runs single-trait", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                            normalize = "zscore", plots = character(0)))
  expect_identical(res$normalize, "zscore")
  agg <- res$subtypes[, .(m = mean(N), s = sd(N)), by = Trait]
  expect_true(all(abs(agg$m) < 1e-8) && all(abs(agg$s - 1) < 1e-8))
})


# ── raw-by-default contract (0.5.4) ───────────────────────────────────────────

test_that("the default analysis rescales nothing and selects on raw effects", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2), plots = character(0)))

  expect_identical(res$normalize,    "none")
  expect_identical(res$value_source, "raw")
  # N is the identity transform of S everywhere it is reported.
  expect_equal(res$subtypes$N, res$subtypes$S)
  expect_equal(res$N$N,        res$S$S)
  # The selection is the raw-effect argmax — i.e. identical to asking for it.
  expl <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2),
            normalize = "none", value_source = "raw", plots = character(0)))
  expect_identical(unique(res$selected[, .(Block, Subtype)]),
                   unique(expl$selected[, .(Block, Subtype)]))
  # And the reported effects are the catalogue's own, in raw units.
  expect_true(max(abs(res$subtypes$S)) > 0)
  expect_equal(sort(res$subtypes$S),
               sort(as.numeric(res$subtypes$Effect)))
})

test_that("value_source = 'normalized' counts as asking, and promotes zscore", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(suppressMessages(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2),
            value_source = "normalized", plots = character(0))))

  expect_identical(res$normalize, "zscore")
  agg <- res$subtypes[, .(m = mean(N), s = sd(N)), by = Trait]
  expect_true(all(abs(agg$m) < 1e-8) && all(abs(agg$s - 1) < 1e-8))
  # It says so rather than promoting silently.
  expect_message(suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1), value_source = "normalized",
            plots = character(0))),
    "normalize = \"zscore\"")
})

test_that("value_source = 'normalized' with an explicit normalize = 'none' errors", {
  fx <- .analysis_fixture()
  expect_error(
    suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
        normalize = "none", value_source = "normalized",
        plots = character(0))),
    "has no meaning")
})

test_that("an explicit normalize is honoured under value_source = 'normalized'", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2),
            normalize = "minmax", value_source = "normalized",
            plots = character(0)))
  expect_identical(res$normalize, "minmax")
  expect_true(all(abs(res$subtypes$N) <= 1 + 1e-12))
})

test_that("a raw result re-plotted as normalized is not labelled a standard score", {
  fx  <- .analysis_fixture()
  raw <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2), plots = character(0)))
  # N holds raw values here, so the subtitle must not claim a standard score.
  p <- plot_block_effect_genome(raw, value_source = "normalized")
  expect_false(grepl("standard score", p$labels$subtitle))

  zs <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2),
            normalize = "zscore", value_source = "normalized",
            plots = character(0)))
  expect_true(grepl("standard score",
                    plot_block_effect_genome(zs)$labels$subtitle))
})

# ── per-trait top blocks ──────────────────────────────────────────────────────

test_that("per-trait top block equals argmax of the selected subtype's value", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2), plots = character(0)))
  for (tr in unique(res$S$Trait)) {
    st  <- res$S[Trait == tr]
    top <- res$top_blocks[Trait == tr]
    expect_true(all(top$S == max(st$S)))
    expect_true(all(st[is_top == TRUE, S] == max(st$S)))
  }
})

# ── data-source resolution ────────────────────────────────────────────────────

test_that("catalogue path runs and returns valid tables", {
  fx   <- .analysis_fixture()
  capt <- suppressWarnings(suppressMessages(
    run_haplotype_catalogue(fx$meta, effect_list = list(Yield = fx$eff1),
                            output_dir = tempfile("cat"), save_profile = FALSE)))
  res  <- suppressWarnings(analyze_block_catalogue(catalogue = capt,
                                                   plots = character(0)))
  expect_s3_class(res$selected, "data.table")
  expect_gt(nrow(res$selected), 0L)
})

test_that("errors when neither catalogue nor geno_meta+effect_list given", {
  expect_error(analyze_block_catalogue(), "catalogue")
})

test_that("selected_blocks restricts the analysis to the chosen blocks", {
  fx   <- .analysis_fixture()
  full <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                   plots = character(0)))
  keep <- head(unique(full$S$Block), 2L)
  res  <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                            selected_blocks = keep, plots = character(0)))
  expect_setequal(unique(res$S$Block), keep)
  expect_setequal(unique(res$subtypes$Block), keep)
})

# ── gated writes ──────────────────────────────────────────────────────────────

test_that("output_dir writes the revised table set (+ PNGs for requested plots)", {
  fx <- .analysis_fixture()
  od <- tempfile("abc")
  suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                           output_dir = od, plots = character(0)))
  expect_setequal(list.files(od),
    c("block_subtypes.txt", "block_selected.txt", "block_S.txt",
      "block_N.txt", "block_top.txt"))

  od2 <- tempfile("abc2")
  suppressWarnings(suppressMessages(analyze_block_catalogue(
    fx$meta, list(Yield = fx$eff1), output_dir = od2, plots = "density")))
  expect_true("block_density.png" %in% list.files(od2))
})

test_that("NULL output_dir writes nothing and no longer warns (0.9.0, D3.4)", {
  fx <- .analysis_fixture()
  d <- withr::local_tempdir(); withr::local_dir(d)
  expect_no_warning(
    res <- suppressMessages(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                    plots = character(0))))
  expect_s3_class(res, "hs_block_analysis")
  expect_identical(list.files(d), character(0))
})

test_that("write_outputs = FALSE is deprecated but still honoured (0.9.0, D3.4)", {
  fx <- .analysis_fixture()
  od <- tempfile("nowrite")
  expect_warning(suppressMessages(
    analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                            output_dir = od, write_outputs = FALSE,
                            plots = character(0))),
    class = "hs_warning_deprecated")
  expect_false(dir.exists(od) && length(list.files(od)) > 0L)
})

# ── plot layer honours value_source (D2) ──────────────────────────────────────

test_that("analyze_block_catalogue builds the requested ggplot figures", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(suppressMessages(analyze_block_catalogue(
    fx$meta, list(Yield = fx$eff1, Protein = fx$eff2),
    plots = c("rectangular", "manhattan", "density", "heatmap"))))
  expect_setequal(names(res$plots),
                  c("rectangular", "manhattan", "density", "heatmap"))
  expect_true(all(vapply(res$plots, inherits, logical(1L), "ggplot")))
})

test_that("each standalone plotter returns a ggplot", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2), plots = character(0)))
  expect_s3_class(plot_block_effect_genome(res), "ggplot")
  expect_s3_class(plot_block_effect_manhattan(res), "ggplot")
  expect_s3_class(plot_block_effect_density(res), "ggplot")
  expect_s3_class(plot_top_block_heatmap(res, top_n = 2L), "ggplot")
})

test_that("built plots fill by the stored value_source, override re-plots the other", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2),
            value_source = "normalized", plots = character(0)))
  gN  <- plot_block_effect_genome(res)                 # stored = normalized
  gR  <- plot_block_effect_genome(res, value_source = "raw")   # override
  expect_equal(sort(gN$data$PlotVal), sort(res$N$N))
  expect_equal(sort(gR$data$PlotVal), sort(res$S$S))
})

test_that("no plotter accepts or uses a 'corrected'/'metric' argument", {
  for (fn in list(plot_block_effect_genome, plot_block_effect_manhattan,
                  plot_block_effect_density, plot_block_effect_circos)) {
    fm <- names(formals(fn))
    expect_false(any(c("metric", "value") %in% fm))
    expect_true("value_source" %in% fm)
  }
})

test_that("plots = character(0) yields no figures", {
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                  plots = character(0)))
  expect_length(res$plots, 0L)
})

test_that("plot_block_effect_circos falls back to a ggplot without circlize", {
  skip_if(requireNamespace("circlize", quietly = TRUE), "circlize is installed")
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta, list(Yield = fx$eff1),
                                                  plots = character(0)))
  expect_message(p <- plot_block_effect_circos(res), "circlize")
  expect_s3_class(p, "ggplot")
})

test_that("plot_block_effect_circos renders a recorded plot with circlize", {
  skip_if_not_installed("circlize")
  fx  <- .analysis_fixture()
  res <- suppressWarnings(analyze_block_catalogue(fx$meta,
            list(Yield = fx$eff1, Protein = fx$eff2), plots = character(0)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_warning(p <- plot_block_effect_circos(res, link_top = TRUE),
                 class = "hs_warning_deprecated")          # row 12 b (0.9.0)
  expect_s3_class(p, "recordedplot")
})

test_that("a plotter errors on a non-analysis object", {
  expect_error(plot_block_effect_genome(list(foo = 1)), "analyze_block_catalogue")
})
