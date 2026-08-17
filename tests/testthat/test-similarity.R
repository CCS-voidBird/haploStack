# =============================================================================
# test-similarity.R
#
# QA tests for the IBD/IBS & marker-signal module (R2):
#   compute_segment_ibs(), compute_block_ibd(), score_marker_signal(),
#   compute_signal_metrics().
# =============================================================================

library(testthat)
library(haploStack)
library(data.table)

# ---------------------------------------------------------------------------
# Hand-verifiable fixture
#
# Family: parents P1, P2; offspring O1. One block B1 = {m1..m4} on chr 1.
#   P1 hap1 = AAAA   P1 hap2 = TTTT
#   P2 hap1 = AAAA   P2 hap2 = TTTT
#   O1 hap1 = AAAA   O1 hap2 = AATT
#
# 4 parental haplotypes x 2 progeny haplotypes = 8 pairwise comparisons.
# Block similarities: P*.1 vs O1.1 = 1.0; P*.1 vs O1.2 = 0.5;
#                     P*.2 vs O1.1 = 0.0; P*.2 vs O1.2 = 0.5.
# ---------------------------------------------------------------------------
.sim_fixture <- local({
  samples <- c("P1", "P2", "O1")
  markers <- c("m1", "m2", "m3", "m4")
  arr <- array(NA_character_, dim = c(3L, 4L, 2L),
               dimnames = list(samples, markers, c("1", "2")))
  arr["P1", , 1] <- c("A","A","A","A"); arr["P1", , 2] <- c("T","T","T","T")
  arr["P2", , 1] <- c("A","A","A","A"); arr["P2", , 2] <- c("T","T","T","T")
  arr["O1", , 1] <- c("A","A","A","A"); arr["O1", , 2] <- c("A","A","T","T")

  map <- data.frame(Chrom = "1", Marker = markers, Pos = c(10, 20, 30, 40),
                    stringsAsFactors = FALSE)
  bmm <- data.frame(Chrom = "1", Marker = markers, Block = "B1",
                    stringsAsFactors = FALSE)

  meta0 <- new_geno_meta(
    geno_arrays           = arr,
    map_info              = map,
    block_marker_map_info = bmm,
    metadata = list(ploidy = 2L, char_dict = c("A" = 1L, "T" = 2L),
                    n_samples = 3L, n_markers = 4L,
                    source_format = "nucleotide", encoding = "nucleotide"))

  fam  <- data.frame(IID = samples,
                     Father = c(NA, NA, "P1"),
                     Mother = c(NA, NA, "P2"), stringsAsFactors = FALSE)
  meta <- attach_family_history(meta0, fam, parent_cols = c("Father", "Mother"))

  list(meta = meta, meta0 = meta0,
       ped = data.frame(Sample = "O1", Parent1 = "P1", Parent2 = "P2",
                        stringsAsFactors = FALSE))
})

# ---------------------------------------------------------------------------
# Selection fixture — the same panel split into TWO blocks, so the per-block
# "top N" rule has something to choose between:
#   B1 = {m1 (Pos 10), m2 (Pos 20)}   B2 = {m3 (Pos 30), m4 (Pos 40)}
# It carries a block_info table (via attach_block_info) for the signal columns
# to be written into, and a family history so score_marker_signal() can run.
# ---------------------------------------------------------------------------
.sel_fixture <- local({
  samples <- c("P1", "P2", "O1")
  markers <- c("m1", "m2", "m3", "m4")
  arr <- array(NA_character_, dim = c(3L, 4L, 2L),
               dimnames = list(samples, markers, c("1", "2")))
  arr["P1", , 1] <- c("A","A","A","A"); arr["P1", , 2] <- c("T","T","T","T")
  arr["P2", , 1] <- c("A","A","A","A"); arr["P2", , 2] <- c("T","T","T","T")
  arr["O1", , 1] <- c("A","A","A","A"); arr["O1", , 2] <- c("A","A","T","T")

  map <- data.frame(Chrom = "1", Marker = markers, Pos = c(10, 20, 30, 40),
                    stringsAsFactors = FALSE)
  bmm <- data.frame(Chrom = "1", Marker = markers,
                    Block = c("B1", "B1", "B2", "B2"),
                    stringsAsFactors = FALSE)

  meta <- new_geno_meta(
    geno_arrays           = arr,
    map_info              = map,
    block_marker_map_info = bmm,
    metadata = list(ploidy = 2L, char_dict = c("A" = 1L, "T" = 2L),
                    n_samples = 3L, n_markers = 4L,
                    source_format = "nucleotide", encoding = "nucleotide"))
  meta <- suppressMessages(attach_block_info(meta, bmm))
  fam  <- data.frame(IID = samples,
                     Father = c(NA, NA, "P1"),
                     Mother = c(NA, NA, "P2"), stringsAsFactors = FALSE)
  meta <- attach_family_history(meta, fam, parent_cols = c("Father", "Mother"))

  # Hand-made metrics: B1 -> m2 wins (0.4 > 0.1); B2 -> m3 wins (0.4 > 0.2).
  metrics <- data.table(Chrom = "1", Marker = markers, Pos = c(10, 20, 30, 40),
                        final_mcc = c(0.1, 0.4, 0.4, 0.2),
                        MAF = c(0.1, 0.2, 0.5, 0.3))
  list(meta = meta, bmm = bmm, metrics = metrics)
})

# ── compute_segment_ibs ───────────────────────────────────────────────────────

test_that("compute_segment_ibs returns the NA-aware match ratio", {
  expect_equal(compute_segment_ibs(c("A","T","G"), c("A","C","G")), 2/3)
  expect_equal(compute_segment_ibs(c("A","T"), c("A","T")), 1)
  expect_equal(compute_segment_ibs(c("A","T"), c("G","C")), 0)
})

test_that("compute_segment_ibs excludes NA / empty positions", {
  expect_equal(compute_segment_ibs(c("A", NA, "G"), c("A", "C", "G")), 1)
  expect_true(is.na(compute_segment_ibs(c(NA, NA), c("A", "T"))))
})

test_that("compute_segment_ibs weights by block length when asked", {
  expect_equal(
    compute_segment_ibs(c("A","T"), c("A","T"), block_length = 1e5,
                        weight_by_length = TRUE), 1e5)
  expect_error(
    compute_segment_ibs(c("A"), c("A"), weight_by_length = TRUE),
    "block_length")
  expect_error(compute_segment_ibs(c("A","T"), "A"), "equal length")
})

# ── compute_block_ibd ─────────────────────────────────────────────────────────

test_that("compute_block_ibd returns one row per pairwise comparison", {
  ibd <- compute_block_ibd(.sim_fixture$meta)
  expect_s3_class(ibd, "data.table")
  expect_true(all(c("Gen","FID","Parent","ParentHap","Progeny","ProgenyHap",
                    "Block","Chrom","Similarity","IBD") %in% names(ibd)))
  expect_equal(nrow(ibd), 8L)                       # 4 parental haps x 2 progeny haps
  expect_setequal(round(ibd$Similarity, 3), c(1, 0.5, 0, 0.5))
})

test_that("compute_block_ibd applies an inclusive threshold", {
  ibd <- compute_block_ibd(.sim_fixture$meta, threshold = 0.5)
  # similarity exactly 0.5 counts as IBD
  expect_true(all(ibd$IBD[abs(ibd$Similarity - 0.5) < 1e-9]))
  # similarity 0 is not IBD
  expect_false(any(ibd$IBD[ibd$Similarity == 0]))
  # a stricter-than-data threshold removes the 0.5 calls
  ibd2 <- compute_block_ibd(.sim_fixture$meta, threshold = 0.75)
  expect_false(any(ibd2$IBD[abs(ibd2$Similarity - 0.5) < 1e-9]))
})

test_that("compute_block_ibd works from a flat pedigree", {
  ibd <- compute_block_ibd(.sim_fixture$meta0, pedigree = .sim_fixture$ped)
  expect_equal(nrow(ibd), 8L)
  expect_true(all(ibd$Gen == "ped"))
})

test_that("compute_block_ibd errors without a family source", {
  expect_error(compute_block_ibd(.sim_fixture$meta0), "family_history")
})

# ── score_marker_signal ───────────────────────────────────────────────────────

test_that("score_marker_signal (IBD mode) builds the standard confusion matrix", {
  res <- score_marker_signal(.sim_fixture$meta, mode = "IBD")
  expect_named(res, c("detail", "marker", "metrics", "signal"))
  expect_true(all(c("Chrom","Marker","TP","FP","TN","FN") %in% names(res$marker)))

  m1 <- res$marker[Marker == "m1"]
  expect_equal(m1$TP, 4); expect_equal(m1$FP, 0)
  expect_equal(m1$FN, 2); expect_equal(m1$TN, 2)
  # every marker is covered by all 8 comparisons
  res$marker[, tot := TP + FP + TN + FN]
  expect_true(all(res$marker$tot == 8))
})

test_that("score_marker_signal (IBS mode) weights cells by similarity", {
  res <- score_marker_signal(.sim_fixture$meta, mode = "IBS")
  m1 <- res$marker[Marker == "m1"]
  expect_equal(m1$TP, 3); expect_equal(m1$FP, 1)
  expect_equal(m1$FN, 1); expect_equal(m1$TN, 3)
  # fractional cells still total the number of comparisons (8) per marker
  res$marker[, tot := TP + FP + TN + FN]
  expect_true(all(abs(res$marker$tot - 8) < 1e-9))
})

test_that("score_marker_signal detail carries the ParentHap dimension", {
  res <- score_marker_signal(.sim_fixture$meta, mode = "IBD")
  expect_true(all(c("Gen","FID","IID_parent","ParentHap","Block","Marker") %in%
                    names(res$detail)))
  expect_setequal(unique(res$detail$ParentHap), c(1L, 2L))
})

test_that("score_marker_signal metrics carry MAF and final scores", {
  res <- score_marker_signal(.sim_fixture$meta, mode = "IBD")
  expect_true(all(c("accuracy","precision","recall","specificity","f1","mcc",
                    "MAF","final_mcc") %in% names(res$metrics)))
  # m1 alleles across all sample-haplotypes: A x4, T x2 -> MAF = 1 - 4/6 = 1/3
  expect_equal(res$metrics[Marker == "m1", MAF], 1/3, tolerance = 1e-6)
})

test_that("score_marker_signal rejects an unknown mode", {
  expect_error(score_marker_signal(.sim_fixture$meta, mode = "foo"))
})

# ── compute_signal_metrics ────────────────────────────────────────────────────

test_that("compute_signal_metrics computes the textbook metrics", {
  tab <- data.table(Chrom = "1", Marker = "m1", TP = 4, FP = 0, TN = 2, FN = 2)
  out <- compute_signal_metrics(tab, maf = c(m1 = 0.5))

  expect_equal(out$accuracy, 0.75)
  expect_equal(out$precision, 1)
  expect_equal(out$recall, 4/6, tolerance = 1e-6)
  expect_equal(out$specificity, 1)
  expect_equal(out$f1, 0.8, tolerance = 1e-6)
  expect_equal(out$mcc, 8 / sqrt(192), tolerance = 1e-6)
  # final = metric * MAF * (1 - MAF) = metric * 0.25
  expect_equal(out$final_accuracy, 0.75 * 0.25, tolerance = 1e-6)
})

test_that("compute_signal_metrics fills zero denominators with 0 (no NA)", {
  tab <- data.table(Chrom = "1", Marker = "mZ", TP = 0, FP = 0, TN = 5, FN = 0)
  out <- compute_signal_metrics(tab, maf = c(mZ = 0.4))
  expect_equal(out$precision, 0)        # TP+FP = 0
  expect_equal(out$recall, 0)           # TP+FN = 0
  expect_equal(out$f1, 0)               # precision+recall = 0
  expect_equal(out$mcc, 0)              # a marginal sum is 0
  expect_equal(out$accuracy, 1)
  metric_cols <- c("accuracy", "precision", "recall", "specificity", "f1", "mcc")
  expect_false(anyNA(out[, c(metric_cols, "MAF",
                             paste0("final_", metric_cols)), with = FALSE]))
})

test_that("compute_signal_metrics zero_division = NA restores NA behaviour", {
  tab <- data.table(Chrom = "1", Marker = "mZ", TP = 0, FP = 0, TN = 5, FN = 0)
  out <- compute_signal_metrics(tab, maf = c(mZ = 0.4), zero_division = NA)
  expect_true(is.na(out$precision))     # TP+FP = 0
  expect_true(is.na(out$recall))        # TP+FN = 0
  expect_equal(out$accuracy, 1)
})

test_that("compute_signal_metrics fills missing MAF with 0 (final_* = 0)", {
  tab <- data.table(Chrom = "1", Marker = "mNoMAF",
                    TP = 4, FP = 0, TN = 2, FN = 2)
  expect_warning(
    out <- compute_signal_metrics(tab, maf = c(other = 0.3)),
    "no MAF")
  expect_equal(out$MAF, 0)
  expect_equal(out$final_mcc, 0)
  expect_false(anyNA(out$final_accuracy))
})

test_that("compute_signal_metrics needs a MAF source", {
  tab <- data.table(Chrom = "1", Marker = "m1", TP = 1, FP = 1, TN = 1, FN = 1)
  expect_error(compute_signal_metrics(tab), "geno_meta")
})

# ── plot_marker_signal ────────────────────────────────────────────────────────

test_that("plot_marker_signal returns a ggplot (positions from geno_meta)", {
  res <- score_marker_signal(.sim_fixture$meta, mode = "IBD")
  p   <- plot_marker_signal(res$metrics, geno_meta = .sim_fixture$meta,
                            score = "final_mcc")
  expect_s3_class(p, "ggplot")
})

test_that("plot_marker_signal errors without a position source", {
  res <- score_marker_signal(.sim_fixture$meta, mode = "IBD")
  expect_error(plot_marker_signal(res$metrics), "geno_meta|map_info")
})

test_that("plot_marker_signal errors on an unknown score column", {
  res <- score_marker_signal(.sim_fixture$meta, mode = "IBD")
  expect_error(
    plot_marker_signal(res$metrics, geno_meta = .sim_fixture$meta,
                       score = "not_a_column"),
    "not found")
})

# Count the text layers of a plot — marker labels add one on top of the
# always-present "max = ..." annotation and the block-strip labels.
.n_text_layers <- function(p)
  sum(vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1L)))

test_that("plot_marker_signal draws no marker labels and no title by default", {
  sel  <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 1)
  meta <- attach_signal_markers(.sel_fixture$meta, sel)
  p    <- plot_marker_signal(.sel_fixture$metrics, geno_meta = meta)

  expect_s3_class(p, "ggplot")
  expect_null(p$labels$title)
  base_text <- .n_text_layers(p)

  p_lab <- plot_marker_signal(.sel_fixture$metrics, geno_meta = meta,
                              show_marker_labels = TRUE, top_n = 2)
  expect_equal(.n_text_layers(p_lab), base_text + 1L)
})

test_that("plot_marker_signal draws a title only when asked", {
  p_off <- plot_marker_signal(.sel_fixture$metrics, geno_meta = .sel_fixture$meta,
                              highlight = "quantile")
  expect_null(p_off$labels$title)

  p_on <- plot_marker_signal(.sel_fixture$metrics, geno_meta = .sel_fixture$meta,
                             highlight = "quantile", show_title = TRUE)
  expect_equal(p_on$labels$title, "Marker signal — final_mcc")

  # an explicit title implies show_title
  p_txt <- plot_marker_signal(.sel_fixture$metrics, geno_meta = .sel_fixture$meta,
                              highlight = "quantile", title = "My run")
  expect_equal(p_txt$labels$title, "My run")
})

test_that("plot_marker_signal highlights the selected signal markers", {
  sel  <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 1)
  meta <- attach_signal_markers(.sel_fixture$meta, sel)

  p <- plot_marker_signal(.sel_fixture$metrics, geno_meta = meta)
  expect_setequal(p$data[candidate == TRUE]$Marker, c("m2", "m3"))

  # the same selection passed directly, without touching the container
  p2 <- plot_marker_signal(.sel_fixture$metrics, geno_meta = .sel_fixture$meta,
                           signal_markers = sel)
  expect_setequal(p2$data[candidate == TRUE]$Marker, c("m2", "m3"))

  # The two rules genuinely differ: give B2 only weak markers, and its best
  # marker (m3) is still that block's representative even though it is nowhere
  # near the genome-wide top quantile — which is the point of the change.
  weak <- data.table::copy(.sel_fixture$metrics)
  weak[Marker == "m3", final_mcc := 0.05]
  weak[Marker == "m4", final_mcc := 0.02]
  meta_weak <- attach_signal_markers(
    .sel_fixture$meta,
    select_signal_markers(weak, geno_meta = .sel_fixture$meta, n = 1))

  p_sig <- plot_marker_signal(weak, geno_meta = meta_weak)
  expect_setequal(p_sig$data[candidate == TRUE]$Marker, c("m2", "m3"))

  p_qnt <- plot_marker_signal(weak, geno_meta = meta_weak,
                              highlight = "quantile", quantile = 0.95)
  expect_setequal(p_qnt$data[candidate == TRUE]$Marker, "m2")

  expect_equal(sum(plot_marker_signal(.sel_fixture$metrics, geno_meta = meta,
                                      highlight = "none")$data$candidate), 0L)
})

test_that("plot_marker_signal falls back to the quantile rule with a message", {
  expect_message(
    p <- plot_marker_signal(.sel_fixture$metrics, geno_meta = .sel_fixture$meta),
    "no signal-marker selection")
  expect_s3_class(p, "ggplot")

  expect_error(
    plot_marker_signal(.sel_fixture$metrics, geno_meta = .sel_fixture$meta,
                       highlight = "signal"),
    "no signal markers")
})

test_that("plot_marker_signal carries a marker x block table", {
  sel  <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 1)
  meta <- attach_signal_markers(.sel_fixture$meta, sel)
  p    <- plot_marker_signal(.sel_fixture$metrics, geno_meta = meta)

  tab <- attr(p, "signal_table")
  expect_s3_class(tab, "data.table")
  expect_equal(nrow(tab), 4L)                      # one row per plotted marker
  expect_true(all(c("Chrom", "Pos", "Marker", "Block", "is_signal",
                    "score_val", "block_score", "block_n_markers") %in%
                    names(tab)))
  expect_setequal(tab[is_signal == TRUE]$Marker, c("m2", "m3"))
  expect_equal(tab[Block == "B1"]$block_n_markers, c(2L, 2L))
  expect_equal(unique(tab[Block == "B1"]$block_score), 0.4)
})

test_that("save_marker_signal writes the table beside the image", {
  skip_if_not_installed("ggplot2")
  sel  <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 1)
  meta <- attach_signal_markers(.sel_fixture$meta, sel)
  p    <- plot_marker_signal(.sel_fixture$metrics, geno_meta = meta)

  png <- file.path(tempdir(), "marker_signal_test.png")
  csv <- file.path(tempdir(), "marker_signal_test.csv")
  on.exit(unlink(c(png, csv)), add = TRUE)

  suppressMessages(save_marker_signal(p, png, width = 6, height = 3, dpi = 72))
  expect_true(file.exists(png))
  expect_true(file.exists(csv))
  expect_equal(nrow(data.table::fread(csv)), 4L)
})

test_that("plot_marker_signal writes the table to table_file on request", {
  csv <- file.path(tempdir(), "marker_signal_direct.csv")
  on.exit(unlink(csv), add = TRUE)
  suppressMessages(
    plot_marker_signal(.sel_fixture$metrics, geno_meta = .sel_fixture$meta,
                       highlight = "quantile", table_file = csv))
  expect_true(file.exists(csv))
})

# ── select_signal_markers / attach / get ──────────────────────────────────────

test_that("select_signal_markers picks the best marker of each block", {
  sel <- select_signal_markers(.sel_fixture$metrics,
                               geno_meta = .sel_fixture$meta, n = 1)
  expect_s3_class(sel, "data.table")
  expect_equal(nrow(sel), 2L)
  expect_setequal(sel$Block, c("B1", "B2"))
  expect_equal(sel[Block == "B1"]$Marker, "m2")
  expect_equal(sel[Block == "B2"]$Marker, "m3")
  expect_true(all(sel$signal_rank == 1L))
  expect_equal(attr(sel, "score"), "final_mcc")
})

test_that("select_signal_markers keeps n per block (and fewer when short)", {
  sel <- select_signal_markers(.sel_fixture$metrics,
                               geno_meta = .sel_fixture$meta, n = 2)
  expect_equal(nrow(sel), 4L)
  expect_equal(sel[Block == "B1"][order(signal_rank)]$Marker, c("m2", "m1"))
  # asking for more than a block holds yields what exists, not NA padding
  sel3 <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 5)
  expect_equal(nrow(sel3), 4L)
  expect_equal(sel3[, .N, by = Block]$N, c(2L, 2L))
})

test_that("select_signal_markers breaks ties deterministically", {
  m <- data.table::copy(.sel_fixture$metrics)
  m[Marker == "m1", final_mcc := 0.4]        # m1 and m2 now tie inside B1

  by_pos <- select_signal_markers(m, geno_meta = .sel_fixture$meta, n = 1,
                                  ties = "position")
  expect_equal(by_pos[Block == "B1"]$Marker, "m1")   # Pos 10 < 20

  by_maf <- select_signal_markers(m, geno_meta = .sel_fixture$meta, n = 1,
                                  ties = "maf")
  expect_equal(by_maf[Block == "B1"]$Marker, "m2")   # MAF 0.2 > 0.1
})

test_that("select_signal_markers honours min_score and validates its inputs", {
  sel <- select_signal_markers(.sel_fixture$metrics,
                               geno_meta = .sel_fixture$meta, n = 1,
                               min_score = 0.35)
  expect_equal(nrow(sel), 2L)                       # both block winners are 0.4
  sel2 <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 1,
                                min_score = 0.45)
  expect_equal(nrow(sel2), 0L)                      # nothing clears the floor
  expect_error(select_signal_markers(.sel_fixture$metrics,
                                     geno_meta = .sel_fixture$meta, n = 0),
               "positive integer")
  expect_error(select_signal_markers(.sel_fixture$metrics,
                                     geno_meta = .sel_fixture$meta,
                                     score = "nope"), "not found")
  expect_error(select_signal_markers(.sel_fixture$metrics), "block")
})

test_that("attach_signal_markers records the selection on block_info", {
  sel  <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 1)
  meta <- attach_signal_markers(.sel_fixture$meta, sel)

  expect_true(is_geno_meta(meta))
  bi <- as.data.table(meta$block_info)
  expect_true(all(c("Signal_Marker", "Signal_Score", "N_Signal") %in% names(bi)))
  expect_equal(bi[Block == "B1"]$Signal_Marker, "m2")
  expect_equal(bi[Block == "B1"]$Signal_Score, 0.4)
  expect_equal(bi[Block == "B2"]$N_Signal, 1L)
  expect_equal(meta$metadata$signal_selection$score, "final_mcc")
  expect_equal(meta$metadata$signal_selection$n_blocks, 2L)
})

test_that("get_signal_markers round-trips a multi-marker selection", {
  sel  <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 2)
  meta <- attach_signal_markers(.sel_fixture$meta, sel)

  long <- get_signal_markers(meta)
  expect_equal(nrow(long), 4L)
  expect_equal(long[Block == "B1"][order(signal_rank)]$Marker, c("m2", "m1"))
  expect_setequal(get_signal_markers(meta, long = FALSE),
                  c("m1", "m2", "m3", "m4"))
  # no selection stored -> empty, not an error
  expect_equal(nrow(get_signal_markers(.sel_fixture$meta)), 0L)
})

test_that("attach_signal_markers derives block_info when it is missing", {
  # .sim_fixture$meta carries a block map but no block_info (the common case
  # for a GenoMeta built straight from a marker->block table).
  expect_null(.sim_fixture$meta$block_info)
  res <- score_marker_signal(.sim_fixture$meta, mode = "IBD", n_signal = 1)
  expect_message(meta <- attach_signal_markers(.sim_fixture$meta, res$signal),
                 "derived it")
  expect_true(is_geno_meta(meta))
  expect_equal(nrow(get_signal_markers(meta)), 1L)   # single block B1
})

test_that("attach_signal_markers refuses to silently overwrite", {
  sel  <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 1)
  meta <- attach_signal_markers(.sel_fixture$meta, sel)
  expect_error(attach_signal_markers(meta, sel, overwrite = FALSE),
               "overwrite")
  expect_silent(suppressMessages(attach_signal_markers(meta, sel)))
})

test_that("score_marker_signal selects signal markers as part of scoring", {
  res <- score_marker_signal(.sel_fixture$meta, mode = "IBD", n_signal = 1)
  expect_s3_class(res$signal, "data.table")
  expect_setequal(res$signal$Block, c("B1", "B2"))
  expect_true(all(res$signal$signal_rank == 1L))

  direct <- select_signal_markers(res$metrics, geno_meta = .sel_fixture$meta,
                                  n = 1)
  expect_equal(res$signal$Marker, direct$Marker)

  res0 <- score_marker_signal(.sel_fixture$meta, mode = "IBD", n_signal = 0)
  expect_equal(nrow(res0$signal), 0L)
})

# ── selection lifecycle across the container ──────────────────────────────────

test_that("filter_geno_meta prunes signal markers that leave the panel", {
  sel  <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 2)
  meta <- attach_signal_markers(.sel_fixture$meta, sel)

  expect_message(
    sub <- filter_geno_meta(meta, keep_markers = c("m1", "m3", "m4")),
    "signal marker")
  bi <- as.data.table(sub$block_info)
  expect_equal(bi[Block == "B1"]$Signal_Marker, "m1")   # m2 was pruned
  expect_equal(bi[Block == "B1"]$N_Signal, 1L)
  expect_true(is.na(bi[Block == "B1"]$Signal_Score))     # rank-1 marker is gone
  expect_equal(bi[Block == "B2"]$N_Signal, 2L)           # untouched
})

test_that("attach_block_info drops a stale selection with a message", {
  sel  <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 1)
  meta <- attach_signal_markers(.sel_fixture$meta, sel)

  expect_message(meta2 <- attach_block_info(meta, .sel_fixture$bmm),
                 "signal-marker selection was dropped")
  expect_false("Signal_Marker" %in% names(meta2$block_info))
})

test_that("build_block_location_info carries the selection through", {
  sel  <- select_signal_markers(.sel_fixture$metrics,
                                geno_meta = .sel_fixture$meta, n = 1)
  meta <- attach_signal_markers(.sel_fixture$meta, sel)

  loc <- build_block_location_info(meta)
  expect_true("Signal_Marker" %in% names(loc$block_info))
  expect_equal(as.data.table(loc$block_info)[Block == "B1"]$Signal_Marker, "m2")
})

# ── plot_signal_precision_recall ──────────────────────────────────────────────

test_that("plot_signal_precision_recall returns a ggplot", {
  res <- score_marker_signal(.sim_fixture$meta, mode = "IBD")
  p   <- plot_signal_precision_recall(res$metrics)
  expect_s3_class(p, "ggplot")
})

test_that("plot_signal_precision_recall needs precision/recall/MAF", {
  expect_error(
    plot_signal_precision_recall(data.table(Marker = "m1", precision = 1)),
    "missing column")
})
