# =============================================================================
# test-catalogue.R
#
# QA tests for the catalogue module:
#   geno_meta_views(), filter_geno_meta(), save/load_geno_meta(),
#   recode_geno_meta(), compute_profile_localgebv(), build_block_location_info(),
#   plus the new sample_info / marker_effects GenoMeta fields.
# =============================================================================

library(testthat)
library(haploStack)

# Shared test data (created once per file via local helper)
.test_data <- local({
  set.seed(42L)
  geno <- sim_geno_matrix(n_samples = 6L, n_markers = 12L,
                           encoding = "nucleotide")
  map  <- sim_map_data(rownames(geno), n_chrom = 2L)
  bmap <- sim_block_map(map, block_size = 3L)
  eff  <- sim_effect_df(rownames(geno), alleles = c("A","T"))
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 12L,
                         n_chrom = 2L, block_size = 3L, seed = 42L)
  list(geno = geno, map = map, bmap = bmap, eff = eff, meta = meta)
})

# Profile views are derived on demand (profile is no longer a separate object)
.views <- function(meta) geno_meta_views(meta)

# ── geno_meta_views (profile views derived from GenoMeta) ──────────────────────

test_that("geno_meta_views returns valid view tables", {
  v <- .views(.test_data$meta)
  expect_true(is.list(v))
  expect_true(is.data.table(v$alleles))
  expect_true(is.data.table(v$marker_meta))
  expect_true(is.data.table(v$block_meta))
})

test_that("geno_meta_views alleles table has correct dimensions", {
  v      <- .views(.test_data$meta)
  n_samp <- ncol(.test_data$geno)
  n_mark <- nrow(.test_data$geno)
  # Diploid encoding → 2 rows per sample
  expect_equal(nrow(v$alleles), n_samp * 2L)
  # Columns: Sample + Haplotype_idx + one column per marker
  expect_equal(ncol(v$alleles), n_mark + 2L)
})

test_that("geno_meta_views detects encoding automatically", {
  v <- .views(.test_data$meta)
  expect_equal(v$encoding, "nucleotide")
})

test_that("geno_meta_views errors when input is not a GenoMeta", {
  expect_error(.views(list()),           "GenoMeta")
  expect_error(.views(.test_data$geno),  "GenoMeta")
})

test_that("geno_meta_views returns empty block_meta when blocks are absent", {
  meta_no_blocks <- sim_geno_meta(n_samples = 4L, n_markers = 8L,
                                   block_size = NULL, seed = 99L)
  v <- .views(meta_no_blocks)
  expect_equal(nrow(v$block_meta), 0L)
})

test_that("geno_meta_views marker_meta is sorted by genomic position", {
  v     <- .views(.test_data$meta)
  pos   <- v$marker_meta$position
  chrom <- v$marker_meta$chromosome
  expect_true(all(diff(as.numeric(pos[chrom == "1"])) >= 0))
})

test_that("geno_meta_views block_meta n_markers matches block_map", {
  v           <- .views(.test_data$meta)
  block_sizes <- table(.test_data$bmap$Block)
  for (blk in v$block_meta$block_id) {
    expect_equal(
      v$block_meta[block_id == blk, n_markers],
      as.integer(block_sizes[blk]),
      label = paste("Block", blk, "marker count")
    )
  }
})

test_that("geno_meta_views works with nucleotide encoding", {
  meta_nuc <- sim_geno_meta(n_samples = 4L, n_markers = 8L,
                              encoding = "nucleotide", block_size = 4L)
  v <- .views(meta_nuc)
  expect_equal(v$encoding, "nucleotide")
})

# ── sample_info (required field) ──────────────────────────────────────────────

test_that("new_geno_meta auto-builds a FAM-structured sample_info", {
  meta <- sim_geno_meta(n_samples = 4L, n_markers = 8L, block_size = 4L)
  expect_true(is.data.frame(meta$sample_info))
  expect_equal(names(meta$sample_info),
               c("Family", "Sample", "Parent1", "Parent2", "Sex"))
  expect_equal(meta$sample_info$Sample, dimnames(meta$geno_arrays)[[1L]])
  # non-Sample columns are dummy-filled with "0"
  expect_true(all(meta$sample_info$Family  == "0"))
  expect_true(all(meta$sample_info$Parent1 == "0"))
  expect_true(all(meta$sample_info$Parent2 == "0"))
  expect_true(all(meta$sample_info$Sex     == "0"))
})

test_that("new_geno_meta reads a positional fam into sample_info", {
  meta  <- sim_geno_meta(n_samples = 4L, n_markers = 8L, block_size = 4L)
  samp  <- dimnames(meta$geno_arrays)[[1L]]
  fam   <- data.frame(fid = c("A","A","B","B"), iid = samp,
                      pa = c("0","0","0","0"), ma = c("0","0","0","0"),
                      sx = c("1","2","1","2"), stringsAsFactors = FALSE)
  meta2 <- new_geno_meta(meta$geno_arrays, meta$map_info,
                         block_info            = meta$block_info,
                         block_marker_map_info = meta$block_marker_map_info,
                         metadata              = meta$metadata,
                         sample_info           = fam)
  expect_equal(names(meta2$sample_info),
               c("Family", "Sample", "Parent1", "Parent2", "Sex"))
  expect_equal(meta2$sample_info$Family, c("A","A","B","B"))
  expect_equal(meta2$sample_info$Sex,    c("1","2","1","2"))
})

test_that("is_geno_meta rejects a misaligned sample_info", {
  meta <- .test_data$meta
  meta$sample_info <- meta$sample_info[-1L, , drop = FALSE]
  expect_false(suppressMessages(is_geno_meta(meta, verbose = FALSE)))
})

test_that("new_geno_meta errors on sample_info not matching genotype samples", {
  meta   <- .test_data$meta
  bad_si <- data.frame(Sample = paste0("X", seq_len(nrow(meta$sample_info))),
                       stringsAsFactors = FALSE)
  expect_error(
    new_geno_meta(meta$geno_arrays, meta$map_info,
                  block_info            = meta$block_info,
                  block_marker_map_info = meta$block_marker_map_info,
                  metadata              = meta$metadata,
                  sample_info           = bad_si),
    "sample set"
  )
})

# ── filter_geno_meta (supersedes filter_profile) ──────────────────────────────

test_that("filter_geno_meta restricts samples correctly", {
  sub <- filter_geno_meta(.test_data$meta, keep_samples = c("S1","S2"))
  expect_equal(sort(dimnames(sub$geno_arrays)[[1L]]), c("S1","S2"))
  expect_equal(sort(sub$sample_info$Sample), c("S1","S2"))
})

test_that("filter_geno_meta restricts blocks correctly", {
  v      <- .views(.test_data$meta)
  target <- v$block_meta$block_id[1L]
  sub    <- filter_geno_meta(.test_data$meta, blocks = target)
  sub_v  <- .views(sub)
  expect_true(all(sub_v$block_meta$block_id == target))
})

test_that("filter_geno_meta restricts chromosomes correctly", {
  sub   <- filter_geno_meta(.test_data$meta,
                            chromosomes = list(list(chr = "1")))
  sub_v <- .views(sub)
  expect_true(all(sub_v$marker_meta$chromosome == "1"))
})

test_that("filter_geno_meta warns on unknown sample (no error)", {
  expect_warning(filter_geno_meta(.test_data$meta, keep_samples = "NOEXIST"),
                 "not present")
})

test_that("filter_geno_meta warns on unknown marker (no error)", {
  expect_warning(filter_geno_meta(.test_data$meta, keep_markers = "NOEXIST"),
                 "not present")
})

test_that("filter_geno_meta returns a copy, not modifying the original", {
  n_orig <- length(dimnames(.test_data$meta$geno_arrays)[[1L]])
  sub    <- filter_geno_meta(.test_data$meta, keep_samples = "S1")
  expect_equal(length(dimnames(.test_data$meta$geno_arrays)[[1L]]), n_orig)
})

# ── save / load GenoMeta ──────────────────────────────────────────────────────

test_that("save_geno_meta and load_geno_meta round-trip correctly", {
  meta <- .test_data$meta
  tmp  <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  save_geno_meta(meta, tmp)
  meta2 <- load_geno_meta(tmp)
  expect_equal(meta$geno_arrays,           meta2$geno_arrays)
  expect_equal(meta$metadata$encoding,     meta2$metadata$encoding)
  expect_equal(meta$sample_info,           meta2$sample_info)
})

test_that("load_geno_meta errors on non-existent file", {
  expect_error(load_geno_meta("/no/such/file.rds"), "not found")
})

# ── recode_geno_meta ──────────────────────────────────────────────────────────

test_that("recode_geno_meta converts selectiontools to nucleotide", {
  rec <- recode_geno_meta(.test_data$meta, to_type = "nucleotide")
  expect_equal(rec$metadata$encoding, "nucleotide")
  expect_true(all(stats::na.omit(as.character(rec$geno_arrays)) %in%
                    c("A","G","C","T","-")))
})

test_that("recode_geno_meta returns input unchanged when encoding matches", {
  rec <- recode_geno_meta(.test_data$meta, to_type = "nucleotide")
  expect_identical(rec, .test_data$meta)
  # the old target is defunct (0.9.0, C-5)
  expect_error(recode_geno_meta(.test_data$meta, to_type = "selectiontools"),
               class = "hs_error_defunct")
})

test_that("nucleotide<->selectiontools codes are alphabetical A=1,C=2,G=3,T=4", {
  expect_equal(unlist(DEFAULT_ENCODING_DICT[c("A", "C", "G", "T", "-")]),
               c(A = 1L, C = 2L, G = 3L, T = 4L, `-` = -1L))
  # inverse: codes 2/3 map back to C/G (the fixed swap)
  expect_equal(unname(REVERSED_ENCODING_DICT[c("1", "2", "3", "4")]),
               c("A", "C", "G", "T"))
  # end-to-end: codes 1/2 and 3/4 declared as SelectionTools translate to
  # A/C and G/T at the portal (deprecated since 0.9.0, C-5; behind the switch, C-6)
  withr::local_options(list(haploStack.st_formats = TRUE))
  expect_warning(
    m <- suppressMessages(build_geno_meta(
      geno_2d = matrix(c("1/2", "3/4"), nrow = 1,
                       dimnames = list("S1", c("M1", "M2"))),
      map     = data.frame(Chrom = "1", Marker = c("M1", "M2"), Pos = c(1e6, 2e6)),
      encoding = "selectiontools")),
    class = "hs_warning_deprecated")
  expect_setequal(as.character(m$geno_arrays["S1", "M1", ]), c("A", "C"))
  expect_setequal(as.character(m$geno_arrays["S1", "M2", ]), c("G", "T"))
})

# ── compute_profile_localgebv ─────────────────────────────────────────────────

test_that("compute_profile_localgebv returns wide data.table with Block column", {
  gebv <- compute_profile_localgebv(.test_data$meta, .test_data$eff)
  expect_true(is.data.table(gebv))
  expect_true("Block" %in% names(gebv))
})

test_that("compute_profile_localgebv dimensions: n_blocks rows, n_samples+1 cols", {
  v         <- .views(.test_data$meta)
  gebv      <- compute_profile_localgebv(.test_data$meta, .test_data$eff)
  n_blocks  <- nrow(v$block_meta)
  n_samples <- length(v$samples)
  expect_equal(nrow(gebv), n_blocks)
  expect_equal(ncol(gebv), n_samples + 1L)  # +1 for Block column
})

test_that("compute_profile_localgebv haplotype level returns 2x sample columns", {
  v         <- .views(.test_data$meta)
  gebv      <- compute_profile_localgebv(.test_data$meta, .test_data$eff,
                                         level = "haplotype")
  n_blocks  <- nrow(v$block_meta)
  n_samples <- length(v$samples)
  expect_equal(nrow(gebv), n_blocks)
  expect_equal(ncol(gebv), n_samples * 2L + 1L)
})

test_that("compute_profile_localgebv fill_na fills missing blocks", {
  empty_eff <- data.frame(Marker = character(0L), Allele = character(0L),
                           Effect = numeric(0L))
  gebv        <- compute_profile_localgebv(.test_data$meta, empty_eff,
                                           fill_na = -999)
  sample_cols <- setdiff(names(gebv), "Block")
  vals        <- as.numeric(unlist(gebv[, sample_cols, with = FALSE]))
  expect_true(all(vals == -999 | is.na(vals) | vals == 0))
})

test_that("compute_profile_localgebv block order matches block_meta", {
  v           <- .views(.test_data$meta)
  gebv        <- compute_profile_localgebv(.test_data$meta, .test_data$eff)
  block_order <- v$block_meta[order(chromosome, start_pos), block_id]
  returned    <- gebv$Block
  expect_equal(returned, block_order[block_order %in% returned])
})

test_that("compute_profile_localgebv errors on invalid GenoMeta", {
  expect_error(compute_profile_localgebv(list(), data.frame()), "GenoMeta")
})

test_that("compute_profile_localgebv falls back to stored marker_effects", {
  meta <- .test_data$meta
  meta$marker_effects <- .test_data$eff
  expect_true(is_geno_meta(meta))
  gebv <- compute_profile_localgebv(meta)        # effect_df = NULL
  expect_true(is.data.table(gebv))
  expect_true("Block" %in% names(gebv))
})

test_that("compute_profile_localgebv errors when no effects are available", {
  expect_error(compute_profile_localgebv(.test_data$meta), "no effects")
})

# ── marker_effects (optional stored copy) ─────────────────────────────────────

test_that("marker_effects is pruned by filter_geno_meta and survives save/load", {
  meta <- .test_data$meta
  meta$marker_effects <- .test_data$eff
  expect_true(is_geno_meta(meta))

  # filter_geno_meta normalises marker_effects to the canonical per-trait list.
  sub <- filter_geno_meta(meta, chromosomes = list(list(chr = "1")))
  expect_true(all(sub$marker_effects[[1L]]$Marker %in% sub$map_info$Marker))

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  save_geno_meta(meta, tmp)
  meta2 <- load_geno_meta(tmp)
  expect_equal(meta$marker_effects, meta2$marker_effects)
})

# ── build_block_location_info ─────────────────────────────────────────────────

test_that("build_block_location_info returns required columns", {
  result <- build_block_location_info(.test_data$meta)
  info   <- result$block_info
  expect_true(all(c("Chrom","Block","Start","End","Pos",
                     "Num_SNPs","Block_length") %in% names(info)))
})

test_that("build_block_location_info Block_start <= Block_end", {
  info <- build_block_location_info(.test_data$meta)$block_info
  expect_true(all(info$Start <= info$End))
})

test_that("build_block_location_info Num_SNPs matches block_map", {
  info     <- build_block_location_info(.test_data$meta)$block_info
  expected <- as.data.table(.test_data$bmap)[, .N, by = Block]
  for (blk in info$Block) {
    expect_equal(
      info[Block == blk, Num_SNPs],
      expected[Block == blk, N]
    )
  }
})

# ── Manhattan plots (regression: .create_genomic_positions name collision) ─────

test_that("plot_manhattan_variance returns a ggplot on multi-chrom data", {
  dt <- data.table(Chrom = c("1","1","2","2"), Pos = c(10,20,30,40),
                   BlockVariance = c(1,2,3,4))
  p  <- plot_manhattan_variance(dt)
  expect_s3_class(p, "ggplot")
})

test_that("plot_manhattan_effect returns a ggplot on multi-chrom data", {
  dt <- data.table(Chrom = c("1","1","2","2"), Pos = c(10,20,30,40),
                   Effect = c(-0.5, 0.2, 0.0, 0.7))
  p  <- plot_manhattan_effect(dt)
  expect_s3_class(p, "ggplot")
})

# ── plot_localgebv_alignment ──────────────────────────────────────────────────

# A GenoMeta carrying both block_info and a traced family history. The pedigree
# makes samp[1],samp[2] founders of an F_1 with offspring samp[3],samp[4].
.lg_meta <- local({
  meta <- sim_geno_meta(n_samples = 8L, n_markers = 12L,
                        n_chrom = 2L, block_size = 3L, seed = 7L)
  samp <- dimnames(meta$geno_arrays)[[1L]]
  fam  <- data.frame(
    IID    = samp,
    Father = c(NA, NA, samp[1], samp[1], samp[3], samp[3], samp[5], samp[5]),
    Mother = c(NA, NA, samp[2], samp[2], samp[4], samp[4], samp[6], samp[6]),
    stringsAsFactors = FALSE)
  attach_family_history(meta, fam, parent_cols = c("Father", "Mother"))
})

test_that("plot_localgebv_alignment returns a nested gen -> family ggplot list", {
  lg <- compute_profile_localgebv(.lg_meta, .test_data$eff)
  p  <- plot_localgebv_alignment(lg, .lg_meta)

  expect_type(p, "list")
  expect_true(all(grepl("^gen_", names(p))))
  expect_true(all(grepl("^F_", names(p[["gen_1"]]))))
  expect_s3_class(p[["gen_1"]][[1L]], "ggplot")
})

test_that("plot_localgebv_alignment puts offspring rows above founder rows", {
  lg  <- compute_profile_localgebv(.lg_meta, .test_data$eff)
  p   <- plot_localgebv_alignment(lg, .lg_meta)

  gen <- as.data.table(.lg_meta$metadata$family_history[["gen_1"]])
  fid <- names(p[["gen_1"]])[1L]
  off <- unique(gen$IID[gen$FID == fid & gen$Label == "O"])
  fnd <- unique(gen$IID[gen$FID == fid & gen$Label == "P"])

  lv  <- levels(p[["gen_1"]][[fid]]$data$y_label)  # level 1 = bottom row
  expect_true(tail(lv, 1L) %in% off)               # top row is an offspring
  expect_true(head(lv, 1L) %in% fnd)               # bottom row is a founder
})

test_that("plot_localgebv_alignment honours selected_blocks", {
  lg     <- compute_profile_localgebv(.lg_meta, .test_data$eff)
  keep   <- lg$Block[1:2]
  p      <- plot_localgebv_alignment(lg, .lg_meta, selected_blocks = keep)
  blocks <- sort(unique(as.character(p[["gen_1"]][[1L]]$data$Block)))
  expect_equal(blocks, sort(keep))
})

test_that("plot_localgebv_alignment auto-detects haplotype-level columns", {
  lg <- compute_profile_localgebv(.lg_meta, .test_data$eff, level = "haplotype")
  p  <- plot_localgebv_alignment(lg, .lg_meta)
  # each genotype sample expands to two haplotype rows ("<Sample>.<idx>")
  expect_true(all(grepl("\\.\\d+$",
                        levels(p[["gen_1"]][[1L]]$data$y_label))))
})

test_that("plot_localgebv_alignment warns when given a bare family-history list", {
  lg <- compute_profile_localgebv(.lg_meta, .test_data$eff)
  expect_warning(
    plot_localgebv_alignment(lg, .lg_meta$metadata$family_history),
    "chromosome labels"
  )
})

test_that("plot_localgebv_alignment errors on a non-localgebv input", {
  expect_error(
    plot_localgebv_alignment(data.frame(x = 1), .lg_meta),
    "Block"
  )
})

# ── compute_block_catalogue ───────────────────────────────────────────────────

# Identification-first column order, with the subtype signature always present.
.catalogue_cols <- c("Chrom", "Pos", "Block", "Sample", "Haplotype", "signature",
                     "Effect", "Count", "Block_start", "Block_end",
                     "Num_SNPs_in_Block", "Block_length", "Block_chr")

test_that("compute_block_catalogue uses id-first order with a signature column", {
  cat_dt <- compute_block_catalogue(.test_data$meta, .test_data$eff)
  expect_true(is.data.table(cat_dt))
  expect_identical(names(cat_dt), .catalogue_cols)
  # identity columns lead; measures come after
  expect_true(which(names(cat_dt) == "Chrom") < which(names(cat_dt) == "Effect"))
  expect_true(which(names(cat_dt) == "Block") < which(names(cat_dt) == "Count"))
  expect_type(cat_dt$signature, "character")
})

test_that("compute_block_catalogue signature = representative carrier's combo", {
  cat_dt <- compute_block_catalogue(.test_data$meta, .test_data$eff)
  sigs   <- haploStack:::.compute_block_signatures(.test_data$meta)
  # each row's signature must equal the signature of its representative carrier
  chk <- merge(cat_dt[, .(Block, Sample, Haplotype, signature)],
               sigs[, .(Block, Sample, Haplotype = Haplotype_idx,
                        ref = signature)],
               by = c("Block", "Sample", "Haplotype"), all.x = TRUE)
  chk[is.na(ref), ref := ""]
  expect_equal(chk$signature, chk$ref)
})

test_that("compute_block_catalogue Count sums to n_samples * ploidy per block", {
  meta   <- .test_data$meta
  cat_dt <- compute_block_catalogue(meta, .test_data$eff)
  nhap   <- dim(meta$geno_arrays)[1L] * (meta$metadata$ploidy %||% 2L)
  per_block <- cat_dt[, .(tot = sum(Count)), by = Block]
  expect_true(all(per_block$tot == nhap))
  # one row per distinct effect state; Haplotype is the representative index
  expect_true(all(cat_dt$Haplotype %in% seq_len(meta$metadata$ploidy %||% 2L)))
  expect_true(all(cat_dt$Num_SNPs_in_Block >= 1L))
})

test_that("compute_block_catalogue errors on a non-GenoMeta", {
  expect_error(compute_block_catalogue(data.frame(x = 1), .test_data$eff),
               "GenoMeta")
})

test_that("compute_block_catalogue group_by='effect' is the default", {
  a <- compute_block_catalogue(.test_data$meta, .test_data$eff)
  b <- compute_block_catalogue(.test_data$meta, .test_data$eff,
                               group_by = "effect")
  expect_identical(a, b)
  expect_true("signature" %in% names(a))   # signature is now always present
})

test_that("compute_block_catalogue group_by='signature' keeps the id-first order", {
  cat_dt <- compute_block_catalogue(.test_data$meta, .test_data$eff,
                                    group_by = "signature")
  expect_true("signature" %in% names(cat_dt))
  expect_type(cat_dt$signature, "character")
  # same identification-first column order as the default path
  expect_identical(names(cat_dt), .catalogue_cols)
})

test_that("compute_block_catalogue signature: Count still sums to n*ploidy", {
  meta   <- .test_data$meta
  cat_dt <- compute_block_catalogue(meta, .test_data$eff, group_by = "signature")
  nhap   <- dim(meta$geno_arrays)[1L] * (meta$metadata$ploidy %||% 2L)
  per_block <- cat_dt[, .(tot = sum(Count)), by = Block]
  expect_true(all(per_block$tot == nhap))
})

test_that("compute_block_catalogue signature is >= as fine-grained as effect", {
  eff_dt <- compute_block_catalogue(.test_data$meta, .test_data$eff,
                                    group_by = "effect")
  sig_dt <- compute_block_catalogue(.test_data$meta, .test_data$eff,
                                    group_by = "signature")
  n_eff <- eff_dt[, .(n_eff = .N), by = Block]
  n_sig <- sig_dt[, .(n_sig = .N), by = Block]
  m     <- merge(n_eff, n_sig, by = "Block")
  # one allele combination -> one effect, so signatures never coarser than effects
  expect_true(all(m$n_sig >= m$n_eff))
  # each signature maps to exactly one Effect value within its block
  per <- sig_dt[, .(n_distinct_eff = uniqueN(round(Effect, 10L))),
                by = .(Block, signature)]
  expect_true(all(per$n_distinct_eff == 1L))
})

test_that("compute_block_catalogue rejects an unknown group_by", {
  expect_error(
    compute_block_catalogue(.test_data$meta, .test_data$eff, group_by = "nope"),
    "should be one of")
})

test_that("run_haplotype_catalogue writes the Haplotype_Catalogue TSV", {
  out <- withr::local_tempdir()
  res <- run_haplotype_catalogue(
    geno_meta   = .test_data$meta,
    effect_list = list(T1 = .test_data$eff),
    output_dir  = out,
    ld_threshold = 0.5, tolerance = 3L,
    save_profile = FALSE)

  f <- file.path(out, "Haplotype_Catalogue_0.5_t3_T1.txt")
  expect_true(file.exists(f))
  hdr <- names(data.table::fread(f, nrows = 0L))
  expect_identical(hdr[1:5], c("Chrom", "Pos", "Block", "Sample", "Haplotype"))
  expect_true("signature" %in% hdr)
  expect_true(is.data.table(res$T1$hap_catalogue))
})

test_that("run_haplotype_catalogue writes the per-block Block_Summary table", {
  out <- withr::local_tempdir()
  res <- run_haplotype_catalogue(
    geno_meta   = .test_data$meta,
    effect_list = list(T1 = .test_data$eff),
    output_dir  = out,
    ld_threshold = 0.5, tolerance = 3L,
    save_profile = FALSE)

  # Returned element carries the summary columns …
  bs <- res$T1$block_summary
  expect_true(is.data.table(bs))
  expect_true(all(c("Block", "BlockVariance", "BlockMax", "BlockMin",
                    "BlockMean", "Chrom", "Pos") %in% names(bs)))
  # identification columns lead the summary
  expect_identical(names(bs)[1:3], c("Chrom", "Pos", "Block"))
  expect_true(all(bs$BlockMax >= bs$BlockMean))
  expect_true(all(bs$BlockMean >= bs$BlockMin))

  # … and the same summary is written to disk, one file per trait.
  f <- file.path(out, "Block_Summary_0.5_t3_T1.txt")
  expect_true(file.exists(f))
  hdr <- names(data.table::fread(f, nrows = 0L))
  expect_true(all(c("BlockVariance", "BlockMax", "BlockMin", "BlockMean")
                  %in% hdr))
})

test_that("plot_manhattan_variance/effect share the genomic-scatter body", {
  # Both wrappers must return a ggplot; effect draws the zero reference line,
  # variance does not (only difference beyond the value column).
  dt_v <- data.table(Chrom = c("1","1","2"), Pos = c(10,20,30),
                     BlockVariance = c(1,2,3))
  dt_e <- data.table(Chrom = c("1","1","2"), Pos = c(10,20,30),
                     Effect = c(-0.5, 0.2, 0.7))
  expect_s3_class(plot_manhattan_variance(dt_v), "ggplot")
  expect_s3_class(plot_manhattan_effect(dt_e),   "ggplot")
})
