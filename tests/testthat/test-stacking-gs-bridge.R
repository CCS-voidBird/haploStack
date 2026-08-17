# =============================================================================
# test-stacking-gs-bridge.R
#
# Tests for the GenoMeta <-> genomicSimulation bridge (R/stacking_gs_bridge.R).
# The file/parse helpers run without genomicSimulation; the live-crossing
# round-trip is skipped unless the package is installed.
# =============================================================================

library(testthat)
library(haploStack)

.gsb <- local({
  samples <- c("P1", "P2", "P3")
  markers <- c("m1", "m2", "m3", "m4")
  arr <- array(NA_character_, dim = c(3L, 4L, 2L),
               dimnames = list(samples, markers, c("1", "2")))
  arr["P1", , 1] <- c("A","A","A","A"); arr["P1", , 2] <- c("T","T","T","T")
  arr["P2", , 1] <- c("A","T","A","T"); arr["P2", , 2] <- c("T","A","T","A")
  arr["P3", , 1] <- c("A","A","T","T"); arr["P3", , 2] <- c("T","T","A","A")

  map <- data.frame(Chrom = "1", Marker = markers, Pos = c(10, 20, 30, 40),
                    stringsAsFactors = FALSE)
  bmm <- data.frame(Chrom = "1", Marker = markers,
                    Block = c("B1","B1","B2","B2"), stringsAsFactors = FALSE)
  eff <- data.frame(Marker = rep(markers, each = 2L),
                    Allele = rep(c("A","T"), 4L),
                    Effect = seq_len(8L) / 10, stringsAsFactors = FALSE)

  meta <- new_geno_meta(
    geno_arrays = arr, map_info = map, block_marker_map_info = bmm,
    metadata = list(ploidy = 2L, char_dict = c(A = 1L, T = 2L),
                    n_samples = 3L, n_markers = 4L,
                    source_format = "nucleotide", encoding = "nucleotide"),
    marker_effects = eff)
  list(meta = meta)
})

# ── file / parse helpers (no genomicSimulation needed) ────────────────────────

test_that(".gs_write_files writes geno/map/eff in genomicSimulation layout", {
  f <- haploStack:::.gs_write_files(.gsb$meta, work_dir = tempfile("gsb"))
  expect_true(file.exists(f[["geno"]]))
  expect_true(file.exists(f[["map"]]))
  expect_true(file.exists(f[["eff"]]))

  gl <- readLines(f[["geno"]])
  expect_length(gl, 5L)                         # header + 4 markers (markers as rows)
  cells <- strsplit(gl[2L], "\t")[[1L]]
  expect_match(cells[2L], "^[ACGT]{2}$")        # two-character genotype

  ml <- strsplit(readLines(f[["map"]])[1L], "\t")[[1L]]
  expect_equal(ml, c("Marker", "Chrom", "Pos"))
})

test_that(".gs_positions_to_cM scales bp maps to cM but leaves cM maps alone", {
  # physical (bp) positions -> cM at 1 cM/Mb
  expect_equal(haploStack:::.gs_positions_to_cM(c(1e6, 3.9e7, 9.86e7)),
               c(1, 39, 98.6))
  # already-cM positions (small) pass through unchanged
  expect_equal(haploStack:::.gs_positions_to_cM(c(10, 20, 30, 40)),
               c(10, 20, 30, 40))
  # rate override
  expect_equal(haploStack:::.gs_positions_to_cM(1e6, cM_per_mb = 2), 2)
})

test_that(".gs_write_files emits cM (non-scientific) positions for a bp map", {
  m <- .gsb$meta
  m$map_info$Pos <- c(1e6, 3.9e7, 5.78e7, 9.86e7)      # bp
  f  <- haploStack:::.gs_write_files(m, work_dir = tempfile("gsb"))
  ml <- read.delim(f[["map"]], stringsAsFactors = FALSE)
  expect_equal(ml$Pos, c(1, 39, 57.8, 98.6))            # bp -> cM
  # written plainly, not as "3.9e+07" (gS's parser needs decimal notation)
  expect_false(any(grepl("e\\+?[0-9]", readLines(f[["map"]]))))
})

test_that(".gs_write_files omits eff when no effects are available", {
  m0 <- .gsb$meta
  m0$marker_effects <- NULL
  f <- haploStack:::.gs_write_files(m0, effect_df = NULL,
                                    work_dir = tempfile("gsb"))
  expect_true(is.na(f[["eff"]]))
})

test_that(".gs_genes_to_array parses 2-character genotypes to a 3-D array", {
  gm <- matrix(c("AT","GG","CA","TT"), nrow = 2L,
               dimnames = list(c("m1","m2"), c("s1","s2")))
  arr <- haploStack:::.gs_genes_to_array(gm)
  expect_equal(dim(arr), c(2L, 2L, 2L))                  # samples x markers x 2
  expect_equal(unname(arr["s1","m1", ]), c("A","T"))
  expect_equal(unname(arr["s2","m2", ]), c("T","T"))
  expect_equal(dimnames(arr)[[1L]], c("s1","s2"))
})

# ── update_geno_meta_from_matrix (no genomicSimulation needed) ────────────────

test_that("update_geno_meta_from_matrix round-trips a gene matrix to a GenoMeta", {
  # markers x samples, 2-char alleles — the see.group.gene.data() layout.
  markers <- c("m1", "m2", "m3", "m4")
  gm <- matrix(c("AT","AT","TA","AA",     # P1
                 "AA","TT","AT","TA",     # P2
                 "TT","AA","AA","TT"),    # P3
               nrow = 4L, ncol = 3L,
               dimnames = list(markers, c("P1", "P2", "P3")))

  out <- update_geno_meta_from_matrix(.gsb$meta, gm)

  expect_true(is_geno_meta(out))
  # samples x markers x 2, names/order preserved from the matrix
  expect_equal(dim(out$geno_arrays), c(3L, 4L, 2L))
  expect_equal(dimnames(out$geno_arrays)[[1L]], c("P1", "P2", "P3"))
  expect_equal(dimnames(out$geno_arrays)[[2L]], markers)
  # alleles land in the array in the right cells
  expect_equal(unname(out$geno_arrays["P1", "m1", ]), c("A", "T"))
  expect_equal(unname(out$geno_arrays["P3", "m4", ]), c("T", "T"))
  # metadata sample/marker counts reflect the new matrix
  expect_equal(out$metadata$n_samples, 3L)
  expect_equal(out$metadata$n_markers, 4L)
})

test_that("update_geno_meta_from_matrix reuses the template's marker metadata", {
  gm <- matrix("AT", nrow = 4L, ncol = 2L,
               dimnames = list(c("m1","m2","m3","m4"), c("s1","s2")))
  out <- update_geno_meta_from_matrix(.gsb$meta, gm)

  # map / block-map / effects come from the template, reordered to the matrix
  expect_equal(out$map_info$Marker, c("m1","m2","m3","m4"))
  expect_equal(nrow(out$block_marker_map_info), 4L)
  expect_identical(out$marker_effects, .gsb$meta$marker_effects)
  expect_equal(out$metadata$encoding, .gsb$meta$metadata$encoding)
})

test_that("update_geno_meta_from_matrix validates its inputs", {
  gm <- matrix("AT", nrow = 4L, ncol = 2L,
               dimnames = list(c("m1","m2","m3","m4"), c("s1","s2")))
  expect_error(update_geno_meta_from_matrix(list(), gm), "GenoMeta")
  # matrix without dimnames is rejected
  expect_error(update_geno_meta_from_matrix(.gsb$meta, matrix("AT", 2, 2)),
               "row \\(marker\\) and column \\(sample\\) names")
})

test_that("update_geno_meta_from_matrix populates offspring FAM from cross_plan", {
  markers <- c("m1","m2","m3","m4")
  # three offspring: two full-sibs of P1xP2, one of P1xP3
  gm <- matrix("AT", nrow = 4L, ncol = 3L,
               dimnames = list(markers, c("F01n01","F01n02","F02n01")))
  plan <- data.frame(
    Parent1 = c("P1","P1","P1"),
    Parent2 = c("P2","P2","P3"),
    CrossID = c("F01n01","F01n02","F02n01"),
    stringsAsFactors = FALSE)

  out <- update_geno_meta_from_matrix(.gsb$meta, gm, cross_plan = plan)
  si  <- out$sample_info

  expect_true(is_geno_meta(out))
  # rows aligned to genotype (offspring) order
  expect_equal(si$Sample, c("F01n01","F01n02","F02n01"))
  expect_equal(si$Parent1, c("P1","P1","P1"))
  expect_equal(si$Parent2, c("P2","P2","P3"))
  # one parent-pair == one family; full-sibs share a Family id, the P1xP3 differs
  expect_equal(si$Family[1], si$Family[2])
  expect_false(si$Family[1] == si$Family[3])
})

test_that("update_geno_meta_from_matrix accepts Cross_ID and reorders to genotypes", {
  markers <- c("m1","m2","m3","m4")
  # matrix sample order differs from plan row order — FAM must follow the matrix
  gm <- matrix("AT", nrow = 4L, ncol = 2L,
               dimnames = list(markers, c("F02n01","F01n01")))
  plan <- data.frame(Parent1 = c("P1","P1"), Parent2 = c("P2","P3"),
                     Cross_ID = c("F01n01","F02n01"), stringsAsFactors = FALSE)

  out <- update_geno_meta_from_matrix(.gsb$meta, gm, cross_plan = plan)
  expect_equal(out$sample_info$Sample,  c("F02n01","F01n01"))
  expect_equal(out$sample_info$Parent2, c("P3","P2"))
})

test_that("update_geno_meta_from_matrix errors on unmatched offspring / bad plan", {
  markers <- c("m1","m2","m3","m4")
  gm <- matrix("AT", nrow = 4L, ncol = 2L,
               dimnames = list(markers, c("F01n01","F99n99")))
  plan <- data.frame(Parent1 = "P1", Parent2 = "P2", CrossID = "F01n01",
                     stringsAsFactors = FALSE)
  expect_error(update_geno_meta_from_matrix(.gsb$meta, gm, cross_plan = plan),
               "not found in 'cross_plan'")
  # plan missing parent columns
  bad <- data.frame(CrossID = c("F01n01","F99n99"), stringsAsFactors = FALSE)
  expect_error(update_geno_meta_from_matrix(.gsb$meta, gm, cross_plan = bad),
               "'Parent1' and 'Parent2'")
})

# ── session helpers ───────────────────────────────────────────────────────────

test_that("gs_close is a safe no-op on NULL / non-handles", {
  expect_null(gs_close(NULL))          # no genomicSimulation needed
  expect_null(gs_close(list()))
})

# ── genomicSimulation-dependent paths ─────────────────────────────────────────

test_that("gs_* functions guard cleanly when genomicSimulation is absent", {
  skip_if(requireNamespace("genomicSimulation", quietly = TRUE),
          "genomicSimulation is installed")
  expect_error(gs_load_geno_meta(.gsb$meta), "genomicSimulation")
})

test_that("persistent gS session crosses across cycles with name sync (FR1+FR3)", {
  skip_if_not_installed("genomicSimulation")
  h <- gs_load_geno_meta(.gsb$meta)

  # cycle 1: cross the founders; the handle now points at the offspring group.
  plan1 <- make_cross_design(c("P1","P2","P3"), method = "all_pairs",
                             n_offspring = 2L)
  cp1 <- gs_cross_plan(h, plan1)
  expect_s3_class(cp1$handle, "gs_handle")

  # rename offspring to compact IDs and push them into the live group (FR3).
  old <- cp1$plans$CrossID
  new <- paste0("C1_", seq_along(old))
  gs_sync_group_names(cp1$handle, old, new)

  # cycle 2: cross using the NEW names as parents — resolves only if the sync
  # took effect on the live group (no re-load between cycles).
  plan2 <- make_cross_design(new[1:3], method = "all_pairs", n_offspring = 1L)
  cp2 <- gs_cross_plan(cp1$handle, plan2)
  expect_true(all(c("CrossID", "F1mean") %in% names(cp2$plans)))

  off <- gs_group_to_geno_meta(cp2$handle, cross_plan = cp2$plans)
  expect_true(is_geno_meta(off))
  # offspring FAM carries the cycle-1 parents used for cycle-2 crosses
  expect_true(all(off$sample_info$Parent1 %in% new))
  expect_null(gs_close(cp2$handle))
})

test_that("genomicSimulation cross round-trip (integration)", {
  skip_if_not_installed("genomicSimulation")
  h    <- gs_load_geno_meta(.gsb$meta)
  plan <- make_cross_design(c("P1","P2","P3"), method = "all_pairs",
                            n_offspring = 2L)
  res  <- gs_cross_plan(h, plan)
  expect_true(all(c("CrossID", "F1mean") %in% names(res$plans)))
  expect_s3_class(res$handle, "gs_handle")

  off <- gs_group_to_geno_meta(res$handle)
  expect_true(is_geno_meta(off))
  expect_equal(dimnames(off$geno_arrays)[[2L]],
               dimnames(.gsb$meta$geno_arrays)[[2L]])   # same marker set
})
