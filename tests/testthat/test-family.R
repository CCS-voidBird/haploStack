# =============================================================================
# test-family.R
#
# QA tests for the family-structure tools:
#   trace_family_history(), attach_family_history(),
#   family_history_to_pedigree(), and the bridge to trace_haplotype_family().
# =============================================================================

library(testthat)
library(haploStack)

# 3-generation pedigree:
#   founders  F1..F4  (no parents)
#   gen 1     O1,O2 = F1xF2 ;  O3 = F3xF4
#   gen 2     G1    = O1xO3
.fam <- data.frame(
  IID    = c("F1","F2","F3","F4","O1","O2","O3","G1"),
  Father = c(NA,  NA,  NA,  NA,  "F1","F1","F3","O1"),
  Mother = c(NA,  NA,  NA,  NA,  "F2","F2","F4","O3"),
  stringsAsFactors = FALSE
)

# ── trace_family_history ──────────────────────────────────────────────────────

test_that("trace_family_history walks oldest -> newest across generations", {
  h <- trace_family_history(.fam, parent_cols = c("Father","Mother"))
  expect_type(h, "list")
  expect_length(h, 2L)                       # gen 1 + gen 2
  expect_named(h, c("gen_1","gen_2"))

  # Generation 1: two families, offspring O1/O2/O3
  g1 <- h[["gen_1"]]
  expect_setequal(unique(g1$IID[g1$Label == "O"]), c("O1","O2","O3"))
  expect_setequal(unique(g1$IID[g1$Label == "P"]), c("F1","F2","F3","F4"))

  # Generation 2: G1 from parents O1, O3
  g2 <- h[["gen_2"]]
  expect_setequal(unique(g2$IID[g2$Label == "O"]), "G1")
  expect_setequal(unique(g2$IID[g2$Label == "P"]), c("O1","O3"))
})

test_that("trace_family_history groups full-sib offspring into one family", {
  h  <- trace_family_history(.fam, parent_cols = c("Father","Mother"))
  g1 <- h[["gen_1"]]
  # O1 and O2 share parents F1,F2 -> same FID
  fid_o1 <- g1$FID[g1$IID == "O1" & g1$Label == "O"]
  fid_o2 <- g1$FID[g1$IID == "O2" & g1$Label == "O"]
  fid_o3 <- g1$FID[g1$IID == "O3" & g1$Label == "O"]
  expect_equal(fid_o1, fid_o2)
  expect_false(identical(fid_o1, fid_o3))
})

test_that("trace_family_history accepts a single parent column", {
  h <- trace_family_history(.fam, parent_cols = "Father")
  expect_length(h, 2L)
  expect_setequal(unique(h[["gen_1"]]$IID[h[["gen_1"]]$Label == "O"]),
                  c("O1","O2","O3"))
})

test_that("trace_family_history accepts parent columns by index", {
  h_name <- trace_family_history(.fam, parent_cols = c("Father","Mother"))
  h_idx  <- trace_family_history(.fam, parent_cols = c(2L, 3L))
  expect_equal(h_idx, h_name)
})

test_that("trace_family_history resolves IDs from the first column (no IID col)", {
  fam2 <- .fam
  names(fam2)[1] <- "Sample"          # no "IID" column -> first-column fallback
  h <- trace_family_history(fam2, parent_cols = c("Father","Mother"))
  expect_length(h, 2L)
})

test_that("trace_family_history returns empty list when no families observable", {
  founders_only <- data.frame(IID = c("A","B"), Father = c(NA, NA),
                              Mother = c(NA, NA), stringsAsFactors = FALSE)
  expect_length(trace_family_history(founders_only,
                                     parent_cols = c("Father","Mother")), 0L)
})

# ── attach_family_history ─────────────────────────────────────────────────────

test_that("attach_family_history stores supplementary metadata and keeps validity", {
  meta <- sim_geno_meta(n_samples = 8L, n_markers = 12L, block_size = 3L)
  # Map the simulated sample IDs onto a pedigree
  samp <- dimnames(meta$geno_arrays)[[1L]]
  fam  <- data.frame(
    IID    = samp,
    Father = c(NA, NA, samp[1], samp[1], samp[3], samp[3], samp[5], samp[5]),
    Mother = c(NA, NA, samp[2], samp[2], samp[4], samp[4], samp[6], samp[6]),
    stringsAsFactors = FALSE
  )
  meta <- attach_family_history(meta, fam, parent_cols = c("Father","Mother"))

  expect_true(is.list(meta$metadata$family_history))
  expect_true(length(meta$metadata$family_history) >= 1L)
  expect_true(is_geno_meta(meta))           # supplementary, not validated away

  # Write-back: offspring samples get a generation-qualified Family ID
  fam_col <- meta$sample_info$Family[meta$sample_info$Sample %in% samp[3:8]]
  expect_true(any(grepl("^gen_", fam_col)))

  # Survives a filter round-trip (carried in metadata)
  sub <- filter_geno_meta(meta, keep_samples = samp[1:6])
  expect_false(is.null(sub$metadata$family_history))
})

test_that("attach_family_history traces the object's own sample_info by default", {
  meta <- sim_geno_meta(n_samples = 6L, n_markers = 10L, block_size = 5L)
  samp <- dimnames(meta$geno_arrays)[[1L]]
  # Populate sample_info parents (cols 3,4) directly: S3..S6 are S1xS2 offspring
  meta$sample_info$Parent1 <- c("0","0", samp[1], samp[1], samp[1], samp[1])
  meta$sample_info$Parent2 <- c("0","0", samp[2], samp[2], samp[2], samp[2])
  expect_true(is_geno_meta(meta))

  meta <- attach_family_history(meta)          # no external fam → uses sample_info
  expect_true(length(meta$metadata$family_history) >= 1L)
  expect_true(is_geno_meta(meta))
  # The four offspring now carry a non-"0" Family ID
  off_fam <- meta$sample_info$Family[meta$sample_info$Sample %in% samp[3:6]]
  expect_true(all(off_fam != "0"))
})

test_that("attach_family_history warns when no individual matches the genotypes", {
  meta <- sim_geno_meta(n_samples = 4L, n_markers = 8L, block_size = 4L)
  expect_warning(
    attach_family_history(meta, .fam, parent_cols = c("Father","Mother")),
    "no family-history individual"
  )
})

# ── family_history_to_pedigree + bridge to trace_haplotype_family ─────────────

test_that("family_history_to_pedigree yields Sample/Parent1/Parent2", {
  h   <- trace_family_history(.fam, parent_cols = c("Father","Mother"))
  ped <- family_history_to_pedigree(h)
  expect_true(all(c("Sample","Parent1","Parent2") %in% names(ped)))
  expect_setequal(ped$Sample, c("O1","O2","O3","G1"))
  expect_setequal(unlist(ped[ped$Sample == "G1", c("Parent1","Parent2")]),
                  c("O1","O3"))
})

test_that("family_history_to_pedigree drives trace_haplotype_family end to end", {
  meta <- sim_geno_meta(n_samples = 8L, n_markers = 12L, block_size = 3L)
  samp <- dimnames(meta$geno_arrays)[[1L]]
  fam  <- data.frame(
    IID    = samp,
    Father = c(NA, NA, samp[1], samp[1], samp[3], samp[3], samp[5], samp[5]),
    Mother = c(NA, NA, samp[2], samp[2], samp[4], samp[4], samp[6], samp[6]),
    stringsAsFactors = FALSE
  )
  meta <- attach_family_history(meta, fam, parent_cols = c("Father","Mother"))
  ped  <- family_history_to_pedigree(meta)
  fam_align <- trace_haplotype_family(meta, ped)
  expect_s3_class(fam_align, "data.table")
  expect_true(all(c("Sample","Block","Parent","n_match","n_total") %in%
                    names(fam_align)))
})
