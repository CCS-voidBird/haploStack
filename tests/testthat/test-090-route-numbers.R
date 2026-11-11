# =============================================================================
# test-090-route-numbers.R — route-dependent numbers (design D3 §3.7, D3.5)
#
#   R1  every .resolve_center() caller defaults per .HS_CENTER_DEFAULT: calling
#       with no centre equals calling with the table's value; "none" for the
#       seven catalogue / selection entry points, "locus_mean" for
#       count_favourable_blocks(); center_by_freq is defunct (0.9.0)
#   R3  selected_blocks restricts hap_catalogue (cc/F05) via
#       compute_block_catalogue(blocks =)
#   T   count_favourable_blocks(tol =): |value| <= tol is Zero; default = 0.8.8 rule
#   C   the centre travels with the result: provenance or attr "hs_center"
# (R2 = centre check with catalogue is in test-090-analysis.R A4; R4–R6 land
#  with D3.7 / D3.8.)
# =============================================================================

.rn_fixture <- function() {
  ds <- sim_dataset(n_samples = 12L, n_markers = 24L, block_size = 4L, seed = 11L)
  m <- ds$geno_meta
  m$marker_effects <- list(T1 = ds$eff_nuc)
  list(meta = m, eff = ds$eff_nuc)
}
.rn_quiet <- function(expr) suppressMessages(suppressWarnings(expr))
.rn_strip <- function(x) { x <- as.data.frame(x); attr(x, "hs_center") <- NULL; x }

test_that("R1: the eight entry points default per .HS_CENTER_DEFAULT", {
  fx <- .rn_fixture(); m <- fx$meta; eff <- fx$eff
  tab <- haploStack:::.HS_CENTER_DEFAULT
  expect_identical(unname(tab[c("run_haplotype_catalogue", "compute_profile_localgebv",
                                "predict_gebv", "compute_block_catalogue",
                                "analyze_block_catalogue", "count_block_subtypes",
                                "count_favourable_subtypes")]), rep("none", 7L))
  expect_identical(unname(tab["count_favourable_blocks"]), "locus_mean")

  same <- function(f, ..., default) {
    a <- .rn_quiet(f(...)); b <- .rn_quiet(f(..., center = default))
    expect_equal(.rn_strip(a), .rn_strip(b))
    a
  }
  same(compute_profile_localgebv, m, effect_df = eff, default = "none")
  same(compute_block_catalogue, m, effect_df = eff, default = "none")
  same(count_block_subtypes, m, trait = "T1", default = "none")
  same(count_favourable_subtypes, m, trait = "T1", default = "none")
  a <- .rn_quiet(count_favourable_blocks(m, trait = "T1", plot = FALSE))
  b <- .rn_quiet(count_favourable_blocks(m, trait = "T1", plot = FALSE, center = "locus_mean"))
  expect_equal(a$table, b$table); expect_equal(a$dummy, b$dummy)
  expect_identical(hs_provenance(.rn_quiet(run_haplotype_catalogue(m, list(T1 = eff))))$center, "none")
  expect_identical(hs_provenance(.rn_quiet(analyze_block_catalogue(m, list(T1 = eff),
                                                                   plots = character(0))))$center, "none")

  # the 0.8.8 numbers stay reachable, and differ from the new default
  old <- .rn_quiet(compute_block_catalogue(m, effect_df = eff, center = "allele_freq"))
  new <- .rn_quiet(compute_block_catalogue(m, effect_df = eff))
  expect_false(isTRUE(all.equal(old$Effect, new$Effect)))

  # center_by_freq is defunct since 0.9.0 (row 26 a): the error names the mode
  expect_error(compute_profile_localgebv(m, effect_df = eff, center_by_freq = TRUE),
               'center = "allele_freq"', class = "hs_error_defunct", fixed = TRUE)
})

test_that("R3: selected_blocks restricts hap_catalogue through compute_block_catalogue(blocks =)", {
  fx <- .rn_fixture(); m <- fx$meta; eff <- fx$eff
  blocks <- unique(as.character(m$block_marker_map_info$Block))
  keep <- blocks[1:2]
  cb <- .rn_quiet(compute_block_catalogue(m, effect_df = eff, blocks = keep))
  expect_setequal(unique(cb$Block), keep)
  full <- .rn_quiet(compute_block_catalogue(m, effect_df = eff))
  expect_equal(.rn_strip(cb), .rn_strip(full[full$Block %in% keep, ]), ignore_attr = TRUE)
  res <- .rn_quiet(run_haplotype_catalogue(m, list(T1 = eff), selected_blocks = data.frame(Block = keep)))
  expect_setequal(unique(res$T1$hap_catalogue$Block), keep)
  expect_setequal(unique(res$T1$gebv_matrix$Block), keep)
})

test_that("T: count_favourable_blocks(tol =) widens the Zero band; default is the strict sign", {
  fx <- .rn_fixture(); m <- fx$meta
  a <- .rn_quiet(count_favourable_blocks(m, trait = "T1", plot = FALSE))
  b <- .rn_quiet(count_favourable_blocks(m, trait = "T1", plot = FALSE, tol = 0))
  expect_equal(a$table, b$table)
  big <- .rn_quiet(count_favourable_blocks(m, trait = "T1", plot = FALSE, tol = 1e6))
  expect_true(all(big$table$nPositive == 0L) && all(big$table$nNegative == 0L))
  expect_true(all(a$table$nPositive + a$table$nNegative >= big$table$nPositive + big$table$nNegative))
  expect_error(count_favourable_blocks(m, trait = "T1", plot = FALSE, tol = -1), class = "hs_error_input")
  expect_error(count_favourable_blocks(m, trait = "T1", plot = FALSE, tol = c(0, 1)), class = "hs_error_input")
})

test_that("C: the centre used travels with the result", {
  fx <- .rn_fixture(); m <- fx$meta; eff <- fx$eff
  expect_identical(attr(.rn_quiet(compute_block_catalogue(m, effect_df = eff)), "hs_center"), "none")
  expect_identical(attr(.rn_quiet(compute_block_catalogue(m, effect_df = eff, center = "allele_freq")),
                        "hs_center"), "allele_freq")
  expect_identical(attr(.rn_quiet(count_block_subtypes(m, trait = "T1")), "hs_center"), "none")
  expect_identical(attr(.rn_quiet(count_favourable_subtypes(m, trait = "T1")), "hs_center"), "none")
  cfb <- .rn_quiet(count_favourable_blocks(m, trait = "T1", plot = FALSE))
  expect_identical(attr(cfb, "hs_center"), "locus_mean")
  expect_identical(hs_provenance(.rn_quiet(run_haplotype_catalogue(m, list(T1 = eff),
                                                                   center = "locus_mean")))$center,
                   "locus_mean")
})

test_that("R5: trace_haplotype_family() reads both parental copies (row 17 a)", {
  # One block of six markers. P1 is heterozygous everywhere (copy 1 all A, copy 2
  # all T); P2 is C/G. The offspring carries P1's SECOND copy and P2's first.
  mk  <- paste0("M", 1:6)
  arr <- array(NA_character_, dim = c(3L, 6L, 2L), dimnames = list(c("P1", "P2", "O"), mk, NULL))
  arr["P1", , 1L] <- "A"; arr["P1", , 2L] <- "T"
  arr["P2", , 1L] <- "C"; arr["P2", , 2L] <- "G"
  arr["O",  , 1L] <- "T"; arr["O",  , 2L] <- "C"
  map <- data.frame(Chrom = "1", Marker = mk, Pos = 1:6 * 100, stringsAsFactors = FALSE)
  bmm <- data.frame(Chrom = "1", Marker = mk, Block = "b1", stringsAsFactors = FALSE)
  m <- suppressMessages(new_geno_meta(arr, map_info = map, block_marker_map_info = bmm,
                                      metadata = list(source_format = "user_supplied",
                                                      encoding = "nucleotide")))
  ped <- data.frame(Sample = "O", Parent1 = "P1", Parent2 = "P2", stringsAsFactors = FALSE)
  tr <- trace_haplotype_family(m, ped)
  expect_identical(names(tr), c("Sample", "Haplotype_idx", "Block", "Parent", "n_match", "n_total", "ParentHap"))
  h1 <- tr[Haplotype_idx == 1L]; h2 <- tr[Haplotype_idx == 2L]
  expect_identical(h1$Parent, "Parent1"); expect_identical(h1$ParentHap, 2L); expect_identical(h1$n_match, 6L)
  expect_identical(h2$Parent, "Parent2"); expect_identical(h2$ParentHap, 1L); expect_identical(h2$n_match, 6L)
  expect_identical(tr$n_total, c(6L, 6L))
  # 0.8.8 read parent copy 1 only: offspring copy 1 ("T") matched neither P1 copy 1 ("A")
  # nor P2 copy 1 ("C") -> "unknown". The fix is what the test pins.
  # unknown parents: 0 matches, "unknown", NA copy
  tr0 <- trace_haplotype_family(m, data.frame(Sample = "O", Parent1 = NA, Parent2 = NA))
  expect_identical(unique(tr0$Parent), "unknown"); expect_true(all(is.na(tr0$ParentHap)))
  expect_identical(tr0$n_match, c(0L, 0L))
  # a parent absent from the panel behaves as unknown for its side
  tr1 <- trace_haplotype_family(m, data.frame(Sample = "O", Parent1 = "P1", Parent2 = "ghost"))
  expect_identical(tr1[Haplotype_idx == 1L]$Parent, "Parent1")
  expect_identical(tr1[Haplotype_idx == 2L]$Parent, "unknown")
  # ambiguous: both parents' best copies tie with > 0 matches
  arr2 <- arr; arr2["P2", , 1L] <- "T"
  m2 <- suppressMessages(new_geno_meta(arr2, map_info = map, block_marker_map_info = bmm,
                                       metadata = list(source_format = "user_supplied",
                                                       encoding = "nucleotide")))
  tr2 <- trace_haplotype_family(m2, ped)
  expect_identical(tr2[Haplotype_idx == 1L]$Parent, "ambiguous")
  expect_true(is.na(tr2[Haplotype_idx == 1L]$ParentHap))
})

test_that("R9: the catalogue core and the signatures built on shared views are identical (cc/F09)", {
  fx <- .rn_fixture(); m <- fx$meta
  v  <- geno_meta_views(m)
  expect_identical(haploStack:::.catalogue_core(m, views = v), haploStack:::.catalogue_core(m))
  expect_identical(haploStack:::.compute_block_signatures(m, views = v),
                   haploStack:::.compute_block_signatures(m))
})

