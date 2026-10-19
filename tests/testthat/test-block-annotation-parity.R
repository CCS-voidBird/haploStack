# =============================================================================
# test-block-annotation-parity.R — BAN (#033) stage-2 acceptance: the golden
# parity contract (Q25.c). At flank_bp = 0, export_block_annotation() must
# reproduce the Python block_annotation 0.3.0 step files line-for-line on the
# 2026-09-20 reference run (GRDC barley AGG1, 1H+2H, 300 lines, HapSelect
# blocks at r² >= 0.5): 1,774 blocks / 228 with >= 1 gene / 3,111 (block, gene)
# records / 86 (marker, gene) pairs.
#
# The reference inputs and outputs live OUTSIDE the package (the repo's
# ready_to_use/workshop_products_0.7.7/ and the GRDC GFF3), so this test is
# skip-gated on their presence: it runs on the project machine and skips in a
# plain R CMD check. Re-cut the reference whenever the Python side changes
# (PRD §6, "two implementations of one rule set").
# =============================================================================

.parity_paths <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "..", "..", ".."),
                        mustWork = FALSE)
  list(
    ref  = file.path(root, "ready_to_use", "workshop_products_0.7.7",
                     "07_annotation"),
    gff3 = file.path(root, "data", "GRDC", "hordeum_vulgare",
                     "GRDC_barley_AGG1", "raw", "Barley_genes.gff3")
  )
}

test_that("flank_bp = 0 reproduces the Python reference step files", {
  p <- .parity_paths()
  skip_if(!dir.exists(p$ref), "Python reference run not on this machine")
  skip_if(!file.exists(p$gff3), "MorexV3 GFF3 not on this machine")

  snp <- data.table::fread(file.path(p$ref, "haplostack_out", "snp_map_ba.txt"),
                           colClasses = "character")
  bmm <- data.table::fread(file.path(p$ref, "haplostack_out",
                                     "workshop_block_marker_map.csv"),
                           colClasses = "character")

  map <- data.frame(Chrom = snp$chr, Marker = snp$name,
                    Pos = as.numeric(snp$pos), stringsAsFactors = FALSE)

  # A minimal container: the annotation route needs only map + blocks.
  n_mk <- nrow(map)
  arr  <- array("A", dim = c(2L, n_mk, 2L),
                dimnames = list(c("S1", "S2"), map$Marker, c("1", "2")))
  meta <- new_geno_meta(
    geno_arrays = arr, map_info = map,
    metadata = list(ploidy = 2L, char_dict = c("A" = 1L), n_samples = 2L,
                    n_markers = n_mk, source_format = "nucleotide"))
  meta <- attach_block_info(meta, as.data.frame(bmm))

  genes <- read_gff3_genes(p$gff3)
  ann <- compute_block_annotation(
    meta, genes,
    chr_map      = c("1" = "1H", "2" = "2H", "3" = "3H", "4" = "4H",
                     "5" = "5H", "6" = "6H", "7" = "7H"),
    flank_bp     = 0,
    browser_base = "https://plants.ensembl.org/Hordeum_vulgare")

  # the headline counts of the reference run
  expect_equal(nrow(ann$blocks), 1774L)
  expect_equal(sum(ann$blocks$gene_count > 0L), 228L)
  expect_equal(nrow(ann$features), 3111L)
  expect_equal(nrow(ann$marker_gene), 86L)
  expect_equal(length(unique(ann$marker_gene$gene_id)), 70L)

  meta <- attach_block_annotation(meta, ann)
  od <- file.path(tempdir(), "hs_parity_export")
  unlink(od, recursive = TRUE)
  paths <- export_block_annotation(meta, od)

  ref_file <- function(nm) file.path(p$ref, "annotation_out", nm)
  for (nm in c("step1_block_coordinates.tsv", "step1_block_marker_map.tsv",
               "step3_local_annotation.tsv",
               "step3_local_annotation_records.tsv",
               "step4_deep_annotation.tsv")) {
    got  <- readLines(file.path(od, nm), warn = FALSE)
    want <- readLines(ref_file(nm), warn = FALSE)
    expect_identical(got, want, label = paste0("lines of ", nm))
  }
})
