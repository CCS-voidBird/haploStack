# =============================================================================
# test-block-annotation-parity.R — BAN (#033) stage-2 acceptance: the golden
# parity contract (Q25.c). At flank_bp = 0, export_block_annotation() must
# reproduce the Python block_annotation 0.3.0 step files line for line.
#
# Two references:
#   P1  the fixture shipped in inst/extdata/parity (0.9.0, decision-form row
#       36 b; design D5-10): synthetic markers laid on the public MorexV3
#       1H + 2H gene rows, the Python reference re-cut once on them
#       (1 195 blocks / 568 with >= 1 gene / 5 397 records / 63 pairs). Runs
#       everywhere, from an installed build too.
#   P2  the 2026-09-20 project run (GRDC barley AGG1, 1H+2H, 300 lines,
#       HapSelect blocks at r2 >= 0.5: 1 774 / 228 / 3 111 / 86), whose inputs
#       live outside the package, so it is skip-gated on their presence.
# =============================================================================

.parity_fixture_dir <- function() {
  d <- system.file("extdata", "parity", package = "haploStack")
  if (!nzchar(d) || !file.exists(file.path(d, "genes_1H2H.gff3.gz")))
    d <- testthat::test_path("..", "..", "inst", "extdata", "parity")
  d
}

.parity_container <- function(map, bmm) {
  n_mk <- nrow(map)
  arr  <- array("A", dim = c(2L, n_mk, 2L),
                dimnames = list(c("S1", "S2"), map$Marker, c("1", "2")))
  meta <- new_geno_meta(
    geno_arrays = arr, map_info = map,
    metadata = list(ploidy = 2L, char_dict = c("A" = 1L), n_samples = 2L,
                    n_markers = n_mk, source_format = "nucleotide"))
  attach_block_info(meta, as.data.frame(bmm))
}

.parity_compare <- function(od, ref_dir, ref_suffix = "") {
  for (nm in c("step1_block_coordinates.tsv", "step1_block_marker_map.tsv",
               "step3_local_annotation.tsv",
               "step3_local_annotation_records.tsv",
               "step4_deep_annotation.tsv")) {
    got  <- readLines(file.path(od, nm), warn = FALSE)
    rf   <- file.path(ref_dir, paste0(nm, ref_suffix))
    want <- if (endsWith(rf, ".gz")) readLines(gzfile(rf), warn = FALSE)
            else readLines(rf, warn = FALSE)
    expect_identical(got, want, label = paste0("lines of ", nm))
  }
}

test_that("P1: the shipped fixture reproduces the Python reference step files (un-gated)", {
  d <- .parity_fixture_dir()
  expect_true(file.exists(file.path(d, "genes_1H2H.gff3.gz")))
  snp <- data.table::fread(file.path(d, "snp_map.txt.gz"), colClasses = "character")
  bmm <- data.table::fread(file.path(d, "block_marker_map.csv.gz"), colClasses = "character")
  map <- data.frame(Chrom = snp$chr, Marker = snp$name, Pos = as.numeric(snp$pos),
                    stringsAsFactors = FALSE)
  meta <- .parity_container(map, bmm)

  genes <- read_gff3_genes(file.path(d, "genes_1H2H.gff3.gz"))
  expect_identical(nrow(genes), 10245L)
  ann <- compute_block_annotation(meta, genes, chr_map = c("1" = "1H", "2" = "2H"),
                                  flank_bp = 0,
                                  browser_base = "https://plants.ensembl.org/Hordeum_vulgare")
  expect_identical(nrow(ann$blocks), 1195L)
  expect_identical(sum(ann$blocks$gene_count > 0L), 568L)
  expect_identical(nrow(ann$features), 5397L)
  expect_identical(nrow(ann$marker_gene), 63L)

  meta <- attach_block_annotation(meta, ann)
  od <- withr::local_tempdir()
  export_block_annotation(meta, od)
  .parity_compare(od, d, ".gz")
})


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

test_that("P2: the 2026-09-20 GRDC run (machine-gated) is reproduced", {
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

  meta <- .parity_container(map, bmm)

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

  .parity_compare(od, file.path(p$ref, "annotation_out"))
})
