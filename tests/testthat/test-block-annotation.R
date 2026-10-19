# =============================================================================
# test-block-annotation.R — BAN (#033): readers, compute semantics, container
# carry-forward. The golden-file parity against the Python reference run lives
# in test-block-annotation-parity.R (skip-gated on the reference paths).
# =============================================================================

# --- fixtures ----------------------------------------------------------------

# A GFF3 exercising: comment lines, an Ensembl-style gene (ID/Name/biotype,
# gene: prefix), a gene with no Name, an NCBI-style gene (gene_id/gene_biotype),
# a non-gene row, and a ##FASTA tail.
.write_toy_gff3 <- function(path, gz = FALSE) {
  lines <- c(
    "##gff-version 3",
    "##sequence-region   1 1 100000000",
    paste("1", "src", "gene", "1000", "2000", ".", "+", ".",
          "ID=gene:G1;Name=GeneOne;biotype=protein_coding", sep = "\t"),
    paste("1", "src", "mRNA", "1000", "2000", ".", "+", ".",
          "ID=transcript:G1.1;Parent=gene:G1", sep = "\t"),
    paste("1", "src", "gene", "5000", "9000", ".", "-", ".",
          "ID=gene:G2;biotype=protein_coding", sep = "\t"),
    paste("2", "src", "gene", "100", "900", ".", "+", ".",
          "gene_id=G3;gene_name=GeneThree;gene_biotype=lncRNA", sep = "\t"),
    "##FASTA",
    ">1",
    "ACGTACGT"
  )
  con <- if (gz) gzfile(path, "wt") else file(path, "wt")
  writeLines(lines, con)
  close(con)
  path
}

.write_toy_xref <- function(path, with_interpro = TRUE) {
  rows <- c(
    "gene_stable_id\txref_id\txref_label\tdescription\tdb_name",
    "G1\tGO:0001\tGO:0001\tfirst process\tGO",
    "G1\tGO:0002\tGO:0002\tsecond process\tGO",
    "G1\tGO:0001\tGO:0001\tfirst process\tGO",             # duplicate — deduped
    if (with_interpro) "G1\tIPR001\tIPR001\tkinase domain\tInterProScan",
    "G1\t12345\tLOC12345\tsome ncbi name\tNCBI gene (formerly EntrezGene)",
    "G1\tP001\tP001\tSwiss description\tUniProtKB/Swiss-Prot",
    "G1\tQ001\tQ001\tTrEMBL description\tUniProtKB/TrEMBL",
    "G2\tQ002\tQ002\tonly TrEMBL\tUniProtKB/TrEMBL"
  )
  writeLines(rows[!vapply(rows, is.null, logical(1L))], path)
  path
}

# A blocked container whose block spans are known, with a GFF3 built against
# them. sim_geno_meta(): 20 markers, 2 chromosomes, 4-marker blocks.
.toy_annotated_inputs <- function() {
  meta <- sim_geno_meta(n_samples = 4L, n_markers = 20L,
                        encoding = "nucleotide", n_chrom = 2L,
                        block_size = 4L, seed = 7L)
  bi <- as.data.frame(build_block_location_info(meta)$block_info)
  bi <- bi[order(bi$Chrom, bi$Start), ]
  b1 <- bi[1L, ]   # first block, chrom 1
  b2 <- bi[2L, ]   # second block, chrom 1

  genes <- data.table::data.table(
    chrom     = c("1", "1", "1", "2"),
    start     = c(b1$Start + 1,                 # inside block 1
                  b1$End - 1,                   # straddles blocks 1 and 2
                  b2$End + 5000,                # 5 kb beyond block 2
                  1),                           # chrom 2, before everything
    end       = c(b1$Start + 2,
                  b2$Start + 1,
                  b2$End + 6000,
                  2),
    strand    = "+",
    gene_id   = c("gene:IN1", "gene:STRAD", "gene:NEAR", "gene:FAR2"),
    gene_name = c("InOne", NA, "NearTwo", "FarTwo"),
    biotype   = "protein_coding"
  )
  list(meta = meta, bi = bi, genes = genes, b1 = b1, b2 = b2)
}


# --- read_gff3_genes ---------------------------------------------------------

test_that("read_gff3_genes parses genes, synonyms and the FASTA tail", {
  gff <- .write_toy_gff3(tempfile(fileext = ".gff3"))
  g <- read_gff3_genes(gff)

  expect_equal(nrow(g), 3L)                       # mRNA and FASTA rows skipped
  expect_equal(g$gene_id, c("gene:G1", "gene:G2", "G3"))  # prefix kept
  expect_equal(g$gene_name, c("GeneOne", NA, "GeneThree"))
  expect_equal(g$biotype, c("protein_coding", "protein_coding", "lncRNA"))
  expect_equal(g$start, c(1000L, 5000L, 100L))
  expect_s3_class(g, "data.table")
  expect_false(is.null(attr(g, "hs_annotation_source")))
})

test_that("read_gff3_genes reads gzipped files and applies chr_map", {
  gff <- .write_toy_gff3(tempfile(fileext = ".gff3.gz"), gz = TRUE)
  g <- read_gff3_genes(gff, chr_map = c("1" = "1H"))
  expect_equal(nrow(g), 3L)
  expect_equal(sort(unique(g$chrom)), c("1H", "2"))
})

test_that("read_gff3_genes errors usefully on a wrong feature type", {
  gff <- .write_toy_gff3(tempfile(fileext = ".gff3"))
  expect_error(read_gff3_genes(gff, feature_type = "exon"),
               "no 'exon' rows")
})


# --- read_gene_xrefs ---------------------------------------------------------

test_that("read_gene_xrefs replicates the load_gene_xrefs collapse rules", {
  xp <- .write_toy_xref(tempfile(fileext = ".tsv"))
  xr <- read_gene_xrefs(xp)

  g1 <- xr[xr$gene_id == "G1", ]
  expect_equal(g1$go_terms,
               "GO:0001:first process; GO:0002:second process")  # dedup, order
  expect_equal(g1$interpro_domains, "IPR001:kinase domain")
  expect_equal(g1$gene_symbol, "LOC12345")
  expect_equal(g1$gene_description, "Swiss description")  # Swiss-Prot > TrEMBL

  g2 <- xr[xr$gene_id == "G2", ]
  expect_equal(g2$gene_description, "only TrEMBL")
  expect_equal(g2$go_terms, "")
})

test_that("read_gene_xrefs warns when the dump has no InterPro", {
  xp <- .write_toy_xref(tempfile(fileext = ".tsv"), with_interpro = FALSE)
  expect_warning(read_gene_xrefs(xp), "no InterPro")
})


# --- compute_block_annotation ------------------------------------------------

test_that("compute_block_annotation: overlap, flank and the marker-block gate", {
  tt <- .toy_annotated_inputs()

  # flank 0 — the parity semantics
  ann0 <- compute_block_annotation(tt$meta, tt$genes, flank_bp = 0)
  bl0  <- as.data.frame(ann0$blocks)
  b1r  <- bl0[bl0$Block == tt$b1$Block, ]
  b2r  <- bl0[bl0$Block == tt$b2$Block, ]

  expect_equal(b1r$gene_count, 2L)            # IN1 + the straddler
  expect_equal(b1r$gene_ids, "gene:IN1; gene:STRAD")
  expect_equal(b1r$gene_names, "InOne; ")     # absent Name joins as ""
  expect_equal(b2r$gene_count, 1L)            # the straddler only
  expect_equal(b2r$gene_ids, "gene:STRAD")

  # the gate: STRAD touches blocks 1 and 2, but each marker-gene pair must lie
  # inside its own block, so no marker of block 1 may be attributed via the
  # block-2 record and vice versa
  mg0 <- as.data.frame(ann0$marker_gene)
  map <- tt$meta$map_info
  for (i in seq_len(nrow(mg0))) {
    blk <- mg0$Block[i]
    bi  <- tt$bi[tt$bi$Block == blk, ]
    pos <- mg0$Pos[i]
    expect_true(pos >= bi$Start && pos <= bi$End)
    expect_true(pos >= mg0$gene_start[i] && pos <= mg0$gene_end[i])
  }

  # flank 10 kb — NEAR (5 kb past block 2) comes in; counts monotone
  ann10 <- compute_block_annotation(tt$meta, tt$genes, flank_bp = 10000)
  bl10  <- as.data.frame(ann10$blocks)
  expect_true(all(bl10$gene_count >= bl0$gene_count))
  expect_true(grepl("gene:NEAR",
                    bl10$gene_ids[bl10$Block == tt$b2$Block], fixed = TRUE))
  # the reported spans never pad
  expect_equal(bl10$Start, bl0$Start)
  expect_equal(bl10$End, bl0$End)

  # provenance records the dial
  expect_equal(ann10$provenance$flank_bp, 10000)
  expect_equal(ann0$provenance$flank_bp, 0)
})

test_that("compute_block_annotation warns on a vocabulary mismatch and chr_map fixes it", {
  tt <- .toy_annotated_inputs()
  g  <- data.table::copy(tt$genes)
  g[, chrom := paste0(chrom, "H")]            # GFF3 says 1H/2H, container 1/2

  expect_warning(compute_block_annotation(tt$meta, g, flank_bp = 0),
                 "no block overlaps any gene")

  ann <- compute_block_annotation(tt$meta, g, flank_bp = 0,
                                  chr_map = c("1" = "1H", "2" = "2H"))
  expect_gt(sum(ann$blocks$gene_count), 0L)
  # stored chromosomes stay in the container vocabulary
  expect_true(all(ann$blocks$Chrom %in% c("1", "2")))
})

test_that("compute_block_annotation: QTL table and browser URLs", {
  tt  <- .toy_annotated_inputs()
  qtl <- data.frame(chrom = "1",
                    start = tt$b1$Start, end = tt$b1$End,
                    name  = "QTL_A", stringsAsFactors = FALSE)
  ann <- compute_block_annotation(tt$meta, tt$genes, qtl = qtl, flank_bp = 0,
                                  browser_base = "https://plants.ensembl.org/X")
  bl <- as.data.frame(ann$blocks)
  expect_true(all(c("qtl_count", "qtl_names") %in% names(bl)))
  expect_equal(bl$qtl_count[bl$Block == tt$b1$Block], 1L)
  expect_equal(bl$qtl_names[bl$Block == tt$b1$Block], "QTL_A")
  expect_equal(sum(bl$qtl_count), 1L)

  mg <- as.data.frame(ann$marker_gene)
  if (nrow(mg) > 0L)
    expect_true(all(grepl("^https://plants\\.ensembl\\.org/X/Location/View\\?r=",
                          mg$ensembl_browser_url)))
})

test_that("compute_block_annotation joins xrefs on the stripped id only", {
  tt <- .toy_annotated_inputs()
  xr <- data.table::data.table(
    gene_id = c("IN1", "STRAD"),              # bare ids, as a dump stores them
    gene_symbol = c("LOCIN1", "LOCSTRAD"),
    gene_description = c("d1", "d2"),
    go_terms = c("GO:1:x", ""),
    interpro_domains = c("", ""))
  ann <- compute_block_annotation(tt$meta, tt$genes, xrefs = xr, flank_bp = 0)
  mg  <- as.data.frame(ann$marker_gene)
  hit <- mg[mg$gene_id == "gene:IN1", ]
  if (nrow(hit) > 0L) {
    expect_equal(unique(hit$gene_symbol), "LOCIN1")
    expect_equal(unique(hit$gene_id), "gene:IN1")   # stored id keeps the prefix
  }
  expect_true(all(mg$gene_symbol[mg$gene_id == "gene:NEAR"] == ""))
})

test_that("compute_block_annotation refuses a non-bp container", {
  tt <- .toy_annotated_inputs()
  tt$meta$metadata$pos_unit <- "mbp"
  expect_error(compute_block_annotation(tt$meta, tt$genes),
               "bp")
})


# --- attach / validate / export ---------------------------------------------

test_that("attach_block_annotation stores tables + provenance and validates keys", {
  tt  <- .toy_annotated_inputs()
  ann <- compute_block_annotation(tt$meta, tt$genes, flank_bp = 0)
  m2  <- attach_block_annotation(tt$meta, ann)

  expect_true(is_geno_meta(m2))
  expect_identical(names(m2$annotation), c("blocks", "features", "marker_gene"))
  expect_identical(m2$metadata$annotation$flank_bp, 0)
  expect_null(tt$meta$annotation)              # input untouched

  # an annotation computed against other blocks is refused
  bad <- ann
  bad$blocks$Block[1L] <- "NOT_A_BLOCK"
  expect_error(attach_block_annotation(tt$meta, bad), "different block")

  # is_geno_meta flags a corrupted stored annotation
  m3 <- m2
  m3$annotation$marker_gene$Marker[1L] <- "NOT_A_MARKER"
  if (nrow(m3$annotation$marker_gene) > 0L)
    expect_false(is_geno_meta(m3, verbose = FALSE))
})

test_that("export_block_annotation writes the Python schema and messages a flank", {
  tt  <- .toy_annotated_inputs()
  ann <- compute_block_annotation(tt$meta, tt$genes, flank_bp = 0)
  m2  <- attach_block_annotation(tt$meta, ann)

  od <- file.path(tempdir(), "hs_ann_export")
  unlink(od, recursive = TRUE)
  paths <- export_block_annotation(m2, od)
  expect_true(all(file.exists(paths)))
  expect_setequal(basename(unname(paths)),
                  c("step1_block_coordinates.tsv", "step1_block_marker_map.tsv",
                    "step3_local_annotation.tsv",
                    "step3_local_annotation_records.tsv",
                    "step4_deep_annotation.tsv"))

  s3 <- data.table::fread(paths[["step3_local_annotation"]], sep = "\t",
                          colClasses = "character")
  expect_identical(names(s3),
                   c("block", "chrom", "start_pos", "end_pos", "n_markers",
                     "length_bp", "gene_count", "gene_ids", "gene_names"))
  expect_equal(nrow(s3), nrow(ann$blocks))

  s4 <- data.table::fread(paths[["step4_deep_annotation"]], sep = "\t",
                          colClasses = "character")
  expect_identical(names(s4),
                   c("marker_name", "marker_chr", "marker_pos", "block",
                     "gene_id", "gene_symbol", "gene_name", "gene_biotype",
                     "gene_start", "gene_end", "gene_strand",
                     "gene_description", "go_terms", "interpro_domains",
                     "ensembl_browser_url"))

  ann10 <- compute_block_annotation(tt$meta, tt$genes)   # default 10 kb
  expect_message(export_block_annotation(ann10, od), "flank_bp = 10,000")
})


# --- carry-forward (Q25.i) ---------------------------------------------------

test_that("filter_geno_meta prunes the stored annotation", {
  tt  <- .toy_annotated_inputs()
  ann <- compute_block_annotation(tt$meta, tt$genes, flank_bp = 0)
  m2  <- attach_block_annotation(tt$meta, ann)

  f <- filter_geno_meta(m2, chromosomes = list(list(chr = "1")))
  expect_false(is.null(f$annotation))
  expect_true(all(f$annotation$blocks$Chrom == "1"))
  expect_true(all(f$annotation$marker_gene$Marker %in% f$map_info$Marker))
  expect_true(is_geno_meta(f))
  # provenance travels in metadata
  expect_identical(f$metadata$annotation$flank_bp, 0)
})

test_that("attach_block_info drops the stored annotation with a message", {
  tt  <- .toy_annotated_inputs()
  ann <- compute_block_annotation(tt$meta, tt$genes, flank_bp = 0)
  m2  <- attach_block_annotation(tt$meta, ann)

  bmap <- sim_block_map(m2$map_info, block_size = 5L)   # new block definitions
  expect_message(m3 <- attach_block_info(m2, bmap), "annotation was dropped")
  expect_null(m3$annotation)
  expect_null(m3$metadata$annotation)
})

test_that("phenotype / signal / recode transforms keep the annotation", {
  tt  <- .toy_annotated_inputs()
  ann <- compute_block_annotation(tt$meta, tt$genes, flank_bp = 0)
  m2  <- attach_block_annotation(tt$meta, ann)

  ph <- data.frame(Sample = dimnames(m2$geno_arrays)[[1L]],
                   Yield = seq_len(dim(m2$geno_arrays)[1L]))
  m_ph <- attach_phenotype(m2, ph)
  expect_identical(m_ph$annotation, m2$annotation)

  m_rc <- recode_geno_meta(m2, to_type = "selectiontools")
  expect_identical(m_rc$annotation, m2$annotation)
})

test_that("a stored annotation survives an .rds round trip identically", {
  tt  <- .toy_annotated_inputs()
  ann <- compute_block_annotation(tt$meta, tt$genes, flank_bp = 0)
  m2  <- attach_block_annotation(tt$meta, ann)

  p <- tempfile(fileext = ".rds")
  saveRDS(m2, p)
  back <- readRDS(p)
  expect_true(is_geno_meta(back))
  expect_identical(back$annotation, m2$annotation)
  expect_identical(back$metadata$annotation, m2$metadata$annotation)
})
