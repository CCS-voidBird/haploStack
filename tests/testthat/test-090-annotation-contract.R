# =============================================================================
# test-090-annotation-contract.R — one locus contract (design D5 §3.6, steps
# D5-7 / D5-8 / D5-9)
#
#   N1  chr_map is container -> annotation on standardize / attach / track;
#       a matched row stores the container's spelling
#   N2  read_gff3_genes(rename_seqid =); the deprecated chr_map warns
#   N3  partial vocabulary: hs_warning_input names the chromosomes, suggests
#       the affix map, fills provenance$chrom_unmatched; strict errors
#   N4  one unit rule: pos_unit present and not bp -> error; absent + small
#       positions -> hs_warning_input
#   N5  .hs_nearest_marker() == the 0.8.8 per-row loop on random maps with
#       ties and duplicate positions
#   N6  the attach shortcut rejects a foreign Marker / Block
#   N7  provenance: md5 + base name, no absolute path / mtime; ISO created
#   N8  QTL BED: track / browser / # lines skipped; start > end refused
#   N9  a GFF3 coordinate of 3e9 reads; feature_type carried
#   N10 exporter manifests: absolute paths, byte counts (three exporters)
# =============================================================================

.n_meta <- function() {
  sim_geno_meta(n_samples = 4L, n_markers = 20L, encoding = "nucleotide",
                n_chrom = 2L, block_size = 4L, seed = 7L)
}

.n_genes <- function(meta, chrom_suffix = "") {
  bi <- as.data.frame(build_block_location_info(meta)$block_info)
  bi <- bi[order(bi$Chrom, bi$Start), ]
  b1 <- bi[1L, ]; b2 <- bi[2L, ]
  data.table::data.table(
    chrom     = paste0(c("1", "1", "1", "2"), chrom_suffix),
    start     = c(b1$Start + 1, b1$End - 1, b2$End + 5000, 1),
    end       = c(b1$Start + 2, b2$Start + 1, b2$End + 6000, 2),
    strand    = "+",
    gene_id   = c("gene:IN1", "gene:STRAD", "gene:NEAR", "gene:FAR2"),
    gene_name = c("InOne", NA, "NearTwo", "FarTwo"),
    biotype   = "protein_coding")
}

.n_gff3 <- function(path, big = FALSE) {
  s <- if (big) "3000000000" else "1000"
  e <- if (big) "3000000500" else "2000"
  writeLines(c("##gff-version 3",
               paste("1", "src", "gene", s, e, ".", "+", ".",
                     "ID=gene:G1;Name=GeneOne;biotype=protein_coding", sep = "\t"),
               paste("1", "src", "Gene", "5000", "9000", ".", "-", ".",
                     "ID=gene:G2;biotype=protein_coding", sep = "\t"),
               paste("2", "src", "mRNA", "100", "900", ".", "+", ".",
                     "ID=transcript:T1", sep = "\t")), path)
  path
}

test_that("N1: chr_map is container -> annotation; a match keeps the container spelling", {
  m <- .n_meta()
  map <- m$map_info
  # the user's table spells the chromosome 1H; the container says 1
  d <- data.frame(Chrom = "1H", Pos = as.numeric(map$Pos[3L]) + 10, Name = "Q")
  expect_warning(s0 <- standardize_marker_annotations(d, meta = m, match_bp = 100),
                 "match no marker")
  expect_true(is.na(s0$Marker))                       # no map: 1H matches nothing
  s1 <- standardize_marker_annotations(d, meta = m, match_bp = 100,
                                       chr_map = c("1" = "1H"))
  expect_identical(s1$Marker, map$Marker[3L])
  expect_identical(s1$Chrom, "1")                     # container spelling stored
  expect_equal(s1$matched_dist_bp, 10)
  # attach passes it through and records it
  a <- suppressMessages(attach_marker_annotations(m, d, match_bp = 100,
                                                  chr_map = c("1" = "1H")))
  expect_identical(a$annotation$custom$Marker, map$Marker[3L])
  expect_identical(a$metadata$annotation$custom$chr_map, c("1" = "1H"))
  expect_identical(a$metadata$annotation$custom$unplaced, character(0L))
  # the label track too
  tr <- track_marker_labels(d, meta = m, match_bp = 100, chr_map = c("1" = "1H"))
  expect_identical(tr$data$Chrom, "1")
})

test_that("N2: read_gff3_genes(rename_seqid =) renames file -> stored; chr_map is the deprecated spelling", {
  gff <- .n_gff3(tempfile(fileext = ".gff3"))
  g <- read_gff3_genes(gff, rename_seqid = c("1" = "chr1"))
  expect_identical(g$chrom, c("chr1", "chr1"))
  expect_warning(g2 <- read_gff3_genes(gff, chr_map = c("1" = "chr1")),
                 "read_gff3_genes(chr_map =) is deprecated", fixed = TRUE,
                 class = "hs_warning_deprecated")
  expect_identical(as.data.frame(g2), as.data.frame(g))
})

test_that("N3: a partial vocabulary warns before the overlap, names chromosomes, suggests the affix map", {
  m <- .n_meta()
  g <- .n_genes(m)
  g[chrom == "2", chrom := "chr2"]                     # chromosome 2 unmatched
  w <- tryCatch(compute_block_annotation(m, g, flank_bp = 0), warning = identity)
  expect_s3_class(w, "hs_warning_input")
  expect_match(conditionMessage(w), "1 container chromosome(s) have no gene", fixed = TRUE)
  expect_match(conditionMessage(w), 'chr_map = c("2" = "chr2")', fixed = TRUE)
  ann <- suppressWarnings(compute_block_annotation(m, g, flank_bp = 0))
  expect_identical(ann$provenance$chrom_unmatched, "2")
  expect_identical(ann$provenance$chrom_vocab$container, c("1", "2"))
  expect_identical(ann$provenance$chrom_vocab$annotation, c("1", "chr2"))
  # the chr_map clears it
  ann2 <- compute_block_annotation(m, g, flank_bp = 0, chr_map = c("2" = "chr2"))
  expect_identical(ann2$provenance$chrom_unmatched, character(0L))
  expect_true(all(ann2$blocks$Chrom %in% c("1", "2")))
  # strict
  e <- tryCatch(compute_block_annotation(m, g, flank_bp = 0, strict = TRUE),
                error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_identical(e$bad, "2")
  # a QTL table with a foreign seqid is checked the same way
  qtl <- data.frame(chrom = "9H", start = 1, end = 2, name = "q")
  expect_warning(compute_block_annotation(m, .n_genes(m), flank_bp = 0, qtl = qtl),
                 "no QTL", class = "hs_warning_input")
})

test_that("N4: one unit rule — declared non-bp errors, undeclared small positions warn", {
  m <- .n_meta()
  m_mb <- m; m_mb$metadata$pos_unit <- "mbp"
  e <- tryCatch(haploStack:::.hs_require_bp(m_mb, "f"), error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "^f: map_info\\$Pos is stored in 'mbp'")
  expect_error(compute_block_annotation(m_mb, .n_genes(m)), class = "hs_error_input")
  expect_error(standardize_marker_annotations(
    data.frame(Chrom = "1", Pos = 1, Name = "x"), meta = m_mb), class = "hs_error_input")
  # declared bp: silent
  m_bp <- m; m_bp$metadata$pos_unit <- "bp"
  expect_silent(haploStack:::.hs_require_bp(m_bp, "f"))
  # undeclared, positions that look like Mb
  m_small <- m; m_small$map_info$Pos <- as.numeric(m$map_info$Pos) / 1e5
  expect_warning(haploStack:::.hs_require_bp(m_small, "f"),
                 "look like Mb or cM", class = "hs_warning_input")
  expect_silent(haploStack:::.hs_require_bp(m, "f"))   # bp-scale positions pass
})

test_that("N5: the vectorised nearest-marker matcher equals the 0.8.8 loop (ties, duplicates)", {
  oracle <- function(map, chrom, pos, max_bp) {
    idx <- rep(NA_integer_, length(chrom)); d_out <- rep(NA_real_, length(chrom))
    for (i in seq_along(chrom)) {
      if (is.na(chrom[i]) || is.na(pos[i])) next     # the 0.8.8 loop ran on located rows only
      cand <- which(map$Chrom == chrom[i])
      if (!length(cand)) next
      d <- abs(map$Pos[cand] - pos[i])
      j <- which.min(d)
      if (d[j] <= max_bp) { idx[i] <- cand[j]; d_out[i] <- d[j] }
    }
    list(index = idx, dist = d_out)
  }
  set.seed(90)
  for (rep in 1:20) {
    n_map <- sample(5:40, 1L)
    map <- data.frame(Chrom = sample(c("1", "2", "3"), n_map, replace = TRUE),
                      Pos = sample(0:30, n_map, replace = TRUE),    # many duplicates
                      stringsAsFactors = FALSE)
    nq <- sample(1:30, 1L)
    chrom <- sample(c("1", "2", "3", "4", NA), nq, replace = TRUE)
    pos   <- sample(c(-3:33, NA), nq, replace = TRUE)
    max_bp <- sample(c(0, 1, 2, 5, 100), 1L)
    o <- oracle(map, chrom, pos, max_bp)
    v <- haploStack:::.hs_nearest_marker(map$Chrom, map$Pos, chrom, pos, max_bp)
    expect_identical(v$index, o$index)
    expect_equal(v$dist, o$dist)
  }
  # the standardizer uses it: identical result to the loop on a real map
  m <- .n_meta()
  map <- m$map_info
  d <- data.frame(Chrom = c("1", "2", "1"), Name = c("a", "b", "c"),
                  Pos = c(as.numeric(map$Pos[2L]) + 5, as.numeric(map$Pos[15L]) - 7, 1))
  expect_warning(s <- standardize_marker_annotations(d, meta = m, match_bp = 10),
                 "1 of 3 row(s) match no marker", fixed = TRUE)
  expect_identical(s$Marker, c(map$Marker[2L], map$Marker[15L], NA))
  expect_equal(s$matched_dist_bp, c(5, 7, NA))
})

test_that("N6: the attach shortcut checks a standardized table's Marker and Block", {
  m <- .n_meta()
  std <- standardize_marker_annotations(
    data.frame(Marker = m$map_info$Marker[1:2], Name = c("a", "b")), meta = m)
  ok <- suppressMessages(attach_marker_annotations(m, std))
  expect_identical(ok$annotation$custom$Marker, m$map_info$Marker[1:2])
  bad <- std; bad$Marker[2L] <- "NOT_A_MARKER"
  e <- tryCatch(attach_marker_annotations(m, bad), error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "not in the container's map (first: b)", fixed = TRUE)
  bad2 <- std; bad2$Block[1L] <- "NOT_A_BLOCK"
  expect_error(attach_marker_annotations(m, bad2), "not in the container's blocks",
               class = "hs_error_input")
  # the custom provenance has the new stamps
  expect_match(ok$metadata$annotation$custom$created,
               "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$")
  expect_null(ok$metadata$annotation$custom$date)
})

test_that("N7: provenance carries md5 + base name, no absolute path, an ISO-8601 created stamp", {
  m <- .n_meta()
  gff <- .n_gff3(tempfile(fileext = ".gff3"))
  g <- read_gff3_genes(gff)
  src <- attr(g, "hs_annotation_source")
  expect_identical(names(src), c("file", "md5", "size", "n"))
  expect_identical(src$file, basename(gff))
  expect_identical(src$md5, unname(tools::md5sum(gff)))
  expect_identical(src$n, 2L)
  ann <- suppressWarnings(compute_block_annotation(m, g, flank_bp = 0))
  expect_identical(ann$provenance$gff3, src)
  expect_null(ann$provenance$gff3$path)
  expect_null(ann$provenance$date)
  expect_match(ann$provenance$created, "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$")
  # an in-memory table has a minimal record
  ann2 <- compute_block_annotation(m, .n_genes(m), flank_bp = 0)
  expect_true(is.na(ann2$provenance$gff3$file))
  expect_identical(ann2$provenance$gff3$n, 4L)
  # a QTL path gets the same record
  qp <- tempfile(fileext = ".bed")
  bi <- as.data.frame(build_block_location_info(m)$block_info)
  writeLines(paste("1", bi$Start[1L], bi$End[1L], "qtl1", sep = "\t"), qp)
  ann3 <- compute_block_annotation(m, .n_genes(m), flank_bp = 0, qtl = qp)
  expect_identical(ann3$provenance$qtl$file, basename(qp))
  expect_identical(ann3$provenance$qtl$n, 1L)
})

test_that("N8: the QTL BED reader skips track / browser / comment lines and refuses start > end", {
  p <- tempfile(fileext = ".bed")
  writeLines(c("browser position chr1:1-100", "track name=qtl", "# a comment",
               paste("1", 10, 20, "q1", sep = "\t"), "",
               paste("2", 5, 5, "q2", sep = "\t")), p)
  q <- haploStack:::.read_qtl_bed(p)
  expect_identical(q$name, c("q1", "q2"))
  expect_equal(q$start, c(10, 5))
  expect_identical(attr(q, "hs_annotation_source")$n, 2L)
  writeLines(paste("1", 30, 20, "rev", sep = "\t"), p)
  expect_error(haploStack:::.read_qtl_bed(p), "start > end", class = "hs_error_format")
  writeLines(paste("1", "x", 20, "nan", sep = "\t"), p)
  expect_error(haploStack:::.read_qtl_bed(p), "non-numeric", class = "hs_error_format")
  writeLines(paste("1", 1, sep = "\t"), p)
  expect_error(haploStack:::.read_qtl_bed(p), "4 tab-separated", class = "hs_error_format")
})

test_that("N9: a coordinate beyond int32 reads; feature_type is carried", {
  gff <- .n_gff3(tempfile(fileext = ".gff3"), big = TRUE)
  g <- read_gff3_genes(gff)
  expect_type(g$start, "double")
  expect_equal(g$start[1L], 3e9)
  expect_identical(g$feature_type, c("gene", "Gene"))   # the file's spelling
  m <- .n_meta()
  ann <- suppressWarnings(compute_block_annotation(m, .n_genes(m)[, feature_type := "locus"],
                                                   flank_bp = 0))
  expect_true(all(ann$features$feature_type == "locus"))
  ann0 <- compute_block_annotation(m, .n_genes(m), flank_bp = 0)
  expect_true(all(ann0$features$feature_type == "gene"))
})

test_that("N11: the annotation counts follow the tables through a prune (ann/F7)", {
  m <- .n_meta()
  ann <- compute_block_annotation(m, .n_genes(m), flank_bp = 0)
  m2 <- attach_block_annotation(m, ann)
  d <- data.frame(Marker = c(m$map_info$Marker[1:2], "nowhere"), Name = c("a", "b", "c"))
  m2 <- suppressWarnings(suppressMessages(attach_marker_annotations(m2, d)))
  expect_identical(m2$metadata$annotation$custom$unplaced, "c")
  f <- filter_geno_meta(m2, chromosomes = list(list(chr = "2")))
  p <- f$metadata$annotation
  expect_identical(p$n_blocks, nrow(f$annotation$blocks))
  expect_identical(p$n_blocks_hit, sum(f$annotation$blocks$gene_count > 0L))
  expect_identical(p$n_features, nrow(f$annotation$features))
  expect_identical(p$n_marker_gene, nrow(f$annotation$marker_gene))
  expect_true(p$n_blocks < ann$provenance$n_blocks)
  expect_identical(p$custom$n, nrow(f$annotation$custom))
  expect_identical(p$custom$n_matched, sum(!is.na(f$annotation$custom$Marker)))
  expect_identical(p$custom$unplaced, "c")      # a position-only row stays
})

test_that("N10: the three exporters return manifests with absolute paths and byte counts", {
  m <- .n_meta()
  d <- withr::local_tempdir()
  is_abs <- function(p) all(grepl("^([A-Za-z]:)?/", p))
  # export_block_annotation(): absolute paths back, manifest + notes attributes
  ann <- compute_block_annotation(m, .n_genes(m), flank_bp = 0)
  m2 <- attach_block_annotation(m, ann)
  paths <- export_block_annotation(m2, file.path(d, "ann"))
  expect_true(is_abs(paths))
  mf <- attr(paths, "hs_files")
  expect_identical(names(mf), c("path", "kind", "object", "bytes"))
  expect_identical(nrow(mf), 5L)
  expect_identical(unname(mf$bytes), unname(as.numeric(file.size(mf$path))))
  expect_identical(attr(paths, "hs_notes"), character(0L))
  ann10 <- compute_block_annotation(m, .n_genes(m))
  p10 <- suppressMessages(export_block_annotation(ann10, file.path(d, "ann10")))
  expect_match(attr(p10, "hs_notes"), "flank_bp = 10,000", fixed = TRUE)
  # stored blocks are in block_info order; the files are in the Python order
  expect_identical(m2$annotation$blocks$Block, as.character(m2$block_info$Block))
  s1 <- data.table::fread(paths[["step1_block_coordinates"]], colClasses = "character")
  expect_identical(s1$block, ann$blocks$Block)
  # export_marker_effects(): `files` as before, manifest in hs_files
  eff <- sim_effect_df(m$map_info$Marker, alleles = c("A", "T"), seed = 2L)
  me <- m; me$marker_effects <- list(Yield = eff)
  out <- suppressMessages(export_marker_effects(me, output_dir = file.path(d, "eff")))
  expect_identical(basename(attr(out, "files")), "marker_effects_Yield.txt")
  expect_true(is_abs(attr(out, "hs_files")$path))
  expect_identical(attr(out, "hs_files")$bytes, as.numeric(file.size(attr(out, "files"))))
  expect_identical(attr(out, "hs_notes"), character(0L))
  me2 <- m; me2$marker_effects <- list("Yield (t/ha)" = eff)
  out2 <- suppressMessages(export_marker_effects(me2, output_dir = file.path(d, "eff2")))
  expect_match(attr(out2, "hs_notes"), "written as 'Yield_t_ha'", fixed = TRUE)
  # build_block_location_info(output_dir =)
  r <- suppressMessages(build_block_location_info(m, output_dir = file.path(d, "blk"), prefix = "p"))
  expect_identical(nrow(attr(r, "hs_files")), 2L)
  expect_true(is_abs(attr(r, "hs_files")$path))
  expect_null(attr(build_block_location_info(m), "hs_files"))
})
