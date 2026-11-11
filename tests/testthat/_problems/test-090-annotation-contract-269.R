# Extracted from test-090-annotation-contract.R:269

# prequel ----------------------------------------------------------------------
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

# test -------------------------------------------------------------------------
m <- .n_meta()
d <- withr::local_tempdir()
is_abs <- function(p) all(grepl("^([A-Za-z]:)?/", p))
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
expect_identical(m2$annotation$blocks$Block, as.character(m2$block_info$Block))
s1 <- data.table::fread(paths[["step1_block_coordinates"]], colClasses = "character")
expect_identical(s1$block, ann$blocks$Block)
eff <- sim_effect_df(m$map_info$Marker, alleles = c("A", "T"), seed = 2L)
me <- attach_marker_effects(m, eff, trait = "Yield")
