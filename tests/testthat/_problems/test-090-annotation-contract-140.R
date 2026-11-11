# Extracted from test-090-annotation-contract.R:140

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
oracle <- function(map, chrom, pos, max_bp) {
    idx <- rep(NA_integer_, length(chrom)); d_out <- rep(NA_real_, length(chrom))
    for (i in seq_along(chrom)) {
      cand <- which(map$Chrom == chrom[i])
      if (!length(cand)) next
      d <- abs(map$Pos[cand] - pos[i])
      j <- which.min(d)
      if (d[j] <= max_bp) { idx[i] <- cand[j]; d_out[i] <- d[j] }
    }
    list(index = idx, dist = d_out)
  }
