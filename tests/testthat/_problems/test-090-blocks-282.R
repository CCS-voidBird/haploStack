# Extracted from test-090-blocks.R:282

# prequel ----------------------------------------------------------------------
.c090_meta <- function() sim_geno_meta(n_samples = 8L, n_markers = 16L,
                                       block_size = 4L, seed = 1L)
.c090_regions <- function(m) {
  bi <- m$block_info
  data.frame(Chrom = bi$Chrom, Block = bi$Block, Start = bi$Start, End = bi$End,
             stringsAsFactors = FALSE)
}
.c090_write <- function(df, path, ...) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, ...)
  path
}

# test -------------------------------------------------------------------------
m <- .c090_meta()
d <- tempfile("m9_")
dir.create(d)
expect_message(paths <- build_block_location_info(m, output_dir = d, prefix = "x"),
                 "Block files written", class = "hs_message")
