# Extracted from test-090-blocks.R:70

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
bmap <- sim_block_map(m$map_info, block_size = 8L)
bmap$Block <- ifelse(bmap$Chrom == "1", "007", "1e5")
d <- tempfile("m1_")
dir.create(d)
f <- .c090_write(stats::setNames(bmap, c("chrom", "marker", "block")),
                   file.path(d, "my_blocks.tsv"))
a <- attach_block_info(m, f)
expect_true(is_geno_meta(a, strict = TRUE))
expect_setequal(unique(a$block_marker_map_info$Block), c("007", "1e5"))
expect_identical(a$block_marker_map_info$Marker, m$map_info$Marker)
p <- a$metadata$blocks
expect_identical(p$source, "user_map")
expect_identical(p$label, "my_blocks")
expect_identical(p$file, "my_blocks.tsv")
expect_identical(p$file_md5, unname(tools::md5sum(f)))
expect_identical(p[c("n_blocks", "n_singletons", "n_matched", "n_unmatched")],
                   list(n_blocks = 2L, n_singletons = 0L, n_matched = 16L, n_unmatched = 0L))
expect_identical(p$params$layout, "block_marker_map")
b <- attach_block_info(m, bmap)
expect_identical(b$block_marker_map_info, a$block_marker_map_info)
expect_identical(b$metadata$blocks$map_md5, p$map_md5)
expect_true(is.na(b$metadata$blocks$file))
w <- suppressMessages(build_block_location_info(a, output_dir = d, prefix = "w"))
