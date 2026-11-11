# =============================================================================
# test-090-blocks.R — the block map (phase C: D2-5 one resolver, D2-6
# provenance, D2-7 one writer; owner steer 2026-10-29: two shapes defined by the
# functions, SelectionTools formats behind one switch)
#
#   M1  a block-marker map file with the canonical header attaches; identifiers
#       stay text; the file and the data.frame give the same map_md5
#   M2  anything but the two shapes is refused (hs_error_block_map naming them):
#       a headerless file, an alias header, a 2-column table; an unnamed
#       3-/4-column data.frame is read positionally
#   M3  a block-region map attaches: markers by position, Region_Start/End
#       beside the marker-derived span, an empty region omitted with a warning;
#       expand_block_info() returns what the 0.8.8 loop returned
#   M4  chromosomes: a consistent rename is matched on Marker and recorded; an
#       inconsistent one, a block on two chromosomes and an unknown region
#       chromosome are errors (with a chrom_map suggestion); chrom_map applies
#   M5  match rate: no marker in the panel is an error; under 50 % a warning
#   M6  overlap: exact duplicates dropped silently; a marker in two blocks is an
#       error with `bad`, or kept-first with a warning under overlap = "first"
#   M7  readers build then attach: block_file takes a data.frame, forwards
#       block_overlap, keeps the FAM; .read_block_file() is gone
#   M8  provenance: user map, LD resolver (kept through attach), pruned_from on
#       filter, replaced on re-blocking
#   M9  one writer: 7 columns, plain integers, round trip through attach;
#       the resolver's attach = FALSE files are the same
#   M10 the SelectionTools switch: off by default (hs_error_input with hint),
#       use_st_formats() returns the previous value, st_formats = TRUE per call,
#       a malformed option is refused, the deprecated names are gated too
#   M11 st_convert(): GenoMeta -> SelectionTools files -> GenoMeta round trip
#       (no SelectionTools needed); SelectionTools' block table -> marker map
# =============================================================================

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

test_that("M1: a block-marker map file attaches; identifiers stay text; one md5", {
  m <- .c090_meta()
  bmap <- sim_block_map(m$map_info, block_size = 8L)
  bmap$Block <- ifelse(bmap$Chrom == "1", "007", "1e5")      # ids that type inference rewrites
  d <- tempfile("m1_"); dir.create(d)
  f <- .c090_write(stats::setNames(bmap, c("chrom", "marker", "block")),
                   file.path(d, "my_blocks.tsv"))             # lower-case header is fine
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
  # the same map as a data.frame: same block set, same identity
  b <- attach_block_info(m, bmap)
  expect_identical(b$block_marker_map_info, a$block_marker_map_info)
  expect_identical(b$metadata$blocks$map_md5, p$map_md5)
  expect_true(is.na(b$metadata$blocks$file))
  # and through the CSV the writer produces
  w <- suppressMessages(build_block_location_info(a, output_dir = d, prefix = "w"))
  c2 <- suppressMessages(attach_block_info(m, file.path(d, "w_block_marker_map.csv")))
  expect_identical(c2$metadata$blocks$map_md5, p$map_md5)
  expect_identical(c2$block_marker_map_info, a$block_marker_map_info)
})

test_that("M2: only the two shapes are accepted", {
  m <- .c090_meta()
  bmap <- sim_block_map(m$map_info, block_size = 8L)
  d <- tempfile("m2_"); dir.create(d)
  # a headerless file: its first line is read as the header and does not fit
  f <- .c090_write(bmap, file.path(d, "noheader.tsv"), col.names = FALSE)
  e <- tryCatch(attach_block_info(m, f), error = identity)
  expect_s3_class(e, "hs_error_block_map")
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "block-marker map with columns Chrom, Marker, Block")
  expect_match(conditionMessage(e), "block-region map with columns Chrom, Block, Start, End")
  # an alias header (what 0.8.8 silently accepted)
  e <- tryCatch(attach_block_info(m, data.frame(Chr = bmap$Chrom, SNP = bmap$Marker,
                                                 Block = bmap$Block)), error = identity)
  expect_s3_class(e, "hs_error_block_map")
  expect_identical(e$columns, c("Chr", "SNP", "Block"))
  # two columns
  expect_error(attach_block_info(m, bmap[, c("Marker", "Block")]), class = "hs_error_block_map")
  # forcing a shape the table does not have
  expect_error(attach_block_info(m, bmap, source_type = "block_region_map"),
               "block-region map", class = "hs_error_block_map")
  # unnamed data.frames are read positionally
  u3 <- stats::setNames(bmap, paste0("V", 1:3))
  expect_identical(attach_block_info(m, u3)$block_marker_map_info,
                   attach_block_info(m, bmap)$block_marker_map_info)
  u4 <- stats::setNames(.c090_regions(m), paste0("V", 1:4))
  expect_identical(attach_block_info(m, u4)$block_marker_map_info$Block,
                   m$block_marker_map_info$Block)
  # not a source at all
  expect_error(attach_block_info(m, 42), class = "hs_error_block_map")
  expect_error(attach_block_info(m, file.path(d, "nope.tsv")), "file not found")
  # a bad Start
  reg <- .c090_regions(m); reg$Start[2L] <- "x"
  expect_error(attach_block_info(m, reg), "Start is not numeric", class = "hs_error_block_map")
  reg <- .c090_regions(m); reg$Start[2L] <- reg$End[2L] + 1
  expect_error(attach_block_info(m, reg), "Start > End", class = "hs_error_block_map")
})

test_that("M3: a block-region map attaches by position; expand_block_info() is the old loop", {
  m <- .c090_meta()
  reg <- .c090_regions(m)
  b <- attach_block_info(m, reg)
  expect_identical(b$block_marker_map_info, m$block_marker_map_info)
  expect_true(all(c("Region_Start", "Region_End") %in% names(b$block_info)))
  expect_identical(b$block_info$Region_Start, b$block_info$Start)
  expect_identical(b$block_info$Region_End, b$block_info$End)
  expect_identical(b$metadata$blocks$params$layout, "block_region_map")
  # source_type = "block_info" still means the region map
  expect_identical(attach_block_info(m, reg, source_type = "block_info")$block_info,
                   b$block_info)
  # a region with no marker is omitted with a warning, and counted
  reg2 <- rbind(reg, data.frame(Chrom = "1", Block = "empty", Start = 1, End = 2))
  expect_warning(b2 <- attach_block_info(m, reg2), "contain no marker")
  expect_false("empty" %in% b2$block_info$Block)
  expect_identical(b2$metadata$blocks$n_unmatched, 1L)
  # a region on a chromosome the panel lacks
  reg3 <- rbind(reg, data.frame(Chrom = "9", Block = "x", Start = 1, End = 2))
  expect_error(attach_block_info(m, reg3), "not in this panel", class = "hs_error_block_map")
  # a wider declared span is kept beside the marker-derived one
  reg4 <- reg; reg4$Start <- reg4$Start - 100; reg4$End <- reg4$End + 100
  b4 <- attach_block_info(m, reg4)
  expect_identical(b4$block_info$Start, b$block_info$Start)
  expect_identical(b4$block_info$Region_Start, b$block_info$Start - 100)
  # the public wrapper returns exactly what the 0.8.8 loop returned
  old_loop <- function(map, bi) {
    rows <- lapply(seq_len(nrow(bi)), function(i) {
      hits <- map[map$Chrom == bi$Chrom[i] & map$Pos >= bi$Start[i] & map$Pos <= bi$End[i], ]
      if (nrow(hits) == 0L) return(NULL)
      data.frame(Chrom = bi$Chrom[i], Marker = hits$Marker, Block = bi$Block[i],
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
    rownames(out) <- NULL
    out
  }
  shuffled <- reg[c(3, 1, 4, 2), ]
  expect_identical(expand_block_info(m$map_info, shuffled), old_loop(m$map_info, shuffled))
  expect_warning(expand_block_info(m$map_info, reg2), "contains no markers")
  ov <- rbind(reg, data.frame(Chrom = "1", Block = "wide", Start = 0, End = 1e9))
  expect_warning(expand_block_info(m$map_info, ov), "more than one block")
  expect_warning(
    expect_warning(expand_block_info(m$map_info, data.frame(Chrom = "9", Block = "x",
                                                            Start = 1, End = 2)),
                   "contains no markers"),
    "no markers matched")
})

test_that("M4: chromosome names are checked, renamed consistently or refused", {
  m <- .c090_meta()
  bmap <- sim_block_map(m$map_info, block_size = 8L)
  # a consistent rename: 1 -> 1H, 2 -> 2H in the map, panel has 1 / 2
  ren <- bmap; ren$Chrom <- paste0(ren$Chrom, "H")
  expect_message(a <- attach_block_info(m, ren), "1H -> 1, 2H -> 2", class = "hs_message")
  expect_identical(unique(a$block_marker_map_info$Chrom), c("1", "2"))   # the panel's spelling
  expect_identical(a$metadata$blocks$params$chrom_renamed, c("1H" = "1", "2H" = "2"))
  # chrom_map does the same without the message
  expect_silent(b <- attach_block_info(m, ren, chrom_map = c("1H" = "1", "2H" = "2")))
  expect_identical(b$block_marker_map_info, a$block_marker_map_info)
  expect_identical(b$metadata$blocks$params$chrom_map, c("1H" = "1", "2H" = "2"))
  expect_error(attach_block_info(m, ren, chrom_map = c("1", "2")), "named character")
  # inconsistent: one map chromosome lands on two panel chromosomes
  bad <- bmap; bad$Chrom <- "1H"
  e <- tryCatch(attach_block_info(m, bad), error = identity)
  expect_s3_class(e, "hs_error_block_map")
  expect_match(conditionMessage(e), "one-to-one")
  # a block whose markers sit on two panel chromosomes
  span <- bmap; span$Block <- "one_block"
  e <- tryCatch(attach_block_info(m, span), error = identity)
  expect_match(conditionMessage(e), "more than one chromosome")
  expect_identical(e$bad, "one_block")
  # a region map on another spelling: refused with the chrom_map to pass
  reg <- .c090_regions(m); reg$Chrom <- paste0(reg$Chrom, "H")
  e <- tryCatch(attach_block_info(m, reg), error = identity)
  expect_s3_class(e, "hs_error_block_map")
  expect_match(conditionMessage(e), 'chrom_map = c("1H" = "1", "2H" = "2")', fixed = TRUE)
  expect_identical(e$chrom_map, c("1H" = "1", "2H" = "2"))
  r <- attach_block_info(m, reg, chrom_map = e$chrom_map)
  expect_identical(r$block_marker_map_info, m$block_marker_map_info)
})

test_that("M5: the match rate is checked", {
  m <- .c090_meta()
  bmap <- sim_block_map(m$map_info, block_size = 8L)
  none <- bmap; none$Marker <- paste0("x", none$Marker)
  e <- tryCatch(attach_block_info(m, none), error = identity)
  expect_s3_class(e, "hs_error_block_map")
  expect_match(conditionMessage(e), "no marker of the block map is in this panel")
  few <- rbind(bmap[1:4, ], data.frame(Chrom = "1", Marker = paste0("y", 1:12), Block = "B9"))
  expect_warning(a <- attach_block_info(m, few), "only 4 of 16")
  expect_identical(a$metadata$blocks$n_unmatched, 12L)
  expect_identical(a$metadata$blocks$n_singletons, 12L)   # the other panel markers
  expect_identical(sum(grepl("^dummy_", a$block_marker_map_info$Block)), 12L)
})

test_that("M6: overlap policy", {
  m <- .c090_meta()
  bmap <- sim_block_map(m$map_info, block_size = 8L)
  dup <- rbind(bmap, bmap[3:5, ])                         # exact duplicates: dropped quietly
  expect_no_warning(a <- attach_block_info(m, dup))
  expect_identical(a$block_marker_map_info, attach_block_info(m, bmap)$block_marker_map_info)
  ov <- rbind(bmap, data.frame(Chrom = "1", Marker = bmap$Marker[1:2], Block = "other"))
  e <- tryCatch(attach_block_info(m, ov), error = identity)
  expect_s3_class(e, "hs_error_block_map")
  expect_identical(e$bad, bmap$Marker[1:2])
  expect_match(conditionMessage(e), 'overlap = "first"', fixed = TRUE)
  expect_warning(f <- attach_block_info(m, ov, overlap = "first"), "first block of each is kept")
  expect_identical(f$block_marker_map_info, attach_block_info(m, bmap)$block_marker_map_info)
  expect_identical(f$metadata$blocks$params$n_overlap, 2L)
  expect_identical(f$metadata$blocks$params$overlap, "first")
})

test_that("M7: readers build, then attach", {
  m <- .c090_meta()
  bmap <- sim_block_map(m$map_info, block_size = 8L)
  d <- tempfile("m7_"); dir.create(d)
  suppressMessages(write_geno_meta(recode_geno_meta(m, "nucleotide"), file.path(d, "x"),
                                   preset = "PL"))
  r <- suppressMessages(suppressWarnings(
    read_plink_ped(file.path(d, "x.ped"), file.path(d, "x.map"), block_file = bmap)))
  expect_identical(r$block_marker_map_info$Block, attach_block_info(m, bmap)$block_marker_map_info$Block)
  expect_identical(r$metadata$blocks$source, "user_map")
  expect_identical(r$sample_info$Sample, m$sample_info$Sample)
  ov <- rbind(bmap, data.frame(Chrom = "1", Marker = bmap$Marker[1], Block = "other"))
  expect_error(suppressMessages(suppressWarnings(
    read_plink_ped(file.path(d, "x.ped"), file.path(d, "x.map"), block_file = ov))),
    class = "hs_error_block_map")
  expect_warning(suppressMessages(
    read_plink_ped(file.path(d, "x.ped"), file.path(d, "x.map"), block_file = ov,
                   block_overlap = "first")), "first block")
  g <- suppressMessages(suppressWarnings(
    read_plink_ped(file.path(d, "x.ped"), file.path(d, "x.map"),
                   block_file = .c090_regions(m))))
  expect_identical(g$metadata$blocks$params$layout, "block_region_map")
  expect_false(exists(".read_block_file", envir = asNamespace("haploStack")))
})

test_that("M8: provenance follows the blocks through the transforms", {
  skip_if_not_installed("HapSelect")
  toy <- st_toy()
  r <- hapS_ld_blocks(toy, ld_threshold = 0.5, tolerance = 3L, attach = FALSE, verbose = FALSE)
  expect_named(r, c("block_marker_map", "block_info", "provenance"))
  p <- r$provenance
  expect_identical(p$source, "hapS_ld_blocks")
  expect_identical(p$label, "hapS_ld_blocks_r0.5_t3")
  expect_identical(p$params$method, "flanking")
  expect_identical(p$n_blocks, length(unique(r$block_marker_map$Block)))
  # the list attaches and keeps the resolver's identity
  a <- suppressMessages(attach_block_info(toy, r))
  expect_identical(a$metadata$blocks$source, "hapS_ld_blocks")
  expect_identical(a$metadata$blocks$map_md5, p$map_md5)
  b <- suppressMessages(hapS_ld_blocks(toy, 0.5, 3L, verbose = FALSE))
  expect_identical(b$metadata$blocks, a$metadata$blocks)
  # a filter records the loss; a re-blocking replaces the record
  keep <- dimnames(a$geno_arrays)[[2L]][1:10]
  f <- suppressMessages(filter_geno_meta(a, keep_markers = keep))
  expect_identical(f$metadata$blocks$pruned_from, p$n_blocks)
  expect_lt(f$metadata$blocks$n_blocks, p$n_blocks)
  u <- suppressMessages(attach_block_info(a, sim_block_map(a$map_info, block_size = 5L)))
  expect_identical(u$metadata$blocks$source, "user_map")
  expect_null(u$metadata$blocks$pruned_from)
})

test_that("M9: one writer, seven columns, plain integers, exact round trip", {
  m <- .c090_meta()
  d <- tempfile("m9_"); dir.create(d)
  expect_message(paths <- build_block_location_info(m, output_dir = d, prefix = "x"),
                 "Block files written", class = "hs_message")
  bi <- readLines(file.path(d, "x_block_info.csv"))
  expect_identical(bi[1L], "Chrom,Block,Start,End,Pos,Num_SNPs,Block_length")
  expect_false(any(grepl("e+", bi, fixed = TRUE)))          # 84900000, not 8.49e+07
  expect_true(any(grepl(",84900000,", bi, fixed = TRUE)))
  bmm <- readLines(file.path(d, "x_block_marker_map.csv"))
  expect_identical(bmm[1L], "Chrom,Marker,Block")
  expect_length(bmm, 17L)
  back <- suppressMessages(attach_block_info(m, file.path(d, "x_block_marker_map.csv")))
  expect_identical(back$block_marker_map_info, m$block_marker_map_info)
  expect_identical(back$block_info, m$block_info)
  # the region file re-expands by position to the same blocks
  back2 <- suppressMessages(attach_block_info(m, file.path(d, "x_block_info.csv")))
  expect_identical(back2$block_marker_map_info, m$block_marker_map_info)
  # a region map's extra columns are written; a selection too
  r <- attach_block_info(m, .c090_regions(m))
  invisible(suppressMessages(build_block_location_info(r, output_dir = d, prefix = "r")))
  expect_identical(readLines(file.path(d, "r_block_info.csv"), n = 1L),
                   "Chrom,Block,Start,End,Pos,Num_SNPs,Block_length,Region_Start,Region_End")
  # the resolver's attach = FALSE route uses the same writer
  skip_if_not_installed("HapSelect")
  toy <- st_toy()
  invisible(hapS_ld_blocks(toy, 0.5, 3L, attach = FALSE, output_dir = d, prefix = "h",
                           verbose = FALSE))
  expect_identical(readLines(file.path(d, "h_block_info.csv"), n = 1L),
                   "Chrom,Block,Start,End,Pos,Num_SNPs,Block_length")
})

test_that("M10: SelectionTools formats are behind one switch", {
  withr::local_options(list(haploStack.st_formats = NULL))
  m <- .c090_meta()
  e <- tryCatch(st_ld_blocks(m), error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "SelectionTools formats are off")
  expect_identical(e$hint, "use_st_formats(TRUE)")
  expect_error(st_convert(m, to = "selectiontools", output_dir = tempfile()),
               "use_st_formats", class = "hs_error_input")
  expect_error(suppressWarnings(st_ld_blocks_files("a", "b"), classes = "hs_warning_deprecated"),
               "SelectionTools formats are off")
  expect_error(define_haplotype_blocks_meta(m), class = "hs_error_defunct")   # defunct first (0.9.0)
  # the setter returns the previous value and validates
  expect_message(old <- use_st_formats(TRUE), "enabled", class = "hs_message")
  expect_null(old)
  expect_true(getOption("haploStack.st_formats"))
  expect_message(prev <- use_st_formats(FALSE), "disabled")
  expect_true(prev)
  expect_error(use_st_formats("yes"), "TRUE or FALSE")
  # per call, and a malformed option is refused rather than obeyed
  d <- tempfile("m10_"); dir.create(d)
  expect_message(st_convert(m, to = "selectiontools", output_dir = d, st_formats = TRUE), "wrote")
  options(haploStack.st_formats = "yes")
  expect_error(st_convert(m, to = "selectiontools", output_dir = d), "TRUE or FALSE")
  # the registry knows it
  expect_true("st_formats" %in% names(haploStack:::.HS_SETTINGS))
})

test_that("M11: st_convert() round trip and the SelectionTools block table", {
  withr::local_options(list(haploStack.st_formats = TRUE))
  m <- .c090_meta()
  d <- tempfile("m11_"); dir.create(d)
  out <- suppressMessages(st_convert(m, to = "selectiontools", output_dir = d, prefix = "p"))
  expect_named(out, c("geno_file", "map_file", "chrom_keys"))
  expect_identical(basename(out$geno_file), "p_geno.txt")
  g <- readLines(out$geno_file)
  expect_identical(strsplit(g[1L], "\t")[[1L]], c("m/i", dimnames(m$geno_arrays)[[1L]]))
  expect_length(g, 17L)
  expect_identical(readLines(out$map_file, n = 1L), "name\tchrom\tpos")
  expect_identical(out$chrom_keys, data.frame(Chrom = c("1", "2"), key = 1:2,
                                              stringsAsFactors = FALSE))
  back <- suppressMessages(st_convert(out$geno_file, map = out$map_file))
  expect_true(is_geno_meta(back, strict = TRUE))
  expect_identical(dimnames(back$geno_arrays)[1:2], dimnames(m$geno_arrays)[1:2])
  # alleles were recoded per marker to 1, 2 in sorted order, and the read back
  # translates the codes to nucleotides (1 -> A, 2 -> C; C-5): the pattern holds
  for (j in seq_len(dim(m$geno_arrays)[2L])) {
    orig <- m$geno_arrays[, j, ]
    exp  <- c("A", "C", "G", "T")[match(orig, sort(unique(as.vector(orig))))]
    dim(exp) <- dim(orig)
    expect_identical(unname(back$geno_arrays[, j, ]), exp, info = j)
  }
  expect_identical(back$metadata$encoding, "nucleotide")
  expect_equal(back$map_info$Pos,
               as.numeric(m$map_info$Pos[match(back$map_info$Marker, m$map_info$Marker)]))
  expect_identical(back$metadata$source_format, "selectiontools")
  # SelectionTools' block table, as a data.frame and as a file
  mk  <- dimnames(m$geno_arrays)[[2L]]
  blk <- data.frame(Name = c("b000001", "b000002"),
                    Markers = c(paste0(paste(mk[1:3], collapse = ";"), ";"),
                                paste0(paste(mk[5:6], collapse = "; "), ";")),
                    stringsAsFactors = FALSE)
  bmm <- st_convert(blk, map = m)
  expect_identical(bmm, data.frame(Chrom = m$map_info$Chrom[c(1:3, 5:6)],
                                   Marker = mk[c(1:3, 5:6)],
                                   Block = rep(c("b000001", "b000002"), c(3L, 2L)),
                                   stringsAsFactors = FALSE))
  f <- file.path(d, "blocks.ld.txt")
  utils::write.table(blk, f, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_identical(st_convert(f, map = out$map_file)$Marker, bmm$Marker)
  a <- suppressMessages(attach_block_info(m, bmm))
  expect_identical(a$metadata$blocks$n_blocks, 2L)
  expect_identical(a$metadata$blocks$n_singletons, 11L)
  # what it refuses
  expect_error(st_convert(blk), "needs 'map'")
  expect_error(st_convert(data.frame(a = 1), map = m), "Name and Markers")
  expect_error(st_convert(out$geno_file), "needs 'map'")
  blk$Markers[1L] <- "nope;"
  expect_error(st_convert(blk, map = m), "not in 'map'")
  expect_error(st_convert(m), "already a GenoMeta")
  expect_error(st_convert(m, to = "selectiontools"), "output_dir")
})
