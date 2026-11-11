# =============================================================================
# test-090-ld-helpers.R — the LD-resolver helpers on toy data, no backend
# (design D5 §3.7, deps/F16; step D5-12)
#
#   LH1  .ld_check_rule(): the rule, the backend's zero note
#   LH2  .ld_placed_markers(): which markers a resolver may place
#   LH3  .hapS_chrom_key(): integer keys, natural order, barley unchanged
#   LH4  .hapS_dosage_frame(): dosage of the least frequent allele, NA kept
#   LH5  .st_write_inputs(): aliases, codes, positions in fixed notation,
#        a > 5-state marker left out
#   LH6  .st_parse_blocks(): members back to container names, unknown alias
#   LH7  .ld_finish(): attach = FALSE tables + provenance; attach = TRUE
#        container; files on request
# =============================================================================

.lh_meta <- function() {
  sim_geno_meta(n_samples = 6L, n_markers = 10L, n_chrom = 2L, block_size = 5L,
                encoding = "nucleotide", seed = 21L)
}

test_that("LH1: .ld_check_rule() refuses what the backends cannot take", {
  f <- haploStack:::.ld_check_rule
  expect_true(f(0.5, 0, "x"))
  expect_true(f(1, 3L, "x"))
  expect_error(f(0, 0, "x"), "'ld_threshold' must be one number in (0, 1]", fixed = TRUE)
  expect_error(f(0, 0, "x", zero_note = haploStack:::.st_zero_note),
               "aborts the R session", class = "hs_error_input")
  expect_error(f(1.5, 0, "x"), "ld_threshold")
  expect_error(f(NA, 0, "x"), "ld_threshold")
  expect_error(f(0.5, -1, "x"), "'tolerance' must be one whole number", fixed = TRUE)
  expect_error(f(0.5, 1.5, "x"), "tolerance")
  expect_error(f(0.5, c(1, 2), "x"), "tolerance")
})

test_that("LH2: .ld_placed_markers() keeps finite positive positions on placed chromosomes", {
  map <- data.frame(Marker = paste0("m", 1:8),
                    Chrom  = c("1H", "1H", "Un", "0", "", NA, "chrUn", "2H"),
                    Pos    = c(10, 0, 5, 5, 5, 5, 5, NA),
                    stringsAsFactors = FALSE)
  expect_identical(haploStack:::.ld_placed_markers(map), "m1")
  map2 <- data.frame(Marker = c("a", "b"), Chrom = c("3", "unplaced"), Pos = c("12", "7"))
  expect_identical(haploStack:::.ld_placed_markers(map2), "a")       # text positions parse
})

test_that("LH3: .hapS_chrom_key() numbers chromosomes the way HapSelect block ids embed them", {
  k <- haploStack:::.hapS_chrom_key(c("2H", "1H", "7H", "1H"))
  expect_identical(k$Chrom, c("1H", "2H", "7H"))
  expect_identical(k$Chromosome, c(1L, 2L, 3L))                   # natural order rank
  k2 <- haploStack:::.hapS_chrom_key(c("3", "1", "2"))
  expect_identical(k2$Chrom, c("1", "2", "3"))
  expect_identical(k2$Chromosome, c(1L, 2L, 3L))                  # numeric names keep their value
  k3 <- haploStack:::.hapS_chrom_key(c("chr10", "chr2", "chr1"))
  expect_identical(k3$Chrom, c("chr1", "chr2", "chr10"))          # 10 after 2 (row 29 a)
})

test_that("LH4: .hapS_dosage_frame() counts the least frequent allele and keeps NA", {
  m <- .lh_meta()
  arr <- m$geno_arrays
  # marker 1: five A/A and one A/T -> T is the counted allele; sample 6 missing
  arr[, 1L, 1L] <- "A"; arr[, 1L, 2L] <- "A"; arr[6L, 1L, 2L] <- "T"
  arr[2L, 2L, ] <- NA
  m$geno_arrays <- arr
  key <- haploStack:::.hapS_chrom_key(m$map_info$Chrom)
  mk  <- m$map_info$Marker[1:3]
  d <- haploStack:::.hapS_dosage_frame(m, mk, key)
  expect_true(is.list(d))
  fr <- d$geno %||% d[[1L]]
  expect_true(is.data.frame(fr))
  expect_identical(names(fr)[1:3], c("SNP", "Chromosome", "Position"))
  expect_identical(as.character(fr$SNP), mk)
  dose1 <- as.numeric(fr[fr$SNP == mk[1L], -(1:3)])
  expect_identical(dose1, c(0, 0, 0, 0, 0, 1))                     # one T in sample 6
  dose2 <- as.numeric(fr[fr$SNP == mk[2L], -(1:3)])
  expect_true(is.na(dose2[2L]))                                    # the missing call stays NA
})

test_that("LH5: .st_write_inputs() writes aliases, codes and fixed-notation positions", {
  m <- .lh_meta()
  m$map_info$Pos[1L] <- 88000000
  d <- withr::local_tempdir()
  mk <- m$map_info$Marker
  io <- haploStack:::.st_write_inputs(m, mk, d, alias = TRUE, prefix = "t")
  expect_true(file.exists(io$geno_file) && file.exists(io$map_file))
  g <- readLines(io$geno_file)
  expect_identical(strsplit(g[1L], "\t")[[1L]], c("m/i", paste0("i", 1:6)))
  expect_identical(length(g), length(mk) + 1L)
  expect_true(all(grepl("^m[0-9]+\t([0-9-]/[0-9-]\t)*[0-9-]/[0-9-]$", g[-1L])))
  mp <- readLines(io$map_file)
  expect_identical(mp[1L], "name\tchrom\tpos")
  expect_true(any(grepl("\t88000000$", mp)))                       # no 8.8e+07
  expect_identical(io$markers$Marker, mk[order(match(m$map_info$Chrom, io$chrom_keys$Chrom),
                                               as.numeric(m$map_info$Pos))])
  expect_identical(io$n_left_out, 0L)
  # real names when alias = FALSE; a six-state marker is left out with a warning
  io2 <- haploStack:::.st_write_inputs(m, mk, d, alias = FALSE, prefix = "r")
  expect_identical(strsplit(readLines(io2$geno_file)[1L], "\t")[[1L]][-1L],
                   dimnames(m$geno_arrays)[[1L]])
  m6 <- m; m6$geno_arrays[, 3L, 1L] <- c("A", "C", "G", "T", "R", "Y")
  expect_warning(io3 <- haploStack:::.st_write_inputs(m6, mk, d, prefix = "s"),
                 "more than 5 allele states")
  expect_identical(io3$n_left_out, 1L)
  expect_false(mk[3L] %in% io3$markers$Marker)
  # dict = "fixed" refuses a non-nucleotide state
  expect_error(haploStack:::.st_write_inputs(m6, mk, d, prefix = "f", dict = "fixed"),
               'dict = "fixed"', fixed = TRUE)
})

test_that("LH6: .st_parse_blocks() maps aliases back and refuses unknown ones", {
  alias_map <- data.frame(alias = c("m1", "m2", "m3"), Marker = c("A1", "B2", "C3"),
                          Chrom = c("1H", "1H", "2H"), stringsAsFactors = FALSE)
  raw <- data.frame(Name = c("1:1", "2:1"), Markers = c("m1; m2", "m3"),
                    stringsAsFactors = FALSE)
  bm <- haploStack:::.st_parse_blocks(raw, alias_map)
  expect_identical(bm, data.frame(Chrom = c("1H", "1H", "2H"), Marker = c("A1", "B2", "C3"),
                                  Block = c("1:1", "1:1", "2:1"), stringsAsFactors = FALSE))
  expect_error(haploStack:::.st_parse_blocks(data.frame(Name = "x", Markers = "m9"), alias_map),
               "the bridge did not write", class = "hs_error_input")
  expect_error(haploStack:::.st_parse_blocks(data.frame(Foo = 1), alias_map),
               "unexpected SelectionTools block table", class = "hs_error_internal")
  expect_error(haploStack:::.st_parse_blocks(data.frame(Name = "x", Markers = "m9"), alias_map,
                                             caller = "st_convert"),
               "not in 'map'", class = "hs_error_input")
})

test_that("LH7: .ld_finish() builds the tables, the provenance and the container", {
  m <- .lh_meta()
  map <- m$map_info
  bm <- data.frame(Chrom = map$Chrom[1:4], Marker = map$Marker[1:4],
                   Block = c("1:1", "1:1", "1:2", "1:2"), stringsAsFactors = FALSE)
  map_df <- data.frame(Marker = map$Marker[1:4], Chrom = map$Chrom[1:4],
                       Pos = as.numeric(map$Pos[1:4]), stringsAsFactors = FALSE)
  rule <- list(backend = "toy", ld_threshold = 0.5, tolerance = 2L, n_unplaced = 6L)
  r <- haploStack:::.ld_finish(m, bm, map_df, attach = FALSE, output_dir = NULL,
                               prefix = "p", caller = "toy", rule = rule)
  expect_named(r, c("block_marker_map", "block_info", "provenance"))
  expect_identical(nrow(r$block_info), 2L)
  expect_identical(r$provenance$source, "toy")
  expect_identical(r$provenance$label, "toy_r0.5_t2")
  expect_identical(r$provenance$n_unmatched, 6L)
  d <- withr::local_tempdir()
  r2 <- suppressMessages(haploStack:::.ld_finish(m, bm, map_df, attach = FALSE, output_dir = d,
                                                 prefix = "p", caller = "toy", rule = rule))
  expect_setequal(list.files(d), c("p_block_info.csv", "p_block_marker_map.csv"))
  m2 <- suppressMessages(haploStack:::.ld_finish(m, bm, map_df, attach = TRUE, output_dir = NULL,
                                                 prefix = "p", caller = "toy", rule = rule))
  expect_true(is_geno_meta(m2))
  expect_identical(m2$metadata$blocks$source, "toy")
  expect_identical(sort(unique(as.character(m2$block_marker_map_info$Block)))[1:2], c("1:1", "1:2"))
})
