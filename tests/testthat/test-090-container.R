# =============================================================================
# test-090-container.R — container mechanical fixes (design D2-1)
#
#   K1  char_dict is always a named list of integer codes: a vector is
#       converted, an absent one is derived from the observed alleles
#   K2  a dummy singleton's block_info row carries Pos / Num_SNPs /
#       Block_length when those columns exist (not NA)
#   K3  one annotation check shared by new_geno_meta(), is_geno_meta() and
#       attach_block_annotation()
#   K4  expand_block_info() takes a map, not a container (help and example)
#   K5  .reshape_wide_alleles() refuses a multi-character allele anywhere
#       (cont/F24), with the format class
# =============================================================================

.k090_arrays <- function() {
  s <- c("S1", "S2", "S3"); m <- c("M1", "M2", "M3", "M4")
  h1 <- matrix(c("A", "A", "G", "T",
                 "T", "A", "A", "A",
                 "A", "T", "T", "A"), nrow = 3, byrow = TRUE, dimnames = list(s, m))
  h2 <- matrix(c("A", "T", "G", "T",
                 "T", "A", "T", "A",
                 "T", "T", "G", "A"), nrow = 3, byrow = TRUE, dimnames = list(s, m))
  arr <- array(c(h1, h2), dim = c(dim(h1), 2L), dimnames = c(dimnames(h1), list(c("1", "2"))))
  map <- data.frame(Chrom = c("1", "1", "2", "2"), Marker = m,
                    Pos = c(1e6, 2e6, 3e6, 4e6), stringsAsFactors = FALSE)
  list(arr = arr, map = map)
}
.k090_meta <- function(char_dict) {
  md <- list(ploidy = 2L, n_samples = 3L, n_markers = 4L,
             source_format = "hand_built", encoding = "nucleotide")
  if (!missing(char_dict)) md$char_dict <- char_dict
  md
}

test_that("K1: char_dict is normalised to a named list; derived when absent", {
  d <- .k090_arrays()
  gm <- suppressMessages(new_geno_meta(d$arr, d$map, metadata = .k090_meta(c(A = 1L, G = 2L, T = 3L))))
  expect_identical(gm$metadata$char_dict, list(A = 1L, G = 2L, T = 3L))
  gm <- suppressMessages(new_geno_meta(d$arr, d$map, metadata = .k090_meta(list(A = 1, T = 2))))
  expect_identical(gm$metadata$char_dict, list(A = 1L, T = 2L))
  gm <- suppressMessages(new_geno_meta(d$arr, d$map, metadata = .k090_meta()))
  expect_identical(gm$metadata$char_dict, list(A = 1L, G = 2L, T = 3L))
  expect_true(is_geno_meta(gm))
  expect_error(suppressMessages(new_geno_meta(d$arr, d$map, metadata = .k090_meta(c(1L, 2L)))),
               class = "hs_error_input")
  # every builder gives a list
  expect_true(is.list(sim_geno_meta(4, 6, seed = 1L)$metadata$char_dict))
})

test_that("K2: dummy singleton rows carry their position and length", {
  d <- .k090_arrays()
  bmm <- data.frame(Chrom = c("1", "1"), Marker = c("M1", "M2"), Block = c("B1", "B1"),
                    stringsAsFactors = FALSE)
  bi  <- data.frame(Chrom = "1", Block = "B1", Start = 1e6, End = 2e6, Pos = 1.5e6,
                    Num_SNPs = 2L, Block_length = 1e6 + 1, stringsAsFactors = FALSE)
  gm <- suppressMessages(new_geno_meta(d$arr, d$map, block_info = bi,
                                       block_marker_map_info = bmm,
                                       metadata = .k090_meta(list(A = 1L, G = 2L, T = 3L))))
  dum <- gm$block_info[grepl("^dummy_", gm$block_info$Block), ]
  expect_identical(dum$Block, c("dummy_M3", "dummy_M4"))
  expect_identical(dum$Start, c(3e6, 4e6))
  expect_identical(dum$End, dum$Start)
  expect_identical(dum$Pos, dum$Start)
  expect_identical(dum$Num_SNPs, c(1L, 1L))
  expect_identical(dum$Block_length, c(1, 1))
  # a column the constructor cannot fill stays NA
  bi$Note <- "x"
  gm2 <- suppressMessages(new_geno_meta(d$arr, d$map, block_info = bi,
                                        block_marker_map_info = bmm,
                                        metadata = .k090_meta(list(A = 1L, G = 2L, T = 3L))))
  expect_true(all(is.na(gm2$block_info$Note[grepl("^dummy_", gm2$block_info$Block)])))
})

test_that("K3: one annotation check for the constructor, the validator and the attach", {
  d <- .k090_arrays()
  bmm <- data.frame(Chrom = c("1", "1", "2", "2"), Marker = c("M1", "M2", "M3", "M4"),
                    Block = c("B1", "B1", "B2", "B2"), stringsAsFactors = FALSE)
  md  <- .k090_meta(list(A = 1L, G = 2L, T = 3L))
  good <- list(blocks = data.frame(Block = "B1", n_genes = 1L),
               features = data.frame(Block = "B1", gene_id = "g1"),
               marker_gene = data.frame(Marker = "M1", gene_id = "g1"))
  gm <- suppressMessages(new_geno_meta(d$arr, d$map, block_marker_map_info = bmm,
                                       metadata = md, annotation = good))
  expect_identical(gm$annotation, good)
  expect_true(is_geno_meta(gm))

  bad <- good; bad$blocks$Block <- "B9"
  e <- tryCatch(suppressMessages(new_geno_meta(d$arr, d$map, block_marker_map_info = bmm,
                                               metadata = md, annotation = bad)),
                error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "^new_geno_meta: 1 annotated block\\(s\\) absent from the container's blocks: B9")
  reason <- haploStack:::.hs_check_annotation(bad, c("B1", "B2"), d$map$Marker)
  expect_identical(conditionMessage(e), paste0("new_geno_meta: ", reason))
  expect_null(haploStack:::.hs_check_annotation(good, c("B1", "B2"), d$map$Marker))
  expect_match(haploStack:::.hs_check_annotation(good["blocks"], "B1", "M1"), "only part of the gene-table trio")
  expect_match(haploStack:::.hs_check_annotation(list(custom = data.frame(Marker = "M9")), "B1", "M1"),
               "^1 custom-annotation marker\\(s\\) absent from map_info: M9")

  gm$annotation <- bad
  expect_false(is_geno_meta(gm, verbose = FALSE))
  expect_message(is_geno_meta(gm), "absent from the container's blocks", class = "hs_message")
  gm$annotation <- NULL
  e <- tryCatch(attach_block_annotation(gm, c(bad, list(provenance = list()))), error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "^attach_block_annotation: 1 annotated block\\(s\\) absent from the container's blocks: B9\\. The annotation was computed against different block definitions\\.$")
})

test_that("K4: expand_block_info() documents and takes a map", {
  d <- .k090_arrays()
  coords <- data.frame(Chrom = "1", Block = "B1", Start = 5e5, End = 2.5e6)
  out <- suppressWarnings(expand_block_info(d$map, coords))
  expect_identical(names(out), c("Chrom", "Marker", "Block"))
  expect_identical(out$Marker, c("M1", "M2"))
  rd <- testthat::test_path("..", "..", "man", "expand_block_info.Rd")
  skip_if_not(file.exists(rd), "man/ not reachable")
  txt <- readLines(rd, warn = FALSE)
  expect_true(any(grepl("\\item{map_info}", txt, fixed = TRUE)))
  expect_false(any(grepl("\\item{geno_meta}", txt, fixed = TRUE)))
  expect_true(any(grepl("expand_block_info(meta$map_info, block_coords)", txt, fixed = TRUE)))
})

test_that("K5: a multi-character allele anywhere in a wide matrix is refused", {
  m <- matrix("A", nrow = 5, ncol = 4, dimnames = list(paste0("S", 1:5), NULL))
  expect_identical(dim(haploStack:::.reshape_wide_alleles(m, 2L, 2L)), c(5L, 2L))
  m[5, 3] <- "AT"                                    # row 4+: passed up to 0.8.8
  e <- tryCatch(haploStack:::.reshape_wide_alleles(m, 2L, 2L), error = identity)
  expect_s3_class(e, "hs_error_format")
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "multi-character")
})


# =============================================================================
# Phase B (design D2-2 / D2-3): one rebuild path and the cheap validator
#
#   C1  carry-forward matrix: every transform keeps FAM, effects, phenotype
#       and annotation, or drops exactly what it declares
#   C2  build_geno_meta(fam =, blocks =) and read_genotype(block_file =) keep
#       the FAM; attach_block_info() and an LD resolver carry the same fields
#   C3  counters and ploidy are derived; a hand-set value is ignored by the
#       constructor and fails is_geno_meta() on a hand-edited object
#   C4  is_geno_meta() does not scan a stamped container; it does under
#       strict = TRUE and on an unstamped one
#   C5  a reversed map fails; filter_geno_meta(keep_markers = <reordered>)
#       keeps map and array aligned; map order follows the array
#   C8  scope rules: signal_selection cleared on re-blocking, its n_blocks
#       refreshed on filter; family_history dropped on offspring, kept on
#       filter; no date in signal_selection
#   C9  invalid annotation -> hs_error_input; filter prunes it; re-blocking
#       keeps custom with Block re-derived
#   C12 row 6: a fully filtered block's row is gone
# =============================================================================

.b090_full <- function(seed = 11L) {
  m <- sim_geno_meta(n_samples = 8L, n_markers = 16L, block_size = 4L, seed = seed)
  m$marker_effects <- list(T1 = sim_effect_df(dimnames(m$geno_arrays)[[2L]],
                                              alleles = c("A", "T"), seed = seed))
  ph <- data.frame(Sample = dimnames(m$geno_arrays)[[1L]], T1 = seq_len(8L))
  fam <- data.frame(Family = rep(c("F1", "F2"), each = 4L),
                    Sample = dimnames(m$geno_arrays)[[1L]],
                    Parent1 = c("0", "0", "S1", "S1", "0", "0", "S5", "S5"),
                    Parent2 = c("0", "0", "S2", "S2", "0", "0", "S6", "S6"),
                    Sex = "0", stringsAsFactors = FALSE)
  m <- suppressMessages(new_geno_meta(m$geno_arrays, m$map_info,
                                      block_info = m$block_info,
                                      block_marker_map_info = m$block_marker_map_info,
                                      metadata = m$metadata, sample_info = fam,
                                      marker_effects = m$marker_effects, phenotype = ph))
  m <- suppressWarnings(attach_family_history(m))
  m
}
.b090_same <- function(a, b, fields = c("sample_info", "marker_effects", "phenotype"))
  vapply(fields, function(f) identical(a[[f]], b[[f]]), logical(1L))

test_that("C1: every transform forwards FAM, effects, phenotype and family history", {
  m <- .b090_full()
  expect_true(is_geno_meta(m))
  expect_false(is.null(m$metadata$family_history))
  # keep everything
  p2 <- suppressWarnings(attach_phenotype(m, m$phenotype))
  expect_true(all(.b090_same(m, p2, c("sample_info", "marker_effects"))))
  expect_identical(p2$metadata$family_history, m$metadata$family_history)
  r2 <- suppressMessages(recode_geno_meta(m, "nucleotide"))
  expect_true(all(.b090_same(m, r2)))
  expect_identical(r2$metadata$family_history, m$metadata$family_history)
  expect_identical(r2$block_marker_map_info, m$block_marker_map_info)
  b2 <- suppressMessages(attach_block_info(m, sim_block_map(m$map_info, block_size = 2L)))
  expect_true(all(.b090_same(m, b2)))
  expect_identical(b2$metadata$family_history, m$metadata$family_history)
  expect_identical(length(unique(b2$block_marker_map_info$Block)), 8L)
  # declared losses: offspring drop phenotype and family history, keep effects
  cd  <- make_cross_design(c("S1", "S2"), n_offspring = 2L)
  off <- haploStack:::.cross_block(m, cd)$geno_meta
  expect_null(off$phenotype)
  expect_null(off$metadata$family_history)
  expect_identical(off$marker_effects, m$marker_effects)
  expect_identical(off$block_marker_map_info, m$block_marker_map_info)
  expect_true(all(off$sample_info$Parent1 == "0"))
  gm <- t(apply(m$geno_arrays[1:2, , ], c(1, 2), paste, collapse = ""))
  colnames(gm) <- c("K1", "K2")
  up <- update_geno_meta_from_matrix(m, gm)
  expect_null(up$phenotype)
  expect_null(up$metadata$family_history)
  expect_identical(up$marker_effects, m$marker_effects)
  expect_identical(dimnames(up$geno_arrays)[[1L]], c("K1", "K2"))
  # the rebuild helper refuses what it does not know
  expect_error(haploStack:::.rebuild_geno_meta(m, nonsense = 1), class = "hs_error_internal")
  expect_error(haploStack:::.rebuild_geno_meta(m, .drop = "nonsense"), class = "hs_error_internal")
})

test_that("C2: the FAM survives block attachment on every route", {
  m <- .b090_full()
  geno <- haploStack:::.meta_to_allele_matrix(m)
  bmap <- sim_block_map(m$map_info, block_size = 4L)
  b <- suppressMessages(build_geno_meta(geno, m$map_info, fam = m$sample_info,
                                        blocks = bmap, split_char = "/",
                                        encoding = "nucleotide"))
  expect_identical(b$sample_info, m$sample_info)
  expect_identical(length(unique(b$block_marker_map_info$Block)), 4L)
  # the same through a PED read with a block file
  d <- tempfile("ped_"); dir.create(d)
  suppressMessages(write_geno_meta(recode_geno_meta(m, "nucleotide"),
                                   file.path(d, "x"), preset = "PL"))
  bf <- file.path(d, "blocks.txt")
  utils::write.table(bmap, bf, sep = "\t", quote = FALSE, row.names = FALSE)
  r <- suppressMessages(suppressWarnings(
    read_plink_ped(file.path(d, "x.ped"), file.path(d, "x.map"), block_file = bf)))
  expect_identical(r$sample_info$Parent1, m$sample_info$Parent1)
  expect_identical(r$sample_info$Family, m$sample_info$Family)
  expect_identical(length(unique(r$block_marker_map_info$Block)), 4L)
})

test_that("C3: counters and ploidy are derived, never trusted", {
  m <- .b090_full()
  md <- m$metadata; md$n_markers <- 999L; md$ploidy <- 7L; md$n_samples <- 1L
  gm <- suppressMessages(new_geno_meta(m$geno_arrays, m$map_info, metadata = md))
  expect_identical(gm$metadata$n_markers, 16L)
  expect_identical(gm$metadata$n_samples, 8L)
  expect_identical(gm$metadata$ploidy, 2L)
  md2 <- list(source_format = "x", encoding = "nucleotide")   # the minimum
  gm2 <- suppressMessages(new_geno_meta(m$geno_arrays, m$map_info, metadata = md2))
  expect_true(is_geno_meta(gm2))
  expect_error(suppressMessages(new_geno_meta(m$geno_arrays, m$map_info, metadata = list())),
               "source_format", class = "hs_error_input")
  # a hand edit fails the cheap check, with the reason
  h <- m; h$metadata$n_markers <- 999L
  expect_false(is_geno_meta(h, verbose = FALSE))
  expect_message(is_geno_meta(h), "metadata\\$n_markers \\(999\\)", class = "hs_message")
  h <- m; h$metadata$ploidy <- NULL
  expect_message(is_geno_meta(h), "metadata\\$ploidy is missing", class = "hs_message")
})

test_that("C4: a stamped container is not scanned; strict and unstamped are", {
  m <- .b090_full()
  calls <- 0L
  testthat::local_mocked_bindings(
    .missing_token_report = function(x, source_format = NULL) { calls <<- calls + 1L; NULL },
    .package = "haploStack")
  expect_true(is_geno_meta(m))
  expect_identical(calls, 0L)
  expect_true(is_geno_meta(m, strict = TRUE))
  expect_identical(calls, 1L)
  u <- m; u$metadata$missing_tokens <- NULL
  expect_true(is_geno_meta(u))
  expect_identical(calls, 2L)
  haploStack:::.hs_check_meta(m, "f")
  expect_identical(calls, 2L)
  haploStack:::.hs_check_meta(u, "f")
  expect_identical(calls, 3L)
})

test_that("C5: map order follows the array; a reversed map fails", {
  m <- .b090_full()
  rev <- m; rev$map_info <- rev$map_info[nrow(rev$map_info):1, ]
  expect_false(is_geno_meta(rev, verbose = FALSE))
  expect_message(is_geno_meta(rev), "array's marker order", class = "hs_message")
  mk <- dimnames(m$geno_arrays)[[2L]]
  keep <- mk[c(12, 3, 7, 1)]                       # a reordered request
  f <- suppressWarnings(filter_geno_meta(m, keep_markers = keep))
  expect_identical(dimnames(f$geno_arrays)[[2L]], mk[c(1, 3, 7, 12)])
  expect_identical(as.character(f$map_info$Marker), dimnames(f$geno_arrays)[[2L]])
  expect_true(is_geno_meta(f, strict = TRUE))
  expect_identical(f$metadata$n_markers, 4L)
})

test_that("C8: block-scoped and sample-scoped metadata follow the scope rules", {
  m <- .b090_full()
  sig <- data.frame(Block = unique(m$block_marker_map_info$Block)[1:2],
                    Marker = dimnames(m$geno_arrays)[[2L]][c(1, 5)],
                    Score = c(0.9, 0.8), stringsAsFactors = FALSE)
  s <- suppressMessages(attach_signal_markers(m, sig))
  expect_identical(s$metadata$signal_selection$n_blocks, 2L)
  expect_null(s$metadata$signal_selection$date)
  expect_identical(s$metadata$family_history, m$metadata$family_history)
  # filter: the selection is pruned and its count refreshed; history kept
  f <- suppressMessages(suppressWarnings(
    filter_geno_meta(s, keep_markers = dimnames(m$geno_arrays)[[2L]][5:16])))
  expect_identical(f$metadata$signal_selection$n_blocks, 1L)
  expect_identical(f$metadata$family_history, s$metadata$family_history)
  # re-blocking clears the selection; history (sample-scoped) stays
  b <- suppressMessages(attach_block_info(s, sim_block_map(s$map_info, block_size = 8L)))
  expect_null(b$metadata$signal_selection)
  expect_false(any(haploStack:::.signal_marker_cols() %in% names(b$block_info)))
  expect_identical(b$metadata$family_history, s$metadata$family_history)
  # identical inputs give identical containers (no timestamps)
  s2 <- suppressMessages(attach_signal_markers(m, sig))
  expect_identical(s, s2)
})

test_that("C9: annotation is validated, pruned on filter, custom kept on re-blocking", {
  m <- .b090_full()
  blocks <- unique(m$block_marker_map_info$Block)
  ann <- list(blocks = data.frame(Block = blocks[1:2], gene_count = 1L),
              features = data.frame(Block = blocks[1:2], gene_id = c("g1", "g2")),
              marker_gene = data.frame(Marker = dimnames(m$geno_arrays)[[2L]][c(1, 5)],
                                       Block = blocks[1:2], gene_id = c("g1", "g2")))
  a <- suppressMessages(new_geno_meta(m$geno_arrays, m$map_info,
                                      block_marker_map_info = m$block_marker_map_info,
                                      metadata = m$metadata, annotation = ann))
  bad <- ann; bad$blocks$Block[1] <- "NOPE"
  expect_error(suppressMessages(new_geno_meta(m$geno_arrays, m$map_info,
                                              block_marker_map_info = m$block_marker_map_info,
                                              metadata = m$metadata, annotation = bad)),
               class = "hs_error_input")
  f <- suppressMessages(suppressWarnings(
    filter_geno_meta(a, keep_markers = dimnames(m$geno_arrays)[[2L]][5:16])))
  expect_identical(f$annotation$blocks$Block, blocks[2])
  expect_identical(f$annotation$marker_gene$gene_id, "g2")
  custom <- data.frame(Marker = dimnames(m$geno_arrays)[[2L]][c(2, 9)], Name = c("L1", "L2"))
  c1 <- suppressMessages(attach_marker_annotations(a, custom))
  expect_setequal(names(c1$annotation), c("blocks", "features", "marker_gene", "custom"))
  nb <- sim_block_map(m$map_info, block_size = 8L)
  expect_message(b <- attach_block_info(c1, nb), "custom marker annotations are kept")
  expect_identical(names(b$annotation), "custom")
  expect_identical(b$annotation$custom$Block, nb$Block[match(b$annotation$custom$Marker, nb$Marker)])
  expect_identical(names(b$metadata$annotation), "custom")
  expect_true(is_geno_meta(b, strict = TRUE))
})

test_that("C12: a fully filtered block loses its block_info row (row 6)", {
  m <- .b090_full()
  blocks <- unique(m$block_marker_map_info$Block)
  keep <- m$block_marker_map_info$Marker[m$block_marker_map_info$Block != blocks[1]]
  f <- suppressWarnings(filter_geno_meta(m, keep_markers = keep))
  expect_false(blocks[1] %in% f$block_info$Block)
  expect_setequal(f$block_info$Block, blocks[-1])
  expect_true(is_geno_meta(f, strict = TRUE))
})
