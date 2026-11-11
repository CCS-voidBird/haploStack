# =============================================================================
# test-090-gs-missing.R — missing vs unknown at the genomicSimulation boundary
# (owner request 2026-10-27; definitions in R/stacking_gs_bridge.R)
#
#   GM1  a gS true-missing allele (the sentinel) on read-back warns once with
#        code "gs_true_missing", carries the count, and is read as NA; the
#        placeholder alleles ("-", "N") are read as NA silently
#   GM2  (gS installed) the one gene-data reader requests the sentinel, so a
#        true missing and a literal "-" placeholder come back distinguishable
# =============================================================================

test_that("GM1: a gS true missing warns and reads as NA; placeholders read as NA silently", {
  gm <- matrix(c("AT", "A?", "??", "TT",
                 "--", "AN", "TT", "AA"), nrow = 2, byrow = TRUE,
               dimnames = list(c("m1", "m2"), c("S1", "S2", "S3", "S4")))
  r <- haploStack:::.hs_with_notes(haploStack:::.gs_genes_to_array(gm), muffle = TRUE)
  arr <- r$value
  expect_identical(dim(arr), c(4L, 2L, 2L))
  expect_true(all(is.na(arr["S3", "m1", ])))          # "??" true missing
  expect_identical(unname(arr["S2", "m1", ]), c("A", NA))   # "A?" half true missing
  expect_true(all(is.na(arr["S1", "m2", ])))          # "--" placeholder
  expect_identical(unname(arr["S2", "m2", ]), c("A", NA))   # legacy "N"
  expect_identical(unname(arr["S4", "m2", ]), c("A", "A"))
  expect_identical(r$notes$code, "gs_true_missing")
  expect_identical(r$notes$stage, "gs_read_back")
  expect_match(r$notes$message, "3 allele\\(s\\)")
  expect_identical(r$notes$n, 3L)
  # without a true missing: no note at all
  gm2 <- gm; gm2[1, 2:3] <- c("AT", "TT")
  r2 <- haploStack:::.hs_with_notes(haploStack:::.gs_genes_to_array(gm2), muffle = TRUE)
  expect_identical(nrow(r2$notes), 0L)
  expect_true(all(is.na(r2$value["S1", "m2", ])))
  # the warning reaches a plain caller as hs_warning_fallback
  expect_warning(haploStack:::.gs_genes_to_array(gm), class = "hs_warning_fallback")
})

test_that("GM2: the gene-data reader keeps a true missing apart from the '-' placeholder", {
  skip_if_not_installed("genomicSimulation")
  d <- withr::local_tempdir()
  map <- file.path(d, "map.txt"); gen <- file.path(d, "geno.txt")
  writeLines(c("marker chr pos", "m1 1 1.0", "m2 1 2.0"), map)
  # G2's last cell is empty: the only input gS loads as a true missing (a
  # trailing empty cell parses; an interior one would shift the columns)
  writeLines(c("name\tG1\tG2", "m1\tAA\tTT", "m2\t--\t"), gen)
  s <- genomicSimulation::load.data(allele.file = gen, map.file = map)
  on.exit(genomicSimulation::clear.simdata(), add = TRUE)
  g <- haploStack:::.gs_group_genes(s$groupNum)
  expect_identical(unname(g["m2", ]), c("--", "??"))     # placeholder vs true missing
  dflt <- genomicSimulation::see.group.gene.data(s$groupNum)
  expect_identical(unname(dflt["m2", ]), c("--", "--"))  # gS's default would merge them
  r <- haploStack:::.hs_with_notes(haploStack:::.gs_genes_to_array(g), muffle = TRUE)
  expect_true(all(is.na(r$value[, "m2", ])))
  expect_identical(r$notes$n, 2L)
})
