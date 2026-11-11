# =============================================================================
# test-090-deprecation.R — the deprecation ladder (design D5-1, framework F7)
#
#   L1  .hs_deprecated(): class vector, text, fields; function and argument
#   L2  .hs_defunct(): class vector, text, fields
#   L3  the two define_haplotype_blocks* names are defunct (row 26 a)
#   L4  hygiene: no .Deprecated() / .Defunct() outside the ladder file
#   L5  center_by_freq defunct at the one resolver site (row 26 a)
#   L6  the row-27 / Q-c rungs (D5-3): output_format, the two selectors,
#       extend, process/transform_genotype, st_ld_blocks_files, legend
# =============================================================================

test_that("L1: .hs_deprecated() warns with the ladder classes and text", {
  w <- tryCatch(haploStack:::.hs_deprecated("old_fn()", "new_fn()", "0.9.0"),
                warning = identity)
  expect_identical(class(w)[1:3],
                   c("hs_warning_deprecated", "deprecatedWarning", "hs_warning"))
  expect_identical(conditionMessage(w),
                   paste0("'old_fn()' is deprecated since haploStack 0.9.0.\n",
                          "Use 'new_fn()' instead.\nSee help(\"haploStack-deprecated\")."))
  expect_identical(w$old, "old_fn()")
  expect_identical(w$since, "0.9.0")
  expect_null(conditionCall(w))

  w <- tryCatch(haploStack:::.hs_deprecated("old_arg", "new_arg", "0.9.0",
                                            what = "argument", fn = "f"),
                warning = identity)
  expect_identical(conditionMessage(w),
                   "f(old_arg =) is deprecated since haploStack 0.9.0; use f(new_arg =).")
  # it warns on every call: no session memory
  expect_warning(haploStack:::.hs_deprecated("a()", "b()", "0.9.0"), class = "hs_warning_deprecated")
  expect_warning(haploStack:::.hs_deprecated("a()", "b()", "0.9.0"), class = "hs_warning_deprecated")
})

test_that("L2: .hs_defunct() errors with the ladder classes and text", {
  e <- tryCatch(haploStack:::.hs_defunct("old_fn()", "new_fn()", "0.9.0", dep_since = "0.7.10"),
                error = identity)
  expect_identical(class(e)[1:4],
                   c("hs_error_defunct", "defunctError", "hs_error_input", "hs_error"))
  expect_identical(conditionMessage(e),
                   paste0("'old_fn()' is defunct since haploStack 0.9.0 (deprecated 0.7.10). ",
                          "Use 'new_fn()'.\nSee help(\"haploStack-defunct\")."))
  expect_identical(e$new, "new_fn()")
  e <- tryCatch(haploStack:::.hs_defunct("x", "y", "0.9.0", what = "argument", fn = "g"),
                error = identity)
  expect_identical(conditionMessage(e),
                   "g(x =) is defunct since haploStack 0.9.0. Use g(y =).\nSee help(\"haploStack-defunct\").")
})

test_that("L3: define_haplotype_blocks*() are defunct (row 26 a) and name the replacement", {
  e <- tryCatch(define_haplotype_blocks("g.txt", "m.txt"), error = identity)
  expect_s3_class(e, "hs_error_defunct")
  expect_s3_class(e, "defunctError")
  expect_match(conditionMessage(e), "defunct since haploStack 0.9.0 (deprecated 0.7.10)", fixed = TRUE)
  expect_match(conditionMessage(e), "st_ld_blocks()", fixed = TRUE)
  e <- tryCatch(define_haplotype_blocks_meta(list()), error = identity)
  expect_s3_class(e, "hs_error_defunct")
  expect_identical(e$new, "st_ld_blocks()")
})

test_that("L5: center_by_freq is defunct at the one resolver site, names the mode", {
  e <- tryCatch(haploStack:::.resolve_center(NULL, TRUE, "f"), error = identity)
  expect_s3_class(e, "hs_error_defunct")
  expect_match(conditionMessage(e), 'f(center_by_freq =) is defunct since haploStack 0.9.0 (deprecated 0.6.0). Use f(center = "allele_freq" =).', fixed = TRUE)
  e <- tryCatch(haploStack:::.resolve_center(NULL, FALSE, "g"), error = identity)
  expect_identical(e$new, 'center = "none"')
  # the argument is still in the eight signatures (one release)
  for (f in c(run_haplotype_catalogue, compute_profile_localgebv, compute_block_catalogue,
              analyze_block_catalogue, count_block_subtypes, count_favourable_subtypes,
              count_favourable_blocks, predict_gebv))
    expect_true("center_by_freq" %in% names(formals(f)))
  m <- sim_geno_meta(4, 6, n_chrom = 1L, seed = 3L)
  eff <- sim_effect_df(m$map_info$Marker, alleles = c("A", "T"), seed = 3L)
  expect_error(compute_profile_localgebv(m, effect_df = eff, center_by_freq = TRUE),
               class = "hs_error_defunct")
})

test_that("L6: the row-27 / Q-c rungs warn with the ladder class and return the replacement's result", {
  m   <- sim_geno_meta(4, 6, n_chrom = 1L, seed = 3L)
  # output_format -> preset (both writers)
  expect_warning(a <- geno_meta_to_matrix(m, output_format = "gs"),
                 "geno_meta_to_matrix(output_format =) is deprecated", fixed = TRUE,
                 class = "hs_warning_deprecated")
  expect_identical(a$matrix, geno_meta_to_matrix(m, preset = "GS")$matrix)
  expect_identical(a$spec, geno_meta_to_matrix(m, preset = "GS")$spec)
  f <- tempfile(fileext = ".csv")
  expect_warning(suppressMessages(write_geno_meta(m, f, output_format = "gs")),
                 class = "hs_warning_deprecated")
  expect_true(file.exists(f))
  # preset wins when both are given; output_format still warns
  expect_warning(b <- geno_meta_to_matrix(m, preset = "GS", output_format = "plink_ped"),
                 class = "hs_warning_deprecated")
  expect_identical(b$matrix, a$matrix)

  # the two selectors
  lg <- data.table::data.table(Block = paste0("B", 1:5),
                               S1 = c(1, -1, 1, 0, 2), S2 = c(-1, 1, 1, 1, 0),
                               S3 = c(1, 1, -1, -1, 1), S4 = c(0, 0, 0, 1, 1))
  expect_warning(g <- greedy_selection(lg, n_parents = 2L), "ohs_selection(method = \"Greedy\")",
                 fixed = TRUE, class = "hs_warning_deprecated")
  expect_identical(g, haploStack:::.greedy_selection(lg, n_parents = 2L))
  expect_warning(genetic_algorithm_selection(lg, n_parents = 2L, pop_size = 6L,
                                             n_generations = 3L),
                 class = "hs_warning_deprecated")

  # extend = on the two effect functions
  eff <- data.frame(Marker = c("M1", "M1"), Allele = c("A", "T"), Effect = c(0.5, -0.5))
  expect_warning(x <- transform_effects(eff, extend = TRUE),
                 "transform_effects(extend =) is deprecated since haploStack 0.9.0 and is ignored.",
                 fixed = TRUE, class = "hs_warning_deprecated")
  expect_identical(x, transform_effects(eff))
  p <- tempfile(fileext = ".txt"); utils::write.table(eff, p, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_warning(read_effect_file(p, extend = FALSE), class = "hs_warning_deprecated")
  expect_identical(suppressWarnings(read_effect_file(p, extend = FALSE)), read_effect_file(p))

  # process_genotype() / transform_genotype()
  geno <- matrix(c("A/A", "A/T", "T/T", "A/A"), nrow = 2,
                 dimnames = list(c("M1", "M2"), c("S1", "S2")))
  expect_warning(pg <- process_genotype(geno, preset = "GS"), "build_geno_meta()", fixed = TRUE,
                 class = "hs_warning_deprecated")
  # (no map in `pg`, so the call stops after the warning; the warning is what is tested)
  expect_warning(try(transform_genotype(pg, output_preset = "PL"), silent = TRUE),
                 "write_geno_meta()", fixed = TRUE, class = "hs_warning_deprecated")

  # st_ld_blocks_files(): warns before the gate stops it
  withr::local_options(list(haploStack.st_formats = NULL))
  expect_warning(try(st_ld_blocks_files("a", "b"), silent = TRUE),
                 "st_convert()", fixed = TRUE, class = "hs_warning_deprecated")

  # plot_meta_circos(legend =)
  skip_if_not_installed("circlize")
  d <- data.frame(Chrom = "1", Name = "Rrs1", Start = 1e6, End = 2e6)
  expect_warning(ob <- plot_meta_circos(track_marker_labels(d), chr_sizes = c(`1` = 1e7),
                                        draw = FALSE, legend = FALSE),
                 "legend_position", class = "hs_warning_deprecated")
  expect_identical(ob$legend_position, "none")
})

test_that("L4: no .Deprecated() / .Defunct() outside the ladder file", {
  rdir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(rdir), "package R/ sources not reachable from the test dir")
  files <- setdiff(list.files(rdir, pattern = "[.]R$", full.names = TRUE),
                   file.path(rdir, "haploStack-deprecated.R"))
  hits <- unlist(lapply(files, function(f) {
    src <- readLines(f, warn = FALSE)
    i <- grep("\\.(Deprecated|Defunct)\\(", src)
    if (length(i)) paste0(basename(f), ":", i)
  }))
  expect_identical(hits, NULL)
})
