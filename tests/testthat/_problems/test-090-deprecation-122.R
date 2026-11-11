# Extracted from test-090-deprecation.R:122

# test -------------------------------------------------------------------------
m   <- sim_geno_meta(4, 6, n_chrom = 1L, seed = 3L)
expect_warning(a <- geno_meta_to_matrix(m, output_format = "gs"),
                 "geno_meta_to_matrix(output_format =) is deprecated", fixed = TRUE,
                 class = "hs_warning_deprecated")
expect_identical(a$matrix, geno_meta_to_matrix(m, preset = "GS")$matrix)
expect_identical(a$spec, geno_meta_to_matrix(m, preset = "GS")$spec)
f <- tempfile(fileext = ".csv")
expect_warning(suppressMessages(write_geno_meta(m, f, output_format = "gs")),
                 class = "hs_warning_deprecated")
expect_true(file.exists(f))
expect_warning(b <- geno_meta_to_matrix(m, preset = "GS", output_format = "plink_ped"),
                 class = "hs_warning_deprecated")
expect_identical(b$matrix, a$matrix)
lg <- data.table::data.table(Block = paste0("B", 1:5),
                               S1 = c(1, -1, 1, 0, 2), S2 = c(-1, 1, 1, 1, 0),
                               S3 = c(1, 1, -1, -1, 1), S4 = c(0, 0, 0, 1, 1))
expect_warning(g <- greedy_selection(lg, n_parents = 2L), "ohs_selection(method = \"Greedy\")",
                 fixed = TRUE, class = "hs_warning_deprecated")
expect_identical(g, haploStack:::.greedy_selection(lg, n_parents = 2L))
expect_warning(genetic_algorithm_selection(lg, n_parents = 2L, pop_size = 6L,
                                             n_generations = 3L),
                 class = "hs_warning_deprecated")
eff <- data.frame(Marker = c("M1", "M1"), Allele = c("A", "T"), Effect = c(0.5, -0.5))
expect_warning(x <- transform_effects(eff, extend = TRUE),
                 "transform_effects(extend =) is deprecated since haploStack 0.9.0 and is ignored.",
                 fixed = TRUE, class = "hs_warning_deprecated")
expect_identical(x, transform_effects(eff))
p <- tempfile(fileext = ".txt")
utils::write.table(eff, p, sep = "\t", quote = FALSE, row.names = FALSE)
expect_warning(read_effect_file(p, extend = FALSE), class = "hs_warning_deprecated")
expect_identical(suppressWarnings(read_effect_file(p, extend = FALSE)), read_effect_file(p))
geno <- matrix(c("A/A", "A/T", "T/T", "A/A"), nrow = 2,
                 dimnames = list(c("M1", "M2"), c("S1", "S2")))
expect_warning(pg <- process_genotype(geno, preset = "GS"), "build_geno_meta()", fixed = TRUE,
                 class = "hs_warning_deprecated")
expect_warning(transform_genotype(pg, output_preset = "PL"), "write_geno_meta()", fixed = TRUE,
                 class = "hs_warning_deprecated")
