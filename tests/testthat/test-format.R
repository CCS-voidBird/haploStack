# =============================================================================
# test-format.R
#
# QA tests for the format module:
#   detect_genotype_encoding(), split_haplotypes(), combine_haplotypes(),
#   translate_alleles(), process_genotype(), transform_effects(),
#   read_dart_genotype()
# =============================================================================

library(testthat)
library(haploStack)

# ── detect_genotype_encoding ──────────────────────────────────────────────────

test_that("detect_genotype_encoding refuses SelectionTools codes (0.9.0, C-5)", {
  geno <- matrix(c("1/1","1/4","4/4","-1/-1"), nrow = 2)
  expect_error(detect_genotype_encoding(geno), class = "hs_error_format")
  expect_equal(haploStack:::.detect_encoding(geno), "selectiontools")
})

test_that("detect_genotype_encoding correctly identifies nucleotide slash format", {
  geno <- matrix(c("A/A","A/T","T/T","G/C"), nrow = 2)
  expect_equal(detect_genotype_encoding(geno), "nucleotide")
})

test_that("detect_genotype_encoding correctly identifies dosage numeric format", {
  geno <- matrix(c("0","1","2","0","1","2"), nrow = 3)
  expect_equal(detect_genotype_encoding(geno), "numeric")
})

test_that("detect_genotype_encoding errors on all-NA input", {
  geno <- matrix(rep(NA_character_, 4), nrow = 2)
  expect_error(detect_genotype_encoding(geno), "No non-NA")
})

test_that("detect_genotype_encoding handles nucleotide no-slash format", {
  geno <- matrix(c("AA","AT","TT","AA"), nrow = 2)
  expect_equal(detect_genotype_encoding(geno), "nucleotide")
})

test_that("detect_genotype_encoding accepts data.frame input", {
  df <- as.data.frame(matrix(c("A/A","T/T"), nrow = 1))
  expect_equal(detect_genotype_encoding(df), "nucleotide")
})

# ── split_haplotypes ──────────────────────────────────────────────────────────

test_that("split_haplotypes returns correct 3D array dimensions", {
  geno <- matrix(c("1/1","1/4","4/4","1/1"), nrow = 2,
                 dimnames = list(c("M1","M2"), c("S1","S2")))
  res <- split_haplotypes(geno, orientation = 2L)
  expect_equal(dim(res$array), c(2L, 2L, 2L))
  expect_equal(res$samples, c("S1","S2"))
  expect_equal(res$markers, c("M1","M2"))
})

test_that("split_haplotypes correctly separates alleles", {
  geno <- matrix("1/4", nrow = 1, ncol = 1,
                 dimnames = list("M1","S1"))
  res <- split_haplotypes(geno, orientation = 2L)
  expect_equal(res$array["S1","M1","H1"], "1")
  expect_equal(res$array["S1","M1","H2"], "4")
})

test_that("split_haplotypes handles nucleotide no-delimiter format", {
  geno <- matrix("AT", nrow = 1, ncol = 1,
                 dimnames = list("M1","S1"))
  res <- split_haplotypes(geno, split_char = "", orientation = 2L)
  expect_equal(res$array["S1","M1","H1"], "A")
  expect_equal(res$array["S1","M1","H2"], "T")
})

test_that("split_haplotypes auto-detects orientation from aspect ratio", {
  # More rows than cols → rows = markers
  geno <- matrix(paste0("1/", c(1,4,4,1,4,4)), nrow = 6, ncol = 2,
                 dimnames = list(paste0("M",1:6), c("S1","S2")))
  res <- split_haplotypes(geno)
  expect_equal(length(res$samples), 2L)
  expect_equal(length(res$markers), 6L)
})

# ── combine_haplotypes ────────────────────────────────────────────────────────

test_that("combine_haplotypes round-trips split_haplotypes", {
  geno <- matrix(c("1/1","1/4","4/4","4/1"), nrow = 2,
                 dimnames = list(c("M1","M2"), c("S1","S2")))
  split_res <- split_haplotypes(geno, orientation = 2L)
  combined  <- combine_haplotypes(split_res$array, split_res$samples)
  # Should have same dimensions as original (transposed back)
  expect_equal(nrow(combined), 2L)  # samples as rows
  expect_equal(ncol(combined), 2L)  # markers as cols
})

test_that("combine_haplotypes respects output_orientation = 2", {
  geno <- matrix("1/4", nrow = 1, dimnames = list("M1","S1"))
  split_res <- split_haplotypes(geno, orientation = 2L)
  out <- combine_haplotypes(split_res$array, "S1", output_orientation = 2L)
  expect_equal(nrow(out), 1L)  # markers as rows
})

# ── translate_alleles ─────────────────────────────────────────────────────────

test_that("translate_alleles converts nucleotide to selectiontools correctly", {
  withr::local_options(list(haploStack.st_formats = TRUE))        # 0.9.0, C-6
  arr <- array(c("A","T","T","A"), dim = c(1,2,2),
               dimnames = list("S1", c("M1","M2"), c("H1","H2")))
  expect_warning(out <- translate_alleles(arr, "nucleotide", "selectiontools"),
                 class = "hs_warning_deprecated")             # 0.9.0, C-5
  expect_equal(as.character(out["S1","M1","H1"]), "1")  # A → 1
  expect_equal(as.character(out["S1","M1","H2"]), "4")  # T → 4
})

test_that("translate_alleles returns unchanged array when encodings match", {
  arr <- array("1", dim = c(1,1,1))
  out <- translate_alleles(arr, "selectiontools", "selectiontools")
  expect_identical(arr, out)
})

test_that("translate_alleles converts selectiontools to nucleotide correctly", {
  withr::local_options(list(haploStack.st_formats = TRUE))        # 0.9.0, C-6
  arr <- array(c("1","4"), dim = c(1,1,2))
  out <- suppressWarnings(translate_alleles(arr, "selectiontools", "nucleotide"))
  expect_equal(as.character(out[1,1,1]), "A")  # 1 → A
  expect_equal(as.character(out[1,1,2]), "T")  # 4 → T
})

# ── the retired "char" / "number" vocabulary (0.5.1) ─────────────────────────

test_that("the removed encoding tokens error, naming their replacement", {
  arr <- array(c("A","T"), dim = c(1,1,2))
  # Removed at every layer that used to accept them, and the message says what
  # to write instead — the old failure mode was a confusing "unsupported
  # encoding" from whichever helper happened to see the token first, or (worse)
  # a GenoMeta that built fine and was rejected later by geno_meta_views().
  expect_error(translate_alleles(arr, "char", "number"),
               'was removed in 0.5.1')
  expect_error(translate_alleles(arr, "char", "nucleotide"), 'use "nucleotide"')
  expect_error(translate_alleles(arr, "nucleotide", "number"),
               'use "selectiontools"')

  geno <- matrix(c("1/1", "1/4", "4/4", "1/1"), nrow = 2L,
                 dimnames = list(c("S1", "S2"), c("M1", "M2")))
  map  <- data.frame(Chrom = c("1", "1"), Marker = c("M1", "M2"),
                     Pos = c(1e6, 2e6), stringsAsFactors = FALSE)
  expect_error(build_geno_meta(geno, map, encoding = "number"),
               'use "selectiontools"')
  expect_error(build_geno_meta(geno, map, encoding = "char"),
               'use "nucleotide"')
  expect_error(suppressWarnings(process_genotype(geno, encode = "number"),
                                classes = "hs_warning_deprecated"),   # row 27 a
               'use "selectiontools"')

  # An unknown (never-supported) token is a different, plainer error.
  expect_error(translate_alleles(arr, "nucleotide", "banana"),
               "unknown encoding")
})

test_that("build_geno_meta stores an encoding geno_meta_views() accepts", {
  # The 0.5.1 fix: metadata$encoding is now always from the canonical
  # vocabulary, so the views layer can consume whatever build_geno_meta emits.
  # Since 0.9.0 (C-5) that vocabulary is nucleotide / numeric; a declared
  # "selectiontools" is translated (with a deprecation warning).
  geno <- matrix(c("A/A", "A/T", "T/T", "A/A"), nrow = 2L,
                 dimnames = list(c("S1", "S2"), c("M1", "M2")))
  map  <- data.frame(Chrom = c("1", "1"), Marker = c("M1", "M2"),
                     Pos = c(1e6, 2e6), stringsAsFactors = FALSE)
  # list(), not c(), so the NULL (auto-detect) case is actually exercised.
  for (enc in list("nucleotide", NULL)) {
    meta <- suppressMessages(build_geno_meta(geno, map, encoding = enc))
    expect_true(meta$metadata$encoding %in% c("nucleotide", "numeric"))
    expect_no_error(geno_meta_views(meta))
  }
  codes <- matrix(c("1/1", "1/4", "4/4", "1/1"), nrow = 2L,
                  dimnames = list(c("S1", "S2"), c("M1", "M2")))
  withr::local_options(list(haploStack.st_formats = TRUE))        # 0.9.0, C-6
  expect_warning(meta <- suppressMessages(build_geno_meta(codes, map, encoding = "selectiontools")),
                 class = "hs_warning_deprecated")
  expect_identical(meta$metadata$encoding, "nucleotide")
  expect_no_error(geno_meta_views(meta))
})

# ── transform_effects ─────────────────────────────────────────────────────────

test_that("transform_effects standardises column names by position", {
  df <- data.frame(x = c("M1","M1"), y = c("A","T"), z = c(0.5,-0.5))
  out <- transform_effects(df)
  expect_true("Marker" %in% names(out))
  expect_true("Allele" %in% names(out))
  expect_true("Effect" %in% names(out))
})

test_that("transform_effects coerces Effect to numeric", {
  df <- data.frame(Marker = "M1", Allele = "A", Effect = "0.5")
  out <- transform_effects(df)
  expect_type(out$Effect, "double")
})

test_that("transform_effects raises error on < 3 columns", {
  df <- data.frame(x = "M1", y = "A")
  expect_error(transform_effects(df), "3 columns")
})

# ── read_dart_genotype ────────────────────────────────────────────────────────

# Helper: build a minimal in-memory DArT CSV and write to tmp file
.dart_tmp <- local({
  set.seed(10L)
  n_mk  <- 20L
  n_smp <- 8L
  # Simulate dosage matrix (0/1/2); first 2 markers low call rate (force filter)
  dos <- matrix(sample(c(0L, 1L, 2L), n_mk * n_smp, replace = TRUE),
                nrow = n_mk, ncol = n_smp)
  # Force call rate filter for marker 1 (set 80% to NA — drops below 0.8 threshold)
  dos[1L, seq_len(ceiling(n_smp * 0.85))] <- NA
  # Force MAF filter for marker 2 (set all to 0 → maf = 0)
  dos[2L, ] <- 0L

  df <- as.data.frame(dos, stringsAsFactors = FALSE)
  colnames(df) <- paste0("Sample", seq_len(n_smp))

  # Add metadata columns
  meta <- data.frame(
    AlleleID = paste0("Mkr", seq_len(n_mk)),
    Chrom    = rep(c("1H","2H"), each = n_mk / 2L),
    ChromPosSnp = seq(1e6L, by = 5e5L, length.out = n_mk),
    stringsAsFactors = FALSE
  )
  full_df <- cbind(meta, df)
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(full_df, tmp, row.names = FALSE)
  list(tmp = tmp, n_mk = n_mk, n_smp = n_smp, dos = dos)
})

test_that("read_dart_genotype returns a valid GenoMeta", {
  gm <- suppressMessages(read_dart_genotype(.dart_tmp$tmp))
  expect_true(is_geno_meta(gm))
  expect_equal(gm$metadata$source_format, "dart_seq")
})

test_that("read_dart_genotype geno_arrays has correct dimensions", {
  gm <- suppressMessages(read_dart_genotype(.dart_tmp$tmp))
  expect_equal(length(dim(gm$geno_arrays)), 3L)
  expect_equal(dim(gm$geno_arrays)[1L], .dart_tmp$n_smp)
  expect_equal(dim(gm$geno_arrays)[3L], 2L)  # diploid
})

test_that("read_dart_genotype correctly splits dosage to A/B haplotypes", {
  # Build a tiny hand-crafted DArT CSV with known dosage values
  tiny_df <- data.frame(
    AlleleID    = c("Mkr_AA","Mkr_AB","Mkr_BB"),
    Chrom       = c("1H","1H","2H"),
    ChromPosSnp = c(100000L, 200000L, 300000L),
    S1 = c(2L, 1L, 0L),
    S2 = c(2L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  tmp2 <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp2), add = TRUE)
  utils::write.csv(tiny_df, tmp2, row.names = FALSE)

  gm  <- suppressMessages(read_dart_genotype(tmp2, maf_threshold = 0,
                                               call_rate_threshold = 0))
  arr <- gm$geno_arrays

  # Mkr_AA (dosage 2): both haplotypes should be "A"
  expect_equal(arr["S1", "Mkr_AA", 1L], "A")
  expect_equal(arr["S1", "Mkr_AA", 2L], "A")

  # Mkr_BB (dosage 0): both haplotypes should be "B"
  expect_equal(arr["S1", "Mkr_BB", 1L], "B")
  expect_equal(arr["S1", "Mkr_BB", 2L], "B")

  # Mkr_AB (dosage 1): hap1 = "A", hap2 = "B"
  expect_equal(arr["S1", "Mkr_AB", 1L], "A")
  expect_equal(arr["S1", "Mkr_AB", 2L], "B")
})

test_that("read_dart_genotype map_info has Chr, Marker, Pos columns", {
  gm <- suppressMessages(read_dart_genotype(.dart_tmp$tmp))
  expect_true(all(c("Chrom","Marker","Pos") %in% names(gm$map_info)))
})

test_that("read_dart_genotype removes low-call-rate markers", {
  n_before <- .dart_tmp$n_mk
  gm <- suppressMessages(read_dart_genotype(.dart_tmp$tmp,
                                             call_rate_threshold = 0.8))
  expect_true(gm$metadata$n_markers < n_before)
})

test_that("read_dart_genotype removes MAF-zero markers", {
  gm_strict <- suppressMessages(
    read_dart_genotype(.dart_tmp$tmp, call_rate_threshold = 0, maf_threshold = 0.01)
  )
  # Marker 2 (all-zero dosage → MAF = 0) must be filtered
  expect_false("Mkr2" %in% gm_strict$map_info$Marker)
})

test_that("read_dart_genotype errors on non-existent file", {
  expect_error(read_dart_genotype("/no/such/file.csv"), "not found")
})

# ── read_vcf_genotype ─────────────────────────────────────────────────────────

# Helper: write a minimal VCF to a temp file
.vcf_tmp <- local({
  vcf_lines <- c(
    "##fileformat=VCFv4.1",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
    # Mkr_AA: all homozygous ref (dosage 0) → MAF=0, filtered by MAF
    "1H\t100000\tMkr_ref\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0",
    # Mkr_AB: heterozygous (dosage 1) → MAF=0.5, passes
    "1H\t200000\tMkr_het\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/1\t1/0\t1/1",
    # Mkr_AA2: all alt (dosage 2) → MAF=0, filtered by MAF
    "2H\t300000\tMkr_alt\tA\tT\t.\tPASS\t.\tGT\t1/1\t1/1\t1/1\t1/1",
    # Mkr_mixed: mixed dosage → MAF > 0, passes
    "2H\t400000\tMkr_mix\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1",
    # Mkr_lowcall: low call rate → filtered
    "3\t500000\tMkr_lc\tA\tT\t.\tPASS\t.\tGT\t./.\t./.\t./.\t0/1",
    # Mkr_multi: multi-allelic → filtered by biallelic_only
    "3\t600000\tMkr_ma\tA\tT,C\t.\tPASS\t.\tGT\t0/1\t1/2\t0/2\t0/0",
    # Mkr_phased: phased GT (use | separator) → should parse correctly
    "1\t700000\tMkr_ph\tA\tT\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0\t0|1"
  )
  tmp <- tempfile(fileext = ".vcf")
  writeLines(vcf_lines, tmp)
  tmp
})

test_that("read_vcf_genotype returns a valid GenoMeta", {
  skip_if_not_installed("vcfR")
  gm <- suppressMessages(read_vcf_genotype(.vcf_tmp))
  expect_true(is_geno_meta(gm))
  expect_equal(gm$metadata$source_format, "vcf")
})

test_that("read_vcf_genotype geno_arrays is a 3D array with 2 haplotype slices", {
  skip_if_not_installed("vcfR")
  gm <- suppressMessages(read_vcf_genotype(.vcf_tmp))
  expect_equal(length(dim(gm$geno_arrays)), 3L)
  expect_equal(dim(gm$geno_arrays)[3L], 2L)
})

test_that("read_vcf_genotype map_info has Chr, Marker, Pos columns", {
  skip_if_not_installed("vcfR")
  gm <- suppressMessages(read_vcf_genotype(.vcf_tmp))
  expect_true(all(c("Chrom","Marker","Pos") %in% names(gm$map_info)))
})

test_that("read_vcf_genotype strips chromosome suffixes", {
  skip_if_not_installed("vcfR")
  gm <- suppressMessages(read_vcf_genotype(.vcf_tmp))
  # 1H and 2H should become "1" and "2"
  expect_false(any(grepl("H$", gm$map_info$Chrom)))
})

test_that("read_vcf_genotype removes multi-allelic sites", {
  skip_if_not_installed("vcfR")
  gm <- suppressMessages(read_vcf_genotype(.vcf_tmp, maf_threshold = 0,
                                             call_rate_threshold = 0))
  expect_false("Mkr_ma" %in% gm$map_info$Marker)
})

test_that("read_vcf_genotype removes low-call-rate markers", {
  skip_if_not_installed("vcfR")
  gm <- suppressMessages(read_vcf_genotype(.vcf_tmp, call_rate_threshold = 0.5))
  expect_false("Mkr_lc" %in% gm$map_info$Marker)
})

test_that("read_vcf_genotype correctly splits dosage to A/B haplotypes", {
  skip_if_not_installed("vcfR")
  # Tiny 3-marker VCF with known dosage: 0, 1, 2
  tiny_vcf <- c(
    "##fileformat=VCFv4.1",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "1\t100\tMkr_BB\tA\tT\t.\t.\t.\tGT\t0/0\t0/0",
    "1\t200\tMkr_AB\tA\tT\t.\t.\t.\tGT\t0/1\t1/0",
    "1\t300\tMkr_AA\tA\tT\t.\t.\t.\tGT\t1/1\t1/1"
  )
  tmp2 <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp2), add = TRUE)
  writeLines(tiny_vcf, tmp2)

  gm  <- suppressMessages(read_vcf_genotype(tmp2, maf_threshold = 0,
                                              call_rate_threshold = 0))
  arr <- gm$geno_arrays

  # Mkr_BB (dosage 0): both haplotypes = "B"
  expect_equal(arr["S1", "Mkr_BB", 1L], "B")
  expect_equal(arr["S1", "Mkr_BB", 2L], "B")

  # Mkr_AB (dosage 1): hap1 = "A", hap2 = "B"
  expect_equal(arr["S1", "Mkr_AB", 1L], "A")
  expect_equal(arr["S1", "Mkr_AB", 2L], "B")

  # Mkr_AA (dosage 2): both haplotypes = "A"
  expect_equal(arr["S1", "Mkr_AA", 1L], "A")
  expect_equal(arr["S1", "Mkr_AA", 2L], "A")
})

test_that("read_vcf_genotype handles phased GT (| separator)", {
  skip_if_not_installed("vcfR")
  phased_vcf <- c(
    "##fileformat=VCFv4.1",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "1\t100\tMkr1\tA\tT\t.\t.\t.\tGT\t0|1\t1|1\t0|0"
  )
  tmp3 <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp3), add = TRUE)
  writeLines(phased_vcf, tmp3)

  gm <- suppressMessages(read_vcf_genotype(tmp3, maf_threshold = 0,
                                             call_rate_threshold = 0))
  expect_true(is_geno_meta(gm))
  expect_equal(dim(gm$geno_arrays)[1L], 3L)  # 3 samples
})

test_that("read_vcf_genotype errors on non-existent file", {
  skip_if_not_installed("vcfR")
  expect_error(read_vcf_genotype("/no/such/file.vcf"), "not found")
})

# ── build_geno_meta automatic marker intersection (geno ∩ map) ────────────────

.bgm_geno <- matrix(
  c("A/A","A/T","T/T","A/A",
    "T/T","A/T","A/A","T/T",
    "A/A","A/A","T/T","A/T"),
  nrow = 3L, byrow = TRUE,
  dimnames = list(c("S1","S2","S3"), c("M1","M2","M3","M4"))
)
.bgm_map <- data.frame(Chrom = c("1","1","2"), Marker = c("M1","M2","M3"),
                       Pos = c(1e6, 2e6, 3e6), stringsAsFactors = FALSE)

test_that("build_geno_meta drops genotype markers absent from the map", {
  # expect_warning() aborts the call at the warning, so assert and capture
  # in separate calls.
  expect_warning(suppressMessages(build_geno_meta(.bgm_geno, .bgm_map)),
                 "dropped")
  gm <- suppressMessages(suppressWarnings(build_geno_meta(.bgm_geno, .bgm_map)))
  expect_true(is_geno_meta(gm))
  expect_equal(gm$metadata$n_markers, 3L)            # M4 filtered out
  expect_false("M4" %in% dimnames(gm$geno_arrays)[[2L]])
  expect_setequal(dimnames(gm$geno_arrays)[[2L]], c("M1","M2","M3"))
  expect_equal(dim(gm$geno_arrays)[1L], 3L)          # samples untouched
})

test_that("build_geno_meta excludes map markers absent from the genotypes", {
  geno2 <- .bgm_geno[, c("M1","M2"), drop = FALSE]
  expect_warning(suppressMessages(build_geno_meta(geno2, .bgm_map)),
                 "absent from array")
  gm <- suppressMessages(suppressWarnings(build_geno_meta(geno2, .bgm_map)))
  expect_true(is_geno_meta(gm))
  expect_equal(gm$metadata$n_markers, 2L)
  expect_setequal(dimnames(gm$geno_arrays)[[2L]], c("M1","M2"))
})

test_that("build_geno_meta errors when no markers overlap the map", {
  geno_x <- .bgm_geno[, 1:3, drop = FALSE]
  colnames(geno_x) <- c("X1","X2","X3")
  expect_error(suppressWarnings(build_geno_meta(geno_x, .bgm_map)),
               "no markers in common")
})

test_that("build_geno_meta leaves a fully matched marker set unchanged", {
  geno_ok <- .bgm_geno[, c("M1","M2","M3"), drop = FALSE]
  gm <- suppressMessages(build_geno_meta(geno_ok, .bgm_map))
  expect_true(is_geno_meta(gm))
  expect_equal(gm$metadata$n_markers, 3L)
  expect_setequal(dimnames(gm$geno_arrays)[[2L]], c("M1","M2","M3"))
})

test_that("build_geno_meta stores a positional fam as FAM sample_info", {
  geno_ok <- .bgm_geno[, c("M1","M2","M3"), drop = FALSE]
  fam <- data.frame(fid = c("A","A","B"), iid = c("S1","S2","S3"),
                    pa = c("0","0","S1"), ma = c("0","0","S2"),
                    sx = c("1","2","1"), stringsAsFactors = FALSE)
  gm <- suppressMessages(build_geno_meta(geno_ok, .bgm_map, fam = fam))
  expect_true(is_geno_meta(gm))
  expect_equal(names(gm$sample_info),
               c("Family", "Sample", "Parent1", "Parent2", "Sex"))
  expect_equal(gm$sample_info$Sample, dimnames(gm$geno_arrays)[[1L]])
  expect_setequal(gm$sample_info$Family, c("A","B"))
})

test_that("read_genotype delegates marker intersection to build_geno_meta", {
  geno_tmp <- tempfile(fileext = ".txt")
  map_tmp  <- tempfile(fileext = ".txt")
  on.exit(unlink(c(geno_tmp, map_tmp)), add = TRUE)

  # Genotype: markers as rows (M1..M4), samples as header
  writeLines(c("ID\tS1\tS2\tS3",
               "M1\t1/1\t1/4\t4/4",
               "M2\t4/4\t1/4\t1/1",
               "M3\t1/1\t1/1\t4/4",
               "M4\t1/4\t4/4\t1/1"), geno_tmp)
  # Map covers only M1..M3 -> M4 must be dropped downstream
  writeLines(c("Marker\tChrom\tPos",
               "M1\t1\t1000000",
               "M2\t1\t2000000",
               "M3\t2\t3000000"), map_tmp)

  # ST preset: markers as rows, slash-separated SelectionTools cells,
  # map column order Marker | Chrom | Pos. Gated since 0.9.0 (C-5).
  gm <- suppressMessages(suppressWarnings(
    read_genotype(geno_tmp, map_file = map_tmp, preset = "ST", st_formats = TRUE)))
  expect_true(is_geno_meta(gm))
  expect_setequal(dimnames(gm$geno_arrays)[[2L]], c("M1","M2","M3"))
  expect_equal(dim(gm$geno_arrays)[1L], 3L)          # 3 samples
})

# ── resolve_geno_format (preset / parameter-set resolver) ─────────────────────

test_that("resolve_geno_format returns the default parameter set for NULL", {
  spec <- resolve_geno_format(NULL)
  expect_equal(spec$geno_orient, "sample")
  expect_equal(spec$fam_cols, 0L)
  expect_equal(spec$map_order, c("Chrom", "Marker", "Pos"))
  expect_equal(spec$encoding, "nucleotide")
})

test_that("resolve_geno_format resolves built-in presets and aliases", {
  st <- resolve_geno_format("ST")
  expect_equal(st$geno_orient, "marker")
  expect_equal(st$encoding, "selectiontools")
  expect_equal(st$map_order, c("Marker", "Chrom", "Pos"))

  pl <- resolve_geno_format("plink_ped")          # alias for "PL"
  expect_equal(pl$fam_cols, 6L)
  expect_equal(pl$geno_orient, "sample")
})

test_that("resolve_geno_format accepts a custom parameter-set list", {
  spec <- resolve_geno_format(list(geno_orient = "marker",
                                   split_char  = "/",
                                   encoding    = "nucleotide"))
  expect_equal(spec$geno_orient, "marker")
  expect_equal(spec$split_char, "/")
  expect_equal(spec$encoding, "nucleotide")
})

test_that("resolve_geno_format errors on unknown preset and bad fields", {
  expect_error(resolve_geno_format("NOPE"), "unknown preset")
  expect_error(resolve_geno_format(list(encoding = "number")),
               "selectiontools")
})

test_that("read_genotype handles a PL preset with an embedded FAM block", {
  geno_tmp <- tempfile(fileext = ".txt")
  map_tmp  <- tempfile(fileext = ".txt")
  on.exit(unlink(c(geno_tmp, map_tmp)), add = TRUE)

  # Sample-as-row table; first six columns are FAM fields, markers from col 7.
  writeLines(c("FID\tIID\tPAT\tMAT\tSEX\tPHENO\tM1\tM2\tM3",
               "F1\tS1\t0\t0\t1\t-9\tAA\tAT\tTT",
               "F1\tS2\tS1\t0\t2\t-9\tAT\tAA\tAT",
               "F2\tS3\t0\t0\t1\t-9\tTT\tAT\tAA"), geno_tmp)
  writeLines(c("Chrom\tMarker\tPos",
               "1\tM1\t1000000",
               "1\tM2\t2000000",
               "2\tM3\t3000000"), map_tmp)

  gm <- suppressMessages(suppressWarnings(
    read_genotype(geno_tmp, map_file = map_tmp, preset = "PL")))
  expect_true(is_geno_meta(gm))
  expect_equal(dim(gm$geno_arrays)[1L], 3L)          # 3 samples
  expect_setequal(dimnames(gm$geno_arrays)[[1L]], c("S1", "S2", "S3"))
  expect_setequal(dimnames(gm$geno_arrays)[[2L]], c("M1", "M2", "M3"))
  expect_equal(names(gm$sample_info),
               c("Family", "Sample", "Parent1", "Parent2", "Sex"))
  expect_setequal(gm$sample_info$Family, c("F1", "F2"))
})

test_that("write_geno_meta round-trips through a preset (ST)", {
  geno_tmp <- tempfile(fileext = ".txt")
  map_tmp  <- tempfile(fileext = ".txt")
  out_tmp  <- tempfile(fileext = ".csv")
  on.exit(unlink(c(geno_tmp, map_tmp, out_tmp)), add = TRUE)

  writeLines(c("ID\tS1\tS2\tS3",
               "M1\t1/1\t1/4\t4/4",
               "M2\t4/4\t1/4\t1/1",
               "M3\t1/1\t1/1\t4/4"), geno_tmp)
  writeLines(c("Marker\tChrom\tPos",
               "M1\t1\t1000000",
               "M2\t1\t2000000",
               "M3\t2\t3000000"), map_tmp)

  gm <- suppressMessages(suppressWarnings(
    read_genotype(geno_tmp, map_file = map_tmp, preset = "ST", st_formats = TRUE)))
  res <- suppressMessages(write_geno_meta(gm, out_tmp, preset = "ST", st_formats = TRUE))
  expect_true(file.exists(out_tmp))
  # ST output is marker-as-row: 3 markers + 1 header row
  back <- read.csv(out_tmp, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(nrow(back), 3L)
  expect_equal(colnames(back)[1L], "Marker")
  expect_setequal(back$Marker, c("M1", "M2", "M3"))
})
