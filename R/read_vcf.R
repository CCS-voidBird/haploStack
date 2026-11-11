# =============================================================================
# read_vcf.R — VCF genotype files (through vcfR, a Suggests)
#
# read_vcf_genotype(): GT fields to per-copy alleles. The read defaults (QC,
# imputation, phasing) are revisited with the 0.9.1 ingest facade (D5 §3.8).
# Split out of format_helper.R in 0.9.0 (D2-9).
# =============================================================================

#' Read a VCF file into a GenoMeta object
#'
#' Parses a standard VCF (Variant Call Format) file into the unified
#' \code{GenoMeta} structure using the \pkg{vcfR} package. Genotype calls in
#' the \code{GT} FORMAT field (phased or unphased) are converted to dosage
#' (0/1/2) and split into two haplotype slices
#' (2 -> \code{"A"}/\code{"A"}, 1 -> \code{"A"}/\code{"B"},
#' 0 -> \code{"B"}/\code{"B"}).
#'
#' @section VCF requirements:
#' \itemize{
#'   \item VCF version 4.0+ with a \code{GT} FORMAT field.
#'   \item Diploid samples only (exactly two allele calls per GT).
#'   \item Multi-allelic sites are filtered to biallelic only when
#'     \code{biallelic_only = TRUE} (default).
#' }
#'
#' @param file               Path to the VCF file (plain or bgzipped).
#' @param block_file         Optional block map: a path or a data.frame in one
#'   of the two shapes \code{\link{attach_block_info}} accepts (block-marker
#'   map \code{Chrom, Marker, Block} or block-region map \code{Chrom, Block,
#'   Start, End}; header required). Attached after the read.
#' @param block_overlap      Overlap policy for the block map, as
#'   \code{attach_block_info(overlap =)}: \code{"error"} (default) or
#'   \code{"first"}.
#' @param pos_unit           Position unit: \code{"bp"} (default) or
#'   \code{"mbp"}.
#' @param call_rate_threshold Numeric (0-1). Markers below this call rate are
#'   removed. Default \code{0.8}.
#' @param maf_threshold      Numeric (0-1). Markers below this minor allele
#'   frequency are removed. Default \code{0.01}.
#' @param biallelic_only     Logical. Drop multi-allelic sites. Default
#'   \code{TRUE}.
#' @param missing_chars      Extra missing-allele characters for this read
#'   (see \code{\link{read_genotype}}); only a character registered for the
#'   \code{"nucleotide"} vocabulary of the converted container applies.
#'
#' @return A valid \code{GenoMeta} list with \code{metadata$source_format =
#'   "vcf"} and A/B haplotype encoding.
#'
#' @seealso \code{\link{read_dart_genotype}}, \code{\link{read_plink_ped}},
#'   \code{\link{filter_geno_meta}}, \code{\link{compute_grm}}
#'
#' @export
read_vcf_genotype <- function(file,
                               block_file          = NULL,
                               pos_unit            = "bp",
                               call_rate_threshold = 0.8,
                               maf_threshold       = 0.01,
                               biallelic_only      = TRUE,
                               block_overlap       = c("error", "first"),
                               missing_chars       = NULL) {

  if (!file.exists(file))
    .hs_abort("read_vcf_genotype: file not found: ", file)
  if (!requireNamespace("vcfR", quietly = TRUE))
    .hs_abort(class = "hs_error_dependency",
         data = list(hint = "install.packages(\"vcfR\")"),
         "read_vcf_genotype requires the 'vcfR' package. ",
         "Install with: install.packages('vcfR')")

  vcf <- vcfR::read.vcfR(file, verbose = FALSE)

  # ---- Extract fixed fields (CHROM, POS, ID, ALT) --------------------------
  fix_df <- as.data.frame(vcfR::getFIX(vcf), stringsAsFactors = FALSE)
  colnames(fix_df) <- toupper(colnames(fix_df))

  if (!all(c("CHROM","POS") %in% names(fix_df)))
    .hs_abort("read_vcf_genotype: VCF is missing CHROM or POS columns.")

  # ---- Biallelic filter (before GT extraction) -----------------------------
  if (biallelic_only) {
    alt_col      <- fix_df[["ALT"]]
    is_biallelic <- !grepl(",", alt_col, fixed = TRUE) &
                    !is.na(alt_col) & alt_col != "."
    n_multi <- sum(!is_biallelic)
    if (n_multi > 0L)
      .hs_inform("read_vcf_genotype: removed ", n_multi, " multi-allelic site(s).")
    vcf    <- vcf[is_biallelic, ]
    fix_df <- fix_df[is_biallelic, , drop = FALSE]
  }

  if (nrow(fix_df) == 0L)
    .hs_abort("read_vcf_genotype: no sites remaining after biallelic filter.")

  # ---- Build marker IDs ----------------------------------------------------
  raw_ids    <- fix_df[["ID"]]
  missing_id <- is.na(raw_ids) | raw_ids == "." | raw_ids == ""
  raw_ids[missing_id] <- paste0(fix_df[["CHROM"]][missing_id], ":",
                                 fix_df[["POS"]][missing_id])
  marker_ids <- make.unique(raw_ids, sep = "_dup")

  # ---- Clean chromosome names ----------------------------------------------
  chrom_raw   <- fix_df[["CHROM"]]
  chrom_clean <- gsub("^[Cc]hr", "", chrom_raw)   # strip chr/Chr prefix
  chrom_clean <- gsub("H$",      "", chrom_clean)  # strip barley H-suffix

  pos_num <- suppressWarnings(as.numeric(fix_df[["POS"]]))
  if (tolower(pos_unit) == "mbp") pos_num <- pos_num * 1e6

  # ---- Extract GT matrix ---------------------------------------------------
  gt_mat <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)
  if (is.null(gt_mat) || nrow(gt_mat) == 0L)
    .hs_abort("read_vcf_genotype: no GT data found in VCF.")

  rownames(gt_mat) <- marker_ids

  # ---- Parse GT strings to dosage (count of alt alleles) ------------------
  .gt_to_dosage <- function(gt_vec) {
    gt_clean <- gsub("|", "/", gt_vec, fixed = TRUE)  # phased to unphased
    parts    <- strsplit(gt_clean, "/", fixed = TRUE)
    vapply(parts, function(x) {
      if (length(x) != 2L || any(x == ".")) return(NA_integer_)
      a1 <- suppressWarnings(as.integer(x[1L]))
      a2 <- suppressWarnings(as.integer(x[2L]))
      if (is.na(a1) || is.na(a2)) return(NA_integer_)
      pmin(a1 + a2, 2L)
    }, integer(1L))
  }

  dos_mat <- matrix(
    .gt_to_dosage(as.vector(gt_mat)),
    nrow     = nrow(gt_mat),
    ncol     = ncol(gt_mat),
    dimnames = list(marker_ids, colnames(gt_mat))
  )

  n_mk  <- nrow(dos_mat)
  n_smp <- ncol(dos_mat)
  if (n_mk  == 0L) .hs_abort("read_vcf_genotype: no markers after GT parsing.")
  if (n_smp == 0L) .hs_abort("read_vcf_genotype: no samples found in VCF.")

  sample_names <- colnames(dos_mat)

  # ---- Call rate filter ----------------------------------------------------
  call_rate <- rowMeans(!is.na(dos_mat))
  keep_cr   <- call_rate >= call_rate_threshold
  n_fail_cr <- sum(!keep_cr)
  if (n_fail_cr > 0L)
    .hs_inform("read_vcf_genotype: removed ", n_fail_cr,
            " marker(s) below call rate ", call_rate_threshold, ".")
  dos_mat     <- dos_mat[keep_cr, , drop = FALSE]
  chrom_clean <- chrom_clean[keep_cr]
  pos_num     <- pos_num[keep_cr]

  # ---- MAF filter ----------------------------------------------------------
  p          <- rowMeans(dos_mat, na.rm = TRUE) / 2
  maf        <- pmin(p, 1 - p)
  keep_maf   <- !is.na(maf) & maf >= maf_threshold
  n_fail_maf <- sum(!keep_maf)
  if (n_fail_maf > 0L)
    .hs_inform("read_vcf_genotype: removed ", n_fail_maf,
            " marker(s) below MAF ", maf_threshold, ".")
  dos_mat     <- dos_mat[keep_maf, , drop = FALSE]
  chrom_clean <- chrom_clean[keep_maf]
  pos_num     <- pos_num[keep_maf]

  n_mk_final <- nrow(dos_mat)
  if (n_mk_final == 0L)
    .hs_abort("read_vcf_genotype: no markers remaining after call rate / MAF filtering.")

  marker_names <- rownames(dos_mat)

  # ---- Mean-impute NA genotypes --------------------------------------------
  na_rows <- which(rowSums(is.na(dos_mat)) > 0L)
  if (length(na_rows) > 0L) {
    .hs_inform("read_vcf_genotype: mean-imputing NA genotypes in ",
            length(na_rows), " marker(s).")
    row_means <- rowMeans(dos_mat, na.rm = TRUE)
    for (r in na_rows) {
      na_cols <- which(is.na(dos_mat[r, ]))
      dos_mat[r, na_cols] <- as.integer(round(row_means[r]))
    }
  }

  # ---- Split dosage to A/B haplotype slices --------------------------------
  arr <- array(
    NA_character_,
    dim      = c(n_smp, n_mk_final, 2L),
    dimnames = list(sample_names, marker_names, c("H1", "H2"))
  )
  arr[, , 1L] <- t(ifelse(dos_mat == 0L, "B", "A"))  # hap1: "B" only if dosage=0
  arr[, , 2L] <- t(ifelse(dos_mat >= 2L, "A", "B"))  # hap2: "A" only if dosage=2

  # ---- Map info & block files ----------------------------------------------
  map_df <- data.frame(
    Chrom  = chrom_clean,
    Marker = marker_names,
    Pos    = pos_num,
    stringsAsFactors = FALSE
  )

  # ---- Assemble GenoMeta ---------------------------------------------------
  meta <- new_geno_meta(
    geno_arrays           = arr,
    map_info              = map_df,
    metadata = list(
      ploidy        = 2L,
      char_dict     = list("A" = 1L, "B" = 2L),
      n_samples     = n_smp,
      n_markers     = n_mk_final,
      source_format = "vcf",
      encoding      = "nucleotide"
    ),
    # As for DArT: GT fields are converted to "A"/"B" alleles before this call,
    # so the container carries an allele vocabulary (FR3).
    source_format = "nucleotide",
    missing_chars = missing_chars
  )

  # Blocks: build, then attach (one block path, 0.9.0)
  if (!is.null(block_file))
    meta <- attach_block_info(meta, block_file, overlap = block_overlap)
  meta
}
