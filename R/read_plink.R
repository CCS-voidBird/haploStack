# =============================================================================
# read_plink.R — PLINK 1.9 PED + MAP (the built-in load preset)
#
# read_plink_ped(): the FAM block (six columns) is kept as sample_info, the
# alleles are read as text, PLINK's "0" is a missing call, the map is the
# 4-column .map; a block file can be attached on the way in. Split out of
# format_helper.R in 0.9.0 (D2-9).
# =============================================================================

#' Read a Plink 1.9 PED + MAP file pair into a GenoMeta object
#'
#' Parses a \code{.ped} file (6-column header + allele pairs per marker) and
#' its companion \code{.map} file into the unified \code{GenoMeta} structure.
#'
#' @param ped_file  Path to the \code{.ped} file.
#' @param map_file  Path to the \code{.map} file (4-column PLINK format:
#'   chr, snp_id, cM, bp).
#' @param block_file Optional block map: a path or a data.frame in one of the
#'   two shapes \code{\link{attach_block_info}} accepts (block-marker map
#'   \code{Chrom, Marker, Block} or block-region map \code{Chrom, Block, Start,
#'   End}; header required). Attached after the read, so the FAM and the block
#'   provenance are kept.
#' @param block_overlap Overlap policy for the block map, as
#'   \code{attach_block_info(overlap =)}: \code{"error"} (default) or
#'   \code{"first"}.
#' @param pos_unit  Position unit in \code{map_file}: \code{"bp"} (default)
#'   or \code{"mbp"}.
#' @param missing_chars Extra missing-allele characters for this read:
#'   \code{NULL} (default) takes the session registry
#'   (\code{\link{use_missing_chars}}); a character vector names them directly
#'   (\code{character(0)} ignores the registry). See
#'   \code{\link{new_geno_meta}}.
#'
#' @return A valid \code{GenoMeta} list.
#'
#' @export
read_plink_ped <- function(ped_file, map_file, block_file = NULL,
                            pos_unit = "bp", missing_chars = NULL,
                            block_overlap = c("error", "first")) {
  ped_file <- as.character(ped_file)
  map_file <- as.character(map_file)

  # Parse MAP file (chr, snp_id, cM, bp) — use process_genotype map path.
  # Read as text and coerce Pos here: a MAP column is an identifier until this
  # function says otherwise (see .read_table_chr).
  raw_map <- .read_table_chr(map_file, header = FALSE, strict = TRUE)
  if (ncol(raw_map) >= 4L) {
    map_df <- data.frame(
      Chrom  = as.character(raw_map[[1L]]),
      Marker = as.character(raw_map[[2L]]),
      Pos    = as.numeric(raw_map[[4L]]),
      stringsAsFactors = FALSE
    )
  } else if (ncol(raw_map) == 3L) {
    map_df <- data.frame(
      Chrom  = as.character(raw_map[[1L]]),
      Marker = as.character(raw_map[[2L]]),
      Pos    = as.numeric(raw_map[[3L]]),
      stringsAsFactors = FALSE
    )
  } else {
    .hs_abort("read_plink_ped: MAP file must have 3 or 4 columns.")
  }

  # Q5: delegate Mbp warning, the preset of map is "PL"
  map_df <- .process_map(map_df, marker_names = map_df$Marker,
                         map_preset = "PL", pos_unit = pos_unit)

  n_markers <- nrow(map_df)

  # Parse PED file. Read as text: an allele column must never be type-inferred
  # (a monomorphic-T column would come back logical — see .read_table_chr).
  ped <- .read_table_chr(ped_file, header = FALSE, strict = TRUE)
  # Columns 1-6: FID, IID, father, mother, sex, phenotype; then ploidy*n_markers allele cols
  n_allele_cols <- ncol(ped) - 6L
  if (n_allele_cols %% n_markers != 0L)
    .hs_abort("read_plink_ped: n_allele_cols (", n_allele_cols,
         ") is not divisible by n_markers (", n_markers,
         "). Cannot determine ploidy.")
  ploidy_detected <- as.integer(n_allele_cols %/% n_markers)

  sample_names <- as.character(ped[[2L]])
  n_samples    <- length(sample_names)

  allele_mat <- as.matrix(ped[, 7L:ncol(ped)])
  mode(allele_mat) <- "character"
  rownames(allele_mat) <- sample_names

  # Collapse ploidy consecutive allele cols → (n_samples × n_markers)
  geno_concat <- .reshape_wide_alleles(allele_mat, n_markers, ploidy_detected)
  colnames(geno_concat) <- map_df$Marker

  split_res <- split_haplotypes(geno_concat,
                                 split_char   = "",
                                 orientation  = 1L,
                                 ploidy       = ploidy_detected,
                                 sample_names = sample_names,
                                 marker_names = map_df$Marker)
  geno_arr <- split_res$array   # (n_samples × n_markers × ploidy_detected)

  # Simple char_dict from observed alleles
  obs <- unique(c(geno_arr))
  obs <- obs[!is.na(obs) & obs != "0"]
  char_dict <- setNames(seq_along(obs), obs)

  # Declare the format rather than letting it be detected (FR3). A PED holds
  # alleles — letters, or PLINK's 1/2 (or 1/2/3/4) codes — and never dosages,
  # and "0" is a missing call. This matters: a 1/2-coded PED splits into
  # single digits, which detection would read as "numeric", and under the
  # dosage rule "0" would then be kept as a real genotype. Declaring the
  # format closes that route. Since 0.9.0 (C-5) the one categorical vocabulary
  # is "nucleotide" — categorical allele labels — and PLINK's digit codes are
  # kept as the labels they are: haploStack does not translate them (the
  # SelectionTools table A = 1 .. T = 4 is not PLINK's 1/2 coding).
  ped_encoding <- "nucleotide"

  # PED columns 1-5 are the FAM fields (FID, IID, father, mother, sex) — feed
  # them straight into sample_info (positional FAM order).
  fam_si <- data.frame(
    Family  = as.character(ped[[1L]]),
    Sample  = as.character(ped[[2L]]),
    Parent1 = as.character(ped[[3L]]),
    Parent2 = as.character(ped[[4L]]),
    Sex     = as.character(ped[[5L]]),
    stringsAsFactors = FALSE
  )

  meta <- new_geno_meta(
    geno_arrays           = geno_arr,
    map_info              = map_df,
    metadata = list(
      ploidy        = ploidy_detected,
      char_dict     = char_dict,
      n_samples     = n_samples,
      n_markers     = n_markers,
      source_format = "plink_ped",
      encoding      = ped_encoding
    ),
    sample_info   = fam_si,
    source_format = ped_encoding,
    # Belt and braces: "0" is already missing under both allele vocabularies,
    # and declaring it here says the reader knows its own sentinel rather than
    # inheriting one.
    extra_na      = "0",
    missing_chars = missing_chars
  )

  # Blocks: build, then attach (one block path, 0.9.0)
  if (!is.null(block_file))
    meta <- attach_block_info(meta, block_file, overlap = block_overlap)
  meta
}
