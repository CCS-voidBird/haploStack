# =============================================================================
# read_dart.R — DArT 1-row / 2-row genotype files
#
# read_dart_genotype(): the DArT layouts as haploStack 0.9.0 reads them. The
# dosage model (which allele a count counts, imputation, the phase flag) is
# 0.9.1's P1 (dev/PRD_dosage_model.md); until then this reader keeps its
# 0.8.8 behaviour and the help says what it does not do. Split out of
# format_helper.R in 0.9.0 (D2-9).
# =============================================================================

#' Read a DArTseq SNP report into a GenoMeta object
#'
#' Parses a DArTseq 1-row SNP report (wide format: one row per marker, one
#' column per sample, dosage values 0/1/2) into the unified \code{GenoMeta}
#' structure. Marker metadata columns are detected automatically, and the
#' dosage is split into two haplotype slices for consistency with the rest of
#' the haploStack pipeline (2 → ["A","A"], 1 → ["A","B"], 0 → ["B","B"]).
#'
#' @section DArT report format:
#' The file must be a CSV or tab-delimited text file. Multiple leading metadata
#' rows (rows before the header line) are skipped with \code{n_meta_rows};
#' they are not detected. Sample columns are identified as the first
#' contiguous block of columns containing only numeric values (0, 1, 2, or NA).
#'
#' @param file               Path to the DArTseq SNP report file.
#' @param n_meta_rows        Integer number of leading rows to skip before the
#'   header line (the row holding the column names), or \code{NULL} (default)
#'   when the first line is the header. There is no auto-detection: a report
#'   with leading metadata rows needs the count.
#' @param marker_col         Column name used as the marker ID. Searched in
#'   order: \code{"AlleleID"}, \code{"TargetID"}, \code{"CloneID"},
#'   \code{"MarkerName"}, \code{"Marker"}. Or supply a column name directly.
#'   Default \code{NULL} = auto.
#' @param chrom_col          Column name for chromosome. Searched in order:
#'   first column matching \code{"[Cc]hrom"} but not \code{"[Pp]os"}.
#'   Default \code{NULL} = auto.
#' @param pos_col            Column name for SNP position. Searched in order:
#'   first column matching \code{"[Pp]os"}, \code{"ChromPos"}, \code{"SnpPos"}.
#'   Default \code{NULL} = auto.
#' @param sample_start_col   Integer index of the first sample column. When
#'   \code{NULL} (default), auto-detected as the first column after all
#'   character-type marker-metadata columns.
#' @param block_file         Optional block map: a path or a data.frame in one
#'   of the two shapes \code{\link{attach_block_info}} accepts (block-marker
#'   map \code{Chrom, Marker, Block} or block-region map \code{Chrom, Block,
#'   Start, End}; header required). Attached after the read.
#' @param block_overlap      Overlap policy for the block map, as
#'   \code{attach_block_info(overlap =)}: \code{"error"} (default) or
#'   \code{"first"}.
#' @param pos_unit           Position unit: \code{"bp"} (default) or
#'   \code{"mbp"}.
#' @param call_rate_threshold Numeric (0–1). Markers with a call rate below
#'   this threshold are removed. Default \code{0.8}.
#' @param maf_threshold      Numeric (0–1). Markers with minor allele
#'   frequency below this threshold are removed. Default \code{0.01}.
#' @param missing_chars      Extra missing-allele characters for this read
#'   (see \code{\link{read_genotype}}). A DArT report is dosage on disk, so
#'   only a character registered for the \code{"nucleotide"} vocabulary of the
#'   converted container applies.
#'
#' @return A valid \code{GenoMeta} list with:
#' \describe{
#'   \item{\code{geno_arrays}}{3-D character array \code{[n_samples × n_markers × 2]}.
#'     Allele encoding: \code{"A"} = reference, \code{"B"} = alternative.
#'     \code{NA} = missing call.}
#'   \item{\code{map_info}}{data.frame \code{(Chr, Marker, Pos)}.}
#'   \item{\code{metadata$source_format}}{\code{"dart_seq"}.}
#' }
#'
#' @seealso \code{\link{read_st_genotype}}, \code{\link{read_plink_ped}},
#'   \code{\link{filter_geno_meta}}, \code{\link{compute_grm}}
#'
#' @export
read_dart_genotype <- function(file,
                                n_meta_rows         = NULL,
                                marker_col          = NULL,
                                chrom_col           = NULL,
                                pos_col             = NULL,
                                sample_start_col    = NULL,
                                block_file          = NULL,
                                block_overlap       = c("error", "first"),
                                pos_unit            = "bp",
                                call_rate_threshold = 0.8,
                                maf_threshold       = 0.01,
                                missing_chars       = NULL) {

  file <- as.character(file)
  if (!file.exists(file))
    .hs_abort("read_dart_genotype: file not found: ", file)

  # ---- Read raw file -------------------------------------------------------
  # Read as text. Every value this function needs as a number is re-derived
  # below with as.numeric(as.character(.)) — dosages, positions, the
  # metadata/sample column heuristic — so nothing depends on inference, and a
  # T/F flag column can no longer arrive as logical (see .read_table_chr).
  sep <- if (grepl("[.]csv$", file, ignore.case = TRUE)) "," else "\t"
  raw <- .read_table_chr(file, header = TRUE, sep = sep)

  if (nrow(raw) == 0L || ncol(raw) < 2L)
    .hs_abort("read_dart_genotype: file appears empty or has fewer than 2 columns.")

  # ---- Handle leading metadata rows ----------------------------------------
  # DArT reports sometimes have extra rows at the top (callrate rows, etc.)
  # before the actual per-marker data. Auto-detect by finding the first row
  # where the majority of sample-candidate columns are numeric 0/1/2.
  if (!is.null(n_meta_rows) && n_meta_rows > 0L) {
    colnames(raw) <- as.character(unlist(raw[n_meta_rows, ]))
    raw <- raw[(n_meta_rows + 1L):nrow(raw), , drop = FALSE]
    rownames(raw) <- NULL
  }

  # ---- Identify marker metadata vs sample columns --------------------------
  # Heuristic: marker metadata columns are character; sample columns are
  # coercible to 0/1/2/NA. Walk from right until we find a non-numeric column.
  is_numeric_col <- vapply(seq_len(ncol(raw)), function(j) {
    v <- suppressWarnings(as.numeric(as.character(raw[[j]])))
    non_na <- v[!is.na(v)]
    length(non_na) > 0L && all(non_na %in% c(0, 1, 2))
  }, logical(1L))

  if (!is.null(sample_start_col)) {
    first_sample <- as.integer(sample_start_col)
  } else {
    # Find the first column that looks like a sample column; require that all
    # subsequent columns also pass the numeric test.
    first_sample <- NULL
    for (j in seq_len(ncol(raw))) {
      if (all(is_numeric_col[j:ncol(raw)])) {
        first_sample <- j
        break
      }
    }
    if (is.null(first_sample))
      .hs_abort("read_dart_genotype: could not auto-detect sample columns. ",
           "Please supply 'sample_start_col'.")
  }

  meta_cols_df <- raw[, seq_len(first_sample - 1L), drop = FALSE]
  geno_df      <- raw[, first_sample:ncol(raw), drop = FALSE]
  sample_names <- colnames(geno_df)

  # ---- Identify marker ID column -------------------------------------------
  candidate_marker_cols <- c("AlleleID","TargetID","CloneID","MarkerName","Marker")
  if (!is.null(marker_col)) {
    mk_col <- marker_col
  } else {
    mk_col <- NULL
    for (cn in candidate_marker_cols) {
      if (cn %in% names(meta_cols_df)) { mk_col <- cn; break }
    }
    if (is.null(mk_col)) mk_col <- names(meta_cols_df)[1L]
  }
  marker_names <- as.character(meta_cols_df[[mk_col]])

  # Deduplicate marker names (rare but possible)
  if (anyDuplicated(marker_names)) {
    marker_names <- make.unique(marker_names, sep = "_")
    .hs_warn("read_dart_genotype: duplicate marker IDs found; made unique.")
  }

  # ---- Identify chromosome and position columns ----------------------------
  all_meta_cols <- names(meta_cols_df)

  if (!is.null(chrom_col)) {
    chr_col_name <- chrom_col
  } else {
    chr_col_name <- NULL
    for (cn in all_meta_cols) {
      if (grepl("[Cc]hrom", cn) && !grepl("[Pp]os", cn)) {
        chr_col_name <- cn; break
      }
    }
  }

  if (!is.null(pos_col)) {
    pos_col_name <- pos_col
  } else {
    pos_col_name <- NULL
    pos_patterns <- c("ChromPosSnp", "SnpPosition", "ChromPos",
                      "[Pp]os", "[Pp]osition")
    for (pat in pos_patterns) {
      for (cn in all_meta_cols) {
        if (grepl(pat, cn)) { pos_col_name <- cn; break }
      }
      if (!is.null(pos_col_name)) break
    }
  }

  chr_vec <- if (!is.null(chr_col_name) && chr_col_name %in% all_meta_cols) {
    # Normalise H-suffix (1H → 1) common in barley
    gsub("H$", "", as.character(meta_cols_df[[chr_col_name]]))
  } else {
    rep(NA_character_, nrow(meta_cols_df))
  }

  pos_vec <- if (!is.null(pos_col_name) && pos_col_name %in% all_meta_cols) {
    suppressWarnings(as.numeric(meta_cols_df[[pos_col_name]]))
  } else {
    rep(NA_real_, nrow(meta_cols_df))
  }

  # ---- Build numeric dosage matrix -----------------------------------------
  n_all_markers <- length(marker_names)
  n_samples     <- length(sample_names)

  dosage <- matrix(NA_real_, nrow = n_all_markers, ncol = n_samples,
                   dimnames = list(marker_names, sample_names))
  for (j in seq_len(n_samples)) {
    v <- suppressWarnings(as.numeric(as.character(geno_df[[j]])))
    # DArT missing: -1, blank strings, or literal NA
    v[!is.na(v) & v < 0L] <- NA_real_
    dosage[, j] <- v
  }

  # ---- Filter by call rate --------------------------------------------------
  call_rate <- rowMeans(!is.na(dosage))
  pass_cr   <- call_rate >= call_rate_threshold
  if (sum(pass_cr) < n_all_markers)
    .hs_inform("read_dart_genotype: ", sum(!pass_cr),
            " marker(s) removed (call rate < ", call_rate_threshold, ").")

  dosage       <- dosage[pass_cr, , drop = FALSE]
  chr_vec      <- chr_vec[pass_cr]
  pos_vec      <- pos_vec[pass_cr]
  marker_names <- rownames(dosage)
  n_markers    <- nrow(dosage)

  # ---- Filter by MAF -------------------------------------------------------
  freq <- rowMeans(dosage, na.rm = TRUE) / 2   # allele freq of alt allele
  maf  <- pmin(freq, 1 - freq)
  pass_maf <- !is.na(maf) & maf >= maf_threshold
  if (sum(pass_maf) < n_markers)
    .hs_inform("read_dart_genotype: ", sum(!pass_maf),
            " marker(s) removed (MAF < ", maf_threshold, ").")

  dosage       <- dosage[pass_maf, , drop = FALSE]
  chr_vec      <- chr_vec[pass_maf]
  pos_vec      <- pos_vec[pass_maf]
  marker_names <- rownames(dosage)
  n_markers    <- nrow(dosage)

  if (n_markers == 0L)
    .hs_abort("read_dart_genotype: no markers remain after filtering.")

  # ---- Split dosage to two haplotype slices --------------------------------
  # 2 (AA) → hap1 = "A", hap2 = "A"
  # 1 (AB) → hap1 = "A", hap2 = "B"
  # 0 (BB) → hap1 = "B", hap2 = "B"
  # NA     → hap1 = NA,  hap2 = NA
  hap1 <- matrix(ifelse(dosage == 0L, "B", "A"),
                 nrow = n_markers, ncol = n_samples)
  hap2 <- matrix(ifelse(dosage >= 2L, "A", "B"),
                 nrow = n_markers, ncol = n_samples)
  hap1[is.na(dosage)] <- NA_character_
  hap2[is.na(dosage)] <- NA_character_

  # [n_samples × n_markers × ploidy=2]; hap1/hap2 are [n_markers × n_samples]
  geno_arr <- array(NA_character_,
                    dim      = c(n_samples, n_markers, 2L),
                    dimnames = list(sample_names, marker_names, c("1","2")))
  geno_arr[, , 1L] <- t(hap1)
  geno_arr[, , 2L] <- t(hap2)

  # ---- Build map_info -------------------------------------------------------
  map_df <- data.frame(
    Chrom  = chr_vec,
    Marker = marker_names,
    Pos    = pos_vec,
    stringsAsFactors = FALSE
  )
  map_df <- .process_map(map_df, marker_names = marker_names,
                         map_preset = NULL, pos_unit = pos_unit)

  # ---- Block files ----------------------------------------------------------
  # ---- Assemble GenoMeta ---------------------------------------------------
  char_dict <- list("A" = 1L, "B" = 2L)

  meta <- new_geno_meta(
    geno_arrays           = geno_arr,
    map_info              = map_df,
    metadata = list(
      ploidy        = 2L,
      char_dict     = char_dict,
      n_samples     = n_samples,
      n_markers     = n_markers,
      source_format = "dart_seq",
      encoding      = "nucleotide"
    ),
    # A DArT report is dosage on disk, but the reader has already converted it
    # to "A"/"B" alleles by this point, so the container is an allele one and
    # its vocabulary is the allele vocabulary (FR3).
    source_format = "nucleotide",
    missing_chars = missing_chars
  )

  # Blocks: build, then attach (one block path, 0.9.0)
  if (!is.null(block_file))
    meta <- attach_block_info(meta, block_file, overlap = block_overlap)
  meta
}
