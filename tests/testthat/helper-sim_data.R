# =============================================================================
# helper-sim_data.R
#
# Simulation dataset generators for haploStack QA tests.
#
# All functions return self-contained, deterministic test objects that do NOT
# depend on any external files or package state.
#
# Data structure conventions
# --------------------------
# Genotype matrix : rows = markers, cols = samples
#                   values = "A/A", "A/T", "T/T" (nucleotide) or
#                            "1/1", "1/4", "4/4" (SelectionTools)
# Map data frame  : cols Marker, Chrom, Pos
# Block map       : cols Chrom, Marker, Block
# Effect data frame: cols Marker, Allele, Effect (no header)
# =============================================================================

# Attach data.table so its functions (is.data.table, as.data.table, etc.)
# are available directly in all test files without library(data.table) per file.
library(data.table)

#' Generate a minimal simulated genotype matrix
#'
#' @param n_samples  Number of samples (columns). Default 10.
#' @param n_markers  Number of markers (rows). Default 20.
#' @param encoding   "nucleotide" or "selectiontools". Default "selectiontools".
#' @param seed       Random seed. Default 42.
#' @return Character matrix (n_markers × n_samples), named rows and cols.
sim_geno_matrix <- function(n_samples  = 10L,
                             n_markers  = 20L,
                             encoding   = "selectiontools",
                             seed       = 42L) {
  set.seed(seed)
  if (encoding == "selectiontools") {
    alleles <- c("1", "4")
  } else {
    alleles <- c("A", "T")
  }
  mat <- matrix(NA_character_, nrow = n_markers, ncol = n_samples)
  for (i in seq_len(n_markers)) {
    for (j in seq_len(n_samples)) {
      a1 <- sample(alleles, 1L)
      a2 <- sample(alleles, 1L)
      mat[i, j] <- paste0(a1, "/", a2)
    }
  }
  rownames(mat) <- paste0("Mkr", seq_len(n_markers))
  colnames(mat) <- paste0("S", seq_len(n_samples))
  mat
}

#' Generate a marker map data frame
#'
#' @param marker_names  Character vector of marker IDs.
#' @param n_chrom       Number of chromosomes. Default 3.
#' @return data.frame with columns Marker, Chrom, Pos.
sim_map_data <- function(marker_names,
                          n_chrom = 3L) {
  n_markers <- length(marker_names)
  chrom_id  <- sort(rep(seq_len(n_chrom), length.out = n_markers))
  pos_vals  <- unlist(tapply(chrom_id, chrom_id, function(chr_idx) {
    n_in_chr <- sum(chrom_id == chr_idx[1L])
    sort(sample(seq(1e6L, 1e8L, by = 1e5L), n_in_chr, replace = FALSE))
  }))

  data.frame(
    Marker = marker_names,
    Chrom  = as.character(chrom_id),
    Pos    = pos_vals,
    stringsAsFactors = FALSE
  )
}

#' Generate a block map data frame
#'
#' Assigns consecutive markers to blocks within each chromosome.
#'
#' @param map_data     data.frame from \code{sim_map_data()}.
#' @param block_size   Number of markers per block. Default 4.
#' @return data.frame with columns Chrom, Marker, Block.
sim_block_map <- function(map_data, block_size = 4L) {
  result <- lapply(unique(map_data$Chrom), function(chr) {
    chr_markers <- map_data$Marker[map_data$Chrom == chr]
    n <- length(chr_markers)
    blk_idx <- ceiling(seq_len(n) / block_size)
    data.frame(
      Chrom  = chr,
      Marker = chr_markers,
      Block  = paste0("Chr", chr, "_B", blk_idx),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}

#' Generate a marker effect data frame
#'
#' @param marker_names Character vector of marker IDs.
#' @param alleles      Character vector of allele states. Default c("A","T").
#' @param seed         Random seed. Default 42.
#' @return data.frame with columns Marker, Allele, Effect.
sim_effect_df <- function(marker_names,
                           alleles = c("A", "T"),
                           seed    = 42L) {
  set.seed(seed)
  pairs <- expand.grid(Marker = marker_names, Allele = alleles,
                       stringsAsFactors = FALSE)
  pairs$Effect <- stats::rnorm(nrow(pairs), mean = 0, sd = 0.5)
  pairs <- pairs[order(pairs$Marker, pairs$Allele), ]
  rownames(pairs) <- NULL
  pairs
}

#' Build a GenoMeta from simulated genotype data
#'
#' Wraps \code{sim_geno_matrix()}, \code{sim_map_data()}, and
#' \code{sim_block_map()} into a single call, splits the 2-D genotype matrix
#' into the 3-D array format expected by \code{new_geno_meta()}, and
#' optionally attaches a block map via \code{attach_block_info()}.
#'
#' @param n_samples  Number of samples. Default 10.
#' @param n_markers  Number of markers. Default 20.
#' @param encoding   \code{"selectiontools"} (default) or \code{"nucleotide"}.
#' @param n_chrom    Number of chromosomes. Default 2.
#' @param block_size Markers per block. \code{NULL} (default) leaves block
#'   fields empty in the returned GenoMeta.
#' @param seed       Random seed. Default 42.
#' @return A valid \code{GenoMeta} list.
sim_geno_meta <- function(n_samples  = 10L,
                           n_markers  = 20L,
                           encoding   = "selectiontools",
                           n_chrom    = 2L,
                           block_size = NULL,
                           seed       = 42L) {

  geno <- sim_geno_matrix(n_samples, n_markers, encoding, seed)
  map  <- sim_map_data(rownames(geno), n_chrom = n_chrom)
  bmap <- if (!is.null(block_size)) sim_block_map(map, block_size) else NULL

  marker_names <- rownames(geno)
  sample_names <- colnames(geno)
  n_mk         <- length(marker_names)
  n_samp       <- length(sample_names)

  # Split "1/4" → ["1","4"] per cell into a 3-D array [n_samples × n_markers × 2]
  arr <- array(NA_character_, dim = c(n_samp, n_mk, 2L),
               dimnames = list(sample_names, marker_names, c("1", "2")))
  for (j in seq_len(n_samp)) {
    parts <- strsplit(geno[, j], "/", fixed = TRUE)
    arr[j, , 1L] <- vapply(parts, `[`, character(1L), 1L)
    arr[j, , 2L] <- vapply(parts, `[`, character(1L), 2L)
  }

  char_dict <- if (encoding == "selectiontools") c("1" = 1L, "4" = 2L)
               else                              c("A" = 1L, "T" = 2L)

  meta <- new_geno_meta(
    geno_arrays = arr,
    map_info    = map,
    metadata    = list(
      ploidy        = 2L,
      char_dict     = char_dict,
      n_samples     = n_samp,
      n_markers     = n_mk,
      source_format = encoding
    )
  )

  if (!is.null(bmap)) meta <- attach_block_info(meta, bmap)
  meta
}


#' Generate a complete simulation dataset
#'
#' Returns a named list with all components needed to run haploStack tests:
#' geno, map, block_map, effect_df_nuc (nucleotide), effect_df_st
#' (SelectionTools numeric), geno_meta.
#'
#' @param n_samples  Number of samples. Default 10.
#' @param n_markers  Number of markers. Default 20.
#' @param block_size Markers per block. Default 4.
#' @param seed       Random seed. Default 42.
#' @return Named list.
sim_dataset <- function(n_samples  = 10L,
                         n_markers  = 20L,
                         block_size = 4L,
                         seed       = 42L) {

  geno_st  <- sim_geno_matrix(n_samples, n_markers, "selectiontools", seed)
  geno_nuc <- sim_geno_matrix(n_samples, n_markers, "nucleotide",     seed + 1L)

  map      <- sim_map_data(rownames(geno_st), n_chrom = 3L)
  bmap     <- sim_block_map(map, block_size)

  eff_nuc  <- sim_effect_df(rownames(geno_st), alleles = c("A","T"), seed)
  eff_st   <- sim_effect_df(rownames(geno_st), alleles = c("1","4"), seed)

  meta     <- sim_geno_meta(n_samples, n_markers, "selectiontools",
                             n_chrom = 3L, block_size = block_size, seed = seed)

  list(
    geno_st   = geno_st,
    geno_nuc  = geno_nuc,
    map       = map,
    block_map = bmap,
    eff_nuc   = eff_nuc,
    eff_st    = eff_st,
    geno_meta = meta
  )
}
