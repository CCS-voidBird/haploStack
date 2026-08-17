# =============================================================================
# simulate_dataset.R
#
# Synthetic dataset generators for haploStack.
#
# These helpers build self-contained, deterministic example datasets (genotype
# matrix, marker map, LD block map, marker effects, and assembled GenoMeta) with
# no dependence on external files or package state. They let users try the full
# build -> catalogue -> stacking-simulation workflow without supplying real data,
# and provide reproducible fixtures for examples and tests.
#
# Data-structure conventions
# --------------------------
# Genotype matrix  : rows = markers, cols = samples
#                    values = "A/A", "A/T", "T/T" (nucleotide) or
#                             "1/1", "1/4", "4/4" (SelectionTools numeric)
# Map data frame   : cols Marker, Chrom, Pos
# Block map        : cols Chrom, Marker, Block
# Effect data frame: cols Marker, Allele, Effect
# =============================================================================

#' Simulate a genotype matrix
#'
#' Generate a random biallelic genotype matrix (markers in rows, samples in
#' columns) in either nucleotide or SelectionTools numeric encoding.
#'
#' @param n_samples Number of samples (columns). Default 10.
#' @param n_markers Number of markers (rows). Default 20.
#' @param encoding  `"selectiontools"` (default) or `"nucleotide"`.
#' @param seed      Random seed for reproducibility. Default 42.
#' @return A character matrix (`n_markers` x `n_samples`) with named rows
#'   (`Mkr1`, `Mkr2`, ...) and columns (`S1`, `S2`, ...); each cell is a phased
#'   genotype such as `"1/4"` or `"A/T"`.
#' @seealso [sim_dataset()], [sim_geno_meta()]
#' @family dataset simulation
#' @examples
#' geno <- sim_geno_matrix(n_samples = 6, n_markers = 8)
#' dim(geno)
#' @export
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

#' Simulate a marker map
#'
#' Build a marker map that spreads the given markers across a number of
#' chromosomes with sorted base-pair positions.
#'
#' @param marker_names Character vector of marker IDs.
#' @param n_chrom      Number of chromosomes. Default 3.
#' @return A data frame with columns `Marker`, `Chrom`, `Pos`.
#' @seealso [sim_block_map()], [sim_dataset()]
#' @family dataset simulation
#' @examples
#' map <- sim_map_data(paste0("Mkr", 1:12), n_chrom = 3)
#' head(map)
#' @export
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

#' Simulate an LD block map
#'
#' Assign consecutive markers to blocks within each chromosome.
#'
#' @param map_data   A map data frame from [sim_map_data()].
#' @param block_size Number of markers per block. Default 4.
#' @return A data frame with columns `Chrom`, `Marker`, `Block`.
#' @seealso [sim_map_data()], [sim_dataset()]
#' @family dataset simulation
#' @examples
#' map  <- sim_map_data(paste0("Mkr", 1:12), n_chrom = 3)
#' bmap <- sim_block_map(map, block_size = 4)
#' head(bmap)
#' @export
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

#' Simulate a marker effect table
#'
#' Draw per-allele additive effects from a normal distribution.
#'
#' @param marker_names Character vector of marker IDs.
#' @param alleles      Character vector of allele states. Default `c("A", "T")`.
#' @param seed         Random seed for reproducibility. Default 42.
#' @return A data frame with columns `Marker`, `Allele`, `Effect`.
#' @seealso [sim_dataset()]
#' @family dataset simulation
#' @examples
#' eff <- sim_effect_df(paste0("Mkr", 1:5), alleles = c("1", "4"))
#' head(eff)
#' @export
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

#' Simulate a GenoMeta object
#'
#' Assemble a valid [`GenoMeta`][new_geno_meta] from simulated genotype and map
#' data. Wraps [sim_geno_matrix()], [sim_map_data()], and [sim_block_map()],
#' splits the 2-D genotype matrix into the 3-D allele array expected by
#' [new_geno_meta()], and optionally attaches a block map with
#' [attach_block_info()].
#'
#' @param n_samples  Number of samples. Default 10.
#' @param n_markers  Number of markers. Default 20.
#' @param encoding   `"selectiontools"` (default) or `"nucleotide"`.
#' @param n_chrom    Number of chromosomes. Default 2.
#' @param block_size Markers per block. `NULL` (default) leaves the block fields
#'   empty in the returned GenoMeta.
#' @param seed       Random seed for reproducibility. Default 42.
#' @return A valid `GenoMeta` list.
#' @seealso [sim_dataset()], [new_geno_meta()], [attach_block_info()]
#' @family dataset simulation
#' @examples
#' meta <- sim_geno_meta(n_samples = 6, n_markers = 8, block_size = 4)
#' is_geno_meta(meta)
#' @export
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

  # Split "1/4" -> ["1","4"] per cell into a 3-D array [n_samples x n_markers x 2]
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

#' Simulate a complete workflow dataset
#'
#' Produce every component needed to exercise the haploStack workflow end to
#' end: genotype matrices (both encodings), a marker map, an LD block map,
#' marker effect tables, and an assembled [`GenoMeta`][new_geno_meta].
#'
#' @param n_samples  Number of samples. Default 10.
#' @param n_markers  Number of markers. Default 20.
#' @param block_size Markers per block. Default 4.
#' @param seed       Random seed for reproducibility. Default 42.
#' @return A named list with elements `geno_st`, `geno_nuc`, `map`,
#'   `block_map`, `eff_nuc`, `eff_st`, and `geno_meta`.
#' @seealso [sim_geno_meta()], [run_haplotype_catalogue()],
#'   [run_stacking_simulation()]
#' @family dataset simulation
#' @examples
#' ds <- sim_dataset(n_samples = 6, n_markers = 8, block_size = 4)
#' names(ds)
#' is_geno_meta(ds$geno_meta)
#' @export
sim_dataset <- function(n_samples  = 10L,
                        n_markers  = 20L,
                        block_size = 4L,
                        seed       = 42L) {

  geno_st  <- sim_geno_matrix(n_samples, n_markers, "selectiontools", seed)
  geno_nuc <- sim_geno_matrix(n_samples, n_markers, "nucleotide",     seed + 1L)

  map      <- sim_map_data(rownames(geno_st), n_chrom = 3L)
  bmap     <- sim_block_map(map, block_size)

  eff_nuc  <- sim_effect_df(rownames(geno_st), alleles = c("A", "T"), seed)
  eff_st   <- sim_effect_df(rownames(geno_st), alleles = c("1", "4"), seed)

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
