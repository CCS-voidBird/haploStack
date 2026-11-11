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
# Every draw runs under .with_seed() (R/state.R): the output for a given seed
# is what set.seed(seed) gave up to 0.8.8, and the caller's own random stream
# is left as it was. The exported generators are thin wrappers around the
# unseeded .sim_*_draw() helpers so that a composite generator can run several
# draws on one seeded stream, in the same order as before.
#
# Data-structure conventions
# --------------------------
# Genotype matrix  : rows = markers, cols = samples
#                    values = "A/A", "A/T", "T/T" (nucleotide, the default
#                    since 0.9.0) or "1/1", "1/4", "4/4" (SelectionTools codes,
#                    deprecated: a GenoMeta stores nucleotide alleles only, so
#                    a container asked for in that encoding is built with the
#                    same draws as A/T — C-5)
# Map data frame   : cols Marker, Chrom, Pos
# Block map        : cols Chrom, Marker, Block
# Effect data frame: cols Marker, Allele, Effect
# =============================================================================

.SIM_ENCODINGS <- c("nucleotide", "selectiontools")

.sim_alleles <- function(encoding) {
  if (encoding == "selectiontools") c("1", "4") else c("A", "T")
}

# Unseeded draws ---------------------------------------------------------------

.sim_geno_draw <- function(n_samples, n_markers, alleles) {
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

.sim_map_draw <- function(marker_names, n_chrom) {
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

.sim_effect_draw <- function(marker_names, alleles) {
  pairs <- expand.grid(Marker = marker_names, Allele = alleles,
                       stringsAsFactors = FALSE)
  pairs$Effect <- stats::rnorm(nrow(pairs), mean = 0, sd = 0.5)
  pairs <- pairs[.hs_order(pairs$Marker, pairs$Allele), ]
  rownames(pairs) <- NULL
  pairs
}

# Genotype matrix and map on one seeded stream (the order sim_geno_meta() has
# always used: genotypes first, then positions).
.sim_geno_and_map <- function(n_samples, n_markers, encoding, n_chrom, seed) {
  alleles <- .sim_alleles(encoding)
  .with_seed(seed, {
    geno <- .sim_geno_draw(n_samples, n_markers, alleles)
    map  <- .sim_map_draw(rownames(geno), n_chrom)
    list(geno = geno, map = map)
  })
}

# Assemble a GenoMeta from a drawn genotype matrix and map. A matrix drawn in
# SelectionTools codes is translated first (1 -> A, 4 -> T: the same draws as
# the nucleotide alleles, since .sim_geno_draw() samples by position), so the
# container is always nucleotide (0.9.0, C-5).
.sim_meta_from <- function(geno, map, encoding, block_size = NULL) {
  if (identical(encoding, "selectiontools")) {
    geno     <- .st_decode_cells(geno, split_char = "/")
    encoding <- "nucleotide"
  }
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

  char_dict <- c("A" = 1L, "T" = 2L)

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

# Exported generators -----------------------------------------------------------

#' Simulate a genotype matrix
#'
#' Generate a random biallelic genotype matrix (markers in rows, samples in
#' columns) in nucleotide encoding (`A/T`). The SelectionTools numeric
#' encoding (`1/4`) is deprecated since 0.9.0 — a `GenoMeta` stores
#' nucleotide alleles only — and still returned here, with a warning, for the
#' same draws (the allele is sampled by position, so `1` is `A` and `4` is
#' `T`).
#'
#' @param n_samples Number of samples (columns). Default 10.
#' @param n_markers Number of markers (rows). Default 20.
#' @param encoding  `"nucleotide"` (default since 0.9.0) or
#'   `"selectiontools"` (deprecated).
#' @param seed      Random seed for reproducibility. Default 42. The caller's
#'   own random stream is left untouched; `NULL` draws from it instead.
#' @return A character matrix (`n_markers` x `n_samples`) with named rows
#'   (`Mkr1`, `Mkr2`, ...) and columns (`S1`, `S2`, ...); each cell is a phased
#'   genotype such as `"A/T"`.
#' @seealso [sim_dataset()], [sim_geno_meta()]
#' @family dataset simulation
#' @examples
#' geno <- sim_geno_matrix(n_samples = 6, n_markers = 8)
#' dim(geno)
#' @export
sim_geno_matrix <- function(n_samples  = 10L,
                            n_markers  = 20L,
                            encoding   = c("nucleotide", "selectiontools"),
                            seed       = 42L) {
  encoding <- match.arg(encoding)
  if (encoding == "selectiontools") {
    .st_gate("sim_geno_matrix", per_call = FALSE)        # C-6
    .st_encoding_deprecated("sim_geno_matrix",
                            then = "the 1/4-coded matrix is returned for now (same draws as A/T)")
  }
  .with_seed(seed, .sim_geno_draw(n_samples, n_markers, .sim_alleles(encoding)))
}

#' Simulate a marker map
#'
#' Build a marker map that spreads the given markers across a number of
#' chromosomes with sorted base-pair positions.
#'
#' @param marker_names Character vector of marker IDs.
#' @param n_chrom      Number of chromosomes. Default 3.
#' @param seed         Random seed. `NULL` (default) draws the positions from
#'   the caller's random stream, as before; a seed makes the map reproducible
#'   on its own and leaves the caller's stream untouched.
#' @return A data frame with columns `Marker`, `Chrom`, `Pos`.
#' @seealso [sim_block_map()], [sim_dataset()]
#' @family dataset simulation
#' @examples
#' map <- sim_map_data(paste0("Mkr", 1:12), n_chrom = 3)
#' head(map)
#' @export
sim_map_data <- function(marker_names,
                         n_chrom = 3L,
                         seed    = NULL) {
  .with_seed(seed, .sim_map_draw(marker_names, n_chrom))
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
#' @param seed         Random seed for reproducibility. Default 42. The
#'   caller's own random stream is left untouched; `NULL` draws from it.
#' @return A data frame with columns `Marker`, `Allele`, `Effect`.
#' @seealso [sim_dataset()]
#' @family dataset simulation
#' @examples
#' eff <- sim_effect_df(paste0("Mkr", 1:5))
#' head(eff)
#' @export
sim_effect_df <- function(marker_names,
                          alleles = c("A", "T"),
                          seed    = 42L) {
  .with_seed(seed, .sim_effect_draw(marker_names, alleles))
}

#' Simulate a GenoMeta object
#'
#' Assemble a valid [`GenoMeta`][new_geno_meta] from simulated genotype and map
#' data. Draws the genotype matrix and the map on one seeded stream (the same
#' draws as [sim_geno_matrix()] followed by [sim_map_data()]), splits the 2-D
#' genotype matrix into the 3-D allele array expected by [new_geno_meta()],
#' and optionally attaches a block map with [attach_block_info()].
#'
#' @param n_samples  Number of samples. Default 10.
#' @param n_markers  Number of markers. Default 20.
#' @param encoding   `"nucleotide"` (default since 0.9.0). `"selectiontools"`
#'   is deprecated: a `GenoMeta` stores nucleotide alleles only, so the
#'   container is built with the same draws in `A/T`, with a warning.
#' @param n_chrom    Number of chromosomes. Default 2.
#' @param block_size Markers per block. `NULL` (default) leaves the block fields
#'   empty in the returned GenoMeta.
#' @param seed       Random seed for reproducibility. Default 42. The caller's
#'   own random stream is left untouched.
#' @return A valid `GenoMeta` list.
#' @seealso [sim_dataset()], [new_geno_meta()], [attach_block_info()]
#' @family dataset simulation
#' @examples
#' meta <- sim_geno_meta(n_samples = 6, n_markers = 8, block_size = 4)
#' is_geno_meta(meta)
#' @export
sim_geno_meta <- function(n_samples  = 10L,
                          n_markers  = 20L,
                          encoding   = c("nucleotide", "selectiontools"),
                          n_chrom    = 2L,
                          block_size = NULL,
                          seed       = 42L) {
  encoding <- match.arg(encoding)
  if (encoding == "selectiontools") {
    .st_gate("sim_geno_meta", per_call = FALSE)          # C-6
    .st_encoding_deprecated("sim_geno_meta",
                            then = "the container is built with encoding = \"nucleotide\" (the same draws, A/T in place of 1/4)")
    encoding <- "nucleotide"
  }
  d <- .sim_geno_and_map(n_samples, n_markers, encoding, n_chrom, seed)
  .sim_meta_from(d$geno, d$map, encoding, block_size)
}

#' Simulate a complete workflow dataset
#'
#' Produce every component needed to exercise the haploStack workflow end to
#' end: genotype matrices, a marker map, an LD block map, marker effect
#' tables, and an assembled [`GenoMeta`][new_geno_meta].
#'
#' The returned `map` is the map inside `geno_meta` (`geno_meta$map_info`
#' holds the same positions), and `block_map` is the block map attached to
#' it. Up to 0.8.8 `map` was drawn separately, so its positions differed from
#' the container's.
#'
#' `geno_st` and `eff_st` are the draw in SelectionTools codes (`1/4`), kept
#' for scripts that still read them; since 0.9.0 the container is built from
#' the same draw in nucleotides (`1` is `A`, `4` is `T`), because a `GenoMeta`
#' stores nucleotide alleles only. `geno_nuc` is a second, independent draw
#' (seed + 1), as before.
#'
#' @param n_samples  Number of samples. Default 10.
#' @param n_markers  Number of markers. Default 20.
#' @param block_size Markers per block. Default 4.
#' @param seed       Random seed for reproducibility. Default 42. The caller's
#'   own random stream is left untouched.
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

  # Same draws as sim_geno_meta(n_samples, n_markers, "selectiontools",
  # n_chrom = 3, seed): the container and `map` share one map.
  d        <- .sim_geno_and_map(n_samples, n_markers, "selectiontools", 3L, seed)
  geno_st  <- d$geno
  map      <- d$map
  geno_nuc <- .with_seed(seed + 1L,
                         .sim_geno_draw(n_samples, n_markers, c("A", "T")))
  bmap     <- sim_block_map(map, block_size)

  eff_nuc  <- .with_seed(seed, .sim_effect_draw(rownames(geno_st), c("A", "T")))
  eff_st   <- .with_seed(seed, .sim_effect_draw(rownames(geno_st), c("1", "4")))

  meta     <- .sim_meta_from(geno_st, map, "selectiontools", block_size)

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
