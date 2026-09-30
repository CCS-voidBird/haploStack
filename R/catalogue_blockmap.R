# =============================================================================
# catalogue_blockmap.R
#
# build_block_location_info() — compute block spans from block map + SNP map
# expand_block_info()         — convert (Chr,Block,Start,End) → (Chr,Marker,Block)
#
# Q3 (2026-06-08):
#   - Canonical internal block-marker format is strict (Chr, Marker, Block).
#   - expand_block_info() added.
#   - build_block_location_info() now returns a named list (block_info,
#     block_marker_map) and writes two files when output_dir is given.
#
# 0.7.10: define_haplotype_blocks(), the SelectionTools file route, moved to
# R/blocks_selectiontools.R as st_ld_blocks_files(); the old name is a
# deprecated alias in R/blocks_deprecated.R. Local LD blocking is a side route:
# see R/blocks_ld_common.R.
# =============================================================================

# .compute_block_location(bm, map)
#
# Pure computation: join a (Chrom, Marker, Block) table to a (Marker, Chrom,
# Pos) map and derive per-block coordinate statistics. Called by both
# build_block_location_info() and attach_block_info() to avoid code duplication.
#
# Returns list(block_info = data.table, block_marker_map = data.table).
.compute_block_location <- function(bm, map) {
  bm  <- as.data.table(bm)
  map <- as.data.table(map)

  if (!all(c("Chrom", "Marker", "Block") %in% names(bm)))
    setnames(bm, names(bm)[seq_len(3L)], c("Chrom", "Marker", "Block"))
  if (!all(c("Marker", "Chrom", "Pos") %in% names(map)))
    setnames(map, names(map)[seq_len(3L)], c("Marker", "Chrom", "Pos"))

  bm[,  Chrom := as.character(Chrom)]
  map[, Chrom := as.character(Chrom)]
  map[, Pos   := as.numeric(Pos)]

  merged <- dplyr::left_join(bm, map, by = c("Chrom", "Marker"))
  merged <- subset(merged, !is.na(Pos))

  block_info <- as.data.table(merged) %>%
    dplyr::group_by(Chrom, Block) %>%
    dplyr::summarise(
      Start        = min(as.numeric(Pos)),
      End          = max(as.numeric(Pos)),
      Pos          = round((min(as.numeric(Pos)) + max(as.numeric(Pos))) / 2),
      Num_SNPs     = dplyr::n_distinct(Marker),
      .groups      = "drop"
    ) %>%
    dplyr::mutate(Block_length = End - Start + 1L) %>%
    as.data.table()

  block_marker_map <- unique(
    as.data.table(merged)[, .(Chrom, Marker, Block)]
  )

  list(block_info = block_info, block_marker_map = block_marker_map)
}


#' Build block location information from a GenoMeta object
#'
#' Computes per-block start position, end position, midpoint, chromosome,
#' number of SNPs, and block length (bp) from the block and map data stored
#' inside a \code{GenoMeta}. The result is used for Manhattan plots and
#' block-level summary statistics.
#'
#' The GenoMeta must already contain block assignments — call
#' \code{\link{attach_block_info}} first if \code{block_marker_map_info} is
#' \code{NULL}.
#'
#' When \code{output_dir} is provided, two files are written:
#' \describe{
#'   \item{\code{<prefix>_block_info.csv}}{Block coordinate table:
#'     \code{Chrom, Block, Start, End} (plus the signal-marker columns when the
#'     GenoMeta carries a selection).}
#'   \item{\code{<prefix>_block_marker_map.csv}}{Canonical marker-level table:
#'     \code{Chrom, Marker, Block}.}
#' }
#'
#' The coordinates are recomputed from marker positions, so any signal-marker
#' selection stored by \code{\link{attach_signal_markers}}
#' (\code{Signal_Marker}, \code{Signal_Score}, \code{N_Signal}) is joined back on
#' \code{Block} rather than being silently dropped.
#'
#' @param geno_meta  A valid \code{GenoMeta} with \code{block_marker_map_info}
#'   populated. Create via \code{\link{attach_block_info}}.
#' @param output_dir Optional directory path for CSV output. Default
#'   \code{NULL}.
#' @param prefix     Filename prefix used when \code{output_dir} is set.
#'   Default \code{"haplostack"}.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{block_info}}{data.table: \code{Chrom}, \code{Block},
#'     \code{Start}, \code{End}, \code{Pos} (midpoint), \code{Num_SNPs},
#'     \code{Block_length}.}
#'   \item{\code{block_marker_map}}{data.table: \code{Chrom}, \code{Marker},
#'     \code{Block}.}
#' }
#'
#' @examples
#' \dontrun{
#' meta <- read_st_genotype("geno.csv")
#' meta <- attach_block_info(meta, read.csv("block_map.csv"))
#' loc  <- build_block_location_info(meta, output_dir = "output/")
#' head(loc$block_info)
#' }
#'
#' @seealso \code{\link{attach_block_info}}, \code{\link{hapS_ld_blocks}} and
#'   \code{\link{st_ld_blocks}} (local LD blocking, a side route),
#'   \code{\link{expand_block_info}}
#'
#' @export
build_block_location_info <- function(geno_meta,
                                       output_dir = NULL,
                                       prefix     = "haplostack") {
  if (!is_geno_meta(geno_meta))
    stop("build_block_location_info: 'geno_meta' must be a valid GenoMeta.")
  if (is.null(geno_meta$block_marker_map_info))
    stop("build_block_location_info: GenoMeta has no block_marker_map_info. ",
         "Call attach_block_info() first to assign markers to blocks.")

  result <- .compute_block_location(geno_meta$block_marker_map_info,
                                    geno_meta$map_info)

  # The coordinates above are recomputed from marker positions, so a stored
  # signal-marker selection (attach_signal_markers()) would be lost; join it
  # back on Block so both the returned table and the CSV carry it.
  sig_cols <- intersect(.signal_marker_cols(), names(geno_meta$block_info))
  if (length(sig_cols) > 0L) {
    sig <- unique(as.data.table(geno_meta$block_info)[
      , c("Block", sig_cols), with = FALSE])
    sig[, Block := as.character(Block)]
    result$block_info[sig, on = "Block",
                      (sig_cols) := mget(paste0("i.", sig_cols))]
  }

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE)
    binfo_file <- file.path(output_dir, paste0(prefix, "_block_info.csv"))
    bmap_file  <- file.path(output_dir, paste0(prefix, "_block_marker_map.csv"))
    fwrite(result$block_info[, c("Chrom", "Block", "Start", "End", sig_cols),
                             with = FALSE],
           file = binfo_file, sep = ",", row.names = FALSE, quote = FALSE)
    fwrite(result$block_marker_map,
           file = bmap_file, sep = ",", row.names = FALSE, quote = FALSE)
    message("Block info written to: ", binfo_file)
    message("Block marker map written to: ", bmap_file)
  }

  result
}


#' Convert external block coordinates to a marker-level block map
#'
#' Takes an external block coordinate table \code{(Chrom, Block, Start, End)}
#' and finds all markers in the \code{GenoMeta}'s \code{map_info} whose
#' position falls within each block's genomic range. Returns the canonical
#' internal \code{(Chrom, Marker, Block)} table.
#'
#' Use this function when block coordinates are defined externally (e.g. from a
#' sliding-window analysis or a QTL interval). Pass the result directly to
#' \code{\link{attach_block_info}} to update the GenoMeta.
#'
#' @note Start/End may differ from actual marker positions — blocks can span
#'   regions wider than the available markers. Markers outside all block ranges
#'   are silently excluded. A \code{warning()} is emitted for blocks that
#'   contain no markers and for markers that fall within more than one block.
#'
#' @param geno_meta      A valid \code{GenoMeta} list. The \code{map_info}
#'   field is used to match markers to block ranges.
#' @param block_info_df  data.frame with columns \code{Chrom}, \code{Block},
#'   \code{Start}, \code{End} (first four columns used positionally if names
#'   differ).
#'
#' @return data.frame with columns \code{Chrom}, \code{Marker}, \code{Block}
#'   (strict canonical order), one row per marker × block match. Pass this
#'   directly to \code{\link{attach_block_info}} with
#'   \code{source_type = "block_marker_map"}.
#'
#' @examples
#' \dontrun{
#' meta <- read_st_genotype("geno.csv")
#'
#' block_coords <- data.frame(
#'   Chrom = "1", Block = c("B1", "B2"),
#'   Start = c(100000, 500000), End = c(300000, 800000)
#' )
#' bmap <- expand_block_info(meta, block_coords)
#' meta <- attach_block_info(meta, bmap)
#' }
#'
#' @seealso \code{\link{attach_block_info}}, \code{\link{build_block_location_info}}
#'
#' @export
expand_block_info <- function(map_info, block_info_df) {
  # if (!is_geno_meta(geno_meta))
  #   stop("expand_block_info: 'geno_meta' must be a valid GenoMeta.")

  bi  <- as.data.frame(block_info_df, stringsAsFactors = FALSE)
  map <- as.data.frame(map_info, stringsAsFactors = FALSE)


  # Normalise column names positionally if canonical names absent
  if (!all(c("Chrom", "Block", "Start", "End") %in% names(bi)))
    names(bi)[seq_len(min(4L, ncol(bi)))] <-
      c("Chrom", "Block", "Start", "End")[seq_len(min(4L, ncol(bi)))]
  if (!all(c("Chrom", "Marker", "Pos") %in% names(map)))
    names(map)[seq_len(min(3L, ncol(map)))] <-
      c("Chrom", "Marker", "Pos")[seq_len(min(3L, ncol(map)))]

  bi$Start <- as.numeric(bi$Start)
  bi$End   <- as.numeric(bi$End)
  map$Pos  <- as.numeric(map$Pos)

  result_rows <- vector("list", nrow(bi))

  for (i in seq_len(nrow(bi))) {
    chr_i   <- bi$Chrom[i]
    block_i <- bi$Block[i]
    start_i <- bi$Start[i]
    end_i   <- bi$End[i]

    hits <- map[map$Chrom == chr_i &
                  map$Pos >= start_i &
                  map$Pos <= end_i, , drop = FALSE]

    if (nrow(hits) == 0L) {
      warning("expand_block_info: block '", block_i,
              "' (", chr_i, ":", start_i, "-", end_i,
              ") contains no markers — omitted from output.")
      next
    }

    result_rows[[i]] <- data.frame(
      Chrom  = chr_i,
      Marker = hits$Marker,
      Block  = block_i,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, result_rows[!vapply(result_rows, is.null, logical(1L))])
  if (is.null(out)) {
    warning("expand_block_info: no markers matched any block.")
    return(data.frame(Chrom  = character(0L), Marker = character(0L),
                      Block  = character(0L), stringsAsFactors = FALSE))
  }

  dup_markers <- out$Marker[duplicated(out$Marker)]
  if (length(dup_markers) > 0L)
    warning("expand_block_info: ", length(unique(dup_markers)),
            " marker(s) fall in more than one block: ",
            paste(head(unique(dup_markers), 5L), collapse = ", "),
            if (length(unique(dup_markers)) > 5L) " ..." else "")

  out[, c("Chrom", "Marker", "Block"), drop = FALSE]
}
