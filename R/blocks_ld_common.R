# =============================================================================
# blocks_ld_common.R — what the local LD-block resolvers share
#
# haploStack does not need LD blocking. The recommended route is a block map
# the user already has, given to read_plink_ped(block_file =),
# build_geno_meta(blocks =) or attach_block_info(); it takes no LD parameter.
# Only when blocks must be built locally does a resolver come in, and each is
# named after its backend so that it cannot be mistaken for the main route or
# for the other resolver:
#
#   hapS_ld_blocks()   HapSelect       R/blocks_hapselect.R      the default
#   st_ld_blocks()     SelectionTools  R/blocks_selectiontools.R experimental
#   st_ld_blocks_files()  SelectionTools on its own file formats (same file)
#
# The two backends do not partition alike at the same threshold / tolerance
# (thread #034): a block set is comparable only with one from the same
# backend. They return the same shapes, built by the helpers below.
# =============================================================================


# The blocking rule, checked before any backend is touched. `zero_note` adds a
# backend's reason to refuse 0 (SelectionTools aborts the R session on it).
.ld_check_rule <- function(ld_threshold, tolerance, caller, zero_note = NULL) {
  if (!is.numeric(ld_threshold) || length(ld_threshold) != 1L ||
      is.na(ld_threshold) || ld_threshold <= 0 || ld_threshold > 1)
    stop(caller, ": 'ld_threshold' must be one number in (0, 1] (an r^2).",
         if (!is.null(zero_note)) paste0(" ", zero_note), call. = FALSE)
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) ||
      tolerance < 0 || tolerance != round(tolerance))
    stop(caller, ": 'tolerance' must be one whole number >= 0 (the number ",
         "of sub-threshold markers a block may absorb).", call. = FALSE)
  invisible(TRUE)
}

# Markers a resolver can place: a finite position > 0 on a named, placed
# chromosome. Both backends read markers in position order, so an unplaced
# marker would be blocked by file order; it is left out and becomes a
# singleton block on attach, as any marker absent from a block map does.
.ld_placed_markers <- function(map_info) {
  chrom <- trimws(as.character(map_info$Chrom))
  pos   <- suppressWarnings(as.numeric(map_info$Pos))
  ok <- !is.na(chrom) & nzchar(chrom) &
        !tolower(chrom) %in% c("0", "un", "unknown", "chrun", "unplaced") &
        is.finite(pos) & pos > 0
  as.character(map_info$Marker)[ok]
}

# "N blocks over M markers (S single-marker, largest L)" for the summaries.
.ld_sizes_note <- function(bm) {
  sizes <- table(bm$Block)
  sprintf("%d blocks over %d markers (%d single-marker, largest %d)",
          length(sizes), nrow(bm), sum(sizes == 1L), max(sizes))
}

# The common tail of every resolver.
#   attach = FALSE: list(block_marker_map, block_info), and with `output_dir`
#                   the two canonical files <prefix>_block_info.csv and
#                   <prefix>_block_marker_map.csv;
#   attach = TRUE:  the container with the blocks attached. attach_block_info()
#                   rebuilds it and drops marker_effects and phenotype and
#                   resets the FAM (INDEX discrepancy 20); the samples and
#                   markers are unchanged, so they are carried across.
# bm:     data.frame(Chrom, Marker, Block) of the blocked markers
# map_df: data.frame(Marker, Chrom, Pos) of the same markers
.ld_finish <- function(geno_meta, bm, map_df, attach, output_dir, prefix,
                       caller) {
  if (!isTRUE(attach)) {
    bi <- .compute_block_location(bm, map_df)$block_info
    if (!is.null(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      # scipen: a round position (88000000) must not be written as 8.8e+07
      data.table::fwrite(data.table::as.data.table(bi),
                         file.path(output_dir, paste0(prefix, "_block_info.csv")),
                         scipen = 100L)
      data.table::fwrite(data.table::as.data.table(bm),
                         file.path(output_dir, paste0(prefix, "_block_marker_map.csv")))
    }
    return(list(block_marker_map = data.table::as.data.table(bm),
                block_info       = data.table::as.data.table(bi)))
  }

  out <- attach_block_info(geno_meta, bm, source_type = "block_marker_map")
  for (f in c("marker_effects", "phenotype", "sample_info"))
    if (!is.null(geno_meta[[f]])) out[[f]] <- geno_meta[[f]]
  if (!is_geno_meta(out, verbose = TRUE))
    stop(caller, ": the container failed validation after attaching blocks.",
         call. = FALSE)
  if (!is.null(output_dir))
    invisible(build_block_location_info(out, output_dir = output_dir,
                                        prefix = prefix))
  out
}
