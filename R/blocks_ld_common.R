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
    .hs_abort(caller, ": 'ld_threshold' must be one number in (0, 1] (an r^2).",
         if (!is.null(zero_note)) paste0(" ", zero_note), call. = FALSE)
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) ||
      tolerance < 0 || tolerance != round(tolerance))
    .hs_abort(caller, ": 'tolerance' must be one whole number >= 0 (the number ",
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
#   attach = FALSE: list(block_marker_map, block_info, provenance), and with
#                   `output_dir` the two canonical files (.write_block_files());
#   attach = TRUE:  the container with the blocks attached (through
#                   attach_block_info(), which rebuilds it and keeps
#                   marker_effects, phenotype and the FAM) and, with
#                   `output_dir`, the same two files.
# bm:     data.frame(Chrom, Marker, Block) of the blocked markers
# map_df: data.frame(Marker, Chrom, Pos) of the same markers
# rule:   list(backend, ld_threshold, tolerance, ..., n_unplaced): the source,
#         label and params of the provenance record (metadata$blocks, D2-6)
.ld_finish <- function(geno_meta, bm, map_df, attach, output_dir, prefix,
                       caller, rule) {
  bi   <- .compute_block_location(bm, map_df)$block_info
  prov <- .block_provenance(
    bm, geno_meta$map_info,
    source      = rule$backend,
    label       = sprintf("%s_r%s_t%d", rule$backend, format(rule$ld_threshold),
                          as.integer(rule$tolerance)),
    params      = rule[setdiff(names(rule), c("backend", "n_unplaced"))],
    n_unmatched = rule$n_unplaced %||% 0L)
  res <- list(block_marker_map = data.table::as.data.table(bm),
              block_info       = data.table::as.data.table(bi),
              provenance       = prov)
  if (!isTRUE(attach)) {
    if (!is.null(output_dir)) .write_block_files(res, output_dir, prefix)
    return(res)
  }

  out <- attach_block_info(geno_meta, res)
  if (!is_geno_meta(out, verbose = TRUE))
    .hs_abort(caller, ": the container failed validation after attaching blocks.",
         call. = FALSE)
  if (!is.null(output_dir))
    invisible(build_block_location_info(out, output_dir = output_dir,
                                        prefix = prefix))
  out
}
