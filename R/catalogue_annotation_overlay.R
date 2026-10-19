# =============================================================================
# catalogue_annotation_overlay.R
#
# BAN (thread #033) stage 6, the Q25.o half — user-supplied marker-level
# annotations / QTLs as a data frame, normalised once and reusable everywhere:
#
#   standardize_marker_annotations() — the single normaliser (Marker route or
#                                      Chrom+Pos route, match_bp matching,
#                                      unplaced rows kept + warned)
#   attach_marker_annotations()      — store the resolved table as
#                                      meta$annotation$custom (+ provenance)
#
# Catalogue-layer file: it consumes the container and the df contract, never
# the annotation compute module. track_marker_labels()
# (catalogue_circos_meta.R) is the circos consumer; the `annotations =`
# overlay arguments on the catalogue plotters are the remaining stage-6 item.
# =============================================================================


#' Normalise a marker-level annotation / QTL data frame
#'
#' The single owner of the Q25.o contract. Accepted inputs: a data.frame with
#' \code{Name} plus either \code{Marker} (resolved against
#' \code{meta$map_info}) or \code{Chrom} + \code{Pos} (bp; matched to the
#' nearest marker on that chromosome within \code{match_bp} when \code{meta}
#' is given). Optional \code{Type} (free text — "QTL", "gene", "candidate" …)
#' and \code{Trait} pass through.
#'
#' Rows that resolve to no marker — an unknown \code{Marker} id, or no marker
#' within \code{match_bp} of the position — are \strong{kept}, with
#' \code{Marker}/\code{Block} \code{NA}, and one warning names them: a named
#' locus is still drawable by position even when no marker tags it.
#'
#' @param x        The annotation data.frame.
#' @param meta     A \code{GenoMeta}; required for the \code{Marker} route and
#'   for position-to-marker matching (omit it and a positional table is
#'   normalised without matching).
#' @param match_bp Matching tolerance in bp (default the shared 10 kb,
#'   \code{.HS_ANNOT_DEFAULT_BP} — the same constant as
#'   \code{compute_block_annotation()}'s \code{flank_bp} default).
#'
#' @return A data.frame: \code{Marker, Chrom, Pos, Block, Name, Type, Trait,
#'   matched_dist_bp} — one row per input row, input order kept.
#'   \code{matched_dist_bp} is 0 on the \code{Marker} route, the
#'   marker-to-position distance on a positional match, \code{NA} unplaced.
#'
#' @seealso \code{\link{attach_marker_annotations}},
#'   \code{\link{track_marker_labels}}
#' @export
standardize_marker_annotations <- function(x, meta = NULL,
                                           match_bp = .HS_ANNOT_DEFAULT_BP) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"Name" %in% names(df))
    stop("standardize_marker_annotations: 'x' needs a Name column.")
  has_marker <- "Marker" %in% names(df)
  has_pos    <- all(c("Chrom", "Pos") %in% names(df))
  if (!has_marker && !has_pos)
    stop("standardize_marker_annotations: 'x' needs Marker, or Chrom + Pos.")
  if (!is.numeric(match_bp) || length(match_bp) != 1L || is.na(match_bp) ||
      match_bp < 0)
    stop("standardize_marker_annotations: 'match_bp' must be a single ",
         "number >= 0.")

  n <- nrow(df)
  out <- data.frame(
    Marker = if (has_marker) as.character(df$Marker) else NA_character_,
    Chrom  = if (has_pos) as.character(df$Chrom) else NA_character_,
    Pos    = if (has_pos) as.numeric(df$Pos) else NA_real_,
    Block  = NA_character_,
    Name   = as.character(df$Name),
    Type   = if ("Type" %in% names(df)) as.character(df$Type)
             else NA_character_,
    Trait  = if ("Trait" %in% names(df)) as.character(df$Trait)
             else NA_character_,
    matched_dist_bp = NA_real_,
    stringsAsFactors = FALSE)

  if (!is.null(meta)) {
    if (!is_geno_meta(meta))
      stop("standardize_marker_annotations: 'meta' is not a valid GenoMeta.")
    map <- as.data.frame(meta$map_info, stringsAsFactors = FALSE)
    map$Pos <- as.numeric(map$Pos)

    # Marker route first: known ids take their mapped position.
    if (has_marker) {
      idx <- match(out$Marker, map$Marker)
      hit <- !is.na(idx)
      out$Chrom[hit] <- map$Chrom[idx[hit]]
      out$Pos[hit]   <- map$Pos[idx[hit]]
      out$matched_dist_bp[hit] <- 0
      out$Marker[!hit] <- NA_character_
    }
    # Positional route for everything still unresolved but located.
    todo <- is.na(out$Marker) & !is.na(out$Chrom) & !is.na(out$Pos)
    for (i in which(todo)) {
      cand <- map[map$Chrom == out$Chrom[i], , drop = FALSE]
      if (nrow(cand) == 0L) next
      d <- abs(cand$Pos - out$Pos[i])
      j <- which.min(d)
      if (d[j] <= match_bp) {
        out$Marker[i]          <- cand$Marker[j]
        out$matched_dist_bp[i] <- d[j]
      }
    }
    # Block of the resolved marker.
    bmm <- meta$block_marker_map_info
    if (!is.null(bmm)) {
      bidx <- match(out$Marker, as.character(bmm$Marker))
      out$Block <- ifelse(is.na(bidx), NA_character_,
                          as.character(bmm$Block)[bidx])
    }
  } else if (has_marker && !has_pos) {
    stop("standardize_marker_annotations: the Marker route needs 'meta' ",
         "to resolve positions.")
  }

  unplaced <- is.na(out$Marker)
  if (!is.null(meta) && any(unplaced))
    warning("standardize_marker_annotations: ", sum(unplaced), " of ", n,
            " row(s) match no marker within ",
            format(match_bp, big.mark = ","), " bp — kept with Marker/Block ",
            "NA (first: ",
            paste(utils::head(out$Name[unplaced], 5L), collapse = ", "),
            if (sum(unplaced) > 5L) " ..." else "", ").")
  out
}


#' Attach marker-level annotations / QTLs to the container
#'
#' Normalises \code{x} through \code{\link{standardize_marker_annotations}}
#' and stores the resolved table as \code{meta$annotation$custom}, with its
#' own provenance under \code{metadata$annotation$custom}. The element is
#' independent of the GFF3 route: it attaches to a container that has no gene
#' annotation (then \code{meta$annotation} holds \code{custom} alone).
#'
#' Carry-forward matches the rest of the annotation element (Q25.i):
#' \code{\link{filter_geno_meta}} prunes rows whose resolved marker left the
#' panel (position-only rows stay); \code{\link{attach_block_info}} drops the
#' whole annotation with a message; other transforms keep it.
#'
#' @param meta     A valid \code{GenoMeta}.
#' @param x        The annotation data.frame (Q25.o contract), or an
#'   already-standardized table.
#' @param match_bp Passed to the standardizer (default the shared 10 kb).
#' @param source   Optional free-text provenance note (e.g. the file the df
#'   came from).
#'
#' @return The \code{GenoMeta} with \code{$annotation$custom} set.
#' @seealso \code{\link{standardize_marker_annotations}},
#'   \code{\link{track_marker_labels}}
#' @export
attach_marker_annotations <- function(meta, x,
                                      match_bp = .HS_ANNOT_DEFAULT_BP,
                                      source = NULL) {
  if (!is_geno_meta(meta))
    stop("attach_marker_annotations: 'meta' is not a valid GenoMeta.")
  std_cols <- c("Marker", "Chrom", "Pos", "Block", "Name", "Type", "Trait",
                "matched_dist_bp")
  std <- if (is.data.frame(x) && all(std_cols %in% names(x)))
    as.data.frame(x, stringsAsFactors = FALSE)[, std_cols]
  else
    standardize_marker_annotations(x, meta = meta, match_bp = match_bp)

  if (is.null(meta$annotation)) meta$annotation <- list()
  meta$annotation$custom <- std
  if (is.null(meta$metadata$annotation)) meta$metadata$annotation <- list()
  meta$metadata$annotation$custom <- list(
    source          = source %||% NA_character_,
    match_bp        = as.numeric(match_bp),
    n               = nrow(std),
    n_matched       = sum(!is.na(std$Marker)),
    n_unplaced      = sum(is.na(std$Marker)),
    date            = format(Sys.Date()),
    package_version = as.character(utils::packageVersion("haploStack")))
  meta
}
