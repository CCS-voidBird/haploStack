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
#' @param chr_map  Optional named character vector, container chromosome ->
#'   the annotation's spelling (the one direction every haploStack
#'   \code{chr_map} uses, e.g. \code{c("1" = "1H")}). Positions are matched
#'   on the mapped names; the stored \code{Chrom} of a matched row is the
#'   container's spelling.
#'
#' @return A data.frame: \code{Marker, Chrom, Pos, Block, Name, Type, Trait,
#'   matched_dist_bp} — one row per input row, input order kept.
#'   \code{matched_dist_bp} is 0 on the \code{Marker} route, the
#'   marker-to-position distance on a positional match, \code{NA} unplaced.
#'   A region table (\code{Chrom} + \code{Start}/\code{End}, no \code{Pos})
#'   is matched at its midpoint.
#'
#' @seealso \code{\link{attach_marker_annotations}},
#'   \code{\link{track_marker_labels}}
#' @export
standardize_marker_annotations <- function(x, meta = NULL,
                                           match_bp = .HS_ANNOT_DEFAULT_BP,
                                           chr_map = NULL) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"Name" %in% names(df))
    .hs_abort("standardize_marker_annotations: 'x' needs a Name column.")
  has_marker <- "Marker" %in% names(df)
  has_span   <- "Chrom" %in% names(df) && all(c("Start", "End") %in% names(df))
  has_pos    <- "Chrom" %in% names(df) && ("Pos" %in% names(df) || has_span)
  if (!has_marker && !has_pos)
    .hs_abort("standardize_marker_annotations: 'x' needs Marker, or Chrom + Pos ",
              "(or Chrom + Start + End).")
  if (has_pos && !"Pos" %in% names(df))      # a region: matched at its midpoint
    df$Pos <- (as.numeric(df$Start) + as.numeric(df$End)) / 2
  if (!is.numeric(match_bp) || length(match_bp) != 1L || is.na(match_bp) ||
      match_bp < 0)
    .hs_abort("standardize_marker_annotations: 'match_bp' must be a single ",
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
    .hs_check_meta(meta, "standardize_marker_annotations", "meta")
    if (has_pos) .hs_require_bp(meta, "standardize_marker_annotations")
    # One Marker -> locus lookup (D5 §3.6 c), in map order.
    loc <- .marker_locus_table(meta, keyed = FALSE)

    # Marker route first: known ids take their mapped position.
    if (has_marker) {
      idx <- match(out$Marker, loc$Marker)
      hit <- !is.na(idx)
      out$Chrom[hit] <- loc$Chrom[idx[hit]]
      out$Pos[hit]   <- loc$Pos[idx[hit]]
      out$matched_dist_bp[hit] <- 0
      out$Marker[!hit] <- NA_character_
    }
    # Positional route for everything still unresolved but located: nearest
    # marker on the (mapped) chromosome, one vectorised pass (D5 §3.6 d); a
    # matched row takes the container's chromosome spelling.
    todo <- is.na(out$Marker) & !is.na(out$Chrom) & !is.na(out$Pos)
    if (any(todo)) {
      nm <- .hs_nearest_marker(.hs_map_chrom(loc$Chrom, chr_map), loc$Pos,
                               out$Chrom[todo], out$Pos[todo], match_bp)
      hit <- !is.na(nm$index)
      w <- which(todo)[hit]
      out$Marker[w]          <- loc$Marker[nm$index[hit]]
      out$Chrom[w]           <- loc$Chrom[nm$index[hit]]
      out$matched_dist_bp[w] <- nm$dist[hit]
    }
    # Block of the resolved marker.
    bidx <- match(out$Marker, loc$Marker)
    out$Block <- ifelse(is.na(bidx), NA_character_, loc$Block[bidx])
  } else if (has_marker && !has_pos) {
    .hs_abort("standardize_marker_annotations: the Marker route needs 'meta' ",
         "to resolve positions.")
  }

  unplaced <- is.na(out$Marker)
  if (!is.null(meta) && any(unplaced))
    .hs_warn("standardize_marker_annotations: ", sum(unplaced), " of ", n,
            " row(s) match no marker within ",
            format(match_bp, big.mark = ","), " bp \u2014 kept with Marker/Block ",
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
#' @param chr_map  Passed to the standardizer: container chromosome -> the
#'   table's spelling.
#'
#' @return The \code{GenoMeta} with \code{$annotation$custom} set. The
#'   provenance under \code{metadata$annotation$custom} records the source,
#'   \code{match_bp}, \code{chr_map}, the counts, the unplaced names and an
#'   ISO-8601 \code{created} stamp.
#' @seealso \code{\link{standardize_marker_annotations}},
#'   \code{\link{track_marker_labels}}
#' @export
attach_marker_annotations <- function(meta, x,
                                      match_bp = .HS_ANNOT_DEFAULT_BP,
                                      source = NULL, chr_map = NULL) {
  .hs_check_meta(meta, "attach_marker_annotations", "meta")
  std_cols <- c("Marker", "Chrom", "Pos", "Block", "Name", "Type", "Trait",
                "matched_dist_bp")
  if (is.data.frame(x) && all(std_cols %in% names(x))) {
    # An already-standardized table: the fast path, but checked (ann/F21) —
    # every resolved Marker must be in the map and every Block in the blocks.
    std <- as.data.frame(x, stringsAsFactors = FALSE)[, std_cols]
    mk  <- as.character(std$Marker)
    bad_m <- which(!is.na(mk) & !mk %in% as.character(meta$map_info$Marker))
    if (length(bad_m))
      .hs_abort("attach_marker_annotations: ", length(bad_m), " row(s) name a ",
                "Marker that is not in the container's map (first: ",
                paste(utils::head(std$Name[bad_m], 5L), collapse = ", "),
                if (length(bad_m) > 5L) " ..." else "", ").", call. = FALSE,
                data = list(bad = bad_m))
    known <- unique(as.character(if (!is.null(meta$block_marker_map_info))
      meta$block_marker_map_info$Block else meta$block_info$Block))
    bk <- as.character(std$Block)
    bad_b <- which(!is.na(bk) & !bk %in% known)
    if (length(bad_b))
      .hs_abort("attach_marker_annotations: ", length(bad_b), " row(s) name a ",
                "Block that is not in the container's blocks (first: ",
                paste(utils::head(std$Name[bad_b], 5L), collapse = ", "),
                if (length(bad_b) > 5L) " ..." else "", ").", call. = FALSE,
                data = list(bad = bad_b))
  } else {
    std <- standardize_marker_annotations(x, meta = meta, match_bp = match_bp,
                                          chr_map = chr_map)
  }

  ann <- meta$annotation %||% list()
  ann$custom <- std
  prov <- meta$metadata$annotation %||% list()
  prov$custom <- list(
    source          = source %||% NA_character_,
    match_bp        = as.numeric(match_bp),
    chr_map         = chr_map,
    n               = nrow(std),
    n_matched       = sum(!is.na(std$Marker)),
    n_unplaced      = sum(is.na(std$Marker)),
    unplaced        = as.character(std$Name[is.na(std$Marker)]),
    created         = .hs_utc_now(),
    package_version = as.character(utils::packageVersion("haploStack")))
  # One rebuild path (0.9.0): the custom table joins the stored annotation,
  # everything else is forwarded and the whole container revalidated.
  .rebuild_geno_meta(meta, annotation = ann, metadata = list(annotation = prov),
                     .caller = "attach_marker_annotations")
}


# One nearest-marker matcher (D5 §3.6 d): per chromosome, findInterval() on
# the sorted positions, comparing the left and right neighbours. Ties keep the
# 0.8.8 choice — the first map row in map order among the markers at the
# minimum distance — so the old per-row loop gives identical results.
#   map_chrom, map_pos  the map (after chr_map), in map order
#   chrom, pos          the queries
#   max_bp              beyond it a query is unplaced (NA)
# Returns list(index = map row or NA, dist = distance or NA), one per query.
.hs_nearest_marker <- function(map_chrom, map_pos, chrom, pos, max_bp) {
  n <- length(chrom)
  index <- rep(NA_integer_, n)
  dist  <- rep(NA_real_, n)
  map_chrom <- as.character(map_chrom); map_pos <- as.numeric(map_pos)
  chrom <- as.character(chrom); pos <- as.numeric(pos)
  todo <- !is.na(chrom) & !is.na(pos)
  for (ch in unique(chrom[todo])) {
    rows <- which(map_chrom == ch & !is.na(map_pos))
    if (length(rows) == 0L) next
    # first map row at each distinct position, positions ascending
    first <- tapply(rows, map_pos[rows], min)
    pu    <- as.numeric(names(first))
    o     <- order(pu)
    pu    <- pu[o]; first <- as.integer(first[o])
    q  <- which(todo & chrom == ch)
    k  <- findInterval(pos[q], pu)
    li <- pmax(k, 1L); ri <- pmin(k + 1L, length(pu))
    dl <- abs(pos[q] - pu[li]); dr <- abs(pos[q] - pu[ri])
    dm <- pmin(dl, dr)
    cl <- ifelse(dl == dm, first[li], NA_integer_)
    cr <- ifelse(dr == dm, first[ri], NA_integer_)
    best <- pmin(cl, cr, na.rm = TRUE)
    ok <- dm <= max_bp
    index[q[ok]] <- best[ok]
    dist[q[ok]]  <- dm[ok]
  }
  list(index = index, dist = dist)
}
