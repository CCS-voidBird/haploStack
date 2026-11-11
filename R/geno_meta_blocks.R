# =============================================================================
# geno_meta_blocks.R — the block map: one resolver, one provenance record,
# one writer (0.9.0, design D2-5 … D2-7, revised on the owner's steer of
# 2026-10-29)
#
# The container knows two block-map shapes and nothing else:
#
#   block-marker map   (Chrom, Marker, Block)        blocking by markers
#   block-region map   (Chrom, Block, Start, End)    blocking by regions
#
# Both are defined here, by the functions, not by any external file format: a
# file must carry exactly these column names in its header (case does not
# matter; extra columns are ignored, except the Signal_* columns a region map
# written by build_block_location_info() carries). No alias vocabulary, no
# headerless guessing. SelectionTools' own formats are not read here:
# st_convert() (R/st_formats.R) turns them into one of the two shapes first,
# and only when use_st_formats(TRUE) has been called.
#
#   .as_block_tables()    path | data.frame | resolver list -> validated tables
#                         + the provenance record (chromosome checks, match
#                         rate, overlap policy; classed errors)
#   .block_provenance()   the metadata$blocks record
#   .write_block_files()  the two canonical files, one writer
#
# The marker-level table is the container's truth: block_info is always
# derived from the markers' positions (.compute_block_location(),
# this file since D2-9, catalogue_blockmap.R before); a region map keeps its declared span as
# Region_Start / Region_End beside it (form row 34 a).
# =============================================================================

.HS_BLOCK_SHAPES <- list(
  block_marker_map = c("Chrom", "Marker", "Block"),
  block_region_map = c("Chrom", "Block", "Start", "End"))

# The sentence every shape error ends with.
.hs_block_shape_msg <- function() {
  paste0("a block-marker map with columns Chrom, Marker, Block, or a ",
         "block-region map with columns Chrom, Block, Start, End ",
         "(header required; case does not matter)")
}

.hs_block_abort <- function(caller, ..., data = list()) {
  .hs_abort(class = c("hs_error_block_map", "hs_error_input"), data = data,
            caller, ": ", ..., call. = FALSE)
}

# Identify the shape of `tab` and rename its columns to the canonical names.
# Names decide; only a data.frame without names (NULL or V1, V2, ...) is read
# positionally, and only when its column count is exactly a shape's.
.hs_block_shape <- function(tab, layout, caller) {
  nm      <- names(tab)
  unnamed <- is.null(nm) || all(grepl("^V[0-9]+$", nm))
  lower   <- if (unnamed) character(0L) else tolower(trimws(nm))
  has <- function(shape) all(tolower(.HS_BLOCK_SHAPES[[shape]]) %in% lower)
  fits <- function(shape) has(shape) ||
    (unnamed && ncol(tab) == length(.HS_BLOCK_SHAPES[[shape]]))

  found <- if (layout != "auto") {
    if (fits(layout)) layout else NA_character_
  } else if (has("block_region_map")) "block_region_map"
    else if (has("block_marker_map")) "block_marker_map"
    else if (unnamed && ncol(tab) == 3L) "block_marker_map"
    else if (unnamed && ncol(tab) == 4L) "block_region_map"
    else NA_character_

  if (is.na(found))
    .hs_block_abort(caller, "the block map must be ",
                    if (layout == "auto") .hs_block_shape_msg()
                    else paste0("a ", sub("_", "-", sub("_map$", " map", layout)),
                                " with columns ",
                                paste(.HS_BLOCK_SHAPES[[layout]], collapse = ", ")),
                    "; got ",
                    if (unnamed) paste0(ncol(tab), " unnamed column(s)")
                    else paste0("columns: ", paste(nm, collapse = ", ")),
                    ".", data = list(columns = nm))

  cols <- .HS_BLOCK_SHAPES[[found]]
  if (unnamed) {
    names(tab)[seq_along(cols)] <- cols
  } else {
    for (cc in cols) names(tab)[match(tolower(cc), lower)] <- cc
  }
  list(layout = found, tab = tab)
}

# Markers of `map` that fall inside the regions (by chromosome and position):
# data.frame(Chrom, Marker, Block, srow, mrow) ordered by region row, then map
# row — exactly what the per-region loop of expand_block_info() returned up
# to 0.8.8 (blk/F14). Rows without a chromosome or position take no part.
.hs_expand_regions <- function(map, regions) {
  empty <- data.frame(Chrom = character(0L), Marker = character(0L),
                      Block = character(0L), srow = integer(0L),
                      mrow = integer(0L), stringsAsFactors = FALSE)
  mp <- data.table::data.table(
    Chrom  = as.character(map$Chrom),
    Marker = as.character(map$Marker),
    Pos    = suppressWarnings(as.numeric(map$Pos)),
    mrow   = seq_len(nrow(map)))
  mp <- mp[!is.na(Chrom) & !is.na(Pos)]
  rg <- data.table::data.table(
    Chrom = as.character(regions$Chrom),
    Block = as.character(regions$Block),
    Start = suppressWarnings(as.numeric(regions$Start)),
    End   = suppressWarnings(as.numeric(regions$End)),
    srow  = seq_len(nrow(regions)))
  rg <- rg[!is.na(Chrom) & !is.na(Start) & !is.na(End) & Start <= End]
  if (nrow(mp) == 0L || nrow(rg) == 0L) return(empty)
  mp[, Pos2 := Pos]
  data.table::setkey(rg, Chrom, Start, End)
  hits <- data.table::foverlaps(mp, rg, by.x = c("Chrom", "Pos", "Pos2"),
                                type = "within", nomatch = NULL)
  if (nrow(hits) == 0L) return(empty)
  data.table::setorder(hits, srow, mrow)
  out <- as.data.frame(hits[, list(Chrom, Marker, Block, srow, mrow)],
                       stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

# Suggest a chrom_map for source chromosome names absent from the panel: the
# same name once "chr" and a trailing letter group are stripped ("1H" -> "1",
# "chr1" -> "1"). Named character vector (names = source spelling), possibly
# empty.
.hs_chrom_suggest <- function(absent, panel) {
  norm <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x <- sub("^chr(om)?[_. ]?", "", x)
    sub("[a-z]+$", "", x)
  }
  pn <- norm(panel)
  an <- norm(absent)
  m  <- match(an, pn)
  ok <- !is.na(m) & nzchar(an) & !duplicated(m)
  if (!any(ok)) return(character(0L))
  stats::setNames(as.character(panel[m[ok]]), absent[ok])
}

# The overlap policy (form row 8 a): a marker in two blocks is an error by
# default; overlap = "first" keeps each marker's first block and warns.
.hs_block_overlap <- function(bmm, overlap, caller) {
  dup  <- duplicated(bmm$Marker)
  bad  <- unique(bmm$Marker[dup])
  n    <- length(bad)
  if (n == 0L) return(list(bmm = bmm, n_overlap = 0L))
  desc <- vapply(utils::head(bad, 5L), function(m)
    paste0(m, ": ", paste(unique(bmm$Block[bmm$Marker == m]), collapse = ", ")),
    character(1L))
  shown <- paste0(paste(desc, collapse = "; "), if (n > 5L) "; ..." else "")
  if (overlap == "error")
    .hs_block_abort(caller, n, " marker(s) fall in more than one block (",
                    shown, "). Blocks must not overlap: fix the map, or pass ",
                    "overlap = \"first\" to keep each marker's first block.",
                    data = list(bad = bad))
  .hs_warn(caller, ": ", n, " marker(s) fall in more than one block; the ",
           "first block of each is kept (overlap = \"first\"): ", shown,
           call. = FALSE)
  list(bmm = bmm[!dup, , drop = FALSE], n_overlap = as.integer(n))
}

# The identity of a block set: md5 of its (Marker, Block) pairs in Marker
# order. The same map from a file, a data.frame or a round-tripped CSV gets
# the same value.
.hs_map_md5 <- function(bmm) {
  mk <- as.character(bmm$Marker)
  bl <- as.character(bmm$Block)
  o  <- .hs_order(mk, bl)
  tf <- tempfile("hs_bmm_")
  on.exit(unlink(tf), add = TRUE)
  writeLines(paste(mk[o], bl[o], sep = "\t"), tf)
  unname(tools::md5sum(tf))
}

# The metadata$blocks record (design D2 §3.5). `bmm` is the marker-level table
# about to be attached; `map` the panel's map (every marker of the container),
# so n_singletons counts the panel markers the map does not cover — they
# become singleton dummy blocks in the constructor.
.block_provenance <- function(bmm, map, source, label, params = list(),
                              file = NA_character_, file_md5 = NA_character_,
                              n_unmatched = 0L) {
  covered <- unique(as.character(bmm$Marker))
  list(source       = as.character(source),
       label        = as.character(label),
       params       = params,
       file         = as.character(file),
       file_md5     = as.character(file_md5),
       map_md5      = .hs_map_md5(bmm),
       n_blocks     = length(unique(as.character(bmm$Block))),
       n_singletons = as.integer(length(setdiff(as.character(map$Marker), covered))),
       n_matched    = length(covered),
       n_unmatched  = as.integer(n_unmatched))
}


# ---------------------------------------------------------------------------
# The resolver
# ---------------------------------------------------------------------------

# Turn a block-map source into the two validated tables and the provenance
# record, against the panel map `map_df` (Marker, Chrom, Pos).
#
#   source     a path (text file with a header), a data.frame in either shape,
#              or the list an LD resolver returns (block_marker_map [+
#              block_info] [+ provenance]).
#   chrom_map  named character: names = chromosome spelling in the source,
#              values = the panel's. Applied before any check.
#   overlap    "error" (default) or "first" — .hs_block_overlap().
#   layout     "auto" (default) or a shape name to force.
#
# Returns list(block_info, block_marker_map, provenance):
#   block_marker_map  data.frame(Chrom, Marker, Block) — Chrom is the PANEL's
#                     spelling (the join is on Marker); only markers of the
#                     panel; one row per marker.
#   block_info        the 7-column marker-derived table (Chrom, Block, Start,
#                     End, Pos, Num_SNPs, Block_length), plus Region_Start /
#                     Region_End for a region map and any Signal_* columns a
#                     region source carried.
#
# Policies (form rows 7 a, 8 a):
#   * shape by canonical names only (.hs_block_shape());
#   * identifiers stay text; Start/End must parse as numbers;
#   * marker map: joined on Marker. A consistent chromosome rename ("1H" -> "1",
#     one-to-one) is reported once and recorded; an inconsistent one, or a
#     block whose markers sit on two panel chromosomes, is an error;
#   * region map: a chromosome absent from the panel is an error naming it
#     and, when a spelling match exists, the chrom_map to pass;
#   * zero matched markers is always an error; under 50 % a warning;
#   * exact duplicate rows are dropped silently; a marker in two blocks goes
#     to the overlap policy.
.as_block_tables <- function(source, map_df, chrom_map = NULL,
                             overlap = c("error", "first"),
                             layout = c("auto", "block_marker_map",
                                        "block_region_map"),
                             caller = "attach_block_info") {
  overlap <- match.arg(overlap)
  layout  <- match.arg(layout)

  # ---- 1. read ---------------------------------------------------------------
  base_prov <- NULL
  file      <- NA_character_
  file_md5  <- NA_character_
  label     <- "user_map"
  if (is.list(source) && !is.data.frame(source) &&
      (!is.null(source$block_marker_map) || !is.null(source$block_marker_map_info))) {
    tab <- source$block_marker_map %||% source$block_marker_map_info
    if (!is.data.frame(tab))
      .hs_block_abort(caller, "the list's block_marker_map must be a data.frame.")
    tab       <- as.data.frame(tab, stringsAsFactors = FALSE)
    base_prov <- source$provenance
    layout    <- "block_marker_map"
  } else if (is.character(source) && length(source) == 1L) {
    if (!file.exists(source))
      .hs_abort(caller, ": file not found: ", source, call. = FALSE)
    # Text in, text out: "007" stays "007", "1e5" stays "1e5" (blk/F1).
    tab      <- .read_table_chr(source, header = TRUE)
    file     <- basename(source)
    file_md5 <- unname(tools::md5sum(source))
    label    <- sub("[.][^.]*$", "", file)
  } else if (is.data.frame(source)) {
    tab <- as.data.frame(source, stringsAsFactors = FALSE)
  } else {
    .hs_block_abort(caller, "'block_source' must be a data.frame, a file path, ",
                    "or the list an LD resolver returns; got ", class(source)[1L],
                    ".")
  }
  if (nrow(tab) == 0L)
    .hs_block_abort(caller, "the block map has no rows.")

  # ---- 2. shape and types ----------------------------------------------------
  sh     <- .hs_block_shape(tab, layout, caller)
  layout <- sh$layout
  tab    <- sh$tab
  cols   <- .HS_BLOCK_SHAPES[[layout]]
  for (cc in intersect(c("Chrom", "Marker", "Block"), cols))
    tab[[cc]] <- trimws(as.character(tab[[cc]]))
  need <- if (layout == "block_marker_map") c("Marker", "Block")
          else c("Chrom", "Block")
  blank <- Reduce(`|`, lapply(tab[need], function(v) is.na(v) | !nzchar(v)))
  if (any(blank))
    .hs_block_abort(caller, sum(blank), " row(s) have an empty ",
                    paste(need, collapse = " or "), " (first: row ",
                    which(blank)[1L], ").", data = list(bad = which(blank)))
  if (layout == "block_region_map") {
    for (cc in c("Start", "End")) {
      raw <- tab[[cc]]
      v   <- suppressWarnings(as.numeric(raw))
      badv <- which(is.na(v) & !is.na(raw))
      if (length(badv))
        .hs_block_abort(caller, cc, " is not numeric in ", length(badv),
                        " row(s) (first: row ", badv[1L], ", '", raw[badv[1L]],
                        "').", data = list(bad = badv))
      if (anyNA(v))
        .hs_block_abort(caller, cc, " is missing in ", sum(is.na(v)), " row(s).",
                        data = list(bad = which(is.na(v))))
      tab[[cc]] <- v
    }
    rev <- which(tab$Start > tab$End)
    if (length(rev))
      .hs_block_abort(caller, length(rev), " row(s) have Start > End (first: ",
                      "block '", tab$Block[rev[1L]], "', ", tab$Chrom[rev[1L]],
                      ":", format(tab$Start[rev[1L]], scientific = FALSE),
                      "-", format(tab$End[rev[1L]], scientific = FALSE), ").",
                      data = list(bad = rev))
  }
  tab <- tab[!duplicated(tab[cols]), , drop = FALSE]
  rownames(tab) <- NULL

  # ---- 3. the panel map ------------------------------------------------------
  map <- as.data.frame(map_df, stringsAsFactors = FALSE)
  map$Chrom  <- as.character(map$Chrom)
  map$Marker <- as.character(map$Marker)
  map$Pos    <- suppressWarnings(as.numeric(map$Pos))
  panel_chroms <- unique(map$Chrom[!is.na(map$Chrom) & nzchar(map$Chrom)])

  # ---- 4. chromosome spelling ------------------------------------------------
  if (!is.null(chrom_map)) {
    if (!is.character(chrom_map) || is.null(names(chrom_map)) ||
        any(!nzchar(names(chrom_map))))
      .hs_abort(caller, ": 'chrom_map' must be a named character vector ",
                "(names = chromosome spelling in the block map, values = the ",
                "panel's).", call. = FALSE)
    hit <- tab$Chrom %in% names(chrom_map)
    tab$Chrom[hit] <- unname(chrom_map[tab$Chrom[hit]])
  }
  chrom_renamed <- NULL
  n_unmatched   <- 0L

  # ---- 5. match --------------------------------------------------------------
  if (layout == "block_region_map") {
    absent <- setdiff(unique(tab$Chrom), panel_chroms)
    if (length(absent)) {
      sug <- .hs_chrom_suggest(absent, panel_chroms)
      .hs_block_abort(
        caller, length(absent), " chromosome(s) of the block map are not in ",
        "this panel: ", paste(utils::head(absent, 5L), collapse = ", "),
        if (length(absent) > 5L) " ..." else "",
        ". The panel has: ", paste(utils::head(panel_chroms, 8L), collapse = ", "),
        if (length(panel_chroms) > 8L) " ..." else "", ".",
        if (length(sug)) paste0(" Try chrom_map = c(",
                                paste0('"', names(sug), '" = "', sug, '"',
                                       collapse = ", "), ").") else "",
        data = list(bad = absent, chrom_map = sug))
    }
    hits  <- .hs_expand_regions(map, tab)
    empty <- setdiff(seq_len(nrow(tab)), unique(hits$srow))
    n_unmatched <- length(empty)
    if (nrow(hits) == 0L)
      .hs_block_abort(caller, "no marker of this panel falls in any of the ",
                      nrow(tab), " block region(s): wrong file, or positions ",
                      "in another unit?")
    if (length(empty)) {
      shown <- sprintf("%s (%s:%s-%s)", tab$Block[empty], tab$Chrom[empty],
                       format(tab$Start[empty], scientific = FALSE, trim = TRUE),
                       format(tab$End[empty], scientific = FALSE, trim = TRUE))
      .hs_warn(caller, ": ", length(empty), " block region(s) contain no marker ",
               "of this panel and were omitted: ",
               paste(utils::head(shown, 5L), collapse = ", "),
               if (length(empty) > 5L) " ..." else "", call. = FALSE)
    }
    bmm <- hits[, c("Chrom", "Marker", "Block"), drop = FALSE]
  } else {
    idx         <- match(tab$Marker, map$Marker)
    n_unmatched <- sum(is.na(idx))
    n_matched   <- sum(!is.na(idx))
    if (n_matched == 0L)
      .hs_block_abort(caller, "no marker of the block map is in this panel (0 of ",
                      nrow(tab), "): wrong file, or marker IDs spelled ",
                      "differently? First markers of the map: ",
                      paste(utils::head(tab$Marker, 3L), collapse = ", "),
                      "; of the panel: ",
                      paste(utils::head(map$Marker, 3L), collapse = ", "), ".")
    if (n_matched < 0.5 * nrow(tab))
      .hs_warn(caller, ": only ", n_matched, " of ", nrow(tab), " marker(s) of ",
               "the block map are in this panel (",
               round(100 * n_matched / nrow(tab)), "%); the rest are ignored.",
               call. = FALSE)
    keep <- !is.na(idx)
    tab  <- tab[keep, , drop = FALSE]
    idx  <- idx[keep]
    src  <- tab$Chrom
    dst  <- map$Chrom[idx]
    known <- !is.na(src) & nzchar(src) & !is.na(dst) & nzchar(dst)
    if (any(known & src != dst)) {
      pairs <- unique(data.frame(src = src[known], dst = dst[known],
                                 stringsAsFactors = FALSE))
      many_dst <- pairs$src[duplicated(pairs$src)]
      many_src <- pairs$dst[duplicated(pairs$dst)]
      if (length(many_dst) || length(many_src)) {
        show <- c(
          vapply(utils::head(unique(many_dst), 3L), function(s)
            paste0("'", s, "' -> ", paste(pairs$dst[pairs$src == s], collapse = ", ")),
            character(1L)),
          vapply(utils::head(unique(many_src), 3L), function(d)
            paste0(paste0("'", pairs$src[pairs$dst == d], "'", collapse = ", "),
                   " -> ", d), character(1L)))
        .hs_block_abort(
          caller, "the block map's chromosome names do not match this panel's ",
          "one-to-one (joined on Marker): ", paste(show, collapse = "; "),
          ". Fix the map's Chrom column, or pass chrom_map = c(<map ",
          "spelling> = <panel spelling>).", data = list(bad = pairs))
      }
      ren <- pairs[pairs$src != pairs$dst, , drop = FALSE]
      chrom_renamed <- stats::setNames(ren$dst, ren$src)
      .hs_inform(caller, ": chromosome names of the block map were matched to ",
                 "the panel's (joined on Marker): ",
                 paste0(ren$src, " -> ", ren$dst, collapse = ", "), ".")
    }
    per_block <- tapply(dst, tab$Block, function(x) length(unique(x[!is.na(x)])))
    split_bl  <- names(per_block)[per_block > 1L]
    if (length(split_bl))
      .hs_block_abort(
        caller, length(split_bl), " block(s) have markers on more than one ",
        "chromosome of this panel: ",
        paste(vapply(utils::head(split_bl, 5L), function(b)
          paste0(b, " (", paste(unique(dst[tab$Block == b]), collapse = ", "), ")"),
          character(1L)), collapse = "; "),
        if (length(split_bl) > 5L) "; ..." else "", ".",
        data = list(bad = split_bl))
    bmm <- data.frame(Chrom = dst, Marker = tab$Marker, Block = tab$Block,
                      stringsAsFactors = FALSE)
  }

  # ---- 6. overlap ------------------------------------------------------------
  ov  <- .hs_block_overlap(bmm, overlap, caller)
  bmm <- ov$bmm
  rownames(bmm) <- NULL

  # ---- 7. block_info ---------------------------------------------------------
  bi <- as.data.frame(.compute_block_location(bmm, map)$block_info,
                      stringsAsFactors = FALSE)
  if (layout == "block_region_map") {
    # A block declared on several rows spans the union of its regions.
    blocks <- unique(as.character(tab$Block))
    reg <- data.frame(
      Block        = blocks,
      Region_Start = vapply(split(tab$Start, tab$Block), min, numeric(1L))[blocks],
      Region_End   = vapply(split(tab$End,   tab$Block), max, numeric(1L))[blocks],
      stringsAsFactors = FALSE)
    bi$Region_Start <- reg$Region_Start[match(bi$Block, reg$Block)]
    bi$Region_End   <- reg$Region_End[match(bi$Block, reg$Block)]
    sig <- intersect(.signal_marker_cols(), names(tab))
    if (length(sig)) {
      first <- tab[!duplicated(tab$Block), c("Block", sig), drop = FALSE]
      for (cc in sig) bi[[cc]] <- first[[cc]][match(bi$Block, first$Block)]
    }
  }
  rownames(bi) <- NULL

  # ---- 8. provenance ---------------------------------------------------------
  params <- if (!is.null(base_prov)) base_prov$params %||% list()
            else list(layout = layout, overlap = overlap)
  if (!is.null(chrom_map))     params$chrom_map     <- chrom_map
  if (!is.null(chrom_renamed)) params$chrom_renamed <- chrom_renamed
  if (ov$n_overlap > 0L)       params$n_overlap     <- ov$n_overlap
  prov <- .block_provenance(
    bmm, map,
    source      = base_prov$source   %||% "user_map",
    label       = base_prov$label    %||% label,
    params      = params,
    file        = base_prov$file     %||% file,
    file_md5    = base_prov$file_md5 %||% file_md5,
    n_unmatched = if (!is.null(base_prov$n_unmatched)) base_prov$n_unmatched
                  else n_unmatched)

  list(block_info = bi, block_marker_map = bmm, provenance = prov)
}


# ---------------------------------------------------------------------------
# One writer (form row 9 a)
# ---------------------------------------------------------------------------

# Write the two canonical block files:
#   <prefix>_block_info.csv        Chrom, Block, Start, End, Pos, Num_SNPs,
#                                  Block_length [, Region_Start, Region_End]
#                                  [, Signal_*]
#   <prefix>_block_marker_map.csv  Chrom, Marker, Block
# Positions are written as plain integers (scipen), never "8.8e+07". The
# marker-map file is the exact round trip; the coordinate file re-expands by
# position. `tables` is list(block_info, block_marker_map) (a resolver's
# result, or the container's two tables under those names). Returns the two
# paths invisibly, after one message.
.write_block_files <- function(tables, output_dir, prefix) {
  if (!dir.exists(output_dir))
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  bi  <- data.table::as.data.table(tables$block_info)
  bmm <- data.table::as.data.table(tables$block_marker_map)
  keep <- c(intersect(c("Chrom", "Block", "Start", "End", "Pos", "Num_SNPs",
                        "Block_length", "Region_Start", "Region_End"), names(bi)),
            intersect(.signal_marker_cols(), names(bi)))
  f_bi  <- file.path(output_dir, paste0(prefix, "_block_info.csv"))
  f_bmm <- file.path(output_dir, paste0(prefix, "_block_marker_map.csv"))
  data.table::fwrite(bi[, keep, with = FALSE], f_bi, sep = ",", quote = FALSE,
                     scipen = 100L, showProgress = FALSE)
  data.table::fwrite(bmm[, c("Chrom", "Marker", "Block"), with = FALSE], f_bmm,
                     sep = ",", quote = FALSE, scipen = 100L, showProgress = FALSE)
  .hs_inform("Block files written: ", f_bi, " and ", f_bmm)
  invisible(c(block_info = f_bi, block_marker_map = f_bmm))
}


# ---- moved here in 0.9.0 (D2-9) from catalogue_blockmap.R, format_helper.R ----

#' Attach or update block information in a GenoMeta object
#'
#' Updates an existing \code{GenoMeta} with block assignment data and returns
#' a new \code{GenoMeta} with both \code{block_info} and
#' \code{block_marker_map_info} populated. The original object is not modified.
#' Any previously stored block data is replaced.
#'
#' A block map you already have is the \strong{recommended} way to give
#' haploStack its blocks, and it needs no LD parameter. Pass it here, or at read
#' time: \code{\link{read_plink_ped}(block_file =)} or
#' \code{\link{build_geno_meta}(blocks =)}, which also keeps the FAM. Markers
#' absent from the map become singleton blocks. Only when no map exists do the
#' local LD resolvers come in: \code{\link{hapS_ld_blocks}} (HapSelect, the
#' default) or \code{\link{st_ld_blocks}} (SelectionTools, experimental).
#'
#' This function rebuilds the container and (since 0.9.0) keeps
#' \code{marker_effects}, \code{phenotype} and the FAM. What a re-blocking
#' invalidates is dropped, with a message: a stored signal-marker selection
#' and the gene annotation (custom marker annotations are kept with their
#' \code{Block} re-derived).
#'
#' @section The two shapes:
#' The block map is defined by this function, not by any file format. Pass
#' one of:
#' \describe{
#'   \item{block-marker map}{columns \code{Chrom}, \code{Marker},
#'     \code{Block}: one row per marker (blocking by markers). The join is on
#'     \code{Marker}; \code{Chrom} is checked against the panel. Markers absent
#'     from the panel are ignored (a warning below 50 \% matched, an error at
#'     none); panel markers absent from the map become singleton blocks.}
#'   \item{block-region map}{columns \code{Chrom}, \code{Block}, \code{Start},
#'     \code{End}: one row per region (blocking by regions). Every panel marker
#'     whose position lies in a region takes its block; a region with no marker
#'     is omitted with a warning. The declared span is kept as
#'     \code{Region_Start} / \code{Region_End} in \code{block_info}, whose
#'     \code{Start} / \code{End} are always marker-derived.}
#' }
#' A file must carry these column names in its header (case does not matter;
#' extra columns are ignored); a data.frame is read by name, or positionally
#' when it has no names and exactly three or four columns. Identifiers are kept
#' as text. Anything else stops with an \code{hs_error_block_map} naming the
#' two shapes. SelectionTools' own block table is converted first by
#' \code{\link{st_convert}}.
#'
#' @section Chromosomes and overlap:
#' A consistent rename between the map and the panel (\code{1H} vs \code{1},
#' one-to-one) is resolved on \code{Marker}, reported once and recorded; an
#' inconsistent one, a block whose markers sit on two panel chromosomes, or a
#' region on a chromosome the panel does not have is an error (the message
#' suggests a \code{chrom_map} when a spelling match exists). Exact duplicate
#' rows are dropped silently; a marker in two blocks is an error unless
#' \code{overlap = "first"}, which keeps each marker's first block and warns.
#'
#' @section Provenance:
#' The container records where its blocks came from in
#' \code{metadata$blocks}: \code{source} (\code{"user_map"} or the LD
#' resolver's name), \code{label}, the file name and md5, \code{map_md5} (the
#' identity of the block set: the same for a file, a data.frame or a
#' round-tripped CSV), \code{n_blocks}, \code{n_singletons}, \code{n_matched},
#' \code{n_unmatched} and the \code{params} applied (layout, overlap,
#' \code{chrom_map}, renames). \code{\link{filter_geno_meta}} adds
#' \code{pruned_from} when blocks are lost.
#'
#' @param geno_meta    A valid \code{GenoMeta} list.
#' @param block_source A data.frame in one of the two shapes, a file path, or
#'   the list an LD resolver returns with \code{attach = FALSE}.
#' @param source_type  Character. \code{"auto"} (default) reads the shape from
#'   the column names; \code{"block_marker_map"} or \code{"block_region_map"}
#'   (alias \code{"block_info"}) forces it; \code{"file"} requires a path.
#' @param chrom_map    Optional named character vector renaming chromosomes
#'   before any check: names are the spelling in the block map, values the
#'   panel's, e.g. \code{c("1H" = "1")}.
#' @param overlap      \code{"error"} (default) or \code{"first"}; see above.
#'
#' @return A new \code{GenoMeta} with \code{block_info} and
#'   \code{block_marker_map_info} populated.
#'
#' @examples
#' \dontrun{
#' # Typical workflow: read data, then attach blocks from a marker-level table
#' meta  <- read_st_genotype("geno.csv")
#' bmap  <- read.csv("block_marker_map.csv")           # Chrom, Marker, Block
#' meta  <- attach_block_info(meta, bmap)
#'
#' # From block coordinates (Start / End); markers are matched by position
#' binfo <- data.frame(Chrom = "1", Block = "B1", Start = 1e5, End = 5e5)
#' meta  <- attach_block_info(meta, binfo, source_type = "block_info")
#'
#' # No block map? Build LD blocks locally (a side route) — these attach the
#' # blocks themselves and keep effects, phenotype and FAM:
#' meta <- hapS_ld_blocks(meta, ld_threshold = 0.5, tolerance = 3)  # HapSelect
#' }
#'
#' @seealso \code{\link{expand_block_info}}, \code{\link{filter_geno_meta}},
#'   \code{\link{new_geno_meta}}; when no block map exists,
#'   \code{\link{hapS_ld_blocks}} (default) or \code{\link{st_ld_blocks}}
#'   (experimental)
#'
#' @export
attach_block_info <- function(geno_meta,
                               block_source,
                               source_type = c("auto", "block_marker_map",
                                               "block_region_map", "block_info",
                                               "file"),
                               chrom_map   = NULL,
                               overlap     = c("error", "first")) {
  .hs_check_meta(geno_meta, "attach_block_info", "geno_meta")
  source_type <- match.arg(source_type)
  overlap     <- match.arg(overlap)
  if (source_type == "file" &&
      !(is.character(block_source) && length(block_source) == 1L))
    .hs_abort("attach_block_info: source_type = \"file\" needs a file path.")
  layout <- switch(source_type,
                   auto = , file = "auto",
                   block_marker_map = "block_marker_map",
                   block_region_map = , block_info = "block_region_map")
  tabs <- .as_block_tables(block_source, geno_meta$map_info,
                           chrom_map = chrom_map, overlap = overlap,
                           layout = layout, caller = "attach_block_info")
  .replace_block_in_meta(geno_meta, tabs$block_marker_map, tabs$block_info,
                         tabs$provenance)
}


# Internal: replace block fields and return updated GenoMeta copy
# Rebuilds via new_geno_meta() so the same block filtering, singleton dummy-block
# creation, and block ordering are applied as in any other GenoMeta construction.
#
# Note on signal markers: the incoming block_info is recomputed from marker
# positions, so any Signal_Marker/Signal_Score/N_Signal columns previously
# written by attach_signal_markers() are lost here. That is correct —
# redefining the blocks invalidates a selection made against the old ones — but
# it must not be silent, hence the message.
.replace_block_in_meta <- function(geno_meta, block_marker_map, block_info,
                                   provenance = NULL) {
  old_sig <- intersect(.signal_marker_cols(), names(geno_meta$block_info %||%
                                                      data.frame()))
  new_sig <- intersect(.signal_marker_cols(), names(block_info))
  if (length(old_sig) > 0L && length(new_sig) == 0L)
    .hs_inform("attach_block_info: block definitions replaced \u2014 the stored ",
            "signal-marker selection was dropped. Re-run ",
            "score_marker_signal()/select_signal_markers() against the new ",
            "blocks.")

  # One rebuild path (0.9.0): FAM, effects, phenotype and the sample-scoped
  # metadata are forwarded. Replacing the blocks clears what was defined
  # against the old ones: the signal selection, the gene annotation (with a
  # message; custom marker annotations are kept with Block re-derived) and
  # writes the new block provenance (metadata$blocks, D2-6).
  .rebuild_geno_meta(
    geno_meta,
    block_info            = as.data.frame(block_info, stringsAsFactors = FALSE),
    block_marker_map_info = as.data.frame(block_marker_map,
                                          stringsAsFactors = FALSE),
    metadata              = list(blocks = provenance),
    .caller = "attach_block_info")
}


# ---------------------------------------------------------------------------
# Export (3D → 2D)
# ---------------------------------------------------------------------------


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

  merged <- .left_join_dt(bm, map, by = c("Chrom", "Marker"))
  merged <- merged[!is.na(Pos)]

  # One row per (Chrom, Block), sorted by both in C-locale order with NA last,
  # as dplyr >= 1.1 group_by() + summarise() returned up to 0.8.8.
  block_info <- merged[, .(
    Start    = min(as.numeric(Pos)),
    End      = max(as.numeric(Pos)),
    Pos      = round((min(as.numeric(Pos)) + max(as.numeric(Pos))) / 2),
    Num_SNPs = uniqueN(Marker)
  ), by = .(Chrom, Block)]
  # 0.9.0 (row 29 a, D2-8): the one natural chromosome order, then block id
  block_info[, .chr_key := .hs_chrom_rank(Chrom)]
  setorderv(block_info, c(".chr_key", "Block"), na.last = TRUE)
  block_info[, .chr_key := NULL]
  block_info[, Block_length := End - Start + 1L]

  block_marker_map <- unique(merged[, .(Chrom, Marker, Block)])

  list(block_info = block_info, block_marker_map = block_marker_map)
}


# Left join in data.table with dplyr::left_join() semantics: every row of `x`
# in its own order (a row with several matches in `y` repeats, in `y` order),
# the columns of `x` first, then the non-key columns of `y`, and a non-key
# name present in both gets the ".x" / ".y" suffixes. Neither input is
# modified.
.left_join_dt <- function(x, y, by) {
  x <- copy(as.data.table(x))
  y <- copy(as.data.table(y))
  shared <- setdiff(intersect(names(x), names(y)), by)
  if (length(shared)) {
    setnames(x, shared, paste0(shared, ".x"))
    setnames(y, shared, paste0(shared, ".y"))
  }
  out <- y[x, on = by, allow.cartesian = TRUE]
  setcolorder(out, c(names(x), setdiff(names(y), by)))
  out[]
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
#'     \code{Chrom, Block, Start, End, Pos, Num_SNPs, Block_length} (plus
#'     \code{Region_Start, Region_End} after a block-region map, and the
#'     signal-marker columns when the GenoMeta carries a selection); positions
#'     as plain integers. One writer serves every block-file output of the
#'     package (0.9.0).}
#'   \item{\code{<prefix>_block_marker_map.csv}}{Canonical marker-level table:
#'     \code{Chrom, Marker, Block}.}
#' }
#'
#' The coordinates are recomputed from marker positions, so any signal-marker
#' selection stored by \code{\link{attach_signal_markers}}
#' (\code{Signal_Marker}, \code{Signal_Score}, \code{N_Signal}) and the declared
#' span of a block-region map (\code{Region_Start}, \code{Region_End}) are joined
#' back on \code{Block} rather than being silently dropped.
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
  .hs_check_meta(geno_meta, "build_block_location_info", "geno_meta")
  if (is.null(geno_meta$block_marker_map_info))
    .hs_abort("build_block_location_info: GenoMeta has no block_marker_map_info. ",
         "Call attach_block_info() first to assign markers to blocks.")

  result <- .compute_block_location(geno_meta$block_marker_map_info,
                                    geno_meta$map_info)

  # The coordinates above are recomputed from marker positions, so a stored
  # signal-marker selection (attach_signal_markers()) and the declared span of
  # a block-region map (Region_Start / Region_End) would be lost; join them
  # back on Block so both the returned table and the CSV carry them.
  sig_cols <- intersect(c("Region_Start", "Region_End", .signal_marker_cols()),
                        names(geno_meta$block_info))
  if (length(sig_cols) > 0L) {
    sig <- unique(as.data.table(geno_meta$block_info)[
      , c("Block", sig_cols), with = FALSE])
    sig[, Block := as.character(Block)]
    result$block_info[sig, on = "Block",
                      (sig_cols) := mget(paste0("i.", sig_cols))]
  }

  if (!is.null(output_dir)) {
    f <- .write_block_files(result, output_dir, prefix)
    nm <- names(f)
    f <- normalizePath(f, winslash = "/", mustWork = TRUE)   # drops the names
    attr(result, "hs_files") <- .hs_manifest(f, "table", nm, file.info(f)$size)
  }

  result
}


#' Convert external block coordinates to a marker-level block map
#'
#' Takes an external block coordinate table \code{(Chrom, Block, Start, End)}
#' and finds all markers in the \code{GenoMeta}'s \code{map_info} whose
#' position falls within each block's genomic range. Returns the canonical
#' internal \code{(Chrom, Marker, Block)} table. Since 0.9.0 this is the engine
#' \code{\link{attach_block_info}} runs for a block-region map
#' (\code{data.table::foverlaps}); rows come out in the region's order, then
#' the map's, as the loop up to 0.8.8 returned them.
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
#' @param map_info       The marker map (\code{Chrom}, \code{Marker},
#'   \code{Pos}), typically \code{geno_meta$map_info}; markers are matched to
#'   the block ranges by chromosome and position.
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
#' bmap <- expand_block_info(meta$map_info, block_coords)
#' meta <- attach_block_info(meta, bmap)
#' }
#'
#' @seealso \code{\link{attach_block_info}}, \code{\link{build_block_location_info}}
#'
#' @export
expand_block_info <- function(map_info, block_info_df) {
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

  hits  <- .hs_expand_regions(map, bi)
  empty <- setdiff(seq_len(nrow(bi)), unique(hits$srow))
  for (i in empty)
    .hs_warn("expand_block_info: block '", bi$Block[i],
             "' (", bi$Chrom[i], ":", bi$Start[i], "-", bi$End[i],
             ") contains no markers \u2014 omitted from output.")

  if (nrow(hits) == 0L) {
    .hs_warn("expand_block_info: no markers matched any block.")
    return(data.frame(Chrom  = character(0L), Marker = character(0L),
                      Block  = character(0L), stringsAsFactors = FALSE))
  }
  out <- hits[, c("Chrom", "Marker", "Block"), drop = FALSE]
  rownames(out) <- NULL

  dup_markers <- out$Marker[duplicated(out$Marker)]
  if (length(dup_markers) > 0L)
    .hs_warn("expand_block_info: ", length(unique(dup_markers)),
            " marker(s) fall in more than one block: ",
            paste(head(unique(dup_markers), 5L), collapse = ", "),
            if (length(unique(dup_markers)) > 5L) " ..." else "")

  out[, c("Chrom", "Marker", "Block"), drop = FALSE]
}
