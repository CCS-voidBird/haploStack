# =============================================================================
# catalogue_profile.R
#
# Profile views derived on demand from a GenoMeta object, plus GenoMeta
# serialisation and allele recoding.
#
# GenoMeta (see new_geno_meta) is the single central data structure. The
# denormalised "profile views" used by the analysis functions are built on
# demand by geno_meta_views(); they are never stored as a separate object.
#
# Profile view structure (list returned by geno_meta_views):
#   $alleles      data.table — wide: Sample, Haplotype_idx, [marker cols]
#   $marker_meta  data.table — marker, chromosome, position, block_id
#   $block_meta   data.table — block_id, chromosome, start_pos, end_pos, n_markers
#   $samples      character  — sorted unique sample names
#   $encoding     character  — "selectiontools" | "nucleotide" | "numeric"
# =============================================================================


#' Derive on-demand profile views from a GenoMeta object
#'
#' Builds the denormalised \emph{profile views} from a
#' \code{\link{new_geno_meta}} container: the wide allele table, marker-level
#' metadata, block-level metadata, the sample vector and the resolved encoding.
#' The \code{GenoMeta} object remains the single source of truth — these views
#' are computed on demand and never stored. The analysis functions
#' (\code{\link{compute_profile_localgebv}}, \code{\link{plot_haplotype_alignment}},
#' \code{\link{trace_haplotype_family}}) call this internally; users can call it
#' directly to inspect the derived tables.
#'
#' For diploid encodings (SelectionTools or nucleotide) each sample contributes
#' two rows to \code{alleles} (one per haplotype). For dosage encoding
#' (\code{"numeric"}) each sample has one row. When
#' \code{geno_meta$block_marker_map_info} is \code{NULL}, \code{block_meta} is
#' empty.
#'
#' @param geno_meta A valid \code{GenoMeta} list.
#' @param encoding  Optional encoding override; otherwise
#'   \code{metadata$encoding} is used, falling back to
#'   \code{\link{detect_genotype_encoding}}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{alleles}}{Wide \code{data.table}: Sample, Haplotype_idx,
#'     one column per marker.}
#'   \item{\code{marker_meta}}{data.table: marker, chromosome, position,
#'     block_id.}
#'   \item{\code{block_meta}}{data.table: block_id, chromosome, start_pos,
#'     end_pos, n_markers.}
#'   \item{\code{samples}}{Sorted character vector of sample IDs.}
#'   \item{\code{encoding}}{Resolved encoding string.}
#' }
#'
#' @seealso \code{\link{new_geno_meta}}, \code{\link{filter_geno_meta}},
#'   \code{\link{compute_profile_localgebv}}
#' @export
geno_meta_views <- function(geno_meta, encoding = NULL) {

  if (!is_geno_meta(geno_meta))
    stop("geno_meta_views: 'geno_meta' must be a valid GenoMeta. ",
         "Use new_geno_meta() or a read_*() function to create one, then ",
         "attach_block_info() to add block assignments.")

  arr          <- geno_meta$geno_arrays
  sample_names <- dimnames(arr)[[1L]]
  marker_names <- dimnames(arr)[[2L]]

  # Reconstruct 2-D character matrix for encoding detection and allele builders
  geno_mat <- .meta_to_allele_matrix(geno_meta)     # [n_samples × n_markers]

  # Resolve encoding: explicit arg > metadata$encoding > auto-detect
  if (is.null(encoding)) encoding <- geno_meta$metadata$encoding
  if (is.null(encoding)) {
    encoding <- detect_genotype_encoding(geno_mat)
  } else {
    encoding <- match.arg(encoding, c("selectiontools", "nucleotide", "numeric"))
  }

  geno_t <- geno_mat   # [n_samples × n_markers] (builders expect samples-as-rows)

  # Build wide alleles table using existing encoding-specific builders
  alleles_wide <- switch(
    encoding,
    selectiontools = .build_alleles_st(geno_t, sample_names, marker_names),
    nucleotide     = .build_alleles_nuc(geno_t, sample_names, marker_names),
    numeric        = .build_alleles_num(geno_t, sample_names, marker_names),
    stop("geno_meta_views: unsupported encoding '", encoding, "'.")
  )

  # map_dt — already normalised inside GenoMeta
  map_dt <- as.data.table(geno_meta$map_info)

  # marker_meta — no normalisation needed; GenoMeta guarantees Chrom/Marker/Pos
  marker_meta <- map_dt[Marker %in% marker_names,
                        .(marker     = Marker,
                          chromosome = as.character(Chrom),
                          position   = as.numeric(Pos))]

  # Block assignments — NULL is valid (block_meta will be empty)
  bmm <- geno_meta$block_marker_map_info
  if (!is.null(bmm)) {
    blk_lookup  <- as.data.table(bmm)[Marker %in% marker_names, .(Marker, Block)]
    marker_meta <- merge(marker_meta, blk_lookup,
                         by.x = "marker", by.y = "Marker", all.x = TRUE)
    setnames(marker_meta, "Block", "block_id")
  } else {
    marker_meta[, block_id := NA_character_]
  }

  marker_meta <- marker_meta[order(chromosome, position)]

  # block_meta — guard min/max against data.table's empty-group dry run
  block_meta <- marker_meta[!is.na(block_id),
                             .(chromosome = chromosome[1L],
                               start_pos  = if (.N) min(position) else NA_real_,
                               end_pos    = if (.N) max(position) else NA_real_,
                               n_markers  = .N),
                             by = block_id]

  # Reorder allele columns to genomic position order
  ordered_markers <- marker_meta$marker
  ordered_markers <- ordered_markers[ordered_markers %in% names(alleles_wide)]
  data.table::setcolorder(alleles_wide,
                          c("Sample", "Haplotype_idx", ordered_markers))

  list(
    alleles     = alleles_wide,
    marker_meta = marker_meta,
    block_meta  = block_meta,
    samples     = sort(unique(sample_names)),
    encoding    = encoding
  )
}


#' Save a GenoMeta object to an RDS file
#'
#' @param geno_meta A valid \code{GenoMeta} list.
#' @param path      File path (should end in \code{.rds}).
#'
#' @return Invisible \code{NULL}.
#'
#' @seealso \code{\link{load_geno_meta}}
#' @export
save_geno_meta <- function(geno_meta, path) {
  if (!is_geno_meta(geno_meta))
    stop("save_geno_meta: 'geno_meta' is not a valid GenoMeta.")
  saveRDS(geno_meta, file = path)
  invisible(NULL)
}


#' Load a GenoMeta object from an RDS file
#'
#' @param path Path to an RDS file saved by \code{\link{save_geno_meta}}.
#'
#' @return A \code{GenoMeta} list.
#'
#' @seealso \code{\link{save_geno_meta}}
#' @export
load_geno_meta <- function(path) {
  if (!file.exists(path))
    stop("load_geno_meta: File not found: ", path)
  geno_meta <- readRDS(path)
  if (!is_geno_meta(geno_meta))
    stop("load_geno_meta: Not a valid GenoMeta.")
  geno_meta
}


#' Recode a GenoMeta object to a different allele encoding
#'
#' Converts the allele values in \code{geno_arrays} between
#' \code{"selectiontools"} (numeric codes) and \code{"nucleotide"} (characters)
#' using \code{\link{DEFAULT_ENCODING_DICT}} /
#' \code{\link{REVERSED_ENCODING_DICT}}, and updates \code{metadata$char_dict}
#' and \code{metadata$encoding}.
#'
#' @param geno_meta A valid \code{GenoMeta} list.
#' @param to_type   Target encoding: \code{"nucleotide"} or
#'   \code{"selectiontools"}.
#'
#' @return A new \code{GenoMeta} with recoded \code{geno_arrays}.
#'
#' @export
recode_geno_meta <- function(geno_meta, to_type = "nucleotide") {
  if (!is_geno_meta(geno_meta))
    stop("recode_geno_meta: 'geno_meta' is not a valid GenoMeta.")
  to_type <- match.arg(to_type, c("nucleotide", "selectiontools"))

  cur <- geno_meta$metadata$encoding
  if (is.null(cur))
    cur <- detect_genotype_encoding(.meta_to_allele_matrix(geno_meta))
  if (identical(cur, to_type)) return(geno_meta)

  arr  <- geno_meta$geno_arrays
  vals <- as.character(arr)

  if (to_type == "nucleotide") {
    recoded <- unname(REVERSED_ENCODING_DICT[vals])
  } else {
    fwd     <- setNames(as.character(unlist(DEFAULT_ENCODING_DICT)),
                        names(DEFAULT_ENCODING_DICT))
    recoded <- unname(fwd[vals])
  }

  arr2 <- array(recoded, dim = dim(arr), dimnames = dimnames(arr))

  # Rebuild char_dict for the new encoding from the alleles actually present
  present <- unique(stats::na.omit(as.character(arr2)))
  if (to_type == "nucleotide") {
    keep      <- intersect(names(DEFAULT_ENCODING_DICT), present)
    char_dict <- unlist(DEFAULT_ENCODING_DICT[keep])
  } else {
    char_dict <- setNames(as.integer(present), present)
  }

  new_meta               <- geno_meta$metadata
  new_meta$char_dict     <- char_dict
  new_meta$encoding      <- to_type

  new_geno_meta(
    geno_arrays           = arr2,
    map_info              = geno_meta$map_info,
    block_info            = geno_meta$block_info,
    block_marker_map_info = geno_meta$block_marker_map_info,
    metadata              = new_meta,
    sample_info           = geno_meta$sample_info,
    marker_effects        = geno_meta$marker_effects
  )
}


# ── Internal builders ─────────────────────────────────────────────────────────

.build_alleles_st <- function(geno_t, sample_names, marker_names) {
  n_samp <- nrow(geno_t); n_mark <- ncol(geno_t)
  hap1 <- matrix(NA_character_, n_samp, n_mark)
  hap2 <- matrix(NA_character_, n_samp, n_mark)
  for (j in seq_len(n_mark)) {
    parts <- strsplit(as.character(geno_t[, j]), "/", fixed = TRUE)
    hap1[, j] <- vapply(parts, function(p) if (length(p) >= 1L) p[1L] else NA_character_, character(1L))
    hap2[, j] <- vapply(parts, function(p) if (length(p) >= 2L) p[2L] else NA_character_, character(1L))
  }
  .bind_haplotypes(hap1, hap2, sample_names, marker_names)
}

.build_alleles_nuc <- function(geno_t, sample_names, marker_names) {
  n_samp <- nrow(geno_t); n_mark <- ncol(geno_t)
  sample_vals <- as.character(geno_t[1L, ])
  has_slash <- any(grepl("/", sample_vals, fixed = TRUE))
  hap1 <- matrix(NA_character_, n_samp, n_mark)
  hap2 <- matrix(NA_character_, n_samp, n_mark)
  for (j in seq_len(n_mark)) {
    vals <- as.character(geno_t[, j])
    if (has_slash) {
      parts <- strsplit(vals, "/", fixed = TRUE)
      hap1[, j] <- vapply(parts, function(p) if (length(p) >= 1L) p[1L] else NA_character_, character(1L))
      hap2[, j] <- vapply(parts, function(p) if (length(p) >= 2L) p[2L] else NA_character_, character(1L))
    } else {
      hap1[, j] <- substr(vals, 1L, 1L)
      hap2[, j] <- substr(vals, 2L, 2L)
      hap2[nchar(hap2[, j]) == 0L, j] <- NA_character_
    }
  }
  .bind_haplotypes(hap1, hap2, sample_names, marker_names)
}

.build_alleles_num <- function(geno_t, sample_names, marker_names) {
  dt <- as.data.table(as.data.frame(matrix(as.character(geno_t),
                                           nrow = nrow(geno_t),
                                           ncol = ncol(geno_t))))
  setnames(dt, marker_names)
  dt[, Sample := sample_names]
  dt[, Haplotype_idx := 1L]
  setcolorder(dt, c("Sample", "Haplotype_idx", marker_names))
  dt
}

.bind_haplotypes <- function(hap1, hap2, sample_names, marker_names) {
  colnames(hap1) <- marker_names; colnames(hap2) <- marker_names
  # check.names = FALSE: preserve non-syntactic marker IDs (e.g. "1H-123")
  # verbatim so they match geno_meta$map_info$Marker downstream.
  dt1 <- as.data.table(as.data.frame(hap1, check.names = FALSE)); dt1[, Sample := sample_names]; dt1[, Haplotype_idx := 1L]
  dt2 <- as.data.table(as.data.frame(hap2, check.names = FALSE)); dt2[, Sample := sample_names]; dt2[, Haplotype_idx := 2L]
  dt  <- rbindlist(list(dt1, dt2))
  setcolorder(dt, c("Sample", "Haplotype_idx", marker_names))
  dt
}

.normalise_map_dt <- function(map_data) {
  dt <- as.data.table(map_data)
  if (!all(c("Marker","Chrom","Pos") %in% names(dt))) {
    if (ncol(dt) >= 3L) setnames(dt, names(dt)[1:3], c("Marker","Chrom","Pos"))
    else stop("map_data needs at least 3 columns: Marker, Chrom, Pos.")
  }
  dt
}

.normalise_block_dt <- function(block_map) {
  dt <- as.data.table(block_map)
  if (!all(c("Chrom","Marker","Block") %in% names(dt))) {
    if (ncol(dt) >= 3L) setnames(dt, names(dt)[1:3], c("Chrom","Marker","Block"))
    else stop("block_map needs at least 3 columns: Chrom, Marker, Block.")
  }
  dt
}
