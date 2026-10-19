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

  # Resolve encoding: explicit arg > metadata$encoding > auto-detect. The
  # collapsed [n_samples x n_markers] genotype matrix is used for DETECTION
  # only: .meta_to_allele_matrix() writes NA for any cell with a missing copy,
  # so it must never be the source of the per-copy allele table (0.7.6,
  # discrepancy 29 — before that fix a half-called cell (A, NA) was re-split
  # into (NA, NA) here and lost its observed allele in every consumer).
  if (is.null(encoding)) encoding <- geno_meta$metadata$encoding
  if (is.null(encoding)) {
    encoding <- detect_genotype_encoding(.meta_to_allele_matrix(geno_meta))
  } else {
    encoding <- match.arg(encoding, c("selectiontools", "nucleotide", "numeric"))
  }

  # Build the wide alleles table. For the diploid (allele) encodings the 3-D
  # array IS the per-copy view: layer k of geno_arrays is haplotype k, one
  # allele or NA per cell, already normalised by new_geno_meta(). The dosage
  # encoding keeps its one-row-per-sample layout from the collapsed matrix.
  alleles_wide <- switch(
    encoding,
    selectiontools = ,
    nucleotide     = .build_alleles_layers(arr, sample_names, marker_names),
    numeric        = .build_alleles_num(.meta_to_allele_matrix(geno_meta),
                                        sample_names, marker_names),
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
#' @section Containers saved before 0.6.0:
#' Since 0.6.0 a missing allele inside a \code{GenoMeta} is \code{NA} and
#' nothing else. A file saved earlier can still hold the source token
#' (PLINK's \code{"0"}, \code{"N"}, \code{"."}, \code{"-"}) as an allele, which
#' shifts marker effects and frequency-centred local GEBVs. Loading one says so
#' and names \code{\link{repair_geno_meta}} — and \strong{does not rewrite the
#' file}. Repairing what you saved without being asked would be the same class
#' of silent decision the invariant exists to remove, so the repair is yours to
#' make and yours to persist:
#'
#' \preformatted{meta <- load_geno_meta("panel.rds")
#' meta <- repair_geno_meta(meta)
#' save_geno_meta(meta, "panel.rds")}
#'
#' @param path Path to an RDS file saved by \code{\link{save_geno_meta}}.
#'
#' @return A \code{GenoMeta} list, exactly as stored.
#'
#' @seealso \code{\link{save_geno_meta}}, \code{\link{repair_geno_meta}}
#' @export
load_geno_meta <- function(path) {
  if (!file.exists(path))
    stop("load_geno_meta: File not found: ", path)
  geno_meta <- readRDS(path)

  # Report the missing-token condition against the file it came from, so muffle
  # the generic is_geno_meta() warning and issue the path-aware one below. Under
  # options(haploStack.strict_missing = TRUE) is_geno_meta() errors instead, and
  # that error is left to propagate.
  ok <- withCallingHandlers(
    is_geno_meta(geno_meta),
    haploStack_missing_token = function(w) invokeRestart("muffleWarning"))
  if (!ok)
    stop("load_geno_meta: Not a valid GenoMeta.")

  rep <- .missing_token_report(geno_meta)
  if (!is.null(rep))
    message("load_geno_meta: '", path, "' holds ", rep$n_cells,
            " allele cell(s) across ", rep$n_markers, " marker(s) as a missing ",
            "token (", paste0("\"", rep$tokens, "\"", collapse = ", "),
            ") rather than NA — it predates the 0.6.0 invariant. The file has ",
            "NOT been modified. Repair the loaded object with ",
            "repair_geno_meta(), then save_geno_meta() it if you want the file ",
            "updated.")

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
#' @section Both ends are declared:
#' A transform needs to know what it is converting \emph{from} as well as
#' \emph{to}. Before 0.6.0 only \code{to_type} was an argument and the source
#' was taken from \code{metadata$encoding} as an early-return guard, so a
#' container recorded as \code{"numeric"} sailed past the guard and had its
#' dosages pushed through the nucleotide dictionary — \code{1 -> "A"},
#' \code{2 -> "C"}, \code{0 -> NA} — with no warning (INDEX discrepancy 16,
#' item 4).
#'
#' \code{from_type} now names the source explicitly. It defaults to the recorded
#' provenance, but the recorded value is \strong{validated against the data}
#' rather than trusted, and a dosage source with an allele target is refused
#' outright: recoding a vocabulary and reinterpreting a scale are different
#' operations, and the defect was one masquerading as the other.
#'
#' @section Rendering missing back out (\code{na_as}):
#' Inside a \code{GenoMeta} a missing allele is \code{NA}. Converting one back
#' to a character is possible, and irreversible in the direction that matters:
#' the character the \emph{source file} wrote there is gone the moment the
#' container was built, so \code{na_as} produces a new rendering, never a
#' recovery. Recoding a PED-derived container with \code{na_as = "-"} does not
#' give you back the \code{"0"} the PED had, and nothing in the package can.
#'
#' \code{na_as = NULL} (the default) keeps \code{NA} and preserves the
#' invariant. Supplying a string writes it in place of every missing allele —
#' useful when the target is a file format with its own absence code
#' (\code{"-"} for the nucleotide dictionary, \code{"-1"} for SelectionTools).
#' The result is then, deliberately, \strong{not} a normalised container:
#' \code{\link{is_geno_meta}} will say so, and \code{\link{repair_geno_meta}}
#' undoes the rendering. Treat it as an export step, not a storage format.
#'
#' @param geno_meta A valid \code{GenoMeta} list.
#' @param to_type   Target encoding: \code{"nucleotide"} or
#'   \code{"selectiontools"}.
#' @param from_type Source encoding. \code{NULL} (default) uses
#'   \code{metadata$encoding}, or detects it when that is absent; either way the
#'   value is checked against the alleles actually present.
#' @param na_as     Optional string to write in place of every missing allele.
#'   \code{NULL} (default) keeps \code{NA}. See \emph{Rendering missing back
#'   out}.
#'
#' @return A new \code{GenoMeta} with recoded \code{geno_arrays}.
#'
#' @seealso \code{\link{repair_geno_meta}}, \code{\link{use_missing_chars}}
#' @export
recode_geno_meta <- function(geno_meta, to_type = "nucleotide",
                             from_type = NULL, na_as = NULL) {
  if (!is_geno_meta(geno_meta))
    stop("recode_geno_meta: 'geno_meta' is not a valid GenoMeta.")
  to_type <- match.arg(to_type, c("nucleotide", "selectiontools"))

  if (!is.null(na_as)) {
    if (!is.character(na_as) || length(na_as) != 1L || is.na(na_as) ||
        !nzchar(na_as))
      stop("recode_geno_meta: 'na_as' must be a single non-empty string, ",
           "or NULL to keep NA.")
  }

  # ---- Resolve and validate the source encoding -----------------------------
  cur <- from_type %||% geno_meta$metadata$encoding
  if (is.null(cur))
    cur <- detect_genotype_encoding(.meta_to_allele_matrix(geno_meta))
  cur <- .encode_arg(cur, "recode_geno_meta",
                     allowed = c("selectiontools", "nucleotide", "numeric",
                                 "allele"),
                     strict  = FALSE)

  # A dosage scale is not an allele vocabulary. Refuse rather than translate.
  if (cur %in% c("numeric", "allele"))
    stop("recode_geno_meta: source encoding '", cur, "' is a dosage scale, ",
         "not an allele vocabulary — its values count copies of an allele, ",
         "they do not name one. Recoding it to '", to_type, "' would read ",
         "those counts as dictionary keys (1 -> \"A\", 2 -> \"C\"), which is a ",
         "reinterpretation, not a recode. Convert the dosage to alleles at ",
         "read time instead.")

  .validate_recode_source(geno_meta$geno_arrays, cur)

  if (identical(cur, to_type) && is.null(na_as)) return(geno_meta)

  arr  <- geno_meta$geno_arrays
  vals <- as.character(arr)

  if (identical(cur, to_type)) {
    recoded <- vals                       # na_as-only pass, no vocabulary change
  } else if (to_type == "nucleotide") {
    recoded <- unname(REVERSED_ENCODING_DICT[vals])
  } else {
    fwd     <- setNames(as.character(unlist(DEFAULT_ENCODING_DICT)),
                        names(DEFAULT_ENCODING_DICT))
    recoded <- unname(fwd[vals])
  }

  # ---- Render missing, if asked ---------------------------------------------
  n_rendered <- 0L
  if (!is.null(na_as)) {
    n_rendered <- sum(is.na(recoded))
    recoded[is.na(recoded)] <- na_as
  }

  arr2 <- array(recoded, dim = dim(arr), dimnames = dimnames(arr))

  # Rebuild char_dict for the new encoding from the alleles actually present
  present <- unique(stats::na.omit(as.character(arr2)))
  if (!is.null(na_as)) present <- setdiff(present, na_as)
  if (to_type == "nucleotide") {
    keep      <- intersect(names(DEFAULT_ENCODING_DICT), present)
    char_dict <- unlist(DEFAULT_ENCODING_DICT[keep])
  } else {
    char_dict <- setNames(as.integer(present), present)
  }

  new_meta               <- geno_meta$metadata
  new_meta$char_dict     <- char_dict
  new_meta$encoding      <- to_type
  # The old stamp described the old vocabulary; drop it so the portal restamps
  # against the target one.
  new_meta$missing_tokens <- NULL

  if (n_rendered > 0L)
    message("recode_geno_meta: rendered ", n_rendered,
            " missing allele(s) as \"", na_as, "\". The result is no longer a ",
            "normalised container — repair_geno_meta() undoes the rendering, ",
            "but the source file's own missing code is not recoverable.")

  out <- new_geno_meta(
    geno_arrays           = arr2,
    map_info              = geno_meta$map_info,
    block_info            = geno_meta$block_info,
    block_marker_map_info = geno_meta$block_marker_map_info,
    metadata              = new_meta,
    sample_info           = geno_meta$sample_info,
    marker_effects        = geno_meta$marker_effects,
    source_format         = to_type
  )

  # A deliberate rendering must survive the portal it just passed through.
  if (!is.null(na_as)) {
    out$geno_arrays <- arr2
    out$metadata$missing_tokens <- NULL
  }
  # Recoding changes the allele vocabulary, not the blocks — a stored block
  # annotation stays valid (#033). (The phenotype, by contrast, is a known
  # casualty of this rebuild; see geno_meta_structure.md.)
  out$annotation <- geno_meta$annotation
  out
}


# Internal: check a declared source encoding against the alleles actually
# present, so a stale or wrong metadata$encoding fails loudly instead of being
# translated through the wrong dictionary. Only impossible combinations are
# refused — "1" and "2" are legal under more than one vocabulary, so ambiguity
# is left alone rather than guessed at.
.validate_recode_source <- function(arr, from_type) {
  present <- unique(stats::na.omit(as.character(arr)))
  present <- setdiff(present, .missing_tokens(from_type))
  if (length(present) == 0L) return(invisible(TRUE))

  alpha <- grepl("[A-Za-z]", present)
  if (identical(from_type, "selectiontools") && any(alpha))
    stop("recode_geno_meta: source declared 'selectiontools' but the container ",
         "holds alphabetic allele(s): ",
         paste0("\"", utils::head(present[alpha], 5L), "\"", collapse = ", "),
         ". Pass from_type = \"nucleotide\", or fix metadata$encoding.")
  if (identical(from_type, "nucleotide") && !any(alpha))
    stop("recode_geno_meta: source declared 'nucleotide' but no allele in the ",
         "container is alphabetic (present: ",
         paste0("\"", utils::head(present, 5L), "\"", collapse = ", "),
         "). Pass from_type = \"selectiontools\", or fix metadata$encoding.")
  invisible(TRUE)
}


# ── Internal builders ─────────────────────────────────────────────────────────

# .build_alleles_layers(arr, sample_names, marker_names)
#
# Per-copy allele table straight from the ploidy layers of geno_arrays: one
# block of rows per layer k (Haplotype_idx = k), one column per marker, cells
# exactly as stored (an allele or NA). This is the only place the views may
# read alleles from, because it is the only route that keeps a half-called
# cell's observed copy: (A, NA) stays A on copy 1 and NA on copy 2. The old
# builders re-split the collapsed genotype string, which is NA whenever ANY
# copy is missing, and so turned (A, NA) into (NA, NA) (found 2026-09-18 in
# the barley bench; 0.7.6). Layer order is preserved; nothing is sorted.
.build_alleles_layers <- function(arr, sample_names, marker_names) {
  n_samp   <- length(sample_names)
  n_mark   <- length(marker_names)
  n_ploidy <- dim(arr)[3L]
  parts <- vector("list", n_ploidy)
  for (k in seq_len(n_ploidy)) {
    # matrix(): a 1-sample or 1-marker slice drops to a vector; refill it in
    # the array's column-major order so the shape is always [n_samp x n_mark].
    hap <- matrix(as.character(arr[, , k]), nrow = n_samp, ncol = n_mark)
    colnames(hap) <- marker_names
    # check.names = FALSE: preserve non-syntactic marker IDs (e.g. "1H-123")
    # verbatim so they match geno_meta$map_info$Marker downstream.
    dt <- as.data.table(as.data.frame(hap, check.names = FALSE,
                                      stringsAsFactors = FALSE))
    dt[, Sample := sample_names]
    dt[, Haplotype_idx := as.integer(k)]
    parts[[k]] <- dt
  }
  dt <- rbindlist(parts)
  setcolorder(dt, c("Sample", "Haplotype_idx", marker_names))
  dt
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
