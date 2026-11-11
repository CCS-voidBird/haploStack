# =============================================================================
# geno_meta_io.R — the container on disk and as matrices
#
# save_geno_meta() / load_geno_meta() (trusted input only; schema upgrade in
# memory), geno_meta_views() (the wide allele table and the marker / block
# tables the catalogue reads), .meta_to_allele_matrix() /
# .meta_to_dosage_matrix() (the two read seams), geno_meta_to_matrix() and
# write_geno_meta() (the writers, by preset). Split out of format_helper.R and
# catalogue_profile.R in 0.9.0 (D2-9).
# =============================================================================

# .meta_to_allele_matrix(geno_meta, sep = "/")
#
# Reconstruct a 2-D character genotype matrix from the 3-D split array.
# Each cell is the allele calls from all ploidy layers joined by `sep`
# (e.g. "A/T" for diploid). NA in any ploidy layer → NA in the output cell.
# No subsetting is performed; call filter_geno_meta() first if needed.
#
# Returns: character matrix [n_samples × n_markers].
.meta_to_allele_matrix <- function(geno_meta, sep = "/") {
  arr          <- geno_meta$geno_arrays
  n_samples    <- dim(arr)[1L]
  n_markers    <- dim(arr)[2L]
  n_ploidy     <- dim(arr)[3L]
  sample_names <- dimnames(arr)[[1L]]
  marker_names <- dimnames(arr)[[2L]]

  result  <- arr[, , 1L, drop = TRUE]          # [n_samples × n_markers] character
  na_mask <- is.na(result)

  for (k in seq(2L, n_ploidy)) {
    layer_k  <- arr[, , k, drop = TRUE]
    na_mask  <- na_mask | is.na(layer_k)
    result   <- matrix(paste(result, layer_k, sep = sep),
                       nrow = n_samples, ncol = n_markers)
  }

  result[na_mask] <- NA_character_
  dimnames(result) <- list(sample_names, marker_names)
  result   # [n_samples × n_markers]
}


# .meta_to_dosage_matrix(geno_meta, target_alleles = NULL)
#
# Count how many ploidy layers carry a specified allele per marker × sample.
#
# target_alleles: character vector, length == n_markers (positionally matched
#   to dim 1 of geno_arrays). No names required.
#   NULL → auto minor allele per marker (least frequent across all layers).
#
# Strict: if length(target_alleles) != n_markers an error is thrown.
# Subsetting must be done by the caller via filter_geno_meta() first.
# NA in any ploidy layer → NA for that marker × sample cell.
#
# Returns: integer matrix [n_samples × n_markers].
.meta_to_dosage_matrix <- function(geno_meta, target_alleles = NULL) {
  arr          <- geno_meta$geno_arrays
  n_samples    <- dim(arr)[1L]
  n_markers    <- dim(arr)[2L]
  n_ploidy     <- dim(arr)[3L]
  sample_names <- dimnames(arr)[[1L]]
  marker_names <- dimnames(arr)[[2L]]

  if (!is.null(target_alleles)) {
    if (length(target_alleles) != n_markers)
      .hs_abort(".meta_to_dosage_matrix: length(target_alleles) (", length(target_alleles),
           ") must equal n_markers (", n_markers, "). ",
           "Call filter_geno_meta() to subset the GenoMeta to the desired markers first.")
    target_alleles <- as.character(target_alleles)
  } else {
    # Auto: minor allele = least frequent OBSERVED allele per marker (ties broken
    # alphabetically, matching names(sort(table(obs)))[1L]). Detection is in R;
    # counting is the marker_allele_counts kernel.
    #
    # `na_vals` is defence-in-depth (FR7): a no-op on a container that came
    # through the portal, and on a legacy one it stops a missing code being
    # picked as the minor allele — which would make the dosage a count of
    # missingness.
    counts <- .marker_allele_counts(
      arr, na_vals = .missing_tokens(geno_meta$metadata$encoding,
                                     .meta_extra_missing(geno_meta)))
    levs   <- colnames(counts)
    target_alleles <- vapply(seq_len(n_markers), function(m) {
      cj <- counts[m, ]
      if (all(cj == 0L)) return(NA_character_)
      cj[cj == 0L] <- .Machine$integer.max         # ignore alleles unobserved here
      levs[which.min(cj)]                          # least frequent; ties -> alpha
    }, character(1L))
  }

  # Dosage = count of target allele per cell; NA if ANY ploidy layer is missing
  # (na_if_any = TRUE preserves the original has_na semantics).
  .dosage_matrix(arr, target_alleles, na_if_any = TRUE)
}


# ---------------------------------------------------------------------------
# Block attachment
# ---------------------------------------------------------------------------


#' Convert a GenoMeta 3-D array to a 2-D genotype matrix
#'
#' Collapses the haplotype dimension of \code{geno_arrays} into genotype calls
#' laid out per a format \code{preset} — the \emph{same} preset/parameter-set
#' mechanism used by \code{\link{read_genotype}}, so reading and writing share
#' one description of each format. The preset's \code{geno_orient} sets the
#' output orientation and \code{split_char} the in-cell allele separator; the
#' map is returned with its columns reordered to the preset's \code{map_order}.
#'
#' @param meta          A valid \code{GenoMeta} list.
#' @param preset        Format preset name, a custom parameter-set list, or
#'   \code{NULL}. See \code{\link{GENO_FORMAT_PRESETS}}.
#' @param output_format Deprecated since 0.9.0 (warns on every call, see
#'   \link{haploStack-deprecated}): the old spelling of \code{preset}, with the
#'   legacy values \code{"plink_ped"}, \code{"st"}, \code{"gs"}. Used only when
#'   \code{preset} is \code{NULL}.
#' @param presets       Preset registry. Default \code{\link{GENO_FORMAT_PRESETS}}.
#' @param split_char    Optional override for the in-cell allele separator.
#' @param st_formats    For the \code{"ST"} preset (SelectionTools' layout):
#'   \code{TRUE} enables it for this call; \code{NULL} (default) takes the
#'   session switch (\code{\link{use_st_formats}}, off by default). The cells
#'   written are the container's nucleotides; \code{\link{st_convert}} writes
#'   SelectionTools' own codes.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{matrix}}{2-D character matrix in the preset's orientation.}
#'   \item{\code{map}}{data.frame with columns ordered per the preset's
#'     \code{map_order}.}
#'   \item{\code{fam}}{The \code{sample_info} FAM table.}
#'   \item{\code{output_format}}{The resolved preset name/spec source.}
#'   \item{\code{spec}}{The resolved format parameter set.}
#' }
#'
#' @seealso \code{\link{resolve_geno_format}}, \code{\link{write_geno_meta}}
#'
#' @export
geno_meta_to_matrix <- function(meta,
                                 preset        = NULL,
                                 output_format = NULL,
                                 presets       = GENO_FORMAT_PRESETS,
                                 split_char    = NULL,
                                 st_formats    = NULL) {
  .hs_check_meta(meta, "geno_meta_to_matrix", "meta")

  # output_format is on the deprecation ladder since 0.9.0 (Q-c a; it never
  # warned before): mapped as before, defunct at 0.9.x.
  if (!is.null(output_format)) {
    .hs_deprecated("output_format", "preset", "0.9.0", what = "argument",
                   fn = "geno_meta_to_matrix")
    if (is.null(preset)) preset <- output_format
  }
  if (is.null(preset)) preset <- "PL"
  spec <- resolve_geno_format(preset, presets = presets,
                              split_char = split_char)
  # The SelectionTools layout is behind use_st_formats() like every other
  # SelectionTools format (C-6). The cells written are the container's
  # nucleotides; st_convert(to = "selectiontools") writes the codes.
  if (identical(spec$encoding, "selectiontools"))
    .st_gate("geno_meta_to_matrix", st_formats)

  # Join ploidy layers with the preset separator → [n_samples × n_markers].
  am <- .meta_to_allele_matrix(meta, sep = spec$split_char)

  # Orient: "sample" keeps samples-as-rows (the native layout); "marker"
  # transposes to markers-as-rows.
  mat <- if (spec$geno_orient == "marker") t(am) else am

  # Map columns reordered to the preset's file order.
  base_map <- meta$map_info[, c("Chrom", "Marker", "Pos"), drop = FALSE]
  map_out  <- base_map[, spec$map_order, drop = FALSE]

  list(matrix        = mat,
       map           = map_out,
       fam           = meta$sample_info,
       output_format = if (is.character(preset)) preset else "custom",
       spec          = spec)
}


#' Write a GenoMeta object to disk in a specified format
#'
#' Resolves a format \code{preset} (the same parameter-set mechanism as
#' \code{\link{read_genotype}}) and writes the genotype to \code{path}. When the
#' preset embeds FAM columns (\code{fam_cols > 0}, e.g. \code{"PL"}), a PLINK
#' \code{.ped}/\code{.map} pair is written: the FAM block is taken from
#' \code{meta$sample_info} and alleles are space-separated. Otherwise a single
#' delimited file is written in the preset's orientation, with the row
#' identifier (marker or sample IDs) in the first column.
#'
#' That single file is a **CSV**. The \code{"ST"} preset gives SelectionTools'
#' \emph{layout} only (markers in rows, slash-separated cells): the alleles are
#' written as the container holds them — nucleotides since 0.9.0 — and
#' SelectionTools' own reader (\code{st.read.marker.data}) splits fields on a
#' tab, a blank or \code{;}, so it would read this file as zero individuals.
#' SelectionTools' own formats, with its \code{1/2/3/4} codes, are written by
#' \code{\link{st_convert}(to = "selectiontools")}; \code{\link{st_ld_blocks}}
#' writes them for its own run (and consider first whether you need LD
#' blocking: a block map you already have goes to \code{\link{attach_block_info}}).
#'
#' @param meta          A valid \code{GenoMeta} list.
#' @param path          Output file path (without extension for PLINK output;
#'   full path including extension otherwise).
#' @param preset        Format preset name, a custom parameter-set list, or
#'   \code{NULL}. See \code{\link{GENO_FORMAT_PRESETS}}.
#' @param output_format Deprecated since 0.9.0 (warns on every call, see
#'   \link{haploStack-deprecated}): the old spelling of \code{preset}, with the
#'   legacy values \code{"plink_ped"}, \code{"st"}, \code{"gs"}. Used only when
#'   \code{preset} is \code{NULL}.
#' @param presets       Preset registry. Default \code{\link{GENO_FORMAT_PRESETS}}.
#' @param split_char    Optional override for the in-cell allele separator.
#' @param st_formats    For the \code{"ST"} preset: \code{TRUE} enables
#'   SelectionTools' layout for this call; \code{NULL} (default) takes the
#'   session switch (\code{\link{use_st_formats}}). See
#'   \code{\link{geno_meta_to_matrix}}.
#' @param pos_unit      Unit for written map positions: \code{"bp"} (default) or
#'   \code{"mbp"}. Positions are stored internally in bp and converted on write.
#'
#' @return Invisibly, the result of \code{\link{geno_meta_to_matrix}}.
#'
#' @seealso \code{\link{resolve_geno_format}}, \code{\link{geno_meta_to_matrix}}
#'
#' @export
write_geno_meta <- function(meta, path,
                            preset        = NULL,
                            output_format = NULL,
                            presets       = GENO_FORMAT_PRESETS,
                            split_char    = NULL,
                            pos_unit      = "bp",
                            st_formats    = NULL) {
  if (!is.null(output_format)) {
    .hs_deprecated("output_format", "preset", "0.9.0", what = "argument",
                   fn = "write_geno_meta")
    if (is.null(preset)) preset <- output_format
  }
  if (is.null(preset)) preset <- "PL"
  result <- geno_meta_to_matrix(meta, preset = preset, presets = presets,
                                split_char = split_char, st_formats = st_formats)
  spec   <- result$spec

  if (spec$fam_cols > 0L) {
    # ---- PLINK-style output: .ped (FAM block + space-separated alleles) ------
    arr     <- meta$geno_arrays
    sn      <- dimnames(arr)[[1L]]
    mn      <- dimnames(arr)[[2L]]
    ploidy  <- dim(arr)[3L]

    allele_block <- matrix(NA_character_, nrow = length(sn), ncol = length(mn),
                           dimnames = list(sn, mn))
    for (j in seq_along(mn)) {
      layers <- vapply(seq_len(ploidy),
                       function(h) as.character(arr[, j, h]),
                       character(length(sn)))
      if (is.null(dim(layers))) layers <- matrix(layers, nrow = length(sn))
      # A missing allele is PLINK's "0". Writing R's NA as the text "NA" made
      # read_plink_ped() read a fully missing call back as the half-call
      # (NA, "A") (0.8.8 and earlier).
      layers[is.na(layers)] <- "0"
      allele_block[, j] <- apply(layers, 1L, function(r) paste(r, collapse = " "))
    }

    si <- meta$sample_info
    ped_df <- cbind(
      si$Family, si$Sample, si$Parent1, si$Parent2, si$Sex, -9L,
      allele_block
    )
    ped_path <- paste0(path, ".ped")
    utils::write.table(ped_df, ped_path, sep = " ",
                       row.names = FALSE, col.names = FALSE, quote = FALSE)

    # Companion .map: chr, marker, cM (0), position (in requested unit).
    map_out <- data.frame(
      Chrom  = meta$map_info$Chrom,
      Marker = meta$map_info$Marker,
      cM     = 0,
      Pos    = .convert_pos_unit(suppressWarnings(as.numeric(meta$map_info$Pos)),
                                 from_unit = "bp", to_unit = pos_unit),
      stringsAsFactors = FALSE
    )
    map_path <- paste0(path, ".map")
    utils::write.table(map_out, map_path, sep = "\t",
                       row.names = FALSE, col.names = FALSE, quote = FALSE)
    .hs_inform("Written: ", ped_path, " and ", map_path)

  } else {
    # ---- Single delimited file in the preset's orientation -------------------
    mat      <- result$matrix
    id_label <- if (spec$geno_orient == "marker") "Marker" else "Sample"
    out_df   <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
    out_df   <- cbind(setNames(
                        data.frame(rownames(mat), stringsAsFactors = FALSE),
                        id_label),
                      out_df)
    utils::write.csv(out_df, path, row.names = FALSE, quote = FALSE)
    .hs_inform("Written: ", path)
  }

  invisible(result)
}


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

  .hs_check_meta(geno_meta, "geno_meta_views", "geno_meta",
                 hint = paste0("Use new_geno_meta() or a read_*() function to create ",
                               "one, then attach_block_info() to add block assignments."))

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
    encoding <- match.arg(encoding, c("nucleotide", "numeric"))
  }

  # Build the wide alleles table. For the diploid (allele) encodings the 3-D
  # array IS the per-copy view: layer k of geno_arrays is haplotype k, one
  # allele or NA per cell, already normalised by new_geno_meta(). The dosage
  # encoding keeps its one-row-per-sample layout from the collapsed matrix.
  alleles_wide <- switch(
    encoding,
    nucleotide     = .build_alleles_layers(arr, sample_names, marker_names),
    numeric        = .build_alleles_num(.meta_to_allele_matrix(geno_meta),
                                        sample_names, marker_names),
    .hs_abort("geno_meta_views: unsupported encoding '", encoding, "'.")
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

  # 0.9.0 (row 29 a): the one natural chromosome order, then position
  marker_meta <- marker_meta[order(.hs_chrom_rank(chromosome), position)]

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
    samples     = .hs_sort(unique(sample_names)),
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
  .hs_check_meta(geno_meta, "save_geno_meta", "geno_meta")
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
#' @section Trusted input only:
#' \code{readRDS()} unserialises arbitrary R objects, which is not a safe
#' boundary for a file you did not write yourself (R core's advice, and a
#' decompression-bomb vector). Before reading, the file is checked to be
#' under \code{max_bytes} and to start like an R serialisation (a gzip,
#' bzip2 or xz header, or the uncompressed \code{X}/\code{A}/\code{B} header);
#' anything else is refused with an \code{hs_error_format} error. That is a
#' plausibility check, not a security boundary: \strong{never call this on a
#' user upload}. A service accepts genotypes through the text readers
#' (\code{\link{read_genotype}}, \code{\link{read_plink_ped}}, ...) and keeps
#' containers in its own store.
#'
#' @section Containers saved by an earlier haploStack:
#' A container carries \code{metadata$schema_version} since 0.9.0. A file
#' written by 0.6.0 to 0.8.8 (schema 1) is upgraded in memory with
#' \code{\link{upgrade_geno_meta}} — metadata only, one message — and the file
#' is left unchanged; \code{\link{save_geno_meta}} persists the upgraded
#' object if you want it. A file from a newer haploStack than this one is
#' refused (\code{hs_error_schema}).
#'
#' @param path Path to an RDS file saved by \code{\link{save_geno_meta}}.
#' @param max_bytes Numeric. Refuse a file larger than this many bytes
#'   (\code{Inf}, the default, accepts any size).
#' @param strict_missing Logical or \code{NULL}. Whether a container that still
#'   carries a missing-allele token is an error (\code{TRUE}) or reported with
#'   a message (\code{FALSE}); \code{NULL} (default) takes
#'   \code{options(haploStack.strict_missing)}.
#'
#' @return A \code{GenoMeta} list at the current schema (the stored object,
#'   upgraded in memory when it was saved by an earlier haploStack).
#'
#' @seealso \code{\link{save_geno_meta}}, \code{\link{upgrade_geno_meta}},
#'   \code{\link{repair_geno_meta}}
#' @export
load_geno_meta <- function(path, max_bytes = Inf, strict_missing = NULL) {
  if (!file.exists(path))
    .hs_abort("load_geno_meta: File not found: ", path)
  .hs_check_rds_file(path, max_bytes = max_bytes, caller = "load_geno_meta")
  geno_meta <- readRDS(path)

  if (!is.list(geno_meta) || !is.list(geno_meta$metadata))
    .hs_abort(class = c("hs_error_schema", "hs_error_input"),
              "load_geno_meta: '", path, "' is not a GenoMeta container (no ",
              "metadata list).")

  # An older schema is upgraded in memory (metadata only); the file stays as
  # it is. The token notice is handled below, against the path.
  muffle <- function(w) invokeRestart("muffleWarning")
  v <- .hs_schema_of(geno_meta)
  if (v > .HS_SCHEMA_VERSION)
    .hs_abort(class = c("hs_error_schema", "hs_error_input"),
              "load_geno_meta: '", path, "' carries container schema ", v,
              ", newer than this haploStack writes (", .HS_SCHEMA_VERSION,
              "); it was saved by haploStack ",
              geno_meta$metadata$built_with %||% "(unknown)",
              ". Upgrade the package to read it.")
  if (v < .HS_SCHEMA_VERSION) {
    geno_meta <- withCallingHandlers(
      upgrade_geno_meta(geno_meta, verbose = FALSE, strict_missing = FALSE),
      haploStack_missing_token = muffle)
    .hs_inform("load_geno_meta: '", path, "' was saved at container schema ",
               v, "; upgraded to schema ", .HS_SCHEMA_VERSION, " in memory ",
               "(metadata only). The file is unchanged; save_geno_meta() ",
               "persists the upgrade.")
  } else if (identical(geno_meta$metadata$encoding, "selectiontools")) {
    # A schema-2 container saved by a 0.9.0 development build before C-5 still
    # holds SelectionTools codes: the recode runs whatever the schema says (the
    # one migration that needs no use_st_formats(); C-6, found by test S12).
    geno_meta <- withCallingHandlers(
      upgrade_geno_meta(geno_meta, verbose = FALSE, strict_missing = FALSE),
      haploStack_missing_token = muffle)
    .hs_inform("load_geno_meta: '", path, "' stores SelectionTools allele codes ",
               "(1/2/3/4); recoded to nucleotides in memory. The file is ",
               "unchanged; save_geno_meta() persists the recode.")
  }

  # Report the missing-token condition against the file it came from, so muffle
  # the generic is_geno_meta() warning and issue the path-aware one below. Under
  # strict_missing = TRUE (argument or option) is_geno_meta() errors instead,
  # and that error is left to propagate.
  ok <- withCallingHandlers(
    is_geno_meta(geno_meta, strict = TRUE, strict_missing = strict_missing),
    haploStack_missing_token = muffle)
  if (!ok)
    .hs_abort("load_geno_meta: Not a valid GenoMeta.")

  rep <- .missing_token_report(geno_meta)
  if (!is.null(rep))
    .hs_inform("load_geno_meta: '", path, "' holds ", rep$n_cells,
            " allele cell(s) across ", rep$n_markers, " marker(s) as a missing ",
            "token (", paste0("\"", rep$tokens, "\"", collapse = ", "),
            ") rather than NA \u2014 it predates the 0.6.0 invariant. The file has ",
            "NOT been modified. Repair the loaded object with ",
            "repair_geno_meta(), then save_geno_meta() it if you want the file ",
            "updated.")

  geno_meta
}


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
