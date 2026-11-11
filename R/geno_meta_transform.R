# =============================================================================
# geno_meta_transform.R — transforms that return a new container
#
# filter_geno_meta() (markers / samples / chromosomes, one rebuild path) and
# recode_geno_meta() (the explicit nucleotide recode of a container that still
# carries SelectionTools codes; to_type = "selectiontools" is defunct). Split
# out of format_helper.R and catalogue_profile.R in 0.9.0 (D2-9).
# =============================================================================

#' Filter a GenoMeta object by samples, markers, blocks, or chromosomes
#'
#' Intersects \code{geno_arrays} dim 1 with \code{map_info$Marker} and
#' optionally restricts to user-supplied samples, markers, haplotype blocks, or
#' chromosome regions. Multiple marker-level criteria (\code{keep_markers},
#' \code{blocks}, \code{chromosomes}) are combined with \strong{AND}
#' (intersection) logic. \code{block_marker_map_info} is pruned to the retained
#' markers and a block whose every marker left loses its \code{block_info} row
#' (the row is kept while at least one marker remains). \code{sample_info},
#' \code{marker_effects}, \code{phenotype} and \code{annotation} are subset to
#' the survivors, \code{map_info} follows the array's marker order, and the
#' metadata counters follow the array (one rebuild path, 0.9.0).
#'
#' @param meta         A \code{GenoMeta} list (from \code{\link{new_geno_meta}}
#'   or any reader function).
#' @param keep_markers Optional character vector of marker names to retain.
#'   Further intersected with the existing marker set.
#' @param keep_samples Optional character vector of sample names to retain.
#' @param blocks       Optional character vector of block IDs to retain
#'   (resolved to markers via \code{block_marker_map_info}).
#' @param chromosomes  Optional list of chromosome region specs, each a list
#'   with a required \code{chr} field and optional \code{start}/\code{end}
#'   fields, e.g. \code{list(list(chr = "1", start = 1e6, end = 5e6))}.
#'
#' @return A new \code{GenoMeta} list (input is not modified).
#'
#' @export
filter_geno_meta <- function(meta, keep_markers = NULL, keep_samples = NULL,
                             blocks = NULL, chromosomes = NULL) {
  .hs_check_meta(meta, "filter_geno_meta", "meta")

  current_samples <- dimnames(meta$geno_arrays)[[1L]]
  current_markers <- dimnames(meta$geno_arrays)[[2L]]

  # Intersect genotype dim2 with map_info
  valid_markers <- intersect(current_markers, meta$map_info$Marker)

  # --- Marker-level criteria (AND logic across keep_markers/blocks/chromosomes)
  marker_sets <- list()

  if (!is.null(keep_markers)) {
    not_found <- setdiff(keep_markers, valid_markers)
    if (length(not_found) > 0L)
      .hs_warn("filter_geno_meta: ", length(not_found),
              " requested marker(s) not present: ",
              paste(head(not_found, 5L), collapse = ", "),
              if (length(not_found) > 5L) " ..." else "")
    marker_sets$markers <- intersect(valid_markers, keep_markers)
  }

  if (!is.null(blocks)) {
    bmm <- meta$block_marker_map_info
    avail_blocks <- if (!is.null(bmm)) unique(bmm$Block) else character(0L)
    miss <- setdiff(blocks, avail_blocks)
    if (length(miss) > 0L)
      .hs_warn("filter_geno_meta: ", length(miss),
              " requested block(s) not present: ",
              paste(head(miss, 5L), collapse = ", "),
              if (length(miss) > 5L) " ..." else "")
    blk_markers <- if (!is.null(bmm)) bmm$Marker[bmm$Block %in% blocks]
                   else character(0L)
    marker_sets$blocks <- intersect(valid_markers, blk_markers)
  }

  if (!is.null(chromosomes)) {
    map <- meta$map_info
    chr_markers <- character(0L)
    for (region in chromosomes) {
      if (is.null(region$chr)) {
        .hs_warn("filter_geno_meta: region missing 'chr' field, skipping.")
        next
      }
      chr_val <- as.character(region$chr)
      hit     <- map$Chrom == chr_val
      if (!is.null(region$start)) hit <- hit & as.numeric(map$Pos) >= region$start
      if (!is.null(region$end))   hit <- hit & as.numeric(map$Pos) <= region$end
      chr_markers <- union(chr_markers, map$Marker[hit])
    }
    marker_sets$chromosomes <- intersect(valid_markers, chr_markers)
  }

  if (length(marker_sets) > 0L)
    valid_markers <- Reduce(intersect, marker_sets)

  valid_samples <- current_samples
  if (!is.null(keep_samples)) {
    not_found_s <- setdiff(keep_samples, current_samples)
    if (length(not_found_s) > 0L)
      .hs_warn("filter_geno_meta: ", length(not_found_s),
              " requested sample(s) not present: ",
              paste(head(not_found_s, 5L), collapse = ", "),
              if (length(not_found_s) > 5L) " ..." else "")
    valid_samples <- intersect(current_samples, keep_samples)
  }

  # Subset array
  new_arr <- meta$geno_arrays[valid_samples, valid_markers, , drop = FALSE]

  # One rebuild path (0.9.0): the map is re-ordered to the array, the block
  # tables, FAM, effects, phenotype and annotation are pruned to the survivors,
  # a stored signal-marker selection loses the markers that left (with a
  # message), and the metadata counters follow the array. A block whose every
  # marker left loses its block_info row (form row 6).
  .rebuild_geno_meta(meta, geno_arrays = new_arr, .caller = "filter_geno_meta")
}


# ---------------------------------------------------------------------------
# Internal: 3-D array → 2-D conversions
# ---------------------------------------------------------------------------


#' Recode a GenoMeta object to a different allele encoding
#'
#' Converts the allele values in \code{geno_arrays} from SelectionTools'
#' numeric codes (\code{"selectiontools"}) to nucleotide characters
#' (\code{"nucleotide"}) using \code{\link{REVERSED_ENCODING_DICT}}, and
#' updates \code{metadata$char_dict} and \code{metadata$encoding}. Since 0.9.0
#' this is the explicit route for a container that still carries the old
#' vocabulary (\code{\link{upgrade_geno_meta}} and \code{\link{load_geno_meta}}
#' do it on their own); \code{to_type = "selectiontools"} is \strong{defunct},
#' because a \code{GenoMeta} stores nucleotide alleles only — write
#' SelectionTools' formats with \code{\link{st_convert}}.
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
#' @param to_type   Target encoding: \code{"nucleotide"} (the only one a
#'   container stores). \code{"selectiontools"} is defunct since 0.9.0 (see
#'   \link{haploStack-defunct}).
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
  .hs_check_meta(geno_meta, "recode_geno_meta", "geno_meta")
  to_type <- match.arg(to_type, c("nucleotide", "selectiontools"))
  # Defunct without a deprecated rung (0.9.0, C-5): the result could not be a
  # container, since a GenoMeta stores nucleotide alleles only.
  if (to_type == "selectiontools")
    .hs_abort(class = c("hs_error_defunct", "defunctError", "hs_error_input"),
              data  = list(old = "recode_geno_meta(to_type = \"selectiontools\")",
                           new = "st_convert(to = \"selectiontools\")",
                           since = "0.9.0"),
              "recode_geno_meta(to_type = \"selectiontools\") is defunct since ",
              "haploStack 0.9.0: a GenoMeta stores nucleotide alleles only, so ",
              "there is no container to return. Write SelectionTools' formats ",
              "with st_convert(to = \"selectiontools\") after use_st_formats(TRUE).",
              "\nSee help(\"haploStack-defunct\").", call. = FALSE)

  if (!is.null(na_as)) {
    if (!is.character(na_as) || length(na_as) != 1L || is.na(na_as) ||
        !nzchar(na_as))
      .hs_abort("recode_geno_meta: 'na_as' must be a single non-empty string, ",
           "or NULL to keep NA.")
  }

  # ---- Resolve and validate the source encoding -----------------------------
  cur <- from_type %||% geno_meta$metadata$encoding
  if (is.null(cur))     # the internal detector: an old container may hold codes
    cur <- .detect_encoding(.meta_to_allele_matrix(geno_meta))
  cur <- .encode_arg(cur, "recode_geno_meta",
                     allowed = c("selectiontools", "nucleotide", "numeric",
                                 "allele"),
                     strict  = FALSE)

  # A dosage scale is not an allele vocabulary. Refuse rather than translate.
  if (cur %in% c("numeric", "allele"))
    .hs_abort("recode_geno_meta: source encoding '", cur, "' is a dosage scale, ",
         "not an allele vocabulary \u2014 its values count copies of an allele, ",
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
  } else {
    # the adapter's recode (C-6): 1/2/3/4 -> A/C/G/T, -1 -> "-", anything
    # else verbatim (a change of vocabulary, not a repair — as upgrade_geno_meta())
    recoded <- .st_to_nucleotide(vals)
  }

  # ---- Render missing, if asked ---------------------------------------------
  n_rendered <- 0L
  if (!is.null(na_as)) {
    n_rendered <- sum(is.na(recoded))
    recoded[is.na(recoded)] <- na_as
  }

  arr2 <- array(recoded, dim = dim(arr), dimnames = dimnames(arr))

  # Rebuild char_dict for the new encoding from the alleles actually present
  present   <- unique(stats::na.omit(as.character(arr2)))
  if (!is.null(na_as)) present <- setdiff(present, na_as)
  keep      <- intersect(names(DEFAULT_ENCODING_DICT), present)
  char_dict <- unlist(DEFAULT_ENCODING_DICT[keep])

  if (n_rendered > 0L)
    .hs_inform("recode_geno_meta: rendered ", n_rendered,
            " missing allele(s) as \"", na_as, "\". The result is no longer a ",
            "normalised container \u2014 repair_geno_meta() undoes the rendering, ",
            "but the source file's own missing code is not recoverable.")

  # One rebuild path (0.9.0): blocks, FAM, effects, phenotype and annotation
  # are forwarded; the vocabulary keys are patched, and the old stamp is
  # dropped so the portal restamps against the target vocabulary.
  out <- .rebuild_geno_meta(
    geno_meta,
    geno_arrays   = arr2,
    metadata      = list(char_dict = char_dict, encoding = to_type,
                         missing_tokens = NULL),
    source_format = to_type,
    .caller       = "recode_geno_meta")

  # A deliberate rendering must survive the portal it just passed through.
  if (!is.null(na_as)) {
    out$geno_arrays <- arr2
    out$metadata$missing_tokens <- NULL
  }
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
    .hs_abort("recode_geno_meta: source declared 'selectiontools' but the container ",
         "holds alphabetic allele(s): ",
         paste0("\"", utils::head(present[alpha], 5L), "\"", collapse = ", "),
         ". Pass from_type = \"nucleotide\", or fix metadata$encoding.")
  if (identical(from_type, "nucleotide") && !any(alpha))
    .hs_abort("recode_geno_meta: source declared 'nucleotide' but no allele in the ",
         "container is alphabetic (present: ",
         paste0("\"", utils::head(present, 5L), "\"", collapse = ", "),
         "). Pass from_type = \"selectiontools\", or fix metadata$encoding.")
  invisible(TRUE)
}


# ── Internal builders ─────────────────────────────────────────────────────────
