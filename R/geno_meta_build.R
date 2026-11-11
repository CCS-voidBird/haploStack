# =============================================================================
# geno_meta_build.R — build_geno_meta(): the one building route
#
# A genotype matrix plus a map, with optional FAM, blocks, effects and
# phenotype, become a GenoMeta through new_geno_meta() (geno_meta_core.R).
# The readers end here; the 0.9.1 facade (P2) will route the converters here
# too. Split out of format_helper.R in 0.9.0 (D2-9).
# =============================================================================

#' Build a GenoMeta object from fixed-shape in-memory variables
#'
#' Assembles a valid \code{GenoMeta} from three primary elements with a
#' \strong{fixed input contract}: a genotype matrix, a map, and an optional fam
#' table. This is the single assembly entry point shared by
#' \code{\link{read_genotype}} and by callers holding in-memory data. The file
#' layout — orientation, headers, embedded FAM columns, map column order — is
#' resolved \emph{before} this point (in the reader); \code{build_geno_meta()}
#' itself takes only the changeable \emph{format} parameters \code{split_char},
#' \code{pos_unit}, and \code{encoding}.
#'
#' @section Fixed input contract:
#' \describe{
#'   \item{\code{geno_2d}}{Genotype matrix in \code{(n_samples × n_markers)}
#'     orientation only — rows are samples, columns are markers. Row names must
#'     be sample IDs and column names must be marker IDs. Orientation is
#'     \strong{not} auto-detected.}
#'   \item{\code{map}}{At least three columns in the order
#'     \code{Chrom, Marker, Pos}.}
#'   \item{\code{fam}}{Columns in the order \code{Family, Sample, Parent1,
#'     Parent2, Sex} (the \code{Sample} column, position 2, identifies samples).}
#' }
#'
#' Marker and sample name intersections are resolved internally with informative
#' warnings when names are dropped.
#'
#' @param geno_2d      Character matrix or data.frame of raw genotype calls in
#'   \code{(n_samples × n_markers)} orientation with named rows (samples) and
#'   columns (markers). Cells may be slash-separated haplotypes (\code{"1/4"}),
#'   two-character concatenations (\code{"AA"}), or single codes.
#' @param map          data.frame whose first three columns are, in order,
#'   \code{Chrom}, \code{Marker}, \code{Pos}.
#' @param fam          Optional fam data.frame with columns in the order
#'   \code{Family}, \code{Sample}, \code{Parent1}, \code{Parent2}, \code{Sex}.
#'   The \code{Sample} column (position 2) supplies the sample IDs; a
#'   single-column fam is treated as the sample list. Samples absent from
#'   \code{fam} are dropped with a warning. \code{NULL} (default) keeps all
#'   samples and builds a dummy fam.
#' @param blocks       Optional block map passed to
#'   \code{\link{attach_block_info}} once the container is built: a path or a
#'   data.frame in one of the two shapes it accepts (block-marker map
#'   \code{Chrom, Marker, Block} or block-region map \code{Chrom, Block, Start,
#'   End}; header required). \code{NULL} (default) skips block attachment.
#' @param block_overlap Overlap policy for \code{blocks}, as
#'   \code{attach_block_info(overlap =)}: \code{"error"} (default) or
#'   \code{"first"}.
#' @param ploidy       Integer ploidy. Default \code{2L}.
#' @param split_char   Allele separator used in genotype cells: \code{"/"}
#'   (default) or \code{""} (one character per allele).
#' @param pos_unit     Position unit of \code{map$Pos}: \code{"bp"} (default) or
#'   \code{"mbp"}. Internally all positions are stored in bp.
#' @param encoding     Allele encoding: \code{"nucleotide"} (categorical allele
#'   labels — \code{A/C/G/T}, or digits used as labels such as PLINK's 1/2
#'   coding) or \code{"numeric"} (dosage). \code{NULL} (default) auto-detects,
#'   except that integer codes are refused with a pointer to the two
#'   declarations: \code{"selectiontools"} (deprecated since 0.9.0; the codes
#'   \code{1/2/3/4} are translated to \code{A/C/G/T} and the container is
#'   stamped \code{"nucleotide"}) or \code{"nucleotide"} (the digits are kept as
#'   labels). SelectionTools' files go through \code{\link{st_convert}}.
#' @param source_format Label stored in \code{metadata$source_format}.
#'   Default \code{"user_supplied"}.
#' @param extra_na     Optional character vector of additional missing
#'   sentinels this file uses, passed through to the normalisation portal in
#'   \code{\link{new_geno_meta}}. This is the per-read escape hatch for a source
#'   whose absence code the package does not know — including a dosage panel
#'   that writes a digit such as \code{"9"} for missing, which
#'   \code{\link{use_missing_chars}} deliberately refuses to register globally.
#' @param missing_chars Extra missing-allele characters for this build:
#'   \code{NULL} (default) takes the session registry
#'   (\code{\link{use_missing_chars}}); a character vector names them for the
#'   resolved encoding directly (\code{character(0)} ignores the registry).
#'   See \code{\link{new_geno_meta}}.
#'
#' @section Building from a dosage source:
#' A dosage matrix (0/1/2 counts, \code{encoding = "numeric"}) records \emph{how
#' many} copies of an allele a sample carries, not \emph{which} copy sits on
#' which chromosome. The container built from one therefore carries \strong{no
#' strand resolution}: per-slot output from it is not a haplotype, and
#' \code{level = "haplotype"} in \code{\link{compute_profile_localgebv}} should
#' not be read as one.
#'
#' Two consequences worth knowing before you build:
#'
#' \itemize{
#'   \item A dosage source built at \code{ploidy = 2} produces a second ploidy
#'     layer that is entirely \code{NA}, so \code{\link{compute_marker_stats}}
#'     reports \code{Call_rate} 0.5 on every marker and
#'     \code{\link{detect_genotype_encoding}} fails on the reconstructed allele
#'     matrix. The local-GEBV path survives it (the absent layer contributes 0
#'     at genotype level and is dropped at haplotype level), but the container
#'     is describing half its data as missing. Tracked as INDEX discrepancy 25;
#'     fixing the shape changes \code{Call_rate} on every dosage panel and is
#'     scoped separately.
#'   \item For LD estimation and block definition on a dosage frame, prefer
#'     \strong{HapSelect}, which works on dosage natively and never claims
#'     strands, over converting to alleles first.
#' }
#'
#' @return A valid \code{GenoMeta} list (see \code{\link{new_geno_meta}}).
#'
#' @details
#' \strong{Automatic marker filtering.} \code{build_geno_meta()} intersects the
#' genotype markers with \code{map$Marker} and keeps only the shared set — no
#' manual pre-filtering is needed. Genotype markers absent from the map are
#' dropped (with a warning), map markers absent from the genotypes are excluded,
#' and a complete mismatch (no shared markers) raises an error. Every reader
#' (\code{\link{read_genotype}} and its wrappers) delegates here, so this is the
#' single place marker intersection happens. An optional \code{fam} table applies
#' the analogous automatic intersection on the sample axis.
#'
#' @examples
#' \dontrun{
#' # geno_2d: rows = samples, cols = markers (fixed orientation)
#' geno <- matrix(c("A/A", "A/T", "T/T",
#'                  "T/T", "A/T", "A/A"),
#'                nrow = 2L, byrow = TRUE,
#'                dimnames = list(c("S1", "S2"), c("M1", "M2", "M3")))
#' # map: columns in order Chrom, Marker, Pos
#' map  <- data.frame(Chrom = c("1", "1", "2"),
#'                    Marker = c("M1", "M2", "M3"),
#'                    Pos = c(1e6, 2e6, 3e6), stringsAsFactors = FALSE)
#' meta <- build_geno_meta(geno, map)
#' is_geno_meta(meta)  # TRUE
#' }
#'
#' @seealso \code{\link{new_geno_meta}}, \code{\link{attach_block_info}},
#'   \code{\link{read_genotype}}, \code{\link{resolve_geno_format}}
#'
#' @export
build_geno_meta <- function(geno_2d,
                             map,
                             fam           = NULL,
                             blocks        = NULL,
                             ploidy        = 2L,
                             split_char    = "/",
                             pos_unit      = "bp",
                             encoding      = NULL,
                             source_format = "user_supplied",
                             extra_na      = NULL,
                             missing_chars = NULL,
                             block_overlap = c("error", "first")) {

  # ---- Validate genotype matrix (FIXED sample × marker orientation) ----------
  if (!is.matrix(geno_2d) && !is.data.frame(geno_2d))
    .hs_abort("build_geno_meta: 'geno_2d' must be a matrix or data.frame in ",
         "(n_samples \u00D7 n_markers) orientation.")
  geno_2d <- as.matrix(geno_2d)
  mode(geno_2d) <- "character"
  if (is.null(rownames(geno_2d)) || is.null(colnames(geno_2d)))
    .hs_abort("build_geno_meta: 'geno_2d' must have row names (sample IDs) and ",
         "column names (marker IDs). Orientation is fixed: rows = samples, ",
         "columns = markers. Re-orient the input before calling.")
  sample_names <- rownames(geno_2d)
  marker_names <- colnames(geno_2d)

  # ---- Validate / normalise map (FIXED column order: Chrom, Marker, Pos) -----
  if (missing(map) || is.null(map))
    .hs_abort("build_geno_meta: 'map' is required (columns in order ",
         "Chrom, Marker, Pos).")
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  if (ncol(map) < 3L)
    .hs_abort("build_geno_meta: 'map' must have at least three columns in the ",
         "order Chrom, Marker, Pos.")
  # Positional contract: the first three columns ARE Chrom, Marker, Pos. Honour
  # an already-canonical naming, otherwise rename the first three by position.
  if (all(c("Chrom", "Marker", "Pos") %in% names(map))) {
    map <- map[, c("Chrom", "Marker", "Pos"), drop = FALSE]
  } else {
    map <- map[, seq_len(3L), drop = FALSE]
    names(map) <- c("Chrom", "Marker", "Pos")
  }
  # Convert position units (→ bp) and emit the Mbp-scale warning (Q5).
  map <- .process_map(map, marker_names = map$Marker, map_preset = NULL,
                      pos_unit = pos_unit)
  map$Chrom  <- as.character(map$Chrom)
  map$Marker <- as.character(map$Marker)

  # ---- Marker intersection (canonical auto-filter: geno ∩ map) ---------------
  # Single source of truth for marker filtering — every reader delegates here,
  # so callers never need to pre-intersect geno_2d and map by hand.
  map_markers <- map$Marker
  extra_geno  <- setdiff(marker_names, map_markers)
  extra_map   <- setdiff(map_markers,  marker_names)
  if (length(extra_geno) > 0L)
    .hs_warn("build_geno_meta: ", length(extra_geno),
            " marker(s) in geno_2d not found in map and will be dropped: ",
            paste(head(extra_geno, 5L), collapse = ", "),
            if (length(extra_geno) > 5L) " ..." else "")
  if (length(extra_map) > 0L)
    .hs_warn("build_geno_meta: ", length(extra_map),
            " marker(s) in map not found in geno_2d (absent from array): ",
            paste(head(extra_map, 5L), collapse = ", "),
            if (length(extra_map) > 5L) " ..." else "")
  keep_mk <- intersect(marker_names, map_markers)
  if (length(keep_mk) == 0L)
    .hs_abort("build_geno_meta: no markers in common between geno_2d and map.")

  # ---- Sample intersection (optional fam) → becomes sample_info --------------
  # The fam is a positional FAM table (Family, Sample, Parent1, Parent2, Sex);
  # the Sample column is read by position (col 2, or col 1 if single-column).
  keep_smp <- sample_names
  fam_si   <- NULL
  if (!is.null(fam)) {
    fam <- as.data.frame(fam, stringsAsFactors = FALSE)
    if (ncol(fam) < 1L)
      .hs_abort("build_geno_meta: 'fam' must have at least one column (Sample).")
    sample_col <- if (ncol(fam) == 1L) 1L else 2L
    fam_ids    <- as.character(fam[[sample_col]])
    extra_smp  <- setdiff(sample_names, fam_ids)
    if (length(extra_smp) > 0L)
      .hs_warn("build_geno_meta: ", length(extra_smp),
              " sample(s) in geno_2d not found in fam and will be dropped: ",
              paste(head(extra_smp, 5L), collapse = ", "),
              if (length(extra_smp) > 5L) " ..." else "")
    keep_smp <- intersect(sample_names, fam_ids)
    if (length(keep_smp) == 0L)
      .hs_abort("build_geno_meta: no samples in common between geno_2d and fam.")
    # Retain fam rows for the kept samples; new_geno_meta() standardises and
    # aligns them to the genotype array order.
    fam_si <- fam[fam_ids %in% keep_smp, , drop = FALSE]
  }

  # ---- Subset and split to 3D array ------------------------------------------
  # Orientation is fixed (rows = samples): always orientation = 1L.
  geno_sub  <- geno_2d[keep_smp, keep_mk, drop = FALSE]
  split_res <- split_haplotypes(geno_sub, split_char = split_char,
                                orientation = 1L, ploidy = as.integer(ploidy),
                                sample_names = keep_smp, marker_names = keep_mk)
  # split_res$array is (n_samples, n_markers, ploidy)

  # ---- Resolve encoding and char_dict ----------------------------------------
  enc <- encoding %||% split_res$encode
  if (is.null(enc)) enc <- "unknown"
  # Integer codes with no declared encoding (0.9.0, C-5): SelectionTools'
  # 1/2/3/4 or PLINK's 1/2 labels — the caller must say which.
  if (is.null(encoding) && (identical(enc, "selectiontools") ||
                            (identical(enc, "unknown") && .st_is_coded(split_res$array))))
    .st_refuse_codes("build_geno_meta")
  # A declared "selectiontools" is behind the same switch as the files (C-6);
  # the portal (new_geno_meta()) would stop too, but this names the caller.
  if (identical(encoding, "selectiontools")) .st_gate("build_geno_meta", per_call = FALSE)
  .hs_inform("build_geno_meta: detected encoding '", enc, "'")
  # Canonical encodings: "nucleotide" (A/T/C/G, or any categorical allele
  # labels) and "numeric" (what detect_genotype_encoding() emits for
  # 0/1/2-looking data, kept as an alias here). A declared "selectiontools" is
  # deprecated and translated at the portal (new_geno_meta()). The legacy
  # "number" / "char" tokens were removed in 0.5.1 — they used to be accepted
  # here and then stored verbatim in metadata$encoding, producing a GenoMeta
  # that geno_meta_views() rejected.
  enc <- .encode_arg(enc, "build_geno_meta",
                     allowed = c("selectiontools", "nucleotide", "numeric"),
                     strict  = FALSE)
  char_dict <- switch(enc,
    selectiontools = ,
    numeric        = as.list(DEFAULT_ENCODING_DICT),
    nucleotide     = setNames(as.list(seq_along(EXTENDED_ALLELE_DICT)),
                               EXTENDED_ALLELE_DICT),
    {
      # Derive from observed alleles
      obs <- as.vector(split_res$array)
      obs <- unique(obs[!is.na(obs)])
      setNames(as.list(seq_along(obs)), obs)
    }
  )

  # ---- Subset and reorder map to match array marker order --------------------
  map_sub <- map[map$Marker %in% keep_mk, , drop = FALSE]
  map_sub <- map_sub[match(split_res$markers, map_sub$Marker), , drop = FALSE]
  rownames(map_sub) <- NULL

  # ---- Build base GenoMeta ---------------------------------------------------
  meta <- new_geno_meta(
    geno_arrays = split_res$array,
    map_info    = map_sub,
    metadata    = list(
      ploidy        = as.integer(ploidy),
      char_dict     = char_dict,
      n_samples     = length(split_res$samples),
      n_markers     = length(split_res$markers),
      source_format = source_format,
      encoding      = enc
    ),
    sample_info   = fam_si,  # NULL → dummy FAM auto-built in new_geno_meta()
    # `enc` is the vocabulary this call just resolved — the caller's
    # declaration at construction time, which is what the portal keys on.
    source_format = enc,
    extra_na      = extra_na,
    missing_chars = missing_chars
  )

  # ---- Attach blocks if provided ---------------------------------------------
  if (!is.null(blocks))
    meta <- attach_block_info(meta, blocks, overlap = block_overlap)

  meta
}


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------
