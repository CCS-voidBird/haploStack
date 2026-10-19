# =============================================================================
# annotation_compute.R
#
# BAN (thread #033) stage 2 + 3 — the R-native local block annotation:
#
#   compute_block_annotation() — genes (+ optional xrefs, QTLs) x the
#                                container's blocks; object in, object out,
#                                writes nothing
#   attach_block_annotation()  — store the result as meta$annotation (+
#                                metadata$annotation provenance)
#   export_block_annotation()  — the Python block_annotation step1 / step3 /
#                                step3_records / step4 TSVs, byte-comparable
#                                at flank_bp = 0 (Q25.c)
#
# The reference semantics are block_annotation 0.3.0's local route
# (annotate_blocks_local(), _marker_gene_overlaps(),
# annotate_features_deep_local()); every deliberate difference is an argument
# recorded in provenance (flank_bp — Q25.d-iv) and messaged at export.
#
# No plotting in this module (owner clarification 2026-10-19): results are
# tables; the catalogue layer consumes them.
# =============================================================================


# The shared default for the annotation matching thresholds (Q25.d-iv, owner
# insertion 2026-10-19): `flank_bp` here and `match_bp` in the catalogue-side
# overlay (0.8.1) both read this constant.
.HS_ANNOT_DEFAULT_BP <- 10000L

# Internal: integer-safe column for export — barley/wheat bp fit in int32; a
# value beyond it stays numeric (fwrite prints both without scientific noise).
.hs_int_col <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x)) || max(abs(x), na.rm = TRUE) < .Machine$integer.max)
    as.integer(round(x)) else round(x)
}

# Internal: `gene:HORVU...` -> `HORVU...` — the xref join key. Identical to the
# Python _strip_feature_prefix(): everything up to the FIRST colon goes.
.hs_strip_prefix <- function(x) {
  has <- grepl(":", x, fixed = TRUE)
  x[has] <- sub("^[^:]*:", "", x[has])
  x
}

# Internal: apply a container-chromosome -> annotation-vocabulary map.
# Unmapped names pass through unchanged (the Python _map_chrom() fallback).
.hs_map_chrom <- function(chrom, chr_map) {
  if (is.null(chr_map)) return(chrom)
  if (is.null(names(chr_map)) || any(!nzchar(names(chr_map))))
    stop("'chr_map' must be a fully named character vector ",
         "(container name -> annotation seqid).")
  hit <- chrom %in% names(chr_map)
  chrom[hit] <- unname(chr_map[as.character(chrom[hit])])
  chrom
}


#' Annotate the container's blocks from local gene / QTL tables
#'
#' The R-native offline annotation route (thread #033): every GFF3 gene (and
#' optionally every QTL interval) is matched to every block of the
#' \code{GenoMeta} by interval overlap, and every (marker, gene)
#' co-localisation inside a block is tabulated — the same three set operations
#' as \code{block_annotation}'s \code{annotate_blocks_local()} +
#' \code{annotate_features_deep_local()}, computed with
#' \code{data.table::foverlaps()} on the tables the container already holds.
#' Object in, object out: nothing is written (see
#' \code{\link{export_block_annotation}}) and the container is not touched
#' (see \code{\link{attach_block_annotation}}).
#'
#' @section Overlap semantics:
#' Coordinates are 1-based inclusive on both sides (GFF3 and PLINK bp). A gene
#' hits a block iff \code{gene_end >= Start - flank_bp} and
#' \code{gene_start <= End + flank_bp}. A (marker, gene) pair requires the
#' marker inside the gene \emph{and} inside its own (padded) block — the gate
#' that stops a gene straddling several blocks from claiming markers it does
#' not contain. \code{flank_bp} widens the overlap window only: the reported
#' block spans are never padded.
#'
#' \strong{The default \code{flank_bp = 10000} differs from the Python route
#' (which has no padding).} It exists because single-marker blocks are 1 bp
#' wide and would otherwise rarely touch a gene; run with \code{flank_bp = 0}
#' for results byte-comparable to \code{annotate_blocks_local()}. The active
#' flank is recorded in provenance and messaged at export.
#'
#' @section Chromosome vocabulary:
#' \code{chr_map} is the only mapping applied (no species table): names are the
#' \emph{container's} chromosome spellings, values the annotation files'
#' seqids (e.g. \code{c("1" = "1H")}). When no block overlaps any gene the
#' warning prints both vocabularies side by side — the most common silent
#' failure of the Python bridge.
#'
#' @param meta      A valid \code{GenoMeta} with blocks attached
#'   (\code{\link{attach_block_info}}). \code{map_info$Pos} must be bp.
#' @param genes     Gene table from \code{\link{read_gff3_genes}} (or any
#'   data.frame with \code{chrom, start, end, strand, gene_id, gene_name,
#'   biotype}).
#' @param xrefs     Optional per-gene annotation from
#'   \code{\link{read_gene_xrefs}}; joined onto the (marker, gene) table.
#' @param qtl       Optional QTL table: \code{chrom, start, end, name} (the
#'   4-column BED layout; a path is read via the internal BED reader).
#'   Coordinates are used exactly as given, as the Python route does.
#' @param chr_map   Named character vector, container chromosome ->
#'   annotation seqid. \code{NULL} when the vocabularies already agree.
#' @param flank_bp  Matching threshold in bp (default \code{10000}; see
#'   \emph{Overlap semantics}). \code{0} reproduces the Python route exactly.
#' @param strip_prefix Strip the GFF3 type prefix (\code{gene:}) from gene ids
#'   \strong{at the xref join only} (default \code{TRUE}, the Python
#'   behaviour). Stored ids always keep the file's spelling.
#' @param browser_base Optional Ensembl browser base URL (e.g.
#'   \code{"https://plants.ensembl.org/Hordeum_vulgare"}); when given, the
#'   (marker, gene) table carries
#'   \code{<base>/Location/View?r=<chrom>\%3A<start>-<end>} links built with
#'   the mapped chromosome, as the Python route does via its species config.
#'   \code{NULL} leaves the column empty.
#'
#' @return A list of class-less tables (haploStack column convention):
#' \describe{
#'   \item{\code{blocks}}{one row per block: \code{Chrom, Block, Start, End,
#'     Num_SNPs, Block_length, gene_count, gene_ids, gene_names} and, with
#'     \code{qtl}, \code{qtl_count, qtl_names}. Joined strings use
#'     \code{"; "}, empty-for-absent, exactly as the Python step3.}
#'   \item{\code{features}}{one row per (block, gene): \code{Block, Chrom,
#'     Block_start, Block_end, Num_SNPs, Block_length, feature_type,
#'     feature_id, feature_name, feature_biotype, feature_start, feature_end,
#'     feature_strand}.}
#'   \item{\code{marker_gene}}{one row per (marker, gene) pair: \code{Marker,
#'     Chrom, Pos, Block, gene_id, gene_symbol, gene_name, gene_biotype,
#'     gene_start, gene_end, gene_strand, gene_description, go_terms,
#'     interpro_domains, ensembl_browser_url}.}
#'   \item{\code{provenance}}{inputs, counts, \code{chr_map}, \code{flank_bp},
#'     \code{strip_prefix}, date and package version — written to
#'     \code{metadata$annotation} by \code{\link{attach_block_annotation}}.}
#' }
#'
#' @seealso \code{\link{read_gff3_genes}}, \code{\link{read_gene_xrefs}},
#'   \code{\link{attach_block_annotation}},
#'   \code{\link{export_block_annotation}}
#' @export
compute_block_annotation <- function(meta, genes, xrefs = NULL, qtl = NULL,
                                     chr_map = NULL,
                                     flank_bp = .HS_ANNOT_DEFAULT_BP,
                                     strip_prefix = TRUE,
                                     browser_base = NULL) {
  if (!is_geno_meta(meta))
    stop("compute_block_annotation: 'meta' is not a valid GenoMeta.")
  if (is.null(meta$block_marker_map_info))
    stop("compute_block_annotation: the GenoMeta has no blocks — call ",
         "attach_block_info() first.")
  pu <- meta$metadata$pos_unit
  if (!is.null(pu) && !identical(tolower(as.character(pu)), "bp"))
    stop("compute_block_annotation: map_info$Pos is stored in '", pu,
         "' but annotation coordinates are bp — convert the container ",
         "positions to bp before annotating.")
  if (!is.numeric(flank_bp) || length(flank_bp) != 1L || is.na(flank_bp) ||
      flank_bp < 0)
    stop("compute_block_annotation: 'flank_bp' must be a single number >= 0.")
  flank_bp <- as.numeric(flank_bp)

  # --- genes table -----------------------------------------------------------
  g <- as.data.table(genes)
  need_g <- c("chrom", "start", "end", "gene_id")
  miss <- setdiff(need_g, names(g))
  if (length(miss) > 0L)
    stop("compute_block_annotation: 'genes' is missing column(s): ",
         paste(miss, collapse = ", "), " — use read_gff3_genes().")
  for (cc in c("gene_name", "biotype", "strand"))
    if (!cc %in% names(g)) g[, (cc) := NA_character_]
  g <- g[, .(chrom = as.character(chrom), start = as.numeric(start),
             end = as.numeric(end), strand = as.character(strand),
             gene_id = as.character(gene_id),
             gene_name = as.character(gene_name),
             biotype = as.character(biotype))]
  g[, gidx := .I]

  # --- block spans (recomputed from marker positions, the step-1 contract) ---
  loc <- build_block_location_info(meta)
  bl  <- as.data.table(loc$block_info)[
    , .(Chrom = as.character(Chrom), Block = as.character(Block),
        Start = as.numeric(Start), End = as.numeric(End),
        Num_SNPs = as.integer(Num_SNPs),
        Block_length = as.numeric(Block_length))]
  # Python step order: sort_values(["chrom", "start_pos"]) — string chrom.
  setorder(bl, Chrom, Start)
  bl[, bidx := .I]
  bl[, MChrom := .hs_map_chrom(Chrom, chr_map)]

  # --- gene x block overlap (flank pads the query window only) ---------------
  gk <- data.table::copy(g)
  setkeyv(gk, c("chrom", "start", "end"))
  qb <- bl[, .(MChrom, qs = Start - flank_bp, qe = End + flank_bp, bidx)]
  ov <- foverlaps(qb, gk, by.x = c("MChrom", "qs", "qe"),
                  type = "any", nomatch = NULL)
  setorder(ov, bidx, gidx)

  if (nrow(ov) == 0L)
    warning("compute_block_annotation: no block overlaps any gene. ",
            "Block chromosomes (mapped): ",
            paste(sort(unique(bl$MChrom)), collapse = ", "),
            " | gene seqids: ",
            paste(utils::head(sort(unique(g$chrom)), 12L), collapse = ", "),
            if (length(unique(g$chrom)) > 12L) " ..." else "",
            ". Pass 'chr_map' to translate the container's chromosome ",
            "names onto the annotation's.")

  # Joined strings follow the Python "; ".join over attrs.get(..., "")
  # exactly: an absent id/name joins as the empty string.
  blank <- function(x) data.table::fifelse(is.na(x), "", x)
  agg <- ov[, .(gene_count = .N,
                gene_ids   = paste(blank(gene_id),   collapse = "; "),
                gene_names = paste(blank(gene_name), collapse = "; ")),
            by = bidx]
  blocks <- data.table::copy(bl)
  blocks[, `:=`(gene_count = 0L, gene_ids = "", gene_names = "")]
  blocks[agg, on = "bidx",
         `:=`(gene_count = i.gene_count, gene_ids = i.gene_ids,
              gene_names = i.gene_names)]

  # --- features: one row per (block, gene) -----------------------------------
  features <- ov[, .(bidx, gidx,
                     Block = bl$Block[bidx], Chrom = bl$Chrom[bidx],
                     MChrom = bl$MChrom[bidx],
                     Block_start = bl$Start[bidx], Block_end = bl$End[bidx],
                     Num_SNPs = bl$Num_SNPs[bidx],
                     Block_length = bl$Block_length[bidx],
                     feature_type = "gene",
                     feature_id = gene_id, feature_name = gene_name,
                     feature_biotype = biotype,
                     feature_start = start, feature_end = end,
                     feature_strand = strand)]
  features[, ridx := .I]

  # --- optional QTL x block overlap ------------------------------------------
  has_qtl <- !is.null(qtl)
  if (has_qtl) {
    q <- if (is.character(qtl) && length(qtl) == 1L) .read_qtl_bed(qtl)
         else as.data.table(qtl)
    miss_q <- setdiff(c("chrom", "start", "end", "name"), names(q))
    if (length(miss_q) > 0L)
      stop("compute_block_annotation: 'qtl' is missing column(s): ",
           paste(miss_q, collapse = ", "))
    q <- q[, .(chrom = as.character(chrom), start = as.numeric(start),
               end = as.numeric(end), name = as.character(name))]
    q[, qidx := .I]
    setkeyv(q, c("chrom", "start", "end"))
    ovq <- foverlaps(qb, q, by.x = c("MChrom", "qs", "qe"),
                     type = "any", nomatch = NULL)
    setorder(ovq, bidx, qidx)
    # The Python qtl_names drops a missing name from the join but still
    # counts the hit (pandas dropna on a real NaN).
    aggq <- ovq[, .(qtl_count = .N,
                    qtl_names = paste(name[!is.na(name)], collapse = "; ")),
                by = bidx]
    blocks[, `:=`(qtl_count = 0L, qtl_names = "")]
    blocks[aggq, on = "bidx",
           `:=`(qtl_count = i.qtl_count, qtl_names = i.qtl_names)]
  }

  # --- (marker, gene) pairs: marker in gene AND in its own padded block ------
  mk <- as.data.table(meta$map_info)[
    , .(Marker = as.character(Marker), Chrom = as.character(Chrom),
        Pos = as.numeric(Pos))]
  mk[, midx := .I]
  mk[, MChrom := .hs_map_chrom(Chrom, chr_map)]
  fk <- features[, .(ridx, MChrom, feature_start, feature_end)]
  setkeyv(fk, c("MChrom", "feature_start", "feature_end"))
  qm <- mk[, .(MChrom, ms = Pos, me = Pos, midx)]
  ovm <- foverlaps(qm, fk, by.x = c("MChrom", "ms", "me"),
                   type = "any", nomatch = NULL)
  pairs <- cbind(mk[ovm$midx,
                    .(Marker, Chrom, Pos, midx)],
                 features[ovm$ridx,
                          .(ridx, Block, Block_start, Block_end,
                            gene_id = feature_id, gene_name = feature_name,
                            gene_biotype = feature_biotype,
                            gene_start = feature_start,
                            gene_end = feature_end,
                            gene_strand = feature_strand)])
  pairs <- pairs[Pos >= Block_start - flank_bp & Pos <= Block_end + flank_bp]
  # Python order: chromosomes ascending as strings (pandas groupby), then the
  # marker map's row order, then the feature-record order.
  pairs[, chr_key := match(Chrom, sort(unique(mk$Chrom)))]
  setorder(pairs, chr_key, midx, ridx)
  pairs[, chr_key := NULL]

  # xref enrichment on the (optionally) stripped id — never on the stored id.
  key <- if (isTRUE(strip_prefix)) .hs_strip_prefix(pairs$gene_id)
         else pairs$gene_id
  ann_cols <- c("gene_symbol", "gene_description", "go_terms",
                "interpro_domains")
  if (!is.null(xrefs)) {
    xr <- as.data.table(xrefs)
    if (!"gene_id" %in% names(xr))
      stop("compute_block_annotation: 'xrefs' needs a gene_id column — ",
           "use read_gene_xrefs().")
    idx <- match(key, xr$gene_id)
    for (cc in ann_cols) {
      v <- if (cc %in% names(xr)) as.character(xr[[cc]])[idx]
           else NA_character_
      pairs[, (cc) := data.table::fifelse(is.na(v), "", v)]
    }
  } else {
    for (cc in ann_cols) pairs[, (cc) := ""]
  }

  # Browser URL (the Python route builds these from its species config; here
  # the base is explicit — Q25.e, no species table).
  pairs[, ensembl_browser_url :=
          if (!is.null(browser_base) && nzchar(browser_base))
            paste0(browser_base, "/Location/View?r=",
                   .hs_map_chrom(Chrom, chr_map), "%3A",
                   .hs_int_col(gene_start), "-", .hs_int_col(gene_end))
          else ""]

  marker_gene <- pairs[, .(Marker, Chrom, Pos, Block, gene_id, gene_symbol,
                           gene_name, gene_biotype, gene_start, gene_end,
                           gene_strand, gene_description, go_terms,
                           interpro_domains, ensembl_browser_url)]

  src_of <- function(x) {
    s <- attr(x, "hs_annotation_source", exact = TRUE)
    if (is.null(s)) list(path = NA_character_, size = NA_real_,
                         mtime = NA_character_, n = nrow(as.data.frame(x)))
    else s
  }
  provenance <- list(
    gff3            = src_of(genes),
    n_genes         = nrow(g),
    xref            = if (is.null(xrefs)) NULL else src_of(xrefs),
    qtl             = if (!has_qtl) NULL
                      else list(n = nrow(q),
                                path = if (is.character(qtl)) qtl
                                       else NA_character_),
    chr_map         = chr_map,
    flank_bp        = flank_bp,
    strip_prefix    = isTRUE(strip_prefix),
    browser_base    = browser_base,
    n_blocks        = nrow(blocks),
    n_blocks_hit    = sum(blocks$gene_count > 0L),
    n_features      = nrow(features),
    n_marker_gene   = nrow(marker_gene),
    date            = format(Sys.Date()),
    package_version = as.character(utils::packageVersion("haploStack"))
  )

  blocks_out <- blocks[, c("Chrom", "Block", "Start", "End", "Num_SNPs",
                           "Block_length", "gene_count", "gene_ids",
                           "gene_names",
                           if (has_qtl) c("qtl_count", "qtl_names")),
                       with = FALSE]
  features_out <- features[, .(Block, Chrom, Block_start, Block_end, Num_SNPs,
                               Block_length, feature_type, feature_id,
                               feature_name, feature_biotype, feature_start,
                               feature_end, feature_strand)]

  list(blocks      = blocks_out,
       features    = features_out,
       marker_gene = marker_gene,
       provenance  = provenance)
}


# Internal: is this a compute_block_annotation() result?
.hs_is_annotation <- function(x) {
  is.list(x) && all(c("blocks", "features", "marker_gene") %in% names(x)) &&
    is.data.frame(x$blocks)
}


#' Store a block annotation on the GenoMeta
#'
#' Writes the \code{\link{compute_block_annotation}} result onto the container
#' as the optional element \code{meta$annotation} (tables \code{blocks},
#' \code{features}, \code{marker_gene}) and its provenance as
#' \code{metadata$annotation}. Every \code{Block} must exist in the
#' container's block vocabulary and every \code{Marker} in \code{map_info} —
#' an annotation computed against different blocks is refused.
#'
#' Carry-forward (Q25.i): \code{\link{filter_geno_meta}} prunes the stored
#' annotation to the surviving blocks and markers;
#' \code{\link{attach_block_info}} \strong{drops} it with a message (new
#' blocks invalidate every overlap); phenotype / signal-marker / recode
#' transforms and the simulation's offspring builders keep it.
#'
#' @param meta       A valid \code{GenoMeta}.
#' @param annotation The list returned by
#'   \code{\link{compute_block_annotation}}.
#'
#' @return The \code{GenoMeta} with \code{$annotation} set. The input object
#'   is not modified.
#'
#' @seealso \code{\link{compute_block_annotation}},
#'   \code{\link{export_block_annotation}}
#' @export
attach_block_annotation <- function(meta, annotation) {
  if (!is_geno_meta(meta))
    stop("attach_block_annotation: 'meta' is not a valid GenoMeta.")
  if (!.hs_is_annotation(annotation))
    stop("attach_block_annotation: 'annotation' must be the list returned by ",
         "compute_block_annotation().")

  known_blocks <- unique(as.character(
    if (!is.null(meta$block_marker_map_info)) meta$block_marker_map_info$Block
    else meta$block_info$Block))
  bad_b <- setdiff(unique(as.character(annotation$blocks$Block)), known_blocks)
  if (length(bad_b) > 0L)
    stop("attach_block_annotation: ", length(bad_b), " annotated block(s) ",
         "not in the container's blocks (first: ",
         paste(utils::head(bad_b, 5L), collapse = ", "),
         ") — the annotation was computed against different block ",
         "definitions.")
  bad_m <- setdiff(unique(as.character(annotation$marker_gene$Marker)),
                   meta$map_info$Marker)
  if (length(bad_m) > 0L)
    stop("attach_block_annotation: ", length(bad_m), " annotated marker(s) ",
         "not in map_info (first: ",
         paste(utils::head(bad_m, 5L), collapse = ", "), ").")

  # Stored as plain data.frames — the container's own convention (map_info,
  # block_info are data.frames), and a data.table's .internal.selfref pointer
  # would break the save/load identical() contract (stage-3 acceptance).
  meta$annotation <- lapply(annotation[c("blocks", "features", "marker_gene")],
                            function(x) as.data.frame(x,
                                                      stringsAsFactors = FALSE))
  meta$metadata$annotation <- annotation$provenance
  meta
}


# Internal: prune a stored annotation to surviving blocks / markers
# (filter_geno_meta()). Returns NULL for NULL in.
.hs_prune_annotation <- function(ann, keep_blocks, keep_markers) {
  if (is.null(ann)) return(NULL)
  ann$blocks   <- ann$blocks[as.character(ann$blocks$Block) %in% keep_blocks, ,
                             drop = FALSE]
  ann$features <- ann$features[as.character(ann$features$Block) %in%
                                 keep_blocks, , drop = FALSE]
  ann$marker_gene <- ann$marker_gene[
    as.character(ann$marker_gene$Marker) %in% keep_markers &
      as.character(ann$marker_gene$Block) %in% keep_blocks, , drop = FALSE]
  ann
}


#' Export a block annotation in the Python step-file schema
#'
#' Writes the \code{block_annotation} pipeline's step files from an R-side
#' annotation, so every figure, join and document written against the Python
#' TSVs works unchanged (Q25.c):
#' \itemize{
#'   \item \code{step1_block_coordinates.tsv} — \code{chrom, block, start_pos,
#'     end_pos, n_markers, length_bp}
#'   \item \code{step1_block_marker_map.tsv} — \code{chrom, marker, block}
#'     (only when \code{x} is a \code{GenoMeta}; an annotation list does not
#'     carry the marker map)
#'   \item \code{step3_local_annotation.tsv} — one row per block
#'   \item \code{step3_local_annotation_records.tsv} — one row per
#'     (block, gene)
#'   \item \code{step4_deep_annotation.tsv} — one row per (marker, gene)
#' }
#'
#' At \code{flank_bp = 0} the output is byte-comparable to
#' \code{annotate_blocks_local()} / \code{annotate_features_deep_local()} on
#' the same inputs; a non-zero flank changes the hit sets, so exporting one
#' prints a message naming it (the schema is unchanged).
#'
#' @param x          A \code{GenoMeta} carrying \code{$annotation}, or the
#'   list returned by \code{\link{compute_block_annotation}}.
#' @param output_dir Directory for the TSVs (created when missing).
#' @param prefix     Optional filename prefix (default \code{NULL} writes the
#'   Python filenames exactly).
#'
#' @return Invisibly, the written paths (named character vector).
#'
#' @seealso \code{\link{compute_block_annotation}},
#'   \code{\link{attach_block_annotation}}
#' @export
export_block_annotation <- function(x, output_dir, prefix = NULL) {
  ann <- if (.hs_is_annotation(x)) x
         else if (is.list(x) && !is.null(x$annotation)) {
           a <- x$annotation
           a$provenance <- x$metadata$annotation
           a
         } else
           stop("export_block_annotation: 'x' must be a ",
                "compute_block_annotation() result or a GenoMeta carrying ",
                "$annotation (attach_block_annotation()).")
  from_meta <- !.hs_is_annotation(x)

  if (missing(output_dir) || !is.character(output_dir))
    stop("export_block_annotation: 'output_dir' is required.")
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  flank <- ann$provenance$flank_bp
  if (!is.null(flank) && is.numeric(flank) && flank != 0)
    message("export_block_annotation: this annotation was computed with ",
            "flank_bp = ", format(flank, big.mark = ","), " — hit counts are ",
            "not comparable to a Python block_annotation run, which has no ",
            "padding. Recompute with flank_bp = 0 for byte-level parity.")

  pfx <- if (is.null(prefix) || !nzchar(prefix)) "" else paste0(prefix, "_")
  path_of <- function(nm) file.path(output_dir, paste0(pfx, nm))
  wr <- function(dt, nm) {
    p <- path_of(nm)
    data.table::fwrite(dt, p, sep = "\t", quote = FALSE, na = "", eol = "\n")
    p
  }

  bl <- as.data.table(ann$blocks)
  out <- c()

  step1 <- bl[, .(chrom = Chrom, block = Block,
                  start_pos = .hs_int_col(Start), end_pos = .hs_int_col(End),
                  n_markers = .hs_int_col(Num_SNPs),
                  length_bp = .hs_int_col(Block_length))]
  out["step1_block_coordinates"] <- wr(step1, "step1_block_coordinates.tsv")

  if (from_meta && !is.null(x$block_marker_map_info)) {
    bmm <- as.data.table(x$block_marker_map_info)
    out["step1_block_marker_map"] <- wr(
      bmm[, .(chrom = as.character(Chrom), marker = as.character(Marker),
              block = as.character(Block))],
      "step1_block_marker_map.tsv")
  } else if (!from_meta) {
    message("export_block_annotation: no GenoMeta given — ",
            "step1_block_marker_map.tsv skipped (the annotation list does ",
            "not carry the marker map).")
  }

  has_qtl <- all(c("qtl_count", "qtl_names") %in% names(bl))
  step3 <- bl[, c(list(block = Block, chrom = Chrom,
                       start_pos = .hs_int_col(Start),
                       end_pos = .hs_int_col(End),
                       n_markers = .hs_int_col(Num_SNPs),
                       length_bp = .hs_int_col(Block_length),
                       gene_count = .hs_int_col(gene_count),
                       gene_ids = gene_ids, gene_names = gene_names),
                  if (has_qtl) list(qtl_count = .hs_int_col(qtl_count),
                                    qtl_names = qtl_names))]
  out["step3_local_annotation"] <- wr(step3, "step3_local_annotation.tsv")

  fe <- as.data.table(ann$features)
  recs <- fe[, .(block = Block, chrom = Chrom,
                 block_start = .hs_int_col(Block_start),
                 block_end = .hs_int_col(Block_end),
                 n_markers = .hs_int_col(Num_SNPs),
                 block_length_bp = .hs_int_col(Block_length),
                 feature_type, feature_id,
                 feature_name = data.table::fifelse(is.na(feature_name), "",
                                                    feature_name),
                 feature_biotype = data.table::fifelse(is.na(feature_biotype),
                                                       "", feature_biotype),
                 feature_start = .hs_int_col(feature_start),
                 feature_end = .hs_int_col(feature_end),
                 feature_strand = data.table::fifelse(is.na(feature_strand),
                                                      "", feature_strand),
                 ensembl_browser_url = "")]
  out["step3_local_annotation_records"] <-
    wr(recs, "step3_local_annotation_records.tsv")

  mg <- as.data.table(ann$marker_gene)
  step4 <- mg[, .(marker_name = Marker, marker_chr = Chrom,
                  marker_pos = .hs_int_col(Pos), block = Block, gene_id,
                  gene_symbol,
                  gene_name = data.table::fifelse(is.na(gene_name), "",
                                                  gene_name),
                  gene_biotype = data.table::fifelse(is.na(gene_biotype), "",
                                                     gene_biotype),
                  gene_start = .hs_int_col(gene_start),
                  gene_end = .hs_int_col(gene_end),
                  gene_strand = data.table::fifelse(is.na(gene_strand), "",
                                                    gene_strand),
                  gene_description, go_terms, interpro_domains,
                  ensembl_browser_url)]
  out["step4_deep_annotation"] <- wr(step4, "step4_deep_annotation.tsv")

  invisible(out)
}
