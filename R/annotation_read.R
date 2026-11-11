# =============================================================================
# annotation_read.R
#
# BAN (thread #033) stage 1 — readers for the R-native local block annotation.
# Three table producers, no plotting anywhere in this module (Q25, owner
# clarification 2026-10-19: the annotation side outputs table/df products only;
# the catalogue layer is the plotting consumer):
#
#   read_gff3_genes()  — gene rows of a (gzipped) GFF3 as a data.table
#   read_gene_xrefs()  — an Ensembl Plants xref dump collapsed to one row per
#                        gene, replicating block_annotation 0.3.0's
#                        load_gene_xrefs() rules exactly (Q25.g)
#   .read_qtl_bed()    — internal: 4-column QTL BED (chrom start end name)
#
# Every reader is a descendant of the .read_table_chr() lineage (thread #032):
# text columns stay text, numeric columns are coerced explicitly, and nothing
# about a file's meaning is left to type inference.
# =============================================================================


# Internal: pull one key's value out of a GFF3 attribute column (column 9).
# Vectorised over rows; returns NA_character_ where the key is absent.
# A value runs to the next ";" — GFF3 escapes literal ";" as %3B, so this is
# exact, and percent-decoding is deliberately NOT applied (the Python reference
# does not decode either; IDs must round-trip byte-identically).
.gff3_attr <- function(attrs, key) {
  pat <- paste0("(?:^|;)\\s*", key, "\\s*=([^;]*)")
  m   <- regexpr(pat, attrs, perl = TRUE)
  out <- rep(NA_character_, length(attrs))
  hit <- m > 0L
  if (any(hit)) {
    cap <- attr(m, "capture.start")[hit, 1L]
    len <- attr(m, "capture.length")[hit, 1L]
    out[hit] <- trimws(substr(attrs[hit], cap, cap + len - 1L))
  }
  out[!is.na(out) & !nzchar(out)] <- NA_character_
  out
}


#' Read gene rows from a GFF3 annotation file
#'
#' A minimal GFF3 parser for the offline annotation route (thread #033): it
#' extracts rows of one feature type (default \code{"gene"}) with their
#' coordinates, strand and the \code{ID=} / \code{Name=} / \code{biotype=}
#' attributes. It replicates \code{block_annotation}'s
#' \code{_parse_gff3_genes()}: comment lines are skipped, a \code{##FASTA}
#' tail ends the scan, no attribute but the coordinates is required, and
#' feature IDs are kept \strong{exactly as the file spells them} (the
#' \code{gene:} prefix included — it is stripped only at the xref join inside
#' \code{\link{compute_block_annotation}}).
#'
#' NCBI-style files write \code{gene_biotype=} and \code{gene_id=} /
#' \code{gene_name=}; both spellings are accepted (Q25.f), with the Ensembl
#' spelling winning when both are present.
#'
#' Gzipped files (\code{.gff3.gz}) are read through a \code{gzfile()}
#' connection; no new dependency.
#'
#' @param path         Path to the GFF3 file (plain or \code{.gz}).
#' @param feature_type GFF3 type (column 3) to keep, case-insensitive.
#'   Default \code{"gene"} — the only type the local annotation route uses.
#' @param rename_seqid Optional named character vector renaming the GFF3
#'   seqids as they are stored (names = seqid as the file spells it, values =
#'   the name to store). Usually left \code{NULL}:
#'   \code{\link{compute_block_annotation}} maps the \emph{container's}
#'   chromosome names onto the GFF3 vocabulary with its \code{chr_map}
#'   instead, which is the direction every haploStack \code{chr_map} uses
#'   (container spelling -> the other table's spelling).
#' @param chr_map      Deprecated since 0.9.0 (warns): the old name of
#'   \code{rename_seqid}, which renames in the opposite direction to every
#'   other \code{chr_map} in the package.
#'
#' @return A \code{data.table} with columns \code{chrom}, \code{start},
#'   \code{end} (numeric, 1-based inclusive), \code{strand}, \code{gene_id},
#'   \code{gene_name}, \code{biotype} (character; \code{NA} where the
#'   attribute is absent) and \code{feature_type} (the file's spelling of the
#'   type kept), in file order. The source is recorded in the
#'   \code{"hs_annotation_source"} attribute for provenance as \code{file}
#'   (the base name), \code{md5}, \code{size} and \code{n}; no absolute path
#'   and no modification time, so two reads of one file give one record.
#'
#' @seealso \code{\link{compute_block_annotation}}, \code{\link{read_gene_xrefs}}
#' @export
read_gff3_genes <- function(path, feature_type = "gene", rename_seqid = NULL,
                            chr_map = NULL) {
  path <- as.character(path)
  if (!file.exists(path))
    .hs_abort("read_gff3_genes: file not found: ", path)
  if (!is.null(chr_map)) {                      # deprecated since 0.9.0 (D5 X1)
    .hs_deprecated("chr_map", "rename_seqid", "0.9.0", what = "argument",
                   fn = "read_gff3_genes")
    if (is.null(rename_seqid)) rename_seqid <- chr_map
  }
  if (!is.character(feature_type) || length(feature_type) != 1L)
    .hs_abort("read_gff3_genes: 'feature_type' must be a single character value.")

  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)

  # A ##FASTA directive ends the annotation section; everything after it is
  # sequence and must not be parsed as rows.
  fasta_at <- match(TRUE, startsWith(lines, "##FASTA"))
  if (!is.na(fasta_at)) lines <- lines[seq_len(fasta_at - 1L)]

  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) == 0L)
    .hs_abort("read_gff3_genes: no annotation rows in ", path)

  # Pre-filter (ann/F12): only lines that can hold the feature type in a
  # tab-delimited column are split. A superset (column 2 could match too), so
  # the exact type test below stays and the result is identical; the peak on
  # MorexV3 drops from about 400 MB to under 100 MB.
  ft_rx <- paste0("\t", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", feature_type), "\t")
  cand   <- lines[grepl(ft_rx, lines, ignore.case = TRUE)]
  fields <- strsplit(cand, "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 9L]
  type   <- vapply(fields, `[[`, character(1L), 3L)
  fields <- fields[tolower(type) == tolower(feature_type)]
  if (length(fields) == 0L) {
    all_f <- strsplit(lines, "\t", fixed = TRUE)
    all_f <- all_f[lengths(all_f) >= 9L]
    if (length(all_f) == 0L)
      .hs_abort("read_gff3_genes: no 9-column GFF3 rows in ", path)
    types <- .hs_sort(unique(tolower(vapply(all_f, `[[`, character(1L), 3L))))
    .hs_abort("read_gff3_genes: no '", feature_type, "' rows in ", path,
         ". Feature types present: ",
         paste(utils::head(types, 10L), collapse = ", "),
         if (length(types) > 10L) " ..." else "", ".")
  }
  type <- vapply(fields, `[[`, character(1L), 3L)

  attrs <- vapply(fields, `[[`, character(1L), 9L)
  out <- data.table::data.table(
    chrom     = vapply(fields, `[[`, character(1L), 1L),
    start     = as.numeric(vapply(fields, `[[`, character(1L), 4L)),
    end       = as.numeric(vapply(fields, `[[`, character(1L), 5L)),
    strand    = vapply(fields, `[[`, character(1L), 7L),
    gene_id   = data.table::fcoalesce(.gff3_attr(attrs, "ID"),
                                      .gff3_attr(attrs, "gene_id")),
    gene_name = data.table::fcoalesce(.gff3_attr(attrs, "Name"),
                                      .gff3_attr(attrs, "gene_name")),
    biotype   = data.table::fcoalesce(.gff3_attr(attrs, "biotype"),
                                      .gff3_attr(attrs, "gene_biotype")),
    feature_type = type
  )

  if (anyNA(out$start) || anyNA(out$end))
    .hs_abort("read_gff3_genes: non-numeric start/end coordinate in ", path)

  if (!is.null(rename_seqid)) {
    if (is.null(names(rename_seqid)) || any(!nzchar(names(rename_seqid))))
      .hs_abort("read_gff3_genes: 'rename_seqid' must be a fully named character vector.")
    hit <- out$chrom %in% names(rename_seqid)
    out$chrom[hit] <- unname(rename_seqid[out$chrom[hit]])
  }

  data.table::setattr(out, "hs_annotation_source", .hs_source_record(path, nrow(out)))
  out[]
}

# The provenance record of an annotation input (ann/F6, F20): base name, md5,
# size and row count. No absolute path (a container travels between machines)
# and no modification time (two reads of one file are one record).
.hs_source_record <- function(path, n) {
  list(file = basename(path),
       md5  = unname(tools::md5sum(path)),
       size = as.numeric(file.size(path)),
       n    = as.integer(n))
}


#' Read an Ensembl xref dump into a per-gene annotation table
#'
#' Collapses an Ensembl Plants FTP \code{*.xref.tsv(.gz)} dump to one row per
#' gene, replicating \code{block_annotation} 0.3.0's \code{load_gene_xrefs()}
#' exactly (Q25.g):
#' \itemize{
#'   \item \code{go_terms} — rows with \code{db_name == "GO"} as
#'     \code{ID:description} pairs, de-duplicated in order, joined by
#'     \code{"; "};
#'   \item \code{interpro_domains} — rows with \code{db_name} INTERPRO or
#'     INTERPROSCAN, likewise;
#'   \item \code{gene_symbol} — \code{"LOC"} + the \code{xref_id} of the first
#'     NCBI GENE row;
#'   \item \code{gene_description} — the first description in priority order
#'     Swiss-Prot, TrEMBL, NCBI GENE, RFAM, tRNAscan-SE.
#' }
#'
#' A dump with no InterPro rows triggers a warning (Ensembl ships protein
#' domains in a separate file), matching the Python route.
#'
#' @param path Path to the xref TSV (plain or \code{.gz}). Required columns:
#'   \code{gene_stable_id}, \code{xref_id}, \code{xref_label},
#'   \code{description}, \code{db_name}.
#'
#' @return A \code{data.table} with columns \code{gene_id} (the \strong{bare}
#'   stable ID, no \code{gene:} prefix), \code{gene_symbol},
#'   \code{gene_description}, \code{go_terms}, \code{interpro_domains}
#'   (character, \code{""} where absent). Source provenance in the
#'   \code{"hs_annotation_source"} attribute.
#'
#' @seealso \code{\link{compute_block_annotation}}, \code{\link{read_gff3_genes}}
#' @export
read_gene_xrefs <- function(path) {
  path <- as.character(path)
  if (!file.exists(path))
    .hs_abort("read_gene_xrefs: file not found: ", path)

  need <- c("gene_stable_id", "xref_id", "xref_label", "description", "db_name")
  xr <- data.table::fread(path, sep = "\t", select = need,
                          colClasses = "character", data.table = TRUE,
                          showProgress = FALSE)
  miss <- setdiff(need, names(xr))
  if (length(miss) > 0L)
    .hs_abort("read_gene_xrefs: required column(s) missing from ", path, ": ",
         paste(miss, collapse = ", "))
  for (j in need) data.table::set(xr, j = j,
                                  value = data.table::fifelse(is.na(xr[[j]]),
                                                              "", xr[[j]]))
  xr <- xr[nzchar(gene_stable_id)]
  if (nrow(xr) == 0L)
    .hs_abort("read_gene_xrefs: no xref rows with a gene_stable_id in ", path)

  dbu <- toupper(xr$db_name)

  # ID:description pairs joined per gene, de-duplicated in file order —
  # the exact collapse of the Python _joined().
  joined <- function(mask, id_col) {
    sub <- xr[mask & nzchar(xr[[id_col]]) & nzchar(xr$description)]
    if (nrow(sub) == 0L)
      return(data.table::data.table(gene_stable_id = character(0L),
                                    val = character(0L)))
    sub[, .(val = paste(unique(paste0(get(id_col), ":", description)),
                        collapse = "; ")),
        by = gene_stable_id]
  }

  go  <- joined(dbu == "GO", "xref_id")
  ipr <- joined(dbu %in% c("INTERPRO", "INTERPROSCAN"), "xref_id")

  ncbi <- xr[startsWith(dbu, "NCBI GENE") & nzchar(xref_id),
             .(val = paste0("LOC", xref_id[1L])), by = gene_stable_id]

  # Description priority: an earlier database wins; later ones only fill genes
  # still undescribed (combine_first semantics).
  desc <- data.table::data.table(gene_stable_id = character(0L),
                                 val = character(0L))
  for (nm in c("UNIPROTKB/SWISS-PROT", "UNIPROTKB/TREMBL",
               "NCBI GENE (FORMERLY ENTREZGENE)", "RFAM", "TRNASCAN-SE")) {
    sub <- xr[dbu == nm & nzchar(description),
              .(val = description[1L]), by = gene_stable_id]
    if (nrow(sub) == 0L) next
    desc <- rbind(desc, sub[!gene_stable_id %in% desc$gene_stable_id])
  }

  genes <- data.table::data.table(gene_id = unique(xr$gene_stable_id))
  grab <- function(tbl) {
    v <- tbl$val[match(genes$gene_id, tbl$gene_stable_id)]
    data.table::fifelse(is.na(v), "", v)
  }
  out <- genes[, `:=`(gene_symbol      = grab(ncbi),
                      gene_description = grab(desc),
                      go_terms         = grab(go),
                      interpro_domains = grab(ipr))]

  if (!any(nzchar(out$interpro_domains)))
    .hs_warn("read_gene_xrefs: no InterPro entries in this xref dump \u2014 ",
            "Ensembl ships protein domains in a separate file, so ",
            "'interpro_domains' stays empty.")

  data.table::setattr(out, "hs_annotation_source", .hs_source_record(path, nrow(out)))
  out[]
}


# Internal: read a 4-column QTL BED (chrom, start, end, name). Comment lines
# ("#..."), "track" and "browser" lines are skipped before parsing (ann/F27).
# Coordinates are used exactly as stored — the Python route applies no
# 0-based/half-open correction, and parity (Q25.c) pins that choice; a caller
# holding true BED semantics shifts start before writing. A non-numeric or
# missing coordinate, or start > end, is refused with hs_error_format.
.read_qtl_bed <- function(path) {
  path <- as.character(path)
  if (!file.exists(path))
    .hs_abort(".read_qtl_bed: file not found: ", path)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#") &
                   !grepl("^(track|browser)([[:space:]]|$)", lines)]
  if (length(lines) == 0L)
    .hs_abort(class = c("hs_error_format", "hs_error_input"),
              ".read_qtl_bed: no QTL rows in ", path, ".")
  fields <- strsplit(lines, "\t", fixed = TRUE)
  if (any(lengths(fields) < 4L))
    .hs_abort(class = c("hs_error_format", "hs_error_input"),
              ".read_qtl_bed: a QTL BED needs 4 tab-separated columns ",
              "(chrom start end name); line ",
              which(lengths(fields) < 4L)[1L], " of ", path, " has ",
              min(lengths(fields)), ".")
  out <- data.table::data.table(chrom = vapply(fields, `[[`, "", 1L),
                                start = suppressWarnings(as.numeric(vapply(fields, `[[`, "", 2L))),
                                end   = suppressWarnings(as.numeric(vapply(fields, `[[`, "", 3L))),
                                name  = vapply(fields, `[[`, "", 4L))
  bad <- which(is.na(out$start) | is.na(out$end))
  if (length(bad))
    .hs_abort(class = c("hs_error_format", "hs_error_input"),
              ".read_qtl_bed: non-numeric start/end on line ", bad[1L], " of ", path, ".")
  rev <- which(out$start > out$end)
  if (length(rev))
    .hs_abort(class = c("hs_error_format", "hs_error_input"),
              ".read_qtl_bed: start > end on line ", rev[1L], " of ", path,
              " (", out$name[rev[1L]], ").")
  data.table::setattr(out, "hs_annotation_source", .hs_source_record(path, nrow(out)))
  out
}
