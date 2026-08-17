# =============================================================================
# catalogue_family.R
#
# Family-structure tools for the haplotype catalogue:
#   trace_family_history()      — walk a pedigree generation by generation
#   attach_family_history()     — store the history as GenoMeta metadata
#   family_history_to_pedigree() — adapter to the trace_haplotype_family() input
#   trace_haplotype_family()    — block-level parent/offspring similarity
# =============================================================================


# ---------------------------------------------------------------------------
# Internal: fam-table normalisation
# ---------------------------------------------------------------------------

# Resolve the individual-ID vector of a fam table, following the package
# convention used by build_geno_meta(): explicit id_col > "IID" column >
# meaningful rownames > first column.
.fam_resolve_ids <- function(fam, id_col = NULL) {
  fam <- as.data.frame(fam, stringsAsFactors = FALSE)
  if (!is.null(id_col)) {
    if (is.character(id_col)) {
      if (!id_col %in% names(fam))
        stop(".fam_resolve_ids: id_col '", id_col, "' not found in fam.")
      return(as.character(fam[[id_col]]))
    }
    if (id_col < 1L || id_col > ncol(fam))
      stop(".fam_resolve_ids: id_col index ", id_col, " out of range.")
    return(as.character(fam[[id_col]]))
  }
  ids <- if ("IID" %in% names(fam)) fam$IID else rownames(fam)
  if (is.null(ids) || all(ids == as.character(seq_len(nrow(fam)))))
    ids <- fam[[1L]]   # fall back to first column
  as.character(ids)
}

# Convert a (wide) fam table to long parent edges: data.table(IID, Parent).
# Each column in parent_cols contributes one (IID, Parent) row per individual.
# NA, "", and "0" (PLINK missing) are treated as founder (Parent = NA).
.fam_to_edges <- function(fam, parent_cols, id_col = NULL) {
  fam <- as.data.frame(fam, stringsAsFactors = FALSE)
  ids <- .fam_resolve_ids(fam, id_col)

  # Resolve parent columns (accept names or indices)
  if (is.character(parent_cols)) {
    miss <- setdiff(parent_cols, names(fam))
    if (length(miss) > 0L)
      stop(".fam_to_edges: parent column(s) not found: ",
           paste(miss, collapse = ", "))
    pcols <- parent_cols
  } else {
    if (any(parent_cols < 1L | parent_cols > ncol(fam)))
      stop(".fam_to_edges: parent_cols index out of range.")
    pcols <- names(fam)[parent_cols]
  }

  edges <- data.table::rbindlist(lapply(pcols, function(pc) {
    data.table::data.table(IID = ids, Parent = as.character(fam[[pc]]))
  }))
  edges[Parent %in% c("", "0", "NA"), Parent := NA_character_]
  unique(edges)
}

# Group one generation: starting from a set of founders, find their children
# and bucket children sharing the same (sorted, unique) parent set into a
# family (FID). Returns data.table(FID, IID, Label) with Label "P"/"O".
.family_grouping <- function(edges, founder_identifiers = NA) {
  empty <- data.table::data.table(FID = character(0L), IID = character(0L),
                                  Label = character(0L))

  # Generation 1 founders are individuals with a missing-parent row; later
  # generations take the previous generation's offspring directly as the new
  # founders (the reference looked one level too deep, dropping the final
  # generation's families).
  if (length(founder_identifiers) == 1L && is.na(founder_identifiers[1L])) {
    founders_iid <- unique(edges$IID[is.na(edges$Parent)])
  } else {
    founders_iid <- unique(as.character(founder_identifiers))
  }
  if (length(founders_iid) == 0L) return(empty)

  # Children of the founders (drop NA-parent founder rows themselves)
  f1 <- edges[!is.na(Parent) & Parent %in% founders_iid]
  if (nrow(f1) == 0L) return(empty)

  # Sorted unique parent set per child
  combs <- lapply(split(f1$Parent, f1$IID), function(p) sort(unique(p)))

  unique_combs <- unique(unname(combs))
  parts <- lapply(seq_along(unique_combs), function(i) {
    key      <- unique_combs[[i]]
    members  <- names(combs)[vapply(combs, function(x) identical(x, key),
                                    logical(1L))]
    fid      <- paste0("F_", i)
    data.table::rbindlist(list(
      data.table::data.table(
        FID = fid, IID = members,
        Label = ifelse(members %in% founders_iid, "P", "O")),
      data.table::data.table(FID = fid, IID = key, Label = "P")
    ))
  })

  unique(data.table::rbindlist(parts))
}


# ---------------------------------------------------------------------------
# Exported: family-history tracing
# ---------------------------------------------------------------------------

#' Trace a pedigree generation by generation
#'
#' Walks a family/pedigree (\dQuote{fam}) table from the oldest generation to
#' the newest, grouping individuals into families by their unique parent
#' combination. Adapted from the project's standalone family-tracing routine.
#'
#' The fam table is interpreted in wide form: one or more columns hold the
#' parent IDs (given by \code{parent_cols}). Each parent column contributes one
#' parent slot; a missing parent (\code{NA}, \code{""}, or \code{"0"}) marks a
#' founder. Individual IDs follow the package convention used by
#' \code{\link{build_geno_meta}}: an explicit \code{id_col}, else an \code{IID}
#' column, else meaningful row names, else the first column.
#'
#' @param fam         A data.frame / data.table pedigree. Holds an individual-ID
#'   column and one or more parent columns.
#' @param parent_cols Column name(s) or index(es) of the parent columns.
#' @param id_col      Optional column name or index for the individual ID.
#'   \code{NULL} (default) uses the fam convention above.
#' @param founder_id  Value identifying the earliest generation's parents.
#'   Default \code{NA} (founders are individuals with no recorded parent).
#'
#' @return A list of per-generation \code{data.table}s, each with columns
#'   \code{FID} (family ID), \code{IID} (individual ID) and \code{Label}
#'   (\code{"P"} parent / \code{"O"} offspring). Element \code{[[1]]} is the
#'   oldest generation. Empty list when no families are observable.
#'
#' @examples
#' \dontrun{
#' fam <- data.frame(IID     = c("F1","F2","O1","O2","G1"),
#'                   Father  = c(NA,  NA,  "F1", "F1", "O1"),
#'                   Mother  = c(NA,  NA,  "F2", "F2", "O2"))
#' hist <- trace_family_history(fam, parent_cols = c("Father","Mother"))
#' length(hist)   # number of generations traced
#' }
#'
#' @seealso \code{\link{attach_family_history}},
#'   \code{\link{family_history_to_pedigree}}, \code{\link{trace_haplotype_family}}
#'
#' @export
trace_family_history <- function(fam, parent_cols, id_col = NULL,
                                 founder_id = NA) {
  edges <- .fam_to_edges(fam, parent_cols, id_col)

  history <- list()
  hid     <- 1L
  fi      <- founder_id
  repeat {
    fm <- .family_grouping(edges, fi)
    if (nrow(fm) == 0L) break
    history[[hid]] <- fm
    fi  <- unique(fm$IID[fm$Label == "O"])
    hid <- hid + 1L
    if (length(fi) == 0L) break
  }

  if (length(history) > 0L)
    names(history) <- paste0("gen_", seq_along(history))
  history
}


#' Attach a family-tracing history to a GenoMeta and update sample_info
#'
#' Computes \code{\link{trace_family_history}}, stores the result at
#' \code{geno_meta$metadata$family_history}, and writes each sample's
#' (generation-qualified) family ID back into the \code{Family} column of
#' \code{geno_meta$sample_info}, keeping the pedigree table and the tracing
#' history in sync.
#'
#' By default the pedigree is read from the object's own FAM-structured
#' \code{sample_info} (individual ID at column 2, parents at columns 3-4). Pass an
#' explicit \code{fam} + \code{parent_cols} to trace a separate pedigree instead.
#' \code{metadata$family_history} is supplementary: \code{\link{is_geno_meta}}
#' does not validate it, and it (like the updated \code{sample_info}) rides
#' through \code{\link{filter_geno_meta}} / \code{\link{recode_geno_meta}}.
#'
#' @param geno_meta  A valid \code{GenoMeta} list.
#' @param fam        Optional pedigree table (see
#'   \code{\link{trace_family_history}}). \code{NULL} (default) traces
#'   \code{geno_meta$sample_info}.
#' @param parent_cols Parent column name(s) or index(es). When tracing
#'   \code{sample_info} (\code{fam = NULL}) defaults to \code{c(3, 4)}
#'   (\code{Parent1}/\code{Parent2}); required when an external \code{fam} is given.
#' @param id_col     Individual-ID column (name or index). Defaults to \code{2}
#'   (the \code{Sample} column) when tracing \code{sample_info}; \code{NULL}
#'   otherwise (resolved by the fam convention in \code{\link{trace_family_history}}).
#' @param founder_id Founder identifier. Default \code{NA}.
#'
#' @return The \code{GenoMeta} with \code{metadata$family_history} populated and
#'   \code{sample_info$Family} updated.
#'
#' @seealso \code{\link{trace_family_history}},
#'   \code{\link{family_history_to_pedigree}}
#'
#' @export
attach_family_history <- function(geno_meta, fam = NULL,
                                  parent_cols = NULL, id_col = NULL,
                                  founder_id = NA) {
  if (!is_geno_meta(geno_meta))
    stop("attach_family_history: 'geno_meta' is not a valid GenoMeta.")

  if (is.null(fam)) {
    fam <- geno_meta$sample_info                 # trace the object's own FAM
    if (is.null(parent_cols)) parent_cols <- c(3L, 4L)
    if (is.null(id_col))      id_col      <- 2L
  } else if (is.null(parent_cols)) {
    stop("attach_family_history: 'parent_cols' is required when 'fam' is given.")
  }
  history <- trace_family_history(fam, parent_cols, id_col, founder_id)

  hist_ids <- unique(unlist(lapply(history, function(g) g$IID)))
  geno_samples <- dimnames(geno_meta$geno_arrays)[[1L]]
  if (length(intersect(hist_ids, geno_samples)) == 0L)
    warning("attach_family_history: no family-history individual matches the ",
            "genotype sample set — history stored but not linkable to genotypes.")

  geno_meta$metadata$family_history <- history

  # Write-back: each sample's offspring-family FID (generation-qualified, e.g.
  # "gen_1.F_3") into sample_info$Family, leaving founders / unmatched as-is.
  fid_map <- .sample_family_map(history)
  si <- geno_meta$sample_info
  hit <- si$Sample %in% names(fid_map)
  si$Family[hit] <- unname(fid_map[si$Sample[hit]])
  geno_meta$sample_info <- si

  geno_meta
}


# Internal: map each offspring sample to its generation-qualified family ID,
# from a trace_family_history() result. A sample keeps the FID of the family in
# which it is an offspring (its birth family).
.sample_family_map <- function(history) {
  out <- character(0L)
  for (g in seq_along(history)) {
    gen   <- as.data.frame(history[[g]], stringsAsFactors = FALSE)
    label <- if (!is.null(names(history))) names(history)[g] else paste0("gen_", g)
    offs  <- gen[gen$Label == "O", , drop = FALSE]
    keep  <- !(offs$IID %in% names(out))          # first (oldest) birth family wins
    add   <- offs[keep, , drop = FALSE]
    if (nrow(add) > 0L)
      out[add$IID] <- paste0(label, ".", add$FID)
  }
  out
}


#' Convert a family history into a Parent1/Parent2 pedigree
#'
#' Flattens the per-generation family history into the
#' \code{Sample}/\code{Parent1}/\code{Parent2} pedigree consumed by
#' \code{\link{trace_haplotype_family}}. One row per offspring; single-parent
#' families pad \code{Parent2 = NA}. Families with more than two parents are
#' truncated to the first two (with a warning).
#'
#' @param x A \code{GenoMeta} carrying \code{metadata$family_history}, or a
#'   family-history list as returned by \code{\link{trace_family_history}}.
#'
#' @return A \code{data.frame} with columns \code{Sample}, \code{Parent1},
#'   \code{Parent2}.
#'
#' @seealso \code{\link{trace_family_history}}, \code{\link{trace_haplotype_family}}
#'
#' @export
family_history_to_pedigree <- function(x) {
  history <- if (is.list(x) && !is.null(x$metadata))
    x$metadata$family_history
  else x
  if (is.null(history) || length(history) == 0L)
    return(data.frame(Sample = character(0L), Parent1 = character(0L),
                      Parent2 = character(0L), stringsAsFactors = FALSE))

  # FIDs ("F_1", ...) are only unique within a generation, so iterate per
  # generation and per family rather than pooling all rows together.
  warned <- FALSE
  rows <- list()
  for (gen in history) {
    gen <- data.table::as.data.table(gen)
    for (fid in unique(gen$FID[gen$Label == "O"])) {
      grp     <- gen[FID == fid]
      parents <- unique(grp$IID[grp$Label == "P"])
      offs    <- unique(grp$IID[grp$Label == "O"])
      if (length(parents) > 2L) {
        if (!warned) {
          warning("family_history_to_pedigree: family '", fid, "' has >2 ",
                  "parents; truncating to the first two.")
          warned <- TRUE
        }
        parents <- parents[1:2]
      }
      rows[[length(rows) + 1L]] <- data.frame(
        Sample  = offs,
        Parent1 = if (length(parents) >= 1L) parents[1L] else NA_character_,
        Parent2 = if (length(parents) >= 2L) parents[2L] else NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  unique(do.call(rbind, rows))
}


#' Trace haplotype inheritance across a pedigree
#'
#' Given a GenoMeta object and a pedigree table, identifies which haplotype
#' alleles in offspring samples are inherited from founder parents. Returns an
#' inheritance assignment table linking each sample-haplotype-marker combination
#' to its inferred parental origin.
#'
#' @param geno_meta A valid \code{GenoMeta} list from
#'   \code{\link{new_geno_meta}} or any \code{read_*} function.
#' @param pedigree data.frame with columns \code{Sample}, \code{Parent1},
#'   \code{Parent2}. Parents must be present in the genotype sample set or
#'   have \code{NA} for founders.
#' @param blocks   Optional character vector of block IDs to restrict the
#'   analysis. Default \code{NULL} (all blocks).
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{\code{Sample}}{Offspring sample ID.}
#'   \item{\code{Haplotype_idx}}{Haplotype index (1 or 2).}
#'   \item{\code{Block}}{Block ID.}
#'   \item{\code{Parent}}{Inferred parental origin (\code{Parent1},
#'     \code{Parent2}, \code{"ambiguous"}, or \code{"unknown"}).}
#'   \item{\code{n_match}}{Number of markers matching the assigned parent.}
#'   \item{\code{n_total}}{Total markers compared in the block.}
#' }
#'
#' @details
#' For each offspring-block combination, the function counts how many markers
#' in each haplotype match each parent's alleles and assigns the parent with
#' the higher count. Ties are labelled \code{"ambiguous"}; offspring with no
#' parent in the profile are labelled \code{"unknown"}.
#'
#' @examples
#' \dontrun{
#' ped <- data.frame(Sample  = c("Offspring1","Offspring2"),
#'                   Parent1 = c("ParentA","ParentA"),
#'                   Parent2 = c("ParentB","ParentC"))
#' fam <- trace_haplotype_family(meta, ped)
#' head(fam)
#' }
#'
#' @export
trace_haplotype_family <- function(geno_meta, pedigree, blocks = NULL) {
  if (!is_geno_meta(geno_meta))
    stop("trace_haplotype_family: 'geno_meta' is not a valid GenoMeta.")

  ped <- as.data.table(pedigree)
  if (!all(c("Sample","Parent1","Parent2") %in% names(ped)))
    stop("trace_haplotype_family: pedigree must have columns Sample, Parent1, Parent2.")

  # Filter to relevant blocks, then build profile views on demand
  gm   <- if (!is.null(blocks)) filter_geno_meta(geno_meta, blocks = blocks)
          else geno_meta
  prof <- geno_meta_views(gm)

  if (nrow(prof$alleles) == 0L)
    return(data.table(Sample = character(0L), Haplotype_idx = integer(0L),
                      Block  = character(0L), Parent = character(0L),
                      n_match = integer(0L),  n_total = integer(0L)))

  alleles     <- prof$alleles
  marker_meta <- prof$marker_meta
  marker_names <- marker_meta[!is.na(block_id), marker]

  # Melt to long: Sample, Haplotype_idx, Marker, Allele
  long_dt <- data.table::melt(
    alleles[, c("Sample","Haplotype_idx", marker_names), with = FALSE],
    id.vars       = c("Sample","Haplotype_idx"),
    measure.vars  = marker_names,
    variable.name = "Marker",
    value.name    = "Allele",
    variable.factor = FALSE
  )

  # Attach Block
  blk_map <- marker_meta[!is.na(block_id), .(Marker = marker, Block = block_id)]
  long_dt[blk_map, on = "Marker", Block := i.Block]
  long_dt <- long_dt[!is.na(Block)]

  # For each offspring in pedigree
  results <- lapply(seq_len(nrow(ped)), function(idx) {
    row      <- ped[idx]
    offspring <- row$Sample
    p1        <- row$Parent1
    p2        <- row$Parent2

    off_dt <- long_dt[Sample == offspring]
    if (nrow(off_dt) == 0L) return(NULL)

    # Get parent alleles per marker (use haplotype 1 as reference)
    p1_alleles <- if (!is.na(p1) && p1 %in% prof$samples)
      long_dt[Sample == p1 & Haplotype_idx == 1L, .(Marker, P1_Allele = Allele)]
    else data.table(Marker = character(0L), P1_Allele = character(0L))

    p2_alleles <- if (!is.na(p2) && p2 %in% prof$samples)
      long_dt[Sample == p2 & Haplotype_idx == 1L, .(Marker, P2_Allele = Allele)]
    else data.table(Marker = character(0L), P2_Allele = character(0L))

    off_dt[p1_alleles, on = "Marker", P1_Allele := i.P1_Allele]
    off_dt[p2_alleles, on = "Marker", P2_Allele := i.P2_Allele]

    # Count matches per (Haplotype_idx, Block)
    block_result <- off_dt[, {
      n_p1 <- sum(Allele == P1_Allele, na.rm = TRUE)
      n_p2 <- sum(Allele == P2_Allele, na.rm = TRUE)
      n_tot <- .N
      parent <- if (is.na(p1) && is.na(p2)) "unknown"
                else if (n_p1 > n_p2)  "Parent1"
                else if (n_p2 > n_p1)  "Parent2"
                else if (n_p1 == n_p2 && n_p1 > 0L) "ambiguous"
                else "unknown"
      .(Parent  = parent,
        n_match = max(n_p1, n_p2),
        n_total = n_tot)
    }, by = .(Sample, Haplotype_idx, Block)]

    block_result
  })

  rbindlist(Filter(Negate(is.null), results))
}
