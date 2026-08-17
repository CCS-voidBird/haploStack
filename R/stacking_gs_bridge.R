# =============================================================================
# stacking_gs_bridge.R
#
# Bridge between the GenoMeta world (crossing-plan making) and the
# genomicSimulation package (linkage-aware crossing simulation).
#
# genomicSimulation only ingests data from FILES and identifies genotypes by an
# internal index ("X") that is reshuffled by delete.group(); these functions
# encapsulate the file round-trip and the name<->index mapping so the rest of
# the stacking workflow can stay GenoMeta-native.
#
# Functions:
#   gs_load_geno_meta()      — write temp files, load.data() -> a "gs handle"
#   gs_cross_plan()          — realise a crossing plan (make.targeted.crosses)
#   update_geno_meta_from_matrix() — gene-data matrix + template -> GenoMeta (gS-free)
#   gs_group_to_geno_meta()  — read the live group's genotypes back to a GenoMeta
#   gs_group_gebv()          — named breeding-value vector for the live group
#   gs_sync_group_names()    — push compact GenoMeta IDs into the live group (FR3)
#   gs_close()               — free a persistent session's live group (FR1)
#
# Internal:
#   .gs_require()            — guard: genomicSimulation must be installed
#   .gs_write_files()        — GenoMeta -> (geno, map, eff) files for load.data()
#   .gs_genes_to_array()     — gene-data matrix -> 3-D allele array
#
# NB: genomicSimulation is in Suggests. Every function that touches the package
# is guarded; the file/parse helpers work without it. The exact load.data() file
# layout may need tuning against the installed genomicSimulation version
# (tracked in dev/PRD_stacking_simulation.md, issue 04 crossing-backend).
# =============================================================================


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.gs_require <- function() {
  if (!requireNamespace("genomicSimulation", quietly = TRUE))
    stop("This function needs the 'genomicSimulation' package (Suggests). ",
         "Install it, or use cross_backend = 'block'.", call. = FALSE)
}

# Convert map positions to the centimorgans genomicSimulation expects. gS reads
# the map's Pos column as cM (recombination distance): the expected number of
# crossovers per gamete on a chromosome is proportional to its length in Morgans.
# A physical map (bp) passed through unscaled therefore becomes astronomically
# long (e.g. 99 Mb -> 99,000 Morgans), which makes meiosis simulation
# catastrophically slow AND genetically meaningless (every interval recombines).
# Positions that look physical (max > 1e4, implausible for a cM map) are scaled
# bp -> cM at `cM_per_mb` (default 1 cM/Mb, the common human/plant rule of thumb);
# positions already in a cM range are passed through unchanged.
.gs_positions_to_cM <- function(pos, cM_per_mb = 1) {
  pos <- as.numeric(pos)
  if (length(pos) == 0L || all(is.na(pos))) return(pos)
  if (max(pos, na.rm = TRUE) > 1e4) pos <- pos / 1e6 * cM_per_mb   # bp -> cM
  pos
}

# Write a GenoMeta (+ effects) as the genotype / map / effect file trio that
# genomicSimulation::load.data() reads. Returns the three paths.
#   geno : markers x samples, two-character alleles, row+col names (GS preset).
#   map  : Marker  Chrom  Pos(cM)   (header).
#   eff  : Marker  Allele Effect (header).
.gs_write_files <- function(geno_meta, effect_df = NULL,
                            work_dir = tempdir(), prefix = "gs",
                            cM_per_mb = 1) {
  if (!is_geno_meta(geno_meta))
    stop(".gs_write_files: 'geno_meta' is not a valid GenoMeta.")
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

  geno_file <- file.path(work_dir, paste0(prefix, "_geno.txt"))
  map_file  <- file.path(work_dir, paste0(prefix, "_map.txt"))
  eff_file  <- file.path(work_dir, paste0(prefix, "_eff.txt"))

  # Genotype matrix in genomicSimulation orientation (markers x samples, "AT").
  mat <- geno_meta_to_matrix(geno_meta, preset = "GS")$matrix
  utils::write.table(mat, geno_file, sep = "\t", quote = FALSE,
                     col.names = NA)          # leading blank header cell

  # Map: Marker Chrom Pos(cM). Positions written non-scientific so gS parses them.
  map_out <- data.frame(
    Marker = as.character(geno_meta$map_info$Marker),
    Chrom  = as.character(geno_meta$map_info$Chrom),
    Pos    = format(.gs_positions_to_cM(geno_meta$map_info$Pos, cM_per_mb),
                    scientific = FALSE, trim = TRUE),
    stringsAsFactors = FALSE)
  utils::write.table(map_out, map_file, sep = "\t", quote = FALSE,
                     row.names = FALSE)

  # Effects: Marker Allele Effect (NULL -> empty file, BV will be 0). The stored
  # marker_effects may be a named per-trait list; collapse to a single table.
  eff <- .coerce_single_effect_df(effect_df %||% geno_meta$marker_effects)
  if (!is.null(eff)) {
    eff <- as.data.frame(eff[, 1:3, drop = FALSE], stringsAsFactors = FALSE)
    names(eff) <- c("Marker", "Allele", "Effect")
    utils::write.table(eff, eff_file, sep = "\t", quote = FALSE,
                       row.names = FALSE)
  } else {
    eff_file <- NULL
  }

  c(geno = geno_file, map = map_file, eff = if (is.null(eff_file)) NA else eff_file)
}

# Convert a genomicSimulation gene-data matrix (markers x individuals, each cell
# a 2-character genotype, e.g. "AT") into a 3-D allele array
# [n_samples x n_markers x 2] matching GenoMeta$geno_arrays.
.gs_genes_to_array <- function(gene_mat) {
  markers <- rownames(gene_mat)
  samples <- colnames(gene_mat)
  h1 <- matrix(substr(gene_mat, 1L, 1L), nrow = nrow(gene_mat),
               dimnames = dimnames(gene_mat))
  h2 <- matrix(substr(gene_mat, 2L, 2L), nrow = nrow(gene_mat),
               dimnames = dimnames(gene_mat))
  arr <- array(NA_character_, dim = c(length(samples), length(markers), 2L),
               dimnames = list(samples, markers, c("1", "2")))
  arr[, , 1L] <- t(h1)   # markers x samples -> samples x markers
  arr[, , 2L] <- t(h2)
  arr
}


# ---------------------------------------------------------------------------
# Exported bridge functions
# ---------------------------------------------------------------------------

#' Load a GenoMeta into genomicSimulation
#'
#' Writes the GenoMeta (and effects) to a temporary geno/map/effect file trio
#' and loads them with \code{genomicSimulation::load.data()}, returning a
#' lightweight \emph{handle} that the other bridge functions consume. The
#' GenoMeta is kept in the handle as a template so map/block/effect metadata can
#' be reattached when offspring are read back.
#'
#' @param geno_meta A valid \code{GenoMeta} list (the founder / parent pool).
#' @param effect_df Optional \code{(Marker, Allele, Effect)} table for breeding
#'   values. Defaults to \code{geno_meta$marker_effects}; when that holds a
#'   named per-trait list, its first table is used.
#' @param work_dir  Directory for the temp files. Default a session tempdir.
#' @param prefix    File-name prefix. Default \code{"gs"}.
#' @param cM_per_mb Recombination rate used to convert physical map positions
#'   (bp) to the centimorgans genomicSimulation expects. Default \code{1}
#'   (1 cM/Mb). Positions that already look like a cM map (max \eqn{\le 10^4}) are
#'   passed through unchanged. Passing a bp map unscaled makes gS meiosis both
#'   catastrophically slow and genetically meaningless, so this conversion is on
#'   by default.
#'
#' @return A \code{gs_handle} list: \code{group} (genomicSimulation group
#'   number), \code{template} (the input GenoMeta), \code{effect_df},
#'   \code{work_dir}, \code{files}.
#'
#' @seealso \code{\link{gs_cross_plan}}, \code{\link{gs_group_to_geno_meta}}
#' @export
gs_load_geno_meta <- function(geno_meta, effect_df = NULL,
                              work_dir = tempdir(), prefix = "gs",
                              cM_per_mb = 1) {
  .gs_require()
  if (!is_geno_meta(geno_meta))
    stop("gs_load_geno_meta: 'geno_meta' is not a valid GenoMeta.")

  eff   <- .coerce_single_effect_df(effect_df %||% geno_meta$marker_effects)
  files <- .gs_write_files(geno_meta, eff, work_dir, prefix, cM_per_mb = cM_per_mb)

  # Specify the genotype-matrix layout explicitly rather than trusting
  # genomicSimulation's auto-detection: in gS 0.3 the auto-detector mis-reads the
  # header .gs_write_files() writes (col.names = NA) as "no header row" and drops
  # the genotype names (all NA), which then breaks the name->index lookup in
  # gs_cross_plan(). Our file is always: header row, markers as rows, phased
  # 2-char allele pairs. The format helper (and load.data's `format` argument)
  # were added together, so pass `format` only when the helper is available —
  # older gS versions keep the original auto-detecting call.
  args <- list(files[["geno"]], files[["map"]],
               if (is.na(files[["eff"]])) NULL else files[["eff"]])
  if (exists("define.matrix.format.details", where = asNamespace("genomicSimulation")))
    args$format <- genomicSimulation::define.matrix.format.details(
      has.header = TRUE, markers.as.rows = TRUE, cell.style = "P")

  loaded <- do.call(genomicSimulation::load.data, args)
  group <- if (is.list(loaded)) loaded$groupNum else loaded

  structure(
    list(group = group, template = geno_meta, effect_df = eff,
         work_dir = work_dir, files = files),
    class = "gs_handle")
}


#' Realise a crossing plan in genomicSimulation
#'
#' Maps a crossing plan (parent \emph{names}) to the live group's current
#' genomicSimulation indices and performs the crosses with
#' \code{make.targeted.crosses()} — one offspring per plan row, with each pair
#' replicated \code{n_offspring} times when that column is present (origin
#' behaviour). Offspring are renamed to their \code{CrossID}.
#'
#' @param handle  A \code{gs_handle} from \code{\link{gs_load_geno_meta}}.
#' @param plans   data.frame with \code{Parent1}, \code{Parent2} (sample names);
#'   optional \code{Cross_ID}/\code{CrossID} and \code{n_offspring} (as produced
#'   by \code{\link{make_cross_design}}).
#' @param offspring_prefix Prefix for generated CrossIDs when none is supplied.
#'   Default \code{"F"}.
#' @param advance Logical. When \code{TRUE} (default) the parent group is deleted
#'   and the returned handle points at the offspring group (recurrent selection).
#'
#' @return A list: \code{handle} (updated, pointing at the offspring group) and
#'   \code{plans} (the expanded plan with \code{CrossID} and \code{F1mean}).
#'
#' @seealso \code{\link{gs_load_geno_meta}}, \code{\link{make_cross_design}}
#' @export
gs_cross_plan <- function(handle, plans, offspring_prefix = "F",
                          advance = TRUE) {
  .gs_require()
  if (!inherits(handle, "gs_handle"))
    stop("gs_cross_plan: 'handle' must be a gs_handle from gs_load_geno_meta().")
  plans <- as.data.frame(plans, stringsAsFactors = FALSE)
  if (!all(c("Parent1", "Parent2") %in% names(plans)))
    stop("gs_cross_plan: 'plans' needs 'Parent1' and 'Parent2' columns.")

  # Resolve CrossID (accept CrossID or Cross_ID; else synthesise)
  cross_id <- if ("CrossID" %in% names(plans)) plans$CrossID
              else if ("Cross_ID" %in% names(plans)) plans$Cross_ID
              else paste0(offspring_prefix, sprintf("%04i", seq_len(nrow(plans))))

  # Expand by n_offspring (one offspring per resulting row)
  reps <- if ("n_offspring" %in% names(plans)) as.integer(plans$n_offspring)
          else rep(1L, nrow(plans))
  idx  <- rep(seq_len(nrow(plans)), times = reps)
  exp_plans <- plans[idx, , drop = FALSE]
  exp_plans$CrossID <- paste0(cross_id[idx], "n",
                              sprintf("%02i", sequence(reps)))

  # Map parent names -> current genomicSimulation X indices (recompute: indexes
  # shift after delete.group, so never cache across cycles).
  cur_names <- genomicSimulation::see.group.data(handle$group, "N")
  cur_x     <- genomicSimulation::see.group.data(handle$group, "X")
  p1x <- cur_x[match(exp_plans$Parent1, cur_names)]
  p2x <- cur_x[match(exp_plans$Parent2, cur_names)]
  if (anyNA(p1x) || anyNA(p2x))
    stop("gs_cross_plan: some plan parents are not in the current group.")

  off_group <- genomicSimulation::make.targeted.crosses(
    first.parents = p1x, second.parents = p2x)
  genomicSimulation::change.names.to.values(exp_plans$CrossID, off_group)
  exp_plans$F1mean <- genomicSimulation::see.group.data(off_group, "B")

  if (isTRUE(advance))
    genomicSimulation::delete.group(handle$group)

  handle$group <- off_group
  list(handle = handle, plans = exp_plans)
}


# Build an offspring FAM sample_info from a crossing plan, aligned to the
# gene matrix's sample (offspring) names. Returns a 5-column FAM data.frame
# (Family, Sample, Parent1, Parent2, Sex) or NULL when no plan is supplied (the
# caller then falls back to new_geno_meta()'s dummy founder FAM).
#
# Offspring are matched to plan rows by their CrossID (accept CrossID/Cross_ID),
# which equals the gene-matrix column names after gs_cross_plan()'s rename. Family
# is the parent-pair group index in first-appearance order — one (Parent1,Parent2)
# pair == one family, matching .stk_rename_generation()'s FamID (.GRP) semantics.
.offspring_fam_from_plan <- function(samples, cross_plan) {
  cp <- as.data.frame(cross_plan, stringsAsFactors = FALSE)
  if (!all(c("Parent1", "Parent2") %in% names(cp)))
    stop("update_geno_meta_from_matrix: 'cross_plan' needs 'Parent1' and ",
         "'Parent2' columns.")
  off_id <- if ("CrossID" %in% names(cp)) as.character(cp$CrossID)
            else if ("Cross_ID" %in% names(cp)) as.character(cp$Cross_ID)
            else stop("update_geno_meta_from_matrix: 'cross_plan' needs a ",
                      "'CrossID' (or 'Cross_ID') column to match offspring.")

  idx <- match(samples, off_id)
  if (anyNA(idx))
    stop("update_geno_meta_from_matrix: offspring not found in 'cross_plan' ",
         "(by CrossID): ", paste(samples[is.na(idx)], collapse = ", "))

  p1 <- as.character(cp$Parent1[idx])
  p2 <- as.character(cp$Parent2[idx])
  pair   <- paste(p1, p2, sep = "\r")             # \r never appears in an ID
  family <- as.character(match(pair, unique(pair)))  # first-appearance FamID

  data.frame(Family = family, Sample = as.character(samples),
             Parent1 = p1, Parent2 = p2, Sex = "0",
             stringsAsFactors = FALSE)
}

#' Rebuild a GenoMeta from a gene-data matrix and a template
#'
#' Converts a genomicSimulation-style gene-data matrix (markers x samples,
#' 2-character alleles — the layout returned by \code{see.group.gene.data()})
#' into a \code{GenoMeta}, reusing the map / block-map / encoding / effects from a
#' \code{template} GenoMeta. Crossing never changes the marker set, so the
#' template's marker metadata applies unchanged to the new generation.
#'
#' This is the genomicSimulation-free core of \code{\link{gs_group_to_geno_meta}}:
#' it takes a plain matrix rather than a live gS group, so it is unit-testable
#' without the genomicSimulation package. It is the reusable primitive behind the
#' persistent-session crossing backend (each cycle pulls genotypes in-memory and
#' rebuilds a GenoMeta with this function).
#'
#' When a \code{cross_plan} is supplied, the rebuilt GenoMeta describes the
#' \emph{progeny} population: each offspring's \code{sample_info} FAM row carries
#' its real \code{Parent1} / \code{Parent2} (from the plan) and a \code{Family}
#' index grouping full-sibs by parent pair. Without a plan the FAM defaults to the
#' founder-style dummy table (parents \code{"0"}), preserving the previous
#' behaviour.
#'
#' @param template A \code{GenoMeta} whose \code{map_info},
#'   \code{block_info} / \code{block_marker_map_info}, \code{metadata} (encoding,
#'   ploidy) and \code{marker_effects} are reused for the rebuilt GenoMeta.
#' @param gene_matrix A character matrix, markers (rows) x samples (columns),
#'   each cell a 2-character genotype (e.g. \code{"AT"}), with marker and sample
#'   names as row/column names.
#' @param cross_plan Optional crossing plan (as returned by
#'   \code{\link{gs_cross_plan}} / \code{\link{make_cross_design}}) with
#'   \code{Parent1}, \code{Parent2} and \code{CrossID} (or \code{Cross_ID})
#'   columns. When given, the offspring FAM (\code{sample_info}) is populated from
#'   it — offspring are matched to plan rows by \code{CrossID} (which equals the
#'   \code{gene_matrix} column names). Default \code{NULL} (dummy founder FAM).
#'
#' @return A \code{GenoMeta} carrying \code{gene_matrix}'s genotypes with the
#'   template's marker metadata reattached, and — when \code{cross_plan} is given —
#'   an offspring FAM describing each progeny's parents.
#'
#' @seealso \code{\link{gs_group_to_geno_meta}}, \code{\link{new_geno_meta}}
#' @export
update_geno_meta_from_matrix <- function(template, gene_matrix,
                                         cross_plan = NULL) {
  if (!is_geno_meta(template))
    stop("update_geno_meta_from_matrix: 'template' is not a valid GenoMeta.")
  if (!is.matrix(gene_matrix) || is.null(rownames(gene_matrix)) ||
      is.null(colnames(gene_matrix)))
    stop("update_geno_meta_from_matrix: 'gene_matrix' must be a markers x ",
         "samples matrix with row (marker) and column (sample) names.")

  arr     <- .gs_genes_to_array(gene_matrix)
  markers <- dimnames(arr)[[2L]]
  samples <- dimnames(arr)[[1L]]

  map_info <- as.data.frame(template$map_info, stringsAsFactors = FALSE)
  map_info <- map_info[match(markers, map_info$Marker), , drop = FALSE]

  bmm <- template$block_marker_map_info
  if (!is.null(bmm))
    bmm <- as.data.frame(bmm)[as.character(bmm$Marker) %in% markers, , drop = FALSE]

  # Offspring FAM from the plan (real parents); NULL -> dummy founder FAM.
  sample_info <- if (!is.null(cross_plan))
    .offspring_fam_from_plan(samples, cross_plan) else NULL

  meta_md <- template$metadata
  meta_md$n_samples <- dim(arr)[1L]
  meta_md$n_markers <- dim(arr)[2L]

  new_geno_meta(
    geno_arrays           = arr,
    map_info              = map_info,
    block_info            = template$block_info,
    block_marker_map_info = bmm,
    sample_info           = sample_info,
    metadata              = meta_md,
    marker_effects        = template$marker_effects)
}


#' Read a genomicSimulation group's genotypes back into a GenoMeta
#'
#' Pulls the live group's gene data with \code{see.group.gene.data()} and rebuilds
#' a \code{GenoMeta}, reattaching map / block-map / encoding / effects from the
#' handle's template (crossing does not change the marker set). A thin
#' genomicSimulation wrapper over \code{\link{update_geno_meta_from_matrix}}.
#'
#' @param handle A \code{gs_handle} (its \code{group} is read).
#' @param cross_plan Optional crossing plan (e.g. the \code{plans} element from
#'   \code{\link{gs_cross_plan}}) used to populate the offspring FAM with real
#'   parents; forwarded to \code{\link{update_geno_meta_from_matrix}}. Default
#'   \code{NULL} (dummy founder FAM).
#'
#' @return A \code{GenoMeta} for the current offspring group.
#'
#' @seealso \code{\link{update_geno_meta_from_matrix}}, \code{\link{gs_cross_plan}}
#' @export
gs_group_to_geno_meta <- function(handle, cross_plan = NULL) {
  .gs_require()
  if (!inherits(handle, "gs_handle"))
    stop("gs_group_to_geno_meta: 'handle' must be a gs_handle.")

  update_geno_meta_from_matrix(
    handle$template,
    genomicSimulation::see.group.gene.data(handle$group),
    cross_plan = cross_plan)
}


#' Breeding values for the live genomicSimulation group
#'
#' @param handle A \code{gs_handle}.
#' @return A named numeric vector of breeding values (names = genotype names).
#' @seealso \code{\link{gs_cross_plan}}
#' @export
gs_group_gebv <- function(handle) {
  .gs_require()
  if (!inherits(handle, "gs_handle"))
    stop("gs_group_gebv: 'handle' must be a gs_handle.")
  stats::setNames(
    genomicSimulation::see.group.data(handle$group, "B"),
    genomicSimulation::see.group.data(handle$group, "N"))
}


#' Rename the live genomicSimulation group to match compact GenoMeta IDs
#'
#' After the driver renames a generation to compact IDs, the persistent gS group
#' must carry the \emph{same} names so the next cycle's parent-name -> \code{"X"}
#' index lookup in \code{\link{gs_cross_plan}} resolves. Maps the group's current
#' names (\code{"N"}, in the group's internal order) through
#' \code{old_ids -> new_ids} and pushes the result with
#' \code{genomicSimulation::change.names.to.values()}. Names not found in
#' \code{old_ids} are left unchanged.
#'
#' @param handle  A \code{gs_handle} (its live \code{group} is renamed in place).
#' @param old_ids Character vector of the group's current names to remap.
#' @param new_ids Character vector of replacement names, parallel to \code{old_ids}.
#'
#' @return The \code{handle}, invisibly (the rename is a side effect on the group).
#' @seealso \code{\link{gs_cross_plan}}, \code{\link{gs_group_to_geno_meta}}
#' @export
gs_sync_group_names <- function(handle, old_ids, new_ids) {
  .gs_require()
  if (!inherits(handle, "gs_handle"))
    stop("gs_sync_group_names: 'handle' must be a gs_handle.")

  cur    <- genomicSimulation::see.group.data(handle$group, "N")
  mapped <- new_ids[match(cur, old_ids)]
  mapped[is.na(mapped)] <- cur[is.na(mapped)]           # keep unmatched as-is
  genomicSimulation::change.names.to.values(as.character(mapped), handle$group)
  invisible(handle)
}


#' Close a persistent genomicSimulation session
#'
#' Frees the handle's live group so a \code{(rep, scheme)} job leaves no groups
#' behind (matters for long / parallel runs). Safe to call on \code{NULL} or a
#' non-handle (no-op), so it can guard \code{on.exit()} unconditionally.
#'
#' @param handle A \code{gs_handle}, or \code{NULL}.
#' @return \code{NULL}, invisibly.
#' @seealso \code{\link{gs_load_geno_meta}}
#' @export
gs_close <- function(handle) {
  if (is.null(handle) || !inherits(handle, "gs_handle"))
    return(invisible(NULL))
  if (requireNamespace("genomicSimulation", quietly = TRUE) &&
      !is.null(handle$group))
    try(genomicSimulation::delete.group(handle$group), silent = TRUE)
  invisible(NULL)
}
