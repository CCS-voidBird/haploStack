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

# The installed genomicSimulation version (a function so a test can mock it).
.gs_version <- function() utils::packageVersion("genomicSimulation")

# genomicSimulation >= 0.3 (0.9.0, io/F21): the bridge writes the file layout
# and calls the matrix-format helper of 0.3; the auto-detecting branch for
# older versions is gone.
.gs_require <- function() {
  if (!requireNamespace("genomicSimulation", quietly = TRUE))
    .hs_abort(class = "hs_error_dependency",
         data = list(hint = "remotes::install_github(\"vllrs/genomicSimulation\")"),
         "This function needs the 'genomicSimulation' package (Suggests). ",
         "Install it, or use cross_backend = 'block'.", call. = FALSE)
  v <- .gs_version()
  if (v < "0.3")
    .hs_abort(class = "hs_error_dependency",
         data = list(hint = "remotes::install_github(\"vllrs/genomicSimulation\")"),
         "haploStack needs genomicSimulation >= 0.3 (installed: ", as.character(v),
         "): its file layout and matrix-format helper are what the bridge writes ",
         "and calls. Update it, or use cross_backend = 'block'.", call. = FALSE)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# The one-session guard (0.9.0, design D1 §3.6, form row 4 a)
#
# genomicSimulation keeps ONE simulation in the R process (its C library holds
# the loaded data; delete.group() reshuffles indices for every handle). Up to
# 0.8.8 a second gs_load_geno_meta() loaded on top of the first, and gs_close()
# on a handle from an earlier session deleted a group of the current one.
# The package now declares "one genomicSimulation job per R process":
#   - a second load while a session is open stops with hs_error_state, unless
#     replace = TRUE closes the open session first (for concurrent gS jobs run
#     one R process per job; the driver's n_cores > 1 does exactly that);
#   - every handle carries the session token it was made under; the four
#     bridge functions that act on the live group check it, so a stale handle
#     stops with hs_error_state instead of touching another job's groups;
#   - gs_close() frees a handle only when it is the live session, and removes
#     the private scratch directory the load created (work_dir = NULL).
# ---------------------------------------------------------------------------
.hs_gs <- new.env(parent = emptyenv())
.hs_gs$token   <- 0L
.hs_gs$open    <- FALSE
.hs_gs$own_dir <- NULL

# Close whatever session is open: every gS group and its data (clear.simdata,
# gS >= 0.3), then the private scratch. Used by replace = TRUE and by the
# session reset. Never errors.
.gs_session_clear <- function() {
  if (isTRUE(.hs_gs$open) && requireNamespace("genomicSimulation", quietly = TRUE))
    try(genomicSimulation::clear.simdata(), silent = TRUE)
  if (!is.null(.hs_gs$own_dir))
    unlink(.hs_gs$own_dir, recursive = TRUE, force = TRUE)
  .hs_gs$open    <- FALSE
  .hs_gs$own_dir <- NULL
  invisible(NULL)
}

# Is this handle the live session? (the check the four bridge functions run)
.gs_check_handle <- function(handle, fn) {
  if (!inherits(handle, "gs_handle"))
    .hs_abort(fn, ": 'handle' must be a gs_handle from gs_load_geno_meta().",
              call. = FALSE)
  live <- isTRUE(.hs_gs$open) && identical(as.integer(handle$token), .hs_gs$token)
  if (!live)
    .hs_abort(class = "hs_error_state",
              fn, ": this handle's genomicSimulation session was closed or ",
              "replaced (one session per R process; see ?gs_load_geno_meta).",
              call. = FALSE)
  invisible(TRUE)
}

# Does a session hold the process right now? (tests, the session reset)
.gs_session_open <- function() isTRUE(.hs_gs$open)

# Convert map positions to the centimorgans genomicSimulation expects. gS reads
# the map's position column as cM (recombination distance): the expected number
# of crossovers per gamete on a chromosome is proportional to its length in
# Morgans. A physical map (bp) passed through unscaled therefore becomes
# astronomically long (e.g. 99 Mb -> 99,000 Morgans), which makes meiosis
# simulation catastrophically slow AND genetically meaningless (every interval
# recombines). Three rules (0.9.0, io/F11; `backend_opts$map_unit`):
#   "auto"  (default, = 0.8.8) positions that look physical (max > 1e4,
#           implausible for a cM map) are scaled bp -> cM at `cM_per_mb`
#           (1 cM/Mb, the common plant rule of thumb), and the scaling is
#           announced as a note (hs_warning_fallback, code gs_map_scaled);
#           positions already in a cM range pass through unchanged;
#   "bp"    always scaled;
#   "cM"    never scaled.
.gs_positions_to_cM <- function(pos, cM_per_mb = 1, map_unit = c("auto", "bp", "cM")) {
  map_unit <- match.arg(map_unit)
  pos <- as.numeric(pos)
  if (length(pos) == 0L || all(is.na(pos))) return(pos)
  scale <- switch(map_unit,
                  auto = max(pos, na.rm = TRUE) > 1e4,
                  bp   = TRUE,
                  cM   = FALSE)
  if (scale) {
    if (map_unit == "auto")
      .hs_note("gs_load", "gs_map_scaled",
               paste0("map positions look physical (max ", format(max(pos, na.rm = TRUE)),
                      "); scaled bp -> cM at ", format(cM_per_mb),
                      " cM/Mb for genomicSimulation (backend_opts$map_unit = \"bp\" ",
                      "or \"cM\" states the unit)."))
    pos <- pos / 1e6 * cM_per_mb   # bp -> cM
  }
  pos
}

# Write a GenoMeta (+ effects) as the genotype / map / effect file trio that
# genomicSimulation::load.data() reads. Returns the three paths.
#   geno : markers x samples, two-character alleles, row+col names; a missing
#          call (either copy NA, the container's NA-any rule) is "--",
#          genomicSimulation's own missing form (form row 22 a; up to 0.8.8 the
#          text was "NA", which gS read as the alleles N and A — the A then
#          segregated).
#   map  : marker  chr  pos(cM)     with the header names gS 0.3 recognises
#   eff  : marker  allele  eff      (0.9.0, I2: "Marker Chrom Pos" was read as
#          data and logged "Failed to parse 1 line(s)" on every load).
.gs_write_files <- function(geno_meta, effect_df = NULL,
                            work_dir = tempdir(), prefix = "gs",
                            cM_per_mb = 1, use_meta_effects = TRUE,
                            map_unit = "auto") {
  .hs_check_meta(geno_meta, ".gs_write_files", "geno_meta")
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

  geno_file <- file.path(work_dir, paste0(prefix, "_geno.txt"))
  map_file  <- file.path(work_dir, paste0(prefix, "_map.txt"))
  eff_file  <- file.path(work_dir, paste0(prefix, "_eff.txt"))

  # Genotype matrix in genomicSimulation orientation (markers x samples, "AT"),
  # built from the two allele layers; a cell with either copy missing is "--"
  # (see the missing-vs-unknown note above .gs_genes_to_array()).
  arr <- geno_meta$geno_arrays
  if (dim(arr)[3L] != 2L)
    .hs_abort(".gs_write_files: genomicSimulation takes diploid genotypes (two ",
              "alleles per cell); this container has ploidy ", dim(arr)[3L], ".",
              call. = FALSE)
  a1 <- arr[, , 1L]; a2 <- arr[, , 2L]
  cell <- ifelse(is.na(a1) | is.na(a2), "--", paste0(a1, a2))
  dim(cell) <- dim(arr)[1:2]
  dimnames(cell) <- dimnames(arr)[1:2]
  mat <- t(cell)                               # markers x samples
  utils::write.table(mat, geno_file, sep = "\t", quote = FALSE,
                     col.names = NA)          # leading blank header cell

  # Map: marker chr pos(cM). Positions written non-scientific so gS parses them.
  map_out <- data.frame(
    marker = as.character(geno_meta$map_info$Marker),
    chr    = as.character(geno_meta$map_info$Chrom),
    pos    = format(.gs_positions_to_cM(geno_meta$map_info$Pos, cM_per_mb, map_unit),
                    scientific = FALSE, trim = TRUE),
    stringsAsFactors = FALSE)
  utils::write.table(map_out, map_file, sep = "\t", quote = FALSE,
                     row.names = FALSE)

  # Effects: marker allele eff (NULL -> no effect file, F1mean is NA). The
  # stored marker_effects may be a named per-trait list; collapse to a single
  # table. gs_load_geno_meta() resolves the trait itself and passes
  # use_meta_effects = FALSE, so a NULL there means "no effects".
  eff <- .coerce_single_effect_df(
    if (isTRUE(use_meta_effects)) effect_df %||% geno_meta$marker_effects else effect_df)
  if (!is.null(eff)) {
    eff <- as.data.frame(eff[, 1:3, drop = FALSE], stringsAsFactors = FALSE)
    names(eff) <- c("marker", "allele", "eff")
    utils::write.table(eff, eff_file, sep = "\t", quote = FALSE,
                       row.names = FALSE)
  } else {
    eff_file <- NULL
  }

  c(geno = geno_file, map = map_file, eff = if (is.null(eff_file)) NA else eff_file)
}

# ---------------------------------------------------------------------------
# Missing vs unknown at the genomicSimulation boundary (owner, 2026-10-27)
#
#   haploStack "missing"   a GenoMeta cell that is NA. It crosses into gS as the
#                          PLACEHOLDER allele pair "--" written by
#                          .gs_write_files() (form row 22 a, D4-7: gS's own
#                          missing form, for a cell with either copy missing).
#                          gS treats the placeholder as a literal allele: it
#                          propagates through crosses and carries no effect, so
#                          it behaves as a missing call and reads back as the
#                          same character.
#   gS "missing" (true)    an allele gS itself could not load (stored as '\0').
#                          It also propagates and carries no effect, but it is
#                          gS's, not ours: the package has no definition for it
#                          yet. On read-back it is requested under the sentinel
#                          .GS_TRUE_MISSING so it cannot be confused with the
#                          placeholder, it is mapped to NA, and the run gets a
#                          note (hs_warning_fallback, code "gs_true_missing").
#                          PLACEHOLDER: real-missing handling (what a gS missing
#                          should mean for the container) is not designed; the
#                          warning marks the risk until it is.
#
# Never read gene data with gS's default unknown.allele = "-": a true missing
# and the "-" placeholder would then be indistinguishable.
# ---------------------------------------------------------------------------
.GS_TRUE_MISSING <- "?"

# The one gene-data reader: every read-back goes through it so the sentinel is
# always in force.
.gs_group_genes <- function(group) {
  genomicSimulation::see.group.gene.data(group, unknown.allele = .GS_TRUE_MISSING)
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
  # gS's own missing alleles (see the definitions above): count, warn, map to NA.
  n_true <- sum(arr == .GS_TRUE_MISSING, na.rm = TRUE)
  if (n_true > 0L) {
    .hs_note("gs_read_back", "gs_true_missing",
             paste0("genomicSimulation reported ", n_true, " allele(s) it holds as ",
                    "missing (its own missing, not haploStack's placeholder); they ",
                    "are read as NA. Real-missing handling is not implemented: ",
                    "check the input the session was loaded from."),
             n = n_true)
    arr[arr == .GS_TRUE_MISSING] <- NA_character_
  }
  # A missing founder allele goes out to genomicSimulation as a placeholder
  # allele pair and is inherited by offspring as a literal allele; inside a
  # GenoMeta a missing allele is NA and nothing else (0.6.0 contract), so restore
  # it here — otherwise the offspring container carries the placeholder as an
  # allele and is_geno_meta() warns on every cycle of a gS-backend run (0.7.5).
  # The placeholder ("-" from D4-7, "N" from earlier builds) is in the nucleotide
  # missing set.
  arr[arr %in% .missing_tokens("nucleotide")] <- NA_character_
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
#' @section One session per R process (0.9.0):
#' genomicSimulation keeps one simulation per R process, so haploStack keeps
#' one session: a second call while a session is open stops with an
#' \code{hs_error_state} error, unless \code{replace = TRUE} closes the open
#' session first. Every handle carries the token of the session it was made
#' under, and \code{\link{gs_cross_plan}}, \code{\link{gs_group_to_geno_meta}},
#' \code{\link{gs_group_gebv}} and \code{\link{gs_sync_group_names}} refuse a
#' handle whose session was closed or replaced. Close a session with
#' \code{\link{gs_close}}. For concurrent genomicSimulation jobs run one R
#' process per job (\code{run_stacking_simulation(n_cores > 1)} does).
#'
#' @param geno_meta A valid \code{GenoMeta} list (the founder / parent pool).
#' @param effect_df Optional \code{(Marker, Allele, Effect)} table for breeding
#'   values. Defaults to \code{geno_meta$marker_effects}; when that holds a
#'   named per-trait list, its first table is used.
#' @param work_dir  Directory for the genotype / map / effect files
#'   genomicSimulation loads. \code{NULL} (default since 0.9.0) creates a
#'   private directory under \code{tempdir()} that \code{\link{gs_close}}
#'   removes; a directory you pass is never deleted.
#' @param prefix    File-name prefix. Default \code{"gs"}.
#' @param cM_per_mb Recombination rate used to convert physical map positions
#'   (bp) to the centimorgans genomicSimulation expects. Default \code{1}
#'   (1 cM/Mb). Positions that already look like a cM map (max \eqn{\le 10^4}) are
#'   passed through unchanged. Passing a bp map unscaled makes gS meiosis both
#'   catastrophically slow and genetically meaningless, so this conversion is on
#'   by default.
#' @param replace   \code{FALSE} (default): a session already open is an error.
#'   \code{TRUE}: close it (all its groups and its private files) and open this
#'   one; handles of the old session become invalid.
#' @param trait     When the effects hold several traits (a named list, or one
#'   long table with a \code{Trait} column), the one to load; with several
#'   traits and \code{trait = NULL} the call stops (0.9.0; up to 0.8.8 the
#'   first trait was taken silently). Default \code{NULL}.
#' @param map_unit  The unit of \code{map_info$Pos}: \code{"auto"} (default,
#'   the 0.8.8 rule: positions above \code{1e4} are taken as bp and scaled,
#'   with a note \code{gs_map_scaled}), \code{"bp"} (always scaled at
#'   \code{cM_per_mb}) or \code{"cM"} (never scaled).
#'
#' @return A \code{gs_handle} list: \code{group} (genomicSimulation group
#'   number), \code{template} (the input GenoMeta), \code{effect_df},
#'   \code{work_dir}, \code{files}, \code{token} (the session), \code{own_dir}
#'   (whether \code{work_dir} is the private directory).
#'
#' @seealso \code{\link{gs_cross_plan}}, \code{\link{gs_group_to_geno_meta}},
#'   \code{\link{gs_close}}
#' @export
gs_load_geno_meta <- function(geno_meta, effect_df = NULL,
                              work_dir = NULL, prefix = "gs",
                              cM_per_mb = 1, replace = FALSE, trait = NULL,
                              map_unit = c("auto", "bp", "cM")) {
  .gs_require()
  map_unit <- match.arg(map_unit)
  .hs_check_meta(geno_meta, "gs_load_geno_meta", "geno_meta")
  if (.gs_session_open()) {
    if (!isTRUE(replace))
      .hs_abort(class = "hs_error_state",
                "gs_load_geno_meta: genomicSimulation holds one session per R ",
                "process and one is open; close it with gs_close() first, or ",
                "pass replace = TRUE. For concurrent genomicSimulation jobs run ",
                "one R process per job.", call. = FALSE)
    .gs_session_clear()
  }
  own_dir <- is.null(work_dir)
  if (own_dir) work_dir <- tempfile("hs_gs_")

  # One trait's table, or none (row 19 a: several traits need `trait`).
  eff   <- .stk_resolve_effects(effect_df, geno_meta, trait, "gs_load_geno_meta",
                                allow_none = TRUE)
  files <- .gs_write_files(geno_meta, eff, work_dir, prefix, cM_per_mb = cM_per_mb,
                           use_meta_effects = FALSE, map_unit = map_unit)

  # Specify the genotype-matrix layout explicitly rather than trusting
  # genomicSimulation's auto-detection: in gS 0.3 the auto-detector mis-reads the
  # header .gs_write_files() writes (col.names = NA) as "no header row" and drops
  # the genotype names (all NA), which then breaks the name->index lookup in
  # gs_cross_plan(). Our file is always: header row, markers as rows, phased
  # 2-char allele pairs. gS >= 0.3 is required (.gs_require()), so the format
  # helper is always there.
  args <- list(files[["geno"]], files[["map"]],
               if (is.na(files[["eff"]])) NULL else files[["eff"]],
               format = genomicSimulation::define.matrix.format.details(
                 has.header = TRUE, markers.as.rows = TRUE, cell.style = "P"))

  loaded <- tryCatch(do.call(genomicSimulation::load.data, args),
                     error = function(e) {
                       if (own_dir) unlink(work_dir, recursive = TRUE, force = TRUE)
                       stop(e)
                     })
  group <- if (is.list(loaded)) loaded$groupNum else loaded

  # The session is open from here: next token, private dir on record.
  .hs_gs$token   <- .hs_gs$token + 1L
  .hs_gs$open    <- TRUE
  .hs_gs$own_dir <- if (own_dir) work_dir else NULL

  structure(
    list(group = group, template = geno_meta, effect_df = eff,
         work_dir = work_dir, files = files,
         token = .hs_gs$token, own_dir = own_dir),
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
  .gs_check_handle(handle, "gs_cross_plan")
  plans <- as.data.frame(plans, stringsAsFactors = FALSE)
  if (!all(c("Parent1", "Parent2") %in% names(plans)))
    .hs_abort("gs_cross_plan: 'plans' needs 'Parent1' and 'Parent2' columns.")

  # One expansion per offspring with the block backend's naming (CrossID /
  # Cross_ID accepted, else synthesised; replicated by n_offspring) — the shared
  # .expand_cross_plan() (io/F19), keeping the plan's own columns.
  exp_plans <- .expand_cross_plan(plans, offspring_prefix, keep_cols = TRUE)

  # Map parent names -> current genomicSimulation X indices (recompute: indexes
  # shift after delete.group, so never cache across cycles).
  cur_names <- genomicSimulation::see.group.data(handle$group, "N")
  cur_x     <- genomicSimulation::see.group.data(handle$group, "X")
  p1x <- cur_x[match(exp_plans$Parent1, cur_names)]
  p2x <- cur_x[match(exp_plans$Parent2, cur_names)]
  if (anyNA(p1x) || anyNA(p2x))
    .hs_abort("gs_cross_plan: some plan parents are not in the current group.")

  off_group <- genomicSimulation::make.targeted.crosses(
    first.parents = p1x, second.parents = p2x)
  genomicSimulation::change.names.to.values(exp_plans$CrossID, off_group)
  # No effect file loaded -> no breeding value: NA, not 0 (io/F10).
  exp_plans$F1mean <- if (is.na(handle$files[["eff"]])) NA_real_
                      else genomicSimulation::see.group.data(off_group, "B")

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
    .hs_abort("update_geno_meta_from_matrix: 'cross_plan' needs 'Parent1' and ",
         "'Parent2' columns.")
  off_id <- if ("CrossID" %in% names(cp)) as.character(cp$CrossID)
            else if ("Cross_ID" %in% names(cp)) as.character(cp$Cross_ID)
            else .hs_abort("update_geno_meta_from_matrix: 'cross_plan' needs a ",
                      "'CrossID' (or 'Cross_ID') column to match offspring.")

  idx <- match(samples, off_id)
  if (anyNA(idx))
    .hs_abort("update_geno_meta_from_matrix: offspring not found in 'cross_plan' ",
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
  .hs_check_meta(template, "update_geno_meta_from_matrix", "template")
  if (!is.matrix(gene_matrix) || is.null(rownames(gene_matrix)) ||
      is.null(colnames(gene_matrix)))
    .hs_abort("update_geno_meta_from_matrix: 'gene_matrix' must be a markers x ",
         "samples matrix with row (marker) and column (sample) names.")

  arr     <- .gs_genes_to_array(gene_matrix)
  markers <- dimnames(arr)[[2L]]
  samples <- dimnames(arr)[[1L]]

  if (!all(markers %in% as.character(template$map_info$Marker)))
    .hs_abort("update_geno_meta_from_matrix: ", sum(!markers %in% template$map_info$Marker),
         " marker(s) of 'gene_matrix' are not in the template's map.")

  # Offspring FAM from the plan (real parents); NULL -> dummy founder FAM.
  sample_info <- if (!is.null(cross_plan))
    .offspring_fam_from_plan(samples, cross_plan) else NULL

  # One rebuild path (0.9.0): the template's map, blocks, effects and
  # annotation carry over (blocks do not change across cycles, #033); the
  # founders' phenotype and family history are declared losses (X7) — same
  # rule as the block-backend offspring builder.
  .rebuild_geno_meta(template, geno_arrays = arr, sample_info = sample_info,
                     .drop = c("family_history", "phenotype"),
                     .caller = "update_geno_meta_from_matrix")
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
  .gs_check_handle(handle, "gs_group_to_geno_meta")

  update_geno_meta_from_matrix(
    handle$template,
    .gs_group_genes(handle$group),
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
  .gs_check_handle(handle, "gs_group_gebv")
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
  .gs_check_handle(handle, "gs_sync_group_names")

  cur    <- genomicSimulation::see.group.data(handle$group, "N")
  mapped <- new_ids[match(cur, old_ids)]
  mapped[is.na(mapped)] <- cur[is.na(mapped)]           # keep unmatched as-is
  genomicSimulation::change.names.to.values(as.character(mapped), handle$group)
  invisible(handle)
}


#' Close a persistent genomicSimulation session
#'
#' Frees the handle's live group so a \code{(rep, scheme)} job leaves no groups
#' behind (matters for long / parallel runs), marks the session closed and
#' removes the private scratch directory \code{\link{gs_load_geno_meta}}
#' created (a \code{work_dir} you passed is kept). Safe to call on \code{NULL},
#' a non-handle or a handle of a session already closed or replaced (no-op),
#' so it can guard \code{on.exit()} unconditionally. Idempotent.
#'
#' @param handle A \code{gs_handle}, or \code{NULL}.
#' @return \code{NULL}, invisibly.
#' @seealso \code{\link{gs_load_geno_meta}}
#' @export
gs_close <- function(handle) {
  if (is.null(handle) || !inherits(handle, "gs_handle"))
    return(invisible(NULL))
  # Only the live session's handle acts (0.9.0): a stale handle must not
  # delete a group of the session that replaced it.
  live <- isTRUE(.hs_gs$open) && identical(as.integer(handle$token), .hs_gs$token)
  if (!live) return(invisible(NULL))
  if (requireNamespace("genomicSimulation", quietly = TRUE) &&
      !is.null(handle$group))
    try(genomicSimulation::delete.group(handle$group), silent = TRUE)
  if (!is.null(.hs_gs$own_dir))
    unlink(.hs_gs$own_dir, recursive = TRUE, force = TRUE)
  .hs_gs$open    <- FALSE
  .hs_gs$own_dir <- NULL
  invisible(NULL)
}
