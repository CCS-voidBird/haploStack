# =============================================================================
# catalogue_main.R
#
# run_haplotype_catalogue() — the one-call catalogue for one or more traits
# (0.9.0, design D3 §3.2–§3.4):
#
#   1. Resolve the centring mode, the effects and the traits
#   2. Build block location info, the trait-invariant compute core and the
#      allele-combination signatures ONCE (an empty core stops here, before
#      any work or write)
#   3. Per trait (optionally in parallel): local GEBV matrix, per-sample GEBV,
#      block counts, the block-haplotype catalogue, the per-block summary, and
#      three figure RECIPES (nothing is drawn)
#   4. Return an `hs_catalogue` result (trait-keyed; notes and provenance as
#      attributes); when `output_dir` is given, write_outputs() reproduces the
#      0.8.8 file set in it and the result is returned invisibly
#
# write_outputs.hs_catalogue() is the writer: the same names and bytes whether
# called at compute time or later on a result kept in memory.
# =============================================================================

#' Run the haplotype catalogue pipeline
#'
#' Orchestrates the complete catalogue workflow for one or more traits. Builds
#' a haplotype profile from the supplied \code{GenoMeta}, computes per-block
#' local GEBVs, summarises counts, builds the block-haplotype catalogue and
#' prepares the per-trait figures. All input data (genotypes, map, block
#' assignments) are read from the \code{GenoMeta}. Marker effects may be passed
#' via \code{effect_list} or detected from \code{geno_meta$marker_effects} (for
#' example after \code{\link{fit_genomic_prediction}}).
#'
#' Since 0.9.0 the function \strong{computes and returns}; it writes files only
#' when \code{output_dir} is given. The result is an
#' \code{\link[=hs_result]{hs_catalogue}}: a list keyed by trait whose elements
#' are the 0.8.8 tables, plus unrendered figure recipes, with notes and
#' provenance reached through \code{\link{hs_notes}} and
#' \code{\link{hs_provenance}}. \code{\link{write_outputs}} writes the same
#' files \code{output_dir} would have, at any later time.
#'
#' The \code{GenoMeta} must already contain block assignments — call
#' \code{\link{attach_block_info}} before this function if
#' \code{block_marker_map_info} is \code{NULL}.
#'
#' Trait processing is embarrassingly parallel: set \code{n_cores > 1} to use
#' \pkg{future.apply} for parallel execution (requires \pkg{future} and
#' \pkg{future.apply}). Falls back to sequential with a warning if those
#' packages are absent.
#'
#' @section Output label:
#' File names carry a label where 0.8.8 stamped \code{<ld>_t<tol>} (by default
#' \code{0.5_t0}, whatever made the blocks). Since 0.9.0 the label comes from
#' the container's block provenance (\code{metadata$blocks}, written by
#' \code{\link{attach_block_info}} and the LD resolvers): LD-built blocks keep
#' their \code{<r2>_t<tol>} name; a user block map gives the map file's stem
#' (or \code{"usermap"} for an in-memory table); a container with no
#' provenance keeps \code{0.5_t0} and the result carries a note
#' (\code{label_legacy_default}). An explicit \code{label}, or the 0.8.8
#' \code{ld_threshold} / \code{tolerance} pair, overrides it.
#'
#' @section Files written:
#' With \code{output_dir} (or later through \code{write_outputs()}), per trait
#' \code{<T>} and label \code{<L>}: \code{Haplotype_GA_Matrix_<L>_<T>.txt},
#' \code{GEBV_<L>_<T>.txt}, \code{Block_Counts_<L>_<T>.txt},
#' \code{Haplotype_Catalogue_<L>_<T>.txt}, \code{Block_Summary_<L>_<T>.txt},
#' and the figures \code{Block_Counts_<T>.png}, \code{Manhattan_Variance_<L>_<T>.png},
#' \code{Manhattan_Effect_<L>_<T>.png}; plus \code{geno_meta_<L>.rds} once when
#' \code{save_profile = TRUE}. Trait names are sanitised for file names (a
#' changed stem is a note); two traits that collide after sanitising are an
#' error.
#'
#' @param geno_meta       A valid \code{GenoMeta} with
#'   \code{block_marker_map_info} populated.
#' @param effect_list     Marker effects to catalogue. Either a named list of
#'   effect data.frames (one per trait, each with columns Marker, Allele,
#'   Effect) or a single data.frame (optionally with a \code{Trait} column to
#'   split on). When \code{NULL} (default), \code{geno_meta$marker_effects} is
#'   used. An explicit \code{effect_list} always overrides the stored copy.
#' @param trait_names     Character vector of trait names to process. Must
#'   match names in the resolved \code{effect_list}. Default (\code{NULL}): all
#'   traits in the resolved effects.
#' @param output_dir      Directory to write the files into. Default
#'   \code{NULL}: nothing is written; the result is returned (visibly). Pass
#'   \code{"catalogue_output"} for the 0.8.8 behaviour.
#' @param selected_blocks Optional data.frame with column \code{Block} (or a
#'   character vector of block ids) restricting the local GEBV matrix to a
#'   subset of blocks. Default \code{NULL} (all blocks).
#' @param label           Output label (see \emph{Output label}). Default
#'   \code{NULL}: from the block provenance.
#' @param ld_threshold,tolerance The 0.8.8 label pair. Default \code{NULL};
#'   when either is given the label is \code{<ld>_t<tol>} (defaults 0.5 and 0
#'   for the missing one), overriding the provenance.
#' @param save_profile    Logical. Save the \code{GenoMeta} object to
#'   \code{output_dir} as \code{geno_meta_<L>.rds}. Default \code{TRUE}; only
#'   consulted when \code{output_dir} is given.
#' @param figures         Logical. Build the per-trait figure recipes (cheap;
#'   nothing is rendered until written). Default \code{TRUE}.
#' @param center        Centring mode applied to the allele effects:
#'   \code{"none"} (default), \code{"locus_mean"} or \code{"allele_freq"}. See
#'   \code{\link{compute_profile_localgebv}} for what each preserves.
#' @param center_by_freq  Deprecated logical form of \code{center}
#'   (\code{TRUE} = \code{"allele_freq"}, \code{FALSE} = \code{"none"}).
#' @param n_cores         Integer. Parallel workers for trait processing.
#'   \code{1L} (default) is sequential. Values \code{> 1} require
#'   \pkg{future} and \pkg{future.apply}.
#' @param progress        \code{NULL} (default) or a function of one argument,
#'   called after each trait with an event list (\code{stage = "catalogue"},
#'   \code{unit = "trait"}, \code{i}, \code{n}, \code{label}, \code{elapsed},
#'   \code{payload = list(trait)}). Returning \code{FALSE} cancels the run with
#'   \code{hs_error_cancelled}, whose \code{partial} field holds the traits
#'   finished so far.
#'
#' @return An \code{hs_catalogue} (see \code{\link{hs_result}}): a named list
#'   with one element per trait, each holding
#' \describe{
#'   \item{\code{gebv_matrix}}{data.table — wide Block × Sample GEBV matrix.}
#'   \item{\code{gebv_summary}}{data.table — per-sample total GEBV.}
#'   \item{\code{block_counts}}{data.table — per-sample positive/zero/negative
#'     block counts (\code{Sample} is character; \code{\link{plot_block_counts}}
#'     orders the bars).}
#'   \item{\code{hap_catalogue}}{data.table — long-format block-haplotype
#'     catalogue (one row per distinct block effect state with its
#'     allele-combination \code{signature}, population \code{Count} and block
#'     location). See \code{\link{compute_block_catalogue}}.}
#'   \item{\code{block_summary}}{data.table — one row per block:
#'     \code{Chrom}, \code{Pos}, \code{Block}, then \code{BlockVariance},
#'     \code{BlockMax}, \code{BlockMin}, \code{BlockMean}.}
#'   \item{\code{figures}}{the three figure recipes (\code{block_counts},
#'     \code{manhattan_variance}, \code{manhattan_effect}); build them with
#'     \code{hs_figures(res, build = TRUE)}.}
#' }
#'   Returned visibly, or invisibly when \code{output_dir} was written. The
#'   per-trait copy of the container that 0.8.8 returned is gone: the caller
#'   has it, and \code{write_outputs(geno_meta =)} writes it.
#'
#' @examples
#' \dontrun{
#' res <- run_haplotype_catalogue(meta, effect_list = list(Yield = yield_effects))
#' hs_notes(res); hs_provenance(res)$label
#' write_outputs(res, "catalogue_output", geno_meta = meta)
#' }
#'
#' @seealso \code{\link{write_outputs}}, \code{\link{attach_block_info}},
#'   \code{\link{compute_profile_localgebv}},
#'   \code{\link{hapS_ld_blocks}} (local LD blocks, when no block map exists)
#'
#' @export
run_haplotype_catalogue <- function(
  geno_meta,
  effect_list     = NULL,
  trait_names     = NULL,
  output_dir      = NULL,
  selected_blocks = NULL,
  label           = NULL,
  ld_threshold    = NULL,
  tolerance       = NULL,
  save_profile    = TRUE,
  figures         = TRUE,
  center          = NULL,
  center_by_freq  = NULL,
  n_cores         = 1L,
  progress        = NULL
) {
  center <- .resolve_center(center, center_by_freq, "run_haplotype_catalogue")
  .hs_check_meta(geno_meta, "run_haplotype_catalogue", "geno_meta",
                 hint = "Use a read_*() function and then attach_block_info() to prepare it.")
  if (is.null(geno_meta$block_marker_map_info))
    .hs_abort("run_haplotype_catalogue: GenoMeta has no block_marker_map_info. ",
         "Call attach_block_info() first.")
  if (!is.null(output_dir) &&
      (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir)))
    .hs_abort("run_haplotype_catalogue: 'output_dir' must be NULL or a single directory path.")

  # Effect source: an explicit effect_list overrides the GenoMeta's stored
  # marker_effects. Both are normalised to a named list of per-trait tables.
  if (is.null(effect_list)) effect_list <- geno_meta$marker_effects
  effect_list <- .as_effect_list(effect_list)
  if (is.null(effect_list) || length(effect_list) == 0L)
    .hs_abort("run_haplotype_catalogue: no marker effects supplied — pass ",
         "'effect_list' or attach 'marker_effects' to the GenoMeta ",
         "(e.g. via fit_genomic_prediction()).")
  if (is.null(trait_names)) trait_names <- names(effect_list)
  trait_names <- as.character(trait_names)
  missing_tr <- setdiff(trait_names, names(effect_list))
  if (length(missing_tr))
    .hs_abort("run_haplotype_catalogue: no effects for trait(s) ",
              paste(sQuote(missing_tr, FALSE), collapse = ", "), ".")

  blocks <- NULL
  if (!is.null(selected_blocks)) {
    blocks <- if (is.data.frame(selected_blocks)) {
      if (!"Block" %in% names(selected_blocks))
        .hs_abort("run_haplotype_catalogue: 'selected_blocks' needs a 'Block' column.")
      as.character(selected_blocks$Block)
    } else as.character(selected_blocks)
  }

  trk <- .hs_progress_tracker(progress, "catalogue", "trait", length(trait_names))

  run <- .hs_with_notes({
    lab   <- .hs_block_label(geno_meta, label, ld_threshold, tolerance)
    stems <- .hs_file_stems(trait_names)
    names(stems) <- trait_names

    # ── 1. Block location info ───────────────────────────────────────────────
    bloc_info_dt <- as.data.table(build_block_location_info(geno_meta)$block_info)

    # ── 2. Trait-invariant compute core (built ONCE for all traits) ──────────
    # The long-format allele table, allele frequencies, block/sample orderings
    # and (when Rcpp is active) the integer-coded matrices do not depend on the
    # marker effects. An empty core (no block-assigned marker in the genotype)
    # stops here, before any trait is computed or any file is written (cc/F16).
    # One geno_meta_views() per run (cc/F09, D3.10): the core and the block
    # signatures read the same views instead of melting the array twice.
    views        <- geno_meta_views(geno_meta)
    compute_core <- .catalogue_core(geno_meta, views = views)
    if (isTRUE(compute_core$empty))
      .hs_abort("run_haplotype_catalogue: no block-assigned marker matches the ",
                "genotype (the catalogue core is empty). Check that the block map's ",
                "marker ids match the container's.", call. = FALSE)
    compute_sigs <- .compute_block_signatures(geno_meta, views = views)
    rm(views)

    # ── 3. Per-trait analysis (tables + figure recipes; no writes) ───────────
    .run_one_trait <- function(trait) {
      .hs_inform("Processing trait: ", trait)
      effect_df <- effect_list[[trait]]
      stem <- stems[[trait]]

      gebv_matrix <- compute_profile_localgebv(
        geno_meta = geno_meta, effect_df = effect_df, blocks = blocks,
        level = "genotype", center = center, .core = compute_core)

      melt_gebv    <- .melt_gebv(gebv_matrix)
      gebv_summary <- melt_gebv[, .(GEBV = sum(Effect, na.rm = TRUE)), by = Sample]

      block_counts <- melt_gebv[!is.na(Effect), .(Count = .N), by = .(Sample, Label)]
      block_counts <- merge(block_counts, gebv_summary, by = "Sample", all.x = TRUE)
      block_counts$Sample <- as.character(block_counts$Sample)       # cc/F24
      block_counts$Label  <- factor(block_counts$Label,
                                    levels = c("Positive", "Zero", "Negative"))

      hap_catalogue <- compute_block_catalogue(
        geno_meta = geno_meta, effect_df = effect_df, block_info = bloc_info_dt,
        center = center, blocks = blocks,                # cc/F05: selected_blocks everywhere
        .core = compute_core, .signatures = compute_sigs)

      block_summary <- melt_gebv[!is.na(Effect), .(
        BlockVariance = stats::var(Effect, na.rm = TRUE),
        BlockMax      = max(Effect, na.rm = TRUE),
        BlockMin      = min(Effect, na.rm = TRUE),
        BlockMean     = mean(Effect, na.rm = TRUE)
      ), by = Block]
      block_summary <- merge(block_summary, bloc_info_dt[, .(Block, Chrom, Pos)],
                             by = "Block", all.x = TRUE)
      .lead_identifier_cols(block_summary)   # Chrom, Pos, Block lead the summary

      out <- list(gebv_matrix = gebv_matrix, gebv_summary = gebv_summary,
                  block_counts = block_counts, hap_catalogue = hap_catalogue,
                  block_summary = block_summary)
      if (isTRUE(figures)) {
        # Legacy render sizes (0.8.8 ggsave inches × 25.4) and file names.
        out$figures <- list(
          block_counts = .hs_figure(
            "Block_Counts", "plot_block_counts",
            list(haplotype_counts = .hs_ref("block_counts", trait = trait),
                 trait_name = trait),
            width_mm = 12 * 25.4, height_mm = 6 * 25.4, dpi = 150,
            file = paste0("Block_Counts_", stem, ".png")),
          manhattan_variance = .hs_figure(
            "Manhattan_Variance", "plot_manhattan_variance",
            list(variance_data = .hs_ref("block_summary", trait = trait),
                 trait_name = trait),
            width_mm = 12 * 25.4, height_mm = 5 * 25.4, dpi = 150,
            file = paste0("Manhattan_Variance_", lab, "_", stem, ".png")),
          manhattan_effect = .hs_figure(
            "Manhattan_Effect", "plot_manhattan_effect",
            list(effect_data = .hs_ref(NULL, trait = trait, derive = "manhattan_effect_data"),
                 trait_name = trait),
            width_mm = 12 * 25.4, height_mm = 5 * 25.4, dpi = 150,
            file = paste0("Manhattan_Effect_", lab, "_", stem, ".png")))
      }
      .hs_inform("  Trait ", trait, " completed.")
      out
    }

    # ── 4. Dispatch ──────────────────────────────────────────────────────────
    use_parallel <- .hs_parallel_ok(n_cores)    # D1-4: one rule, one plan owner

    n_tr <- length(trait_names)
    results_list <- vector("list", n_tr)
    names(results_list) <- trait_names
    if (use_parallel) {
      got <- .hs_with_plan(n_cores,
        future.apply::future_lapply(trait_names, .run_one_trait, future.seed = TRUE))
      for (i in seq_len(n_tr)) {
        results_list[[i]] <- got[[i]]
        trk$tick(i, trait_names[i], payload = list(trait = trait_names[i]),
                 partial = function() results_list[seq_len(i)])
      }
    } else {
      for (i in seq_len(n_tr)) {
        results_list[[i]] <- .run_one_trait(trait_names[i])
        trk$tick(i, trait_names[i], payload = list(trait = trait_names[i]),
                 partial = function() results_list[seq_len(i)])
      }
    }
    results_list
  })

  use_rcpp <- tryCatch(.hs_setting("use_rcpp", NULL), error = function(e) NULL)
  prov <- .hs_provenance(
    "run_haplotype_catalogue",
    center     = center,
    traits     = trait_names,
    file_stems = stems,
    label      = lab,
    blocks     = geno_meta$metadata$blocks,
    engine     = if (isTRUE(use_rcpp)) "rcpp" else if (isFALSE(use_rcpp)) "r" else NA_character_,
    settings   = list(use_rcpp = use_rcpp))
  res <- .hs_result(run$value, "hs_catalogue", notes = run$notes, provenance = prov)

  if (!is.null(output_dir)) {
    w <- .hs_with_notes(write_outputs(res, output_dir,
                                      geno_meta = if (isTRUE(save_profile)) geno_meta,
                                      overwrite = TRUE, create_dir = TRUE))
    if (NROW(w$notes)) {
      all_notes <- data.table::rbindlist(list(hs_notes(res), w$notes), fill = TRUE)
      attr(res, "hs_notes") <- .hs_collapse_notes(as.data.frame(all_notes))
    }
    .hs_inform("Catalogue written to: ", normalizePath(output_dir, winslash = "/"))
    return(invisible(res))
  }
  res
}

# Long form of a wide Block × Sample GEBV matrix with the sign label (the 0.8.8
# idiom, shared by the run and by the Manhattan-effect figure derivation).
.melt_gebv <- function(gebv_matrix) {
  m <- data.table::melt(as.data.table(gebv_matrix), id.vars = "Block",
                        variable.name = "Sample", value.name = "Effect")
  m[, Label := data.table::fcase(
    Effect > 0,  "Positive",
    Effect < 0,  "Negative",
    Effect == 0, "Zero"
  )]
  m
}

# The output label (D3 §3.2, form row 13 a). Explicit `label` wins; then the
# 0.8.8 `ld_threshold` / `tolerance` pair; then the block provenance: an LD
# resolver's `<r2>_t<tol>`, a user map's file stem, "usermap" for a table; no
# provenance at all keeps "0.5_t0" and leaves a note so old scripts on old
# containers write the same names.
.hs_block_label <- function(geno_meta, label = NULL, ld_threshold = NULL,
                            tolerance = NULL) {
  if (!is.null(label)) {
    label <- as.character(label)[1L]
    if (is.na(label) || !nzchar(label))
      .hs_abort("run_haplotype_catalogue: 'label' must be a non-empty string.",
                call. = FALSE)
    return(label)
  }
  if (!is.null(ld_threshold) || !is.null(tolerance))
    return(paste0(format(ld_threshold %||% 0.5), "_t", as.integer(tolerance %||% 0)))

  b <- geno_meta$metadata$blocks
  if (is.null(b) || is.null(b$source) || identical(as.character(b$source), "unknown")) {
    .hs_note("run_haplotype_catalogue", "label_legacy_default",
             paste0("the container carries no block provenance; the 0.8.8 output ",
                    "label '0.5_t0' is kept (pass label = to change it)"))
    return("0.5_t0")
  }
  p <- b$params
  if (!is.null(p$ld_threshold) && !is.null(p$tolerance))
    return(paste0(format(p$ld_threshold), "_t", as.integer(p$tolerance)))
  f <- b$file
  if (length(f) == 1L && !is.na(f) && nzchar(f))
    return(.hs_file_stem(tools::file_path_sans_ext(basename(f)), what = "label"))
  l <- b$label
  if (length(l) == 1L && !is.na(l) && nzchar(l) && !identical(l, "data.frame"))
    return(.hs_file_stem(tools::file_path_sans_ext(basename(l)), what = "label"))
  "usermap"
}

# ---- the writer --------------------------------------------------------------

#' @param what Which kinds of output to write: any of \code{"tables"},
#'   \code{"figures"}, \code{"container"} (default all three).
#' @param overwrite Replace existing files (default \code{FALSE} when called
#'   directly: an existing target is an error that names the files; the entry
#'   points pass \code{TRUE}).
#' @param geno_meta For an \code{hs_catalogue}: the container to write as
#'   \code{geno_meta_<label>.rds}; \code{NULL} (default) writes none.
#' @param create_dir Create \code{dir} when it does not exist (default
#'   \code{TRUE}).
#' @rdname write_outputs
#' @export
write_outputs.hs_catalogue <- function(x, dir, what = c("tables", "figures", "container"),
                                       overwrite = FALSE, geno_meta = NULL,
                                       create_dir = TRUE, ...) {
  what <- match.arg(what, several.ok = TRUE)
  dir  <- .hs_dir_ready(dir, create_dir = create_dir)
  prov <- hs_provenance(x)
  lab  <- prov$label %||% "0.5_t0"
  traits <- names(x)
  stems  <- prov$file_stems
  if (is.null(stems) || !all(traits %in% names(stems))) {
    stems <- .hs_file_stems(traits)
    names(stems) <- traits
  }
  tab_prefix <- c(gebv_matrix   = "Haplotype_GA_Matrix",
                  gebv_summary  = "GEBV",
                  block_counts  = "Block_Counts",
                  hap_catalogue = "Haplotype_Catalogue",
                  block_summary = "Block_Summary")

  # Plan every target first so overwrite = FALSE refuses before anything is written.
  plan <- list()
  for (tr in traits) {
    el <- x[[tr]]
    if (!is.list(el)) next
    if ("tables" %in% what) {
      for (nm in names(tab_prefix)) if (!is.null(el[[nm]]))
        plan[[length(plan) + 1L]] <- list(
          kind = "table", trait = tr, name = nm,
          path = file.path(dir, paste0(tab_prefix[[nm]], "_", lab, "_", stems[[tr]], ".txt")))
    }
    if ("figures" %in% what && length(el$figures)) {
      fig_names <- names(el$figures) %||% as.character(seq_along(el$figures))
      for (k in seq_along(el$figures)) {
        fig <- el$figures[[k]]
        fname <- if (inherits(fig, "hs_figure") && length(fig$file) == 1L &&
                     !is.na(fig$file) && nzchar(fig$file)) fig$file
                 else paste0(.hs_file_stem(fig_names[k], what = "figure"), "_",
                             stems[[tr]], ".png")
        plan[[length(plan) + 1L]] <- list(kind = "figure", trait = tr,
                                          name = fig_names[k], path = file.path(dir, fname))
      }
    }
  }
  if ("container" %in% what && !is.null(geno_meta)) {
    .hs_check_meta(geno_meta, "write_outputs", "geno_meta")
    plan[[length(plan) + 1L]] <- list(kind = "container", trait = NA_character_,
                                      name = "geno_meta",
                                      path = file.path(dir, paste0("geno_meta_", lab, ".rds")))
  }
  if (!length(plan)) return(invisible(.hs_manifest()))
  .hs_check_overwrite(vapply(plan, `[[`, character(1L), "path"), overwrite)

  rows <- vector("list", length(plan))
  for (k in seq_along(plan)) {
    p <- plan[[k]]
    object <- if (is.na(p$trait)) p$name else paste0(p$trait, "/", p$name)
    rows[[k]] <- switch(p$kind,
      table = .hs_write_table(x[[p$trait]][[p$name]], p$path, overwrite = TRUE,
                              object = object),
      figure = {
        fig <- x[[p$trait]]$figures[[p$name]]
        bytes <- tryCatch({
          obj <- if (inherits(fig, "hs_figure")) hs_build_figure(fig, result = x) else fig
          if (inherits(fig, "hs_figure"))
            hs_render_figure(obj, p$path, width_mm = fig$size$width_mm,
                             height_mm = fig$size$height_mm, dpi = fig$size$dpi,
                             overwrite = TRUE)
          else hs_render_figure(obj, p$path, overwrite = TRUE)
          as.numeric(file.info(p$path)$size)
        }, error = function(e) {
          .hs_note("write_outputs", "figure_failed",
                   paste0("figure ", object, " not written: ", conditionMessage(e)))
          NA_real_
        })
        .hs_manifest(normalizePath(p$path, winslash = "/", mustWork = FALSE), "figure",
                     object, bytes)
      },
      container = {
        save_geno_meta(geno_meta, p$path)
        .hs_manifest(normalizePath(p$path, winslash = "/", mustWork = TRUE), "container",
                     object, as.numeric(file.info(p$path)$size))
      })
  }
  invisible(do.call(rbind, rows))
}
