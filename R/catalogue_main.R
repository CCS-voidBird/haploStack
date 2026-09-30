# =============================================================================
# catalogue_main.R
#
# run_haplotype_catalogue() — high-level entry point that orchestrates the
# full catalogue pipeline for one or more traits:
#
#   1. Load or build block map
#   2. Build block location info
#   3. Read and load genotype data
#   4. Create haplotype profile
#   5. For each trait (optionally in parallel):
#      a. Compute local GEBV matrix
#      b. Compute per-sample GEBV and block counts
#      c. Produce Manhattan plots and count plots
#   6. Return results invisibly
# =============================================================================

#' Run the haplotype catalogue pipeline
#'
#' Orchestrates the complete catalogue workflow for one or more traits. Builds
#' a haplotype profile from the supplied \code{GenoMeta}, computes per-block
#' local GEBVs, summarises counts, and generates visualisation plots for each
#' trait. All input data (genotypes, map, block assignments) are read directly
#' from the \code{GenoMeta} fields. Marker effects may be passed explicitly via
#' \code{effect_list} or detected from \code{geno_meta$marker_effects} (for
#' example after \code{\link{fit_genomic_prediction}}); effects are required to
#' run, but never to construct a valid \code{GenoMeta}.
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
#' @param geno_meta       A valid \code{GenoMeta} with
#'   \code{block_marker_map_info} populated. Use \code{\link{read_st_genotype}}
#'   (or another reader) followed by \code{\link{attach_block_info}} to prepare
#'   the input.
#' @param effect_list     Marker effects to catalogue. Either a named list of
#'   effect data.frames (one per trait, each with columns Marker, Allele,
#'   Effect) or a single data.frame (optionally with a \code{Trait} column to
#'   split on). When \code{NULL} (default), \code{geno_meta$marker_effects} is
#'   used — so effects attached by \code{\link{fit_genomic_prediction}} are
#'   picked up automatically. An explicit \code{effect_list} always overrides
#'   the stored copy.
#' @param trait_names     Character vector of trait names to process. Must
#'   match names in the resolved \code{effect_list}. Default (\code{NULL}): all
#'   traits in the resolved effects.
#' @param output_dir      Directory for output files and plots. Created if it
#'   does not exist. Default \code{"catalogue_output"}.
#' @param selected_blocks Optional data.frame with column \code{Block}
#'   restricting analysis to a subset of blocks. Default \code{NULL} (all
#'   blocks).
#' @param ld_threshold    LD threshold label used in output file names only.
#'   Default \code{0.5}.
#' @param tolerance       Tolerance label used in output file names only.
#'   Default \code{0}.
#' @param save_profile    Logical. Save the \code{GenoMeta} object to
#'   \code{output_dir} as an RDS file. Default \code{TRUE}.
#' @param center        Centring mode applied to the allele effects:
#'   \code{"none"}, \code{"locus_mean"} or \code{"allele_freq"}. See
#'   \code{\link{compute_profile_localgebv}} for what each preserves.
#' @param center_by_freq  Logical. Mean-centre allele effects by population
#'   frequency (catalogue-style). Default \code{TRUE}.
#' @param n_cores         Integer. Parallel workers for trait processing.
#'   \code{1L} (default) is sequential. Values \code{> 1} require
#'   \pkg{future} and \pkg{future.apply}.
#'
#' @return Invisibly, a named list with one element per trait containing:
#' \describe{
#'   \item{\code{gebv_matrix}}{data.table — wide Block × Sample GEBV matrix.}
#'   \item{\code{gebv_summary}}{data.table — per-sample total GEBV.}
#'   \item{\code{block_counts}}{data.table — per-sample positive/zero/negative
#'     block counts.}
#'   \item{\code{hap_catalogue}}{data.table — long-format block-haplotype
#'     catalogue (one row per distinct block effect state with its allele-combination
#'     \code{signature}, population \code{Count} and block location; identification
#'     columns lead), written as
#'     \code{Haplotype_Catalogue_<ld>_t<tol>_<trait>.txt}. See
#'     \code{\link{compute_block_catalogue}}.}
#'   \item{\code{block_summary}}{data.table — one row per block summarising its
#'     across-sample local-GEBV effect, identification columns first
#'     (\code{Chrom}, \code{Pos}, \code{Block}) then \code{BlockVariance},
#'     \code{BlockMax}, \code{BlockMin}, \code{BlockMean}, written as
#'     \code{Block_Summary_<ld>_t<tol>_<trait>.txt}. Drives the per-trait
#'     variance scatter (\code{\link{plot_manhattan_variance}}).}
#'   \item{\code{geno_meta}}{The \code{GenoMeta} object analysed.}
#' }
#'
#' @examples
#' \dontrun{
#' # Prepare GenoMeta with block info
#' meta <- read_st_genotype("geno.csv")
#' meta <- attach_block_info(meta, read.csv("block_marker_map.csv"))
#'
#' results <- run_haplotype_catalogue(
#'   geno_meta   = meta,
#'   effect_list = list(Yield = yield_effects, Protein = prot_effects),
#'   output_dir  = "catalogue_output",
#'   n_cores     = 2L
#' )
#' }
#'
#' @seealso \code{\link{attach_block_info}},
#'   \code{\link{compute_profile_localgebv}},
#'   \code{\link{hapS_ld_blocks}} (local LD blocks, when no block map exists)
#'
#' @export
run_haplotype_catalogue <- function(
  geno_meta,
  effect_list     = NULL,
  trait_names     = NULL,
  output_dir      = "catalogue_output",
  selected_blocks = NULL,
  ld_threshold    = 0.5,
  tolerance       = 0,
  save_profile    = TRUE,
  center          = NULL,
  center_by_freq  = NULL,
  n_cores         = 1L
) {
  # CEN / FR10: the default centring is "none" from 0.6.0, matching
  # compute_profile_localgebv(). Until now the two entry points to the same
  # engine disagreed about what the default answer was. This CHANGES PUBLISHED
  # BLOCK VARIANCES (median factor ~3 on the barley panel);
  # center = "allele_freq" reproduces the 0.5.5 numbers.
  center <- .resolve_center(center, center_by_freq,
                            "run_haplotype_catalogue", default = "none")
  if (!is_geno_meta(geno_meta))
    stop("run_haplotype_catalogue: 'geno_meta' must be a valid GenoMeta. ",
         "Use a read_*() function and then attach_block_info() to prepare it.")
  if (is.null(geno_meta$block_marker_map_info))
    stop("run_haplotype_catalogue: GenoMeta has no block_marker_map_info. ",
         "Call attach_block_info() first.")

  # Effect source: an explicit effect_list overrides the GenoMeta's stored
  # marker_effects. Both are normalised to a named list of per-trait tables, so
  # effects attached by fit_genomic_prediction() are picked up automatically.
  if (is.null(effect_list)) effect_list <- geno_meta$marker_effects
  effect_list <- .as_effect_list(effect_list)
  if (is.null(effect_list) || length(effect_list) == 0L)
    stop("run_haplotype_catalogue: no marker effects supplied — pass ",
         "'effect_list' or attach 'marker_effects' to the GenoMeta ",
         "(e.g. via fit_genomic_prediction()).")
  if (is.null(trait_names)) trait_names <- names(effect_list)

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  # ── 1. Block location info ─────────────────────────────────────────────────
  bloc_result  <- build_block_location_info(geno_meta)
  bloc_info_dt <- as.data.table(bloc_result$block_info)

  # ── 2. Selected blocks ─────────────────────────────────────────────────────
  if (is.null(selected_blocks)) {
    selected_blocks <- data.table(
      Block    = unique(geno_meta$block_marker_map_info$Block),
      Resource = "All"
    )
  }

  # ── 2b. Trait-invariant compute core (built ONCE for all traits) ───────────
  # The long-format allele table, allele frequencies, block/sample orderings and
  # (when Rcpp is active) the integer-coded matrices do not depend on the marker
  # effects, so they are computed once here and shared by every trait — instead
  # of being rebuilt 2x per trait inside compute_profile_localgebv() and
  # compute_block_catalogue(). Captured by the .run_one_trait() closure below
  # (and serialised once to each future worker).
  compute_core <- .catalogue_core(geno_meta)

  # Allele-combination signatures are trait-invariant too — compute once and
  # share across traits so compute_block_catalogue() can label each row's subtype
  # without recomputing per trait.
  compute_sigs <- .compute_block_signatures(geno_meta)

  # ── 3. Persist the GenoMeta object ─────────────────────────────────────────
  if (save_profile) {
    gm_file <- file.path(output_dir,
      paste0("geno_meta_", ld_threshold, "_t", tolerance, ".rds"))
    save_geno_meta(geno_meta, gm_file)
    message("GenoMeta saved to: ", gm_file)
  }

  # ── 4. Per-trait analysis ──────────────────────────────────────────────────
  # Inner function captures geno_meta, bloc_info_dt, effect_list, output_dir,
  # ld_threshold, tolerance, selected_blocks, center_by_freq by closure.
  .run_one_trait <- function(trait) {
    message("Processing trait: ", trait)

    effect_df <- effect_list[[trait]]
    if (is.null(effect_df))
      stop("run_haplotype_catalogue: No effects for trait '", trait, "'.")

    # Compute local GEBV matrix
    gebv_matrix <- compute_profile_localgebv(
      geno_meta      = geno_meta,
      effect_df      = effect_df,
      blocks         = if (!is.null(selected_blocks)) selected_blocks$Block else NULL,
      level          = "genotype",
      center         = center,
      .core          = compute_core
    )

    # Save matrix
    gebv_file <- file.path(output_dir,
      paste0("Haplotype_GA_Matrix_", ld_threshold, "_t", tolerance,
             "_", trait, ".txt"))
    fwrite(gebv_matrix, file = gebv_file, sep = "\t")

    # Melt for summaries
    melt_gebv <- data.table::melt(
      gebv_matrix,
      id.vars       = "Block",
      variable.name = "Sample",
      value.name    = "Effect"
    )

    melt_gebv[, Label := data.table::fcase(
      Effect > 0,  "Positive",
      Effect < 0,  "Negative",
      Effect == 0, "Zero"
    )]

    # Per-sample total GEBV
    gebv_summary <- melt_gebv[, .(GEBV = sum(Effect, na.rm = TRUE)),
                               by = Sample]
    fwrite(gebv_summary,
           file = file.path(output_dir,
             paste0("GEBV_", ld_threshold, "_t", tolerance, "_", trait, ".txt")),
           sep = "\t")

    # Block counts
    block_counts <- melt_gebv[!is.na(Effect),
                               .(Count = .N), by = .(Sample, Label)]
    block_counts <- merge(block_counts, gebv_summary, by = "Sample", all.x = TRUE)

    # Order samples by positive block count
    pos_order <- block_counts[Label == "Positive"][order(Count), Sample]
    pos_order <- as.character(pos_order)
    for (s in unique(block_counts$Sample)) {
      if (!s %in% pos_order) pos_order <- c(s, pos_order)
    }
    block_counts$Sample <- factor(block_counts$Sample, levels = pos_order)
    block_counts$Label  <- factor(block_counts$Label,
                                   levels = c("Positive","Zero","Negative"))

    fwrite(block_counts,
           file = file.path(output_dir,
             paste0("Block_Counts_", ld_threshold, "_t", tolerance,
                    "_", trait, ".txt")),
           sep = "\t")

    # ── Plots ─────────────────────────────────────────────────────────────────
    p_counts <- plot_block_counts(block_counts, trait_name = trait)
    ggplot2::ggsave(
      file.path(output_dir, paste0("Block_Counts_", trait, ".png")),
      p_counts, width = 12, height = 6, dpi = 150
    )

    # message plotting progress and trait completion
    message(
      "  Trait ", trait, " completed. GEBV matrix, summary, and block counts saved."
    )

    # ── Block-haplotype catalogue (long, origin-compatible 12-column TSV) ──────
    hap_catalogue <- compute_block_catalogue(
      geno_meta      = geno_meta,
      effect_df      = effect_df,
      block_info     = bloc_info_dt,
      center         = center,
      .core          = compute_core,
      .signatures    = compute_sigs
    )
    fwrite(hap_catalogue,
           file = file.path(output_dir,
             paste0("Haplotype_Catalogue_", ld_threshold, "_t", tolerance,
                    "_", trait, ".txt")),
           sep = "\t")

    # ── Per-block effect summary (across samples) ─────────────────────────────
    # One row per block: the spread (BlockVariance), extremes (BlockMax /
    # BlockMin) and centre (BlockMean) of its local-GEBV effect across samples,
    # with block location. Written as a table and used for the variance scatter.
    block_summary <- melt_gebv[!is.na(Effect), .(
      BlockVariance = stats::var(Effect, na.rm = TRUE),
      BlockMax      = max(Effect, na.rm = TRUE),
      BlockMin      = min(Effect, na.rm = TRUE),
      BlockMean     = mean(Effect, na.rm = TRUE)
    ), by = Block]
    block_summary <- merge(block_summary,
                           bloc_info_dt[, .(Block, Chrom, Pos)],
                           by = "Block", all.x = TRUE)
    .lead_identifier_cols(block_summary)   # Chrom, Pos, Block lead the summary
    fwrite(block_summary,
           file = file.path(output_dir,
             paste0("Block_Summary_", ld_threshold, "_t", tolerance,
                    "_", trait, ".txt")),
           sep = "\t")

    # Variance scatter (genomic scatter of BlockVariance, one plot per trait)
    p_var <- plot_manhattan_variance(block_summary, trait_name = trait)
    ggplot2::ggsave(
      file.path(output_dir,
        paste0("Manhattan_Variance_", ld_threshold, "_t", tolerance,
               "_", trait, ".png")),
      p_var, width = 12, height = 5, dpi = 150
    )

    # Manhattan — effect
    effect_plot_data <- merge(melt_gebv,
                               bloc_info_dt[, .(Block, Chrom, Pos)],
                               by = "Block", all.x = TRUE)
    p_eff <- plot_manhattan_effect(effect_plot_data, trait_name = trait)
    ggplot2::ggsave(
      file.path(output_dir,
        paste0("Manhattan_Effect_", ld_threshold, "_t", tolerance,
               "_", trait, ".png")),
      p_eff, width = 12, height = 5, dpi = 150
    )

    message("  Trait ", trait, " completed.")
    list(
      gebv_matrix   = gebv_matrix,
      gebv_summary  = gebv_summary,
      block_counts  = block_counts,
      hap_catalogue = hap_catalogue,
      block_summary = block_summary,
      geno_meta     = geno_meta
    )
  }

  # ── Parallel dispatch ──────────────────────────────────────────────────────
  use_parallel <- n_cores > 1L &&
    requireNamespace("future",       quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)

  if (n_cores > 1L && !use_parallel)
    warning("'future' and 'future.apply' required for n_cores > 1. ",
            "Falling back to sequential execution.")

  if (use_parallel) {
    oplan <- future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(oplan), add = TRUE)
    results_list <- future.apply::future_lapply(trait_names, .run_one_trait,
                                                future.seed = TRUE)
  } else {
    results_list <- lapply(trait_names, .run_one_trait)
  }

  names(results_list) <- trait_names
  invisible(results_list)
}
