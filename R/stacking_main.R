# =============================================================================
# stacking_main.R
#
# run_stacking_simulation() — the parameterised recurrent-selection driver for
# the haplotype stacking breeding simulation (PRD #013, issue 03 / Phase 4).
#
# One call replaces the origin `init_simulation.R` config script. The driver
# runs the route
#
#     select -> design -> cross -> advance
#
# over replications x schemes x cycles, entirely GenoMeta-native. Every cycle:
#   1. compute_profile_localgebv(level = "haplotype")  (the shared scheme input)
#   2. scheme dispatch (TS / OHS / OCS / GBS) on the haplotype-level local GEBV
#   3. make_cross_design() (all selected-parent pairs x num_offspring)
#   4. .cross_backend_dispatch() ("block" default, or "genomicSimulation")
#   5. the offspring GenoMeta becomes the next cycle's population.
#
# `scheme` accepts a vector -> the schemes are compared in one run. Decision
# cycles run the chosen scheme; intervening `selfing_cycles` generations advance
# under truncation (GS), and `ga_start_by` delays the chosen scheme until that
# cycle (truncation before it). The whole timeline therefore spans
# `select_cycles * (selfing_cycles + 1)` generations.
#
# Output is origin-compatible (see FR5): `plans`, `summary`, `ideal_genotypes`,
# `final_populations`, plus a `trajectories` compat view for the existing
# plot/monitor helpers.
# =============================================================================


# ---------------------------------------------------------------------------
# Internal — per-sample total local GEBV from a haplotype-level table
# ---------------------------------------------------------------------------

# Collapse the wide "Block" + "<Sample>.<idx>" table to genotype level and sum
# over blocks -> a named numeric vector of each sample's total local GEBV.
.stk_total_gebv <- function(localgebv) {
  geno <- .collapse_hap_to_geno(localgebv)
  cols <- setdiff(names(geno), "Block")
  colSums(as.matrix(geno[, cols, with = FALSE]), na.rm = TRUE)
}


# ---------------------------------------------------------------------------
# Internal — ideal haplotype local GEBV (theoretical best genotype ceiling)
# ---------------------------------------------------------------------------

# The theoretical best genotype assembled from the material: at each block take
# the top-`ploidy` best haplotype segment values present in the pool and sum them,
# then sum over blocks. A cheap, D2Ideal-free ceiling on a wide "Block" +
# "<Sample>.<idx>" haplotype-level GEBV table. Measurement only.
.stk_ideal_hap_gebv <- function(localgebv, ploidy = 2L) {
  dt   <- as.data.table(localgebv)
  cols <- setdiff(names(dt), "Block")
  if (!length(cols) || nrow(dt) == 0L) return(NA_real_)
  m <- as.matrix(dt[, cols, with = FALSE])          # [blocks x copies]
  per_block <- apply(m, 1L, function(r) {
    r <- r[!is.na(r)]
    if (!length(r)) return(NA_real_)
    sum(sort(r, decreasing = TRUE)[seq_len(min(ploidy, length(r)))])
  })
  sum(per_block, na.rm = TRUE)
}


# ---------------------------------------------------------------------------
# Internal — block restriction for selection (selected_blocks / select_top_blocks)
# ---------------------------------------------------------------------------

# Restrict the haplotype-level GEBV table to the blocks selection should focus
# on: an explicit `selected_blocks` set and/or the top `top_n` blocks by summed
# |Effect|. Crossing still transmits every block; this only narrows what the
# selection objective sees.
.stk_limit_blocks <- function(localgebv, selected_blocks = NULL, top_n = NULL) {
  dt <- as.data.table(localgebv)
  if (!is.null(selected_blocks))
    dt <- dt[Block %in% selected_blocks]
  if (!is.null(top_n) && nrow(dt) > top_n) {
    sample_cols <- setdiff(names(dt), "Block")
    score <- rowSums(abs(as.matrix(dt[, sample_cols, with = FALSE])), na.rm = TRUE)
    dt <- dt[order(score, decreasing = TRUE)][seq_len(top_n)]
  }
  dt[]
}


# ---------------------------------------------------------------------------
# Internal — truncation (TS / GS) selection on total local GEBV
# ---------------------------------------------------------------------------

.stk_select_truncation <- function(localgebv, n_parents, select_priority = "high") {
  v   <- .stk_total_gebv(localgebv)
  dec <- select_priority == "high"
  names(sort(v, decreasing = dec))[seq_len(min(n_parents, length(v)))]
}


# ---------------------------------------------------------------------------
# Internal — scheme dispatch
# ---------------------------------------------------------------------------

# Select `n_parents` from the current pool under one scheme. `opts` carries the
# scheme-specific knobs (ohs_method, gbs_method, select_priority, max_kinship,
# ...). Always returns >= 2 distinct parents (falls back to truncation when a
# scheme yields too few, e.g. OCS stopping early on the kinship cap).
.stk_select_parents <- function(scheme, localgebv, geno_meta, n_parents,
                                ideal_genotypes, opts) {
  sel <- switch(
    scheme,
    TS  = .stk_select_truncation(localgebv, n_parents, opts$select_priority),
    OHS = ohs_selection(localgebv, n_parents = n_parents,
                        method = opts$ohs_method, seed = opts$seed),
    OCS = {
      v   <- .stk_total_gebv(localgebv)
      grm <- tryCatch(suppressMessages(compute_grm(geno_meta)),
                      error = function(e) NULL)
      ocs_selection(v, kinship_matrix = grm, n_parents = n_parents,
                    max_kinship = opts$max_kinship,
                    freq_penalty_weight = opts$freq_penalty_weight,
                    seed = opts$seed)
    },
    GBS = gbs_selection(localgebv, ideal_genotypes = ideal_genotypes,
                        n_parents = n_parents, gbs_method = opts$gbs_method,
                        gbs_objective = opts$gbs_objective,
                        gbs_engine = opts$gbs_engine,
                        penalty_mode = opts$penalty_mode,
                        gbs_penalty = opts$gbs_penalty,
                        penalty_scale = opts$penalty_scale,
                        gbs_effect_weight = opts$gbs_effect_weight,
                        select_priority = opts$select_priority,
                        seed = opts$seed),
    stop(".stk_select_parents: unknown scheme '", scheme, "'.")
  )

  sel <- unique(sel)
  if (length(sel) < 2L)
    sel <- .stk_select_truncation(localgebv, max(n_parents, 2L),
                                  opts$select_priority)
  sel
}


# ---------------------------------------------------------------------------
# Internal — coancestry among the selected parents
# ---------------------------------------------------------------------------

# Mean off-diagonal (non-self) block-IBS coancestry among the selected parents
# (FR5 `init_coancestry`), computed by .stk_block_ibs_coancestry() in
# stacking_monitor.R. This is the genuine block-level relatedness measure (no
# pedigree required), superseding the Phase-4 VanRaden-GRM stand-in. `blocks`
# scopes it to the blocks selection focuses on, matching what the schemes see.
.stk_coancestry <- function(geno_meta, parents, blocks = NULL) {
  .stk_block_ibs_coancestry(geno_meta, parents, blocks = blocks)
}


# ---------------------------------------------------------------------------
# Internal — compact generational naming
# ---------------------------------------------------------------------------

# Rename a generation of offspring with compact positional IDs so identifiers
# stop encoding the full crossing history (which otherwise compounds without
# bound across cycles). The lineage is preserved in the Parent1/Parent2 columns
# + the cycle, not in the string.
#
#   sample name : C<cycle>F<family>_<child>   (e.g. C03F045_2)
#   cross id    : C<cycle>F<family>           (e.g. C03F045)  one pair = one family
#
# `cross_res` is a backend result list(geno_meta, plans) where plans has one row
# per offspring (Parent1, Parent2, CrossID, F1mean) in the GenoMeta's sample
# order. Returns list(off = renamed GenoMeta, bp = per-offspring data.table with
# the new ids + FamID, cross_plans = per-cross data.table CrossID/Parent1/
# Parent2/F1count/F1mean).
.stk_rename_generation <- function(cross_res, cycle, opts) {
  bp <- as.data.table(cross_res$plans)
  bp[, OldID := as.character(CrossID)]
  # Families in first-appearance order == cross-design order (one pair = one
  # family; all_pairs yields each pair once).
  bp[, FamID := .GRP, by = .(Parent1, Parent2)]
  bp[, Child := seq_len(.N), by = FamID]
  bp[, NewID := sprintf(opts$id_fmt_sample, cycle, FamID, Child)]

  src       <- cross_res$geno_meta
  arr       <- src$geno_arrays
  old_names <- dimnames(arr)[[1L]]
  idx       <- match(old_names, bp$OldID)
  new_names <- ifelse(is.na(idx), old_names, bp$NewID[idx])  # keep unmatched as-is
  dimnames(arr)[[1L]] <- new_names

  # Offspring FAM: each progeny's real lineage (Family = FamID, its Parent1/
  # Parent2), aligned to the genotype sample order — so the returned GenoMeta is a
  # true progeny population for BOTH backends (the block/gS backends otherwise
  # leave a dummy founder FAM). Unmatched samples fall back to dummy "0".
  mi  <- match(new_names, bp$NewID)
  fam <- data.frame(
    Family  = ifelse(is.na(mi), "0", as.character(bp$FamID[mi])),
    Sample  = new_names,
    Parent1 = ifelse(is.na(mi), "0", as.character(bp$Parent1[mi])),
    Parent2 = ifelse(is.na(mi), "0", as.character(bp$Parent2[mi])),
    Sex     = "0",
    stringsAsFactors = FALSE)

  # Rebuild so the canonical sample_info regenerates against the new names (a
  # bare dimnames swap would leave sample_info$Sample stale and fail validation).
  md <- src$metadata
  md$n_samples <- length(new_names)
  off <- new_geno_meta(
    geno_arrays           = arr,
    map_info              = src$map_info,
    block_info            = src$block_info,
    block_marker_map_info = src$block_marker_map_info,
    sample_info           = fam,
    metadata              = md,
    marker_effects        = src$marker_effects)

  cross_plans <- bp[, .(F1count = .N, F1mean = mean(F1mean, na.rm = TRUE)),
                    by = .(FamID, Parent1, Parent2)]
  cross_plans[, CrossID := sprintf(opts$id_fmt_cross, cycle, FamID)]

  bp[, CrossID := NewID]
  list(off = off, bp = bp, cross_plans = cross_plans)
}


# ---------------------------------------------------------------------------
# Internal — advance one generation under a scheme
# ---------------------------------------------------------------------------

# select -> design -> cross -> rename -> score. Returns the (compactly renamed)
# offspring GenoMeta, the per-cross plan, the per-offspring F1mean vector and the
# cycle metrics (D2Ideal, coancestry, the ideal-genotype target and its
# ideal_phenotype, marker count). Under `opts$ideal_scope = "founder"` (default)
# the ideal genotype is `opts$founder_ideal`, built once from the founders
# (Kemper 2012: targets fixed for the whole run, a lost state stays a lost
# target); under `"cycle"` it is recomputed here from the current pool (0.6.0 /
# origin behaviour: a lost state drops out and the remaining targets are
# re-normalised). `cycle` is the whole-cycle index used to build the compact
# offspring IDs (see .stk_rename_generation).
.stk_advance <- function(current, scheme, effect_df, opts, cycle, session = NULL) {
  want <- function(m) m %in% opts$metrics

  # Per-step wall-clock timing (surfaced in verbose mode so the user can see which
  # step — selection, crossing, offspring GEBV, metrics — dominates the run).
  timings <- c(parent_gebv = 0, ideal = 0, selection = 0, crossing = 0,
               offspring_gebv = 0, metrics = 0)
  .tk <- proc.time()[["elapsed"]]
  lap <- function(step) {
    now <- proc.time()[["elapsed"]]
    timings[[step]] <<- timings[[step]] + (now - .tk)
    .tk <<- now
  }

  lg <- compute_profile_localgebv(current, effect_df, level = "haplotype")
  # A founder-fixed ideal defines the block set of GBS selection for the whole
  # run (the founders' `select_top_blocks`), so a target block cannot silently
  # drop out of a later cycle's top-N. Only GBS selection is pinned: TS/OHS/OCS
  # keep the per-cycle top-N, and the d2ideal metric is measured on the full
  # offspring table, so neither is affected by the scope.
  use_founder <- identical(opts$ideal_scope, "founder") &&
    !is.null(opts$founder_ideal)
  lg <- if (use_founder && scheme == "GBS")
    .stk_limit_blocks(lg, opts$founder_blocks, NULL)
  else
    .stk_limit_blocks(lg, opts$selected_blocks, opts$select_top_blocks)
  if (nrow(lg) == 0L)
    stop(".stk_advance: no blocks remain for selection (check 'selected_blocks' ",
         "/ block assignments).")
  lap("parent_gebv")

  # Ideal genotype has two INDEPENDENT roles: (a) the SELECTION target for GBS
  # (built whenever GBS runs, never influenced by 'metrics'), and (b) the input to
  # the d2ideal MEASUREMENT. Built once and shared when both apply, so metrics can
  # never change GBS selection.
  need_ideal <- scheme == "GBS" || want("d2ideal")
  ideal <- if (!need_ideal) NULL
  else if (use_founder) opts$founder_ideal
  else construct_ideal_genotype(lg, method = opts$construct_method,
                                k_penalty = opts$k_penalty,
                                select_priority = opts$select_priority)
  lap("ideal")

  parents <- .stk_select_parents(scheme, lg, current, opts$top_founders,
                                 ideal, opts)
  lap("selection")

  cross_design <- make_cross_design(parents, method = "all_pairs",
                                    avoid_selfing = opts$avoid_selfing,
                                    n_offspring = opts$num_offspring)

  res <- .cross_backend_dispatch(opts$cross_backend, current, cross_design,
                                 effect_df = effect_df, work_dir = opts$work_dir,
                                 session = session)

  # Compact, history-free IDs for the new generation (stops the cross-id /
  # sample-name blow-up; lineage stays in Parent1/Parent2 + cycle).
  gen <- .stk_rename_generation(res, cycle, opts)
  off <- gen$off

  # FR3 (persistent gS session): push the SAME compact names into the live group
  # so next cycle's parent-name -> "X" lookup in gs_cross_plan() resolves. Ordering
  # is load-bearing: cross -> rename -> sync -> next cycle. No-op for the block
  # backend (session is NULL).
  session <- res$session
  if (!is.null(session))
    session <- gs_sync_group_names(session, gen$bp$OldID, gen$bp$NewID)
  lap("crossing")

  # ── Optional metrics — all measured on the OFFSPRING population ──────────────
  # Offspring local GEBV recompute is the expensive extra pass; do it only when a
  # metric needs it. Everything below fills NA unless its token is in opts$metrics.
  need_offlg <- want("d2ideal") || want("ideal_hap_gebv")
  off_lg <- if (need_offlg)
    tryCatch(compute_profile_localgebv(off, effect_df, level = "haplotype"),
             error = function(e) NULL)
  else NULL
  lap("offspring_gebv")

  d2 <- if (want("d2ideal") && !is.null(ideal) &&
            !is.null(off_lg) && nrow(off_lg) > 0L)
    tryCatch(.distance_to_ideal(off_lg, ideal, gbs_method = opts$gbs_method),
             error = function(e) NA_real_)
  else NA_real_

  ideal_hap <- if (want("ideal_hap_gebv") && !is.null(off_lg) && nrow(off_lg) > 0L)
    tryCatch(.stk_ideal_hap_gebv(off_lg, ploidy = opts$ploidy),
             error = function(e) NA_real_)
  else NA_real_

  # Nei's (1973) gene diversity of the offspring population (marker-level).
  gdiv <- if (want("gene_diversity"))
    tryCatch(as.numeric(compute_gene_diversity(off, level = "marker")),
             error = function(e) NA_real_)
  else NA_real_

  # Coancestry = mean block-IBS relatedness among the OFFSPRING (a diversity
  # measure of the produced population, distinct from Nei gene diversity).
  coan <- if (want("coancestry"))
    tryCatch(.stk_coancestry(off, dimnames(off$geno_arrays)[[1L]],
                             blocks = unique(as.character(lg$Block))),
             error = function(e) NA_real_)
  else NA_real_
  lap("metrics")

  list(off             = off,
       cross_plans     = gen$cross_plans,
       f1              = as.numeric(gen$bp$F1mean),
       ideal_phenotype = if (is.null(ideal)) NA_real_
                         else as.numeric(attr(ideal, "ideal_phenotype")),
       d2              = d2,
       ideal_hap       = ideal_hap,
       coan            = coan,
       gene_div        = gdiv,
       parents         = parents,
       n_marker        = dim(off$geno_arrays)[2L],
       session         = session,
       timings         = timings)
}


# ---------------------------------------------------------------------------
# Internal — verbose per-step timing breakdown
# ---------------------------------------------------------------------------

# Labels for the timed steps of .stk_advance() (shared by the progress line and
# the end-of-job breakdown).
.STK_STEP_LABELS <- c(
  parent_gebv    = "parent local GEBV",
  ideal          = "ideal genotype",
  selection      = "parent selection",
  crossing       = "crossing (simulation)",
  offspring_gebv = "offspring local GEBV",
  metrics        = "offspring metrics")

# Human-readable per-step timing breakdown for a (scheme, rep) job: steps sorted
# by elapsed time with percentages, dominant step called out. Verbose only.
.stk_timing_summary <- function(scheme, rep, step_time) {
  tot <- sum(step_time)
  if (!is.finite(tot) || tot <= 0) return(NULL)
  st    <- sort(step_time, decreasing = TRUE)
  lines <- sprintf("      %-22s %8.2fs  (%4.1f%%)",
                   .STK_STEP_LABELS[names(st)], st, 100 * st / tot)
  paste0(sprintf("  [%s rep %d] step time — total %.1fs; dominant: %s\n",
                 scheme, rep, tot,
                 .STK_STEP_LABELS[[names(which.max(step_time))]]),
         paste(lines, collapse = "\n"))
}


# ---------------------------------------------------------------------------
# Internal — origin-shaped per-cycle records
# ---------------------------------------------------------------------------

# Per-cross `plans` rows (FR5 column order): the pre-built per-cross plan
# (compact CrossID + parents) augmented with the cycle's identifiers and metrics.
.stk_plans_record <- function(adv, select_type, cycle, rep, sel_scheme) {
  cp <- data.table::copy(adv$cross_plans)
  cp[, FamID := NULL]
  cp[, `:=`(select_type     = select_type,
            cycle           = cycle,
            rep             = rep,
            sel_scheme      = sel_scheme,
            D2Ideal         = adv$d2,
            coancestry      = adv$coan,
            ideal_phenotype = adv$ideal_phenotype)]
  setcolorder(cp, c("CrossID", "Parent1", "Parent2", "F1count", "F1mean",
                    "select_type", "cycle", "rep", "sel_scheme",
                    "D2Ideal", "coancestry", "ideal_phenotype"))
  cp[]
}

# One `summary` row (FR5 columns) for a cycle.
.stk_summary_record <- function(adv, select_type, cycle, rep, sel_scheme) {
  f1 <- adv$f1
  data.table(cycle                = cycle,
             select_type          = select_type,
             rep                  = rep,
             sel_scheme           = sel_scheme,
             mean_F1mean          = mean(f1, na.rm = TRUE),
             var_F1mean           = stats::var(f1, na.rm = TRUE),
             mean_D2Ideal         = adv$d2,
             mean_coancestry      = adv$coan,
             mean_gene_diversity  = adv$gene_div,
             mean_ideal_phenotype = adv$ideal_phenotype,
             ideal_hap_gebv       = adv$ideal_hap,
             marker_size          = adv$n_marker)
}


# ---------------------------------------------------------------------------
# run_stacking_simulation()
# ---------------------------------------------------------------------------

#' Run the haplotype stacking breeding simulation
#'
#' Parameterised recurrent-selection driver: from a founder \code{GenoMeta} and
#' a marker-effect table it runs multi-cycle breeding (select -> design -> cross
#' -> advance) across replications and one or more selection schemes, scoring
#' each cycle on haplotype-level local GEBV. One call replaces the origin
#' \code{init_simulation.R} config script.
#'
#' @section Selection schemes:
#' \code{scheme} chooses one or more of (a vector compares them in one run):
#' \describe{
#'   \item{\code{"TS"}}{Truncation (genomic) selection: top \code{top_founders}
#'     by total local GEBV.}
#'   \item{\code{"OHS"}}{Optimal Haplotype Selection — maximise next-cross gain
#'     via \code{\link{ohs_selection}} (\code{ohs_method}: \code{"GA_pairwise"},
#'     \code{"GA_coverage"}, \code{"Greedy"}).}
#'   \item{\code{"OCS"}}{Optimal Contribution Selection — kinship-constrained
#'     (\code{\link{ocs_selection}} on a VanRaden GRM, \code{max_kinship}).}
#'   \item{\code{"GBS"}}{Genotype Building Strategy — drive block-haplotype
#'     frequencies toward an ideal genotype via \code{\link{gbs_selection}}
#'     (\code{gbs_method}: \code{"ARCSIN"}/\code{"LOG"}).}
#' }
#' \code{ga_start_by} delays the chosen scheme until that cycle (truncation
#' before it); \code{selfing_cycles} interleaves that many truncation-advanced
#' selfing generations between decision cycles.
#'
#' @section Crossing backends:
#' \describe{
#'   \item{\code{"block"} (default)}{Pure-R, dependency-free. The haplotype
#'     blocks are the crossing loci and segregate as indivisible units: each
#'     offspring inherits whole parental blocks with \emph{no intra-block
#'     recombination}, the natural model for the upper bound of stacking gain.}
#'   \item{\code{"genomicSimulation"}}{Opt-in (\pkg{genomicSimulation}, in
#'     Suggests). Marker-level crossing with real linkage and recombination.}
#' }
#' Both backends are interchangeable (internal \code{.cross_backend_dispatch()}).
#'
#' Replications x schemes are embarrassingly parallel: set \code{n_cores > 1} to
#' use \pkg{future.apply}. Each replication uses a deterministic seed
#' (\code{seed + rep}) so results reproduce regardless of parallelism.
#'
#' @param geno_meta Founder \code{GenoMeta}. For back-compatibility a
#'   \code{catalogue_results} list (see \code{\link{run_haplotype_catalogue}})
#'   may be passed here or via \code{catalogue_results}.
#' @param effect_df data.frame of marker effects (Marker, Allele, Effect). When
#'   \code{NULL}, the \code{GenoMeta}'s attached \code{marker_effects} are used;
#'   if those hold a named per-trait list, the first trait's table is taken (this
#'   pipeline is single-trait).
#' @param selected_blocks Optional character vector restricting selection to
#'   these blocks. Default \code{NULL} (all blocks).
#' @param replications Integer. Independent simulation runs. Default \code{10}.
#' @param select_cycles Integer. Decision cycles per replication. Default
#'   \code{100}.
#' @param selfing_cycles Integer. Truncation-advanced selfing generations
#'   between decision cycles. Default \code{0}.
#' @param ga_start_by Integer. Cycle from which the chosen scheme is applied
#'   (truncation before it). Default \code{1}.
#' @param num_offspring Integer. Offspring per cross. Default \code{4}.
#' @param top_founders Integer. Parents selected per cycle. Default \code{16}.
#' @param avoid_selfing Logical. Exclude self-crosses from the design. Default
#'   \code{FALSE}.
#' @param scheme Character vector of schemes to run/compare: any of \code{"TS"},
#'   \code{"OHS"}, \code{"OCS"}, \code{"GBS"}.
#' @param ohs_method OHS engine (see \code{\link{ohs_selection}}).
#' @param gbs_method GBS deviation form, \code{"ARCSIN"} (default) or
#'   \code{"LOG"} — Kemper's two modes (see \code{\link{gbs_selection}}).
#' @param gbs_preset \code{"haploStack"} (default) or \code{"kemper"} (the
#'   published method: fixed \code{-100} penalty). Fills \code{gbs_objective},
#'   \code{gbs_engine}, \code{penalty_mode}, \code{gbs_penalty},
#'   \code{penalty_scale} and \code{gbs_effect_weight} when they are \code{NULL};
#'   an explicit value always wins. See \code{\link{gbs_selection}}.
#' @param gbs_objective GBS objective, \code{"candidate"} or \code{"set"};
#'   \code{NULL} (default) takes the preset value (\code{"candidate"}).
#' @param gbs_engine GBS engine, \code{"greedy"} or \code{"GA"}; \code{NULL}
#'   (default) takes the preset value (\code{"greedy"}).
#' @param penalty_mode GBS absent-state penalty scaling, \code{"dynamic"} or
#'   \code{"fixed"}; \code{NULL} (default) takes the preset value
#'   (\code{"dynamic"} for \code{"haploStack"}, \code{"fixed"} for \code{"kemper"}).
#' @param gbs_penalty GBS base penalty (\code{<= 0}); \code{NULL} (default) takes
#'   the preset value \code{-10}.
#' @param penalty_scale GBS penalty scale or cap (\code{>= 1}); \code{NULL}
#'   (default) takes the preset value \code{10}.
#' @param gbs_effect_weight GBS effect-magnitude weight; \code{NULL} (default)
#'   takes the preset value \code{0}.
#' @param ideal_scope \code{"founder"} (default) builds the ideal genotype once
#'   from the founders and keeps it for every cycle — Kemper's design, and the
#'   scope under which the \code{D2Ideal} trajectory is comparable across cycles;
#'   GBS selection is then also pinned to the founders' block set.
#'   \code{"cycle"} rebuilds it from the current pool every cycle (the 0.6.0 /
#'   origin behaviour: a lost state drops out and the remaining targets are
#'   re-normalised).
#' @param select_priority \code{"high"} (default) or \code{"low"} — selection
#'   direction (shared by TS/OHS-ideal/GBS and the ideal genotype).
#' @param construct_method Ideal-genotype construction (see
#'   \code{\link{construct_ideal_genotype}}).
#' @param k_penalty Diversity penalty for \code{construct_method = "custom"}.
#'   Default \code{0.5}.
#' @param max_kinship OCS mean-kinship cap. Default \code{0.125}.
#' @param freq_penalty_weight OCS frequency-deviation weight. Default \code{1}.
#' @param select_top_blocks Integer. Cap on blocks selection considers (top by
#'   summed \code{|Effect|}). Default \code{200}.
#' @param cross_backend \code{"block"} (default) or \code{"genomicSimulation"}.
#' @param metrics Character vector of per-cycle metrics to compute on the
#'   \strong{offspring} population. Default \code{"gebv"} (mean/var of F1mean only).
#'   Append any of \code{"d2ideal"} (distance to the ideal genotype),
#'   \code{"gene_diversity"} (Nei marker-level), \code{"coancestry"} (mean block-IBS
#'   relatedness among offspring), \code{"ideal_hap_gebv"} (theoretical best genotype
#'   = per-block top-\code{ploidy} segment sum), or \code{"all"} for every metric.
#'   \strong{Metrics are measurement-only: they never affect parental selection or
#'   the offspring produced} — for a fixed seed the populations are identical for any
#'   \code{metrics}. Unrequested metric columns are present but \code{NA}. (GBS still
#'   builds its ideal-genotype selection target regardless of this argument.)
#' @param seed Base random seed (replication \code{r} uses \code{seed + r}).
#'   Default \code{42}.
#' @param output_dir Directory for output files/plots. Created if absent.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#' @param monitor Logical. Live monitor: when \code{TRUE} \strong{and}
#'   \code{output_dir} is set, the accumulating record and figures are refreshed
#'   after \emph{each} \code{(rep, scheme)} job completes (rather than only at the
#'   end), so a long run shows progress mid-flight — \code{live_summary.txt},
#'   \code{live_plans.txt}, \code{live_metrics.png}
#'   (\code{\link{plot_simulation_metrics}}) and \code{live_trajectories.png}
#'   (\code{\link{plot_simulation_trajectories}}) are overwritten in
#'   \code{output_dir}. Works in both sequential and parallel runs (jobs are
#'   harvested as they resolve). \strong{Observation-only}: for a fixed seed the
#'   returned \code{plans}/\code{summary}/\code{final_populations} are identical to
#'   a \code{monitor = FALSE} run. Default \code{FALSE}.
#' @param n_cores Integer. Parallel workers (\code{> 1} needs \pkg{future} and
#'   \pkg{future.apply}). Default \code{1}.
#' @param catalogue_results,trait_name,select_method Deprecated back-compat
#'   inputs. \code{catalogue_results}/\code{trait_name} supply the founder
#'   \code{GenoMeta}; \code{select_method} (\code{"GA"}/\code{"GS"}/\code{"OCS"})
#'   maps to \code{scheme} (\code{GA -> OHS(GA_coverage)}, \code{GS -> TS},
#'   \code{OCS -> OCS}).
#'
#' @return A named list:
#' \describe{
#'   \item{\code{plans}}{Per-cross data.table: \code{CrossID, Parent1, Parent2,
#'     F1count, F1mean, select_type, cycle, rep, sel_scheme, D2Ideal,
#'     coancestry, ideal_phenotype}.}
#'   \item{\code{summary}}{Per-cycle data.table: \code{cycle, select_type, rep,
#'     sel_scheme, mean_F1mean, var_F1mean, mean_D2Ideal, mean_coancestry,
#'     mean_gene_diversity, mean_ideal_phenotype, marker_size} (with a
#'     \code{cycle = 0} founder row). \code{mean_gene_diversity} is Nei's (1973)
#'     gene diversity of the cycle's population (see
#'     \code{\link{compute_gene_diversity}}).}
#'   \item{\code{ideal_genotypes}}{The founder ideal-genotype target:
#'     \code{Block, Effect, CFreq, TargetFreq}.}
#'   \item{\code{final_populations}}{Named list (\code{"<scheme>_rep<r>"} ->
#'     final offspring \code{GenoMeta}).}
#'   \item{\code{trajectories}}{Compat view for
#'     \code{\link{plot_simulation_trajectories}}: \code{Rep, Cycle, Method,
#'     MeanGEBV, SDev, BlockCoverage}.}
#' }
#'
#' @seealso \code{\link{ohs_selection}}, \code{\link{ocs_selection}},
#'   \code{\link{gbs_selection}}, \code{\link{construct_ideal_genotype}},
#'   \code{\link{make_cross_design}}, \code{\link{plot_simulation_trajectories}}
#'
#' @export
run_stacking_simulation <- function(
  geno_meta           = NULL,
  effect_df           = NULL,
  selected_blocks     = NULL,
  replications        = 10L,
  select_cycles       = 100L,
  selfing_cycles      = 0L,
  ga_start_by         = 1L,
  num_offspring       = 4L,
  top_founders        = 16L,
  avoid_selfing       = FALSE,
  scheme              = c("TS", "OHS", "OCS", "GBS"),
  ohs_method          = c("GA_pairwise", "GA_coverage", "Greedy"),
  gbs_method          = c("ARCSIN", "LOG"),
  gbs_preset          = c("haploStack", "kemper"),
  gbs_objective       = NULL,
  gbs_engine          = NULL,
  penalty_mode        = NULL,
  gbs_penalty         = NULL,
  penalty_scale       = NULL,
  gbs_effect_weight   = NULL,
  ideal_scope         = c("founder", "cycle"),
  select_priority     = c("high", "low"),
  construct_method    = c("custom", "argmax", "ratio"),
  k_penalty           = 0.5,
  max_kinship         = 0.125,
  freq_penalty_weight = 1,
  select_top_blocks   = 200L,
  cross_backend       = c("block", "genomicSimulation"),
  metrics             = "gebv",
  seed                = 42L,
  output_dir          = "stacking_output",
  verbose             = TRUE,
  monitor             = FALSE,
  n_cores             = 1L,
  # ── deprecated back-compat ───────────────────────────────────────────────
  catalogue_results   = NULL,
  trait_name          = NULL,
  select_method       = NULL
) {

  # ── Back-compat: catalogue_results / select_method ────────────────────────
  if (is.null(catalogue_results) && !is.null(geno_meta) && !is_geno_meta(geno_meta))
    catalogue_results <- geno_meta            # catalogue list passed positionally

  if (!is.null(catalogue_results)) {
    trait_name <- trait_name %||% names(catalogue_results)[1L]
    if (!trait_name %in% names(catalogue_results))
      stop("run_stacking_simulation: trait '", trait_name,
           "' not in catalogue_results.")
    geno_meta <- catalogue_results[[trait_name]]$geno_meta
  }

  if (!is.null(select_method)) {              # GA/GS/OCS -> scheme
    select_method <- match.arg(select_method, c("GA", "GS", "OCS"))
    scheme <- switch(select_method, GA = "OHS", GS = "TS", OCS = "OCS")
    if (select_method == "GA") ohs_method <- "GA_coverage"
  }

  if (!is_geno_meta(geno_meta))
    stop("run_stacking_simulation: 'geno_meta' must be a valid GenoMeta ",
         "(or supply 'catalogue_results').")

  scheme           <- match.arg(scheme, c("TS", "OHS", "OCS", "GBS"),
                                several.ok = TRUE)
  ohs_method       <- match.arg(ohs_method)
  gbs_method       <- match.arg(gbs_method)
  gbs_preset       <- match.arg(gbs_preset)
  ideal_scope      <- match.arg(ideal_scope)
  select_priority  <- match.arg(select_priority)
  construct_method <- match.arg(construct_method)
  cross_backend    <- match.arg(cross_backend)

  # GBS knobs: the preset fills whatever was left NULL; explicit values win.
  gbs_spec <- .gbs_apply_preset(gbs_preset, list(
    gbs_objective = if (is.null(gbs_objective)) NULL
                    else match.arg(gbs_objective, c("candidate", "set")),
    gbs_engine    = if (is.null(gbs_engine)) NULL
                    else match.arg(gbs_engine, c("greedy", "GA")),
    penalty_mode  = if (is.null(penalty_mode)) NULL
                    else match.arg(penalty_mode, c("dynamic", "fixed")),
    gbs_penalty       = gbs_penalty,
    penalty_scale     = penalty_scale,
    gbs_effect_weight = gbs_effect_weight))
  .gbs_check_penalty_args(gbs_spec$gbs_penalty, gbs_spec$penalty_mode,
                          gbs_spec$penalty_scale)
  if (gbs_spec$gbs_engine == "GA" && gbs_spec$gbs_objective != "set")
    stop("run_stacking_simulation: gbs_engine = \"GA\" requires ",
         "gbs_objective = \"set\".")

  # ── Metrics vector (measurement-only; never affects selection) ─────────────
  # "gebv" (F1mean) is the always-on base; append any optional offspring metrics,
  # or "all" for every one. Every value is measured on the offspring population.
  all_metrics <- c("gebv", "d2ideal", "gene_diversity", "coancestry",
                   "ideal_hap_gebv")
  if (any(metrics == "all")) metrics <- all_metrics
  metrics <- match.arg(metrics, all_metrics, several.ok = TRUE)
  metrics <- union("gebv", metrics)          # F1mean cannot be disabled

  # Single-trait pipeline: collapse a stored per-trait marker_effects list to one
  # effect table (uses the first trait when several are present).
  if (is.null(effect_df)) effect_df <- geno_meta$marker_effects
  effect_df <- .coerce_single_effect_df(effect_df)
  if (is.null(effect_df))
    stop("run_stacking_simulation: no effects supplied — pass 'effect_df' or ",
         "attach 'marker_effects' to the GenoMeta.")

  if (!is.null(output_dir) && !dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  whole_cycles <- select_cycles * (selfing_cycles + 1L)

  # Fixed-width fields for the compact generational IDs (C<cycle>F<family>_<child>):
  # widths are sized to the run so IDs stay short, uniform and sortable.
  fam_max <- max(2L, choose(top_founders, 2L))
  cyc_w   <- max(2L, nchar(as.character(whole_cycles)))
  fam_w   <- max(2L, nchar(as.character(fam_max)))
  ch_w    <- max(1L, nchar(as.character(num_offspring)))

  # Scheme-specific knobs captured for the cycle helpers.
  opts <- list(
    selected_blocks   = selected_blocks,
    select_top_blocks = select_top_blocks,
    top_founders      = top_founders,
    num_offspring     = num_offspring,
    avoid_selfing     = avoid_selfing,
    ohs_method        = ohs_method,
    gbs_method        = gbs_method,
    gbs_objective     = gbs_spec$gbs_objective,
    gbs_engine        = gbs_spec$gbs_engine,
    penalty_mode      = gbs_spec$penalty_mode,
    gbs_penalty       = gbs_spec$gbs_penalty,
    penalty_scale     = gbs_spec$penalty_scale,
    gbs_effect_weight = gbs_spec$gbs_effect_weight,
    ideal_scope       = ideal_scope,
    founder_ideal     = NULL,          # filled below once the founder ideal exists
    founder_blocks    = NULL,
    select_priority   = select_priority,
    construct_method  = construct_method,
    k_penalty         = k_penalty,
    max_kinship       = max_kinship,
    freq_penalty_weight = freq_penalty_weight,
    cross_backend     = cross_backend,
    metrics           = metrics,
    ploidy            = geno_meta$metadata$ploidy %||% 2L,
    seed              = seed,
    # Base dir for the genomicSimulation backend's transient scratch files. Kept
    # under tempdir() (never output_dir, which may be on a synced/slow drive) and
    # made unique per (rep, scheme) job in .run_one() so parallel workers don't
    # clobber each other's gs_geno/map/eff trio. Unused by the block backend.
    work_dir          = file.path(tempdir(), "stacking_gs_scratch"),
    id_fmt_sample     = paste0("C%0", cyc_w, "dF%0", fam_w, "d_%0", ch_w, "d"),
    id_fmt_cross      = paste0("C%0", cyc_w, "dF%0", fam_w, "d")
  )

  # Founder haplotype-level GEBV (baseline; drives the cycle-0 record).
  founder_lg <- .stk_limit_blocks(
    compute_profile_localgebv(geno_meta, effect_df, level = "haplotype"),
    selected_blocks, select_top_blocks)
  # Founder ideal-genotype target returned to the caller: only meaningful when GBS
  # will use it for selection or d2ideal is measured; else NULL (metrics-only).
  ideal_genotypes <- if ("GBS" %in% scheme || "d2ideal" %in% metrics)
    construct_ideal_genotype(founder_lg, method = construct_method,
                             k_penalty = k_penalty,
                             select_priority = select_priority)
  else NULL
  # Founder-fixed scope: every cycle selects and measures against this table,
  # and GBS selection is pinned to the founders' block set (see .stk_advance).
  if (ideal_scope == "founder" && !is.null(ideal_genotypes)) {
    opts$founder_ideal  <- ideal_genotypes
    opts$founder_blocks <- unique(as.character(founder_lg$Block))
  }

  # ── One (rep, scheme) run of the full route ───────────────────────────────
  .run_one <- function(rep, sch) {
    if (verbose) message("Scheme ", sch, " — replication ", rep, " / ", replications)
    set.seed(seed + rep)

    # Isolate this job's genomicSimulation scratch files (no cross-worker clobber).
    opts$work_dir <- file.path(opts$work_dir, paste0(sch, "_rep", rep))
    if (opts$cross_backend == "genomicSimulation")
      dir.create(opts$work_dir, showWarnings = FALSE, recursive = TRUE)

    current     <- geno_meta
    plan_recs   <- vector("list", whole_cycles + 1L)
    summ_recs   <- vector("list", whole_cycles + 1L)

    # Persistent genomicSimulation session for this (rep, scheme) job (FR1): opened
    # lazily on the first gS cross, carried across all cycles + selfing generations,
    # and freed on exit so no live group leaks (even if a cycle errors). NULL and a
    # no-op for the block backend.
    session <- NULL
    on.exit(gs_close(session), add = TRUE)

    # cycle-0 founder record (start_fam_record). Optional metrics measure the
    # founder population; each is filled only when its token is in `metrics`.
    wantm   <- function(m) m %in% metrics
    f_ideal <- if (wantm("d2ideal"))
      construct_ideal_genotype(founder_lg, method = construct_method,
                               k_penalty = k_penalty,
                               select_priority = select_priority)
    else NULL
    v0 <- .stk_total_gebv(founder_lg)
    summ_recs[[1L]] <- data.table(
      cycle = 0L, select_type = "founder", rep = rep, sel_scheme = sch,
      mean_F1mean = mean(v0, na.rm = TRUE), var_F1mean = stats::var(v0, na.rm = TRUE),
      mean_D2Ideal = if (wantm("d2ideal"))
        .distance_to_ideal(founder_lg, f_ideal, gbs_method = gbs_method)
        else NA_real_,
      mean_coancestry = if (wantm("coancestry"))
        tryCatch(.stk_coancestry(geno_meta, dimnames(geno_meta$geno_arrays)[[1L]],
                                 blocks = unique(as.character(founder_lg$Block))),
                 error = function(e) NA_real_)
        else NA_real_,
      mean_gene_diversity = if (wantm("gene_diversity"))
        tryCatch(as.numeric(compute_gene_diversity(geno_meta, level = "marker")),
                 error = function(e) NA_real_)
        else NA_real_,
      mean_ideal_phenotype = if (is.null(f_ideal)) NA_real_
        else as.numeric(attr(f_ideal, "ideal_phenotype")),
      ideal_hap_gebv = if (wantm("ideal_hap_gebv"))
        .stk_ideal_hap_gebv(founder_lg, ploidy = opts$ploidy) else NA_real_,
      marker_size = dim(geno_meta$geno_arrays)[2L])

    slot <- 1L
    wc   <- 0L
    step_time <- c(parent_gebv = 0, ideal = 0, selection = 0, crossing = 0,
                   offspring_gebv = 0, metrics = 0)
    for (cyc in seq_len(select_cycles)) {
      # decision cycle: chosen scheme once ga_start_by reached, else truncation
      sch_this <- if (cyc >= ga_start_by) sch else "TS"
      wc   <- wc + 1L
      slot <- slot + 1L
      adv  <- .stk_advance(current, sch_this, effect_df, opts, wc, session = session)
      step_time <- step_time + adv$timings
      plan_recs[[slot]] <- .stk_plans_record(adv, sch_this, wc, rep, sch)
      summ_recs[[slot]] <- .stk_summary_record(adv, sch_this, wc, rep, sch)
      current <- adv$off
      session <- adv$session

      # intervening selfing generations advance under truncation (GS)
      for (sf in seq_len(selfing_cycles)) {
        wc   <- wc + 1L
        slot <- slot + 1L
        adv  <- .stk_advance(current, "TS", effect_df, opts, wc, session = session)
        step_time <- step_time + adv$timings
        plan_recs[[slot]] <- .stk_plans_record(adv, "selfing", wc, rep, sch)
        summ_recs[[slot]] <- .stk_summary_record(adv, "selfing", wc, rep, sch)
        current <- adv$off
        session <- adv$session
      }

      # Progress line: cycle counter + the step consuming the most time so far, so
      # a long run tells the user WHERE the time is going without waiting for the end.
      if (verbose && cyc %% 10L == 0L) {
        dom <- .STK_STEP_LABELS[[names(which.max(step_time))]]
        message(sprintf("  %s rep %d: cycle %d / %d — %.0fs elapsed, dominant step: %s",
                        sch, rep, cyc, select_cycles, sum(step_time), dom))
      }
    }

    # End-of-job breakdown: which step took the majority of this (scheme, rep) run.
    if (verbose) {
      summ <- .stk_timing_summary(sch, rep, step_time)
      if (!is.null(summ)) message(summ)
    }

    list(plans   = rbindlist(plan_recs),
         summary = rbindlist(summ_recs),
         final   = current,
         key     = paste0(sch, "_rep", rep))
  }

  # ── Dispatch over the (rep x scheme) grid ─────────────────────────────────
  grid <- expand.grid(rep = seq_len(replications), scheme = scheme,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  n_jobs <- nrow(grid)

  use_parallel <- n_cores > 1L &&
    requireNamespace("future",       quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)
  if (n_cores > 1L && !use_parallel)
    warning("'future' and 'future.apply' required for n_cores > 1. ",
            "Falling back to sequential execution.")

  run_k <- function(k) .run_one(grid$rep[k], grid$scheme[k])

  # Live monitor needs an output_dir to write into; announce and disable if absent.
  live_on <- isTRUE(monitor)
  if (live_on && is.null(output_dir)) {
    message("run_stacking_simulation: 'monitor = TRUE' needs 'output_dir'; ",
            "live figures/tables are skipped.")
    live_on <- FALSE
  }
  live_tag <- paste(unique(scheme), collapse = "-")

  job_results <- vector("list", n_jobs)

  if (live_on) {
    # Harvest jobs as they resolve and refresh the record/figures after each one.
    # The accumulated per-job outputs are exactly those used in the final combine
    # below, so the return value is identical to a monitor = FALSE run.
    acc_summ  <- vector("list", n_jobs)
    acc_plans <- vector("list", n_jobs)
    n_done    <- 0L
    collect <- function(k, res) {
      job_results[[k]] <<- res
      acc_summ[[k]]    <<- res$summary
      acc_plans[[k]]   <<- res$plans
      n_done <<- n_done + 1L
      if (verbose)
        message(sprintf("[%d/%d] %s rep %d done - mean_F1mean=%.4g",
                        n_done, n_jobs, grid$scheme[k], grid$rep[k],
                        mean(res$summary[cycle > 0L]$mean_F1mean, na.rm = TRUE)))
      keep <- !vapply(acc_summ, is.null, logical(1L))
      .stk_render_live(data.table::rbindlist(acc_summ[keep]),
                       data.table::rbindlist(acc_plans[keep]),
                       output_dir, live_tag)
    }

    if (use_parallel) {
      oplan <- future::plan(future::multisession, workers = n_cores)
      on.exit(future::plan(oplan), add = TRUE)
      futs    <- lapply(seq_len(n_jobs),
                        function(k) future::future(run_k(k), seed = NULL))
      pending <- seq_len(n_jobs)
      while (length(pending)) {
        ready <- pending[vapply(futs[pending], future::resolved, logical(1L))]
        if (length(ready)) {
          for (k in ready) collect(k, future::value(futs[[k]]))
          pending <- setdiff(pending, ready)
        } else {
          Sys.sleep(0.2)
        }
      }
    } else {
      for (k in seq_len(n_jobs)) collect(k, run_k(k))
    }
  } else if (use_parallel) {
    oplan <- future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(oplan), add = TRUE)
    job_results <- future.apply::future_lapply(seq_len(n_jobs), run_k,
                                               future.seed = NULL)
  } else {
    job_results <- lapply(seq_len(n_jobs), run_k)
  }

  # ── Combine ────────────────────────────────────────────────────────────────
  plans   <- rbindlist(lapply(job_results, `[[`, "plans"))
  summary <- rbindlist(lapply(job_results, `[[`, "summary"))
  final_populations <- stats::setNames(
    lapply(job_results, `[[`, "final"),
    vapply(job_results, `[[`, character(1L), "key"))

  # trajectories compat view (Rep/Cycle/Method/MeanGEBV/SDev/BlockCoverage)
  trajectories <- summary[cycle > 0L, .(
    Rep           = rep,
    Cycle         = cycle,
    Method        = sel_scheme,
    MeanGEBV      = mean_F1mean,
    SDev          = sqrt(var_F1mean),
    BlockCoverage = NA_real_)]

  # ── Write output ────────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    tag <- paste(unique(scheme), collapse = "-")
    fwrite(plans,   file.path(output_dir, paste0("stacking_plans_",   tag, ".txt")),
           sep = "\t")
    fwrite(summary, file.path(output_dir, paste0("stacking_summary_", tag, ".txt")),
           sep = "\t")
    p <- tryCatch(
      plot_simulation_trajectories(trajectories, trait_name = tag),
      error = function(e) NULL)
    if (!is.null(p))
      tryCatch(
        ggplot2::ggsave(file.path(output_dir, paste0("Trajectories_", tag, ".png")),
                        p, width = 10, height = 5, dpi = 150),
        error = function(e) NULL)
  }

  list(
    plans             = plans,
    summary           = summary,
    ideal_genotypes   = ideal_genotypes,
    final_populations = final_populations,
    trajectories      = trajectories
  )
}
