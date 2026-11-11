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
# Steps 2-3 are the parents-first protocol (`selection_unit = "parents"`).
# Under `selection_unit = "pairs"` (direct pair selection, thread #031) the
# scheme returns the crosses themselves: .stk_select_pairs() dispatches to the
# family scorers of stacking_pairs.R (TS mid-parent / usefulness, OHS M[i, j],
# GBS on the family's expected offspring) under the budget resolved by
# .stk_selection_budget() — the same number of families as the traversal, or
# the same number of distinct parents — and make_cross_design("custom") turns
# the family table into the cycle's design.
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
# Internal — the haplotype column names of a wide local-GEBV table (sc/F24)
# ---------------------------------------------------------------------------

# The one parser of "<Sample>.<haplotype_idx>" column names (0.9.0, D4-1): the
# eight hand-written `sub("\\.[0-9]+$", "", ...)` copies of 0.8.8 went through
# this function, with the same regex, so nothing moves. `cols` are the sample
# columns of the table (everything but "Block"). Returns
#   cols     the input
#   base     the sample name of every column
#   idx      the haplotype index of every column (integer; NA when absent)
#   samples  unique(base), in first-appearance order
#   ploidy   `ploidy` when given, else the largest number of columns one
#            sample has (the 0.8.8 inference)
.hap_cols <- function(cols, ploidy = NULL) {
  cols <- as.character(cols)
  base <- sub("\\.[0-9]+$", "", cols)
  idx  <- suppressWarnings(as.integer(sub("^.*\\.([0-9]+)$", "\\1", cols)))
  idx[!grepl("\\.[0-9]+$", cols)] <- NA_integer_
  samples <- unique(base)
  if (is.null(ploidy))
    ploidy <- if (length(cols)) max(tabulate(match(base, samples))) else 0L
  list(cols = cols, base = base, idx = idx, samples = samples,
       ploidy = as.integer(ploidy))
}


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
# Internal — the scheme registry (0.9.0, D4 §3.5; sc/F07)
# ---------------------------------------------------------------------------

# One entry per selection scheme. Up to 0.8.8 the scheme name was matched as a
# string at twelve sites of this file; every one of them now reads the registry:
#   parents             selector under selection_unit = "parents":
#                       function(lg, geno_meta, n, ideal, opts) -> parent names
#   pairs               family scorer under selection_unit = "pairs":
#                       function(lg, geno_meta, ideal, opts) -> pair_selection,
#                       or NULL when the scheme has no family form (OCS, until
#                       the KCT rename)
#   needs_ideal         the scheme selects against the ideal genotype (so the
#                       ideal is built whatever `metrics` says)
#   direction           selection direction the scheme honours: "both",
#                       "high", or "argmax" (= "high" unless the ideal builder
#                       is argmax; form row 19, D4-6)
#   pin_founder_blocks  under ideal_scope = "founder" the scheme is scored on
#                       the founders' block set
#   in_default          the scheme is in the driver's default `scheme` vector
# The selector bodies are the 0.8.8 switch() arms, moved verbatim.
# `.STK_BASELINE` is the scheme of selfing generations and of the decision
# cycles before `ga_start_by`.

.STK_BASELINE <- "TS"

.stk_sel_ts <- function(lg, geno_meta, n, ideal, opts)
  .stk_select_truncation(lg, n, opts$select_priority)

# The GA controls of the spec (`ga_control`, sc/F30): the engines' own defaults
# when the list is empty, so nothing moves without it.
.stk_ga <- function(opts) {
  ga <- opts$ga_control
  list(pop_size = ga$pop_size %||% 100L,
       n_generations = ga$n_generations %||% 200L,
       mutation_rate = ga$mutation_rate %||% 0.1)
}

.stk_sel_ohs <- function(lg, geno_meta, n, ideal, opts) {
  ga <- .stk_ga(opts)
  ohs_selection(lg, n_parents = n, method = opts$ohs_method,
                elite_ratio = opts$ohs_elite_ratio,
                elite_level = opts$ohs_elite_level,
                pop_size = ga$pop_size, n_generations = ga$n_generations,
                mutation_rate = ga$mutation_rate,
                seed = opts$seed)
}

.stk_sel_ocs <- function(lg, geno_meta, n, ideal, opts) {
  v   <- .stk_total_gebv(lg)
  grm <- tryCatch(suppressMessages(compute_grm(geno_meta)),
                  error = function(e) {
                    .hs_note("selection", "ocs_grm_failed",
                             paste0("compute_grm() failed (", conditionMessage(e),
                                    "); OCS ran as truncation on GEBV."))
                    NULL
                  })
  ocs_selection(v, kinship_matrix = grm, n_parents = n,
                max_kinship = opts$max_kinship,
                freq_penalty_weight = opts$freq_penalty_weight,
                seed = opts$seed)
}

.stk_sel_gbs <- function(lg, geno_meta, n, ideal, opts)
  gbs_selection(lg, ideal_genotypes = ideal, n_parents = n,
                gbs_method = opts$gbs_method,
                gbs_objective = opts$gbs_objective,
                gbs_engine = opts$gbs_engine,
                penalty_mode = opts$penalty_mode,
                gbs_penalty = opts$gbs_penalty,
                penalty_scale = opts$penalty_scale,
                penalty_rank_scale = opts$penalty_rank_scale,
                gbs_effect_weight = opts$gbs_effect_weight,
                select_priority = opts$select_priority,
                ga_control = if (is.null(opts$ga_control)) list() else opts$ga_control,
                seed = opts$seed)

.stk_pairs_ts <- function(lg, geno_meta, ideal, opts)
  .stk_pairs_truncation(lg, opts$pair_budget, opts)

.stk_pairs_ohs <- function(lg, geno_meta, ideal, opts) {
  b  <- opts$pair_budget
  ga <- .stk_ga(opts)
  ohs_pair_selection(lg, pair_budget = b$pair_budget,
                     top_founders = b$n_parents_ref,
                     allow_selfing = b$allow_selfing,
                     ohs_method = opts$ohs_method,
                     elite_ratio = opts$ohs_elite_ratio,
                     elite_level = opts$ohs_elite_level,
                     shortlist_families = opts$ohs_shortlist_families,
                     max_use = opts$max_use, overshoot = opts$overshoot,
                     pop_size = ga$pop_size, n_generations = ga$n_generations,
                     mutation_rate = ga$mutation_rate,
                     seed = opts$seed)
}

.stk_pairs_gbs <- function(lg, geno_meta, ideal, opts) {
  b <- opts$pair_budget
  gbs_pair_selection(lg, ideal_genotypes = ideal,
                     pair_budget = b$pair_budget,
                     top_founders = b$n_parents_ref,
                     num_offspring = opts$num_offspring,
                     allow_selfing = b$allow_selfing,
                     gbs_method = opts$gbs_method,
                     gbs_objective = opts$gbs_objective,
                     penalty_mode = opts$penalty_mode,
                     gbs_penalty = opts$gbs_penalty,
                     penalty_scale = opts$penalty_scale,
                     penalty_rank_scale = opts$penalty_rank_scale,
                     gbs_effect_weight = opts$gbs_effect_weight,
                     gbs_pair_level = opts$gbs_pair_level,
                     pair_freq = opts$pair_freq,
                     pair_penalty = opts$pair_penalty,
                     pair_pool = opts$pair_pool,
                     shortlist_k = opts$shortlist_k,
                     max_use = opts$max_use, overshoot = opts$overshoot,
                     select_priority = opts$select_priority)
}

.STK_SCHEMES <- list(
  TS  = list(parents = .stk_sel_ts,  pairs = .stk_pairs_ts,  needs_ideal = FALSE,
             direction = "both",   pin_founder_blocks = FALSE, in_default = TRUE),
  OHS = list(parents = .stk_sel_ohs, pairs = .stk_pairs_ohs, needs_ideal = FALSE,
             direction = "high",   pin_founder_blocks = FALSE, in_default = TRUE),
  OCS = list(parents = .stk_sel_ocs, pairs = NULL,           needs_ideal = FALSE,
             direction = "high",   pin_founder_blocks = FALSE, in_default = FALSE),
  GBS = list(parents = .stk_sel_gbs, pairs = .stk_pairs_gbs, needs_ideal = TRUE,
             direction = "argmax", pin_founder_blocks = TRUE,  in_default = TRUE))

# The registry entry of a scheme, or an internal error for an unknown name
# (the entry points validate `scheme` against names(.STK_SCHEMES) first).
.stk_scheme <- function(scheme) {
  s <- .STK_SCHEMES[[scheme]]
  if (is.null(s))
    .hs_abort(class = "hs_error_internal", "unknown scheme '", scheme, "'.",
              call. = FALSE)
  s
}

# A property of several schemes at once (e.g. .stk_scheme_flag(scheme, "needs_ideal")).
.stk_scheme_flag <- function(schemes, what) {
  vapply(schemes, function(s) isTRUE(.stk_scheme(s)[[what]]), logical(1L),
         USE.NAMES = FALSE)
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
  sel <- .stk_scheme(scheme)$parents(localgebv, geno_meta, n_parents,
                                     ideal_genotypes, opts)
  sel <- unique(sel)
  if (length(sel) < 2L) {
    .hs_note("selection", "fallback_truncation",
             paste0("scheme ", scheme, " returned ", length(sel), " parent(s); ",
                    "truncation selection filled the parent set."))
    sel <- .stk_select_truncation(localgebv, max(n_parents, 2L),
                                  opts$select_priority)
  }
  sel
}


# ---------------------------------------------------------------------------
# Internal — selection unit: parents-first vs direct pairs (budget + dispatch)
# ---------------------------------------------------------------------------

# Resolve the per-cycle budget the driver runs under. Parents-first
# (`selection_unit = "parents"`) crosses `top_founders` parents in a full
# traversal, so the family count is fixed at .stk_n_families(). Direct pairing
# ("pairs") must be tied to that protocol by one of the two conflicting budgets
# of pair_selection_budget(): the same number of families, or the same number
# of distinct parents. `allow_selfing` is what the design actually does: under
# parents-first, make_cross_design("all_pairs") enumerates combn(parents, 2)
# and never contains a self-cross (INDEX discrepancy 27 / Q24 c), so the driver
# passes FALSE; under pairs the candidate pool, this budget and the design all
# take `!avoid_selfing`, so the three agree.
.stk_selection_budget <- function(selection_unit, pair_budget, top_founders,
                                  allow_selfing = FALSE) {
  n <- as.integer(top_founders)
  if (selection_unit == "parents") {
    full <- .stk_n_families(n, allow_selfing)
    return(list(selection_unit = "parents", pair_budget = NA_character_,
                n_parents_ref = n, allow_selfing = isTRUE(allow_selfing),
                n_families = full, n_families_min = full, n_families_max = full,
                n_parents_max = n))
  }
  if (is.null(pair_budget))
    .hs_abort("run_stacking_simulation: selection_unit = \"pairs\" needs ",
         "'pair_budget' = \"families\" (same offspring population size as the ",
         "parents-first traversal) or \"parents\" (same number of distinct ",
         "parents) \u2014 the two conflict, so one must be chosen.")
  pair_budget <- match.arg(pair_budget, c("families", "parents"))
  c(list(selection_unit = "pairs"),
    pair_selection_budget(n, pair_budget, allow_selfing))
}

# One-line human-readable statement of a resolved budget (verbose header and
# the output check of .stk_select_pairs()).
.stk_budget_text <- function(b) {
  selfing <- if (isTRUE(b$allow_selfing)) "self-crosses counted"
             else "self-crosses not in the all-pairs design"
  if (b$selection_unit == "parents")
    sprintf("parents-first: %d parents -> %d families per cycle (%s)",
            b$n_parents_ref, b$n_families, selfing)
  else if (b$pair_budget == "families")
    sprintf(paste0("direct pairing, pair_budget = \"families\": %d families per ",
                   "cycle (= the traversal of %d parents; %s), distinct parents free"),
            b$n_families, b$n_parents_ref, selfing)
  else
    sprintf(paste0("direct pairing, pair_budget = \"parents\": add pairs until %d ",
                   "distinct parents are in use (%d..%d families per cycle; %s)"),
            b$n_parents_max, b$n_families_min, b$n_families_max, selfing)
}

# Direct pair selection (thread #031): the scheme returns the families of the
# cycle as a "pair_selection" table (Parent1, Parent2, score, rank) under the
# budget in `opts$pair_budget`. Scorers live in stacking_pairs.R; this is the
# dispatch (the pairs counterpart of .stk_select_parents()) plus the output
# check. A scheme that returns no family errors with its name — no truncation
# fallback for this unit. OCS has no family form until the KCT rename (plan
# #030 Phase 0 / pair plan Phase C).
.stk_select_pairs <- function(scheme, localgebv, geno_meta, ideal_genotypes,
                              opts) {
  b <- opts$pair_budget
  f <- .stk_scheme(scheme)$pairs
  if (is.null(f))
    .hs_abort("run_stacking_simulation: scheme \"", scheme, "\" has no direct-pair ",
              "form yet (it follows the KCT rename; pair plan Phase C). Use ",
              "selection_unit = \"parents\" for ", scheme, ".", call. = FALSE)
  pairs <- f(localgebv, geno_meta, ideal_genotypes, opts)
  .stk_check_pairs(pairs, dimnames(geno_meta$geno_arrays)[[1L]], b, scheme)
  pairs
}

# Output contract of a family scorer: >= 1 family, parents in the pool, no
# duplicate family (in either orientation), a self only when allowed.
.stk_check_pairs <- function(pairs, pool, budget, scheme) {
  if (is.null(pairs) || nrow(pairs) == 0L)
    .hs_abort("run_stacking_simulation: scheme \"", scheme, "\" returned no family ",
         "under selection_unit = \"pairs\" (", .stk_budget_text(budget), ").",
         call. = FALSE)
  miss <- setdiff(unique(c(pairs$Parent1, pairs$Parent2)), pool)
  if (length(miss))
    .hs_abort(".stk_select_pairs: parent(s) not in the pool: ",
         paste(miss, collapse = ", "))
  key <- paste(pmin(pairs$Parent1, pairs$Parent2), pmax(pairs$Parent1, pairs$Parent2))
  if (anyDuplicated(key))
    .hs_abort(".stk_select_pairs: scheme \"", scheme, "\" returned a duplicate family.")
  if (!isTRUE(budget$allow_selfing) && any(pairs$Parent1 == pairs$Parent2))
    .hs_abort(".stk_select_pairs: scheme \"", scheme, "\" returned a self-cross with ",
         "avoid_selfing = TRUE.")
  invisible(TRUE)
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
  # bare dimnames swap would leave sample_info$Sample stale and fail
  # validation). One rebuild path (0.9.0): blocks, effects and annotation carry
  # to the renamed offspring (#033); the founders' phenotype and family history
  # are declared losses (X7).
  off <- .rebuild_geno_meta(src, geno_arrays = arr, sample_info = fam,
                            .drop = c("family_history", "phenotype"),
                            .caller = ".advance_generation")

  bp[, CrossID := NewID]
  list(off = off, bp = bp)
}

# The per-cross plan of a generation from its per-offspring table (one row per
# family: F1count, the mean F1mean of its offspring, the compact CrossID). Built
# after the offspring F1mean is known — for the block backend that is the
# shared local-GEBV pass of .stk_advance() (0.9.0, D4-2), not the crossing.
.stk_cross_plans <- function(bp, cycle, opts) {
  cross_plans <- bp[, .(F1count = .N, F1mean = mean(F1mean, na.rm = TRUE)),
                    by = .(FamID, Parent1, Parent2)]
  cross_plans[, CrossID := sprintf(opts$id_fmt_cross, cycle, FamID)]
  cross_plans
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
#
# One local-GEBV pass per generation (0.9.0, form row 20 a / D4-2): the
# offspring are scored ONCE at haplotype level, on the renamed container, and
# that table serves three uses that were three passes up to 0.8.8 — the block
# backend's F1mean (its genotype-level pass; the haplotype sum is identical,
# probe D4 §3.4), the d2ideal / ideal_hap_gebv metrics, and the next
# generation's parent table, handed back as `off_lg` and passed in as `lg`.
# `lg = NULL` computes the parent table here (first generation, or after a
# failed pass). Under genomicSimulation F1mean stays gS's native breeding
# value and the pass runs only when a metric needs it.
.stk_advance <- function(current, scheme, effect_df, opts, cycle, session = NULL,
                         lg = NULL) {
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

  if (is.null(lg))
    lg <- compute_profile_localgebv(current, effect_df, level = "haplotype")
  # A founder-fixed ideal defines the block set of GBS selection for the whole
  # run (the founders' `select_top_blocks`), so a target block cannot silently
  # drop out of a later cycle's top-N. Only GBS selection is pinned: TS/OHS/OCS
  # keep the per-cycle top-N, and the d2ideal metric is measured on the full
  # offspring table, so neither is affected by the scope.
  use_founder <- identical(opts$ideal_scope, "founder") &&
    !is.null(opts$founder_ideal)
  sch <- .stk_scheme(scheme)
  lg <- if (use_founder && isTRUE(sch$pin_founder_blocks))
    .stk_limit_blocks(lg, opts$founder_blocks, NULL)
  else
    .stk_limit_blocks(lg, opts$selected_blocks, opts$select_top_blocks)
  if (nrow(lg) == 0L)
    .hs_abort(".stk_advance: no blocks remain for selection (check 'selected_blocks' ",
         "/ block assignments).")
  lap("parent_gebv")

  # Ideal genotype has two INDEPENDENT roles: (a) the SELECTION target for GBS
  # (built whenever GBS runs, never influenced by 'metrics'), and (b) the input to
  # the d2ideal MEASUREMENT. Built once and shared when both apply, so metrics can
  # never change GBS selection.
  need_ideal <- isTRUE(sch$needs_ideal) || want("d2ideal")
  ideal <- if (!need_ideal) NULL
  else if (use_founder) opts$founder_ideal
  else construct_ideal_genotype(lg, method = opts$construct_method,
                                k_penalty = opts$k_penalty,
                                select_priority = opts$select_priority)
  lap("ideal")

  # Selection unit. Parents-first: a parent subgroup, then its full traversal.
  # Direct pairing: the scheme returns the crosses themselves under the
  # resolved budget; the family table (a fresh copy — make_cross_design()
  # adds columns by reference) becomes the custom design, self-crosses kept
  # only when the budget allowed them as candidates.
  pairs <- NULL
  if (identical(opts$selection_unit, "pairs")) {
    pairs   <- .stk_select_pairs(scheme, lg, current, ideal, opts)
    parents <- unique(c(pairs$Parent1, pairs$Parent2))
    lap("selection")
    cross_design <- make_cross_design(
      parents, method = "custom",
      cross_table   = data.table(Parent1 = pairs$Parent1, Parent2 = pairs$Parent2),
      avoid_selfing = !isTRUE(opts$pair_budget$allow_selfing),
      n_offspring   = opts$num_offspring)
  } else {
    parents <- .stk_select_parents(scheme, lg, current, opts$top_founders,
                                   ideal, opts)
    lap("selection")
    cross_design <- make_cross_design(parents, method = "all_pairs",
                                      avoid_selfing = opts$avoid_selfing,
                                      n_offspring = opts$num_offspring)
  }

  # The driver rebuilds the offspring FAM itself in .stk_rename_generation(),
  # so the genomicSimulation read-back skips its own (io/F25: one rebuild less
  # per generation; the block backend never built one).
  # The block backend's F1mean comes from the shared pass below (f1mean = FALSE
  # skips its own genotype-level pass); genomicSimulation returns its native BV.
  res <- .cross_backend_dispatch(opts$cross_backend, current, cross_design,
                                 effect_df = effect_df, work_dir = opts$work_dir,
                                 session = session, offspring_fam = FALSE,
                                 f1mean = FALSE)

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

  # ── The one offspring local-GEBV pass (row 20 a) ─────────────────────────────
  # Always for the block backend (it is the F1mean); for genomicSimulation only
  # when a metric needs it. On error the generation keeps 0.8.8's behaviour —
  # F1mean NA, metrics NA, the next generation recomputes its parent table —
  # and says so (note code f1mean_failed).
  is_block   <- identical(opts$cross_backend, "block")
  need_offlg <- is_block || want("d2ideal") || want("ideal_hap_gebv")
  off_lg <- if (need_offlg)
    tryCatch(compute_profile_localgebv(off, effect_df, level = "haplotype"),
             error = function(e) {
               .hs_note("simulation", "f1mean_failed",
                        paste0("offspring local GEBV pass failed (",
                               conditionMessage(e), "); ",
                               if (is_block) "F1mean and " else "",
                               "the offspring metrics are NA for the generation."))
               NULL
             })
  else NULL
  if (is_block) {
    f1 <- if (is.null(off_lg) || nrow(off_lg) == 0L) rep(NA_real_, nrow(gen$bp))
          else unname(.stk_total_gebv(off_lg)[gen$bp$NewID])
    gen$bp[, F1mean := f1]
  }
  lap("offspring_gebv")

  # Per-cross plan (F1count / F1mean per family), the per-family selection score
  # (direct pairing; NA under parents-first) and the cycle's family /
  # distinct-parent counts — both units, so the two protocols are compared on
  # the same columns.
  cross_plans <- .stk_cross_plans(gen$bp, cycle, opts)
  cross_plans[, pair_score := NA_real_]
  if (!is.null(pairs))
    cross_plans[data.table(Parent1 = pairs$Parent1, Parent2 = pairs$Parent2,
                           s = pairs$score),
                on = .(Parent1, Parent2), pair_score := i.s]
  n_families     <- nrow(cross_plans)
  n_parents_used <- length(unique(c(cross_plans$Parent1, cross_plans$Parent2)))

  # A metric that fails is NA for the generation and says so (note code
  # metric_failed, with the metric's name and the error; sc/F19).
  metric <- function(name, expr) tryCatch(expr, error = function(e) {
    .hs_note("metrics", "metric_failed",
             paste0(name, " failed (", conditionMessage(e), "); NA for the generation."))
    NA_real_
  })
  d2 <- if (want("d2ideal") && !is.null(ideal) &&
            !is.null(off_lg) && nrow(off_lg) > 0L)
    metric("d2ideal", .distance_to_ideal(off_lg, ideal, gbs_method = opts$gbs_method))
  else NA_real_

  ideal_hap <- if (want("ideal_hap_gebv") && !is.null(off_lg) && nrow(off_lg) > 0L)
    metric("ideal_hap_gebv", .stk_ideal_hap_gebv(off_lg, ploidy = opts$ploidy))
  else NA_real_

  # Nei's (1973) gene diversity of the offspring population (marker-level).
  gdiv <- if (want("gene_diversity"))
    metric("gene_diversity", as.numeric(compute_gene_diversity(off, level = "marker")))
  else NA_real_

  # Coancestry = mean block-IBS relatedness among the OFFSPRING (a diversity
  # measure of the produced population, distinct from Nei gene diversity).
  coan <- if (want("coancestry"))
    metric("coancestry", .stk_coancestry(off, dimnames(off$geno_arrays)[[1L]],
                                         blocks = unique(as.character(lg$Block))))
  else NA_real_
  lap("metrics")

  list(off             = off,
       off_lg          = off_lg,
       cross_plans     = cross_plans,
       n_families      = n_families,
       n_parents_used  = n_parents_used,
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
  paste0(sprintf("  [%s rep %d] step time \u2014 total %.1fs; dominant: %s\n",
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
                    "D2Ideal", "coancestry", "ideal_phenotype", "pair_score"))
  cp[]
}

# One `summary` row (FR5 columns) for a cycle. The founder row (cycle 0) goes
# through the same function from the `adv`-shaped list .stk_founder_measures()
# builds (sc/F26), so the two rows can never drift apart in their arithmetic.
.stk_summary_record <- function(adv, select_type, cycle, rep, sel_scheme) {
  f1 <- adv$f1
  data.table(cycle                = cycle,
             select_type          = select_type,
             rep                  = rep,
             sel_scheme           = sel_scheme,
             mean_F1mean          = mean(f1, na.rm = TRUE),
             var_F1mean           = stats::var(f1, na.rm = TRUE),
             max_F1mean           = if (any(!is.na(f1))) max(f1, na.rm = TRUE) else NA_real_,
             mean_D2Ideal         = adv$d2,
             mean_coancestry      = adv$coan,
             mean_gene_diversity  = adv$gene_div,
             mean_ideal_phenotype = adv$ideal_phenotype,
             ideal_hap_gebv       = adv$ideal_hap,
             marker_size          = adv$n_marker,
             n_families           = as.integer(adv$n_families),
             n_parents_used       = as.integer(adv$n_parents_used))
}

# The founder population measured like a generation (the cycle-0 record of
# 0.8.8, `start_fam_record`): each optional metric is filled only when its
# token is in `opts$metrics`; the ideal is the founder ideal already built by
# the entry point (identical input, no random draw — sc/F26), never rebuilt.
.stk_founder_measures <- function(geno_meta, founder_lg, founder_ideal, opts) {
  want <- function(m) m %in% opts$metrics
  v0 <- .stk_total_gebv(founder_lg)
  f_ideal <- if (want("d2ideal")) founder_ideal else NULL
  list(f1              = v0,
       d2              = if (want("d2ideal"))
                           .distance_to_ideal(founder_lg, f_ideal,
                                              gbs_method = opts$gbs_method)
                         else NA_real_,
       coan            = if (want("coancestry"))
                           tryCatch(.stk_coancestry(
                             geno_meta, dimnames(geno_meta$geno_arrays)[[1L]],
                             blocks = unique(as.character(founder_lg$Block))),
                             error = function(e) NA_real_)
                         else NA_real_,
       gene_div        = if (want("gene_diversity"))
                           tryCatch(as.numeric(compute_gene_diversity(geno_meta, level = "marker")),
                                    error = function(e) NA_real_)
                         else NA_real_,
       ideal_phenotype = if (is.null(f_ideal)) NA_real_
                         else as.numeric(attr(f_ideal, "ideal_phenotype")),
       ideal_hap       = if (want("ideal_hap_gebv"))
                           .stk_ideal_hap_gebv(founder_lg, ploidy = opts$ploidy)
                         else NA_real_,
       n_marker        = dim(geno_meta$geno_arrays)[2L],
       n_families      = NA_integer_,
       n_parents_used  = NA_integer_)
}

# The `trajectories` compat view of a summary table (Rep / Cycle / Method /
# MeanGEBV / SDev / BlockCoverage), one definition for the final combine and
# the live monitor (sc/F23).
.stk_trajectories_view <- function(summary) {
  s <- as.data.table(summary)
  # The legend name is method + mode (owner's naming rule, 2026-11-08): the
  # pairs mode is "cross design", suffix -CD; OHS-block x cross is BECD.
  # Calls keep `scheme` and `selection_unit` separate; only the label joins
  # them. A summary without the unit column (0.8.8 tables) keeps the scheme.
  label <- s$sel_scheme
  if ("selection_unit" %in% names(s)) {
    cd <- !is.na(s$selection_unit) & s$selection_unit == "pairs"
    be <- cd & s$sel_scheme %in% "OHS" & "ohs_elite_level" %in% names(s) &
          !is.na(s$ohs_elite_level) & s$ohs_elite_level == "block"
    label <- ifelse(be, "BECD", ifelse(cd, paste0(s$sel_scheme, "-CD"), s$sel_scheme))
  }
  s[, Method := label]
  s[cycle > 0L, .(
    Rep           = rep,
    Cycle         = cycle,
    Method        = Method,
    MeanGEBV      = mean_F1mean,
    SDev          = sqrt(var_F1mean),
    BlockCoverage = NA_real_)]
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
#'     \code{"Greedy"}; \code{ohs_elite_ratio} /
#'     \code{ohs_elite_level} restrict the objective to the set's elite
#'     families — the crossing design is unchanged).}
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
#' @section Selection unit:
#' Under \code{selection_unit = "parents"} (default) every scheme selects a
#' \emph{parent subgroup} of \code{top_founders}, which is then crossed in a
#' full traversal (\code{choose(top_founders, 2)} families per cycle). Under
#' \code{selection_unit = "pairs"} (direct pair selection) every scheme scores
#' \emph{families} — the expected F1 of one cross — and one greedy adds families
#' until the budget stops it: \code{pair_budget = "families"} makes as many
#' families as the traversal (distinct parents free), \code{"parents"} adds
#' families until \code{top_founders} distinct parents are in use (families
#' free; see \code{\link{pair_selection_budget}}). What each scheme scores:
#' \describe{
#'   \item{\code{"TS"}}{the mid-parent value (the expected F1 mean), or with
#'     \code{ts_pair_score = "usefulness"} the mid-parent plus
#'     \code{usefulness_k} F1 standard deviations (exact under the block model).}
#'   \item{\code{"OHS"}}{the family's stacking ceiling \code{M[i, j]}
#'     (\code{\link{ohs_pair_selection}}; under \code{"parents"} the best
#'     parent set and its traversal). With \code{ohs_elite_level = "block"} and a
#'     ratio below 1 the
#'     family set's per-block elite is scored instead
#'     (\code{ohs_shortlist_families} candidate families).}
#'   \item{\code{"GBS"}}{the ideal-genotype objective on the family's expected
#'     offspring (\code{\link{gbs_pair_selection}}: \code{gbs_pair_level},
#'     \code{pair_freq}, \code{pair_penalty}, \code{pair_pool}; the default
#'     objective is \code{"set"}).}
#'   \item{\code{"OCS"}}{no family form yet (follows the KCT rename).}
#' }
#' \code{summary} reports the cycle's \code{n_families} and
#' \code{n_parents_used} for both units, and \code{plans} the per-family
#' \code{pair_score}.
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
#' @param effect_df Marker effects: one \code{(Marker, Allele, Effect)} table, a
#'   named per-trait list of them, or one long table with a \code{Trait}
#'   column. When \code{NULL}, the \code{GenoMeta}'s attached
#'   \code{marker_effects} are used. This pipeline breeds on one trait: with
#'   several traits pass \code{trait = "<name>"} (since 0.9.0 the call stops
#'   otherwise, naming the traits; up to 0.8.8 the first trait was taken
#'   silently).
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
#'   Must not exceed the population (since 0.9.0 that is an error; up to 0.8.8
#'   the set was capped silently).
#' @param avoid_selfing Logical. Exclude self-crosses. Default \code{TRUE}.
#'   The \code{"all_pairs"} design of the parents-first protocol never contains
#'   a self-cross in any case (it enumerates \code{combn(parents, 2)}), so the
#'   flag acts only under \code{selection_unit = "pairs"}, where
#'   \code{avoid_selfing = FALSE} offers the selfs as candidate families and
#'   counts them in the budget.
#' @param selection_unit What a cycle selects. \code{"parents"} (default): a
#'   subgroup of \code{top_founders} parents, then every pair of them is
#'   crossed. \code{"pairs"}: the crosses themselves are chosen (direct
#'   pairing; see the \emph{Selection unit} section). Requires
#'   \code{pair_budget}.
#' @param pair_budget How a direct-pairing cycle is matched to the parents-first
#'   protocol; required when \code{selection_unit = "pairs"}, ignored otherwise.
#'   The two options conflict, so there is no default. \code{"families"}: the
#'   same offspring population size — as many families as the traversal of
#'   \code{top_founders} parents makes (\code{choose(top_founders, 2)}), with the
#'   number of distinct parents free. \code{"parents"}: the same number of
#'   parents — pairs are added until \code{top_founders} distinct parent names
#'   are in use, with the number of families free. See
#'   \code{\link{pair_selection_budget}} for the numbers a given
#'   \code{top_founders} implies.
#' @param gbs_pair_level,pair_freq,pair_penalty,pair_pool,shortlist_k Direct
#'   pairing, GBS: the level the family is scored at (\code{"population"} /
#'   \code{"family"}), the frequency projection (\code{NULL} = by budget),
#'   the absent-state penalty form (\code{"loss"} / \code{"hard"}), the
#'   candidate space and its shortlist size — see \code{\link{gbs_pair_selection}}.
#' @param max_use Direct pairing: per-parent usage cap (families a parent may
#'   appear in). Default \code{Inf}.
#' @param overshoot Direct pairing, \code{pair_budget = "parents"}: a family
#'   that would exceed the parent limit is \code{"skip"}ped (default; the search
#'   goes on among the parents in use, up to the traversal) or \code{"stop"}s
#'   the search.
#' @param ts_pair_score,usefulness_k Direct pairing, TS: \code{"midparent"}
#'   (default) or \code{"usefulness"} = mid-parent \code{+ usefulness_k *
#'   sd(F1)}. Default \code{usefulness_k = 1}.
#' @param scheme Character vector of schemes to run/compare: any of \code{"TS"},
#'   \code{"OHS"}, \code{"OCS"}, \code{"GBS"}. Default \code{c("TS", "OHS",
#'   "GBS")} since 0.9.0 (OCS stays available by name; up to 0.8.8 the default
#'   ran all four).
#' @param ohs_method OHS engine (see \code{\link{ohs_selection}}).
#' @param ohs_elite_ratio Fraction in \code{(0, 1]} of a set's families the OHS
#'   objective counts (see \code{\link{ohs_selection}}). Default \code{1} =
#'   every family; the crossing design has the same size at every ratio.
#' @param ohs_elite_level \code{"pair"} (default) or \code{"block"}: the elite
#'   families are chosen once per set, or once per block (the ratio's small
#'   limit is then the best two haplotypes of every block — Kemper's ultimate
#'   genotype — rather than the best single cross).
#' @param ohs_shortlist_families Direct pairing, OHS block level: candidate
#'   families held with their per-block values; \code{NULL} = the larger of
#'   \code{4 * n_families} and \code{100} (see \code{\link{ohs_pair_selection}}).
#' @param gbs_method GBS deviation form, \code{"ARCSIN"} (default) or
#'   \code{"LOG"} — Kemper's two modes (see \code{\link{gbs_selection}}).
#' @param gbs_preset \code{"haploStack"} (default) or \code{"kemper"} (the
#'   published method: fixed \code{-100} penalty). Fills \code{gbs_objective},
#'   \code{gbs_engine}, \code{penalty_mode}, \code{gbs_penalty},
#'   \code{penalty_scale} and \code{gbs_effect_weight} when they are \code{NULL};
#'   an explicit value always wins. See \code{\link{gbs_selection}}.
#' @param gbs_objective GBS objective, \code{"candidate"} or \code{"set"};
#'   \code{NULL} (default) takes the preset value (\code{"candidate"} under
#'   parents-first; \code{"set"} under direct pairing).
#' @param gbs_engine GBS engine, \code{"greedy"} or \code{"GA"}; \code{NULL}
#'   (default) takes the preset value (\code{"greedy"}).
#' @param penalty_mode GBS absent-state penalty scaling, \code{"dynamic"} or
#'   \code{"fixed"}; \code{NULL} (default) takes the preset value
#'   (\code{"dynamic"} for \code{"haploStack"}, \code{"fixed"} for \code{"kemper"}).
#' @param gbs_penalty GBS base penalty (\code{<= 0}); \code{NULL} (default) takes
#'   the preset value \code{-10}.
#' @param penalty_scale GBS penalty scale or cap (\code{>= 1}); \code{NULL}
#'   (default) takes the preset value \code{10}.
#' @param penalty_rank_scale GBS per-rank multipliers of the absent-state
#'   penalty (not preset-backed; default \code{1} = the behaviour before
#'   0.7.11). See the rank-scaled absence penalty section of
#'   \code{\link{gbs_selection}}.
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
#'   direction. \code{"low"} is honoured by TS (both selection units) and by
#'   GBS with \code{construct_method = "argmax"}; OHS, OCS and GBS with the
#'   \code{"custom"} or \code{"ratio"} builder select for the high tail
#'   whatever the priority says, so \code{"low"} is refused with them (0.9.0;
#'   up to 0.8.8 it was ignored silently).
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
#' @param output_dir \code{NULL} (default since 0.9.0): the result is returned
#'   in memory and nothing is written. A directory (created if absent) receives
#'   the files 0.8.8 wrote, through \code{\link{write_outputs}}:
#'   \code{stacking_plans_<tag>.txt}, \code{stacking_summary_<tag>.txt} and
#'   \code{Trajectories_<tag>.png} (\code{<tag>} = the schemes joined by
#'   \code{"-"}); the result is then returned invisibly. Up to 0.8.8 the default
#'   was \code{"stacking_output"} in the working directory.
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
#' @param trait The trait to breed on when \code{effect_df} (or the
#'   container's \code{marker_effects}) holds several traits. Default
#'   \code{NULL}.
#' @param keep_populations \code{"final"} (default): keep every job's final
#'   population in \code{final_populations}; \code{"none"}: an empty list.
#' @param backend_opts Options of the crossing backend (genomicSimulation):
#'   a list with any of \code{cM_per_mb}, \code{map_unit}, \code{keep_files};
#'   see \code{\link{stk_spec}}.
#' @param ga_control \code{list(pop_size, n_generations, mutation_rate)} for
#'   the GA engines; see \code{\link{stk_spec}}.
#' @param progress \code{NULL} (default) or a function of one argument, the
#'   package's progress contract: with \code{n_cores = 1} it is called after
#'   every generation (\code{unit = "cycle"}, see \code{\link{stk_step}}), with
#'   \code{n_cores > 1} after every job (\code{unit = "job"}). Returning
#'   \code{FALSE} cancels the run with an \code{hs_error_cancelled} error whose
#'   \code{partial} field holds the result of the generations and jobs done so
#'   far.
#' @param catalogue_results,trait_name,select_method Deprecated since 0.9.0
#'   (each warns on every call, see \link{haploStack-deprecated}) back-compat
#'   inputs. \code{catalogue_results}/\code{trait_name} supply the founder
#'   \code{GenoMeta}; \code{select_method} (\code{"GA"}/\code{"GS"}/\code{"OCS"})
#'   maps to \code{scheme} (\code{GA -> OHS(GA_pairwise)}, \code{GS -> TS},
#'   \code{OCS -> OCS}).
#'
#' @return An \code{hs_simulation} result (class \code{c("hs_simulation",
#'   "hs_result", "list")}, see \code{\link{hs_result}}): the five 0.8.8
#'   elements below, plus \code{notes} (every silent fallback of the run, one
#'   row per job and code: \code{fallback_truncation}, \code{ocs_grm_failed},
#'   \code{ocs_kinship_cap}, \code{metric_failed}, \code{f1mean_failed},
#'   \code{monitor_write_failed}, \code{figure_failed}; one warning per code
#'   is also raised at the end of the run), \code{provenance} (package, status,
#'   the resolved \code{spec}, trait, seed, RNG kind, input sizes, block
#'   provenance, backend, the founder ideal's summary, the file tag, one row per
#'   job) and \code{figures} (empty: the trajectories figure is built by
#'   \code{write_outputs()} on demand). The result survives \code{saveRDS()}.
#' \describe{
#'   \item{\code{plans}}{Per-cross data.table: \code{CrossID, Parent1, Parent2,
#'     F1count, F1mean, select_type, cycle, rep, sel_scheme, D2Ideal,
#'     coancestry, ideal_phenotype, pair_score} (the last is the family's
#'     selection score under direct pairing, \code{NA} otherwise).}
#'   \item{\code{summary}}{Per-cycle data.table: \code{cycle, select_type, rep,
#'     sel_scheme, mean_F1mean, var_F1mean, max_F1mean, mean_D2Ideal,
#'     mean_coancestry, mean_gene_diversity, mean_ideal_phenotype,
#'     ideal_hap_gebv, marker_size, n_families, n_parents_used} (with a
#'     \code{cycle = 0} founder row). \code{max_F1mean} is the breeding value
#'     of the cycle's best single offspring (the top elite — the realised
#'     counterpart of a set-level "ultimate genotype" objective), against the
#'     population mean in \code{mean_F1mean}. \code{mean_gene_diversity} is
#'     Nei's (1973) gene diversity of the cycle's population (see
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
#'   \code{\link{gbs_selection}}, \code{\link{ohs_pair_selection}},
#'   \code{\link{gbs_pair_selection}}, \code{\link{pair_selection_budget}},
#'   \code{\link{construct_ideal_genotype}}, \code{\link{make_cross_design}},
#'   \code{\link{plot_simulation_trajectories}}
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
  avoid_selfing       = TRUE,
  selection_unit      = c("parents", "pairs"),
  pair_budget         = NULL,
  gbs_pair_level      = c("population", "family"),
  pair_freq           = NULL,
  pair_penalty        = c("loss", "hard"),
  pair_pool           = c("shortlist", "all"),
  shortlist_k         = NULL,
  max_use             = Inf,
  overshoot           = c("skip", "stop"),
  ts_pair_score       = c("midparent", "usefulness"),
  usefulness_k        = 1,
  scheme              = c("TS", "OHS", "GBS"),
  ohs_method          = c("GA_pairwise", "Greedy"),
  ohs_elite_ratio     = 1,
  ohs_elite_level     = c("pair", "block"),
  ohs_shortlist_families = NULL,
  gbs_method          = c("ARCSIN", "LOG"),
  gbs_preset          = c("haploStack", "kemper"),
  gbs_objective       = NULL,
  gbs_engine          = NULL,
  penalty_mode        = NULL,
  gbs_penalty         = NULL,
  penalty_scale       = NULL,
  penalty_rank_scale  = 1,
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
  output_dir          = NULL,
  verbose             = TRUE,
  monitor             = FALSE,
  n_cores             = 1L,
  trait               = NULL,
  keep_populations    = c("final", "none"),
  backend_opts        = list(),
  ga_control          = list(),
  progress            = NULL,
  # ── deprecated back-compat ───────────────────────────────────────────────
  catalogue_results   = NULL,
  trait_name          = NULL,
  select_method       = NULL
) {

  # ── Back-compat: catalogue_results / select_method ────────────────────────
  # The three inputs are on the deprecation ladder since 0.9.0 (sc/F31; they
  # never warned before): hs_warning_deprecated, mapped as before.
  if (is.null(catalogue_results) && !is.null(geno_meta) && !is_geno_meta(geno_meta)) {
    .hs_deprecated("a catalogue list passed as 'geno_meta'", "geno_meta = <GenoMeta>",
                   "0.9.0", what = "argument", fn = "run_stacking_simulation")
    catalogue_results <- geno_meta            # catalogue list passed positionally
  } else if (!is.null(catalogue_results)) {
    .hs_deprecated("catalogue_results", "geno_meta = <GenoMeta>", "0.9.0",
                   what = "argument", fn = "run_stacking_simulation")
  }
  if (!is.null(trait_name))
    .hs_deprecated("trait_name", "trait", "0.9.0", what = "argument",
                   fn = "run_stacking_simulation")

  if (!is.null(catalogue_results)) {
    trait_name <- trait_name %||% names(catalogue_results)[1L]
    if (!trait_name %in% names(catalogue_results))
      .hs_abort("run_stacking_simulation: trait '", trait_name,
           "' not in catalogue_results.")
    geno_meta <- catalogue_results[[trait_name]]$geno_meta
  }

  if (!is.null(select_method)) {              # GA/GS/OCS -> scheme
    .hs_deprecated("select_method", "scheme", "0.9.0", what = "argument",
                   fn = "run_stacking_simulation")
    select_method <- match.arg(select_method, c("GA", "GS", "OCS"))
    scheme <- switch(select_method, GA = "OHS", GS = "TS", OCS = "OCS")
    if (select_method == "GA") ohs_method <- "GA_pairwise"   # GA_coverage removed 0.7.4
  }

  .hs_check_meta(geno_meta, "run_stacking_simulation", "geno_meta",
                 hint = "(Or supply 'catalogue_results'.)")
  keep_populations <- match.arg(keep_populations)
  if (!is.null(output_dir) &&
      (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir)))
    .hs_abort("run_stacking_simulation: 'output_dir' must be NULL or a single ",
              "directory path.", call. = FALSE)

  # ── The spec: every science argument validated once (stk_spec, D4-3) ──────
  spec <- stk_spec(
    scheme = scheme, replications = replications, select_cycles = select_cycles,
    selfing_cycles = selfing_cycles, ga_start_by = ga_start_by,
    num_offspring = num_offspring, top_founders = top_founders,
    avoid_selfing = avoid_selfing, selection_unit = selection_unit,
    pair_budget = pair_budget, gbs_pair_level = gbs_pair_level,
    pair_freq = pair_freq, pair_penalty = pair_penalty, pair_pool = pair_pool,
    shortlist_k = shortlist_k, max_use = max_use, overshoot = overshoot,
    ts_pair_score = ts_pair_score, usefulness_k = usefulness_k,
    ohs_method = ohs_method, ohs_elite_ratio = ohs_elite_ratio,
    ohs_elite_level = ohs_elite_level,
    ohs_shortlist_families = ohs_shortlist_families, gbs_method = gbs_method,
    gbs_preset = gbs_preset, gbs_objective = gbs_objective,
    gbs_engine = gbs_engine, penalty_mode = penalty_mode,
    gbs_penalty = gbs_penalty, penalty_scale = penalty_scale,
    penalty_rank_scale = penalty_rank_scale,
    gbs_effect_weight = gbs_effect_weight, ideal_scope = ideal_scope,
    select_priority = select_priority, construct_method = construct_method,
    k_penalty = k_penalty, max_kinship = max_kinship,
    freq_penalty_weight = freq_penalty_weight,
    select_top_blocks = select_top_blocks, selected_blocks = selected_blocks,
    cross_backend = cross_backend, metrics = metrics, trait = trait,
    seed = seed, backend_opts = backend_opts, ga_control = ga_control)

  # What every job shares: the effect table, the founders' table, the founder
  # ideal, the opts list. Computed once, no random draw; its notes (row 19's
  # d2ideal_direction_blind) join the result's.
  prep_run <- .hs_with_notes(
    .stk_prepare(spec, geno_meta, effect_df, caller = "run_stacking_simulation"),
    muffle = TRUE)
  prep <- prep_run$value

  if (!is.null(output_dir) && !dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  jobs   <- stk_jobs(spec)
  n_jobs <- nrow(jobs)

  # The quick calculation of what `top_founders` implies per cycle.
  if (verbose) .hs_inform("Selection unit - ", .stk_budget_text(spec$budget))

  # Progress (X5): per generation when the jobs run here; per job when they
  # run in workers (events fire as each job's record is collected).
  use_parallel   <- .hs_parallel_ok(n_cores)
  cycle_progress <- if (use_parallel) NULL else progress
  job_trk <- .hs_progress_tracker(if (use_parallel) progress else NULL,
                                  "simulation", "job", n_jobs)

  # ── One (rep, scheme) job: init -> step until done -> record ──────────────
  run_job <- function(k) {
    st <- stk_init(spec, geno_meta, job = k, .prep = prep)
    on.exit(stk_close(st), add = TRUE)
    if (verbose)
      .hs_inform("Scheme ", st$scheme, " — replication ", st$rep, " / ",
                 spec$replications)
    while (identical(st$status, "running")) {
      st <- stk_step(st, cycle_progress)
      if (verbose) .stk_verbose_line(st)
    }
    if (verbose) {
      summ <- .stk_timing_summary(st$scheme, st$rep, st$step_time)
      if (!is.null(summ)) .hs_inform(summ)
    }
    .stk_job_record(st, announce = FALSE)     # the wrapper announces once per code
  }

  # ── Harvest: the live monitor, the job progress event, the cancel ─────────
  live_on <- isTRUE(monitor)
  if (live_on && is.null(output_dir)) {
    .hs_inform("run_stacking_simulation: 'monitor = TRUE' needs 'output_dir'; ",
            "live figures/tables are skipped.")
    live_on <- FALSE
  }
  live_tag  <- paste(unique(spec$scheme), collapse = "-")
  records   <- vector("list", n_jobs)
  n_done    <- 0L
  done_so_far <- function()
    stk_result(records[!vapply(records, is.null, logical(1L))], keep_populations)
  on_job <- function(k, rec) {
    records[[k]] <<- rec
    n_done <<- n_done + 1L
    if (live_on) {
      if (verbose)
        .hs_inform(sprintf("[%d/%d] %s rep %d done - mean_F1mean=%.4g",
                        n_done, n_jobs, rec$scheme, rec$rep,
                        mean(rec$summary[cycle > 0L]$mean_F1mean, na.rm = TRUE)))
      keep <- !vapply(records, is.null, logical(1L))
      .stk_render_live(data.table::rbindlist(lapply(records[keep], `[[`, "summary")),
                       data.table::rbindlist(lapply(records[keep], `[[`, "plans")),
                       output_dir, live_tag)
    }
    if (identical(rec$status, "cancelled"))
      .hs_abort("run_stacking_simulation: cancelled by the progress callback in ",
                rec$scheme, " rep ", rec$rep, " after generation ", rec$generation,
                " of ", rec$n_generations, ".", call. = FALSE,
                class = "hs_error_cancelled",
                data = list(partial = done_so_far(),
                            at = list(unit = "cycle", i = rec$generation,
                                      n = rec$n_generations, job = k)))
    job_trk$tick(k, label = rec$key,
                 payload = list(job = k, scheme = rec$scheme, rep = rec$rep,
                                status = rec$status, elapsed = rec$elapsed),
                 partial = done_so_far)
  }

  # The caller's stream is left as it was (every job seeds itself, row 3). The
  # jobs keep their own notes; whatever is raised outside them (the live
  # monitor) is collected here (X4).
  run <- .hs_with_notes(
    .hs_keep_rng(.stk_dispatch(n_jobs, run_job, if (use_parallel) n_cores else 1L,
                               on_job = on_job)),
    muffle = TRUE)

  res <- stk_result(records, keep_populations)
  res$provenance$fun <- "run_stacking_simulation"
  res$notes <- .stk_merge_notes(list(prep_run$notes, res$notes, run$notes))

  # ── Write output (row 10 a): only when a directory is named ───────────────
  if (!is.null(output_dir)) {
    w <- .hs_with_notes(write_outputs(res, output_dir, overwrite = TRUE), muffle = TRUE)
    res$notes <- .stk_merge_notes(list(res$notes, w$notes))
  }

  # One warning per note code for the whole run; the detail is hs_notes(res).
  .stk_announce_notes(res$notes)

  if (is.null(output_dir)) res else invisible(res)
}

# The 0.8.8 verbose progress line, once every 10 decision cycles (after that
# cycle's selfing generations): the cycle counter and the step consuming the
# most time so far, so a long run tells the user WHERE the time is going.
.stk_verbose_line <- function(st) {
  per <- st$spec$selfing_cycles + 1L
  if (st$generation %% per != 0L) return(invisible(NULL))
  cyc <- st$generation %/% per
  if (cyc %% 10L != 0L) return(invisible(NULL))
  dom <- .STK_STEP_LABELS[[names(which.max(st$step_time))]]
  .hs_inform(sprintf("  %s rep %d: cycle %d / %d — %.0fs elapsed, dominant step: %s",
                     st$scheme, st$rep, cyc, st$spec$select_cycles,
                     sum(st$step_time), dom))
  invisible(NULL)
}
