# =============================================================================
# stacking_step.R — the step API of the stacking simulation (0.9.0, design D4
# §3.1–3.2, form row 18 a)
#
#   stk_spec()    every science argument of run_stacking_simulation(), validated
#                 once and held once (the defaults live in this signature)
#   stk_jobs()    the (rep, scheme) grid of a spec, in the 0.8.8 order
#   stk_init()    one job's state: founders, founder record, the job's own
#                 random stream, the genomicSimulation session when that
#                 backend is used
#   stk_step()    advance the state by ONE generation under the job's stream;
#                 fires the progress event; never throws on cancel
#   stk_result()  one state or several -> the simulation result
#   stk_close()   free the gS session and the scratch; idempotent
#
# run_stacking_simulation() (stacking_main.R) is a wrapper over these five
# calls with identical numbers: a job's random draws happen in the same order
# as in 0.8.8 — set.seed(seed + rep) at the start, then per generation the
# selection draws and the crossing draws — only now the job's stream lives in
# the state (`rng`) and is swapped in for each step, so what happens between
# steps (a callback that draws, another job stepped in between, a save and a
# reload) cannot change a job's numbers, and the caller's own stream is left
# as it was.
# =============================================================================


# ---------------------------------------------------------------------------
# backend options (D4 §3.8; the map-unit and keep-files knobs land with D4-7)
# ---------------------------------------------------------------------------

.STK_BACKEND_OPTS <- list(cM_per_mb = 1, map_unit = "auto", keep_files = FALSE)

.stk_backend_opts <- function(x, caller = "stk_spec") {
  if (is.null(x)) x <- list()
  if (!is.list(x))
    .hs_abort(caller, ": 'backend_opts' must be a list.", call. = FALSE)
  bad <- setdiff(names(x), names(.STK_BACKEND_OPTS))
  if (length(bad))
    .hs_abort(caller, ": unknown backend_opts: ", paste(bad, collapse = ", "),
              " (known: ", paste(names(.STK_BACKEND_OPTS), collapse = ", "), ").",
              call. = FALSE)
  out <- utils::modifyList(.STK_BACKEND_OPTS, x)
  if (!is.numeric(out$cM_per_mb) || length(out$cM_per_mb) != 1L ||
      !is.finite(out$cM_per_mb) || out$cM_per_mb <= 0)
    .hs_abort(caller, ": backend_opts$cM_per_mb must be one positive number.",
              call. = FALSE)
  out$map_unit  <- match.arg(out$map_unit, c("auto", "bp", "cM"))
  out$keep_files <- isTRUE(out$keep_files)
  out
}


# ---------------------------------------------------------------------------
# stk_spec()
# ---------------------------------------------------------------------------

#' The specification of a stacking simulation (step API)
#'
#' The step API of \code{\link{run_stacking_simulation}} (0.9.0): a run is a
#' \emph{spec} (every science argument, validated once), a grid of
#' \emph{jobs} (\code{replications x scheme}), and per job a \emph{state} that
#' \code{\link{stk_step}} advances one generation at a time. The wrapper
#' \code{run_stacking_simulation()} is exactly
#' \preformatted{
#'   spec <- stk_spec(...)
#'   for (k in seq_len(nrow(stk_jobs(spec)))) {
#'     st <- stk_init(spec, geno_meta, effect_df, job = k)
#'     while (st$status == "running") st <- stk_step(st)
#'     stk_close(st)
#'   }
#'   stk_result(states)
#' }
#' and gives the same numbers: a job's random draws happen in the same order
#' as before, from \code{set.seed(seed + rep)} on, only the stream now lives in
#' the state and is swapped in for each step. So a website or a scheduler can
#' step a job, save it (\code{saveRDS(state)} works for the \code{"block"}
#' backend), resume it later, interleave jobs, or cancel one, and get the
#' numbers the one-call run gives.
#'
#' @inheritParams run_stacking_simulation
#' @param trait The trait whose effects the run breeds on when the container
#'   holds several (see \code{\link{run_stacking_simulation}}). Default
#'   \code{NULL}.
#' @param seed Base seed; replication \code{r} runs on \code{seed + r}.
#'   \code{NULL} draws one integer from the caller's stream and records it in
#'   the spec (\code{spec$seed}). Default \code{42}.
#' @param backend_opts Options of the crossing backend, a list with any of
#'   \code{cM_per_mb} (default \code{1}), \code{map_unit} (\code{"auto"},
#'   \code{"bp"}, \code{"cM"}) and \code{keep_files} (default \code{FALSE}:
#'   the genomicSimulation scratch files of a job are removed when the job is
#'   closed). Used by the genomicSimulation backend only.
#' @param ga_control A list with any of \code{pop_size}, \code{n_generations},
#'   \code{mutation_rate} for the GA engines (OHS \code{"GA_pairwise"} and GBS
#'   \code{gbs_engine = "GA"}). Default \code{list()}: the engines' defaults.
#'
#' @return \code{stk_spec()}: a list of class \code{"stk_spec"} holding every
#'   resolved argument plus \code{budget} (the per-cycle budget of
#'   \code{\link{pair_selection_budget}}), \code{gbs} (the resolved GBS
#'   specification), \code{n_generations} (\code{select_cycles *
#'   (selfing_cycles + 1)}) and the compact-ID formats. It holds no data, so it
#'   serialises and can be reused across containers.
#'
#' @seealso \code{\link{run_stacking_simulation}}, \code{\link{stk_step}}
#' @export
stk_spec <- function(
  scheme              = c("TS", "OHS", "GBS"),
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
  selected_blocks     = NULL,
  cross_backend       = c("block", "genomicSimulation"),
  metrics             = "gebv",
  trait               = NULL,
  seed                = 42L,
  backend_opts        = list(),
  ga_control          = list()
) {
  caller <- "stk_spec"
  scheme           <- match.arg(scheme, names(.STK_SCHEMES), several.ok = TRUE)
  ohs_method       <- match.arg(ohs_method)
  ohs_elite_level  <- match.arg(ohs_elite_level)
  .ohs_check_elite(ohs_elite_ratio)
  gbs_method       <- match.arg(gbs_method)
  gbs_preset       <- match.arg(gbs_preset)
  ideal_scope      <- match.arg(ideal_scope)
  select_priority  <- match.arg(select_priority)
  construct_method <- match.arg(construct_method)
  cross_backend    <- match.arg(cross_backend)
  selection_unit   <- match.arg(selection_unit)
  gbs_pair_level   <- match.arg(gbs_pair_level)
  pair_penalty     <- match.arg(pair_penalty)
  pair_pool        <- match.arg(pair_pool)
  overshoot        <- match.arg(overshoot)
  ts_pair_score    <- match.arg(ts_pair_score)
  if (!is.null(pair_freq)) pair_freq <- match.arg(pair_freq, c("projected", "running"))

  for (nm in c("replications", "select_cycles", "selfing_cycles", "ga_start_by",
               "num_offspring", "top_founders", "select_top_blocks")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0 || v != round(v))
      .hs_abort(caller, ": '", nm, "' must be one whole number >= 0.", call. = FALSE)
    assign(nm, as.integer(v))
  }
  if (replications < 1L || select_cycles < 1L || num_offspring < 1L ||
      top_founders < 2L || ga_start_by < 1L)
    .hs_abort(caller, ": 'replications', 'select_cycles', 'num_offspring' must be ",
              ">= 1, 'ga_start_by' >= 1 and 'top_founders' >= 2.", call. = FALSE)
  if (!is.null(trait) && (!is.character(trait) || length(trait) != 1L || is.na(trait)))
    .hs_abort(caller, ": 'trait' must be NULL or one trait name.", call. = FALSE)
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  } else if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    .hs_abort(caller, ": 'seed' must be one number or NULL.", call. = FALSE)
  }
  if (!is.list(ga_control) ||
      length(setdiff(names(ga_control), c("pop_size", "n_generations", "mutation_rate"))))
    .hs_abort(caller, ": 'ga_control' must be a list with any of pop_size, ",
              "n_generations, mutation_rate.", call. = FALSE)
  backend_opts <- .stk_backend_opts(backend_opts, caller)

  # Row 19 a (intake 07): select_priority = "low" is refused where a scheme
  # cannot honour it — OHS and OCS select for the high tail only, and GBS's
  # "custom" / "ratio" ideal builders target the high states whatever the
  # priority says (only the argmax builder takes the direction, N2).
  if (select_priority == "low") {
    dir <- vapply(scheme, function(s) .stk_scheme(s)$direction, character(1L))
    bad <- scheme[dir == "high" | (dir == "argmax" & construct_method != "argmax")]
    if (length(bad))
      .hs_abort(caller, ": select_priority = \"low\" is not supported by scheme ",
                paste(sQuote(bad, FALSE), collapse = ", "),
                if (any(dir[scheme %in% bad] == "argmax"))
                  " (GBS honours the direction only with construct_method = \"argmax\")"
                else "",
                "; TS does, in both selection units.", call. = FALSE)
  }

  # Selection unit and its per-cycle budget (parents-first: the traversal of
  # `top_founders`; direct pairing: the chosen `pair_budget`). Resolved up front
  # so a missing `pair_budget` fails before anything runs, and so the family-ID
  # width follows the family count. Selfing: the all-pairs design never selfs,
  # so parents-first passes FALSE; under pairs the flag is honoured by the
  # candidate pool, the budget and the design alike.
  if (selection_unit == "parents" && !is.null(pair_budget))
    .hs_warn(caller, ": 'pair_budget' is ignored unless selection_unit = \"pairs\".",
             call. = FALSE)
  budget <- .stk_selection_budget(
    selection_unit, pair_budget, top_founders,
    allow_selfing = if (selection_unit == "pairs") !isTRUE(avoid_selfing) else FALSE)
  if (selection_unit == "pairs") {
    no_pairs <- scheme[vapply(scheme, function(s) is.null(.stk_scheme(s)$pairs),
                              logical(1L))]
    if (length(no_pairs))
      .hs_abort(caller, ": scheme \"", no_pairs[1L], "\" has no direct-pair form yet ",
                "(pair plan Phase C); run it under selection_unit = \"parents\".",
                call. = FALSE)
    if (!is.null(gbs_engine) && gbs_engine == "GA")
      .hs_abort(caller, ": gbs_engine = \"GA\" over family sets is not implemented ",
                "(pair plan Phase D); the pairs engine is the greedy.", call. = FALSE)
    # The pairs default objective is "set" (the candidate score is separable
    # over a family); an explicit value still wins.
    if (is.null(gbs_objective)) gbs_objective <- "set"
  }

  # GBS knobs: the preset fills whatever was left NULL; explicit values win
  # (one resolver for the driver and the two GBS entry points, sc/F25).
  gbs <- .gbs_spec(gbs_preset, gbs_objective = gbs_objective,
                   gbs_engine = gbs_engine, penalty_mode = penalty_mode,
                   gbs_penalty = gbs_penalty, penalty_scale = penalty_scale,
                   gbs_effect_weight = gbs_effect_weight,
                   penalty_rank_scale = penalty_rank_scale, caller = caller)

  # Metrics (measurement-only; never affect selection). "gebv" (F1mean) is the
  # always-on base; "all" = every one.
  all_metrics <- c("gebv", "d2ideal", "gene_diversity", "coancestry",
                   "ideal_hap_gebv")
  if (any(metrics == "all")) metrics <- all_metrics
  metrics <- match.arg(metrics, all_metrics, several.ok = TRUE)
  metrics <- union("gebv", metrics)

  n_generations <- select_cycles * (selfing_cycles + 1L)

  # Fixed-width fields for the compact generational IDs (C<cycle>F<family>_<child>):
  # widths are sized to the run so IDs stay short, uniform and sortable. The
  # family width follows the most families a cycle can make under the budget
  # (parents-first: choose(top_founders, 2), unchanged).
  fam_max <- max(2L, budget$n_families_max)
  cyc_w   <- max(2L, nchar(as.character(n_generations)))
  fam_w   <- max(2L, nchar(as.character(fam_max)))
  ch_w    <- max(1L, nchar(as.character(num_offspring)))

  structure(list(
    scheme = scheme, replications = replications, select_cycles = select_cycles,
    selfing_cycles = selfing_cycles, ga_start_by = ga_start_by,
    num_offspring = num_offspring, top_founders = top_founders,
    avoid_selfing = isTRUE(avoid_selfing), selection_unit = selection_unit,
    pair_budget = budget$pair_budget, budget = budget,
    gbs_pair_level = gbs_pair_level, pair_freq = pair_freq,
    pair_penalty = pair_penalty, pair_pool = pair_pool, shortlist_k = shortlist_k,
    max_use = max_use, overshoot = overshoot, ts_pair_score = ts_pair_score,
    usefulness_k = usefulness_k, ohs_method = ohs_method,
    ohs_elite_ratio = ohs_elite_ratio, ohs_elite_level = ohs_elite_level,
    ohs_shortlist_families = ohs_shortlist_families, gbs_method = gbs_method,
    gbs_preset = gbs_preset, gbs = gbs, penalty_rank_scale = penalty_rank_scale,
    ideal_scope = ideal_scope, select_priority = select_priority,
    construct_method = construct_method, k_penalty = k_penalty,
    max_kinship = max_kinship, freq_penalty_weight = freq_penalty_weight,
    select_top_blocks = select_top_blocks, selected_blocks = selected_blocks,
    cross_backend = cross_backend, metrics = metrics, trait = trait,
    seed = seed, backend_opts = backend_opts, ga_control = ga_control,
    n_generations = n_generations,
    id_fmt_sample = paste0("C%0", cyc_w, "dF%0", fam_w, "d_%0", ch_w, "d"),
    id_fmt_cross  = paste0("C%0", cyc_w, "dF%0", fam_w, "d")),
    class = "stk_spec")
}

#' @export
print.stk_spec <- function(x, ...) {
  cat("<stk_spec> ", paste(x$scheme, collapse = "/"), " x ", x$replications,
      " rep(s) x ", x$select_cycles, " cycle(s)",
      if (x$selfing_cycles > 0L) paste0(" (+", x$selfing_cycles, " selfing)") else "",
      " = ", x$n_generations, " generation(s); ", x$selection_unit, "; ",
      x$cross_backend, "; seed ", x$seed, "\n", sep = "")
  invisible(x)
}

.stk_check_spec <- function(spec, fn) {
  if (!inherits(spec, "stk_spec"))
    .hs_abort(fn, ": 'spec' must come from stk_spec().", call. = FALSE)
  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# stk_jobs()
# ---------------------------------------------------------------------------

#' @rdname stk_spec
#' @param spec An \code{stk_spec}.
#' @return \code{stk_jobs()}: a data.frame \code{job, rep, scheme} — one row
#'   per (replication, scheme), replications varying fastest (the 0.8.8 grid).
#' @export
stk_jobs <- function(spec) {
  .stk_check_spec(spec, "stk_jobs")
  g <- expand.grid(rep = seq_len(spec$replications), scheme = spec$scheme,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  data.frame(job = seq_len(nrow(g)), rep = as.integer(g$rep), scheme = g$scheme,
             stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------------
# Internal — what every job of a run shares (computed once, D4 §3.2)
# ---------------------------------------------------------------------------

# The resolved effect table, the founders' haplotype-level table, the founder
# ideal (when a scheme selects against it or d2ideal is measured), and the
# `opts` list the cycle helpers read. No random draw happens here.
.stk_prepare <- function(spec, geno_meta, effect_df = NULL, caller = "stk_init") {
  .hs_check_meta(geno_meta, caller, "geno_meta")
  effect_df <- .stk_resolve_effects(effect_df, geno_meta, spec$trait, caller)

  # Row 19 a: a parent set larger than the population was capped silently.
  n_pool <- dim(geno_meta$geno_arrays)[1L]
  if (spec$top_founders > n_pool)
    .hs_abort(caller, ": top_founders (", spec$top_founders, ") exceeds the ",
              "population (", n_pool, " samples).", call. = FALSE,
              data = list(bad = spec$top_founders))
  # The d2ideal metric measures distance to the ideal built by the "custom"
  # or "ratio" builder, which targets the high states: under
  # select_priority = "low" with a direction-blind scheme (TS) the metric
  # still measures distance to the HIGH ideal. Said once per run.
  if (spec$select_priority == "low" && "d2ideal" %in% spec$metrics &&
      spec$construct_method != "argmax")
    .hs_note("prepare", "d2ideal_direction_blind",
             paste0("d2ideal is measured against the ideal of construct_method = \"",
                    spec$construct_method, "\", which targets the high states; ",
                    "under select_priority = \"low\" it is a distance to the ",
                    "high ideal (use construct_method = \"argmax\" for a low ideal)."))

  founder_lg <- .stk_limit_blocks(
    compute_profile_localgebv(geno_meta, effect_df, level = "haplotype"),
    spec$selected_blocks, spec$select_top_blocks)
  if (nrow(founder_lg) == 0L)
    .hs_abort(caller, ": no blocks remain for selection (check 'selected_blocks' ",
              "/ block assignments).", call. = FALSE)

  # Founder ideal-genotype target returned to the caller: only meaningful when a
  # scheme selects against it or d2ideal is measured; else NULL (metrics-only).
  ideal <- if (any(.stk_scheme_flag(spec$scheme, "needs_ideal")) ||
               "d2ideal" %in% spec$metrics)
    construct_ideal_genotype(founder_lg, method = spec$construct_method,
                             k_penalty = spec$k_penalty,
                             select_priority = spec$select_priority)
  else NULL

  gbs <- spec$gbs
  opts <- list(
    selected_blocks   = spec$selected_blocks,
    select_top_blocks = spec$select_top_blocks,
    top_founders      = spec$top_founders,
    num_offspring     = spec$num_offspring,
    avoid_selfing     = spec$avoid_selfing,
    selection_unit    = spec$selection_unit,
    pair_budget       = spec$budget,       # resolved budget (see .stk_selection_budget)
    gbs_pair_level    = spec$gbs_pair_level,
    pair_freq         = spec$pair_freq,
    pair_penalty      = spec$pair_penalty,
    pair_pool         = spec$pair_pool,
    shortlist_k       = spec$shortlist_k,
    max_use           = spec$max_use,
    overshoot         = spec$overshoot,
    ts_pair_score     = spec$ts_pair_score,
    usefulness_k      = spec$usefulness_k,
    ohs_method        = spec$ohs_method,
    ohs_elite_ratio   = spec$ohs_elite_ratio,
    ohs_elite_level   = spec$ohs_elite_level,
    ohs_shortlist_families = spec$ohs_shortlist_families,
    gbs_method        = spec$gbs_method,
    gbs_objective     = gbs$gbs_objective,
    gbs_engine        = gbs$gbs_engine,
    penalty_mode      = gbs$penalty_mode,
    gbs_penalty       = gbs$gbs_penalty,
    penalty_scale     = gbs$penalty_scale,
    penalty_rank_scale = spec$penalty_rank_scale,
    gbs_effect_weight = gbs$gbs_effect_weight,
    ideal_scope       = spec$ideal_scope,
    founder_ideal     = NULL,
    founder_blocks    = NULL,
    select_priority   = spec$select_priority,
    construct_method  = spec$construct_method,
    k_penalty         = spec$k_penalty,
    max_kinship       = spec$max_kinship,
    freq_penalty_weight = spec$freq_penalty_weight,
    cross_backend     = spec$cross_backend,
    metrics           = spec$metrics,
    ploidy            = geno_meta$metadata$ploidy %||% 2L,
    seed              = spec$seed,
    ga_control        = spec$ga_control,
    work_dir          = NULL,              # per job (stk_init)
    id_fmt_sample     = spec$id_fmt_sample,
    id_fmt_cross      = spec$id_fmt_cross)
  # Founder-fixed scope: every cycle selects and measures against this table,
  # and GBS selection is pinned to the founders' block set (see .stk_advance).
  if (spec$ideal_scope == "founder" && !is.null(ideal)) {
    opts$founder_ideal  <- ideal
    opts$founder_blocks <- unique(as.character(founder_lg$Block))
  }
  list(effect_df = effect_df, founder_lg = founder_lg, ideal = ideal, opts = opts)
}

# The one effect table a run breeds on (row 19 a, io/F15, N3). `effect_df`
# (or the container's marker_effects) may be one table, a named per-trait
# list, or one long table with a `Trait` column (.as_effect_list() reads all
# three). One trait is taken as it is; several need `trait`, which must name
# one of them — up to 0.8.8 the first trait was taken silently, and a `trait`
# that matched nothing fell back to the first table too. `allow_none = TRUE`
# (the genomicSimulation loader) returns NULL when there are no effects.
.stk_resolve_effects <- function(effect_df, geno_meta, trait = NULL,
                                 caller = "run_stacking_simulation",
                                 allow_none = FALSE) {
  if (is.null(effect_df)) effect_df <- geno_meta$marker_effects
  tabs <- .as_effect_list(effect_df)
  if (is.null(tabs) || !length(tabs)) {
    if (isTRUE(allow_none)) return(NULL)
    .hs_abort(caller, ": no effects supplied — pass 'effect_df' or ",
              "attach 'marker_effects' to the GenoMeta.", call. = FALSE)
  }
  if (!is.null(trait)) {
    if (!trait %in% names(tabs))
      .hs_abort(caller, ": trait '", trait, "' is not among the effects (",
                paste(sQuote(names(tabs), FALSE), collapse = ", "), ").",
                call. = FALSE, data = list(bad = trait))
    return(tabs[[trait]])
  }
  if (length(tabs) > 1L)
    .hs_abort(caller, ": the effects hold ", length(tabs), " traits (",
              paste(sQuote(names(tabs), FALSE), collapse = ", "),
              "); this run breeds on one — pass trait = \"<name>\" (or one ",
              "effect table, e.g. effect_df = meta$marker_effects[[\"", names(tabs)[1L],
              "\"]]).", call. = FALSE, data = list(bad = names(tabs)))
  tabs[[1L]]
}


# ---------------------------------------------------------------------------
# stk_init()
# ---------------------------------------------------------------------------

#' Initialise, step, close and collect one simulation job (step API)
#'
#' See \code{\link{stk_spec}} for the whole picture. \code{stk_init()} builds
#' the state of one \code{(rep, scheme)} job: the founders and their cycle-0
#' record, the job's random stream (\code{set.seed(seed + rep)}, captured),
#' and — under \code{cross_backend = "genomicSimulation"} — the live session,
#' opened on the founders (one per R process, see
#' \code{\link{gs_load_geno_meta}}). \code{stk_step()} advances it by one
#' generation: the decision generation of a cycle runs the job's scheme (or
#' truncation before \code{ga_start_by}), the \code{selfing_cycles}
#' generations after it run truncation (\code{select_type = "selfing"}).
#'
#' @section The state:
#' A plain list (class \code{"stk_state"}): \code{spec, job, rep, scheme,
#' generation, n_generations, current} (the current population, a
#' \code{GenoMeta}), \code{lg} (its haplotype-level local GEBV table, the next
#' generation's parent table), \code{rng} (the job's stream), \code{session}
#' (the genomicSimulation handle or \code{NULL}), \code{plan_recs},
#' \code{summ_recs} (the per-generation records), \code{notes}, \code{status}
#' (\code{"running"}, \code{"done"}, \code{"cancelled"}). A \code{"block"}
#' state survives \code{saveRDS()} / \code{readRDS()} and resumes identically;
#' a genomicSimulation state holds a live C-level group and cannot be resumed
#' from a file (\code{stk_step()} on one stops with \code{hs_error_state}).
#'
#' @section Progress and cancel:
#' \code{progress} is the package's one callback contract: after each
#' generation it receives \code{list(stage = "simulation", unit = "cycle", i,
#' n, label, elapsed, payload)} where \code{payload} holds \code{job, scheme,
#' rep, generation, cycle, select_type, summary} (that generation's summary
#' row), \code{elapsed} (the job's wall time) and \code{timings}. Returning
#' \code{FALSE} cancels: the state comes back with \code{status =
#' "cancelled"} and every generation so far on record — the step API never
#' throws on a cancel; the wrapper \code{run_stacking_simulation()} turns it
#' into an \code{hs_error_cancelled} error carrying the partial result.
#'
#' @param spec An \code{\link{stk_spec}}.
#' @param geno_meta Founder \code{GenoMeta}.
#' @param effect_df Marker effects (one table, or a per-trait list with
#'   \code{spec$trait} naming the one to use); \code{NULL} takes the
#'   container's \code{marker_effects}.
#' @param job Which row of \code{\link{stk_jobs}}\code{(spec)} to run. Default
#'   \code{1}.
#' @param state An \code{stk_state} from \code{stk_init()} / \code{stk_step()}.
#' @param progress \code{NULL} (default) or a function of one argument (the
#'   event; see the section above).
#' @param .prep Internal: what \code{run_stacking_simulation()} computes once
#'   for all jobs. Leave \code{NULL}.
#'
#' @return \code{stk_init()} and \code{stk_step()}: the state. \code{stk_close()}:
#'   the state, invisibly, with its session freed (idempotent; a state can be
#'   closed any number of times). \code{stk_result()}: see
#'   \code{\link{run_stacking_simulation}}.
#'
#' @seealso \code{\link{stk_spec}}, \code{\link{run_stacking_simulation}}
#' @name stk_state
#' @export
stk_init <- function(spec, geno_meta, effect_df = NULL, job = 1L, .prep = NULL) {
  .stk_check_spec(spec, "stk_init")
  jobs <- stk_jobs(spec)
  if (!is.numeric(job) || length(job) != 1L || is.na(job) || job < 1L ||
      job > nrow(jobs) || job != round(job))
    .hs_abort("stk_init: 'job' must be a row number of stk_jobs(spec) (1..",
              nrow(jobs), ").", call. = FALSE)
  job   <- as.integer(job)
  prep  <- if (is.null(.prep)) .stk_prepare(spec, geno_meta, effect_df) else .prep
  rep   <- jobs$rep[job]
  sch   <- jobs$scheme[job]
  opts  <- prep$opts

  # The job's own stream, from exactly where 0.8.8's loop began:
  # set.seed(seed + rep), captured; the caller's stream is put back.
  rng <- .with_rng_state(spec$seed + rep, NULL)$state

  # genomicSimulation: the session is opened here, on the founders, so the
  # state always holds the live handle (sc/F21: no cycle can leak one). Its
  # scratch is private to the job unless backend_opts$keep_files keeps it.
  session  <- NULL
  work_dir <- NULL
  if (identical(spec$cross_backend, "genomicSimulation")) {
    bo <- spec$backend_opts
    if (isTRUE(bo$keep_files)) {
      work_dir <- tempfile("hs_gs_")
      dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    }
    session  <- gs_load_geno_meta(geno_meta, prep$effect_df, work_dir = work_dir,
                                  prefix = "gs", cM_per_mb = bo$cM_per_mb,
                                  map_unit = bo$map_unit)
    work_dir <- session$work_dir
  }
  opts$work_dir <- work_dir

  n_gen <- spec$n_generations
  plan_recs <- vector("list", n_gen)
  summ_recs <- vector("list", n_gen + 1L)
  summ_recs[[1L]] <- .stk_summary_record(
    .stk_founder_measures(geno_meta, prep$founder_lg, prep$ideal, opts),
    "founder", 0L, rep, sch)

  d <- dim(geno_meta$geno_arrays)
  structure(list(
    spec = spec, job = job, rep = rep, scheme = sch,
    generation = 0L, n_generations = n_gen,
    current = geno_meta, lg = NULL, rng = rng, session = session,
    work_dir = work_dir, effect_df = prep$effect_df, opts = opts,
    ideal = prep$ideal,
    inputs = list(n_samples = d[1L], n_markers = d[2L],
                  n_blocks = length(unique(as.character(
                    geno_meta$block_marker_map_info$Block)))),
    blocks = geno_meta$metadata$blocks,
    plan_recs = plan_recs, summ_recs = summ_recs,
    step_time = c(parent_gebv = 0, ideal = 0, selection = 0, crossing = 0,
                  offspring_gebv = 0, metrics = 0),
    notes = .stk_notes_df(), status = "running",
    started = as.numeric(Sys.time()), elapsed = 0),
    class = "stk_state")
}

.stk_check_state <- function(state, fn) {
  if (!inherits(state, "stk_state"))
    .hs_abort(fn, ": 'state' must come from stk_init() / stk_step().", call. = FALSE)
  invisible(TRUE)
}

# The notes table of a job: the framework columns plus job and first_cycle.
.stk_notes_df <- function() {
  cbind(.hs_notes_df(), job = integer(), first_cycle = integer())
}

# Fold the notes one generation raised into the job's table: same (stage, code,
# message) -> n summed, first_cycle kept.
.stk_add_notes <- function(state, notes, generation) {
  if (!NROW(notes)) return(state)
  notes <- notes[, c("stage", "code", "message", "n"), drop = FALSE]
  notes$job <- state$job
  notes$first_cycle <- as.integer(generation)
  cur <- state$notes
  key_new <- paste(notes$stage, notes$code, notes$message, sep = "\r")
  key_cur <- paste(cur$stage, cur$code, cur$message, sep = "\r")
  hit <- match(key_new, key_cur)
  if (any(!is.na(hit)))
    cur$n[hit[!is.na(hit)]] <- cur$n[hit[!is.na(hit)]] + notes$n[!is.na(hit)]
  cur <- rbind(cur, notes[is.na(hit), , drop = FALSE])
  rownames(cur) <- NULL
  state$notes <- cur
  state
}

# The schedule is a pure function of the generation number (D4 §3.1): `g` is a
# decision generation when (g - 1) %% (selfing_cycles + 1) == 0, with decision
# index cyc; the scheme used is the baseline on selfing generations and before
# ga_start_by, else the job's scheme.
.stk_schedule <- function(spec, scheme, g) {
  per <- spec$selfing_cycles + 1L
  decision <- (g - 1L) %% per == 0L
  cyc      <- (g - 1L) %/% per + 1L
  used <- if (!decision || cyc < spec$ga_start_by) .STK_BASELINE else scheme
  list(cycle = cyc, decision = decision, scheme_used = used,
       select_type = if (decision) used else "selfing")
}


# ---------------------------------------------------------------------------
# stk_step()
# ---------------------------------------------------------------------------

#' @rdname stk_state
#' @export
stk_step <- function(state, progress = NULL) {
  .stk_check_state(state, "stk_step")
  if (!identical(state$status, "running")) return(state)
  trk <- .hs_progress_tracker(progress, "simulation", "cycle", state$n_generations)
  if (!is.null(state$session) && !.gs_handle_live(state$session))
    .hs_abort(class = "hs_error_state",
              "stk_step: a genomicSimulation session cannot be resumed from a saved ",
              "state (it holds a live C-level group); rerun the job.", call. = FALSE)

  g  <- state$generation + 1L
  sc <- .stk_schedule(state$spec, state$scheme, g)

  # The generation runs on the job's own stream; every haploStack warning it
  # raises becomes a note of the job (X4) and is muffled here — the wrapper
  # re-raises one per code at the end of the job.
  run <- .hs_with_notes(
    .with_rng_state(state$rng,
      .stk_advance(state$current, sc$scheme_used, state$effect_df, state$opts, g,
                   session = state$session, lg = state$lg)),
    muffle = TRUE)
  adv       <- run$value$value
  state$rng <- run$value$state
  state     <- .stk_add_notes(state, run$notes, g)

  state$plan_recs[[g]]      <- .stk_plans_record(adv, sc$select_type, g, state$rep, state$scheme)
  state$summ_recs[[g + 1L]] <- .stk_summary_record(adv, sc$select_type, g, state$rep, state$scheme)
  state$current    <- adv$off
  state$lg         <- adv$off_lg
  state$session    <- adv$session
  state$step_time  <- state$step_time + adv$timings
  state$generation <- g
  state$elapsed    <- as.numeric(Sys.time()) - state$started
  if (g >= state$n_generations) {
    state$status <- "done"
    state <- stk_close(state)          # nothing more to cross: free the session
  }

  # The progress event fires outside the job's stream. FALSE = cancel: the
  # state stays valid, with this generation on record, and never throws here.
  if (!is.null(progress)) {
    payload <- list(job = state$job, scheme = state$scheme, rep = state$rep,
                    generation = g, n_generations = state$n_generations,
                    cycle = sc$cycle, select_type = sc$select_type,
                    summary = state$summ_recs[[g + 1L]],
                    elapsed = state$elapsed, timings = adv$timings)
    cancelled <- tryCatch({
      trk$tick(g, label = sprintf("%s rep %d", state$scheme, state$rep),
               payload = payload)
      FALSE
    }, hs_error_cancelled = function(e) TRUE)
    if (cancelled && identical(state$status, "running")) {
      state$status <- "cancelled"
      state <- stk_close(state)
    }
  }
  state
}

# Is this handle the live genomicSimulation session? (a deserialised state's
# handle is not: the token no longer matches, or no session is open)
.gs_handle_live <- function(handle) {
  inherits(handle, "gs_handle") && isTRUE(.hs_gs$open) &&
    identical(as.integer(handle$token), .hs_gs$token)
}


# ---------------------------------------------------------------------------
# stk_close()
# ---------------------------------------------------------------------------

#' @rdname stk_state
#' @export
stk_close <- function(state) {
  .stk_check_state(state, "stk_close")
  if (!is.null(state$session)) {
    gs_close(state$session)                 # frees the group and its private dir
    state$session <- NULL
  }
  if (!is.null(state$work_dir) && !isTRUE(state$spec$backend_opts$keep_files) &&
      dir.exists(state$work_dir))
    unlink(state$work_dir, recursive = TRUE, force = TRUE)
  invisible(state)
}


# ---------------------------------------------------------------------------
# stk_result()
# ---------------------------------------------------------------------------

# The serialisable record of a finished (or cancelled) job: the tables, the
# final population, the notes, no live handle. With `announce = TRUE` (the
# step API) one hs_warning_fallback per note code is raised here ("n times"),
# so a caller hears what the job muffled; the wrapper announces once per code
# for the whole run instead (X4).
.stk_job_record <- function(state, announce = TRUE) {
  .stk_check_state(state, ".stk_job_record")
  keep <- !vapply(state$plan_recs, is.null, logical(1L))
  notes <- state$notes
  if (isTRUE(announce) && NROW(notes)) {
    for (i in seq_len(nrow(notes)))
      .hs_note(notes$stage[i], notes$code[i],
               paste0(notes$message[i],
                      if (notes$n[i] > 1L) sprintf(" (%d times)", notes$n[i]) else "",
                      " [", state$scheme, " rep ", state$rep, ", first in generation ",
                      notes$first_cycle[i], "]"),
               job = state$job, first_cycle = notes$first_cycle[i])
  }
  list(job = state$job, rep = state$rep, scheme = state$scheme,
       key = paste0(state$scheme, "_rep", state$rep),
       plans = rbindlist(state$plan_recs[keep]),
       summary = rbindlist(state$summ_recs[!vapply(state$summ_recs, is.null, logical(1L))]),
       final = state$current, ideal = state$ideal, notes = notes,
       status = state$status, generation = state$generation,
       n_generations = state$n_generations, step_time = state$step_time,
       elapsed = state$elapsed, work_dir = state$work_dir,
       kept_files = isTRUE(state$spec$backend_opts$keep_files) && !is.null(state$work_dir),
       spec = state$spec, inputs = state$inputs, blocks = state$blocks)
}

.stk_is_record <- function(x) is.list(x) && !is.null(x$key) && !is.null(x$summary)

# Several notes tables (per job, plus the wrapper's own) -> one, with the
# framework columns first; rows ordered by job then first generation.
.stk_merge_notes <- function(tables) {
  tables <- Filter(function(t) !is.null(t) && NROW(t) > 0L, tables)
  if (!length(tables)) return(.stk_notes_df())
  cols <- c("stage", "code", "message", "n", "job", "first_cycle")
  tables <- lapply(tables, function(t) {
    t <- as.data.frame(t, stringsAsFactors = FALSE)
    for (cn in cols) if (!cn %in% names(t)) t[[cn]] <- NA
    t$n <- as.integer(t$n); t$job <- as.integer(t$job)
    t$first_cycle <- as.integer(t$first_cycle)
    extra <- setdiff(names(t), cols)
    t[, c(cols, extra), drop = FALSE]
  })
  all_cols <- unique(unlist(lapply(tables, names)))
  tables <- lapply(tables, function(t) { t[setdiff(all_cols, names(t))] <- NA; t[all_cols] })
  out <- do.call(rbind, tables)
  out <- out[order(out$job, out$first_cycle, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# One hs_warning_fallback per distinct note code for a whole run.
.stk_announce_notes <- function(notes) {
  if (!NROW(notes)) return(invisible(NULL))
  for (code in unique(notes$code)) {
    rows <- notes[notes$code == code, , drop = FALSE]
    n_tot <- sum(rows$n, na.rm = TRUE)
    n_job <- length(unique(rows$job[!is.na(rows$job)]))
    .hs_warn(paste0(code, ": ", rows$message[1L],
                    if (n_tot > 1L) sprintf(" (%d times%s; see hs_notes())", n_tot,
                                            if (n_job > 1L) sprintf(" in %d jobs", n_job) else "")
                    else ""),
             call. = FALSE, class = "hs_warning_fallback",
             data = list(stage = rows$stage[1L], code = code, n = n_tot))
  }
  invisible(NULL)
}

#' @rdname stk_state
#' @param states One \code{stk_state}, or a list of states (or of the job
#'   records \code{run_stacking_simulation()} collects).
#' @param keep_populations \code{"final"} (default): the final population of
#'   every job in \code{final_populations}; \code{"none"}: an empty list (the
#'   largest object of a result; a website's choice).
#' @export
stk_result <- function(states, keep_populations = c("final", "none")) {
  keep_populations <- match.arg(keep_populations)
  if (inherits(states, "stk_state")) states <- list(states)
  if (!is.list(states) || !length(states))
    .hs_abort("stk_result: 'states' must be an stk_state or a list of them.",
              call. = FALSE)
  recs <- lapply(states, function(s) {
    if (inherits(s, "stk_state")) .stk_job_record(s, announce = TRUE)
    else if (.stk_is_record(s)) s
    else .hs_abort("stk_result: every element must be an stk_state.", call. = FALSE)
  })
  spec    <- recs[[1L]]$spec
  plans   <- rbindlist(lapply(recs, `[[`, "plans"))
  summary <- rbindlist(lapply(recs, `[[`, "summary"))
  keys    <- vapply(recs, `[[`, character(1L), "key")
  final_populations <- if (keep_populations == "final")
    stats::setNames(lapply(recs, `[[`, "final"), keys)
  else stats::setNames(list(), character())
  ideal  <- recs[[1L]]$ideal
  status <- vapply(recs, `[[`, character(1L), "status")
  tag    <- paste(unique(spec$scheme), collapse = "-")
  gs_ver <- if (identical(spec$cross_backend, "genomicSimulation") &&
                requireNamespace("genomicSimulation", quietly = TRUE))
    as.character(utils::packageVersion("genomicSimulation")) else NA_character_
  kept <- unlist(lapply(recs, function(r) if (isTRUE(r$kept_files)) r$work_dir))
  prov <- .hs_provenance(
    "stk_result",
    status   = if (any(status == "cancelled")) "cancelled" else "done",
    spec     = unclass(spec),
    trait    = spec$trait,
    seed     = spec$seed,
    rng_kind = RNGkind()[1L],
    inputs   = recs[[1L]]$inputs,
    blocks   = recs[[1L]]$blocks,
    backend  = list(name = spec$cross_backend, version = gs_ver),
    ideal    = if (is.null(ideal)) NULL else
      list(ideal_phenotype = as.numeric(attr(ideal, "ideal_phenotype")),
           ploidy = as.integer(attr(ideal, "ploidy")),
           method = spec$construct_method, k_penalty = spec$k_penalty),
    tag      = tag,
    jobs     = data.frame(job = vapply(recs, `[[`, integer(1L), "job"),
                          rep = vapply(recs, `[[`, integer(1L), "rep"),
                          scheme = vapply(recs, `[[`, character(1L), "scheme"),
                          status = status,
                          generation = vapply(recs, `[[`, integer(1L), "generation"),
                          elapsed = vapply(recs, `[[`, numeric(1L), "elapsed"),
                          stringsAsFactors = FALSE),
    kept_files = if (length(kept)) as.character(kept) else character())
  .hs_result(
    list(plans             = plans,
         summary           = summary,
         ideal_genotypes   = ideal,
         final_populations = final_populations,
         trajectories      = .stk_trajectories_view(summary)),
    "hs_simulation",
    notes = .stk_merge_notes(lapply(recs, `[[`, "notes")),
    provenance = prov)
}

#' @export
print.hs_simulation <- function(x, ...) {
  prov <- hs_provenance(x)
  s <- as.data.table(x$summary)
  cat("<haploStack result: hs_simulation>\n")
  if (nrow(s)) {
    last <- max(s$cycle)
    cat("  schemes: ", paste(unique(s$sel_scheme), collapse = ", "), " x ",
        length(unique(s$rep)), " replication(s) x ", last, " generation(s)\n", sep = "")
    fin <- s[cycle == last, .(m = mean(mean_F1mean, na.rm = TRUE)), by = sel_scheme]
    for (i in seq_len(nrow(fin)))
      cat(sprintf("  %-6s final mean F1mean %s\n", fin$sel_scheme[i], format(fin$m[i], digits = 4)))
  }
  cat("  status : ", prov$status %||% NA_character_, "   notes: ", NROW(hs_notes(x)),
      "   backend: ", prov$backend$name %||% NA_character_, "\n", sep = "")
  cat("  Use write_outputs(x, dir) to write the plans, the summary and the trajectories figure.\n")
  invisible(x)
}

#' @rdname write_outputs
#' @export
write_outputs.hs_simulation <- function(x, dir, what = c("tables", "figures"),
                                        overwrite = FALSE, create_dir = TRUE, ...) {
  what <- match.arg(what, several.ok = TRUE)
  dir  <- .hs_dir_ready(dir, create_dir = create_dir)
  tag  <- hs_provenance(x)$tag %||% "sim"
  plan <- list()
  if ("tables" %in% what) {
    plan[[length(plan) + 1L]] <- list(kind = "table", name = "plans",
                                      path = file.path(dir, paste0("stacking_plans_", tag, ".txt")))
    plan[[length(plan) + 1L]] <- list(kind = "table", name = "summary",
                                      path = file.path(dir, paste0("stacking_summary_", tag, ".txt")))
  }
  if ("figures" %in% what)
    plan[[length(plan) + 1L]] <- list(kind = "figure", name = "trajectories",
                                      path = file.path(dir, paste0("Trajectories_", tag, ".png")))
  if (!length(plan)) return(invisible(.hs_manifest()))
  .hs_check_overwrite(vapply(plan, `[[`, character(1L), "path"), overwrite)
  rows <- lapply(plan, function(p) {
    if (p$kind == "table")
      return(.hs_write_table(x[[p$name]], p$path, overwrite = TRUE, object = p$name))
    bytes <- tryCatch({
      pt <- plot_simulation_trajectories(x$trajectories, trait_name = tag)
      # one replicate has no ribbon: ggplot2's "Removed n rows" note on the
      # empty ribbon is cosmetic and is not a note of the run
      withCallingHandlers(
        hs_render_figure(pt, p$path, width_mm = 254, height_mm = 127, dpi = 150,
                         overwrite = TRUE),
        warning = function(w) {
          if (grepl("Removed .* row", conditionMessage(w))) invokeRestart("muffleWarning")
        })
      as.numeric(file.info(p$path)$size)
    }, error = function(e) {
      .hs_note("write_outputs", "figure_failed",
               paste0("figure '", p$name, "' could not be rendered to ", p$path,
                      ": ", conditionMessage(e)))
      NA_real_
    })
    .hs_manifest(p$path, "figure", p$name, bytes)
  })
  invisible(do.call(rbind, rows))
}


# ---------------------------------------------------------------------------
# Internal — the dispatcher (D4 §3.7; bounded queue, sc/F08)
# ---------------------------------------------------------------------------

# Run `run_job(k)` for k in 1..n_jobs and return the records in job order.
# Sequential: one after the other. Parallel (n_cores > 1 with future +
# future.apply): at most `n_cores` futures alive at a time, refilled as each
# resolves (0.8.8 created every future up front, which blocked once the
# workers were full). `on_job(k, record)` fires in this process as each job's
# record arrives (the live monitor). Each job seeds its own stream and leaves
# the worker's stream as it found it, so the futures take seed = NULL as in
# 0.8.8 (the generator kind in a worker then stays the default, which keeps
# parallel and sequential runs identical).
.stk_dispatch <- function(n_jobs, run_job, n_cores, on_job = NULL) {
  out <- vector("list", n_jobs)
  if (!.hs_parallel_ok(n_cores)) {
    for (k in seq_len(n_jobs)) {
      out[[k]] <- run_job(k)
      if (!is.null(on_job)) on_job(k, out[[k]])
    }
    return(out)
  }
  .hs_with_plan(n_cores, {
    mk <- function(k) future::future(run_job(k), seed = NULL)
    futs <- vector("list", n_jobs)
    pending <- integer(0L)
    nxt <- 1L
    while (nxt <= n_jobs || length(pending)) {
      while (nxt <= n_jobs && length(pending) < n_cores) {
        futs[[nxt]] <- mk(nxt)
        pending <- c(pending, nxt)
        nxt <- nxt + 1L
      }
      ready <- pending[vapply(futs[pending], future::resolved, logical(1L))]
      if (length(ready)) {
        for (k in ready) {
          out[[k]] <- future::value(futs[[k]])
          futs[k]  <- list(NULL)
          if (!is.null(on_job)) on_job(k, out[[k]])
        }
        pending <- setdiff(pending, ready)
      } else {
        Sys.sleep(0.2)
      }
    }
  })
  out
}
