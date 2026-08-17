# =============================================================================
# stacking_cross.R
#
# Crossing backends for the stacking simulation. A single dispatcher,
# .cross_backend_dispatch(), routes a crossing design to one of two clearly
# separated, interchangeable models — both GenoMeta -> GenoMeta, so a cycle can
# swap backends without touching the rest of the driver:
#
#   "block"             — pure-R, dependency-free. The haplotype blocks are the
#                         crossing loci: each offspring gamete inherits, per
#                         block, one whole parental haplotype (the parent's real
#                         marker alleles across that block), with free
#                         recombination BETWEEN blocks and none WITHIN a block.
#                         Implemented by .cross_block().
#   "genomicSimulation" — opt-in (Suggests). Marker-level crossing with real
#                         linkage/recombination via the Phase-0 bridge
#                         (gs_load_geno_meta / gs_cross_plan / gs_group_to_geno_meta).
#
# Both return list(geno_meta = <offspring GenoMeta>, plans = <one row per
# offspring: Parent1, Parent2, CrossID, F1mean>). F1mean is each offspring's
# breeding value: genomicSimulation's native BV, or — for the block backend —
# the offspring's total local GEBV (summed over blocks) when effects are
# available (NA only when no effects are supplied). Offspring are named
# "<CrossID>n<NN>", identical to gs_cross_plan(), so the two backends produce
# directly comparable output.
#
# NB: .cross_block() superseded the legacy genotype-level
# .simulate_offspring_gebv(); the driver (run_stacking_simulation(), issue 03)
# now routes every cross through .cross_backend_dispatch().
# =============================================================================


# ---------------------------------------------------------------------------
# Internal — expand a crossing design to one row per offspring
# ---------------------------------------------------------------------------

# Mirrors gs_cross_plan()'s expansion/naming so the block and genomicSimulation
# backends generate the same offspring identifiers. Returns a data.table with
# Parent1, Parent2, CrossID.
.expand_cross_plan <- function(cross_design, offspring_prefix = "F") {
  cd <- as.data.frame(cross_design, stringsAsFactors = FALSE)
  if (!all(c("Parent1", "Parent2") %in% names(cd)))
    stop(".expand_cross_plan: 'cross_design' needs 'Parent1' and 'Parent2'.")

  cross_id <- if ("CrossID" %in% names(cd)) cd$CrossID
              else if ("Cross_ID" %in% names(cd)) cd$Cross_ID
              else paste0(offspring_prefix, sprintf("%04i", seq_len(nrow(cd))))

  reps <- if ("n_offspring" %in% names(cd)) as.integer(cd$n_offspring)
          else rep(1L, nrow(cd))
  idx  <- rep(seq_len(nrow(cd)), times = reps)

  data.table(
    Parent1 = cd$Parent1[idx],
    Parent2 = cd$Parent2[idx],
    CrossID = paste0(cross_id[idx], "n", sprintf("%02i", sequence(reps)))
  )
}


# ---------------------------------------------------------------------------
# Internal — per-marker inheritance unit ids (one block = one indivisible unit)
# ---------------------------------------------------------------------------

# Returns a character vector aligned to the marker dimension: the block id of
# each marker, with unblocked markers each forming their own singleton unit
# (free recombination). `markers` is dimnames(geno_arrays)[[2]].
.block_units <- function(geno_meta, markers) {
  bmm <- geno_meta$block_marker_map_info
  if (is.null(bmm)) return(markers)                 # every marker its own unit
  unit <- as.character(bmm$Block)[match(markers, as.character(bmm$Marker))]
  na   <- is.na(unit)
  unit[na] <- markers[na]                           # unblocked -> singleton unit
  unit
}


# Draw one gamete from parent `p`: for each inheritance unit pick one of the
# parent's two haplotypes, then copy that haplotype's alleles for the unit's
# markers. Vectorised over markers. Returns a length-n_markers allele vector.
.block_gamete <- function(arr_p, unit_id, unique_units) {
  hap_for_unit <- sample(c(1L, 2L), length(unique_units), replace = TRUE)
  hv <- hap_for_unit[match(unit_id, unique_units)]
  arr_p[cbind(seq_len(nrow(arr_p)), hv)]
}


# ---------------------------------------------------------------------------
# Internal — pure-R block crossing (GenoMeta -> GenoMeta)
# ---------------------------------------------------------------------------

# Per-offspring breeding value (F1mean) = each offspring's total local GEBV
# (summed over blocks), the block-backend analogue of genomicSimulation's BV.
# Returns NA for every offspring when no effects are available.
.block_f1mean <- function(off_meta, eff, off_names) {
  if (is.null(eff)) return(rep(NA_real_, length(off_names)))
  gebv <- tryCatch(
    compute_profile_localgebv(off_meta, effect_df = eff, level = "genotype"),
    error = function(e) NULL)
  if (is.null(gebv) || nrow(gebv) == 0L)
    return(rep(NA_real_, length(off_names)))
  cols <- setdiff(names(gebv), "Block")
  tot  <- colSums(as.matrix(gebv[, cols, with = FALSE]), na.rm = TRUE)
  unname(tot[off_names])
}


# Realise `cross_design` under the block model and return a fresh offspring
# GenoMeta (+ a per-offspring plan). Diploid only (gamete model). When effects
# are available (`effect_df` or the GenoMeta's `marker_effects`), each
# offspring's F1mean (total local GEBV) is computed; otherwise F1mean is NA.
.cross_block <- function(geno_meta, cross_design, offspring_prefix = "F",
                         effect_df = NULL) {
  if (!is_geno_meta(geno_meta))
    stop(".cross_block: 'geno_meta' is not a valid GenoMeta.")
  ploidy <- geno_meta$metadata$ploidy %||% 2L
  if (ploidy != 2L)
    stop(".cross_block: the block backend models a diploid gamete (ploidy 2); ",
         "got ploidy ", ploidy, ". Use cross_backend = 'genomicSimulation'.")

  arr     <- geno_meta$geno_arrays
  samples <- dimnames(arr)[[1L]]
  markers <- dimnames(arr)[[2L]]
  n_mk    <- length(markers)

  unit_id      <- .block_units(geno_meta, markers)
  unique_units <- unique(unit_id)

  exp   <- .expand_cross_plan(cross_design, offspring_prefix)
  n_off <- nrow(exp)
  if (n_off == 0L)
    stop(".cross_block: 'cross_design' expands to zero offspring.")

  miss <- setdiff(unique(c(exp$Parent1, exp$Parent2)), samples)
  if (length(miss) > 0L)
    stop(".cross_block: parent(s) not in the pool: ",
         paste(miss, collapse = ", "))

  off_names <- make.unique(exp$CrossID, sep = "_")
  out_arr <- array(NA_character_, dim = c(n_off, n_mk, 2L),
                   dimnames = list(off_names, markers, c("1", "2")))

  for (i in seq_len(n_off)) {
    arr_p1 <- arr[exp$Parent1[i], , ]            # [n_markers x 2]
    arr_p2 <- arr[exp$Parent2[i], , ]
    out_arr[i, , 1L] <- .block_gamete(arr_p1, unit_id, unique_units)
    out_arr[i, , 2L] <- .block_gamete(arr_p2, unit_id, unique_units)
  }

  meta_md <- geno_meta$metadata
  meta_md$n_samples <- n_off
  meta_md$n_markers <- n_mk

  off_meta <- new_geno_meta(
    geno_arrays           = out_arr,
    map_info              = geno_meta$map_info,
    block_info            = geno_meta$block_info,
    block_marker_map_info = geno_meta$block_marker_map_info,
    metadata              = meta_md,
    marker_effects        = geno_meta$marker_effects)

  eff    <- .coerce_single_effect_df(effect_df %||% geno_meta$marker_effects)
  f1mean <- .block_f1mean(off_meta, eff, off_names)
  plans  <- data.table(Parent1 = exp$Parent1, Parent2 = exp$Parent2,
                       CrossID = off_names, F1mean = f1mean)

  list(geno_meta = off_meta, plans = plans)
}


# ---------------------------------------------------------------------------
# Internal — crossing-backend dispatcher
# ---------------------------------------------------------------------------

# Route a crossing design to the chosen backend. Both backends return the same
# shape: list(geno_meta = <offspring GenoMeta>, plans = data.table(Parent1,
# Parent2, CrossID, F1mean), session = <gs_handle or NULL>).
#
# `session` implements the persistent genomicSimulation session (FR1): NULL for
# the block backend, and for gS either NULL (first cross of a job — open one from
# `geno_meta`, the single founder write) or a live handle to reuse (later cycles
# cross in-memory, no re-write / re-load). The returned `session` points at the
# offspring group and is threaded to the next cross by .stk_advance()/.run_one().
.cross_backend_dispatch <- function(cross_backend = c("block", "genomicSimulation"),
                                    geno_meta, cross_design,
                                    effect_df        = NULL,
                                    offspring_prefix = "F",
                                    work_dir         = tempdir(),
                                    prefix           = "gs",
                                    session          = NULL) {
  cross_backend <- match.arg(cross_backend)

  if (cross_backend == "block") {
    res <- .cross_block(geno_meta, cross_design, offspring_prefix,
                        effect_df = effect_df)
    return(list(geno_meta = res$geno_meta, plans = res$plans, session = NULL))
  }

  # genomicSimulation — persistent session: open once at the founders, then reuse
  # the live handle each cycle (no per-cycle genotype write + load.data()).
  handle <- session %||% gs_load_geno_meta(geno_meta, effect_df,
                                           work_dir = work_dir, prefix = prefix)
  cp     <- gs_cross_plan(handle, cross_design,
                          offspring_prefix = offspring_prefix, advance = TRUE)
  # Offspring FAM from the plan (real parents) — driver re-derives compact FAM in
  # .stk_rename_generation(), but this keeps the bridge-level GenoMeta correct too.
  off    <- gs_group_to_geno_meta(cp$handle, cross_plan = cp$plans)

  cpd <- as.data.table(cp$plans)
  plans <- data.table(
    Parent1 = as.character(cpd$Parent1),
    Parent2 = as.character(cpd$Parent2),
    CrossID = as.character(cpd$CrossID),
    F1mean  = as.numeric(cpd$F1mean))

  list(geno_meta = off, plans = plans, session = cp$handle)
}
