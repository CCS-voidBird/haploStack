# =============================================================================
# Evaluate_stacking_blocks.R  —  user-level pipeline, no package code touched
#
# Counts how many of a SELECTED haplotype catalogue's target subtypes each line
# in a GenoMeta carries, and draws the result as a dummy alignment merged with a
# per-sample bar panel.
#
#   1. .hs_catalogue_targets()    catalogue (any package shape) -> (Block, Subtype)
#      resolve_search_subtypes()  fallback when the table carries no signature,
#                                 only (Block, Trait, Effect): finds the exact
#                                 subtype by matching this panel's catalogue
#   2. count_selected_subtypes()  targets -> per-sample counts + dummy matrix
#   3. plot_selected_alignment()  dummy alignment, offspring and parents each
#                                 ordered small -> large WITHIN their own group
#      plot_selected_meta()       the same panel merged with the hit-count bars
#                                 via plot_meta_alignment()
#   evaluate_selected_subtypes()  runs all of it and returns everything
#
# Counting rules (deliberate, and the whole point of the "dummy" wording):
#   * count_by = "haplotype" — the ACTUAL number of copies of the target subtype
#     in that sample x block (0..ploidy).
#   * count_by = "genotype"  — mere OCCURRENCE: does that sample x block carry
#     the subtype at all (0/1).
#   * Every true shot is worth 1. The subtype's effect NEVER enters the score —
#     it is carried through to the tables for reference only.
#   * The dummy GEBV of a sample is the sum of its counts, i.e. exactly the
#     column sum of the dummy Block x Sample matrix.
# =============================================================================

# USAGE
#
# source("Evaluate_stacking_blocks.R")
#
# meta <- load_geno_meta("my_geno_meta.rds")
# meta <- attach_family_history(meta)        # needed for GID/FID/Label + figures
#
# ana <- analyze_block_catalogue(meta)       # or run_haplotype_catalogue(), or
# res <- evaluate_selected_subtypes(meta, ana, count_by = "haplotype")
#
# res$by_sample                 # GID FID IID Label nTarget nCarried Frac GEBV Rank
# res$dummy                     # Block x Sample count matrix (colSums == GEBV)
# res$detail                    # one row per sample x target, zeros included
# res$plots[["gen_1"]][["F_1"]] # dummy alignment
# res$meta [["gen_1"]][["F_1"]] # alignment + dummy-GEBV bars, one figure

library(haploStack)
library(data.table)


# ---------------------------------------------------------------------------
# Small local helpers
# ---------------------------------------------------------------------------

# Pick the effect table for one trait out of geno_meta$marker_effects. Mirrors
# haploStack:::.resolve_effect_df(), but returns the table so it can be handed
# to compute_block_catalogue(effect_df=) even when the GenoMeta stores a bare
# (unnamed) effects table and the search table still names a trait.
.hs_effect_table <- function(geno_meta, trait) {
  me <- geno_meta$marker_effects
  if (is.null(me))
    stop("evaluate_selected_subtypes: this GenoMeta carries no marker_effects ",
         "(attach them first, e.g. via fit_genomic_prediction()).")
  if (is.data.frame(me)) return(me)
  if (length(me) == 0L)
    stop("evaluate_selected_subtypes: geno_meta$marker_effects is empty.")
  if (!is.na(trait) && nzchar(trait)) {
    if (trait %in% names(me)) return(me[[trait]])
    stop("evaluate_selected_subtypes: trait '", trait, "' is not in ",
         "geno_meta$marker_effects. Available: ",
         paste(names(me), collapse = ", "), ".")
  }
  if (length(me) == 1L) return(me[[1L]])
  stop("evaluate_selected_subtypes: the search table has no Trait column but ",
       "the GenoMeta stores ", length(me), " traits (",
       paste(names(me), collapse = ", "), "); add a Trait column.")
}

# (IID, GID, FID, Label) from the family history, following the package's own
# convention (see .sample_generation_map() in R/stacking_evaluate.R): the oldest
# generation in which a sample is an offspring wins, else its oldest appearance.
.hs_pedigree_map <- function(geno_meta) {
  empty <- data.table(IID = character(0L), GID = character(0L),
                      FID = character(0L), Label = character(0L))
  h <- geno_meta$metadata$family_history
  if (is.null(h) || length(h) == 0L) return(empty)
  gn <- names(h)
  if (is.null(gn)) gn <- paste0("gen_", seq_along(h))
  long <- rbindlist(lapply(seq_along(h), function(g) {
    d <- as.data.table(h[[g]])
    if (nrow(d) == 0L || !all(c("FID", "IID", "Label") %in% names(d)))
      return(NULL)
    data.table(IID = as.character(d$IID), GID = gn[g],
               FID = as.character(d$FID), Label = as.character(d$Label),
               gen_i = g)
  }))
  if (nrow(long) == 0L) return(empty)
  long[, is_off := Label == "O"]
  setorder(long, IID, -is_off, gen_i)
  unique(long, by = "IID")[, .(IID, GID, FID, Label)]
}

# A drawn panel's row order, bottom -> top. Same introspection the package uses
# in .discrete_axis_levels(); repeated here so this script calls nothing internal.
.hs_panel_levels <- function(p, which = "y") {
  lv <- tryCatch({
    pp <- ggplot2::ggplot_build(p)$layout$panel_params[[1L]]
    as.character(pp[[which]]$get_labels())
  }, error = function(e) character(0L))
  lv[!is.na(lv)]
}

# Normalise an external search table: a data.frame or a file path, columns
# matched case-insensitively. Block is required; Effect is required unless a
# ready-made signature column is present.
.hs_read_search_table <- function(x) {
  if (is.character(x) && length(x) == 1L) x <- data.table::fread(x)
  if (!is.data.frame(x))
    stop("evaluate_selected_subtypes: 'search_tab' must be a data.frame or a ",
         "path to a delimited file.")
  d <- as.data.table(x)
  pick <- function(cands) {
    i <- which(tolower(names(d)) %in% tolower(cands))[1L]
    if (is.na(i)) NULL else names(d)[i]
  }
  cb <- pick("Block"); ce <- pick("Effect"); ct <- pick("Trait")
  cs <- pick(c("Subtype", "signature", "Signature"))
  if (is.null(cb))
    stop("evaluate_selected_subtypes: 'search_tab' needs a 'Block' column.")
  if (is.null(ce) && is.null(cs))
    stop("evaluate_selected_subtypes: 'search_tab' needs an 'Effect' column ",
         "(or an already-resolved 'Subtype'/'signature' column).")

  out <- data.table(
    .row         = seq_len(nrow(d)),
    Trait        = if (is.null(ct)) NA_character_ else as.character(d[[ct]]),
    Block        = as.character(d[[cb]]),
    EffectTarget = if (is.null(ce)) NA_real_ else as.numeric(d[[ce]]),
    SubtypeGiven = if (is.null(cs)) NA_character_ else as.character(d[[cs]]))
  out[is.na(Trait) | !nzchar(Trait), Trait := NA_character_]
  out[!is.na(SubtypeGiven) & !nzchar(SubtypeGiven), SubtypeGiven := NA_character_]
  out[, TraitKey := fifelse(is.na(Trait), "", Trait)]
  out[]
}


# ---------------------------------------------------------------------------
# STEP 1 — a selected haplotype catalogue -> the (Block, Subtype) target list
# ---------------------------------------------------------------------------

# Accepts every shape the package produces, and takes EVERY row as a target —
# no fitness selection, no one-per-block reduction, so a catalogue holding
# several subtypes for the same block is honoured as several targets:
#   * analyze_block_catalogue() result -> $selected (or $top_blocks / $subtypes
#     via `use`)
#   * run_haplotype_catalogue()-shaped named list of per-trait $hap_catalogue
#   * a catalogue data.frame: Block + signature/Subtype (+ Trait, Effect, Count)
#   * a path to such a table
# Returns data.table(Trait, Block, Subtype, Effect) or NULL when the table has
# no signature column at all (the caller then falls back to effect matching).
.hs_catalogue_targets <- function(catalogue,
                                  use = c("selected", "top_blocks", "subtypes")) {
  use <- match.arg(use)
  if (is.character(catalogue) && length(catalogue) == 1L)
    catalogue <- data.table::fread(catalogue)

  # analyze_block_catalogue() result
  if (is.list(catalogue) && !is.data.frame(catalogue) &&
      any(c("selected", "top_blocks", "subtypes") %in% names(catalogue))) {
    tb <- catalogue[[use]]
    if (is.null(tb))
      stop(".hs_catalogue_targets: the analysis result has no '", use, "' table.")
    long <- as.data.table(tb)

  # run_haplotype_catalogue()-shaped named list (one element per trait)
  } else if (is.list(catalogue) && !is.data.frame(catalogue)) {
    nm <- names(catalogue)
    if (is.null(nm) || any(!nzchar(nm)))
      stop(".hs_catalogue_targets: a list catalogue must be named, one element ",
           "per trait.")
    long <- rbindlist(lapply(nm, function(tr) {
      el <- catalogue[[tr]]
      hc <- if (is.data.frame(el)) el else el$hap_catalogue
      if (is.null(hc))
        stop(".hs_catalogue_targets: element '", tr, "' has no 'hap_catalogue'.")
      d <- as.data.table(hc)
      d[, Trait := tr]
      d
    }), fill = TRUE)

  } else if (is.data.frame(catalogue)) {
    long <- as.data.table(catalogue)
  } else {
    stop(".hs_catalogue_targets: 'catalogue' must be a package catalogue ",
         "result, a data.frame, or a path to one.")
  }

  if (!"Block" %in% names(long))
    stop(".hs_catalogue_targets: the catalogue needs a 'Block' column.")
  scol <- intersect(c("Subtype", "signature", "Signature"), names(long))
  if (length(scol) == 0L) return(NULL)          # -> effect-matching fallback

  out <- data.table(
    Trait   = if ("Trait" %in% names(long)) as.character(long$Trait)
              else NA_character_,
    Block   = as.character(long$Block),
    Subtype = as.character(long[[scol[1L]]]),
    Effect  = if ("Effect" %in% names(long)) as.numeric(long$Effect)
              else NA_real_)
  out <- out[!is.na(Block) & nzchar(Block) & !is.na(Subtype)]
  if (nrow(out) == 0L)
    stop(".hs_catalogue_targets: the catalogue holds no usable (Block, ",
         "Subtype) rows.")
  unique(out)[]
}


# ---------------------------------------------------------------------------
# STEP 1b — fallback: (Block, Trait, Effect) -> the exact Subtype signature
# ---------------------------------------------------------------------------

#' @param search_tab     External table: Block, Trait, Effect (+ optional
#'                       Subtype/signature, which is then honoured as given).
#' @param tol,rel_tol    A candidate matches when
#'                       |Effect_panel - Effect_target| <= max(tol,
#'                       rel_tol * |Effect_target|).
#' @param center_by_freq Must match how the search table's effects were made.
#' @param match_by       "effect" (default) or "rank": within each (Trait,
#'                       Block) the search rows are ranked by descending Effect
#'                       and paired with the panel's subtypes ranked the same
#'                       way — use when only the ordering is comparable.
#' @param on_ambiguous   Several subtypes inside the tolerance (real: distinct
#'                       allele combinations can share an effect):
#'                       "all" keeps every one of them (default, honest),
#'                       "best_count" keeps the commonest, "first" the nearest.
#' @return list(report, targets, candidates)
resolve_search_subtypes <- function(geno_meta, search_tab,
                                    tol            = 1e-8,
                                    rel_tol        = 1e-6,
                                    center_by_freq = TRUE,
                                    match_by       = c("effect", "rank"),
                                    on_ambiguous   = c("all", "best_count",
                                                       "first"),
                                    prefilter      = TRUE,
                                    verbose        = TRUE) {

  match_by     <- match.arg(match_by)
  on_ambiguous <- match.arg(on_ambiguous)
  if (!is_geno_meta(geno_meta))
    stop("resolve_search_subtypes: 'geno_meta' is not a valid GenoMeta.")

  want  <- .hs_read_search_table(search_tab)
  known <- unique(as.character(geno_meta$block_marker_map_info$Block))

  # --- this panel's subtype catalogue, once per trait ----------------------
  cand <- rbindlist(lapply(unique(want$TraitKey), function(k) {
    tr  <- if (nzchar(k)) k else NA_character_
    blk <- intersect(unique(want[TraitKey == k, Block]), known)
    if (length(blk) == 0L) return(NULL)
    eff <- .hs_effect_table(geno_meta, tr)
    gm  <- if (isTRUE(prefilter))
             filter_geno_meta(geno_meta, blocks = blk) else geno_meta
    cd  <- as.data.table(compute_block_catalogue(
             gm, effect_df = eff, center_by_freq = center_by_freq,
             group_by = "signature"))
    cd <- cd[Block %in% blk, .(TraitKey = k,
                               Block    = as.character(Block),
                               Subtype  = as.character(signature),
                               Effect   = as.numeric(Effect),
                               Count    = as.integer(Count))]
    cd[, Freq := Count / sum(Count), by = Block]
    cd[]
  }))
  if (is.null(cand) || nrow(cand) == 0L)
    cand <- data.table(TraitKey = character(0L), Block = character(0L),
                       Subtype = character(0L), Effect = numeric(0L),
                       Count = integer(0L), Freq = numeric(0L))

  # --- rank mode: pair by within-block effect rank -------------------------
  if (match_by == "rank") {
    want[, RankIn := frank(-EffectTarget, ties.method = "first"),
         by = .(TraitKey, Block)]
    cnd <- copy(cand)[, RankIn := frank(-Effect, ties.method = "first"),
                      by = .(TraitKey, Block)]
    rep <- merge(want[, .(.row, Trait, TraitKey, Block, EffectTarget, RankIn)],
                 cnd, by = c("TraitKey", "Block", "RankIn"), all.x = TRUE)
    rep[, `:=`(Delta      = Effect - EffectTarget,
               nCandidate = NA_integer_, nWithinTol = NA_integer_,
               Status     = fifelse(is.na(Subtype),
                                    fifelse(Block %in% known,
                                            "no_such_rank", "block_absent"),
                                    "matched_by_rank"))]
    report <- rep[, .(.row, Trait, Block, EffectTarget, Subtype, Effect, Delta,
                      Count, Freq, nCandidate, nWithinTol, Status)]
    setorder(report, .row)
  } else {

    # --- effect mode: nearest candidate inside the tolerance ---------------
    report <- rbindlist(lapply(seq_len(nrow(want)), function(i) {
      w  <- want[i]
      ci <- cand[TraitKey == w$TraitKey & Block == w$Block]
      row <- function(sel, status, nwt) data.table(
        .row = w$.row, Trait = w$Trait, Block = w$Block,
        EffectTarget = w$EffectTarget,
        Subtype = ci$Subtype[sel], Effect = ci$Effect[sel],
        Delta = ci$Effect[sel] - w$EffectTarget,
        Count = ci$Count[sel], Freq = ci$Freq[sel],
        nCandidate = nrow(ci), nWithinTol = nwt, Status = status)

      if (nrow(ci) == 0L)
        return(data.table(.row = w$.row, Trait = w$Trait, Block = w$Block,
                          EffectTarget = w$EffectTarget,
                          Subtype = NA_character_, Effect = NA_real_,
                          Delta = NA_real_, Count = NA_integer_,
                          Freq = NA_real_, nCandidate = 0L, nWithinTol = 0L,
                          Status = if (w$Block %in% known)
                                     "no_subtypes" else "block_absent"))

      # A signature supplied in the search table is authoritative.
      if (!is.na(w$SubtypeGiven)) {
        j <- which(ci$Subtype == w$SubtypeGiven)
        if (length(j) == 0L)
          return(data.table(.row = w$.row, Trait = w$Trait, Block = w$Block,
                            EffectTarget = w$EffectTarget,
                            Subtype = w$SubtypeGiven, Effect = NA_real_,
                            Delta = NA_real_, Count = 0L, Freq = 0,
                            nCandidate = nrow(ci), nWithinTol = 0L,
                            Status = "given_absent"))
        return(row(j[1L], "given", 1L))
      }

      d      <- abs(ci$Effect - w$EffectTarget)
      tol_i  <- max(tol, rel_tol * abs(w$EffectTarget), na.rm = TRUE)
      within <- which(d <= tol_i)
      if (length(within) == 0L) return(row(which.min(d), "unmatched", 0L))
      if (length(within) == 1L) return(row(within, "matched", 1L))
      sel <- switch(on_ambiguous,
                    all        = within,
                    best_count = within[which.max(ci$Count[within])],
                    first      = within[which.min(d[within])])
      row(sel, "ambiguous", length(within))
    }))
  }

  ok      <- c("matched", "ambiguous", "given", "matched_by_rank")
  targets <- report[Status %in% ok & !is.na(Subtype)]

  if (isTRUE(verbose)) {
    message("resolve_search_subtypes: ",
            uniqueN(report$.row[report$Status %in% ok]), "/", nrow(want),
            " search rows resolved -> ", nrow(unique(targets[, .(Block, Subtype)])),
            " distinct (Block, Subtype) targets over ",
            uniqueN(targets$Block), " block(s).")
    bad <- report[!Status %in% ok]
    if (nrow(bad) > 0L) {
      message("  unresolved: ",
              paste(sprintf("%s=%d", names(table(bad$Status)),
                            as.integer(table(bad$Status))), collapse = ", "))
      um <- bad[Status == "unmatched" & !is.na(Delta)]
      if (nrow(um) > 0L)
        message("  nearest-candidate |Delta| median = ",
                signif(stats::median(abs(um$Delta)), 3),
                " — a large, systematic value means the search table's effects ",
                "were built on another scale (check center_by_freq / the ",
                "source panel), not that the haplotypes are absent.")
    }
  }

  list(report     = report[],
       targets    = targets[],
       candidates = cand[])
}


# ---------------------------------------------------------------------------
# STEP 2 — count the target subtypes across the panel
# ---------------------------------------------------------------------------

#' @param targets   data.table with at least Block + Subtype (Trait / Effect are
#'                  carried through for reference only — they never enter the
#'                  score). The `targets` element of .hs_catalogue_targets() or
#'                  resolve_search_subtypes() fits as-is.
#' @param count_by  "haplotype" (default) — the actual number of copies of the
#'                  subtype in that sample x block (0..ploidy); "genotype" —
#'                  occurrence only, 0/1 per sample x block.
#' @return list(by_sample, detail, dummy, targets, summary, count_by)
count_selected_subtypes <- function(geno_meta, targets,
                                    count_by      = c("haplotype", "genotype"),
                                    output_prefix = NULL) {

  count_by <- match.arg(count_by)
  tg <- as.data.table(copy(targets))
  if (!all(c("Block", "Subtype") %in% names(tg)))
    stop("count_selected_subtypes: 'targets' needs 'Block' and 'Subtype'.")
  if (!"Effect" %in% names(tg)) tg[, Effect := NA_real_]
  if (!"Trait"  %in% names(tg)) tg[, Trait  := NA_character_]
  tg <- unique(tg[!is.na(Subtype), .(Trait, Block, Subtype, Effect)])

  # One keyed pass over every (Block, Subtype) target — several blocks and
  # several subtypes per block are exactly what this function is for.
  key <- unique(tg[, .(Block, Subtype)])
  car <- count_subtype_carriage(geno_meta, key, count_by = count_by)

  # Every true shot is worth 1: the score is the COUNT of hits, never the
  # effect. Effect/Trait ride along for reference only.
  detail <- merge(car$detail, unique(tg[, .(Block, Subtype, Trait, Effect)]),
                  by = c("Block", "Subtype"), all.x = TRUE,
                  allow.cartesian = TRUE)
  setorder(detail, Block, Subtype, Sample)

  # Dummy Block x Sample matrix: cell = counts of target hits in that block.
  # Same shape as compute_profile_localgebv(level = "genotype"), so
  # plot_localgebv_alignment() takes it unchanged.
  dummy <- dcast(car$detail, Block ~ Sample, value.var = "Count",
                 fun.aggregate = sum, fill = 0L)
  bi <- tryCatch(as.data.table(build_block_location_info(geno_meta)$block_info),
                 error = function(e) NULL)
  if (!is.null(bi)) {
    ord <- bi$Block[bi$Block %in% dummy$Block]
    if (length(ord) == nrow(dummy)) dummy <- dummy[match(ord, Block)]
  }

  # The dummy GEBV IS the column sum of that matrix — computed from it, so the
  # table and the figure can never drift apart.
  scol <- setdiff(names(dummy), "Block")
  gebv <- data.table(IID  = scol,
                     GEBV = as.numeric(colSums(
                       as.matrix(dummy[, ..scol]), na.rm = TRUE)))

  by_sample <- copy(car$by_sample)
  setnames(by_sample, "Sample", "IID")
  by_sample[gebv, on = "IID", GEBV := i.GEBV]
  by_sample[car$detail[Carried == TRUE, .(nBlockHit = uniqueN(Block)),
                       by = .(IID = Sample)],
            on = "IID", nBlockHit := i.nBlockHit]
  by_sample[is.na(nBlockHit), nBlockHit := 0L]

  ped <- .hs_pedigree_map(geno_meta)
  if (nrow(ped) > 0L)
    by_sample[ped, on = "IID", `:=`(GID = i.GID, FID = i.FID, Label = i.Label)]
  else
    by_sample[, `:=`(GID = NA_character_, FID = NA_character_,
                     Label = NA_character_)]

  setcolorder(by_sample, c("GID", "FID", "IID", "Label", "nTarget", "nCarried",
                           "nBlockHit", "Frac", "GEBV"))
  by_sample[, Total := NULL]                 # == GEBV, kept once under one name
  setorder(by_sample, -GEBV, -nCarried, IID)
  by_sample[, Rank := seq_len(.N)]

  summary <- data.table(
    nTarget     = car$n_target,
    nTargetBlock = uniqueN(tg$Block),
    nSample     = car$n_samples,
    count_by    = count_by,
    nCarried    = car$n_carried,
    TotalHits   = sum(by_sample$GEBV),
    MeanGEBV    = mean(by_sample$GEBV),
    MeanFrac    = mean(by_sample$Frac))

  if (!is.null(output_prefix)) {
    fwrite(by_sample, paste0(output_prefix, "_by_sample.txt"), sep = "\t")
    fwrite(detail,    paste0(output_prefix, "_detail.txt"),    sep = "\t")
    fwrite(dummy,     paste0(output_prefix, "_dummy.txt"),     sep = "\t")
    fwrite(tg,        paste0(output_prefix, "_targets.txt"),   sep = "\t")
    message("count_selected_subtypes: written ", output_prefix,
            "_{by_sample,detail,dummy,targets}.txt")
  }

  list(by_sample = by_sample[], detail = detail[], dummy = dummy[],
       targets = tg[], summary = summary[], count_by = count_by)
}


# ---------------------------------------------------------------------------
# STEP 3 — the dummy alignment, ordered small -> large within each group
# ---------------------------------------------------------------------------

# plot_localgebv_alignment() applies a `sample_order` INSIDE each group (see
# .hs_order_within() / the offspring-founder split in R/catalogue_viz.R), so
# handing it the whole panel ranked by dummy GEBV orders offspring and parents
# from small to large SEPARATELY, and the founder boundary line stays put —
# offspring always occupy the upper block, founders the lower one.
#
# `ascending = TRUE` (default) reads UP the y axis: within each group the
# smallest dummy GEBV sits at the bottom and the largest at the top. The vector
# is therefore handed over in DESCENDING order, because the alignment draws the
# first element of `sample_order` at the TOP of its group
# (`y_levels <- rev(c(off, fnd))`, R/catalogue_viz.R).
plot_selected_alignment <- function(res, geno_meta, ascending = TRUE,
                                    fill_name = "Target hits", ...) {
  ord <- copy(res$by_sample)
  setorderv(ord, c("GEBV", "IID"), c(if (ascending) -1L else 1L, 1L))
  mx  <- suppressWarnings(max(as.matrix(res$dummy[, -1L, with = FALSE]),
                              na.rm = TRUE))
  if (!is.finite(mx) || mx <= 0) mx <- 1
  tryCatch(
    plot_localgebv_alignment(res$dummy, geno_meta,
                             sample_order = as.character(ord$IID),
                             fill_limits  = c(0, mx),
                             fill_name    = fill_name, ...),
    error = function(e) {
      warning("plot_selected_alignment: ", conditionMessage(e),
              " — returning NULL (run attach_family_history() first).")
      NULL
    })
}


# ---------------------------------------------------------------------------
# STEP 3b — the meta figure: alignment + per-sample dummy-GEBV bars
# ---------------------------------------------------------------------------

# One plot_meta_alignment() figure per family: the dummy alignment on the left,
# the sample's dummy GEBV (sum of its counts) as horizontal bars on the right,
# both locked to the same row order. Needs the 'cowplot' package.
plot_selected_meta <- function(res, geno_meta, rel_widths = c(3, 1),
                               ascending = TRUE, legend = "bottom",
                               bar_name = "Dummy GEBV", alignment = NULL, ...) {

  al <- if (is.null(alignment))
          plot_selected_alignment(res, geno_meta, ascending = ascending, ...)
        else alignment
  if (is.null(al) || length(al) == 0L) return(NULL)

  bars_dt <- data.frame(Sample = as.character(res$by_sample$IID),
                        Value  = as.numeric(res$by_sample$GEBV),
                        stringsAsFactors = FALSE)
  names(bars_dt)[2L] <- bar_name

  lapply(al, function(gen) lapply(gen, function(p) {
    lv <- .hs_panel_levels(p, "y")
    if (length(lv) == 0L) return(p)
    bars <- plot_phenotype_bars(bars_dt[bars_dt$Sample %in% lv, , drop = FALSE],
                                orientation   = "horizontal",
                                sample_levels = lv,
                                palette       = c("#b2182b", "#4C72B0")) +
      ggplot2::labs(x = bar_name)     # the bars are counts, not a phenotype
    plot_meta_alignment(list(p, bars), rel_widths = rel_widths, legend = legend)
  }))
}


# ---------------------------------------------------------------------------
# One-call wrapper
# ---------------------------------------------------------------------------

#' @param catalogue A selected haplotype catalogue holding the target subtypes:
#'   an analyze_block_catalogue() result, a run_haplotype_catalogue()-shaped
#'   named list, a catalogue data.frame, or a path to one. EVERY row is a
#'   target, so several subtypes per block are counted as several targets.
#'   A table with no signature column falls back to matching (Block, Trait,
#'   Effect) against this panel's own catalogue — see resolve_search_subtypes().
#' @param use Which table of an analyze_block_catalogue() result supplies the
#'   targets: "selected" (default), "top_blocks" or "subtypes".
evaluate_selected_subtypes <- function(geno_meta, catalogue,
                                       count_by       = c("haplotype",
                                                          "genotype"),
                                       use            = c("selected",
                                                          "top_blocks",
                                                          "subtypes"),
                                       plot           = TRUE,
                                       meta_plot      = TRUE,
                                       ascending      = TRUE,
                                       rel_widths     = c(3, 1),
                                       output_prefix  = NULL,
                                       verbose        = TRUE,
                                       ...) {

  count_by <- match.arg(count_by)
  if (!is_geno_meta(geno_meta))
    stop("evaluate_selected_subtypes: 'geno_meta' is not a valid GenoMeta.")

  tg  <- .hs_catalogue_targets(catalogue, use = match.arg(use))
  rep <- NULL
  if (is.null(tg)) {                 # no signature column: resolve by effect
    if (isTRUE(verbose))
      message("evaluate_selected_subtypes: the catalogue carries no subtype ",
              "signature — resolving targets from (Block, Trait, Effect).")
    rs  <- resolve_search_subtypes(geno_meta, catalogue, verbose = verbose, ...)
    tg  <- rs$targets[, .(Trait, Block, Subtype, Effect)]
    rep <- rs$report
    if (nrow(tg) == 0L)
      stop("evaluate_selected_subtypes: no catalogue row could be resolved to ",
           "a subtype — inspect $resolution (Status / Delta).")
  } else if (isTRUE(verbose)) {
    message("evaluate_selected_subtypes: ", nrow(unique(tg[, .(Block, Subtype)])),
            " target (Block, Subtype) pair(s) over ", uniqueN(tg$Block),
            " block(s), counted by ", count_by, ".")
  }

  out <- count_selected_subtypes(geno_meta, tg, count_by = count_by,
                                 output_prefix = output_prefix)
  out$resolution <- rep
  out$plots <- if (isTRUE(plot))
    plot_selected_alignment(out, geno_meta, ascending = ascending) else NULL
  out$meta  <- if (isTRUE(plot) && isTRUE(meta_plot) && !is.null(out$plots))
    plot_selected_meta(out, geno_meta, rel_widths = rel_widths,
                       alignment = out$plots) else NULL

  if (!is.null(output_prefix) && !is.null(rep))
    fwrite(rep, paste0(output_prefix, "_resolution.txt"), sep = "\t")
  out
}
