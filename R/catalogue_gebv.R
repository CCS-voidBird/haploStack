# =============================================================================
# catalogue_gebv.R — local GEBVs, the block catalogue and the allele sheet
#
#   compute_profile_localgebv()   — exported: per-block local GEBV matrix at the
#                                   genotype or haplotype level (one shared
#                                   compute core, R/catalogue_engine.R)
#   compute_block_catalogue()     — exported: one row per (block, subtype) with
#                                   its allele-combination signature and count
#   write_block_allele_sheet()    — exported: per-block allele counts / freq
#   .compute_block_signatures()   — internal: per-copy block signatures
#   .validate_effect_encoding()   — internal: effect alleles vs the container
#
# Centring is one table (.HS_CENTER_DEFAULT, R/catalogue_engine.R; 0.9.0 row
# 15 c); the arithmetic lives in .hap_effects() / .local_gebv_engine().
# (Header rewritten 2026-11-09, cc/F23: the previous one described functions
# that no longer exist here.)
# =============================================================================

# ---------------------------------------------------------------------------
# Internal: resolve a single-trait effect table
# ---------------------------------------------------------------------------

#' Resolve a single-trait effect table for the per-block GEBV functions
#'
#' An explicit \code{effect_df} always wins. Otherwise the effects are taken
#' from \code{geno_meta$marker_effects}, which may be the canonical named list of
#' per-trait tables or a bare data.frame. \code{trait} selects an element of the
#' list; when it is \code{NULL} and exactly one trait is stored that one is used,
#' otherwise an error asks the caller to disambiguate. Detecting effects from the
#' GenoMeta is a convenience — they are never required to construct one.
#'
#' @param geno_meta A valid \code{GenoMeta} list.
#' @param effect_df Explicit effect table, or \code{NULL} to fall back to
#'   \code{geno_meta$marker_effects}.
#' @param trait     Optional trait name selecting one element of a stored
#'   per-trait list.
#' @param caller    Calling-function name, used in error messages.
#' @return A single effect data.frame.
#' @keywords internal
#' @noRd
.resolve_effect_df <- function(geno_meta, effect_df, trait = NULL, caller) {
  if (!is.null(effect_df)) return(effect_df)

  me <- geno_meta$marker_effects
  if (is.null(me))
    .hs_abort(caller, ": no effects supplied \u2014 pass 'effect_df' or attach ",
         "'marker_effects' to the GenoMeta (e.g. via fit_genomic_prediction()).")

  if (is.data.frame(me)) return(me)        # bare table (back-compat)

  # canonical named list of per-trait tables
  if (!is.null(trait)) {
    if (!trait %in% names(me))
      .hs_abort(caller, ": trait '", trait, "' not in geno_meta$marker_effects. ",
           "Available: ", paste(names(me), collapse = ", "), ".")
    return(me[[trait]])
  }
  if (length(me) == 1L) return(me[[1L]])
  .hs_abort(caller, ": geno_meta$marker_effects holds multiple traits (",
       paste(names(me), collapse = ", "), "); pass 'effect_df' or set 'trait' ",
       "to choose one.")
}


# ---------------------------------------------------------------------------
# Internal: effect computation
# ---------------------------------------------------------------------------

#' Attach effects and compute haplotype-level GEBV per block
#'
#' Internal helper consumed by \code{compute_profile_localgebv()}.
#' Joins the allele effect table to the long-format allele data, optionally
#' mean-centres by population allele frequency, and sums to per-block
#' haplotype effects.
#'
#' @param long_dt     Long-format \code{data.table} with columns
#'   \code{Sample}, \code{Haplotype_idx}, \code{Marker}, \code{Allele},
#'   \code{Block} (no NA blocks).
#' @param eff         Standardised effect \code{data.table} with columns
#'   \code{Marker}, \code{Allele}, \code{Effect}.
#' @param center        Centring mode applied to the allele effects:
#'   \code{"none"}, \code{"locus_mean"} or \code{"allele_freq"}. See
#'   \code{\link{compute_profile_localgebv}} for what each preserves.
#' @param center_by_freq Logical; centre effects by population allele frequency
#'   before summing.
#' @param freq_dt Optional pre-computed allele-frequency table (Marker, Allele,
#'   AlleleFreq). When \code{NULL} (default) it is derived from \code{long_dt}, as
#'   before. The catalogue engine injects the core's cached table so the
#'   frequency pass is not repeated per trait; the values are identical because
#'   the cached table is built from the same long-format data.
#' @return \code{data.table} with columns \code{Sample}, \code{Haplotype_idx},
#'   \code{Block}, \code{HapEffect}.
#' @keywords internal
#' @noRd
.compute_effect_size <- function(long_dt, eff, center = "none",
                                 freq_dt = NULL) {
  center <- match.arg(as.character(center)[1L], HS_CENTER_MODES)

  # Attach Effect via UPDATE JOIN (avoids column-structure drift)
  long_dt[, Effect := NA_real_]
  long_dt[eff, on = .(Marker, Allele), Effect := i.Effect]
  long_dt[is.na(Effect), Effect := 0]

  if (center != "none" && is.null(freq_dt)) {
    freq_dt <- long_dt[, .(AlleleCount = .N), by = .(Marker, Allele)]
    tot_dt  <- long_dt[, .(MarkerTotal = .N), by = Marker]
    freq_dt <- freq_dt[tot_dt, on = "Marker"][
      , AlleleFreq := AlleleCount / MarkerTotal
    ][, c("AlleleCount", "MarkerTotal") := NULL]
  }

  if (center == "allele_freq") {
    # Each allele centred by ITS OWN frequency — a reweighting, so the A-vs-B
    # contrast moves from e_A - e_B to q*e_A - p*e_B and block variance changes.
    long_dt[, AlleleFreq := NA_real_]
    long_dt[freq_dt, on = .(Marker, Allele), AlleleFreq := i.AlleleFreq]
    long_dt[is.na(AlleleFreq), AlleleFreq := 0]

    hap_dt <- long_dt[,
      .(HapEffect = sum(Effect, na.rm = TRUE) -
                    sum(AlleleFreq * Effect, na.rm = TRUE)),
      by = .(Sample, Haplotype_idx, Block)]
  } else if (center == "locus_mean") {
    # ONE mean per marker, subtracted from every allele there — a pure location
    # shift. The contrast is preserved and block variance is unchanged, which is
    # what makes this the mode that agrees with HapSelect's mean_adjust value
    # for value rather than merely on spread.
    long_dt[, LocusMean := NA_real_]
    long_dt[.locus_means(freq_dt, eff), on = "Marker", LocusMean := i.LocusMean]
    long_dt[is.na(LocusMean), LocusMean := 0]

    hap_dt <- long_dt[,
      .(HapEffect = sum(Effect - LocusMean, na.rm = TRUE)),
      by = .(Sample, Haplotype_idx, Block)]
  } else {
    hap_dt <- long_dt[,
      .(HapEffect = sum(Effect, na.rm = TRUE)),
      by = .(Sample, Haplotype_idx, Block)]
  }

  hap_dt
}


# ---------------------------------------------------------------------------
# Exported: compute_profile_localgebv
# ---------------------------------------------------------------------------

#' Compute block-level local GEBVs from a GenoMeta object
#'
#' Calculates per-block genomic estimated breeding values (local GEBVs) for
#' each sample by summing marker allele effects within each haplotype block.
#' Returns a wide-format \code{data.table} with blocks as rows and samples
#' (or sample-haplotype pairs) as columns.
#'
#' Encoding validation: allele values in \code{effect_df} are checked against
#' the GenoMeta encoding. Mismatched alleles trigger a \code{warning()} and are
#' silently excluded from the GEBV sum (they do not contribute any effect).
#' No auto-conversion is performed — ensure \code{effect_df} uses the same
#' encoding as the genotype data before calling this function.
#'
#' @section Centring:
#' The three modes are different operations, not settings of one. Writing
#' \code{e_a} for an allele's effect and \code{p_a} for its population
#' frequency, and \code{mean_i = sum_a p_a * e_a} for the locus mean of marker
#' \code{i}:
#'
#' \tabular{llll}{
#'   \strong{mode} \tab \strong{contribution} \tab \strong{allele contrast} \tab \strong{block variance} \cr
#'   \code{"none"} \tab \code{e_a} \tab \code{e_A - e_B} \tab reference \cr
#'   \code{"locus_mean"} \tab \code{e_a - mean_i} \tab \strong{preserved} \tab \strong{unchanged} \cr
#'   \code{"allele_freq"} \tab \code{e_a * (1 - p_a)} \tab becomes \code{q*e_A - p*e_B} \tab changed
#' }
#'
#' \code{"locus_mean"} subtracts \strong{one} value per marker from every allele
#' at it, so it is a pure location shift: the A-vs-B contrast survives and every
#' block variance is identical to the uncentred one. It is also invariant to
#' adding a constant to both alleles' effects at a marker — the property 0.5.3
#' wanted from centring — and it reproduces HapSelect's
#' \code{mean_adjust = TRUE} \strong{value for value}, because
#' \code{d*e_A + (ploidy - d)*e_B - ploidy*mean_i} is exactly
#' \code{(d - ploidy*p_A) * (e_A - e_B)}.
#'
#' \code{"allele_freq"} centres each allele by \strong{its own} frequency, so
#' the two alleles of a marker are weighted differently. That is a reweighting
#' rather than a shift: block variances change (median factor ~3 on the barley
#' panel) and rankings can move with them. It was the catalogue default up to
#' 0.5.5 and is what reproduces those numbers.
#'
#' Use \code{"locus_mean"} for selection and for any cross-package comparison;
#' \code{"none"} is the default because it is the raw quantity and because it is
#' what both entry points to this engine now agree on.
#'
#' @section What a haplotype slot is:
#' \code{level = "haplotype"} reports the genotype matrix's copies in the order
#' the \strong{reader} wrote them. It is a \strong{slot assignment, not phase}.
#'
#' A genotype field in an unphased source — a PED, a VCF without phasing — is
#' unordered by definition, so a heterozygote's two alleles are written in one
#' fixed order by whatever produced the file. On the GRDC barley panel every
#' heterozygote is written \code{AT} and none \code{TA} (7,584,722 vs 0), which
#' makes slot 1 the A-strand panel-wide. Following the genotype matrix and
#' assigning its copies to slots is exactly what this function should do — but
#' per-slot output is not a chromosome, and nothing in the package checks
#' whether the source was phased. Do not read \code{level = "haplotype"} columns
#' as haplotypes of a real chromosome unless the source itself is phased.
#'
#' \code{level = "genotype"} sums over every ploidy copy, so it is the level to
#' compare against a genotype-level figure from another package. Across levels
#' the factor is between \code{ploidy} and \code{ploidy^2} — a population
#' property that cannot be divided out.
#'
#' @param geno_meta     A valid \code{GenoMeta} list from
#'   \code{\link{new_geno_meta}} or any \code{read_*} function.
#' @param effect_df     data.frame or data.table with \strong{at least three
#'   columns read by position}: col 1 = Marker ID, col 2 = Allele state, col
#'   3 = Effect size. Additional columns are ignored. Column names are ignored.
#'   Allele values must use the \strong{same encoding} as the genotype data.
#'   When \code{NULL} (default), \code{geno_meta$marker_effects} is used as the
#'   effect source (see \code{trait}); an explicit \code{effect_df} always
#'   overrides it.
#' @param trait         Optional trait name selecting one table from a stored
#'   per-trait \code{geno_meta$marker_effects} list when \code{effect_df} is
#'   \code{NULL}. With a single stored trait it is picked automatically; with
#'   several, \code{trait} (or an explicit \code{effect_df}) is required.
#' @param samples       Optional character vector of sample names to include.
#' @param blocks        Optional character vector of block IDs to include.
#' @param chromosomes   Optional list of chromosome region specs (see
#'   \code{\link{filter_geno_meta}}).
#' @param markers       Optional character vector of marker names to include.
#' @param level         Output resolution. \code{"genotype"} (default) sums
#'   \strong{every ploidy copy} per sample (2 for a diploid, 4 for a tetraploid);
#'   \code{"haplotype"} keeps them separate with column names
#'   \code{"<Sample>.<haplotype_idx>"}. See \emph{What a haplotype slot is}.
#' @param center        Centring mode applied to the allele effects before they
#'   are summed: \code{"none"} (default), \code{"locus_mean"} or
#'   \code{"allele_freq"}. See \emph{Centring}.
#' @param center_by_freq \strong{Deprecated} logical, kept for one release:
#'   \code{TRUE} maps to \code{center = "allele_freq"} and \code{FALSE} to
#'   \code{center = "none"}. Supplying both this and \code{center} is an error.
#' @param fill_na       Numeric. Value for (Block, Sample) cells with no
#'   effect information. Default \code{0}.
#' @param .core         Internal. A pre-built trait-invariant compute core
#'   (\code{.catalogue_core()}) shared by \code{\link{run_haplotype_catalogue}}
#'   across traits to avoid rebuilding the long-format allele table per trait.
#'   Not for direct use; cannot be combined with the
#'   samples/markers/chromosomes filters. Default \code{NULL} (build internally).
#'
#' @return A \code{data.table} in wide format:
#' \describe{
#'   \item{\code{Block}}{Block IDs ordered by genomic position.}
#'   \item{Sample columns}{One column per sample (\code{level = "genotype"})
#'     or one per sample × haplotype (\code{level = "haplotype"}).}
#' }
#'
#' @details
#' Effect computation is delegated to the internal
#' \code{.compute_effect_size()} function. Plot functions receive pre-computed
#' GEBV vectors and do not call this function directly.
#'
#' This function uses \strong{update joins} (not right-join reassignments) when
#' attaching Effect and Block columns to the long-format allele table. This
#' prevents the column-structure drift bug present in an earlier version of the
#' catalogue pipeline where \code{Sample} drifted to an unexpected position.
#'
#' @examples
#' \dontrun{
#' meta     <- load_geno_meta("my_geno_meta.rds")
#' eff      <- read_effect_file("effects/trait1.GS.Eff")
#' gebv_mat <- compute_profile_localgebv(meta, eff)
#' dim(gebv_mat)   # n_blocks × (n_samples + 1)
#' }
#'
#' @seealso \code{\link{new_geno_meta}}, \code{\link{filter_geno_meta}},
#'   \code{\link{run_haplotype_catalogue}}
#'
#' @export
compute_profile_localgebv <- function(
  geno_meta,
  effect_df      = NULL,
  trait          = NULL,
  samples        = NULL,
  blocks         = NULL,
  chromosomes    = NULL,
  markers        = NULL,
  level          = c("genotype", "haplotype"),
  center         = NULL,
  center_by_freq = NULL,
  fill_na        = 0,
  .core          = NULL
) {
  center <- .resolve_center(center, center_by_freq, "compute_profile_localgebv")
  level <- match.arg(level)
  .hs_check_meta(geno_meta, "compute_profile_localgebv", "geno_meta")

  # Effect source: explicit effect_df overrides the stored marker_effects; a
  # per-trait list in the GenoMeta is resolved via `trait`.
  effect_df <- .resolve_effect_df(geno_meta, effect_df, trait,
                                  "compute_profile_localgebv")
  eff <- .standardise_effect_dt(effect_df, "compute_profile_localgebv")

  # Compute core: a caller (run_haplotype_catalogue) may pass a pre-built,
  # trait-invariant core via `.core` so the long-format allele table is built
  # ONCE per run rather than once per trait. A shared core may not be combined
  # with per-sample / per-marker / per-chromosome filters (those change the data
  # the core was built from); a `blocks` subset is applied as keep_blocks.
  any_subset_filter <- !is.null(samples) || !is.null(chromosomes) ||
    !is.null(markers)
  if (!is.null(.core)) {
    if (any_subset_filter)
      .hs_abort("compute_profile_localgebv: '.core' cannot be combined with ",
           "samples/markers/chromosomes filters.")
    core        <- .core
    keep_blocks <- blocks
  } else {
    any_filter <- any_subset_filter || !is.null(blocks)
    gm <- if (any_filter)
      filter_geno_meta(geno_meta, keep_samples = samples, keep_markers = markers,
                       blocks = blocks, chromosomes = chromosomes)
    else geno_meta
    core        <- .catalogue_core(gm)
    keep_blocks <- NULL                 # filtering already baked into the core
  }

  # Validate allele encoding against the resolved encoding (warns + drops on
  # mismatch, exactly as before).
  .validate_effect_encoding(eff, core$encoding, observed = core$freq_core$Allele)

  if (isTRUE(core$empty)) {
    .hs_warn("compute_profile_localgebv: empty genotype data after filtering, ",
            "or no markers with block assignments remain.")
    return(data.table(Block = character(0L)))
  }

  .local_gebv_engine(core, eff, center, level = level,
                     fill_na = fill_na, keep_blocks = keep_blocks)
}


# ---------------------------------------------------------------------------
# Block-haplotype catalogue (long, origin-compatible)
# ---------------------------------------------------------------------------

#' Build a block-haplotype effect catalogue
#'
#' Enumerates the distinct haplotype \emph{effect states} per block and their
#' population frequency, reproducing the long-format catalogue of the origin
#' GRDC workflow (e.g. \code{Haplotype_Catalogue_0.5_t3_test1.txt}). Each row is
#' one distinct block-haplotype local-GEBV value: its \code{Effect}, a
#' representative carrier (\code{Sample}, \code{Haplotype}), that carrier's
#' allele-combination \code{signature} (so the subtype's member-marker alleles
#' are readable straight from the row), the number of haplotype copies in the
#' population carrying that state (\code{Count}), and the block's genomic
#' location. \code{Count} sums to \code{n_samples * ploidy} within each block.
#'
#' Columns are ordered \strong{identification-first} — genomic / carrier
#' identity (\code{Chrom}, \code{Pos}, \code{Block}, \code{Sample},
#' \code{Haplotype}, \code{signature}) leads, with the measures and block-span
#' detail following.
#'
#' The effects are the haplotype-level local GEBVs from
#' \code{\link{compute_profile_localgebv}} (\code{level = "haplotype"}); states
#' are grouped by their local-GEBV value (haplotypes with an identical block
#' effect are one catalogue entry, including the zero/absent state).
#'
#' @param geno_meta A valid \code{GenoMeta} with \code{block_marker_map_info}
#'   populated (see \code{\link{attach_block_info}}).
#' @param effect_df data.frame of marker effects (Marker, Allele, Effect). When
#'   \code{NULL}, the GenoMeta's attached \code{marker_effects} are used (see
#'   \code{trait}).
#' @param trait     Optional trait name selecting one table from a stored
#'   per-trait \code{geno_meta$marker_effects} list when \code{effect_df} is
#'   \code{NULL}.
#' @param block_info Optional pre-computed block-location table (as returned by
#'   \code{\link{build_block_location_info}}'s \code{block_info}: \code{Chrom},
#'   \code{Block}, \code{Start}, \code{End}, \code{Pos}, \code{Num_SNPs},
#'   \code{Block_length}). Recomputed from \code{geno_meta} when \code{NULL}.
#' @param center        Centring mode applied to the allele effects:
#'   \code{"none"}, \code{"locus_mean"} or \code{"allele_freq"}. See
#'   \code{\link{compute_profile_localgebv}} for what each preserves.
#' @param center_by_freq \strong{Deprecated} logical form of \code{center}
#'   (\code{TRUE} = \code{"allele_freq"}, \code{FALSE} = \code{"none"}); warns
#'   on every call. The default centring is \code{"none"} since 0.9.0 (it was
#'   \code{"allele_freq"} up to 0.8.8; pass \code{center = "allele_freq"} to
#'   reproduce those numbers).
#' @param blocks Optional character vector of block ids restricting the
#'   catalogue to those blocks (0.9.0, cc/F05; \code{\link{run_haplotype_catalogue}}
#'   passes its \code{selected_blocks}). Default \code{NULL}: every block.
#' @param group_by How distinct subtypes are enumerated within a block.
#'   \code{"effect"} (default) groups haplotype states by their local-GEBV
#'   value (origin granularity); each row carries the \code{signature} of its
#'   representative carrier. \code{"signature"} groups by the
#'   \strong{allele-combination signature} — the ordered alleles across the
#'   block's member markers, read from \code{geno_meta$geno_arrays} — so two
#'   distinct allele combinations sharing an effect become separate rows. Either
#'   way a \code{signature} column is present. Signatures are content-based and
#'   trait-invariant (independent of \code{effect_df}).
#' @param .core Internal. A pre-built trait-invariant compute core
#'   (\code{.catalogue_core()}) shared across traits by
#'   \code{\link{run_haplotype_catalogue}} / \code{\link{analyze_block_catalogue}}.
#'   Not for direct use. Default \code{NULL} (build internally).
#' @param .signatures Internal. A pre-built per-copy block-signature table
#'   (\code{.compute_block_signatures()}); reused across traits so it is computed
#'   once (used on both \code{group_by} paths). Default \code{NULL}.
#'
#' @return A \code{data.table}, one row per distinct block-haplotype state,
#'   ordered by block then descending \code{Count}, with identification columns
#'   first: \code{Chrom, Pos, Block, Sample, Haplotype, signature, Effect, Count,
#'   Block_start, Block_end, Num_SNPs_in_Block, Block_length, Block_chr}. The
#'   \code{signature} is the representative carrier's allele combination over the
#'   block's member markers (under \code{group_by = "signature"} it is also the
#'   grouping key). Effect values and \code{Count} follow the identity columns.
#'
#' @seealso \code{\link{compute_profile_localgebv}},
#'   \code{\link{build_block_location_info}}, \code{\link{run_haplotype_catalogue}}
#' @export
compute_block_catalogue <- function(geno_meta,
                                    effect_df      = NULL,
                                    trait          = NULL,
                                    block_info     = NULL,
                                    center         = NULL,
                                    center_by_freq = NULL,
                                    group_by       = c("effect", "signature"),
                                    blocks         = NULL,
                                    .core          = NULL,
                                    .signatures    = NULL) {
  group_by <- match.arg(group_by)
  center   <- .resolve_center(center, center_by_freq, "compute_block_catalogue")
  if (!is.null(blocks)) blocks <- as.character(blocks)
  .hs_check_meta(geno_meta, "compute_block_catalogue", "geno_meta")

  empty <- data.table(
    Block = character(0L), Effect = numeric(0L), Sample = character(0L),
    Haplotype = integer(0L), Count = integer(0L), Chrom = character(0L),
    Block_start = numeric(0L), Block_end = numeric(0L), Pos = numeric(0L),
    Block_chr = character(0L), Num_SNPs_in_Block = integer(0L),
    Block_length = numeric(0L), signature = character(0L))
  .lead_identifier_cols(empty)

  # Haplotype-level local GEBV (wide: Block + "<Sample>.<haplotype_idx>").
  # Reuse the run-level trait-invariant core when supplied, else build one here.
  effect_df <- .resolve_effect_df(geno_meta, effect_df, trait,
                                  "compute_block_catalogue")
  eff  <- .standardise_effect_dt(effect_df, "compute_block_catalogue")
  core <- if (!is.null(.core)) .core else .catalogue_core(geno_meta)
  .validate_effect_encoding(eff, core$encoding, observed = core$freq_core$Allele)
  lg <- .local_gebv_engine(core, eff, center, level = "haplotype",
                           fill_na = 0, keep_blocks = blocks)
  lg <- as.data.table(lg)
  sample_cols <- setdiff(names(lg), "Block")
  if (nrow(lg) == 0L || length(sample_cols) == 0L) return(empty[])

  # Long form: one row per (Sample, Haplotype, Block) haplotype copy.
  long <- data.table::melt(lg, id.vars = "Block", measure.vars = sample_cols,
                           variable.name = "SampleHap", value.name = "Effect",
                           variable.factor = FALSE)
  long[, Sample    := sub("\\.[0-9]+$", "", SampleHap)]
  long[, Haplotype := as.integer(sub("^.*\\.([0-9]+)$", "\\1", SampleHap))]
  long[, Effect    := as.numeric(Effect)]

  # The allele-combination signature is trait-invariant; a caller running many
  # traits (analyze_block_catalogue / run_haplotype_catalogue) computes it once
  # and injects it via .signatures. It labels each catalogue row's subtype so the
  # member-marker alleles are readable straight from the table.
  sig_dt <- if (!is.null(.signatures)) .signatures
            else .compute_block_signatures(geno_meta)

  # Collapse to distinct subtypes: copy count + representative carrier (smallest
  # Sample/Haplotype). "effect" groups by local-GEBV value (origin behaviour);
  # "signature" groups by the allele-combination content over the block markers.
  if (group_by == "signature") {
    long[sig_dt, on = c("Block", "Sample", "Haplotype" = "Haplotype_idx"),
         signature := i.signature]
    long[is.na(signature), signature := ""]   # absent/missing state
    data.table::setorder(long, Block, signature, Sample, Haplotype)
    cat_dt <- long[, .(Effect    = Effect[1L],
                       Sample    = Sample[1L],
                       Haplotype = Haplotype[1L],
                       Count     = .N),
                   by = .(Block, signature)]
  } else {
    # Round guards floating-point grouping jitter.
    data.table::setorder(long, Block, Sample, Haplotype)
    long[, EffKey := round(Effect, 10L)]
    cat_dt <- long[, .(Effect    = Effect[1L],
                       Sample    = Sample[1L],
                       Haplotype = Haplotype[1L],
                       Count     = .N),
                   by = .(Block, EffKey)]
    cat_dt[, EffKey := NULL]
    # Label each effect-state row with its representative carrier's signature.
    cat_dt[sig_dt, on = c("Block", "Sample", "Haplotype" = "Haplotype_idx"),
           signature := i.signature]
    cat_dt[is.na(signature), signature := ""]
  }

  # Attach block location.
  if (is.null(block_info))
    block_info <- build_block_location_info(geno_meta)$block_info
  bi <- as.data.table(block_info)
  loc <- bi[, .(Block,
                Chrom             = as.character(Chrom),
                Block_start       = as.numeric(Start),
                Block_end         = as.numeric(End),
                Pos               = as.numeric(Pos),
                Num_SNPs_in_Block = as.integer(Num_SNPs),
                Block_length      = as.numeric(Block_length))]
  cat_dt <- merge(cat_dt, loc, by = "Block", all.x = TRUE)
  cat_dt[, Block_chr := Chrom]

  # Identification columns lead (Chrom, Pos, Block, Sample, Haplotype,
  # signature, …); measures and block-span detail follow.
  .lead_identifier_cols(cat_dt)
  data.table::setorder(cat_dt, Block, -Count)
  data.table::setattr(cat_dt, "hs_center", center)   # D3.5: the centre travels with the table
  cat_dt[]
}


# ---------------------------------------------------------------------------
# Internal: allele-combination signatures
# ---------------------------------------------------------------------------

#' Build per-haplotype-copy allele signatures for each block
#'
#' Internal helper for \code{compute_block_catalogue(group_by = "signature")}.
#' For every haplotype copy \code{(Sample, Haplotype_idx)} and every block, it
#' concatenates the alleles across the block's member markers in canonical
#' genomic order, producing a content-based \code{signature} string. The
#' signature is \strong{trait-invariant} (it depends only on genotype, not on
#' marker effects), so it can be computed once and reused across traits.
#'
#' An unobserved allele is written as \code{fill_na} (default \code{"-"}, see
#' \code{\link{HS_FILL_NA_DEFAULT}}) rather than dropped, so every signature is
#' exactly as wide as its block and each character keeps its marker's position.
#' This is what makes the signature a faithful identifier: dropping the gap
#' instead would collapse distinct combinations onto one another — over markers
#' \code{(M1, M2, M3)} the copies \code{(A, T, -)} and \code{(A, -, T)} would both
#' paste to \code{"AT"}. A copy whose alleles are all unobserved therefore still
#' produces a row (an all-\code{fill_na} signature), so the copy count per block
#' is always \code{n_samples x ploidy}.
#'
#' Cost is linear in the genotype matrix: a single melt + one grouped paste over
#' \code{(Block, Sample, Haplotype_idx)}; distinct subtypes are later found by a
#' hashed group-by on the signature (never pairwise string comparison).
#'
#' @param geno_meta A valid \code{GenoMeta} with block assignments.
#' @param fill_na Placeholder for an unobserved allele. Default \code{NULL},
#'   meaning the \code{haploStack.fill_na} option (itself defaulting to
#'   \code{"-"}). Must be a single non-empty string.
#' @param blocks Optional character vector of block IDs. When supplied only those
#'   blocks' member markers enter the melt, which is what makes a single-block
#'   lookup cheap (\code{\link{count_subtype_occurrence}}); the remaining
#'   per-call cost is \code{geno_meta_views()}, so callers doing many lookups
#'   should build this table once and pass it on. Default \code{NULL} (all
#'   blocks).
#' @return \code{data.table} with columns \code{Block}, \code{Sample},
#'   \code{Haplotype_idx} (integer), \code{signature} (character).
#' @keywords internal
#' @noRd
.compute_block_signatures <- function(geno_meta, fill_na = NULL, blocks = NULL,
                                      views = NULL) {
  fill_na <- .hs_fill_na(fill_na)
  prof      <- if (is.null(views)) geno_meta_views(geno_meta) else views   # cc/F09
  meta_cols <- c("Sample", "Haplotype_idx")
  all_mks   <- setdiff(names(prof$alleles), meta_cols)

  empty <- data.table(Block = character(0L), Sample = character(0L),
                      Haplotype_idx = integer(0L), signature = character(0L))

  bm <- prof$marker_meta[!is.na(block_id),
                         .(Marker     = as.character(marker),
                           Block      = as.character(block_id),
                           position   = as.numeric(position))]
  # Restricting to the requested blocks shrinks the melt to their member markers
  # (a handful of columns instead of the whole panel).
  if (!is.null(blocks)) bm <- bm[Block %in% unique(as.character(blocks))]
  marker_cols <- intersect(all_mks, bm$Marker)
  if (length(marker_cols) == 0L) return(empty[])

  long <- data.table::melt(
    prof$alleles[, c(meta_cols, marker_cols), with = FALSE],
    id.vars         = meta_cols,
    measure.vars    = marker_cols,
    variable.name   = "Marker",
    value.name      = "Allele",
    variable.factor = FALSE)
  setDT(long)
  long[, Marker := as.character(Marker)]
  long[, Allele := as.character(Allele)]
  # Fill, never drop: the signature is positional, so a missing allele must keep
  # its slot or distinct combinations collapse onto one another.
  long[is.na(Allele) | !nzchar(trimws(Allele)), Allele := fill_na]

  long[bm, on = "Marker", `:=`(Block = i.Block, position = i.position)]
  long <- long[!is.na(Block)]
  if (nrow(long) == 0L) return(empty[])

  # Canonical order within a copy's block: by genomic position then marker name.
  data.table::setorder(long, Block, Sample, Haplotype_idx, position, Marker)
  sig <- long[, .(signature = paste(Allele, collapse = "")),
              by = .(Block, Sample, Haplotype_idx)]
  sig[, Block         := as.character(Block)]
  sig[, Sample        := as.character(Sample)]
  sig[, Haplotype_idx := as.integer(Haplotype_idx)]
  sig[]
}


# ---------------------------------------------------------------------------
# Exported: haplotype-block marker-allele combination sheet
# ---------------------------------------------------------------------------

#' Write a haplotype-block marker-allele combination sheet
#'
#' Enumerates, for every block, the distinct \strong{marker-allele combinations}
#' (haplotype subtypes) observed in the population and their frequency, as a single
#' tidy table. A combination is the tuple of alleles a haplotype copy carries across
#' the block's member markers, taken in canonical genomic order. Because blocks hold
#' different member markers, the marker set is stored \emph{inside} \code{Markers} /
#' \code{Alleles} cells (semicolon-separated) rather than as per-marker columns, so
#' every block's rows share the same shape and the sheet scales to whole genomes.
#'
#' The combination content is derived from genotype only (independent of any marker
#' effects); it mirrors the signature construction of the catalogue's internal
#' \code{group_by = "signature"} path, and uses the same \code{fill_na}
#' placeholder, so \code{Signature} here and the catalogue's \code{signature} /
#' \code{\link{count_block_subtypes}}'s \code{Subtype} are the same string and can
#' be joined directly.
#'
#' @param geno_meta A valid \code{GenoMeta} with block assignments (see
#'   \code{\link{attach_block_info}}).
#' @param output_file Optional path. When non-\code{NULL}, the table is written as a
#'   tab-separated file via \code{data.table::fwrite()}; the table is returned either
#'   way. Default \code{NULL} (return only).
#' @param blocks Optional character vector of block IDs to restrict to. Default
#'   \code{NULL} (all blocks).
#' @param count_by How occurrences are counted. \code{"haplotype"} (default) counts
#'   haplotype copies, so \code{Count} sums to \code{n_samples x ploidy} within a
#'   block. \code{"genotype"} counts distinct samples carrying the combination on at
#'   least one copy.
#' @param fill_na Placeholder written for an unobserved allele in both
#'   \code{Signature} and \code{Alleles}. Default \code{NULL}, meaning the
#'   \code{haploStack.fill_na} option (itself defaulting to \code{"-"}, see
#'   \code{\link{HS_FILL_NA_DEFAULT}}). Must be a single non-empty string — an
#'   empty placeholder would shorten \code{Signature} and let distinct
#'   combinations collide.
#'
#' @return A \code{data.table}, ordered by block then descending \code{Count}, with
#'   columns:
#' \describe{
#'   \item{\code{Block}}{Block ID.}
#'   \item{\code{Chrom}}{Chromosome of the block.}
#'   \item{\code{Signature}}{The combination's alleles concatenated in canonical
#'     genomic order (no separator), with \code{fill_na} in place of an unobserved
#'     allele. Always exactly as wide as the block's marker count.}
#'   \item{\code{Markers}}{The block's member marker IDs, genomic order, joined by
#'     \code{";"} (identical for every row of a block).}
#'   \item{\code{Alleles}}{This combination's allele per member marker, joined by
#'     \code{";"} and index-aligned with \code{Markers}. Unobserved alleles appear
#'     as \code{fill_na}.}
#'   \item{\code{Count}}{Occurrence count per \code{count_by}.}
#'   \item{\code{Freq}}{\code{Count} divided by the block's total count (sums to 1
#'     within a block).}
#' }
#'
#' @examples
#' \dontrun{
#' meta <- load_geno_meta("my_geno_meta.rds")
#' sheet <- write_block_allele_sheet(meta, output_file = "block_allele_sheet.txt")
#' head(sheet)
#' }
#'
#' @seealso \code{\link{compute_block_catalogue}}, \code{\link{geno_meta_views}},
#'   \code{\link{compute_marker_stats}}
#' @export
write_block_allele_sheet <- function(geno_meta, output_file = NULL,
                                     blocks = NULL,
                                     count_by = c("haplotype", "genotype"),
                                     fill_na = NULL) {
  count_by <- match.arg(count_by)
  fill_na  <- .hs_fill_na(fill_na)
  .hs_check_meta(geno_meta, "write_block_allele_sheet", "geno_meta")

  views <- geno_meta_views(geno_meta)
  mm    <- views$marker_meta[!is.na(block_id)]        # already genomic-ordered
  if (!is.null(blocks)) mm <- mm[block_id %in% as.character(blocks)]

  empty <- data.table(Block = character(0L), Chrom = character(0L),
                      Signature = character(0L), Markers = character(0L),
                      Alleles = character(0L), Count = integer(0L),
                      Freq = numeric(0L))
  if (nrow(mm) == 0L) {
    if (!is.null(output_file)) data.table::fwrite(empty, output_file, sep = "\t")
    return(empty[])
  }

  alleles   <- views$alleles
  block_ids <- unique(mm$block_id)

  per_block <- lapply(block_ids, function(b) {
    bm    <- mm[block_id == b]
    mk    <- bm$marker[bm$marker %in% names(alleles)]  # member markers, in order
    if (length(mk) == 0L) return(NULL)
    chrom <- as.character(bm$chromosome[1L])

    # One row per haplotype copy; missing coded as fill_na so a copy is never
    # dropped, Count sums to n_samples x ploidy, and Signature keeps one
    # character per member marker (an empty fill would let distinct
    # combinations collide — see .compute_block_signatures()).
    sub <- data.table::copy(alleles[, c("Sample", mk), with = FALSE])
    for (col in mk) {
      v <- as.character(sub[[col]])
      v[is.na(v) | !nzchar(trimws(v))] <- fill_na
      data.table::set(sub, j = col, value = v)
    }

    if (count_by == "genotype")
      sub <- unique(sub, by = c("Sample", mk))         # a sample counts once/combo

    combos <- sub[, .(Count = .N), by = mk]
    combos[, Alleles   := do.call(paste, c(.SD, sep = ";")), .SDcols = mk]
    combos[, Signature := do.call(paste0, .SD),               .SDcols = mk]
    combos[, `:=`(Block   = as.character(b),
                  Chrom   = chrom,
                  Markers = paste(mk, collapse = ";"),
                  Freq    = Count / sum(Count))]
    combos[, .(Block, Chrom, Signature, Markers, Alleles, Count, Freq)]
  })

  out <- data.table::rbindlist(per_block, use.names = TRUE)
  if (nrow(out) == 0L) out <- empty
  else data.table::setorder(out, Block, -Count)

  if (!is.null(output_file))
    data.table::fwrite(out, output_file, sep = "\t")
  out[]
}


# ---------------------------------------------------------------------------
# Internal: encoding validation
# ---------------------------------------------------------------------------

#' Warn when effect_df alleles are incompatible with the profile encoding
#'
#' Q1 decision: strictly use the encoding specified by the profile; warn (not
#' error) when \code{effect_df} contains allele values outside the expected
#' set. The downstream join silently drops unmatched alleles; this warning
#' makes that silent exclusion visible.
#'
#' @param eff      Standardised effect data.table (Marker, Allele, Effect).
#' @param encoding Profile encoding string ("nucleotide", "numeric").
#' @param observed Alleles observed in the container (the core's frequency
#'   table). The nucleotide vocabulary is "categorical allele labels": a panel
#'   read from a PLINK 1/2-coded PED keeps its digit labels, so an effect
#'   table on those labels is valid — the join below is on the label.
#' @keywords internal
#' @noRd
.validate_effect_encoding <- function(eff, encoding, observed = NULL) {
  valid_alleles <- switch(encoding,
    "nucleotide"     = unique(c(names(DEFAULT_ENCODING_DICT),
                                as.character(observed))),
    "numeric"        = c("0", "1", "2"),
    NULL
  )
  if (is.null(valid_alleles)) return(invisible(NULL))

  bad <- setdiff(unique(eff$Allele), valid_alleles)
  if (length(bad) > 0L) {
    .hs_warn(
      "compute_profile_localgebv: effect_df contains allele(s) not in the ",
      "'", encoding, "' encoding dict: ",
      paste(bad, collapse = ", "), ". ",
      "These rows will not contribute to local GEBV. ",
      "Ensure effect_df uses the same encoding as the profile."
    )
  }
  invisible(NULL)
}
