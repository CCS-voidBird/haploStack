# =============================================================================
# block_subtypes.R
#
# count_block_subtypes()      — exported: ENUMERATE the subtypes present in one
#                               or more blocks (+ optional per-trait Effect).
# count_subtype_occurrence()  — exported: elemental LOOKUP — how often one named
#                               subtype occurs in one named block (0 if absent).
# count_subtype_carriage()    — exported: carriage of a SET of (block, subtype)
#                               targets, per sample plus pooled totals.
# count_favourable_subtypes() — exported: pipeline — derive the top subtype per
#                               block from a catalogue, then count its carriage.
# .as_subtype_targets()       — internal: normalise the target specification.
# .resolve_top_subtypes()     — internal: catalogue -> one target per block.
#
# The counting reuses .compute_block_signatures(); effects come wholly from
# compute_block_catalogue(); the "which subtype wins a block" rule is
# .select_block_subtype() (R/catalogue_analysis.R), so "favourable" means exactly
# what it means in analyze_block_catalogue(). No signature, allele-frequency,
# local-GEBV or selection logic is reimplemented here.
#
# Layering: count_subtype_occurrence() is the primitive and defines the contract;
# count_subtype_carriage() answers the same question for many targets in ONE
# keyed join rather than by looping it (a test pins the two to agree), because
# the per-call cost is geno_meta_views() and looping would pay it per target.
# =============================================================================

#' Count the haplotype subtypes carried in one or more blocks
#'
#' Searches a GenoMeta for the distinct \strong{allele combinations}
#' ("subtypes") observed in the named blocks and counts how often each occurs in
#' the population, optionally attaching the subtype's block effect for a chosen
#' trait. It answers the question "which haplotypes exist in this block, how
#' common are they, and are they good?" for a handful of blocks without running
#' a genome-wide catalogue.
#'
#' A subtype is identified by its \code{Subtype} signature: the alleles a
#' haplotype copy carries across the block's member markers, concatenated in
#' canonical genomic order. The signature is derived from genotype alone, so
#' \code{Block}, \code{Subtype}, \code{Count} and \code{Freq} never depend on
#' \code{trait}.
#'
#' @section Effects and cost:
#' \code{trait = NULL} (the default) leaves \code{Trait} and \code{Effect} as
#' \code{NA} and does no effect computation at all. Naming a \code{trait} calls
#' \code{\link{compute_block_catalogue}}, which has no block filter and so runs a
#' \strong{genome-wide} local-GEBV pass even when one block was requested. For
#' repeated interrogation of the same GenoMeta, either subset it first with
#' \code{\link{filter_geno_meta}} or reuse a pre-built signature table via
#' \code{.signatures}.
#'
#' @section Missing alleles:
#' An unobserved allele is written as \code{fill_na} (default \code{"-"}, see
#' \code{\link{HS_FILL_NA_DEFAULT}}) and keeps its slot, so every \code{Subtype}
#' is exactly as wide as its block and every haplotype copy is counted:
#' \code{sum(Count)} within a block is always \code{n_samples * ploidy} under
#' \code{count_by = "haplotype"}. A copy with nothing observed appears as an
#' all-\code{fill_na} subtype rather than vanishing.
#'
#' Filling rather than dropping is what keeps the signature a faithful
#' identifier: over markers \code{(M1, M2, M3)} the copies \code{(A, T, -)} and
#' \code{(A, -, T)} would both paste to \code{"AT"} if the gap were dropped, and
#' would be merged into one subtype. They are now \code{"AT-"} and \code{"A-T"}.
#'
#' Note that the container uncalls a genotype cell whenever \emph{any} copy of it
#' is missing, so on a diploid a single missing copy shows as \code{fill_na} on
#' both copies of that sample at that marker.
#' \code{\link{write_block_allele_sheet}} uses the same placeholder, so its
#' \code{Signature} column and this function's \code{Subtype} are the same string
#' and can be joined directly; use it when you also want the marker-by-marker
#' breakdown.
#'
#' @param geno_meta A valid \code{GenoMeta} with block assignments (see
#'   \code{\link{attach_block_info}}).
#' @param blocks Character vector of one or more block IDs to search. Block IDs
#'   not present raise a warning and are skipped. Default \code{NULL} (every
#'   block).
#' @param trait Optional trait name selecting the marker effects used for the
#'   \code{Effect} column. Default \code{NULL} — no effects are attached.
#' @param count_by How occurrences are counted. \code{"haplotype"} (default)
#'   counts haplotype copies, so a heterozygous carrier contributes to two
#'   subtypes and \code{sum(Count)} within a block is \code{n_samples * ploidy}.
#'   \code{"genotype"} counts distinct samples carrying the subtype on at least
#'   one copy, so a homozygote counts once.
#' @param center_by_freq Passed to \code{\link{compute_block_catalogue}} when
#'   \code{trait} is given; centres allele effects by population frequency before
#'   summing. Ignored when \code{trait} is \code{NULL}. Default \code{TRUE}.
#' @param output_file Optional path. When non-\code{NULL} the table is written
#'   there as a tab-separated file via \code{data.table::fwrite()}; the table is
#'   returned either way. Default \code{NULL}.
#' @param fill_na Placeholder written for an unobserved allele in
#'   \code{Subtype}. Default \code{NULL}, meaning the \code{haploStack.fill_na}
#'   option (itself defaulting to \code{"-"}). Ignored when \code{.signatures} is
#'   supplied, since that table is already built. Must be a single non-empty
#'   string.
#' @param .signatures Internal. A pre-built per-copy block-signature table
#'   (\code{.compute_block_signatures()}) to reuse across repeated calls. Default
#'   \code{NULL} (build internally).
#'
#' @return A \code{data.table}, one row per \code{(Block, Subtype)}, ordered by
#'   block then descending \code{Count}. The six columns are always present, so
#'   results of several calls bind together cleanly:
#' \describe{
#'   \item{\code{Trait}}{Trait whose effects were attached, or \code{NA} when
#'     \code{trait} was \code{NULL}.}
#'   \item{\code{Block}}{Block ID.}
#'   \item{\code{Subtype}}{The allele-combination signature — the copy's alleles
#'     over the block's member markers, concatenated in canonical genomic order
#'     with no separator.}
#'   \item{\code{Count}}{Occurrences, per \code{count_by}.}
#'   \item{\code{Freq}}{\code{Count} divided by the block's total count (sums to
#'     1 within a block). Under \code{count_by = "genotype"} the denominator is
#'     the number of sample-subtype pairs, not the number of samples.}
#'   \item{\code{Effect}}{The subtype's block local GEBV, taken from its
#'     representative carrier; \code{NA} when \code{trait} was \code{NULL}.}
#' }
#'
#' @examples
#' \dontrun{
#' meta <- load_geno_meta("my_geno_meta.rds")
#'
#' # Which haplotypes live in these two blocks, and how common are they?
#' count_block_subtypes(meta, blocks = c("Chr1_B3", "Chr2_B7"))
#'
#' # Same, with the yield effect of each subtype, counted per sample
#' count_block_subtypes(meta, blocks = "Chr1_B3", trait = "Yield",
#'                      count_by = "genotype")
#'
#' # Whole genome, written to file
#' count_block_subtypes(meta, output_file = "block_subtypes.txt")
#' }
#'
#' @seealso \code{\link{write_block_allele_sheet}},
#'   \code{\link{compute_block_catalogue}}, \code{\link{geno_meta_views}}
#' @export
count_block_subtypes <- function(geno_meta,
                                 blocks         = NULL,
                                 trait          = NULL,
                                 count_by       = c("haplotype", "genotype"),
                                 center_by_freq = TRUE,
                                 output_file    = NULL,
                                 fill_na        = NULL,
                                 .signatures    = NULL) {

  count_by <- match.arg(count_by)
  if (!is_geno_meta(geno_meta))
    stop("count_block_subtypes: 'geno_meta' is not a valid GenoMeta.")

  empty <- data.table(Trait   = character(0L), Block = character(0L),
                      Subtype = character(0L), Count = integer(0L),
                      Freq    = numeric(0L),   Effect = numeric(0L))

  # Trait-invariant signatures: the single source for both the counts and the
  # key the effects are joined on (see the "Missing alleles" section).
  # With no trait the signatures are all we need, so restrict the melt to the
  # requested blocks. With a trait, compute_block_catalogue() below needs a table
  # covering every block (it labels them all), so build the full one.
  sigs <- .signatures
  if (is.null(sigs))
    sigs <- .compute_block_signatures(geno_meta, fill_na,
                                      blocks = if (is.null(trait)) blocks else NULL)

  # The absent-block warning is checked against the GenoMeta's own block
  # universe, NOT against `sigs`: when every requested block is missing the
  # restricted signature table comes back empty, and keying the warning off it
  # would silently skip the very case the caller most needs to hear about.
  req <- if (is.null(blocks)) NULL else unique(as.character(blocks))
  if (!is.null(req)) {
    known  <- unique(as.character(geno_meta$block_marker_map_info$Block))
    absent <- setdiff(req, known)
    if (length(absent) > 0L)
      warning("count_block_subtypes: ", length(absent),
              " block(s) not found in this GenoMeta: ",
              paste(utils::head(absent, 5L), collapse = ", "),
              if (length(absent) > 5L) ", ..." else "")
  }

  sel <- if (is.null(req)) sigs else sigs[Block %in% req]
  if (nrow(sel) == 0L) {
    if (!is.null(output_file)) data.table::fwrite(empty, output_file, sep = "\t")
    return(empty[])
  }

  cnt <- if (count_by == "haplotype")
           sel[, .(Count = .N), by = .(Block, Subtype = signature)]
         else
           sel[, .(Count = data.table::uniqueN(Sample)),
               by = .(Block, Subtype = signature)]
  cnt[, Count := as.integer(Count)]
  cnt[, Freq  := Count / sum(Count), by = Block]

  if (is.null(trait)) {
    cnt[, `:=`(Trait = NA_character_, Effect = NA_real_)]
    me <- geno_meta$marker_effects
    if (!is.null(me)) {
      avail <- names(.as_effect_list(me))
      message("count_block_subtypes: effects available for trait(s) ",
              paste(avail, collapse = ", "),
              " — pass trait= to attach an Effect column.")
    }
  } else {
    # The UNFILTERED signature table is injected: compute_block_catalogue()
    # labels every block, and a pre-filtered table would leave the rest NA.
    cat_dt <- compute_block_catalogue(geno_meta, trait = trait,
                                      center_by_freq = center_by_freq,
                                      group_by       = "signature",
                                      .signatures    = sigs)
    cnt[, Effect := NA_real_]
    cnt[cat_dt, on = c("Block", "Subtype" = "signature"), Effect := i.Effect]
    cnt[, Trait := as.character(trait)]
  }

  .lead_identifier_cols(cnt)
  data.table::setcolorder(cnt,
    c("Trait", "Block", "Subtype", "Count", "Freq", "Effect"))
  data.table::setorder(cnt, Block, -Count)

  if (!is.null(output_file))
    data.table::fwrite(cnt, output_file, sep = "\t")
  cnt[]
}


# ---------------------------------------------------------------------------
# Exported: elemental lookup — one block, one subtype
# ---------------------------------------------------------------------------

#' How often does one named subtype occur in one named block?
#'
#' The elemental subtype lookup: given a single block ID and a single signature
#' string, count how many times that exact haplotype occurs in the GenoMeta. A
#' subtype that is simply not present returns \code{0} rather than raising an
#' error, so the function is safe to call over candidate signatures that may or
#' may not segregate in the panel.
#'
#' Signatures must be written the way the package builds them — alleles in
#' canonical genomic order with \code{fill_na} (default \code{"-"}) standing in
#' for an unobserved allele. Take them from \code{\link{count_block_subtypes}},
#' \code{\link{write_block_allele_sheet}} or
#' \code{\link{compute_block_catalogue}}; a hand-written string of the wrong
#' width will legitimately count \code{0}.
#'
#' @section Cost:
#' Only the block's member markers are read, but every call still rebuilds the
#' GenoMeta's allele view, and that dominates. For more than a couple of lookups
#' use \code{\link{count_subtype_carriage}}, which answers the same question for
#' many targets in one pass, or build the signature table once with
#' \code{haploStack:::.compute_block_signatures()} and pass it as
#' \code{.signatures}.
#'
#' @param geno_meta A valid \code{GenoMeta} with block assignments.
#' @param block A single block ID. A block absent from the GenoMeta warns and
#'   returns \code{0}.
#' @param subtype A single signature string.
#' @param count_by \code{"haplotype"} (default) counts haplotype copies;
#'   \code{"genotype"} counts distinct samples carrying the subtype on at least
#'   one copy.
#' @param fill_na Placeholder for an unobserved allele when signatures are built.
#'   Default \code{NULL} (the \code{haploStack.fill_na} option, \code{"-"}).
#'   Ignored when \code{.signatures} is supplied.
#' @param .signatures Optional pre-built per-copy signature table from
#'   \code{.compute_block_signatures()}, reused across calls. Default
#'   \code{NULL} (build one restricted to this block).
#'
#' @return A single non-negative \code{integer}. \code{0} when the subtype does
#'   not occur, and \code{0} with a warning when \code{block} is not in the
#'   GenoMeta.
#'
#' @examples
#' \dontrun{
#' meta <- load_geno_meta("my_geno_meta.rds")
#'
#' # What does this block hold, and how common is its commonest subtype?
#' tab <- count_block_subtypes(meta, blocks = "Chr1_B3")
#' count_subtype_occurrence(meta, "Chr1_B3", tab$Subtype[1])
#'
#' # A subtype that does not segregate here simply counts 0
#' count_subtype_occurrence(meta, "Chr1_B3", "ACGT")
#' }
#'
#' @seealso \code{\link{count_subtype_carriage}} for a set of targets,
#'   \code{\link{count_block_subtypes}} to enumerate what a block holds,
#'   \code{\link{count_favourable_subtypes}} for the catalogue-driven pipeline
#' @export
count_subtype_occurrence <- function(geno_meta, block, subtype,
                                     count_by    = c("haplotype", "genotype"),
                                     fill_na     = NULL,
                                     .signatures = NULL) {

  count_by <- match.arg(count_by)
  if (!is_geno_meta(geno_meta))
    stop("count_subtype_occurrence: 'geno_meta' is not a valid GenoMeta.")

  block   <- as.character(block)
  subtype <- as.character(subtype)
  if (length(block) != 1L || is.na(block) || !nzchar(block))
    stop("count_subtype_occurrence: 'block' must be a single non-empty block ID; ",
         "use count_subtype_carriage() for several blocks.")
  if (length(subtype) != 1L || is.na(subtype))
    stop("count_subtype_occurrence: 'subtype' must be a single non-NA signature ",
         "string; use count_subtype_carriage() for several subtypes.")

  sigs <- .signatures
  if (is.null(sigs))
    sigs <- .compute_block_signatures(geno_meta, fill_na, blocks = block)

  in_block <- sigs[Block == block]
  if (nrow(in_block) == 0L) {
    warning("count_subtype_occurrence: block '", block,
            "' not found in this GenoMeta; returning 0.")
    return(0L)
  }

  hit <- in_block[signature == subtype]
  if (nrow(hit) == 0L) return(0L)          # not observed — a valid answer
  if (count_by == "haplotype") nrow(hit) else data.table::uniqueN(hit$Sample)
}


# ---------------------------------------------------------------------------
# Internal: normalise a (Block, Subtype) target specification
# ---------------------------------------------------------------------------

# Accepts a data.frame carrying `Block` plus one of `Subtype` / `signature` /
# `Signature` — so a table straight out of count_block_subtypes(),
# write_block_allele_sheet() or compute_block_catalogue() can be handed over
# unchanged — or a named character vector whose names are the block IDs.
# Returns data.table(Block, Subtype); duplicates are dropped with a warning.
# @keywords internal
.as_subtype_targets <- function(subtypes, caller) {
  if (is.data.frame(subtypes)) {
    dt <- data.table::as.data.table(subtypes)
    if (!"Block" %in% names(dt))
      stop(caller, ": 'subtypes' data.frame needs a 'Block' column.")
    scol <- intersect(c("Subtype", "signature", "Signature"), names(dt))
    if (length(scol) == 0L)
      stop(caller, ": 'subtypes' data.frame needs one of 'Subtype', ",
           "'signature' or 'Signature'.")
    out <- data.table::data.table(Block   = as.character(dt[["Block"]]),
                                  Subtype = as.character(dt[[scol[1L]]]))
  } else if (is.character(subtypes) || is.factor(subtypes)) {
    nm <- names(subtypes)
    if (is.null(nm) || any(!nzchar(nm)))
      stop(caller, ": a vector of 'subtypes' must be fully named, with the ",
           "block ID as each element's name.")
    out <- data.table::data.table(Block = nm, Subtype = as.character(subtypes))
  } else {
    stop(caller, ": 'subtypes' must be a data.frame with 'Block' plus ",
         "'Subtype'/'signature'/'Signature', or a named character vector ",
         "(names = block IDs).")
  }

  out <- out[!is.na(Block) & nzchar(Block) & !is.na(Subtype)]
  if (nrow(out) == 0L)
    stop(caller, ": no usable (Block, Subtype) targets were supplied.")

  n_dup <- sum(duplicated(out))
  if (n_dup > 0L) {
    warning(caller, ": ", n_dup,
            " duplicate (Block, Subtype) target(s) dropped.")
    out <- unique(out)
  }
  out[]
}


# ---------------------------------------------------------------------------
# Exported: carriage of a set of targets
# ---------------------------------------------------------------------------

#' Count carriage of a set of block subtypes
#'
#' Given a set of \code{(Block, Subtype)} targets, report how many of them each
#' sample carries. This is \code{\link{count_subtype_occurrence}} asked of many
#' targets at once — the usual question being "how many of these favourable
#' haplotypes does this line have?", so it is typically pointed at a sub-GenoMeta
#' holding a single sample. A multi-sample GenoMeta returns one row per sample,
#' which ranks lines by carriage without a wrapper loop.
#'
#' A target whose subtype does not occur contributes a \code{Count} of \code{0}
#' rather than dropping out, so \code{detail} always holds exactly
#' \code{n_target * n_samples} rows and the denominator of any fraction stays
#' explicit.
#'
#' @param geno_meta A valid \code{GenoMeta} with block assignments.
#' @param subtypes The targets: a data.frame with a \code{Block} column plus one
#'   of \code{Subtype} / \code{signature} / \code{Signature} (so the output of
#'   \code{\link{count_block_subtypes}}, \code{\link{write_block_allele_sheet}} or
#'   \code{\link{compute_block_catalogue}} passes straight in), or a named
#'   character vector whose names are block IDs. Duplicates are dropped with a
#'   warning; target blocks absent from the GenoMeta warn and count \code{0}.
#' @param count_by \code{"haplotype"} (default) counts matching haplotype copies
#'   per sample, so \code{Count} runs \code{0}..\code{ploidy}.
#'   \code{"genotype"} reduces each sample to carried / not carried (\code{Count}
#'   is \code{0} or \code{1}), so the pooled total is the number of carrying
#'   samples.
#' @param fill_na Placeholder for an unobserved allele when signatures are built.
#'   Default \code{NULL} (the \code{haploStack.fill_na} option). Ignored when
#'   \code{.signatures} is supplied.
#' @param .signatures Optional pre-built signature table to reuse. Default
#'   \code{NULL}.
#'
#' @return A named \code{list}:
#' \describe{
#'   \item{\code{n_target}}{Number of distinct \code{(Block, Subtype)} targets.}
#'   \item{\code{n_samples}}{Number of samples in the GenoMeta.}
#'   \item{\code{n_carried}}{Targets carried, pooled over samples — the headline
#'     number, and for a single-sample GenoMeta the answer to "how many of these
#'     does this line carry?".}
#'   \item{\code{total}}{Sum of \code{Count} over every sample and target
#'     (matching copies under \code{"haplotype"}).}
#'   \item{\code{by_sample}}{\code{data.table(Sample, nTarget, nCarried, Total,
#'     Frac)}, ordered by descending \code{nCarried}; \code{Frac} is
#'     \code{nCarried / nTarget}.}
#'   \item{\code{detail}}{\code{data.table(Block, Sample, Subtype, Count,
#'     Carried)} — one row per sample x target, zeros included.}
#'   \item{\code{count_by}}{The counting mode used.}
#' }
#'
#' @examples
#' \dontrun{
#' meta <- load_geno_meta("my_geno_meta.rds")
#'
#' # Targets as a named vector ...
#' res <- count_subtype_carriage(meta, c(Chr1_B3 = "4441", Chr2_B7 = "1114"))
#' res$n_carried
#'
#' # ... or any table carrying Block plus a signature column
#' top <- count_block_subtypes(meta, trait = "Yield")[, .SD[1L], by = Block]
#' count_subtype_carriage(meta, top)$by_sample
#' }
#'
#' @seealso \code{\link{count_subtype_occurrence}} (the primitive),
#'   \code{\link{count_favourable_subtypes}} (derives the targets from a
#'   catalogue), \code{\link{count_block_subtypes}}
#' @export
count_subtype_carriage <- function(geno_meta, subtypes,
                                   count_by    = c("haplotype", "genotype"),
                                   fill_na     = NULL,
                                   .signatures = NULL) {

  count_by <- match.arg(count_by)
  if (!is_geno_meta(geno_meta))
    stop("count_subtype_carriage: 'geno_meta' is not a valid GenoMeta.")

  tgt     <- .as_subtype_targets(subtypes, "count_subtype_carriage")
  samples <- dimnames(geno_meta$geno_arrays)[[1L]]

  sigs <- .signatures
  if (is.null(sigs))
    sigs <- .compute_block_signatures(geno_meta, fill_na,
                                      blocks = unique(tgt$Block))

  absent <- setdiff(unique(tgt$Block), unique(sigs$Block))
  if (length(absent) > 0L)
    warning("count_subtype_carriage: ", length(absent),
            " target block(s) not found in this GenoMeta (counted as 0): ",
            paste(utils::head(absent, 5L), collapse = ", "),
            if (length(absent) > 5L) ", ..." else "")

  # Full sample x target grid, so an uncarried target stays visible as a zero.
  grid <- data.table::CJ(Sample = samples, Target = seq_len(nrow(tgt)),
                         sorted = FALSE)
  grid[, `:=`(Block = tgt$Block[Target], Subtype = tgt$Subtype[Target])]
  grid[, Target := NULL]

  # One pass over the signatures: matching copies per (Sample, Block, Subtype).
  obs <- sigs[signature %in% unique(tgt$Subtype),
              .(Copies = .N), by = .(Sample, Block, Subtype = signature)]

  detail <- merge(grid, obs, by = c("Sample", "Block", "Subtype"), all.x = TRUE)
  detail[is.na(Copies), Copies := 0L]
  detail[, Count   := if (count_by == "haplotype") as.integer(Copies)
                      else as.integer(Copies > 0L)]
  detail[, Carried := Count > 0L]
  detail[, Copies  := NULL]

  by_sample <- detail[, .(nTarget  = .N,
                          nCarried = sum(Carried),
                          Total    = sum(Count)), by = Sample]
  by_sample[, Frac := nCarried / nTarget]
  data.table::setorder(by_sample, -nCarried, Sample)

  .lead_identifier_cols(detail)
  data.table::setorder(detail, Block, Sample)

  list(n_target  = nrow(tgt),
       n_samples = length(samples),
       n_carried = sum(detail$Carried),
       total     = sum(detail$Count),
       by_sample = by_sample[],
       detail    = detail[],
       count_by  = count_by)
}


# ---------------------------------------------------------------------------
# Internal: catalogue -> one favourable subtype per block
# ---------------------------------------------------------------------------

# Resolve a catalogue, in any of the shapes the package produces, to
# data.table(Block, Subtype, Effect) with exactly one row per block. Accepts:
#   * an analyze_block_catalogue() result  — `selected` already holds the picks,
#     so its choice is honoured verbatim and Effect is summed over traits from `S`
#   * a run_haplotype_catalogue()-shaped named list (per-trait `$hap_catalogue`)
#   * a single catalogue data.frame (Block + signature/Subtype + Effect [+ Count])
#   * NULL — compute one from geno_meta + trait
# The winner-per-block rule is always .select_block_subtype() (catalogue_analysis.R),
# so "favourable" means precisely what it means in analyze_block_catalogue().
# @keywords internal
.resolve_top_subtypes <- function(catalogue, geno_meta, trait, fitness_fn,
                                  min_count, center_by_freq, fill_na,
                                  caller = "count_favourable_subtypes") {

  # ── An analysis result: the blocks are already chosen. ────────────────────
  if (is.list(catalogue) && !is.data.frame(catalogue) &&
      !is.null(catalogue$selected)) {
    sel <- data.table::as.data.table(catalogue$selected)
    if (!all(c("Block", "Subtype") %in% names(sel)))
      stop(caller, ": the analysis result's 'selected' table needs 'Block' and ",
           "'Subtype' columns.")
    ref <- unique(sel[, .(Block = as.character(Block),
                          Subtype = as.character(Subtype))], by = "Block")
    if (!is.null(catalogue$S)) {
      s <- data.table::as.data.table(catalogue$S)
      eff <- s[, .(Effect = sum(as.numeric(S), na.rm = TRUE)),
               by = .(Block = as.character(Block), Subtype = as.character(Subtype))]
      ref <- merge(ref, eff, by = c("Block", "Subtype"), all.x = TRUE)
    } else {
      ref[, Effect := NA_real_]
    }
    return(ref[])
  }

  # ── Otherwise assemble a long (Block, Subtype, Count, Trait, Effect) table. ─
  if (is.data.frame(catalogue)) {
    long <- data.table::as.data.table(catalogue)
    if (!"Trait" %in% names(long)) long[, Trait := "Trait"]
  } else if (is.list(catalogue) && !is.null(catalogue)) {
    nm <- names(catalogue)
    if (is.null(nm) || any(!nzchar(nm)))
      stop(caller, ": a list 'catalogue' must be named (one element per trait).")
    long <- data.table::rbindlist(lapply(nm, function(tr) {
      el <- catalogue[[tr]]
      hc <- if (is.data.frame(el)) el else el$hap_catalogue
      if (is.null(hc))
        stop(caller, ": catalogue element '", tr, "' has no 'hap_catalogue'.")
      d <- data.table::as.data.table(hc)
      d[, Trait := tr]
      d
    }), fill = TRUE)
  } else {
    if (is.null(geno_meta))
      stop(caller, ": supply a 'catalogue', or a 'geno_meta' carrying effects.")
    long <- compute_block_catalogue(geno_meta, trait = trait,
                                    center_by_freq = center_by_freq,
                                    group_by = "signature")
    long <- data.table::as.data.table(long)
    long[, Trait := as.character(trait %||% "Trait")]
  }

  scol <- intersect(c("Subtype", "signature", "Signature"), names(long))
  if (length(scol) == 0L)
    stop(caller, ": the catalogue needs a subtype column ('Subtype', ",
         "'signature' or 'Signature').")
  if (!"Effect" %in% names(long))
    stop(caller, ": the catalogue needs an 'Effect' column to rank subtypes by.")
  if (!"Count" %in% names(long)) long[, Count := 1L]

  # .select_block_subtype() expects Block/Subtype/Count/Trait plus Vbase, and a
  # named per-trait weight vector. Equal unit weights: a constant scaling of
  # every weight leaves the argmax unchanged.
  st <- long[, .(Block   = as.character(Block),
                 Subtype = as.character(get(scol[1L])),
                 Count   = as.numeric(Count),
                 Trait   = as.character(Trait),
                 Vbase   = as.numeric(Effect))]
  st <- st[!is.na(Block) & !is.na(Subtype)]
  if (nrow(st) == 0L) stop(caller, ": the catalogue holds no usable subtypes.")

  traits <- unique(st$Trait)
  w      <- stats::setNames(rep(1, length(traits)), traits)
  pick   <- .select_block_subtype(st, fitness_fn, min_count, w)

  eff <- st[, .(Effect = sum(Vbase, na.rm = TRUE)), by = .(Block, Subtype)]
  ref <- merge(pick, eff, by = c("Block", "Subtype"), all.x = TRUE)
  unique(ref, by = "Block")[]
}


# ---------------------------------------------------------------------------
# Exported: favourable-subtype carriage pipeline
# ---------------------------------------------------------------------------

#' Count how many of a catalogue's favourable block subtypes a panel carries
#'
#' A small pipeline over the two functions above. It takes a haplotype catalogue,
#' picks the \strong{one most favourable subtype per block} (by default the
#' highest block effect), and then counts how many of those a GenoMeta's samples
#' actually carry. Pointed at a single-sample sub-GenoMeta it answers "how many
#' favourable haplotypes does this line already have?"; pointed at a panel it
#' ranks the lines.
#'
#' The winner-per-block rule is \code{\link{analyze_block_catalogue}}'s own
#' selector, so \code{fitness_fn} takes the same values and \dQuote{favourable}
#' means the same thing in both places. When the catalogue covers several traits
#' the per-trait effects of a subtype are summed with equal weight before the
#' block's winner is chosen; pass a pre-computed
#' \code{\link{analyze_block_catalogue}} result to control that with
#' \code{trait_weights}, and its \code{selected} picks are then honoured verbatim.
#'
#' @param geno_meta A valid \code{GenoMeta} with block assignments — the panel
#'   whose carriage is counted (often a one-sample subset).
#' @param catalogue The reference. Any of: a catalogue table from
#'   \code{\link{compute_block_catalogue}(group_by = "signature")}; a named list
#'   of per-trait catalogues as \code{\link{run_haplotype_catalogue}} returns; an
#'   \code{\link{analyze_block_catalogue}} result, whose \code{selected} picks are
#'   used as they stand. \code{NULL} (default) computes a catalogue from
#'   \code{geno_meta} and \code{trait}.
#' @param trait Trait name, used only when \code{catalogue} is \code{NULL}.
#' @param fitness_fn How each block's winner is chosen, passed to the shared
#'   selector: \code{"max_effect"} (default — highest effect),
#'   \code{"frequent_best"} (best among subtypes with \code{Count >= min_count}),
#'   \code{"expected"} (frequency-weighted effect), \code{"min_effect"}
#'   (lowest, for counting unfavourable haplotypes), or a function taking a
#'   block's subtype rows and returning exactly one.
#' @param min_count Minimum subtype count for \code{"frequent_best"}. Default 1.
#' @param top_n Optionally keep only the \code{top_n} strongest blocks, ranked by
#'   the selected subtype's effect (ascending under \code{"min_effect"}). Default
#'   \code{NULL} (every block in the catalogue).
#' @param count_by \code{"haplotype"} (default) counts matching copies per
#'   sample; \code{"genotype"} counts each sample once per target carried. See
#'   \code{\link{count_subtype_carriage}}.
#' @param center_by_freq Passed to \code{\link{compute_block_catalogue}} when the
#'   catalogue is computed here. Default \code{TRUE}.
#' @param fill_na Placeholder for an unobserved allele. Default \code{NULL} (the
#'   \code{haploStack.fill_na} option).
#' @param output_file Optional path; writes the \code{by_sample} table as a
#'   tab-separated file. The full result is returned either way. Default
#'   \code{NULL}.
#'
#' @return The \code{\link{count_subtype_carriage}} list, plus:
#' \describe{
#'   \item{\code{reference}}{\code{data.table(Block, Subtype, Effect)} — the
#'     favourable target chosen for each block, ordered by descending
#'     \code{Effect}.}
#'   \item{\code{detail}}{as \code{\link{count_subtype_carriage}}, with the
#'     target's \code{Effect} joined on.}
#'   \item{\code{fitness_fn}}{The rule used, as a string (\code{"custom"} for a
#'     function).}
#' }
#'
#' @examples
#' \dontrun{
#' meta <- load_geno_meta("my_geno_meta.rds")
#' cat_dt <- compute_block_catalogue(meta, trait = "Yield",
#'                                   group_by = "signature")
#'
#' # How many favourable haplotypes does one line carry?
#' line <- filter_geno_meta(meta, keep_samples = "L0021")
#' res  <- count_favourable_subtypes(line, cat_dt)
#' res$n_carried
#' res$detail[Carried == TRUE]
#'
#' # Rank the whole panel on the 50 strongest blocks
#' count_favourable_subtypes(meta, cat_dt, top_n = 50)$by_sample
#' }
#'
#' @seealso \code{\link{count_subtype_carriage}},
#'   \code{\link{count_subtype_occurrence}},
#'   \code{\link{analyze_block_catalogue}},
#'   \code{\link{compute_block_catalogue}}
#' @export
count_favourable_subtypes <- function(geno_meta,
                                      catalogue      = NULL,
                                      trait          = NULL,
                                      fitness_fn     = c("max_effect",
                                                         "frequent_best",
                                                         "expected",
                                                         "min_effect"),
                                      min_count      = 1,
                                      top_n          = NULL,
                                      count_by       = c("haplotype", "genotype"),
                                      center_by_freq = TRUE,
                                      fill_na        = NULL,
                                      output_file    = NULL) {

  if (!is.function(fitness_fn)) fitness_fn <- match.arg(fitness_fn)
  count_by <- match.arg(count_by)
  if (!is_geno_meta(geno_meta))
    stop("count_favourable_subtypes: 'geno_meta' is not a valid GenoMeta.")

  ref <- .resolve_top_subtypes(catalogue, geno_meta, trait, fitness_fn,
                               min_count, center_by_freq, fill_na)

  # Rank blocks by the selected subtype's effect, following the fitness
  # direction, so top_n keeps the strongest under either sense of "strong".
  if (identical(fitness_fn, "min_effect")) data.table::setorder(ref, Effect,  Block)
  else                                     data.table::setorder(ref, -Effect, Block)

  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)[1L]
    if (is.na(top_n) || top_n < 1L)
      stop("count_favourable_subtypes: 'top_n' must be a positive integer.")
    if (top_n < nrow(ref)) ref <- ref[seq_len(top_n)]
  }

  res <- count_subtype_carriage(geno_meta, ref, count_by = count_by,
                               fill_na = fill_na)

  # Carry each target's effect into the per-row detail.
  res$detail <- merge(res$detail, ref[, .(Block, Subtype, Effect)],
                      by = c("Block", "Subtype"), all.x = TRUE)
  .lead_identifier_cols(res$detail)
  data.table::setorder(res$detail, Block, Sample)

  res$reference  <- ref[]
  res$fitness_fn <- if (is.function(fitness_fn)) "custom" else fitness_fn

  if (!is.null(output_file)) {
    data.table::fwrite(res$by_sample, output_file, sep = "\t")
    message("Favourable-subtype carriage written to: ", output_file)
  }
  res
}
