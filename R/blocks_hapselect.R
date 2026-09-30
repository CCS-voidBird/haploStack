# =============================================================================
# blocks_hapselect.R — hapS_ld_blocks(): LD haplotype blocks via HapSelect
#
# The default local resolver (R/blocks_ld_common.R says when a resolver is
# needed at all: only when no block map exists). Promoted in 0.7.10 from the
# dev-tree bridge dev/hapselect_bridge/hapselect_bridge.R (thread #028), whose
# LD half it reproduces call for call:
#
#   1. one chromosome key, character Chrom <-> the integer HapSelect requires
#      (its check_ld_matrix() refuses a character column, and letting it coerce
#      by itself numbers the geno and map frames independently);
#   2. the dosage of one allele per marker (the minor one — r^2 does not depend
#      on the choice), NA when any ploidy copy is missing;
#   3. markers sorted by (chromosome, position): pairwise_ld() takes ROW ORDER
#      as genomic order and def_blocks() calls Locus2 == Locus1 + 1 adjacent,
#      so unsorted input does not error, it blocks the wrong pairs;
#   4. HapSelect::pairwise_ld() + def_blocks() + block_obj_to_df(), and the
#      one-row-per-block table ("m1;m2;...", Block_ID "<key>:<n>") mapped back
#      to (Chrom, Marker, Block) with the container's chromosome names.
#
# On the same markers it returns the bridge's partition and block IDs; unlike
# the bridge it leaves unplaced markers out (they would be blocked by file
# order) and keeps effects / phenotype / FAM on attach.
#
# HapSelect (github.com/wrshf7/HapSelect, GPL-3) is in Suggests: haploStack
# installs and runs without it. Two of its habits are contained here: with
# parallelize = TRUE pairwise_ld() replaces the user's future plan and resets
# it to sequential, and both functions leave progressr's handlers option set.
# =============================================================================


# Fail early, with the install route and the alternatives.
.hapS_require <- function(caller) {
  if (!requireNamespace("HapSelect", quietly = TRUE))
    stop(caller, " requires the HapSelect package, which is not installed:\n",
         "  remotes::install_github(\"wrshf7/HapSelect\")   ",
         "# builds C++: Rtools on Windows\n",
         "Or skip local LD blocking — the recommended route: give a block map ",
         "to read_plink_ped(block_file =), build_geno_meta(blocks =) or ",
         "attach_block_info().", call. = FALSE)
  invisible(TRUE)
}

# Character chromosome -> integer key, as the bridge builds it (so block IDs
# match every HapSelect result made through it). Chromosomes that parse as
# whole numbers keep their value; otherwise they are numbered by the digits
# they carry ("1H" -> 1), falling back to plain sort order when some carry none.
# returns data.frame(Chrom, Chromosome), ordered by Chromosome.
.hapS_chrom_key <- function(chroms) {
  chroms <- unique(as.character(chroms))
  num_direct <- suppressWarnings(as.numeric(chroms))
  if (!anyNA(num_direct) && all(num_direct == round(num_direct))) {
    codes <- as.integer(num_direct)
  } else {
    digits <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", chroms)))
    ord    <- if (anyNA(digits)) order(chroms) else order(digits, chroms)
    codes  <- integer(length(chroms))
    codes[ord] <- seq_along(chroms)
  }
  key <- data.frame(Chrom = chroms, Chromosome = as.integer(codes),
                    stringsAsFactors = FALSE)
  key[order(key$Chromosome), , drop = FALSE]
}

# The genotype frame HapSelect reads: SNP, Chromosome, Position, then one
# dosage column per sample; rows in genomic order. `markers` are the placed
# markers. Also returns how many alleles each marker carries.
.hapS_dosage_frame <- function(geno_meta, markers, chrom_key) {
  arr     <- geno_meta$geno_arrays[, markers, , drop = FALSE]
  samples <- dimnames(arr)[[1L]]
  n_s <- dim(arr)[1L]; n_m <- dim(arr)[2L]; ploidy <- dim(arr)[3L]
  miss <- .missing_tokens(geno_meta$metadata$encoding %||%
                            geno_meta$metadata$source_format)

  # Allele counts per marker; the counted allele is the least frequent one
  # observed, ties broken alphabetically (so it never equals the other allele
  # of a perfectly balanced marker).
  levs   <- sort(setdiff(unique(as.vector(arr)), c(NA_character_, miss)))
  counts <- matrix(0L, nrow = n_m, ncol = length(levs))
  for (k in seq_len(ploidy)) {
    layer <- matrix(arr[, , k], nrow = n_s, ncol = n_m)
    for (li in seq_along(levs))
      counts[, li] <- counts[, li] + colSums(layer == levs[li], na.rm = TRUE)
  }
  counted <- rep(NA_character_, n_m)
  n_all   <- integer(n_m)
  for (j in seq_len(n_m)) {
    obs <- which(counts[j, ] > 0L)
    n_all[j] <- length(obs)
    if (length(obs) == 0L) next
    ranked <- levs[obs][order(-counts[j, obs])]
    counted[j] <- ranked[length(ranked)]
  }

  dose  <- matrix(0, nrow = n_s, ncol = n_m)
  na_ct <- matrix(FALSE, nrow = n_s, ncol = n_m)
  for (k in seq_len(ploidy)) {
    layer    <- matrix(arr[, , k], nrow = n_s, ncol = n_m)
    layer_na <- is.na(layer) | layer %in% miss
    na_ct    <- na_ct | layer_na
    hit      <- !layer_na & layer == counted[col(layer)]
    hit[is.na(hit)] <- FALSE          # a marker with no observed allele
    dose     <- dose + hit
  }
  dose[na_ct] <- NA_real_

  map <- geno_meta$map_info
  mi  <- match(markers, map$Marker)
  geno <- data.frame(
    SNP        = markers,
    Chromosome = chrom_key$Chromosome[match(as.character(map$Chrom[mi]),
                                            chrom_key$Chrom)],
    Position   = as.numeric(map$Pos[mi]),
    stringsAsFactors = FALSE)
  dose_df <- as.data.frame(t(dose), stringsAsFactors = FALSE)
  names(dose_df) <- samples
  geno <- cbind(geno, dose_df)
  ord  <- order(geno$Chromosome, geno$Position)
  geno <- geno[ord, , drop = FALSE]
  rownames(geno) <- geno$SNP
  list(geno = geno, n_alleles = n_all[ord])
}


#' LD haplotype blocks with HapSelect (the default local resolver)
#'
#' \strong{A side route.} haploStack does not need LD blocking: the recommended
#' input is a block map you already have, given to
#' \code{\link{read_plink_ped}} (\code{block_file =}),
#' \code{\link{build_geno_meta}} (\code{blocks =}) or
#' \code{\link{attach_block_info}} — no LD parameter involved. Use this
#' function only when blocks must be built locally. It is the default of the two
#' resolvers; \code{\link{st_ld_blocks}} (SelectionTools) is the experimental
#' alternative, and \strong{the two do not partition alike} at the same
#' settings — compare block sets only within one backend.
#'
#' It turns the container into HapSelect's inputs (one allele's dosage per
#' marker, markers in genomic order, an integer chromosome key), runs
#' \code{HapSelect::pairwise_ld()} and \code{HapSelect::def_blocks()}, maps the
#' blocks back to the container's marker IDs and chromosome names, and by
#' default attaches them. It reproduces, call for call, the LD half of the
#' dev-tree bridge \code{dev/hapselect_bridge/hapselect_bridge.R} that the
#' workshop and the repository's dataset segmentations used before 0.7.10.
#'
#' @section HapSelect's defaults are not used:
#' \code{HapSelect::def_blocks()} defaults to \code{threshold = 0.7,
#' tolerance = 1}. This function passes the whole rule itself, with defaults of
#' \eqn{r^2 \ge 0.5}, tolerance 3, \code{tol_reset = TRUE}, the
#' \code{"flanking"} method and \code{start = "LD"}: the settings of the
#' repository's segmented datasets and of the workshop.
#'
#' @section Markers left out, and markers without variation:
#' Markers without a usable position (non-finite or \eqn{\le} 0) or on an
#' unplaced chromosome (\code{0}, \code{UN}, \code{unknown}) are not blocked,
#' because HapSelect reads row order as genomic order; on attach they become
#' singleton blocks, as any marker absent from a block map does.
#'
#' A marker whose genotype never varies in the panel has a constant dosage:
#' monomorphic, or every line heterozygous (in an inbred panel usually a
#' collapsed paralog, a QC finding). Its \eqn{r^2} is \code{NA} with every
#' other marker, and \code{pairwise_ld()} warns \emph{"the standard deviation is
#' zero"} once per such pair — thousands of times on a real panel. Those
#' warnings are muffled here and the markers counted in the summary. HapSelect
#' does not make them singletons: a run of them becomes a block of its own, or,
#' with enough tolerance, part of the neighbouring block.
#'
#' @section Cost:
#' \code{pairwise_ld()} is pure R and quadratic in the markers per chromosome:
#' about 1.6 min for 1 000 markers x 300 lines, about 14 min for the 9 856-marker
#' GRDC barley panel. Rehearse on one chromosome. \code{parallelize = TRUE} uses
#' HapSelect's own \pkg{future} parallelism; the caller's \pkg{future} plan is
#' restored afterwards.
#'
#' @param geno_meta    A valid \code{GenoMeta} with at least 2 samples.
#' @param ld_threshold \eqn{r^2} at which flanking markers stay in one block;
#'   one number in \eqn{(0, 1]}. Default \code{0.5}.
#' @param tolerance    Whole number \eqn{\ge 0}: how many sub-threshold markers
#'   a block may absorb. Default \code{3}.
#' @param tol_reset    Logical: HapSelect's \code{tol_reset} (reset the
#'   tolerance count after a marker in LD). Default \code{TRUE}.
#' @param method       HapSelect's blocking method, \code{"flanking"} (default)
#'   or \code{"average"}.
#' @param start        Where a block starts, \code{"LD"} (default) or
#'   \code{"beginning"}; passed to \code{def_blocks()}.
#' @param attach       \code{TRUE} (default): return \code{geno_meta} with the
#'   blocks attached, its \code{marker_effects}, \code{phenotype} and
#'   \code{sample_info} carried across (\code{\link{attach_block_info}} alone
#'   would drop them). \code{FALSE}: return the block tables only.
#' @param output_dir   Optional directory for the two canonical block files,
#'   \code{<prefix>_block_info.csv} and \code{<prefix>_block_marker_map.csv}.
#' @param prefix       File-name prefix for \code{output_dir}. Default
#'   \code{"hapS_blocks"}.
#' @param parallelize  Passed to \code{HapSelect::pairwise_ld()}. Default
#'   \code{FALSE}.
#' @param verbose      \code{TRUE} (default) prints a summary (and HapSelect's
#'   progress bar where \pkg{progressr} shows one); \code{FALSE} is silent.
#'
#' @return With \code{attach = TRUE}, the \code{GenoMeta} with
#'   \code{block_marker_map_info} / \code{block_info} set. With
#'   \code{attach = FALSE}, a list with \code{block_marker_map}
#'   (\code{Chrom, Marker, Block}) and \code{block_info}
#'   (\code{Chrom, Block, Start, End, Pos, Num_SNPs, Block_length}), covering
#'   the blocked markers only. Block IDs are HapSelect's own,
#'   \code{"<chromosome key>:<n>"} (\code{"1:1"}, \code{"1:2"}, ...).
#'
#' @seealso \code{\link{attach_block_info}} (the recommended route: a block map
#'   you already have), \code{\link{st_ld_blocks}} (the SelectionTools
#'   resolver, experimental), \code{\link{build_block_location_info}},
#'   \code{\link{run_haplotype_catalogue}}
#'
#' @examples
#' \dontrun{
#' # The recommended route needs no LD blocking at all:
#' meta <- read_plink_ped("geno.ped", "geno.map", block_file = "blocks.tsv")
#'
#' # No block map? Build LD blocks locally with HapSelect:
#' meta <- read_plink_ped("geno.ped", "geno.map")
#' meta <- hapS_ld_blocks(meta)                        # r2 >= 0.5, tolerance 3
#' blk  <- hapS_ld_blocks(meta, ld_threshold = 0.6, tolerance = 1,
#'                        attach = FALSE, output_dir = "blocks")
#' }
#'
#' @export
hapS_ld_blocks <- function(geno_meta,
                           ld_threshold = 0.5,
                           tolerance    = 3L,
                           tol_reset    = TRUE,
                           method       = c("flanking", "average"),
                           start        = c("LD", "beginning"),
                           attach       = TRUE,
                           output_dir   = NULL,
                           prefix       = "hapS_blocks",
                           parallelize  = FALSE,
                           verbose      = TRUE) {
  caller <- "hapS_ld_blocks"
  if (!is_geno_meta(geno_meta, verbose = FALSE))
    stop(caller, ": 'geno_meta' must be a valid GenoMeta.", call. = FALSE)
  .ld_check_rule(ld_threshold, tolerance, caller)
  if (!is.logical(tol_reset) || length(tol_reset) != 1L || is.na(tol_reset))
    stop(caller, ": 'tol_reset' must be TRUE or FALSE.", call. = FALSE)
  method <- match.arg(method)
  start  <- match.arg(start)
  .hapS_require(caller)
  d <- dim(geno_meta$geno_arrays)
  if (d[1L] < 2L)
    stop(caller, ": LD needs at least 2 samples.", call. = FALSE)

  placed <- .ld_placed_markers(geno_meta$map_info)
  n_unplaced <- d[2L] - length(placed)
  if (length(placed) < 2L)
    stop(caller, ": fewer than 2 markers carry a usable position.",
         call. = FALSE)

  # Contain HapSelect's global side effects (see the file header).
  old_handlers <- getOption("progressr.handlers")
  on.exit(options(progressr.handlers = old_handlers), add = TRUE)
  if (!isTRUE(verbose)) {
    old_enable <- options(progressr.enable = FALSE)
    on.exit(options(old_enable), add = TRUE)
  }
  if (isTRUE(parallelize)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
  }

  t0 <- Sys.time()
  chrom_key <- .hapS_chrom_key(geno_meta$map_info$Chrom[
    match(placed, geno_meta$map_info$Marker)])
  fr   <- .hapS_dosage_frame(geno_meta, placed, chrom_key)
  geno <- fr$geno
  multi <- geno$SNP[fr$n_alleles > 2L]
  if (length(multi))
    warning(caller, ": ", length(multi), " marker(s) carry more than two ",
            "alleles; their LD is computed on the dosage of one allele ",
            "against all others (e.g. ", paste(utils::head(multi, 3L),
            collapse = ", "), ").", call. = FALSE)
  map <- geno[, c("SNP", "Chromosome", "Position")]
  rownames(map) <- NULL

  dose <- as.matrix(geno[, -(1:3), drop = FALSE])
  flat <- geno$SNP[apply(dose, 1L, function(x)
    length(unique(x[!is.na(x)])) <= 1L)]
  rm(dose)

  n_zero_sd <- 0L
  ld <- withCallingHandlers(
    HapSelect::pairwise_ld(geno, parallelize = isTRUE(parallelize)),
    warning = function(w) {
      if (grepl("standard deviation is zero", conditionMessage(w),
                fixed = TRUE)) {
        n_zero_sd <<- n_zero_sd + 1L
        invokeRestart("muffleWarning")
      }
    })
  blocks <- HapSelect::def_blocks(ld, map, method = method,
                                  tolerance = as.integer(tolerance),
                                  tol_reset = tol_reset,
                                  threshold = ld_threshold, start = start)
  hb <- HapSelect::block_obj_to_df(blocks, map)

  members <- strsplit(as.character(hb$Block), ";", fixed = TRUE)
  n_each  <- lengths(members)
  bm <- data.frame(
    Chrom  = chrom_key$Chrom[match(rep(hb$Chrom, n_each),
                                   chrom_key$Chromosome)],
    Marker = unlist(members, use.names = FALSE),
    Block  = rep(as.character(hb$Block_ID), n_each),
    stringsAsFactors = FALSE)
  if (anyDuplicated(bm$Marker) || !setequal(bm$Marker, geno$SNP))
    stop(caller, ": HapSelect's blocks do not cover each placed marker ",
         "exactly once (", length(unique(bm$Marker)), " of ", nrow(geno),
         " covered).", call. = FALSE)
  bm <- bm[order(match(bm$Marker, geno$SNP)), , drop = FALSE]
  rownames(bm) <- NULL

  if (isTRUE(verbose))
    message(sprintf(paste0(
      "%s: HapSelect blocks (r^2 >= %s, tolerance %d%s, %s): %s; %d unplaced ",
      "marker(s) left unblocked%s; %.1f s"),
      caller, format(ld_threshold), as.integer(tolerance),
      if (tol_reset) ", tol_reset" else "", method, .ld_sizes_note(bm),
      n_unplaced,
      if (length(flat)) sprintf(paste0(
        "; %d marker(s) without genotype variation (no LD; %d 'standard ",
        "deviation is zero' warnings muffled)"), length(flat), n_zero_sd)
      else "",
      as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  .ld_finish(geno_meta, bm,
             data.frame(Marker = geno$SNP,
                        Chrom  = chrom_key$Chrom[match(geno$Chromosome,
                                                       chrom_key$Chromosome)],
                        Pos    = geno$Position, stringsAsFactors = FALSE),
             attach, output_dir, prefix, caller)
}
