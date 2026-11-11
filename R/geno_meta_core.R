# =============================================================================
# geno_meta_core.R — the container's schema, its one rebuild path, and its
# one validator core (0.9.0, design D2 §3.2–3.4, framework F6)
#
#   .HS_SCHEMA_VERSION    what this build writes into metadata$schema_version
#   .HS_META_SCOPE        which metadata keys are derived, vocabulary, sample-
#                         scoped or block-scoped (the carry-forward rules)
#   .rebuild_geno_meta()  the ONE way to derive a container from another:
#                         every transform goes through it, so what is kept
#                         and what is dropped is decided in one place
#   .hs_meta_reason()     the validator core; is_geno_meta() and
#                         .hs_check_meta() are its two faces
#   .hs_check_meta()      the entry-point guard (replaces ~60 copies of
#                         `if (!is_geno_meta(x)) .hs_abort(...)`)
#   upgrade_geno_meta()   schema 1 (0.6.0–0.8.8) -> 2, metadata only
#   .hs_check_rds_file()  the pre-read check load_geno_meta() applies
#
# The constructor itself, new_geno_meta(), moved here in D2-9 (it stayed in format_helper.R until the
# file split (D2-9).
# =============================================================================

.HS_SCHEMA_VERSION <- 2L

# Metadata keys by scope. `derived` keys are always recomputed by the
# constructor; `vocabulary` keys are kept unless a transform patches them;
# the `sample` key is kept only while every sample of the new container was a
# sample of the old one (offspring drop it); `block` keys are cleared when the
# block definitions are replaced and refreshed when they are merely pruned.
.HS_META_SCOPE <- list(
  derived    = c("n_samples", "n_markers", "ploidy", "schema_version", "built_with"),
  vocabulary = c("encoding", "char_dict", "missing_tokens", "source_format"),
  sample     = "family_history",
  block      = c("signal_selection", "blocks", "annotation"))

# The schema a container carries: absent means 1 (0.6.0–0.8.8; an unstamped
# missing-token vocabulary on top of that means pre-0.6.0).
.hs_schema_of <- function(x) {
  v <- x$metadata$schema_version
  if (is.null(v) || length(v) == 0L || is.na(v[1L])) 1L else as.integer(v[1L])
}

.hs_version_cache <- new.env(parent = emptyenv())
.hs_pkg_version <- function() {
  v <- .hs_version_cache$v
  if (is.null(v)) {
    v <- tryCatch(as.character(utils::packageVersion("haploStack")),
                  error = function(e) NA_character_)
    .hs_version_cache$v <- v
  }
  v
}

# Stamp the derived metadata keys from the array. The caller's values for
# them are not trusted (D2 §3.1): the array is the fact.
.hs_stamp_metadata <- function(metadata, geno_arrays) {
  d <- dim(geno_arrays)
  metadata$n_samples      <- as.integer(d[1L])
  metadata$n_markers      <- as.integer(d[2L])
  metadata$ploidy         <- as.integer(d[3L])
  metadata$schema_version <- .HS_SCHEMA_VERSION
  metadata$built_with     <- .hs_pkg_version()
  metadata
}


# ---------------------------------------------------------------------------
# One rebuild path
# ---------------------------------------------------------------------------

# Derive a container from `old`. Named arguments in `...` replace elements:
#   geno_arrays, map_info, block_info, block_marker_map_info, sample_info,
#   marker_effects, phenotype, annotation, source_format, extra_na,
#   missing_chars, and `metadata`, which is a PATCH merged into old$metadata
#   (a NULL value deletes the key).
# `.drop` names elements (or metadata keys) to leave out — a declared loss,
# e.g. c("family_history", "phenotype") for offspring. Everything else is
# forwarded, then aligned to the new array:
#   * map_info, when not given, is the old map subset to the new markers in
#     array order (C-D2); block_marker_map_info is pruned to them;
#   * sample_info, when not given, is the old FAM when every new sample is an
#     old one, else NULL (the constructor builds the dummy FAM);
#   * marker_effects and phenotype are pruned to the surviving markers /
#     samples; annotation is pruned with .hs_prune_annotation();
#   * metadata follows .HS_META_SCOPE: replacing the blocks (a new
#     block_marker_map_info) clears signal_selection, the gene annotation and
#     the block provenance; pruning refreshes their counts; new samples drop
#     family_history.
# `.caller` names the transform in messages.
.rebuild_geno_meta <- function(old, ..., .drop = character(0L),
                               .caller = ".rebuild_geno_meta") {
  repl <- list(...)
  allowed <- c("geno_arrays", "map_info", "block_info", "block_marker_map_info",
               "sample_info", "marker_effects", "phenotype", "annotation",
               "source_format", "extra_na", "missing_chars", "metadata")
  bad <- setdiff(names(repl), allowed)
  if (length(bad) || (length(repl) && (is.null(names(repl)) || any(!nzchar(names(repl))))))
    .hs_abort(class = "hs_error_internal", ".rebuild_geno_meta: unknown element(s) ",
              paste(bad, collapse = ", "), ".", call. = FALSE)
  elements <- setdiff(allowed, c("source_format", "extra_na", "missing_chars", "metadata"))
  bad_drop <- setdiff(.drop, c(elements, unlist(.HS_META_SCOPE), names(old$metadata)))
  if (length(bad_drop))
    .hs_abort(class = "hs_error_internal", ".rebuild_geno_meta: cannot drop ",
              paste(bad_drop, collapse = ", "), ".", call. = FALSE)

  # A v1 container is brought to the current schema first (metadata only).
  if (.hs_schema_of(old) < .HS_SCHEMA_VERSION)
    old <- .hs_upgrade_metadata(old, verbose = FALSE)

  given <- function(f) f %in% names(repl)
  take  <- function(f) {
    if (f %in% .drop) return(NULL)
    if (given(f)) repl[[f]] else old[[f]]
  }

  arr   <- take("geno_arrays")
  if (!is.array(arr) || length(dim(arr)) != 3L)
    .hs_abort(.caller, ": 'geno_arrays' must be a 3-D array.", call. = FALSE)
  new_s <- dimnames(arr)[[1L]]
  new_m <- dimnames(arr)[[2L]]
  old_s <- dimnames(old$geno_arrays)[[1L]]
  old_m <- dimnames(old$geno_arrays)[[2L]]
  markers_changed <- !identical(new_m, old_m)
  samples_new     <- !all(new_s %in% old_s)
  blocks_replaced <- given("block_marker_map_info") &&
    !"block_marker_map_info" %in% .drop

  # ---- map ------------------------------------------------------------------
  map_info <- take("map_info")
  if (!given("map_info") && markers_changed) {
    idx <- match(new_m, as.character(map_info$Marker))
    if (anyNA(idx))
      .hs_abort(class = "hs_error_internal", .caller, ": ", sum(is.na(idx)),
                " marker(s) of the new array are not in the old map.", call. = FALSE)
    map_info <- map_info[idx, , drop = FALSE]
    rownames(map_info) <- NULL
  }

  # ---- blocks ---------------------------------------------------------------
  bmm <- take("block_marker_map_info")
  bi  <- take("block_info")
  if (!blocks_replaced && !is.null(bmm) && markers_changed) {
    bmm <- as.data.frame(bmm, stringsAsFactors = FALSE)
    bmm <- bmm[as.character(bmm$Marker) %in% new_m, , drop = FALSE]
    rownames(bmm) <- NULL
  }
  if (!given("block_info") && !is.null(bi) && markers_changed)
    bi <- .prune_signal_markers(bi, new_m, caller = .caller)

  # ---- samples --------------------------------------------------------------
  sample_info <- if (given("sample_info") || "sample_info" %in% .drop) take("sample_info")
                 else if (!samples_new && !is.null(old$sample_info))
                   old$sample_info[match(new_s, as.character(old$sample_info$Sample)), ,
                                   drop = FALSE]
                 else NULL
  if (!is.null(sample_info)) rownames(sample_info) <- NULL

  marker_effects <- take("marker_effects")
  if (!given("marker_effects") && !is.null(marker_effects) && markers_changed) {
    marker_effects <- if (is.data.frame(marker_effects))
      marker_effects[marker_effects$Marker %in% new_m, , drop = FALSE]
    else lapply(marker_effects, function(df) df[df$Marker %in% new_m, , drop = FALSE])
  }

  phenotype <- take("phenotype")
  if (!given("phenotype") && !is.null(phenotype) && !identical(new_s, old_s))
    phenotype <- phenotype[as.character(phenotype$Sample) %in% new_s, , drop = FALSE]

  # ---- metadata (scope rules) ----------------------------------------------
  md <- old$metadata
  patch <- if (given("metadata")) repl$metadata else list()
  if (!is.list(patch))
    .hs_abort(class = "hs_error_internal", .caller,
              ": the metadata patch must be a list.", call. = FALSE)
  for (k in intersect(.drop, c(unlist(.HS_META_SCOPE), names(md)))) md[[k]] <- NULL
  if (samples_new && !"family_history" %in% names(patch))
    md$family_history <- NULL

  # ---- annotation -----------------------------------------------------------
  annotation <- take("annotation")
  custom_rederive <- FALSE
  if (!given("annotation") && !"annotation" %in% .drop) {
    if (blocks_replaced) {
      had_any <- !is.null(annotation) || !is.null(md$annotation)
      custom  <- annotation$custom
      if (had_any) {
        .hs_inform(.caller, ": block definitions replaced — the stored ",
                   "block annotation was dropped. Re-run ",
                   "compute_block_annotation()/attach_block_annotation() against ",
                   "the new blocks.",
                   if (!is.null(custom))
                     " The custom marker annotations are kept; their Block is re-derived."
                   else "")
      }
      annotation <- if (!is.null(custom)) list(custom = custom) else NULL
      md$annotation <- if (!is.null(custom) && !is.null(md$annotation$custom))
        list(custom = md$annotation$custom) else NULL
      custom_rederive <- !is.null(custom)
    } else if (!is.null(annotation) && markers_changed) {
      keep_blocks <- if (!is.null(bmm)) unique(as.character(bmm$Block))
                     else if (!is.null(bi)) unique(as.character(bi$Block))
                     else character(0L)
      annotation <- .hs_prune_annotation(annotation, keep_blocks, new_m)
      md$annotation <- .hs_annotation_counts(annotation, md$annotation)   # ann/F7
    }
  }
  if (blocks_replaced) {
    if (!"signal_selection" %in% names(patch)) md$signal_selection <- NULL
    if (!"blocks" %in% names(patch)) md$blocks <- NULL
  }
  for (k in names(patch)) {
    if (is.null(patch[[k]])) md[[k]] <- NULL else md[[k]] <- patch[[k]]
  }

  out <- new_geno_meta(
    geno_arrays           = arr,
    map_info              = map_info,
    block_info            = bi,
    block_marker_map_info = bmm,
    metadata              = md,
    sample_info           = sample_info,
    marker_effects        = marker_effects,
    phenotype             = phenotype,
    annotation            = annotation,
    source_format         = repl$source_format,
    extra_na              = repl$extra_na,
    # A derived container keeps the vocabulary its source was built with: the
    # session registry is not consulted again unless the transform says so.
    missing_chars         = if (given("missing_chars")) repl$missing_chars
                            else character(0L))

  # Custom marker annotations survive a re-blocking with their Block re-derived
  # against the blocks the container now holds (ann/F23).
  if (custom_rederive && !is.null(out$annotation$custom))
    out$annotation$custom <- .hs_rederive_custom_block(out$annotation$custom,
                                                      out$block_marker_map_info)

  # Pruned blocks: refresh the counts the block-scoped keys carry.
  if (!blocks_replaced && markers_changed) {
    if (!is.null(out$metadata$signal_selection) && !is.null(out$block_info$N_Signal))
      out$metadata$signal_selection$n_blocks <-
        as.integer(sum(out$block_info$N_Signal > 0L, na.rm = TRUE))
    if (!is.null(out$metadata$blocks) && !is.null(out$block_marker_map_info)) {
      n_new <- length(unique(as.character(out$block_marker_map_info$Block)))
      n_old <- out$metadata$blocks$n_blocks
      if (!is.null(n_old) && !identical(as.integer(n_old), n_new)) {
        out$metadata$blocks$pruned_from <- out$metadata$blocks$pruned_from %||% as.integer(n_old)
        out$metadata$blocks$n_blocks    <- n_new
      }
    }
  }
  out
}

# Re-derive the Block column of a custom marker-annotation table against a
# block marker map: a row whose marker is in the map takes that block, a
# position-only row (Marker NA) or an unmapped marker gets NA.
.hs_rederive_custom_block <- function(custom, bmm) {
  if (is.null(custom) || !"Block" %in% names(custom)) return(custom)
  if (is.null(bmm)) { custom$Block <- NA_character_; return(custom) }
  idx <- match(as.character(custom$Marker), as.character(bmm$Marker))
  custom$Block <- as.character(bmm$Block)[idx]
  custom
}


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

# The validator core: NULL when `x` is a valid GenoMeta, else one sentence
# naming the first failing condition. Every check is O(markers + samples);
# the array itself is never scanned here (the missing-token scan is the
# caller's decision, see is_geno_meta()).
#   strict = TRUE adds: sample_info rows in array order (not just the same
#   set) and no block_info row without a marker (orphan rows).
.hs_meta_reason <- function(x, strict = FALSE) {
  few <- function(v) paste0(paste(utils::head(v, 5L), collapse = ", "),
                            if (length(v) > 5L) " ..." else "")
  if (!is.list(x))
    return("'x' is not a list.")

  required <- c("geno_arrays", "map_info", "block_info",
                "block_marker_map_info", "sample_info", "metadata")
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0L)
    return(paste0("missing required field(s): ",
                  paste(missing_fields, collapse = ", "), "."))

  arr <- x$geno_arrays
  if (!is.array(arr) || length(dim(arr)) != 3L)
    return("'geno_arrays' must be a 3-D array.")
  d <- dim(arr)

  if (!is.data.frame(x$map_info))
    return("'map_info' must be a data.frame.")
  if (d[2L] != nrow(x$map_info))
    return(paste0("dim(geno_arrays)[2] (", d[2L],
                  ") does not equal nrow(map_info) (", nrow(x$map_info), ")."))
  markers <- dimnames(arr)[[2L]]
  if (!is.null(markers) && !identical(markers, as.character(x$map_info$Marker)))
    return(paste0("map_info$Marker is not in the array's marker order ",
                  "(dimnames(geno_arrays)[[2]]); rebuild the container with ",
                  "new_geno_meta() or filter_geno_meta()."))

  md <- x$metadata
  if (!is.list(md))
    return("'metadata' must be a list.")
  v <- .hs_schema_of(x)
  if (v > .HS_SCHEMA_VERSION)
    return(paste0("metadata$schema_version ", v, " is newer than this ",
                  "haploStack writes (", .HS_SCHEMA_VERSION, "); it was saved by ",
                  "haploStack ", md$built_with %||% "(unknown)",
                  ". Upgrade the package to read it."))
  for (k in c("n_samples", "n_markers", "ploidy")) {
    i <- match(k, c("n_samples", "n_markers", "ploidy"))
    val <- md[[k]]
    if (is.null(val))
      return(paste0("metadata$", k, " is missing."))
    if (length(val) != 1L || is.na(val) || as.integer(val) != d[i])
      return(paste0("metadata$", k, " (", paste(val, collapse = ","),
                    ") does not equal dim(geno_arrays)[", i, "] (", d[i], ")."))
  }

  if (!is.null(x$block_marker_map_info)) {
    bmm <- x$block_marker_map_info
    if (!is.data.frame(bmm))
      return("'block_marker_map_info' must be a data.frame or NULL.")
    bad <- setdiff(as.character(bmm$Marker), as.character(x$map_info$Marker))
    if (length(bad) > 0L)
      return(paste0(length(bad),
                    " marker(s) in block_marker_map_info absent from map_info: ",
                    few(bad), "."))
    dup <- as.character(bmm$Marker)[duplicated(as.character(bmm$Marker))]
    if (length(dup) > 0L)
      return(paste0(length(unique(dup)), " marker(s) assigned to more than ",
                    "one block in block_marker_map_info: ", few(unique(dup)),
                    ". Blocks must not overlap."))
    if (isTRUE(strict) && !is.null(x$block_info) && is.data.frame(x$block_info)) {
      orphan <- setdiff(as.character(x$block_info$Block), as.character(bmm$Block))
      if (length(orphan) > 0L)
        return(paste0(length(orphan), " block_info row(s) without a marker in ",
                      "block_marker_map_info: ", few(orphan), "."))
    }
  }

  # sample_info: required, must align with the genotype sample set
  samples <- dimnames(arr)[[1L]]
  if (!is.data.frame(x$sample_info))
    return("'sample_info' must be a data.frame.")
  if (!"Sample" %in% names(x$sample_info))
    return("'sample_info' must contain a 'Sample' column.")
  si <- as.character(x$sample_info$Sample)
  if (!setequal(si, samples))
    return(paste0("'sample_info$Sample' does not match the genotype ",
                  "sample set (dimnames(geno_arrays)[[1]])."))
  if (isTRUE(strict) && !identical(si, samples))
    return("'sample_info' rows are not in the array's sample order.")

  # marker_effects: always optional (never required for a valid GenoMeta). The
  # canonical form is a named list of per-trait effect tables; a bare data.frame
  # is also accepted. When present, every table must carry Marker/Allele/Effect.
  if (!is.null(x$marker_effects)) {
    me_list <- if (is.data.frame(x$marker_effects)) list(x$marker_effects)
               else x$marker_effects
    if (!is.list(me_list))
      return(paste0("'marker_effects' must be a data.frame, a named list ",
                    "of data.frames, or NULL."))
    need <- c("Marker", "Allele", "Effect")
    for (me in me_list) {
      if (!is.data.frame(me))
        return("each 'marker_effects' table must be a data.frame.")
      miss <- setdiff(need, names(me))
      if (length(miss) > 0L)
        return(paste0("'marker_effects' table missing column(s): ",
                      paste(miss, collapse = ", "), "."))
    }
  }

  # phenotype: optional (absent in legacy objects); when present must be a
  # data.frame with a Sample column, a subset of the genotype sample set.
  if ("phenotype" %in% names(x) && !is.null(x$phenotype)) {
    if (!is.data.frame(x$phenotype))
      return("'phenotype' must be a data.frame or NULL.")
    if (!"Sample" %in% names(x$phenotype))
      return("'phenotype' must contain a 'Sample' column.")
    bad_p <- setdiff(as.character(x$phenotype$Sample), samples)
    if (length(bad_p) > 0L)
      return(paste0(length(bad_p),
                    " phenotype sample(s) absent from the genotype set: ",
                    few(bad_p), "."))
  }

  # annotation: optional (#033); every annotated Block must exist in the
  # container's block vocabulary and every Marker in map_info.
  if ("annotation" %in% names(x) && !is.null(x$annotation)) {
    known_blocks <- if (!is.null(x$block_marker_map_info)) x$block_marker_map_info$Block
                    else if (!is.null(x$block_info)) x$block_info$Block
                    else character(0L)
    reason <- .hs_check_annotation(x$annotation, known_blocks, x$map_info$Marker)
    if (!is.null(reason)) return(reason)
  }
  NULL
}

# The entry-point guard. Aborts with hs_error_input when `x` is not a valid
# GenoMeta, naming the caller, the argument and the reason (plus `hint`, a
# sentence on how to get a valid one, also carried as the condition's `hint`
# field); the condition's call is the caller's own call. On an unstamped
# (pre-0.6.0 or hand-built) container the missing-token notice is raised as
# is_geno_meta() raises it.
.hs_check_meta <- function(x, caller, arg = "geno_meta", strict = FALSE,
                           hint = NULL) {
  reason <- .hs_meta_reason(x, strict = strict)
  if (!is.null(reason)) {
    msg <- paste0(caller, ": '", arg, "' is not a valid GenoMeta: ", reason,
                  if (!is.null(hint)) paste0(" ", hint) else "")
    stop(.hs_condition(msg, sys.call(-1L), c("hs_error_input", "hs_error"),
                       if (is.null(hint)) list() else list(hint = hint),
                       errorCondition))
  }
  if (isTRUE(strict) || is.null(x$metadata$missing_tokens))
    .notify_missing_tokens(.missing_token_report(x), caller)
  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# Schema upgrade
# ---------------------------------------------------------------------------

# Metadata-only v1 -> v2 patch, without validation (the upgrade's first half;
# .rebuild_geno_meta() applies it silently to a v1 input).
.hs_upgrade_metadata <- function(x, verbose = TRUE) {
  md  <- x$metadata
  arr <- x$geno_arrays
  if (is.array(arr) && length(dim(arr)) == 3L) {
    d <- dim(arr)
    changed <- character(0L)
    for (i in 1:3) {
      k <- c("n_samples", "n_markers", "ploidy")[i]
      if (is.null(md[[k]]) || length(md[[k]]) != 1L || is.na(md[[k]]) ||
          as.integer(md[[k]]) != d[i]) changed <- c(changed, k)
      md[[k]] <- as.integer(d[i])
    }
    if (length(changed) && isTRUE(verbose))
      .hs_inform("upgrade_geno_meta: metadata$", paste(changed, collapse = ", $"),
                 " re-derived from the array.")
    md$char_dict <- tryCatch(.hs_char_dict(md$char_dict, arr),
                             error = function(e) .hs_char_dict(NULL, arr))
  }
  bmm <- x$block_marker_map_info
  if (!is.null(bmm) && is.data.frame(bmm) && is.null(md$blocks))
    md$blocks <- list(source = "unknown",
                      n_blocks = length(unique(as.character(bmm$Block))))
  if (!is.null(md$signal_selection) && is.data.frame(x$block_info) &&
      "N_Signal" %in% names(x$block_info))
    md$signal_selection$n_blocks <-
      as.integer(sum(x$block_info$N_Signal > 0L, na.rm = TRUE))
  md$schema_version <- .HS_SCHEMA_VERSION
  md$built_with     <- .hs_pkg_version()
  x$metadata <- md
  x
}


#' Bring a GenoMeta saved by an earlier haploStack to the current schema
#'
#' Containers carry \code{metadata$schema_version} since 0.9.0 (the current
#' schema is 2; a container without the key is schema 1, written by 0.6.0 to
#' 0.8.8). \code{\link{load_geno_meta}} upgrades on read; call this directly
#' on a container you hold in memory, for instance one restored with
#' \code{readRDS()} yourself.
#'
#' The upgrade touches \strong{metadata only}: the counters
#' \code{n_samples} / \code{n_markers} / \code{ploidy} are re-derived from the
#' array, \code{char_dict} becomes a named list, a container with blocks gets
#' a block provenance record (\code{metadata$blocks}, source \code{"unknown"}),
#' a signal-marker selection has its block count refreshed, and the schema
#' and package version are stamped. Missing tokens are never repaired —
#' \code{\link{repair_geno_meta}} is a separate, explicit step. The one change
#' to the genotype array is the SelectionTools recode: a container saved with
#' \code{metadata$encoding = "selectiontools"} (1/2/3/4 codes) has its alleles
#' translated to nucleotides and is stamped \code{"nucleotide"}, since 0.9.0
#' stores no other categorical vocabulary (see \code{\link{use_st_formats}});
#' that step runs whatever the schema says.
#'
#' A container that cannot pass the strict check after the upgrade (two
#' blocks sharing a marker, an annotation pointing at blocks the container
#' does not have) is refused with a classed error naming the fix.
#'
#' @param x       A \code{GenoMeta} list (any schema).
#' @param verbose Logical. Report what changed (default \code{TRUE}).
#' @param strict_missing Logical or \code{NULL}. Whether a container that still
#'   carries a missing-allele token is an error (\code{TRUE}) or a warning
#'   (\code{FALSE}); \code{NULL} (default) takes
#'   \code{options(haploStack.strict_missing)}.
#'
#' @return The container at schema \code{2}. A container already at the
#'   current schema is returned unchanged, unless it carries the
#'   SelectionTools vocabulary (then it is recoded).
#'
#' @examples
#' meta <- sim_geno_meta(6, 12, block_size = 4, seed = 1)
#' meta$metadata$schema_version   # 2
#' old <- meta; old$metadata$schema_version <- NULL   # as 0.8.8 wrote it
#' identical(upgrade_geno_meta(old, verbose = FALSE)$metadata$schema_version, 2L)
#'
#' @seealso \code{\link{load_geno_meta}}, \code{\link{is_geno_meta}},
#'   \code{\link{repair_geno_meta}}
#' @export
upgrade_geno_meta <- function(x, verbose = TRUE, strict_missing = NULL) {
  if (!is.list(x) || !is.list(x$metadata))
    .hs_abort(class = c("hs_error_schema", "hs_error_input"),
              "upgrade_geno_meta: 'x' is not a GenoMeta container (no metadata ",
              "list).")
  v <- .hs_schema_of(x)
  if (v > .HS_SCHEMA_VERSION)
    .hs_abort(class = c("hs_error_schema", "hs_error_input"),
              "upgrade_geno_meta: the container carries schema ", v,
              ", newer than this haploStack writes (", .HS_SCHEMA_VERSION,
              "); it was saved by haploStack ",
              x$metadata$built_with %||% "(unknown)", ". Upgrade the package.")
  # The SelectionTools recode (0.9.0, C-5) runs whatever the schema says: a
  # container built by a 0.9.0 development build is already schema 2 and may
  # still carry the old vocabulary.
  x <- .st_upgrade_container(x, verbose = verbose)
  if (v == .HS_SCHEMA_VERSION) return(x)

  x <- .hs_upgrade_metadata(x, verbose = verbose)

  bmm <- x$block_marker_map_info
  if (!is.null(bmm) && is.data.frame(bmm)) {
    dup <- unique(as.character(bmm$Marker)[duplicated(as.character(bmm$Marker))])
    if (length(dup) > 0L)
      .hs_abort(class = c("hs_error_block_map", "hs_error_input"),
                data = list(bad = dup),
                "upgrade_geno_meta: ", length(dup), " marker(s) belong to more ",
                "than one block (", paste(utils::head(dup, 5L), collapse = ", "),
                if (length(dup) > 5L) " ..." else "",
                "). Such a container cannot run a catalogue; re-attach the ",
                "blocks with attach_block_info(overlap = \"first\"), or fix ",
                "the map first.")
  }
  if (!is.null(x$annotation)) {
    known_blocks <- if (!is.null(bmm)) bmm$Block
                    else if (!is.null(x$block_info)) x$block_info$Block
                    else character(0L)
    reason <- .hs_check_annotation(x$annotation, known_blocks, x$map_info$Marker)
    if (!is.null(reason))
      .hs_abort(class = c("hs_error_schema", "hs_error_input"),
                "upgrade_geno_meta: the stored annotation does not fit the ",
                "container (", reason, ") Recompute it with ",
                "compute_block_annotation().")
  }
  reason <- .hs_meta_reason(x, strict = TRUE)
  if (!is.null(reason))
    .hs_abort(class = c("hs_error_schema", "hs_error_input"),
              "upgrade_geno_meta: the container cannot be upgraded: ", reason)
  .notify_missing_tokens(.missing_token_report(x), "upgrade_geno_meta",
                         strict_missing = strict_missing)
  if (isTRUE(verbose))
    .hs_inform("upgrade_geno_meta: schema ", v, " -> ", .HS_SCHEMA_VERSION,
               " (metadata only; the genotype array is unchanged).")
  x
}


# ---------------------------------------------------------------------------
# The pre-read check for an RDS file (D1-9, form row 32)
# ---------------------------------------------------------------------------

# readRDS() of untrusted bytes is not a safe boundary (R core's own advice;
# it is also a decompression-bomb vector). Before unserialising, check that
# the file is under `max_bytes` and starts like an R serialisation: a gzip,
# bzip2 or xz magic, or the uncompressed "X\n" / "A\n" / "B\n" header.
# Anything else is refused with hs_error_format. This is a plausibility
# check, not a security boundary: load_geno_meta() reads trusted files only.
.hs_check_rds_file <- function(path, max_bytes = Inf, caller = "load_geno_meta") {
  sz <- file.size(path)
  if (is.finite(max_bytes) && !is.na(sz) && sz > max_bytes)
    .hs_abort(class = c("hs_error_format", "hs_error_input"),
              caller, ": '", path, "' is ", sz, " bytes, over max_bytes = ",
              format(max_bytes, scientific = FALSE), ".")
  magic <- tryCatch(readBin(path, what = "raw", n = 6L), error = function(e) raw(0))
  starts <- function(bytes) length(magic) >= length(bytes) &&
    identical(magic[seq_along(bytes)], as.raw(bytes))
  ok <- starts(c(0x1f, 0x8b)) ||                               # gzip
        starts(charToRaw("BZh")) ||                             # bzip2
        starts(c(0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00)) ||        # xz
        starts(charToRaw("X\n")) || starts(charToRaw("A\n")) || # uncompressed
        starts(charToRaw("B\n"))
  if (!ok)
    .hs_abort(class = c("hs_error_format", "hs_error_input"),
              caller, ": '", path, "' does not start like an R serialisation ",
              "(RDS): no gzip/bzip2/xz header and no uncompressed X/A/B ",
              "header. Only a file written by save_geno_meta() / saveRDS() is ",
              "read here, and only from a trusted source.")
  invisible(TRUE)
}


# ---- moved here in 0.9.0 (D2-9) from format_helper.R ----

#' Create and validate a GenoMeta internal metadata object
#'
#' A \strong{GenoMeta} is a plain named list that consolidates all genotype
#' and map data into a single structure consumed by the haploStack pipeline.
#' No formal class system is used — any field can be read or modified directly.
#' Strict dimension validation is applied only here and in transform functions.
#'
#' @section Structure:
#' \describe{
#'   \item{\code{geno_arrays}}{3-D integer/character array
#'     \code{[n_samples × n_markers × ploidy]}. Values are encoded per
#'     \code{metadata$char_dict}.}
#'   \item{\code{map_info}}{data.frame \code{(Chrom, Marker, Pos)}.
#'     \code{nrow(map_info) == dim(geno_arrays)[2]}.}
#'   \item{\code{block_info}}{data.frame \code{(Chr, Block, Start, End)} or
#'     \code{NULL}. A block keeps its row only while at least one of its
#'     markers is in the container (a fully filtered block is pruned).}
#'   \item{\code{block_marker_map_info}}{data.frame \code{(Chr, Marker, Block)}
#'     or \code{NULL}. \code{Marker} column must be a subset of
#'     \code{map_info$Marker}.}
#'   \item{\code{sample_info}}{Required FAM-structured data.frame with positional
#'     columns \code{(Family, Sample, Parent1, Parent2, Sex)}. Column \emph{index}
#'     is the contract — an input fam is read by position and stored under these
#'     canonical names (so \code{$Sample} is stable). \code{Sample} (position 2)
#'     must match \code{dimnames(geno_arrays)[[1]]}. When \code{NULL}, a dummy
#'     table is built from the genotype sample names (other columns \code{"0"}).
#'     See \code{docs/geno_meta_structure.md}.}
#'   \item{\code{marker_effects}}{Optional allele effects, or \code{NULL}.
#'     Stored canonically as a \strong{named list of per-trait data.frames},
#'     each \code{(Marker, Allele, Effect)} with the list name giving the trait
#'     — the same shape as the \code{effect_list} consumed by
#'     \code{\link{run_haplotype_catalogue}} and produced by
#'     \code{\link{fit_genomic_prediction}}. A bare data.frame (optionally with a
#'     \code{Trait} column to split on) is accepted and normalised to this list.
#'     Used as the default effect source by
#'     \code{\link{compute_profile_localgebv}} and
#'     \code{\link{run_haplotype_catalogue}}; an explicit effect argument always
#'     overrides it. Never required for a valid GenoMeta.}
#'   \item{\code{phenotype}}{Optional data.frame of sample-level trait values, or
#'     \code{NULL}. First column \code{Sample} (matching
#'     \code{dimnames(geno_arrays)[[1]]}); remaining columns are numeric trait
#'     measurements. Keyed on \code{Sample} the way \code{marker_effects} is keyed
#'     on \code{Marker}: standardised on construction, pruned to the surviving
#'     samples by \code{\link{filter_geno_meta}}. Attach via
#'     \code{\link{attach_phenotype}} / \code{\link{read_phenotype}}.}
#'   \item{\code{metadata}}{Named list:
#'     \code{ploidy}, \code{n_samples}, \code{n_markers} (integers, always
#'     derived from the array), \code{schema_version} (\code{2L} since 0.9.0;
#'     absent means 1, see \code{\link{upgrade_geno_meta}}), \code{built_with}
#'     (the package version that built it), \code{char_dict} (named list of
#'     integer codes: the observed allele vocabulary, informational, derived
#'     when absent), \code{source_format} (character), and \code{encoding}
#'     (auto-detected when absent). Optionally \code{missing_tokens}
#'     (character): the missing-allele vocabulary this container was
#'     normalised against, written wherever normalisation ran (see
#'     \code{\link{repair_geno_meta}}, \code{\link{use_missing_chars}}); absent
#'     on containers built before 0.6.0, and on any container that has not
#'     been normalised. Transforms may add \code{family_history},
#'     \code{signal_selection}, \code{blocks} and \code{annotation}.}
#' }
#'
#' @param geno_arrays         3-D array \code{[n_samples × n_markers × ploidy]}.
#' @param map_info            data.frame with \code{Chr}, \code{Marker},
#'   \code{Pos}.
#' @param block_info          data.frame with \code{Chr}, \code{Block},
#'   \code{Start}, \code{End}, or \code{NULL}.
#' @param block_marker_map_info data.frame with \code{Chr}, \code{Marker},
#'   \code{Block}, or \code{NULL}.
#' @param metadata            Named list with at minimum \code{source_format}.
#'   \code{ploidy}, \code{n_samples} and \code{n_markers} are derived from the
#'   array (a caller value is ignored); \code{schema_version} and
#'   \code{built_with} are stamped. An optional \code{char_dict} (allele ->
#'   code) is normalised to a named list of integer codes; when absent it is
#'   derived from the observed alleles.
#' @param sample_info         FAM-structured data.frame read by position
#'   (\code{Family}, \code{Sample}, \code{Parent1}, \code{Parent2}, \code{Sex}),
#'   or \code{NULL} (default) to auto-build a dummy table from the genotype
#'   sample names. A one-column input is treated as the \code{Sample} list.
#' @param marker_effects      Optional allele effects, or \code{NULL} (default).
#'   Either a named list of per-trait data.frames (\code{Marker}, \code{Allele},
#'   \code{Effect}) or a single data.frame (optionally with a \code{Trait}
#'   column); normalised on construction to the named-list form.
#' @param phenotype          Optional data.frame of sample-level trait values
#'   (\code{Sample} + numeric trait columns), or \code{NULL} (default). Aligned to
#'   the genotype sample set on construction (see \code{\link{attach_phenotype}}).
#' @param annotation         Optional annotation list (the gene tables
#'   \code{blocks}/\code{features}/\code{marker_gene} and/or \code{custom}),
#'   validated against the block and marker vocabularies exactly as
#'   \code{\link{is_geno_meta}} validates it; see
#'   \code{\link{attach_block_annotation}}.
#' @param source_format      Allele vocabulary of the source being parsed
#'   (\code{"nucleotide"}, \code{"selectiontools"}, \code{"numeric"}), used to
#'   normalise the missing alleles. Defaults to \code{metadata$encoding}.
#'   Declare it whenever the format is known rather than detected: a 1/2-coded
#'   PED, for instance, looks like dosage once split, and under the dosage rule
#'   its \code{"0"} would then be kept as a real genotype.
#' @param extra_na           Optional character vector of additional missing
#'   sentinels this particular source uses. For a rule that should apply to
#'   every read in a session, use \code{\link{use_missing_chars}} instead.
#' @param missing_chars      Extra missing-allele characters for this build:
#'   \code{NULL} (default) takes the session registry
#'   (\code{\link{use_missing_chars}}); a character vector names them for the
#'   source format directly (\code{character(0)} ignores the registry); a
#'   per-mode list is read like the registry. Resolved here and nowhere later:
#'   the vocabulary used is stamped in \code{metadata$missing_tokens}.
#'
#' @section The normalisation portal:
#' Since 0.6.0 this function is where a missing allele becomes \code{NA}, and it
#' is the \strong{only} such place. That is a property of the call graph rather
#' than a convention: \code{\link{read_plink_ped}},
#' \code{\link{read_dart_genotype}} and \code{\link{read_vcf_genotype}} bypass
#' \code{\link{build_geno_meta}} and construct here directly, as do the internal
#' transforms, so every route from a source file to a container passes through
#' this one function. Normalisation runs after the encoding is resolved and
#' before anything is returned, so no consumer can observe a container that has
#' not been through it.
#'
#' The vocabulary used is recorded in \code{metadata$missing_tokens}. A
#' container that already carries a stamp covering the vocabulary in force is
#' left alone rather than re-scanned, which is what keeps a derived container (a
#' filter, a block swap, a crossing generation) cheap.
#'
#' \strong{Normalisation is one-way.} Once a container exists, the source
#' token that a given \code{NA} came from is gone — \code{NA} is the only
#' representation left, and no function can tell which character the file wrote
#' there. \code{\link{recode_geno_meta}} can render missing back out as a
#' character of your choosing, but that is a new rendering, not a recovery of
#' the original file.
#'
#' @export
new_geno_meta <- function(geno_arrays,
                           map_info,
                           block_info              = NULL,
                           block_marker_map_info   = NULL,
                           metadata,
                           sample_info             = NULL,
                           marker_effects          = NULL,
                           phenotype               = NULL,
                           annotation              = NULL,
                           source_format           = NULL,
                           extra_na                = NULL,
                           missing_chars           = NULL) {

  # ---- Validate required fields ----------------------------------------
  if (!is.array(geno_arrays) || length(dim(geno_arrays)) != 3L)
    .hs_abort("new_geno_meta: 'geno_arrays' must be a 3-D array.")
  if (!is.data.frame(map_info))
    .hs_abort("new_geno_meta: 'map_info' must be a data.frame.")
  required_map <- c("Chrom", "Marker", "Pos")
  missing_map  <- setdiff(required_map, names(map_info))
  if (length(missing_map) > 0L)
    .hs_abort("new_geno_meta: 'map_info' missing columns: ",
         paste(missing_map, collapse = ", "))
  if (dim(geno_arrays)[2L] != nrow(map_info))
    .hs_abort("new_geno_meta: dim(geno_arrays)[2] (", dim(geno_arrays)[2L],
         ") must equal nrow(map_info) (", nrow(map_info), ").")

  if (!is.null(block_marker_map_info)) {
    if (!is.data.frame(block_marker_map_info))
      .hs_abort("new_geno_meta: 'block_marker_map_info' must be a data.frame or NULL.")

    # Filter: keep only markers present in map_info
    bad_markers <- setdiff(block_marker_map_info$Marker, map_info$Marker)
    if (length(bad_markers) > 0L) {
      .hs_warn("new_geno_meta: ", length(bad_markers),
              " marker(s) in block_marker_map_info absent from map_info \u2014 dropped: ",
              paste(head(bad_markers, 5L), collapse = ", "),
              if (length(bad_markers) > 5L) " ..." else "")
      block_marker_map_info <- block_marker_map_info[
        block_marker_map_info$Marker %in% map_info$Marker, , drop = FALSE
      ]
    }

    # Filter block_info to only blocks still represented after marker pruning
    if (!is.null(block_info)) {
      retained_blocks <- unique(block_marker_map_info$Block)
      block_info <- block_info[
        block_info$Block %in% retained_blocks, , drop = FALSE
      ]
    }

    # Coerce identifier columns now so rbind() below stays type-consistent
    block_marker_map_info$Chrom  <- as.character(block_marker_map_info$Chrom)
    block_marker_map_info$Marker <- as.character(block_marker_map_info$Marker)
    block_marker_map_info$Block  <- as.character(block_marker_map_info$Block)

    # Dummy blocks: every map_info marker absent from block_marker_map_info
    # becomes its own singleton block (one marker per dummy block).
    absent_markers <- setdiff(map_info$Marker, block_marker_map_info$Marker)
    if (length(absent_markers) > 0L) {
      .hs_inform("new_geno_meta: ", length(absent_markers),
              " map marker(s) absent from block_marker_map_info \u2014 ",
              "created singleton dummy block(s).")

      map_idx     <- match(absent_markers, map_info$Marker)
      dummy_chrom <- as.character(map_info$Chrom[map_idx])
      dummy_pos   <- as.numeric(map_info$Pos[map_idx])
      dummy_block <- paste0("dummy_", absent_markers)

      block_marker_map_info <- rbind(
        block_marker_map_info[, c("Chrom", "Marker", "Block"), drop = FALSE],
        data.frame(Chrom  = dummy_chrom,
                   Marker = as.character(absent_markers),
                   Block  = dummy_block,
                   stringsAsFactors = FALSE)
      )

      # Update block_info with one row per dummy block (Start = End = Pos)
      if (!is.null(block_info)) {
        block_info$Chrom <- as.character(block_info$Chrom)
        block_info$Block <- as.character(block_info$Block)
        dummy_bi <- data.frame(Chrom = dummy_chrom,
                               Block = dummy_block,
                               Start = dummy_pos,
                               End   = dummy_pos,
                               stringsAsFactors = FALSE)
        # Fill the extra block_info columns so rbind aligns. A singleton has a
        # real position and length, so the descriptive columns get their true
        # values rather than NA (0.9.0); anything else is NA.
        extra_cols <- setdiff(names(block_info), names(dummy_bi))
        for (cc in extra_cols) {
          dummy_bi[[cc]] <- switch(cc,
            Pos          = dummy_pos,
            Num_SNPs     = rep(1L, length(dummy_pos)),
            Block_length = rep(1, length(dummy_pos)),
            NA)
        }
        dummy_bi <- dummy_bi[, names(block_info), drop = FALSE]
        block_info <- rbind(block_info, dummy_bi)
      }
    }

    # Reorder block_marker_map_info to follow map_info marker order
    block_marker_map_info <- block_marker_map_info[
      order(match(block_marker_map_info$Marker, map_info$Marker)), , drop = FALSE
    ]
    rownames(block_marker_map_info) <- NULL

    # Reorder block_info by genomic position (Chrom, then Start)
    if (!is.null(block_info)) {
      block_info <- block_info[
        .hs_order(block_info$Chrom, block_info$Start), , drop = FALSE
      ]
      rownames(block_info) <- NULL
    }
  }

  if (!is.list(metadata) || !"source_format" %in% names(metadata))
    .hs_abort("new_geno_meta: 'metadata' must be a named list with at least ",
         "'source_format'.")
  # ---- SelectionTools allele codes (0.9.0, C-5) -------------------------------
  # A GenoMeta stores nucleotide alleles only. A declared "selectiontools"
  # encoding is deprecated: the array is translated here, at the portal, and
  # the container is stamped "nucleotide"; integer codes with no declared
  # encoding are refused (they may be SelectionTools codes or PLINK's 1/2
  # labels, and the caller must say which). See R/st_formats.R.
  st <- .st_portal(geno_arrays, metadata, "new_geno_meta")
  geno_arrays <- st$geno_arrays
  metadata    <- st$metadata
  # Derived, not trusted (0.9.0): n_samples, n_markers and ploidy come from
  # the array; the schema version and the building package are stamped.
  metadata <- .hs_stamp_metadata(metadata, geno_arrays)
  # char_dict is the observed allele vocabulary, informational: always stored
  # as a named list of integer codes (0.9.0), derived from the array when the
  # caller did not give one.
  metadata$char_dict <- .hs_char_dict(metadata$char_dict, geno_arrays)

  # ---- sample_info: required FAM-structured field, aligned to genotypes ------
  geno_samples <- dimnames(geno_arrays)[[1L]]
  sample_info  <- .standardise_sample_info(sample_info, geno_samples,
                                           caller = "new_geno_meta")

  # ---- marker_effects: optional; normalise to the canonical named list of
  #      per-trait effect tables (Marker, Allele, Effect). A bare data.frame
  #      (optionally carrying a Trait column) is accepted for convenience. ------
  if (!is.null(marker_effects)) {
    marker_effects <- .as_effect_list(marker_effects)
  }

  # ---- phenotype: optional, align to the genotype sample set -----------------
  if (!is.null(phenotype)) {
    phenotype <- .standardise_phenotype(phenotype, geno_samples,
                                        caller = "new_geno_meta")
  }

  # ---- annotation: optional; the same check is_geno_meta() applies ----------
  if (!is.null(annotation)) {
    known_blocks <- unique(as.character(
      if (!is.null(block_marker_map_info)) block_marker_map_info$Block
      else if (!is.null(block_info)) block_info$Block
      else character(0L)))
    reason <- .hs_check_annotation(annotation, known_blocks,
                                   as.character(map_info$Marker))
    if (!is.null(reason)) .hs_abort("new_geno_meta: ", reason)
  }

  # ---- encoding: store in metadata when absent (best-effort detection) -------
  # The internal detector: SelectionTools codes were already refused above, so
  # the verdict here is "nucleotide", "numeric" or nothing.
  if (is.null(metadata$encoding)) {
    metadata$encoding <- tryCatch(
      .detect_encoding(.meta_to_allele_matrix(list(geno_arrays = geno_arrays))),
      error = function(e) NULL
    )
    if (identical(metadata$encoding, "selectiontools")) metadata$encoding <- NULL
  }

  # ---- Normalise missing alleles — THE portal (FR3) --------------------------
  # Every route from a source file to a GenoMeta ends here: build_geno_meta(),
  # read_plink_ped(), read_dart_genotype() and read_vcf_genotype() all construct
  # through this function, and so does every internal transform. Putting the
  # normalisation anywhere else would leave a way in that skips it.
  #
  # It runs AFTER the encoding is resolved and BEFORE anything is returned, so
  # the vocabulary is known when the mapping is applied and no consumer can
  # observe a container that has not been through it.
  norm <- .normalise_at_portal(geno_arrays, metadata, source_format, extra_na,
                               missing_chars = missing_chars,
                               caller = "new_geno_meta")
  geno_arrays <- norm$array
  metadata    <- norm$metadata

  # ---- Coerce Chrom / Marker / Block columns to character --------------------
  map_info$Chrom  <- as.character(map_info$Chrom)
  map_info$Marker <- as.character(map_info$Marker)

  if (!is.null(block_marker_map_info)) {
    block_marker_map_info$Chrom  <- as.character(block_marker_map_info$Chrom)
    block_marker_map_info$Marker <- as.character(block_marker_map_info$Marker)
    block_marker_map_info$Block  <- as.character(block_marker_map_info$Block)
  }

  if (!is.null(block_info)) {
    block_info$Chrom <- as.character(block_info$Chrom)
    block_info$Block <- as.character(block_info$Block)
  }

  # Return plain named list — no class attribute
  out <- list(
    geno_arrays           = geno_arrays,
    map_info              = map_info,
    block_info            = block_info,
    block_marker_map_info = block_marker_map_info,
    sample_info           = sample_info,
    marker_effects        = marker_effects,
    phenotype             = phenotype,
    metadata              = metadata
  )
  if (!is.null(annotation)) out$annotation <- annotation
  out
}


# Internal: the container's char_dict as a named list of integer codes.
# Accepts the list the readers build, the named integer vector older callers
# (and 0.8.8 containers) carry, or NULL, which derives the vocabulary from the
# alleles observed in the array (missing tokens excluded, NA ignored).
.hs_char_dict <- function(x, geno_arrays) {
  if (is.null(x) || length(x) == 0L) {
    obs <- .hs_sort(unique(as.vector(geno_arrays)))
    obs <- obs[!is.na(obs) & nzchar(obs)]
    return(stats::setNames(as.list(seq_along(obs)), obs))
  }
  if (is.null(names(x)) || any(!nzchar(names(x))))
    .hs_abort("new_geno_meta: 'metadata$char_dict' must be named (allele -> code).")
  vals <- unlist(x, use.names = FALSE)
  if (!is.numeric(vals) || length(vals) != length(x))
    .hs_abort("new_geno_meta: 'metadata$char_dict' must map each allele to one ",
              "integer code.")
  stats::setNames(as.list(as.integer(vals)), names(x))
}


# Internal: the one annotation check, shared by new_geno_meta(),
# is_geno_meta() and attach_block_annotation(). `ann` is the container's
# annotation list; `known_blocks` / `markers` the vocabularies it must point
# into. Returns NULL when valid, else the reason as one string (the caller
# decides whether that is a FALSE or an error).
.hs_check_annotation <- function(ann, known_blocks, markers) {
  gene3 <- c("blocks", "features", "marker_gene")
  if (!is.list(ann) || length(ann) == 0L ||
      !all(names(ann) %in% c(gene3, "custom")))
    return(paste0("'annotation' must be a list of the tables 'blocks'/",
                  "'features'/'marker_gene' and/or 'custom' (see ",
                  "attach_block_annotation() / attach_marker_annotations())."))
  if (any(gene3 %in% names(ann)) && !all(gene3 %in% names(ann)))
    return(paste0("'annotation' carries only part of the gene-table ",
                  "trio (blocks/features/marker_gene)."))
  few <- function(v) paste0(paste(utils::head(v, 5L), collapse = ", "),
                            if (length(v) > 5L) " ..." else "")
  known_blocks <- unique(as.character(known_blocks))
  markers      <- as.character(markers)
  if ("blocks" %in% names(ann)) {
    bad_b <- setdiff(unique(as.character(ann$blocks$Block)), known_blocks)
    if (length(bad_b) > 0L)
      return(paste0(length(bad_b), " annotated block(s) absent from ",
                    "the container's blocks: ", few(bad_b), "."))
    bad_m <- setdiff(unique(as.character(ann$marker_gene$Marker)), markers)
    if (length(bad_m) > 0L)
      return(paste0(length(bad_m), " annotated marker(s) absent from ",
                    "map_info: ", few(bad_m), "."))
  }
  if ("custom" %in% names(ann)) {
    cm <- ann$custom$Marker
    bad_c <- setdiff(unique(as.character(cm[!is.na(cm)])), markers)
    if (length(bad_c) > 0L)
      return(paste0(length(bad_c), " custom-annotation marker(s) ",
                    "absent from map_info: ", few(bad_c), "."))
  }
  NULL
}


# Internal: normalise any marker-effects object to the canonical named list of
# per-trait effect tables. Each element is a data.frame (Marker, Allele, Effect)
# and the element name is the trait. Accepts:
#   * NULL                                    -> NULL
#   * a named list of data.frames             -> standardised element-wise
#   * a single data.frame with a Trait column -> split into one table per trait
#   * a single data.frame without Trait       -> wrapped under `default_name`
.as_effect_list <- function(x, default_name = "Trait") {
  if (is.null(x)) return(NULL)

  if (is.data.frame(x)) {
    if ("Trait" %in% names(x)) {
      keep <- setdiff(names(x), "Trait")
      sp   <- split(x[, keep, drop = FALSE], as.character(x$Trait))
      return(lapply(sp, .standardise_effect_df))
    }
    return(stats::setNames(list(.standardise_effect_df(x)), default_name))
  }

  if (is.list(x)) {
    if (length(x) == 0L) return(list())
    if (is.null(names(x)) || any(!nzchar(names(x))))
      .hs_abort("marker effects given as a list must be fully named ",
           "(one element name per trait).")
    return(lapply(x, .standardise_effect_df))
  }

  .hs_abort("marker effects must be a data.frame or a named list of data.frames.")
}


# Internal: standardise a single effect table to (Marker, Allele, Effect),
# reading the first three columns by position.
.standardise_effect_df <- function(d) {
  if (!is.data.frame(d) && !inherits(d, "data.table"))
    .hs_abort("each marker-effect table must be a data.frame.")
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  if (ncol(d) < 3L)
    .hs_abort("each marker-effect table needs >= 3 columns (Marker, Allele, Effect).")
  data.frame(
    Marker = as.character(d[[1L]]),
    Allele = as.character(d[[2L]]),
    Effect = as.numeric(d[[3L]]),
    stringsAsFactors = FALSE
  )
}


# Internal: collapse a marker-effects object (named list of per-trait tables or
# a bare data.frame) to a single effect table, or NULL when empty. Used where a
# single-trait table is required (e.g. founder-pool breeding values). `trait`
# selects a named element; otherwise the first table is used.
.coerce_single_effect_df <- function(x, trait = NULL) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    if (length(x) == 0L) return(NULL)
    if (!is.null(trait) && trait %in% names(x)) return(x[[trait]])
    return(x[[1L]])
  }
  NULL
}


# Internal: standardise a phenotype table to (Sample, <trait...>) aligned to the
# genotype sample set.
#
#   ph           data.frame with a Sample column (or a first column treated as
#                Sample) plus one or more numeric trait columns
#   geno_samples character vector dimnames(geno_arrays)[[1]]
#
# The Sample column is matched against geno_samples: rows for unknown samples are
# dropped with a warning, genotype samples lacking a phenotype row are added with
# NA trait values (left-join semantics), and the result is reordered to the
# genotype sample order. Trait columns are coerced to numeric.
.standardise_phenotype <- function(ph, geno_samples, caller = "new_geno_meta") {
  if (!is.data.frame(ph))
    .hs_abort(caller, ": 'phenotype' must be a data.frame or NULL.")
  ph <- as.data.frame(ph, stringsAsFactors = FALSE)
  if (ncol(ph) < 2L)
    .hs_abort(caller, ": 'phenotype' must have a Sample column plus at least one ",
         "trait column.")

  # Resolve the Sample column: prefer an explicit "Sample" name, else position 1.
  sample_col <- if ("Sample" %in% names(ph)) "Sample" else names(ph)[1L]
  sample_ids <- as.character(ph[[sample_col]])
  trait_cols <- setdiff(names(ph), sample_col)

  std <- data.frame(Sample = sample_ids, stringsAsFactors = FALSE)
  for (tc in trait_cols)
    std[[tc]] <- suppressWarnings(as.numeric(ph[[tc]]))

  geno_samples <- as.character(geno_samples)

  # Drop phenotype rows for samples absent from the genotype set.
  extra <- setdiff(std$Sample, geno_samples)
  if (length(extra) > 0L) {
    .hs_warn(caller, ": ", length(extra),
            " phenotype sample(s) absent from the genotype set \u2014 dropped: ",
            paste(utils::head(extra, 5L), collapse = ", "),
            if (length(extra) > 5L) " ..." else "")
    std <- std[std$Sample %in% geno_samples, , drop = FALSE]
  }
  if (anyDuplicated(std$Sample))
    .hs_abort(caller, ": 'phenotype' has duplicate Sample IDs.")

  # Left-join onto the genotype sample order: genotype samples without a
  # phenotype row get NA traits.
  missing_pheno <- setdiff(geno_samples, std$Sample)
  if (length(missing_pheno) > 0L) {
    .hs_warn(caller, ": ", length(missing_pheno),
            " genotype sample(s) lack a phenotype row \u2014 filled with NA: ",
            paste(utils::head(missing_pheno, 5L), collapse = ", "),
            if (length(missing_pheno) > 5L) " ..." else "")
    pad <- std[seq_len(length(missing_pheno)), , drop = FALSE]
    pad$Sample <- missing_pheno
    for (tc in setdiff(names(std), "Sample")) pad[[tc]] <- NA_real_
    std <- rbind(std, pad)
  }

  std <- std[match(geno_samples, std$Sample), , drop = FALSE]
  rownames(std) <- NULL
  std
}


# Canonical FAM column names for sample_info (positional contract; see
# docs/geno_meta_structure.md).
.SAMPLE_INFO_COLS <- c("Family", "Sample", "Parent1", "Parent2", "Sex")


# Internal: coerce any fam-like input to the canonical 5-column FAM sample_info,
# aligned (by the Sample column = position 2) to the genotype sample set.
#
#   x            data.frame fam, or NULL (build a dummy table from geno_samples)
#   geno_samples character vector dimnames(geno_arrays)[[1]]
#
# Columns are read BY POSITION (1=Family, 2=Sample, 3=Parent1, 4=Parent2,
# 5=Sex); a one-column input is treated as Sample for back-compatibility.
# Missing/short columns are dummy-filled with "0". The result always carries the
# canonical names so sample_info$Sample stays stable.
.standardise_sample_info <- function(x, geno_samples, caller = "new_geno_meta") {
  n <- length(geno_samples)

  if (is.null(x)) {
    si <- data.frame(
      Family  = rep("0", n),
      Sample  = as.character(geno_samples),
      Parent1 = rep("0", n),
      Parent2 = rep("0", n),
      Sex     = rep("0", n),
      stringsAsFactors = FALSE
    )
    return(si)
  }

  if (!is.data.frame(x))
    .hs_abort(caller, ": 'sample_info' must be a data.frame or NULL.")
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  if (ncol(x) == 1L) {
    # Back-compat: a single column is the Sample/IID list.
    si <- data.frame(
      Family  = rep("0", nrow(x)),
      Sample  = as.character(x[[1L]]),
      Parent1 = rep("0", nrow(x)),
      Parent2 = rep("0", nrow(x)),
      Sex     = rep("0", nrow(x)),
      stringsAsFactors = FALSE
    )
  } else {
    # Read up to five columns by POSITION; pad any missing with "0".
    col <- function(i) if (ncol(x) >= i) as.character(x[[i]]) else rep("0", nrow(x))
    si <- data.frame(
      Family  = col(1L),
      Sample  = col(2L),
      Parent1 = col(3L),
      Parent2 = col(4L),
      Sex     = col(5L),
      stringsAsFactors = FALSE
    )
  }

  if (!setequal(si$Sample, geno_samples))
    .hs_abort(caller, ": 'sample_info' Sample column (position 2) must match the ",
         "genotype sample set exactly (dimnames(geno_arrays)[[1]]).")

  # Reorder rows to genotype sample order
  si <- si[match(geno_samples, si$Sample), , drop = FALSE]
  rownames(si) <- NULL
  si
}


#' Check whether an object is a valid GenoMeta
#'
#' Structural validation of a container, cheap by default: every check is
#' O(markers + samples) and the genotype array itself is not scanned. The
#' checks are: the required elements; a 3-D \code{geno_arrays}; a
#' \code{map_info} with one row per marker \strong{in the array's marker
#' order}; the metadata counters \code{n_samples} / \code{n_markers} /
#' \code{ploidy} equal to the array's dimensions; a schema version this
#' haploStack can read; no marker in two blocks; \code{sample_info} matching
#' the sample set; well-formed \code{marker_effects}, \code{phenotype} and
#' \code{annotation}.
#'
#' When the check fails and \code{verbose = TRUE}, a message describing the
#' first failing condition is emitted before \code{FALSE} is returned. Inside
#' the package the same core runs as the entry-point guard, whose error
#' carries the reason.
#'
#' @section The missing-allele invariant:
#' Inside a \code{GenoMeta} a missing allele is \code{NA} and nothing else.
#' A container built by 0.6.0 or later records the vocabulary it was
#' normalised against in \code{metadata$missing_tokens}, and that stamp is
#' trusted: the array is not re-scanned. A container \strong{without} the
#' stamp (built before 0.6.0, or by hand) is scanned, and if it still carries
#' a source token as an allele (PLINK's \code{"0"}, \code{"N"}, \code{"."},
#' \code{"-"}) it is \strong{structurally valid} — the function still returns
#' \code{TRUE} — but a warning names the token, the number of markers
#' affected, and \code{\link{repair_geno_meta}}. \code{strict = TRUE} scans
#' regardless of the stamp. \code{strict_missing = TRUE} (or the option
#' \code{haploStack.strict_missing}) turns that warning into an error.
#'
#' @param x       Any R object.
#' @param verbose Logical. When \code{TRUE} (default), emit a diagnostic
#'   message explaining why \code{x} is not a valid GenoMeta.
#' @param strict  Logical. \code{TRUE} adds the full checks: the missing-token
#'   scan regardless of the stamp, \code{sample_info} rows in the array's
#'   sample order (not only the same set), and no \code{block_info} row
#'   without a marker. Used by \code{\link{load_geno_meta}} and
#'   \code{\link{upgrade_geno_meta}}; \code{FALSE} (default) everywhere else.
#' @param strict_missing Logical or \code{NULL}. Whether a missing token still
#'   in the array is an error (\code{TRUE}) or a warning (\code{FALSE});
#'   \code{NULL} (default) takes \code{options(haploStack.strict_missing)},
#'   itself \code{FALSE} by default.
#' @return \code{TRUE} if \code{x} appears to be a valid GenoMeta;
#'   \code{FALSE} otherwise.
#'
#' @seealso \code{\link{repair_geno_meta}}, \code{\link{upgrade_geno_meta}}
#' @export
is_geno_meta <- function(x, verbose = TRUE, strict = FALSE, strict_missing = NULL) {
  reason <- .hs_meta_reason(x, strict = strict)
  if (!is.null(reason)) {
    if (verbose) .hs_inform("is_geno_meta: ", reason)
    return(FALSE)
  }
  # The container is structurally valid. The one content invariant (a missing
  # allele is NA and nothing else, FR4) is checked by scanning the array only
  # when the container carries no stamp — every container the package builds
  # does — or under strict = TRUE. A pre-0.6.0 container can still carry the
  # source token, so this warns and names the repair rather than refusing;
  # strict_missing turns it into an error.
  if (isTRUE(strict) || is.null(x$metadata$missing_tokens))
    .notify_missing_tokens(.missing_token_report(x), "is_geno_meta",
                           strict_missing = strict_missing)
  TRUE
}


# ---------------------------------------------------------------------------
# Format readers
# ---------------------------------------------------------------------------
