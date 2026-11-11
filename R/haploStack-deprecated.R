# =============================================================================
# haploStack-deprecated.R — the deprecation ladder (0.9.0, framework F7)
#
# Every name or argument the package retires climbs the same ladder:
#
#   deprecated  .hs_deprecated(): a warning on every call, class
#               c("hs_warning_deprecated", "deprecatedWarning", "hs_warning");
#               the old name keeps working and forwards to the replacement
#   defunct     .hs_defunct(): an error, class c("hs_error_defunct",
#               "defunctError", "hs_error_input", "hs_error"); a defunct
#               *argument* stays in the signature for one release so the
#               caller is told the replacement instead of "unused argument"
#   removed     at the next 0.x
#
# at least one minor release apart. No other file calls .Deprecated() or
# .Defunct() (hygiene test). The two Rd pages below (`haploStack-deprecated`,
# `haploStack-defunct`) carry the table of every name on the ladder; each
# rung also gets a CHANGELOG "Deprecated" / "Defunct" line.
# =============================================================================

.hs_ladder_text <- function(old, new, since, what, fn) {
  if (what == "argument") {
    if (is.null(fn)) fn <- "the function"
    c(sprintf("%s(%s =)", fn, old),
      if (is.null(new)) NA_character_ else sprintf("%s(%s =)", fn, new))
  } else {
    c(sprintf("'%s'", old), sprintf("'%s'", new))
  }
}

# Warn that `old` is deprecated since `since` in favour of `new`. `what`
# says whether `old` is a function or an argument (of `fn`); an argument
# retired without a replacement passes `new = NULL` ("... and is ignored").
# Warns on every call: there is no session memory to leak between jobs.
.hs_deprecated <- function(old, new, since, what = c("function", "argument"),
                           fn = NULL) {
  what <- match.arg(what)
  lab  <- .hs_ladder_text(old, new, since, what, fn)
  msg  <- if (what == "argument" && is.na(lab[2L]))
    sprintf("%s is deprecated since haploStack %s and is ignored.", lab[1L], since)
  else if (what == "argument")
    sprintf("%s is deprecated since haploStack %s; use %s.", lab[1L], since, lab[2L])
  else
    sprintf(paste0("%s is deprecated since haploStack %s.\nUse %s instead.\n",
                   "See help(\"haploStack-deprecated\")."), lab[1L], since, lab[2L])
  .hs_warn(msg, call. = FALSE,
           class = c("hs_warning_deprecated", "deprecatedWarning"),
           data = list(old = old, new = new, since = since))
}

# Stop because `old` is defunct since `since` (deprecated since `dep_since`).
.hs_defunct <- function(old, new, since, dep_since = NULL,
                        what = c("function", "argument"), fn = NULL) {
  what <- match.arg(what)
  lab  <- .hs_ladder_text(old, new, since, what, fn)
  dep  <- if (is.null(dep_since)) "" else sprintf(" (deprecated %s)", dep_since)
  msg  <- sprintf(paste0("%s is defunct since haploStack %s%s. Use %s.\n",
                         "See help(\"haploStack-defunct\")."),
                  lab[1L], since, dep, lab[2L])
  .hs_abort(msg, call. = FALSE,
            class = c("hs_error_defunct", "defunctError", "hs_error_input"),
            data = list(old = old, new = new, since = since))
}


#' Deprecated names and arguments
#'
#' Names on the first rung of haploStack's deprecation ladder: each still
#' works, forwards to its replacement and warns on every call (warning class
#' \code{"hs_warning_deprecated"}, also \code{"deprecatedWarning"}). One minor
#' release later it becomes defunct (\link{haploStack-defunct}) and errors with
#' the replacement in the message; it is removed at the next 0.x.
#'
#' \tabular{llll}{
#'   \strong{Deprecated} \tab \strong{Kind} \tab \strong{Since} \tab \strong{Use instead} \cr
#'   \code{greedy_selection()} \tab function \tab 0.9.0 \tab \code{\link{ohs_selection}(method = "Greedy")}, or \code{\link{run_stacking_simulation}(scheme = "OHS", ohs_method = "Greedy")} \cr
#'   \code{genetic_algorithm_selection()} \tab function \tab 0.9.0 \tab \code{\link{ohs_selection}()} (\code{method = "GA_pairwise"}), or \code{run_stacking_simulation(scheme = "OHS")} \cr
#'   \code{st_ld_blocks_files()} \tab function \tab 0.9.0 \tab \code{\link{st_convert}()} then \code{\link{st_ld_blocks}()} \cr
#'   \code{process_genotype()} \tab function \tab 0.9.0 \tab \code{\link{build_geno_meta}()} (or a reader) \cr
#'   \code{transform_genotype()} \tab function \tab 0.9.0 \tab \code{\link{write_geno_meta}()} / \code{\link{geno_meta_to_matrix}()} \cr
#'   \code{output_format =} on \code{\link{geno_meta_to_matrix}}, \code{\link{write_geno_meta}} \tab argument \tab 0.9.0 \tab \code{preset = "PL"} / \code{"ST"} / \code{"GS"} \cr
#'   \code{catalogue_results =}, \code{trait_name =}, \code{select_method =} and a catalogue list passed as \code{geno_meta} on \code{\link{run_stacking_simulation}} \tab argument \tab 0.9.0 \tab \code{geno_meta = <GenoMeta>}, \code{trait =}, \code{scheme =} \cr
#'   \code{extend =} on \code{\link{transform_effects}}, \code{\link{read_effect_file}} \tab argument \tab 0.9.0 \tab nothing: the rows it promised were never added \cr
#'   \code{legend =} on \code{\link{plot_meta_circos}} \tab argument \tab 0.9.0 \tab \code{legend_position = "none"} \cr
#'   \code{write_outputs =} on \code{\link{analyze_block_catalogue}} \tab argument \tab 0.9.0 \tab \code{output_dir =} alone decides \cr
#'   \code{table_file =} on \code{\link{plot_marker_signal}}, \code{\link{plot_marker_effects}} \tab argument \tab 0.9.0 \tab \code{\link{marker_signal_table}()}, \code{\link{save_marker_signal}()} \cr
#'   \code{\link{plot_block_effect_circos}()} \tab function \tab 0.9.0 \tab \code{plot_meta_circos(lapply(traits, track_block_effect, analysis = x))} \cr
#'   \code{encoding = "selectiontools"} in \code{\link{build_geno_meta}}, \code{\link{new_geno_meta}}, \code{\link{sim_geno_meta}}, \code{\link{sim_geno_matrix}} \tab argument \tab 0.9.0 \tab \code{\link{st_convert}()} for SelectionTools' files; the alleles are translated to nucleotides \cr
#'   \code{"selectiontools"} as a side of \code{\link{translate_alleles}()}, \code{\link{transform_effects}()} \tab argument \tab 0.9.0 \tab \code{\link{st_convert}()} \cr
#'   \code{use_missing_chars(mode = "selectiontools")} \tab argument \tab 0.9.0 \tab \code{st_convert(missing_chars =)} \cr
#' }
#'
#' \code{greedy_selection()} and \code{genetic_algorithm_selection()} were a
#' second implementation of the coverage selection the driver runs as the
#' \code{OHS} scheme (\code{\link{ohs_selection}}); they go in 0.10.
#' \code{process_genotype()} / \code{transform_genotype()} are the pre-GenoMeta
#' layer: a \code{GenoMeta} is built by \code{\link{build_geno_meta}} or a
#' reader and written by \code{\link{write_geno_meta}}. \code{st_ld_blocks_files()}
#' is \code{st_convert()} followed by \code{st_ld_blocks()}; all three
#' SelectionTools items are also behind \code{\link{use_st_formats}}.
#'
#' The \code{"selectiontools"} rows are one change: since 0.9.0 a
#' \code{GenoMeta} stores nucleotide alleles only, and SelectionTools' code
#' table (\code{A = 1, C = 2, G = 3, T = 4}, missing \code{-1}) is applied at
#' the boundary by \code{\link{st_convert}}. All three rows are behind
#' \code{\link{use_st_formats}}: while the switch is off they stop with an
#' \code{hs_error_input}; with it on they warn and translate.
#'
#' @return What the replacement returns.
#'
#' @name haploStack-deprecated
#' @seealso \link{haploStack-defunct}
#' @keywords internal
NULL

#' Defunct names and arguments
#'
#' Names on the second rung of haploStack's deprecation ladder: each errors
#' (class \code{"hs_error_defunct"}, also \code{"defunctError"}) with the
#' replacement in the message. A defunct argument stays in the signature for
#' one release so that a call naming it is told what to use instead of
#' failing with "unused argument". Removal follows at the next 0.x.
#'
#' \tabular{lllll}{
#'   \strong{Defunct} \tab \strong{Kind} \tab \strong{Deprecated} \tab \strong{Defunct} \tab \strong{Use instead} \cr
#'   \code{define_haplotype_blocks()} \tab function \tab 0.7.10 \tab 0.9.0 \tab \code{\link{st_convert}()} then \code{\link{st_ld_blocks}()} (or \code{\link{hapS_ld_blocks}()}, the default resolver) \cr
#'   \code{define_haplotype_blocks_meta()} \tab function \tab 0.7.10 \tab 0.9.0 \tab \code{\link{st_ld_blocks}()} \cr
#'   \code{center_by_freq =} on the eight centring entry points \tab argument \tab 0.6.0 \tab 0.9.0 \tab \code{center = "allele_freq"} (\code{TRUE}) / \code{center = "none"} (\code{FALSE}) \cr
#'   \code{recode_geno_meta(to_type = "selectiontools")} \tab argument \tab -- \tab 0.9.0 \tab \code{st_convert(to = "selectiontools")} \cr
#' }
#'
#' The two \code{define_haplotype_blocks*} names ran SelectionTools, now the
#' experimental resolver; the default one is \code{\link{hapS_ld_blocks}}
#' (HapSelect). Neither is needed when you already have a block map: give it to
#' \code{\link{attach_block_info}}, or to the reader. \code{center_by_freq}
#' stays in the eight signatures for this release so that a call naming it is
#' told the replacement; the entry points are \code{\link{run_haplotype_catalogue}},
#' \code{\link{compute_profile_localgebv}}, \code{\link{compute_block_catalogue}},
#' \code{\link{analyze_block_catalogue}}, \code{\link{count_block_subtypes}},
#' \code{\link{count_favourable_subtypes}}, \code{\link{count_favourable_blocks}}
#' and \code{\link{predict_gebv}}.
#'
#' \code{recode_geno_meta(to_type = "selectiontools")} skipped the deprecated
#' rung: since 0.9.0 a \code{GenoMeta} stores nucleotide alleles only, so
#' there was no container it could have returned. The other names climb the
#' ladder from \link{haploStack-deprecated}.
#'
#' @param geno_data_file,geno_map_file,ld_threshold,tolerance,output_dir,geno_meta,attach,prefix,keep_files,verbose
#'   Kept so that an old call is told the replacement instead of failing with
#'   "unused argument"; ignored.
#'
#' @name haploStack-defunct
#' @seealso \link{haploStack-deprecated}
#' @keywords internal
NULL

# Defunct since 0.9.0 (decision-form row 26 a: they warned since 0.7.10). The
# stubs keep their signatures and export until 0.10 so an old call is told the
# replacement.
#' @rdname haploStack-defunct
#' @export
define_haplotype_blocks <- function(geno_data_file,
                                    geno_map_file,
                                    ld_threshold = 0.5,
                                    tolerance    = 0,
                                    output_dir   = "haplotype") {
  .hs_defunct("define_haplotype_blocks()", "st_convert() + st_ld_blocks()",
              "0.9.0", dep_since = "0.7.10")
}

#' @rdname haploStack-defunct
#' @export
define_haplotype_blocks_meta <- function(geno_meta,
                                         ld_threshold = 0.5,
                                         tolerance    = 0L,
                                         attach       = TRUE,
                                         output_dir   = NULL,
                                         prefix       = "st_blocks",
                                         keep_files   = FALSE,
                                         verbose      = TRUE) {
  .hs_defunct("define_haplotype_blocks_meta()", "st_ld_blocks()", "0.9.0",
              dep_since = "0.7.10")
}
