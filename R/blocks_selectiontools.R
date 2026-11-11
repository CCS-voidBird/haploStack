# =============================================================================
# blocks_selectiontools.R — LD haplotype blocks via SelectionTools (EXPERIMENTAL)
#
# The SelectionTools resolver, one of the two local LD-block backends
# (R/blocks_ld_common.R says when to use a resolver at all; HapSelect,
# R/blocks_hapselect.R, is the default):
#
#   st_ld_blocks()        a GenoMeta in (named define_haplotype_blocks_meta()
#                         in 0.7.9)
#   st_ld_blocks_files()  SelectionTools' own file formats in (named
#                         define_haplotype_blocks() up to 0.7.9)
#
# Both are behind use_st_formats() (0.9.0: SelectionTools is being retired),
# and the exchange of formats lives in R/st_formats.R (st_convert() and the
# helpers .st_write_inputs() / .st_parse_blocks() the resolver uses).
#
# Experimental: SelectionTools needed every guard below to run safely from a
# package (thread #034), and its partition differs from HapSelect's at the same
# settings. The container route exists because the file route could not be
# driven from a container (INDEX discrepancy 34):
#
#   * write_geno_meta(output_format = "st") writes haploStack's ST *encoding*
#     as CSV; SelectionTools' matrix reader separates on tab, blank or ";" only,
#     so it read that file as 0 individuals;
#   * SelectionTools keeps the leading integer of a chromosome name ("1H" -> 1),
#     so joining its output back to a named map loses every position;
#   * st.calc.ld / st.def.hblocks / st.recode.hil write their side files
#     ("ld", "blocks.hmap", "blocks.hall") into the working directory, and a
#     read into data set "default" replaces whatever the user kept there;
#   * ld.threshold = 0 ends the R session inside SelectionTools' C code.
#
# The bridge therefore owns the whole exchange: it writes SelectionTools' own
# matrix ("m") and map ("mcp") formats under ALIASES (m1.., i1.., integer
# chromosome keys) into a private working directory, runs SelectionTools there
# on a dedicated data set, parses st.def.hblocks()'s one-row-per-block table
# ("Markers" = "m1;m2;..."), translates aliases back to the container's marker
# IDs and chromosome names, and deletes the directory. Aliases also make marker
# and sample names immaterial to SelectionTools' tokenizer.
# =============================================================================


# Fail early, with the alternatives, when SelectionTools is absent.
.st_require <- function(caller) {
  if (!requireNamespace("SelectionTools", quietly = TRUE))
    .hs_abort(class = "hs_error_dependency",
         data = list(hint = "SelectionTools is distributed outside CRAN; see the haploStack README"),
         caller, " requires the SelectionTools package, which is not ",
         "installed. Alternatives: a block map you already have (the ",
         "recommended route: attach_block_info(), or block_file = / blocks = ",
         "when reading), or hapS_ld_blocks(), the default LD resolver.",
         call. = FALSE)
  invisible(TRUE)
}

# ld.threshold = 0 kills the R session inside st.def.hblocks() (C level, no
# error), so .ld_check_rule() refuses it with this reason.
.st_zero_note <- "0 is refused: SelectionTools aborts the R session on it."

# SelectionTools resolves every file name against the variables st.input.dir
# and st.output.dir. It holds neither itself: its .onLoad creates both, as "",
# in the GLOBAL environment (where its manual has users set them), and its
# functions find them there. Since 0.9.0 (D1 §3.8) the bridge points both at
# its private working directory for the call and restores whatever the user
# had — present or absent — on exit, so the working directory is never
# changed and a user's own setting is neither refused nor lost. After
# rm(list = ls()) the variables are gone and SelectionTools alone would stop
# with "object 'st.input.dir' not found"; the bridge supplies them and removes
# them again.
.st_with_io_dirs <- function(work_dir, code) {
  vars <- c("st.input.dir", "st.output.dir")
  ge   <- globalenv()
  had  <- vapply(vars, exists, logical(1L), envir = ge, inherits = FALSE)
  old  <- lapply(vars, function(v)
    if (exists(v, envir = ge, inherits = FALSE)) get(v, envir = ge, inherits = FALSE)
    else NULL)
  on.exit({
    for (i in seq_along(vars)) {
      if (had[i]) assign(vars[i], old[[i]], envir = ge)
      else suppressWarnings(rm(list = vars[i], envir = ge))
    }
  }, add = TRUE)
  dir <- paste0(normalizePath(work_dir, winslash = "/", mustWork = TRUE), "/")
  for (v in vars) assign(v, dir, envir = ge)
  force(code)
}

# Run SelectionTools' LD blocking on two files already in the current working
# directory. Returns st.def.hblocks()'s table (one row per block; Name =
# block ID, Markers = "a;b;c;").
.st_run_hblocks <- function(geno_file, map_file, ld_threshold, tolerance,
                            data_set = "hs_ld_blocks", verbose = TRUE) {
  run <- function() {
    SelectionTools::st.read.marker.data(geno_file, format = "m",
                                        data.set = data_set)
    SelectionTools::st.read.map(map_file, m.stretch = 1, format = "mcp",
                                skip = 1, data.set = data_set)
    # auxfiles = FALSE: the LD table stays in SelectionTools' memory, where
    # st.def.hblocks() reads it; writing it out only to read it back into R
    # costs minutes on a genome-wide panel and returns nothing we use.
    SelectionTools::st.calc.ld(ld.measure = "r2", auxfiles = FALSE,
                               data.set = data_set)
    SelectionTools::st.def.hblocks(ld.threshold = ld_threshold,
                                   ld.criterion = "flanking",
                                   tolerance    = as.integer(tolerance),
                                   data.set     = data_set)
  }
  if (isTRUE(verbose)) return(run())
  res <- NULL
  utils::capture.output(res <- run(), type = "output")
  res
}

#' LD haplotype blocks with SelectionTools (experimental)
#'
#' \strong{A side route, and experimental.} haploStack does not need LD
#' blocking: the recommended input is a block map you already have, given to
#' \code{\link{read_plink_ped}} (\code{block_file =}),
#' \code{\link{build_geno_meta}} (\code{blocks =}) or
#' \code{\link{attach_block_info}} — no LD parameter involved. When blocks must
#' be built locally, \code{\link{hapS_ld_blocks}} (HapSelect) is the default
#' resolver; this function is the SelectionTools alternative. It is marked
#' experimental because SelectionTools needed several guards to run safely from
#' a package (see \emph{Known SelectionTools behaviour}), and because
#' \strong{its partition is not HapSelect's}: at the same \code{ld_threshold} /
#' \code{tolerance} the two give very different blocks, so compare block sets
#' only within one backend.
#'
#' Give it a container, an r\eqn{^2} threshold and a tolerance; it writes the
#' container into SelectionTools' own input formats in a private temporary
#' directory, runs \code{st.calc.ld()} + \code{st.def.hblocks(ld.criterion =
#' "flanking")} there, reads the blocks back in the container's own marker IDs
#' and chromosome names, deletes the directory, and (by default) attaches the
#' blocks. Named \code{define_haplotype_blocks_meta()} in 0.7.9; that name still
#' works, with a deprecation warning. SelectionTools is being retired from
#' haploStack: this function is behind \code{\link{use_st_formats}} (off by
#' default), and the exchange of formats goes through \code{\link{st_convert}}.
#'
#' @section What SelectionTools sees:
#' Markers are written under aliases (\code{m1}, \code{m2}, ...), samples as
#' \code{i1}, \code{i2}, ..., chromosomes as integer keys in natural order, and
#' each marker's alleles recoded to \code{1, 2, ...} (missing = \code{-}). LD
#' (\eqn{r^2}) does not depend on allele labels, and aliases make marker and
#' sample names — spaces, \code{|}, very long IDs — irrelevant to
#' SelectionTools' tokenizer. Everything is translated back before return.
#'
#' Markers without a usable position (non-finite or \eqn{\le} 0) or on an
#' unplaced chromosome (\code{0}, \code{UN}, \code{unknown}) are not blocked:
#' SelectionTools orders markers by position. On attach they become singleton
#' blocks, as any marker absent from a block map does. Markers with more than
#' five allele states cannot be written in SelectionTools' allele codes and are
#' treated the same way, with a warning.
#'
#' SelectionTools' side files (\code{ld}, \code{blocks.hmap}, the \code{.srt}
#' map) are written into the private directory, never the working directory,
#' and its data are loaded into the dedicated data set \code{"hs_ld_blocks"},
#' so a SelectionTools session's \code{"default"} data set is left alone.
#'
#' @section SelectionTools' directory variables:
#' SelectionTools resolves file names against \code{st.input.dir} and
#' \code{st.output.dir}, two variables its \code{.onLoad} creates as \code{""}
#' in the global environment. For the call both are pointed at a private
#' temporary directory and restored afterwards (present or absent, whatever
#' you had); the working directory is not changed. Process state SelectionTools
#' cannot release (its thread count, the resident LD data set and those
#' globals) is documented here because the package cannot undo it: run
#' SelectionTools in its own R process when that matters.
#'
#' @section Known SelectionTools behaviour:
#' \itemize{
#'   \item \code{ld_threshold = 0} ends the R session inside SelectionTools' C
#'     code, so it is refused.
#'   \item At \code{ld_threshold = 1} the last marker pair of the map still
#'     merges into one block.
#'   \item Its tolerance is far more permissive than HapSelect's. On 300 barley
#'     lines (1H + 2H, 2 677 markers) at \eqn{r^2 \ge 0.5} it gives 1 418 / 836
#'     / 294 blocks at tolerance 0 / 1 / 3, where \code{hapS_ld_blocks()} gives
#'     1 774 / 1 425 at tolerance 1 / 3. Hence this function's default tolerance
#'     of 0, against 3 for \code{hapS_ld_blocks()}; the partitions still differ.
#'   \item Diploid genotypes only.
#' }
#'
#' @param geno_meta    A valid, diploid \code{GenoMeta}.
#' @param ld_threshold \eqn{r^2} above which flanking markers stay in one
#'   block; one number in \eqn{(0, 1]}. Default \code{0.5}. \code{0} is refused:
#'   SelectionTools aborts the R session on it.
#' @param tolerance    Whole number \eqn{\ge 0}: how many markers below the
#'   threshold a block may absorb (SelectionTools' \code{tolerance}). Default
#'   \code{0}.
#' @param attach       \code{TRUE} (default): return \code{geno_meta} with the
#'   blocks attached. Its \code{marker_effects}, \code{phenotype} and
#'   \code{sample_info} are carried across — \code{\link{attach_block_info}}
#'   alone would drop them. \code{FALSE}: return the block tables only.
#' @param output_dir   Optional directory. When given, the two canonical block
#'   files are written there: \code{<prefix>_block_info.csv} and
#'   \code{<prefix>_block_marker_map.csv}.
#' @param prefix       File-name prefix for \code{output_dir}. Default
#'   \code{"st_blocks"}.
#' @param keep_files   \code{FALSE} (default) deletes the private directory;
#'   \code{TRUE} keeps it (the SelectionTools input files, for inspection) and
#'   reports its path.
#' @param verbose      \code{TRUE} (default) shows SelectionTools' console
#'   output and a summary; \code{FALSE} silences both.
#' @param st_formats   \code{TRUE} enables SelectionTools formats for this call;
#'   \code{NULL} (default) takes the session switch
#'   (\code{\link{use_st_formats}}, off by default: SelectionTools is being
#'   retired from haploStack).
#' @param progress     \code{NULL} (default) or a function of one argument (the
#'   package's progress contract, 0.9.0). SelectionTools' C code cannot be
#'   interrupted, so there are two events only: \code{unit = "run"},
#'   \code{i = 1} (\code{"start"}) before the run and \code{i = 2}
#'   (\code{"done"}, \code{payload$n_blocks}) after it; returning \code{FALSE}
#'   at \code{"start"} cancels with an \code{hs_error_cancelled} error (no
#'   partial result).
#'
#' @return With \code{attach = TRUE}, the \code{GenoMeta} with
#'   \code{block_marker_map_info} / \code{block_info} set. With
#'   \code{attach = FALSE}, a list with \code{block_marker_map}
#'   (\code{Chrom, Marker, Block}) and \code{block_info}
#'   (\code{Chrom, Block, Start, End, Pos, Num_SNPs, Block_length}), covering
#'   the blocked markers only. Block IDs are SelectionTools' own
#'   (\code{b000000}, \code{b000001}, ...), numbered along the genome.
#'
#' @seealso \code{\link{hapS_ld_blocks}} (the default resolver),
#'   \code{\link{attach_block_info}} (the recommended route: a block map you
#'   already have), \code{\link{st_ld_blocks_files}} (SelectionTools files in),
#'   \code{\link{build_block_location_info}}
#'
#' @examples
#' \dontrun{
#' # The recommended route needs no LD blocking at all:
#' meta <- read_plink_ped("geno.ped", "geno.map", block_file = "blocks.tsv")
#'
#' # Building blocks locally with SelectionTools instead:
#' meta <- read_plink_ped("geno.ped", "geno.map")
#' meta <- st_ld_blocks(meta, ld_threshold = 0.5, tolerance = 0)
#'
#' # the tables only, plus the canonical files
#' blk <- st_ld_blocks(meta, 0.5, 0, attach = FALSE, output_dir = "blocks")
#' }
#'
#' @export
st_ld_blocks <- function(geno_meta,
                         ld_threshold = 0.5,
                         tolerance    = 0L,
                         attach       = TRUE,
                         output_dir   = NULL,
                         prefix       = "st_blocks",
                         keep_files   = FALSE,
                         verbose      = TRUE,
                         st_formats   = NULL,
                         progress     = NULL) {
  caller <- "st_ld_blocks"
  .st_gate(caller, st_formats)
  trk <- .hs_progress_tracker(progress, "ld_blocks", "run", 2L)
  .hs_check_meta(geno_meta, caller, "geno_meta")
  .ld_check_rule(ld_threshold, tolerance, caller, zero_note = .st_zero_note)
  .st_require(caller)
  d <- dim(geno_meta$geno_arrays)
  if (d[3L] != 2L)
    .hs_abort(caller, ": SelectionTools reads diploid genotypes; this container ",
         "has ploidy ", d[3L], ".", call. = FALSE)
  if (d[1L] < 2L)
    .hs_abort(caller, ": LD needs at least 2 samples.", call. = FALSE)

  placed <- .ld_placed_markers(geno_meta$map_info)
  n_unplaced <- d[2L] - length(placed)
  if (length(placed) < 2L)
    .hs_abort(caller, ": fewer than 2 markers carry a usable position.",
         call. = FALSE)

  t0 <- Sys.time()
  work_dir <- tempfile("hs_st_")
  dir.create(work_dir, recursive = TRUE)
  if (isTRUE(keep_files)) {
    .hs_inform(caller, ": SelectionTools files kept in ", work_dir)
  } else {
    on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }

  io <- .st_write_inputs(geno_meta, placed, work_dir, caller = caller)

  # SelectionTools reads and writes through st.input.dir / st.output.dir,
  # which point at the private directory for this call (0.9.0: the working
  # directory is not changed).
  trk$tick(1L, label = "start")
  raw <- .st_with_io_dirs(work_dir,
                          .st_run_hblocks(basename(io$geno_file), basename(io$map_file),
                                          ld_threshold, tolerance, verbose = verbose))

  bm <- .st_parse_blocks(raw, io$markers, caller = caller)
  trk$tick(2L, label = "done",
           payload = list(n_blocks = length(unique(bm$Block))))

  if (isTRUE(verbose))
    .hs_inform(sprintf(paste0(
      "%s: SelectionTools blocks (r^2 >= %s, tolerance %d): %s; %d unplaced ",
      "marker(s) left unblocked%s; %.1f s"),
      caller, format(ld_threshold), as.integer(tolerance), .ld_sizes_note(bm),
      n_unplaced,
      if (io$n_left_out) sprintf(", %d with > 5 alleles", io$n_left_out) else "",
      as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  .ld_finish(geno_meta, bm,
             data.frame(Marker = io$markers$Marker, Chrom = io$markers$Chrom,
                        Pos = io$markers$Pos, stringsAsFactors = FALSE),
             attach, output_dir, prefix, caller,
             rule = list(backend = "st_ld_blocks", ld_threshold = ld_threshold,
                         tolerance = tolerance,
                         n_unplaced = n_unplaced + io$n_left_out))
}


#' LD haplotype blocks from SelectionTools files (experimental)
#'
#' The file route of the SelectionTools resolver, for genotype and map files
#' already in SelectionTools' own formats: \code{\link{st_convert}} reads them
#' into a \code{GenoMeta}, \code{\link{st_ld_blocks}} blocks it, and the two
#' canonical block files are written to \code{output_dir}. Named
#' \code{define_haplotype_blocks()} up to 0.7.9; that name still works, with a
#' deprecation warning. SelectionTools is being retired from haploStack:
#' \strong{this function is deprecated since 0.9.0} (it warns on every call
#' and goes in 0.10; it is \code{\link{st_convert}} followed by
#' \code{\link{st_ld_blocks}}), it is behind \code{\link{use_st_formats}}
#' (off by default), and
#' SelectionTools' own side files (\code{<tag>.ld.txt},
#' \code{<tag>.haplotype_map.txt}, written up to 0.8.8) are no longer produced.
#' Before using it, check whether you need LD blocking at all: the recommended
#' route is a block map you already have (\code{\link{attach_block_info}}), and
#' the default local resolver is \code{\link{hapS_ld_blocks}}.
#'
#' Two files are written to \code{output_dir}, named after the genotype file
#' and the rule (\code{<tag>} = \code{<geno>_<ld_threshold>_t<tolerance>}):
#' \code{<tag>_block_info.csv} and \code{<tag>_block_marker_map.csv}, the
#' canonical tables (see \code{\link{build_block_location_info}}).
#'
#' @param geno_data_file Path to a genotype file in SelectionTools' own matrix
#'   format (\code{st.read.marker.data(format = "m")}): markers in rows, the
#'   first row the individual names, fields separated by a tab, a blank or
#'   \code{;} (not a comma). Starting from a \code{GenoMeta}, use
#'   \code{\link{st_ld_blocks}}.
#' @param geno_map_file  Path to the map, SelectionTools format \code{"mcp"}
#'   with a header line: columns \code{name, chrom, pos}. The returned tables
#'   take \code{Chrom} (e.g. \code{"1H"}) from this file.
#' @inheritParams st_ld_blocks
#' @param output_dir     Directory where the two block files are written.
#'
#' @return A named list with \code{block_marker_map} (\code{Chrom, Marker,
#'   Block}), \code{block_info} (\code{Chrom, Block, Start, End, Pos, Num_SNPs,
#'   Block_length}) and \code{provenance}, as \code{st_ld_blocks(attach =
#'   FALSE)} returns it.
#'
#' @details
#' Requires the \pkg{SelectionTools} package. If not installed, an informative
#' error is raised so the user can supply a pre-computed block map instead.
#'
#' @seealso \code{\link{st_ld_blocks}} (start from a \code{GenoMeta}),
#'   \code{\link{st_convert}} (the format converter),
#'   \code{\link{hapS_ld_blocks}} (the default resolver),
#'   \code{\link{attach_block_info}}, \code{\link{build_block_location_info}}
#'
#' @export
st_ld_blocks_files <- function(geno_data_file,
                               geno_map_file,
                               ld_threshold = 0.5,
                               tolerance    = 0,
                               output_dir   = "haplotype",
                               st_formats   = NULL,
                               progress     = NULL) {
  # Deprecated since 0.9.0 (row 27 a): st_convert() then st_ld_blocks().
  .hs_deprecated("st_ld_blocks_files()", "st_convert() followed by st_ld_blocks()",
                 "0.9.0")
  caller <- "st_ld_blocks_files"
  .st_gate(caller, st_formats)
  .ld_check_rule(ld_threshold, tolerance, caller, zero_note = .st_zero_note)
  .st_require(caller)
  tag  <- paste0(gsub("\\..*$", "", basename(geno_data_file)), "_",
                 ld_threshold, "_t", tolerance)
  meta <- st_convert(geno_data_file, map = geno_map_file, st_formats = TRUE)
  st_ld_blocks(meta, ld_threshold = ld_threshold, tolerance = tolerance,
               attach = FALSE, output_dir = output_dir, prefix = tag,
               verbose = FALSE, st_formats = TRUE, progress = progress)
}
