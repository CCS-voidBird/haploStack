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
    stop(caller, " requires the SelectionTools package, which is not ",
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
# functions find them there. The bridge passes bare names inside its private
# working directory, so both must be "".
#   * Set to a directory: refused. SelectionTools would read and write there,
#     outside the private directory.
#   * Gone: every st.* reader stops with "object 'st.input.dir' not found".
#     That happens after rm(list = ls()) once SelectionTools has loaded, which
#     is what a script re-run in the same session does. The bridge supplies ""
#     for the call and removes it again on exit (.st_drop_io_dirs()), so the
#     global environment is left as it was found.
# Returns the names it supplied.
.st_check_io_dirs <- function(caller) {
  ns   <- asNamespace("SelectionTools")
  vars <- c("st.input.dir", "st.output.dir")
  have <- vapply(vars, exists, logical(1L), envir = ns)   # SelectionTools' lookup
  dirs <- vapply(vars, function(v)
    if (exists(v, envir = ns)) paste(as.character(get(v, envir = ns)), collapse = "")
    else "", "")
  if (any(nzchar(dirs)))
    stop(caller, ": SelectionTools' ", paste(vars[nzchar(dirs)], collapse = " / "),
         " is set (", paste(dirs[nzchar(dirs)], collapse = ", "), "). The bridge ",
         "keeps SelectionTools' files in a private directory and needs the ",
         "default: run ", paste0(vars[nzchar(dirs)], " <- \"\"", collapse = "; "),
         " and call again.", call. = FALSE)
  supplied <- vars[!have]
  for (v in supplied) assign(v, "", envir = globalenv())
  invisible(supplied)
}

# Undo what .st_check_io_dirs() supplied (on exit).
.st_drop_io_dirs <- function(supplied) {
  if (length(supplied))
    suppressWarnings(rm(list = supplied, envir = globalenv()))
  invisible(NULL)
}

# Write the two SelectionTools input files for `markers` (placed, any order)
# into `dir`, under aliases. Returns the files and the alias tables.
#
#   hs_geno.txt  matrix format "m": header "m/i  i1 i2 ...", one line per
#                marker "m7  1/2 2/2 -/1 ...", tab-separated. Alleles are
#                recoded per marker to 1, 2, ... in sorted order — r^2 does not
#                depend on which allele is called what, and SelectionTools
#                reads 1..5 whatever the source vocabulary was (nucleotides,
#                1/2/3/4 codes, DArT A/B). Missing = "-".
#   hs_map.txt   map format "mcp" with a header: "name chrom pos", chrom as an
#                integer key in natural order, pos in bp.
.st_write_inputs <- function(geno_meta, markers, dir) {
  map <- as.data.frame(geno_meta$map_info, stringsAsFactors = FALSE)
  map <- map[match(markers, map$Marker), c("Chrom", "Marker", "Pos"), drop = FALSE]
  map$Chrom <- as.character(map$Chrom)
  map$Pos   <- as.numeric(map$Pos)

  chrom_levels <- .natural_chrom_levels(map$Chrom)
  map$key <- match(map$Chrom, chrom_levels)
  map <- map[order(map$key, map$Pos, seq_len(nrow(map))), , drop = FALSE]

  arr <- geno_meta$geno_arrays[, map$Marker, , drop = FALSE]
  a1  <- matrix(as.character(arr[, , 1L]), nrow = dim(arr)[1L])
  a2  <- matrix(as.character(arr[, , 2L]), nrow = dim(arr)[1L])
  c1  <- matrix(NA_integer_, nrow(a1), ncol(a1))
  c2  <- c1
  n_all <- integer(ncol(a1))
  for (j in seq_len(ncol(a1))) {
    lv <- sort(unique(c(a1[, j], a2[, j])))
    lv <- lv[!is.na(lv)]
    n_all[j] <- length(lv)
    c1[, j] <- match(a1[, j], lv)
    c2[, j] <- match(a2[, j], lv)
  }
  # SelectionTools codes alleles 1..5; a marker with more states cannot be
  # written faithfully and is left out (it becomes a singleton on attach).
  too_many <- n_all > 5L
  if (any(too_many)) {
    warning("st_ld_blocks: ", sum(too_many), " marker(s) carry ",
            "more than 5 allele states and are left out of SelectionTools' ",
            "blocking (e.g. ", paste(head(map$Marker[too_many], 3L),
            collapse = ", "), ").", call. = FALSE)
    keep <- !too_many
    map <- map[keep, , drop = FALSE]
    c1 <- c1[, keep, drop = FALSE]
    c2 <- c2[, keep, drop = FALSE]
  }

  m_alias <- paste0("m", seq_len(nrow(map)))
  i_alias <- paste0("i", seq_len(nrow(c1)))
  cell <- function(x) { x <- as.character(x); x[is.na(x)] <- "-"; x }
  cells <- matrix(paste0(cell(c1), "/", cell(c2)), nrow = nrow(c1))

  geno_file <- file.path(dir, "hs_geno.txt")
  map_file  <- file.path(dir, "hs_map.txt")
  body <- data.table::as.data.table(t(cells))
  data.table::setnames(body, i_alias)
  body <- cbind(data.table::data.table(`m/i` = m_alias), body)
  data.table::fwrite(body, geno_file, sep = "\t", quote = FALSE,
                     col.names = TRUE, showProgress = FALSE)
  # scipen: positions in fixed notation — SelectionTools' own map parser must
  # never be handed "8.8e+07" for a round position.
  data.table::fwrite(data.table::data.table(name = m_alias, chrom = map$key,
                                            pos = map$Pos),
                     map_file, sep = "\t", quote = FALSE, col.names = TRUE,
                     scipen = 100L, showProgress = FALSE)

  list(geno_file = geno_file, map_file = map_file,
       markers   = data.frame(alias = m_alias, Marker = map$Marker,
                              Chrom = map$Chrom, Pos = map$Pos,
                              stringsAsFactors = FALSE),
       n_left_out = sum(too_many))
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

# st.def.hblocks() table -> (Chrom, Marker, Block) in container IDs and names.
.st_parse_blocks <- function(raw, alias_map) {
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  if (!all(c("Name", "Markers") %in% names(raw)))
    stop("st_ld_blocks: unexpected SelectionTools output ",
         "(columns: ", paste(names(raw), collapse = ", "), ").", call. = FALSE)
  members <- strsplit(gsub("[[:space:]]", "", raw$Markers), ";", fixed = TRUE)
  members <- lapply(members, function(x) x[nzchar(x)])
  bm <- data.frame(alias = unlist(members, use.names = FALSE),
                   Block = rep(as.character(raw$Name), lengths(members)),
                   stringsAsFactors = FALSE)
  idx <- match(bm$alias, alias_map$alias)
  if (anyNA(idx))
    stop("st_ld_blocks: SelectionTools returned marker(s) the ",
         "bridge did not write (", paste(head(bm$alias[is.na(idx)], 3L),
         collapse = ", "), ").", call. = FALSE)
  out <- data.frame(Chrom  = alias_map$Chrom[idx],
                    Marker = alias_map$Marker[idx],
                    Block  = bm$Block,
                    stringsAsFactors = FALSE)
  out[order(match(out$Marker, alias_map$Marker)), , drop = FALSE]
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
#' works, with a deprecation warning.
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
#' \code{st.output.dir}. It creates both, as \code{""}, in the global
#' environment when it loads. Both must be \code{""} here: a directory set there
#' is refused, and the error names the remedy. If they are gone — after
#' \code{rm(list = ls())}, for instance — \code{""} is supplied for the call and
#' removed again afterwards, so the global environment is left as it was found.
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
                         verbose      = TRUE) {
  caller <- "st_ld_blocks"
  if (!is_geno_meta(geno_meta, verbose = FALSE))
    stop(caller, ": 'geno_meta' must be a valid GenoMeta.", call. = FALSE)
  .ld_check_rule(ld_threshold, tolerance, caller, zero_note = .st_zero_note)
  .st_require(caller)
  st_dirs <- .st_check_io_dirs(caller)
  on.exit(.st_drop_io_dirs(st_dirs), add = TRUE)
  d <- dim(geno_meta$geno_arrays)
  if (d[3L] != 2L)
    stop(caller, ": SelectionTools reads diploid genotypes; this container ",
         "has ploidy ", d[3L], ".", call. = FALSE)
  if (d[1L] < 2L)
    stop(caller, ": LD needs at least 2 samples.", call. = FALSE)

  placed <- .ld_placed_markers(geno_meta$map_info)
  n_unplaced <- d[2L] - length(placed)
  if (length(placed) < 2L)
    stop(caller, ": fewer than 2 markers carry a usable position.",
         call. = FALSE)

  t0 <- Sys.time()
  work_dir <- tempfile("hs_st_")
  dir.create(work_dir, recursive = TRUE)
  if (isTRUE(keep_files)) {
    message(caller, ": SelectionTools files kept in ", work_dir)
  } else {
    on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }

  io <- .st_write_inputs(geno_meta, placed, work_dir)

  # SelectionTools writes its side files into the working directory.
  old_wd <- setwd(work_dir)
  on.exit(setwd(old_wd), add = TRUE, after = FALSE)
  raw <- .st_run_hblocks(basename(io$geno_file), basename(io$map_file),
                         ld_threshold, tolerance, verbose = verbose)
  setwd(old_wd)

  bm <- .st_parse_blocks(raw, io$markers)

  if (isTRUE(verbose))
    message(sprintf(paste0(
      "%s: SelectionTools blocks (r^2 >= %s, tolerance %d): %s; %d unplaced ",
      "marker(s) left unblocked%s; %.1f s"),
      caller, format(ld_threshold), as.integer(tolerance), .ld_sizes_note(bm),
      n_unplaced,
      if (io$n_left_out) sprintf(", %d with > 5 alleles", io$n_left_out) else "",
      as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  .ld_finish(geno_meta, bm,
             data.frame(Marker = io$markers$Marker, Chrom = io$markers$Chrom,
                        Pos = io$markers$Pos, stringsAsFactors = FALSE),
             attach, output_dir, prefix, caller)
}


#' LD haplotype blocks from SelectionTools files (experimental)
#'
#' The file route of the SelectionTools resolver, for genotype and map files
#' already in SelectionTools' own formats. Named
#' \code{define_haplotype_blocks()} up to 0.7.9; that name still works, with a
#' deprecation warning. For data in a \code{GenoMeta} use
#' \code{\link{st_ld_blocks}}. Before either, check whether you need LD
#' blocking at all: the recommended route is a block map you already have
#' (\code{\link{attach_block_info}}), and the default local resolver is
#' \code{\link{hapS_ld_blocks}}. Experimental, for the reasons
#' \code{\link{st_ld_blocks}} gives.
#'
#' Four files are written to \code{output_dir}, named after the genotype file
#' and the rule (\code{<tag>} = \code{<geno>_<ld_threshold>_t<tolerance>}):
#' \describe{
#'   \item{\code{<tag>_block_info.csv}}{Block coordinate table:
#'     \code{Chrom, Block, Start, End}.}
#'   \item{\code{<tag>_block_marker_map.csv}}{Canonical marker-level table:
#'     \code{Chrom, Marker, Block}.}
#'   \item{\code{<tag>.ld.txt}}{SelectionTools' block table, as
#'     \code{st.def.hblocks()} returns it.}
#'   \item{\code{<tag>.haplotype_map.txt}}{Per-marker haplotype-allele map
#'     from \code{st.recode.hil()}.}
#' }
#'
#' @param geno_data_file Path to a genotype file in SelectionTools' own matrix
#'   format (\code{st.read.marker.data(format = "m")}): markers in rows, the
#'   first row the individual names, fields separated by a tab, a blank or
#'   \code{;} — **not a comma**, so the CSV that
#'   \code{write_geno_meta(output_format = "st")} writes is not one. Starting
#'   from a \code{GenoMeta}, use \code{\link{st_ld_blocks}}.
#' @param geno_map_file  Path to the map, SelectionTools format \code{"mcp"}
#'   with a header line: columns \code{name, chrom, pos}. SelectionTools keeps
#'   only the leading integer of a chromosome name; the returned tables take
#'   \code{Chrom} (e.g. \code{"1H"}) from this file.
#' @param ld_threshold   Minimum LD (\code{r^2}) to merge adjacent markers; one
#'   number in \eqn{(0, 1]}. Default \code{0.5}. \code{0} is refused:
#'   SelectionTools aborts the R session on it.
#' @param tolerance      Tolerance for \code{st.def.hblocks()}; a whole number
#'   \eqn{\ge 0}. Default \code{0}.
#' @param output_dir     Directory where block-map files are written.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{block_info}}{data.table: \code{Chrom, Block, Start, End}.}
#'   \item{\code{block_marker_map}}{data.table: \code{Chrom, Marker, Block}.}
#' }
#'
#' @details
#' Requires the \pkg{SelectionTools} package. If not installed, an informative
#' error is raised so the user can supply a pre-computed block map instead.
#'
#' SelectionTools runs in a private temporary directory (its side files —
#' \code{ld}, \code{blocks.hmap}, \code{blocks.hall}, the \code{.srt} map — never
#' touch the working directory or the input folder) and on the data set
#' \code{"ld.blocks"}; its \code{"default"} data set is left alone (0.7.9).
#' SelectionTools' \code{st.input.dir} / \code{st.output.dir} are handled as
#' \code{\link{st_ld_blocks}} describes.
#'
#' After calling this function, pass the \code{block_marker_map} element of
#' the returned list to \code{\link{attach_block_info}} to register the block
#' assignments in a \code{GenoMeta} object:
#' \preformatted{
#'   blocks <- st_ld_blocks_files("geno.txt", "map.txt")
#'   meta   <- attach_block_info(meta, blocks$block_marker_map)
#' }
#'
#' @seealso \code{\link{st_ld_blocks}} (start from a \code{GenoMeta}),
#'   \code{\link{hapS_ld_blocks}} (the default resolver),
#'   \code{\link{attach_block_info}}, \code{\link{run_haplotype_catalogue}},
#'   \code{\link{build_block_location_info}}, \code{\link{expand_block_info}}
#'
#' @export
st_ld_blocks_files <- function(geno_data_file,
                               geno_map_file,
                               ld_threshold = 0.5,
                               tolerance    = 0,
                               output_dir   = "haplotype") {
  caller <- "st_ld_blocks_files"
  .ld_check_rule(ld_threshold, tolerance, caller, zero_note = .st_zero_note)
  .st_require(caller)
  st_dirs <- .st_check_io_dirs(caller)
  on.exit(.st_drop_io_dirs(st_dirs), add = TRUE)

  geno_data_file <- normalizePath(geno_data_file, winslash = "/", mustWork = TRUE)
  geno_map_file  <- normalizePath(geno_map_file,  winslash = "/", mustWork = TRUE)
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  geno_basename <- gsub("\\..*$", "", basename(geno_data_file))
  tag           <- paste0(geno_basename, "_", ld_threshold, "_t", tolerance)

  # SelectionTools writes side files into the working directory ("ld",
  # "blocks.hmap", "blocks.hall") and next to the map it reads ("<map>.srt").
  # Run it in a private directory holding its own copy of the map, so neither
  # the user's working directory nor the input folder is touched, and two runs
  # cannot overwrite each other's files.
  work_dir <- tempfile("hs_st_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  file.copy(geno_map_file, file.path(work_dir, "hs_map.txt"))
  old_wd <- setwd(work_dir)
  on.exit(setwd(old_wd), add = TRUE, after = FALSE)

  # Load straight into the dedicated data set: reading into "default" (and
  # copying) replaced whatever a SelectionTools session kept there.
  ds <- "ld.blocks"
  SelectionTools::st.read.marker.data(geno_data_file, format = "m",
                                 data.set = ds)
  SelectionTools::st.read.map("hs_map.txt",
                               m.stretch  = 1,
                               format     = "mcp",
                               skip       = 1,
                               data.set   = ds)

  message("Computing pairwise LD (r²)...")
  # auxfiles = FALSE: st.def.hblocks() reads the LD from SelectionTools'
  # memory; writing the whole pairwise table to disk only to read it back into
  # R (and discard it) was the slowest step on a genome-wide panel.
  SelectionTools::st.calc.ld(ld.measure = "r2",
                              data.set   = ds,
                              auxfiles   = FALSE)

  message("Defining haplotype blocks (LD threshold = ", ld_threshold,
          ", tolerance = ", tolerance, ")...")
  block_info_raw <- SelectionTools::st.def.hblocks(
    ld.threshold = ld_threshold,
    ld.criterion = "flanking",
    tolerance    = as.integer(tolerance),
    data.set     = ds
  )
  colnames(block_info_raw)[3L] <- "Block"

  # Write raw LD block file
  utils::write.table(block_info_raw,
    file = file.path(output_dir, paste0(tag, ".ld.txt")),
    sep = " ", row.names = FALSE, quote = FALSE)

  # Build full haplotype map
  hrec <- SelectionTools::st.recode.hil(data.set = ds)
  setwd(old_wd)
  hrec$Block_id <- paste0(hrec$Block, "_", hrec$AlleleNr)
  block_info_full <- dplyr::left_join(hrec, block_info_raw, by = "Block")

  setDT(block_info_full)
  block_info_map <- block_info_full[, .(
    Marker     = unlist(strsplit(gsub(" ", "", Markers), ";")),
    MarkerName = paste0(
      unlist(strsplit(gsub(" ", "", Markers), ";")), ".",
      unlist(strsplit(gsub(" ", "", AlleleDef), ";"))),
    Pos   = rep(Pos, lengths(strsplit(gsub(" ", "", Markers), ";"))),
    Block = rep(Block, lengths(strsplit(gsub(" ", "", Markers), ";")))
  ), by = .(Chrom, Block_id)][, c("Chrom","Marker","MarkerName",
                                   "Block","Block_id","Pos")]

  # Read map file for precise positions — and for the chromosome NAMES.
  # SelectionTools keeps only the leading integer of a chromosome ("1H" -> 1),
  # so joining its Chrom back to a named map matched nothing and every
  # position came out NA (up to 0.7.8). Join on Marker alone and take Chrom
  # from the map. Marker and Chrom are read as text (a numeric-looking marker
  # ID must stay the string SelectionTools returned); Pos keeps fread's type,
  # so the written files format positions exactly as before.
  map_info <- fread(geno_map_file, header = TRUE,
                    colClasses = list(character = 1:2))
  map_info <- map_info[, 1:3]
  colnames(map_info) <- c("Marker", "Chrom", "Pos")

  block_info_map <- dplyr::left_join(
    block_info_map[, c("Marker","MarkerName","Block","Block_id")],
    map_info, by = "Marker", relationship = "many-to-many"
  )[, c("Chrom","Marker","MarkerName","Block","Block_id","Pos")]

  fwrite(as.data.table(block_info_map),
         file = file.path(output_dir, paste0(tag, ".haplotype_map.txt")),
         sep = " ", row.names = FALSE, quote = FALSE)

  # Build canonical block_marker_map: (Chrom, Marker, Block)
  block_marker_map <- unique(
    as.data.table(block_info_map)[, .(Chrom, Marker, Block)]
  )

  # Build block_info: (Chrom, Block, Start, End) from marker positions
  block_info <- block_marker_map[
    as.data.table(map_info)[, .(Marker, Pos)], on = "Marker", nomatch = 0
  ][, .(Start = min(Pos), End = max(Pos)), by = .(Chrom, Block)]

  # Write both canonical files
  binfo_file  <- file.path(output_dir, paste0(tag, "_block_info.csv"))
  bmap_file   <- file.path(output_dir, paste0(tag, "_block_marker_map.csv"))
  fwrite(block_info,      file = binfo_file,  sep = ",", row.names = FALSE, quote = FALSE)
  fwrite(block_marker_map, file = bmap_file,  sep = ",", row.names = FALSE, quote = FALSE)
  message("Block files written to: ", output_dir)

  list(block_info = block_info, block_marker_map = block_marker_map)
}
