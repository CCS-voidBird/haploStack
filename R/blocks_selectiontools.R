# =============================================================================
# blocks_selectiontools.R — LD haplotype blocks for a GenoMeta via SelectionTools
#
# define_haplotype_blocks_meta() is the GenoMeta-in route to SelectionTools'
# LD blocking (st.calc.ld + st.def.hblocks, "flanking" criterion). It exists
# because the file route could not be driven from a container (INDEX
# discrepancy 34):
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


# Fail early, with the install route, when SelectionTools is absent.
.st_require <- function(caller) {
  if (!requireNamespace("SelectionTools", quietly = TRUE))
    stop(caller, " requires the SelectionTools package, which is not ",
         "installed. Alternatives: HapSelect through ",
         "dev/hapselect_bridge/hapselect_bridge.R, or a pre-computed block map ",
         "via attach_block_info().", call. = FALSE)
  invisible(TRUE)
}

# The blocking rule SelectionTools can run. ld.threshold = 0 kills the R
# session inside st.def.hblocks() (C level, no error), so it is refused here.
.st_check_rule <- function(ld_threshold, tolerance, caller) {
  if (!is.numeric(ld_threshold) || length(ld_threshold) != 1L ||
      is.na(ld_threshold) || ld_threshold <= 0 || ld_threshold > 1)
    stop(caller, ": 'ld_threshold' must be one number in (0, 1] ",
         "(an r^2). 0 is refused: SelectionTools aborts the R session on it.",
         call. = FALSE)
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) ||
      tolerance < 0 || tolerance != round(tolerance))
    stop(caller, ": 'tolerance' must be one whole number >= 0 (the number ",
         "of sub-threshold markers a block may absorb).", call. = FALSE)
  invisible(TRUE)
}

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

# Markers SelectionTools can place: a finite position > 0 on a named,
# placed chromosome. It orders markers by map position, so an unplaced marker
# would be blocked by file order.
.st_placed_markers <- function(map_info) {
  chrom <- trimws(as.character(map_info$Chrom))
  pos   <- suppressWarnings(as.numeric(map_info$Pos))
  ok <- !is.na(chrom) & nzchar(chrom) &
        !tolower(chrom) %in% c("0", "un", "unknown", "chrun", "unplaced") &
        is.finite(pos) & pos > 0
  as.character(map_info$Marker)[ok]
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
    warning("define_haplotype_blocks_meta: ", sum(too_many), " marker(s) carry ",
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
    stop("define_haplotype_blocks_meta: unexpected SelectionTools output ",
         "(columns: ", paste(names(raw), collapse = ", "), ").", call. = FALSE)
  members <- strsplit(gsub("[[:space:]]", "", raw$Markers), ";", fixed = TRUE)
  members <- lapply(members, function(x) x[nzchar(x)])
  bm <- data.frame(alias = unlist(members, use.names = FALSE),
                   Block = rep(as.character(raw$Name), lengths(members)),
                   stringsAsFactors = FALSE)
  idx <- match(bm$alias, alias_map$alias)
  if (anyNA(idx))
    stop("define_haplotype_blocks_meta: SelectionTools returned marker(s) the ",
         "bridge did not write (", paste(head(bm$alias[is.na(idx)], 3L),
         collapse = ", "), ").", call. = FALSE)
  out <- data.frame(Chrom  = alias_map$Chrom[idx],
                    Marker = alias_map$Marker[idx],
                    Block  = bm$Block,
                    stringsAsFactors = FALSE)
  out[order(match(out$Marker, alias_map$Marker)), , drop = FALSE]
}


#' LD haplotype blocks for a GenoMeta, with SelectionTools
#'
#' The GenoMeta-in route to SelectionTools' LD blocking. Give it a container,
#' an r\eqn{^2} threshold and a tolerance; it writes the container into
#' SelectionTools' own input formats in a private temporary directory, runs
#' \code{st.calc.ld()} + \code{st.def.hblocks(ld.criterion = "flanking")}
#' there, reads the blocks back in the container's own marker IDs and
#' chromosome names, deletes the directory, and (by default) attaches the
#' blocks to the container.
#'
#' Use this rather than \code{\link{define_haplotype_blocks}} when the data is
#' already a \code{GenoMeta}: that function takes SelectionTools files, and
#' haploStack's own "ST" export (\code{write_geno_meta(output_format = "st")})
#' is a CSV that SelectionTools cannot read.
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
#' @seealso \code{\link{define_haplotype_blocks}} (the file route),
#'   \code{\link{attach_block_info}}, \code{\link{build_block_location_info}}
#'
#' @examples
#' \dontrun{
#' meta <- read_plink_ped("geno.ped", "geno.map")
#' meta <- define_haplotype_blocks_meta(meta, ld_threshold = 0.5, tolerance = 3)
#' head(meta$block_marker_map_info)
#'
#' # the tables only, plus the canonical files
#' blk <- define_haplotype_blocks_meta(meta, 0.5, 3, attach = FALSE,
#'                                     output_dir = "blocks")
#' }
#'
#' @export
define_haplotype_blocks_meta <- function(geno_meta,
                                         ld_threshold = 0.5,
                                         tolerance    = 0L,
                                         attach       = TRUE,
                                         output_dir   = NULL,
                                         prefix       = "st_blocks",
                                         keep_files   = FALSE,
                                         verbose      = TRUE) {
  caller <- "define_haplotype_blocks_meta"
  if (!is_geno_meta(geno_meta, verbose = FALSE))
    stop(caller, ": 'geno_meta' must be a valid GenoMeta.", call. = FALSE)
  .st_check_rule(ld_threshold, tolerance, caller)
  .st_require(caller)
  st_dirs <- .st_check_io_dirs(caller)
  on.exit(.st_drop_io_dirs(st_dirs), add = TRUE)
  d <- dim(geno_meta$geno_arrays)
  if (d[3L] != 2L)
    stop(caller, ": SelectionTools reads diploid genotypes; this container ",
         "has ploidy ", d[3L], ".", call. = FALSE)
  if (d[1L] < 2L)
    stop(caller, ": LD needs at least 2 samples.", call. = FALSE)

  placed <- .st_placed_markers(geno_meta$map_info)
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
  loc <- .compute_block_location(
    bm, data.frame(Marker = io$markers$Marker, Chrom = io$markers$Chrom,
                   Pos = io$markers$Pos, stringsAsFactors = FALSE))
  bi <- loc$block_info

  if (isTRUE(verbose)) {
    sizes <- table(bm$Block)
    message(sprintf(paste0(
      "%s: SelectionTools blocks (r^2 >= %s, tolerance %d): %d blocks over %d ",
      "markers (%d single-marker, largest %d); %d unplaced marker(s) left ",
      "unblocked%s; %.1f s"),
      caller, format(ld_threshold), as.integer(tolerance), length(sizes),
      nrow(bm), sum(sizes == 1L), max(sizes), n_unplaced,
      if (io$n_left_out) sprintf(", %d with > 5 alleles", io$n_left_out) else "",
      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  if (!isTRUE(attach)) {
    if (!is.null(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      # scipen: a round position (88000000) must not be written as 8.8e+07
      data.table::fwrite(data.table::as.data.table(bi),
                         file.path(output_dir, paste0(prefix, "_block_info.csv")),
                         scipen = 100L)
      data.table::fwrite(data.table::as.data.table(bm),
                         file.path(output_dir, paste0(prefix, "_block_marker_map.csv")))
    }
    return(list(block_marker_map = data.table::as.data.table(bm),
                block_info       = data.table::as.data.table(bi)))
  }

  out <- attach_block_info(geno_meta, bm, source_type = "block_marker_map")
  # attach_block_info() rebuilds the container and drops these (INDEX
  # discrepancy 20); the samples and markers are unchanged, so carry them over.
  for (f in c("marker_effects", "phenotype", "sample_info"))
    if (!is.null(geno_meta[[f]])) out[[f]] <- geno_meta[[f]]
  if (!is_geno_meta(out, verbose = TRUE))
    stop(caller, ": the container failed validation after attaching blocks.",
         call. = FALSE)
  if (!is.null(output_dir))
    invisible(build_block_location_info(out, output_dir = output_dir,
                                        prefix = prefix))
  out
}
