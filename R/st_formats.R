# =============================================================================
# st_formats.R — SelectionTools' own file formats: one converter, one switch
# (0.9.0, owner steer 2026-10-29)
#
# SelectionTools is being retired from haploStack. Its formats are still worth
# reading and writing — users hold data in them — but the core never sees
# them: st_convert() turns them into the shapes the core defines (a GenoMeta,
# or a block-marker map for attach_block_info()) and writes a GenoMeta out in
# them, and nothing else in the package parses or emits a SelectionTools
# format. Every entry point that processes one — st_convert(), st_ld_blocks(),
# st_ld_blocks_files() and the deprecated names that call them — is behind one
# global switch, off by default:
#
#   use_st_formats(TRUE)              this session
#   options(haploStack.st_formats =)  the same, by hand (validated on read)
#   st_formats = TRUE                 one call (framework F1: argument > option)
#
# The formats (SelectionTools' manual):
#   matrix "m"   genotype: header "m/i <ind> <ind> ...", one line per marker
#                "<marker> 1/2 2/2 -/1 ...", fields separated by a tab, a
#                blank or ";" (not a comma); alleles 1..5, "-" missing
#   map "mcp"    "name chrom pos" with a header line; SelectionTools keeps
#                only the leading integer of a chromosome name ("1H" -> 1)
#   blocks       st.def.hblocks()'s one-row-per-block table: Name = block ID,
#                Markers = "m1;m2;..."
#
# The "ST" PRESET of read_genotype() (markers in rows, 1/4 codes) is a generic
# delimited layout, but its cells are SelectionTools' allele codes, so it is
# gated too and its cells are translated on the way in (C-5, below).
#
# The allele encoding (C-5, owner clarification 2026-10-29). The real
# difficulty of the SelectionTools format is not its layout but its
# numeric-categorical allele encoding with a slash separator — A = 1, C = 2,
# G = 3, T = 4, missing = -1, so the diploid genotype "A-" is written "1/-1".
# Up to 0.8.8 that was a vocabulary the container could carry
# (metadata$encoding "selectiontools"). Since 0.9.0 a GenoMeta stores
# nucleotide alleles only (plus dosage): the code table lives in this file and
# is applied at the boundary — st_convert() and the "ST" preset translate on
# the way in, the writer recodes on the way out — and the portal
# (new_geno_meta()) turns a declared "selectiontools" encoding into
# nucleotides with a deprecation warning, refuses undeclared codes with a
# pointer here, and upgrade_geno_meta() recodes a container saved with the
# old vocabulary.
#
# The switch covers the token as well as the files (C-6, owner 2026-11-03):
# encoding = "selectiontools" on build_geno_meta() / new_geno_meta() / the
# sim_* generators, a "selectiontools" side of translate_alleles() /
# transform_effects(), use_missing_chars(mode = "selectiontools") and the "ST"
# writer preset (geno_meta_to_matrix(), write_geno_meta()) all stop with the
# use_st_formats(TRUE) hint while it is off. The one deliberate exception is
# upgrade_geno_meta() (so load_geno_meta() and recode_geno_meta() on a container
# saved with the old vocabulary): a one-way migration that announces itself and
# must work without a switch, or an 0.8.8 .rds could not be opened. Every
# recode goes through .st_translate() with the two tables below.
# =============================================================================


#' Switch SelectionTools formats on or off for the session
#'
#' SelectionTools is being retired from haploStack. Everything SelectionTools-
#' specific is behind this one switch, which is \strong{off} by default: while
#' it is off, such a call stops with an \code{hs_error_input} naming this
#' function (field \code{hint}). The switch covers
#' \itemize{
#'   \item its \strong{file formats}: \code{\link{st_convert}},
#'     \code{\link{st_ld_blocks}}, \code{\link{st_ld_blocks_files}} and the
#'     deprecated names that call them; the \code{"ST"} preset of
#'     \code{\link{read_genotype}} / \code{\link{read_st_genotype}} (its cells
#'     are SelectionTools' allele codes) and of \code{\link{geno_meta_to_matrix}}
#'     / \code{\link{write_geno_meta}}. Each of these also takes
#'     \code{st_formats = TRUE} for one call.
#'   \item its \strong{allele encoding token}: \code{encoding = "selectiontools"}
#'     on \code{\link{build_geno_meta}}, \code{\link{new_geno_meta}},
#'     \code{\link{sim_geno_meta}} and \code{\link{sim_geno_matrix}}; a
#'     \code{"selectiontools"} side of \code{\link{translate_alleles}} or
#'     \code{\link{transform_effects}}; and
#'     \code{\link{use_missing_chars}(mode = "selectiontools")}. These are
#'     deprecated (see below) and take no per-call argument: the session switch
#'     is the only way in.
#' }
#' The one exception is \code{\link{upgrade_geno_meta}} — so
#' \code{\link{load_geno_meta}} and \code{\link{recode_geno_meta}} on a
#' container saved with the old vocabulary — which recodes without the switch:
#' it is a one-way migration that announces itself, and a 0.8.8 \code{.rds}
#' must open.
#'
#' @section The SelectionTools allele encoding:
#' SelectionTools codes alleles as integers — \code{A = 1, C = 2, G = 3,
#' T = 4}, missing \code{-1} — and writes a diploid genotype as \code{1/-1}.
#' That code table is the real reason the formats are being retired: since
#' 0.9.0 a \code{GenoMeta} stores nucleotide alleles only (\code{A/C/G/T},
#' plus the dosage vocabulary), and the table is applied at the boundary.
#' \code{\link{st_convert}} and the \code{"ST"} preset translate on the way
#' in; \code{st_convert(to = "selectiontools")} recodes on the way out;
#' \code{encoding = "selectiontools"} passed to \code{\link{build_geno_meta}},
#' \code{\link{new_geno_meta}} or the \code{sim_*} generators is deprecated
#' (with the switch on, the alleles are translated and the container is
#' stamped \code{"nucleotide"}); integer codes with no declared encoding are
#' refused with an \code{hs_error_format} naming this page; and
#' \code{\link{upgrade_geno_meta}} (so \code{\link{load_geno_meta}}) recodes a
#' container saved with the old vocabulary.
#'
#' @param enable \code{TRUE} (default) or \code{FALSE}.
#' @return The previous value of the option, invisibly (\code{NULL} when it
#'   was never set), so it can be put back with
#'   \code{options(haploStack.st_formats = )}.
#' @examples
#' old <- use_st_formats(TRUE)
#' options(haploStack.st_formats = old)
#' @seealso \code{\link{st_convert}}
#' @export
use_st_formats <- function(enable = TRUE) {
  enable <- .check_flag(enable, "argument", "st_formats")
  old <- .hs_option_raw("st_formats")
  options(haploStack.st_formats = enable)
  .hs_inform("use_st_formats: SelectionTools formats ",
             if (enable) "enabled" else "disabled", " for this session.")
  invisible(old)
}

# The gate every SelectionTools entry point passes first: the file formats and
# the "selectiontools" encoding token alike (C-6, owner 2026-11-03: the adapter
# is inert unless the switch is on). `per_call = FALSE` for the deprecated
# token on builders, generators, the legacy layer and the registry, which take
# no st_formats argument: the session switch is the only way in.
.st_gate <- function(caller, st_formats = NULL, per_call = TRUE) {
  if (!isTRUE(.hs_setting("st_formats", st_formats)))
    .hs_abort(class = "hs_error_input",
              data = list(hint = "use_st_formats(TRUE)"),
              caller, ": SelectionTools formats are off (SelectionTools is ",
              "being retired from haploStack). Call use_st_formats(TRUE) for ",
              "this session",
              if (per_call) ", or pass st_formats = TRUE to this call." else ".",
              call. = FALSE)
  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# The SelectionTools allele encoding (C-5)
# ---------------------------------------------------------------------------

#' The SelectionTools allele code table
#'
#' SelectionTools codes the nucleotides in alphabetical order, \code{A = 1,
#' C = 2, G = 3, T = 4}, and a missing allele as \code{-1}.
#' \code{DEFAULT_ENCODING_DICT} is that table as a named list (character to
#' integer); \code{REVERSED_ENCODING_DICT} is its inverse (code, as a string,
#' to character). Since 0.9.0 the table is used only at the SelectionTools
#' boundary (\code{\link{st_convert}}, the \code{"ST"} preset, the
#' \code{\link{upgrade_geno_meta}} recode): a \code{GenoMeta} stores
#' nucleotide alleles, never these codes. The legacy
#' \code{\link{translate_alleles}} / \code{\link{transform_effects}} layer
#' still reads it.
#'
#' @format \code{DEFAULT_ENCODING_DICT}: a named list of five integers.
#'   \code{REVERSED_ENCODING_DICT}: a named character vector of five.
#' @seealso \code{\link{use_st_formats}}, \code{\link{EXTENDED_ALLELE_DICT}}
#' @export
DEFAULT_ENCODING_DICT <- list(
  "A" = 1L,
  "C" = 2L,
  "G" = 3L,
  "T" = 4L,
  "-" = -1L
)

#' @rdname DEFAULT_ENCODING_DICT
#' @export
REVERSED_ENCODING_DICT <- setNames(
  names(DEFAULT_ENCODING_DICT),
  as.character(unname(unlist(DEFAULT_ENCODING_DICT)))
)

# The four allele codes, nucleotide -> code (as strings).
.ST_CODES <- c(A = "1", C = "2", G = "3", T = "4")

# The one translation table, both ways (C-6: every recode in the package reads
# these two vectors; nothing else indexes DEFAULT_ENCODING_DICT for a recode).
# SelectionTools' missing code -1 pairs with the nucleotide missing character.
.ST_TO_NUCLEOTIDE   <- c(stats::setNames(names(.ST_CODES), .ST_CODES), "-1" = "-")
.ST_FROM_NUCLEOTIDE <- c(.ST_CODES, "-" = "-1")

# Translate SelectionTools allele tokens to nucleotides, element-wise, keeping
# dim and dimnames: 1/2/3/4 -> A/C/G/T, and SelectionTools' own missing code
# -1 -> "-" (the nucleotide missing character). Every other token — "0", "-",
# a registered character, anything unknown — is kept verbatim: this is a
# change of vocabulary, not a repair, so what the portal normalises (or
# .missing_token_report() flags on an old container) is unchanged.
.st_to_nucleotide <- function(x) .st_translate(x, .ST_TO_NUCLEOTIDE)

# The inverse, for the legacy translate_alleles() / transform_effects() layer:
# A/C/G/T -> 1/2/3/4, "-" -> -1, anything else verbatim.
.st_from_nucleotide <- function(x) .st_translate(x, .ST_FROM_NUCLEOTIDE)

.st_translate <- function(x, lut) {
  v    <- trimws(as.character(x))
  out  <- unname(lut[v])
  keep <- is.na(out) & !is.na(v)
  out[keep] <- as.character(x)[keep]
  attributes(out) <- attributes(x)
  out
}

# Translate a 2-D matrix of genotype cells ("1/4", "1/-1", "14") to nucleotide
# cells ("A/T", "A/-", "AT"); a missing allele becomes "-", which the
# nucleotide vocabulary reads as missing at the portal. Dimnames are kept.
# Whole digit tokens only: a digit that is part of a longer token (say "12"
# in a delimited cell) is not a code and is left alone.
.st_decode_cells <- function(cells, split_char = "/", missing_chars = NULL) {
  cells <- as.matrix(cells)
  v <- as.character(cells)
  if (nzchar(split_char)) {
    v <- gsub("(?<![0-9])-1(?![0-9])", "-", v, perl = TRUE)
    for (tk in .extra_missing_chars("selectiontools", missing_chars))
      v <- gsub(paste0("(?<![^", split_char, "])", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", tk),
                       "(?![^", split_char, "])"), "-", v, perl = TRUE)
    v <- gsub("(?<![0-9])0(?![0-9])", "-", v, perl = TRUE)
    for (nt in names(.ST_CODES))
      v <- gsub(paste0("(?<![0-9])", .ST_CODES[[nt]], "(?![0-9])"), nt, v, perl = TRUE)
  } else {
    v <- chartr("01234", "-ACGT", v)
  }
  attributes(v) <- attributes(cells)
  v
}

# Are these allele tokens SelectionTools' integer codes? TRUE when every
# non-missing token is one of 1/2/3/4 (so a PLINK 1/2-coded array says TRUE
# too, which is why the portal asks for a declaration rather than guessing).
.st_is_coded <- function(x) {
  u <- unique(trimws(as.character(x)))
  u <- u[!is.na(u) & !u %in% c("", "NA", "N", "n", ".", "-", "0", "-1")]
  length(u) > 0L && all(u %in% .ST_CODES)
}

# Integer allele codes reached a container builder with no declared encoding.
.st_refuse_codes <- function(caller) {
  .hs_abort(class = c("hs_error_format", "hs_error_input"),
            data = list(hint = "st_convert()"),
            caller, ": the alleles are integer codes (1/2/3/4) with no declared ",
            "encoding. haploStack 0.9.0 stores nucleotide alleles only. If they ",
            "are SelectionTools codes (A = 1, C = 2, G = 3, T = 4, missing -1), ",
            "read the file with st_convert() after use_st_formats(TRUE), or pass ",
            "encoding = \"selectiontools\" to have them translated (deprecated). ",
            "If the digits are allele labels of your own (PLINK's 1/2 coding, ",
            "for instance), pass encoding = \"nucleotide\" to keep them as labels.",
            call. = FALSE)
}

# The deprecation warning for the "selectiontools" encoding token (ladder
# rung 1, F7): `fn(what)` names the call, `then` says what happens instead.
.st_encoding_deprecated <- function(fn, what = "encoding = \"selectiontools\"",
                                    then = paste0("the alleles are translated to ",
                                                  "nucleotides (A/C/G/T)")) {
  .hs_warn(sprintf(paste0("%s(%s) is deprecated since haploStack 0.9.0: a GenoMeta ",
                          "stores nucleotide alleles only, so %s. Use st_convert() ",
                          "for SelectionTools' own formats.\n",
                          "See help(\"haploStack-deprecated\")."),
                   fn, what, then),
           call. = FALSE,
           class = c("hs_warning_deprecated", "deprecatedWarning"),
           data  = list(old = paste0(fn, "(", what, ")"), new = "st_convert()",
                        since = "0.9.0"))
}

# The char_dict of a translated array: the observed nucleotide labels, the
# missing tokens left out (the translation runs before the portal normalises,
# so "0" and "-" may still be in the array).
.st_char_dict <- function(arr) {
  obs <- .hs_sort(unique(as.vector(arr)))
  obs <- obs[!is.na(obs) & nzchar(obs) & !obs %in% .missing_tokens("nucleotide")]
  stats::setNames(as.list(seq_along(obs)), obs)
}

# The portal's half (new_geno_meta()): a declared "selectiontools" encoding
# is deprecated — the array is translated, the container is stamped
# "nucleotide" and char_dict / missing_tokens are re-derived; integer codes
# with no declared encoding are refused.
.st_portal <- function(geno_arrays, metadata, caller = "new_geno_meta") {
  enc <- metadata$encoding
  if (identical(enc, "selectiontools")) {
    .st_gate(caller, per_call = FALSE)               # C-6: the token is gated too
    .st_encoding_deprecated(caller)
    geno_arrays             <- .st_to_nucleotide(geno_arrays)
    metadata$encoding       <- "nucleotide"
    metadata$char_dict      <- .st_char_dict(geno_arrays)
    metadata$missing_tokens <- NULL
  } else if (is.null(enc) && is.character(geno_arrays) && .st_is_coded(geno_arrays)) {
    .st_refuse_codes(caller)
  }
  list(geno_arrays = geno_arrays, metadata = metadata)
}

# upgrade_geno_meta()'s half: a container saved with metadata$encoding
# "selectiontools" (0.8.8 and earlier, or a 0.9.0 development build) is
# recoded in place — the array, char_dict and the Allele column of any stored
# marker effects — and stamped "nucleotide". Runs whatever the schema says.
.st_upgrade_container <- function(x, verbose = TRUE) {
  if (!identical(x$metadata$encoding, "selectiontools")) return(x)
  x$geno_arrays <- .st_to_nucleotide(x$geno_arrays)
  md <- x$metadata
  md$encoding  <- "nucleotide"
  md$char_dict <- .st_char_dict(x$geno_arrays)
  x$metadata   <- md
  if (is.list(x$marker_effects) && length(x$marker_effects))
    x$marker_effects <- lapply(x$marker_effects, function(d) {
      if (is.data.frame(d) && "Allele" %in% names(d))
        d$Allele <- .st_to_nucleotide(as.character(d$Allele))
      d
    })
  if (isTRUE(verbose))
    .hs_inform("upgrade_geno_meta: SelectionTools allele codes (1/2/3/4) recoded ",
               "to nucleotides (A/C/G/T) — haploStack 0.9.0 stores nucleotide ",
               "alleles only.")
  x
}


#' Convert between SelectionTools' formats and haploStack's
#'
#' The one place SelectionTools' own file formats are read or written. It is
#' off until \code{\link{use_st_formats}} is called (or \code{st_formats =
#' TRUE} is passed), because SelectionTools is being retired. Three
#' conversions, chosen by what \code{x} is:
#' \describe{
#'   \item{SelectionTools genotype + map files in}{\code{x} = path to a
#'     genotype file in SelectionTools' matrix format (\code{"m"}: header
#'     \code{m/i <ind> ...}, one line per marker, fields separated by a tab, a
#'     blank or \code{;}); \code{map} = its map (\code{"mcp"}: \code{name chrom
#'     pos} with a header line; a path, a data.frame or a \code{GenoMeta}).
#'     Returns a \code{GenoMeta} built by \code{\link{build_geno_meta}}: the
#'     allele codes are translated to nucleotides (\code{1/2/3/4} to
#'     \code{A/C/G/T}; \code{-1}, \code{0} and \code{-} are missing), the
#'     container is stamped \code{encoding = "nucleotide"} with
#'     \code{source_format = "selectiontools"} as its provenance, and markers
#'     absent from the map are dropped as in any read.}
#'   \item{SelectionTools block table in}{\code{x} = the table
#'     \code{st.def.hblocks()} returns (a data.frame, or the path of one written
#'     with \code{write.table}), recognised by its \code{Name} and
#'     \code{Markers} columns; \code{map} = the panel's map (a \code{GenoMeta},
#'     its \code{map_info}, or the path of the \code{"mcp"} file). Returns the
#'     block-marker map \code{(Chrom, Marker, Block)} — \code{Chrom} taken from
#'     the map, since SelectionTools keeps only the leading integer of a
#'     chromosome name — ready for \code{\link{attach_block_info}}.}
#'   \item{GenoMeta out}{\code{x} = a \code{GenoMeta}, \code{to =
#'     "selectiontools"}. Writes \code{<prefix>_geno.txt} (matrix \code{"m"},
#'     tab-separated, missing \code{-}; with \code{dict = "marker"}, the
#'     default, each marker's alleles are recoded to \code{1, 2, ...} in
#'     sorted order, which is all SelectionTools' LD needs; with \code{dict =
#'     "fixed"} the fixed table \code{A = 1, C = 2, G = 3, T = 4} is written,
#'     and a marker carrying any other allele is an error) and \code{<prefix>_map.txt}
#'     (\code{"mcp"}; chromosomes as integer keys in natural order, positions
#'     in bp) into \code{output_dir}. Markers without a usable position, and
#'     markers with more than five allele states, are left out (SelectionTools
#'     orders by position and codes alleles 1..5). Returns the two paths and
#'     the chromosome key table, invisibly. \code{\link{st_ld_blocks}} writes
#'     the same files for its own run, under aliases.}
#' }
#'
#' @param x A file path, a data.frame or a \code{GenoMeta}; see above.
#' @param map The map that goes with \code{x} (files in): a path, a data.frame
#'   (\code{Chrom, Marker, Pos}, or SelectionTools' \code{name, chrom, pos}) or
#'   a \code{GenoMeta}. Ignored for \code{to = "selectiontools"}.
#' @param to \code{"haploStack"} (default): read SelectionTools formats in.
#'   \code{"selectiontools"}: write a \code{GenoMeta} out.
#' @param output_dir Directory for \code{to = "selectiontools"} (created).
#' @param prefix File-name prefix for \code{to = "selectiontools"}. Default
#'   \code{"hs"}.
#' @param dict Allele codes for \code{to = "selectiontools"}: \code{"marker"}
#'   (default; \code{1, 2, ...} per marker) or \code{"fixed"} (\code{A = 1,
#'   C = 2, G = 3, T = 4}); see above.
#' @param missing_chars Extra missing-allele characters for a genotype read:
#'   characters that mean missing in the SelectionTools file, beyond
#'   \code{-1}, \code{0} and \code{-} (see \code{\link{build_geno_meta}}).
#' @param st_formats \code{TRUE} enables SelectionTools formats for this call;
#'   \code{NULL} (default) takes the session switch
#'   (\code{\link{use_st_formats}}).
#' @return A \code{GenoMeta}, a block-marker map data.frame, or (invisibly) a
#'   list \code{geno_file, map_file, chrom_keys}; see above.
#' @examples
#' \dontrun{
#' use_st_formats(TRUE)
#' meta <- st_convert("geno.txt", map = "map.txt")               # files in
#' st_convert(meta, to = "selectiontools", output_dir = "st")    # files out
#' blk  <- read.table("blocks.ld.txt", header = TRUE, stringsAsFactors = FALSE)
#' meta <- attach_block_info(meta, st_convert(blk, map = meta))  # blocks in
#' }
#' @seealso \code{\link{use_st_formats}}, \code{\link{st_ld_blocks}},
#'   \code{\link{attach_block_info}}, \code{\link{build_geno_meta}}
#' @export
st_convert <- function(x, map = NULL, to = c("haploStack", "selectiontools"),
                       output_dir = NULL, prefix = "hs",
                       dict = c("marker", "fixed"), missing_chars = NULL,
                       st_formats = NULL) {
  caller <- "st_convert"
  .st_gate(caller, st_formats)
  to   <- match.arg(to)
  dict <- match.arg(dict)

  # ---- GenoMeta out ----------------------------------------------------------
  if (to == "selectiontools") {
    .hs_check_meta(x, caller, "x")
    if (is.null(output_dir))
      .hs_abort(caller, ": 'output_dir' is needed to write SelectionTools files.",
                call. = FALSE)
    if (dim(x$geno_arrays)[3L] != 2L)
      .hs_abort(caller, ": SelectionTools reads diploid genotypes; this ",
                "container has ploidy ", dim(x$geno_arrays)[3L], ".", call. = FALSE)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    placed <- .ld_placed_markers(x$map_info)
    if (!length(placed))
      .hs_abort(caller, ": no marker carries a usable position.", call. = FALSE)
    io <- .st_write_inputs(x, placed, output_dir, alias = FALSE, prefix = prefix,
                           caller = caller, dict = dict)
    .hs_inform(caller, ": wrote ", io$geno_file, " and ", io$map_file, " (",
               nrow(io$markers), " marker(s), ", dim(x$geno_arrays)[1L],
               " individual(s)).")
    return(invisible(list(geno_file = io$geno_file, map_file = io$map_file,
                          chrom_keys = io$chrom_keys)))
  }

  # ---- files in --------------------------------------------------------------
  if (is.list(x) && !is.data.frame(x) && !is.null(x$geno_arrays))
    .hs_abort(caller, ": 'x' is already a GenoMeta; to = \"selectiontools\" ",
              "writes it out.", call. = FALSE)
  tab <- NULL
  if (is.data.frame(x)) {
    tab <- as.data.frame(x, stringsAsFactors = FALSE)
    if (!all(c("Name", "Markers") %in% names(tab)))
      .hs_abort(caller, ": a data.frame 'x' must be SelectionTools' block table ",
                "(columns Name and Markers); got ", paste(names(tab), collapse = ", "),
                ". A genotype comes in as a file path.", call. = FALSE)
  } else if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x))
      .hs_abort(caller, ": file not found: ", x, call. = FALSE)
    peek <- .read_table_chr(x, header = TRUE)
    if (all(c("Name", "Markers") %in% names(peek))) tab <- peek
  } else {
    .hs_abort(caller, ": 'x' must be a file path, a data.frame or a GenoMeta.",
              call. = FALSE)
  }

  if (!is.null(tab)) {
    if (is.null(map))
      .hs_abort(caller, ": a SelectionTools block table needs 'map' (the panel's ",
                "map, or its GenoMeta) for the chromosome names.", call. = FALSE)
    mp <- .st_as_map(map, caller)
    alias_map <- data.frame(alias = mp$Marker, Marker = mp$Marker,
                            Chrom = mp$Chrom, stringsAsFactors = FALSE)
    return(.st_parse_blocks(tab, alias_map, caller = caller))
  }

  if (is.null(map))
    .hs_abort(caller, ": a SelectionTools genotype file needs 'map' (its \"mcp\" ",
              "map file).", call. = FALSE)
  io <- .st_read_inputs(x, caller)
  mp <- .st_as_map(map, caller)
  # The code table is applied here, at the boundary: the builder sees
  # nucleotides and the container never holds a SelectionTools code.
  cells <- .st_decode_cells(io$geno, split_char = "/", missing_chars = missing_chars)
  build_geno_meta(cells, mp, split_char = "/", encoding = "nucleotide",
                  source_format = "selectiontools", missing_chars = missing_chars)
}


# ---------------------------------------------------------------------------
# The three formats
# ---------------------------------------------------------------------------

# A map from a path (SelectionTools "mcp": name chrom pos, header), a
# data.frame (canonical or SelectionTools names, else positional) or a
# GenoMeta -> data.frame(Chrom, Marker, Pos), identifiers as text.
.st_as_map <- function(map, caller) {
  if (is.list(map) && !is.data.frame(map) && !is.null(map$map_info))
    map <- map$map_info
  if (is.character(map) && length(map) == 1L) {
    if (!file.exists(map))
      .hs_abort(caller, ": map file not found: ", map, call. = FALSE)
    m <- .read_table_chr(map, header = TRUE)
    if (ncol(m) < 3L)
      .hs_abort(caller, ": a SelectionTools map file has three columns (name, ",
                "chrom, pos); ", basename(map), " has ", ncol(m), ".", call. = FALSE)
    return(data.frame(Chrom  = as.character(m[[2L]]),
                      Marker = as.character(m[[1L]]),
                      Pos    = suppressWarnings(as.numeric(m[[3L]])),
                      stringsAsFactors = FALSE))
  }
  if (!is.data.frame(map))
    .hs_abort(caller, ": 'map' must be a file path, a data.frame or a GenoMeta.",
              call. = FALSE)
  m  <- as.data.frame(map, stringsAsFactors = FALSE)
  nm <- tolower(names(m))
  pick <- function(canon, st) {
    i <- match(c(canon, st), nm)
    i <- i[!is.na(i)]
    if (length(i)) m[[i[1L]]] else NULL
  }
  chrom  <- pick("chrom", "chr")
  marker <- pick("marker", "name")
  pos    <- pick("pos", "position")
  if (is.null(chrom) || is.null(marker) || is.null(pos)) {
    if (ncol(m) < 3L)
      .hs_abort(caller, ": 'map' needs Chrom, Marker and Pos (or SelectionTools' ",
                "name, chrom, pos).", call. = FALSE)
    marker <- m[[1L]]; chrom <- m[[2L]]; pos <- m[[3L]]   # SelectionTools order
  }
  data.frame(Chrom  = as.character(chrom),
             Marker = as.character(marker),
             Pos    = suppressWarnings(as.numeric(pos)),
             stringsAsFactors = FALSE)
}

# SelectionTools matrix format "m" -> sample x marker matrix of "a/b" cells.
.st_read_inputs <- function(geno_file, caller) {
  g <- .read_table_chr(geno_file, header = TRUE)
  if (ncol(g) < 2L || nrow(g) < 1L)
    .hs_abort(caller, ": ", basename(geno_file), " does not look like ",
              "SelectionTools' matrix format (header 'm/i <ind> ...', one line ",
              "per marker).", call. = FALSE)
  markers <- as.character(g[[1L]])
  samples <- names(g)[-1L]
  cells   <- t(as.matrix(g[, -1L, drop = FALSE]))
  dimnames(cells) <- list(samples, markers)
  list(geno = cells, markers = markers, samples = samples)
}

# Write the two SelectionTools input files for `markers` (placed, any order)
# into `dir`. Returns the files, the marker table and the chromosome keys.
#
#   <prefix>_geno.txt  matrix format "m": header "m/i  i1 i2 ...", one line per
#                      marker "m7  1/2 2/2 -/1 ...", tab-separated. Alleles are
#                      recoded per marker to 1, 2, ... in sorted order — r^2
#                      does not depend on which allele is called what, and
#                      SelectionTools reads 1..5 whatever the source vocabulary
#                      was (nucleotides, 1/2/3/4 codes, DArT A/B). Missing = "-".
#   <prefix>_map.txt   map format "mcp" with a header: "name chrom pos", chrom
#                      as an integer key in natural order, pos in bp.
#
# alias = TRUE (the resolver) writes markers as m1, m2, ... and individuals as
# i1, i2, ..., so that marker and sample names — spaces, "|", very long IDs —
# are immaterial to SelectionTools' tokenizer; the `markers` table translates
# back. alias = FALSE (st_convert()) writes the real names.
.st_write_inputs <- function(geno_meta, markers, dir, alias = TRUE,
                             prefix = "hs", caller = "st_ld_blocks",
                             dict = c("marker", "fixed")) {
  dict <- match.arg(dict)
  map <- as.data.frame(geno_meta$map_info, stringsAsFactors = FALSE)
  map <- map[match(markers, map$Marker), c("Chrom", "Marker", "Pos"), drop = FALSE]
  map$Chrom <- as.character(map$Chrom)
  map$Pos   <- as.numeric(map$Pos)

  chrom_levels <- .natural_chrom_levels(map$Chrom)
  map$key <- match(map$Chrom, chrom_levels)
  map <- map[.hs_order(map$key, map$Pos, seq_len(nrow(map))), , drop = FALSE]

  arr <- geno_meta$geno_arrays[, map$Marker, , drop = FALSE]
  a1  <- matrix(as.character(arr[, , 1L]), nrow = dim(arr)[1L])
  a2  <- matrix(as.character(arr[, , 2L]), nrow = dim(arr)[1L])
  c1  <- matrix(NA_integer_, nrow(a1), ncol(a1))
  c2  <- c1
  n_all <- integer(ncol(a1))
  for (j in seq_len(ncol(a1))) {
    lv <- .hs_sort(unique(c(a1[, j], a2[, j])))
    lv <- lv[!is.na(lv)]
    n_all[j] <- length(lv)
    if (dict == "fixed") {
      bad <- setdiff(lv, names(.ST_CODES))
      if (length(bad))
        .hs_abort(caller, ": dict = \"fixed\" writes A = 1, C = 2, G = 3, T = 4 ",
                  "only; marker ", map$Marker[j], " carries ",
                  paste0("\"", bad, "\"", collapse = ", "), ".", call. = FALSE)
      c1[, j] <- as.integer(.ST_CODES[a1[, j]])
      c2[, j] <- as.integer(.ST_CODES[a2[, j]])
    } else {
      c1[, j] <- match(a1[, j], lv)
      c2[, j] <- match(a2[, j], lv)
    }
  }
  # SelectionTools codes alleles 1..5; a marker with more states cannot be
  # written faithfully and is left out (it becomes a singleton on attach).
  too_many <- n_all > 5L
  if (any(too_many)) {
    .hs_warn(caller, ": ", sum(too_many), " marker(s) carry ",
             "more than 5 allele states and are left out of SelectionTools' ",
             "files (e.g. ", paste(utils::head(map$Marker[too_many], 3L),
             collapse = ", "), ").", call. = FALSE)
    keep <- !too_many
    map <- map[keep, , drop = FALSE]
    c1 <- c1[, keep, drop = FALSE]
    c2 <- c2[, keep, drop = FALSE]
  }

  m_alias <- if (isTRUE(alias)) paste0("m", seq_len(nrow(map))) else map$Marker
  i_alias <- if (isTRUE(alias)) paste0("i", seq_len(nrow(c1)))
             else dimnames(geno_meta$geno_arrays)[[1L]]
  cell <- function(x) { x <- as.character(x); x[is.na(x)] <- "-"; x }
  cells <- matrix(paste0(cell(c1), "/", cell(c2)), nrow = nrow(c1))

  geno_file <- file.path(dir, paste0(prefix, "_geno.txt"))
  map_file  <- file.path(dir, paste0(prefix, "_map.txt"))
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

  list(geno_file  = geno_file, map_file = map_file,
       markers    = data.frame(alias = m_alias, Marker = map$Marker,
                               Chrom = map$Chrom, Pos = map$Pos,
                               stringsAsFactors = FALSE),
       chrom_keys = data.frame(Chrom = chrom_levels, key = seq_along(chrom_levels),
                               stringsAsFactors = FALSE),
       n_left_out = sum(too_many))
}

# st.def.hblocks() table -> (Chrom, Marker, Block) in the container's marker
# IDs and chromosome names. `alias_map` has alias, Marker, Chrom (alias ==
# Marker when the table was built on real names).
.st_parse_blocks <- function(raw, alias_map, caller = "st_ld_blocks") {
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  if (!all(c("Name", "Markers") %in% names(raw)))
    .hs_abort(class = if (identical(caller, "st_convert")) "hs_error_input"
                      else "hs_error_internal",
              caller, ": unexpected SelectionTools block table ",
              "(columns: ", paste(names(raw), collapse = ", "),
              "; expected Name and Markers).", call. = FALSE)
  members <- strsplit(gsub("[[:space:]]", "", as.character(raw$Markers)), ";",
                      fixed = TRUE)
  members <- lapply(members, function(x) x[nzchar(x)])
  bm <- data.frame(alias = unlist(members, use.names = FALSE),
                   Block = rep(as.character(raw$Name), lengths(members)),
                   stringsAsFactors = FALSE)
  idx <- match(bm$alias, alias_map$alias)
  if (anyNA(idx))
    .hs_abort(caller, ": the SelectionTools block table names marker(s) ",
              if (identical(caller, "st_convert")) "not in 'map'"
              else "the bridge did not write",
              " (", paste(utils::head(bm$alias[is.na(idx)], 3L), collapse = ", "),
              ").", call. = FALSE)
  out <- data.frame(Chrom  = alias_map$Chrom[idx],
                    Marker = alias_map$Marker[idx],
                    Block  = bm$Block,
                    stringsAsFactors = FALSE)
  out <- out[order(match(out$Marker, alias_map$Marker)), , drop = FALSE]
  rownames(out) <- NULL
  out
}
