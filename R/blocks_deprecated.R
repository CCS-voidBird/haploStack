# =============================================================================
# blocks_deprecated.R — the pre-0.7.10 names of the SelectionTools resolver
#
# 0.7.10 named every local LD-block function after its backend (hapS_ /
# st_), so that none reads as haploStack's way of getting blocks — that is a
# block map the user already has (R/blocks_ld_common.R). The old names still
# work, with a deprecation warning, and forward their arguments unchanged.
# =============================================================================

#' Deprecated LD-block function names
#'
#' Renamed in 0.7.10 so that each local LD-block function names its backend:
#' \tabular{ll}{
#'   \strong{Old name (≤ 0.7.9)} \tab \strong{Use instead} \cr
#'   \code{define_haplotype_blocks()} \tab \code{\link{st_ld_blocks_files}()} \cr
#'   \code{define_haplotype_blocks_meta()} \tab \code{\link{st_ld_blocks}()} \cr
#' }
#' Both old names ran SelectionTools, which is now the experimental resolver;
#' the default one is \code{\link{hapS_ld_blocks}} (HapSelect). Neither is
#' needed when you already have a block map — give it to
#' \code{\link{attach_block_info}}, or to the reader. Each old name warns once
#' per call and returns exactly what its replacement returns.
#'
#' @param geno_data_file,geno_map_file,ld_threshold,tolerance,output_dir See
#'   \code{\link{st_ld_blocks_files}}.
#' @param geno_meta,attach,prefix,keep_files,verbose See
#'   \code{\link{st_ld_blocks}}.
#'
#' @return What the replacement returns.
#'
#' @name haploStack-deprecated
#' @keywords internal
NULL

#' @rdname haploStack-deprecated
#' @export
define_haplotype_blocks <- function(geno_data_file,
                                    geno_map_file,
                                    ld_threshold = 0.5,
                                    tolerance    = 0,
                                    output_dir   = "haplotype") {
  .Deprecated("st_ld_blocks_files", package = "haploStack", msg = paste0(
    "define_haplotype_blocks() is deprecated since haploStack 0.7.10: it is ",
    "the SelectionTools file route, now st_ld_blocks_files(). With a block ",
    "map you already have, no LD blocking is needed: attach_block_info()."))
  st_ld_blocks_files(geno_data_file, geno_map_file,
                     ld_threshold = ld_threshold, tolerance = tolerance,
                     output_dir = output_dir)
}

#' @rdname haploStack-deprecated
#' @export
define_haplotype_blocks_meta <- function(geno_meta,
                                         ld_threshold = 0.5,
                                         tolerance    = 0L,
                                         attach       = TRUE,
                                         output_dir   = NULL,
                                         prefix       = "st_blocks",
                                         keep_files   = FALSE,
                                         verbose      = TRUE) {
  .Deprecated("st_ld_blocks", package = "haploStack", msg = paste0(
    "define_haplotype_blocks_meta() is deprecated since haploStack 0.7.10: it ",
    "is the SelectionTools resolver, now st_ld_blocks(). The default resolver ",
    "is hapS_ld_blocks() (HapSelect); with a block map you already have, no LD ",
    "blocking is needed: attach_block_info()."))
  st_ld_blocks(geno_meta, ld_threshold = ld_threshold, tolerance = tolerance,
               attach = attach, output_dir = output_dir, prefix = prefix,
               keep_files = keep_files, verbose = verbose)
}
