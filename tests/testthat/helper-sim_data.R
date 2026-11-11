# =============================================================================
# helper-sim_data.R
#
# Attach data.table so its functions (is.data.table, as.data.table, etc.)
# are available directly in all test files without library(data.table) per
# file.
#
# Up to 0.8.8 this file also held private copies of the six simulated-data
# generators (sim_geno_matrix(), sim_map_data(), sim_block_map(),
# sim_effect_df(), sim_geno_meta(), sim_dataset()). testthat sources helpers
# into the package environment, so every test that called a generator ran the
# copy, not the export. The copies are gone (0.9.0): tests use the package's
# own generators, which draw under a local seed (R/state.R).
# =============================================================================

library(data.table)


# Parallel workers (future multisession) load the INSTALLED haploStack, not a
# namespace attached by pkgload::load_all(); a test that runs package code in a
# worker can only be meaningful from an installed build (phase-gate run from a
# throwaway library). Skip it under a dev namespace (0.9.0, phase E).
.skip_if_dev_namespace <- function() {
  if (requireNamespace("pkgload", quietly = TRUE) &&
      isTRUE(pkgload::is_dev_package("haploStack")))
    testthat::skip("parallel workers need an installed haploStack (dev namespace loaded)")
  invisible(TRUE)
}
