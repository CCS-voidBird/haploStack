# =============================================================================
# test-090-hygiene.R — process-state hygiene of the source tree (design D1-3)
#
# A haploStack call must leave the R process as it found it (framework F1/F2).
# These tests grep R/ for the constructs that touch process state and allow
# each only in the file that owns it. Each allowlist names the build step that
# removes the remaining sites; shrinking a list is the only edit expected here.
#
#   H1  getOption() only in R/state.R (.hs_setting / .hs_options_get) and R/zzz.R
#   H2  set.seed() only in R/state.R
#   H3  setwd() nowhere (D1-8 removed the SelectionTools site, F-5)
#   H4  assign(..., globalenv()) only in the SelectionTools io-dir shim (D1-8)
#   H5  future::plan() only in R/state.R (.hs_with_plan / .hs_keep_plan, D1-4)
#
# Skipped when R/ is not reachable from the test directory (installed package).
# =============================================================================

.h090_r_dir <- function() {
  d <- testthat::test_path("..", "..", "R")
  if (!dir.exists(d)) skip("R/ not reachable from the test directory")
  d
}

# Code lines of every R/ file, comments and roxygen stripped, as a data.frame
# (file, line, text). A line whose first non-blank character is '#' is a
# comment; trailing comments are not stripped (none of the patterns below is
# ever legitimately inside one).
.h090_code_lines <- function() {
  d <- .h090_r_dir()
  files <- sort(list.files(d, pattern = "[.]R$", full.names = TRUE))
  out <- lapply(files, function(f) {
    txt <- readLines(f, warn = FALSE, encoding = "UTF-8")
    keep <- !grepl("^\\s*#", txt) & nzchar(trimws(txt))
    data.frame(file = basename(f), line = which(keep), text = txt[keep],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.h090_sites <- function(lines, pattern) {
  hit <- lines[grepl(pattern, lines$text, perl = TRUE), , drop = FALSE]
  sprintf("%s:%d", hit$file, hit$line)
}

test_that("H1: options are read only through .hs_setting() (state.R) and zzz.R", {
  lines <- .h090_code_lines()
  # no literal haploStack.* option read anywhere but the load hook
  lit <- .h090_sites(lines, 'getOption\\("haploStack\\.')
  expect_true(all(sub(":.*$", "", lit) %in% "zzz.R"),
              label = paste("literal haploStack option reads outside zzz.R:",
                            paste(lit, collapse = ", ")))
  # and getOption() itself only in state.R (resolver + option snapshot) and
  # the load hook (D1-5 moved the progressr containment into state.R)
  any_opt <- .h090_sites(lines, "getOption\\(")
  allowed <- c("state.R", "zzz.R")
  expect_true(all(sub(":.*$", "", any_opt) %in% allowed),
              label = paste("getOption() outside the allowlist:",
                            paste(any_opt, collapse = ", ")))
  expect_true(any(grepl("^state\\.R:", any_opt)))     # the resolver is there
})

test_that("H2: set.seed() lives in state.R only (D4-3 and D4-5 removed the two driver sites)", {
  lines <- .h090_code_lines()
  sites <- .h090_sites(lines, "(^|[^._[:alnum:]])set\\.seed\\(")
  files <- unique(sub(":.*$", "", sites))
  # D4-3 (job boundary, stacking_main.R) and D4-5 (OCS, stacking_selection.R)
  # replace these with .with_rng_state() / .with_seed(); drop them here then.
  allowed <- c("state.R")
  expect_true(all(files %in% allowed),
              label = paste("set.seed() outside the allowlist:",
                            paste(sites[!sub(":.*$", "", sites) %in% allowed],
                                  collapse = ", ")))
})

test_that("H3/H4: working directory and global environment (D1-8: no setwd() anywhere, the ST shim only)", {
  lines <- .h090_code_lines()
  wd <- .h090_sites(lines, "(^|[^._[:alnum:]])setwd\\(")
  expect_identical(wd, character(0L))           # D1-8 (F-5) removed the last site
  ge <- .h090_sites(lines, "assign\\(.*(globalenv\\(\\)|\\.GlobalEnv)")
  expect_true(all(sub(":.*$", "", ge) %in% "blocks_selectiontools.R"),
              label = paste("assign() into the global environment outside the ST shim:",
                            paste(ge, collapse = ", ")))
})

test_that("H5: future::plan() only in state.R (.hs_with_plan / .hs_keep_plan, D1-4)", {
  lines <- .h090_code_lines()
  pl <- .h090_sites(lines, "future::plan\\(")
  allowed <- c("state.R")
  expect_true(all(sub(":.*$", "", pl) %in% allowed),
              label = paste("future::plan() outside the allowlist:",
                            paste(pl, collapse = ", ")))
})
