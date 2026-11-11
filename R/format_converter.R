# =============================================================================
# format_converter.R — launcher for the Shiny genotype format converter
#
# The app itself (UI + server) lives in inst/shiny/format_converter/app.R since
# 0.9.0 (decision-form row 28 a): it is shiny-only code and ships as an app
# directory, so it can be run with shiny::runApp() on its own and the package
# namespace carries no Shiny closures. This file is the one exported entry.
# =============================================================================

#' Launch the interactive genotype format converter
#'
#' Opens a Shiny application that wraps the readers, allowing users to upload
#' a genotype file, select input and output formats, optionally attach a block
#' file, preview the result, and download the converted files, all without
#' writing R code.
#'
#' The app covers the three supported format pairs:
#' \itemize{
#'   \item Plink 1.9 PED + MAP (\code{\link{read_plink_ped}})
#'   \item SelectionTools layout (\code{\link{read_st_genotype}}; the app
#'     enables SelectionTools' formats for its own calls, see
#'     \code{\link{use_st_formats}})
#'   \item genomicSimulation CSV (\code{\link{read_gs_genotype}})
#' }
#' Output is written by \code{\link{write_geno_meta}} with the preset the user
#' picks (\code{"PL"}, \code{"ST"} or \code{"GS"}). A warning raised while
#' reading goes to the app's Log tab and the read completes (up to 0.8.8 the
#' first warning dropped the read).
#'
#' The app directory is \code{system.file("shiny", "format_converter", package
#' = "haploStack")}; \code{shiny::runApp()} on that path is the same thing.
#'
#' @param launch_browser Logical. Pass to \code{shiny::runApp()}. Default
#'   \code{TRUE} (opens a browser tab).
#' @param port Integer. Port for the local server. \code{NULL} uses Shiny's
#'   default auto-select. Default \code{NULL}.
#' @param ... Additional arguments forwarded to \code{shiny::runApp()}.
#'
#' @return Invisibly \code{NULL}. Function blocks until the app is closed.
#'
#' @examples
#' \dontrun{
#' run_format_converter_app()
#' }
#'
#' @export
run_format_converter_app <- function(launch_browser = TRUE, port = NULL, ...) {
  if (!requireNamespace("shiny", quietly = TRUE))
    .hs_abort(class = "hs_error_dependency",
      data = list(hint = "install.packages(\"shiny\")"),
      "run_format_converter_app() requires the 'shiny' package.\n",
      "Install it with: install.packages(\"shiny\")"
    )
  app_dir <- .hs_converter_dir()
  run_args <- list(app_dir, launch.browser = launch_browser, ...)
  if (!is.null(port)) run_args$port <- port
  do.call(shiny::runApp, run_args)
  invisible(NULL)
}

# The app directory: installed under inst/shiny (system.file resolves inst/
# under pkgload::load_all too). Missing = a broken install, not an input.
.hs_converter_dir <- function() {
  d <- system.file("shiny", "format_converter", package = "haploStack")
  if (!nzchar(d) || !file.exists(file.path(d, "app.R")))
    .hs_abort(class = "hs_error_state",
              "run_format_converter_app: the app directory ",
              "'shiny/format_converter' is missing from the installed package.",
              call. = FALSE)
  d
}
