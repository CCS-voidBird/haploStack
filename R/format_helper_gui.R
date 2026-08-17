# =============================================================================
# format_helper_gui.R
#
# run_format_converter_app() — Shiny UI for GenoMeta format conversion.
#
# Q9 (2026-06-08): Shiny confirmed; lives in Suggests.
# Requires: shiny (Suggests field of DESCRIPTION).
# =============================================================================

#' Launch the interactive genotype format converter
#'
#' Opens a Shiny application that wraps the \code{\link{format_helper}} reader
#' functions, allowing users to upload a genotype file, select input and output
#' formats, optionally attach a block file, preview the result, and download the
#' converted files — all without writing R code.
#'
#' The app covers the three supported format pairs:
#' \itemize{
#'   \item Plink 1.9 PED + MAP (\code{\link{read_plink_ped}})
#'   \item SelectionTools CSV (\code{\link{read_st_genotype}})
#'   \item genomicSimulation CSV (\code{\link{read_gs_genotype}})
#' }
#' Output can be written to any of the same three formats via
#' \code{\link{write_geno_meta}}.
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
    stop(
      "run_format_converter_app() requires the 'shiny' package.\n",
      "Install it with: install.packages(\"shiny\")"
    )

  app <- shiny::shinyApp(
    ui     = .fca_ui(),
    server = .fca_server
  )

  run_args <- list(app, launch.browser = launch_browser, ...)
  if (!is.null(port)) run_args$port <- port

  do.call(shiny::runApp, run_args)
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

.fca_ui <- function() {
  shiny::fluidPage(
    shiny::titlePanel("haploStack — Genotype Format Converter"),

    shiny::sidebarLayout(

      shiny::sidebarPanel(
        width = 3,

        shiny::h4("Input"),
        shiny::selectInput("input_format", "Input format",
          choices = c(
            "Plink PED + MAP" = "plink_ped",
            "SelectionTools CSV" = "st",
            "genomicSimulation CSV" = "gs"
          )
        ),

        # PED needs two file uploads; ST/GS only one
        shiny::conditionalPanel(
          condition = "input.input_format == 'plink_ped'",
          shiny::fileInput("ped_file", "PED file (.ped)", accept = ".ped"),
          shiny::fileInput("map_file", "MAP file (.map)", accept = ".map")
        ),
        shiny::conditionalPanel(
          condition = "input.input_format != 'plink_ped'",
          shiny::fileInput("geno_file", "Genotype file (.csv / .txt)")
        ),

        shiny::fileInput("map_file_extra", "Map file (optional, CSV/tab)",
                         accept = c(".csv",".txt",".tsv")),

        shiny::selectInput("pos_unit", "Position unit in map",
          choices = c("bp", "mbp"), selected = "bp"),

        shiny::fileInput("block_file", "Block file (optional)",
                         accept = c(".csv",".txt",".tsv")),

        shiny::hr(),
        shiny::h4("Output"),
        shiny::selectInput("output_format", "Output format",
          choices = c(
            "Plink PED + MAP" = "plink_ped",
            "SelectionTools CSV" = "st",
            "genomicSimulation CSV" = "gs"
          ), selected = "st"
        ),
        shiny::textInput("out_prefix", "Output file prefix", value = "converted"),

        shiny::hr(),
        shiny::actionButton("btn_run", "Convert", class = "btn-primary"),
        shiny::br(), shiny::br(),
        shiny::downloadButton("dl_files", "Download result")
      ),

      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(

          shiny::tabPanel("Summary",
            shiny::br(),
            shiny::verbatimTextOutput("summary_text")
          ),

          shiny::tabPanel("Map preview",
            shiny::br(),
            shiny::tableOutput("map_preview")
          ),

          shiny::tabPanel("Genotype preview",
            shiny::br(),
            shiny::p(shiny::em("First 10 rows × 10 columns of the converted matrix.")),
            shiny::tableOutput("geno_preview")
          ),

          shiny::tabPanel("Block info",
            shiny::br(),
            shiny::tableOutput("block_info_preview")
          ),

          shiny::tabPanel("Log",
            shiny::br(),
            shiny::verbatimTextOutput("log_text")
          )
        )
      )
    )
  )
}


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

.fca_server <- function(input, output, session) {

  # Reactive: read GenoMeta from uploaded files
  geno_meta_r <- shiny::eventReactive(input$btn_run, {
    log_msgs <- character(0L)

    tryCatch({

      map_path_extra <- if (!is.null(input$map_file_extra))
        input$map_file_extra$datapath else NULL
      block_path     <- if (!is.null(input$block_file))
        input$block_file$datapath else NULL

      meta <- if (input$input_format == "plink_ped") {
        shiny::req(input$ped_file, input$map_file)
        log_msgs <<- c(log_msgs, "Reading Plink PED + MAP...")
        read_plink_ped(
          ped_file   = input$ped_file$datapath,
          map_file   = input$map_file$datapath,
          block_file = block_path,
          pos_unit   = input$pos_unit
        )
      } else if (input$input_format == "st") {
        shiny::req(input$geno_file)
        log_msgs <<- c(log_msgs, "Reading SelectionTools genotype file...")
        read_st_genotype(
          file       = input$geno_file$datapath,
          map_file   = map_path_extra,
          block_file = block_path,
          pos_unit   = input$pos_unit
        )
      } else {
        shiny::req(input$geno_file)
        log_msgs <<- c(log_msgs, "Reading genomicSimulation genotype file...")
        read_gs_genotype(
          file       = input$geno_file$datapath,
          map_file   = map_path_extra,
          block_file = block_path,
          pos_unit   = input$pos_unit
        )
      }

      log_msgs <<- c(log_msgs, paste0(
        "Read complete: ",
        meta$metadata$n_markers, " markers × ",
        meta$metadata$n_samples, " samples (ploidy = ",
        meta$metadata$ploidy, ")"
      ))

      list(meta = meta, log = log_msgs)

    }, warning = function(w) {
      log_msgs <<- c(log_msgs, paste0("[WARNING] ", conditionMessage(w)))
      list(meta = NULL, log = log_msgs)
    }, error = function(e) {
      log_msgs <<- c(log_msgs, paste0("[ERROR] ", conditionMessage(e)))
      list(meta = NULL, log = log_msgs)
    })
  })


  # Reactive: convert GenoMeta to output matrix
  converted_r <- shiny::reactive({
    res <- geno_meta_r()
    shiny::req(!is.null(res$meta))
    geno_meta_to_matrix(res$meta, output_format = input$output_format)
  })


  # Summary tab
  output$summary_text <- shiny::renderText({
    res <- geno_meta_r()
    if (is.null(res$meta))
      return("No data loaded. Check the Log tab for details.")
    m <- res$meta$metadata
    paste0(
      "Source format:  ", m$source_format, "\n",
      "Markers:        ", m$n_markers, "\n",
      "Samples:        ", m$n_samples, "\n",
      "Ploidy:         ", m$ploidy, "\n",
      "Block info:     ",
        if (!is.null(res$meta$block_info)) nrow(res$meta$block_info) else 0,
        " blocks\n",
      "Block-marker:   ",
        if (!is.null(res$meta$block_marker_map_info))
          nrow(res$meta$block_marker_map_info) else 0,
        " rows\n",
      "\nOutput format:  ", input$output_format
    )
  })


  # Map preview tab
  output$map_preview <- shiny::renderTable({
    res <- geno_meta_r()
    shiny::req(!is.null(res$meta))
    head(res$meta$map_info, 20L)
  })


  # Genotype preview tab
  output$geno_preview <- shiny::renderTable({
    cv <- converted_r()
    mat <- cv$matrix
    as.data.frame(mat[seq_len(min(10L, nrow(mat))),
                      seq_len(min(10L, ncol(mat))), drop = FALSE],
                  check.names = FALSE)
  })


  # Block info preview tab
  output$block_info_preview <- shiny::renderTable({
    res <- geno_meta_r()
    shiny::req(!is.null(res$meta))
    bi <- res$meta$block_info
    if (is.null(bi))
      return(data.frame(message = "No block information loaded."))
    head(bi, 20L)
  })


  # Log tab
  output$log_text <- shiny::renderText({
    res <- geno_meta_r()
    paste(res$log, collapse = "\n")
  })


  # Download handler — write converted file(s) to a temp dir and zip
  output$dl_files <- shiny::downloadHandler(
    filename = function() paste0(input$out_prefix, "_converted.zip"),
    content  = function(file) {
      res <- geno_meta_r()
      shiny::req(!is.null(res$meta))

      tmp_dir <- tempfile()
      dir.create(tmp_dir)
      out_path <- file.path(tmp_dir, input$out_prefix)

      write_geno_meta(res$meta, path = out_path,
                      output_format = input$output_format)

      # Collect all written files
      written <- list.files(tmp_dir, full.names = TRUE)

      # Write block files too if available
      if (!is.null(res$meta$block_info)) {
        bi_path <- file.path(tmp_dir, paste0(input$out_prefix, "_block_info.csv"))
        utils::write.csv(res$meta$block_info, bi_path, row.names = FALSE)
        written <- c(written, bi_path)
      }
      if (!is.null(res$meta$block_marker_map_info)) {
        bm_path <- file.path(tmp_dir, paste0(input$out_prefix, "_block_marker_map.csv"))
        utils::write.csv(res$meta$block_marker_map_info, bm_path, row.names = FALSE)
        written <- c(written, bm_path)
      }

      zip::zip(file, files = unique(written), mode = "cherry-pick")
    }
  )
}
