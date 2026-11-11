# haploStack — genotype format converter (Shiny app directory)
#
# Launched by haploStack::run_format_converter_app(), or directly with
# shiny::runApp(system.file("shiny", "format_converter", package = "haploStack")).
# Moved out of the package namespace in 0.9.0 (decision-form row 28 a); it only
# calls exported haploStack functions. SelectionTools' layout is enabled per
# call here (st_formats = TRUE), so the user's session switch is not touched.

ui <- shiny::fluidPage(
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
                       accept = c(".csv", ".txt", ".tsv")),

      shiny::selectInput("pos_unit", "Position unit in map",
        choices = c("bp", "mbp"), selected = "bp"),

      shiny::fileInput("block_file", "Block file (optional)",
                       accept = c(".csv", ".txt", ".tsv")),

      shiny::hr(),
      shiny::h4("Output"),
      shiny::selectInput("output_format", "Output format",
        choices = c(
          "Plink PED + MAP" = "PL",
          "SelectionTools CSV" = "ST",
          "genomicSimulation CSV" = "GS"
        ), selected = "PL"
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

server <- function(input, output, session) {

  # Reactive: read GenoMeta from uploaded files. A warning is logged and the
  # read completes (0.9.0; up to 0.8.8 the first warning dropped the read).
  geno_meta_r <- shiny::eventReactive(input$btn_run, {
    log_msgs <- character(0L)
    note <- function(...) log_msgs <<- c(log_msgs, paste0(...))

    read_one <- function() {
      map_path_extra <- if (!is.null(input$map_file_extra))
        input$map_file_extra$datapath else NULL
      block_path     <- if (!is.null(input$block_file))
        input$block_file$datapath else NULL

      if (input$input_format == "plink_ped") {
        shiny::req(input$ped_file, input$map_file)
        note("Reading Plink PED + MAP...")
        haploStack::read_plink_ped(
          ped_file   = input$ped_file$datapath,
          map_file   = input$map_file$datapath,
          block_file = block_path,
          pos_unit   = input$pos_unit
        )
      } else if (input$input_format == "st") {
        shiny::req(input$geno_file)
        note("Reading SelectionTools genotype file...")
        haploStack::read_st_genotype(
          file       = input$geno_file$datapath,
          map_file   = map_path_extra,
          block_file = block_path,
          pos_unit   = input$pos_unit,
          st_formats = TRUE
        )
      } else {
        shiny::req(input$geno_file)
        note("Reading genomicSimulation genotype file...")
        haploStack::read_gs_genotype(
          file       = input$geno_file$datapath,
          map_file   = map_path_extra,
          block_file = block_path,
          pos_unit   = input$pos_unit
        )
      }
    }

    tryCatch({
      meta <- withCallingHandlers(
        read_one(),
        warning = function(w) {
          note("[WARNING] ", conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          note(trimws(conditionMessage(m)))
          invokeRestart("muffleMessage")
        })
      note("Read complete: ",
           meta$metadata$n_markers, " markers × ",
           meta$metadata$n_samples, " samples (ploidy = ",
           meta$metadata$ploidy, ")")
      list(meta = meta, log = log_msgs)
    }, error = function(e) {
      note("[ERROR] ", conditionMessage(e))
      list(meta = NULL, log = log_msgs)
    })
  })

  # Reactive: lay the GenoMeta out in the chosen preset
  converted_r <- shiny::reactive({
    res <- geno_meta_r()
    shiny::req(!is.null(res$meta))
    haploStack::geno_meta_to_matrix(res$meta, preset = input$output_format,
                                    st_formats = TRUE)
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

  # Download handler: write the converted file(s) to a temp dir and zip them
  output$dl_files <- shiny::downloadHandler(
    filename = function() paste0(input$out_prefix, "_converted.zip"),
    content  = function(file) {
      res <- geno_meta_r()
      shiny::req(!is.null(res$meta))

      tmp_dir <- tempfile("hs_converter_")
      dir.create(tmp_dir)
      on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
      out_path <- file.path(tmp_dir, input$out_prefix)

      haploStack::write_geno_meta(res$meta, path = out_path,
                                  preset = input$output_format,
                                  st_formats = TRUE)

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

      # zip is a Suggests: fall back to the system zip through utils when
      # it is absent (-j stores the files without their temp-dir paths).
      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zip(file, files = unique(written), mode = "cherry-pick")
      } else {
        utils::zip(file, files = unique(written), flags = "-j")
      }
    }
  )
}

shiny::shinyApp(ui = ui, server = server)
