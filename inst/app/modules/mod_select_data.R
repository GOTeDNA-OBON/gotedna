# Select data and show them on the map
mod_select_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML(sprintf("
          #%s .shiny-file-input-btn:hover {
            background-color: white !important;
            color: #2241a7 !important;
            border-color: #2241a7 !important;
          }
        ", ns("download_files")))),
    div(
      id = "data_request",
      div(
        class = "section_header",
        div(
          # title and top buttons
          class = "title-container",
          h1("Data request",
               icon("info-circle", class = "definition",
                    title = "Processing speed is impacted if data request is broad. Please select specific taxonomic level or species for increased speed."),
             ),
          div(
            class = "buttons-container",
            actionButton(ns("reset"), "Reset",
              title = "Reset selection to default values"
            ),
            actionButton(ns("hide_fields"), "Show/Hide fields",
              title = "Hide or show fields"
            )
          )
        ),
        # fields
        div(
          id = ns("data_request_fields"),
          class = "inputs-container",
          ## top fields
          div(
            id = "data_request_top_fields",
            fluidRow(
              column(
                3,
                selectInput(ns("datasource"),
                  label = "Data source",
                  choices = list(
                    "GOTeDNA" = "gotedna",
                    "Upload data" = "external_data",
                    "Download OBIS Data" = "download_from_obis"
                  ),
                  selected = "gotedna"
                )
              ),
              column(
                3,
                uiOutput(ns("data_button_placeholder"))
              ),
              column(
                3,
                selectInput(ns("data_type"),
                  label = "Type of data",
                  choices = list(
                    "Multi-species (metabarcoding)" = "metabarcoding"
                    #"Species specific (qPCR)" = "qPCR"
                  ),
                  selected = "metabarcoding"
                )
              )
            )
          ),
          ## bottom fields
          div(
            id = "data_request_bottom_fields",
            fluidRow(
              column(
                3,
                selectInput(
                  ns("taxo_lvl"),
                  "Taxonomy selection",
                  choices = taxonomic_ranks
                )
              ),
              column(
                3, selectInput(ns("taxo_id"),
                  "Taxa",
                  choices = NULL
                )
              ),
              column(
                3,
                selectizeInput(ns("slc_spe"),
                  "Species",
                  choices = "All"
                )
              ),
              column(
                3,
                # https://stackoverflow.com/questions/50218614/shiny-selectinput-to-select-all-from-dropdown
                htmltools::tagQuery(shinyWidgets::pickerInput(ns("primer"),
                  div(
                    "Primer set",
                    icon("info-circle", class = "definition", onclick = "fakeClick('primer-info')",
                         title = "Primer details"),
                  ),
                  choices = "All",
                  options = list(`actions-box` = TRUE),
                  multiple = TRUE
                ))$find(".btn")$removeAttrs("data-toggle")$addAttrs(`data-bs-toggle` = "dropdown")$allTags()
              )
            )
          )
        )
      ),
      div(
        class = "section_footer",
        id = ns("section_footer_data_request"),
        column(12, uiOutput(outputId = ns("n_smpl_data")))
      )
    ),
    div(
      id = "area_selection",
      div(
        class = "section_header",
        div(
          class = "title-container",
          h1(
            "Area selection",
            span(
              id = ns("lock"),
              class = "lock", icon("lock"),
              title = "Map view locked"
            ),
            span(
              id = ns("restrict"),
              class = "restrict", icon("warning"),
              title = "Selection restricted to spatial area"
            )
          ),
          div(
            class = "buttons-container",
            id = "button_map",
            actionButton(ns("confirm"), "Confirm",
              title = "Confirm spatial selection"
            ),
            actionButton(ns("lock"), "Lock view",
              title = "Set and lock map bounds"
            ),
            actionButton(ns("clear_area"), "Clear area",
              title = "Clear current spatial selection"
            ),
            actionButton(ns("hide_map"), "Show/Hide map",
              title = "Hide or show map"
            )
          )
        ),
        uiOutput(ns("map_warning"))
      ),
      div(
        class = "leaflet_container",
        style = "display: grid; justify-content: center; position: relative;",
        id = ns("map_container"),
        mapedit::editModUI(ns("map-select"), height = "75vh", width = "65vw"),
      ),
      div(
        class = "section_footer",
        id = ns("section_footer_area_selection"),
        div(
          actionButton(ns("show_map_info"), "How to select",
            title = "Display information about how to use the map below"
          )
        ),
        div(
          uiOutput(outputId = ns("n_smpl_map"))
        )
      )
    )
  )
}


mod_select_data_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns


    # Generate map
    sf_edits <<- callModule(
      mapedit::editMod,
      leafmap = basemap(),
      id = "map-select"
    )

    ## hide and show
    observeEvent(input$hide_fields, {
      shinyjs::toggle("data_request_fields")
      shinyjs::toggle("section_footer_data_request")
    })

    observeEvent(input$hide_map, {
      shinyjs::toggle("map_container")
      shinyjs::toggle("section_footer_area_selection")
    })

    observe({
      if (is.null(r$geom_slc)) {
        shinyjs::hide("restrict")
      } else {
        shinyjs::show("restrict")
      }
    })

    shinyjs::hide("lock")
    observeEvent(input$lock, {
      shinyjs::toggle("lock")
      r$lock_view <- !(r$lock_view)
    })

    output$map_warning <- renderUI({
      if (r$taxon_selected < 4) {
        div(
          style = "
        margin-bottom: 8px;
        padding: 8px 12px;
        background-color: #fff3cd;
        border: 1px solid #ffeeba;
        border-radius: 4px;
        color: #856404;
        font-weight: 500;
        display: flex;
        justify-content: center;
      ",
          icon("exclamation-triangle"),
          "Map is available after taxon selection."
        )
      } else {
        NULL
      }
    })

    # Dynamically render button based on datasource
    output$data_button_placeholder <- renderUI({
      req(input$datasource)

      if (input$datasource == "external_data") {
        div(
          id = ns("external_files"),
          class = "file_input-container",
          fileInput(
            ns("external_file"),
            "Upload your files",
            multiple = TRUE,
            buttonLabel = "Browse...",
            placeholder = "No file selected"
          ),
          tags$small(
            style = "display:block; margin-top:-12px; margin-bottom:12px; color:#8B0000;",
            HTML('Warning: uploaded data must use the exact column names and formatting shown in the
              <a href="gotedna_data_template.xlsx" download>GOTeDNA template</a>.')
          )
        )
      } else if (input$datasource == "download_from_obis") {
        div(
          id = ns("download_files"),
          class = "file_input-container",
          actionButton(
            ns("download_button"),
            "Download File From OBIS",
            class = "shiny-file-input-btn",
            style = "
              width: 100%;
              margin-top: 32px;
              background-color: #2241a7;
              color: white;
              border-color: #2241a7;
              display: flex;
              align-items: center;
              justify-content: center;
            "
                    ),
          tags$small(
            style = "display:block; margin-top:0px; margin-bottom:12px; color:#8B0000;",
            "Click to download all OBIS GOTeDNA data as one file. This may take hours."
          )
        )
      } else {
        NULL  # No button for gotedna
      }
    })

    observeEvent(input$download_button, {
      showModal(
        modalDialog(
          title = "Confirm download",

          div(
            style = "
          padding: 16px 20px;
          border: 1px solid #e0e0e0;
          border-radius: 6px;
          background-color: #f9fafb;
        ",

            tags$p(
              "This will download a single file containing all data from OBIS that are suitable for GOTeDNA."
            ),

            tags$p(
              style = "margin-top: 8px;",
              "After that you can upload that file into the app and explore the data in GOTeDNA."
            ),

            tags$p(
              style = "margin-top: 8px;",
              paste0("The data currently in GOTeDNA were pulled from OBIS at ", last_obis_download_ts, ". A fresh download is only useful if usable data have been uploaded to OBIS since then.")
            ),

            tags$p(
              style = "margin-top: 8px;",
              "The request may take hours depending on server load."
            ),

            tags$div(
              style = "
            margin-top: 16px;
            padding: 12px;
            background-color: #fff3cd;
            border: 1px solid #ffeeba;
            border-radius: 4px;
            color: #856404;
          ",
              icon("warning"),
              " This is a potentially long-running operation and cannot be stopped from the app (must be cancelled on the server or in rStudio)."
            )
          ),

          footer = tagList(
            modalButton("Cancel"),
            downloadButton(
              ns("confirm_download"),
              "OK, download",
              class = "btn-primary",
              style = "
                display: inline-flex;
                align-items: center;
                justify-content: center;
              "
            )
          ),

          easyClose = TRUE
        )
      )
    })


    # File upload handling (existing logic)
    observeEvent(input$external_file, {
      req(input$external_file)
      df <- read_uploaded_file(input$external_file[1, ])
      showModal(modalDialog(
        title = "Please wait",
        "Processing coordinates into station clusters. This may take several minutes for larger files (e.g. >20MB).",
        footer = NULL,
        easyClose = FALSE
      ))
      on.exit(removeModal(), add = TRUE)

      df_with_assigned_stations <- update_station_variable(df)
      r$upload_data <- df_with_assigned_stations
      r$upload_stations <- get_station(df_with_assigned_stations)
      r$upload_primers <- newprob_mb <- calc_det_prob(r$upload_data)
      scaledprobs_mb <- scale_newprob(r$upload_data, newprob_mb)

      upload_gotedna_primer <- list()
      for (i in c("kingdom", "phylum", "class", "order", "family", "genus", "scientificName")) {
        upload_gotedna_primer[[i]] <- primer_sort(i, scaledprobs_mb) |>
          mutate(text = paste0(primer, " (", detects, "/", total, " ", perc, "%)"))
      }

      r$upload_primers <- upload_gotedna_primer
      r$cur_data <- r$upload_data
      r$data_station <- r$upload_stations
      r$protocol_ID <- paste0(r$cur_data$protocol_ID)
    })


    output$confirm_download <- downloadHandler(
      filename = function() {
        paste0(
          "obis_dataset_",
          Sys.Date(),
          ".csv"
        )
      },

      content = function(file) {

        removeModal()  # close confirmation modal

        showModal(
          modalDialog(
            title = "Downloading data",

            div(
              style = "
        padding: 16px;
        background-color: #f9fafb;
        border: 1px solid #e0e0e0;
        border-radius: 6px;
      ",
              p("Downloading all GOTeDNA data from OBIS."),
              p("This may take hours, depending on data amount and bandwidth."),
              tags$div(
                style = "
          margin-top: 12px;
          padding: 10px;
          background-color: #eef2ff;
          border-left: 4px solid #2241a7;
        ",
                icon("info-circle"),
                " You may close this dialog. The download will continue in the background. If you would like to stop the download, you need to stop the Shiny app (e.g. interrupt it in rStudio)."
              )
            ),

            footer = tagList(
              modalButton("Close"),
              span(style = "margin-left: 8px; font-size: 0.9em; color: #666;",
                   "Closing will not stop the download")
            ),

            easyClose = TRUE
          )
        )

        on.exit(removeModal(), add = TRUE)

        # --- ACTUAL DOWNLOAD ---
        obis_data <- tryCatch(
          big_OBIS_data_pull(),
          error = function(e) {
            showNotification(
              paste("Download failed:", e$message),
              type = "error"
            )
            stop(e)
          }
        )

        write.csv(obis_data, file, row.names = FALSE)
      }
    )


    ## load data
    observe({
      if (input$datasource == "gotedna") {
        shinyjs::hide("external_files")
        # Keep GOTeDNA data as before
        gotedna_data <- gotedna_data0
        gotedna_station <- gotedna_station0
        r$data_type <- input$data_type
        r$cur_data <- gotedna_data[[input$data_type]]
        r$data_station <- gotedna_station[[input$data_type]]
        r$protocol_ID <- paste0(
          r$cur_data$protocol_ID
        )

        r$primer_choices_all <- get_primer_selection(
          r$taxon_lvl_slc,
          filter_taxon(
            isolate(r$cur_data),
            r$taxon_lvl_slc,
            r$taxon_id_slc,
            r$scientificName
          ),
          gotedna_primer
        )
        # Disable fileInput and visually gray it out
        shinyjs::disable("external_file")

      } else {
        shinyjs::show("external_files")
        gotedna_data <- gotedna_data0
        gotedna_station <- gotedna_station0
        shinyjs::enable("external_file")
        req(r$upload_data)

        r$cur_data <- r$upload_data
        r$data_station <- r$upload_stations

        r$primer_choices_all <- get_primer_selection(
          r$taxon_lvl_slc,
          filter_taxon(
            isolate(r$cur_data),
            r$taxon_lvl_slc,
            r$taxon_id_slc,
            r$scientificName
          ),
          r$upload_primers
        )
        r$protocol_ID <- paste0(
          r$cur_data$protocol_ID
        )
      }
    })


    observeEvent(input$data_type, {
      r$data_type <- input$data_type
      r$cur_data <- gotedna_data[[input$data_type]]
      r$data_station <- gotedna_station[[input$data_type]]
      r$protocol_ID <- paste0(
        r$cur_data$protocol_ID
      )
    })

    observe(
      if (!is.null(r$station_slc)) {
        r$cur_data_sta_slc <- r$cur_data |>
          dplyr::filter(station %in% r$station_slc)
      } else {
        r$cur_data_sta_slc <- r$cur_data
      }
    )

    ## update Taxon data
    observe({
      updateSelectInput(
        session,
        "taxo_id",
        choices = c(
          r$cur_data[[input$taxo_lvl]] |>
            unique() |>
            sort()
        ),
        selected = NULL
      )
    })

    observeEvent(input$taxo_id, {
      if (input$taxo_id != "All") {
        updateSelectizeInput(
          session,
          "slc_spe",
          choices = c(
            "All",
            r$cur_data %>%
              filter(
                .data[[input$taxo_lvl]] == input$taxo_id,
                str_detect(scientificName, " ")  # only names with a space
              ) %>%
              pull(scientificName) %>%
              unique() %>%
              sort()
          ),
          selected = "All"
        )
      } else {
        updateSelectizeInput(
          session,
          "slc_spe",
          choices = c(
            "All",
            r$cur_data_sta_slc$scientificName |>
              unique() |>
              sort()
          ),
          selected = "All",
        )
      }
    })

    observeEvent(input$reset, {
      updateSelectInput(
        session,
        "taxo_lvl",
        selected = "kingdom"
      )
      updateSelectizeInput(
        session,
        "slc_spe",
        selected = "All"
      )
      r$reset <- r$reset + 1
      r$data_ready <- NULL
    })


    observe({
      r$primer_choices_all <- get_primer_selection(
        r$taxon_lvl_slc,
        filter_taxon(
          isolate(r$cur_data_sta_slc),
          r$taxon_lvl_slc,
          r$taxon_id_slc,
          r$scientificName
        ),
      )
      shinyWidgets::updatePickerInput(
        session,
        "primer",
        choices = r$primer_choices_all,
        selected = r$primer_choices_all
      )
    })

    observeEvent(input$primer, {
      r$primer <- input$primer
    })

    # sample number
    output$n_smpl_map <- renderUI({
      tagList(
        div(
          class = "sample_selected",
          p(
            "Total number of samples: ",
            span(
              class = "sample_selected_map",
              format(r$n_sample, big.mark = ",")
            )
          )
        )
      )
    })


    listenMapData <- reactive({
      list(
        input$taxo_id,
        input$slc_spe,
        input$data_type,
        input$primer
      )
    })

    observeEvent(listenMapData() |> debounce(100), {
      r$scientificName <- input$slc_spe
      r$taxon_id_slc <- input$taxo_id
      r$taxon_selected <- r$taxon_selected + 1
      if (!is.null(input$slc_spe) && input$slc_spe != "All") {
        r$taxon_lvl_slc <- "scientificName"
      } else {
        r$taxon_lvl_slc <- input$taxo_lvl
      }
      r$fig_ready <- FALSE
      # count data
      r$geom <- filter_station(r)
      r$reload_map <- r$reload_map + 1
    })

    observeEvent(input$confirm, {
      if (!is.null(sf_edits()$all)) {
        cli::cli_alert_info("Using geom(s) drawn to select region")
        if (!is.null(r$geom)) {
          id_slc <- sf::st_contains(sf_edits()$all, r$geom, sparse = FALSE) |>
            apply(2, any)
          if (sum(id_slc)) {
            r$geom <- r$geom[id_slc, ]
            r$geom_slc <- sf_edits()$all
            r$station_slc <- r$geom$station

            geom_coords <- sf::st_bbox(r$geom_slc)
          } else {
            showNotification("No station selected", type = "warning")
          }
        } else {
          cli::cli_alert_info("`r$geom` is null")
        }
      } else {
        cli::cli_alert_info("`sf_edits()$all` is null")
        showNotification("Empty spatial selection", type = "warning")
      }
      r$reload_map <- r$reload_map + 1
    })


    observeEvent(input$show_map_info, r$show_map_info <- TRUE)

    observeEvent(r$geom, {
      r$n_sample <- sum(r$geom$count)
    })

    observeEvent(input$clear_area, {
      r$geom_slc <- r$station_slc <- NULL
      r$primer_choices_all <- get_primer_selection(
        r$taxon_lvl_slc,
        filter_taxon(
          isolate(r$cur_data),
          r$taxon_lvl_slc,
          r$taxon_id_slc,
          r$scientificName
        ),
        gotedna_primer
      )
      r$geom <- filter_station(r)
      # Reload module
      sf_edits <<- callModule(
        mapedit::editMod,
        leafmap = basemap(),
        id = "map-select"
      )
      # shinyWidgets::updatePickerInput(
      #   session,
      #   "primer",
      #   selected = r$primer
      # )
      r$reload_map <- r$reload_map + 1
    })

    observeEvent(r$reload_map, {
      prx <- leaflet::leafletProxy("map-select")
      prx$id <- "slc_data-map-select-map"
      update_map(prx, isolate(r$geom), isolate(r$geom_slc), lock_view = r$lock_view)
    })
  })
}


# filter via inner join and used to count samples
filter_station <- function(r) {
  if (length(r$station_slc)) {
    sta <- r$data_station |>
      dplyr::filter(station %in% r$station_slc)
  } else {
    sta <- r$data_station
  }
  dff <- filter_taxon(
    r$cur_data, r$taxon_lvl_slc, r$taxon_id_slc, r$scientificName,
    r$primer
  ) |>
    dplyr::filter(primer %in% r$primer)
  sta |>
    dplyr::inner_join(
      dff |>
        dplyr::group_by(station, samp_name) |>
        dplyr::summarise(
          count = dplyr::n_distinct(samp_name),
          success = sum(detected)
        ),
      # dplyr::group_by(station, primer, species) |>
      # dplyr::summarise(
      #   success = sum(detected),
      #   count = dplyr::n_distinct(samp_name)),
      join_by(station)
    )
}


filter_taxon <- function(data, taxon_lvl, taxon_id, scientificName, primer) {
  out <- data
  if (!is.null(taxon_lvl)) {
    if (taxon_lvl == "scientificName") {
      out <- out |>
        dplyr::filter(scientificName == {{ scientificName }}, primer == primer)
    } else {
      if (taxon_id != "All") {
        out <- out |>
          dplyr::filter(!!dplyr::ensym(taxon_lvl) == taxon_id, primer == primer)
      }
    }
  }
  out
}

update_map <- function(proxy, geom, geom_slc, lock_view = FALSE) {
  proxy |>
    leaflet::clearGroup("station") |>
    leaflet::clearGroup("select_polygon") |>
    leaflet::addMarkers(
      data = geom,
      clusterOptions = leaflet::markerClusterOptions(),
      label = ~ paste(success, "observations"),
      group = "station"
    )
  if (!is.null(geom_slc)) {
    geom_type <- geom_slc |>
      sf::st_geometry_type() |>
      as.character()
    ind <- geom_type == "POLYGON"
    if (length(ind)) {
      proxy |>
        leaflet::addPolygons(
          data = geom_slc[ind, ],
          color = "#75f9c6",
          fillOpacity = 0.1,
          group = "select_polygon"
        )
    } else {
      showNotification(
        "Only rectangles and polygons can be used",
        type = "warning"
      )
    }
  }
  if (!lock_view) {
    bb <- geom |>
      sf::st_bbox() |>
      as.vector()
    proxy |>
      leaflet::fitBounds(bb[1], bb[2], bb[3], bb[4])
  }
}

read_uploaded_file <- function(fileinfo) {
  stopifnot(nrow(fileinfo) == 1)

  path <- fileinfo$datapath
  name <- fileinfo$name
  ext  <- tolower(tools::file_ext(name))

  if (ext == "csv") {
    read.csv(path, stringsAsFactors = FALSE)
  } else if (ext %in% c("xls", "xlsx")) {
    readxl::read_excel(path) |> as.data.frame()
  } else {
    stop("Unsupported file type: ", ext)
  }
}


