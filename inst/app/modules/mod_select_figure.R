# Select data and show them on the map
mod_select_figure_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(
      id = "figure_selection",
      div(
        class = "section_header",
        div(
          class = "title-container",
          h1("Figure selection"),
          div(
            class = "buttons-container",
            actionButton(ns("hide_figs"), "Hide/Show figures",
              title = "Hide or show fields"
            )
          )
        )
      ),
      div(
        id = ns("figure_selection_main"),
        div(
          class = "d-flex justify-content-center buttons-container",
          id = "select_all_figures",
          actionButton(
            ns("select_all"),
            "Select all figures",
            title = "Select all figures"
          ),
          actionButton(
            ns("deselect_all"),
            "Deselect all figures",
            title = "Deselect all figures"
          )
        ),
        div(
          id = "thumbnail_container",
          fluidRow(
            class = "justify-content-center",
            add_figure_selection(
              ns("fig_detect"),
              "Monthly eDNA detection probability",
              "img/thumbnails/tn_thresh.svg"
            ),
            add_figure_selection(
              ns("fig_effort"),
              "Guidance on sampling effort",
              "img/thumbnails/tn_effort.svg"
            ),
            add_figure_selection(
              ns("fig_heatmap"),
              "Species detection heatmap",
              "img/thumbnails/tn_heatmap.svg"
            ),
            add_figure_selection(
              ns("fig_samples"),
              "Data variation",
              "img/thumbnails/tn_samples.svg"
            )
          )
        ),
        div(
          class = "d-flex justify-content-center",
          id = "confirm_figures_selection",
          actionButton(
            ns("compute_and_visualize"),
            label = "Compute & Visualize",
            title = "Compute optimal detection window",
            class = "primary-button"
          )
        ),
      )
    ),
    div(
      id = "observation_protocol",
      div(
        id = "protocol_panel",
        class = "title-container",
        h1("Protocol Selection"),
        div(
          class = "buttons-container",
          actionButton(ns("hide_prots"), "Hide/Show Protocols",
                       title = "Hide or show fields"
          )
        )
      ),

      fluidRow(
        id = ns("protocol_section"),
        class = "",
        column(
          3,
          class = "control-panel",
          div(
            id = "fig_left_panel",
            shinyWidgets::pickerInput(
              inputId = ns("prot_id"),
              label = tagList(
                "Protocol ID "
              ),
              choices = "Not available",
              multiple = TRUE,
              options = list(
                `actions-box` = TRUE,      # select all / deselect all
                `live-search` = FALSE       # search bar (nice if we have lots of protocols)
              )
            )
          )
        ),
        column(
          width = 8,
          uiOutput(ns("protocol_details"))
        )
      )
    ),

    div(
      id = "observation",
      fluidRow(
        class = "panels-container",
        column(
          3,
          class = "control-panel",
          div(
            id = "obs_panel",
            class = "section_header",
            h1("Observation")
          ),
          div(
            id = "fig_left_panel",
            selectInput(
              ns("threshold"),
              "Threshold",
              choices = ls_threshold,
              selected = 75
            ),
            div(
              id = "fig_sampling_info",
              h4("Guidance"),
              div(
                class = "sampling_info-item",
                h6("Optimal timing: "),
                uiOutput(ns("opt_sampl"), class = "fig_text_output")
              ),
              div(
                class = "sampling_info-item",
                h6("Confidence: "),
                uiOutput(ns("conf"), class = "fig_text_output")
              ),
              div(
                class = "sampling_info-item",
                h6("Consistency among years: "),
                uiOutput(ns("var_year"), class = "fig_text_output"),
              )
            ),
            # downloadButton(
            #   ns("export_pdf"),
            #   "Export to PDF",
            #   title = "Export figures to PDF",
            #    class = "primary-button"
            # )
          )
        ),
        column(
          9,
          class = "show-panels",
          div(
            class = "fig_main_container",
            div(
              class = "fig_main_container-fig",
              ui_fig_smooth("fig_smooth", "Monthly eDNA Detection Probability", c("sample_size.html", "detection.html"), ns),
              ui_fig_detect_bottom("fig_detect", ns),
              ui_fig_effort("fig_effort", "Guidance on sampling effort", "sample_size.html", ns),
              ui_fig_hm("fig_heatmap", "Species detection heatmap", "heatmap.html", ns),
              ui_fig_samples("fig_samples", "Data variation", "field_sample.html", ns)
            )
          )
        )
      )
    ),
    div(
      id = "reference_data_authorship",
      div(
        class = "table_title-container",
        h1("Reference data authorship")
      ),
      DT::DTOutput(ns("data_authorship"))
    ),
    # )
    # )
  )
}


mod_select_figure_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$hide_figs, {
      shinyjs::toggle("figure_selection_main")
    })

    observeEvent(input$hide_prots, {
      shinyjs::toggle("protocol_section")
    })

    observeEvent(input$select_all, {
      for (i in c("fig_smooth", "fig_detect", "fig_effort", "fig_heatmap", "fig_samples")) {
        show_fig(i)
        r$fig_slc[[i]] <- TRUE
      }
    })

    observeEvent(input$deselect_all, {
      for (i in c("fig_smooth", "fig_detect", "fig_effort", "fig_heatmap", "fig_samples")) {
        hide_fig(i)
        r$fig_slc[[i]] <- FALSE
      }
    })

    observeEvent(r$reset,
      {
        for (i in c("fig_smooth", "fig_detect", "fig_effort", "fig_heatmap", "fig_samples")) {
          hide_fig(i)
          r$fig_slc[[i]] <- FALSE
        }
      },
      ignoreInit = TRUE
    )

    hide_fig("fig_smooth")
    hide_fig("fig_detect")
    observeEvent(input$fig_detect, {
      toggle_legends("thresh_legend")
      toggle_fig("fig_detect")
      toggle_fig("fig_smooth")
      r$fig_slc$fig_smooth <- !r$fig_slc$fig_smooth
      r$fig_slc$fig_detect <- !r$fig_slc$fig_detect
    })

    hide_fig("fig_effort")
    observeEvent(input$fig_effort, {
      toggle_fig("fig_effort")
      r$fig_slc$fig_effort <- !r$fig_slc$fig_effort
    })

    hide_fig("fig_heatmap")
    observeEvent(input$fig_heatmap, {
      toggle_fig("fig_heatmap")
      r$fig_slc$fig_heatmap <- !r$fig_slc$fig_heatmap
    })

    hide_fig("fig_samples")
    observeEvent(input$fig_samples, {
      toggle_fig("fig_samples")
      r$fig_slc$fig_samples <- !r$fig_slc$fig_samples
    })

    # THIS RUNS WHEN WE CLICK COMPUTE & VISUALIZE
    observeEvent(
      ignoreInit = TRUE,
      list(input$compute_and_visualize),
      {
        req(input$compute_and_visualize)
        r$frozen_selected_taxon_level <- isolate(r$taxon_lvl_slc)
        r$frozen_selected_taxon_id <- isolate(r$taxon_id_slc)

        update_data_active()
        update_protocol_menu()
        req(input$prot_id)
        validate(
          need(input$prot_id != "Not available", "Protocol not selected yet")
        )

        r$data_ready <- prepare_data(r) |>
          filter(protocol_ID %in% input$prot_id)

        if (nrow(r$data_ready)) {
          showNotification(
            paste0(
              "Computing detection window with threshold set to ",
              input$threshold, "%",
              ifelse(
                nrow(r$data_ready) > 1e4,
                paste0(" (", nrow(r$data_ready), " observations, this may take some time)"),
                ""
              )
            ),
            type = "message",
            duration = NULL,
            id = "notif_calc_win"
          )

          # Compute detection probability
          r$newprob <- calc_det_prob(r$data_ready, r$frozen_selected_taxon_level, r$frozen_selected_taxon_id, pool_primers = TRUE)
          # Safely scale probabilities
          r$scaledprobs <- tryCatch({
            if (length(r$newprob$newP_agg) == 0 && length(r$newprob$newP_yr) == 0) {
              cat("calc_det_prob returned empty")
              NULL
            } else {
              scale_newprob(r$data_ready, r$newprob, r$frozen_selected_taxon_level)
            }
          }, error = function(e) {
            cat("scale_newprob failed:", conditionMessage(e), "\n")
            NULL
          })


          # Do the same for scientificName if level == genus
          if(r$frozen_selected_taxon_level == "genus") {
            r$newprob_by_scientificName <- calc_det_prob(r$data_ready, "scientificName", "All", pool_primers = TRUE)

            r$scaledprobs_by_scientificName <- tryCatch({
              if (length(r$newprob_by_scientificName$newP_agg) == 0 && length(r$newprob_by_scientificName$newP_yr) == 0) {
                cat("calc_det_prob returned empty for level: scientificName\n")
                NULL
              } else {
                scale_newprob(r$data_ready, r$newprob_by_scientificName, "scientificName")
              }
            }, error = function(e) {
              cat("scale_newprob failed for scientificName:", conditionMessage(e), "\n")
              NULL
            })
          }


          cli::cli_alert_info("Computing optimal detection window")

          thresh_slc <- input$threshold
          win <- calc_window(
            threshold = input$threshold,
            scaledprobs = r$scaledprobs
          )

          j.sim <- jaccard_test(
            r$scaledprobs,
            input$threshold
          )

          if (is.null(win)) {
            #  showNotification("No optimal detection window", type = "warning")
            output$opt_sampl <- renderUI("No single detection window")
            output$conf <- renderUI("NA")
            output$var_year <- renderUI(
              paste(j.sim)
            )
          } else {
            output$opt_sampl <- renderUI(
              paste(win$opt_sampling$period)
            )
            output$conf <- renderUI(paste(win$fshTest$confidence))
            output$var_year <- renderUI(
              paste(j.sim)
            )
          }


          r$fig_ready <- TRUE

          session$onFlushed(function() {
            removeNotification(id = "notif_calc_win")
          }, once = TRUE)

          # create protocol vector

          v_prot <- r$scaledprobs$protocol_ID |> unique()
        } else {
          # showNotification("Data selection is empty", type = "warning")
        }

        shinyscreenshot::screenshot(
          selector = "#data_request_top_fields",
          filename = "data_top",
          download = FALSE, server_dir = tempdir()
        )

        shinyscreenshot::screenshot(
          selector = "#data_request_bottom_fields",
          filename = "data_btm",
          download = FALSE, server_dir = tempdir()
        )

        shinyscreenshot::screenshot(
          selector = "#reference_data_authorship",
          filename = "dat_auth",
          download = FALSE, server_dir = tempdir()
        )
      }
    )

    #THIS RUNS WHEN WE CHANGE PROTOCOL ID OR THRESHOLD
    observeEvent(
      ignoreInit = TRUE,
      list(input$threshold, input$prot_id),
      {
        req(input$compute_and_visualize)
        r$frozen_selected_taxon_level <- isolate(r$taxon_lvl_slc)
        r$frozen_selected_taxon_id <- isolate(r$taxon_id_slc)

        req(input$prot_id)
        validate(
          need(input$prot_id != "Not available", "Protocol not selected yet")
        )

        r$data_ready <- prepare_data(r) |>
          filter(protocol_ID == input$prot_id)

        if (nrow(r$data_ready)) {
          showNotification(
            paste0(
              "Computing detection window with threshold set to ",
              input$threshold, "%",
              ifelse(
                nrow(r$data_ready) > 1e4,
                paste0(" (", nrow(r$data_ready), " observations, this may take some time)"),
                ""
              )
            ),
            type = "message",
            duration = NULL,
            id = "notif_calc_win"
          )

          # Compute detection probability
          r$newprob <- calc_det_prob(r$data_ready, r$frozen_selected_taxon_level, r$frozen_selected_taxon_id, pool_primers = TRUE)
          # Safely scale probabilities
          r$scaledprobs <- tryCatch({
            if (length(r$newprob$newP_agg) == 0 && length(r$newprob$newP_yr) == 0) {
              cat("calc_det_prob returned empty")
              NULL
            } else {
              scale_newprob(r$data_ready, r$newprob, r$frozen_selected_taxon_level)
            }
          }, error = function(e) {
            cat("scale_newprob failed:", conditionMessage(e), "\n")
            NULL
          })


          # Do the same for scientificName if level == genus
          if(r$frozen_selected_taxon_level == "genus") {
            r$newprob_by_scientificName <- calc_det_prob(r$data_ready, "scientificName", "All", pool_primers = TRUE)

            r$scaledprobs_by_scientificName <- tryCatch({
              if (length(r$newprob_by_scientificName$newP_agg) == 0 && length(r$newprob_by_scientificName$newP_yr) == 0) {
                cat("calc_det_prob returned empty for level: scientificName\n")
                NULL
              } else {
                scale_newprob(r$data_ready, r$newprob_by_scientificName, "scientificName")
              }
            }, error = function(e) {
              cat("scale_newprob failed for scientificName:", conditionMessage(e), "\n")
              NULL
            })
          }


          cli::cli_alert_info("Computing optimal detection window")

          thresh_slc <- input$threshold
          win <- calc_window(
            threshold = input$threshold,
            scaledprobs = r$scaledprobs
          )

          j.sim <- jaccard_test(
            r$scaledprobs,
            input$threshold
          )

          if (is.null(win)) {
            #  showNotification("No optimal detection window", type = "warning")
            output$opt_sampl <- renderUI("No single detection window")
            output$conf <- renderUI("NA")
            output$var_year <- renderUI(
              paste(j.sim)
            )
          } else {
            output$opt_sampl <- renderUI(
              paste(win$opt_sampling$period)
            )
            output$conf <- renderUI(paste(win$fshTest$confidence))
            output$var_year <- renderUI(
              paste(j.sim)
            )
          }


          r$fig_ready <- TRUE

          session$onFlushed(function() {
            removeNotification(id = "notif_calc_win")
          }, once = TRUE)

          # create protocol vector

          v_prot <- r$scaledprobs$protocol_ID |> unique()
        } else {
          # showNotification("Data selection is empty", type = "warning")
        }

        shinyscreenshot::screenshot(
          selector = "#data_request_top_fields",
          filename = "data_top",
          download = FALSE, server_dir = tempdir()
        )

        shinyscreenshot::screenshot(
          selector = "#data_request_bottom_fields",
          filename = "data_btm",
          download = FALSE, server_dir = tempdir()
        )

        shinyscreenshot::screenshot(
          selector = "#reference_data_authorship",
          filename = "dat_auth",
          download = FALSE, server_dir = tempdir()
        )
      }
    )

    update_data_active <- function() {
      out <- r$cur_data_sta_slc
      if (!is.null(r$frozen_selected_taxon_level)) {
        if (r$frozen_selected_taxon_level == "scientificName") {
          out <- out |>
            dplyr::filter(scientificName == r$scientificName)
        } else {
          if (r$frozen_selected_taxon_id != "All") {
            out <- out[
              out[[r$frozen_selected_taxon_level]] == r$frozen_selected_taxon_id,
            ]
          }
        }
      }

      # primer-based subset
      r$data_active <- out |>
        dplyr::filter(primer %in% r$primer)
      # prevent computation when new data are selected
      r$fig_ready <- FALSE
    }

    #observe({ update_data_active() })

    update_protocol_menu <- function() {
      if (!is.null(r$data_active) && nrow(r$data_active) > 0) {
        protocols_by_sample <- r$data_active %>%
          dplyr::distinct(protocol_ID, samp_name)  # unique protocol × sample

        v_prot <- protocols_by_sample$protocol_ID |> table()
        # Sort protocols by decreasing count
        sorted_prot <- sort(v_prot, decreasing = TRUE)
        # Create choices: display name -> value
        l_prot <- as.list(names(sorted_prot))
        names(l_prot) <- paste0(
          "Protocol ",
          names(sorted_prot),
          " (", sorted_prot, " samples)"
        )

        selected_prot <- names(sorted_prot)[1]  # value with most observations

        shinyWidgets::updatePickerInput(
          session,
          "prot_id",
          choices = l_prot,
          selected = selected_prot
        )

      } else {
        # No data → empty the menu
        shinyWidgets::updatePickerInput(
          session,
          "prot_id",
          choices = list(),
          selected = NULL
        )
      }
    }

    observe({ update_protocol_menu() })

    selected_protocol_rows <- reactive({
      req(input$prot_id)

      protocol_info %>%
        dplyr::filter(
          protocol_ID == input$prot_id
        )
    })

    # --- base protocol per Location (the "starting choice") ---
    base_protocol <- reactiveVal(4)

    # live details card
    output$protocol_details <- renderUI({
      df <- selected_protocol_rows()

      if (nrow(df) == 0) {
        return(tags$div(style="margin-top:10px;", em("No rows found for this ProtocolID.")))
      }

      method_groups <- list(
        "Field Methods" = c("samp_size","size_frac","filter_material","samp_mat_process",
                            "minimumDepthInMeters","maximumDepthInMeters"),
        "Storage Methods" = c("samp_store_temp","samp_store_sol"),
        "Lab Methods" = c("target_gene","pcr_primer_forward","pcr_primer_reverse","nucl_acid_ext_kit"),
        "Library Preparation" = c("platform","instrument","seq_kit"),
        "Bioinformatic Methods" = c("otu_db","tax_assign_cat","otu_seq_comp_appr")
      )

      pick_display <- function(x) {
        x <- as.character(x)
        x <- x[!is.na(x) & nzchar(trimws(x))]
        if (length(x) == 0) return(NA_character_)
        ux <- unique(x)
        if (length(ux) == 1) ux else paste(ux, collapse = " | ")
      }

      # Build base df (for comparison)
      bp <- base_protocol()
      df_base <- NULL
      if (!is.null(bp) && nzchar(bp)) {
        df_base <- protocol_info %>%
          dplyr::filter(protocol_ID == bp)
        if (nrow(df_base) == 0) df_base <- NULL
      }

      # UI

      tags$div(
        class = "protocol-details-wrap",
        tagList(
          lapply(names(method_groups), function(group_name) {
            fields <- method_groups[[group_name]]

            # Keep only fields that actually exist in df
            fields <- fields[fields %in% names(df)]
            if (length(fields) == 0) return(NULL)  # hide empty groups

            tags$div(
              tags$h5(class = "protocol-group-title", group_name),

              tags$div(
                class="protocol-grid",

                lapply(fields, function(f) {
                  cur <- pick_display(df[[f]])
                  bas <- if (!is.null(df_base) && f %in% names(df_base)) pick_display(df_base[[f]]) else NA_character_

                  cur2 <- ifelse(is.na(cur), "", trimws(as.character(cur)))
                  bas2 <- ifelse(is.na(bas), "", trimws(as.character(bas)))

                  changed <- nzchar(cur2) && nzchar(bas2) && !identical(cur2, bas2)
                  is_na   <- !nzchar(cur2)

                  display_value <- if (nchar(cur2) > 60) {
                    paste0(substr(cur2, 1, 60), "…")
                  } else {
                    cur2
                  }
                  has_url <- grepl("https?://", cur2)

                  display_value <- if (has_url) {
                    # Extract URL
                    url <- sub(".*(https?://[^ )]+).*", "\\1", cur2)

                    # Extract label before URL
                    label <- trimws(gsub("\\(https?://.*\\)", "", cur2))

                    # Fallback if label is empty
                    if (!nzchar(label)) label <- "View link"

                    tags$a(
                      href = url,
                      target = "_blank",
                      label
                    )
                  } else {
                    cur2
                  }


                  tags$div(
                    tags$div(
                      class="protocol-field-title",
                      protocol_labels[[f]] %||% f
                    ),

                    tags$div(
                      class = paste("protocol-card", if (changed) "changed", if (is_na) "na"),
                      if (is_na) "—" else display_value
                    )
                  )
                })
              ),
              tags$hr(style="border: 0; border-bottom: 1px solid #444; margin: 10px 0;")
            )
          })
        )
      )
    })







    # FIGURES
    ## DETECTION part 1
    output$fig_smooth_plot_output <- renderPlot({
      if (req(r$fig_ready, cancelOutput = TRUE)) {
        draw_fig_smooth(isolate(r), r$fig_ready && r$fig_slc$fig_detect,
          id = input$prot_id # using prot_id as string
        )
      }
    })
    ## DETECTION part 2
    output$fig_detect_plot_output <- renderPlot({
      if (req(r$fig_ready, cancelOutput = TRUE)) {
        draw_fig_detect(r, r$fig_ready && r$fig_slc$fig_detect, input$threshold)
      }
    })

    ## EFFORT NEEDED
    output$fig_effort_plot_output <- plotly::renderPlotly({
      req(r$fig_ready, cancelOutput = TRUE)
      if (r$frozen_selected_taxon_level == 'genus') {
        effort_needed_fig(r$scaledprobs_by_scientificName, selected_taxon_level = "scientificName")
      } else {
        effort_needed_fig(r$scaledprobs, selected_taxon_level = r$frozen_selected_taxon_level)
      }
    })


    ## HEATMAP
    output$fig_heatmap_plot_output <- plotly::renderPlotly({
      if (req(r$fig_ready, cancelOutput = TRUE)) {
        plt_ready <- r$fig_ready && r$fig_slc$fig_heatmap
        if (r$frozen_selected_taxon_level == 'genus') {
          hm_fig(r$scaledprobs_by_scientificName, selected_taxon_level = "scientificName")
        } else {
          hm_fig(r$scaledprobs, selected_taxon_level = r$frozen_selected_taxon_level)
        }
      }
    })

    ## DATA VARIATION
    output$fig_samples_plot_output <- plotly::renderPlotly({
      req(input$year_selected)

      # Keep your original logic
      plt_ready <- r$fig_ready && r$fig_slc$fig_samples
      num_of_species <- r$data_ready$scientificName |> unique() |> length()

      if (num_of_species > 26) {
        plotly::plot_ly(
          type = "scatter",
          mode = "text",
          text = "Too many species to plot for this taxon. Please restrict your search to less than 26 species.",
          x = 0, y = 0,
          textfont = list(size = 12)
        ) %>%
          plotly::layout(
            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
          )
      } else {
        # Wrap draw_fig_samples in tryCatch to catch runtime errors without breaking the app
        tryCatch({
          ggp <- draw_fig_samples(r$newprob_by_scientificName, plt_ready, year = input$year_selected)
        }, error = function(e) {
          message("*** ERROR in draw_fig_samples: ", conditionMessage(e))
          return(plotly::plot_ly(type = "scatter", mode = "text",
                                 text = "Error rendering figure", x = 0, y = 0))
        })
      }
    })

    output$fig_samples_controls <- renderUI({
      req(r$data_ready)

      yrs <- sort(unique(r$data_ready$year))

      div(
        style = "max-width: 150px; margin-bottom: 10px;",  # narrower + small margin
        selectInput(
          ns("year_selected"),
          "Year",
          choices = yrs[order(yrs, decreasing = TRUE)],
          selected = max(yrs)
        )
      )
    })

    # output$debug_field_sample <- renderPrint({
    #   req(input$year_selected)
    #
    #   # Filter data exactly as for the figure
    #   fig_data <- r$data_ready |>
    #     filter(protocol_ID == input$prot_id, year == input$year_selected)
    #
    #   # Show structure and some summaries
    #   cat("=== field_sample_fig debug ===\n")
    #   cat("Rows:", nrow(fig_data), "\n")
    #   cat("Cols:", paste(names(fig_data), collapse=", "), "\n\n")
    #
    #   # Show a brief summary without printing everything
    #   str(fig_data)
    #
    #   # Optional: show first few rows
    #   head(fig_data, 10)
    # })


    # DATA AUTHORSHIP TABLE
    output$data_authorship <- DT::renderDT({
      req(!is.null(r$cur_data), nrow(r$cur_data) > 0)
      req(!is.null(input$prot_id))
      # req(input$calc_window)
      authorship_data <- r$cur_data
      authorship_data <- authorship_data[rowSums(is.na(authorship_data)) != ncol(authorship_data), ]
      if (length(r$station_slc)) {
        authorship_data <- authorship_data |> dplyr::filter(station %in% r$station_slc)
      }
      if (!is.null(r$taxon_lvl_slc)) {
        if (r$taxon_lvl_slc == "scientificName") {
          authorship_data <- authorship_data |> dplyr::filter(scientificName == r$scientificName)
        } else {
          if (r$taxon_id_slc != "All") {
            authorship_data <- authorship_data[
              authorship_data[[r$taxon_lvl_slc]] == r$taxon_id_slc,
            ]
          }
        }
      }

      # Remove all-NA rows
      all_na_rows <- apply(authorship_data, 1, function(x) all(is.na(x)))
      authorship_data <- authorship_data[!all_na_rows, ]

      authorship_data <- authorship_data %>%
        dplyr::mutate(
          LClabel = as.character(LClabel),
          LClabel = dplyr::na_if(LClabel, "NA")
        )

      dt_data <- authorship_data %>%
        dplyr::ungroup() %>%
        dplyr::group_by(protocol_ID, protocolVersion, bibliographicCitation) %>%
        dplyr::summarise(
          `# of Samples` = dplyr::n_distinct(samp_name, na.rm = TRUE),
          `# of Locations` = dplyr::n_distinct(station, na.rm = TRUE),
          Contact = paste(unique(ownerContact), collapse = "; "),
          # collapse LClabels with datasetID_obis prepended as links
          LClabel = {
            valid <- !is.na(LClabel)

            if (!any(valid)) {
              NA_character_
            } else {
              ul <- unique(LClabel[valid])
              out <- vapply(ul, function(lbl) {

                if (is.na(lbl) || !nzchar(lbl)) return(NA_character_)

                ids <- unique(datasetID_obis[valid & LClabel == lbl])
                ids <- ids[!is.na(ids)]

                if (length(ids) == 0) return(NA_character_)

                links <- paste0(
                  '<a href="https://obis.org/dataset/', ids,
                  '" target="_blank">', ids, '</a>'
                )

                paste0(
                  lbl,
                  '<div style="margin-top:10px; font-size: 0.9em;">',
                  '<strong>Data available on the Ocean Biodiversity Information System (OBIS):</strong><br>',
                  paste(links, collapse = ", "),
                  '</div>'
                )
              }, character(1))


              out <- out[!is.na(out)]

              if (length(out) == 0) NA_character_ else paste(out, collapse = "<br><br>")
            }
          }


        ) %>%
        dplyr::ungroup()


      # Helper to generate shiny action buttons in DT
      shinyInput <- function(FUN, len, id, ...) {
        vapply(seq_len(len), function(i) {
          as.character(FUN(paste0(id, "_", i), ...))
        }, character(1))
      }

      dt_data <- dt_data %>%
        dplyr::mutate(
          `Indigenous Contributions` = NA_character_  # default: no button
        )

      # Only assign buttons to rows where LClabel is not NULL
      rows_with_label <- which(!is.na(dt_data$LClabel) & nzchar(dt_data$LClabel))

      dt_data$`Indigenous Contributions`[rows_with_label] <- sprintf(
        '<button id="%s" style="border:0; background:transparent;" title="Indigenous Contributions">
           <img src="img/fn_logo.png" height="25"/>
         </button>',
        ns(paste0("fn_btn_", rows_with_label))
      )

      dt_data <- dt_data %>%
        dplyr::rename(
          "Protocol ID" = "protocol_ID",
          "Protocol Version" = "protocolVersion",
          "Publication" = "bibliographicCitation"
        ) |>
        dplyr::relocate(
          `Protocol ID`, `Protocol Version`, `# of Samples`, `# of Locations`, `Indigenous Contributions`, Publication
        )

      r$dt_data(dt_data)  # save the current table for observers

      DT::datatable(
        dt_data,
        escape = FALSE,
        rownames = FALSE,
        options = list(
          columnDefs = list(
            list(className = "dt-center", targets = "_all"),
            list(visible = FALSE, targets = which(names(dt_data) == "LClabel") - 1)
        ),
        drawCallback = JS(sprintf(
            "function(settings) {
        for (var i = 0; i < %d; i++) {
          var btn = document.getElementById('%s' + (i+1));
          if (btn) {
            btn.onclick = function(e) {
              Shiny.setInputValue('%s', this.id, {priority: 'event'});
            };
          }
        }
      }",
            nrow(dt_data),
            ns("fn_btn_"),
            ns("fn_btn_clicked")
          ))
        )
      )
    })

    observeEvent(input$fn_btn_clicked, {
      dt_data <- r$dt_data()
      req(dt_data)

      row_idx <- as.numeric(sub(
        paste0(session$ns("fn_btn_")),
        "",
        input$fn_btn_clicked
      ))

      text_to_show <- dt_data$LClabel[row_idx]

      showModal(
        modalDialog(
          title = "Indigenous Contributions",
          size = "xl",
          easyClose = TRUE,
          footer = modalButton("Close"),
          div(
            class = "indig_container",
            div(
              class = "indig_img",
              tags$img(
                src = "img/fn_logo.png",
                alt = "Indigenous icon"
              )
            ),
            div(
              class = "indig_txt",
              HTML(text_to_show)
            )
          ),
          class = "indig_modal"
        )
      )
    })


    # EXPORT PDF
    # output$export_pdf <- downloadHandler(
    #   filename = function() {
    #     paste0("GOTeDNA_report_", Sys.Date(), ".pdf")
    #   },
    #   contentType = "application/pdf",
    #   content = function(file) {
    #     # Data Request
    #     #   datSrc <- input$datasource
    #     #  datTyp <- input$data_type
    #     # primers <- r$primer

    #     # Area Selection
    #     if (!is.null(r$geom_slc)) {
    #       geom_coords <- sf::st_bbox(r$geom_slc) |>
    #         as.matrix() |>
    #         t() |>
    #         as.data.frame()
    #     } else {
    #       geom_coords <- c("Area selection not confirmed")
    #     }

    #     mapDL <- leaflet::addMarkers(
    #       basemap(),
    #       data = r$geom,
    #       clusterOptions = leaflet::markerClusterOptions(),
    #       label = ~ paste(success, "observations"),
    #       group = "station"
    #     ) |>
    #       leaflet::addScaleBar()

    #     mapview::mapviewOptions(fgb = FALSE)

    #     mapview::mapshot(
    #       mapDL,
    #       file = file.path(tempdir(), "mapDL.png")
    #     )

    #     # Observation
    #     thresh_slc <- input$threshold
    #     protID <- input$prot_id
    #     sampWin <- calc_window(
    #       threshold = thresh_slc,
    #       scaledprobs = r$scaledprobs
    #     )

    #     cons <- jaccard_test(
    #       r$scaledprobs,
    #       thresh_slc
    #     )


    #     # Reference Data Authorship
    #     FNdata <- r$data_ready$LClabel

    #     tempReport <- file.path(tempdir(), "report.Rmd")
    #     tempDFOlogo <- file.path(tempdir(), "DFOlogo.png")
    #     tempGOTlogo <- file.path(tempdir(), "GOTeDNAlogo.png")
    #     temphmLeg <- file.path(tempdir(), "hmLegend.png")
    #     tempthreshAx <- file.path(tempdir(), "threshAxis.png")
    #     tempthreshLeg <- file.path(tempdir(), "threshLegend.png")

    #     file.copy("Report.rmd", tempReport, overwrite = TRUE)
    #     file.copy("DFOlogo.png", tempDFOlogo, overwrite = TRUE)
    #     file.copy("GOTeDNAlogo.png", tempGOTlogo, overwrite = TRUE)
    #     file.copy("hm_legend.png", temphmLeg, overwrite = TRUE)
    #     file.copy("thresh_axis.png", tempthreshAx, overwrite = TRUE)
    #     file.copy("thresh_legend.png", tempthreshLeg, overwrite = TRUE)

    #     # params <- list(win, j.sim)

    #     out <- rmarkdown::render(tempReport) # ,
    #     # params = params)
    #     file.rename(out, file)
    #   }
    # )
  })
}


#------- INTERNALS

prepare_data <- function(r) {
  out <- r$cur_data
  if (length(r$station_slc)) {
    out <- out |>
      dplyr::filter(station %in% r$station_slc)
  }
  if (!is.null(r$frozen_selected_taxon_level)) {
    if (r$frozen_selected_taxon_level == "scientificName") {
      out <- out |>
        dplyr::filter(scientificName == r$scientificName)
    } else {
      if (r$frozen_selected_taxon_id != "All") {
        out <- out[
          out[[r$frozen_selected_taxon_level]] == r$frozen_selected_taxon_id,
        ]
      }
    }
  }
  # primer-based subset
  out |>
    dplyr::filter(primer %in% r$primer)
}

## FIG HELPERS

plotText <- function(txt, size = 6) {
  data.frame(x = 0.5, y = 0.5, txt = txt) |>
    ggplot2::ggplot(ggplot2::aes(x, y, label = txt)) +
    ggplot2::geom_text(size = size) +
    ggplot2::theme_void()
}

plotNotAvailable <- function() {
  plotText(stringr::str_wrap("Plot not available. Click 'Compute & visualize'"))
}

plotNotAvailableError <- function() {
  plotText(stringr::str_wrap(
    "An error occured while rendering the plot (sample size may be too small).",
    width = 60
  ))
}

## Detection part 1
draw_fig_smooth <- function(r, ready, id) {
  if (ready) {
    plt <- try(
      smooth_fig(
        r$scaledprobs
      )
    )
    if (inherits(plt, "try-error")) {
      plotNotAvailableError()
      # there are lots of erros due to predLoess, this is better than a misleading
      # message on the plot
    } else {
      plt
    }
  } else {
    plotNotAvailable()
  }
}

## Detection part 2
draw_fig_detect <- function(r, ready, threshold) {
  if (ready) {
    plt <- thresh_fig(
      threshold,
      r$scaledprobs
    )
    plt
  } else {
    plotNotAvailable()
  }
}


## Data variation
draw_fig_samples <- function(newprob, ready, year) {
  if (ready) {
    global_max_n <- max(unlist(lapply(newprob$newP_yr, `[[`, "n")), na.rm = TRUE)
    p <- field_sample_fig(
      newprob, year, global_max_n
    )
    p
  } else {
    plotNotAvailable()
  }
}

# Top plot: smooth figure
ui_fig_smooth <- function(fig_id, title, caption_files, ns) {
  div(
    id = paste0(ns(fig_id), "_fig_container"),
    class = "fig_container",
    style = "padding-bottom: 0rem;",
    h4(title),
    div(
      class = "fig_caption-container",
      div(
        class = "fig_caption",
        includeHTML(file.path("www", "doc", "caption", caption_files[[1]]))
      ),
      div(
        class = "fig_caption",
        includeHTML(file.path("www", "doc", "caption", caption_files[[2]]))
      )
    ),
    div(
      class = "fig_panel_container",
      style = "padding-bottom: 0rem;",
      div(
        class = "fig_panel",
        style = "padding-bottom: 0rem;",
        bslib::layout_columns(
          # Left legend
          bslib::card_image(
            file = "www/img/fixed-legends/thresh_axis.png",
            fill = FALSE,
            style = "max-width:120px; width:100%; height:auto;"
          ),
          # Middle plot
          bslib::card_body(
            plotOutput(ns("fig_smooth_plot_output"), width = "100%", height = "auto"),
            fillable = TRUE
          ),
          col_widths = c(2, 10)
        )
      )
    )
  )
}

# Bottom plot: detect figure
ui_fig_detect_bottom <- function(fig_id, ns) {
  div(
    id = paste0(ns(fig_id), "_fig_container"),
    class = "fig_container",
    style = "padding-top: 0rem;",
    div(
      class = "fig_panel_container",
      style = "padding-top: 0rem;",
      div(
        class = "fig_panel",
        style = "padding-top: 0rem;",
        bslib::layout_columns(
          # Left empty spacer
          div(),
          # Middle plot
          bslib::card_body(
            plotOutput(ns("fig_detect_plot_output"), width = "100%", height = "100%"),
            style = "min-height:450px;",
            fillable = TRUE
          ),
          # Right legend
          div(
            style = "display:flex; align-items:center; height:100%;",
            bslib::card_image(
              file = "www/img/fixed-legends/thresh_legend.png",
              fill = FALSE,
              style = "max-width:620px; width:100%; height:auto;"
            )
          ),
          col_widths = c(2, 7, 3)
        )
      )
    )
  )
}


ui_fig_hm <- function(fig_id, title, caption_file, ns) {
  div(
    id = paste0(ns(fig_id), "_fig_container"),
    class = "fig_container",
    h4(title),
    div(
      class = "fig_caption-container",
      div(
        class = "fig_caption",
        includeHTML(file.path("www", "doc", "caption", caption_file))
      )
    ),
    div(
      class = "fig_panel_container",
      div(
        class = "fig_panel",
        bslib::layout_columns(
          bslib::card_body(
            plotly::plotlyOutput(paste0(ns(fig_id), "_plot_output"),
              height = "auto"
            ),
            fillable = TRUE,
          ),
          bslib::card_image(
            file = "www/img/fixed-legends/hm_legend.png",
            fill = FALSE,
            style = "min-width:100px; max-width:200px;"
          ),
          col_widths = bslib::breakpoints(
            sm = c(9, 3),
            md = c(10, 2)
          )
        )
      )
    )
  )
}

ui_fig_effort <- function(fig_id, title, caption_file, ns) {
  div(
    id = paste0(ns(fig_id), "_fig_container"),
    class = "fig_container",
    style = "width:100%;",
    h4(title),
    div(
      class = "fig_caption-container",
      style = "width:100%;",
      div(
        class = "fig_caption",
        includeHTML(file.path("www", "doc", "caption", caption_file))
      )
    ),
    div(
      class = "fig_panel_container",
      div(
        class = "fig_panel",
        bslib::layout_columns(
          bslib::card_body(
            plotly::plotlyOutput(paste0(ns(fig_id), "_plot_output"),
                                 height = "auto"
            ),
            fillable = TRUE,
          ),
          bslib::card_image(
            file = "www/img/fixed-legends/effort_legend.png",
            fill = FALSE,
            style = "min-width:50px; max-width:100px;"
          ),
          col_widths = bslib::breakpoints(
            sm = c(9, 3),
            md = c(10, 2)
          )
        )
      )
    )
  )
}

ui_fig_samples <- function(fig_id, title, caption_file, ns) {
  div(
    id = paste0(ns(fig_id), "_fig_container"),
    class = "fig_container",

    h4(title),

    div(
      class = "fig_caption-container",
      div(
        class = "fig_caption",
        includeHTML(file.path("www", "doc", "caption", caption_file))
      )
    ),
    # somewhere near your plot output
    verbatimTextOutput(ns("debug_field_sample")),


    uiOutput(ns("fig_samples_controls")),

    div(
      class = "fig_panel_container",
      div(
        class = "fig_panel",
        plotly::plotlyOutput(
          paste0(ns(fig_id), "_plot_output"),
          height = "auto"
        )
      )
    )
  )
}


## Additional helpers

add_thumbnail_button <- function(id, src, title, alt = "Figure thumbnail") {
  # https://stackoverflow.com/questions/44841346/adding-an-image-to-shiny-action-button
  tagList(
    div(
      class = "thumnail_button_container",
      tags$button(
        id = id,
        class = "btn action-button",
        tags$img(
          src = src,
          alt = alt,
          id = "fig-thumbnail"
        ),
        title = title
      ),
      div(
        id = paste0(id, "_thumbnail_selected"),
        class = "thumbnail_selected",
        p("selected")
      )
    ),
  )
}

add_figure_selection <- function(id, title, scr = NULL, info = title) {
  column(
    2,
    div(
      class = "thumbnail",
      add_thumbnail_button(id, scr, info),
      h5(title)
    )
  )
}

show_fig <- function(fig_id) {
  shinyjs::show(paste0(fig_id, "_thumbnail_selected"))
  shinyjs::show(paste0(fig_id, "_fig_container"))
}

hide_fig <- function(fig_id) {
  shinyjs::hide(paste0(fig_id, "_thumbnail_selected"))
  shinyjs::hide(paste0(fig_id, "_fig_container"))
}

toggle_legends <- function(legend_id) {
  shinyjs::toggle(legend_id)
}

toggle_fig <- function(fig_id) {
  shinyjs::toggle(paste0(fig_id, "_thumbnail_selected"))
  shinyjs::toggle(paste0(fig_id, "_fig_container"))
}


























