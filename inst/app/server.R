server <- function(input, output, session) {

  # -----------------------------
  # Existing reactive values
  # -----------------------------
  r <- reactiveValues(
    geom = NULL,
    geom_slc = NULL,
    station_slc = NULL,
    protocol_slc = NULL,
    protocol_info = default_protocol_info,
    protocol_summary = NULL,
    upload_data = NULL,
    upload_stations = NULL,
    data_selected = NULL,
    upload_primers = NULL,
    active_primers = NULL,
    confirmed_prot = NULL,
    confirmed_thresh = NULL,
    scientificName = NULL,
    taxon_lvl_slc = NULL,
    taxon_selected = 0,
    dt_data = reactiveVal(NULL),
    taxon_id_slc = NULL,
    show_map_info = FALSE,
    reload_map = 0,
    fig_ready = FALSE,
    fig_slc = list(
      fig_heatmap = FALSE,
      fig_effort = FALSE,
      fig_samples = FALSE,
      fig_detect = FALSE,
      fig_smooth = FALSE
    ),
    current_fig = "fig1",
    lock_view = FALSE,
    reset = 0
  )

  mod_dialog_disclaimers_server("show_dialog", r)
  observe({
    if (is.null(app_choice()) && is.null(r$disclaimer_shown)) {
      r$show_dialog <- TRUE
      r$disclaimer_shown <- TRUE
    }
  })
  # -----------------------------
  # 1️⃣ App choice module
  # -----------------------------
  app_choice <- mod_app_choice_server("app_choice")

  # -----------------------------
  # 2️⃣ Render choose page
  # -----------------------------
  output$choose_ui <- renderUI({
    choice <- app_choice()
    if (is.null(choice)) {
      mod_app_choice_ui("app_choice")
    } else {
      NULL
    }
  })

  app_b_started <- reactiveVal(TRUE)

  output$app_b_ui <- renderUI({
      app_b_ui()
  })

  outputOptions(output, "app_b_ui", suspendWhenHidden = FALSE)

  observe({
    choice <- app_choice()

    if (!is.null(choice) && choice == "A") {
      shinyjs::show("gotedna_app")
      shinyjs::hide("new_app")

      mod_select_data_server("slc_data", r)
      # mod_dialog_disclaimers_server("show_dialog", r)
      observeEvent(input$show_dialog, r$show_dialog <- TRUE)
      observeEvent(input$show_help, r$show_help <- TRUE)
      mod_dialog_map_info_server("show_map_info", r)
      mod_glossary_server("glossary")
      mod_primers_server("primer_seq")
      observeEvent(input$show_source, r$show_source <- TRUE)
      observeEvent(input$reset, { shinyjs::reset("data_authorship") })
      mod_select_figure_server("slc_fig", r)

    } else if (!is.null(choice) && choice == "B") {
      shinyjs::hide("choice_page")
      shinyjs::hide("gotedna_app")
      shinyjs::hide("app_b_loading")

      shinyjs::show("new_app")
      shinyjs::show("app_b_footer")
    }
  })


  # Make sure shinyjs::useShinyjs() is in your UI
  shinyjs::onclick("logo_gotedna", {
    app_choice(NULL)            # reset choice to show choice page
    shinyjs::hide("gotedna_app")
    shinyjs::hide("new_app")
  })

  shinyjs::onclick("logo_mpa", {
    app_choice(NULL)
    shinyjs::show("choice_page")
    shinyjs::hide("gotedna_app")
    shinyjs::hide("new_app")
    shinyjs::hide("app_b_loading")
  })
}

