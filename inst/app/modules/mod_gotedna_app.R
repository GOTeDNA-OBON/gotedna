# -------------------------
# GOTeDNA Module
# -------------------------
mod_gotedna_ui <- function(id) {
  ns <- NS(id)

  # Wrap the existing GOTeDNA UI
  div(
    id = "gotedna_app",

    navbarPage(id = "navbar",
               position = "fixed-top",
               img(
                 src = "img/logo/GOTeDNA_logo_white_got.svg",
                 alt = "GOTeDNA logo",
                 title = "GOTeDNA logo",
                 id = "logo_gotedna"
               ),
               tabPanel(
                 "Home",
                 mod_select_data_ui(ns("slc_data")),
                 mod_select_figure_ui(ns("slc_fig")),
                 tags$script(type = "text/javascript", src = "js/definitionEvents.js")
               ),
               tabPanel(
                 "Data Structure",
                 div(
                   class = "standalone_container",
                   div(
                     class = "standalone_60",
                     h1("Data Structure"),
                     includeHTML(file.path("www","doc","structure.html"))
                   )
                 )
               ),
               tabPanel(
                 "Interpretation Guide",
                 value = "interp-guide",
                 div(
                   class = "standalone_container",
                   div(
                     class = "standalone_60",
                     h1("Interpretation Guide"),
                     includeHTML(file.path("www", "doc", "interp_guide.html"))
                   )
                 )
               ),
               tabPanel(
                 title = "Primers",
                 value = "primer-info",
                 div(
                   class = "standalone_container",
                   div(
                     class = "standalone_80",
                     mod_primers_ui(ns("primer_seq"))
                   )
                 )
               ),
               tabPanel(
                 "Contact",
                 value = "contact",
                 div(
                   class = "standalone_container",
                   div(
                     class = "standalone_60",
                     h1("Contact"),
                     includeHTML(file.path("www", "doc", "contact.html"))
                   )
                 )
               ),
               tabPanel(
                 "Disclaimers",
                 value = "disc",
                 div(
                   class = "standalone_container",
                   div(
                     class = "standalone_60",
                     h1("Disclaimers"),
                     includeHTML(file.path("www", "doc", "disclaimer.html"))
                   )
                 )
               ),
               tabPanel(
                 "Glossary",
                 value = "glossary",
                 div(
                   class = "standalone_container",
                   div(
                     class = "standalone_80",
                     mod_glossary_ui(ns("glossary"))
                   )
                 )
               )
    )
  )
}


mod_gotedna_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {

    # Exact same server modules as current GOTeDNA
    mod_select_data_server("slc_data", r)

    mod_dialog_disclaimers_server("show_dialog", r)
    observeEvent(input$show_dialog, r$show_dialog <- TRUE)
    observeEvent(input$show_help, r$show_help <- TRUE)
    mod_dialog_map_info_server("show_map_info", r)
    mod_glossary_server("glossary")

    mod_primers_server("primer_seq")
    observeEvent(input$show_source, r$show_source <- TRUE)

    observeEvent(input$reset, {
      shinyjs::reset("data_authorship")
    })

    mod_select_figure_server("slc_fig", r)
  })
}
