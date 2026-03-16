mod_app_choice_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "standalone_container",
    div(
      class = "standalone_60",
      h1("Choose Application"),
      p("Please select which tool you want to launch."),

      actionButton(ns("app_a"), "Explore eDNA Detection Rates"),
      br(), br(),
      actionButton(ns("app_b"), "Explore Marine Protected Areas")
    )
  )
}

mod_app_choice_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    choice <- reactiveVal(NULL)

    observeEvent(input$app_a, { choice("A") })
    observeEvent(input$app_b, { choice("B") })

    choice
  })
}
