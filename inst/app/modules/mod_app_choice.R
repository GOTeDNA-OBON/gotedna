mod_app_choice_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "standalone_container choice-container",
    div(
      class = "standalone_60",
      h1("Choose GOTeDNA Application"),
      p("Please select which tool you want to launch."),

      div(class = "choice-btn-row",
          actionButton(ns("app_a"), "Explore eDNA Detection Rates", class = "choice-btn"),
          actionButton(ns("app_b"), "Explore Marine Protected Areas", class = "choice-btn")
      )
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
