# -------------------------
# New App Placeholder Module
# -------------------------

mod_new_app_ui <- function(id) {
  ns <- NS(id)

  div(
    id = "new_app",
    navbarPage(id = ns("navbar"),
               position = "fixed-top",
               img(
                 src = "img/logo/GOTeDNA_logo_white_got.svg",  # temporary logo
                 alt = "New App logo",
                 title = "New App logo",
                 id = "logo_new_app"
               ),
               tabPanel(
                 "Home",
                 div(
                   class = "standalone_container",
                   div(
                     class = "standalone_60",
                     h1("New App Placeholder"),
                     p("This will become the new application with its own navbar and modules.")
                   )
                 )
               ),
               tabPanel(
                 "Documentation",
                 div(
                   class = "standalone_container",
                   div(
                     class = "standalone_60",
                     h1("Documentation Placeholder"),
                     p("Include HTML docs or instructions here.")
                   )
                 )
               )
    )
  )
}

mod_new_app_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # placeholder server logic
    # You can add reactiveValues or modules here later
  })
}
