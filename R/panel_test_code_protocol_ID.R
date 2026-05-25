#Tester code for Protocol ID and target gene/primer panel

library(shiny)
library(leaflet)

ui <- fluidPage(

  tags$head(

    # ----- CSS -----
    tags$style(HTML("

      html, body {
        height: 100%;
        margin: 0;
        padding: 0;
        overflow: hidden;
      }

      .container-fluid {
        padding: 0;
        height: 100%;
      }

      #map_wrap {
        position: relative;
        width: 100%;
        height: 100vh;
      }

      #map {
        width: 100%;
        height: 100%;
      }

      /* ----- Sliding panel ----- */

      #side_panel {
        position: absolute;
        top: 0;
        left: calc(-25% + 42px);
        width: 25%;
        height: 100%;
        background: white;
        z-index: 1000;
        transition: left 0.3s ease;
        box-shadow: 2px 0 10px rgba(0,0,0,0.3);
        overflow-y: auto;
        padding: 20px;
      }

      #side_panel.open {
        left: 0;
      }

      /* ----- Clickable tab ----- */

      #panel_tab {
        position: absolute;
        top: 100px;
        right: -42px;

        width: 42px;
        height: 140px;

        background: white;

        border: 1px solid #ccc;
        border-left: none;

        cursor: pointer;

        display: flex;
        align-items: center;
        justify-content: center;

        writing-mode: vertical-rl;
        text-orientation: mixed;

        font-weight: bold;
        font-size: 16px;

        box-shadow: 2px 2px 6px rgba(0,0,0,0.2);

        z-index: 1001;
      }

    ")),

    # ----- Javascript -----
    tags$script(HTML("

      $(document).on('click', '#panel_tab', function() {
        $('#side_panel').toggleClass('open');
      });

    "))
  ),

  # ----- Map Wrapper -----
  div(
    id = "map_wrap",

    leafletOutput("map"),

    # ----- Side Panel -----
    div(
      id = "side_panel",

      # Tab
      div(
        id = "panel_tab",
        "FILTERS"
      ),

      h2("Map Controls"),

      br(),

      selectInput(
        "year",
        "Select Year",
        choices = c("All", 2021, 2022, 2023)
      ),

      checkboxInput(
        "show_points",
        "Show Sampling Points",
        TRUE
      ),

      selectInput(
        "group",
        "Taxonomic Group",
        choices = c(
          "All",
          "Fish",
          "Mammals",
          "Birds",
          "Invertebrates"
        )
      ),

      sliderInput(
        "depth",
        "Depth Range",
        min = 0,
        max = 500,
        value = c(0, 500)
      ),

      br(),

      actionButton(
        "confirm",
        "Confirm Filters"
      )
    )
  )
)

server <- function(input, output, session) {

  # ----- Base map -----
  output$map <- renderLeaflet({

    leaflet() %>%

      addProviderTiles(
        providers$CartoDB.Positron
      ) %>%

      setView(
        lng = -63.6,
        lat = 44.7,
        zoom = 6
      )
  })

  # ----- Example points -----
  observeEvent(input$confirm, {

    proxy <- leafletProxy("map")

    proxy %>%
      clearMarkers()

    if (isTRUE(input$show_points)) {

      proxy %>%

        addCircleMarkers(
          lng = c(-63.6, -64.2, -62.9),
          lat = c(44.7, 45.1, 44.3),

          radius = 8,

          fillOpacity = 0.9,

          label = c(
            "Site 1",
            "Site 2",
            "Site 3"
          )
        )
    }
  })
}

shinyApp(ui, server)
