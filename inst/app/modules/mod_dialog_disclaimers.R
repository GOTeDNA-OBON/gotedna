mod_dialog_disclaimers_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    query_modal <- modalDialog(
      title = NULL,
      footer = NULL,
      class = "disclaimer_modal",
      tagList(
        fluidRow(
          column(
            6,
            div(
              id = "disclaimer_left",
              img(
                src = "img/logo/logo_disclaimer.svg",
                alt = "Logo GOTeDNA",
                id = "logo_disclaimer"
              ),
              includeHTML(file.path("www","doc","welcome.html"))
            )
          ),
          column(
            6,
            div(
              id = "disclaimer_right",
              includeHTML(file.path("www","doc", "disclaimer.html")),
              fluidRow(
                column(
                  10,
                  checkboxInput(ns("agreed"),
                                "I have read and understood the above disclaimers.",
                                value = FALSE, width = "90%")
                ),
                column(2, actionButton(ns("dismiss"), "OK", class = "btn-blue"))
              )
            )
          )
        )
      ),
      easyClose = FALSE,
      size = "xl"
    )

    # Show modal when r$show_dialog is TRUE
    observeEvent(r$show_dialog, {
      if (isTRUE(r$show_dialog)) {
        showModal(query_modal)
        # Prevent body scrolling while modal is open
        shinyjs::runjs("document.body.style.overflow = 'hidden';")
      }
    })

    # Handle dismiss button
    observeEvent(input$dismiss, {
      if (isTRUE(input$agreed)) {
        removeModal()
        r$show_dialog <- FALSE

        # ✅ Restore scrolling after modal closes
        shinyjs::runjs("
          document.body.style.overflow = 'auto';
          window.dispatchEvent(new Event('resize'));
        ")
      } else {
        showNotification("You must check the box!", type = "warning")
      }
    })
  })
}
