ui <- fluidPage(
  theme = bslib::bs_theme(version = 5),
  shinyjs::useShinyjs(),

  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = paste0("app_a_styles.css?v=", Sys.time())),
    tags$link(rel = "stylesheet", type = "text/css", href = "fonts.css"),
    tags$link(rel = "stylesheet", type = "text/css", href = paste0("choice_page_styles.css?v=", Sys.time())),
    tags$script(type = "text/javascript", src = "js/scrollPage.js"),
    tags$script(type = "text/javascript", src = "js/fakeClick.js")
  ),

  # =========================
  # CHOICE PAGE ONLY
  # =========================

  div(
    id = "choice_page",

    uiOutput("choose_ui"),

    div(
      id = "choice_footer",

      fluidRow(
        class = "align-items-center",

        column(
          3,
          a(
            img(
              title = "Fisheries and Oceans Canada",
              src = "img/logo_partners/DFO_logo_sq.svg",
              alt = "DFO Logo",
              id = "logo_dfo"
            ),
            href = "https://www.dfo-mpo.gc.ca/index-eng.html",
            target = "_blank"
          )
        ),

        column(
          3,
          a(
            img(
              title = "Maine-eDNA",
              src = "img/logo_partners/logo_Maine_eDNA_nbg_w.png",
              alt = "Maine eDNA Logo",
              id = "logo_mswc"
            ),
            href = "https://umaine.edu/edna/",
            target = "_blank"
          )
        ),

        column(
          3,
          a(
            img(
              title = "UN Ocean Decade",
              src = "img/logo_partners/logo_undossd.svg",
              alt = "UN Logo",
              id = "logo_undossd"
            ),
            href = "https://oceandecade.org/",
            target = "_blank"
          )
        ),

        column(
          3,
          a(
            img(
              title = "OBON",
              src = "img/logo_partners/logo_obon.svg",
              alt = "OBON Logo",
              id = "logo_obon"
            ),
            href = "https://obon-ocean.org/",
            target = "_blank"
          )
        )
      )
    )
  ),

  # =========================
  # APP A - Original GOTeDNA App
  # =========================
  tags$button(
    id = "scroll-top",
    type = "button",
    class = "btn btn-default",
    onclick = "topFunction()",
    style = "display:none; position:fixed; bottom:40px; right:40px; z-index:1000;",
    "^ Top"
  ),

  div(
    id = "gotedna_app",
    class = "app_A",
    style = "display:none;",

    navbarPage(id = "navbar",
               position = "fixed-top",
               # Wrap the logo in an <a> for clickable action
               tags$a(
                 href = "#",
                 id = "logo_gotedna",
                 img(
                   src = "img/logo/GOTeDNA_logo_white_got.svg",
                   alt = "GOTeDNA logo",
                   title = "Back to App Selection",
                   style = "height:40px;"   # adjust height as needed
                 )
               ),
               tabPanel(
                 "Home",
                 mod_select_data_ui("slc_data"),
                 mod_select_figure_ui("slc_fig"),
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
                 "Dataflow",
                 value = "interp-guide",
                 div(
                   class = "standalone_container",
                   div(
                     class = "standalone_60",
                     h1("Dataflow"),
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
                     mod_primers_ui("primer_seq")
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
                     mod_glossary_ui("glossary")
                   )
                 )
               )
    )
  ),

  # =========================
  # APP B - GOTeDNA-MPA
  # =========================
  div(
    id = "app_b_loading",
    style = "display:none; padding-top:120px; text-align:center;",
    h3("Loading GOTeDNA-MPA..."),
    p("Preparing map and analysis modules.")
  ),

  div(
    id = "new_app",
    class = "app_B",
    style = "display:none;",

    uiOutput("app_b_ui"),

  div(
    id = "app_b_footer",
    fluidRow(
      class = "align-items-center",

      column(
        3,
        a(
          img(
            title = "Fisheries and Oceans Canada",
            src = "img/logo_partners/DFO_logo_sq.svg",
            alt = "DFO Logo",
            id = "logo_dfo"
          ),
          href = "https://www.dfo-mpo.gc.ca/index-eng.html",
          target = "_blank"
        )
      ),

      column(
        3,
        a(
          img(
            title = "Maine-eDNA",
            src = "img/logo_partners/logo_Maine_eDNA_nbg_w.png",
            alt = "Maine eDNA Logo",
            id = "logo_mswc"
          ),
          href = "https://umaine.edu/edna/",
          target = "_blank"
        )
      ),

      column(
        3,
        a(
          img(
            title = "UN Ocean Decade",
            src = "img/logo_partners/logo_undossd.svg",
            alt = "UN Logo",
            id = "logo_undossd"
          ),
          href = "https://oceandecade.org/",
          target = "_blank"
        )
      ),

      column(
        3,
        a(
          img(
            title = "OBON",
            src = "img/logo_partners/logo_obon.svg",
            alt = "OBON Logo",
            id = "logo_obon"
          ),
          href = "https://obon-ocean.org/",
          target = "_blank"
        )
       )
      )
    )
  )
)
