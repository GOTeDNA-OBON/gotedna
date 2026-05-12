# -------------------------
# New App Placeholder Module
# -------------------------

app_b_ui <- function() {

  fluidPage(
    shinyjs::useShinyjs(),
    theme = bslib::bs_theme(version = 5),

    tags$head(

      tags$link(
        href = "https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap",
        rel = "stylesheet"
      ),

      # jQuery UI (needed for draggable/resizable)
      tags$link(
        rel  = "stylesheet",
        href = "https://code.jquery.com/ui/1.13.2/themes/base/jquery-ui.css"
      ),
      tags$script(src = "https://code.jquery.com/ui/1.13.2/jquery-ui.min.js"),

      tags$link(
        rel = "stylesheet",
        type = "text/css",
        href = paste0("app_b_styles.css?v=", Sys.time())
      ),

      tags$script(HTML("
$(function(){

  // Smooth scroll navbar
  $(document).on('click', 'a.nav-scroll', function(e){
    e.preventDefault();
    var target = $(this).data('target');
    var el = document.getElementById(target);
    if(el){
      el.scrollIntoView({behavior:'smooth', block:'start'});
    }
  });

    // Ensure Bootstrap 3 collapse is initialized on our target (prevents silent no-op)
  function ensureCollapseInit(){
    var $body = $('#floating_body');
    if(!$body.length) return;

    if (typeof $body.collapse === 'function') {
      // Initialize plugin without toggling
      $body.collapse({ toggle: false });
    }
  }

  function clampFloatingToMap(){
    var $p = $('#floating_panel');
    var $wrap = $('#map_wrap');
    if(!$p.length || !$wrap.length) return;

    var PAD = 14;

    var wrapW = $wrap.innerWidth();
    var wrapH = $wrap.innerHeight();

    $p.css({
      'max-width':  (wrapW - PAD*2) + 'px',
      'max-height': (wrapH - PAD*2) + 'px'
    });

    var pos = $p.position();
    var pW  = $p.outerWidth();
    var pH  = $p.outerHeight();

    var left = pos.left;
    var top  = pos.top;

    var minL = PAD;
    var minT = PAD;
    var maxL = wrapW - pW - PAD;
    var maxT = wrapH - pH - PAD;

    if(maxL < minL) maxL = minL;
    if(maxT < minT) maxT = minT;

    left = Math.min(Math.max(left, minL), maxL);
    top  = Math.min(Math.max(top,  minT), maxT);

    $p.css({ left: left + 'px', top: top + 'px' });
  }

  function setFloatingCollapsedUI(isCollapsed){
    var $p = $('#floating_panel');
    var $toggle = $('#floating_toggle');
    if(!$p.length || !$toggle.length) return;

    if(isCollapsed){
      if($('#floating_body').hasClass('in')){
        $p.data('open_h', $p.outerHeight());
      }

      $p.addClass('is-collapsed');

      var headerH = Math.max($toggle.outerHeight(true) || 0, 45) + 1;
      if($p.hasClass('ui-resizable')){
        $p.resizable('option', 'minHeight', headerH);
      }
      $p.css({ height: headerH + 'px' });

    } else {
      $p.removeClass('is-collapsed');

      if($p.hasClass('ui-resizable')){
        $p.resizable('option', 'minHeight', 180);
      }

      // If we've never been opened before, pick a sane default height
      var oh = $p.data('open_h');
      if(!oh){
        var $wrap = $('#map_wrap');
        var wrapH = $wrap.length ? $wrap.innerHeight() : 700;
        var isLaptop = window.innerWidth <= 1440 || window.innerHeight <= 900;
        oh = isLaptop ? wrapH * 0.72 : wrapH * 0.99;
        $p.data('open_h', oh);
      }
      $p.css({ height: oh + 'px' });
    }

    setTimeout(clampFloatingToMap, 0);
  }

  function enableFloatingResize(){
    var $p = $('#floating_panel');
    if(!$p.length) return;

    // IMPORTANT: we control draggable here. So set absolutePanel(draggable=FALSE).
    if($p.hasClass('ui-draggable')){
      $p.draggable('option', 'containment', '#map_wrap');
      $p.draggable('option', 'handle', '#floating_toggle');
    } else {
      $p.draggable({ containment: '#map_wrap', handle: '#floating_toggle' });
    }

    if(!$p.hasClass('ui-resizable')){
      $p.resizable({
        handles: 'e,s,se',
        minWidth: 320,
        minHeight: 200,
        containment: '#map_wrap'
      });
    } else {
      $p.resizable('option', 'containment', '#map_wrap');
    }

    $p.off('dragstop.clamp').on('dragstop.clamp', clampFloatingToMap);
    $p.off('resizestop.clamp').on('resizestop.clamp', clampFloatingToMap);

    $(window).off('scroll.clampFloating').on('scroll.clampFloating', clampFloatingToMap);
    $(window).off('resize.clampFloating').on('resize.clampFloating', clampFloatingToMap);

    setTimeout(clampFloatingToMap, 0);
    setTimeout(clampFloatingToMap, 250);
  }

  function initFloatingPanel(){
    ensureCollapseInit();
    enableFloatingResize();
    setFloatingCollapsedUI(true); // start collapsed
  }

  $(document).on('shiny:connected', initFloatingPanel);
  initFloatingPanel();

  // Bootstrap 3 collapse events
  $(document).on('hide.bs.collapse', '#floating_body', function(){
    setFloatingCollapsedUI(true);
  });

  $(document).on('hidden.bs.collapse', '#floating_body', function(){
    enableFloatingResize();
    setFloatingCollapsedUI(true);
    clampFloatingToMap();
  });

  $(document).on('shown.bs.collapse', '#floating_body', function(){
    enableFloatingResize();
    setFloatingCollapsedUI(false);
    clampFloatingToMap();
  });

  // Custom message to open the panel
  Shiny.addCustomMessageHandler('openFloating', function(message){
    var $p    = $('#floating_panel');
    var $body = $('#floating_body');
    var $tog  = $('#floating_toggle');
    if(!$p.length || !$body.length || !$tog.length) return;

     ensureCollapseInit();  // <<< add this too (safe)

    if($body.hasClass('in') || $body.hasClass('show')){
      $tog.attr('aria-expanded','true');
      setFloatingCollapsedUI(false);
      enableFloatingResize();
      clampFloatingToMap();
      return;
    }

    $body.one('shown.bs.collapse.openFloating', function(){
      $tog.attr('aria-expanded','true');
      requestAnimationFrame(function(){
        enableFloatingResize();
        setFloatingCollapsedUI(false);
        clampFloatingToMap();
      });
    });

    if(typeof $body.collapse === 'function'){
      $body.collapse('show');
    } else {
      $body.addClass('in').css('display','block');
      $body.trigger('shown.bs.collapse');
    }
  });
});

    "))
    ),

    div(
      class = "app_b",

      # ---- REAL NAVBAR ----
      tags$nav(
        class = "navbar navbar-light fixed-top",
        tags$div(
          class = "container-fluid",
          tags$a(
            class = "navbar-brand nav-scroll",
            href = "#",
            id = "logo_mpa",       # ID for shinyjs click
            img(
              src = "img/logo/GOTeDNA_logo_white_got.svg",
              alt = "GOTeDNA-MPA",
              title = "Back to App Selection",
              style = "height:40px;"  # adjust height as needed
            )
          ),
          tags$ul(
            class = "nav navbar-nav",
            tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_map",    "Interactive Map")),
            tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_sara",   "Detection Details")),
            tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_method", "Explore Protocols")),
            tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_datsel", "Data Selection and Download")),
            tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_div",    "Diversity Metrics")),
            tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_pie",    "Taxonomic Pie Chart"))
            #tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_dwnld",    "Download Data File")),
            #tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_refdat",    "Reference Data Authorship")),
          )
        )
      ),

      # ---- MAP SECTION ----
      div(
        id = "sec_map",
        class = "scroll-section",

        div(
          id = "map_wrap",

          leafletOutput(
            "map",
            height = "calc(100vh - 120px)"
          ),

          div(
            id = "map_loading_overlay",
            div(class = "loader")
          ),

          div(
            id = "monthly_plot_control",
            class = "leaflet-control",
            div(id = "monthly_plot_title", "Monthly Number of Samples Collected"),
            div(id = "monthly_plot_subtitle", textOutput("monthly_plot_subtitle", inline = TRUE)),
            plotOutput("monthly_circular_plot", height = "250px", width = "100%")
          ),

          absolutePanel(
            id = "floating_panel",
            fixed = FALSE, draggable = FALSE,
            top = 10, left = 70, width = 360,

            tags$button(
              id = "floating_toggle",
              type = "button",
              class = "btn btn-default btn-secondary",
              `data-bs-toggle` = "collapse",
              `data-bs-target` = "#floating_body",
              `aria-expanded` = "false",
              `aria-controls` = "floating_body",
              tagList(
                tags$span("Select a Site"),
                tags$span(class = "caret-icon", HTML("&#9662;"))
              )
            ),

            div(
              id = "floating_body",
              class = "collapse",
              div(
                class = "panel-body",
                h4("Filter"),
                selectInput(
                  "sel_year",
                  "Year",
                  choices = c("All"), #use this code if you don't want to be able to select more than one year at a time
                  selected = "All"),

                #choices  = "All",  #use this code to be able to select more than one year at a time
                #selected = "All",
                #multiple = TRUE
                #),

                #h4("Group"),

                # --- Row 1: 4 across ---
                div(
                  class = "filter-btn-grid-4",
                  actionButton(
                    "total_fish",
                    label = tags$img(
                      src = "img/species_buttons/fish_centred.png",
                      alt = "Fish",
                      style = "height:40px;"  # adjust as needed
                    ),
                    title = "Fish",
                    class = "btn btn-default btn-secondary filter-btn"
                  ),
                  actionButton(
                    "total_sharks",
                    label = tags$img(
                      src = "img/species_buttons/shark.png",
                      alt = "Sharks & Rays",
                      style = "height:40px;"  # adjust as needed
                    ),
                    title = "Sharks & Rays",
                    class = "btn btn-default btn-secondary filter-btn"
                  ),
                  actionButton(
                    "total_mammals",
                    label = tags$img(
                      src = "img/species_buttons/whale2.png",
                      alt = "Mammals",
                      style = "height:40px;"  # adjust as needed
                    ),
                    title = "Mammals",
                    class = "btn btn-default btn-secondary filter-btn"
                  ),
                  actionButton(
                    "total_reptiles",
                    label = tags$img(
                      src = "img/species_buttons/turtle.png",
                      alt = "Turtles",
                      style = "height:40px;"  # adjust as needed
                    ),
                    title = "Turtles",
                    class = "btn btn-default btn-secondary filter-btn"
                  )
                ),

                tags$div(style="height:8px;"),  # optional spacing between rows

                # --- Row 2: 4 across ---
                div(
                  class = "filter-btn-grid-4",
                  actionButton(
                    "total_birds",
                    label = tags$img(
                      src = "img/species_buttons/bird.png",
                      alt = "Birds",
                      style = "height:40px;"  # adjust as needed
                    ),
                    title = "Birds",
                    class = "btn btn-default btn-secondary filter-btn"
                  ),
                  actionButton(
                    "total_molluscs",
                    label = tags$img(
                      src = "img/species_buttons/mollusc.png",
                      alt = "Molluscs",
                      style = "height:40px;"  # adjust as needed
                    ),
                    title = "Molluscs",
                    class = "btn btn-default btn-secondary filter-btn"
                  ),
                  actionButton(
                    "total_arthropods",
                    label = tags$img(
                      src = "img/species_buttons/lobster.png",
                      alt = "Arthropods",
                      style = "height:40px;"  # adjust as needed
                    ),
                    title = "Arthropods",
                    class = "btn btn-default btn-secondary filter-btn"
                  ),
                  actionButton(
                    "total_plants",
                    label = tags$img(
                      src = "img/species_buttons/plant.png",
                      alt = "Plants & Algae",
                      style = "height:40px;"  # adjust as needed
                    ),
                    title = "Plants & Algae",
                    class = "btn btn-default btn-secondary filter-btn"
                  ),
                ),

                tags$div(style="height:8px;"),

                # --- Row 3: 2 across ---
                div(
                  class = "filter-btn-grid-3",
                  actionButton("IUCN", "IUCN", title="IUCN Red List", class = "btn btn-default btn-secondary filter-btn filter-btn-short"),
                  actionButton("SARA", "SARA", title="Species at Risk", class = "btn btn-default btn-secondary filter-btn filter-btn-short"),
                  actionButton("AIS",  "AIS",  title="Aquatic Invasive Species", class = "btn btn-default btn-secondary filter-btn filter-btn-short")
                ),

                hr(),
                h4("Species List"),
                uiOutput("species_panel")
              )
            )
          )
        )
      ),

      # ---- DETECTION TABLE SECTION ----
      div(
        id = "sec_sara", class = "scroll-section",
        tabsetPanel(
          tabPanel("Detection Details",                                DT::DTOutput("detections_tbl")),
          tabPanel("IUCN Red List Details",                             DT::DTOutput("iucn_details")),
          tabPanel("Species at Risk Act (SARA): Schedule 1-3 Details", DT::DTOutput("sara_details")),
          tabPanel("Aquatic Invasive Species (AIS) Details",           DT::DTOutput("ais_details"))
        )
      ),

      # ---- Explore Protocols SECTION ----
      div(
        id = "sec_method", class = "scroll-section",
        h3("Explore Protocols"),

        # ---- Row 1: target gene + primer menus ----
        div(
          class = "data-select-grid",

          div(
            class = "data-select-item",
            selectizeInput(
              "div_target_gene",
              "Target gene",
              choices = NULL,
              selected = NULL,
              multiple = FALSE,
              options = list(
                placeholder = "Select target gene(s)"
              )
            )
          ),

          div(
            class = "data-select-item",
            selectizeInput(
              "div_primer",
              "Primer",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list(
                plugins = list("remove_button"),
                placeholder = "Select primer(s)"
              )
            ),
            div(
              class = "primer-btn-row",
              actionButton("div_primer_all", "Select all", class = "btn btn-default btn-secondary btn-sm"),
              actionButton("div_primer_none", "Deselect all", class = "btn btn-default btn-secondary btn-sm")
            )
          )
        ),

        # ---- Row 2: protocol dropdown + cards/plots ----
        div(
          id = "data_request_wrap",
          class = "method-comparison-wrap",
          h3(" "),

          # h4("Data Request"),

          fluidRow(
            column(
              width = 2,
              selectInput(
                "req_protocol",
                "Protocol ID",
                choices = NULL,
                selected = NULL,
                selectize = TRUE
              )
            ),

            column(
              width = 4,
              uiOutput("protocol_details"),
            ),

            column(
              width = 6,
              plotly::plotlyOutput("protocol_nmds_plot", height = "400px"),
              plotly::plotlyOutput("protocol_barplot", height = "350px")
            )
          )
        )
      ),

      # ---- DATA SELECTION SECTION ----
      div(
        id = "sec_datsel", class = "scroll-section",
        h3("Data Selection and Download"),

        div(
          class = "data-select-grid",

          div(
            class = "data-select-item",
            shinyWidgets::pickerInput(
              inputId = "prot_id",
              label = "Select Protocol IDs",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = shinyWidgets::pickerOptions(
                actionsBox = TRUE,
                liveSearch = TRUE,
                noneSelectedText = "Select protocol ID(s)"
              )
            )
          ),

          div(
            class = "data-select-item",
            selectizeInput(
              "div_compare_polygons",
              "Polygon Comparison",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list(
                plugins = list("remove_button"),
                placeholder = "Select polygon(s)"
              )
            ),
            div(
              style = "margin-top: 6px; font-size: 13px; color: #666;",
              "Used for diversity plots and statistics after clicking Confirm."
            )
          ),

          div(
            class = "data-select-item rarefaction-selection",
            numericInput(
              "rarefaction_depth",
              "Rarefaction depth",
              value = 5000,
              min = 0,
              max = 1000000000,
              step = 100
            )
          ),

          div(
            class = "data-select-item confirm-slot",
            div(
              class = "confirm-btn-row",
              actionButton(
                "div_apply",
                "Confirm",
                class = "btn btn-primary"
              ),
              shinyjs::disabled(
                downloadButton(
                  "downloadData",
                  "Download Data",
                  class = "btn-download-got"
                )
              )
            ),
            div(
              class = "download-hint-wrap",
              textOutput("download_hint", inline = TRUE)
            )
          )
        )
      ),

      # ---- DIVERSITY METRICS SECTION ----
      div(
        id = "sec_div", class = "scroll-section",
        h3("Diversity Metrics"),

        # ---- Top controls ----
        div(
          class = "data-select-grid",

          div(
            class = "data-select-item",
            selectInput(
              "tax_rank",
              "Taxonomy selection",
              choices = c(
                "Kingdom" = "kingdom",
                "Phylum"  = "phylum",
                "Class"   = "class",
                "Order"   = "order",
                "Family"  = "family",
                "Genus"   = "genus",
                "Species" = "scientificName"
              ),
              selected = "scientificName"
            )
          ),

          div(
            class = "data-select-item",
            selectInput(
              "alpha_metric",
              "Alpha Diversity",
              choices = c(
                "Observed Richness" = "observed",
                "Shannon" = "shannon",
                "Simpson" = "simpson",
                "Inverse Simpson" = "invsimpson",
                "ACE" = "ace",
                "Pielou's Evenness" = "pielou"
              ),
              selected = "observed"
            )
          ),

          div(class = "data-select-item"),
          div(class = "data-select-item")
        ),

        # ---- Alpha plot ----
        fluidRow(
          column(
            width = 8,
            offset = 2,
            div(
              class = "custom-loader-wrap",

              plotly::plotlyOutput("alpha_boxplot", height = "700px"),

              div(
                id = "alpha_loading_overlay",
                class = "custom-loading-overlay custom-placeholder-overlay",
                div(
                  class = "custom-placeholder-text",
                  "Select the Confirm button to view alpha diversity"
                )
              )
            )
          )
        ),

        # ---- Alpha warning ----
        fluidRow(
          column(
            width = 8,
            offset = 2,
            tags$div(
              style = "margin-top: 10px; font-size: 15px; color: #8a6d6d;",
              textOutput("alpha_warning_text", inline = TRUE)
            )
          )
        ),

        # ---- Alpha overall stats ----
        fluidRow(
          column(
            width = 8,
            offset = 2,
            tags$div(
              style = "margin-top: 10px; font-size: 16px;",
              strong(textOutput("alpha_stats_text", inline = TRUE))
            )
          )
        ),

        # ---- Alpha tabs below plot ----

        fluidRow(
          column(
            width = 8,
            offset = 2,
            div(
              style = "margin-top: 16px;",
              tabsetPanel(
                id = "alpha_results_tabs",
                type = "tabs",

                tabPanel(
                  "Summary statistics",
                  br(),
                  DT::DTOutput("alpha_summary_tbl")
                ),

                tabPanel(
                  "Pairwise p-values",
                  br(),
                  DT::DTOutput("alpha_pairwise_tbl")
                )
              )
            ),
            tags$div(
              style = "margin-top: 8px; font-size: 13px; color: #555;",
              HTML("* p ≤ 0.05&nbsp;&nbsp;&nbsp;** p ≤ 0.01&nbsp;&nbsp;&nbsp;*** p ≤ 0.001&nbsp;&nbsp;&nbsp;. p ≤ 0.1")
            )
          )
        ),

        # ---- Beta controls ----
        div(
          class = "data-select-grid",
          style = "margin-top: 28px;",

          div(class = "data-select-item"),

          div(
            class = "data-select-item",
            selectInput(
              "beta_metric",
              "Beta Diversity",
              choices = c(
                "Bray-Curtis"      = "bray",
                "Jaccard"          = "jaccard",
                "Euclidean"        = "euclidean",
                "Aitchison"        = "aitchison",
                "Robust Aitchison" = "robust.aitchison"
              ),
              selected = "bray"
            )
          ),

          div(class = "data-select-item"),
          div(class = "data-select-item")
        ),

        # ---- Beta plot ----
        fluidRow(
          column(
            width = 8,
            offset = 2,
            div(
              class = "custom-loader-wrap",

              plotly::plotlyOutput("beta_pcoa", height = "700px"),

              div(
                id = "beta_loading_overlay",
                class = "custom-loading-overlay custom-placeholder-overlay",
                div(
                  class = "custom-placeholder-text",
                  "Select the Confirm button to view beta diversity"
                )
              )
            )
          )
        ),

        fluidRow(
          column(
            width = 8,
            offset = 2,
            tags$div(
              style = "margin-top: 10px; font-size: 15px; color: #8a6d6d;",
              textOutput("beta_warning_text", inline = TRUE)
            )
          )
        ),

        fluidRow(
          column(
            width = 8,
            offset = 2,
            tags$div(
              style = "margin-top: 10px; font-size: 16px;",
              strong(textOutput("beta_stats_text", inline = TRUE))
            ),
            tags$div(
              style = "margin-top: 6px; font-size: 15px;",
              textOutput("beta_dispersion_text", inline = TRUE)
            )
          )
        )
      ),


      # ---- Taxonomic Pie Chart ----
      div(
        id = "sec_pie", class = "scroll-section",
        h3("Taxonomic Pie Chart"),

        div(
          class = "custom-loader-wrap",

          taxplore::KronaChartOutput("tax_krona", height = "750px"),

          div(
            id = "tax_loading_overlay",
            class = "custom-loading-overlay custom-placeholder-overlay",
            div(
              class = "custom-placeholder-text",
              "Select the Confirm button to view taxonomy"
            )
          )
        )
      )
    )
  )
}

assign_protocol_ID <- function(df,
                               protocol_columns,
                               protocol_sheet = NULL) {

  # Remove existing protocol_ID if present
  df <- df %>%
    select(-any_of("protocol_ID"))

  # Get distinct protocol definitions from incoming data
  new_protocol_combos <- df %>%
    select(all_of(protocol_columns)) %>%
    distinct()

  # --------------------------------------------------------------
  # CASE 1: No existing protocol_sheet → build from scratch
  # --------------------------------------------------------------
  if (is.null(protocol_sheet) || nrow(protocol_sheet) == 0) {

    protocol_sheet <- new_protocol_combos %>%
      mutate(protocol_ID = row_number())

  } else {

    # Ensure protocol_sheet has required structure
    required_cols <- c(protocol_columns, "protocol_ID")
    missing_cols <- setdiff(required_cols, names(protocol_sheet))

    if (length(missing_cols) > 0) {
      stop("protocol_sheet is missing required columns: ",
           paste(missing_cols, collapse = ", "))
    }

    # --------------------------------------------------------------
    # Add new protocol_IDs for unseen combinations
    # --------------------------------------------------------------

    unseen_protocols <- anti_join(
      new_protocol_combos,
      protocol_sheet %>% select(all_of(protocol_columns)),
      by = protocol_columns
    )

    if (nrow(unseen_protocols) > 0) {

      max_id <- max(protocol_sheet$protocol_ID, na.rm = TRUE)

      unseen_protocols <- unseen_protocols %>%
        mutate(protocol_ID = row_number() + max_id)

      protocol_sheet <- bind_rows(protocol_sheet, unseen_protocols)
    }
  }

  # --------------------------------------------------------------
  # Assign protocol_ID back to df
  # --------------------------------------------------------------

  df_with_ids <- df %>%
    left_join(protocol_sheet, by = protocol_columns)

  return(list(
    data = df_with_ids,
    protocol_sheet = protocol_sheet
  ))
}

app_b_server <- function(input, output, session){

  protocol_source <- reactive({
    df <- selection_panel_df()

    shiny::validate(
      shiny::need(!is.null(df), "Select a site/cell/polygon to view protocols."),
      shiny::need(nrow(df) > 0, "No detections available for this filtered selection.")
    )

    live_filters <- list(
      target_gene = input$div_target_gene %||% character(0),
      primers     = input$div_primer %||% character(0)
    )

    df <- apply_diversity_dropdown_filters(df, live_filters)

    shiny::validate(
      shiny::need(nrow(df) > 0, "No detections available for the selected target gene/primer filters.")
    )

    protocol_columns <- c(
      "nucl_acid_ext_kit",
      "platform",
      "instrument",
      "seq_kit",
      "otu_db",
      "tax_assign_cat",
      "otu_seq_comp_appr",
      "min_depth_floor",
      "max_depth_floor",
      "samp_size_mid",
      "size_frac",
      "filter_material",
      "samp_mat_process",
      "samp_store_temp",
      "samp_store_sol"
    )

    protocol_columns <- protocol_columns[protocol_columns %in% names(df)]

    shiny::validate(
      shiny::need(length(protocol_columns) > 0, "No protocol columns were found.")
    )

    assign_protocol_ID(
      df = df,
      protocol_columns = protocol_columns,
      protocol_sheet = NULL
    )$data
  })

  meta_all <- reactive({
    protocol_source()
  })

  protocol_data <- reactive({
    df <- meta_all()

    if (!"detected" %in% names(df)) {
      df <- df %>%
        dplyr::mutate(
          detected = dplyr::if_else(
            !is.na(organismQuantity) & organismQuantity > 0,
            1L, 0L
          )
        )
    }

    df
  })

  protocol_info <- reactive({
    protocol_data() %>%
      dplyr::group_by(protocol_ID) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()
  })

  protocol_summary <- reactive({
    protocol_data() %>%
      dplyr::group_by(protocol_ID, samp_name) %>%
      dplyr::summarise(
        detected_sample = any(detected == 1, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::group_by(protocol_ID) %>%
      dplyr::summarise(
        total_samples = dplyr::n(),
        total_detections = sum(detected_sample),
        detection_rate = 100 * total_detections / total_samples,
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(total_detections))
  })

  observeEvent(protocol_summary(), {
    ps <- protocol_summary()

    shiny::validate(
      shiny::need(nrow(ps) > 0, "No protocols found.")
    )

    choices <- as.list(ps$protocol_ID)

    names(choices) <- paste0(
      "Protocol ", ps$protocol_ID,
      " (", ps$total_detections, " detections | ",
      sprintf("%.1f", ps$detection_rate), "%)"
    )

    # ---- Single protocol dropdown ----
    session$onFlushed(function() {

      current <- suppressWarnings(as.numeric(isolate(input$req_protocol)))

      selected_id <- if (
        length(current) == 1 &&
        !is.na(current) &&
        current %in% ps$protocol_ID
      ) {
        current
      } else {
        ps$protocol_ID[1]
      }

      updateSelectInput(
        session,
        "req_protocol",
        choices = choices,
        selected = selected_id
      )

    }, once = TRUE)

    # ---- Multi-protocol picker ----
    shinyWidgets::updatePickerInput(
      session,
      inputId = "prot_id",
      choices = choices,
      selected = NULL
    )

  }, ignoreInit = FALSE)

  # helper: pick a single display value (unique or collapse)
  pick_display <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(trimws(x))]
    if (length(x) == 0) return(NA_character_)
    ux <- unique(x)
    if (length(ux) == 1) ux else paste(ux, collapse = " | ")
  }

  # rows matching current ProtocolID
  selected_protocol_rows <- reactive({
    ps <- protocol_summary()

    selected_id <- suppressWarnings(as.numeric(input$req_protocol))

    if (is.na(selected_id) && nrow(ps) > 0) {
      selected_id <- ps$protocol_ID[1]
    }

    protocol_info() %>%
      dplyr::filter(protocol_ID == selected_id)
  })

  div_unlocked <- reactiveVal(FALSE)

  # --- base protocol per Location (the "starting choice") ---
  base_protocol <- reactiveVal(NULL)

  observeEvent(protocol_summary(), {
    ps <- protocol_summary()
    base_protocol(if (nrow(ps) > 0) ps$protocol_ID[1] else NULL)
  }, ignoreInit = FALSE)

  # live details card
  output$protocol_details <- renderUI({
    df <- selected_protocol_rows()

    if (nrow(df) == 0) {
      return(tags$div(style="margin-top:10px;", em("No rows found for this Location + ProtocolID.")))
    }

    method_groups <- list(
      "Field Methods" = c("samp_size","size_frac","filter_material","samp_mat_process",
                          "minimumDepthInMeters","maximumDepthInMeters"),
      "Storage and Extraction Methods" = c("samp_store_temp","samp_store_sol", "nucl_acid_ext_kit"),
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
      df_base <- protocol_info() %>%
        dplyr::filter(protocol_ID == as.numeric(bp))
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

                tags$div(
                  tags$div(class="protocol-field-title", f),
                  tags$div(
                    class = paste("protocol-card", if (changed) "changed", if (is_na) "na"),
                    if (is_na) "—" else cur2
                  )
                )
              })
            )
          )
        })
      )
    )
  })

  output$protocol_nmds_plot <- plotly::renderPlotly({
    ps <- protocol_summary()

    shiny::validate(
      shiny::need(nrow(ps) > 2, "NMDS requires more than two protocols.")
    )

    top_ids <- ps %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::pull(protocol_ID)

    filtered_protocol_sheet <- protocol_info() %>%
      dplyr::filter(protocol_ID %in% top_ids) %>%
      dplyr::mutate(
        dplyr::across(
          where(~ inherits(.x, c("Date", "POSIXct", "POSIXlt"))),
          as.character
        )
      )

    protocol_nmds(filtered_protocol_sheet)
  })

  output$protocol_barplot <- plotly::renderPlotly({
    ps <- protocol_summary()

    shiny::validate(
      shiny::need(nrow(ps) > 0, "No protocol summary available.")
    )

    top_protocols <- ps %>%
      dplyr::slice_max(order_by = total_detections, n = 10, with_ties = FALSE) %>%
      dplyr::arrange(protocol_ID)

    protocol_bargraph(top_protocols)
  })

  add_primer_combo_vec <- function(fwd, rev) {
    fwd <- trimws(fwd)
    rev <- trimws(rev)

    fwd[fwd == ""] <- NA_character_
    rev[rev == ""] <- NA_character_

    dplyr::case_when(
      !is.na(fwd) & !is.na(rev) ~ paste(fwd, rev, sep = " | "),
      !is.na(fwd) &  is.na(rev) ~ fwd,
      is.na(fwd) & !is.na(rev) ~ rev,
      TRUE ~ NA_character_
    )
  }

  safe_chr <- function(data, col) {
    if (col %in% names(data)) as.character(data[[col]]) else NA_character_
  }

  safe_num <- function(data, col) {
    if (col %in% names(data)) suppressWarnings(as.numeric(data[[col]])) else NA_real_
  }

  selection_map_df <- reactive({
    det <- selected_detections()

    if (is.null(det) || nrow(det) == 0) return(det)

    det %>%
      sf::st_drop_geometry() %>%
      dplyr::transmute(
        id = safe_chr(., "id"),
        occurrenceID = safe_chr(., "occurrenceID"),
        samp_name = safe_chr(., "samp_name"),
        scientificName = safe_chr(., "scientificName"),
        year = safe_chr(., "year"),
        month = safe_chr(., "month"),
        target_gene = safe_chr(., "target_gene"),
        organismQuantity = safe_num(., "organismQuantity"),
        kingdom = safe_chr(., "kingdom"),
        phylum = safe_chr(., "phylum"),
        class = safe_chr(., "class"),
        order = safe_chr(., "order"),
        family = safe_chr(., "family"),
        genus = safe_chr(., "genus"),
        pcr_primer_name_forward = safe_chr(., "pcr_primer_name_forward"),
        pcr_primer_name_reverse = safe_chr(., "pcr_primer_name_reverse"),
        eventDate_clean = safe_chr(., "eventDate_clean"),
        samp_size = safe_num(., "samp_size"),
        size_frac = safe_chr(., "size_frac"),
        filter_material = safe_chr(., "filter_material"),
        samp_mat_process = safe_chr(., "samp_mat_process"),
        samp_store_temp = safe_chr(., "samp_store_temp"),
        samp_store_sol = safe_chr(., "samp_store_sol"),
        nucl_acid_ext_kit = safe_chr(., "nucl_acid_ext_kit"),
        platform = safe_chr(., "platform"),
        instrument = safe_chr(., "instrument"),
        seq_kit = safe_chr(., "seq_kit"),
        otu_db = safe_chr(., "otu_db"),
        tax_assign_cat = safe_chr(., "tax_assign_cat"),
        otu_seq_comp_appr = safe_chr(., "otu_seq_comp_appr"),
        category = safe_chr(., "category"),
        minimumDepthInMeters = safe_num(., "minimumDepthInMeters"),
        maximumDepthInMeters = safe_num(., "maximumDepthInMeters"),
        min_depth_floor = safe_num(., "min_depth_floor"),
        max_depth_floor = safe_num(., "max_depth_floor"),
        samp_size_mid = safe_num(., "samp_size_mid"),
        datasetID_obis = safe_chr(., "datasetID_obis"),
        ownerContact = safe_chr(., "ownerContact"),
        eventDate = safe_chr(., "eventDate"),
        bibliographicCitation = safe_chr(., "bibliographicCitation"),
        max_depth_bin = safe_chr(., "max_depth_bin"),
        min_depth_bin = safe_chr(., "min_depth_bin"),
        samp_size_bin = safe_chr(., "samp_size_bin")
      )
  })

  selection_panel_df <- reactive({
    df <- selection_map_df()
    if (is.null(df) || nrow(df) == 0) return(df)

    live_filters <- list(
      year = sel_year_chr(),
      groups = active_groups(),
      iucn_on = filter_iucn_on(),
      sara_on = filter_sara_on(),
      ais_on = filter_ais_on()
    )

    apply_species_filters(df, live_filters)
  })

  # ---- confirmed diversity controls ----
  div_filters <- reactiveVal(
    list(
      target_gene = character(0),
      primers     = character(0),
      polygons      = character(0)
    )
  )

  selection_selection_df <- reactive({
    req(confirmed_map_filters())

    df <- selection_map_df()

    df <- apply_species_filters(
      df,
      confirmed_map_filters()
    )

    if (is.null(df) || nrow(df) == 0) return(df)

    live_filters <- list(
      target_gene = input$div_target_gene %||% character(0),
      primers     = input$div_primer %||% character(0),
      polygons    = input$div_compare_polygons %||% character(0)
    )

    apply_diversity_dropdown_filters(df, live_filters)
  })

  panel_ids <- reactive({
    df <- selection_panel_df()
    if (is.null(df) || nrow(df) == 0) return(character(0))
    unique(df$id)
  })

  selection_ids <- reactive({
    df <- selection_selection_df()
    if (is.null(df) || nrow(df) == 0) return(character(0))
    unique(df$id)
  })


  mpa_membership_base_df <- reactive({

    if (is.null(diversity_mpa_df)) {
      stop("mpa_membership_base_df: diversity_mpa_df is NULL")
    }
    diversity_mpa_df %>%
      dplyr::filter(
        !is.na(site_name), trimws(site_name) != "",
        !is.na(site_type), trimws(site_type) != ""
      )
  })

  mpa_membership_panel_df <- reactive({
    df <- mpa_membership_base_df()

    if (is.null(df) || nrow(df) == 0) {
      return(df)
    }

    yr <- sel_year_chr()
    if (yr != "All") {
      df <- df %>%
        dplyr::filter(as.character(year) == yr)
    }

    df <- apply_species_filters(df)

    df
  })

  mpa_membership_selection_df <- reactive({
    df <- mpa_membership_panel_df()

    if (is.null(df) || nrow(df) == 0) {
      return(df)
    }

    live_filters <- list(
      target_gene = input$div_target_gene %||% character(0),
      primers     = input$div_primer %||% character(0),
      polygons    = input$div_compare_polygons %||% character(0)
    )

    df <- apply_diversity_dropdown_filters(df, live_filters)
    apply_compare_polygon_filter(df, live_filters$polygons)
  })

  polygon_membership_base_df <- reactive({
    pts <- diversity_beta_df
    polys <- drawn_polys()

    if (is.null(pts) || nrow(pts) == 0) {
      return(pts)
    }

    if (is.null(polys) || nrow(polys) == 0) {
      return(sf::st_drop_geometry(pts[0, , drop = FALSE]))
    }

    drawn_list <- lapply(seq_len(nrow(polys)), function(i) {
      g_i <- sf::st_geometry(polys[i, , drop = FALSE])
      lab <- polys$draw_label[i]

      inside_i <- pts[within_any(pts, g_i), , drop = FALSE]
      if (nrow(inside_i) == 0) return(NULL)

      inside_i %>%
        dplyr::mutate(
          site_name = lab,
          site_type = "User"
        )
    })

    out <- dplyr::bind_rows(drawn_list)

    if (is.null(out) || nrow(out) == 0) {
      return(pts[0, , drop = FALSE] %>% sf::st_drop_geometry())
    }

    out <- out %>%
      sf::st_drop_geometry() %>%
      dplyr::arrange(occurrenceID) %>%
      dplyr::group_by(
        occurrenceID, point_key, samp_name, scientificName, year, target_gene,
        site_name, site_type
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    out
  })

  polygon_membership_panel_df <- reactive({
    df <- polygon_membership_base_df()

    if (is.null(df) || nrow(df) == 0) {
      return(df)
    }

    yr <- sel_year_chr()
    if (yr != "All") {
      df <- df %>% dplyr::filter(as.character(year) == yr)
    }

    df <- apply_species_filters(df)

    df
  })

  polygon_membership_selection_df <- reactive({
    df <- polygon_membership_panel_df()

    if (is.null(df) || nrow(df) == 0) {
      return(df)
    }

    #df <- apply_diversity_dropdown_filters(df, div_filters())
    live_filters <- list(
      target_gene = input$div_target_gene %||% character(0),
      primers     = input$div_primer %||% character(0),
      polygons    = input$div_compare_polygons %||% character(0)
    )

    df <- apply_diversity_dropdown_filters(df, live_filters)
    df <- apply_compare_polygon_filter(df, live_filters$polygons)

    if (nrow(df) == 0) {
      return(df)
    }

    df %>%
      dplyr::arrange(occurrenceID) %>%
      dplyr::group_by(
        occurrenceID, point_key, samp_name, scientificName, year, target_gene,
        site_name, site_type
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
   })

  observeEvent(input$div_apply, {

    div_unlocked(TRUE)

    shinyjs::html(
      "alpha_loading_overlay",
      '<div class="loader"></div>'
    )
    shinyjs::removeClass("alpha_loading_overlay", "hidden")

    shinyjs::html(
      "beta_loading_overlay",
      '<div class="loader"></div>'
    )
    shinyjs::removeClass("beta_loading_overlay", "hidden")

    shinyjs::html(
      "tax_loading_overlay",
      '<div class="loader"></div>'
    )
    shinyjs::removeClass("tax_loading_overlay", "hidden")

  }, ignoreInit = TRUE)

  ############################################################


  #NEW CODE ABOVE HERE



  ############################################################







  ############################################################



  #OLD CODE STARTS HERE


  ############################################################
























































  prune_cache <- function(cache_list, max_n = 100) {
    nms <- names(cache_list)
    if (length(nms) <= max_n) return(cache_list)

    keep <- tail(nms, max_n)
    cache_list[keep]
  }

  `%||%` <- function(x, y) if (is.null(x)) y else x

  make_sample_id <- function(df) {
    df %>%
      dplyr::mutate(
        samp_name  = as.character(samp_name),
        year       = if ("year" %in% names(.)) as.character(year) else NA_character_,
        eventDate  = if ("eventDate" %in% names(.)) as.character(eventDate) else NA_character_,
        point_key  = if ("point_key" %in% names(.)) as.character(point_key) else NA_character_,
        sample_id  = dplyr::case_when(
          !is.na(point_key) & point_key != "" ~ point_key,
          !is.na(samp_name) & samp_name != "" & !is.na(eventDate) & eventDate != "" ~
            paste(samp_name, eventDate, sep = " || "),
          !is.na(samp_name) & samp_name != "" & !is.na(year) & year != "" ~
            paste(samp_name, year, sep = " || "),
          !is.na(samp_name) & samp_name != "" ~
            samp_name,
          TRUE ~ NA_character_
        )
      )
  }

  species_cache <- reactiveValues(data = list())

  make_species_cache_key <- function(selected_key, year, groups, iucn_on, sara_on, ais_on) {
    paste(
      selected_key,
      year,
      paste(sort(groups), collapse = "|"),
      paste0("iucn=", iucn_on),
      paste0("sara=", sara_on),
      paste0("ais=", ais_on),
      sep = "~~"
    )
  }

  get_species_by_layer_cached <- function(selected_key, det_sf, year, groups, iucn_on, sara_on, ais_on) {

    cache_key <- make_species_cache_key(
      selected_key = selected_key,
      year         = year,
      groups       = groups,
      iucn_on      = iucn_on,
      sara_on      = sara_on,
      ais_on       = ais_on
    )

    if (!is.null(species_cache$data[[cache_key]])) {
      return(species_cache$data[[cache_key]])
    }


    det_f <- apply_species_filters(det_sf) %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        scientificName = as.character(scientificName),
        target_gene    = as.character(target_gene)
      )

    out <- list(
      All = det_f %>%
        dplyr::pull(scientificName) %>%
        unique() %>%
        stats::na.omit() %>%
        sort()
    )

    for (g in c("12S", "COI", "16S", "18S")) {
      out[[g]] <- det_f %>%
        dplyr::filter(target_gene == g) %>%
        dplyr::pull(scientificName) %>%
        unique() %>%
        stats::na.omit() %>%
        sort()
    }

    species_cache$data[[cache_key]] <- out
    species_cache$data <- prune_cache(species_cache$data, max_n = 100)
    out
  }

  inside_cache <- reactiveValues(data = list())

  make_inside_cache_key <- function(selected_key, year) {
    paste(selected_key, year, sep = "~~")
  }

  get_inside_cached <- function(selected_key, pts, geom, year) {
    cache_key <- make_inside_cache_key(selected_key, year)

    if (!is.null(inside_cache$data[[cache_key]])) {
      return(inside_cache$data[[cache_key]])
    }

    inside <- pts[within_any(pts, geom), , drop = FALSE]
    inside_cache$data[[cache_key]] <- inside
    inside_cache$data <- prune_cache(inside_cache$data, max_n = 100)
    inside
  }

  protocol_bargraph <- function(df) {
    color_vec <- c("#00A08A", "#446455", "#Fdd262", "#5BBCD6", "#046c9a", "#ABDDDE", "#d3dddc")

    df <- df %>%
      dplyr::mutate(
        color = color_vec[(seq_len(dplyr::n()) - 1) %% length(color_vec) + 1],
        protocol_ID = factor(protocol_ID, levels = sort(unique(protocol_ID)))
      )

    plotly::plot_ly(
      data = df,
      x = ~protocol_ID,
      y = ~detection_rate,
      type = "bar",
      marker = list(color = df$color),
      hoverinfo = "text",
      name = "Detection Rate"
    ) %>%
      plotly::layout(
        xaxis = list(title = "Protocol ID"),
        yaxis = list(title = "Detection Rate (%)", showgrid = FALSE),
        margin = list(l = 60, r = 20, t = 50, b = 60),
        legend = list(title = list(text = ""))
      ) %>%
      plotly::config(
        displayModeBar = TRUE,
        modeBarButtonsToAdd = c("resetScale2d")
      )
  }

  get_forward_primer_col <- function(df) {
    cand <- c("pcr_primer_name_forward", "pcr_primer_forward")
    hit <- intersect(cand, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }

  get_reverse_primer_col <- function(df) {
    cand <- c("pcr_primer_name_reverse", "pcr_primer_reverse")
    hit <- intersect(cand, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }

  add_primer_combo <- function(df) {
    fwd_col <- get_forward_primer_col(df)
    rev_col <- get_reverse_primer_col(df)

    if (is.null(fwd_col) && is.null(rev_col)) {
      df$primer_combo <- NA_character_
      return(df)
    }

    fwd <- if (!is.null(fwd_col)) as.character(df[[fwd_col]]) else rep(NA_character_, nrow(df))
    rev <- if (!is.null(rev_col)) as.character(df[[rev_col]]) else rep(NA_character_, nrow(df))

    fwd <- trimws(fwd)
    rev <- trimws(rev)

    fwd[fwd == ""] <- NA_character_
    rev[rev == ""] <- NA_character_

    df$primer_combo <- dplyr::case_when(
      !is.na(fwd) & !is.na(rev) ~ paste(fwd, rev, sep = " | "),
      !is.na(fwd) &  is.na(rev) ~ fwd,
      is.na(fwd) & !is.na(rev) ~ rev,
      TRUE ~ NA_character_
    )

    df
  }

  apply_diversity_dropdown_filters <- function(df, filters) {
    df <- add_primer_combo(df)

    if (length(filters$target_gene) > 0) {
      df <- df %>%
        dplyr::filter(as.character(target_gene) %in% filters$target_gene)
    }

    if (length(filters$primers) > 0) {
      df <- df %>%
        dplyr::filter(primer_combo %in% filters$primers)
    }

    df
  }

  apply_compare_polygon_filter <- function(df, polygons) {
    polygons <- as.character(polygons %||% character(0))

    if (length(polygons) == 0 || is.null(df) || nrow(df) == 0) {
      return(df)
    }

    if (!"site_name" %in% names(df)) {
      return(df[0, , drop = FALSE])
    }

    df %>%
      dplyr::filter(as.character(site_name) %in% polygons)
  }

  within_any <- function(x_sf, geom) {
    if (is.null(geom) || is.null(x_sf) || nrow(x_sf) == 0) {
      return(rep(FALSE, if (is.null(x_sf)) 0 else nrow(x_sf)))
    }

    if (inherits(geom, "sfc")) {
      geom <- sf::st_sf(geometry = geom)
    }

    if (sf::st_crs(x_sf) != sf::st_crs(geom)) {
      geom <- sf::st_transform(geom, sf::st_crs(x_sf))
    }

    lengths(sf::st_intersects(x_sf, geom)) > 0
  }

  # store ALL drawn polygons (one row per polygon), with a stable id
  empty_drawn_sf <- sf::st_sf(
    draw_id    = character(0),
    draw_label = character(0),
    geometry   = sf::st_sfc(crs = 4326)
  )

  drawn_polys <- reactiveVal(empty_drawn_sf)
  selected_draw_id <- reactiveVal(NULL)


  # ---- helper: convert leaflet.draw feature -> sf polygon (EPSG:4326) ----
  feature_to_sf <- function(feature) {
    req(feature$geometry$type)
    type <- feature$geometry$type
    coords <- feature$geometry$coordinates

    close_ring <- function(mat) {
      if (!all(mat[1, ] == mat[nrow(mat), ])) {
        mat <- rbind(mat, mat[1, , drop = FALSE])
      }
      mat
    }

    if (type == "Polygon") {
      rings <- lapply(coords, function(ring) {
        mat <- do.call(rbind, lapply(ring, function(x) c(x[[1]], x[[2]])))
        close_ring(mat)
      })

      poly <- sf::st_polygon(rings)
      out  <- sf::st_sf(geometry = sf::st_sfc(poly, crs = 4326))
      return(sf::st_make_valid(out))
    }

    if (type == "MultiPolygon") {
      mp <- lapply(coords, function(poly_i) {
        lapply(poly_i, function(ring) {
          mat <- do.call(rbind, lapply(ring, function(x) c(x[[1]], x[[2]])))
          close_ring(mat)
        })
      })

      geom <- sf::st_multipolygon(mp)
      out  <- sf::st_sf(geometry = sf::st_sfc(geom, crs = 4326))
      return(sf::st_make_valid(out))
    }

    stop("Drawn feature type not supported: ", type)
  }

  # ---- selection geometry (drawn polygon OR clicked polygon OR clicked grid cell) ----
  selection_geom <- reactive({
    sel_id <- selected_draw_id()
    if (!is.null(sel_id)) {
      polys <- drawn_polys()
      hit <- polys %>% dplyr::filter(draw_id == sel_id)
      if (nrow(hit) > 0) return(sf::st_geometry(hit))
    }

    click <- input$map_shape_click
    if (is.null(click) || is.null(click$id)) return(NULL)

    if (grepl("\\|\\|", click$id)) {
      parts <- strsplit(click$id, "\\|\\|")[[1]]
      p_type <- parts[1]
      p_name <- parts[2]
      poly_sel <- all_polys_click %>% dplyr::filter(site_type == p_type, site_name == p_name)
      if (nrow(poly_sel) == 0) return(NULL)
      return(sf::st_geometry(poly_sel))
    }

    cid <- suppressWarnings(as.integer(click$id))
    if (is.na(cid)) return(NULL)

    cell_poly <- grid_clip %>% dplyr::filter(cell_id == cid)
    if (nrow(cell_poly) == 0) return(NULL)

    sf::st_geometry(cell_poly)
  })

  observe({
    click <- input$map_shape_click
    proxy <- leafletProxy("map")

    proxy %>% clearGroup("Selection outline")

    if (is.null(click) || is.null(click$id)) return()

    id <- as.character(click$id)

    # only draw white outline for clicked grid cells
    cid <- suppressWarnings(as.integer(id))
    if (is.na(cid)) return()

    sel_sf <- grid_clip %>%
      dplyr::filter(cell_id == cid)

    if (nrow(sel_sf) == 0) return()

    proxy %>%
      addPolygons(
        data        = sel_sf,
        group       = "Selection outline",
        layerId     = ~paste0("selected_cell_", cell_id),
        fill        = FALSE,
        color       = "white",
        weight      = 3,
        opacity     = 1,
        options     = pathOptions(
          pane = "pane_selected_top",
          interactive = FALSE
        )
      )
  })

  # ---- Year selection (as character or "All") ----
  sel_year_chr <- reactive({
    yr <- input$sel_year %||% "All"
    as.character(yr)
  })


  # --- Is the Sampling points layer currently visible? ---
  sampling_points_layer_on <- reactive({
    groups_on <- input$map_groups %||% character(0)
    "Sampling Points" %in% groups_on
  })

  # Which richness layers are currently ON (including "All")
  active_richness_layers <- reactive({
    groups_on <- input$map_groups %||% character(0)
    intersect(groups_on, c("All","12S","COI","16S","18S"))
  })

  make_list <- function(vec, max_h = 220) {
    if (length(vec) == 0) return(em("No species detected for this layer."))
    tags$div(
      style = paste0("max-height:", max_h, "px; overflow-y:auto; padding-left: 10px;"),
      tags$ul(lapply(vec, tags$li))
    )
  }

  #remove cache of deleted polygons
  observeEvent(input$map_draw_deleted_features, {
    del <- input$map_draw_deleted_features
    if (is.null(del) || is.null(del$features) || length(del$features) == 0) return()

    del_ids <- vapply(del$features, function(f) {
      if (!is.null(f$properties) && !is.null(f$properties$`_leaflet_id`)) {
        as.character(f$properties$`_leaflet_id`)
      } else {
        NA_character_
      }
    }, character(1))

    del_ids <- del_ids[!is.na(del_ids) & nzchar(del_ids)]
    if (length(del_ids) == 0) return()

    old_keys <- names(species_cache$data)
    if (length(old_keys) > 0) {
      keep_keys <- old_keys[!vapply(old_keys, function(k) {
        any(vapply(del_ids, function(id) grepl(paste0("draw:", id), k, fixed = TRUE), logical(1)))
      }, logical(1))]

      species_cache$data <- species_cache$data[keep_keys]
    }

    old_inside_keys <- names(inside_cache$data)
    if (length(old_inside_keys) > 0) {
      keep_inside_keys <- old_inside_keys[!vapply(old_inside_keys, function(k) {
        any(vapply(del_ids, function(id) grepl(paste0("draw:", id), k, fixed = TRUE), logical(1)))
      }, logical(1))]

      inside_cache$data <- inside_cache$data[keep_inside_keys]
    }

  })

  # ---- year dropdown: only show years present in current selection geom ----
  observeEvent(selection_geom(), {
    pts <- species_sf_min
    shiny::req(pts)
    shiny::req(nrow(pts) > 0)

    g <- selection_geom()

    sel_id <- selected_draw_id()
    click  <- input$map_shape_click

    selected_key <- NULL
    if (!is.null(sel_id)) {
      selected_key <- paste0("draw:", sel_id)
    } else if (!is.null(click$id)) {
      selected_key <- paste0("click:", click$id)
    }

    inside <- if (!is.null(g) && !is.null(selected_key)) {
      get_inside_cached(
        selected_key = selected_key,
        pts          = pts,
        geom         = g,
        year         = "AllYearsForDropdown"
      )
    } else {
      pts
    }

    yrs <- inside %>%
      sf::st_drop_geometry() %>%
      dplyr::pull(year) %>%
      as.character() %>%
      unique() %>%
      stats::na.omit() %>%
      sort()

    cur <- as.character(input$sel_year %||% "All")
    new_choices  <- c("All", yrs)
    new_selected <- if (cur %in% new_choices) cur else "All"

    updateSelectInput(session, "sel_year", choices = new_choices, selected = new_selected)
  }, ignoreInit = FALSE)

  observeEvent(input$map_draw_new_feature, {
    feat <- input$map_draw_new_feature
    req(feat)

    poly_sf <- tryCatch(
      feature_to_sf(feat),
      error = function(e) {
        showNotification(
          paste("Could not read drawn polygon:", e$message),
          type = "error"
        )
        return(NULL)
      }
    )

    req(!is.null(poly_sf), nrow(poly_sf) > 0)

    fid <- NULL
    if (!is.null(feat$properties) && !is.null(feat$properties$`_leaflet_id`)) {
      fid <- as.character(feat$properties$`_leaflet_id`)
    }
    if (is.null(fid) || !nzchar(fid)) {
      fid <- paste0("draw_", as.integer(Sys.time()), "_", sample.int(1e6, 1))
    }

    cur <- drawn_polys()

    # keep label numbering stable by draw order
    next_num <- if (nrow(cur) == 0) 1L else {
      old_nums <- suppressWarnings(as.integer(sub("^Polygon\\s+", "", cur$draw_label)))
      max(old_nums, na.rm = TRUE) + 1L
    }
    if (!is.finite(next_num)) next_num <- 1L

    poly_sf$draw_id    <- fid
    poly_sf$draw_label <- paste0("Polygon ", next_num)

    if (nrow(cur) > 0) {
      cur <- cur[!(cur$draw_id %in% fid), , drop = FALSE]
    }

    drawn_polys(dplyr::bind_rows(cur, poly_sf))
    selected_draw_id(fid)

    session$sendCustomMessage("openFloating", list(id = fid))
  }, ignoreInit = TRUE)

  observeEvent(input$map_draw_deleted_features, {
    del <- input$map_draw_deleted_features
    if (is.null(del) || is.null(del$features) || length(del$features) == 0) return()

    del_ids <- vapply(del$features, function(f){
      if (!is.null(f$properties) && !is.null(f$properties$`_leaflet_id`)) {
        as.character(f$properties$`_leaflet_id`)
      } else {
        NA_character_
      }
    }, character(1))

    del_ids <- del_ids[!is.na(del_ids) & nzchar(del_ids)]
    if (length(del_ids) == 0) return()

    cur <- drawn_polys()
    cur <- cur[!(cur$draw_id %in% del_ids), , drop = FALSE]
    drawn_polys(cur)

    sid <- selected_draw_id()
    if (!is.null(sid) && sid %in% del_ids) {
      if (nrow(cur) > 0) {
        selected_draw_id(cur$draw_id[nrow(cur)])
      } else {
        selected_draw_id(NULL)
      }
    }
  })

  compare_polygon_choices <- reactive({

    live_filters <- list(
      target_gene = input$div_target_gene %||% character(0),
      primers     = input$div_primer %||% character(0)
    )

    base_df <- mpa_membership_panel_df()

    if (!is.null(base_df) && nrow(base_df) > 0) {
      base_df <- apply_diversity_dropdown_filters(base_df, live_filters)
    }

    draw_df <- polygon_membership_panel_df()

    if (!is.null(draw_df) && nrow(draw_df) > 0) {
      draw_df <- apply_diversity_dropdown_filters(draw_df, live_filters)
    }

    get_site_choices <- function(df, user_only = FALSE) {
      if (is.null(df) || nrow(df) == 0) return(character(0))

      if (isTRUE(user_only) && "site_type" %in% names(df)) {
        df <- df %>% dplyr::filter(site_type == "User")
      }

      df %>%
        dplyr::filter(!is.na(site_name), site_name != "") %>%
        dplyr::pull(site_name) %>%
        as.character() %>%
        unique()
    }

    base_groups  <- get_site_choices(base_df)
    drawn_groups <- get_site_choices(draw_df, user_only = TRUE)

    sort(unique(c(base_groups, drawn_groups)))
  })


  #Monthly sampling
  monthly_sample_counts <- reactive({
    det <- selection_map_df()
    req(det)

    det0 <- det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        samp_name = as.character(samp_name),
        month_chr = as.character(month)
      )

    shiny::validate(
      shiny::need(nrow(det0) > 0, "No detections in current selection.")
    )

    month_levels <- month.abb

    counts <- det0 %>%
      dplyr::filter(
        !is.na(samp_name), samp_name != "",
        !is.na(month_chr), month_chr %in% month_levels
      ) %>%
      dplyr::distinct(samp_name, month_chr) %>%
      dplyr::count(month_chr, name = "n_samples") %>%
      dplyr::mutate(month_chr = factor(month_chr, levels = month_levels))

    data.frame(
      month_chr = factor(month_levels, levels = month_levels)
    ) %>%
      dplyr::left_join(counts, by = "month_chr") %>%
      dplyr::mutate(n_samples = dplyr::coalesce(n_samples, 0L))
  })

  output$monthly_plot_subtitle <- renderText({
    click <- input$map_shape_click
    sel_id <- selected_draw_id()
    yr <- sel_year_chr()

    if (!is.null(sel_id)) {
      polys <- drawn_polys()
      hit <- polys %>% dplyr::filter(draw_id == sel_id)
      lab <- if (nrow(hit) > 0) hit$draw_label[1] else "Drawn polygon"
      return(paste0(lab, " | Year: ", yr))
    }

    if (!is.null(click$id)) {
      if (grepl("\\|\\|", click$id)) {
        parts <- strsplit(click$id, "\\|\\|")[[1]]
        return(paste0(parts[2], " | Year: ", yr))
      }

      cid <- suppressWarnings(as.integer(click$id))
      if (!is.na(cid)) {
        return(paste0("Grid cell ", cid, " | Year: ", yr))
      }
    }

    paste0("No selection | Year: ", yr)
  })

  output$monthly_circular_plot <- renderPlot({

    groups_on <- input$map_groups %||% character(0)
    monthly_on <- "Monthly Sampling Plot" %in% groups_on

    # Only render once the monthly plot layer is actually turned on
    req(monthly_on)

    det <- selection_map_df()

    # ---- default placeholder before any cell/polygon is selected ----
    if (is.null(det) || nrow(det) == 0) {
      p_empty <- ggplot2::ggplot(
        data.frame(x = 0.5, y = 0.5, lab = "Select or draw a polygon"),
        ggplot2::aes(x, y)
      ) +
        ggplot2::geom_text(
          ggplot2::aes(label = lab),
          size = 6,
          colour = "grey35"
        ) +
        ggplot2::xlim(0, 1) +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_void() +
        ggplot2::theme(
          plot.background  = ggplot2::element_rect(fill = NA, colour = NA),
          panel.background = ggplot2::element_rect(fill = NA, colour = NA),
          plot.margin = ggplot2::margin(10, 10, 10, 10)
        )

      print(p_empty)
      return(invisible(NULL))
    }

    dat <- monthly_sample_counts()

    shiny::validate(
      shiny::need(nrow(dat) > 0, "No monthly data available."),
      shiny::need(sum(dat$n_samples, na.rm = TRUE) > 0, "No samples available for this selection.")
    )

    dat <- dat %>%
      dplyr::mutate(
        month_chr = factor(month_chr, levels = month.abb),
        fill_group = ifelse(n_samples == 0, "zero", as.character(month_chr))
      )

    ymax <- max(dat$n_samples, na.rm = TRUE)

    outer_max    <- max(1, ceiling(ymax * 1.10))
    inner_offset <- outer_max * 0.42
    label_radius <- outer_max * 1.28
    top_pad      <- outer_max * 0.06

    ring_vals <- c(0, 0.25, 0.50, 0.75, 1.00) * outer_max + inner_offset

    month_cols <- setNames(
      grDevices::colorRampPalette(cool)(12),
      month.abb
    )

    fill_vals <- c(
      zero = "white",
      month_cols
    )

    ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x = month_chr,
        y = n_samples + inner_offset,
        fill = fill_group
      )
    ) +
      ggplot2::geom_hline(
        yintercept = ring_vals,
        colour = "grey80",
        linewidth = 0.5
      ) +
      ggplot2::geom_col(
        width = 0.88,
        colour = NA
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          y = n_samples + inner_offset + outer_max * 0.10,
          label = n_samples
        ),
        size = 3.0,
        color = "black"
      ) +
      ggplot2::geom_text(
        data = dat,
        ggplot2::aes(
          x = month_chr,
          y = inner_offset + label_radius,
          label = month_chr
        ),
        inherit.aes = FALSE,
        size = 4.2,
        color = "black"
      ) +
      ggplot2::coord_polar(start = -pi / 12) +
      ggplot2::scale_y_continuous(
        limits = c(0, inner_offset + label_radius + top_pad),
        expand = c(0, 0)
      ) +
      ggplot2::scale_fill_manual(
        values = fill_vals,
        guide = "none"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_rect(fill = NA, colour = NA),
        panel.background = ggplot2::element_rect(fill = NA, colour = NA),
        plot.margin = ggplot2::margin(2, 6, 2, 6)
      )
  }, res = 110)

  # your existing outputs here:
  output$cell_summary   <- renderUI({ tags$div("...") })

  # helper: pick a single display value (unique or collapse)
  pick_display <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(trimws(x))]
    if (length(x) == 0) return(NA_character_)
    ux <- unique(x)
    if (length(ux) == 1) ux else paste(ux, collapse = " | ")
  }

  selected_detections <- reactive({
    click <- input$map_shape_click
    sel_id <- selected_draw_id()

    # ---- drawn polygon: still needs spatial filtering ----
    if (!is.null(sel_id)) {

      pts <- species_sf_all_with_poly

      g <- selection_geom()

      if (is.null(g) || nrow(pts) == 0) {
        return(NULL)
      }

      keep <- tryCatch(
        within_any(pts, g),
        error = function(e) {
          showNotification(
            paste("Polygon selection failed:", e$message),
            type = "error"
          )
          return(rep(FALSE, nrow(pts)))
        }
      )

      return(pts[keep, , drop = FALSE])
    }

    if (is.null(click) || is.null(click$id)) {
      return(NULL)
    }

    id <- as.character(click$id)

    # ---- clicked MPA/AOI polygon: use precomputed lookup ----
    if (grepl("\\|\\|", id)) {
      parts <- strsplit(id, "\\|\\|")[[1]]
      p_type <- parts[1]
      p_name <- parts[2]

      pts <- species_sf_all_with_poly

      return(
        pts %>%
          dplyr::filter(site_type == p_type, site_name == p_name)
      )
    }

    # ---- clicked grid cell: use precomputed lookup ----
    cid <- suppressWarnings(as.integer(id))
    if (!is.na(cid)) {
      pts <- species_sf_all_with_cell

      return(
        pts %>%
          dplyr::filter(cell_id == cid)
      )
    }

    NULL
  })


  # ---- group filter helper (multi-select; union across selected groups) ----
  apply_group_filter <- function(occ_all, groups) {

    group_map <- list(
      "Fishes"        = list(col = "class",   vals = c("Teleostei")),
      "Sharks & Rays" = list(col = "class",   vals = c("Elasmobranchii")),
      "Mammals"       = list(col = "class",   vals = c("Mammalia")),
      "Turtles"      = list(col = "order",   vals = c("Testudines")),
      "Birds"         = list(col = "class",   vals = c("Aves")),
      "Molluscs"      = list(col = "phylum",  vals = c("Mollusca")),
      "Arthropods"    = list(col = "phylum",  vals = c("Arthropoda")),
      "Plants & Algae"        = list(
        list(col = "kingdom", vals = c("Plantae")),
        list(col = "class",   vals = c("Phaeophyceae"))
      )
    )

    groups <- as.character(groups %||% character(0))
    groups <- intersect(groups, names(group_map))
    if (length(groups) == 0) return(occ_all)

    cols_needed <- unique(vapply(group_map[groups], `[[`, character(1), "col"))
    missing_cols <- setdiff(cols_needed, names(occ_all))

    if (length(missing_cols) > 0) {
      # choose ONE behavior:
      # return(occ_all[0, , drop = FALSE])  # loud fail
      return(occ_all)                       # silent safe
    }

    keep <- rep(FALSE, nrow(occ_all))
    for (g in groups) {
      spec <- group_map[[g]]
      xcol <- tolower(as.character(occ_all[[spec$col]]))
      keep <- keep | (!is.na(xcol) & xcol %in% tolower(spec$vals))
    }

    occ_all[keep, , drop = FALSE]
  }

  # ---- active groups (multi-select) ----
  active_groups <- reactiveVal(character(0))

  toggle_group <- function(g) {
    cur <- active_groups()
    if (g %in% cur) active_groups(setdiff(cur, g)) else active_groups(c(cur, g))
  }

  observeEvent(input$total_fish,       { toggle_group("Fishes") }, ignoreInit = TRUE)
  observeEvent(input$total_sharks,     { toggle_group("Sharks & Rays") }, ignoreInit = TRUE)
  observeEvent(input$total_mammals,    { toggle_group("Mammals") }, ignoreInit = TRUE)
  observeEvent(input$total_reptiles,   { toggle_group("Turtles") }, ignoreInit = TRUE)
  observeEvent(input$total_birds,      { toggle_group("Birds") }, ignoreInit = TRUE)
  observeEvent(input$total_molluscs,   { toggle_group("Molluscs") }, ignoreInit = TRUE)
  observeEvent(input$total_arthropods, { toggle_group("Arthropods") }, ignoreInit = TRUE)
  observeEvent(input$total_plants,     { toggle_group("Plants & Algae") }, ignoreInit = TRUE)

  sync_group_button_classes <- function() {
    cur <- active_groups()

    btn_map <- c(
      total_fish       = "Fishes",
      total_sharks       = "Sharks & Rays",
      total_mammals    = "Mammals",
      total_reptiles   = "Turtles",
      total_birds      = "Birds",
      total_molluscs   = "Molluscs",
      total_arthropods = "Arthropods",
      total_plants     = "Plants & Algae"
    )

    for (btn_id in names(btn_map)) {
      g <- btn_map[[btn_id]]
      if (g %in% cur) shinyjs::addClass(btn_id, "btn-group-on") else shinyjs::removeClass(btn_id, "btn-group-on")
    }
  }

  # IMPORTANT: keep button visuals synced with state
  observeEvent(active_groups(), {
    sync_group_button_classes()
  }, ignoreInit = FALSE)


  # helper: current group label for display/debug if needed
  active_groups_label <- reactive({
    cur <- active_groups()
    if (length(cur) == 0) "All" else paste(cur, collapse = " + ")
  })

  # ---- IUCN/SARA/AIS sets + toggles (assumes IUCN, SARA & AIS exist globally) ----
  iucn_set <- reactive({
    df <- selection_map_df()

    if (is.null(df) || nrow(df) == 0 || !"category" %in% names(df)) {
      return(character(0))
    }

    df %>%
      dplyr::filter(!is.na(category), trimws(as.character(category)) != "") %>%
      dplyr::pull(scientificName) %>%
      unique() %>%
      na.omit() %>%
      as.character()
  })

  sara_set <- reactive(unique(na.omit(SARA$Scientific.Name)))
  ais_set  <- reactive(unique(na.omit(AIS$Scientific.Name)))  # adjust column if needed

  filter_iucn_on <- reactiveVal(FALSE)
  filter_sara_on <- reactiveVal(FALSE)
  filter_ais_on  <- reactiveVal(FALSE)

  confirmed_map_filters <- eventReactive(input$div_apply, {
    list(
      year = sel_year_chr(),
      groups = active_groups(),
      iucn_on = filter_iucn_on(),
      sara_on = filter_sara_on(),
      ais_on = filter_ais_on()
    )
  }, ignoreNULL = FALSE)

  apply_species_filters <- function(occ_all, filters = NULL) {

    if (is.null(filters)) {
      groups   <- active_groups()
      iucn_on  <- filter_iucn_on()
      sara_on  <- filter_sara_on()
      ais_on   <- filter_ais_on()
    } else {
      groups   <- filters$groups
      iucn_on  <- filters$iucn_on
      sara_on  <- filters$sara_on
      ais_on   <- filters$ais_on
    }

    if (!is.null(filters)) {
      yr <- filters$year

      if (!is.null(yr) && yr != "All" && "year" %in% names(occ_all)) {
        occ_all <- occ_all %>%
          dplyr::filter(as.character(year) == yr)
      }
    }

    occ_all <- apply_group_filter(occ_all, groups)

    if ("scientificName" %in% names(occ_all)) {

      spp_keep <- character(0)

      if (isTRUE(iucn_on)) spp_keep <- union(spp_keep, iucn_set())
      if (isTRUE(sara_on)) spp_keep <- union(spp_keep, sara_set())
      if (isTRUE(ais_on))  spp_keep <- union(spp_keep, ais_set())

      if (length(spp_keep) > 0) {
        occ_all <- occ_all %>%
          dplyr::filter(scientificName %in% spp_keep)
      }
    }

    occ_all
  }

  #IUCN icon for marker
  iucn_icon <- leaflet::makeIcon(
    iconUrl = "img/species_buttons/iucn-red-list-logo-red.png",
    iconWidth = 20,
    iconHeight = 20,
    iconAnchorX = 14,
    iconAnchorY = 14,
  )

  #SARA icon for marker
  sara_icon <- leaflet::makeIcon(
    iconUrl = "img/species_buttons/SARA_icon.png",
    iconWidth = 20,
    iconHeight = 20,
    iconAnchorX = 14,
    iconAnchorY = 14
  )

  #AIS icon for marker
  ais_icon <- leaflet::makeIcon(
    iconUrl = "img/species_buttons/invasive_symbol.png",
    iconWidth = 28,
    iconHeight = 28,
    iconAnchorX = 14,
    iconAnchorY = 14
  )

  # redraw points when year changes
  observeEvent(
    list(sel_year_chr(), sampling_points_layer_on(), filter_iucn_on(), filter_sara_on(), filter_ais_on()),
    {
      proxy <- leafletProxy("map")

      if (!isTRUE(sampling_points_layer_on())) {
        proxy %>% clearGroup("Sampling Points")
        return(NULL)
      }

      yr <- sel_year_chr()
      pts <- sampling_pts

      # Year filter ONLY
      if (yr != "All") {
        pts <- pts %>% dplyr::filter(as.character(year) == yr)
      }

      proxy %>% clearGroup("Sampling Points")

      if (isTRUE(filter_iucn_on())) {
        pts_iucn <- pts %>%
          dplyr::filter(!is.na(category), trimws(as.character(category)) != "")

        proxy %>%
          addMarkers(
            data = pts_iucn,
            icon = iucn_icon,
            group = "Sampling Points",
            options = markerOptions(pane = "pane_points"),
            label = ~paste0(
              scientificName,
              ifelse(is.na(category), "", paste0(" | IUCN: ", category)),
              ifelse(is.na(year), "", paste0(" | Year: ", year)),
              ifelse(is.na(samp_name), "", paste0(" | Sample: ", samp_name))
            )
          )
      }

      if (isTRUE(filter_sara_on())) {
        pts_sara <- pts %>%
          dplyr::filter(scientificName %in% sara_set())

        proxy %>%
          addMarkers(
            data = pts_sara,
            icon = sara_icon,
            group = "Sampling Points",
            options = markerOptions(pane = "pane_points"),
            label = ~paste0(
              scientificName,
              ifelse(is.na(year), "", paste0(" | Year: ", year)),
              ifelse(is.na(samp_name), "", paste0(" | Sample: ", samp_name))
            )
          )
      }

      if (isTRUE(filter_ais_on())) {
        pts_ais <- pts %>%
          dplyr::filter(scientificName %in% ais_set())

        proxy %>%
          addMarkers(
            data = pts_ais,
            icon = ais_icon,
            group = "Sampling Points",
            options = markerOptions(pane = "pane_points"),
            label = ~paste0(
              scientificName,
              ifelse(is.na(year), "", paste0(" | Year: ", year)),
              ifelse(is.na(samp_name), "", paste0(" | Sample: ", samp_name))
            )
          )
      }

      if (!isTRUE(filter_iucn_on()) && !isTRUE(filter_sara_on()) && !isTRUE(filter_ais_on())) {
        proxy %>%
          addCircleMarkers(
            data = pts,
            group = "Sampling Points",
            radius = 2,
            stroke = TRUE,
            weight = 1,
            opacity = 1,
            fillOpacity = 0.8,
            options = pathOptions(pane = "pane_points"),
            label = ~paste0(
              "Marker: ", target_gene,
              ifelse(is.na(year), "", paste0(" | Year: ", year)),
              ifelse(is.na(samp_name), "", paste0(" | Sample: ", samp_name))
            )
          )
      }
    },
    ignoreInit = TRUE
  )

  diversity_dropdown_data <- reactive({
    df <- dplyr::bind_rows(
      mpa_membership_panel_df(),
      polygon_membership_panel_df()
    )

    if (is.null(df) || nrow(df) == 0) {
      return(df)
    }

    df %>%
      add_primer_combo() %>%
      dplyr::mutate(
        target_gene  = as.character(target_gene),
        primer_combo = as.character(primer_combo)
      )
  })

  # ---- source data for diversity dropdowns ----
  observeEvent(diversity_dropdown_data(), {
    dd <- diversity_dropdown_data()
    if (is.null(dd) || nrow(dd) == 0) {
      updateSelectizeInput(
        session,
        "div_target_gene",
        choices = character(0),
        selected = character(0),
        server = TRUE
      )
      return()
    }
    gene_choices <- dd %>%
      dplyr::filter(!is.na(target_gene), target_gene != "") %>%
      dplyr::pull(target_gene) %>%
      unique() %>%
      sort()


    freezeReactiveValue(input, "div_target_gene")
    updateSelectizeInput(
      session  = session,
      inputId  = "div_target_gene",
      choices  = gene_choices,
      selected = gene_choices[1],
      server   = TRUE
    )
  }, ignoreInit = FALSE)

  primer_choices_reactive <- reactive({
    dd <- diversity_dropdown_data()
    if (is.null(dd) || nrow(dd) == 0) return(character(0))
    genes_selected <- input$div_target_gene %||% character(0)

    if (length(genes_selected) > 0) {
      dd <- dd %>% dplyr::filter(target_gene %in% genes_selected)
    } else {
      return(character(0))
    }

    dd %>%
      dplyr::filter(!is.na(primer_combo), primer_combo != "") %>%
      dplyr::pull(primer_combo) %>%
      unique() %>%
      sort()
  })

  # Primer choices
  observeEvent(
    list(primer_choices_reactive(), input$div_target_gene),
    {
      primer_choices <- primer_choices_reactive()



      freezeReactiveValue(input, "div_primer")
      updateSelectizeInput(
        session  = session,
        inputId  = "div_primer",
        choices  = primer_choices,
        selected = primer_choices,
        server   = TRUE
      )
    },
    ignoreInit = FALSE
  )

  observeEvent(input$div_primer_all, {
    primer_choices <- isolate(primer_choices_reactive())

    updateSelectizeInput(
      session  = session,
      inputId  = "div_primer",
      choices  = primer_choices,
      selected = primer_choices,
      server   = TRUE
    )
  })

  observeEvent(input$div_primer_none, {
    primer_choices <- isolate(primer_choices_reactive())

    updateSelectizeInput(
      session  = session,
      inputId  = "div_primer",
      choices  = primer_choices,
      selected = character(0),
      server   = TRUE
    )
  })

  observeEvent(compare_polygon_choices(), {
    choices <- compare_polygon_choices()

    cur_sel <- isolate(input$div_compare_polygons %||% character(0))
    sel_keep <- intersect(cur_sel, choices)

    freezeReactiveValue(input, "div_compare_polygons")

    updateSelectizeInput(
      session  = session,
      inputId  = "div_compare_polygons",
      choices  = choices,
      selected = sel_keep,
      server   = TRUE
    )
  }, ignoreInit = FALSE)

  #taxonomic selection reactive #NEW
  active_tax_rank <- reactive({
    req(div_unlocked())
    input$tax_rank %||% "scientificName"
  })

  observeEvent(input$div_apply, {
    div_filters(
      list(
        target_gene = input$div_target_gene %||% character(0),
        primers     = input$div_primer %||% character(0),
        polygons      = input$div_compare_polygons %||% character(0)
      )
    )
  }, ignoreInit = FALSE)

  add_polygon_selection <- function(pts_sf) {
    if (is.null(pts_sf) || nrow(pts_sf) == 0) {
      pts_sf$polygon_selection <- character(0)
      return(pts_sf)
    }

    poly_labels <- rep(NA_character_, nrow(pts_sf))

    # existing MPA/AOI polygons
    if (!is.null(all_polys_click) && nrow(all_polys_click) > 0) {
      mpa_hits <- sf::st_intersects(pts_sf, all_polys_click)

      mpa_labels <- vapply(seq_along(mpa_hits), function(i) {
        hit <- mpa_hits[[i]]
        if (length(hit) == 0) return(NA_character_)

        labs <- all_polys_click$site_name[hit]

        paste(unique(labs), collapse = " | ")
      }, character(1))

      poly_labels <- mpa_labels
    }

    # drawn polygons
    polys_drawn <- drawn_polys()
    if (!is.null(polys_drawn) && nrow(polys_drawn) > 0) {
      drawn_hits <- sf::st_intersects(pts_sf, polys_drawn)

      drawn_labels <- vapply(seq_along(drawn_hits), function(i) {
        hit <- drawn_hits[[i]]
        if (length(hit) == 0) return(NA_character_)

        labs <- polys_drawn$draw_label[hit]
        paste(unique(labs), collapse = " | ")
      }, character(1))

      poly_labels <- ifelse(
        !is.na(poly_labels) & !is.na(drawn_labels),
        paste(poly_labels, drawn_labels, sep = " | "),
        dplyr::coalesce(poly_labels, drawn_labels)
      )
    }

    pts_sf$polygon_selection <- poly_labels
    pts_sf
  }

  species_list_occ_all <- reactive({
    occ_all <- selection_panel_df()
    req(occ_all)

    occ_all %>%
      dplyr::distinct(scientificName, .keep_all = TRUE) %>%
      dplyr::arrange(dplyr::coalesce(worms_valid_name, scientificName))
  })


  observeEvent(input$IUCN, {
    new_state <- !isTRUE(filter_iucn_on())
    filter_iucn_on(new_state)
    if (new_state) shinyjs::addClass("IUCN", "btn-iucn-on") else shinyjs::removeClass("IUCN", "btn-iucn-on")
  })

  observeEvent(input$SARA, {
    new_state <- !isTRUE(filter_sara_on())
    filter_sara_on(new_state)
    if (new_state) shinyjs::addClass("SARA", "btn-sara-on") else shinyjs::removeClass("SARA", "btn-sara-on")
  })

  observeEvent(input$AIS, {
    new_state <- !isTRUE(filter_ais_on())
    filter_ais_on(new_state)
    if (new_state) shinyjs::addClass("AIS", "btn-ais-on") else shinyjs::removeClass("AIS", "btn-ais-on")
  })

  download_ready <- reactiveVal(FALSE)

  observeEvent(input$div_apply, {
    download_ready(TRUE)
    shinyjs::enable("downloadData")
  }, ignoreInit = TRUE)

  observeEvent(
    active_groups(),
    {
      download_ready(FALSE)
      shinyjs::disable("downloadData")
    },
    ignoreInit = TRUE
  )

  # UI

  apply_interest_filter <- function(spp_vec) {
    if (!isTRUE(filter_iucn_on()) && !isTRUE(filter_sara_on()) && !isTRUE(filter_ais_on())) {
      return(spp_vec)
    }

    keep <- character(0)

    if (isTRUE(filter_iucn_on())) keep <- union(keep, iucn_set())
    if (isTRUE(filter_sara_on())) keep <- union(keep, sara_set())
    if (isTRUE(filter_ais_on()))  keep <- union(keep, ais_set())

    intersect(spp_vec, keep)
  }

  active_filters_label <- reactive({
    labs <- c()
    if (isTRUE(filter_iucn_on())) labs <- c(labs, "IUCN")
    if (isTRUE(filter_sara_on())) labs <- c(labs, "SARA")
    if (isTRUE(filter_ais_on()))  labs <- c(labs, "AIS")
    if (length(labs) == 0) "None" else paste(labs, collapse = " + ")
  })

  output$iucn_details <- DT::renderDT({
    if (!isTRUE(filter_iucn_on())) {
      return(DT::datatable(
        data.frame(Message = "Select the “IUCN” button to view species at risk details."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det <- selection_panel_df()
    if (is.null(det) || nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections in the current selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det <- det %>%
      dplyr::mutate(scientificName = as.character(scientificName)) %>%
      dplyr::filter(!is.na(category), trimws(as.character(category)) != "")

    if (nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No IUCN detections for this selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    out <- det %>%
      dplyr::group_by(scientificName, category) %>%
      dplyr::summarise(
        # Common.Name  = paste(sort(unique(na.omit(Common.Name))), collapse = " |OR| "),
        n_detections = dplyr::n(),
        n_samples    = dplyr::n_distinct(samp_name),
        samples      = paste(sort(unique(samp_name)), collapse = ", "),
        years        = paste(sort(unique(na.omit(year))), collapse = ", "),
        months       = paste(sort(unique(na.omit(month))), collapse = ", "),
        markers      = paste(sort(unique(na.omit(target_gene))), collapse = ", "),
        primers      = paste(
          sort(unique(na.omit(
            paste(
              trimws(pcr_primer_name_forward),
              trimws(pcr_primer_name_reverse),
              sep = " | "
            )
          ))),
          collapse = ", "
        ),
        .groups = "drop"
      ) %>%
      dplyr::select(
        scientificName,
        # Common.Name,
        category,
        n_detections,
        n_samples,
        samples,
        years,
        months,
        markers,
        primers
      ) %>%
      dplyr::arrange(category, scientificName)

    DT::datatable(
      out,
      rownames = FALSE,
      colnames = c(
        "Species",
        # "Common Name",
        "IUCN Red List Category",
        "Number of Detections",
        "Number of Samples",
        "Samples",
        "Years",
        "Months",
        "Target Genes",
        "Primers"
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = FALSE,
        scrollCollapse = TRUE
      ),
      class = "nowrap"
    )
  })

  output$sara_details <- DT::renderDT({
    if (!isTRUE(filter_sara_on())) {
      return(DT::datatable(
        data.frame(Message = "Select the “SARA” button to view species at risk details."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det <- selection_panel_df()
    if (is.null(det) || nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections in the current selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det <- det %>%
      dplyr::mutate(scientificName = as.character(scientificName)) %>%
      dplyr::filter(scientificName %in% sara_set())

    if (nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No SARA detections for this selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det2 <- det %>%
      dplyr::left_join(
        SARA %>% dplyr::select(
          Scientific.Name,
          Common.Name,
          Rating
        ),
        by = c("scientificName" = "Scientific.Name")
      )

    out <- det2 %>%
      dplyr::group_by(scientificName, Rating) %>%
      dplyr::summarise(
        Common.Name  = paste(sort(unique(na.omit(Common.Name))), collapse = " |OR| "),
        n_detections = dplyr::n(),
        n_samples    = dplyr::n_distinct(samp_name),
        samples      = paste(sort(unique(samp_name)), collapse = ", "),
        years        = paste(sort(unique(na.omit(year))), collapse = ", "),
        months       = paste(sort(unique(na.omit(month))), collapse = ", "),
        markers      = paste(sort(unique(na.omit(target_gene))), collapse = ", "),
        primers      = paste(
          sort(unique(na.omit(
            paste(
              trimws(pcr_primer_name_forward),
              trimws(pcr_primer_name_reverse),
              sep = " | "
            )
          ))),
          collapse = ", "
        ),
        .groups = "drop"
      ) %>%
      dplyr::select(
        scientificName,
        Common.Name,
        Rating,
        n_detections,
        n_samples,
        samples,
        years,
        months,
        markers,
        primers
      ) %>%
      dplyr::arrange(Rating, scientificName)

    DT::datatable(
      out,
      rownames = FALSE,
      colnames = c(
        "Species",
        "Common Name",
        "SARA Rating",
        "Number of Detections",
        "Number of Samples",
        "Samples",
        "Years",
        "Months",
        "Target Genes",
        "Primers"
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = FALSE,
        scrollCollapse = TRUE
      ),
      class = "nowrap"
    )
  })

  output$ais_details <- DT::renderDT({
    if (!isTRUE(filter_ais_on())) {
      return(DT::datatable(
        data.frame(Message = "Select the “AIS” button to view invasive species details."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det <- selection_panel_df()
    if (is.null(det) || nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections in the current selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det <- det %>%
      dplyr::mutate(scientificName = as.character(scientificName)) %>%
      dplyr::filter(scientificName %in% ais_set())

    if (nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No AIS detections for this selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det2 <- det %>%
      dplyr::left_join(
        AIS %>% dplyr::select(Scientific.Name, Common.Name, Type),
        by = c("scientificName" = "Scientific.Name")
      )

    out <- det2 %>%
      dplyr::group_by(scientificName) %>%
      dplyr::summarise(
        Common.Name   = paste(sort(unique(na.omit(Common.Name))), collapse = " |OR| "),
        Type          = paste(sort(unique(na.omit(Type))), collapse = " |OR| "),
        n_detections  = dplyr::n(),
        n_samples     = dplyr::n_distinct(samp_name),
        samples       = paste(sort(unique(samp_name)), collapse = ", "),
        years         = paste(sort(unique(na.omit(year))), collapse = ", "),
        months        = paste(sort(unique(na.omit(month))), collapse = ", "),
        markers       = paste(sort(unique(na.omit(target_gene))), collapse = ", "),
        primers       = paste(
          sort(unique(na.omit(
            paste(
              trimws(pcr_primer_name_forward),
              trimws(pcr_primer_name_reverse),
              sep = " | "
            )
          ))),
          collapse = ", "
        ),
        .groups = "drop"
      ) %>%
      dplyr::select(
        scientificName,
        Common.Name,
        Type,
        n_detections,
        n_samples,
        samples,
        years,
        months,
        markers,
        primers
      ) %>%
      dplyr::arrange(scientificName)

    DT::datatable(
      out,
      rownames = FALSE,
      colnames = c(
        "Species",
        "Common Name",
        "Type",
        "Number of Detections",
        "Number of Samples",
        "Samples",
        "Years",
        "Months",
        "Target Genes",
        "Primers"
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = FALSE,
        scrollCollapse = TRUE
      ),
      class = "nowrap"
    )
  })

  output$download_hint <- renderText({
    if (!download_ready()) {
      "Click Confirm before downloading and to update diversity plots."
    } else {
      "Filters confirmed. Ready to download."
    }
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("GOTeDNA_filtered_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      tryCatch({
        message("Download started")

        filters <- div_filters()
        message("div_filters read")

        yr <- sel_year_chr()
        message("Year: ", yr)

        dat <- species_sf_all
        message("species_sf_all rows: ", nrow(dat))

        if (yr != "All") {
          dat <- dat %>%
            dplyr::filter(as.character(year) == yr)
          message("After year filter rows: ", nrow(dat))
        }

        dat <- apply_species_filters(dat)
        message("After floating panel filters rows: ", nrow(dat))

        dat <- apply_diversity_dropdown_filters(dat, filters)
        message("After diversity filters rows: ", nrow(dat))

        message("Adding polygon selection")
        dat <- add_polygon_selection(dat)
        message("Polygon selection added")

        if (length(filters$polygons) > 0) {
          dat <- dat %>%
            dplyr::filter(!is.na(polygon_selection)) %>%
            dplyr::filter(
              vapply(
                strsplit(as.character(polygon_selection), " \\| "),
                function(x) any(x %in% filters$polygons),
                logical(1)
              )
            )

          message("After polygon comparison filter rows: ", nrow(dat))
        }

        if (is.null(dat) || nrow(dat) == 0) {
          message("No rows left; writing empty CSV")
          utils::write.csv(data.frame(), file, row.names = FALSE, na = "")
          return()
        }

        message("Extracting coordinates")
        coords <- sf::st_coordinates(dat)
        dat$decimalLongitude <- coords[, 1]
        dat$decimalLatitude  <- coords[, 2]
        message("Coordinates added")

        dat <- dat %>%
          sf::st_drop_geometry()
        message("Geometry dropped")

        dat <- add_primer_combo(dat)
        message("Primer combo added")

        preferred_cols <- c(
          "occurrenceID", "scientificName",
          "kingdom", "phylum", "class", "order", "family", "genus",
          "samp_name", "year", "month", "eventDate",
          "target_gene", "primer_combo",
          "organismQuantity", "organismQuantityType",
          "decimalLatitude", "decimalLongitude", "polygon_selection",
          "minimumDepthInMeters", "maximumDepthInMeters",
          "samp_size", "samp_size_unit"
        )

        keep_cols <- c(
          intersect(preferred_cols, names(dat)),
          setdiff(names(dat), preferred_cols)
        )
        message("Selecting columns")

        dat <- dat %>%
          dplyr::select(dplyr::all_of(keep_cols))
        message("Columns selected")

        utils::write.csv(dat, file, row.names = FALSE, na = "")
        message("Download written successfully")

      }, error = function(e) {
        message("DOWNLOAD ERROR CLASS: ", paste(class(e), collapse = ", "))
        message("DOWNLOAD ERROR MESSAGE: ", conditionMessage(e))
        traceback(2)
        stop(e)
      })
    }
  )

  poly_bbox <- sf::st_bbox(sf::st_buffer(sf::st_union(all_polys_click), dist = 0.2))

  # ---- 1) Render the leaflet map ONCE ----
  output$map <- renderLeaflet({

    # ---- choose initial year + initial layers safely ----
    yrs <- sort(unique(na.omit(as.character(KEY_TBL$year))))

    init_12S <- if (default_year == "All") RICHNESS_GENE_ALL[["12S"]] else RICHNESS_BY_KEY[[paste0("12S_", default_year)]]
    init_COI <- if (default_year == "All") RICHNESS_GENE_ALL[["COI"]] else RICHNESS_BY_KEY[[paste0("COI_", default_year)]]
    init_16S <- if (default_year == "All") RICHNESS_GENE_ALL[["16S"]] else RICHNESS_BY_KEY[[paste0("16S_", default_year)]]
    init_18S <- if (default_year == "All") RICHNESS_GENE_ALL[["18S"]] else RICHNESS_BY_KEY[[paste0("18S_", default_year)]]

    init_ALL <- {
      if (default_year == "All") {
        RICHNESS_ALL
      } else if (default_year %in% names(RICHNESS_ALL_BY_YEAR)) {
        RICHNESS_ALL_BY_YEAR[[default_year]]
      } else {
        RICHNESS_ALL
      }
    }

    # ---- defensive: ensure has_sampling exists where you use it ----
    if (!is.null(init_12S) && !"has_sampling" %in% names(init_12S)) init_12S$has_sampling <- FALSE
    if (!is.null(init_COI) && !"has_sampling" %in% names(init_COI)) init_COI$has_sampling <- FALSE
    if (!is.null(init_16S) && !"has_sampling" %in% names(init_16S)) init_16S$has_sampling <- FALSE
    if (!is.null(init_18S) && !"has_sampling" %in% names(init_18S)) init_18S$has_sampling <- FALSE
    if (!is.null(init_ALL) && !"has_sampling" %in% names(init_ALL)) init_ALL$has_sampling <- FALSE

    m <- leaflet() %>%
      addMapPane("pane_polys",      zIndex = 400) %>%
      addMapPane("pane_zones",      zIndex = 410) %>%
      addMapPane("pane_poly_total", zIndex = 420) %>%
      addMapPane("pane_points",     zIndex = 430) %>%
      addMapPane("pane_drawn_top",    zIndex = 900) %>%   #new code
      addMapPane("pane_selected_top", zIndex = 950) %>%   #new code
      addProviderTiles(providers$CartoDB.Positron, group = "CartoDB Positron") %>%
      addProviderTiles(providers$Esri.OceanBasemap, group = "Esri Ocean Basemap") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery") %>%
      fitBounds(
        lng1 = -65, lat1 = 41,
        lng2 = -59, lat2 = 47
      )

    # ---- add initial richness layers if they exist ----
    if (!is.null(init_12S) && nrow(init_12S) > 0) {
      m <- m %>% addPolygons(
        data        = init_12S,
        group       = "12S",
        layerId     = ~cell_id,
        fillColor   = ~ifelse(has_sampling, pal_rich(n_species), "#D6F4FF"),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling, paste("12S richness:", n_species), "No sampling in this cell"),
        options     = pathOptions(pane = "pane_polys")
      )
    }

    if (!is.null(init_COI) && nrow(init_COI) > 0) {
      m <- m %>% addPolygons(
        data        = init_COI,
        group       = "COI",
        layerId     = ~cell_id,
        fillColor   = ~ifelse(has_sampling, pal_rich(n_species), "#D6F4FF"),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling, paste("COI richness:", n_species), "No sampling in this cell"),
        options     = pathOptions(pane = "pane_polys")
      )
    }

    if (!is.null(init_16S) && nrow(init_16S) > 0) {
      m <- m %>% addPolygons(
        data        = init_16S,
        group       = "16S",
        layerId     = ~cell_id,
        fillColor   = ~ifelse(has_sampling, pal_rich(n_species), "#D6F4FF"),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling, paste("16S richness:", n_species), "No sampling in this cell"),
        options     = pathOptions(pane = "pane_polys")
      )
    }

    if (!is.null(init_18S) && nrow(init_18S) > 0) {
      m <- m %>% addPolygons(
        data        = init_18S,
        group       = "18S",
        layerId     = ~cell_id,
        fillColor   = ~ifelse(has_sampling, pal_rich(n_species), "#D6F4FF"),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling, paste("18S richness:", n_species), "No sampling in this cell"),
        options     = pathOptions(pane = "pane_polys")
      )
    }
    # Always add ALL
    m <- m %>% addPolygons(
      data        = init_ALL,
      group       = "All",
      layerId     = ~cell_id,
      fillColor   = ~pal_rich(n_species_total),
      fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
      color       = NA,
      label       = ~ifelse(has_sampling, paste("Total richness:", n_species_total), "No sampling in this cell"),
      options     = pathOptions(pane = "pane_polys"),
      highlightOptions = highlightOptions(weight = 2, bringToFront = TRUE)
    )

    # ---- rest of your map layers ----
    m <- m %>%
      addPolygons(
        data = depth_layers[["All"]],
        group = "Sampling Depth",
        layerId = ~cell_id,
        fillColor   = ~final_fill,
        fillOpacity = ~final_alpha,
        opacity     = 1,
        color       = NA,
        label       = ~paste0(
          "Depth (m, median): ", ifelse(is.na(depth_med), "NA", round(depth_med, 1)), "<br>",
          "Min: ", ifelse(is.na(depth_min), "NA", round(depth_min, 1)), " | ",
          "Max: ", ifelse(is.na(depth_max), "NA", round(depth_max, 1)), "<br>"
        ) %>% lapply(htmltools::HTML),
        highlightOptions = leaflet::highlightOptions(weight = 2, bringToFront = TRUE),
        options = pathOptions(pane = "pane_polys")
      ) %>%
      addCircleMarkers(
        data        = sampling_pts,
        group       = "Sampling Points",
        radius      = 2,
        stroke      = TRUE,
        weight      = 1,
        opacity     = 1,
        fillOpacity = 0.8,
        options = pathOptions(pane = "pane_points"),
        label       = ~paste0(
          "Marker: ", target_gene,
          ifelse(is.na(year), "", paste0(" | Year: ", year)),
          #ifelse(is.na(source_file), "", paste0(" | File: ", source_file)),
          ifelse(is.na(samp_name), "", paste0(" | Sample: ", samp_name))
        )
      ) %>%
      addPolygons(
        data        = all_polys_zones,
        group       = "MPA/AOI Zone Boundaries",
        fillOpacity = 0,
        color       = "black",
        weight      = 1,
        opacity     = 0.8,
        options     = pathOptions(clickable = FALSE, pane = "pane_zones")
      ) %>%
      addPolygons(
        data        = all_polys_click,
        group       = "MPA/AOI total species richness",
        layerId     = ~paste(site_type, site_name, sep="||"),
        fillOpacity = 0.05,
        color       = "black",
        weight      = 2,
        opacity     = 1,
        popup       = ~site_name,
        options     = pathOptions(pane = "pane_poly_total"),
        highlightOptions = highlightOptions(weight = 3, bringToFront = TRUE)
      )

    # ---- legends + controls ----
    m <- m %>%
      leaflet::addLegend(
        position = "bottomright",
        pal      = pal_rich,
        values   = init_ALL$n_species_total %||% numeric(0),
        title    = "Species richness from eDNA",
        opacity  = 1,
        className = "legend-base legend-richness-box"
      ) %>%
      leaflet::addLegend(
        position  = "bottomright",
        colors    = c(depth_legend_cols, "#feb24c"),
        labels    = c(depth_legend_labs, "Mixed depth (orange overlay)"),
        title     = "Sampling depth",
        opacity   = 1,
        className = "legend-base legend-depth-box"
      ) %>%
      addCircleMarkers(
        lng = 0, lat = 0,
        radius = 1,
        opacity = 0,
        fillOpacity = 0,
        stroke = FALSE,
        group = "Monthly Sampling Plot"
      ) %>%
      addDrawToolbar(
        targetGroup = "drawn",
        polygonOptions = drawPolygonOptions(showArea = TRUE),
        rectangleOptions = drawRectangleOptions(),
        polylineOptions = FALSE,
        markerOptions = FALSE,
        circleOptions        = FALSE,
        circleMarkerOptions  = FALSE,
        editOptions = editToolbarOptions()
      ) %>%
      addLayersControl(
        baseGroups = c("CartoDB Positron", "Esri Ocean Basemap", "Esri World Imagery"),
        overlayGroups = c(
          "MPA/AOI total species richness",
          "All", "12S", "COI", "16S", "18S",
          "MPA/AOI Zone Boundaries",
          "Sampling Points",
          "Sampling Depth",
          "Monthly Sampling Plot"
        ),
        options = layersControlOptions(collapsed = FALSE)
        # )
        # addLayersControl(
        #   baseGroups = c("CartoDB Positron", "Esri Ocean Basemap", "Esri World Imagery"),
        #   overlayGroups = c(
        #     "MPA/AOI total species richness",
        #     "All", "12S", "COI", "16S", "18S",
        #     "MPA/AOI Zone Boundaries",
        #     "Sampling Points",
        #     "Sampling Depth"
        #   ),
        #   options = layersControlOptions(collapsed = FALSE)
      ) %>%
      htmlwidgets::onRender("
                      function(el, x){

                        const richness = new Set([
                          'All',
                          '12S',
                          'COI',
                          '16S',
                          '18S'
                        ]);

                        const context = new Set([
                          'MPA/AOI Zone Boundaries',
                          'Sampling Points',
                          'Sampling Depth',
                          'Monthly Sampling Plot'
                        ]);

                        const depthName = 'Sampling Depth';

                        const monthlyPlotName = 'Monthly Sampling Plot';

                        function getOverlayRows(){
                          const ctrl = el.querySelector('.leaflet-control-layers');
                          if(!ctrl) return [];
                          const rows = ctrl.querySelectorAll('.leaflet-control-layers-overlays label');
                          return Array.from(rows);
                        }

                        function labelText(row){
                          return row.textContent.replace(/\\s+/g,' ').trim();
                        }

                        function inputOf(row){
                          return row.querySelector('input[type=checkbox]');
                        }

                        function clickOffByName(name){
                          const rows = getOverlayRows();
                          rows.forEach(r => {
                            if(labelText(r) === name){
                              const cb = inputOf(r);
                              if(cb && cb.checked) cb.click();
                            }
                          });
                        }

                        function clickOffSet(nameSet){
                          const rows = getOverlayRows();
                          rows.forEach(r => {
                            const nm = labelText(r);
                            if(nameSet.has(nm)){
                              const cb = inputOf(r);
                              if(cb && cb.checked) cb.click();
                            }
                          });
                        }

                        function anyChecked(nameSet){
                          const rows = getOverlayRows();
                          for(const r of rows){
                            const nm = labelText(r);
                            const cb = inputOf(r);
                            if(cb && cb.checked && nameSet.has(nm)) return true;
                          }
                          return false;
                        }

                        function isChecked(name){
                          const rows = getOverlayRows();
                          for(const r of rows){
                            if(labelText(r) === name){
                              const cb = inputOf(r);
                              return cb ? cb.checked : false;
                            }
                          }
                          return false;
                        }

                       function updateMonthlyPlotVisibility(){
  const plotBox = document.getElementById('monthly_plot_control');
  if(!plotBox) return;

  const monthlyOn = isChecked(monthlyPlotName);

  const richLegend  = el.querySelector('.legend-richness-box');
  const depthLegend = el.querySelector('.legend-depth-box');

  let activeLegend = null;
  if(richLegend && !richLegend.classList.contains('legend-hidden')){
    activeLegend = richLegend;
  } else if(depthLegend && !depthLegend.classList.contains('legend-hidden')){
    activeLegend = depthLegend;
  }

  const pad = 14;
  const gap = 14;
  const defaultLegendRight = 10;
  const defaultLegendBottom = 10;

  // Monthly plot OFF
  if(!monthlyOn){
    plotBox.style.display = 'none';
    plotBox.style.left = 'auto';
    plotBox.style.right = 'auto';
    plotBox.style.top = 'auto';
    plotBox.style.bottom = 'auto';

    if(richLegend){
      richLegend.style.left = 'auto';
      richLegend.style.right = defaultLegendRight + 'px';
      richLegend.style.top = 'auto';
      richLegend.style.bottom = defaultLegendBottom + 'px';
      richLegend.style.margin = '0px';
    }

    if(depthLegend){
      depthLegend.style.left = 'auto';
      depthLegend.style.right = defaultLegendRight + 'px';
      depthLegend.style.top = 'auto';
      depthLegend.style.bottom = defaultLegendBottom + 'px';
      depthLegend.style.margin = '0px';
    }

    return;
  }

  // Monthly plot ON
  plotBox.style.display = 'block';
  plotBox.style.left = 'auto';
  plotBox.style.right = pad + 'px';
  plotBox.style.top = 'auto';
  plotBox.style.bottom = pad + 'px';
  plotBox.style.margin = '0px';

  if(activeLegend){
    const plotW = plotBox.offsetWidth || 320;

    activeLegend.style.left = 'auto';
    activeLegend.style.right = (plotW + gap + pad) + 'px';
    activeLegend.style.top = 'auto';
    activeLegend.style.bottom = pad + 'px';
    activeLegend.style.margin = '0px';
  }
}

                        function insertHeadings(){
  const ctrl = el.querySelector('.leaflet-control-layers');
  if(!ctrl) return;

  const overlayBox = ctrl.querySelector('.leaflet-control-layers-overlays');
  if(!overlayBox) return;

  const rows = getOverlayRows();
  if(rows.length === 0) return;

  const richnessOrder = ['All', '12S', '16S', 'COI', '18S'];
  const richnessRows = {};
  let firstCtx = null;

  rows.forEach(r => {
    const nm = labelText(r);
    if(richness.has(nm)){
      richnessRows[nm] = r;
    }
    if(!firstCtx && context.has(nm)){
      firstCtx = r;
    }
  });

  // If our custom layout already exists, do nothing
  if(overlayBox.querySelector('.richness-grid')){
    return;
  }

  const h1 = document.createElement('div');
  h1.className = 'layers-heading';
  h1.textContent = 'Species richness by gene region (grid cells)';

  const grid = document.createElement('div');
  grid.className = 'richness-grid';

  if(richnessRows['All']) grid.appendChild(richnessRows['All']);

  const spacer = document.createElement('div');
  spacer.className = 'richness-spacer';
  grid.appendChild(spacer);

  if(richnessRows['12S']) grid.appendChild(richnessRows['12S']);
  if(richnessRows['16S']) grid.appendChild(richnessRows['16S']);
  if(richnessRows['COI']) grid.appendChild(richnessRows['COI']);
  if(richnessRows['18S']) grid.appendChild(richnessRows['18S']);

  overlayBox.insertBefore(h1, firstCtx || null);
  overlayBox.insertBefore(grid, firstCtx || null);

  if(firstCtx && !overlayBox.querySelector('.layers-sep')){
    const sep = document.createElement('div');
    sep.className = 'layers-sep';
    overlayBox.insertBefore(sep, firstCtx);

    const h2 = document.createElement('div');
    h2.className = 'layers-heading';
    h2.textContent = 'Accessory Layers';
    overlayBox.insertBefore(h2, firstCtx);
  }
}

                        // ---- LEGEND SWAP ---- //

                        function updateLegends(){
                        const richLegend  = el.querySelector('.legend-richness-box');
                        const depthLegend = el.querySelector('.legend-depth-box');

                        const depthOn = isChecked(depthName);
                        const anyRichOn = anyChecked(richness);

                        if(depthLegend){
                        depthLegend.classList.toggle('legend-hidden', !depthOn);
                        }

                        if(richLegend){
                        richLegend.classList.toggle('legend-hidden', depthOn || !anyRichOn);
                        }

                        updateMonthlyPlotVisibility();
                        }

                        function wireExclusivity(){
                          const rows = getOverlayRows();
                          rows.forEach(r => {
                            const nm = labelText(r);
                            const cb = inputOf(r);
                            if(!cb) return;

                            // prevent duplicate listeners if MutationObserver fires
                            if(cb.dataset.wired === '1') return;
                            cb.dataset.wired = '1';

                            cb.addEventListener('change', function(){

                              // Richness ON -> Depth OFF
                              if(this.checked && richness.has(nm)){
                                if(isChecked(depthName)) clickOffByName(depthName);
                              }

                              // Depth ON -> all Richness OFF
                              if(this.checked && nm === depthName){
                                clickOffSet(richness);
                              }

                              updateLegends();
                            });
                          });
                        }

                        // ---- FORCE DEPTH + MONTHLY SAMPLING OFF AT STARTUP ----

                        function forceStartupOverlayState(){
                          clickOffByName(depthName);
                          clickOffByName(monthlyPlotName);
                          updateLegends();
                          updateMonthlyPlotVisibility();
                          }

                        insertHeadings();
                        wireExclusivity();

                        // ensure initial state after the control fully exists
                        setTimeout(forceStartupOverlayState, 0);


                       let observerBusy = false;

const obs = new MutationObserver(() => {
  if(observerBusy) return;

  observerBusy = true;

  try{
    insertHeadings();
    wireExclusivity();
    updateLegends();
  } finally {
    observerBusy = false;
  }
});
                        const ctrl = el.querySelector('.leaflet-control-layers');
                        if(ctrl){
                        obs.observe(ctrl, { childList: true, subtree: true });
                        }

                        setTimeout(function(){
                        $('#map_loading_overlay').fadeOut(250);
                        }, 300);
                      }
                      ")

    m

  })

  observe({
    polys <- drawn_polys()
    proxy <- leafletProxy("map")

    proxy %>% clearGroup("Drawn polygons")
    proxy %>% clearGroup("Selected drawn polygon")

    if (nrow(polys) == 0) return()

    # all drawn polygons
    proxy %>%
      addPolygons(
        data        = polys,
        group       = "Drawn polygons",
        layerId     = ~draw_id,
        fillColor   = "#2241a7",
        fillOpacity = 0.10,
        color       = "#2241a7",
        weight      = 2,
        opacity     = 1,
        label       = ~draw_label,
        options     = pathOptions(pane = "pane_drawn_top")
      )

    # selected polygon highlighted separately
    sid <- selected_draw_id()
    if (!is.null(sid) && sid %in% polys$draw_id) {
      sel <- polys %>% dplyr::filter(draw_id == sid)

      proxy %>%
        addPolygons(
          data        = sel,
          group       = "Selected drawn polygon",
          layerId     = ~draw_id,
          fillColor   = "#Fdd262",
          fillOpacity = 0.18,
          color       = "#Fdd262",
          weight      = 4,
          opacity     = 1,
          label       = ~draw_label,
          options     = pathOptions(pane = "pane_selected_top")
        )
    }
  })

  # ---- 2) When the year changes, swap the richness layers ----
  observeEvent(sel_year_chr(), {
    yr <- sel_year_chr()
    proxy <- leafletProxy("map")

    add_grid_layer <- function(data_sf, group_name, value_col, label_prefix) {
      if (is.null(data_sf) || nrow(data_sf) == 0) {
        proxy %>% clearGroup(group_name)
        return(invisible(NULL))
      }
      if (!"has_sampling" %in% names(data_sf)) data_sf$has_sampling <- FALSE
      data_sf$has_sampling <- as.logical(data_sf$has_sampling)
      if (!value_col %in% names(data_sf)) data_sf[[value_col]] <- NA_real_

      data_sf$.val   <- data_sf[[value_col]]
      data_sf$.label <- ifelse(
        data_sf$has_sampling,
        paste0(label_prefix, ": ", data_sf$.val),
        "No sampling in this cell"
      )

      proxy %>%
        clearGroup(group_name) %>%
        addPolygons(
          data        = data_sf,
          group       = group_name,
          layerId     = ~cell_id,
          fillColor   = ~pal_rich(.val),
          fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
          color       = NA,
          label       = ~.label,
          options     = pathOptions(pane = "pane_polys"),
          highlightOptions = highlightOptions(weight = 2, bringToFront = TRUE)
        )
    }

    grid12  <- if (yr == "All") RICHNESS_GENE_ALL[["12S"]] else RICHNESS_BY_KEY[[paste0("12S_", yr)]]
    gridCOI <- if (yr == "All") RICHNESS_GENE_ALL[["COI"]] else RICHNESS_BY_KEY[[paste0("COI_", yr)]]
    grid16  <- if (yr == "All") RICHNESS_GENE_ALL[["16S"]] else RICHNESS_BY_KEY[[paste0("16S_", yr)]]
    grid18  <- if (yr == "All") RICHNESS_GENE_ALL[["18S"]] else RICHNESS_BY_KEY[[paste0("18S_", yr)]]

    add_grid_layer(grid12,  "12S", "n_species", "12S richness")
    add_grid_layer(gridCOI, "COI", "n_species", "COI richness")
    add_grid_layer(grid16,  "16S", "n_species", "16S richness")
    add_grid_layer(grid18,  "18S", "n_species", "18S richness")

    gridALL <- if (yr == "All") RICHNESS_ALL else RICHNESS_ALL_BY_YEAR[[yr]] %||% RICHNESS_ALL
    add_grid_layer(gridALL, "All", "n_species_total", "Total richness")
  }, ignoreInit = TRUE)

  #Depth toggle
  selected_depth_layer <- function(year) {
    if (identical(year, "All")) return("All")
    as.character(year)
  }

  observe({
    req(input$sel_year)

    yr_key <- selected_depth_layer(input$sel_year)
    if (!yr_key %in% names(depth_layers)) return()

    occ_all <- depth_layers[[yr_key]]

    leafletProxy("map") %>%
      clearGroup("Sampling Depth") %>%
      addPolygons(
        data        = occ_all,
        group       = "Sampling Depth",
        layerId     = ~cell_id,
        fill        = TRUE,
        fillColor   = ~final_fill,
        fillOpacity = ~final_alpha,
        opacity     = 0,
        color       = NA,
        label       = ~paste0(
          "Depth (m, median): ", ifelse(is.na(depth_med), "NA", round(depth_med, 1)), "<br>",
          "Min: ", ifelse(is.na(depth_min), "NA", round(depth_min, 1)), " | ",
          "Max: ", ifelse(is.na(depth_max), "NA", round(depth_max, 1))
        ) %>% lapply(htmltools::HTML)
      )
  })

  # ---- clicked cell helper ----
  clicked_cell <- reactive({
    click <- input$map_shape_click
    if (is.null(click) || is.null(click$id)) return(NA_integer_)
    suppressWarnings(as.integer(click$id))
  })

  # ---- diversity detections for all MPA/AOI polygons ----
  diversity_detections_mpa <- reactive({
    mpa_membership_selection_df()
  })

  #Diversity plots

  # --- build sample x taxon matrix from current selection ---
  comm_mat_mpa <- reactive({
    det <- diversity_detections_mpa()
    req(det)

    occ_all <- det %>%
      sf::st_drop_geometry() %>%
      make_sample_id()

    shiny::validate(
      shiny::need(nrow(occ_all) > 0, "No detections available for the current filters."),
      shiny::need("organismQuantity" %in% names(occ_all),
                  "organismQuantity column is missing.")
    )

    rank_col <- active_tax_rank()

    occ_all2 <- occ_all %>%
      dplyr::mutate(
        site_name = as.character(site_name),
        taxon     = as.character(.data[[rank_col]]),
        value     = as.numeric(organismQuantity)
      ) %>%
      dplyr::filter(
        !is.na(sample_id), sample_id != "",
        !is.na(taxon), taxon != "",
        !is.na(value)
      ) %>%
      dplyr::group_by(sample_id, taxon) %>%
      dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

    shiny::validate(
      shiny::need(nrow(occ_all2) > 0, "No abundance data available after filtering.")
    )

    mat_wide <- occ_all2 %>%
      dplyr::select(sample_id, taxon, value) %>%
      tidyr::pivot_wider(
        names_from  = taxon,
        values_from = value,
        values_fill = 0
      )

    mat <- mat_wide %>%
      dplyr::select(-sample_id) %>%
      as.data.frame()

    rownames(mat) <- mat_wide$sample_id
    mat
  })


  # ---- library size diagnostics from all data ----
  library_sizes <- reactive({
    pts <- diversity_detections_beta()
    req(pts)

    occ_all <- pts %>%
      make_sample_id()

    shiny::validate(
      shiny::need(nrow(occ_all) > 0, "No detections available to calculate rarefaction depth."),
      shiny::need("organismQuantity" %in% names(occ_all), "organismQuantity column is missing.")
    )

    occ_all %>%
      dplyr::mutate(
        value = as.numeric(organismQuantity)
      ) %>%
      dplyr::filter(
        !is.na(sample_id), sample_id != "",
        !is.na(value)
      ) %>%
      dplyr::group_by(sample_id) %>%
      dplyr::summarise(
        library_size = sum(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(library_size)
  })

  confirmed_rarefaction_depth <- reactiveVal(NULL)

  observeEvent(input$div_apply, {
    depth <- input$rarefaction_depth

    shiny::validate(
      shiny::need(!is.null(depth), "Enter a rarefaction depth."),
      shiny::need(!is.na(depth), "Enter a valid rarefaction depth."),
      shiny::need(depth >= 0, "Rarefaction depth must be 0 or greater.")
    )

    confirmed_rarefaction_depth(as.numeric(depth))

    # Force alpha rarefaction summary to calculate/print first
    rr <- isolate(rarefaction_drop_summary())

    message("Rarefaction depth used for alpha diversity: ", rr$depth)
    message("Total samples in current alpha matrix: ", rr$total_samples)
    message("Kept after rarefaction filter: ", rr$kept_samples)
    message("Dropped before rarefaction: ", rr$dropped_samples)

    if (nrow(rr$dropped_table) > 0) {
      print(rr$dropped_table)
    } else {
      message("No samples dropped at this rarefaction depth.")
    }

  }, ignoreInit = TRUE, priority = 100)

  rarefaction_depth <- reactive({
    req(confirmed_rarefaction_depth())
    confirmed_rarefaction_depth()
  })

  rarefaction_drop_summary <- reactive({
    depth <- rarefaction_depth()
    req(depth)

    mat_alpha <- comm_mat_mpa()
    req(mat_alpha)

    lib_sizes <- rowSums(mat_alpha, na.rm = TRUE)

    out <- data.frame(
      sample_id    = rownames(mat_alpha),
      library_size = as.numeric(lib_sizes),
      kept         = lib_sizes >= depth,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::arrange(library_size)

    list(
      depth        = depth,
      total_samples = nrow(out),
      kept_samples  = sum(out$kept, na.rm = TRUE),
      dropped_samples = sum(!out$kept, na.rm = TRUE),
      dropped_table = out %>% dplyr::filter(!kept),
      full_table    = out
    )
  })

  # observe({
  #   rr <- rarefaction_drop_summary()
  #   req(rr)
  #
  #   message("Rarefaction depth: ", rr$depth)
  #   message("Total samples in current alpha matrix: ", rr$total_samples)
  #   message("Kept after rarefaction filter: ", rr$kept_samples)
  #   message("Dropped before rarefaction: ", rr$dropped_samples)
  # })

  # #Check which samples were dropped
  # observe({
  #   rr <- rarefaction_drop_summary()
  #   req(rr)
  #
  #   if (nrow(rr$dropped_table) > 0) {
  #     print(rr$dropped_table)
  #   } else {
  #     message("No samples dropped at this rarefaction depth.")
  #   }
  # })


  #Create rarefied matrix for alpha diversity
  comm_mat_mpa_rarefied <- reactive({
    req(input$div_apply > 0)

    mat <- comm_mat_mpa()
    depth <- rarefaction_depth()

    mat <- round(as.matrix(mat))
    storage.mode(mat) <- "integer"

    # depth = 0 means do NOT rarefy
    if (isTRUE(depth == 0)) {
      return(mat)
    }

    lib_sizes <- rowSums(mat, na.rm = TRUE)
    keep <- lib_sizes >= depth

    mat <- mat[keep, , drop = FALSE]

    shiny::validate(
      shiny::need(
        nrow(mat) > 0,
        "No filtered samples have enough reads to be rarefied at this depth."
      )
    )

    set.seed(123)
    vegan::rrarefy(mat, sample = depth)
  })

  beta_terminal_printed_key <- reactiveVal(NULL)

  sample_meta_mpa <- reactive({
    det <- diversity_detections_mpa()
    req(det)

    meta <- det %>%
      sf::st_drop_geometry() %>%
      make_sample_id() %>%
      dplyr::filter(!is.na(site_name), !is.na(site_type)) %>%
      dplyr::mutate(
        group_label = as.character(site_name)
      ) %>%
      dplyr::distinct(sample_id, group_label, .keep_all = TRUE)

    meta
  })

  # ---- base detections used for beta ordination ----
  # Includes:
  #   - all detections inside MPA/AOI polygons
  #   - plus detections inside any drawn polygons
  # Floating-panel filters and diversity target_gene filter still apply.
  diversity_detections_beta <- reactive({
    df <- dplyr::bind_rows(
      mpa_membership_selection_df(),
      polygon_membership_selection_df()
    )

    if (is.null(df) || nrow(df) == 0) {
      return(df)
    }

    apply_compare_polygon_filter(df, div_filters()$polygons)
  })

  # ---- CLR helper for beta diversity ----
  clr_transform <- function(mat, pseudocount = 1) {
    mat <- as.matrix(mat)

    if (any(mat < 0, na.rm = TRUE)) {
      stop("CLR transformation requires non-negative values.")
    }

    mat_pc <- mat + pseudocount

    log_mat <- log(mat_pc)

    gm <- rowMeans(log_mat, na.rm = TRUE)

    sweep(log_mat, 1, gm, FUN = "-")
  }

  # ---- unique sample x taxon matrix for ONE shared ordination ---- *NEW
  beta_comm_mat <- reactive({
    det <- diversity_detections_beta()
    req(det)

    occ_all <- det %>%
      sf::st_drop_geometry() %>%
      make_sample_id()

    shiny::validate(
      shiny::need(nrow(occ_all) > 0, "No detections available for the current filters."),
      shiny::need("organismQuantity" %in% names(occ_all), "organismQuantity column is missing.")
    )

    rank_col <- active_tax_rank()

    occ_all2 <- occ_all %>%
      dplyr::mutate(
        taxon = as.character(.data[[rank_col]]),
        value = as.numeric(organismQuantity)
      ) %>%
      dplyr::filter(
        !is.na(sample_id), sample_id != "",
        !is.na(taxon), taxon != "",
        !is.na(value)
      ) %>%
      dplyr::group_by(sample_id, taxon) %>%
      dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

    shiny::validate(
      shiny::need(nrow(occ_all2) > 0, "No abundance data available after filtering.")
    )

    mat_wide <- occ_all2 %>%
      tidyr::pivot_wider(
        names_from  = taxon,
        values_from = value,
        values_fill = 0
      )

    mat <- mat_wide %>%
      dplyr::select(-sample_id) %>%
      as.data.frame()

    rownames(mat) <- mat_wide$sample_id

    lib_sizes <- rowSums(mat, na.rm = TRUE)
    keep <- lib_sizes > 0
    mat <- mat[keep, , drop = FALSE]

    shiny::validate(
      shiny::need(nrow(mat) > 0,
                  "No samples with non-zero abundance remain for beta diversity.")
    )

    mat
  })

  # ---- metric-specific beta preprocessing ---- *NEW*
  beta_mat_processed <- reactive({
    mat <- beta_comm_mat()
    req(mat)

    beta_method <- input$beta_metric %||% "bray"
    mat <- as.matrix(mat)

    shiny::validate(
      shiny::need(nrow(mat) > 0, "No samples available for beta diversity."),
      shiny::need(ncol(mat) > 0, "No taxa available for beta diversity.")
    )

    if (beta_method %in% c("bray", "euclidean")) {
      # Hellinger transformation
      mat_out <- vegan::decostand(mat, method = "hellinger")

    } else if (beta_method == "jaccard") {
      # Presence-absence
      mat_out <- (mat > 0) * 1

    } else if (beta_method == "aitchison") {
      # CLR transformation
      mat_out <- clr_transform(mat, pseudocount = 1)

    } else if (beta_method == "robust.aitchison") {
      # Keep raw non-negative matrix; vegan::vegdist will handle this directly
      shiny::validate(
        shiny::need(all(mat >= 0, na.rm = TRUE),
                    "Robust Aitchison distance requires non-negative values.")
      )
      mat_out <- mat

    } else {
      mat_out <- mat
    }

    keep <- rowSums(is.na(mat_out)) == 0
    mat_out <- mat_out[keep, , drop = FALSE]

    shiny::validate(
      shiny::need(nrow(mat_out) > 1, "Not enough valid samples remain for beta diversity.")
    )

    mat_out
  })

  beta_distance <- reactive({     #NEW
    req(div_unlocked())
    mat_use <- beta_mat_processed()
    beta_method <- input$beta_metric %||% "bray"

    shiny::validate(
      shiny::need(nrow(mat_use) > 2, "At least 3 samples are required for beta diversity statistics.")
    )

    if (beta_method == "bray") {
      d <- vegan::vegdist(mat_use, method = "bray")
    } else if (beta_method == "jaccard") {
      d <- vegan::vegdist(mat_use, method = "jaccard", binary = TRUE)
    } else if (beta_method == "euclidean") {
      d <- stats::dist(mat_use, method = "euclidean")
    } else if (beta_method == "aitchison") {
      d <- stats::dist(mat_use, method = "euclidean")
    } else if (beta_method == "robust.aitchison") {
      d <- vegan::vegdist(mat_use, method = "robust.aitchison")
    } else {
      shiny::validate(
        shiny::need(FALSE, "Unsupported beta diversity metric selected.")
      )
    }

    shiny::validate(
      shiny::need(all(is.finite(as.vector(d))),
                  "Distance matrix contains non-finite values for the current filters/metric.")
    )

    d
  })

  #Polygon overlap stats warnings - display once *NEW
  alpha_overlap_warning <- reactive({
    alpha <- alpha_boxplot_occ_all()

    if (is.null(alpha) || nrow(alpha) == 0) {
      return(NULL)
    }

    dup_ids <- alpha %>%
      dplyr::filter(!is.na(sample_id), sample_id != "") %>%
      dplyr::distinct(sample_id, group_label) %>%
      dplyr::count(sample_id, name = "n_groups") %>%
      dplyr::filter(n_groups > 1)

    if (nrow(dup_ids) > 0) {
      return("Some samples belong to more than one polygon. Alpha statistics require independent group membership.")
    }

    NULL
  })

  beta_overlap_warning <- reactive({
    det <- diversity_detections_beta()

    if (is.null(det) || nrow(det) == 0) {
      return(NULL)
    }

    meta <- det %>%
      sf::st_drop_geometry() %>%
      make_sample_id() %>%
      dplyr::mutate(
        site_name   = as.character(site_name),
        group_label = dplyr::case_when(
          !is.na(site_name) & site_name != "" ~ site_name,
          TRUE ~ "NO_SITE"
        )
      ) %>%
      dplyr::filter(
        !is.na(sample_id), sample_id != "",
        !is.na(group_label), group_label != ""
      ) %>%
      dplyr::distinct(sample_id, group_label)

    dup_ids <- meta %>%
      dplyr::count(sample_id, name = "n_groups") %>%
      dplyr::filter(n_groups > 1)

    if (nrow(dup_ids) > 0) {
      return("Some samples belong to more than one polygon. Beta diversity statistics require each sample to belong to only one group.")
    }

    NULL
  })

  #Beta diagnostics - Checks for ordination (beta)  *NEW
  observeEvent(input$div_apply, {
    req(input$div_apply > 0)

    mat_raw <- beta_comm_mat()
    mat_proc <- beta_mat_processed()

    req(mat_raw, mat_proc)

    lib_sizes <- rowSums(mat_raw, na.rm = TRUE)
    zero_sum <- lib_sizes == 0

    beta_key <- paste(
      input$div_apply,
      input$beta_metric,
      input$tax_rank,
      paste(input$div_target_gene %||% character(0), collapse = "|"),
      paste(input$div_primer %||% character(0), collapse = "|"),
      paste(input$div_compare_polygons %||% character(0), collapse = "|"),
      sep = "~~"
    )

    if (identical(beta_terminal_printed_key(), beta_key)) {
      return(NULL)
    }

    beta_terminal_printed_key(beta_key)

    message("BETA raw samples: ", nrow(mat_raw))
    message("BETA zero-sum samples: ", sum(zero_sum, na.rm = TRUE))
    message("BETA min library size: ", min(lib_sizes, na.rm = TRUE))
    message("BETA median library size: ", stats::median(lib_sizes, na.rm = TRUE))
    message("BETA max library size: ", max(lib_sizes, na.rm = TRUE))
    message("BETA processed samples: ", nrow(mat_proc))
    message("BETA processed taxa: ", ncol(mat_proc))
    message("BETA metric: ", input$beta_metric %||% "bray")
  }, ignoreInit = TRUE, priority = 50)

  # ---- metadata for plotting: samples can belong to multiple groups ----
  beta_plot_meta <- reactive({
    det <- diversity_detections_beta()
    req(det)

    det %>%
      sf::st_drop_geometry() %>%
      make_sample_id() %>%
      dplyr::mutate(
        samp_name   = as.character(samp_name),
        site_name   = as.character(site_name),
        site_type   = as.character(site_type),
        year        = if ("year" %in% names(.)) as.character(year) else NA_character_,
        group_label = dplyr::case_when(
          !is.na(site_name) & site_name != "" ~ site_name,
          TRUE ~ "NO_SITE"
        )
      ) %>%
      dplyr::filter(
        !is.na(sample_id), sample_id != "",
        !is.na(group_label), group_label != ""
      ) %>%
      dplyr::distinct(sample_id, samp_name, group_label, site_type, year)
  })

  beta_stats_data <- reactive({    #NEW
    req(div_unlocked())
    if (!is.null(beta_overlap_warning())) {
      return(NULL)
    }

    d <- beta_distance()
    meta <- beta_plot_meta()

    sample_ids <- attr(d, "Labels")

    shiny::validate(
      shiny::need(!is.null(sample_ids), "Distance labels are missing.")
    )

    meta2 <- meta %>%
      dplyr::filter(sample_id %in% sample_ids) %>%
      dplyr::distinct(sample_id, .keep_all = TRUE)

    meta2 <- meta2[match(sample_ids, meta2$sample_id), , drop = FALSE]

    shiny::validate(
      shiny::need(nrow(meta2) == length(sample_ids), "Metadata do not match beta diversity samples."),
      shiny::need(dplyr::n_distinct(meta2$group_label) > 1, "At least 2 polygons are needed for beta comparison.")
    )

    list(
      dist = d,
      meta = meta2
    )
   })

  beta_stats <- reactive({
    dat <- beta_stats_data()
    if (is.null(dat)) {
      return(NULL)
    }

    d   <- dat$dist
    md  <- dat$meta

    md$group_label <- as.factor(md$group_label)

    shiny::validate(
      shiny::need(nlevels(md$group_label) > 1, "At least 2 polygons are needed for beta comparison.")
    )

    permanova <- vegan::adonis2(
      d ~ group_label,
      data = md,
      permutations = 999
    )

    anosim_res <- vegan::anosim(
      x = d,
      grouping = md$group_label,
      permutations = 999
    )

    list(
      permanova = permanova,
      anosim = anosim_res
    )
  })

  p_with_stars <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return(NA_character_)

    star_lab <- dplyr::case_when(
      p <= 0.001 ~ "***",
      p <= 0.01  ~ "**",
      p <= 0.05  ~ "*",
      p <= 0.1   ~ ".",
      TRUE       ~ ""
    )

    paste0(signif(p, 3), ifelse(star_lab == "", "", paste0(" ", star_lab)))
  }

  output$beta_stats_text <- renderText({
    if (!is.null(beta_overlap_warning())) {
      return("")
    }

    st <- beta_stats()
    if (is.null(st)) {
      return("")
    }

    per_tab <- as.data.frame(st$permanova)
    per_row <- per_tab[1, , drop = FALSE]

    paste0(
      "PERMANOVA: F = ",
      round(per_row$F[1], 3),
      ", R2 = ",
      round(per_row$R2[1], 3),
      ", p = ",
      p_with_stars(per_row$`Pr(>F)`[1]),
      " | ANOSIM: R = ",
      round(st$anosim$statistic, 3),
      ", p = ",
      p_with_stars(st$anosim$signif)
    )
  })

  beta_dispersion <- reactive({
    dat <- beta_stats_data()
    if (is.null(dat)) {
      return(NULL)
    }

    d   <- dat$dist
    md  <- dat$meta

    md$group_label <- as.factor(md$group_label)

    bd <- vegan::betadisper(d, md$group_label)
    an <- anova(bd)

    list(
      betadisper = bd,
      anova = an
    )
  })

  output$beta_dispersion_text <- renderText({
    if (!is.null(beta_overlap_warning())) {
      return("")
    }

    bd <- beta_dispersion()
    if (is.null(bd)) {
      return("")
    }

    paste0(
      "Beta dispersion: F = ",
      round(bd$anova$`F value`[1], 3),
      ", p = ",
      p_with_stars(bd$anova$`Pr(>F)`[1])
    )
  })

  output$alpha_warning_text <- renderText({  #NEW
    alpha_overlap_warning() %||% ""
  })

  output$beta_warning_text <- renderText({
    beta_overlap_warning() %||% ""
  })

  #Alpha diversity
  alpha_metric_vec <- reactive({
    req(div_unlocked())
    req(comm_mat_mpa_rarefied())
    mat <- comm_mat_mpa_rarefied()

    metric <- input$alpha_metric %||% "observed"

    vals <- switch(
      metric,
      observed   = vegan::specnumber(mat),
      shannon    = vegan::diversity(mat, index = "shannon"),
      simpson    = vegan::diversity(mat, index = "simpson"),
      invsimpson = vegan::diversity(mat, index = "invsimpson"),
      ace        = vegan::estimateR(mat)["S.ACE", ],
      pielou = {
        shannon_vals <- vegan::diversity(mat, index = "shannon")
        richness_vals <- vegan::specnumber(mat)

        evenness <- ifelse(
          richness_vals > 1,
          shannon_vals / log(richness_vals),
          NA_real_
        )

        evenness
      },
      vegan::specnumber(mat)
    )

    data.frame(
      sample_id = names(vals),
      alpha_val = as.numeric(vals),
      stringsAsFactors = FALSE
    )
  })

  #library size summary
  library_sizes_mpa <- reactive({
    library_sizes()
  })

  #Console check to verify depth
  observeEvent(input$div_apply, {

    df <- library_sizes_mpa()

    req(nrow(df) > 0)

    message("Min library size: ", min(df$library_size, na.rm = TRUE))
    message("Median library size: ", median(df$library_size, na.rm = TRUE))
    message("Max library size: ", max(df$library_size, na.rm = TRUE))

  }, ignoreInit = TRUE)

  observeEvent(input$div_apply, {

    depth <- rarefaction_depth()

    req(depth)

    message("Rarefaction depth used for alpha diversity: ", depth)

  }, ignoreInit = TRUE)

  # --- Alpha boxplot data: selected area + optional drawn polygons ---
  alpha_boxplot_occ_all <- reactive({
    req(div_unlocked())
    selected_polygons <- div_filters()$polygons %||% character(0)

    alpha <- alpha_metric_vec()
    meta  <- sample_meta_mpa()

    # ---- Built-in MPA/AOI alpha values ----
    alpha_mpa <- alpha %>%
      dplyr::left_join(
        meta %>% dplyr::distinct(sample_id, .keep_all = TRUE),
        by = "sample_id"
      ) %>%
      dplyr::filter(!is.na(site_name), site_name != "") %>%
      dplyr::mutate(group_label = as.character(group_label))

    if (length(selected_polygons) > 0) {
      alpha_mpa <- alpha_mpa %>%
        dplyr::filter(as.character(group_label) %in% selected_polygons)
    }

    # ---- Drawn polygon alpha values, using the new filtered polygon pipeline ----
    poly_df <- polygon_membership_selection_df()

    if (!is.null(poly_df) && nrow(poly_df) > 0) {

      rank_col <- active_tax_rank()
      metric   <- input$alpha_metric %||% "observed"
      depth_poly <- rarefaction_depth()

      if (is.finite(depth_poly) && depth_poly > 0 && rank_col %in% names(poly_df)) {

        poly_long <- poly_df %>%
          make_sample_id() %>%
          dplyr::mutate(
            site_name   = as.character(site_name),
            site_type   = as.character(site_type),
            group_label = as.character(site_name),
            samp_name   = as.character(samp_name),
            taxon       = as.character(.data[[rank_col]]),
            value       = as.numeric(organismQuantity)
          ) %>%
          dplyr::filter(
            site_type == "User",
            !is.na(site_name), site_name != "",
            !is.na(sample_id), sample_id != "",
            !is.na(taxon), taxon != "",
            !is.na(value)
          ) %>%
          dplyr::group_by(site_name, site_type, group_label, sample_id, samp_name, taxon) %>%
          dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

        if (length(selected_polygons) > 0) {
          poly_long <- poly_long %>%
            dplyr::filter(group_label %in% selected_polygons)
        }

        alpha_poly <- poly_long %>%
          dplyr::group_split(site_name, site_type, group_label) %>%
          purrr::map_dfr(function(dat) {

            if (nrow(dat) == 0) return(NULL)

            site_name_i   <- dat$site_name[1]
            site_type_i   <- dat$site_type[1]
            group_label_i <- dat$group_label[1]

            mat_poly <- dat %>%
              tidyr::pivot_wider(
                id_cols     = c(sample_id, samp_name),
                names_from  = taxon,
                values_from = value,
                values_fill = 0
              ) %>%
              as.data.frame()

            rownames(mat_poly) <- mat_poly$sample_id

            sample_ids_poly <- data.frame(
              sample_id = mat_poly$sample_id,
              samp_name = mat_poly$samp_name,
              stringsAsFactors = FALSE
            )

            mat_poly$sample_id <- NULL
            mat_poly$samp_name <- NULL

            lib_sizes_poly <- rowSums(mat_poly, na.rm = TRUE)
            keep_poly <- lib_sizes_poly >= depth_poly

            mat_poly <- mat_poly[keep_poly, , drop = FALSE]
            sample_ids_poly <- sample_ids_poly[keep_poly, , drop = FALSE]

            if (nrow(mat_poly) == 0 || ncol(mat_poly) == 0) return(NULL)

            mat_poly <- round(as.matrix(mat_poly))
            storage.mode(mat_poly) <- "integer"

            set.seed(123)
            mat_poly_rarefied <- vegan::rrarefy(mat_poly, sample = depth_poly)

            vals_poly <- switch(
              metric,
              observed   = vegan::specnumber(mat_poly_rarefied),
              shannon    = vegan::diversity(mat_poly_rarefied, index = "shannon"),
              simpson    = vegan::diversity(mat_poly_rarefied, index = "simpson"),
              invsimpson = vegan::diversity(mat_poly_rarefied, index = "invsimpson"),
              ace        = vegan::estimateR(mat_poly_rarefied)["S.ACE", ],
              pielou     = {
                richness <- vegan::specnumber(mat_poly_rarefied)
                shannon  <- vegan::diversity(mat_poly_rarefied, index = "shannon")
                dplyr::if_else(richness > 1, shannon / log(richness), NA_real_)
              },
              vegan::specnumber(mat_poly_rarefied)
            )

            data.frame(
              sample_id   = rownames(mat_poly_rarefied),
              alpha_val   = as.numeric(vals_poly),
              samp_name   = sample_ids_poly$samp_name[
                match(rownames(mat_poly_rarefied), sample_ids_poly$sample_id)
              ],
              site_name   = site_name_i,
              site_type   = site_type_i,
              year        = NA_character_,
              group_label = group_label_i,
              stringsAsFactors = FALSE
            )
          })

        if (!is.null(alpha_poly) && nrow(alpha_poly) > 0) {
          alpha_mpa <- dplyr::bind_rows(alpha_mpa, alpha_poly)
        }
      }
    }

    base_lvls <- unique(alpha_mpa$group_label[alpha_mpa$site_type != "User"])
    draw_lvls <- unique(alpha_mpa$group_label[alpha_mpa$site_type == "User"])

    alpha_mpa %>%
      dplyr::mutate(
        group_label = factor(group_label, levels = c(base_lvls, draw_lvls))
      )
  })

  #-----------------------------------------------------------------------------
  alpha_stats <- reactive({
    if (!is.null(alpha_overlap_warning())) {
      return(NULL)
    }

    alpha <- alpha_boxplot_occ_all()

    shiny::validate(
      shiny::need(nrow(alpha) > 0, "No alpha diversity data available.")
    )

    alpha <- alpha %>%
      dplyr::filter(!is.na(alpha_val), !is.na(group_label), !is.na(sample_id)) %>%
      dplyr::mutate(group_label = droplevels(as.factor(group_label)))

    n_groups <- dplyr::n_distinct(alpha$group_label)

    shiny::validate(
      shiny::need(n_groups > 1, "At least 2 polygons are needed for alpha comparison.")
    )

    if (n_groups == 2) {
      overall <- wilcox.test(alpha_val ~ group_label, data = alpha, exact = FALSE)

      pairwise_df <- data.frame(
        group1 = levels(alpha$group_label)[1],
        group2 = levels(alpha$group_label)[2],
        p_adj  = overall$p.value,
        stringsAsFactors = FALSE
      )

      return(list(
        method   = "wilcox",
        overall  = overall,
        pairwise = pairwise_df
      ))
    }

    overall <- kruskal.test(alpha_val ~ group_label, data = alpha)

    pairwise <- pairwise.wilcox.test(
      x = alpha$alpha_val,
      g = alpha$group_label,
      p.adjust.method = "BH",
      exact = FALSE
    )

    pairwise_df <- as.data.frame(as.table(pairwise$p.value), stringsAsFactors = FALSE)
    names(pairwise_df) <- c("group1", "group2", "p_adj")

    pairwise_df <- pairwise_df %>%
      dplyr::filter(!is.na(p_adj)) %>%
      dplyr::arrange(p_adj)

    list(
      method   = "kruskal",
      overall  = overall,
      pairwise = pairwise_df
    )
  })

  output$alpha_summary_tbl <- DT::renderDT({
    alpha <- alpha_boxplot_occ_all()

    shiny::validate(
      shiny::need(nrow(alpha) > 0, "No alpha diversity data available.")
    )

    out <- alpha %>%
      dplyr::filter(!is.na(group_label), !is.na(alpha_val)) %>%
      dplyr::group_by(group_label) %>%
      dplyr::summarise(
        n_samples = dplyr::n(),
        mean = round(mean(alpha_val, na.rm = TRUE), 2),
        median = round(stats::median(alpha_val, na.rm = TRUE), 2),
        sd = round(stats::sd(alpha_val, na.rm = TRUE), 2),
        min = round(min(alpha_val, na.rm = TRUE), 2),
        q1 = round(as.numeric(stats::quantile(alpha_val, 0.25, na.rm = TRUE, type = 7)), 2),
        q3 = round(as.numeric(stats::quantile(alpha_val, 0.75, na.rm = TRUE, type = 7)), 2),
        max = round(max(alpha_val, na.rm = TRUE), 2),
        .groups = "drop"
      ) %>%
      dplyr::rename(Polygon = group_label)

    DT::datatable(
      out,
      rownames = FALSE,
      colnames = c(
        "Polygon",
        "Number of Samples",
        "Mean",
        "Median",
        "Standard Deviation",
        "Minimum",
        "Q1",
        "Q3",
        "Maximum"
      ),
      options = list(
        pageLength = 10,
        dom = "tip",
        scrollX = TRUE,
        autoWidth = FALSE,
        scrollCollapse = TRUE
      ),
      class = "nowrap"
    )
  })


  output$alpha_stats_text <- renderText({
    if (!is.null(alpha_overlap_warning())) {
      return("")
    }

    st <- alpha_stats()
    if (is.null(st)) {
      return("")
    }

    if (st$method == "wilcox") {
      return(
        paste0(
          "Wilcoxon rank-sum test: W = ",
          round(unname(st$overall$statistic), 3),
          ", p = ",
          p_with_stars(st$overall$p.value)
        )
      )
    }

    paste0(
      "Kruskal-Wallis: chi-squared = ",
      round(unname(st$overall$statistic), 3),
      ", df = ",
      unname(st$overall$parameter),
      ", p = ",
      p_with_stars(st$overall$p.value)
    )
  })

  output$alpha_pairwise_tbl <- DT::renderDT({
    st <- alpha_stats()

    shiny::validate(
      shiny::need(!is.null(st), "")
    )

    df <- st$pairwise

    shiny::validate(
      shiny::need(nrow(df) > 0, "No pairwise comparisons available.")
    )

    df <- df %>%
      dplyr::mutate(
        stars = dplyr::case_when(
          p_adj <= 0.001 ~ "***",
          p_adj <= 0.01  ~ "**",
          p_adj <= 0.05  ~ "*",
          p_adj <= 0.1   ~ ".",
          TRUE           ~ ""
        )
      )

    p_label <- if (st$method == "wilcox") {
      "p-value"
    } else {
      "Adjusted p-value"
    }

    df <- df %>%
      dplyr::mutate(
        p_display = ifelse(
          stars == "",
          as.character(signif(p_adj, 3)),
          paste0(signif(p_adj, 3), " ", stars)
        )
      ) %>%
      dplyr::select(group1, group2, p_display)

    DT::datatable(
      df,
      rownames = FALSE,
      colnames = c("Polygon 1", "Polygon 2", p_label),
      options = list(
        pageLength = 10,
        dom = "tip",
        scrollX = TRUE
      )
    )
  })

  make_cld_letters <- function(alpha_df, st, alpha = 0.05) {
    groups <- levels(droplevels(as.factor(alpha_df$group_label)))

    # If only 1 group, return that group with "A"
    if (length(groups) == 1) {
      return(data.frame(
        group_label = groups,
        letters = "A",
        stringsAsFactors = FALSE
      ))
    }

    # Start with matrix of "not significant"
    pmat <- matrix(
      1,
      nrow = length(groups),
      ncol = length(groups),
      dimnames = list(groups, groups)
    )

    diag(pmat) <- 1

    # Fill matrix from pairwise results
    if (!is.null(st$pairwise) && nrow(st$pairwise) > 0) {
      for (i in seq_len(nrow(st$pairwise))) {
        g1 <- as.character(st$pairwise$group1[i])
        g2 <- as.character(st$pairwise$group2[i])
        p  <- as.numeric(st$pairwise$p_adj[i])

        if (!is.na(g1) && !is.na(g2) && g1 %in% groups && g2 %in% groups) {
          pmat[g1, g2] <- p
          pmat[g2, g1] <- p
        }
      }
    }

    # multcompLetters expects logical matrix:
    # TRUE means groups are NOT significantly different
    letter_obj <- multcompView::multcompLetters(pmat < alpha)

    data.frame(
      group_label = names(letter_obj$Letters),
      letters = unname(letter_obj$Letters),
      stringsAsFactors = FALSE
    )
  }

  alpha_plot_annotations <- reactive({
    st <- alpha_stats()
    alpha <- alpha_boxplot_occ_all()

    req(st, alpha)
    req(nrow(alpha) > 0)

    alpha2 <- alpha %>%
      dplyr::filter(!is.na(group_label), !is.na(alpha_val)) %>%
      dplyr::mutate(group_label = droplevels(as.factor(group_label)))

    cld <- make_cld_letters(alpha2, st, alpha = 0.05)

    y_pos <- alpha2 %>%
      dplyr::group_by(group_label) %>%
      dplyr::summarise(
        y = max(alpha_val, na.rm = TRUE) * 1.08,
        .groups = "drop"
      ) %>%
      dplyr::mutate(group_label = as.character(group_label))

    cld %>%
      dplyr::mutate(group_label = as.character(group_label)) %>%
      dplyr::left_join(y_pos, by = "group_label") %>%
      dplyr::mutate(
        x = group_label
      )
  })

  #----------------------------------------------------------------
  get_group_palette <- function(group_levels) {
    group_levels <- unique(as.character(group_levels))

    base_cols <- c(
      "Eastern Shore Islands AOI" = "#046c9a",
      "St. Anns Bank Marine Protected Area" = "#5BBCD6",
      "Musquash Estuary Marine Protected Area" = "#ABDDDE",
      "Fundian Channel–Browns Bank AOI" = "#446455",
      "Gully Marine Protected Area" = "#FDD262"
    )

    pal <- setNames(rep(NA_character_, length(group_levels)), group_levels)

    known <- intersect(group_levels, names(base_cols))
    pal[known] <- base_cols[known]

    unknown <- setdiff(group_levels, names(base_cols))

    if (length(unknown) > 0) {
      fallback_cols <- grDevices::hcl.colors(
        n = length(unknown),
        palette = "Temps"
      )
      pal[unknown] <- fallback_cols
    }

    pal
  }

  # --- Alpha diversity boxplot (Plotly) ---
  #color_vec <- c("#046c9a", "#5BBCD6", "#ABDDDE", "#446455", "#00A08A","#Fdd262")

  output$alpha_boxplot <- plotly::renderPlotly({

    alpha <- alpha_boxplot_occ_all()
    ann   <- alpha_plot_annotations()

    shiny::validate(
      shiny::need(nrow(alpha) > 0, "No samples available for the current selection/year.")
    )

    y_max <- max(alpha$alpha_val, na.rm = TRUE)
    ann_top <- if (!is.null(ann) && nrow(ann) > 0) {
      max(ann$y, na.rm = TRUE)
    } else {
      y_max * 1.10
    }
    y_upper <- max(y_max * 1.15, ann_top * 1.10)

    metric_label <- switch(
      input$alpha_metric %||% "observed",
      observed   = "Richness (count of taxa)",
      shannon    = "Shannon diversity",
      simpson    = "Simpson diversity",
      invsimpson = "Inverse Simpson diversity",
      ace        = "ACE estimated richness",
      pielou     = "Pielou's evenness",
      "Alpha diversity"
    )

    rank_label   <- tools::toTitleCase(gsub("_", " ", active_tax_rank()))
    gene_label   <- div_filters()$target_gene
    primer_label <- div_filters()$primers

    gene_suffix <- if (length(gene_label) > 0) {
      paste0(" | Genes: ", paste(gene_label, collapse = ", "))
    } else {
      ""
    }

    primer_suffix <- if (length(primer_label) > 0) {
      paste0(" | Primers: ", paste(primer_label, collapse = "; "))
    } else {
      ""
    }

    metric_label <- paste0(
      metric_label, " at ", rank_label, " level",
      gene_suffix, primer_suffix
    )

    group_levels <- unique(as.character(alpha$group_label))
    alpha$group_label <- factor(alpha$group_label, levels = group_levels)

    pal <- get_group_palette(group_levels)

    if (!is.null(ann) && nrow(ann) > 0) {
      ann <- ann %>%
        dplyr::mutate(
          x = as.character(x),
          group_label = factor(x, levels = group_levels)
        )
    }

    p <- plotly::plot_ly()

    for (grp in group_levels) {
      df_grp <- alpha[as.character(alpha$group_label) == grp, , drop = FALSE]
      grp_col <- unname(pal[grp])
      if (is.na(grp_col) || is.null(grp_col)) grp_col <- "#333333"

      p <- p %>%
        plotly::add_trace(
          data = df_grp,
          x = ~group_label,
          y = ~alpha_val,
          type = "box",
          name = grp,
          legendgroup = grp,
          showlegend = TRUE,
          marker = list(
            color = grp_col,
            opacity = 0.8
          ),
          line = list(color = grp_col),
          fillcolor = grDevices::adjustcolor(grp_col, alpha.f = 0.35),
          boxpoints = "all",
          jitter = 0.3,
          pointpos = 0,
          customdata = ~samp_name,
          hovertemplate = paste(
            "<b>Sample:</b> %{customdata}<br>",
            "<b>Group:</b> %{x}<br>",
            "<b>", metric_label, ":</b> %{y}<extra></extra>"
          ),
          inherit = FALSE
        )

      if (!is.null(ann) && nrow(ann) > 0) {
        ann_grp <- ann[ann$x == grp, , drop = FALSE]

        if (nrow(ann_grp) > 0) {
          p <- p %>%
            plotly::add_text(
              data = ann_grp,
              x = ~group_label,
              y = ~y,
              text = ~letters,
              textposition = "top center",
              textfont = list(size = 18, color = "black"),
              showlegend = FALSE,
              legendgroup = grp,
              hoverinfo = "skip",
              inherit = FALSE
            )
        }
      }
    }

    p %>%
      plotly::layout(
        font = list(size = 18),
        xaxis = list(
          title = list(text = "Location (Polygon)", font = list(size = 20)),
          tickfont = list(size = 16),
          showgrid = FALSE,
          showline = TRUE,
          linecolor = "black",
          rangemode = "tozero"
        ),
        yaxis = list(
          title = list(text = metric_label, font = list(size = 20)),
          tickfont = list(size = 16),
          showgrid = FALSE,
          showline = TRUE,
          linecolor = "black",
          range = c(0, y_upper),
          fixedrange = FALSE
        ),
        legend = list(
          title = list(text = "Polygon"),
          font = list(size = 14),
          itemclick = "toggle",
          itemdoubleclick = "toggleothers",
          groupclick = "togglegroup"
        ),
        margin = list(l = 80, r = 30, t = 40, b = 120),
        showlegend = TRUE
      ) %>%
      htmlwidgets::onRender("
        function(el,x){
          setTimeout(function(){
            $('#beta_loading_overlay').addClass('hidden');
          },300);
        }
      ")
  })

  make_ellipse <- function(df, conf = 0.95, npoints = 100) {
    if (nrow(df) < 3) return(NULL)

    center <- c(mean(df$PC1, na.rm = TRUE), mean(df$PC2, na.rm = TRUE))
    cov_mat <- stats::cov(df[, c("PC1", "PC2")], use = "complete.obs")

    if (any(!is.finite(cov_mat)) || det(cov_mat) <= 0) return(NULL)

    angles <- seq(0, 2 * pi, length.out = npoints)
    circle <- cbind(cos(angles), sin(angles))

    radius <- sqrt(stats::qchisq(conf, df = 2))
    ellipse <- t(center + radius * t(circle %*% chol(cov_mat)))

    out <- data.frame(
      PC1 = ellipse[, 1],
      PC2 = ellipse[, 2]
    )
    out
  }

  #Beta diversity *NEW
  output$beta_pcoa <- plotly::renderPlotly({
    if (!is.null(beta_overlap_warning())) {
      return(NULL)
    }

    dat  <- beta_stats_data()
    d    <- dat$dist
    meta <- dat$meta

    sample_ids <- attr(d, "Labels")

    shiny::validate(
      shiny::need(!is.null(sample_ids), "Sample names are missing from the beta diversity distance object."),
      shiny::need(length(sample_ids) > 2, "At least 3 samples are required to compute a PCoA.")
    )

    beta_method <- input$beta_metric %||% "bray"

    shiny::validate(
      shiny::need(
        all(is.finite(as.vector(d))),
        "Distance matrix contains non-finite values for the current filters/metric."
      )
    )

    ord <- tryCatch(
      stats::cmdscale(d, k = 2, eig = TRUE),
      error = function(e) NULL
    )

    shiny::validate(
      shiny::need(!is.null(ord), "PCoA could not be computed for the current beta diversity selection."),
      shiny::need(!is.null(ord$points), "PCoA returned no coordinates.")
    )

    shiny::validate(
      shiny::need(length(sample_ids) == nrow(ord$points), "Mismatch between sample names and ordination points.")
    )

    scores <- data.frame(
      sample_id = sample_ids,
      PC1 = ord$points[, 1],
      PC2 = ord$points[, 2],
      stringsAsFactors = FALSE
    )

    plot_df <- scores %>%
      dplyr::left_join(meta, by = "sample_id") %>%
      dplyr::mutate(
        group_label = as.character(group_label),
        site_type   = as.character(site_type)
      ) %>%
      dplyr::filter(!is.na(group_label), nzchar(group_label)) %>%
      dplyr::distinct(sample_id, group_label, site_type, .keep_all = TRUE)

    shiny::validate(
      shiny::need(nrow(plot_df) > 0, "No plotting metadata available for the current filters.")
    )

    eig <- ord$eig
    eig_pos <- eig[eig > 0]
    pct1 <- if (length(eig_pos) >= 1) round(100 * eig_pos[1] / sum(eig_pos), 1) else NA_real_
    pct2 <- if (length(eig_pos) >= 2) round(100 * eig_pos[2] / sum(eig_pos), 1) else NA_real_

    method_label <- switch(
      beta_method,
      bray      = "Bray-Curtis (Hellinger transformed)",
      jaccard   = "Jaccard (presence-absence)",
      euclidean = "Euclidean (Hellinger transformed)",
      aitchison = "Aitchison (CLR transformed)",
      `robust.aitchison` = "Robust Aitchison",
      beta_method
    )

    rank_label   <- tools::toTitleCase(gsub("_", " ", active_tax_rank()))
    gene_label   <- div_filters()$target_gene
    primer_label <- div_filters()$primers

    gene_suffix <- if (length(gene_label) > 0) {
      paste0(" | Genes: ", paste(gene_label, collapse = ", "))
    } else {
      ""
    }

    primer_suffix <- if (length(primer_label) > 0) {
      paste0(" | Primers: ", paste(primer_label, collapse = "; "))
    } else {
      ""
    }

    plot_title <- paste0(
      "Beta diversity (PCoA; ", method_label, ") at ",
      rank_label, " level",
      gene_suffix, primer_suffix
    )

    plot_df <- plot_df %>%
      dplyr::mutate(
        hover_text = paste0(
          "<b>Sample:</b> ", samp_name,
          "<br><b>Group:</b> ", group_label,
          "<br><b>Type:</b> ", site_type,
          "<br><b>PC1:</b> ", sprintf("%.3f", PC1),
          "<br><b>PC2:</b> ", sprintf("%.3f", PC2)
        )
      )

    make_ellipse <- function(df, conf = 0.95, npoints = 100) {
      if (nrow(df) < 3) return(NULL)

      xy <- df[, c("PC1", "PC2"), drop = FALSE]
      xy <- xy[stats::complete.cases(xy), , drop = FALSE]

      if (nrow(xy) < 3) return(NULL)

      cov_mat <- tryCatch(stats::cov(xy), error = function(e) NULL)
      if (is.null(cov_mat)) return(NULL)
      if (any(!is.finite(cov_mat))) return(NULL)
      if (det(cov_mat) <= 0) return(NULL)

      center <- c(mean(xy$PC1), mean(xy$PC2))
      angles <- seq(0, 2 * pi, length.out = npoints)
      unit_circle <- cbind(cos(angles), sin(angles))

      radius <- sqrt(stats::qchisq(conf, df = 2))

      chol_decomp <- tryCatch(chol(cov_mat), error = function(e) NULL)
      if (is.null(chol_decomp)) return(NULL)

      ellipse_coords <- radius * unit_circle %*% chol_decomp
      ellipse_coords <- sweep(ellipse_coords, 2, center, FUN = "+")

      data.frame(
        PC1 = ellipse_coords[, 1],
        PC2 = ellipse_coords[, 2],
        stringsAsFactors = FALSE
      )
    }

    group_levels <- unique(as.character(plot_df$group_label))
    pal <- get_group_palette(group_levels)

    ellipse_list <- lapply(group_levels, function(grp) {
      df_grp <- plot_df %>% dplyr::filter(group_label == grp)
      ell <- make_ellipse(df_grp, conf = 0.95, npoints = 100)
      if (is.null(ell)) return(NULL)
      ell$group_label <- grp
      ell
    })
    names(ellipse_list) <- group_levels

    p <- plotly::plot_ly()

    for (grp in group_levels) {
      df_grp <- plot_df %>% dplyr::filter(group_label == grp)
      ell_df <- ellipse_list[[grp]]
      grp_col <- unname(pal[grp])
      if (is.na(grp_col) || is.null(grp_col)) grp_col <- "#333333"

      if (!is.null(ell_df) && nrow(ell_df) > 0) {
        p <- p %>%
          plotly::add_trace(
            data = ell_df,
            x = ~PC1,
            y = ~PC2,
            type = "scatter",
            mode = "lines",
            fill = "toself",
            fillcolor = grDevices::adjustcolor(grp_col, alpha.f = 0.20),
            line = list(width = 1, color = grp_col),
            opacity = 1,
            hoverinfo = "skip",
            legendgroup = grp,
            name = grp,
            showlegend = FALSE,
            inherit = FALSE
          )
      }

      p <- p %>%
        plotly::add_trace(
          data = df_grp,
          x = ~PC1,
          y = ~PC2,
          type = "scatter",
          mode = "markers",
          marker = list(
            size = 9,
            color = grp_col
          ),
          text = ~hover_text,
          hoverinfo = "text",
          legendgroup = grp,
          name = grp,
          showlegend = TRUE,
          inherit = FALSE
        )
    }

    p %>%
      plotly::layout(
        title = list(
          text = plot_title,
          font = list(size = 20)
        ),
        font = list(size = 18),
        xaxis = list(
          title = list(
            text = paste0("PC1", if (!is.na(pct1)) paste0(" (", pct1, "%)") else ""),
            font = list(size = 20)
          ),
          tickfont = list(size = 16),
          showgrid = FALSE,
          showline = TRUE,
          linecolor = "black"
        ),
        yaxis = list(
          title = list(
            text = paste0("PC2", if (!is.na(pct2)) paste0(" (", pct2, "%)") else ""),
            font = list(size = 20)
          ),
          tickfont = list(size = 16),
          showgrid = FALSE,
          showline = TRUE,
          linecolor = "black"
        ),
        legend = list(
          title = list(text = "Polygon"),
          font = list(size = 14),
          groupclick = "togglegroup",
          itemdoubleclick = "toggleothers"
        ),
        margin = list(l = 80, r = 30, t = 70, b = 80),
        showlegend = TRUE
      ) %>%
      htmlwidgets::onRender("
    function(el,x){
      setTimeout(function(){
        $('#alpha_loading_overlay').addClass('hidden');
      },300);
    }
  ")
  })



  #Krona plot

  output$tax_krona <- taxplore::renderKronaChart({
    req(div_unlocked())
    ids <- selection_ids()

    shiny::validate(
      shiny::need(
        !is.null(ids) && length(ids) > 0,
        "Select a cell/polygon (or draw a polygon) to display a Krona chart."
      )
    )

    det0 <- selection_selection_df()

    shiny::validate(
      shiny::need(nrow(det0) > 0, "No data available for this selection."),
      shiny::need("organismQuantity" %in% names(det0),
                  "organismQuantity column is missing.")
    )

    tax_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")

    krona_occ_all <- det0 %>%
      dplyr::mutate(
        species = as.character(scientificName),
        qty     = as.numeric(organismQuantity)
      ) %>%
      dplyr::filter(
        !is.na(species), species != "",
        !is.na(qty), qty > 0
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(tax_ranks, "species")))) %>%
      dplyr::summarise(
        magnitude = sum(qty, na.rm = TRUE),
        .groups   = "drop"
      )

    shiny::validate(
      shiny::need(
        nrow(krona_occ_all) > 0,
        "No detections with organismQuantity > 0 after the current target gene and primer filters."
      )
    )

    tax_occ_all <- krona_occ_all %>%
      dplyr::select(dplyr::all_of(c(tax_ranks, "species")))

    taxplore::plot_krona(
      tax_occ_all,
      magnitude   = krona_occ_all$magnitude,
      total_label = "Sum organismQuantity"
    ) %>%
      htmlwidgets::onRender("
        function(el,x){
          setTimeout(function(){
            $('#tax_loading_overlay').addClass('hidden');
          },300);
        }
      ")
  })

  #Helper to generate species vectors
  get_species_vec <- function(det_sf, gene = NULL) {
    req(det_sf)

    occ_all <- det_sf

    if (!is.null(gene)) {
      occ_all <- occ_all %>% dplyr::filter(as.character(target_gene) == as.character(gene))
    }

    occ_all <- apply_species_filters(occ_all)

    occ_all %>%
      sf::st_drop_geometry() %>%
      dplyr::pull(scientificName) %>%
      as.character() %>%
      unique() %>%
      stats::na.omit() %>%
      sort()
  }

  # ---- species panel (drawn polygon OR click) ----
  output$species_panel <- renderUI({

    active_groups()
    filter_sara_on()
    filter_ais_on()

    yr <- sel_year_chr()

    make_list_local <- function(vec, max_h = 220) {
      if (length(vec) == 0) return(em("No species detected for this layer."))
      tags$div(
        style = paste0("max-height: ", max_h, "px; overflow-y: auto; padding-left: 10px;"),
        tags$ul(lapply(vec, tags$li))
      )
    }

    render_by_layer_ui <- function(by_layer, header = NULL, max_h = 220) {
      if (is.null(by_layer) || length(by_layer) == 0) {
        return(em("Turn ON a richness layer (All / 12S / COI / 16S / 18S) to view species lists."))
      }

      ord <- c("All", "12S", "COI", "16S", "18S")
      nm  <- intersect(ord, names(by_layer))
      if (length(nm) == 0) nm <- names(by_layer)

      panels <- tagList()
      for (k in nm) {
        vec <- by_layer[[k]]
        panels <- tagAppendChildren(
          panels,
          tags$details(
            open = TRUE,
            tags$summary(strong(paste0(k, " (", length(vec), ")"))),
            make_list_local(vec, max_h = max_h),
            tags$hr()
          )
        )
      }

      if (!is.null(header)) tagList(header, tags$hr(), panels) else panels
    }

    layers_on <- active_richness_layers()
    if (length(layers_on) == 0) {
      return(em("Turn ON a richness layer (All / 12S / COI / 16S / 18S) to view species lists."))
    }

    sel_id <- selected_draw_id()
    click  <- input$map_shape_click

    # drawn polygon branch
    if (!is.null(sel_id)) {
      polys <- drawn_polys()
      hit <- polys %>% dplyr::filter(draw_id == sel_id)
      if (nrow(hit) == 0) return(em("No species detected for the current selection."))

      pts <- species_sf_by_year[[yr]]
      if (is.null(pts)) pts <- species_sf_min[0, ]

      inside <- get_inside_cached(
        selected_key = paste0("draw:", sel_id),
        pts          = pts,
        geom         = sf::st_geometry(hit),
        year         = yr
      )

      by_layer_all <- get_species_by_layer_cached(
        selected_key = paste0("draw:", sel_id),
        det_sf       = inside,
        year         = yr,
        groups       = active_groups(),
        iucn_on      = isTRUE(filter_iucn_on()),
        sara_on      = isTRUE(filter_sara_on()),
        ais_on       = isTRUE(filter_ais_on())
      )

      by_layer <- by_layer_all[intersect(names(by_layer_all), layers_on)]
      poly_lab <- hit$draw_label[1] %||% sel_id

      return(render_by_layer_ui(by_layer, strong(paste0("Drawn polygon: ", poly_lab)), 180))
    }

    if (is.null(click) || is.null(click$id)) {
      return(em("Click a grid cell or an MPA/AOI outline, or draw a polygon."))
    }

    if (grepl("\\|\\|", click$id)) {
      groups_on <- input$map_groups %||% character(0)
      show_poly_total <- "MPA/AOI total species richness" %in% groups_on
      if (!show_poly_total) {
        return(em("Turn ON “MPA/AOI total species richness” to view polygon species lists."))
      }

      parts  <- strsplit(click$id, "\\|\\|")[[1]]
      p_type <- parts[1]
      p_name <- parts[2]

      poly_sel <- all_polys_click %>%
        dplyr::filter(site_type == p_type, site_name == p_name)

      if (nrow(poly_sel) == 0) return(em("Polygon not found."))

      pts <- species_sf_by_year[[yr]]
      if (is.null(pts)) pts <- species_sf_min[0, ]

      inside <- get_inside_cached(
        selected_key = paste0("click:", click$id),
        pts          = pts,
        geom         = sf::st_geometry(poly_sel),
        year         = yr
      )

      by_layer_all <- get_species_by_layer_cached(
        selected_key = paste0("click:", click$id),
        det_sf       = inside,
        year         = yr,
        groups       = active_groups(),
        iucn_on      = isTRUE(filter_iucn_on()),
        sara_on      = isTRUE(filter_sara_on()),
        ais_on       = isTRUE(filter_ais_on())
      )

      by_layer <- by_layer_all[intersect(names(by_layer_all), layers_on)]
      return(render_by_layer_ui(by_layer, strong(paste0(p_type, ": ", p_name)), 220))
    }

    cid <- suppressWarnings(as.integer(click$id))
    if (is.na(cid)) {
      return(em("Click a grid cell or an MPA/AOI outline."))
    }

    pts <- species_sf_by_year[[yr]]
    if (is.null(pts)) pts <- species_sf_min[0, ]

    cell_poly <- grid_clip %>% dplyr::filter(cell_id == cid)

    det <- get_inside_cached(
      selected_key = paste0("click:", click$id),
      pts          = pts,
      geom         = sf::st_geometry(cell_poly),
      year         = yr
    )

    if (is.null(det) || nrow(det) == 0) {
      return(em("No species detected for the current selection."))
    }

    by_layer_all <- get_species_by_layer_cached(
      selected_key = paste0("click:", click$id),
      det_sf       = det,
      year         = yr,
      groups       = active_groups(),
      iucn_on      = isTRUE(filter_iucn_on()),
      sara_on      = isTRUE(filter_sara_on()),
      ais_on       = isTRUE(filter_ais_on())
    )

    by_layer <- by_layer_all[intersect(names(by_layer_all), layers_on)]
    render_by_layer_ui(by_layer, strong(paste0("Grid cell: ", cid)), 220)
  })

  selection_details_df <- reactive({
    det <- selected_detections()

    if (is.null(det) || nrow(det) == 0) {
      return(det)
    }

    coords <- sf::st_coordinates(det)

    det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        decimalLongitude = coords[, 1],
        decimalLatitude  = coords[, 2]
      )
  })

  selection_details_panel_df <- reactive({
    df <- selection_details_df()

    if (is.null(df) || nrow(df) == 0) {
      return(df)
    }

    apply_species_filters(df)
  })

  output$detections_tbl <- DT::renderDT({
    det_raw <- selection_details_df()

    if (is.null(det_raw) || nrow(det_raw) == 0) {
      return(DT::datatable(
        data.frame(Message = "Click an MPA/AOI polygon/grid cell, or draw a polygon, to view detection details"),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det_show <- selection_details_panel_df()

    if (is.null(det_show) || nrow(det_show) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections match the current Year and/or Group filters."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }


    needed_cols <- c(
      "kingdom", "phylum", "class", "order", "family", "genus",
      "scientificName", "samp_name", "target_gene", "eventDate",
      "pcr_primer_name_forward", "pcr_primer_name_reverse",
      "pcr_primer_forward", "pcr_primer_reverse",
      "organismQuantity", "project_contact", "LClabel",
      "samp_size", "samp_size_unit", "occurrenceID",
      "minimumDepthInMeters", "maximumDepthInMeters",
      "bathymetry", "associatedSequences",
      "decimalLatitude", "decimalLongitude",
      "flags", "dataset_id"
    )

    missing_cols <- setdiff(needed_cols, names(det_show))
    for (nm in missing_cols) {
      det_show[[nm]] <- NA_character_
    }

    fwd_col <- get_forward_primer_col(det_show)
    rev_col <- get_reverse_primer_col(det_show)

    fwd_vals <- if (!is.null(fwd_col)) as.character(det_show[[fwd_col]]) else rep(NA_character_, nrow(det_show))
    rev_vals <- if (!is.null(rev_col)) as.character(det_show[[rev_col]]) else rep(NA_character_, nrow(det_show))

    fwd_vals <- trimws(fwd_vals)
    rev_vals <- trimws(rev_vals)
    fwd_vals[fwd_vals == ""] <- NA_character_
    rev_vals[rev_vals == ""] <- NA_character_

    det_show$Primers <- dplyr::case_when(
      !is.na(fwd_vals) & !is.na(rev_vals) ~ paste(fwd_vals, rev_vals, sep = " | "),
      !is.na(fwd_vals) ~ fwd_vals,
      !is.na(rev_vals) ~ rev_vals,
      TRUE ~ NA_character_
    )

    det_show <- det_show %>%
      dplyr::mutate(
        across(c(
          scientificName, samp_name, target_gene,
          bathymetry, flags, dataset_id, eventDate
        ), ~ as.character(.)),
        Volume = dplyr::case_when(
          !is.na(samp_size) & !is.na(samp_size_unit) ~ paste(samp_size, samp_size_unit),
          !is.na(samp_size) ~ as.character(samp_size),
          !is.na(samp_size_unit) ~ as.character(samp_size_unit),
          TRUE ~ NA_character_
        ),
        `Depth Range (m)` = dplyr::case_when(
          !is.na(minimumDepthInMeters) & !is.na(maximumDepthInMeters) ~
            paste(minimumDepthInMeters, maximumDepthInMeters, sep = " | "),
          !is.na(minimumDepthInMeters) ~ as.character(minimumDepthInMeters),
          !is.na(maximumDepthInMeters) ~ as.character(maximumDepthInMeters),
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::arrange(scientificName, samp_name) %>%
      dplyr::select(
        kingdom,
        phylum,
        class,
        order,
        family,
        genus,
        scientificName,
        samp_name,
        eventDate,
        target_gene,
        Primers,
        organismQuantity,
        Volume,
        `Depth Range (m)`,
        bathymetry,
        decimalLatitude,
        decimalLongitude,
        flags,
        occurrenceID,
        dataset_id,
        associatedSequences,
        project_contact,
        LClabel
      )

    DT::datatable(
      det_show,
      rownames = FALSE,
      colnames = c(
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Sample Name",
        "Date",
        "Target Gene",
        "Primers",
        "# of Sequence Reads",
        "Volume",
        "Minimum/Maximum Depth (m)",
        "Bathymetry",
        "Latitude",
        "Longitude",
        "OBIS Data Flags",
        "OBIS Occurrence ID",
        "OBIS Dataset Identifier",
        "Associated Sequences",
        "Project Contact",
        "Indigenous Acknowledgement & Contributions"
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = FALSE,
        scrollCollapse = TRUE
      ),
      class = "nowrap"
    )
  })

  # highlight selected cell (outline on top)
  observeEvent(input$map_shape_click, {
    click <- input$map_shape_click
    req(click, click$id)

    id <- as.character(click$id)

    polys <- drawn_polys()

    # clicked one of the stored drawn polygons
    if (nrow(polys) > 0 && id %in% polys$draw_id) {
      selected_draw_id(id)
      session$sendCustomMessage("openFloating", list(id = id))
      return()
    }

    # otherwise switch away from drawn polygons
    selected_draw_id(NULL)

    is_mpa  <- grepl("\\|\\|", id)
    is_cell <- !is.na(suppressWarnings(as.integer(id)))

    if (is_mpa || is_cell) {
      session$sendCustomMessage("openFloating", list(id = id))
    }
  }, ignoreInit = TRUE)

}
