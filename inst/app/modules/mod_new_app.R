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
            #tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_method", "Method Comparison")),
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
        id = "sec_map", class = "scroll-section",
        div(
          id = "map_wrap",
          leafletOutput("map"),

          div(
            id = "monthly_plot_control",
            class = "leaflet-control",
            div(id = "monthly_plot_title", "Monthly Number of Samples Collected"),
            div(id = "monthly_plot_subtitle", textOutput("monthly_plot_subtitle", inline = TRUE)),
            plotOutput("monthly_circular_plot", height = "220px", width = "100%")
          ),

          absolutePanel(
            id = "floating_panel",
            fixed = FALSE, draggable = FALSE,
            top = 10, left = 65, width = 360,

            tags$button(
              id = "floating_toggle",
              type = "button",
              class = "btn btn-default btn-secondary",
              `data-bs-toggle`  = "collapse",
              `data-bs-target`  = "#floating_body",
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
                      alt = "Plants",
                      style = "height:40px;"  # adjust as needed
                    ),
                    title = "Plants",
                    class = "btn btn-default btn-secondary filter-btn"
                  ),
                ),

                tags$div(style="height:8px;"),

                # --- Row 3: 2 across ---
                div(
                  class = "filter-btn-grid-2",
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
          tabPanel("Detection Details",              DT::DTOutput("detections_tbl")),
          tabPanel("Species at Risk Act (SARA): Schedule 1 Details", DT::DTOutput("sara_details")),
          tabPanel("Aquatic Invasive Species (AIS) Details",             DT::DTOutput("ais_details"))
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
            selectizeInput(
              "div_target_gene",
              "Target gene",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list(
                plugins = list("remove_button"),
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
                "Observed (Richness)" = "observed",
                "Shannon"             = "shannon",
                "Simpson"             = "simpson",
                "InvSimpson"          = "invsimpson",
                "ACE"                 = "ace"
              ),
              selected = "observed"
            )
          ),

          div(class = "data-select-item"),
          div(class = "data-select-item"),
          div(class = "data-select-item")
        )
      ),

      # ---- Alpha plot (full width row) ----
      fluidRow(
        column(
          width = 8,
          offset = 2,
          shinycssloaders::withSpinner(
            plotly::plotlyOutput("alpha_boxplot", height = "700px"),
            type = 4
          )
        )
      ),

      # ---- Beta plot controls ----
      div(
        class = "data-select-grid",

        div(class = "data-select-item"),   # empty first column

        div(
          class = "data-select-item",
          selectInput(
            "beta_metric",
            "Beta Diversity",
            choices = c(
              "Bray-Curtis"      = "bray",
              "Jaccard"          = "jaccard",
              "Euclidean"        = "euclidean",
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
          width = 9,
          offset = 2,
          shinycssloaders::withSpinner(
            plotly::plotlyOutput("beta_pcoa", height = "700px"),
            type = 4
          )
        )
      ),

      # ---- Taxonomic Pie Chart ----
      div(
        id = "sec_pie", class = "scroll-section",
        h3("Taxonomic Pie Chart"),

        shinycssloaders::withSpinner(
          taxplore::KronaChartOutput("tax_krona", height = "900px"),
          type = 4
        )
      )
    )
  )
}

app_b_server <- function(input, output, session) {

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
        samp_name = as.character(samp_name),
        year      = if ("year" %in% names(.)) as.character(year) else NA_character_,
        eventDate = if ("eventDate" %in% names(.)) as.character(eventDate) else NA_character_,
        sample_id = dplyr::case_when(
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

  make_species_cache_key <- function(selected_key, year, groups, sara_on, ais_on) {
    paste(
      selected_key,
      year,
      paste(sort(groups), collapse = "|"),
      paste0("sara=", sara_on),
      paste0("ais=", ais_on),
      sep = "~~"
    )
  }

  get_species_by_layer_cached <- function(selected_key, det_sf, year, groups, sara_on, ais_on) {

    cache_key <- make_species_cache_key(
      selected_key = selected_key,
      year         = year,
      groups       = groups,
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

  detections_filtered <- reactive({
    det <- selected_detections()
    if (is.null(det) || nrow(det) == 0) return(NULL)

    det %>%
      apply_species_filters()
  })

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

  species_gene_summary <- function(det_sf) {
    det_sf %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        scientificName = as.character(scientificName),
        target_gene    = as.character(target_gene),
        samp_name      = as.character(samp_name)
      ) %>%
      dplyr::filter(!is.na(scientificName), scientificName != "") %>%
      dplyr::group_by(scientificName) %>%
      dplyr::summarise(
        genes       = paste(sort(unique(na.omit(target_gene))), collapse = ", "),
        n_detections = dplyr::n(),
        n_samples    = dplyr::n_distinct(na.omit(samp_name)),
        .groups = "drop"
      ) %>%
      dplyr::arrange(scientificName)
  }

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

  # Return a named list: each name is a layer ("12S", "COI", ... or "All"),
  # each value is the species vector for that layer
  species_by_active_layers <- function(det_sf, layers_on, apply_filters_fn) {
    det_sf <- apply_filters_fn(det_sf)

    out <- list()

    # "All" means no gene filter
    if ("All" %in% layers_on) {
      spp_all <- det_sf %>%
        sf::st_drop_geometry() %>%
        dplyr::pull(scientificName) %>%
        as.character() %>%
        unique() %>%
        stats::na.omit() %>%
        sort()

      out[["All"]] <- spp_all
    }

    # gene-specific lists
    genes <- setdiff(layers_on, "All")
    for (g in genes) {
      spp_g <- det_sf %>%
        dplyr::filter(as.character(target_gene) == g) %>%
        sf::st_drop_geometry() %>%
        dplyr::pull(scientificName) %>%
        as.character() %>%
        unique() %>%
        stats::na.omit() %>%
        sort()

      out[[g]] <- spp_g
    }

    out
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

  div_gene_initialized <- reactiveVal(FALSE)
  div_primer_initialized <- reactiveVal(FALSE)

  # Detections that fall inside any MPA/AOI polygon (drops outside points)
  detections_in_mpa <- reactive({
    yr <- sel_year_chr()

    pts <- species_sf_all_with_poly

    if (yr != "All") {
      pts <- pts %>% dplyr::filter(as.character(year) == yr)
    }

    pts %>%
      dplyr::filter(!is.na(site_name), !is.na(site_type)) %>%
      dplyr::arrange(occurrenceID) %>%
      dplyr::group_by(occurrenceID, samp_name, scientificName, year, target_gene) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  })

  # redraw points when year changes
  observeEvent(
    list(sel_year_chr(), sampling_points_layer_on()),
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

      proxy %>%
        clearGroup("Sampling Points") %>%
        addCircleMarkers(
          data        = pts,
          group       = "Sampling Points",
          radius      = 2,
          stroke      = TRUE,
          weight      = 1,
          opacity     = 1,
          fillOpacity = 0.8,
          options     = pathOptions(pane = "pane_points"),
          label       = ~paste0(
            "Marker: ", target_gene,
            ifelse(is.na(year), "", paste0(" | Year: ", year)),
            ifelse(is.na(samp_name), "", paste0(" | Sample: ", samp_name))
          )
        )
    },
    ignoreInit = TRUE
  )

  #Monthly sampling
  monthly_sample_counts <- reactive({
    det <- selected_detections()
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

    det <- selected_detections()

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
    yr <- sel_year_chr()
    click <- input$map_shape_click
    sel_id <- selected_draw_id()

    # ---- drawn polygon: still needs spatial filtering ----
    if (!is.null(sel_id)) {
      pts <- if (yr == "All") {
        species_sf_all
      } else {
        species_sf_by_year[[yr]] %||% species_sf_all[0, ]
      }

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
      if (yr != "All") {
        pts <- pts %>% dplyr::filter(as.character(year) == yr)
      }

      return(
        pts %>%
          dplyr::filter(site_type == p_type, site_name == p_name)
      )
    }

    # ---- clicked grid cell: use precomputed lookup ----
    cid <- suppressWarnings(as.integer(id))
    if (!is.na(cid)) {
      pts <- species_sf_all_with_cell
      if (yr != "All") {
        pts <- pts %>% dplyr::filter(as.character(year) == yr)
      }

      return(
        pts %>%
          dplyr::filter(cell_id == cid)
      )
    }

    NULL
  })

  selected_detections_min <- reactive({
    yr <- sel_year_chr()
    click <- input$map_shape_click
    sel_id <- selected_draw_id()

    # drawn polygon still spatial
    if (!is.null(sel_id)) {
      pts <- species_sf_by_year[[yr]]
      if (is.null(pts)) pts <- species_sf_min[0, ]

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

    if (is.null(click) || is.null(click$id)) return(NULL)

    id <- as.character(click$id)

    pts <- species_sf_min
    if (yr != "All") {
      pts <- pts %>% dplyr::filter(as.character(year) == yr)
    }

    if (grepl("\\|\\|", id)) {
      parts <- strsplit(id, "\\|\\|")[[1]]
      p_type <- parts[1]
      p_name <- parts[2]

      return(
        pts %>%
          dplyr::left_join(point_poly_lookup, by = "occurrenceID") %>%
          dplyr::filter(site_type == p_type, site_name == p_name)
      )
    }

    cid <- suppressWarnings(as.integer(id))
    if (!is.na(cid)) {
      return(
        pts %>%
          dplyr::left_join(point_cell_lookup, by = "occurrenceID") %>%
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
      "Plants"        = list(col = "kingdom", vals = c("Plantae"))
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
  observeEvent(input$total_plants,     { toggle_group("Plants") }, ignoreInit = TRUE)

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
      total_plants     = "Plants"
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

  # ---- SARA/AIS sets + toggles (assumes SARA & AIS exist globally) ----
  sara_set <- reactive(unique(na.omit(SARA$Scientific.Name)))
  ais_set  <- reactive(unique(na.omit(AIS$Scientific.Name)))  # adjust column if needed

  filter_sara_on <- reactiveVal(FALSE)
  filter_ais_on  <- reactiveVal(FALSE)

  apply_species_filters <- function(occ_all) {

    # 1) group buttons
    occ_all <- apply_group_filter(occ_all, active_groups())

    # 2) SARA/AIS union logic
    if ("scientificName" %in% names(occ_all)) {
      spp_keep <- character(0)

      if (isTRUE(filter_sara_on())) spp_keep <- union(spp_keep, sara_set())
      if (isTRUE(filter_ais_on()))  spp_keep <- union(spp_keep, ais_set())

      if (length(spp_keep) > 0) {
        occ_all <- occ_all %>% dplyr::filter(scientificName %in% spp_keep)
      }
    }

    occ_all
  }

  diversity_dropdown_data <- reactive({
    pts <- species_sf_all

    yr <- sel_year_chr()
    if (yr != "All") {
      pts <- pts %>% dplyr::filter(as.character(year) == yr)
    }

    pts <- apply_species_filters(pts)

    pts %>%
      sf::st_drop_geometry() %>%
      add_primer_combo() %>%
      dplyr::mutate(
        target_gene  = as.character(target_gene),
        primer_combo = as.character(primer_combo)
      )
  })

  # ---- source data for diversity dropdowns ----
  observeEvent(diversity_dropdown_data(), {
    dd <- diversity_dropdown_data()

    gene_choices <- dd %>%
      dplyr::filter(!is.na(target_gene), target_gene != "") %>%
      dplyr::pull(target_gene) %>%
      unique() %>%
      sort()

    cur_gene <- input$div_target_gene %||% character(0)
    sel_gene <- intersect(cur_gene, gene_choices)

    # first load = select all genes
    if (!div_gene_initialized()) {
      sel_gene <- gene_choices
      div_gene_initialized(TRUE)
    }

    freezeReactiveValue(input, "div_target_gene")
    updateSelectizeInput(
      session  = session,
      inputId  = "div_target_gene",
      choices  = gene_choices,
      selected = sel_gene,
      server   = TRUE
    )
  }, ignoreInit = FALSE)

  primer_choices_reactive <- reactive({
    dd <- diversity_dropdown_data()

    genes_selected <- input$div_target_gene %||% character(0)

    if (length(genes_selected) > 0) {
      dd <- dd %>%
        dplyr::filter(target_gene %in% genes_selected)
    } else {
      dd <- dd[0, , drop = FALSE]
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

      cur_primer <- input$div_primer %||% character(0)
      sel_primer <- intersect(cur_primer, primer_choices)

      # first load = select all available primers
      if (!div_primer_initialized()) {
        sel_primer <- primer_choices
        div_primer_initialized(TRUE)
      }

      freezeReactiveValue(input, "div_primer")
      updateSelectizeInput(
        session  = session,
        inputId  = "div_primer",
        choices  = primer_choices,
        selected = sel_primer,
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

  # ---- confirmed diversity controls ----
  div_filters <- reactiveVal(
    list(
      tax_rank    = "scientificName",
      target_gene = character(0),
      primers     = character(0)
    )
  )

  observeEvent(input$div_apply, {
    div_filters(
      list(
        tax_rank    = input$tax_rank %||% "scientificName",
        target_gene = input$div_target_gene %||% character(0),
        primers     = input$div_primer %||% character(0)
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

        labs <- paste0(
          all_polys_click$site_type[hit],
          ": ",
          all_polys_click$site_name[hit]
        )

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
    occ_all <- selected_detections()   # or occ_all_filtered(), etc.
    req(occ_all)

    occ_all <- apply_species_filters(occ_all)

    occ_all %>% dplyr::distinct(scientificName, .keep_all = TRUE) %>% dplyr::arrange(dplyr::coalesce(worms_valid_name, scientificName))
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
    if (!isTRUE(filter_sara_on()) && !isTRUE(filter_ais_on())) {
      return(spp_vec)
    }

    keep <- character(0)

    if (isTRUE(filter_sara_on())) keep <- union(keep, sara_set())
    if (isTRUE(filter_ais_on()))  keep <- union(keep, ais_set())

    intersect(spp_vec, keep)
  }

  active_filters_label <- reactive({
    labs <- c()
    if (isTRUE(filter_sara_on())) labs <- c(labs, "SARA")
    if (isTRUE(filter_ais_on()))  labs <- c(labs, "AIS")
    if (length(labs) == 0) "None" else paste(labs, collapse = " + ")
  })

  output$sara_details <- DT::renderDT({
    if (!isTRUE(filter_sara_on())) {
      return(DT::datatable(
        data.frame(Message = "Select the “SARA” button to view species at risk details."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det <- selected_detections()
    if (is.null(det) || nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections in the current selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det <- det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(scientificName = as.character(scientificName)) %>%
      dplyr::filter(scientificName %in% sara_set())

    if (nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No SARA Schedule 1 detections for this selection."),
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

    det <- selected_detections()
    if (is.null(det) || nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections in the current selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det <- det %>%
      sf::st_drop_geometry() %>%
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
      "Click Confirm before downloading."
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

        if (is.null(dat) || nrow(dat) == 0) {
          message("No rows left; writing empty CSV")
          utils::write.csv(data.frame(), file, row.names = FALSE, na = "")
          return()
        }

        message("Adding polygon selection")
        dat <- add_polygon_selection(dat)
        message("Polygon selection added")

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

        rank_col <- filters$tax_rank %||% "scientificName"
        message("Rank column: ", rank_col)


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

  # ---- detections used by Diversity section only ----
  # Floating panel filters still apply first.
  # tax_rank and div_target_gene only affect diversity calculations.
  # diversity_detections <- reactive({
  #   det <- selected_detections()
  #   req(det)
  #
  #   det <- apply_species_filters(det)
  #
  #   gene_sel <- input$div_target_gene %||% character(0)
  #
  #   # If user deselects everything, keep all genes
  #   if (length(gene_sel) > 0) {
  #     det <- det %>%
  #       dplyr::filter(as.character(target_gene) %in% gene_sel)
  #   }
  #
  #   det
  # })


  # ---- diversity detections for all MPA/AOI polygons ----
  diversity_detections_mpa <- reactive({
    yr <- sel_year_chr()
    filters <- div_filters()

    pts <- species_sf_all_with_poly

    if (yr != "All") {
      pts <- pts %>% dplyr::filter(as.character(year) == yr)
    }

    if (is.null(pts) || nrow(pts) == 0) {
      return(pts)
    }

    pts <- apply_species_filters(pts)
    pts <- apply_diversity_dropdown_filters(pts, filters)

    if (nrow(pts) == 0) {
      return(pts)
    }

    pts %>%
      dplyr::filter(!is.na(site_name), !is.na(site_type)) %>%
      dplyr::arrange(occurrenceID) %>%
      dplyr::group_by(
        occurrenceID, samp_name, scientificName, year, target_gene,
        site_name, site_type
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
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

    rank_col <- div_filters()$tax_rank

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
    pts <- species_sf_all
    req(pts)

    occ_all <- pts %>%
      sf::st_drop_geometry() %>%
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

  rarefaction_depth <- reactive({
    5000
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

  observe({
    rr <- rarefaction_drop_summary()
    req(rr)
    print("running rarefaction_drop_summary")
    message("Rarefaction depth: ", rr$depth)
    message("Total samples in current alpha matrix: ", rr$total_samples)
    message("Kept after rarefaction filter: ", rr$kept_samples)
    message("Dropped before rarefaction: ", rr$dropped_samples)
  })

  #Check which samples were dropped
  observe({
    rr <- rarefaction_drop_summary()
    req(rr)

    if (nrow(rr$dropped_table) > 0) {
      print(rr$dropped_table)
    } else {
      message("No samples dropped at this rarefaction depth.")
    }
  })

  #Checks for ordination (beta)
  observe({
    mat <- beta_comm_mat()
    req(mat)

    lib_sizes <- rowSums(mat, na.rm = TRUE)

    message("BETA total samples: ", nrow(mat))
    message("BETA zero-sum samples: ", sum(lib_sizes == 0, na.rm = TRUE))
    message("BETA min library size: ", min(lib_sizes, na.rm = TRUE))
    message("BETA median library size: ", median(lib_sizes, na.rm = TRUE))
    message("BETA max library size: ", max(lib_sizes, na.rm = TRUE))
  })

  #Create rarefied matrix for alpha diversity
  comm_mat_mpa_rarefied <- reactive({
    mat <- comm_mat_mpa()
    req(mat)

    depth <- rarefaction_depth()

    shiny::validate(
      shiny::need(nrow(mat) > 0, "No samples available for rarefaction."),
      shiny::need(depth > 0, "Rarefaction depth must be greater than 0.")
    )

    lib_sizes <- rowSums(mat, na.rm = TRUE)

    keep <- lib_sizes >= depth
    mat  <- mat[keep, , drop = FALSE]

    shiny::validate(
      shiny::need(nrow(mat) > 0,
                  "No filtered samples have enough reads to be rarefied at the depth.")
    )

    mat <- round(as.matrix(mat))
    storage.mode(mat) <- "integer"

    set.seed(123)
    vegan::rrarefy(mat, sample = depth)
  })


  sample_meta_mpa <- reactive({
    det <- diversity_detections_mpa()
    req(det)

    keep_ids <- rownames(comm_mat_mpa_rarefied())

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
      dplyr::filter(sample_id %in% keep_ids) %>%
      dplyr::distinct(sample_id, samp_name, site_name, site_type, year, group_label)
  })

  # ---- base detections used for beta ordination ----
  # Includes:
  #   - all detections inside MPA/AOI polygons
  #   - plus detections inside any drawn polygons
  # Floating-panel filters and diversity target_gene filter still apply.
  diversity_detections_beta <- reactive({
    yr <- sel_year_chr()
    filters <- div_filters()

    pts <- species_sf_all_with_poly

    if (yr != "All") {
      pts <- pts %>% dplyr::filter(as.character(year) == yr)
    }

    if (is.null(pts) || nrow(pts) == 0) {
      return(pts)
    }

    pts <- apply_species_filters(pts)
    pts <- apply_diversity_dropdown_filters(pts, filters)

    if (nrow(pts) == 0) {
      return(pts)
    }

    # ensure columns exist
    if (!"site_name" %in% names(pts)) pts$site_name <- NA_character_
    if (!"site_type" %in% names(pts)) pts$site_type <- NA_character_

    pts <- pts %>%
      dplyr::mutate(
        site_name = as.character(site_name),
        site_type = as.character(site_type)
      )

    in_mpa <- pts %>%
      dplyr::filter(!is.na(site_name), site_name != "",
                    !is.na(site_type), site_type != "")

    if (nrow(in_mpa) > 0) {
      in_mpa <- in_mpa %>%
        dplyr::arrange(occurrenceID) %>%
        dplyr::group_by(
          occurrenceID, samp_name, scientificName, year, target_gene,
          site_name, site_type
        ) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
    }

    polys <- drawn_polys()

    if (is.null(polys) || nrow(polys) == 0) {
      return(in_mpa)
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

    in_drawn <- dplyr::bind_rows(drawn_list)

    dplyr::bind_rows(in_mpa, in_drawn)
  })

  # ---- unique sample x taxon matrix for ONE shared ordination ----
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

    rank_col <- div_filters()$tax_rank

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

  # --- sample metadata for grouping/hover (Location etc.) ---
  # sample_meta <- reactive({
  #   det <- selected_detections()
  #   req(det)
  #
  #   occ_all <- det %>% sf::st_drop_geometry()
  #
  #   out <- occ_all %>%
  #     dplyr::mutate(
  #       samp_name = as.character(samp_name),
  #       site_name  = if ("site_name" %in% names(.)) as.character(site_name) else NA_character_,
  #       year      = if ("year" %in% names(.)) as.character(year) else NA_character_
  #     ) %>%
  #     dplyr::distinct(samp_name, site_name, year)
  #
  #   out
  # })

  #Alpha diversity
  alpha_metric_vec <- reactive({
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
  observe({
    df <- library_sizes_mpa()
    req(nrow(df) > 0)

    message("Min library size: ", min(df$library_size, na.rm = TRUE))
    message("Median library size: ", median(df$library_size, na.rm = TRUE))
    message("Max library size: ", max(df$library_size, na.rm = TRUE))
  })

  observe({
    depth <- rarefaction_depth()
    req(depth)

    message("Rarefaction depth used for alpha diversity: ", depth)
  })

  # --- Alpha boxplot data: selected area + optional drawn polygons ---
  alpha_boxplot_occ_all <- reactive({
    yr <- sel_year_chr()

    alpha <- alpha_metric_vec()
    meta  <- sample_meta_mpa()

    alpha_mpa <- alpha %>%
      dplyr::left_join(meta, by = "sample_id") %>%
      dplyr::filter(!is.na(site_name), site_name != "") %>%
      dplyr::mutate(group_label = as.character(group_label))

    polys <- drawn_polys()

    if (!is.null(polys) && nrow(polys) > 0) {

      pts_analysis <- species_sf_all
      if (yr != "All") {
        pts_analysis <- pts_analysis %>%
          dplyr::filter(as.character(year) == yr)
      }

      pts_analysis <- pts_analysis %>%
        apply_species_filters() %>%
        apply_diversity_dropdown_filters(div_filters())

      rank_col <- div_filters()$tax_rank
      metric   <- input$alpha_metric %||% "observed"

      alpha_poly_list <- lapply(seq_len(nrow(polys)), function(i) {
        g_i <- sf::st_geometry(polys[i, , drop = FALSE])

        depth_poly <- rarefaction_depth()

        if (!is.finite(depth_poly) || depth_poly <= 0) return(NULL)

        inside <- pts_analysis[within_any(pts_analysis, g_i), , drop = FALSE] %>%
          sf::st_drop_geometry() %>%
          make_sample_id() %>%
          dplyr::mutate(
            samp_name = as.character(samp_name),
            taxon     = as.character(.data[[rank_col]]),
            value     = as.numeric(organismQuantity)
          ) %>%
          dplyr::filter(
            !is.na(sample_id), sample_id != "",
            !is.na(taxon), taxon != "",
            !is.na(value)
          ) %>%
          dplyr::group_by(sample_id, samp_name, taxon) %>%
          dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

        if (nrow(inside) == 0) return(NULL)

        mat_poly <- inside %>%
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

        if (nrow(mat_poly) == 0) return(NULL)

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
          vegan::specnumber(mat_poly_rarefied)
        )

        poly_label <- polys$draw_label[i]

        data.frame(
          sample_id   = rownames(mat_poly_rarefied),
          alpha_val   = as.numeric(vals_poly),
          samp_name   = sample_ids_poly$samp_name[
            match(rownames(mat_poly_rarefied), sample_ids_poly$sample_id)
          ],
          site_name   = poly_label,
          site_type   = "User",
          year        = NA_character_,
          group_label = poly_label,
          stringsAsFactors = FALSE
        )
      })

      alpha_poly <- dplyr::bind_rows(alpha_poly_list)

      if (nrow(alpha_poly) > 0) {
        alpha_mpa <- dplyr::bind_rows(alpha_mpa, alpha_poly)
      }
    }

    base_lvls <- unique(alpha_mpa$group_label[alpha_mpa$site_type != "User"])
    draw_lvls <- unique(alpha_mpa$group_label[alpha_mpa$site_type == "User"])

    alpha_mpa %>%
      dplyr::mutate(group_label = factor(group_label, levels = c(base_lvls, draw_lvls)))
  })

  # --- Alpha diversity boxplot (Plotly) ---
  color_vec <- c("#046c9a", "#5BBCD6", "#ABDDDE", "#446455", "#00A08A","#Fdd262")

  output$alpha_boxplot <- plotly::renderPlotly({
    alpha <- alpha_boxplot_occ_all()

    shiny::validate(
      shiny::need(nrow(alpha) > 0, "No samples available for the current selection/year.")
    )

    metric_label <- switch(
      input$alpha_metric %||% "observed",
      observed   = "Richness (count of taxa)",
      shannon    = "Shannon diversity",
      simpson    = "Simpson diversity",
      invsimpson = "Inverse Simpson diversity",
      ace        = "ACE estimated richness",
      "Alpha diversity"
    )

    rank_label <- tools::toTitleCase(gsub("_", " ", div_filters()$tax_rank))
    gene_label <- div_filters()$target_gene
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


    plotly::plot_ly(
      data = alpha,
      x = ~group_label,
      y = ~alpha_val,
      type = "box",
      color = ~group_label,
      colors = color_vec,
      boxpoints = "all",
      jitter = 0.3,
      pointpos = 0,
      hovertemplate = paste(
        "<b>Sample:</b> %{customdata[0]}<br>",
        "<b>Group:</b> %{x}<br>",
        "<b>", metric_label, ":</b> %{y}<extra></extra>"
      ),
      customdata = ~cbind(samp_name)
    ) %>%
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
          range = c(0, max(alpha$alpha_val, na.rm = TRUE) * 1.05),
          fixedrange = FALSE
        ),
        margin = list(l = 80, r = 30, t = 40, b = 120),
        showlegend = FALSE
      )
  })

  #Beta diversity
  output$beta_pcoa <- plotly::renderPlotly({
    mat  <- beta_comm_mat()
    meta <- beta_plot_meta()

    shiny::validate(
      shiny::need(nrow(mat) > 2, "At least 3 samples are required to compute a PCoA.")
    )


    beta_method <- input$beta_metric %||% "bray"
    mat_use <- as.matrix(mat)

    if (beta_method %in% c("bray", "jaccard")) {
      rs <- rowSums(mat_use, na.rm = TRUE)
      mat_rel <- mat_use
      mat_rel[rs > 0, ] <- mat_use[rs > 0, , drop = FALSE] / rs[rs > 0]

      if (beta_method == "jaccard") {
        mat_dist_input <- (mat_rel > 0) * 1
      } else {
        mat_dist_input <- mat_rel
      }

      d <- vegan::vegdist(mat_dist_input, method = beta_method)

    } else if (beta_method == "euclidean") {
      d <- stats::dist(mat_use, method = "euclidean")

    } else if (beta_method == "robust.aitchison") {
      shiny::validate(
        shiny::need(all(mat_use >= 0, na.rm = TRUE),
                    "Robust Aitchison distance requires non-negative values.")
      )

      d <- vegan::vegdist(mat_use, method = "robust.aitchison")
    }


    ord <- stats::cmdscale(d, k = 2, eig = TRUE)

    sample_ids <- rownames(mat_use)

    shiny::validate(
      shiny::need(!is.null(sample_ids), "Sample names are missing from the beta diversity matrix."),
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
      bray             = "Bray-Curtis",
      jaccard          = "Jaccard",
      euclidean        = "Euclidean",
      `robust.aitchison` = "Robust Aitchison",
      beta_method
    )

    rank_label <- tools::toTitleCase(gsub("_", " ", div_filters()$tax_rank))
    gene_label <- div_filters()$target_gene
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

    # Build hover text explicitly so it always matches plot_df row count
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

    plotly::plot_ly(
      data = plot_df,
      x = ~PC1,
      y = ~PC2,
      type = "scatter",
      mode = "markers",
      color = ~group_label,
      colors = color_vec,
      text = ~hover_text,
      hoverinfo = "text"
    ) %>%
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
          font = list(size = 16)
        ),
        margin = list(l = 80, r = 30, t = 70, b = 80),
        showlegend = TRUE
      )
  })


  #Krona plot
  output$tax_krona <- taxplore::renderKronaChart({
    det <- selected_detections()
    req(det)

    # floating-panel filters
    det <- apply_species_filters(det)

    # data-selection filters: target gene + primer only
    det <- apply_diversity_dropdown_filters(
      det,
      list(
        target_gene = div_filters()$target_gene,
        primers     = div_filters()$primers
      )
    )

    shiny::validate(
      shiny::need(!is.null(det) && nrow(det) > 0,
                  "Select a cell/polygon (or draw a polygon) to display a Krona chart.")
    )

    det0 <- det %>%
      sf::st_drop_geometry()

    shiny::validate(
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
    )
  })


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
      sara_on      = isTRUE(filter_sara_on()),
      ais_on       = isTRUE(filter_ais_on())
    )

    by_layer <- by_layer_all[intersect(names(by_layer_all), layers_on)]
    render_by_layer_ui(by_layer, strong(paste0("Grid cell: ", cid)), 220)
  })

  output$detections_tbl <- DT::renderDT({
    det <- selected_detections()

    if (is.null(det) || nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "Click an MPA/AOI polygon/grid cell, or draw a polygon, to view detection details"),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det_f <- apply_species_filters(det)

    if (nrow(det_f) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections match the current Year and/or Group filters."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det_show <- det_f %>%
      dplyr::mutate(
        decimalLongitude = sf::st_coordinates(.)[, 1],
        decimalLatitude  = sf::st_coordinates(.)[, 2]
      ) %>%
      sf::st_drop_geometry()

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
        ), ~as.character(.)),
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
