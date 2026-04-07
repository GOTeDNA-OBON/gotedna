library(GOTeDNA)
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)
library(plotly)
library(DT)
library(stringr)

# library(leaflet)
# library(sf)
# library(shiny)
# library(shinyjs)
# library(bslib)
# library(plotly)
cli::cli_alert_info("Packages loaded")

list.files("modules", full.names = TRUE) |>
  lapply(source)
cli::cli_alert_info("Modules loaded")

options(shiny.maxRequestSize = 100 * 1024^2)
# Ensure that the user installs these if not already installed
# Necessary to generate PDF report


# if (!tinytex::is_tinytex() == TRUE) {
#   tinytex::install_tinytex()
# }


# Generic helpers
trans_letters <- function(x, pos = 1, fun = toupper) {
  strsplit(x, split = "") |>
    lapply(\(y) {
      y[pos] <- fun(y[pos])
      paste(y, collapse = "")
    }) |>
    unlist()
}

# Data
## import glossary
gloss <- read.csv("data/glossary.csv")
gloss$Term <- paste0('<p align ="right"><b>', trimws(gloss$Term), "</b></p>")
gloss$Definition <- trimws(gloss$Definition)

## import GOTeDNA data
gotedna_data <- gotedna_data0 <- readRDS("data/gotedna_data.rds")
last_obis_download_ts <- as.POSIXct(readLines("data/last_obis_download_ts.txt"), tz = Sys.timezone())


gotedna_station <- gotedna_station0 <- readRDS("data/gotedna_station.rds")
gotedna_primer <- readRDS("data/gotedna_primer.rds")

#taxonomic_ranks <- list("kingdom", "phylum", "family", "order", "class", "genus")
#taxonomic_ranks <- list("domain", "kingdom", "phylum", "class", "order", "family", "genus")
taxonomic_ranks <- list("kingdom", "phylum", "class", "order", "family", "genus")
names(taxonomic_ranks) <- trans_letters(taxonomic_ranks |> unlist())

ls_threshold <- as.list(seq(50, 95, 5))
names(ls_threshold) <- paste0(seq(50, 95, 5), "%")


# function

basemap <- function() {
  leaflet::leaflet() |>
    leafem::addMouseCoordinates() |>
    #  leaflet::addProviderTiles("Esri.OceanBasemap", group = "OceaBasemap") |>
    leaflet::addProviderTiles("OpenStreetMap", group = "OpenStreetMap") |>
    #  leaflet::addProviderTiles("Stadia.StamenTonerLite") |>
    # leaflet::addProviderTiles("Esri.WorldGrayCanvas") |>

    # leaflet::addLayersControl(
    #  baseGroups = c("OpenStreetMap", "Ocean Basemap"),
    #  position = "bottomleft"
    #  ) |>
    leaflet::addScaleBar(
      position = c("bottomright"),
      options = leaflet::scaleBarOptions(maxWidth = 200)
    ) |>
    leaflet.extras::addDrawToolbar(
      polylineOptions = FALSE,
      circleOptions = FALSE,
      polygonOptions = TRUE,
      rectangleOptions = TRUE,
      circleMarkerOptions = FALSE,
      markerOptions = FALSE
    )
}


get_primer_selection <- function(lvl, data, primer_sheet = gotedna_primer) {
  if (is.null(lvl)) {
    return("Not available")
  } else {
    if (lvl == "scientificName") {
      tx_col <- "scientificName"
    } else {
      tx_col <- lvl
    }
    data_available <- data |>
      select({{ tx_col }}, primer) |>
      distinct()
    if (!is.null(data_available) && nrow(data_available)) {
      tmp <- primer_sheet[[lvl]] |>
        inner_join(
          data_available,
          join_by(primer == primer, {{ lvl }} == {{ tx_col }})
        ) |>
        mutate(
          text = paste0(primer, " (", detects, "/", total, " ; ", perc, "%)")
        )
      out <- as.list(tmp$primer)
      names(out) <- tmp$text
      if (is.null(out)) {
        return("Not available")
      } else {
        return(out)
      }
    } else {
      return("Not available")
    }
  }
}

get_station <- function(x) {
  x |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(decimalLongitude)) |>
    dplyr::filter(!is.na(phylum)) |>
    dplyr::select(
      c(decimalLongitude, decimalLatitude, station)
    ) |>
    dplyr::distinct() |>
    dplyr::group_by(station) |>
    dplyr::summarise(
      decimalLongitude = mean(as.numeric(decimalLongitude)),
      decimalLatitude = mean(as.numeric(decimalLatitude))
    ) |>
    dplyr::ungroup() |>
    as.data.frame() |>
    sf::st_as_sf(
      coords = c("decimalLongitude", "decimalLatitude"),
      crs = sf::st_crs(4326)
    )
}


# Primer information for primer tab
## import glossary
primer_seqs <- read.csv("data/primers.csv") |>
  dplyr::rename(
    "Primer set" = "PrimerSet",
    "Type of data" = "Data",
    "Original primer name" = "OriginalName",
    "Sequence (5'-3')" = "Seq",
    "Fragment length (bp)" = "bp"
  )


big_OBIS_data_pull <- function(dataset_ids = NULL) {
  D_mb <- read_data(
    dataset_ids    = dataset_ids,
    scientificname = NULL,
    worms_id       = NULL,
    areaid         = NULL,
    join_by        = c("auto", "occurrenceID", "id"),
    require_absences = TRUE
  )

  D_mb_clean <- add_detected_column(D_mb)
  D_mb_clean
}


default_protocol_info <- readRDS("data/protocol_sheet.rds")

protocol_labels <- c(
  samp_size = "Sample Volume (L)",
  samp_size_mid = "Sample Volume (L)",
  size_frac = "Filter Pore Size",
  filter_material = "Filter Material",
  samp_mat_process = "Sample Processing Method",
  minimumDepthInMeters = "Minimum Depth (m)",
  maximumDepthInMeters = "Maximum Depth (m)",
  min_depth_floor = "Minimum Depth (m)",
  max_depth_floor = "Maximum Depth (m)",

  samp_store_temp = "Storage Temperature (°C)",
  samp_store_sol = "Storage Solution",

  target_gene = "Target Gene",
  pcr_primer_forward = "Forward Primer",
  pcr_primer_reverse = "Reverse Primer",
  nucl_acid_ext_kit = "DNA Extraction Kit",

  platform = "Sequencing Platform",
  instrument = "Instrument",
  seq_kit = "Sequencing Kit",

  otu_db = "Reference Database",
  tax_assign_cat = "Taxonomic Assignment Method",
  otu_seq_comp_appr = "OTU/ASV Approach"
)

protocol_columns <- c(
  'nucl_acid_ext_kit',
  'platform',
  'instrument',
  'seq_kit',
  'otu_db',
  'tax_assign_cat',
  'otu_seq_comp_appr',
  'min_depth_floor',
  'max_depth_floor',
  'samp_size_mid',
  'size_frac',
  'filter_material',
  'samp_mat_process',
  'samp_store_temp',
  'samp_store_sol'
)



required_cols <- c(
  "samp_name",
  "target_gene",
  "pcr_primer_name_forward",
  "pcr_primer_name_reverse",
  "scientificName",
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "eventDate_clean",
  "decimalLatitude",
  "decimalLongitude",
  "organismQuantity"
)

optional_columns <- c(
  'samp_size',
  'size_frac',
  'filter_material',
  'samp_mat_process',
  'samp_store_temp',
  'samp_store_sol',
  'project_contact',
  'nucl_acid_ext_kit',
  'platform',
  'LClabel',
  'instrument',
  'month',
  'year',
  'seq_kit',
  'otu_db',
  'tax_assign_cat',
  'otu_seq_comp_appr',
  'minimumDepthInMeters',
  'maximumDepthInMeters',
  'category',
  'hab'
)

