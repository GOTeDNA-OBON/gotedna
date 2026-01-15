library(GOTeDNA)
library(dplyr)
library(ggplot2)
library(patchwork)
library(plotly)
library(DT)
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


#Code for testing LClabel text ASDF
gotedna_data$metabarcoding <- gotedna_data$metabarcoding %>%
  # create a numeric label for each unique datasetID_obis
  group_by(datasetID_obis) %>%
  mutate(dataset_idx = cur_group_id()) %>%
  ungroup() %>%
  # create the LClabel text without wrapping <div>
  mutate(LClabel = paste0(
    "<p>", dataset_idx, ": We acknowledge that this research was conducted on the unceded and unsurrendered traditional territories of the
            Mi'kmaq, Passamaquoddy, Wabanaki (Bay of Fundy ecodistrict); Mi’kmaq, Wabanaki (Magdalen Shallows, Scotian Shelf);
            Métis peoples (Churchill Estuary-Hudson Bay); and the Inuit homelands of Inuit Nunangat including: Nunatsiavut
            (NL-Labrador Shelves), Nunavik (Southern Hudson Strait), and Nunavut (Baffin Bay/Davis Strait, North Baffin Fjords).</p>
    <p>Arctic data collections were funded by ArcticNet, Polar Knowledge Canada, DFO (Aquatic Invasive Species Monitoring
            Program, Strategic Program for Ecosystem-Based Research and Advice [SPERA], Ocean Protection Plan [OPP] Coastal
            Environmental Baseline Program, Arctic Science Funds, GRDI), Nunavut Wildlife Management Board, Nunavik Marine
            Region Wildlife Board, and World Wildlife Funds to KH and ALR. Field accommodations and/or logistic support during
            Arctic field campaigns were provided by Churchill Northern Studies Centre, Glencore-Raglan, Baffinland Iron and Vale
            Mines, NRCan (Polar Continental Shelf program), Nunatsiavut Government, Environment and Climate Change Canada, and
            Government of Nunavut MV Nuliajuk.</p>"
  )) %>%
  select(-dataset_idx)

gotedna_data0 <- gotedna_data

#gotedna_data$metabarcoding <- readRDS("data/test_obis_animalia.rds")

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

  D_mb_msct <- D_mb %>%
    dplyr::mutate(msct = case_when(
      organismQuantity == 0 ~ TRUE,
      organismQuantity > 10 ~ TRUE
    )) |>
    tidyr::drop_na(msct)

  D_mb_nodetect <- D_mb_msct %>%
    dplyr::group_by(
      protocol_ID, protocolVersion, scientificName, primer, station) %>%
    dplyr::summarise(num_detected = sum(detected)) %>%
    dplyr::filter(num_detected == 0)

  D_mb_clean <- dplyr::anti_join(D_mb_msct, D_mb_nodetect,
                                 by = c("protocol_ID","protocolVersion","scientificName",
                                        "primer", "station"))
  D_mb_clean
}
