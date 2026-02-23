library(GOTeDNA)
library(dplyr)
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


# #Code for testing LClabel text ASDF
# gotedna_data$metabarcoding <- gotedna_data$metabarcoding %>%
#   # create a numeric label for each unique datasetID_obis
#   group_by(datasetID_obis) %>%
#   mutate(dataset_idx = cur_group_id()) %>%
#   ungroup() %>%
#   # create the LClabel text without wrapping <div>
#   mutate(LClabel =
#     "<p>We acknowledge that this research was conducted on the unceded and unsurrendered traditional territories of the
#             Mi'kmaq, Passamaquoddy, Wabanaki (Bay of Fundy ecodistrict); Mi’kmaq, Wabanaki (Magdalen Shallows, Scotian Shelf);
#             Métis peoples (Churchill Estuary-Hudson Bay); and the Inuit homelands of Inuit Nunangat including: Nunatsiavut
#             (NL-Labrador Shelves), Nunavik (Southern Hudson Strait), and Nunavut (Baffin Bay/Davis Strait, North Baffin Fjords).</p>
#     <p>Arctic data collections were funded by ArcticNet, Polar Knowledge Canada, DFO (Aquatic Invasive Species Monitoring
#             Program, Strategic Program for Ecosystem-Based Research and Advice [SPERA], Ocean Protection Plan [OPP] Coastal
#             Environmental Baseline Program, Arctic Science Funds, GRDI), Nunavut Wildlife Management Board, Nunavik Marine
#             Region Wildlife Board, and World Wildlife Funds to KH and ALR. Field accommodations and/or logistic support during
#             Arctic field campaigns were provided by Churchill Northern Studies Centre, Glencore-Raglan, Baffinland Iron and Vale
#             Mines, NRCan (Polar Continental Shelf program), Nunatsiavut Government, Environment and Climate Change Canada, and
#             Government of Nunavut MV Nuliajuk.</p>"
#   ) %>%
#   select(-dataset_idx)
#
# gotedna_data0 <- gotedna_data

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


big_OBIS_data_pull <- function(dataset_ids = c("74b70871-91bd-4b74-91ce-34e9611ce27d", "858f3eb9-7fee-4764-bf7f-04098922f162")) {
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
  'max_depth_floor'
)

version_columns <- c(
  'samp_size_mid',
  'size_frac',
  'filter_material',
  'samp_mat_process',
  'samp_store_temp',
  'samp_store_sol'
)


#
# library(tibble)
#
# raw_protocol_test_sheet <- tibble(
#
#   # --------------------------
#   # PROTOCOL 1 (3 versions)
#   # --------------------------
#   nucl_acid_ext_kit = c(rep("Qiagen", 9),
#                         rep("Zymo", 8),
#                         rep("Macherey-Nagel", 7),
#                         rep("Invitrogen", 6),
#                         rep("Promega", 8)),
#
#   platform = NA_character_,
#
#   instrument = c(rep("MiSeq", 9),
#                  rep("NextSeq", 8),
#                  rep("MiSeq", 7),
#                  rep("NovaSeq", 6),
#                  rep("NextSeq", 8)),
#
#   seq_kit = c(rep("KitA", 9),
#               rep("KitB", 8),
#               rep("KitC", 7),
#               rep("KitD", 6),
#               rep("KitE", 8)),
#
#   otu_db = c(rep("DB1", 9),
#              rep("DB2", 8),
#              rep("DB3", 7),
#              rep("DB4", 6),
#              rep("DB5", 8)),
#
#   tax_assign_cat = c(rep("SINTAX", 9),
#                      rep("RDP", 8),
#                      rep("SINTAX", 7),
#                      rep("BLAST", 6),
#                      rep("RDP", 8)),
#
#   otu_seq_comp_appr = c(rep("vsearch", 9),
#                         rep("blastn", 8),
#                         rep("usearch", 7),
#                         rep("vsearch", 6),
#                         rep("blastn", 8)),
#
#   minimumDepthInMeters = c(rep(5, 9),
#                            rep(10, 8),
#                            rep(20, 7),
#                            rep(15, 6),
#                            rep(30, 8)),
#
#   maximumDepthInMeters = c(rep(15, 9),
#                            rep(20, 8),
#                            rep(30, 7),
#                            rep(25, 6),
#                            rep(40, 8)),
#
#   # --------------------------
#   # VERSION COLUMNS
#   # --------------------------
#
#   # Protocol 1 → 3 versions
#   # Protocol 2 → 2 versions
#   # Protocol 3 → 4 versions
#   # Protocol 4 → 2 versions
#   # Protocol 5 → 3 versions
#
#   samp_size = c(
#     # Protocol 1 (3 versions)
#     rep(0.25, 3), rep(0.50, 3), rep(1.00, 3),
#
#     # Protocol 2 (2 versions)
#     rep(0.25, 4), rep(0.75, 4),
#
#     # Protocol 3 (4 versions)
#     rep(0.10, 2), rep(0.25, 2),
#     rep(0.50, 2), rep(1.00, 1),
#
#     # Protocol 4 (2 versions)
#     rep(0.50, 3), rep(1.50, 3),
#
#     # Protocol 5 (3 versions)
#     rep(0.25, 3), rep(0.75, 3), rep(1.25, 2)
#   ),
#
#   size_frac = sample(c("0.2 µm", "0.45 µm", "0.8 µm"), 38, replace = TRUE),
#   filter_material = sample(c("GF/F", "Cellulose", "Polycarbonate"), 38, replace = TRUE),
#   samp_mat_process = sample(c("Frozen", "Ethanol", "Lyophilized"), 38, replace = TRUE),
#   samp_store_temp = sample(c("-20C", "-80C", "4C"), 38, replace = TRUE),
#   samp_store_sol = sample(c("None", "RNAlater", "Buffer"), 38, replace = TRUE)
# )

#
# assign_protocol_ID <- function(df,
#                                protocol_columns,
#                                version_columns,
#                                protocol_sheet = NULL) {
#
#   # Remove existing IDs if present
#   df <- df %>%
#     select(-any_of(c("protocol_ID", "protocol_version")))
#
#   # Distinct combinations at full granularity
#   new_full_combos <- df %>%
#     select(all_of(c(protocol_columns, version_columns))) %>%
#     distinct()
#
#   # ------------------------------------------------------------------
#   # CASE 1: No existing protocol_sheet → build from scratch
#   # ------------------------------------------------------------------
#   if (is.null(protocol_sheet) || nrow(protocol_sheet) == 0) {
#
#     protocol_sheet <- new_full_combos %>%
#       group_by(across(all_of(protocol_columns))) %>%
#       mutate(
#         protocol_ID = cur_group_id(),
#         protocol_version = row_number()
#       ) %>%
#       ungroup()
#
#   } else {
#
#     # Ensure sheet has required structure
#     required_cols <- c(protocol_columns, version_columns,
#                        "protocol_ID", "protocol_version")
#
#     missing_cols <- setdiff(required_cols, names(protocol_sheet))
#     if (length(missing_cols) > 0) {
#       stop("protocol_sheet is missing required columns: ",
#            paste(missing_cols, collapse = ", "))
#     }
#
#     # --------------------------------------------------------------
#     # STEP 1: Add new protocol_IDs if needed
#     # --------------------------------------------------------------
#
#     new_protocols <- anti_join(
#       new_full_combos %>% select(all_of(protocol_columns)) %>% distinct(),
#       protocol_sheet %>% select(all_of(protocol_columns)) %>% distinct(),
#       by = protocol_columns
#     )
#
#     if (nrow(new_protocols) > 0) {
#
#       max_id <- max(protocol_sheet$protocol_ID, na.rm = TRUE)
#
#       new_protocols <- new_protocols %>%
#         mutate(protocol_ID = row_number() + max_id)
#
#       # give them version 1 initially (will expand below if needed)
#       new_protocols <- new_protocols %>%
#         left_join(new_full_combos, by = protocol_columns) %>%
#         group_by(protocol_ID) %>%
#         mutate(protocol_version = row_number()) %>%
#         ungroup()
#
#       protocol_sheet <- bind_rows(protocol_sheet, new_protocols)
#     }
#
#     # --------------------------------------------------------------
#     # STEP 2: Handle new versions within existing protocol_IDs
#     # --------------------------------------------------------------
#
#     # attach protocol_ID to incoming combos
#     new_full_combos_with_id <- new_full_combos %>%
#       left_join(
#         protocol_sheet %>%
#           select(all_of(protocol_columns), protocol_ID) %>%
#           distinct(),
#         by = protocol_columns
#       )
#
#     # find unseen full combinations
#     unseen_versions <- anti_join(
#       new_full_combos_with_id,
#       protocol_sheet %>%
#         select(all_of(c(protocol_columns,
#                         version_columns,
#                         "protocol_ID"))),
#       by = c(protocol_columns, version_columns, "protocol_ID")
#     )
#
#     if (nrow(unseen_versions) > 0) {
#
#       unseen_versions <- unseen_versions %>%
#         group_by(protocol_ID) %>%
#         mutate(
#           protocol_version =
#             row_number() +
#             max(protocol_sheet$protocol_version[
#               protocol_sheet$protocol_ID == first(protocol_ID)
#             ])
#         ) %>%
#         ungroup()
#
#       protocol_sheet <- bind_rows(protocol_sheet, unseen_versions)
#     }
#   }
#
#   # ------------------------------------------------------------------
#   # FINAL: Assign IDs + versions back to df
#   # ------------------------------------------------------------------
#
#   df_with_ids <- df %>%
#     left_join(
#       protocol_sheet,
#       by = c(protocol_columns, version_columns)
#     )
#
#   return(list(
#     data = df_with_ids,
#     protocol_sheet = protocol_sheet
#   ))
# }
#
# protocol_columns <- c(
#   'nucl_acid_ext_kit',
#   'platform',
#   'instrument',
#   'seq_kit',
#   'otu_db',
#   'tax_assign_cat',
#   'otu_seq_comp_appr',
#   'minimumDepthInMeters',
#   'maximumDepthInMeters'
# )
#
# version_columns <- c(
#   'samp_size',
#   'size_frac',
#   'filter_material',
#   'samp_mat_process',
#   'samp_store_temp',
#   'samp_store_sol'
# )
#
# result <- assign_protocol_ID(
#   df = raw_protocol_test_sheet,
#   protocol_columns = protocol_columns,
#   version_columns = version_columns
# )
#
# protocol_dummy_data <- result$data
# protocol_test_sheet <- result$protocol_sheet
