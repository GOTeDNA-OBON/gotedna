

# Full file paths
existing_files <- list.files(
  "inst/app/data/raw_OBIS",
  pattern = "^dataset-.*\\.rds$",
  full.names = TRUE
)

# Extract dataset ID from filename: "dataset-XXX.rds" → "XXX"
dataset_ids <- sub("^dataset-(.*)\\.rds$", "\\1", basename(existing_files))

# Read files into a named list
datasets <- lapply(existing_files, readRDS)

# Name the list by the extracted IDs
names(datasets) <- dataset_ids

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





processed_datasets <- lapply(names(datasets), function(ds) {
  core_and_extensions <- datasets[[ds]]
  core_and_extensions <- calculate_and_enforce_columns(core_and_extensions, ds)
  core_and_extensions
})

clean_datasets <- processed_datasets[!sapply(processed_datasets, is.null)]

combined_data <- dplyr::bind_rows(clean_datasets)


################################################
#ADD STATIONS TO THE ONE BIG DATASET
################################################

combined_data <- update_location_clusters(combined_data)

################################################
#ADD PROTOCOL_ID
################################################

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


combined_data <- add_quantitative_bins_for_protocol_cols(combined_data)
protocol_result <- assign_protocol_ID(combined_data, protocol_columns)
combined_data <- protocol_result$data

saveRDS(protocol_result$protocol_sheet, 'inst/app/data/protocol_sheet.rds')


################################################
#DETECTION CALCULATION AND FILTERS
################################################

combined_data <- add_detected_column(combined_data)

################################################
#REMOVE UNNEEDED COLUMNS
################################################

added_columns <- c(
  "protocol_ID",
  "detected",
  "datasetID_obis",
  "year",
  "month",
  "stationLabel",
  "ownerContact",
  "eventDate",
  "bibliographicCitation",
  "ownerContact",
  "max_depth_bin",
  "min_depth_bin",
  "samp_size_bin",
  "station",
  "stationLabel",
  "primer"
)

all_cols <- c(required_cols, optional_columns, protocol_columns, added_columns)
final_cols <- intersect(all_cols, names(combined_data))

D_mb_clean <- combined_data[, final_cols, drop = FALSE]


################################################################################################
#STORE IN THE METABARCODING PART OF GOTEDNA_DATA AND RECORD TIMESTAMP OF OBIS PULL
################################################################################################

gotedna_data <- gotedna_data0 <- readRDS("./inst/app/data/gotedna_data.rds")

gotedna_data$metabarcoding <- D_mb_clean

writeLines(
  format(round(Sys.time(), "mins"), "%Y-%m-%d %H:%M %Z"),
  "inst/app/data/last_obis_download_ts.txt"
)



#####################################
#FILTER OUT SPECIES THAT WERE NOT DETECTED WITHIN X METRES
#####################################

#first filter duplicates and undetected
gotedna_data$metabarcoding <- remove_duplicates_and_undetected(gotedna_data$metabarcoding)

#Then use distance function to remove combinations that were not detected within x metres
print(paste0("rows in data before nondetection distance filter: ", nrow(gotedna_data$metabarcoding)))
#removing rows for species that have never been detected at that station
gotedna_data$metabarcoding <- filter_nondetections_all(
  gotedna_data$metabarcoding,
  distance = 500
)
print(paste0("rows in data after nondetection distance filter: ", nrow(gotedna_data$metabarcoding)))



################################################
#ADD STATION FILE FOR MAPPING
################################################


gotedna_station <- list(
  metabarcoding = gotedna_data$metabarcoding |> get_station(),
  qPCR = gotedna_data$qPCR |> get_station()
)
saveRDS(gotedna_station, "inst/app/data/gotedna_station.rds")


################################################
#CHECK PRIMERS FOR TYPOS AND MAKE/STORE PRIMER SHEET
################################################

library(stringdist)
library(dplyr)
library(stringr)

valid_primer_df <- read.csv("inst/app/data/primers.csv")
valid_primer_df <- valid_primer_df %>% filter(Data == "Metabarcoding")

valid_primer_df <- valid_primer_df %>%
  mutate(
    # split OriginalName on <br>
    primers = str_split(OriginalName, "\\s*<br>\\s*"),

    forward_primer = purrr::map_chr(primers, 1),
    reverse_primer = purrr::map_chr(primers, 2),

    primer_display_name = paste0(
      Locus, " | ",
      forward_primer, " / ", reverse_primer
    )
  ) %>%
  select(-primers)

valid_primers <- valid_primer_df$primer_display_name

map_primers <- function(
    observed,
    valid,
    method = "osa",
    max_dist = 2
) {

  # distance matrix: rows = observed, cols = valid
  dist_mat <- stringdistmatrix(
    observed,
    valid,
    method = method
  )

  best_match_idx  <- apply(dist_mat, 1, which.min)
  best_match_dist <- apply(dist_mat, 1, min)

  tibble(
    observed_primer = observed,
    matched_primer  = valid[best_match_idx],
    distance        = best_match_dist,
    status = case_when(
      best_match_dist == 0              ~ "exact",
      best_match_dist <= max_dist       ~ "corrected",
      TRUE                              ~ "new_or_unmatched"
    ),
    final_primer = if_else(
      best_match_dist <= max_dist,
      valid[best_match_idx],
      observed
    )
  )
}

primer_map <- map_primers(
  observed = unique(gotedna_data$metabarcoding$primer),
  valid    = valid_primers,
  max_dist = 10
)

View(primer_map)

gotedna_data$metabarcoding <- gotedna_data$metabarcoding %>%
  left_join(
    primer_map %>% select(observed_primer, final_primer),
    by = c("primer" = "observed_primer")
  ) %>%
  mutate(primer = coalesce(final_primer, primer)) %>%
  select(-final_primer)


#### Need to do for both qPCR and metabarcoding
# Prepare primer data
# gotedna_data <- readRDS("inst/app/data/gotedna_data.rds")

saveRDS(gotedna_data, "inst/app/data/gotedna_data.rds")

taxon_levels <- c(
  "kingdom", "phylum", "class",
  "order", "family", "genus", "scientificName"
)

gotedna_primer <- setNames(
  lapply(taxon_levels, function(level) {
    make_primer_sheet(gotedna_data$metabarcoding, level)
  }),
  taxon_levels
)



saveRDS(gotedna_primer, "inst/app/data/gotedna_primer.rds")

