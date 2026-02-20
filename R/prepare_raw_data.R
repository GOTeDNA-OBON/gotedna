

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
  "organismQuantity",
  "concentration",
  "pcr_primer_lod",
  "datasetID_obis"
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

enforce_schema <- function(df, required, optional) {
  # Check required columns exist
  missing_required <- setdiff(required, names(df))
  if (length(missing_required)) {
    message("Missing required columns: ", paste(missing_required, collapse = ", "))
    return(NULL)
  }

  # Only add optional columns that aren't already in df
  optional_to_add <- setdiff(optional, names(df))
  if (length(optional_to_add)) {
    df[optional_to_add] <- NA
  }

  # Strict allow-list: only include each column once
  # df <- df[, c(required, optional), drop = FALSE]

  df
}


processed_datasets <- lapply(names(datasets), function(ds) {

  core_and_extensions <- datasets[[ds]]

  # Skip if dataset is NULL
  if (is.null(core_and_extensions)) return(NULL)

  # --- Core transformations ---
  core_and_extensions <- core_and_extensions %>%
    filter(!is.na(samp_name)) %>%
    mutate(
      datasetID_obis   = ds,  # now comes from filename
      eventDate_chr    = as.character(eventDate),
      eventDate_clean  = suppressWarnings(ymd(substr(eventDate_chr, 1, 10))),
      year             = year(eventDate_clean),
      month            = month(eventDate_clean),
      decimalLatitude  = suppressWarnings(as.numeric(decimalLatitude)),
      decimalLongitude = suppressWarnings(as.numeric(decimalLongitude))
    ) %>%
    mutate(
      station               = if ("samplingStation" %in% names(.)) samplingStation else NA_character_,
      ownerContact          = project_contact,
      bibliographicCitation = if ("bibliographicCitation" %in% names(.)) bibliographicCitation else NA_character_,
      organismQuantity      = suppressWarnings(as.numeric(organismQuantity))
    ) %>%
    mutate(
      primer     = sprintf(
        "%s | %s / %s",
        target_gene,
        pcr_primer_name_forward,
        pcr_primer_name_reverse
      ),
      eventDate  = eventDate_clean
    )

  # --- Enforce schema ---
  core_and_extensions <- enforce_schema(
    core_and_extensions,
    required = required_cols,
    optional = optional_columns
  )

  core_and_extensions
})





clean_datasets <- processed_datasets[!sapply(processed_datasets, is.null)]

combined_data <- dplyr::bind_rows(clean_datasets)


################################################
#ADD STATIONS TO THE ONE BIG DATASET
################################################


# ----------------------------
# Complete-linkage helper
# ----------------------------
run_complete_linkage <- function(df, threshold_m, cluster_prefix = "C") {
  if (nrow(df) == 1) return(tibble(cluster = paste0(cluster_prefix, "_1")))

  coords <- df[, c("lon", "lat")]
  dist_mat <- geosphere::distm(coords)

  hc <- hclust(as.dist(dist_mat), method = "complete")
  clusters <- cutree(hc, h = threshold_m)

  tibble(cluster = paste0(cluster_prefix, "_", clusters))
}

# ----------------------------
# Recursive hybrid clustering
# ----------------------------
hybrid_cluster_unique <- function(df_unique, threshold_m, max_hc_size, cluster_prefix = "S") {

  if (nrow(df_unique) <= max_hc_size) {
    return(run_complete_linkage(df_unique, threshold_m, cluster_prefix))
  }

  # Compute approximate max distance for coarse DBSCAN
  coords <- df_unique[, c("lon", "lat")]
  dist_mat <- geosphere::distm(coords)
  max_dist <- max(dist_mat)
  eps <- max_dist / 2

  db <- dbscan(coords, eps = eps, minPts = 1)
  df_unique$coarse_cluster <- db$cluster

  result_clusters <- df_unique %>%
    group_by(coarse_cluster) %>%
    group_modify(~{
      if (nrow(.x) > max_hc_size) {
        hybrid_cluster_unique(.x %>% select(-coarse_cluster), threshold_m, max_hc_size,
                              cluster_prefix = paste0(cluster_prefix, "_", unique(.x$coarse_cluster)))
      } else {
        run_complete_linkage(.x %>% select(-coarse_cluster), threshold_m,
                             cluster_prefix = paste0(cluster_prefix, "_", unique(.x$coarse_cluster)))
      }
    }) %>%
    ungroup()

  return(result_clusters)
}

# ----------------------------
# Main function
# ----------------------------
update_location_clusters <- function(df,
                                     lat_col = "decimalLatitude",
                                     long_col = "decimalLongitude",
                                     distance_threshold = 50,
                                     max_hc_size = 500) {

  # ---- 1. Ensure stationLabel exists ----
  if (!"stationLabel" %in% names(df)) {
    if ("station" %in% names(df)) {
      df$stationLabel <- df$station
    } else {
      df$stationLabel <- NA_character_
    }
  }

  # ---- 2. Keep only valid coordinates ----
  valid_df <- df %>%
    filter(!is.na(.data[[lat_col]]), !is.na(.data[[long_col]]))

  if (nrow(valid_df) == 0) {
    df$station <- NA_character_
    return(df)
  }

  # ---- 3. Reduce to unique coordinates ----
  coords_unique <- valid_df %>%
    distinct(.data[[lat_col]], .data[[long_col]]) %>%
    rename(lat = all_of(lat_col), lon = all_of(long_col))

  # ---- 4. Run hybrid clustering on unique coordinates ----
  station_map <- hybrid_cluster_unique(coords_unique,
                                       threshold_m = distance_threshold,
                                       max_hc_size = max_hc_size,
                                       cluster_prefix = "Location")

  # ---- 5. Map station IDs back to full dataset ----
  coords_unique$station <- as.integer(factor(station_map$cluster))

  df$station <- NULL
  df <- df %>%
    left_join(coords_unique, by = setNames(c("lat", "lon"), c(lat_col, long_col)))

  df
}



# df: dataset with lon/lat and station column
# lon_col / lat_col: coordinate columns
# station_col: station IDs
# threshold_m: maximum allowed distance
check_station_distances_unique <- function(df,
                                           lon_col = "decimalLongitude",
                                           lat_col = "decimalLatitude",
                                           station_col = "station",
                                           threshold_m = 500) {

  violations <- df %>%
    group_by(.data[[station_col]]) %>%
    group_modify(~{
      # ---- 1. Reduce to unique coordinates for efficiency ----
      pts <- .x %>%
        distinct(.data[[lon_col]], .data[[lat_col]]) %>%
        select(all_of(c(lon_col, lat_col)))

      n <- nrow(pts)
      if (n <= 1) return(tibble())

      # ---- 2. Compute pairwise distance matrix ----
      dmat <- geosphere::distm(pts)

      # ---- 3. Check for distances exceeding threshold ----
      idx <- which(dmat > threshold_m, arr.ind = TRUE)
      idx <- idx[idx[,1] < idx[,2], , drop = FALSE]  # upper triangle only
      if (nrow(idx) == 0) return(tibble())

      tibble(
        station = unique(.x[[station_col]]),
        point1 = idx[,1],
        point2 = idx[,2],
        distance_m = dmat[idx]
      )
    }) %>%
    ungroup()

  if (nrow(violations) == 0) {
    message("✅ All stations pass: no pair of points exceeds ", threshold_m, " meters.")
  } else {
    warning("⚠️ Found ", nrow(violations), " pairs exceeding ", threshold_m, " meters.")
  }

  return(violations)
}

combined_data <- update_location_clusters(combined_data)


################################################
#ADD PROTOCOL_ID AND PROTOCOLVERSION
################################################


add_quantitative_bins_for_protocol_cols <- function(df) {

  df <- df %>%
    mutate(
      # ---- Coerce samp_size safely ----
      samp_size_num = suppressWarnings(
        as.numeric(gsub("[^0-9.]", "", samp_size))
      ),

      # ---- Depth bins ----
      min_depth_floor = floor(round(minimumDepthInMeters, 6) / 5) * 5,
      max_depth_floor = floor(round(maximumDepthInMeters, 6) / 5) * 5,

      min_depth_bin = paste0(min_depth_floor, "-", min_depth_floor + 5, "m"),
      max_depth_bin = paste0(max_depth_floor, "-", max_depth_floor + 5, "m"),

      # ---- Sample size bins ----
      samp_size_floor = if_else(
        is.na(samp_size_num),
        NA_real_,
        pmax(
          0,
          floor((samp_size_num - 0.125) / 0.25) * 0.25 + 0.125
        )
      ),

      samp_size_upper = samp_size_floor + 0.25,

      samp_size_bin = if_else(
        is.na(samp_size_floor),
        NA_character_,
        sprintf(
          "%.3f-%.3fL",
          round(samp_size_floor, 3),
          round(samp_size_upper, 3)
        )
      )
    ) %>%
    select(-samp_size_upper)

  df
}


assign_protocol_ID <- function(df,
                               protocol_columns,
                               version_columns,
                               protocol_sheet = NULL) {

  # Remove existing IDs if present
  df <- df %>%
    select(-any_of(c("protocol_ID", "protocolVersion")))

  # Distinct combinations at full granularity
  new_full_combos <- df %>%
    select(all_of(c(protocol_columns, version_columns))) %>%
    distinct()

  # ------------------------------------------------------------------
  # CASE 1: No existing protocol_sheet → build from scratch
  # ------------------------------------------------------------------
  if (is.null(protocol_sheet) || nrow(protocol_sheet) == 0) {

    protocol_sheet <- new_full_combos %>%
      group_by(across(all_of(protocol_columns))) %>%
      mutate(
        protocol_ID = cur_group_id(),
        protocolVersion = row_number()
      ) %>%
      ungroup()

  } else {

    # Ensure sheet has required structure
    required_cols <- c(protocol_columns, version_columns,
                       "protocol_ID", "protocolVersion")

    missing_cols <- setdiff(required_cols, names(protocol_sheet))
    if (length(missing_cols) > 0) {
      stop("protocol_sheet is missing required columns: ",
           paste(missing_cols, collapse = ", "))
    }

    # --------------------------------------------------------------
    # STEP 1: Add new protocol_IDs if needed
    # --------------------------------------------------------------

    new_protocols <- anti_join(
      new_full_combos %>% select(all_of(protocol_columns)) %>% distinct(),
      protocol_sheet %>% select(all_of(protocol_columns)) %>% distinct(),
      by = protocol_columns
    )

    if (nrow(new_protocols) > 0) {

      max_id <- max(protocol_sheet$protocol_ID, na.rm = TRUE)

      new_protocols <- new_protocols %>%
        mutate(protocol_ID = row_number() + max_id)

      # give them version 1 initially (will expand below if needed)
      new_protocols <- new_protocols %>%
        left_join(new_full_combos, by = protocol_columns) %>%
        group_by(protocol_ID) %>%
        mutate(protocolVersion = row_number()) %>%
        ungroup()

      protocol_sheet <- bind_rows(protocol_sheet, new_protocols)
    }

    # --------------------------------------------------------------
    # STEP 2: Handle new versions within existing protocol_IDs
    # --------------------------------------------------------------

    # attach protocol_ID to incoming combos
    new_full_combos_with_id <- new_full_combos %>%
      left_join(
        protocol_sheet %>%
          select(all_of(protocol_columns), protocol_ID) %>%
          distinct(),
        by = protocol_columns
      )

    # find unseen full combinations
    unseen_versions <- anti_join(
      new_full_combos_with_id,
      protocol_sheet %>%
        select(all_of(c(protocol_columns,
                        version_columns,
                        "protocol_ID"))),
      by = c(protocol_columns, version_columns, "protocol_ID")
    )

    if (nrow(unseen_versions) > 0) {

      unseen_versions <- unseen_versions %>%
        group_by(protocol_ID) %>%
        mutate(
          protocolVersion =
            row_number() +
            max(protocol_sheet$protocolVersion[
              protocol_sheet$protocol_ID == first(protocol_ID)
            ])
        ) %>%
        ungroup()

      protocol_sheet <- bind_rows(protocol_sheet, unseen_versions)
    }
  }

  # ------------------------------------------------------------------
  # FINAL: Assign IDs + versions back to df
  # ------------------------------------------------------------------

  df_with_ids <- df %>%
    left_join(
      protocol_sheet,
      by = c(protocol_columns, version_columns)
    )

  return(list(
    data = df_with_ids,
    protocol_sheet = protocol_sheet
  ))
}


combined_data <- add_quantitative_bins_for_protocol_cols(combined_data)
protocol_result <- assign_protocol_ID(combined_data)
combined_data <- protocol_result$data

saveRDS(protocol_result$protocol_sheet, 'inst/app/data/protocol_sheet.rds')


################################################
#DETECTION CALCULATION AND FILTERS
################################################

#THIS REMOVES ALL 1-10 READ DATA
combined_data <- combined_data %>%
  dplyr::mutate(msct = case_when(
    organismQuantity == 0 ~ TRUE,
    organismQuantity > 10 ~ TRUE
  )) |>
  tidyr::drop_na(msct)

#NOW DETECTED IS ADDED
combined_data <- combined_data %>%
  mutate(
    detected = dplyr::case_when(
      !is.na(organismQuantity) & organismQuantity > 0 ~ 1L,
      TRUE ~ 0L
    )
  )

#THIS IS REMOVING ROWS, FOR OLD FUNCTIONALITY, REMOVE!
D_mb_nodetect <- combined_data %>%
  dplyr::group_by(
    protocol_ID, protocolVersion, scientificName, primer, station) %>%
  dplyr::summarise(num_detected = sum(detected)) %>%
  dplyr::filter(num_detected == 0)

D_mb_clean <- dplyr::anti_join(combined_data, D_mb_nodetect,
                               by = c("protocol_ID","protocolVersion","scientificName",
                                      "primer", "station"))


################################################
#REMOVE UNNEEDED COLUMNS
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
  'max_depth_floor'
)

version_columns <- c(
  'samp_size_floor',
  'size_frac',
  'filter_material',
  'samp_mat_process',
  'samp_store_temp',
  'samp_store_sol'
)

added_columns <- c(
  "protocol_ID",
  "protocolVersion",
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

)

all_cols <- c(required_cols, optional_columns, protocol_columns, version_columns, added_columns)
final_cols <- intersect(all_cols, names(D_mb_clean))

D_mb_clean <- D_mb_clean[, final_cols, drop = FALSE]


################################################################################################
#STORE IN THE METABARCODING PART OF GOTEDNA_DATA AND RECORD TIMESTAMP OF OBIS PULL
################################################################################################

gotedna_data <- gotedna_data0 <- readRDS("./inst/app/data/gotedna_data.rds")

gotedna_data$metabarcoding <- D_mb_clean

writeLines(
  format(round(Sys.time(), "mins"), "%Y-%m-%d %H:%M %Z"),
  "inst/app/data/last_obis_download_ts.txt"
)

saveRDS(gotedna_data, "inst/app/data/gotedna_data.rds")
################################################
#ADD STATION FILE FOR MAPPING
################################################

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

# View(primer_map)

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

newprob_mb <- calc_det_prob(gotedna_data$metabarcoding)
scaledprobs_mb <- scale_newprob(gotedna_data$metabarcoding, newprob_mb)
scaledprobs_q <- NULL

gotedna_primer <- list()

# this needs to be based on the area selection
for (i in c("kingdom", "phylum", "class", "order", "family", "genus", "scientificName")) {
  gotedna_primer[[i]] <- primer_sort(i, dplyr::bind_rows(scaledprobs_mb, scaledprobs_q)) |>
    mutate(text = paste0(primer, " (", detects, "/", total, " ", perc, "%)"))
}

saveRDS(gotedna_primer, "inst/app/data/gotedna_primer.rds")




















#
# test_bf <- readRDS("inst/app/data/raw_OBIS/dataset-e596c238-1aa1-4d56-a26c-aac622c0c246.rds")
#
# test_sm <-readRDS("inst/app/data/dataset-e596c238-1aa1-4d56-a26c-aac622c0c246.rds")
#






