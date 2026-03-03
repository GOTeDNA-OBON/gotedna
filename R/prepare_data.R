
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

calculate_and_enforce_columns <- function(core_and_extensions, ds = NULL) {

  # Skip if dataset is NULL
  if (is.null(core_and_extensions)) return(NULL)

  core_and_extensions <- core_and_extensions %>%
    filter(!is.na(samp_name))

  if (!is.null(ds)) {
    core_and_extensions <- core_and_extensions %>%
      mutate(
        datasetID_obis   = ds
      )
  }

  # --- Core transformations ---
  core_and_extensions <- core_and_extensions %>%
    mutate(
      eventDate_chr    = as.character(eventDate),
      eventDate_clean  = suppressWarnings(lubridate::ymd(substr(eventDate_chr, 1, 10))),
      year             = lubridate::year(eventDate_clean),
      month            = lubridate::month(eventDate_clean),
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

  # --- Enforce schema ---
  core_and_extensions <- enforce_schema(
    core_and_extensions,
    required = required_cols,
    optional = optional_columns
  )

  core_and_extensions
}

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

################################################
#ADD PROTOCOL_ID
################################################

add_quantitative_bins_for_protocol_cols <- function(df) {

  df <- df %>%
    mutate(
      # ---- Coerce samp_size safely ----
      samp_size_num = suppressWarnings(
        as.numeric(gsub("[^0-9.]", "", samp_size))
      ),
      samp_size_num = dplyr::if_else(
        samp_size_unit == "mL",
        samp_size_num / 1000,
        samp_size_num
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

      samp_size_mid = if_else(
        is.na(samp_size_floor),
        NA_real_,
        round(samp_size_floor + 0.125, 3)
      ),

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


################################################
#DETECTION CALCULATION AND FILTERS
################################################

add_detected_column <- function(df) {
  #THIS REMOVES ALL 1-10 READ DATA
  df <- df %>%
    dplyr::mutate(msct = case_when(
      organismQuantity == 0 ~ TRUE,
      organismQuantity > 10 ~ TRUE
    )) |>
    tidyr::drop_na(msct)

  #NOW DETECTED IS ADDED
  df <- df %>%
    mutate(
      detected = dplyr::case_when(
        !is.na(organismQuantity) & organismQuantity > 0 ~ 1L,
        TRUE ~ 0L
      )
    )
  df
}


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


filter_nondetections_specifics <- function(df_subset, distance = 500) {
  sf::sf_use_s2(TRUE)

  if (nrow(df_subset) == 0) return(df_subset)

  # 1) Coordinate-level detection flag (keep ALL rows after join)
  coord_detection_summary <- df_subset %>%
    group_by(decimalLongitude, decimalLatitude) %>%
    summarise(
      coord_has_detection = any(detected == 1),
      .groups = "drop"
    )

  df_subset <- df_subset %>%
    left_join(coord_detection_summary,
              by = c("decimalLongitude", "decimalLatitude"))

  # 2) Unique coordinates with coord_id
  coords <- df_subset %>%
    distinct(decimalLongitude, decimalLatitude, coord_has_detection) %>%
    mutate(coord_id = row_number())

  df_subset <- df_subset %>%
    left_join(coords, by = c("decimalLongitude", "decimalLatitude", "coord_has_detection"))

  # If only one coordinate, the only way to be "within distance of a detection"
  # is if that coordinate itself has a detection.
  if (nrow(coords) == 1) {
    return(df_subset %>% filter(coord_has_detection))
  }

  # 3) sf points in lon/lat (global-safe)
  coords_sf <- sf::st_as_sf(
    coords,
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326
  )

  # 4) Neighbor list within distance (meters)
  nbrs <- sf::st_is_within_distance(coords_sf, coords_sf, dist = distance)

  # 5) Propagate coord_has_detection through neighbor lists
  coord_detected_within_dist <- vapply(
    seq_along(nbrs),
    function(i) any(coords$coord_has_detection[nbrs[[i]]]),
    logical(1)
  )

  # 6) Attach to all rows + final filter rule:
  # keep rows if:
  #   - taxon was ever detected at that coordinate (coord_has_detection)
  #   OR
  #   - taxon detected within distance of that coordinate (detected_within_dist)
  df_subset %>%
    mutate(detected_within_dist = coord_detected_within_dist[coord_id]) %>%
    filter(coord_has_detection | detected_within_dist)
}


filter_nondetections_all <- function(df,
                                     distance = 500,
                                     selected_taxon_level = "scientificName",
                                     selected_taxon_id = "All") {

  stopifnot(selected_taxon_level %in% names(df))
  stopifnot(all(c("protocol_ID", "decimalLongitude", "decimalLatitude", "detected") %in% names(df)))

  # # Optional filter to one taxon
  # if (!identical(selected_taxon_id, "All")) {
  #   df <- df %>% filter(.data[[selected_taxon_level]] == selected_taxon_id)
  # }
  print(nrow(df))
  if (nrow(df) == 0) return(df)

  df %>%
    group_by(protocol_ID, primer, .data[[selected_taxon_level]]) %>%
    group_modify(~ filter_nondetections_specifics(.x, distance = distance)) %>%
    ungroup()
}
