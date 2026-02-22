#' *New proposal of this page*
#' Read and format metabarcoding metadata and data from OBIS
#'
#' @description Function that reads for any data in OBIS that has a DNA Derived Data
#' extension. Will compile occurrence (core) and DNA Derived Data (extension) files
#' into a single dataframe. Data frames are then formatted appropriately for use in
#' subsequence analysis and visualizations (e.g., date reformatting, merging metadata
#' and data).
#' * Note that template column names are in the format of Darwin Core Archive
#' (DwC-A) using Darwin Core (DwC) data standards where possible.
#'
#' @return A tibble with 25 columns:
#' * `protocol_ID`
#' * `protocolVersion`
#' * `samp_name`
#' * `eventID`
#' * `primer`
#' * `species`
#' * `domain`
#' * `kingdom`
#' * `phylum`
#' * `class`
#' * `order`
#' * `family`
#' * `genus`
#' * `concentration`: provided when choose.method = "qPCR"
#' * `pcr_primer_lod` : provided when choose.method = "qPCR"
#' * `organismQuantity`: provided when choose.method = "metabarcoding"
#' * `date`
#' * ecodistrict`
#' * `LClabel` : Local Contexts label to denote First Nations data sovereignty
#' * `decimalLatitude`
#' * `decimalLongitude`
#' * `station`
#' * `year`
#' * `month`
#' * `detected`
#' * `msct` :logical, where minimum sequence copy threshold = 10
#' * `ownerContact` : email of data owner/steward
#' * `bibliographicCitation` : DOI reference, if applicable
#'
#' @author Anais Lacoursiere-Roussel \email{Anais.Lacoursiere@@dfo-mpo.gc.ca}
#' @rdname read_data
#' @export
#' @examples
#' \dontrun{
#' D_mb <- read_data(
#'  choose.method = "metabarcoding", path.folder = "./inst/testdata"
#' )
#' }


read_raw_data <- function(
    dataset_ids    = NULL,
    scientificname = NULL,
    worms_id       = NULL,
    areaid         = NULL,
    join_by        = c("auto", "occurrenceID", "id"),
    require_absences = TRUE
) {
  library(robis)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(purrr)
  library(dbscan)
  library(geosphere)
  library(sf)

  occurrence_cols <- c(
    "recordedBy",
    "bibliographicCitation",
    "materialSampleID",
    "organismQuantity",
    "organismQuantityType",
    "sampleSizeValue",
    "sampleSizeUnit",
    "associatedSequences",
    "minimumDepthInMeters",
    "maximumDepthInMeters",
    "month",
    "year",
    "scientificNameID",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus"
  )

  dna_cols <- c(
    "id",
    "dna_sequence",
    "target_gene",
    "pcr_primer_forward",
    "pcr_primer_forward",   # appears twice in your list
    "samp_name",
    "env_broad_scale",
    "env_local_scale",
    "env_medium",
    "samp_mat_process",
    "size_frac",
    "samp_size",
    "samp_size_unit",
    "otu_db",
    "seq_kit",
    "otu_seq_comp_appr",
    "pcr_primer_name_forward",
    "pcr_primer_name_reverse",
    "pcr_primer_reference",
    "occurrenceID"
  )


  mof_cols <- c(
    "id",
    "seq_id",
    "samp_category",
    "checkls_ver",
    "assay_name",
    "assay_type",
    "targetTaxonomicAssay",
    "geo_loc_name",
    "technical_rep_id",
    "project_contact",
    "seq_run_id",
    "lib_id",
    "project_id",
    "pcr_0_1",
    "samp_store_sol",
    "samp_store_temp",
    "platform",
    "instrument",
    "tax_assign_cat",
    "LClabel",
    "occurrenceID",
    "nucl_acid_ext",
    "nucl_acid_ext_kit",
    "filter_material"
  )

  added_cols <- c("category", "flags")

  mandatory_obis <- c(
    "occurrenceID",
    "eventDate",
    "decimalLongitude",
    "decimalLatitude",
    "scientificName",
    "occurrenceStatus",
    "basisOfRecord"
  )

  cols_included_from_OBIS <- unique(c(occurrence_cols, dna_cols, mof_cols, added_cols, mandatory_obis))


  join_by <- match.arg(join_by)
  ## 0. If dataset_ids is NULL, discover them via robis::dataset() ----
  if (is.null(dataset_ids)) {
    message("Discovering datasets with DNADerivedData and MeasurementOrFact extensions ...")
    ds_tbl <- robis::dataset(
      scientificname = scientificname,
      areaid         = areaid,
      taxonid        = worms_id,
      hasextensions  = c("DNADerivedData", "MeasurementOrFact")
    )

    ds_tbl <- ds_tbl |>
      dplyr::filter(statistics$absence != 0)

    if (nrow(ds_tbl) == 0L) {
      warning("No datasets found that match scientificname/areaid AND have DNADerivedData.")
      return(tibble::tibble())
    }

    # dataset() typically returns a column called 'id' for dataset id
    if ("id" %in% names(ds_tbl)) {
      dataset_ids <- unique(ds_tbl$id)
    } else if ("datasetid" %in% names(ds_tbl)) {
      dataset_ids <- unique(ds_tbl$datasetid)
    } else {
      stop("Could not find dataset id column in dataset() output.")
    }
    message("Found ", length(dataset_ids), " dataset(s) with DNADerivedData.")
    dataset_ids <- as.character(dataset_ids)
    existing_files <- list.files("inst/app/data/raw_OBIS", pattern = "^dataset-.*\\.rds$")

    saved_ds <- sub("^dataset-(.*)\\.rds$", "\\1", existing_files)

    dataset_ids <- setdiff(dataset_ids, saved_ds)
  }
  dataset_ids <- as.character(dataset_ids)
  print("About to start pulling these datasets: ")
  print(dataset_ids)
  ## 1. Loop over datasets and pull DNADerivedData occurrences ----
  for (ds in dataset_ids) {
    Sys.sleep(3)
    message("Pulling OBIS dataset: ", ds)
    # ---- 1a. Check dataset has DNADerivedData extension (defensive) ----
    ds_meta <- robis::dataset(datasetid = ds)
    exts <- tolower(unlist(ds_meta$extensions))
    if (!"dnaderiveddata" %in% exts) {
      warning("Dataset ", ds, " has no DNADerivedData extension; skipping.")
      return(NULL)
    }

    if (!"measurementorfact" %in% exts) {
      warning("Dataset ", ds, " has no MeasurementOrFact extension; skipping.")
      return(NULL)
    }
    exclude_list <- c("NO_COORD",
      "ZERO_COORD",
      "LON_OUT_OF_RANGE",
      "LAT_OUT_OF_RANGE",
      "NO_MATCH"
    )
    # ---- 1b. Pull occurrence records with DNADerivedData + filters ----
    rec <- tryCatch(
      {
        robis::occurrence(
          datasetid      = ds,
          scientificname = scientificname,
          taxonid        = worms_id,
          areaid         = areaid,
          absence        = "include",
          extensions     = c("DNADerivedData", "MeasurementOrFact"),
          hasextensions  = c("DNADerivedData", "MeasurementOrFact"),
          dropped = "include",
          exclude = exclude_list
        )
      },
      error = function(e) {
        warning("Failed to fetch dataset ", ds, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    # skip iteration if rec is NULL
    if (is.null(rec)) {
      return(NULL)
    }

    if (nrow(rec) == 0L) {
      warning("No occurrence records returned for dataset ", ds, " with these filters.")
      return(NULL)
    }



    # ---- 1c. Build core_occ and filter on occurrenceStatus ----
    core_occ <- rec %>%
      distinct(occurrenceID, .keep_all = TRUE)
    if (!"occurrenceStatus" %in% names(core_occ)) {
      warning("Dataset ", ds, " has no occurrenceStatus column; skipping.")
      return(NULL)
    }
    # Unique non-NA status values
    status_vals <- unique(na.omit(core_occ$occurrenceStatus))
    # Require BOTH "present" and "absent"
    if (require_absences) {
      if (!all(c("present", "absent") %in% status_vals)) {
        warning(
          "Dataset ", ds,
          " does not contain both 'present' and 'absent' in occurrenceStatus; skipping."
        )
        return(NULL)
      }
    }
    # Keep also an id-based core for joining if needed
    core_id <- rec %>%
      distinct(id, .keep_all = TRUE)
    # DNADerivedData extension (includes `id` by default)
    dna_only <- robis::unnest_extension(rec, "DNADerivedData")

    shared_dna_cols <- intersect(cols_included_from_OBIS, names(dna_only))
    print("dna shared cols: ")
    print(shared_dna_cols)
    dna_only <- dna_only %>% select(shared_dna_cols)
    #MeasurementOfFact extension
    mof_only <- unnest_extension(rec, "MeasurementOrFact")

    mof_only <- mof_only %>%
      group_by(occurrenceID, measurementType) %>%
      slice(1) %>%
      ungroup(1)
    # ---- 2. Decide how to join core + extension ----
    wide_mof <- mof_only %>%
      pivot_wider(
        id_cols = c(occurrenceID, id),
        names_from = measurementType,
        values_from = measurementValue
      )

    shared_mof_cols <- intersect(cols_included_from_OBIS, names(mof_only))
    print("mof shared cols: ")
    print(shared_mof_cols)
    mof_only <- mof_only %>% select(shared_mof_cols)
    mof_and_dna <- wide_mof %>% left_join(dna_only, by = "id")

    shared_cols <- intersect(cols_included_from_OBIS, names(rec))
    print("rec shared cols: ")
    print(shared_cols)
    rec <- rec %>% select(shared_cols)
    join_choice <- join_by
    if (join_choice == "auto") {
      # avoid vector-recycling warning by checking non-NA separately
      can_occ <- "occurrenceID" %in% names(core_occ) &&
        "occurrenceID" %in% names(mof_and_dna) &&
        any(!is.na(core_occ$occurrenceID)) &&
        any(!is.na(mof_and_dna$occurrenceID))
      can_id  <- "id" %in% names(core_id) &&
        "id" %in% names(mof_and_dna) &&
        any(!is.na(core_id$id)) &&
        any(!is.na(mof_and_dna$id))
      if (can_occ) {
        join_choice <- "occurrenceID"
      } else if (can_id) {
        join_choice <- "id"
      } else {
        stop("Neither occurrenceID nor id can be used to join core and DNADerivedData for dataset ",
             ds, ".")
      }
    }
    if (join_choice == "occurrenceID") {
      core_and_extensions <- core_occ %>% left_join(mof_and_dna, by = "occurrenceID")
    } else {  # "id"
      core_and_extensions <- core_id %>%
        left_join(mof_and_dna, by = "id")
    }
    # ---- 3. Basic cleaning & derived fields ----

    if (!has_required_cols(core_and_extensions, required_cols)) {
      return(NULL)
    }

    if (is.null(core_and_extensions)) {
      return(NULL)
    }
    dup_names <- names(core_and_extensions)[duplicated(names(core_and_extensions))]
    all_shared_cols <- intersect(core_and_extensions, cols_included_from_OBIS)
    core_and_extensions <- select(all_shared_cols)
    if (length(dup_names) > 0) {
      message("Duplicate column names detected:")
      print(unique(dup_names))
    }
    # ---- 5. Return in GOTeDNA_df-like shape ----
    out <- core_and_extensions %>%
      mutate(
        samp_name = as.character(samp_name),
      )
    saveRDS(out, paste0("inst/app/data/raw_OBIS/dataset-", ds, ".rds"))

  }
}


###########################
#Code below here is for adding locations to the data based on efficient clustering
###########################

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
  "eventDate",
  "decimalLatitude",
  "decimalLongitude",
  "organismQuantity",
  "concentration",
  "pcr_primer_lod"
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


has_required_cols <- function(df, required_cols) {
  missing <- setdiff(required_cols, colnames(df))

  if (length(missing) > 0) {
    message(
      "Skipping dataset — missing columns: ",
      paste(missing, collapse = ", ")
    )
    return(FALSE)
  }

  TRUE
}

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
  df <- df[, c(required, optional), drop = FALSE]

  df
}


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

#columns from the flow chart
#Occurrence


