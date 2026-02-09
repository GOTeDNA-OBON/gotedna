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

read_data <- function(
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
  }
  dataset_ids <- as.character(dataset_ids)
  ## 1. Loop over datasets and pull DNADerivedData occurrences ----
  obis_list <- purrr::map(dataset_ids, function(ds) {
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
    # ---- 1b. Pull occurrence records with DNADerivedData + filters ----
    rec <- robis::occurrence(
      datasetid      = ds,
      scientificname = scientificname,
      taxonid        = worms_id,
      areaid         = areaid,
      absence        = "include",
      extensions     = c("DNADerivedData", "MeasurementOrFact"),
      hasextensions  = c("DNADerivedData", "MeasurementOrFact")
    )

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

    mof_and_dna <- wide_mof %>% left_join(dna_only, by = "id")

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
    core_and_extensions <- core_and_extensions %>%
      filter(!is.na(samp_name)) %>%  # keep only real samples
      mutate(
        datasetID_obis = ds,
        # eventDate usually ISO; strip time & parse
        eventDate_chr   = as.character(eventDate),
        eventDate_clean = suppressWarnings(
          ymd(substr(eventDate_chr, 1, 10))
        ),
        year  = year(eventDate_clean),
        month = month(eventDate_clean),
        decimalLatitude  = suppressWarnings(as.numeric(decimalLatitude)),
        decimalLongitude = suppressWarnings(as.numeric(decimalLongitude))
      )
    # Safe versions if fields are missing
    core_and_extensions <- core_and_extensions %>%
      mutate(
        station = if ("samplingStation" %in% names(.)) samplingStation else NA_character_,
        ownerContact = if ("ownerInstitutionCode" %in% names(.)) ownerInstitutionCode else NA_character_,
        bibliographicCitation = if ("bibliographicCitation" %in% names(.)) bibliographicCitation else NA_character_
      )
    # ---- 4. Metabarcoding detection + primer ----
    core_and_extensions <- core_and_extensions %>%
      mutate(
        organismQuantity = suppressWarnings(as.numeric(organismQuantity)),
        detected = dplyr::case_when(
          !is.na(organismQuantity) & organismQuantity > 0 ~ 1L,
          TRUE ~ 0L
        ),
        primer = dplyr::coalesce(target_subfragment, target_gene)
      )
    if (!check_required_cols(core_and_extensions, old_cols, ds)) {
      return(NULL)
    }
    # ---- 5. Return in GOTeDNA_df-like shape ----
    out <- core_and_extensions %>%
      mutate(
        samp_name = as.character(samp_name),
        sprintf(
          "%s | %s / %s",
          target_gene,
          pcr_primer_name_forward,
          pcr_primer_name_reverse
        )
      )
    out <- core_and_extensions %>%
      transmute(
        protocol_ID           = protocol_ID,
        protocolVersion       = protocolVersion,
        samp_name             = as.character(samp_name),
        primer                = sprintf("%s | %s / %s", target_gene, pcr_primer_name_forward, pcr_primer_name_reverse),
        scientificName        = scientificName,
        kingdom               = kingdom,
        phylum                = phylum,
        class                 = class,
        order                 = order,
        family                = family,
        genus                 = genus,
        eventDate             = eventDate_clean,
        LClabel               = LClabel,
        decimalLatitude       = decimalLatitude,
        decimalLongitude      = decimalLongitude,
        station               = station,
        year                  = year,
        month                 = month,
        organismQuantity      = organismQuantity,
        concentration         = concentration,
        pcr_primer_lod        = pcr_primer_lod,
        detected              = detected,
        ownerContact          = project_contact,
        bibliographicCitation = bibliographicCitation,
        datasetID_obis        = datasetID_obis
      )
    out
  })
  obis_list <- Filter(Negate(is.null), obis_list)
  ## 6. Bind everything together ----
  obis_list <- purrr::compact(obis_list)
  if (length(obis_list) == 0L) {
    warning("No OBIS datasets with DNADerivedData and both present/absent occurrenceStatus returned any records for these filters.")
    return(tibble::tibble())
  }
  GOTeDNA_df <- dplyr::bind_rows(obis_list)
  rownames(GOTeDNA_df) <- NULL
  # saveRDS(GOTeDNA_df, "inst/app/data/temp_obis_data.rds")
  GOTeDNA_df_with_assigned_stations <- update_location_clusters(GOTeDNA_df)

  GOTeDNA_df_with_assigned_stations

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




##############NEW LIST SANDBOX HERE


old_cols <- c(
  "protocol_ID",
  "protocolVersion",
  "samp_name",
  "primer",
  "scientificName",
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "eventDate",
  "LClabel",
  "decimalLatitude",
  "decimalLongitude",
  "station",
  "year",
  "month",
  "organismQuantity",
  "concentration",
  "pcr_primer_lod",
  "detected",
  "project_contact",
  "bibliographicCitation",
  "datasetID_obis"
)

# protocol_ID           = protocol_ID,
# protocolVersion       = protocolVersion,
# samp_name             = as.character(samp_name),
# primer                = sprintf("%s | %s / %s", target_gene, pcr_primer_name_forward, pcr_primer_name_reverse),
# scientificName        = scientificName,
# kingdom               = kingdom,
# phylum                = phylum,
# class                 = class,
# order                 = order,
# family                = family,
# genus                 = genus,
# eventDate             = eventDate_clean,
# LClabel               = LClabel,
# decimalLatitude       = decimalLatitude,
# decimalLongitude      = decimalLongitude,
# station               = station,
# year                  = year,
# month                 = month,
# organismQuantity      = organismQuantity,
# concentration         = concentration,
# pcr_primer_lod        = pcr_primer_lod,
# detected              = detected,
# ownerContact          = project_contact,
# bibliographicCitation = bibliographicCitation,
# datasetID_obis        = datasetID_obis

protocol_columns <- c(
  'samp_size',
  'size_frac',
  'filter_material',
  'samp_mat_process',
  'samp_store_temp',
  'samp_store_sol',
  'target_gene',
  'pcr_primer_forward',
  'pcr_primer_reverse',
  'nucl_acid_ext_kit',
  'platform',
  'instrument',
  'seq_kit',
  'otu_db',
  'tax_assign_cat',
  'otu_seq_comp_appr'
  )

not_found <- setdiff(new_columns, names(core_and_extensions_debug))
#"platform"       "instrument"     "tax_assign_cat"


check_required_cols <- function(df, req_cols, dataset_id) {
  missing <- setdiff(req_cols, names(df))

  if (length(missing) > 0) {
    message(
      sprintf(
        'dataset_id: "%s" did not have required column(s): %s',
        dataset_id,
        paste(missing, collapse = ", ")
      )
    )
    return(FALSE)
  }

  TRUE
}

