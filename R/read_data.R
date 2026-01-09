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
    # ---- 5. Return in GOTeDNA_df-like shape ----
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
        ownerContact          = ownerContact,
        bibliographicCitation = bibliographicCitation,
        datasetID_obis        = datasetID_obis
      )
    out
  })
  ## 6. Bind everything together ----
  obis_list <- purrr::compact(obis_list)
  if (length(obis_list) == 0L) {
    warning("No OBIS datasets with DNADerivedData and both present/absent occurrenceStatus returned any records for these filters.")
    return(tibble::tibble())
  }
  GOTeDNA_df <- dplyr::bind_rows(obis_list)
  rownames(GOTeDNA_df) <- NULL

  GOTeDNA_df_with_assigned_stations <- update_station_variable(GOTeDNA_df)

  GOTeDNA_df_with_assigned_stations

}

update_station_variable <- function(df,
                                    lat_col  = "decimalLatitude",
                                    long_col = "decimalLongitude",
                                    eps_km   = 0.5,
                                    minPts  = 2) {

  # ---- 1. Ensure stationLabel exists, but do NOT overwrite it ----
  if (!"stationLabel" %in% names(df)) {
    if ("station" %in% names(df)) {
      df$stationLabel <- df$station
    } else {
      df$stationLabel <- NA_character_
    }
  }

  # ---- 2. Valid coordinate rows ----
  valid_idx <- which(
    !is.na(df[[lat_col]]) &
      !is.na(df[[long_col]])
  )

  if (length(valid_idx) < minPts) {
    df$station <- NA_character_
    return(df)
  }

  # ---- 3. Project coordinates to meters ----
  sf_pts <- sf::st_as_sf(
    df[valid_idx, ],
    coords = c(long_col, lat_col),
    crs = 4326
  )

  sf_pts <- sf::st_transform(sf_pts, 3857)
  coords_m <- sf::st_coordinates(sf_pts)

  # ---- 4. DBSCAN clustering ----
  eps_m <- eps_km * 1000
  db <- dbscan::dbscan(coords_m, eps = eps_m, minPts = minPts)

  # ---- 5. Assign clusters ----
  df$station <- NA_character_

  clusters <- db$cluster

  # Noise points (0) → unique IDs
  if (any(clusters == 0)) {
    noise_ids <- which(clusters == 0)
    max_cluster <- max(clusters)

    clusters[noise_ids] <- seq(
      from = max_cluster + 1,
      length.out = length(noise_ids)
    )
  }

  df$station[valid_idx] <- as.character(clusters)

  df
}

