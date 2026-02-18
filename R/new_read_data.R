read_data_test <- function(
    dataset_ids    = NULL,
    scientificname = NULL,
    worms_id       = NULL,
    areaid         = 34,                               #Canada: North Atlantic
    geometry       = "POLYGON ((-67.72 40.614, -56.821 40.614, -56.821 47.279, -67.72 47.279, -67.72 40.614))",
    join_by        = c("auto", "occurrenceID", "id"),
    require_absences = TRUE,
    require_dna_sequence = FALSE
) {
  library(robis)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(purrr)
  library(dbscan)
  library(geosphere)
  library(sf)


  cols_not_needed_at_pull <- c(
    "eventDate_clean", "year","month","datasetID_obis", "primer"
  )

  extra_cols_needed_at_pull <- c("dna", "mof", "absence", "occurrenceStatus", "materialSampleID", "eventDate")
  need_cols <- c(
    "samp_name","target_gene","pcr_primer_name_forward","pcr_primer_name_reverse",
    "occurrenceID", "id", "occurrenceStatus","basisOfRecord","scientificName","scientificNameID",
    "kingdom","phylum","class","order","family","genus",
    "eventDate_clean","LClabel","decimalLatitude","decimalLongitude","year","month",
    "organismQuantity","samp_size","size_frac","filter_material","samp_mat_process",
    "samp_store_temp","samp_store_sol","nucl_acid_ext_kit",
    "platform","instrument","seq_kit","otu_db","tax_assign_cat","otu_seq_comp_appr",
    "minimumDepthInMeters","maximumDepthInMeters","project_contact","bibliographicCitation",
    "datasetID_obis","category","hab", "samp_name", "target_subfragment", "target_gene", "primer"
  )

  keep_at_pull_cols <- setdiff(need_cols, cols_not_needed_at_pull)
  keep_at_pull_cols <- c(keep_at_pull_cols, extra_cols_needed_at_pull)

  join_by <- match.arg(join_by)
  ## 0. If dataset_ids is NULL, discover them via robis::dataset() ----
  if (is.null(dataset_ids)) {
    message("Discovering datasets with DNADerivedData and MeasurementOrFact extensions ...")
    ds_tbl <- robis::dataset(
      scientificname = scientificname,
      areaid         = areaid,
      geometry       = geometry,
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
    rec_with_absences <- robis::occurrence(
      datasetid      = ds,
      scientificname = scientificname,
      taxonid        = worms_id,
      areaid         = areaid,
      geometry       = geometry,
      absence        = "include"
    ) %>%
      select(any_of(keep_at_pull_cols))

    if (nrow(rec_with_absences) == 0L) {
      warning("No occurrence records returned for dataset ", ds, " with these filters.")
      return(NULL)
    }
    # Unique non-NA status values
    status_vals <- unique(na.omit(rec_with_absences$occurrenceStatus))
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
    rec_with_absences <- rec_with_absences %>% filter(absence == TRUE)


    rec_with_extensions <- robis::occurrence(
      datasetid      = ds,
      scientificname = scientificname,
      taxonid        = worms_id,
      areaid         = areaid,
      geometry       = geometry,
      extensions  = c("DNADerivedData", "MeasurementOrFact"),
      hasextensions  = c("DNADerivedData", "MeasurementOrFact")
    ) %>%
      select(any_of(keep_at_pull_cols))
    if (nrow(rec_with_extensions) == 0L) {
      warning("No extension records returned for dataset ", ds, " with these filters.")
      return(NULL)
    }

    # ---- 1c. Build core_occ and filter on occurrenceStatus ----
    core_occ <- rec_with_extensions %>%
      distinct(occurrenceID, .keep_all = TRUE)
    if (!"occurrenceStatus" %in% names(core_occ)) {
      warning("Dataset ", ds, " has no occurrenceStatus column; skipping.")
      return(NULL)
    }

    # Keep also an id-based core for joining if needed
    core_id <- rec_with_extensions %>%
      distinct(id, .keep_all = TRUE)
    core_id <- core_id %>% select(-mof, -dna)

    # DNADerivedData extension (includes `id` by default)
    dna_only <- robis::unnest_extension(core_occ, "DNADerivedData")

    dna_only <- dna_only %>%
      select(any_of(keep_at_pull_cols))

    #MeasurementOfFact extension
    mof_only <- robis::unnest_extension(core_occ, "MeasurementOrFact")

    core_occ <- core_occ %>% select(-mof, -dna)
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
      )%>%
      select(any_of(keep_at_pull_cols))

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
      core_and_extensions <- core_occ %>% left_join(mof_and_dna %>% select(-id), by = "occurrenceID")
    } else {  # "id"
      core_and_extensions <- core_id %>%
        left_join(mof_and_dna %>% select(-occurrenceID), by = "id")
    }

    core_and_extensions <- core_and_extensions %>%
      mutate(organismQuantity = as.numeric(organismQuantity))
    rec_with_absences <- rec_with_absences %>%
      mutate(organismQuantity = as.numeric(organismQuantity))
    if (!require_dna_sequence) {
      group_cols <- c(
        "materialSampleID",
        "scientificName",
        "target_gene",
        "pcr_primer_name_forward",
        "pcr_primer_name_reverse")

      # all other columns that will take first() — includes occurrenceID and id
      first_cols <- setdiff(names(core_and_extensions), c(group_cols, "organismQuantity"))

      core_and_extensions <- core_and_extensions %>%
        group_by(across(all_of(group_cols))) %>%
        summarise(
          organismQuantity = sum(as.numeric(organismQuantity), na.rm = TRUE),
          across(all_of(first_cols), first),
          .groups = "drop"
        )
    }



    #HERE WE JOIN THE ABSENCES TO THE WIDER EXTENSIONS DF THAT HAS ONLY POSITIVE DETECTIONS
    extra_cols <- setdiff(names(core_and_extensions), names(rec_with_absences))

    # Exclude columns that are per-row (DNA, organismQuantity, id, etc.)
    metadata_cols <- setdiff(extra_cols, c("DNA_sequence"))

    metadata_ref <- core_and_extensions %>%
      select(all_of(group_cols), all_of(metadata_cols)) %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(across(everything(), first), .groups = "drop")
    with_absence_df_filled <- rec_with_absences %>%
      left_join(metadata_ref, by = c("materialSampleID", "scientificName"))

    dna_cols <- setdiff(names(core_and_extensions), names(with_absence_df_filled))
    if (length(dna_cols) > 0) {
      with_absence_df_filled[dna_cols] <- NA
    }
    core_and_extensions <- bind_rows(core_and_extensions, with_absence_df_filled)

    # ---- 3. Basic cleaning & derived fields ----
    core_and_extensions <- core_and_extensions %>%
      filter(!is.na(samp_name)) %>%  # keep only real samples
      mutate(
        datasetID_obis = ds,
        # eventDate usually ISO; strip time & parse
        eventDate_clean = suppressWarnings(
          ymd(substr(as.character(eventDate), 1, 10))
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
    # Helper: makes sure data will missing columns will still be pulled, adds NA values
    ensure_cols <- function(df, cols) {
      missing <- setdiff(cols, names(df))
      if (length(missing)) {
        df[missing] <- NA
      }
      df
    }

    core_and_extensions <- ensure_cols(core_and_extensions, need_cols)

    out <- core_and_extensions %>%
      mutate(
        samp_name = as.character(samp_name)
      ) %>%
      select(all_of(need_cols)) %>%
      rename(
        eventDate = eventDate_clean,
        ownerContact = project_contact
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
  # saveRDS(GOTeDNA_df, "inst/app/data/temp_obis_data.rds")
  # GOTeDNA_df_with_assigned_stations <- update_location_clusters(GOTeDNA_df)

  # GOTeDNA_df_with_assigned_stations
  GOTeDNA_df
}



#does not work
#7ede12fd-c420-4226-8326-86c307c9f80e

#does work
#1952bb10-d56f-40c6-abd3-33dc1362bda9


find_variation <- function(df, group_cols = c("samp_name", "scientificName")) {

  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      across(
        everything(),
        ~ n_distinct(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = -all_of(group_cols),
      names_to = "column",
      values_to = "n_distinct"
    ) %>%
    filter(n_distinct > 1) %>%
    arrange(desc(n_distinct))
}


