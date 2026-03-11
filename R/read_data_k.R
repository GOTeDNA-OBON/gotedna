
read_data_k <- function(
    dataset_ids       = NULL,
    scientificname    = NULL,
    worms_id          = NULL,
    areaid            = NULL, #Canada: North Atlantic
    #geometry          = "POLYGON ((-67.72 40.614, -56.821 40.614, -56.821 47.279, -67.72 47.279, -67.72 40.614))",
    join_by           = c("auto", "occurrenceID", "id"),
    require_absences  = TRUE
) {

  join_by <- match.arg(join_by)

  # ---- columns you want back ----
  occurrence_cols <- c(
    "recordedBy","bibliographicCitation","materialSampleID",
    "organismQuantity","organismQuantityType",
    "samp_size", "samp_size_unit", "decimalLatitude", "decimalLongitude",
    "minimumDepthInMeters","maximumDepthInMeters","month","year",
    "scientificNameID","kingdom","phylum","class","order","family","genus",
    "dataset_id", "bathymetry", "associatedSequences", "bibliographicCitation"
  )

  dna_cols <- c(
    "id","dna_sequence","target_gene","pcr_primer_forward", "pcr_primer_reverse",
    "samp_name", "env_broad_scale","env_local_scale","env_medium","samp_mat_process",
    "size_frac","samp_size","samp_size_unit","otu_db","seq_kit", "otu_seq_comp_appr",
    "pcr_primer_name_forward","pcr_primer_name_reverse", "pcr_primer_reference",
    "occurrenceID"
  )

  mof_cols <- c(
    "id","seq_id","samp_category","checkls_ver","assay_name","assay_type",
    "targetTaxonomicAssay","geo_loc_name","technical_rep_id","project_contact",
    "seq_run_id","lib_id","project_id","pcr_0_1","samp_store_sol","samp_store_temp",
    "platform","instrument","tax_assign_cat","LClabel","occurrenceID",
    "nucl_acid_ext","nucl_acid_ext_kit","filter_material"
  )

  required_ext_cols <- c(
    "samp_name",
    "target_gene",
    "pcr_primer_name_forward",
    "pcr_primer_name_reverse"
  )

  added_cols <- c("category", "flags")

  mandatory_obis <- c(
    "occurrenceID","eventDate","decimalLongitude","decimalLatitude",
    "scientificName","occurrenceStatus","basisOfRecord"
  )

  cols_included_from_OBIS <- unique(c(
    occurrence_cols, dna_cols, mof_cols, added_cols, mandatory_obis
  ))

  # ---- helper: enforce columns (no external function needed) ----
  enforce_cols <- function(occ_all, cols) {
    # add missing columns as NA
    missing <- setdiff(cols, names(occ_all))
    if (length(missing) > 0) {
      for (m in missing) occ_all[[m]] <- NA
    }
    # keep only requested columns in a consistent order
    occ_all <- occ_all[, intersect(cols, names(occ_all)), drop = FALSE]
    occ_all
  }

  # ---- helper: build extension-wide table (MOF wide + DNA) ----
  join_extensions <- function(rec_ext, cols_included_from_OBIS, join_by) {
    if (is.null(rec_ext) || nrow(rec_ext) == 0L) return(NULL)

    # base cores
    core_occ_ext <- dplyr::distinct(rec_ext, occurrenceID, .keep_all = TRUE)
    core_id_ext  <- dplyr::distinct(rec_ext, id, .keep_all = TRUE)

    # DNADerivedData extension
    dna_only <- robis::unnest_extension(rec_ext, "DNADerivedData")
    shared_dna_cols <- intersect(cols_included_from_OBIS, names(dna_only))
    dna_only <- dplyr::select(dna_only, dplyr::any_of(shared_dna_cols))

    # MeasurementOrFact extension -> wide
    mof_only <- robis::unnest_extension(rec_ext, "MeasurementOrFact") %>%
      dplyr::group_by(occurrenceID, measurementType) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    wide_mof <- tidyr::pivot_wider(
      mof_only,
      id_cols    = c(occurrenceID, id),
      names_from = measurementType,
      values_from = measurementValue
    )

    shared_mof_cols <- intersect(cols_included_from_OBIS, names(wide_mof))
    wide_mof <- dplyr::select(wide_mof, dplyr::any_of(shared_mof_cols))

    mof_and_dna <- dplyr::left_join(wide_mof, dna_only, by = "id")

    # choose join key
    join_choice <- join_by
    if (join_choice == "auto") {
      can_occ <- "occurrenceID" %in% names(core_occ_ext) &&
        "occurrenceID" %in% names(mof_and_dna) &&
        any(!is.na(core_occ_ext$occurrenceID)) &&
        any(!is.na(mof_and_dna$occurrenceID))
      can_id  <- "id" %in% names(core_id_ext) &&
        "id" %in% names(mof_and_dna) &&
        any(!is.na(core_id_ext$id)) &&
        any(!is.na(mof_and_dna$id))
      if (can_occ) join_choice <- "occurrenceID"
      else if (can_id) join_choice <- "id"
      else stop("Neither occurrenceID nor id can be used to join core and extensions.")
    }

    out_ext <- if (join_choice == "occurrenceID") {
      dplyr::left_join(core_occ_ext, mof_and_dna, by = "occurrenceID")
    } else {
      dplyr::left_join(core_id_ext, mof_and_dna, by = "id")
    }

    out_ext <- dplyr::select(out_ext, dplyr::any_of(cols_included_from_OBIS))
    out_ext
  }

  # ---- 0) discover dataset_ids if NULL ----
  if (is.null(dataset_ids)) {
    message("Discovering datasets with DNADerivedData and MeasurementOrFact extensions ...")
    ds_tbl <- robis::dataset(
      scientificname = scientificname,
      areaid         = areaid,
      #geometry       = geometry,
      taxonid        = worms_id,
      hasextensions  = c("DNADerivedData", "MeasurementOrFact")
    ) %>%
      dplyr::filter(statistics$absence != 0)

    if (nrow(ds_tbl) == 0L) {
      warning("No datasets found that match scientificname/areaid AND have DNADerivedData + absences.")
      return(tibble::tibble())
    }

    if ("id" %in% names(ds_tbl)) dataset_ids <- unique(ds_tbl$id)
    else if ("datasetid" %in% names(ds_tbl)) dataset_ids <- unique(ds_tbl$datasetid)
    else stop("Could not find dataset id column in dataset() output.")

    message("Found ", length(dataset_ids), " dataset(s).")
  }

  dataset_ids <- as.character(dataset_ids)

  exclude_list <- c("NO_COORD","ZERO_COORD","LON_OUT_OF_RANGE","LAT_OUT_OF_RANGE","NO_MATCH")

  # ---- 1) loop datasets ----
  obis_list <- purrr::map(dataset_ids, function(ds) {

    Sys.sleep(3)
    message("Pulling OBIS dataset: ", ds)

    # defensive: check extensions exist
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

    # 1a) FULL CORE with absences (NO extensions)
    core_all <- tryCatch(
      robis::occurrence(
        datasetid      = ds,
        scientificname = scientificname,
        taxonid        = worms_id,
        areaid         = areaid,
        #geometry       = geometry,
        absence        = "include",
        dropped        = "include",
        exclude        = exclude_list
      ),
      error = function(e) {
        warning("Failed to fetch CORE (absences) for dataset ", ds, ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(core_all) || nrow(core_all) == 0L) {
      warning("No core occurrence records returned for dataset ", ds, ".")
      return(NULL)
    }

    core_all <- dplyr::distinct(core_all, occurrenceID, .keep_all = TRUE)

    if (!"occurrenceStatus" %in% names(core_all)) {
      warning("Dataset ", ds, " has no occurrenceStatus column; skipping.")
      return(NULL)
    }

    status_vals <- unique(na.omit(core_all$occurrenceStatus))
    if (require_absences && !all(c("present","absent") %in% status_vals)) {
      warning("Dataset ", ds, " does not contain both 'present' and 'absent'; skipping.")
      return(NULL)
    }

    core_all <- dplyr::select(core_all, dplyr::any_of(cols_included_from_OBIS))
    core_all <- enforce_cols(core_all, cols_included_from_OBIS)

    # 1b) EXTENSIONS (present-only; do NOT request absence="include")
    rec_ext <- tryCatch(
      robis::occurrence(
        datasetid      = ds,
        scientificname = scientificname,
        taxonid        = worms_id,
        areaid         = areaid,
        #geometry       = geometry,
        extensions     = c("DNADerivedData", "MeasurementOrFact"),
        hasextensions  = c("DNADerivedData", "MeasurementOrFact"),
        dropped        = "include",
        exclude        = exclude_list
      ),
      error = function(e) {
        warning("Failed to fetch EXTENSIONS for dataset ", ds, ": ", conditionMessage(e))
        NULL
      }
    )

    # If extensions fetch fails, return core (absences preserved)
    if (is.null(rec_ext) || nrow(rec_ext) == 0L) {
      core_all$samp_name <- as.character(dplyr::coalesce(core_all$samp_name, core_all$materialSampleID))
      return(core_all)
    }

    ext_joined <- join_extensions(rec_ext, cols_included_from_OBIS, join_by)

    ext_joined <- ext_joined %>%
      dplyr::distinct(.data$occurrenceID, .keep_all = TRUE)

    if (is.null(ext_joined) || nrow(ext_joined) == 0L) {
      core_all$samp_name <- as.character(dplyr::coalesce(core_all$samp_name, core_all$materialSampleID))
      return(core_all)
    }

    missing_required <- setdiff(required_ext_cols, names(ext_joined))
    if (length(missing_required) > 0) {
      warning(
        "Dataset ", ds,
        " missing required column(s): ",
        paste(missing_required, collapse = ", "),
        " ; skipping."
      )
      return(NULL)
    }

    # 1c) merge extensions onto full core (absences keep NA extension fields)
    out <- dplyr::left_join(core_all, ext_joined, by = "occurrenceID", suffix = c("", ".ext"))

    # coalesce duplicated columns (prefer extension values)
    dup_cols <- intersect(names(core_all), names(ext_joined))
    dup_cols <- setdiff(dup_cols, "occurrenceID")
    for (nm in dup_cols) {
      ext_nm <- paste0(nm, ".ext")
      if (ext_nm %in% names(out)) {
        out[[nm]] <- dplyr::coalesce(out[[ext_nm]], out[[nm]])
        out[[ext_nm]] <- NULL
      }
    }

    out <- dplyr::select(out, dplyr::any_of(cols_included_from_OBIS))
    out <- enforce_cols(out, cols_included_from_OBIS)
    out$samp_name <- as.character(dplyr::coalesce(out$samp_name, out$materialSampleID))
    out
  })

  obis_list <- purrr::compact(obis_list)

  if (length(obis_list) == 0L) {
    warning("No OBIS datasets returned any records for these filters.")
    return(tibble::tibble())
  }

  GOTeDNA_occ_all <- dplyr::bind_rows(obis_list)
  rownames(GOTeDNA_occ_all) <- NULL
  GOTeDNA_occ_all
}

