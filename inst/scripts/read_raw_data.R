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
  library(dplyr)

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
    "pcr_primer_reverse",
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

    # dataset_ids <- setdiff(dataset_ids, saved_ds)
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
      next
    }

    if (!"measurementorfact" %in% exts) {
      warning("Dataset ", ds, " has no MeasurementOrFact extension; skipping.")
      next
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
      next
    }

    if (nrow(rec) == 0L) {
      warning("No occurrence records returned for dataset ", ds, " with these filters.")
      next
    }



    # ---- 1c. Build core_occ and filter on occurrenceStatus ----
    core_occ <- rec %>%
      distinct(occurrenceID, .keep_all = TRUE)
    if (!"occurrenceStatus" %in% names(core_occ)) {
      warning("Dataset ", ds, " has no occurrenceStatus column; skipping.")
      next
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
        next
      }
    }
    # Keep also an id-based core for joining if needed
    core_id <- rec %>%
      distinct(id, .keep_all = TRUE)
    # DNADerivedData extension (includes `id` by default)
    dna_only <- robis::unnest_extension(rec, "DNADerivedData")

    shared_dna_cols <- intersect(cols_included_from_OBIS, names(dna_only))

    dna_only <- dna_only %>% select(shared_dna_cols)
    #MeasurementOfFact extension
    mof_only <- robis::unnest_extension(rec, "MeasurementOrFact")

    mof_only <- mof_only %>%
      group_by(occurrenceID, measurementType) %>%
      slice(1) %>%
      ungroup()
    # ---- 2. Decide how to join core + extension ----
    wide_mof <- mof_only %>%
      tidyr::pivot_wider(
        id_cols = c(occurrenceID, id),
        names_from = measurementType,
        values_from = measurementValue
      )

    shared_mof_cols <- intersect(cols_included_from_OBIS, names(wide_mof))

    wide_mof <- wide_mof %>% select(shared_mof_cols)
    mof_and_dna <- wide_mof %>% left_join(dna_only, by = "id")

    shared_cols <- intersect(cols_included_from_OBIS, names(rec))

    rec <- rec %>% select(shared_cols)

    shared_cols <- intersect(cols_included_from_OBIS, names(core_id))

    core_id <- core_id %>% select(shared_cols)

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

    # if (!has_required_cols(core_and_extensions, required_cols)) {
    #   return(NULL)
    # }

    if (is.null(core_and_extensions)) {
      message("core_and_extensions is NULL...moving to next")
      next
    }
    all_shared_cols <- intersect(names(core_and_extensions), cols_included_from_OBIS)

    core_and_extensions <- core_and_extensions %>% select(all_shared_cols)
    dup_names <- names(core_and_extensions)[duplicated(names(core_and_extensions))]
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

