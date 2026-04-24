#Canadian MPA Network Shapefiles and planning regions extracted from the Canadian Database of Protected and Conserved Areas. This data was filtered for the three MPAs targeted in the Maritimes region.

#Function to find St. Anns Bank, Musquash, Gully (MPA) and Eastern Shore Islands network, Fundian Channel-BrownsBank (AOI) polygons + transform to WGS84 for leaflet
read_poly_wgs84 <- function(...) {
  st_read(file.path("data", "polygons", ...), quiet = TRUE) %>%
    st_transform(4326)
}

mpa_targets <- read_poly_wgs84("mpa_targets.shp")
esi_poly    <- read_poly_wgs84("EasternShoreIslands_networksite.shp")
fcbb_poly   <- read_poly_wgs84("FCBB_Proposed_MPA_Boundary_zones.shp")

#Combine all polygons into one sf object
#Add a common name/type to each polygon set
mpa_polys <- mpa_targets %>%
  st_geometry() %>%                     # keep just geometry
  st_as_sf() %>%                        # convert back to sf
  mutate(
    site_name = mpa_targets$NAME_E,
    site_type = "MPA"
  )

esi_polys <- esi_poly %>%
  st_geometry() %>%
  st_as_sf() %>%
  mutate(
    site_name = "Eastern Shore Islands AOI",
    site_type = "AOI"
  )

fcbb_polys <- fcbb_poly %>%
  st_geometry() %>%
  st_as_sf() %>%
  mutate(
    site_name = "Fundian Channel–Browns Bank AOI",
    site_type = "AOI"
  )

#Bind everything together
all_polys <- bind_rows(
  mpa_polys,
  esi_polys,
  fcbb_polys
)

#Ensure it's still sf
all_polys <- st_as_sf(all_polys)

# ---- Clean once ----
all_polys_clean <- all_polys %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON") %>%
  (\(x) x[!st_is_empty(x), ])() %>%
  st_as_sf()

# ---- A) zone outlines (keep every piece) ----
all_polys_zones <- all_polys_clean

# ---- B) clickable dissolved outlines (one feature per site) ----
all_polys_click <- all_polys_clean %>%
  group_by(site_type, site_name) %>%
  summarise(geometry = st_union(x), .groups = "drop") %>%
  st_as_sf()


#-----------------------------------------------------------------------------
read_data <- function(
    dataset_ids       = NULL,
    scientificname    = NULL,
    worms_id          = NULL,
    areaid            = "34", #Canada: North Atlantic
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

STORED_DATA <- read_data()
saveRDS(STORED_DATA, "./data/OBIS_data.rds")
STORED_DATA <- readRDS("./data/OBIS_data.rds")

################DO NOT CHANGE ABOVE CODE



#Connect code below to stored data object

# ---- standardize types early ----
occ_all <- STORED_DATA %>%
  dplyr::mutate(
    year             = as.character(year),
    samp_name         = as.character(samp_name),
    occurrenceStatus  = tolower(as.character(occurrenceStatus)),
    decimalLatitude   = suppressWarnings(as.numeric(decimalLatitude)),
    decimalLongitude  = suppressWarnings(as.numeric(decimalLongitude)),
    target_gene       = dplyr::case_when(
      stringr::str_detect(tolower(target_gene), "12s") ~ "12S",
      stringr::str_detect(tolower(target_gene), "coi") ~ "COI",
      stringr::str_detect(tolower(target_gene), "16s") ~ "16S",
      stringr::str_detect(tolower(target_gene), "18s") ~ "18S",
      TRUE ~ as.character(target_gene)
    )
  )

# ---- Build available gene-year keys dynamically ----
KEY_TBL <- occ_all %>%
  dplyr::filter(!is.na(year), year != "", !is.na(target_gene), target_gene != "") %>%
  dplyr::distinct(target_gene, year) %>%
  dplyr::mutate(
    target_gene = as.character(target_gene),
    year        = as.character(year),
    key         = paste(target_gene, year, sep = "_")
  ) %>%
  dplyr::arrange(target_gene, year)

KEYS <- KEY_TBL$key

# ---- Split the *data* by key ----
DATA_BY_KEY <- occ_all %>%
  filter(!is.na(year), year != "", !is.na(target_gene), target_gene != "") %>%
  mutate(
    target_gene = as.character(target_gene),
    year        = as.character(year),
    key         = paste(target_gene, year, sep = "_")
  ) %>%
  group_by(key) %>%
  group_split(.keep = TRUE)

names(DATA_BY_KEY) <- KEYS

# ---- Convert each group to SF points (PRESENT only) ----
points_sf_from_occ_all <- function(occ_all) {
  occ_all %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
           tolower(as.character(occurrenceStatus)) == "present") %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
}

SPECIES_SF_BY_KEY <- purrr::imap(DATA_BY_KEY, ~{
  parts <- strsplit(.y, "_")[[1]]
  gene <- parts[1]
  yr   <- parts[2]

  points_sf_from_occ_all(.x) %>%
    mutate(target_gene = gene, year = yr)
})

species_sf_all <- dplyr::bind_rows(SPECIES_SF_BY_KEY)

#############################################################
#Turn species data into sf points  for the spatial join


# ---- ONE join: species points -> clickable MPA/AOI polygons ----
species_in_polys_all <- sf::st_join(
  species_sf_all,
  all_polys_click,          # dissolved clickable polygons
  join = sf::st_within,
  left = FALSE
) %>%
  dplyr::mutate(
    year = as.character(year),
    target_gene = as.character(target_gene),
    scientificName = as.character(scientificName)
  ) %>%
  dplyr::filter(!is.na(scientificName), scientificName != "")

# year-aware species list per polygon
poly_species_year <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(site_name, site_type, year, scientificName) %>%
  dplyr::group_by(site_name, site_type, year) %>%
  dplyr::summarise(species = list(sort(unique(scientificName))), .groups = "drop")

# all-years species list per polygon
poly_species_all <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(site_name, site_type, scientificName) %>%
  dplyr::group_by(site_name, site_type) %>%
  dplyr::summarise(species = list(sort(unique(scientificName))), .groups = "drop")

total_species_year <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(site_name, site_type, year) %>%
  dplyr::summarise(n_species_total = dplyr::n_distinct(scientificName), .groups = "drop")

total_species_all <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(site_name, site_type) %>%
  dplyr::summarise(n_species_total = dplyr::n_distinct(scientificName), .groups = "drop")

species_by_class_year <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(site_name, site_type, year, class) %>%
  dplyr::summarise(n_species = dplyr::n_distinct(scientificName), .groups = "drop")

species_by_class_wide_year <- species_by_class_year %>%
  tidyr::pivot_wider(names_from = class, values_from = n_species, values_fill = 0)

species_by_class_all <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(site_name, site_type, class) %>%
  dplyr::summarise(n_species = dplyr::n_distinct(scientificName), .groups = "drop")

species_by_class_wide_all <- species_by_class_all %>%
  tidyr::pivot_wider(names_from = class, values_from = n_species, values_fill = 0)






# optional convenience table for "All years combined"

species_in_polys_all %>%
  st_drop_geometry() %>%
  count(site_name, site_type, year, target_gene)

poly_species_all_from_year <- poly_species_year %>%
  group_by(site_name, site_type) %>%
  summarise(species = list(sort(unique(unlist(species)))), .groups = "drop")


#Summary report per polygon




##Species Richness Polygons

# 1) Make a grid over polygons
# 1) Work in a projected CRS (Canada Lambert is a good default)
##Species Richness Polygons (Option A: projected CRS to avoid s2 errors)

##Species Richness Polygons (Option A: build grid in EPSG:4326 so it stays upright in leaflet)

# 0) Work in lon/lat (leaflet native)
crs_ll <- 4326

# 1) Clean + dissolve polygons in 4326
poly_union_ll <- all_polys_click %>%
  sf::st_make_valid() %>%
  sf::st_transform(crs_ll) %>%
  sf::st_union() %>%
  sf::st_as_sf()

# 2) Choose an "approx 2000 m" grid size expressed in degrees at your latitude
cell_m <- 2000

# centroid latitude (used to approximate meters->degrees)
cent <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(poly_union_ll)))
lat0 <- mean(cent[,2], na.rm = TRUE)

deg_per_m_lat <- 1 / 111320
deg_per_m_lon <- 1 / (111320 * cos(lat0 * pi/180))

cellsize_deg <- c(cell_m * deg_per_m_lon, cell_m * deg_per_m_lat)  # c(lon_deg, lat_deg)

# 3) Build grid in 4326 (upright in leaflet)
grid_ll <- sf::st_make_grid(
  poly_union_ll,
  cellsize = cellsize_deg,
  square   = TRUE
) %>%
  sf::st_as_sf() %>%
  dplyr::mutate(cell_id = dplyr::row_number())

# 4) Clip grid ONCE (still 4326)
grid_clip_ll <- sf::st_intersection(grid_ll, poly_union_ll) %>%
  sf::st_make_valid() %>%
  sf::st_collection_extract("POLYGON") %>%
  (\(x) x[!sf::st_is_empty(x), ])() %>%
  sf::st_as_sf()

# (Leaflet uses 4326 anyway)
grid_clip <- grid_clip_ll

# 5) Points stay in 4326 too
SPECIES_SF_BY_KEY_ll <- SPECIES_SF_BY_KEY %>%
  purrr::map(~sf::st_transform(.x, crs_ll))

#Make an "all points" sf in 4326 (PRESENT only)
species_sf_all_ll <- dplyr::bind_rows(SPECIES_SF_BY_KEY_ll)

# 6) Richness in 4326 (no CRS mismatch)
make_richness_layer_fast <- function(grid_sf, pts_sf) {
  stopifnot(inherits(grid_sf, "sf"), inherits(pts_sf, "sf"))
  if (!"scientificName" %in% names(pts_sf)) {
    stop("pts_sf must contain a 'scientificName' column.")
  }

  idx <- sf::st_intersects(grid_sf, pts_sf)

  n_sp <- vapply(idx, function(i) {
    if (length(i) == 0) return(0L)
    length(unique(as.character(pts_sf$scientificName[i])))
  }, integer(1))

  grid_sf %>%
    dplyr::mutate(
      n_species    = n_sp,
      has_sampling = n_sp > 0
    )
}

RICHNESS_BY_KEY <- purrr::map(
  SPECIES_SF_BY_KEY_ll,
  ~make_richness_layer_fast(grid_clip_ll, .x)
)

# optional: all-years combined richness
RICHNESS_ALL <- make_richness_layer_fast(
  grid_clip_ll,
  dplyr::bind_rows(SPECIES_SF_BY_KEY_ll)
)

# species_sf_all: sf POINTS with at least cell_id, target_gene, scientificName (and year)
# grid_clip: sf POLYGONS with cell_id + geometry  (your analysis grid)

build_gene_all_grid <- function(gene, grid_sf, pts_sf, tax_col = "scientificName") {

  pts_g <- pts_sf %>%
    dplyr::filter(as.character(target_gene) == gene) %>%
    dplyr::mutate(
      cell_id = as.integer(cell_id),
      taxon   = as.character(.data[[tax_col]])
    ) %>%
    dplyr::filter(!is.na(cell_id), !is.na(taxon), taxon != "")

  # unique taxa per cell across ALL years
  rich_tbl <- pts_g %>%
    sf::st_drop_geometry() %>%
    dplyr::distinct(cell_id, taxon) %>%
    dplyr::count(cell_id, name = "n_species")

  out <- grid_sf %>%
    dplyr::mutate(cell_id = as.integer(cell_id)) %>%
    dplyr::left_join(rich_tbl, by = "cell_id") %>%
    dplyr::mutate(
      n_species    = tidyr::replace_na(n_species, 0L),
      has_sampling = n_species > 0
    )

  out
}

# Use the lon/lat versions you created for the grid + points
# grid_clip_ll and species_sf_all_ll are both EPSG:4326

species_sf_all_cell <- sf::st_join(
  species_sf_all_ll,
  grid_clip_ll %>% dplyr::select(cell_id),
  join = sf::st_within,
  left = FALSE
)

# Now this will work because pts_sf has cell_id
RICHNESS_GENE_ALL <- list(
  "12S" = build_gene_all_grid("12S", grid_clip_ll, species_sf_all_cell, tax_col = "scientificName"),
  "COI" = build_gene_all_grid("COI", grid_clip_ll, species_sf_all_cell, tax_col = "scientificName"),
  "16S" = build_gene_all_grid("16S", grid_clip_ll, species_sf_all_cell, tax_col = "scientificName"),
  "18S" = build_gene_all_grid("18S", grid_clip_ll, species_sf_all_cell, tax_col = "scientificName")
)

# optional: keep dissolved polygon union for leaflet
poly_union <- poly_union_ll


## ===== Unified richness palettes + leaflet map (Option A: 4326-only) =====

# Assumes these already exist from your Option A block:
# - grid_clip_ll (sf, EPSG:4326) with column cell_id
# - SPECIES_SF_BY_KEY_ll (named list of sf points, EPSG:4326) with scientificName, year, target_gene
# - RICHNESS_BY_KEY (named list of sf grids, EPSG:4326) with n_species
# - RICHNESS_ALL (sf grid, EPSG:4326) with n_species (from make_richness_layer_fast)
# - poly_union_ll or all_polys_click for outlines

# 0) Ensure grid has cell_id
if (!"cell_id" %in% names(grid_clip_ll)) {
  grid_clip_ll <- grid_clip_ll %>% dplyr::mutate(cell_id = dplyr::row_number())
}

# 1) Make an "all points" sf in 4326 (present only; should already be present-only, but keep safe)
species_sf_all_ll <- dplyr::bind_rows(SPECIES_SF_BY_KEY_ll) %>%
  dplyr::mutate(
    year             = as.character(year),
    target_gene      = as.character(target_gene),
    scientificName   = as.character(scientificName),
    occurrenceStatus = tolower(as.character(occurrenceStatus))
  ) %>%
  dplyr::filter(
    is.na(occurrenceStatus) | occurrenceStatus == "present",
    !is.na(scientificName), scientificName != ""
  )

# 2) Rename ALL grid column consistently for mapping (“total across whatever you used to compute it”)
# Your RICHNESS_ALL currently has n_species from make_richness_layer_fast().
RICHNESS_ALL <- RICHNESS_ALL %>%
  dplyr::rename(n_species_total = n_species)

# (Optional) If you want per-year ALL-markers layers in 4326:
years_all <- sort(unique(na.omit(species_sf_all_ll$year)))

RICHNESS_ALL_BY_YEAR <- setNames(
  lapply(years_all, function(yr) {
    pts_y <- species_sf_all_ll %>% dplyr::filter(year == yr)
    make_richness_layer_fast(grid_clip_ll, pts_y) %>%
      dplyr::rename(n_species_total = n_species)
  }),
  years_all
)

# 3) Shared richness palette domain across ALL layers
max_rich <- max(
  RICHNESS_ALL$n_species_total,
  unlist(purrr::map(RICHNESS_BY_KEY, ~ .x$n_species)),
  unlist(purrr::map(RICHNESS_ALL_BY_YEAR, ~ .x$n_species_total)),
  na.rm = TRUE
)

rich_domain <- c(0, max_rich)

wes_cont <- function(name, n = 100) {
  grDevices::colorRampPalette(
    wesanderson::wes_palette(name, type = "continuous")
  )(n)
}

pal_vec <- wes_cont("Zissou1", 100)

pal_rich <- leaflet::colorNumeric(
  palette  = pal_vec,
  domain   = rich_domain,
  na.color = "transparent"
)

# 4) Cell -> species list tables (ALL + by key) using *grid_clip_ll* (not projected!)
idx_all <- sf::st_intersects(grid_clip_ll, species_sf_all_ll)

CELL_SPECIES_ALL <- tibble::tibble(cell_id = grid_clip_ll$cell_id) %>%
  dplyr::mutate(
    spp = purrr::map(idx_all, \(i) sort(unique(species_sf_all_ll$scientificName[i])))
  )

CELL_SPECIES_BY_KEY <- purrr::imap(SPECIES_SF_BY_KEY_ll, ~{
  idx <- sf::st_intersects(grid_clip_ll, .x)
  tibble::tibble(cell_id = grid_clip_ll$cell_id) %>%
    dplyr::mutate(
      spp = purrr::map(idx, \(i) sort(unique(.x$scientificName[i])))
    )
})

# Long table: one row per (cell_id, scientificName, target_gene, year)
cell_species_all <- purrr::imap_dfr(SPECIES_SF_BY_KEY_ll, ~{
  key <- .y
  parts <- strsplit(key, "_")[[1]]
  gene <- parts[1]
  yr   <- parts[2]

  idx <- sf::st_intersects(grid_clip_ll, .x)

  tibble::tibble(cell_id = grid_clip_ll$cell_id) %>%
    dplyr::mutate(scientificName = purrr::map(idx, \(i) unique(.x$scientificName[i]))) %>%
    tidyr::unnest(scientificName) %>%
    dplyr::mutate(target_gene = gene, year = yr) %>%
    dplyr::filter(!is.na(scientificName), scientificName != "")
})

cell_species_total <- cell_species_all %>%
  dplyr::group_by(cell_id) %>%
  dplyr::summarise(n_species_total = dplyr::n_distinct(scientificName), .groups = "drop")

# 5) Convenience: split RICHNESS_BY_KEY into gene buckets
grid_12S_by_year <- RICHNESS_BY_KEY[grep("^12S_", names(RICHNESS_BY_KEY))]
grid_COI_by_year <- RICHNESS_BY_KEY[grep("^COI_", names(RICHNESS_BY_KEY))]
grid_16S_by_year <- RICHNESS_BY_KEY[grep("^16S_", names(RICHNESS_BY_KEY))]
grid_18s_by_year <- RICHNESS_BY_KEY[grep("^18S_", names(RICHNESS_BY_KEY))]

# 6) Helper: pick which richness layer to display
# gene: "All" / "12S" / "COI" / "16S" / "18S"
# year: "All" or specific year string
get_richness_layer <- function(gene = "All", year = "All") {
  gene <- as.character(gene)
  year <- as.character(year)

  if (gene == "All") {
    if (year == "All") return(RICHNESS_ALL)
    if (year %in% names(RICHNESS_ALL_BY_YEAR)) return(RICHNESS_ALL_BY_YEAR[[year]])
    return(RICHNESS_ALL)
  }

  # gene-specific layers are only defined for specific years (no "All-years" gene grid)
  if (year == "All") return(NULL)

  key <- paste(gene, year, sep = "_")
  if (key %in% names(RICHNESS_BY_KEY)) return(RICHNESS_BY_KEY[[key]])

  NULL
}

# 7) Defaults
default_gene <- "All"   # or "12S", "COI"
default_year <- "All"   # or "2024", etc.

init_layer <- get_richness_layer(default_gene, default_year)

# Optionally: a safe init if gene-specific All-years returns NULL
if (is.null(init_layer)) init_layer <- RICHNESS_ALL



#Define polygon layers

# Richness layers keyed by gene+year
selected_rich_key <- function(gene, year) {
  gene <- as.character(gene); year <- as.character(year)

  if (gene == "All") return("All")        # handle via RICHNESS_ALL / RICHNESS_ALL_BY_YEAR
  if (year == "All") return(NULL)         # no gene-specific layer for All-years
  paste(gene, year, sep = "_")            # e.g., "12S_2024"
}


#Organize per year layers into names lists

# --- Richness grids by year ---
grid_12S_by_year <- RICHNESS_BY_KEY[grep("^12S_", names(RICHNESS_BY_KEY))]
grid_COI_by_year <- RICHNESS_BY_KEY[grep("^COI_", names(RICHNESS_BY_KEY))]
grid_16S_by_year <- RICHNESS_BY_KEY[grep("^16S_", names(RICHNESS_BY_KEY))]
grid_18S_by_year <- RICHNESS_BY_KEY[grep("^18S_", names(RICHNESS_BY_KEY))]

grid_ALL_static <- RICHNESS_ALL


## ===== Shared richness colour scale across ALL layers =====

# 1) Collect maxima safely (handles NULLs / missing cols)
max_from_key_layers <- function(lst, col) {
  vals <- unlist(purrr::map(lst, ~{
    if (is.null(.x) || !col %in% names(.x)) return(NA_real_)
    suppressWarnings(as.numeric(.x[[col]]))
  }))
  suppressWarnings(max(vals, na.rm = TRUE))
}

max_all_markers <- max(
  suppressWarnings(max(as.numeric(RICHNESS_ALL$n_species_total), na.rm = TRUE)),
  max_from_key_layers(RICHNESS_ALL_BY_YEAR, "n_species_total"),
  na.rm = TRUE
)

max_gene_layers <- max_from_key_layers(RICHNESS_BY_KEY, "n_species")

max_rich <- max(max_all_markers, max_gene_layers, na.rm = TRUE)

# 2) Build ONE palette with ONE domain
rich_domain <- c(0, max_rich)


##NOTE: Right now the SARA Schedule 1 filter is matching based on the data inputted into the app. Once linking the code to OBIS data, edit so that it matches based on WoRMS AphiaID


#Code for the SARA Schedule 1 list (Cleanup and WoRMS linkage)

#Call in file from OBIS_Prep folder
#SARA <- read.xlsx(file.path("data", "sara_ais", "SARA_Clean_Schedule1_specieslist.xlsx"))

SARA <- read.xlsx(file.path("data", "sara_ais", "SARA_Clean_Schedule1_specieslist_noPacific.xlsx"))
AIS  <- read.xlsx(file.path("data", "sara_ais", "Target_AIS_List_Claudio.xlsx"))

# helper: safe lookup for one name
worms_lookup_one <- function(x) {
  if (is.na(x) || !nzchar(x)) return(tibble(AphiaID = NA_integer_, worms_status = NA_character_, worms_valid_name = NA_character_))

  res <- tryCatch(
    worrms::wm_records_name(name = x, fuzzy = TRUE, marine_only = TRUE),
    error = function(e) NULL
  )

  if (is.null(res) || nrow(res) == 0) {
    tibble(AphiaID = NA_integer_, worms_status = NA_character_, worms_valid_name = NA_character_)
  } else {
    tibble(
      AphiaID         = res$AphiaID[[1]],
      worms_status    = res$status[[1]],      # e.g., "accepted", "unaccepted"
      worms_valid_name= res$valid_name[[1]]   # accepted name WoRMS points to
    )
  }
}

##Clean up columns

worms_tbl_AIS <- AIS %>%
  distinct(Scientific.Name) %>%
  mutate(worms = map(Scientific.Name, worms_lookup_one)) %>%
  tidyr::unnest(worms)

AIS <- AIS %>%
  left_join(worms_tbl_AIS, by = "Scientific.Name") %>%
  mutate(
    worms_match = !is.na(AphiaID)
  )

sample_tag <- function(occ_all) {  #NEW
  occ_all %>%
    dplyr::mutate(
      yr = dplyr::if_else(is.na(year), "", as.character(year)),
      mk = dplyr::if_else(is.na(target_gene), "", as.character(target_gene)),
      tag = paste0(
        samp_name,
        dplyr::if_else(
          yr != "" | mk != "",
          paste0(" (", paste(c(yr, mk)[c(yr, mk) != ""], collapse = ", "), ")"),
          ""
        )
      )
    ) %>%
    dplyr::pull(tag) %>%
    unique() %>%
    sort()
}


##Fix grid heatmap by year for all target_genes

grid_ALL_by_year <- RICHNESS_ALL_BY_YEAR   # already in 4326
grid_ALL_static  <- RICHNESS_ALL           # all years combined, already in 4326


# ----------------------------
# Depth colour mapping (3-zone ramps)
# ----------------------------

cool <- c(
  "#ccebc5", "#a8ddb5", "#7fcdbb",
  "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"
)

depth_breaks <- c(0, 25, 50, 75, 100, 500, 1000, Inf)
depth_labels <- c(
  "0–25 m", "26–50 m", "51–75 m", "76–100 m",
  "101–500 m", "501–1000 m", "1000+ m"
)

orange3 <- c("#fed976", "#feb24c", "#fd8d3c")                   # mixed depth overlay

orange_pal    <- grDevices::colorRampPalette(orange3)


get_depth_num <- function(occ_all) {
  depth_col <- dplyr::case_when(
    "minimumDepthInMeters" %in% names(occ_all) ~ "minimumDepthInMeters",
    "Depth" %in% names(occ_all) ~ "Depth",
    TRUE ~ NA_character_
  )
  if (is.na(depth_col)) return(occ_all %>% mutate(depth_m = NA_real_))

  occ_all %>% mutate(depth_m = suppressWarnings(as.numeric(.data[[depth_col]])))
}


# map depth (m) -> binned
depth_to_col <- function(d) {
  out <- rep(NA_character_, length(d))
  ok  <- !is.na(d)

  # findInterval returns 1..7 for these breaks
  idx <- findInterval(d[ok], vec = depth_breaks, rightmost.closed = TRUE)

  # idx could be 0 if d < 0; clamp just in case
  idx <- pmax(1, pmin(length(cool), idx))

  out[ok] <- cool[idx]
  out
}


# mixed depth -> orange shade (more orange with bigger within-cell range)
mixed_to_orange <- function(depth_range, cap_range = 200) {
  out <- rep(NA_character_, length(depth_range))
  i <- !is.na(depth_range)
  if (!any(i)) return(out)

  n <- 80
  rcap <- pmin(depth_range[i], cap_range)
  idx <- floor(rcap / cap_range * (n - 1)) + 1
  idx <- pmax(1, pmin(n, idx))
  out[i] <- orange_pal(n)[idx]
  out
}

DATA_BY_YEAR_GENE <- DATA_BY_KEY


# --- sample points by dataset (ONE per materialSampleID) ---
SAMPLE_PTS_BY_KEY <- purrr::imap(DATA_BY_YEAR_GENE, ~{
  parts <- strsplit(.y, "_")[[1]]
  gene <- parts[1]
  year <- parts[2]

  .x %>%
    get_depth_num() %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
    distinct(samp_name, decimalLatitude, decimalLongitude, depth_m, .keep_all = TRUE) %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
    mutate(target_gene = gene, year = year)
})

sample_pts_all <- bind_rows(SAMPLE_PTS_BY_KEY) %>%
  mutate(year = as.character(year))

#Compute per-cell depth stats(min,max,median + mixed depth flag)
depth_on_grid_stats <- function(grid_sf, sample_pts_sf, depth_col = "depth_m",
                                fun_med = stats::median,
                                mixed_delta_m = 1.0) {  # <- threshold: min/max must differ by >= 1 m to count as "mixed"
  stopifnot("cell_id" %in% names(grid_sf))
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(inherits(sample_pts_sf, "sf"))

  if (sf::st_crs(sample_pts_sf) != sf::st_crs(grid_sf)) {
    sample_pts_sf <- sf::st_transform(sample_pts_sf, sf::st_crs(grid_sf))
  }
  if (sf::st_crs(all_polys) != sf::st_crs(grid_sf)) {
    all_polys2 <- sf::st_transform(all_polys, sf::st_crs(grid_sf))
  } else {
    all_polys2 <- all_polys
  }

  haspt <- lengths(sf::st_intersects(grid_sf, sample_pts_sf)) > 0
  grid_with <- sf::st_join(grid_sf, sample_pts_sf, join = sf::st_intersects, left = TRUE)

  stats_tbl <- grid_with %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      depth_min = if (all(is.na(.data[[depth_col]]))) NA_real_ else min(.data[[depth_col]], na.rm = TRUE),
      depth_max = if (all(is.na(.data[[depth_col]]))) NA_real_ else max(.data[[depth_col]], na.rm = TRUE),
      depth_med = if (all(is.na(.data[[depth_col]]))) NA_real_ else fun_med(.data[[depth_col]], na.rm = TRUE),
      .groups   = "drop"
    )

  out <- grid_sf %>%
    dplyr::left_join(stats_tbl, by = "cell_id") %>%
    dplyr::mutate(
      has_sampling = haspt,
      depth_range  = depth_max - depth_min,
      mixed_depth  = dplyr::if_else(has_sampling & !is.na(depth_range) & depth_range >= mixed_delta_m, TRUE, FALSE)
    ) %>%
    sf::st_as_sf()

  sf::st_intersection(out, all_polys2) %>%
    sf::st_collection_extract("POLYGON") %>%
    (\(x) x[!sf::st_is_empty(x), ])()
}

#Compute depth grids by year only
# years available in sample points

# use the clipped grid that definitely has cell_id
YEARS_DEPTH <- sort(unique(na.omit(as.character(sample_pts_all$year))))

grid_for_depth <- grid_clip

depth_layers <- c(
  list("All" = depth_on_grid_stats(grid_for_depth, sample_pts_all)),
  setNames(
    purrr::map(YEARS_DEPTH, \(yr) depth_on_grid_stats(grid_for_depth, dplyr::filter(sample_pts_all, year == yr))),
    YEARS_DEPTH
  )
)

depth_layers <- lapply(depth_layers, function(occ_all) {
  occ_all %>%
    dplyr::mutate(
      mixed_depth = as.logical(mixed_depth),
      depth_fill  = depth_to_col(depth_med),
      mixed_fill  = mixed_to_orange(depth_range),
      final_fill  = dplyr::if_else(mixed_depth, mixed_fill, depth_fill, missing = depth_fill),
      final_fill  = dplyr::if_else(!has_sampling | is.na(depth_med), "#7bccc4", final_fill),
      final_alpha = dplyr::if_else(!has_sampling | is.na(depth_med), 0.05, 1)
    )
})


depth_legend_cols <- cool
depth_legend_labs <- depth_labels

#For later
#choices_gene <- sort(unique(KEY_TBL$target_gene))
#choices_year <- sort(unique(KEY_TBL$year))
#choices_key  <- KEY_TBL$key

#layer_sf <- RICHNESS_BY_KEY[[input$key]]
#-------------------------------------------------------------------NEW
standardize_month_col <- function(df) {
  df %>%
    dplyr::mutate(
      month = dplyr::case_when(
        !is.na(suppressWarnings(as.integer(as.character(month)))) &
          suppressWarnings(as.integer(as.character(month))) %in% 1:12 ~
          month.abb[suppressWarnings(as.integer(as.character(month)))],
        as.character(month) %in% month.name ~ substr(as.character(month), 1, 3),
        as.character(month) %in% month.abb ~ as.character(month),
        TRUE ~ as.character(month)
      )
    )
}

make_point_key <- function(df) {
  coords <- sf::st_coordinates(df)

  df %>%
    dplyr::mutate(
      samp_name = as.character(samp_name),
      eventDate = as.character(eventDate),
      lon = coords[, 1],
      lat = coords[, 2],
      lat5 = ifelse(!is.na(lat), sprintf("%.5f", round(lat, 5)), NA_character_),
      lon5 = ifelse(!is.na(lon), sprintf("%.5f", round(lon, 5)), NA_character_),
      point_key = dplyr::case_when(
        !is.na(samp_name) & samp_name != "" &
          !is.na(eventDate) & eventDate != "" &
          !is.na(lat5) & !is.na(lon5) ~
          paste(samp_name, eventDate, lat5, lon5, sep = " || "),

        !is.na(samp_name) & samp_name != "" &
          !is.na(eventDate) & eventDate != "" ~
          paste(samp_name, eventDate, sep = " || "),

        !is.na(samp_name) & samp_name != "" &
          !is.na(lat5) & !is.na(lon5) ~
          paste(samp_name, lat5, lon5, sep = " || "),

        !is.na(samp_name) & samp_name != "" ~
          samp_name,

        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-lon, -lat, -lat5, -lon5)
}

occ_all <- standardize_month_col(occ_all)

species_sf_all <- species_sf_all %>%
  standardize_month_col() %>%
  make_point_key()

# --- Sampling points layer for leaflet (keeps year + source_file if present) ---

sampling_pts <- species_sf_all %>%                       #NEW
  mutate(
    year = if ("year" %in% names(.)) as.character(year) else NA_character_
  ) %>%
  distinct(target_gene, year, samp_name, geometry, .keep_all = TRUE)


species_sf_min <- species_sf_all %>%
  dplyr::select(
    point_key,
    occurrenceID,
    samp_name,
    eventDate,
    scientificName,
    target_gene,
    year,
    kingdom,
    phylum,
    class,
    order,
    geometry
  )

species_sf_by_year <- split(species_sf_min, as.character(species_sf_min$year))
species_sf_by_year[["All"]] <- species_sf_min

point_poly_lookup <- sf::st_join(
  species_sf_all,
  all_polys_click %>% dplyr::select(site_type, site_name),
  join = sf::st_within,
  left = FALSE
) %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(point_key, site_type, site_name)

species_sf_all_with_poly <- species_sf_all %>%
  dplyr::left_join(point_poly_lookup, by = "point_key")

point_cell_lookup <- sf::st_join(
  species_sf_all,
  grid_clip %>% dplyr::select(cell_id),
  join = sf::st_within,
  left = FALSE
) %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(point_key, cell_id)

species_sf_all_with_cell <- species_sf_all %>%
  dplyr::left_join(point_cell_lookup, by = "point_key")

#-------------------------------------------------------------

#Bundle the outputs in one list:
APP_DATA <- list(
  occ_all = occ_all,
  KEY_TBL = KEY_TBL,
  sampling_pts = sampling_pts,
  species_sf_all = species_sf_all,
  species_sf_min = species_sf_min,
  species_sf_by_year = species_sf_by_year,
  point_cell_lookup = point_cell_lookup,
  species_sf_all_with_cell = species_sf_all_with_cell,
  point_poly_lookup = point_poly_lookup,
  species_sf_all_with_poly = species_sf_all_with_poly,
  grid_clip = grid_clip,
  RICHNESS_BY_KEY = RICHNESS_BY_KEY,
  RICHNESS_ALL = RICHNESS_ALL,
  RICHNESS_ALL_BY_YEAR = RICHNESS_ALL_BY_YEAR,
  depth_layers = depth_layers,
  all_polys_click = all_polys_click,
  all_polys_zones = all_polys_zones,
  pal_rich = pal_rich
)

#Undo hashtags when ready to load and save data

APP_DATA <- mget(ls()) #stores all variables from the environment into APP_DATA
saveRDS(APP_DATA, "./inst/app/data/v2_essential_app_data_20260424.rds") #saves that as a file
#APP_DATA <- readRDS("./inst/app/data/v2_essential_app_data_20260424.rds") #loads that file
#list2env(APP_DATA, .GlobalEnv) #pulls everything out of APP_DATA into their original names

ggplot2::theme_set(
  ggplot2::theme_minimal(base_family = "sans")
)
