#' Calculate detection probabilities from data that were previously imported
#' with read_data()
#'
#' @description Function that calculates non-parametric probability of detection
#' for each selected taxon and primer for the selected data. Probabilities are
#' calculated both (1) monthly, across all years; and (2) monthly with each year
#' separate. Outputs are used in subsequent functions.
#'
#' @param data (required, data.frame) Data.frame read in with [read_data()].
#'
#' @return Two lists, each with distinct elements of the selected taxon; primer ID
#' containing 5-7 columns:
#' * `month`
#' * `n` total number of samples per month
#' * `nd` number of positive detections per month
#' * `p` calculated monthly detection probability
#' * `s` standard deviation
#' * `protocol_ID`
#' * `year`
#' * `yr.mo` year;month concatenated
#'
#' @author Anais Lacoursiere-Roussel \email{Anais.Lacoursiere@@dfo-mpo.gc.ca}
#' @rdname calc_det_prob
#' @export
#' @examples
#' \dontrun{
#'  calc_det_prob(data = D_mb)
#' }
calc_det_prob <- function(data, selected_taxon_level = "scientificName", selected_taxon_id = "All", pool_primers = FALSE) {

  oop <- options("dplyr.summarise.inform")
  options(dplyr.summarise.inform = FALSE)
  # reset option on exit
  on.exit(options(dplyr.summarise.inform = oop))


  print(paste0("rows in data before nondetection distance filter: ", nrow(data)))
  #removing rows for species that have never been detected at that station
  data <- filter_nondetections_all(
       data,
       distance = 500,
       selected_taxon_level,
       selected_taxon_id
    )
  print(paste0("rows in data after nondetection distance filter: ", nrow(data)))



  data <- data %>%
    dplyr::group_by(
      protocol_ID,
      primer,
      month,
      year,
      .data[[selected_taxon_level]],
      samp_name
    ) %>%
    dplyr::summarise(
      detected = as.integer(any(detected == 1)),
      .groups = "drop"
    )

  bad <- !data$detected %in% c(0, 1)

  if (any(bad)) {
    stop(
      "calc_det_prob(): detected must be 0/1 after collapsing.\n",
      "Found ", sum(bad), " invalid values.\n",
      "Unique invalid values: ",
      paste(unique(data$detected[bad]), collapse = ", ")
    )
  }

  if (pool_primers) {
  data %<>%
    dplyr::mutate(.,
                  id = paste0(protocol_ID, ";", .data[[selected_taxon_level]]),
                  id.yr = paste0(protocol_ID, ";", .data[[selected_taxon_level]], ";ALLPRIMERS;", year)
    )
  } else {
  data %<>%
    dplyr::mutate(.,
      id = paste0(protocol_ID, ";", .data[[selected_taxon_level]], ";", primer),
      id.yr = paste0(protocol_ID, ";", .data[[selected_taxon_level]], ";", primer, ";", year)
    )
  }

  # Create a variable so detection probability is calculated separately for each
  # protocol ID, version, selected_taxon_level, and primer

  # create new list variables to store outputs
  lnd <- length(unique(data$id))
  newP <- vector("list", lnd)
  names(newP) <- unique(data$id)
  SUM <- COM <- comps <- vector("list", length(unique(data$id)))
  names(COM) <- unique(data$id)

  # calculate detection probabilities - year aggregated
  for (occurrence in unique(data$id)) {
    SDF <- data %>%
      dplyr::filter(id == occurrence) %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(
        n = dplyr::n(),
        nd = sum(detected),
        p = nd/n,
        s = sqrt(p * (1 - p)/n)
      ) %>%
      as.data.frame()

    if (any(SDF$n > 1)) {
      newP[[occurrence]] <- SDF
    }
  }
  newP_agg <- newP[lengths(newP) != 0]


  # calculate monthly detection probability for each year
  lny <- length(unique(data$id.yr))
  newP <- vector("list", lny)
  names(newP) <- unique(data$id.yr)
  SUM <- COM <- comps <- vector("list", length(unique(data$id.yr)))
  names(COM) <- unique(data$id.yr)

  for (occurrence in unique(data$id.yr)) {
    SDF <- data %>%
      dplyr::filter(id.yr == occurrence) %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(
        n = dplyr::n(),
        nd = sum(detected),
        p = nd/n,
        s = sqrt(p * (1 - p)/n)
      ) %>%
      as.data.frame()

    if (any(SDF$n > 1)) {
      newP[[occurrence]] <- SDF
    }
  }
  newP <- newP[lengths(newP) != 0]
  newP_yr <- Map(cbind, newP, id.yr = names(newP))

  newP_yr <- lapply(newP_yr, function(x) {
    dplyr::mutate(x,
      protocol_ID = stringr::word(x$id.yr, 1, sep = stringr::fixed(";")),
      year = stringr::word(x$id.yr, -1, sep = stringr::fixed(";")),
      yr.mo = paste0(year, ";", month),
      id.yr = NULL
    )
  }) # to obtain variation among years

  list(newP_agg = newP_agg, newP_yr = newP_yr)
}

drop_all_zero_taxa <- function(df, taxon_col) {
  df %>%
    group_by(
      station,
      protocol_ID,
      .data[[taxon_col]]
    ) %>%
    filter(any(detected == 1)) %>%
    ungroup()
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
