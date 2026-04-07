#' Remove duplicate rows and redundant non-detections
#'
#' @description
#' Cleans a detection dataset by:
#' \enumerate{
#'   \item Removing exact duplicate rows.
#'   \item Removing non-detection records (detected == 0) when a detection
#'         (detected == 1) exists for the same combination of
#'         \code{samp_name}, \code{primer}, \code{protocol_ID},
#'         \code{protocolVersion}, and \code{scientificName}.
#'   \item Removing species that have no positive detections anywhere in
#'         the dataset.
#'   \item Removing species/primer/protocol combinations that contain no
#'         positive detections.
#' }
#'
#' This ensures that retained records represent meaningful detections
#' and eliminates redundant or entirely negative combinations.
#'
#' @param df A data.frame or tibble containing at minimum the columns:
#'   \code{samp_name}, \code{primer}, \code{protocol_ID},
#'   \code{protocolVersion}, \code{scientificName}, and
#'   \code{detected} (0/1).
#'
#' @return A cleaned tibble with duplicates and redundant
#'   non-detections removed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cleaned_data <- remove_duplicates_and_undetected(gotedna_data$metabarcoding)
#' }
remove_duplicates_and_undetected <- function(df) {

  #remove duplicate rows
  df <- dplyr::distinct(df, samp_name, primer, protocol_ID, scientificName, detected, .keep_all = TRUE)

  #remove nondetections if there is a detection for that species, primer, protocol_ID, and samp_name
  df <- df %>%
    group_by(samp_name, primer, protocol_ID, scientificName) %>%
    filter(
      !(any(detected == 1, na.rm = TRUE) &
          any(detected == 0, na.rm = TRUE) &
          detected == 0)
    ) %>%
    ungroup()

  #First remove species and primer/protocol combinations that were not detected anywhere
  df <- df %>%
    # Step 1: remove species with no positive detections anywhere
    group_by(scientificName) %>%
    filter(any(detected == 1, na.rm = TRUE)) %>%
    ungroup() %>%

    # Step 2: remove species/primer/protocol combos with no positives
    group_by(scientificName, primer, protocol_ID) %>%
    filter(any(detected == 1, na.rm = TRUE)) %>%
    ungroup()
  df
}

