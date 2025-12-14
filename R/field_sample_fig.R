#' Display actual sampling effort of taxa within specified taxonomic group.
#'
#' @description This function displays the actual number of eDNA samples taken
#' for all species of a specified taxonomic group per year.
#'
#' @param data (required, data.frame): Data.frame read in with [read_data()]
#' @author Anais Lacoursiere-Roussel \email{Anais.Lacoursiere@@dfo-mpo.gc.ca}
#' @rdname field_sample_fig
#' @export
#' @examples
#' \dontrun{
#' field_sample_fig(
#'   data = gotedna_data$metabarcoding
#' )
#' }
prepare_plotly_samples <- function(data) {
  # Only run if 'detected' exists; otherwise assume already aggregated
  if (!"detected" %in% names(data)) {
    message("Data already aggregated, skipping prepare_plotly_samples()")
    return(data)
  }

  stopifnot("detected" %in% names(data))
  data %>%
    dplyr::group_by(scientificName, month) %>%
    dplyr::summarise(
      n = dplyr::n(),
      nd = sum(detected, na.rm = TRUE),
      freq_det = nd / n,
      .groups = "drop"
    ) %>%
    dplyr::rename(
      `Detection rate` = freq_det,
      `Sample size` = n
    )
}

field_sample_fig <- function(data) {
  cat("\n=== field_sample_fig ===\n")
  cat("Incoming rows:", nrow(data), "\n")
  cat("Incoming cols:", paste(names(data), collapse = ", "), "\n")

  data <- prepare_plotly_samples(data)

  cat("After prepare_plotly_samples()\n")
  cat("Rows:", nrow(data), "\n")
  cat("Cols:", paste(names(data), collapse = ", "), "\n")

  if (nrow(data) == 0) {
    cat("⚠️ ZERO ROWS AFTER SUMMARISE\n")
    return(
      plotly::plot_ly(
        type = "scatter",
        mode = "text",
        text = "No data available for this selection",
        x = 0, y = 0
      )
    )
  }

  # Color palette matching palette("Alphabet")
  species_levels <- sort(unique(data$scientificName))
  alphabet_colors <- RColorBrewer::brewer.pal(8, "Set1")
  colors <- rep(alphabet_colors, length.out = length(species_levels))
  names(colors) <- species_levels

  max_size <- max(data$`Sample size`, na.rm = TRUE)
  sizeref <- 40 * max_size / 100^2

  # Optional: x-jitter
  data$month_jittered <- data$month + runif(nrow(data), -0.25, 0.25)

  plotly::plot_ly(
    data,
    x = ~month_jittered,
    y = ~`Detection rate`,
    type = "scatter",
    mode = "markers",
    color = ~scientificName,
    colors = colors,
    size = ~`Sample size`,
    marker = list(sizemode = 'area', sizeref = sizeref, sizemin = 4),
    text = ~paste("Species:", scientificName, "<br>Month:", month, "<br>Sample size:", `Sample size`),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      xaxis = list(
        title = list(
          text = "Month",
          font = list(size = 18, color = "black")
        ),
        tickmode = "array",
        tickvals = 1:12,
        tickangle = -45,
        ticktext = month.abb,
        tickfont = list(size = 12, color = "black"),
        showticklabels = TRUE,
        range = c(0.5, 12.5),
        ticks = "outside",
        tickwidth = 1,
        tickcolor = "black"
      ),
      yaxis = list(
        title = list(
          text = "Detection rate",
          font = list(size = 18, color = "black")
        ),
        range = c(0, 1),
        tickfont = list(size = 16, color = "#939598")
      ),
      showlegend = TRUE,
      legend = list(
        font = list(size = 14, color = "#939598"),
        orientation = "v",
        x = 1.02, y = 1
      ),
      margin = list(t = 60, b = 50, r = 150)
    )
}


