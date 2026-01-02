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
field_sample_fig <- function(data) {


  ## ---- aggregation (always happens, once) ----
  data <- data %>%
    dplyr::group_by(scientificName, month) %>%
    dplyr::summarise(
      `Sample size` = dplyr::n(),
      nd = sum(detected, na.rm = TRUE),
      `Detection rate` = nd / `Sample size`,
      .groups = "drop"
    )

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

  ## ---- aesthetics prep ----
  species_levels <- sort(unique(data$scientificName))
  alphabet_colors <- RColorBrewer::brewer.pal(8, "Set1")
  colors <- rep(alphabet_colors, length.out = length(species_levels))
  names(colors) <- species_levels

  max_size <- max(data$`Sample size`, na.rm = TRUE)
  sizeref <- 40 * max_size / 100^2

  data$month_jittered <- data$month + runif(nrow(data), -0.25, 0.25)

  ## ---- plot ----
  plotly::plot_ly(
    data,
    x = ~month_jittered,
    y = ~`Detection rate`,
    type = "scatter",
    mode = "markers",
    color = ~scientificName,
    colors = colors,
    size = ~`Sample size`,
    marker = list(sizemode = "area", sizeref = sizeref, sizemin = 4),
    text = ~paste(
      "Species:", scientificName,
      "<br>Month:", month,
      "<br>Detection Rate:", round(`Detection rate`, 2),
      "<br>Sample size:", `Sample size`
    ),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      xaxis = list(
        title = list(text = "Month", font = list(size = 18)),
        tickmode = "array",
        tickvals = 1:12,
        ticktext = month.abb,
        range = c(0.5, 12.5)
      ),
      yaxis = list(
        title = list(text = "Detection rate", font = list(size = 18)),
        range = c(0, 1)
      ),
      legend = list(
        orientation = "v",
        x = 1.02, y = 1
      ),
      margin = list(t = 60, b = 50, r = 150)
    )
}



