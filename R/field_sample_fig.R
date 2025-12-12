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
  data %>%
    dplyr::group_by(scientificName, month) %>%  # remove year grouping
    dplyr::summarise(
      n = dplyr::n(),
      nd = sum(detected, na.rm = TRUE),
      freq_det = nd / n,
      .groups = "drop"
    ) %>%
    dplyr::rename(
      Month = month,
      `Detection rate` = freq_det,
      `Sample size` = n
    )
}

field_sample_fig <- function(data) {
  data <- prepare_plotly_samples(data)

  # Color palette matching palette("Alphabet")
  species_levels <- sort(unique(data$scientificName))
  alphabet_colors <- RColorBrewer::brewer.pal(8, "Set1")
  colors <- rep(alphabet_colors, length.out = length(species_levels))
  names(colors) <- species_levels

  max_size <- max(data$`Sample size`, na.rm = TRUE)
  sizeref <- 40 * max_size / 100^2

  # Optional: x-jitter
  data$Month_jittered <- data$Month + runif(nrow(data), -0.25, 0.25)

  plotly::plot_ly(
    data,
    x = ~Month_jittered,
    y = ~`Detection rate`,
    type = "scatter",
    mode = "markers",
    color = ~scientificName,
    colors = colors,
    size = ~`Sample size`,
    marker = list(sizemode = 'area', sizeref = sizeref, sizemin = 4),
    text = ~paste("Species:", scientificName, "<br>Month:", Month, "<br>Sample size:", `Sample size`),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      xaxis = list(
        title = list(
          text = "Month",
          font = list(size = 18, color = "black")  # <-- increase size here
        ),
        tickmode = "array",
        tickvals = 1:12,
        tickangle = -45,
        ticktext = month.abb,
        tickfont = list(size = 12, color = "black"),
        showticklabels = TRUE,
        range = c(0.5, 12.5),
        ticks = "outside",    # show small ticks outside
        tickwidth = 1,        # line thickness
        tickcolor = "black" # line color
      ),
      yaxis = list(
        title = list(
          text = "Detection rate",
          font = list(size = 18, color = "black")  # increase y-axis title size
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

