protocol_bargraph <- function(df) {
  color_vec <- c("#00A08A", "#446455", "#Fdd262", "#5BBCD6", "#046c9a", "#ABDDDE")

  # calculate non-detections
  df <- df %>%
    dplyr::mutate(
      non_detections = total_samples - total_detections,
      # assign a color from color_vec cycling through rows
      color = color_vec[(seq_len(n()) - 1) %% length(color_vec) + 1]
    )

  plot_ly() %>%
    # bottom: detections with cycling colors
    add_trace(
      data = df,
      x = ~factor(protocol_ID),
      y = ~total_detections,
      type = 'bar',
      name = "Detections",
      marker = list(color = df$color),  # <-- pass vector directly
      hoverinfo = "y+name"
    ) %>%
    # top: non-detections (keep gray)
    add_trace(
      data = df,
      x = ~factor(protocol_ID),
      y = ~non_detections,
      type = 'bar',
      name = "Non-Detections",
      marker = list(color = '#d3dddc'),
      hoverinfo = "y+name"
    ) %>%
    layout(
      barmode = 'stack',
      xaxis = list(title = "Protocol ID"),
      yaxis = list(title = "Total Samples", showgrid = FALSE),
      margin = list(l = 60, r = 20, t = 50, b = 60),
      legend = list(title = list(text = "Status")),
      shapes = list(
        list(
          type = "line",
          x0 = -0.5, x1 = -0.5,       # just before the first bar
          y0 = 0, y1 = max(df$total_samples), # full height of y-axis
          line = list(color = "black", width = 1)
        )
      )
    ) %>%
    config(
      displayModeBar = TRUE,
      modeBarButtonsToAdd = c("resetScale2d")
    )
}


