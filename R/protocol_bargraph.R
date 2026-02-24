protocol_bargraph <- function(df) {
  color_vec <- c("#00A08A", "#446455", "#Fdd262", "#5BBCD6", "#046c9a", "#ABDDDE", "#d3dddc")

  # assign a color cycling through bars
  df <- df %>%
    dplyr::mutate(
      color = color_vec[(seq_len(n()) - 1) %% length(color_vec) + 1],
      protocol_ID = factor(protocol_ID, levels = sort(unique(protocol_ID)))  # ensure left-to-right numeric order
    )

  plot_ly(
    data = df,
    x = ~protocol_ID,
    y = ~detection_rate,
    type = 'bar',
    marker = list(color = df$color),
    hoverinfo = "text",
    name = "Detection Rate"
  ) %>%
    layout(
      xaxis = list(title = "Protocol ID"),
      yaxis = list(title = "Detection Rate (%)", showgrid = FALSE),
      margin = list(l = 60, r = 20, t = 50, b = 60),
      legend = list(title = list(text = "")),
      shapes = list(
        list(
          type = "line",
          x0 = -0.5, x1 = -0.5,
          y0 = 0, y1 = 100,  # detection rate is 0–100%
          line = list(color = "black", width = 1)
        )
      )
    ) %>%
    config(
      displayModeBar = TRUE,
      modeBarButtonsToAdd = c("resetScale2d")
    )
}


