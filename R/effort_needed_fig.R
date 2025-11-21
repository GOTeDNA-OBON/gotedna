#' Calculate sampling effort needed to obtain different detection thresholds.
#'
#' @description This function calculates number of samples needed to obtain
#' the selected taxon detection at different thresholds by using scaled and interpolated
#' data produced with `[scale_newprob()]`.
#'
#' @param scaledprobs  (required, data.frame) Normalized detection
#' probabilities as returned by the element `month` of the list returned by
#' [scale_newprob()] for one taxon and one primer.
#'
#' @author Anais Lacoursiere-Roussel \email{Anais.Lacoursiere@@dfo-mpo.gc.ca}
#' @rdname effort_needed_fig
#' @export
#' @examples
#' \dontrun{
#' newprob <- calc_det_prob(gotedna_data$metabarcoding)
#' scaledprobs <- scale_newprob(gotedna_data$metabarcoding, newprob)
#' effort_needed_fig(
#'   scaledprobs
#'   )
#' )
#' }

effort_needed_fig <- function(scaledprobs, height_per_species = 400, dot_size = 12, selected_taxon_level = "species") {
  # Filter NA years and prepare month labels
  df <- scaledprobs %>%
    dplyr::filter(is.na(year)) %>%
    dplyr::mutate(
      month = factor(
        month,
        levels = 1:12,
        labels = c("Jan","Feb","Mar","Apr","May",
                   "Jun","Jul","Aug","Sep","Oct","Nov","Dec")
      )
    )

  taxa_list <- unique(df[[selected_taxon_level]])
  n_taxa <- length(taxa_list)

  # Make a cyclic month color palette (so Dec ~ Jan)
  cyclic_palette <- colorRampPalette(
    c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725", "#7b3f53")
  )(12)

  # Generate plots for each species
  plots <- lapply(seq_along(taxa_list), function(i) {
    sp <- taxa_list[i]
    df_sp <- df[df[[selected_taxon_level]] == sp, ]

    DF_sp <- expand.grid(
      p = df_sp$fill,
      `Samples needed` = seq_len(25),
      `Detection rate` = NA
    ) %>%
      merge(data.frame(
        p = df_sp$fill,
        Month = df_sp$month,
        Taxon = df_sp[[selected_taxon_level]]
      ))

    for (j in seq_len(nrow(DF_sp))) {
      DF_sp$`Detection rate`[j] <- 1 - dbinom(
        0,
        size = DF_sp$`Samples needed`[j],
        prob = DF_sp$p[j]
      )
    }

    show_legend <- i == 1

    plotly::plot_ly(
      DF_sp,
      x = ~`Samples needed`,
      y = ~`Detection rate`,
      type = "scatter",
      mode = "markers",
      color = ~Month,
      colors = cyclic_palette,
      marker = list(size = dot_size),
      showlegend = show_legend
    ) %>%
      plotly::layout(
        yaxis = list(range = c(0, 1), title = "Detection probability"),
        xaxis = list(title = "Number of samples")
      )
  })

  # Combine all subplots vertically
  subplot_obj <- plotly::subplot(
    plots,
    nrows = n_taxa,
    shareX = TRUE,
    titleY = TRUE,
    heights = rep(1 / n_taxa, n_taxa),
    margin = 0.01
  )

  # Add per-taxon titles via annotations
  heights <- rep(1 / n_taxa, n_taxa)
  padding <- 0.02  # distance from bottom of subplot

  annotations <- lapply(seq_along(taxa_list), function(i) {
    y_bottom <- 1 - sum(heights[1:i])
    subplot_height <- heights[i]
    list(
      text = taxa_list[i],
      x = 0.98,                        # right edge
      y = y_bottom + subplot_height * padding,
      xref = "paper",
      yref = "paper",
      xanchor = "right",
      yanchor = "bottom",
      showarrow = FALSE,
      font = list(size = 22, color = "#333")
    )
  })

  subplot_obj %>%
    plotly::layout(
      height = n_taxa * height_per_species,
      margin = list(t = 80, r = 180, b = 60),
      legend = list(
        title = list(text = "Month", font = list(size = 16)),  # legend title size
        font = list(size = 12),                                 # legend text size
        orientation = "v",
        y = 1,          # top
        x = 1.02,       # right
        xanchor = "left",
        yanchor = "top"
      ),
      annotations = annotations
    )
}

