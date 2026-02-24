library(forcats)


set.seed(123)  # for reproducibility

# Mock protocol sheet
protocol_test_sheet <- tibble::tibble(
  protocol_ID = 1:10,
  nucl_acid_ext_kit = c(
    rep("DNeasy® Blood & Tissue Kit (Qiagen)", 5),
    rep("Zymo Quick-DNA Kit", 5)
  ),
  platform = NA_character_,  # all NA to test missing column handling
  instrument = sample(c("MiSeq", "NextSeq", NA), 10, replace = TRUE),
  seq_kit = c(
    rep("MiSeq Reagent Kit v3 (Illumina)", 6),
    rep("NextSeq 2000 Kit", 4)
  ),
  otu_db = c(
    rep("barque v1.7.2 (Mathon et al. 2021)", 7),
    rep("custom_db v2.0", 3)
  ),
  tax_assign_cat = sample(c("SINTAX", "RDP", NA), 10, replace = TRUE),
  otu_seq_comp_appr = sample(c(
    "blastn (NCBI)", "vsearch", "usearch", NA
  ), 10, replace = TRUE),
  minimumDepthInMeters = sample(1:5, 10, replace = TRUE),
  maximumDepthInMeters = sample(2:6, 10, replace = TRUE)
)


plot_nmds_interactive <- function(nmds_points, stress_val) {

  nmds_points <- nmds_points %>%
    mutate(protocol_id = as.character(protocol_id))

  shape_vec <- c("circle", "square", "diamond", "triangle-up", "triangle-down")
  color_vec <- c("#00A08A", "#446455", "#Fdd262", "#5BBCD6", "#046c9a", "#ABDDDE", "#d3dddc")

  p <- plotly::plot_ly()
  prots_present <- sort(unique(nmds_points$protocol_id))

  for (i in seq_along(prots_present)) {
    prot <- prots_present[i]
    df_sub <- nmds_points[nmds_points$protocol_id == prot, ]
    if (nrow(df_sub) == 0) next

    shape_idx <- ((i - 1) %% length(shape_vec)) + 1
    color_idx <- ((i - 1) %% length(color_vec)) + 1
    fill_col   <- color_vec[color_idx]
    border_col <- grDevices::adjustcolor(fill_col, offset = c(-0.2, -0.2, -0.2, 0))

    p <- p %>%
      plotly::add_trace(
        data = df_sub,
        x = ~NMDS1,
        y = ~NMDS2,
        type = "scatter",
        mode = "markers+text",
        text = ~protocol_id,
        textposition = "top center",
        name = paste("Protocol", prot),
        marker = list(
          size = 18,
          symbol = shape_vec[shape_idx],
          color = fill_col,
          opacity = 0.9,
          line = list(color = border_col, width = 2)
        ),
        showlegend = TRUE
      )
  }

  p %>%
    plotly::layout(
      xaxis = list(showgrid = FALSE, showline = TRUE, linecolor = "black"),
      yaxis = list(showgrid = FALSE, showline = TRUE, linecolor = "black"),
      legend = list(title = list(text = "Protocol")),
      margin = list(l = 60, r = 140, t = 50, b = 60),
      annotations = list(
        list(
          x = 0, y = 1.01,
          xref = "paper", yref = "paper",
          text = paste0("Stress = ", round(stress_val, 3)),
          showarrow = FALSE,
          xanchor = "left",
          yanchor = "top",
          font = list(size = 14)
        )
      )
    )
}



protocol_nmds <- function(df) {
  # Make a copy for labeling
  df_labels <- df

  # Columns to exclude from NMDS (just the label)
  exclude_cols <- c("protocol_ID", "protocolVersion")

  # Fill NAs and convert characters to factors for NMDS
  df_nmds <- df %>%
    dplyr::select(-all_of(exclude_cols)) %>%  # remove protocol_ID
    mutate(across(where(is.character),
                  ~ forcats::fct_na_value_to_level(as.factor(.), level = "Missing"))) %>%
    mutate(across(where(is.character), as.factor))

  # Compute Gower distance
  gower_dist <- cluster::daisy(df_nmds, metric = "gower")

  # Run NMDS
  nmds_res <- vegan::metaMDS(gower_dist, k = 2, trymax = 100)

  # Extract NMDS coordinates
  nmds_points <- as.data.frame(vegan::scores(nmds_res))

  # Add labels from original dataframe
  if ("protocol_ID" %in% colnames(df_labels)) {
    nmds_points$protocol_id <- df_labels$protocol_ID
  } else {
    nmds_points$protocol_id <- seq_len(nrow(df_labels))
  }

  # Stress
  stress_val <- round(nmds_res$stress, 3)

  plot_nmds_interactive(nmds_points, stress_val)

}



# protocol_nmds(protocol_test_sheet)

