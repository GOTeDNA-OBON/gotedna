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

  # Create hover text for each point
  nmds_points <- nmds_points %>%
    mutate(
      hover_text = paste0(
        "Protocol: ", protocol_id, "<br>",
        "NMDS1: ", round(NMDS1, 3), "<br>",
        "NMDS2: ", round(NMDS2, 3), "<br>",
        "Stress: ", round(stress_val, 3)
      )
    )

  # Create interactive plot
  plot_ly(
    data = nmds_points,
    x = ~NMDS1,
    y = ~NMDS2,
    type = 'scatter',
    mode = 'markers+text',
    text = ~protocol_id,
    textposition = 'top center',
    hoverinfo = 'text',
    hovertext = ~hover_text,
    marker = list(
      size = 8,
      color = '#2241a7',
      line = list(width = 1, color = 'black')
    )
  ) %>%
    layout(
      title = list(
        text = paste0("Stress = ", round(stress_val, 3)),
        x = 0.05,
        xanchor = "left"
        ),
      xaxis = list(title = "NMDS1", zeroline = FALSE, showgrid = TRUE, automargin = TRUE),
      yaxis = list(title = "NMDS2", zeroline = FALSE, showgrid = TRUE, automargin = TRUE),
      margin = list(l = 60, r = 20, t = 50, b = 60)  # padding around plot
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

