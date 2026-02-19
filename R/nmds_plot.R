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








protocol_nmds <- function(df) {
  # Make a copy for labeling
  df_labels <- df

  # Columns to exclude from NMDS (just the label)
  exclude_cols <- "protocol_ID"

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

  # Plot
  ggplot(nmds_points, aes(NMDS1, NMDS2, label = protocol_id)) +
    geom_point(size = 4, color = "#440154") +
    geom_text(vjust = -0.7) +
    theme_minimal() +
    labs(
      title = "NMDS of Protocols (Gower distance)",
      subtitle = paste("Stress =", stress_val)
    )
}


# protocol_nmds(protocol_test_sheet)

