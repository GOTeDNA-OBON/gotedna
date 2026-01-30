
#Code for removing duplicate data
# n_before <- nrow(data_df)
#
# df_dedup <- data_df[!duplicated(data_df), ]
#
# n_after <- nrow(df_dedup)
#
# n_duplicates <- n_before - n_after
# n_duplicates



#PULLS ALL SUITABLE DATA FROM OBIS
D_mb <- read_data(
  dataset_ids    = NULL,
  scientificname = NULL,
  worms_id       = NULL,
  areaid         = NULL,
  join_by        = c("auto", "occurrenceID", "id"),
  require_absences = TRUE
)

#Reduce false positives (metabarcoding detections < 10 reads/station/primer excluded)
D_mb_msct <- D_mb %>%
  dplyr::mutate(msct = case_when(
    organismQuantity == 0 ~ TRUE,
    organismQuantity > 10 ~ TRUE
  )) |>
  tidyr::drop_na(msct)

D_mb_nodetect <- D_mb_msct %>%
  dplyr::group_by(
    protocol_ID, protocolVersion, scientificName, primer, station) %>%
  dplyr::summarise(num_detected = sum(detected)) %>%
  dplyr::filter(num_detected == 0)

D_mb_clean <- dplyr::anti_join(D_mb_msct, D_mb_nodetect,
                               by = c("protocol_ID","protocolVersion","scientificName",
                                      "primer", "station"))

# make a list of two data frames
gotedna_data <- list(
  metabarcoding = D_mb_clean |>
    dplyr::filter(!is.na(decimalLongitude),
                  !class %in% c("Aves","Insecta", "Hexapoda"),
                  !order %in% c("Primates","Artiodactyla","Perissodactyla","Rodentia"),
                  !family %in% c("Felidae","Canidae","Procyonidae")) |>
    dplyr::ungroup() |>
    as.data.frame(),
  qPCR = D_qPCR |>
    dplyr::filter(!is.na(decimalLongitude), !is.na(phylum)) |>
    dplyr::ungroup() |>
    as.data.frame()
)

# Data
## import glossary
gloss <- read.csv("./inst/app/data/glossary.csv")
gloss$Term <- paste0('<p align ="right"><b>', trimws(gloss$Term), "</b></p>")
gloss$Definition <- trimws(gloss$Definition)

## import GOTeDNA data
gotedna_data <- gotedna_data0 <- readRDS("./inst/app/data/gotedna_data.rds")

gotedna_data$metabarcoding <- readRDS("./inst/app/data/test_obis_animalia.rds")



writeLines(
  format(round(Sys.time(), "mins"), "%Y-%m-%d %H:%M %Z"),
  "inst/app/data/last_obis_download_ts.txt"
)

# for performances sake, we use a separate object for station to only display
# a few points on map


gotedna_station <- list(
  metabarcoding = gotedna_data$metabarcoding |> get_station(),
  qPCR = gotedna_data$qPCR |> get_station()
)
saveRDS(gotedna_station, "inst/app/data/gotedna_station.rds")

###################
#Check Primers for typos
###################
library(stringdist)
library(dplyr)
library(stringr)

valid_primer_df <- read.csv("inst/app/data/primers.csv")
valid_primer_df <- valid_primer_df %>% filter(Data == "Metabarcoding")

valid_primer_df <- valid_primer_df %>%
  mutate(
    # split OriginalName on <br>
    primers = str_split(OriginalName, "\\s*<br>\\s*"),

    forward_primer = map_chr(primers, 1),
    reverse_primer = map_chr(primers, 2),

    primer_display_name = paste0(
      Locus, " | ",
      forward_primer, " / ", reverse_primer
    )
  ) %>%
  select(-primers)

valid_primers <- valid_primer_df$primer_display_name

map_primers <- function(
    observed,
    valid,
    method = "osa",
    max_dist = 2
) {

  # distance matrix: rows = observed, cols = valid
  dist_mat <- stringdistmatrix(
    observed,
    valid,
    method = method
  )

  best_match_idx  <- apply(dist_mat, 1, which.min)
  best_match_dist <- apply(dist_mat, 1, min)

  tibble(
    observed_primer = observed,
    matched_primer  = valid[best_match_idx],
    distance        = best_match_dist,
    status = case_when(
      best_match_dist == 0              ~ "exact",
      best_match_dist <= max_dist       ~ "corrected",
      TRUE                              ~ "new_or_unmatched"
    ),
    final_primer = if_else(
      best_match_dist <= max_dist,
      valid[best_match_idx],
      observed
    )
  )
}

primer_map <- map_primers(
  observed = unique(gotedna_data$metabarcoding$primer),
  valid    = valid_primers,
  max_dist = 10
)

View(primer_map)

gotedna_data$metabarcoding <- gotedna_data$metabarcoding %>%
  left_join(
    primer_map %>% select(observed_primer, final_primer),
    by = c("primer" = "observed_primer")
  ) %>%
  mutate(primer = coalesce(final_primer, primer)) %>%
  select(-final_primer)

saveRDS(gotedna_data, "inst/app/data/gotedna_data.rds")

#### Need to do for both qPCR and metabarcoding
# Prepare primer data
# gotedna_data <- readRDS("inst/app/data/gotedna_data.rds")

newprob_mb <- calc_det_prob(gotedna_data$metabarcoding)
newprob_q <- calc_det_prob(gotedna_data$qPCR)

scaledprobs_mb <- scale_newprob(gotedna_data$metabarcoding, newprob_mb)
scaledprobs_q <- scale_newprob(gotedna_data$qPCR, newprob_q)
scaledprobs_q <- NULL

gotedna_primer <- list()

# this needs to be based on the area selection
for (i in c("kingdom", "phylum", "class", "order", "family", "genus", "scientificName")) {
  gotedna_primer[[i]] <- primer_sort(i, dplyr::bind_rows(scaledprobs_mb, scaledprobs_q)) |>
    mutate(text = paste0(primer, " (", detects, "/", total, " ", perc, "%)"))
}

saveRDS(gotedna_primer, "inst/app/data/gotedna_primer.rds")
