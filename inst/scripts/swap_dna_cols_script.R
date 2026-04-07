#the following ids are known to be entered incorrectly and should have dna columns switched if using
#read_raw_data() and prepare_raw_data() to prepare the data for the app.
#There may also be other files that need swapping, so feel free to run
#the code below on all existing_files rather than just these ids.
#Conversely, if you are using the read_data() function to prepare data,
#it already swaps columns where necessary for each dataset and returns one bulk dataframe.

ids <- c(
  "5b6b4895-8b97-4fe0-a8f9-e085cb849fba",
  "5fb6a5e8-a2c4-4a75-aff9-3b4ebb2cdc74",
  "88f84d06-da70-486c-8944-ac1601977551",
  "e5337db0-95bd-4b1f-bddd-6202622d5402"
)

# Full file paths
existing_files <- list.files(
  "inst/app/data/raw_OBIS",
  pattern = "^dataset-.*\\.rds$",
  full.names = TRUE
)

pattern_ids <- paste(ids, collapse = "|")

filtered_files <- existing_files[
  grepl(pattern_ids, basename(existing_files))
]

lapply(filtered_files, function(f) {
  df <- readRDS(f)

  before <- sum(is_dna_seq(df$pcr_primer_name_forward), na.rm = TRUE)

  df <- fix_primer_columns(df)

  after <- sum(is_dna_seq(df$pcr_primer_name_forward), na.rm = TRUE)

  message(basename(f), ": swapped rows = ", before - after)

  saveRDS(df, f)
})


