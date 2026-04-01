
is_dna_seq <- function(x, min_len = 10) {
  grepl(paste0("^[ACGTURYKMSWBDHVNI]{", min_len, ",}$"),
        x,
        ignore.case = TRUE)
}

fix_swapped_columns <- function(df, col_correct, col_dna) {
  correct_vals <- df[[col_correct]]
  dna_vals     <- df[[col_dna]]

  # Detect DNA-like values
  correct_is_dna <- is_dna_seq(correct_vals)
  dna_is_dna     <- is_dna_seq(dna_vals)

  # Rows where they are swapped:
  # correct column has DNA AND dna column does NOT
  swapped <- correct_is_dna & !dna_is_dna

  # Swap only those rows
  tmp <- correct_vals[swapped]
  correct_vals[swapped] <- dna_vals[swapped]
  dna_vals[swapped]     <- tmp

  # Put back into df
  df[[col_correct]] <- correct_vals
  df[[col_dna]]     <- dna_vals

  return(df)
}

fix_primer_columns <- function(df) {

  df <- fix_swapped_columns(
    df,
    col_correct = "pcr_primer_name_forward",
    col_dna     = "pcr_primer_forward"
  )

  df <- fix_swapped_columns(
    df,
    col_correct = "pcr_primer_name_reverse",
    col_dna     = "pcr_primer_reverse"
  )

  df
}


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


