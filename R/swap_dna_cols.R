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
