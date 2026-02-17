#!/usr/bin/env Rscript

# Combine all profile run CSVs into a single file, keeping metadata about
# the profiled variable, mode, and original source file.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tools)
})

agg_dir <- "results/aggregated"
pattern <- "^combined_results_.*\\.csv$"
csv_files <- list.files(agg_dir, pattern = pattern, full.names = TRUE)

extract_profile_meta <- function(path) {
  base <- file_path_sans_ext(basename(path))
  mode_match <- regexpr("_(fit|refined|refined_second)$", base, perl = TRUE)
  if (mode_match == -1) {
    return(NULL)
  }
  mode <- sub("^_", "", regmatches(base, mode_match))
  prefix <- substr(base, 1, mode_match - 1)
  prefix <- sub("^combined_results_", "", prefix)
  var_part <- sub("^[^_]+_[^_]+_", "", prefix)
  if (identical(var_part, prefix) || nchar(var_part) == 0) {
    return(NULL) # no third token -> not a profile result
  }
  list(profile_var = var_part, mode = mode)
}

profiles <- lapply(csv_files, function(path) {
  meta <- extract_profile_meta(path)
  if (is.null(meta)) {
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE) |>
    mutate(
      profile_var = meta$profile_var,
      mode = meta$mode,
      source_file = basename(path)
    )
})

profiles <- Filter(Negate(is.null), profiles)

if (length(profiles) == 0) {
  stop("No profile CSVs found in ", agg_dir)
}

combined_profiles <- bind_rows(profiles)
out_file <- file.path(agg_dir, "all_profiles.csv")
write_csv(combined_profiles, out_file)
cat("Wrote", nrow(combined_profiles), "rows to", out_file, "\n")
