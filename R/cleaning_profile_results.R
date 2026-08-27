#' Combine per-parameter profile-likelihood CSVs into one cleaned data frame
#'
#' Reads a set of per-parameter profile-likelihood result files (one CSV per
#' profiled parameter, all sharing a common `prefix`/`sufix` naming scheme),
#' drops rows `readr::read_csv()` flagged as parsing problems, coerces every
#' column to numeric, and keeps only fits above the `max_loglik` threshold.
#' The set of profiled parameter names is inferred from the columns of the
#' first file (for `"b1"`), so all sibling files are expected to share the
#' same non-parameter columns (`loglik`, `loglik.se`, `sample`, `flag`,
#' `delta`).
#'
#' @param prefix Character scalar. Path prefix shared by all per-parameter
#'   CSV files, e.g. `"aggregated_results/profile_Ahmedabad_"`.
#' @param sufix Character scalar. Path suffix shared by all per-parameter
#'   CSV files, e.g. `"_1997_2011.csv"`.
#' @param max_loglik Numeric scalar. Rows with `loglik >= max_loglik` are
#'   dropped (default `-500`).
#'
#' @return A combined data frame with one row per accepted parameter set
#'   across all profiled parameters.
#' @export
clean_profile <- function(prefix, sufix, max_loglik = -500) {
  name <- "b1"
  csv_name <- paste0(prefix, name, sufix)
  result <- read_csv(csv_name, show_col_types = FALSE)

  names_param <- result |>
    dplyr::select(-any_of(c("loglik", "loglik.se", "sample", "flag", "delta"))) |>
    colnames()

  all_res <- c()
  for (name in names_param) {
    csv_name <- paste0(prefix, name, sufix)
    result <- read_csv(csv_name, show_col_types = FALSE)
    result <- result |> filter(!(row_number() %in% (problems(result) |> pull("row"))))
    all_res <- rbind(all_res, result)
  }
  all_res <- all_res |>
    mutate(across(everything(), as.numeric)) |>
    filter(loglik < max_loglik)
  return(all_res)
}



#' Clean a single profile-likelihood result file
#'
#' Single-file counterpart to [clean_profile()]: reads one CSV of
#' profile-likelihood results, drops rows `readr::read_csv()` flagged as
#' parsing problems, coerces every column to numeric, and keeps only fits
#' above the `max_loglik` threshold.
#'
#' @param file_name Character scalar. Path to the profile-likelihood CSV.
#' @param max_loglik Numeric scalar. Rows with `loglik >= max_loglik` are
#'   dropped (default `-500`).
#'
#' @return A cleaned data frame of accepted parameter sets.
#' @export
clean_file_result <- function(file_name, max_loglik = -500) {
  result <- read_csv(file_name, show_col_types = FALSE)

  names_param <- result |>
    dplyr::select(-any_of(c("loglik", "loglik.se", "sample", "flag", "delta"))) |>
    colnames()

  result <- result |> filter(!(row_number() %in% (problems(result) |> pull("row"))))
  all_res <- result |>
    mutate(across(everything(), as.numeric)) |>
    filter(loglik < max_loglik)
  return(all_res)
}
