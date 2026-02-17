source("R/simulation_functions.R")


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
