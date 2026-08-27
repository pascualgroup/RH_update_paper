# Script to run the fitting procedure on a SLURM cluster
# Expects command line arguments:
# 1: city (Ahmedabad or Surat)
# 2: model (const, linear, exp, poly, linear_inf, PET, PET_v2)
# 3: scale (monthly, daily or hourly)
# 4: covariate (name of the covariate column in the VC dataset)
# 5: VC_file (filename of the VC dataset, relative to the city data folder)
# 6+: (optional) extra_cov_dat (filename of the extra covariate dataset, relative to the city data folder)
#      extra_covs (names of extra covariate columns in the extra_cov_dat, can be multiple)

# assumes the following environment variables are set by SLURM:
# SLURM_ARRAY_TASK_ID: the index of the current array job (1-based)
# SLURM_NTASKS: number of tasks (cores) allocated to this job
# SLURM_JOB_ID: the job ID (for manifest filenames)
# SLURM_PROCID: the process ID within the job (for manifest filenames)
# SLURM_LOCALID: the local task ID on the node (for manifest filenames)

# assumes the rhmalaria package is installed (see README.md for cluster
# install instructions); current working directory does not matter
library(rhmalaria)

run_id <- Sys.getenv("RUN_ID")
normalize_extra <- Sys.getenv("NORMALIZE_EXTRA", unset = 1)
normalize_extra <- as.logical(as.integer(normalize_extra))
normalize_covariate <- Sys.getenv("NORMALIZE_EXTRA", unset = 1)
normalize_covariate <- as.logical(as.integer(normalize_covariate))
end_year <- as.numeric(Sys.getenv("END_YEAR", unset = 2011))

arguments <- commandArgs(trailingOnly = TRUE)

cat(length(arguments))
cat("\n")
cat(arguments)
cat("\n")
city <- arguments[1]
model <- arguments[2]
covariate <- arguments[3]
window_start <- as.numeric(arguments[4])
window_end <- as.numeric(arguments[5])
dataset <- arguments[6]

# city <- "Ahmedabad"
# model <- "inf_exponent"
# covariate <- "hadISD"
# window_start <- 8
# window_end <- 9
# dataset <- "dataset_Ahmedabad_2023.csv"


n_cores <- as.numeric(Sys.getenv("SLURM_NTASKS", unset = 10))
# n_cores <-

if (city == "Surat") {
  dataset <- paste0("data/Surat/", dataset)
} else {
  dataset <- paste0("data/Ahmedabad/", dataset)
}

# if (length(arguments) > 5) {
#   if (city == "Surat") {
#     extra_cov_dat <- paste0("data/Surat/", arguments[6])
#     extra_covs <- arguments[7:length(arguments)]
#   } else {
#     extra_cov_dat <- paste0("data/Ahmedabad/", arguments[6])
#     extra_covs <- arguments[7:length(arguments)]
#   }
# } else {
extra_cov_dat <- NULL
extra_covs <- NULL
# }
format_optional_value <- function(value, empty_label = "None") {
  if (is.null(value) || length(value) == 0) {
    return(empty_label)
  }
  paste(value, collapse = ";")
}

summarize_run_parameters <- function() {
  paste(
    sprintf("run_id=%s", if (nzchar(run_id)) run_id else "None"),
    sprintf("city=%s", city),
    sprintf("model=%s", model),
    sprintf("scale=%s", scale),
    sprintf("covariate=%s", covariate),
    sprintf("VC_file=%s", VC_file),
    sprintf("extra_cov_dat=%s", format_optional_value(extra_cov_dat)),
    sprintf("extra_covariates=%s", format_optional_value(extra_covs)),
    sprintf("param_index=%s", param_index),
    sprintf("group_size=%s", group_size),
    sprintf("start_index=%s", start_index),
    sprintf("end_index=%s", end_index),
    sprintf("n_cores=%s", n_cores),
    sprintf("normalize_covariate=%s", normalize_covariate),
    sprintf("normalize_extra=%s", normalize_extra),
    sep = ", "
  )
}

stop_if_empty <- function(result_obj, stage_label) {
  res_tbl <- result_obj[[1]]
  if (is.null(res_tbl) || nrow(res_tbl) == 0) {
    stop(sprintf(
      "fit_model returned no rows during %s stage. Parameters: %s",
      stage_label,
      summarize_run_parameters()
    ))
  }
}

param_index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = 1))
midway_max_jobs <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT",
  unset = Sys.getenv("SLURM_ARRAY_TASK_MAX",
    unset = Sys.getenv("SLURM_ARRAY_SIZE",
      unset = "2500"
    )
  )
))

param_grid_files <- list(
  inf_exponent = "param_grid_inf_exponent.csv",
  baseline = "param_grid_baseline.csv"
)

if (!model %in% names(param_grid_files)) {
  stop(paste("Unknown model:", model))
}

s1 <- read.csv(file.path("param_grids", param_grid_files[[model]]))


group_size <- floor(nrow(s1) / midway_max_jobs)
start_index <- (param_index - 1) * group_size + 1
end_index <- param_index * group_size

parameters <- s1[start_index:end_index, ]

result_fit <- fit_model(
  parameters,
  city = city,
  model = model,
  covariate = covariate,
  start_year = 1997, start_month = 1,
  end_year = end_year, end_month = 12,
  window_start = window_start,
  window_end = window_end,
  dataset = dataset,
  mode = "fit", allow_parallel = TRUE,
  n_cores = n_cores,
  seed_num = as.integer(Sys.time()),
  start_index = start_index, output_to_file = FALSE,
  output_format = "data.frame"
)


stop_if_empty(result_fit, "fit")

parameters <- result_fit[[1]] |>
  arrange(desc(loglik)) |>
  select(-any_of(c("sample", "loglik", "loglik.se", "flag")))

result_refined <- fit_model(
  parameters,
  city = city,
  model = model,
  covariate = covariate,
  start_year = 1997, start_month = 1,
  end_year = end_year, end_month = 12,
  window_start = window_start,
  window_end = window_end,
  dataset = dataset,
  mode = "refined", allow_parallel = TRUE,
  n_cores = n_cores,
  seed_num = as.integer(Sys.time()),
  start_index = start_index, output_to_file = FALSE,
  output_format = "data.frame"
)

stop_if_empty(result_refined, "refined")


parameters <- result_refined[[1]] |>
  arrange(desc(loglik)) |>
  select(-any_of(c("sample", "loglik", "loglik.se", "flag")))

result_refined_second <- fit_model(
  parameters,
  city = city,
  model = model,
  covariate = covariate,
  start_year = 1997, start_month = 1,
  end_year = end_year, end_month = 12,
  window_start = window_start,
  window_end = window_end,
  dataset = dataset,
  mode = "refined_second", allow_parallel = TRUE,
  n_cores = n_cores,
  seed_num = as.integer(Sys.time()),
  start_index = start_index, output_to_file = FALSE,
  output_format = "data.frame"
)

stop_if_empty(result_refined_second, "refined_second")

## only array id matters for suffix
a_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
a_id <- paste0("A", a_id)

## files are saved in the pattern output_<run_id>_<mode>_A<array_id>.csv
## manifest ends in json
# write the three main result tables (if present) and create matching manifests
res1 <- tryCatch(save_results_safe(result_fit[[1]], dir = "results", prefix = paste("output", run_id, "fit", sep = "_"), suffix = a_id), error = function(e) {
  warning(conditionMessage(e))
  NULL
})
if (!is.null(res1)) {
  try(write_manifest_for_result(result_fit[[1]], res1, "fit", manifest_info = result_fit[[2]], dir = "manifests"), silent = TRUE)
}

res2 <- tryCatch(save_results_safe(result_refined[[1]], dir = "results", prefix = paste("output", run_id, "refined", sep = "_"), suffix = a_id), error = function(e) {
  warning(conditionMessage(e))
  NULL
})
if (!is.null(res2)) {
  try(write_manifest_for_result(result_refined[[1]], res2, "refined", manifest_info = result_refined[[2]], dir = "manifests"), silent = TRUE)
}

res3 <- tryCatch(save_results_safe(result_refined_second[[1]], dir = "results", prefix = paste("output", run_id, "refined_second", sep = "_"), suffix = a_id), error = function(e) {
  warning(conditionMessage(e))
  NULL
})
if (!is.null(res3)) {
  try(write_manifest_for_result(result_refined_second[[1]], res3, "refined_second", manifest_info = result_refined_second[[2]], dir = "manifests"), silent = TRUE)
}
