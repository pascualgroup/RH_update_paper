## Helpers to write an atomic JSON run manifest for fit runs

#' Write a run manifest to disk atomically
#'
#' This function writes a run manifest (an R list) to a JSON file specified by
#' `path`. The write is performed atomically by writing to a temporary file in
#' the same directory and then renaming it into place. When `append = TRUE` and
#' the target file already exists, the function attempts to read the existing
#' JSON and append the new `info` entry. Existing single-object manifests are
#' converted to an array before appending.
#'
#' @param path Character scalar. Path to the target JSON file to create or
#'   append to.
#' @param info List. A manifest information structure (typically created with
#'   `create_manifest_info()`) that will be written as JSON.
#' @param pretty Logical scalar. If `TRUE` the JSON is pretty-printed.
#' @param append Logical scalar. If `TRUE` and `path` exists, attempt to read
#'   the existing JSON and append `info` to it. If `FALSE` the file will be
#'   overwritten.
#'
#' @return Invisibly returns `path` (the file path written).
#' @export
write_run_manifest <- function(path, info, pretty = FALSE, append = FALSE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to write JSON manifests. Please install it first.")
  }
  dir <- dirname(path)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  tmp <- tempfile(pattern = "manifest", tmpdir = dir, fileext = ".json")

  if (isTRUE(append) && file.exists(path)) {
    # read existing manifest and try to append
    existing <- tryCatch(jsonlite::read_json(path, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(existing)) {
      # if we couldn't read the existing file, fall back to creating a new one
      jsonlite::write_json(info, tmp, pretty = pretty, auto_unbox = TRUE, na = "null")
    } else if (is.list(existing) && is.null(names(existing))) {
      # existing is a JSON array -> append
      new <- c(existing, list(info))
      jsonlite::write_json(new, tmp, pretty = pretty, auto_unbox = TRUE, na = "null")
    } else if (is.list(existing) && !is.null(names(existing))) {
      # existing is a single JSON object -> convert to array of two objects
      new <- list(existing, info)
      jsonlite::write_json(new, tmp, pretty = pretty, auto_unbox = TRUE, na = "null")
    } else {
      # fallback: wrap existing and new into an array
      new <- list(existing, info)
      jsonlite::write_json(new, tmp, pretty = pretty, auto_unbox = TRUE, na = "null")
    }
  } else {
    # create new (overwrite)
    jsonlite::write_json(info, tmp, pretty = pretty, auto_unbox = TRUE, na = "null")
  }

  # atomic move
  if (!file.rename(tmp, path)) {
    # fallback: try copying then unlink
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

##' Summarize a set of model parameters
#'
#' Summarise parameter grids or named numeric vectors used for fitting. For
#' data frames / tibbles the function returns a list with the type, number of
#' starts, column names and a min/max summary for each parameter. For named
#' numeric vectors a compact representation with names and values is returned.
#'
#' @param parameters A data frame / tibble (grid of parameter starts) or a
#'   named numeric vector.
#'
#' @return A list describing the parameters. The structure depends on the input
#'   type (see details).
#' @export
summarize_parameters <- function(parameters) {
  if (is.data.frame(parameters) || tibble::is_tibble(parameters)) {
    min_vals <- sapply(parameters, min, na.rm = TRUE)
    max_vals <- sapply(parameters, max, na.rm = TRUE)
    param_summary <- setNames(
      lapply(seq_along(min_vals), function(i) c(min_vals[i], max_vals[i])),
      names(min_vals)
    )
    return(list(
      type = "grid",
      n_starts = nrow(parameters),
      colnames = colnames(parameters),
      summary = param_summary
    ))
  }
  if (is.numeric(parameters) && !is.null(names(parameters))) {
    return(list(type = "vector", names = names(parameters), values = as.list(parameters)))
  }
  # fallback
  return(list(type = class(parameters)[1]))
}

##' Create a manifest information skeleton for a run
#'
#' Build a manifest information list containing metadata about a single run:
#' start time, user, host, R version and (if available) Git commit/branch.
#' The `args` slot is intended to capture the command-line / function arguments
#' that produced the run (for reproducibility).
#'
#' @param run_id Optional character. A unique run identifier. If `NULL` a
#'   timestamp-derived id will be generated.
#' @param start_time POSIXct. The run start time. Defaults to `Sys.time()`.
#' @param user Character scalar. Username for provenance. Defaults to system
#'   user.
#' @param args List. Arbitrary named list of run arguments to record.
#'
#' @return A list suitable for writing to a run manifest (see `write_run_manifest`).
#' @export
create_manifest_info <- function(run_id = NULL, start_time = Sys.time(), user = Sys.info()["user"], args = list()) {
  if (is.null(run_id)) run_id <- paste0(gsub("[: -]", "_", format(start_time, tz = "UTC")), "_", as.integer(Sys.time()) %% 100000)
  git_sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) NULL)
  git_branch <- tryCatch(system2("git", c("rev-parse", "--abbrev-ref", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) NULL)
  info <- list(
    run_id = run_id,
    start_time = format(start_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    user = user,
    hostname = tryCatch(Sys.info()["nodename"], error = function(e) NA_character_),
    r_version = R.version.string,
    git = list(sha = if (length(git_sha)) git_sha[[1]] else NULL, branch = if (length(git_branch)) git_branch[[1]] else NULL),
    args = args,
    status = "started"
  )
  info
}
