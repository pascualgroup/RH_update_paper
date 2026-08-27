## Safe write helpers for SLURM array tasks
##
## This file provides small helpers to atomically write per-task results and
## per-task manifest JSON files. Intended to be sourced by worker R scripts in
## SLURM jobs. Files are written to a temp file in the same directory and then
## renamed to the final name to make writes atomic on the same filesystem.

#' Null-coalescing operator
#'
#' Returns `a` unless it is `NULL`, in which case `b` is returned. Used by the
#' manifest-merging helpers below.
#'
#' @param a Value to prefer.
#' @param b Fallback value used when `a` is `NULL`.
#' @return `a`, or `b` if `a` is `NULL`.
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Build a unique per-task filename suffix from SLURM environment variables
#'
#' Construct a unique suffix from SLURM environment variables (job id, array
#' task id, process id, local id, hostname), falling back to the current PID
#' and hostname outside of a SLURM job. Used by [safe_write_csv()] and
#' [safe_write_json()] so concurrent array tasks don't collide when writing
#' per-task output files.
#'
#' @return A character scalar suffix, e.g. `"12345_A3_P0_node07"`.
#' @keywords internal
get_slurm_suffix <- function() {
  jid <- Sys.getenv("SLURM_JOB_ID", unset = Sys.getenv("JOB_ID", "local"))
  a_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  proc <- Sys.getenv("SLURM_PROCID", unset = as.character(Sys.getpid()))
  local <- Sys.getenv("SLURM_LOCALID", unset = "")
  host <- Sys.info()[["nodename"]]
  parts <- c(jid)
  if (nzchar(a_id)) parts <- c(parts, paste0("A", a_id))
  if (nzchar(local)) parts <- c(parts, paste0("L", local))
  parts <- c(parts, paste0("P", proc), host)
  paste(parts, collapse = "_")
}

#' Ensure a directory exists
#'
#' Create a directory (and parents) if it does not already exist. This is a
#' convenience wrapper used by the safe write helpers to ensure the target
#' directory is present before attempting atomic writes.
#'
#' @param dir Character scalar: path to the directory to ensure exists.
#' @return Invisibly returns the input `dir`.
#' @keywords internal
ensure_dir <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(dir)
}

#' Atomically write a data.frame to CSV
#'
#' Write a data.frame or tibble to a CSV file atomically by first writing to a
#' temporary file in the same directory and then renaming it to the final
#' filename. The final filename is constructed from `prefix` and a unique
#' suffix that by default encodes SLURM environment information (so array
#' tasks don't collide).
#'
#' @param df A data.frame or tibble to be written.
#' @param dir Directory where the file will be written (created if missing).
#' @param prefix File name prefix. The final filename will be `paste0(prefix, "_", suffix, ".csv")`.
#' @param suffix Optional suffix to use instead of the SLURM-derived suffix.
#' @param row.names Passed to `write.csv()`; default FALSE.
#' @param ... Additional arguments passed to `utils::write.csv()`.
#' @return Invisibly returns the path to the final CSV file written.
#' @export
safe_write_csv <- function(df, dir = "results", prefix = "run", suffix = NULL, row.names = FALSE, ...) {
  if (!is.data.frame(df)) stop("df must be a data.frame or tibble")
  ensure_dir(dir)
  if (is.null(suffix)) suffix <- get_slurm_suffix()
  final <- file.path(dir, sprintf("%s_%s.csv", prefix, suffix))
  tmp <- file.path(dir, sprintf(".%s_%s.csv.tmp", prefix, paste0(suffix, "_", sample.int(1e6, 1))))
  # write to tmp then rename
  # prefer readr::write_csv when available for performance; fall back to utils::write.csv
  if (isTRUE(row.names)) {
    # if requested, materialize row names as a column (name it .rownames unless that exists)
    rn <- rownames(df)
    if (!is.null(rn)) {
      colname <- ".rownames"
      if (colname %in% names(df)) {
        # find a free name
        i <- 1L
        while (paste0(colname, i) %in% names(df)) i <- i + 1L
        colname <- paste0(colname, i)
      }
      df <- cbind(stats::setNames(data.frame(rn, stringsAsFactors = FALSE), colname), df)
      rownames(df) <- NULL
    }
  }
  if (requireNamespace("readr", quietly = TRUE)) {
    # readr::write_csv does not support row.names; we've already handled row names above
    readr::write_csv(df, tmp, ...)
  } else {
    utils::write.csv(df, tmp, row.names = FALSE, ...)
  }
  if (!file.rename(tmp, final)) {
    # fallback: copy then unlink (less atomic if across filesystems)
    file.copy(tmp, final, overwrite = FALSE)
    unlink(tmp)
  }
  invisible(final)
}

#' Atomically write a list/object manifest as JSON
#'
#' Serialize an R object to JSON and write it atomically to the target
#' directory. Uses `jsonlite::write_json()`; jsonlite must be available.
#'
#' @param obj R object (usually a list) that will be converted to JSON.
#' @param dir Directory where the JSON will be written.
#' @param prefix File name prefix. Final file is `paste0(prefix, "_", suffix, ".json")`.
#' @param suffix Optional suffix; if NULL the SLURM-derived suffix is used.
#' @param pretty Logical, passed to `jsonlite::write_json()`.
#' @param auto_unbox Logical, passed to `jsonlite::write_json()`.
#' @return Invisibly returns the path to the final JSON file written.
#' @export
safe_write_json <- function(obj, dir = "manifests", prefix = "manifest", suffix = NULL, pretty = FALSE, auto_unbox = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite required for safe_write_json")
  ensure_dir(dir)
  if (is.null(suffix)) suffix <- get_slurm_suffix()
  final <- file.path(dir, sprintf("%s_%s.json", prefix, suffix))
  tmp <- file.path(dir, sprintf(".%s_%s.json.tmp", prefix, paste0(suffix, "_", sample.int(1e6, 1))))
  jsonlite::write_json(obj, tmp, pretty = pretty, auto_unbox = auto_unbox)
  if (!file.rename(tmp, final)) {
    file.copy(tmp, final, overwrite = FALSE)
    unlink(tmp)
  }
  invisible(final)
}

#' List per-task manifest files
#'
#' Return the list of per-task manifest files created by workers. This is a
#' small convenience used by aggregator scripts.
#'
#' @param dir Directory to search (default: "manifests").
#' @param pattern Regex pattern used to match manifest filenames.
#' @return Character vector of full paths (possibly empty).
#' @export
list_manifests <- function(dir = "manifests", pattern = "^manifest_.*\\.json$") {
  if (!dir.exists(dir)) {
    return(character(0))
  }
  list.files(dir, pattern = pattern, full.names = TRUE)
}

#' Read a manifest safely
#'
#' Safely read a JSON manifest file and return it as an R object. On error
#' (missing file or malformed JSON) this returns NULL and emits a warning.
#'
#' @param path Path to the JSON manifest file.
#' @return The parsed manifest (as a list) or NULL on failure.
#' @export
safe_read_manifest <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite required for safe_read_manifest")
  tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) {
    warning("failed to read manifest: ", path, " -> ", conditionMessage(e))
    NULL
  })
}

# save results

# Save results safely: prefer `safe_write_csv()` from R/save_safe.R when present
# Falls back to a timestamped write.csv() if helper is not available.
#' Save results with safe write semantics
#'
#' Attempt to write results using `safe_write_csv()` if available (atomic
#' write). If the helper is not present this falls back to a timestamped
#' `write.csv()` file. This is intended for use within SLURM tasks where a
#' per-task file is preferred.
#'
#' @param df A data.frame or tibble with results to save.
#' @param dir Directory where results should be written.
#' @param prefix Filename prefix for the written file.
#' @return Invisibly returns the path to the written CSV file, or NULL if `df` is NULL.
#' @export
save_results_safe <- function(df, dir = "results", prefix = "output", suffix = NULL) {
  if (is.null(df)) {
    return(NULL)
  }
  if (exists("safe_write_csv", mode = "function")) {
    tryCatch(
      {
        path <- safe_write_csv(df, dir = dir, prefix = prefix, suffix = suffix)
        message("wrote:", path)
        return(path)
      },
      error = function(e) {
        warning("safe_write_csv failed: ", conditionMessage(e))
      }
    )
  }
  # fallback: timestamp + pid filename
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  suf <- format(Sys.time(), "%Y%m%dT%H%M%S")
  pid <- Sys.getpid()
  final <- file.path(dir, sprintf("%s_%s_pid%s.csv", prefix, suf, pid))
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(df, final)
  } else {
    utils::write.csv(df, final, row.names = FALSE)
  }
  message("wrote (fallback):", final)
  invisible(final)
}

## helper: write a manifest JSON next to a result CSV atomically
#' Write manifest JSON next to a result CSV
#'
#' Write the provided manifest object as JSON in the same directory as
#' `result_path` and with the same base filename (but with a `.json`
#' extension). The write is atomic when possible (write to a temporary
#' file then rename). If `jsonlite` is not available the function falls back
#' to `safe_write_json()` (which writes into `manifests/`).
#'
#' @param manifest_obj A list representing the manifest to write.
#' @param result_path Path to the result CSV file; the manifest will be written next to it.
#' @return Invisibly returns the path to the written manifest file, or NULL on failure.
#' @export
safe_write_manifest_next_to <- function(manifest_obj, result_path, dir = NULL) {
  if (is.null(result_path) || !nzchar(result_path)) {
    return(NULL)
  }
  if (is.null(dir)) {
    dir <- dirname(result_path)
  }
  # ensure the target directory exists
  ensure_dir(dir)
  base <- tools::file_path_sans_ext(basename(result_path))
  final <- file.path(dir, paste0(base, ".json"))
  tmp <- file.path(dir, paste0(".", base, "_", sample.int(1e6, 1), ".json.tmp"))

  # prefer jsonlite for direct write next to CSV
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    # write compact JSON by default to reduce line breaks and file size
    jsonlite::write_json(manifest_obj, tmp, pretty = FALSE, auto_unbox = TRUE)
    if (!file.rename(tmp, final)) {
      file.copy(tmp, final, overwrite = FALSE)
      unlink(tmp)
    }
    return(final)
  }

  # fallback to package helper if present (writes into manifests/)
  if (exists("safe_write_json", mode = "function")) {
    return(safe_write_json(manifest_obj, dir = "manifests", prefix = base))
  }

  warning("no JSON writer available to write manifest for ", result_path)
  NULL
}

#' Write a run manifest describing a result file
#'
#' Build (or update) a manifest object for a fitting/profiling result and
#' write it next to `result_path` via [safe_write_manifest_next_to()]. If
#' `manifest_info` is supplied (e.g. the info object returned alongside a fit)
#' it is reused and annotated with the result file path, status, and
#' `mode_label`; otherwise a minimal manifest is constructed from scratch.
#'
#' @param df The result data frame written to `result_path`; used only to
#'   determine `status` (`"finished"` if it has rows, `"empty"` otherwise).
#' @param result_path Path to the result CSV file the manifest describes.
#' @param mode_label Character scalar identifying which fitting mode produced
#'   the result (recorded under `parameters$mode`).
#' @param manifest_info Optional list (e.g. from [create_manifest_info()]) to
#'   reuse and annotate instead of building a manifest from scratch.
#' @param dir Optional directory to write the manifest into; defaults to
#'   `dirname(result_path)`.
#'
#' @return Invisibly returns the path to the written manifest file, or `NULL`
#'   if `df` or `result_path` is `NULL`.
#' @export
write_manifest_for_result <- function(df, result_path, mode_label, manifest_info = NULL, dir = NULL) {
  if (is.null(df) || is.null(result_path)) {
    return(NULL)
  }

  # If the fitting function returned a manifest/info object (second element),
  # prefer that and only update/add fields that point to the local result file.
  if (!is.null(manifest_info) && is.list(manifest_info)) {
    # ensure it's a list we can modify safely
    manifest_obj <- manifest_info
    # override or add the result file path so filename matches
    manifest_obj$result_file <- result_path
    # update status based on presence of rows
    manifest_obj$status <- if (nrow(df) > 0) "finished" else "empty"
    # annotate which mode produced this file (helpful when combining manifests)
    if (is.null(manifest_obj$parameters)) manifest_obj$parameters <- list()
    manifest_obj$parameters$mode <- mode_label
  } else {
    # fallback: construct a minimal manifest object
    manifest_obj <- list(
      run_id = paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid()),
      start_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
      user = Sys.info()[["user"]],
      hostname = Sys.info()[["nodename"]],
      r_version = R.version$version.string,
      args = list(
        city = get0("city", ifnotfound = NULL),
        model = get0("model", ifnotfound = NULL),
        scale = get0("scale", ifnotfound = NULL),
        covariate = get0("covariate", ifnotfound = NULL),
        extra_cov = get0("extra_covs", ifnotfound = NULL)
      ),
      status = if (nrow(df) > 0) "finished" else "empty",
      parameters = list(n_rows = nrow(df), mode = mode_label),
      result_file = result_path
    )
  }

  safe_write_manifest_next_to(manifest_obj, result_path, dir)
}


NULL
