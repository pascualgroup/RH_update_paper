#!/usr/bin/env Rscript
# Combine per-array result CSVs and JSON manifests into single files.
#
# Usage:
#   Rscript scripts/combine_array_outputs.R --run_id RUN_ID [--array ARRAYID] [--results-dir results] [--manifests-dir manifests] [--out-dir results/aggregated]
#
# Description:
#   This script collects CSV result files and JSON manifest files produced by
#   individual tasks (written with the safe writers) and combines them into a
#   single CSV and a single JSON manifest. Files are selected by matching the
#   SLURM job id and (optionally) the array index. Writes are atomic (tmp file
#   then rename).

library(jsonlite)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(argv) {
  out <- list(run_id = NULL, array = NULL, results_dir = "results", manifests_dir = "manifests", out_dir = file.path("results", "aggregated"))
  for (a in argv) {
    if (grepl("^--run_id=", a)) {
      out$run_id <- sub("^--run_id=", "", a)
    } else if (grepl("^--array=", a)) {
      out$array <- sub("^--array=", "", a)
    } else if (grepl("^--results-dir=", a)) {
      out$results_dir <- sub("^--results-dir=", "", a)
    } else if (grepl("^--manifests-dir=", a)) {
      out$manifests_dir <- sub("^--manifests-dir=", "", a)
    } else if (grepl("^--out-dir=", a)) out$out_dir <- sub("^--out-dir=", "", a)
  }
  out
}

opts <- parse_args(args)
if (is.null(opts$run_id) || nchar(opts$run_id) == 0) {
  stop("Please provide --run_id RUN_ID (the run id used in filenames)")
}

ensure_dir <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  normalizePath(d, mustWork = FALSE)
}
opts$results_dir <- ensure_dir(opts$results_dir)
opts$manifests_dir <- ensure_dir(opts$manifests_dir)
opts$out_dir <- ensure_dir(opts$out_dir)

pattern_match <- function(files, job, array = NULL) {
  # match files containing the job id, and if array provided, the pattern _A<array>
  m <- grepl(job, files, fixed = TRUE)
  if (!is.null(array) && nzchar(array)) {
    m <- m & grepl(paste0("_A", array), files, fixed = TRUE)
  }
  files[m]
}

cat(sprintf("Combining files for run_id=%s array=%s\n", opts$run_id, ifelse(is.null(opts$array), "<any>", opts$array)))

## Combine per-mode CSVs and manifests
csv_files_all <- list.files(opts$results_dir, pattern = "\\.csv$", full.names = TRUE)
csv_files <- pattern_match(csv_files_all, opts$run_id, opts$array)

if (length(csv_files) == 0) {
  cat("No CSV result files found for the requested run_id/array. Exiting with status 0.\n")
} else {
  cat(sprintf("Found %d CSV files, combining per-mode...\n", length(csv_files)))
  csv_files <- sort(csv_files)

  modes <- c("fit", "refined", "refined_second")
  # track combined files we create so we don't accidentally delete them below
  combined_csv_written <- character(0)
  combined_json_written <- character(0)
  for (mode in modes) {
    b <- basename(csv_files)
    if (mode == "refined") {
      files_for_mode <- csv_files[grepl("_refined_", b, fixed = TRUE) &
        !grepl("_refined_second", b, fixed = TRUE)]
    } else {
      files_for_mode <- csv_files[grepl(paste0("_", mode, "_"), b, fixed = TRUE)]
    }

    if (length(files_for_mode) == 0) {
      cat(sprintf("No CSVs for mode '%s', skipping.\n", mode))
      next
    }
    cat(sprintf("Mode '%s': combining %d CSV files\n", mode, length(files_for_mode)))
    df_list <- lapply(sort(files_for_mode), function(f) {
      tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) {
        warning(sprintf("Failed to read %s: %s", f, e$message))
        NULL
      })
    })
    df_list <- Filter(Negate(is.null), df_list)
    if (length(df_list) == 0) {
      warning(sprintf("No readable CSVs for mode '%s'; skipping.", mode))
      next
    }
    combined <- do.call(dplyr::bind_rows, df_list)
    out_csv <- file.path(opts$out_dir, sprintf("combined_results_%s_%s%s.csv", opts$run_id, mode, ifelse(is.null(opts$array), "", paste0("_A", opts$array))))
    tmp <- tempfile(pattern = paste0("combined_results_", mode), tmpdir = opts$out_dir, fileext = ".csv")
    readr::write_csv(combined, tmp)
    if (!file.rename(tmp, out_csv)) {
      file.copy(tmp, out_csv, overwrite = TRUE)
      unlink(tmp)
    }
    combined_csv_written <- c(combined_csv_written, normalizePath(out_csv, winslash = "/", mustWork = FALSE))
    cat(sprintf("Wrote combined CSV for mode '%s' to %s\n", mode, out_csv))

    # Now combine manifests for this mode
    manifest_files_all <- list.files(opts$manifests_dir, pattern = "\\.json$", full.names = TRUE)
    manifest_files <- pattern_match(manifest_files_all, opts$run_id, opts$array)
    files_for_mode_m <- manifest_files[grepl(paste0("_", mode, "_"), basename(manifest_files), fixed = TRUE)]
    if (length(files_for_mode_m) == 0) {
      cat(sprintf("No manifests for mode '%s', skipping manifest combine.\n", mode))
      next
    }
    cat(sprintf("Mode '%s': combining %d manifest files\n", mode, length(files_for_mode_m)))
    manifests <- lapply(sort(files_for_mode_m), function(f) {
      tryCatch(jsonlite::read_json(f, simplifyVector = FALSE), error = function(e) {
        warning(sprintf("Failed to read JSON %s: %s", f, e$message))
        NULL
      })
    })
    manifests <- Filter(Negate(is.null), manifests)
    combined_manifest <- list(
      combined_from = basename(files_for_mode_m),
      n = length(manifests),
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      manifests = manifests
    )
    out_json <- file.path(opts$out_dir, sprintf("combined_manifest_%s_%s%s.json", opts$run_id, mode, ifelse(is.null(opts$array), "", paste0("_A", opts$array))))
    tmpj <- tempfile(pattern = paste0("combined_manifest_", mode), tmpdir = opts$out_dir, fileext = ".json")
    # write compact JSON (few line breaks) to save space and make files denser
    jsonlite::write_json(combined_manifest, tmpj, pretty = FALSE, auto_unbox = TRUE, na = "null")
    if (!file.rename(tmpj, out_json)) {
      file.copy(tmpj, out_json, overwrite = TRUE)
      unlink(tmpj)
    }
    combined_json_written <- c(combined_json_written, normalizePath(out_json, winslash = "/", mustWork = FALSE))
    cat(sprintf("Wrote combined manifest for mode '%s' to %s\n", mode, out_json))
  }
}

# Clean up: remove per-task input files that belong to this run_id (and optional array)
# but be careful not to remove any combined outputs.
cat("Cleaning up per-task files for run_id=", opts$run_id, " array=", ifelse(is.null(opts$array), "<any>", opts$array), "\n")

# results
all_results <- list.files(opts$results_dir, pattern = "\\.csv$", full.names = TRUE)
candidates_results <- pattern_match(all_results, opts$run_id, opts$array)
if (length(candidates_results) > 0) {
  # exclude any combined files we wrote (and any file whose basename starts with combined_results_)
  is_combined <- basename(candidates_results) %in% basename(combined_csv_written) | grepl("^combined_results_", basename(candidates_results))
  to_delete_results <- candidates_results[!is_combined]
  if (length(to_delete_results) > 0) {
    cat(sprintf("Removing %d per-task result files...\n", length(to_delete_results)))
    sapply(to_delete_results, function(f) {
      if (file.exists(f)) {
        ok <- file.remove(f)
        cat(sprintf("  %s -> %s\n", basename(f), ifelse(ok, "deleted", "failed")))
      }
    })
  } else {
    cat("No per-task result files to remove (all matched combined outputs).\n")
  }
} else {
  cat("No per-task result files found for cleanup.\n")
}

# manifests
all_manifs <- list.files(opts$manifests_dir, pattern = "\\.json$", full.names = TRUE)
candidates_manifs <- pattern_match(all_manifs, opts$run_id, opts$array)
if (length(candidates_manifs) > 0) {
  is_combined_m <- basename(candidates_manifs) %in% basename(combined_json_written) | grepl("^combined_manifest_", basename(candidates_manifs))
  to_delete_manifs <- candidates_manifs[!is_combined_m]
  if (length(to_delete_manifs) > 0) {
    cat(sprintf("Removing %d per-task manifest files...\n", length(to_delete_manifs)))
    sapply(to_delete_manifs, function(f) {
      if (file.exists(f)) {
        ok <- file.remove(f)
        cat(sprintf("  %s -> %s\n", basename(f), ifelse(ok, "deleted", "failed")))
      }
    })
  } else {
    cat("No per-task manifest files to remove (all matched combined outputs).\n")
  }
} else {
  cat("No per-task manifest files found for cleanup.\n")
}

cat("Done.\n")
