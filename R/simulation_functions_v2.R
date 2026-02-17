library(tidyverse)
library(lubridate)
library(pomp)
library(doParallel)
library(doRNG)
library(ggtext)
library(doFuture)
library(circumstance)

# centralized helper that checks either explicit arg or an option
warn_if <- function(msg, verbose = NULL, call. = FALSE) {
  # prefer explicit verbose; otherwise check package-level option
  if (is.null(verbose)) verbose <- getOption("vectorial.verbose", FALSE)
  if (isTRUE(verbose)) warning(msg, call. = call.)
  invisible(NULL)
}


# Helper to generate output/result filenames used by `fit_model`
generate_output_names <- function(city, model, scale, covariate, start_year, start_month, end_year, end_month, mode, profile_var = NULL) {
  if (!is.null(profile_var)) {
    output_result <- paste("output_profile", profile_var, city, model,
      scale, covariate, start_year, start_month,
      end_year, end_month, mode,
      sep = "_"
    ) |> paste0(".csv")
    output_log <- paste("log_profile", profile_var, city, model,
      scale, covariate, start_year, start_month,
      end_year, end_month, mode,
      sep = "_"
    ) |> paste0(".log")
  } else {
    output_result <- paste("output", city, model, scale, covariate,
      start_year, start_month,
      end_year, end_month, mode,
      sep = "_"
    ) |> paste0(".csv")
    output_log <- paste("log", city, model, scale, covariate,
      start_year, start_month,
      end_year, end_month, mode,
      sep = "_"
    ) |> paste0(".log")
  }
  list(output_result = output_result, output_log = output_log)
}

yhour <- function(date_hour) {
  day_of_year <- yday(date_hour)
  # Get the hour of the day
  hour_of_day <- hour(date_hour)
  hour_of_year <- (day_of_year - 1) * 24 + hour_of_day
  return(hour_of_year)
}

days_in_year <- function(year) {
  result <- if_else(leap_year(year), 366, 365)
  return(result)
}

hours_in_year <- function(year) {
  result <- days_in_year(year) * 24
  return(result)
}

mins_in_year <- function(year) {
  result <- hours_in_year(year) * 60
  return(result)
}

find_datetime_column <- function(df) {
  accepted <- c("Date", "POSIXct", "POSIXt")
  for (nm in names(df)) {
    cls <- class(df[[nm]])[1]
    if (cls %in% accepted) {
      return(nm)
    }
  }
  return(NULL)
}

filter_by_period <- function(df, start_year, start_month, end_year, end_month) {
  df |>
    filter(year >= start_year, year <= end_year) |>
    filter(year != start_year | (year == start_year & month >= start_month)) |>
    filter(year != end_year | (year == end_year & month <= end_month))
}

normalize_minmax_safe <- function(x) {
  x_num <- as.numeric(x)
  minx <- min(x_num, na.rm = TRUE)
  maxx <- max(x_num, na.rm = TRUE)
  if (is.infinite(minx) || is.infinite(maxx) || (maxx - minx) <= .Machine$double.eps) {
    # constant or all-NA -> return NA vector (explicit)
    return(rep(NA_real_, length(x_num)))
  }
  (x_num - minx) / (maxx - minx)
}

safe_pop_smooth <- function(time, pop, method = "smooth.spline", params = list(spar = 0.5)) {
  n <- length(time)
  if (n < 5 || length(unique(pop[!is.na(pop)])) < 3) {
    warning("Too few points to smooth population reliably; setting dpopdt = NA")
    return(list(pop = pop, dpopdt = rep(NA_real_, n)))
  }
  if (method == "smooth.spline") {
    spl_both <- loess(pop ~ time)$fitted
    spl_both <- smooth.spline(x = time, spl_both)
    pop <- predict(spl_both, deriv = 0, x = time)$y
    dpopdt <- predict(spl_both, deriv = 1, x = time)$y
    return(list(pop = pop, dpopdt = dpopdt))
  } else {
    stop("method not supported")
  }
}


#' Preprocess covariate and malaria case data for POMP modeling
#'
#' This function prepares covariate (climate/VC) and malaria case datasets for
#' use with the POMP models in this repository. It:
#' - selects and aggregates the requested covariate at the specified `scale` ("monthly", "daily", or "hourly");
#' - finds and expands date/time components (year, month, day, hour where applicable);
#' - optionally merges additional covariates (`extra_cov`) and renames them to `cov1`, `cov2`, ...;
#' - normalizes extra covariates with a safe min-max transform (constant or all-NA columns become NA);
#' - fits a smooth population curve and computes `dpopdt` (time derivative) used by the models.
#'
#' The function returns a list with two elements: `dat` (the covariate table, including `time`, `pop`, `dpopdt`) and
#' `dat_mal` (malaria observations filtered to the selected period). Many downstream functions expect `time` to be a
#' numeric fractional-year index.
#'
#' @param dat A data frame containing covariate/time information. Must include a date or POSIX column.
#' @param dat_mal A data frame with malaria case information (must contain `year`, `month`, and `pop` for smoothing).
#' @param scale Character. One of `"monthly"`, `"daily"`, or `"hourly"` indicating aggregation scale.
#' @param covariate Character. Column name in `dat` to use as the primary covariate.
#' @param start_year Integer. First year to include.
#' @param end_year Integer. Last year to include.
#' @param start_month Integer. Starting month in `start_year`.
#' @param end_month Integer. Final month in `end_year`.
#' @param extra_cov Optional character vector of additional covariate column names to include.
#' @param extra_cov_dataset Optional data frame containing the `extra_cov` columns (if not present in `dat`).
#' @param normalize_extra Logical. If `TRUE`, extra covariates are normalized using a safe min-max transform.
#'
#' @return A list with elements `dat` and `dat_mal`. `dat` includes `time`, `pop`, `dpopdt`, and covariates (named `cov1`, `cov2`, ...).
#'
#' @details
#' - Normalization policy: min-max scaling is applied with `na.rm = TRUE`. If a column is constant or all NA,
#'   the function returns an explicit `NA` vector for that covariate to avoid NaN/Inf results.
#' - Date detection: the function searches for the first column of class `Date`, `POSIXct` or `POSIXt` in the provided
#'   data frames and extracts `year`, `month`, `day`, and `hour` (when present).
#'
preprocess_data <- function(dat, dat_mal, scale, covariate, start_year, end_year,
                            start_month, end_month, extra_cov = NULL, extra_cov_dataset = NULL,
                            normalize_covariate = TRUE, normalize_extra = TRUE) {
  dat <- dat |>
    mutate(covariate = dat[[covariate]]) |>
    drop_na(covariate)

  if (!is.null(extra_cov)) {
    if (!is.null(extra_cov_dataset)) {
      ## check if extra_cov is part of extra_cov_dataset
      if (all(extra_cov %in% colnames(extra_cov_dataset))) {
        ## all of extra_cov are in the dataset
        extra_cov_dataset <- extra_cov_dataset |>
          select(any_of(c(
            "date", "year", "month",
            "day", "hour", "minute", extra_cov
          )))
      } else {
        ## at least one of extra_cov is not in the dataset, throw error
        stop("error when preprocessing data, at least one of extra_cov is not in extra_cov_dataset")
      }
    } else {
      ## supposes extra_cov is in dat
      if (all(extra_cov %in% colnames(dat))) {
        ## all of extra_cov are in the dataset
        extra_cov_dataset <- dat |>
          select(any_of(c(
            "date", "year", "month",
            "day", "hour", "minute", extra_cov
          )))
      } else {
        ## at least one of extra_cov is not in the dataset, throw error
        stop("error when preprocessing data, at least one of extra_cov is not in the covariate dataset")
      }
    }
  }

  column_names <- find_datetime_column(dat)
  if (is.null(column_names)) {
    stop("No date or POSIX column found in covariate dataset")
  }

  if (class(dat[[column_names]])[1] == "Date") {
    dat <- dat |>
      mutate(
        year = year(!!sym(column_names)),
        month = month(!!sym(column_names)),
        day = day(!!sym(column_names))
      )
  } else if (class(dat[[column_names]])[1] %in% c("POSIXct", "POSIXt")) {
    dat <- dat |>
      mutate(
        year = year(!!sym(column_names)),
        month = month(!!sym(column_names)),
        day = day(!!sym(column_names)),
        hour = hour(!!sym(column_names))
      )
  }

  if (!is.null(extra_cov_dataset)) {
    column_names <- find_datetime_column(extra_cov_dataset)
    if (is.null(column_names)) {
      stop("No date or POSIX column found in extra covariate dataset")
    }
    if (class(extra_cov_dataset[[column_names]])[1] == "Date") {
      extra_cov_dataset <- extra_cov_dataset |>
        mutate(
          year = year(!!sym(column_names)),
          month = month(!!sym(column_names)),
          day = day(!!sym(column_names))
        )
    } else if (class(extra_cov_dataset[[column_names]])[1] %in% c("POSIXct", "POSIXt")) {
      extra_cov_dataset <- extra_cov_dataset |>
        mutate(
          year = year(!!sym(column_names)),
          month = month(!!sym(column_names)),
          day = day(!!sym(column_names)),
          hour = hour(!!sym(column_names))
        )
    }
  }



  dat <- filter_by_period(dat, start_year, start_month, end_year, end_month)
  if (!is.null(extra_cov_dataset)) {
    new_names <- setNames(extra_cov, paste0("cov", seq_along(extra_cov)))
    extra_cov_dataset <- filter_by_period(extra_cov_dataset, start_year, start_month, end_year, end_month)
  }
  if (scale == "hourly") {
    # Code for hourly scale
    warning("Caution! Not efficient!\n
    Hourly mode takes a while to compute population variation.\n
    May throw error if covariates are not on hourly scale.")
    dat <- dat |>
      group_by(year, month, day, hour) |>
      summarize(covariate = mean(covariate, na.rm = TRUE)) |>
      ungroup() |>
      mutate(time = year +
        yhour(make_datetime(year, month, day, hour)) /
          hours_in_year(year))
    if (!is.null(extra_cov_dataset)) {
      extra_cov_dataset <- extra_cov_dataset |>
        group_by(year, month, day, hour) |>
        summarize(across(all_of(extra_cov), ~ mean(.x, na.rm = TRUE))) |>
        ungroup() |>
        mutate(time = year +
          yhour(make_datetime(year, month, day, hour)) /
            hours_in_year(year)) |>
        rename(!!!new_names)
    }
  } else if (scale == "daily") {
    # Code for daily scale
    dat <- dat |>
      group_by(year, month, day) |>
      summarize(covariate = mean(covariate, na.rm = TRUE)) |>
      ungroup() |>
      mutate(time = year + yday(make_datetime(year, month, day)) /
        days_in_year(year))
    if (!is.null(extra_cov_dataset)) {
      extra_cov_dataset <- extra_cov_dataset |>
        group_by(year, month, day) |>
        summarize(across(all_of(extra_cov), ~ mean(.x, na.rm = TRUE))) |>
        ungroup() |>
        mutate(time = year + yday(make_datetime(year, month, day)) /
          days_in_year(year)) |>
        rename(!!!new_names)
    }
  } else if (scale == "monthly") {
    # Code for monthly scale
    dat <- dat |>
      group_by(year, month) |>
      summarize(covariate = mean(covariate, na.rm = TRUE)) |>
      ungroup() |>
      mutate(time = year + month / 12)
    if (!is.null(extra_cov_dataset)) {
      extra_cov_dataset <- extra_cov_dataset |>
        group_by(year, month) |>
        summarize(across(all_of(extra_cov), ~ mean(.x, na.rm = TRUE))) |>
        ungroup() |>
        mutate(time = year + month / 12) |>
        rename(!!!new_names)
    }
  } else {
    stop("Invalid scale. Please use 'hourly', 'daily', or 'monthly'.")
  }
  if (!is.null(extra_cov_dataset)) {
    dat <- dat |> left_join(extra_cov_dataset)
  }

  scalers <- NULL
  if (normalize_covariate || normalize_extra) {
    cols_to_scale <- c(
      if (normalize_covariate) "covariate" else NULL,
      if (normalize_extra) paste0("cov", seq_along(extra_cov)) else NULL
    )
    if (length(cols_to_scale)) {
      scalers <- caret::preProcess(dat[, cols_to_scale, drop = FALSE], method = "range")
      dat[, cols_to_scale] <- predict(scalers, dat[, cols_to_scale, drop = FALSE])
    }
  }

  dat <- dat |>
    left_join(dat_mal |> select(year, month, pop))
  spl_res <- safe_pop_smooth(dat$time, dat$pop, method = "smooth.spline")
  dat <- dat |>
    mutate(
      pop = spl_res$pop,
      dpopdt = spl_res$dpopdt
    )
  dat_mal <- dat_mal |>
    filter_by_period(start_year, start_month, end_year, end_month)
  result <- list(dat = dat, dat_mal = dat_mal)
  return(result)
}

#' Prepare a pomp object from a model definition file
#'
#' This helper loads a model definition file `pomp_object_<model>.R` (or a provided
#' `model_file`) into a private environment, validates the presence of required
#' symbols, verifies the covariate / observation data contract, and constructs
#' a `pomp` object ready for fitting or simulation. The function avoids
#' polluting the global environment by sourcing the model file into an internal
#' environment and referencing objects explicitly from there.
#'
#' @param model Character. Short model name used to find `pomp_object_<model>.R`
#'   inside `model_dir`. Ignored if `model_file` is provided.
#' @param param Named numeric vector of parameter starting values passed to
#'   `pomp::pomp`.
#' @param dat_cov Data frame. Covariates table produced by `preprocess_data`.
#'   Must contain `time` (numeric), `pop`, and `dpopdt` columns and at least two
#'   timepoints.
#' @param dat_mal Data frame. Malaria observations. Must contain `time` and
#'   `PF` columns.
#' @param model_dir Character. Directory where `pomp_object_<model>.R` is
#'   looked up (default: `dynamical_model/fitting`).
#' @param model_file Optional character. Absolute path to a model file. If
#'   provided, `model`/`model_dir` are ignored.
#' @param validate_only Logical. If `TRUE` perform validations and return a
#'   small list with `status = 'ok'` rather than building the full `pomp`
#'   object. Useful for CI and dry runs.
#'
#' @return A `pomp` object when `validate_only = FALSE`, or a list with
#'   `status = 'ok'` when `validate_only = TRUE`.
#'
#' @details The model file must define the following objects: `par_names`,
#' `vp_names`, `simul_pf`, `rmeas_pf`, `dmeas_pf`, `initlz_pf`, `state_names`,
#' and the transformation helpers `log_transf`, `logit_transf`, and
#' `barycentric_transf`.
#'
prepare_pomp_model <- function(
    model, param, dat_cov, dat_mal,
    model_dir = "dynamical_model/fitting",
    model_file = NULL,
    validate_only = FALSE) {
  # Resolve model file path
  if (is.null(model_file)) {
    model_file <- file.path(model_dir, paste0("pomp_object_", model, ".R"))
  }
  if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file, ". Provide correct 'model' or absolute 'model_file'.")
  }

  # Source into private environment to avoid global pollution
  model_env <- new.env(parent = baseenv())
  sys.source(model_file, envir = model_env)

  # Required symbols
  required_symbols <- c(
    "par_names", "vp_names",
    "simul_pf", "rmeas_pf", "dmeas_pf", "initlz_pf",
    "state_names",
    "log_transf", "logit_transf", "barycentric_transf"
  )
  missing_syms <- required_symbols[!vapply(required_symbols, exists, logical(1), envir = model_env, inherits = FALSE)]
  if (length(missing_syms) > 0) {
    stop("Model file '", model_file, "' is missing required symbols: ", paste(missing_syms, collapse = ", "))
  }

  # Validate dat_cov contract
  req_cov_cols <- c("time", "pop", "dpopdt")
  if (!all(req_cov_cols %in% colnames(dat_cov))) {
    stop("`dat_cov` must contain columns: ", paste(req_cov_cols, collapse = ", "), " .")
  }
  if (!is.numeric(dat_cov$time)) {
    stop("`dat_cov$time` must be numeric fractional-year time index.")
  }
  if (length(dat_cov$time) < 2) {
    stop("`dat_cov` must contain at least two timepoints for covariate lookup/extrapolation.")
  }
  # Ensure increasing order of time
  if (is.unsorted(dat_cov$time)) {
    warning("`dat_cov$time` appears unsorted: sorting by time.")
    dat_cov <- dat_cov[order(dat_cov$time), , drop = FALSE]
  }
  if (any(duplicated(dat_cov$time))) {
    stop("`dat_cov$time` contains duplicated time entries. Resolve duplicates before calling prepare_pomp_model.")
  }

  # Validate dat_mal contract
  if (!all(c("time", "PF") %in% colnames(dat_mal))) {
    stop("`dat_mal` must contain columns 'time' and 'PF'.")
  }

  # Validate params against model parameter names
  model_paramnames <- c(get("par_names", envir = model_env), get("vp_names", envir = model_env))
  if (!is.null(names(param))) {
    missing_params <- setdiff(model_paramnames, names(param))
    if (length(missing_params) > 0) {
      warning(
        "Provided `param` does not include some model parameters: ", paste(missing_params, collapse = ", "),
        ". pomp() may still accept but fitting could fail."
      )
    }
  }

  if (validate_only) {
    return(list(status = "ok", model_file = normalizePath(model_file)))
  }

  # Build covariate table for pomp
  covariates <- covariate_table(
    dat_cov |> select(time, pop, dpopdt, contains("cov")),
    times = "time"
  )

  # create extrapolation times vector; safe because dat_cov has >=2 sorted timepoints
  t_extrap <- with(dat_cov, c(2 * time[1] - time[2], time))
  covariates <- repair_lookup_table(covariates, t_extrap)

  # assemble pomp object while explicitly referencing functions/objects in model_env
  po <- pomp(
    data = dat_mal |> select(time, PF),
    params = param,
    times = "time",
    t0 = with(dat_cov, 2 * time[1] - time[2]),
    covar = covariates,
    rprocess = euler(step.fun = model_env$simul_pf, delta.t = 1 / (365)),
    rmeasure = model_env$rmeas_pf,
    dmeasure = model_env$dmeas_pf,
    rinit = model_env$initlz_pf,
    paramnames = model_paramnames,
    partrans = parameter_trans(
      log = model_env$log_transf,
      logit = model_env$logit_transf,
      barycentric = model_env$barycentric_transf
    ),
    accumvars = c("cases"),
    statenames = model_env$state_names
  )

  # attach the model environment for debugging
  attr(po, "model_env") <- model_env

  return(po)
}

#' Generate a `pomp::rw_sd` specification from a model file
#'
#' Discover model parameter names from a `pomp_object_<model>.R` file (or an
#' explicit `model_file`) and construct a `pomp::rw_sd` object programmatically.
#'
#' The function creates regular random-walk standard-deviation entries for
#' estimated parameters (`par_names`) and `ivp(...)` entries for initial
#' condition parameters (`vp_names`). Call objects (language) are used for
#' `ivp(...)` entries so no `eval(parse(...))` is required.
#'
#' @param reg_rw Numeric. Default/random-walk standard deviation to apply to
#'   regular parameters (e.g. 0.03).
#' @param ivp_rw Numeric. Standard deviation to apply to initial condition
#'   parameters (used inside an `ivp(...)` call, e.g. 0.1).
#' @param model Optional character. Short model name used to find
#'   `pomp_object_<model>.R` in `model_dir`. Ignored if `model_file` is
#'   provided.
#' @param model_dir Character. Directory where model files live (default:
#'   `dynamical_model/fitting`).
#' @param model_file Optional character. Absolute path to a model file. If
#'   provided this is used instead of `model`/`model_dir`.
#' @param extra_rw Optional named vector or list of overrides. Names are
#'   parameter names. Values may be numeric (regular rw sd) or language/call
#'   objects (for example `call("ivp", 0.05)`). Setting a value to `NA`
#'   removes the corresponding rw entry.
#' @param profile_var Optional character or vector of parameter names to
#'   exclude from the returned `rw.sd` (useful when profiling parameters).
#' @param remove_delta Logical. If `TRUE` remove `delta` from the rw list (we
#'   generally do not fit human mortality in these models). Default: `TRUE`.
#'
#' @return An object of class `rw.sd` (as produced by `pomp::rw_sd`).
#'
#' @export
generate_rw_sd <- function(reg_rw, ivp_rw, model = NULL, model_dir = "dynamical_model/fitting", model_file = NULL,
                           extra_rw = NULL, profile_var = NULL, remove_delta = TRUE) {
  # Resolve model file path
  if (is.null(model_file)) {
    if (is.null(model)) stop("Either 'model' or 'model_file' must be provided to generate_rw_sd()")
    model_file <- file.path(model_dir, paste0("pomp_object_", model, ".R"))
  }
  if (!file.exists(model_file)) stop("Model file not found: ", model_file)

  # Source into private environment and discover parameter names
  model_env <- new.env(parent = baseenv())
  sys.source(model_file, envir = model_env)

  if (!exists("par_names", envir = model_env) || !exists("vp_names", envir = model_env)) {
    stop("Model file must define 'par_names' and 'vp_names' vectors for discovery")
  }

  par_names <- get("par_names", envir = model_env)
  vp_names <- get("vp_names", envir = model_env)

  # Build default rw.sd entries: regular rw for par_names, ivp(...) for vp_names
  entries <- list()
  for (p in par_names) {
    entries[[p]] <- as.numeric(reg_rw)
  }
  for (v in vp_names) {
    # create a call object ivp(ivp_rw) so we avoid parsing strings later
    entries[[v]] <- as.call(list(as.name("ivp"), as.numeric(ivp_rw)))
  }

  # Apply extra_rw overrides (named vector or list); allow adding or overriding entries
  if (!is.null(extra_rw)) {
    if (is.null(names(extra_rw)) || any(names(extra_rw) == "")) {
      stop("'extra_rw' must be a named vector or list of overrides: names -> sd values")
    }
    for (nm in names(extra_rw)) {
      val <- extra_rw[[nm]]
      if (is.na(val)) {
        # explicit NA means remove this rw from entries
        entries[[nm]] <- NULL
      } else if (is.numeric(val)) {
        entries[[nm]] <- as.numeric(val)
      } else if (is.language(val) || is.call(val)) {
        entries[[nm]] <- val
      } else {
        stop("'extra_rw' values must be numeric or language/call objects (e.g. call('ivp', 0.05))")
      }
    }
  }
  if (remove_delta) {
    # Remove delta if present, we dont fit human mortality rate
    entries[["delta"]] <- NULL
  }

  # Remove profiled variable if requested
  if (!is.null(profile_var)) {
    entries[profile_var] <- NULL
  }

  if (length(entries) == 0) {
    stop("No rw.sd entries remain after applying overrides/profile exclusion")
  }

  # call pomp::rw_sd using do.call with numeric values and call objects (no parse/eval)
  args_list <- entries
  new_rdd <- do.call(pomp::rw_sd, args_list)

  return(new_rdd)
}

setClass("mif2_list", contains = "list")
setClass("pfilter_list", contains = "list")

# mif2_try ----
setGeneric(
  "mif2_try",
  function(data, starts, ...) {
    standardGeneric("mif2_try")
  },
  valueClass = "mif2_list"
)

setMethod(
  "mif2_try",
  signature = signature(data = "ANY", starts = "data.frame"),
  definition = function(data, starts, ...,
                        seed = TRUE, chunk.size = NULL, scheduling = 1, verbose = FALSE) {
    foreach(
      iter_i = seq_len(nrow(starts)),
      .combine = c,
      .options.future = list(
        seed = seed,
        chunk.size = chunk.size,
        scheduling = scheduling
      )
    ) %dofuture% {
      res2 <- tryCatch(
        {
          pomp::mif2(data, params = starts[iter_i, ], ...)
        },
        error = function(e) {
          if (verbose) {
            warning(sprintf("mif2 failed for iter %d: %s", iter_i, conditionMessage(e)))
          }
          NULL
        }
      )
      list(res2)
    } -> res

    names(res) <- row.names(starts)

    # determine which mif samples actually completed (non-NULL entries)
    success_flags <- vapply(res, function(x) !is.null(x), logical(1))

    if (all(!success_flags)) {
      warning("All mif2 samples failed")
    }

    attr(res, "success") <- success_flags
    class(res) <- "mif2_list"
    return(res)
  }
)

setClassUnion("pomp_list", c("pfilter_list", "mif2_list"))


setGeneric(
  "pfilter_try",
  function(data, Nrep, ...) {
    standardGeneric("pfilter_try")
  }
)

# pfilter_try ----
setMethod(
  "pfilter_try",
  signature = signature(data = "pomp_list", Nrep = "numeric"),
  definition = function(data, Nrep, ...,
                        seed = TRUE, chunk.size = NULL, scheduling = 1, verbose = FALSE) {
    npo <- length(data)
    njobs <- Nrep * npo

    foreach(
      iter_i = seq_len(njobs),
      .combine = c,
      .options.future = list(
        seed = seed,
        chunk.size = chunk.size,
        scheduling = scheduling
      )
    ) %dofuture% {
      ipo <- (iter_i - 1) %% npo + 1
      rep <- (iter_i - 1) %/% npo + 1
      if (is.null(data[[ipo]])) {
        if (verbose) {
          warning(sprintf("Skipping pfilter for iter %d rep %d: NULL pomp object", ipo, rep))
        }
        pf_res <- NULL
      } else {
        pf_res <- tryCatch(
          {
            pomp::pfilter(data[[ipo]], ...)
          },
          error = function(e) {
            if (verbose) {
              warning(sprintf("pfilter failed for iter %d rep %d: %s", ipo, rep, conditionMessage(e)))
            }
            NULL
          }
        )
      }
      list(pf_res)
    } -> res
    # split into list of lists by original data index
    res <- split(res, rep(seq_len(npo), Nrep, length.out = length(res)))

    class(res) <- "pfilter_list"
    return(res)
  }
)

# pfilter_try ----
setMethod(
  "pfilter_try",
  signature = signature(data = "pomp", Nrep = "numeric"),
  definition = function(data, Nrep, ...,
                        seed = TRUE, chunk.size = NULL, scheduling = 1, verbose = FALSE) {
    njobs <- Nrep

    foreach(
      iter_i = seq_len(njobs),
      .combine = c,
      .options.future = list(
        seed = seed,
        chunk.size = chunk.size,
        scheduling = scheduling
      )
    ) %dofuture% {
      if (is.null(data)) {
        if (verbose) {
          warning(sprintf("Skipping pfilter for rep %d: NULL pomp object", iter_i))
        }
        pf_res <- NULL
      } else {
        pf_res <- tryCatch(
          {
            pomp::pfilter(data, ...)
          },
          error = function(e) {
            if (verbose) {
              warning(sprintf("pfilter failed for rep %d: %s", iter_i, conditionMessage(e)))
            }
            NULL
          }
        )
      }
      list(pf_res)
    } -> res
    # split into list of lists by original data index
    # res <- split(res, rep(seq_len(npo), Nrep, length.out = length(res)))

    class(res) <- "pfilter_list"
    return(res)
  }
)


# logLik_try ----
setGeneric(
  "logLik_try",
  function(data, ...) {
    standardGeneric("logLik_try")
  }
)

setMethod(
  "logLik_try",
  signature = signature(data = "pomp_list"),
  definition = function(data) {
    n_outer <- length(data)
    loglik_list <- vector("list", n_outer)
    for (iter_i in seq_len(n_outer)) {
      if (is.null(data[[iter_i]])) {
        loglik_list[[iter_i]] <- NA_real_
        next
      }
      loglik_vals <- vapply(data[[iter_i]], function(pf) {
        if (is.null(pf)) {
          NA_real_
        } else {
          tryCatch(
            {
              ll <- logLik(pf)
              as.numeric(ll)
            },
            error = function(e) {
              warning(sprintf("logLik failed for iter %d: %s", iter_i, conditionMessage(e)))
              NA_real_
            }
          )
        }
      }, numeric(1))
      loglik_list[[iter_i]] <- loglik_vals
      loglik_list_mat <- do.call(rbind, loglik_list)

      attr(data[[iter_i]], "loglik") <- loglik_vals
    }
    return(loglik_list_mat)
  }
)

setGeneric(
  "coef_try",
  function(data, ...) {
    standardGeneric("coef_try")
  }
)

setMethod(
  "coef_try",
  signature = signature(data = "pomp_list"),
  definition = function(data) {
    # If data is a mif2_list: row per element of list
    if (inherits(data, "mif2_list")) {
      n <- length(data)
      param_names <- NULL
      # Find first non-NULL entry for param names
      for (i in seq_len(n)) {
        if (!is.null(data[[i]])) {
          param_names <- names(coef(data[[i]]))
          break
        }
      }
      if (is.null(param_names)) {
        warning("No valid coefficients found in mif2_list")
        return(matrix(NA_real_, nrow = n, ncol = 0))
      }
      coef_mat <- matrix(NA_real_, nrow = n, ncol = length(param_names))
      colnames(coef_mat) <- param_names
      for (i in seq_len(n)) {
        if (!is.null(data[[i]])) {
          vals <- tryCatch(
            {
              cf <- coef(data[[i]])
              out <- rep(NA_real_, length(param_names))
              names(out) <- param_names
              out[names(cf)] <- as.numeric(cf)
              out
            },
            error = function(e) {
              warning(sprintf("coef failed for mif2_list element %d: %s", i, conditionMessage(e)))
              rep(NA_real_, length(param_names))
            }
          )
          coef_mat[i, ] <- vals
        }
      }
      return(coef_mat)
    }
    # If data is a pfilter_list: get coef from first element of nested list
    if (inherits(data, "pfilter_list")) {
      n_outer <- length(data)
      param_names <- NULL
      # Find first non-NULL entry for param names
      for (i in seq_len(n_outer)) {
        if (length(data[[i]]) > 0 && !is.null(data[[i]][[1]])) {
          param_names <- names(coef(data[[i]][[1]]))
          break
        }
      }
      if (is.null(param_names)) {
        warning("No valid coefficients found in pfilter_list")
        return(matrix(NA_real_, nrow = n_outer, ncol = 0))
      }
      coef_mat <- matrix(NA_real_, nrow = n_outer, ncol = length(param_names))
      colnames(coef_mat) <- param_names
      for (i in seq_len(n_outer)) {
        if (length(data[[i]]) > 0 && !is.null(data[[i]][[1]])) {
          vals <- tryCatch(
            {
              cf <- coef(data[[i]][[1]])
              out <- rep(NA_real_, length(param_names))
              names(out) <- param_names
              out[names(cf)] <- as.numeric(cf)
              out
            },
            error = function(e) {
              warning(sprintf("coef failed for pfilter_list element %d: %s", i, conditionMessage(e)))
              rep(NA_real_, length(param_names))
            }
          )
          coef_mat[i, ] <- vals
        }
      }
      return(coef_mat)
    }
    stop("Unsupported data type for coef_try: must be 'mif2_list' or 'pfilter_list'")
  }
)





#' Run MIF2 fitting and post-process results
#'
#' Low-level worker that executes `circumstance::mif2` on a prepared `pomp`
#' object and post-processes per-sample `pfilter` + `coef` summaries. The
#' function preserves and restores the user's `future` plan when switching
#' between parallel and sequential execution.
#'
#' @param po A `pomp` object prepared by `prepare_pomp_model()`.
#' @param n_cores Integer number of worker processes to use when `allow_parallel = TRUE`.
#' @param parameters A list of starting parameter vectors (format expected by `circumstance::mif2`).
#' @param seed_num Integer RNG seed (passed externally for reproducibility).
#' @param rdd A `pomp::rw_sd` object defining parameter random-walk sizes.
#' @param allow_parallel Logical; if `TRUE` will set a parallel future plan (multicore).
#' @param start_index Integer offset added to sample indices (useful for batched grids).
#' @param output_to_file Logical; when `TRUE` results are appended to `output_result`.
#' @param output_result Path to the CSV file to append results.
#' @param output_log Path to a log file (currently unused but reserved for future logging).
#' @param output_format Character; currently `"data.frame"` is supported.
#'
#' @return A data.frame with fitted parameter estimates, log-likelihood summaries and a `flag` column indicating successful evaluations (1) or failures (0).
#'
#' @details The function internally calls `circumstance::mif2` and then runs
#' `circumstance::pfilter` for each MIF sample to estimate log-likelihoods. It
#' uses `tryCatch` around heavy operations so single-sample failures become NA
#' rows instead of aborting the whole batch.
#'
#' @export
run_fitting <- function(po, n_cores, parameters,
                        seed_num, rdd, allow_parallel,
                        start_index,
                        output_to_file, output_result,
                        output_log, output_format,
                        mif_control = list(Np = 3000, Nmif = 50, cooling.type = "geometric", cooling.fraction.50 = 0.5),
                        pfilter_control = list(Nrep = 10, Np = 10000),
                        verbose = FALSE) {
  # preserve and restore the caller's future plan even on error/exit
  old_plan <- tryCatch(
    {
      future::plan()
    },
    error = function(e) NULL
  )
  on.exit(
    {
      if (!is.null(old_plan)) {
        tryCatch(
          {
            future::plan(old_plan)
          },
          error = function(e) NULL
        )
      }
    },
    add = TRUE
  )

  if (is.null(mif_control)) mif_control <- list(Np = 3000, Nmif = 50, cooling.type = "geometric", cooling.fraction.50 = 0.5)
  if (is.null(pfilter_control)) pfilter_control <- list(Nrep = 10, Np = 10000)

  # choose execution plan
  if (allow_parallel) {
    ## if parallel register cluster
    plan(multicore, workers = n_cores)
  } else {
    plan(sequential)
  }

  if (!is.null(seed_num)) {
    mif_args <- c(
      list(data = po), mif_control,
      list(starts = parameters, rw.sd = rdd, seed = seed_num, verbose = verbose)
    )
  } else {
    mif_args <- c(
      list(data = po), mif_control,
      list(starts = parameters, rw.sd = rdd, seed = TRUE, verbose = verbose)
    )
  }

  mifout <- tryCatch(
    {
      do.call(mif2_try, mif_args)
    },
    error = function(e) {
      if (verbose || is.null(verbose)) {
        warning("mif2 failed: ", conditionMessage(e))
      }
      return(NULL)
    }
  )

  # extract success flags from mifout attribute
  success_flags <- attr(mifout, "success")
  if (!any(success_flags)) {
    if (verbose || is.null(verbose)) {
      warning("All mif2 samples failed")
    }
    return(data.frame())
  }

  # If mifout is NULL or empty, nothing to process
  if (is.null(mifout) || length(mifout) == 0) {
    warning("No successful mif2 outputs to process")
    return(data.frame())
  }

  pfilter_args <- c(
    list(data = mifout), pfilter_control, list(seed = TRUE, verbose = verbose)
  )

  pfout <- tryCatch(
    {
      do.call(pfilter_try, pfilter_args)
    },
    error = function(e) {
      if (verbose || is.null(verbose)) {
        warning("pfilter failed: ", conditionMessage(e))
      }
      return(NULL)
    }
  )

  if (is.null(pfout)) {
    warning("No successful pfilter outputs to process")
    return(data.frame())
  }

  # print(pfout)
  loglik_matrix <- logLik_try(pfout)
  coef_matrix <- coef_try(mifout)
  # compute per-outer-loglik summaries (mean and se). Handle single-row result from apply()
  loglik_raw <- apply(loglik_matrix, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) > 0) {
      pomp::logmeanexp(x, se = TRUE)
    } else {
      c(NA_real_, NA_real_)
    }
  })

  if (is.null(dim(loglik_raw))) {
    # single outer iteration -> loglik_raw is a named numeric vector (mean, se)
    nm <- names(loglik_raw)
    if (is.null(nm)) nm <- c("mean", "se")
    loglik <- matrix(loglik_raw, nrow = length(loglik_raw), ncol = 1)
    rownames(loglik) <- nm
  } else {
    loglik <- loglik_raw
    # ensure rownames exist (some older R versions may drop names)
    if (is.null(rownames(loglik))) {
      rn <- names(loglik[, 1])
      if (!is.null(rn)) rownames(loglik) <- rn
    }
  }
  result <- dplyr::as_tibble(coef_matrix)
  result$sample <- seq_len(nrow(result)) + start_index - 1
  result$loglik <- as.numeric(loglik[1, ])
  result$loglik.se <- as.numeric(loglik[2, ])
  result$flag <- as.integer(success_flags)




  return(result)
}

## helper that encapsulates data loading + pomp preparation
#' Prepare data and a POMP model from CSV inputs
#'
#' Convenience wrapper that reads malaria and covariate CSV files, runs
#' `preprocess_data()` to produce covariate and observation tables, coerces
#' `parameters` into a named numeric vector, and calls
#' `prepare_pomp_model()` to build the final `pomp` object.
#'
#' This helper is intended for interactive and programmatic workflows where
#' the raw inputs are CSV files (paths) rather than already-loaded data
#' frames.
#'
#' @param parameters Named numeric vector or single-row data.frame/tibble with
#'   starting parameter values (when a data.frame the first row is used).
#' @param dataset_mal Path to CSV file with malaria case data (will be read
#'   using `readr::read_csv`). The resulting data must contain `time` or
#'   `time_end` and `PF` after preprocessing.
#' @param dataset_VC Path to CSV file with covariate data (VC/climate).
#' @param extra_cov_dataset Optional path to a CSV containing additional
#'   covariates (if not present in `dataset_VC`).
#' @param extra_cov Optional character vector of additional covariate column
#'   names to include from `extra_cov_dataset` or `dataset_VC`.
#' @param scale Character. Aggregation scale: one of `'monthly'`, `'daily'`, or
#'   `'hourly'` (passed to `preprocess_data`).
#' @param covariate Character. Column name in the covariate CSV to use as the
#'   primary covariate.
#' @param start_year,end_year,start_month,end_month Integers specifying the
#'   inclusive time window to prepare.
#' @param normalize_extra Logical. Whether to normalize extra covariates with
#'   a safe min-max transform.
#' @param model Character short name used to locate `pomp_object_<model>.R`
#'   inside `model_dir` (ignored when `model_file` is provided).
#' @param model_dir Directory where model definition files live.
#' @param model_file Optional absolute path to a `pomp_object_*.R` file to
#'   source instead of using `model` + `model_dir`.
#' @param validate_only Logical; if `TRUE` performs validations and returns a
#'   small list with `status = 'ok'` rather than building the full `pomp`
#'   object. Useful for CI or dry-run checks.
#'
#' @return A list with elements:
#' - `po`: a `pomp` object (or a small validation list when `validate_only=TRUE`)
#' - `dat`: preprocessed covariate table (includes `time`, `pop`, `dpopdt`)
#' - `dat_mal`: preprocessed malaria observation table with columns `time` and `PF`.
#'
#' @seealso `preprocess_data`, `prepare_pomp_model`, `fit_model`
prepare_model <- function(parameters,
                          dataset_mal, dataset_VC,
                          extra_cov_dataset = NULL, extra_cov = NULL,
                          scale, covariate,
                          start_year, end_year, start_month, end_month,
                          normalize_covariate = TRUE,
                          normalize_extra = TRUE,
                          model, model_dir = "dynamical_model/fitting",
                          model_file = NULL, validate_only = FALSE) {
  data_mal <- readr::read_csv(dataset_mal, show_col_types = FALSE)
  data_cov <- readr::read_csv(dataset_VC, show_col_types = FALSE)
  if (!is.null(extra_cov_dataset)) {
    data_extra_cov <- readr::read_csv(extra_cov_dataset, show_col_types = FALSE)
  } else {
    data_extra_cov <- NULL
  }

  all_dat <- preprocess_data(
    data_cov, data_mal, scale, covariate, start_year,
    end_year, start_month, end_month, extra_cov, data_extra_cov,
    normalize_covariate, normalize_extra
  )

  dat <- all_dat[[1]]
  dat_mal_raw <- all_dat[[2]]

  # robust time column selection (accept time or time_end)
  if ("time_end" %in% colnames(dat_mal_raw)) {
    dat_mal <- dat_mal_raw |>
      dplyr::select(time_end, PF) |>
      dplyr::rename(time = time_end)
  } else if ("time" %in% colnames(dat_mal_raw)) {
    dat_mal <- dat_mal_raw |>
      dplyr::select(time, PF)
  } else {
    stop("`dat_mal` produced by preprocess_data must contain 'time' or 'time_end' and 'PF'.")
  }

  ### prepare pomp model - coerce parameters to named numeric vector
  if (tibble::is_tibble(parameters) || is.data.frame(parameters)) {
    param <- as.numeric(parameters[1, ])
    param_names <- colnames(parameters)
    names(param) <- param_names
  } else { ## expects named vector
    param <- parameters
  }

  po <- prepare_pomp_model(model, param, dat, dat_mal, model_dir = model_dir, model_file = model_file, validate_only = validate_only)

  list(po = po, dat = dat, dat_mal = dat_mal)
}


#' Fit a POMP model over a parameter grid or starting parameter vector
#'
#' High-level wrapper that prepares data, constructs a `pomp` model from a
#' model definition file, creates random-walk specifications, and runs the
#' fitting routine. This function provides convenience wiring for experiments
#' and supports fast "validate-only" checks useful for CI.
#'
#' @param parameters Named numeric vector or single-row data.frame/tibble with starting parameter values.
#' @param city Character. City name used to generate default output filenames.
#' @param model Character. Short model name (used to locate `pomp_object_<model>.R` in `model_dir`).
#' @param scale Character. Aggregation scale: `'monthly'`, `'daily'`, or `'hourly'`.
#' @param covariate Character. Covariate column name to use (e.g. `'hadISD_VC'`).
#' @param start_year Integer. First year to include.
#' @param start_month Integer. Starting month in `start_year`.
#' @param end_year Integer. Last year to include.
#' @param end_month Integer. Final month to include.
#' @param dataset_mal Path to CSV with malaria case data.
#' @param dataset_VC Path to CSV with covariate data (VC/climate).
#' @param mode Character. One of `'fit'`, `'refined'`, `'refined_second'` controlling default rw sizes.
#' @param allow_parallel Logical. Whether to allow parallel execution (via futures).
#' @param n_cores Integer. Number of cores/workers to request when parallel.
#' @param seed_num Integer. RNG seed for runs.
#' @param start_index Integer. Index offset for sample numbering (useful when running batched grids).
#' @param output_to_file Logical. Whether to append results to disk.
#' @param output_result Path to result CSV file. If missing a default is generated.
#' @param output_log Path to log file. If missing a default is generated.
#' @param output_format Character. `'data.frame'` or `'pomps'` for simulation outputs.
#' @param profile_var Optional character. Parameter name(s) to exclude from `rw.sd` (for profiling experiments).
#' @param extra_cov Optional character vector of additional covariate names to merge into the covariate dataset.
#' @param extra_cov_dataset Optional path to a CSV containing `extra_cov` columns (if not present in `dataset_VC`).
#' @param normalize_extra Logical. Whether to min-max normalize extra covariates.
#' @param model_file Optional character. Absolute path to a `pomp_object_*.R` model file (overrides `model` + `model_dir`).
#' @param extra_rw Optional named list or vector of rw overrides passed to `generate_rw_sd()` (names -> sd or call objects). Set a value to `NA` to remove that entry.
#' @param validate_only Logical. If `TRUE` perform validations and return early with `list(status = 'ok')` instead of constructing the full `pomp` object and fitting.
#' @param reg_rw Optional numeric. Random-walk SD to use for regular (non-ivp) parameters. If `NULL` defaults are chosen based on `mode`.
#' @param ivp_rw Optional numeric. Random-walk SD to use for initial-condition parameters. If `NULL` defaults are chosen based on `mode`.
#' @param mif_control Optional list of control arguments forwarded to `circumstance::mif2` (e.g. `list(Np=1000, Nmif=50)`). If `NULL` defaults are used.
#' @param pfilter_control Optional list of control arguments forwarded to `circumstance::pfilter` (e.g. `list(Nrep=10, Np=10000)`). If `NULL` defaults are used.
#' @param manifest_append Logical. If `TRUE` and a manifest file already exists at the target path, the new run info will be appended into an array; if `FALSE` the manifest file will be (re)created/overwritten. Default: `FALSE`.
#'
#' @return A data.frame (or list) with fitted parameter estimates and log-likelihood information when `validate_only = FALSE`.
#' When `validate_only = TRUE` a small list with `status = 'ok'` and `model_file` is returned.
#'
#' @export
fit_model <- function(
    parameters, city,
    model = "const", scale = "monthly", covariate = "hadISD_VC",
    start_year = 1997, start_month = 1,
    end_year = 2014, end_month = 12,
    dataset_mal, dataset_VC, mode = "fit", allow_parallel = TRUE,
    n_cores = 1, seed_num = as.integer(Sys.time()),
    start_index = 1, output_to_file = TRUE,
    output_result = NULL, output_log = NULL, output_format = "data.frame",
    profile_var = NULL, extra_cov = NULL,
    extra_cov_dataset = NULL,
    normalize_covariate = TRUE,
    normalize_extra = TRUE,
    model_file = NULL,
    model_dir = "dynamical_model/fitting",
    extra_rw = NULL, validate_only = FALSE,
    reg_rw = NULL, ivp_rw = NULL,
    mif_control = NULL, pfilter_control = NULL,
    manifest_append = FALSE, manifest_path = NULL, verbose = FALSE) {
  # call the helper to load data and prepare the pomp object
  prep <- prepare_model(
    parameters = parameters,
    dataset_mal = dataset_mal,
    dataset_VC = dataset_VC,
    extra_cov_dataset = extra_cov_dataset,
    extra_cov = extra_cov,
    scale = scale,
    covariate = covariate,
    start_year = start_year,
    end_year = end_year,
    start_month = start_month,
    end_month = end_month,
    normalize_covariate = normalize_covariate,
    normalize_extra = normalize_extra,
    model = model,
    model_dir = model_dir,
    model_file = model_file,
    validate_only = validate_only
  )

  po <- prep$po
  dat <- prep$dat
  dat_mal <- prep$dat_mal


  if (validate_only) {
    return(po)
  }

  # create run manifest (start)
  if (is.null(manifest_path)) {
    run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", sample(1e6, 1))
    manifest_path <- file.path("results", paste0("run_", run_id, ".json"))
  } else {
    run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", sample(1e6, 1))
    manifest_path <- file.path(manifest_path)
  }

  if (output_to_file) {
    # names_files <- generate_output_names(city, model, scale, covariate, start_year, start_month, end_year, end_month, mode, profile_var)
    if (rlang::is_missing(output_result) || is.null(output_result)) output_result <- file.path("results", paste0("run_", run_id, ".csv"))
    if (rlang::is_missing(output_log) || is.null(output_log)) output_log <- file.path("results", paste0("run_", run_id, ".log"))
  }

  manifest_args <- list(
    city = city, model = model, scale = scale, covariate = covariate,
    start_year = start_year, start_month = start_month, end_year = end_year, end_month = end_month,
    mode = mode, n_cores = n_cores, seed_num = seed_num, start_index = start_index,
    output_to_file = output_to_file, output_result = output_result, output_log = output_log,
    profile_var = profile_var, extra_cov = extra_cov
  )
  manifest_info <- create_manifest_info(run_id = run_id, start_time = Sys.time(), user = Sys.info()["user"], args = manifest_args)
  manifest_info$parameters <- summarize_parameters(parameters)
  # write initial manifest
  # tryCatch(write_run_manifest(manifest_path, manifest_info, append = manifest_append), error = function(e) warning("Failed to write run manifest: ", conditionMessage(e)))


  if (is.null(reg_rw)) {
    reg_rw <- switch(mode,
      fit = 0.03,
      refined = 0.015,
      refined_second = 0.001
    )
  }
  if (is.null(ivp_rw)) {
    ivp_rw <- switch(mode,
      fit = 0.1,
      refined = 0.05,
      refined_second = 0.01
    )
  }
  rdd <- generate_rw_sd(reg_rw, ivp_rw, model = model, model_dir = model_dir, model_file = model_file, extra_rw = extra_rw, profile_var = profile_var)

  # ### fit model
  result <- run_fitting(
    po, n_cores, parameters,
    seed_num, rdd, allow_parallel,
    start_index,
    output_to_file, output_result,
    output_log, output_format,
    mif_control = mif_control, pfilter_control = pfilter_control
  )

  result <- result |> filter(flag == 1)


  # update manifest with results summary
  try(
    {
      min_vals <- sapply(result |> select(-any_of(c("loglik", "loglik.se", "sample", "flag"))), min, na.rm = TRUE)
      max_vals <- sapply(result |> select(-any_of(c("loglik", "loglik.se", "sample", "flag"))), max, na.rm = TRUE)
      param_summary <- setNames(
        lapply(seq_along(min_vals), function(i) c(min_vals[i], max_vals[i])),
        names(min_vals)
      )

      manifest_info$end_time <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
      manifest_info$status <- if (is.data.frame(result) && nrow(result) > 0) "finished" else "failed"
      manifest_info$summary <- list(
        n_results = if (is.data.frame(result)) nrow(result) else 0,
        param_ranges = param_summary
      )
      # write_run_manifest(manifest_path, manifest_info, append = manifest_append)
    },
    silent = TRUE
  )



  if (output_to_file && !is.null(output_result) && is.data.frame(result) && nrow(result) > 0) {
    readr::write_csv(result, output_result, append = file.exists(output_result))
  }

  return(list(result = result, manifest = manifest_info))
}

compute_logliks <- function(
    po, parameters, pfilter_control = NULL,
    allow_parallel = FALSE, n_cores = 1,
    seed = TRUE, chunk.size = NULL,
    scheduling = 1, verbose = FALSE) {
  # preserve and restore the caller's future plan even on error/exit
  old_plan <- tryCatch(
    {
      future::plan()
    },
    error = function(e) NULL
  )
  on.exit(
    {
      if (!is.null(old_plan)) {
        tryCatch(
          {
            future::plan(old_plan)
          },
          error = function(e) NULL
        )
      }
    },
    add = TRUE
  )

  if (is.null(pfilter_control)) pfilter_control <- list(Nrep = 10, Np = 10000)

  # choose execution plan
  if (allow_parallel) {
    ## if parallel register cluster
    plan(multicore, workers = n_cores)
  } else {
    plan(sequential)
  }
  if (is.null(seed)) {
    seed_num <- TRUE
  } else {
    seed_num <- seed
  }

  npo <- nrow(parameters)
  njobs <- pfilter_control$Nrep * npo

  res <- foreach(
    iter_i = seq_len(njobs),
    .combine = c,
    .options.future = list(
      seed = seed,
      chunk.size = chunk.size,
      scheduling = scheduling
    )
  ) %dofuture% {
    ipo <- (iter_i - 1) %% npo + 1
    rep <- (iter_i - 1) %/% npo + 1
    # assemble parameter vector for this job (ipo)
    param_row <- NULL
    if (is.data.frame(parameters) || tibble::is_tibble(parameters)) {
      param_row <- as.numeric(parameters[ipo, , drop = TRUE])
      names(param_row) <- colnames(parameters)
    } else if (is.list(parameters) && length(parameters) >= ipo) {
      param_row <- unlist(parameters[[ipo]])
    } else if (is.numeric(parameters) && !is.null(names(parameters))) {
      param_row <- parameters
    } else {
      stop("Unsupported 'parameters' format in compute_logliks")
    }

    pfilter_args <- list(
      data = po, params = param_row,
      Np = pfilter_control$Np, verbose = verbose
    )


    if (is.null(po)) {
      if (verbose) {
        warning(sprintf("Skipping pfilter for iter %d rep %d: NULL pomp object", ipo, rep))
      }
      pf_res <- NA_real_
    } else {
      pf_res <- tryCatch(
        {
          pfiltered_po <- do.call(pomp::pfilter, pfilter_args)
          logLik(pfiltered_po)
        },
        error = function(e) {
          if (verbose) {
            warning(sprintf("pfilter failed for iter %d rep %d: %s", ipo, rep, conditionMessage(e)))
          }
          NA_real_
        }
      )
    }
    pf_res
  }

  loglik_matrix <- matrix(res, ncol = pfilter_control$Nrep, byrow = FALSE)
  # compute per-outer-loglik summaries (mean and se). Handle single-row result from apply()
  loglik_raw <- apply(loglik_matrix, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) > 0) {
      pomp::logmeanexp(x, se = TRUE)
    } else {
      c(NA_real_, NA_real_)
    }
  })

  if (is.null(dim(loglik_raw))) {
    # single outer iteration -> loglik_raw is a named numeric vector (mean, se)
    nm <- names(loglik_raw)
    if (is.null(nm)) nm <- c("mean", "se")
    loglik <- matrix(loglik_raw, nrow = length(loglik_raw), ncol = 1)
    rownames(loglik) <- nm
  } else {
    loglik <- loglik_raw
    # ensure rownames exist (some older R versions may drop names)
    if (is.null(rownames(loglik))) {
      rn <- names(loglik[, 1])
      if (!is.null(rn)) rownames(loglik) <- rn
    }
  }
  result <- parameters
  result$loglik <- as.numeric(loglik[1, ])
  result$loglik.se <- as.numeric(loglik[2, ])
  return(result)
}


#' Simulate from a prepared POMP model
#'
#' Convenience wrapper to build a `pomp` object from covariate and case
#' datasets (via `preprocess_data()` and `prepare_pomp_model()`), then run
#' stochastic simulations using `pomp::simulate()`.
#'
#' This function accepts the same `model_dir` / `model_file` options as
#' `fit_model()` and supports a `validate_only` mode to just validate the
#' model file and data inputs without running costly simulations.
#'
#' @param parameters Named numeric vector or single-row data.frame/tibble with parameter values used for simulation.
#' @param city Character. City name (used only for consistency with `fit_model`).
#' @param n_sims Integer. Number of simulation trajectories to generate.
#' @param model Character. Short model name used to locate `pomp_object_<model>.R` in `model_dir`.
#' @param scale Character. Aggregation scale: `'monthly'`, `'daily'`, or `'hourly'`.
#' @param covariate Character. Covariate column name to use in `dataset_VC`.
#' @param start_year,start_month,end_year,end_month Integers selecting time window to preprocess.
#' @param dataset_mal Path to CSV with malaria case data.
#' @param dataset_VC Path to CSV with covariate data (VC/climate).
#' @param output_format Character. One of `'data.frame'` or `'pomps'` returned by `pomp::simulate()`.
#' @param extra_cov Optional character vector of additional covariate names to merge into the covariate dataset.
#' @param extra_cov_dataset Optional path to a CSV containing `extra_cov` columns.
#' @param normalize_extra Logical. Whether to min-max normalize extra covariates.
#' @param model_dir Character. Directory where `pomp_object_<model>.R` files are located. Default: `'dynamical_model/fitting'`.
#' @param model_file Optional character. Absolute path to a `pomp_object_*.R` model file (overrides `model` + `model_dir`).
#' @param validate_only Logical. If `TRUE` perform validations and return early with `list(status = 'ok')`.
#'
#' @return Simulation results in the format requested by `output_format` (matching `pomp::simulate()` behavior) or validation list when `validate_only = TRUE`.
#'
#' @export
simulate_model <- function(
    parameters, city, n_sims = 1000,
    model = "const", scale = "monthly", covariate = "VC",
    start_year = 1997, start_month = 1,
    end_year = 2014, end_month = 12,
    dataset_mal, dataset_VC, output_format = c("data.frame", "pomps"),
    extra_cov = NULL, extra_cov_dataset = NULL,
    normalize_covariate = TRUE, normalize_extra = TRUE,
    model_dir = "dynamical_model/fitting", model_file = NULL, validate_only = FALSE) {
  output_format <- match.arg(output_format)

  data_mal <- read_csv(dataset_mal, show_col_types = FALSE)
  data_cov <- read_csv(dataset_VC, show_col_types = FALSE)
  if (!is.null(extra_cov_dataset)) {
    data_extra_cov <- read_csv(extra_cov_dataset, show_col_types = FALSE)
  } else {
    data_extra_cov <- NULL
  }


  all_dat <- preprocess_data(
    data_cov, data_mal, scale, covariate, start_year,
    end_year, start_month, end_month, extra_cov, data_extra_cov,
    normalize_covariate, normalize_extra
  )

  dat <- all_dat[[1]]
  dat_mal_raw <- all_dat[[2]]
  # robust time column selection (accept time or time_end)

  if ("time_end" %in% colnames(dat_mal_raw)) {
    dat_mal <- dat_mal_raw |>
      select(time_end, PF) |>
      rename(time = time_end)
  } else if ("time" %in% colnames(dat_mal_raw)) {
    dat_mal <- dat_mal_raw |> select(time, PF)
  } else {
    stop("`dat_mal` produced by preprocess_data must contain 'time' or 'time_end' and 'PF'.")
  }

  ### prepare pomp model
  if (tibble::is_tibble(parameters) | is.data.frame(parameters)) {
    param <- as.numeric(parameters[1, ])
    param_names <- colnames(parameters)
    names(param) <- param_names
  } else { ## expects named vector
    param <- parameters
  }

  po <- prepare_pomp_model(model, param, dat, dat_mal, model_dir = model_dir, model_file = model_file, validate_only = validate_only)
  if (validate_only) {
    return(po)
  }

  # Validate that supplied 'parameters' contain the model's expected parameter names
  model_env_attr <- tryCatch(attr(po, "model_env"), error = function(e) NULL)
  if (!is.null(model_env_attr)) {
    model_paramnames <- c(get("par_names", envir = model_env_attr), get("vp_names", envir = model_env_attr))
  } else {
    model_paramnames <- NULL
  }

  supplied_param_names <- NULL
  if (tibble::is_tibble(parameters) || is.data.frame(parameters)) {
    supplied_param_names <- colnames(parameters)
  } else if (!is.null(names(parameters))) {
    supplied_param_names <- names(parameters)
  }

  if (!is.null(model_paramnames) && !is.null(supplied_param_names)) {
    missing_params <- setdiff(model_paramnames, supplied_param_names)
    extra_params <- setdiff(supplied_param_names, model_paramnames)
    if (length(missing_params) > 0) {
      stop(
        "Parameter mismatch: supplied 'parameters' are missing model parameters: ", paste(missing_params, collapse = ", "),
        ".\nExpected parameter names (", length(model_paramnames), "): ", paste(model_paramnames, collapse = ", "),
        ".\nSupplied names (", length(supplied_param_names), "): ", paste(supplied_param_names, collapse = ", ")
      )
    }
    if (length(extra_params) > 0) {
      warning("Supplied 'parameters' contain extra names not present in the model: ", paste(extra_params, collapse = ", "))
    }
  }

  result <- po |> simulate(nsim = n_sims, include.data = TRUE, format = output_format)

  return(result)
}


get_top_sample <- function(dataset, k = 1) {
  dataset2 <- dataset |>
    arrange(desc(loglik)) |>
    select(-one_of("sample", "loglik", "loglik.se", "flag"))
  res <- dataset2[k, ] |>
    unlist() |>
    as.numeric()
  names(res) <- colnames(dataset2)

  return(res)
}


plot_simulation <- function(dataset, mode = c("ribbon", "traj"), CI = 0.80, font_size = 12) {
  mode <- match.arg(mode)

  dataset <- dataset |> mutate(type = if_else(.id == "data", "Data", "Simulation"))
  res <- dataset |> filter(type == "Data")
  head(res)
  max_res <- max(res$PF)

  colors <- c("#000000", "#7570b3")
  colors2 <- c("#FFFFFF", "#D5CCE6")
  if (mode == "ribbon") {
    p <- dataset |> ggplot(aes(
      x = time, y = PF, fill = type, color = type
    )) +
      ggdist::stat_lineribbon(.width = CI, linewidth = 1) +
      theme_bw(font_size) +
      theme(legend.position = "bottom", axis.title.y = element_markdown()) +
      labs(
        x = "Time", y = "Observed *P. falciparum* cases", color = "Type",
        fill = "Type"
      ) +
      theme(legend.position = "bottom") +
      scale_fill_manual(values = colors2) +
      scale_color_manual(values = colors) +
      coord_cartesian(ylim = c(0, 1.5 * max_res))
  } else if (mode == "traj") {
    p <- dataset |> ggplot(aes(
      x = time, y = PF, color = type, group = .id
    )) +
      geom_line(alpha = 0.5) +
      theme_bw(font_size) +
      theme(legend.position = "bottom", axis.title.y = element_markdown()) +
      labs(
        x = "Time", y = "Observed *P. falciparum* cases", color = "Type",
        fill = "Type"
      ) +
      theme(legend.position = "bottom") +
      scale_color_manual(values = colors) +
      coord_cartesian(ylim = c(0, 2 * max_res))
  }
  return(p)
}
