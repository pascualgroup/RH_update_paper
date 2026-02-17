library(tidyverse)
library(pomp)
library(tsibble)
library(feasts)
library(here)

source(here("entomological_model/fitting", "fit_model.r"))

construct_datetime_column <- function(df,
                                      tz = "UTC",
                                      out_col = "date",
                                      existing_date_cols = c("date", "Date", "datetime", "Datetime")) {
  # Skip if an existing date/datetime column is present
  if (any(names(df) %in% existing_date_cols)) {
    candidates <- c("date", "Date", "datetime", "Datetime", "DATE")
    if (!"date" %in% names(df)) {
      src <- intersect(candidates, names(df))
      if (length(src)) names(df)[names(df) == src[1]] <- "date"
    }
    return(df)
  }

  # Case-insensitive matching for components
  nm <- names(df)
  lower <- tolower(nm)
  has <- function(cols) all(cols %in% lower)
  get <- function(name) df[[nm[match(name, lower)]]]

  # Require lubridate for safe constructors
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    stop("Please install.packages('lubridate') to construct datetime.")
  }

  if (all(c("year", "month", "day", "hour", "minute", "second"))) {
    dt <- df |> mutate(lubridate::make_datetime(
      year = get("year"), month = get("month"), day = get("day"),
      hour = get("hour"), min = get("minute"), sec = get("second"), tz = tz
    ))
  } else if (has(c("year", "month", "day", "hour", "minute"))) {
    dt <- lubridate::make_datetime(
      year = get("year"), month = get("month"), day = get("day"),
      hour = get("hour"), min = get("minute"), sec = 0, tz = tz
    )
  } else if (has(c("year", "month", "day", "hour"))) {
    dt <- lubridate::make_datetime(
      year = get("year"), month = get("month"), day = get("day"),
      hour = get("hour"), min = 0, sec = 0, tz = tz
    )
  } else if (has(c("year", "month", "day"))) {
    dt <- as.POSIXct(lubridate::make_date(
      year = get("year"), month = get("month"), day = get("day")
    ), tz = tz)
  } else if (has(c("year", "month"))) {
    dt <- as.POSIXct(lubridate::make_date(
      year = get("year"), month = get("month"), day = 1
    ), tz = tz)
  } else if (has(c("year"))) {
    dt <- as.POSIXct(lubridate::make_date(
      year = get("year"), month = 1, day = 1
    ), tz = tz)
  } else {
    stop("No suitable columns to construct a datetime: need at least 'year'.")
  }

  df[[out_col]] <- dt
  df
}

# Detect scale and pick MSTL period
detect_scale_period <- function(time_vec) {
  t <- sort(as.POSIXct(time_vec, tz = "UTC"))
  dt <- diff(as.numeric(t)) # seconds
  step <- median(dt, na.rm = TRUE)

  hr <- 3600
  day <- 24 * hr

  # simple tolerances
  is_hourly <- abs(step - hr) <= 0.25 * hr
  is_daily <- abs(step - day) <= 0.25 * day
  # monthly ≈ 28–31 days; detect by month jumps
  is_monthly <- !is_hourly && !is_daily &&
    median(diff(year(t) * 12 + month(t)), na.rm = TRUE) == 1

  if (is_hourly) {
    return(list(scale = "hourly", period = c(24, 24 * 365), drop_feb29 = TRUE))
  } # daily seasonality
  else if (is_daily) {
    return(list(scale = "daily", period = 365, drop_feb29 = TRUE))
  } # annual seasonality
  else if (is_monthly) {
    return(list(scale = "monthly", period = 12, drop_feb29 = FALSE))
  } # annual seasonality

  stop("Unsupported sampling scale. Expect hourly, daily, or monthly.")
}


decompose_with_forecast_mstl <- function(df, column, opts) {
  scale <- opts$scale
  if (scale == "hourly") {
    if (opts$drop_feb29) {
      # remove Feb 29 data points
      df <- df |>
        filter(!(month(date) == 2 & day(date) == 29))
    }
    tsx <- forecast::msts(
      df[[column]],
      seasonal.periods = opts$period
    )
    mstl_v <- forecast::mstl(tsx,
      lambda = NULL,
      s.window = "periodic", robust = TRUE
    ) |>
      as_tibble() |>
      mutate(date = df$date, Seasonal = Seasonal24 + Seasonal8760 + mean(Trend, na.rm = TRUE)) |>
      select(date, Data, Seasonal, Trend, Remainder)
  } else if (scale == "daily") {
    if (opts$drop_feb29) {
      # remove Feb 29 data points
      df <- df |>
        filter(!(month(date) == 2 & day(date) == 29))
    }
    tsx <- forecast::msts(
      df[[column]],
      seasonal.periods = opts$period
    )
    mstl_v <- forecast::mstl(tsx,
      lambda = NULL,
      s.window = "periodic", robust = TRUE
    ) |>
      as_tibble() |>
      mutate(date = df$date, Seasonal = Seasonal365 + mean(Trend, na.rm = TRUE)) |>
      select(date, Data, Seasonal, Trend, Remainder)
  } else if (scale == "monthly") {
    if (opts$drop_feb29) {
      # remove Feb 29 data points
      df <- df |>
        filter(!(month(date) == 2 & day(date) == 29))
    }
    tsx <- forecast::msts(
      df[[column]],
      seasonal.periods = opts$period
    )
    mstl_v <- forecast::mstl(tsx,
      lambda = NULL,
      s.window = "periodic", robust = TRUE
    ) |>
      as_tibble() |>
      mutate(date = df$date, Seasonal = Seasonal12 + mean(Trend, na.rm = TRUE)) |>
      select(date, Data, Seasonal, Trend, Remainder)
  } else {
    stop("Unsupported scale: ", scale)
  }
  return(mstl_v)
}


# sum all Seasonal* columns from forecast::mstl()
.total_seasonal <- function(mstl_df) {
  ss <- grep("^Seasonal", names(mstl_df))
  if (length(ss) == 0) stop("No Seasonal* columns in mstl output.")
  as.numeric(rowSums(mstl_df[, ss, drop = FALSE]))
}

# Decompose a *single* numeric vector with MSTL, return tibble with date, key, season_total, trend, remainder
.mstl_one <- function(x, dates, periods, drop_feb29 = FALSE,
                      s.window = "periodic", t.window = NULL, robust = TRUE) {
  df <- tibble::tibble(date = as.POSIXct(dates, tz = "UTC"), value = x)
  if (drop_feb29) df <- dplyr::filter(df, !(lubridate::month(date) == 2 & lubridate::mday(date) == 29))
  tsx <- forecast::msts(df$value, seasonal.periods = periods)
  fit <- forecast::mstl(tsx, s.window = s.window, t.window = t.window, robust = robust)
  out <- as_tibble(fit)
  tibble::tibble(
    date       = df$date,
    season_tot = .total_seasonal(out),
    trend      = out$Trend,
    remainder  = out$Remainder
  )
}

# Decompose *multiple* variables with MSTL (vars: character vector), then (optionally) aggregate to monthly means
decompose_with_forecast_mstl_multi <- function(df, date_col = "date", vars,
                                               periods, drop_feb29,
                                               monthly_mean_after = TRUE,
                                               s.window = "periodic", t.window = NULL) {
  stopifnot(all(c(date_col, vars) %in% names(df)))
  dates <- df[[date_col]]
  pieces <- lapply(vars, function(v) {
    comp <- .mstl_one(df[[v]], dates, periods, drop_feb29, s.window, t.window)
    comp |>
      dplyr::mutate(key = v) |>
      dplyr::select(date, key, season_tot, trend, remainder)
  })
  long <- dplyr::bind_rows(pieces)

  if (!monthly_mean_after) {
    return(long)
  }

  # aggregate to monthly means (linear → preserves: original = clim + anomaly at monthly scale)
  long |>
    dplyr::mutate(year = lubridate::year(date), month = lubridate::month(date)) |>
    dplyr::group_by(key, year, month) |>
    dplyr::summarize(
      season_tot = mean(season_tot, na.rm = TRUE),
      trend = mean(trend, na.rm = TRUE),
      remainder = mean(remainder, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(date = as.Date(sprintf("%04d-%02d-01", year, month))) |>
    dplyr::select(date, key, season_tot, trend, remainder)
}


adjust_traits_by_temperature <- function(
    Temp,
    trait_data,
    VC_name = "VC") {
  model_params <- tibble(
    parameter = c("a", "pdr", "mdr", "efd", "bc.succ", "e2a", "mu"),
    model = c("briere", "briere", "briere", "quad", "quad", "quad", "quad_mu")
  )

  temp_ts <- trait_data |>
    left_join(model_params, by = join_by(parameter)) |>
    split(~parameter) |>
    lapply(df_to_vector_param_2) |>
    lapply(compute_trait_curve_v2,
      Temp = Temp
    ) |>
    bind_cols() |>
    mutate(!!sym(VC_name) := compute_VC(a, pdr, mdr, efd, bc.succ, e2a, mu)) |>
    select(!!sym(VC_name))

  return(temp_ts)
}


ablate_variables <- function(df, date_col = "date", vars) {
  if (length(vars) == 0) {
    return(list(
      original = tibble::tibble(date = as.Date(NA))[0, ],
      clim     = tibble::tibble(date = as.Date(NA))[0, ],
      flat     = tibble::tibble(date = as.Date(NA))[0, ]
    ))
  }

  # Detect scale from the first var's time index (assume same index for all)
  info <- detect_scale_period(df[[date_col]])

  long <- decompose_with_forecast_mstl_multi(df,
    date_col = date_col, vars = vars,
    periods = info$period, drop_feb29 = info$drop_feb29,
    monthly_mean_after = FALSE,
    s.window = "periodic"
  )

  orig <- long |>
    dplyr::group_by(key) |>
    dplyr::mutate(orig = (if ("season_year" %in% names(long)) season_year else season_tot) +
      trend + remainder) |>
    dplyr::ungroup() |>
    dplyr::select(date, key, orig) |>
    tidyr::pivot_wider(names_from = key, values_from = orig)

  # Turn long into wide products
  clim <- long |>
    dplyr::group_by(key) |>
    dplyr::mutate(clim = (if ("season_year" %in% names(long)) season_year else season_tot) +
      mean(trend, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::select(date, key, clim) |>
    tidyr::pivot_wider(names_from = key, values_from = clim)



  flat <- long |>
    dplyr::group_by(key) |>
    dplyr::mutate(val = mean(trend, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::select(date, key, val) |>
    tidyr::pivot_wider(names_from = key, values_from = val)



  # # Monthly original = monthly mean of raw (linear)
  #   dplyr::mutate(ym = as.Date(sprintf(
  #     "%04d-%02d-01", lubridate::year(.data[[date_col]]),
  #     lubridate::month(.data[[date_col]])
  #   ))) |>
  #   dplyr::group_by(ym) |>
  #   dplyr::summarize(dplyr::across(dplyr::all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  #   dplyr::rename(date = ym)

  list(original = orig, clim = clim, flat = flat)
}


prepare_data_for_covariate_importance <- function(
    dataset_mal, dataset_VC,
    scale = "monthly", covariate = "VC",
    temperature_covariate = NULL,
    start_year = 1997, start_month = 1,
    end_year = 2011, end_month = 12,
    extra_cov = NULL, extra_cov_dataset = NULL,
    normalize_covariate = TRUE,
    normalize_extra = TRUE,
    trait_file = here("entomological_model", "data", "point_estimate_params.csv")) {
  if (scale != "monthly") {
    stop("Currently only 'monthly' scale is supported")
  }
  # accept file paths or data.frames/tibbles for inputs
  if (is.character(dataset_mal) && length(dataset_mal) == 1) {
    data_mal <- readr::read_csv(dataset_mal, show_col_types = FALSE)
  } else if (is.data.frame(dataset_mal) || tibble::is_tibble(dataset_mal)) {
    data_mal <- dataset_mal
  } else {
    stop("`dataset_mal` must be a file path or a data.frame/tibble")
  }

  if (is.character(dataset_VC) && length(dataset_VC) == 1) {
    data_cov <- readr::read_csv(dataset_VC, show_col_types = FALSE)
  } else if (is.data.frame(dataset_VC) || tibble::is_tibble(dataset_VC)) {
    data_cov <- dataset_VC
  } else {
    stop("`dataset_VC` must be a file path or a data.frame/tibble")
  }

  if (!is.null(extra_cov_dataset)) {
    if (is.character(extra_cov_dataset) && length(extra_cov_dataset) == 1) {
      data_extra_cov <- readr::read_csv(extra_cov_dataset, show_col_types = FALSE)
    } else if (is.data.frame(extra_cov_dataset) || tibble::is_tibble(extra_cov_dataset)) {
      data_extra_cov <- extra_cov_dataset
    } else {
      stop("`extra_cov_dataset` must be a file path or a data.frame/tibble when provided")
    }
  } else {
    data_extra_cov <- NULL
  }

  if (!is.null(temperature_covariate)) {
    if (!(temperature_covariate %in% names(data_cov))) {
      stop(paste0("Temperature covariate '", temperature_covariate, "' not found in dataset_VC"))
    }
    if (is.null(trait_file)) {
      stop("`trait_file` must be provided when `temperature_covariate` is specified")
    } else {
      trait_data <- readr::read_csv(trait_file, show_col_types = FALSE)
      ### recompute trait values based on temperature covariate
      data_cov <- construct_datetime_column(data_cov, out_col = "date")
      data_cov <- data_cov |>
        dplyr::filter(
          lubridate::year(date) >= start_year,
          lubridate::year(date) <= end_year
        ) |>
        dplyr::filter(!(lubridate::year(date) == start_year & lubridate::month(date) < start_month)) |>
        dplyr::filter(!(lubridate::year(date) == end_year & lubridate::month(date) > end_month))

      VC_df <- ablate_variables(
        df = data_cov,
        date_col = "date",
        vars = c(temperature_covariate)
      ) |> lapply(function(temp_df) {
        adjust_traits_by_temperature(
          Temp = temp_df[[temperature_covariate]],
          trait_data = trait_data,
          VC_name = covariate
        ) |>
          mutate(date = temp_df$date) |>
          select(date, !!sym(covariate))
      })
      data_extra_cov <- construct_datetime_column(data_extra_cov, out_col = "date")
      extra_vals <- ablate_variables(
        df = data_extra_cov, date_col = "date",
        vars = c(extra_cov)
      )

      joined_dfs <- list(
        original = list(cov = VC_df$original, extra = extra_vals$original),
        clim = list(cov = VC_df$clim, extra = extra_vals$clim),
        flat = list(cov = VC_df$flat, extra = extra_vals$flat)
      )

      processed_data <- joined_dfs |> lapply(function(x) {
        res <- preprocess_data(
          x$cov, data_mal, scale, covariate, start_year,
          end_year, start_month, end_month, extra_cov,
          x$extra,
          normalize_covariate = FALSE,
          normalize_extra = FALSE
        )$dat
        new_names <- setNames(
          c(
            "covariate",
            paste0("cov", seq_along(extra_cov))
          ),
          c(covariate, extra_cov)
        )

        res <- res |>
          rename(!!!new_names)
        return(res)
      })

      malaria_data <- preprocess_data(
        joined_dfs$original$cov, data_mal,
        scale, covariate, start_year,
        end_year, start_month, end_month, extra_cov,
        joined_dfs$original$extra,
        normalize_covariate = FALSE,
        normalize_extra = FALSE
      )$dat_mal

      if ("time_end" %in% colnames(malaria_data)) {
        dat_mal <- malaria_data |>
          dplyr::select(time_end, PF) |>
          dplyr::rename(time = time_end)
      } else if ("time" %in% colnames(malaria_data)) {
        dat_mal <- malaria_data |>
          dplyr::select(time, PF)
      } else {
        stop("`dat_mal` produced by preprocess_data must contain 'time' or 'time_end' and 'PF'.")
      }

      original <- processed_data$original
      clim <- processed_data$clim
      flat <- processed_data$flat

      scalers <- NULL
      if (normalize_covariate || normalize_extra) {
        cols_to_scale <- c(
          if (normalize_covariate) covariate else NULL,
          if (normalize_extra) extra_cov else NULL
        )
        if (length(cols_to_scale)) {
          scalers <- caret::preProcess(original[, cols_to_scale, drop = FALSE], method = "range")
          original[, cols_to_scale] <- predict(scalers, original[, cols_to_scale, drop = FALSE])
          clim[, cols_to_scale] <- predict(scalers, clim[, cols_to_scale, drop = FALSE])
          flat[, cols_to_scale] <- predict(scalers, flat[, cols_to_scale, drop = FALSE])
        }
      }

      return(list(
        original = original,
        clim = clim,
        flat = flat,
        scalers = scalers,
        malaria_data = dat_mal
      ))
    }
  } else {
    data_cov <- construct_datetime_column(data_cov, out_col = "date")
    data_cov <- data_cov |>
      dplyr::filter(
        lubridate::year(date) >= start_year,
        lubridate::year(date) <= end_year
      ) |>
      dplyr::filter(!(lubridate::year(date) == start_year & lubridate::month(date) < start_month)) |>
      dplyr::filter(!(lubridate::year(date) == end_year & lubridate::month(date) > end_month))

    VC_df <- ablate_variables(
      df = data_cov,
      date_col = "date",
      vars = c(covariate)
    )
    data_extra_cov <- construct_datetime_column(data_extra_cov, out_col = "date")
    extra_vals <- ablate_variables(
      df = data_extra_cov, date_col = "date",
      vars = c(extra_cov)
    )

    joined_dfs <- list(
      original = list(cov = VC_df$original, extra = extra_vals$original),
      clim = list(cov = VC_df$clim, extra = extra_vals$clim),
      flat = list(cov = VC_df$flat, extra = extra_vals$flat)
    )

    processed_data <- joined_dfs |> lapply(function(x) {
      res <- preprocess_data(
        x$cov, data_mal, scale, covariate, start_year,
        end_year, start_month, end_month, extra_cov,
        x$extra,
        normalize_covariate = FALSE,
        normalize_extra = FALSE
      )$dat
      new_names <- setNames(
        c(
          "covariate",
          paste0("cov", seq_along(extra_cov))
        ),
        c(covariate, extra_cov)
      )

      res <- res |>
        rename(!!!new_names)
      return(res)
    })

    malaria_data <- preprocess_data(
      joined_dfs$original$cov, data_mal,
      scale, covariate, start_year,
      end_year, start_month, end_month, extra_cov,
      joined_dfs$original$extra,
      normalize_covariate = FALSE,
      normalize_extra = FALSE
    )$dat_mal

    if ("time_end" %in% colnames(malaria_data)) {
      dat_mal <- malaria_data |>
        dplyr::select(time_end, PF) |>
        dplyr::rename(time = time_end)
    } else if ("time" %in% colnames(malaria_data)) {
      dat_mal <- malaria_data |>
        dplyr::select(time, PF)
    } else {
      stop("`dat_mal` produced by preprocess_data must contain 'time' or 'time_end' and 'PF'.")
    }

    original <- processed_data$original
    clim <- processed_data$clim
    flat <- processed_data$flat

    scalers <- NULL
    if (normalize_covariate || normalize_extra) {
      cols_to_scale <- c(
        if (normalize_covariate) covariate else NULL,
        if (normalize_extra) extra_cov else NULL
      )
      if (length(cols_to_scale)) {
        scalers <- caret::preProcess(original[, cols_to_scale, drop = FALSE], method = "range")
        original[, cols_to_scale] <- predict(scalers, original[, cols_to_scale, drop = FALSE])
        clim[, cols_to_scale] <- predict(scalers, clim[, cols_to_scale, drop = FALSE])
        flat[, cols_to_scale] <- predict(scalers, flat[, cols_to_scale, drop = FALSE])
      }
    }

    return(list(
      original = original,
      clim = clim,
      flat = flat,
      scalers = scalers,
      malaria_data = dat_mal
    ))
  }
}

# ============================================================
# Build all policy covariate combinations for ablation tests
# ============================================================
build_policy_covariate_combinations <- function(prep_data,
                                                covariates,
                                                policy = c("flat", "clim"),
                                                max_order = length(covariates)) {
  policy <- match.arg(policy)

  # Extract preprocessed tables
  cov_base <- prep_data$original
  cov_flat <- prep_data$flat
  cov_clim <- prep_data$clim

  # Identify numeric columns to combine (covariates only)
  cov_numeric <- intersect(covariates, names(cov_base))

  # Generate all non-empty subsets up to max_order
  combo_list <- unlist(lapply(1:max_order, function(k) {
    combn(cov_numeric, k, simplify = FALSE)
  }), recursive = FALSE)

  # Build function to substitute one combination
  make_policy_df <- function(S) {
    out <- cov_base
    for (v in S) {
      policy_df <- if (policy == "flat") cov_flat else cov_clim
      out[[v]] <- policy_df[[v]]
    }
    out
  }

  # Apply over all combinations
  policy_dfs <- lapply(combo_list, make_policy_df)

  # Add metadata for clarity
  names(policy_dfs) <- sapply(combo_list, function(S) paste(S, collapse = "+"))

  tibble::tibble(
    combination = names(policy_dfs),
    covariates = combo_list,
    policy = policy,
    data = policy_dfs
  )
}

compute_logliks <- function(
    model, parameters, dat, dat_mal,
    model_dir, model_file,
    cov_combo, covariate,
    extra_cov,
    policy,
    combination,
    Np, Nfilter_reps,
    start_year, end_year,
    start_month, end_month,
    allow_parallel = FALSE,
    n_cores = NULL,
    verbose = FALSE,
    seed = NULL) {
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
  new_names <- setNames(
    c(covariate, extra_cov),
    c(
      "covariate",
      paste0("cov", seq_along(extra_cov))
    )
  )
  # choose execution plan
  if (allow_parallel) {
    ## if parallel register cluster
    plan(multicore, workers = n_cores)
  } else {
    plan(sequential)
  }
  dat <- dat |>
    rename(!!!new_names)

  po <- prepare_pomp_model(
    model = model, param = parameters, dat_cov = dat, dat_mal = dat_mal,
    model_dir = model_dir,
    model_file = model_file
  )

  # po <- prepare_model(parameters,
  #   dataset_mal = dat_mal, dataset_VC = dat,
  #   extra_cov_dataset = dat, extra_cov = extra_cov,
  #   scale = "monthly", covariate = covariate,
  #   start_year, end_year, start_month, end_month,
  #   normalize_covariate = FALSE,
  #   normalize_extra = FALSE,
  #   model = model, model_dir = model_dir,
  #   model_file = model_file
  # )$po

  pfilter_args <- c(
    list(data = po),
    Np = Np, Nrep = Nfilter_reps, list(seed = if (!is.null(seed)) seed else TRUE, verbose = verbose)
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

  cond_logliks <- lapply(pfout, function(x) {
    x@cond.logLik |> as_tibble_col()
  }) |>
    bind_cols() |>
    rowMeans(na.rm = TRUE)


  return(cond_logliks)
}

compute_scores <- function(all_results) {
  library(dplyr)
  library(purrr)

  # Identify baseline conditional log-likelihood vector
  base <- all_results |>
    purrr::keep(~ .x$policy == "baseline" || length(.x$combination) == 0) |>
    purrr::pluck(1) # first baseline entry
  base_ll <- base$cond.loglik

  # Combine results into a tibble for easy manipulation
  res_df <- purrr::map_dfr(all_results, function(r) {
    tibble(
      policy = r$policy,
      combination = paste(r$combination, collapse = "+"),
      cond.loglik = list(r$cond.loglik)
    )
  })

  # Compute total and delta scores
  res_df <- res_df %>%
    mutate(
      total_loglik = map_dbl(cond.loglik, sum, na.rm = TRUE),
      delta_loglik = map_dbl(cond.loglik, ~ sum(.x - base_ll, na.rm = TRUE)),
      mean_delta = map_dbl(cond.loglik, ~ mean(.x - base_ll, na.rm = TRUE))
    )

  res_df
}



estimate_covariate_importance <- function(
    parameters,
    Np = 1000, Nfilter_reps = 10,
    NBlocks = 1000, Block_size = 12,
    model = "const", scale = "monthly",
    covariate = "VC",
    temperature_covariate = NULL,
    start_year = 1997, start_month = 1,
    end_year = 2014, end_month = 12,
    dataset_mal, dataset_VC,
    extra_cov = NULL, extra_cov_dataset = NULL,
    normalize_covariate = TRUE,
    normalize_extra = TRUE,
    model_dir = "dynamical_model/fitting",
    model_file = NULL,
    trait_file = NULL,
    max_order = min(2, length(c(covariate, extra_cov))),
    allow_parallel = FALSE,
    n_cores = NULL) {
  # prepare data
  prep_data <- prepare_data_for_covariate_importance(
    dataset_mal = dataset_mal, dataset_VC = dataset_VC,
    scale = scale, covariate = covariate,
    temperature_covariate = temperature_covariate,
    start_year = start_year, start_month = start_month,
    end_year = end_year, end_month = end_month,
    extra_cov = extra_cov, extra_cov_dataset = extra_cov_dataset,
    normalize_covariate = normalize_covariate,
    normalize_extra = normalize_extra,
    trait_file = trait_file
  )

  # prepare policy covariate combinations
  policy_combos <- rbind(
    tibble(
      combination = "BASE",
      covariates  = list(character(0)),
      policy      = "baseline",
      data        = list(prep_data$original)
    ),
    build_policy_covariate_combinations(
      prep_data = prep_data,
      covariates = c(covariate, extra_cov),
      policy = "flat",
      max_order = max_order
    ),
    build_policy_covariate_combinations(
      prep_data = prep_data,
      covariates = c(covariate, extra_cov),
      policy = "clim",
      max_order = max_order
    )
  )
  dat_mal <- prep_data$malaria_data
  all_results <- list(nrow(policy_combos))

  seed_num <- as.integer(123)

  for (i in seq_len(nrow(policy_combos))) {
    cov_combo <- policy_combos$covariates[[i]]
    policy <- policy_combos$policy[i]
    combination <- policy_combos$combination[i]
    dat <- policy_combos$data[[i]]
    # |>
    #   mutate(time = (year + (month) / 12)) |>
    #   select(time, everything(), -year, -month)
    result <- compute_logliks(model, parameters, dat, dat_mal,
      model_dir = model_dir, model_file = model_file,
      cov_combo = cov_combo, covariate = covariate,
      extra_cov = extra_cov,
      policy = policy,
      combination = combination,
      Np = Np, Nfilter_reps = Nfilter_reps,
      start_year, end_year,
      start_month, end_month,
      seed = seed_num
    )
    all_results[[i]] <- list(
      policy = policy,
      combination = combination,
      cond.loglik = result,
      cov_combo = cov_combo
    )
  }
  # return(all_results)
  # all_results <- unlist(all_results, recursive = FALSE)

  result <- compute_scores(all_results)
  return(result)
}

canon_split <- function(key_str) {
  if (is.na(key_str) || key_str == "" || key_str == "BASE" || key_str == "<BASE>") {
    return(character(0))
  }
  parts <- unlist(str_split(key_str, "\\+"))
  parts <- str_trim(parts)
  parts <- parts[nzchar(parts)]
  sort(unique(parts))
}
canon_join <- function(parts) if (length(parts) == 0) "<BASE>" else paste(sort(parts), collapse = "+")
canon_key <- function(key_str) canon_join(canon_split(key_str)) # combo string -> canonical key
key_of <- function(S) canon_join(S) # subset vector -> canonical key

# Safe combn for k=0..n
subsets_k <- function(items, k) if (k == 0) list(character(0)) else as.list(combn(items, k, simplify = FALSE))
w_shapley <- function(k, p) factorial(k) * factorial(p - k - 1) / factorial(p)


shapley_from_policy <- function(df_policy, base_val = 0) {
  dfp <- df_policy %>%
    mutate(key = vapply(combination, canon_key, character(1))) %>%
    group_by(key) %>%
    summarize(total_loglik = dplyr::first(total_loglik), .groups = "drop")

  vars <- dfp$key[dfp$key != "<BASE>"] %>%
    purrr::map(canon_split) %>%
    unlist() %>%
    unique() %>%
    sort()
  p <- length(vars)
  if (p == 0) {
    return(tibble(
      covariate = character(), shapley = numeric(),
      total_gain = numeric(), share = numeric()
    ))
  }
  # Coalition value map (canonical keys only); define empty coalition
  vmap <- dfp %>%
    transmute(key, v = total_loglik) %>%
    deframe()
  if (!("<BASE>" %in% names(vmap))) vmap[["<BASE>"]] <- base_val

  # Compute exact Shapley
  phi <- setNames(numeric(p), vars)
  missing <- character(0)

  for (i in vars) {
    contrib <- 0
    others <- setdiff(vars, i)
    for (k in 0:(p - 1)) {
      for (S in subsets_k(others, k)) {
        S_key <- key_of(S) # uses the SAME canonicalization
        Si_key <- key_of(c(S, i))
        if (!(S_key %in% names(vmap)) || !(Si_key %in% names(vmap))) {
          missing <- c(missing, sprintf("('%s','%s')", S_key, Si_key))
          next
        }
        contrib <- contrib + (vmap[[Si_key]] - vmap[[S_key]]) * w_shapley(k, p)
      }
    }
    phi[i] <- contrib
  }

  if (length(missing)) {
    warning(sprintf(
      "Missing %d coalition pairs (exact Shapley needs all 2^p). Examples: %s",
      length(missing), paste(utils::head(unique(missing), 3), collapse = ", ")
    ))
  }

  total_gain <- sum(phi) # equals v(N)-v(∅); no BASE required
  res <- tibble(
    covariate = vars,
    shapley = phi,
    total_gain = total_gain
  ) |>
    mutate(share = if_else(total_gain == 0, NA_real_, shapley / total_gain))
  return(res)
}

# wrapper over all policies
compute_shapley_all_policies <- function(df) {
  base_val <- df |>
    filter(policy == "baseline") |>
    pull(total_loglik)
  base_val <- base_val[1]

  df %>%
    group_by(policy) |>
    group_modify(~ shapley_from_policy(.x, base_val = base_val)) |>
    ungroup()
}


# shap_tab: columns policy, covariate, shapley  (one row per policy×covariate)
decompose_scales <- function(shap_tab, make_positive_harm = FALSE) {
  wide <- shap_tab %>%
    select(policy, covariate, shapley) %>%
    mutate(phi = shapley) %>%
    select(-shapley) %>%
    pivot_wider(names_from = policy, values_from = phi)

  stopifnot(all(c("clim", "flat") %in% names(wide)))

  out <- wide %>%
    mutate(
      interannual = clim,
      total       = flat,
      seasonal    = flat - clim
    ) %>%
    select(covariate, interannual, seasonal, total)

  if (make_positive_harm) {
    out <- out %>% mutate(
      interannual = -interannual,
      seasonal    = -seasonal,
      total       = -total
    )
  }


  # Optional: shares that sum to 1 within each scale
  shares <- out %>%
    pivot_longer(-covariate, names_to = "scale", values_to = "phi") %>%
    group_by(scale) %>%
    mutate(share = if (sum(phi) == 0) NA_real_ else phi / sum(phi)) %>%
    arrange(scale, covariate) |>
    ungroup()

  return(shares)
}
