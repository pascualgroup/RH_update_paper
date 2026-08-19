library(dplyr)
library(ggplot2)

mcap <- function(lp, parameter, confidence = 0.95, lambda = 0.75, Ngrid = 1000) {
  smooth_fit <- loess(lp ~ parameter, span = lambda)
  parameter_grid <- seq(min(parameter, na.rm = T), max(parameter, na.rm = T), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit, newdata = parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter - smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda * length(dist))]
  maxdist <- max(dist[included], na.rm = T)
  weight <- rep(0, length(parameter))
  weight[included] <- (1 - (dist[included] / maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b,
    weight = weight,
    data = data.frame(lp = lp, b = parameter, a = -parameter^2)
  )
  b <- unname(coef(quadratic_fit)["b"])
  a <- unname(coef(quadratic_fit)["a"])
  m <- vcov(quadratic_fit)
  var_b <- m["b", "b"]
  var_a <- m["a", "a"]
  cov_ab <- m["a", "b"]
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b / a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1 / (2 * a)
  se_total_squared <- se_mc_squared + se_stat_squared
  delta <- qchisq(confidence, df = 1) * (a * se_mc_squared + 0.5)
  loglik_diff <- max(smoothed_loglik, na.rm = T) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  list(
    lp = lp, parameter = parameter, confidence = confidence,
    quadratic_fit = quadratic_fit, quadratic_max = b / (2 * a),
    smooth_fit = smooth_fit,
    fit = data.frame(
      parameter = parameter_grid,
      smoothed = smoothed_loglik,
      quadratic = predict(quadratic_fit, list(b = parameter_grid, a = -parameter_grid^2))
    ),
    mle = smooth_arg_max, ci = ci, delta = delta,
    se_stat = sqrt(se_stat_squared), se_mc = sqrt(se_mc_squared), se = sqrt(se_total_squared)
  )
}

simple_mcap <- function(lp, parameter, confidence = 0.95, lambda = 0.75, Ngrid = 1000) {
  smooth_fit <- loess(lp ~ parameter, span = lambda, control = loess.control(surface = "direct"))
  parameter_grid <- seq(min(parameter, na.rm = T), max(parameter, na.rm = T), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit, newdata = parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter - smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda * length(dist))]
  maxdist <- max(dist[included], na.rm = T)
  weight <- rep(0, length(parameter))
  weight[included] <- (1 - (dist[included] / maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b,
    weight = weight,
    data = data.frame(lp = lp, b = parameter, a = -parameter^2)
  )
  b <- unname(coef(quadratic_fit)["b"])
  a <- unname(coef(quadratic_fit)["a"])
  m <- vcov(quadratic_fit)
  var_b <- m["b", "b"]
  var_a <- m["a", "a"]
  cov_ab <- m["a", "b"]
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b / a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1 / (2 * a)
  se_total_squared <- se_mc_squared + se_stat_squared
  delta <- qchisq(confidence, df = 1) / 2
  loglik_diff <- max(smoothed_loglik, na.rm = T) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  list(
    lp = lp, parameter = parameter, confidence = confidence,
    quadratic_fit = quadratic_fit, quadratic_max = b / (2 * a),
    smooth_fit = smooth_fit,
    fit = data.frame(
      parameter = parameter_grid,
      smoothed = smoothed_loglik,
      quadratic = predict(quadratic_fit, list(b = parameter_grid, a = -parameter_grid^2))
    ),
    mle = smooth_arg_max, ci = ci, delta = delta,
    se_stat = sqrt(se_stat_squared), se_mc = sqrt(se_mc_squared), se = sqrt(se_total_squared)
  )
}

mcap_plot <- function(mcap_output, font_size = 16) {
  p <- with(mcap_output, {
    data_points <- data.frame(
      parameter = parameter,
      loglik = lp
    )
    if (is.infinite(ci[1])) {
      ci[1] <- min(data_points$parameter)
    }
    if (is.infinite(ci[2])) {
      ci[2] <- max(data_points$parameter)
    }
    p <- fit |> ggplot(aes(x = parameter)) +
      geom_line(aes(y = smoothed), color = "blue", linewidth = 2) +
      # geom_line(aes(y = quadratic), color = "red") +
      geom_point(aes(x = parameter, y = loglik), data = data_points) +
      geom_vline(xintercept = ci, linetype = "dashed") +
      # geom_vline(xintercept = mle) +
      geom_hline(yintercept = max(fit$smoothed, na.rm = T) - delta) +
      geom_point(x = data_points$parameter[which.max(data_points$loglik)], y = data_points$loglik[which.max(data_points$loglik)], color = "green", size = 3) +
      theme_bw(font_size) +
      labs(title = paste(
        "MLE:", round(data_points$parameter[which.max(data_points$loglik)], digits = 4),
        paste0("CI: (", paste0(round(ci[1], digits = 4), ",", round(ci[2], digits = 4), ")"))
      )) +
      ylim(min(c(data_points$loglik, max(fit$smoothed, na.rm = T) - delta - 1)), max(data_points$loglik))
  })

  return(p)
}

bin_by_quantile <- function(dataset, column_name, n_quantiles = 20) {
  cut_borders <- function(x) {
    pattern <- "(\\(|\\[)(-*[0-9]+\\.*[0-9]*),(-*[0-9]+\\.*[0-9]*)(\\)|\\])"

    start <- as.numeric(gsub(pattern, "\\2", x))
    end <- as.numeric(gsub(pattern, "\\3", x))

    data.frame(start, end)
  }
  quant_values <- seq(0, 1, length.out = n_quantiles + 1)
  x.quant <- quantile(dataset |> pull(column_name), probs = quant_values)
  dataset <- dataset |> mutate(quantile = cut(!!sym(column_name), breaks = x.quant))
  dataset <- cbind(dataset, cut_borders(dataset$quantile))
  return(dataset)
}

plot_mcap_binned <- function(dataset, column_name, n_quantiles = 20, max_per_bin = 20, diff_LL = NULL, remove_outlier = TRUE, zscore = TRUE) {
  dataset <- dataset[complete.cases(dataset), ]
  if (!is.null(diff_LL)) {
    dataset <- dataset |> filter(loglik >= max(loglik) - diff_LL)
  }

  dataset <- dataset |>
    bin_by_quantile(column_name, n_quantiles) |>
    group_by(quantile) |>
    slice_max(loglik, n = max_per_bin) |>
    ungroup() |>
    mutate(mid_val = (start + end) / 2) |>
    # group_by(mid_val) |>
    # summarize(loglik = mean(loglik)) |>
    # ungroup() |>
    remove_missing()

  if (remove_outlier) {
    if (zscore) {
      dataset <- dataset |>
        mutate(zscore = (!!sym(column_name) - mean(!!sym(column_name))) / (sd(!!sym(column_name)))) |>
        filter(abs(zscore) <= 2)
    } else {
      col_vals <- dataset |> pull(column_name)
      iqr <- IQR(col_vals)
      q_min <- quantile(col_vals, 0.25) - 1.5 * iqr
      q_max <- quantile(col_vals, 0.75) + 1.5 * iqr
      dataset <- dataset |> filter(!!sym(column_name) >= q_min, !!sym(column_name) <= q_max)
    }
  }
  # dataset <- dataset |>
  #   bin_by_quantile(column_name, n_quantiles) |>
  #   group_by(quantile) |>
  #   slice_max(loglik, n = max_per_bin) |>
  #   ungroup() |>
  #   mutate(mid_val = (start + end) / 2) |>
  #   group_by(mid_val) |>
  #   summarize(loglik = mean(loglik)) |>
  #   ungroup() |>
  #   remove_missing()


  parameter <- dataset |> pull(column_name)
  lp <- dataset |> pull("loglik")
  mcap_res <- mcap(lp, parameter, lambda = 0.75)
  p <- mcap_plot(mcap_res) + labs(x = column_name, y = "loglikelihood")
  return(p)
}

plot_simple_CI <- function(dataset, column_name, n_quantiles = 20, max_per_bin = 20, diff_LL = NULL, remove_outlier = TRUE, zscore = TRUE) {
  dataset <- dataset[complete.cases(dataset), ]

  if (!is.null(diff_LL)) {
    dataset <- dataset |> filter(loglik >= max(loglik) - diff_LL)
  }
  if (remove_outlier) {
    if (zscore) {
      dataset <- dataset |>
        mutate(zscore = (!!sym(column_name) - mean(!!sym(column_name))) / (sd(!!sym(column_name)))) |>
        filter(abs(zscore) <= 2)
    } else {
      col_vals <- dataset |> pull(column_name)
      iqr <- IQR(col_vals)
      q_min <- quantile(col_vals, 0.25) - 1.5 * iqr
      q_max <- quantile(col_vals, 0.75) + 1.5 * iqr
      dataset <- dataset |> filter(!!sym(column_name) >= q_min, !!sym(column_name) <= q_max)
    }
  }

  dataset <- dataset |>
    bin_by_quantile(column_name, n_quantiles) |>
    group_by(quantile) |>
    slice_max(loglik, n = max_per_bin) |>
    ungroup() |>
    mutate(mid_val = (start + end) / 2) |>
    # group_by(mid_val) |>
    # summarize(loglik = mean(loglik)) |>
    # ungroup() |>
    remove_missing()
  parameter <- dataset |> pull(column_name)
  lp <- dataset |> pull("loglik")
  mcap_res <- simple_mcap(lp, parameter, lambda = 0.75)
  mcap_res$delta <- qchisq(0.95, df = 1) / 2
  p <- mcap_plot(mcap_res) + labs(y = "loglikelihood", x = column_name)
  return(p)
}

mcap_checked <- function(
  lp,
  parameter,
  confidence = 0.95,
  lambda = 0.75,
  Ngrid = 1000,
  external_mle = NULL,
  edge_fraction = 0.05,
  min_local_unique = 5
) {
  ok <- is.finite(lp) & is.finite(parameter)
  lp <- lp[ok]
  parameter <- parameter[ok]

  ## ---- LOESS smoother ----

  smooth_fit <- loess(
    lp ~ parameter,
    span = lambda,
    degree = 2,
    control = loess.control(surface = "direct")
  )

  parameter_grid <- seq(
    min(parameter),
    max(parameter),
    length.out = Ngrid
  )

  smoothed_loglik <- predict(
    smooth_fit,
    newdata = data.frame(parameter = parameter_grid)
  )

  valid_grid <- is.finite(smoothed_loglik)

  parameter_grid <- parameter_grid[valid_grid]
  smoothed_loglik <- smoothed_loglik[valid_grid]

  i_max <- which.max(smoothed_loglik)

  smooth_arg_max <- parameter_grid[i_max]
  smooth_ll_max <- smoothed_loglik[i_max]

  ## Explored range actually supported by LOESS predictions
  explored_range <- range(parameter_grid, na.rm = TRUE)

  ## ---- External MLE diagnostic ----
  ##
  ## Calculate this before any possible early return so the
  ## returned object always contains mle_gap.

  mle_gap <- NA_real_

  if (!is.null(external_mle)) {
    ll_external <- predict(
      smooth_fit,
      newdata = data.frame(
        parameter = external_mle
      )
    )

    if (length(ll_external) == 1L && is.finite(ll_external)) {
      mle_gap <- smooth_ll_max - ll_external
    }
  }

  ## ---- Check edge maximum ----

  parameter_range <- range(parameter, na.rm = TRUE)

  edge_width <- edge_fraction * diff(parameter_range)

  peak_at_edge <-
    smooth_arg_max <= parameter_range[1] + edge_width ||
      smooth_arg_max >= parameter_range[2] - edge_width

  ## ---- LOESS-equivalent local weights ----

  dist <- abs(parameter - smooth_arg_max)

  n_local <- max(
    3L,
    ceiling(lambda * length(parameter))
  )

  # Avoid indexing beyond the available observations
  n_local <- min(n_local, length(parameter))

  cutoff_dist <- sort(dist)[n_local]

  included <- dist <= cutoff_dist

  n_unique_local <- length(
    unique(parameter[included])
  )

  maxdist <- max(dist[included], na.rm = TRUE)

  weight <- numeric(length(parameter))

  if (is.finite(maxdist) && maxdist > 0) {
    weight[included] <- (
      1 -
        (dist[included] / maxdist)^3
    )^3
  } else {
    # Degenerate case: all included focal values are identical
    weight[included] <- 1
  }

  ## ---- Local quadratic metamodel ----
  ##
  ## Parameterization:
  ##
  ##   lp = c - a * parameter^2 + b * parameter
  ##
  ## Therefore an upside-down parabola requires a > 0.

  quadratic_fit <- lm(
    lp ~ a + b,
    weights = weight,
    data = data.frame(
      lp = lp,
      b = parameter,
      a = -parameter^2
    )
  )

  a <- unname(coef(quadratic_fit)["a"])
  b <- unname(coef(quadratic_fit)["b"])

  concave <- is.finite(a) && a > 0

  quadratic_stationary <- if (
    is.finite(a) &&
      is.finite(b) &&
      a != 0
  ) {
    b / (2 * a)
  } else {
    NA_real_
  }

  quadratic_max <- if (concave) {
    quadratic_stationary
  } else {
    NA_real_
  }

  quadratic_pred <- predict(
    quadratic_fit,
    newdata = data.frame(
      b = parameter_grid,
      a = -parameter_grid^2
    )
  )

  ## ---- Initial diagnostics ----

  status <- character()

  if (!concave) {
    status <- c(status, "non_concave")
  }

  if (peak_at_edge) {
    status <- c(status, "peak_at_edge")
  }

  if (n_unique_local < min_local_unique) {
    status <- c(status, "few_local_x")
  }

  ## ============================================================
  ## Non-concave case
  ## ============================================================
  ##
  ## The MCAP variance and cutoff formulas are not valid when
  ## a <= 0. However, return the explored parameter range so that
  ## downstream plotting/table code still has finite limits.
  ##
  ## IMPORTANT:
  ## ci here is NOT a 95% MCAP CI. It is the explored profile range.
  ## ============================================================

  if (!concave) {
    return(list(
      lp = lp,
      parameter = parameter,
      confidence = confidence,
      quadratic_fit = quadratic_fit,
      quadratic_stationary = quadratic_stationary,
      quadratic_max = quadratic_max,
      smooth_fit = smooth_fit,
      fit = data.frame(
        parameter = parameter_grid,
        smoothed = smoothed_loglik,
        quadratic = quadratic_pred
      ),
      mle = smooth_arg_max,

      ## Return finite explored limits rather than NA
      ci = explored_range,

      ## Both sides are range-limited, not inferential MCAP bounds
      ci_truncated = c(
        lower = TRUE,
        upper = TRUE
      ),

      ## Explicitly indicate this is not a valid MCAP CI
      ci_valid = FALSE,
      ci_type = "explored_range",
      explored_range = explored_range,
      delta = NA_real_,
      se_stat = NA_real_,
      se_mc = NA_real_,
      se = NA_real_,
      a = a,
      b = b,
      concave = FALSE,
      peak_at_edge = peak_at_edge,
      n_unique_local = n_unique_local,
      mle_gap = mle_gap,
      status = unique(status)
    ))
  }

  ## ============================================================
  ## MC uncertainty
  ## ============================================================

  m <- vcov(quadratic_fit)

  var_b <- m["b", "b"]
  var_a <- m["a", "a"]
  cov_ab <- m["a", "b"]

  se_mc_squared <-
    (1 / (4 * a^2)) *
      (
        var_b -
          (2 * b / a) * cov_ab +
          (b^2 / a^2) * var_a
      )

  se_stat_squared <- 1 / (2 * a)

  se_total_squared <-
    se_mc_squared +
    se_stat_squared

  delta <-
    qchisq(confidence, df = 1) *
      (a * se_mc_squared + 0.5)

  ## ---- CI support ----

  loglik_diff <-
    smooth_ll_max -
    smoothed_loglik

  inside <- loglik_diff < delta

  ## Check whether we actually crossed the threshold
  ## before reaching either edge of the evaluated profile.

  left_supported <- if (i_max > 1) {
    any(!inside[seq_len(i_max - 1)])
  } else {
    FALSE
  }

  right_supported <- if (i_max < length(inside)) {
    any(!inside[(i_max + 1):length(inside)])
  } else {
    FALSE
  }

  if (!left_supported) {
    status <- c(status, "left_unbounded")
  }

  if (!right_supported) {
    status <- c(status, "right_unbounded")
  }

  ## ---- CI endpoints ----
  ##
  ## If the threshold crossing exists, use the MCAP-supported
  ## endpoint.
  ##
  ## If not, use the edge of the explored profile range and
  ## mark that side as truncated.

  if (any(inside)) {
    inside_range <- range(
      parameter_grid[inside],
      na.rm = TRUE
    )
  } else {
    inside_range <- c(
      smooth_arg_max,
      smooth_arg_max
    )
  }

  ci_lower <- if (left_supported) {
    inside_range[1]
  } else {
    explored_range[1]
  }

  ci_upper <- if (right_supported) {
    inside_range[2]
  } else {
    explored_range[2]
  }

  ci <- c(
    lower = ci_lower,
    upper = ci_upper
  )

  ci_truncated <- c(
    lower = !left_supported,
    upper = !right_supported
  )

  ## CI interpretation
  ##
  ## If at least one side is truncated, the MCAP confidence
  ## set extends beyond the explored parameter range.

  ci_type <- if (
    left_supported &&
      right_supported
  ) {
    "mcap"
  } else if (
    !left_supported &&
      !right_supported
  ) {
    "mcap_beyond_both_ranges"
  } else if (!left_supported) {
    "mcap_beyond_lower_range"
  } else {
    "mcap_beyond_upper_range"
  }

  ## A concave fit gives a valid MCAP construction.
  ## "Valid" here does not mean both limits were numerically
  ## observed within the explored domain.
  ci_valid <- TRUE

  if (length(status) == 0) {
    status <- "ok"
  }

  ## ---- Return ----

  list(
    lp = lp,
    parameter = parameter,
    confidence = confidence,
    quadratic_fit = quadratic_fit,
    quadratic_stationary = quadratic_stationary,
    quadratic_max = quadratic_max,
    smooth_fit = smooth_fit,
    fit = data.frame(
      parameter = parameter_grid,
      smoothed = smoothed_loglik,
      quadratic = quadratic_pred
    ),
    mle = smooth_arg_max,
    ci = ci,
    ci_truncated = ci_truncated,
    ci_valid = ci_valid,
    ci_type = ci_type,
    explored_range = explored_range,
    delta = delta,
    a = a,
    b = b,
    concave = TRUE,
    se_stat = sqrt(se_stat_squared),
    se_mc = sqrt(se_mc_squared),
    se = sqrt(se_total_squared),
    peak_at_edge = peak_at_edge,
    n_unique_local = n_unique_local,
    mle_gap = mle_gap,
    status = unique(status)
  )
}
