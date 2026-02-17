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
