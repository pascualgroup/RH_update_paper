#' Compute a Monte Carlo Adjusted Profile (MCAP) confidence interval
#'
#' Implements the MCAP method (Ionides et al.) for turning a noisy
#' Monte-Carlo profile log-likelihood into a smoothed profile and a
#' confidence interval for one parameter. The profile is smoothed with
#' loess, a weighted local quadratic is fit near the smoothed maximum, and
#' the resulting curvature is used to combine statistical and Monte-Carlo
#' (simulation) uncertainty into an inflated log-likelihood cutoff.
#'
#' @param lp Numeric vector of profile log-likelihood values.
#' @param parameter Numeric vector of parameter values corresponding to
#'   `lp` (same length).
#' @param confidence Confidence level for the interval.
#' @param lambda Loess span, and fraction of points (nearest to the
#'   smoothed maximum) used to weight the local quadratic fit.
#' @param Ngrid Number of points in the grid used to evaluate the loess
#'   smooth and locate the confidence interval.
#'
#' @return A list with elements: `lp`, `parameter`, `confidence` (the
#'   inputs); `quadratic_fit` (the weighted local `lm` fit); `quadratic_max`
#'   (vertex of the quadratic fit); `smooth_fit` (the `loess` object);
#'   `fit` (data frame with columns `parameter`, `smoothed`, `quadratic`
#'   giving the loess and quadratic predictions on the grid); `mle`
#'   (parameter value maximizing the smoothed profile); `ci` (two-element
#'   vector, the MCAP confidence interval); `delta` (the log-likelihood
#'   drop used as the CI cutoff); `se_stat`, `se_mc`, `se` (statistical,
#'   Monte Carlo, and total standard errors of the estimate).
#'
#' @details The local quadratic is parameterized as
#'   `lp = c - a * parameter^2 + b * parameter`, so its curvature `a` and
#'   the sampling variance of `a` and `b` are used to estimate both the
#'   Monte Carlo variance (`se_mc`) coming from simulation noise in `lp`
#'   and the usual statistical variance (`se_stat`). `delta` inflates the
#'   classical `qchisq(confidence, df = 1) / 2` cutoff by a term
#'   proportional to `se_mc` so the resulting `ci` accounts for both
#'   sources of uncertainty.
#' @export
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

#' Compute an MCAP profile using a classical, non-Monte-Carlo-adjusted cutoff
#'
#' A variant of [mcap()] that smooths and fits the local quadratic the same
#' way, but ignores the Monte Carlo variance term when defining the
#' confidence interval. Use this when the profile is believed to be
#' essentially noise-free (e.g. profile points already well-converged) and
#' a plain likelihood-ratio cutoff is preferred over the MC-inflated one.
#'
#' @param lp Numeric vector of profile log-likelihood values.
#' @param parameter Numeric vector of parameter values corresponding to
#'   `lp` (same length).
#' @param confidence Confidence level for the interval.
#' @param lambda Loess span, and fraction of points (nearest to the
#'   smoothed maximum) used to weight the local quadratic fit.
#' @param Ngrid Number of points in the grid used to evaluate the loess
#'   smooth and locate the confidence interval.
#'
#' @return A list with the same elements as [mcap()]: `lp`, `parameter`,
#'   `confidence`, `quadratic_fit`, `quadratic_max`, `smooth_fit`, `fit`,
#'   `mle`, `ci`, `delta`, `se_stat`, `se_mc`, `se`.
#'
#' @details Differs from [mcap()] in two ways. First, the loess smoother is
#'   fit with `loess.control(surface = "direct")`, computing exact
#'   predictions instead of the (faster, approximate) default interpolation
#'   surface. Second, and more importantly, `delta` is set to the plain
#'   `qchisq(confidence, df = 1) / 2` cutoff rather than [mcap()]'s
#'   MC-inflated cutoff — `se_mc` and `se_stat` are still computed and
#'   returned, but `se_mc` is not folded into `delta`, so the resulting
#'   `ci` is a classical profile-likelihood interval that does not account
#'   for Monte Carlo noise in `lp`.
#' @export
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

#' Plot an MCAP profile with its confidence interval
#'
#' Builds a diagnostic plot of a profile log-likelihood together with its
#' MCAP smoothing and confidence interval, as returned by [mcap()],
#' [simple_mcap()], or the concave case of [mcap_checked()].
#'
#' @param mcap_output A list as returned by [mcap()] / [simple_mcap()] /
#'   [mcap_checked()], containing at least `parameter`, `lp`, `fit`, `ci`,
#'   and `delta`.
#' @param font_size Base font size passed to `theme_bw()`.
#'
#' @return A ggplot object showing the loess-smoothed profile (blue line),
#'   the raw profile log-likelihood points, the maximizing point
#'   (highlighted in green), the confidence-interval bounds (vertical
#'   dashed lines), and the log-likelihood cutoff used to define them
#'   (horizontal line). The plot title reports the MLE and CI bounds.
#'
#' @details Infinite CI bounds (as can occur when the cutoff is never
#'   crossed within the explored grid) are replaced with the min/max of the
#'   observed parameter values so the plot always has finite limits.
#' @export
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

#' Bin a dataset into quantile groups of a column
#'
#' Cuts a numeric column into `n_quantiles` equal-probability quantile bins
#' and adds the bin membership and numeric bin edges as new columns. Used by
#' [plot_mcap_binned()] and [plot_simple_CI()] to thin dense profile output
#' before fitting MCAP.
#'
#' @param dataset A data frame or tibble.
#' @param column_name Character scalar naming the numeric column to bin.
#' @param n_quantiles Number of quantile bins to create.
#'
#' @return `dataset` with additional columns: `quantile` (a factor giving
#'   the bin interval label from `cut()`), and `start`, `end` (numeric
#'   lower/upper edges of that interval, parsed out of the `quantile`
#'   label).
#' @export
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

#' Build an MCAP profile-likelihood plot from binned, thinned profile output
#'
#' Thins a raw profile-likelihood dataset (typically many Monte Carlo
#' replicates per parameter value) by binning `column_name` into quantiles
#' and keeping only the highest-loglik rows per bin, then fits [mcap()] to
#' the thinned data and plots the result with [mcap_plot()].
#'
#' @param dataset A data frame with (at least) columns `loglik` and
#'   `column_name`.
#' @param column_name Character scalar naming the parameter column to
#'   profile over.
#' @param n_quantiles Number of quantile bins passed to
#'   [bin_by_quantile()].
#' @param max_per_bin Maximum number of highest-`loglik` rows kept per
#'   quantile bin.
#' @param diff_LL Optional numeric; if supplied, rows with `loglik` more
#'   than `diff_LL` below the maximum `loglik` are dropped before binning.
#' @param remove_outlier Logical; whether to filter outlying values of
#'   `column_name` after binning/thinning.
#' @param zscore Logical; if `TRUE` outliers are removed using an absolute
#'   z-score threshold of 2, otherwise using 1.5*IQR fences.
#'
#' @return A ggplot object (from [mcap_plot()]) with `column_name` as the
#'   x-axis label and `"loglikelihood"` as the y-axis label.
#'
#' @details Pipeline: drop incomplete rows, optionally drop rows far below
#'   the maximum `loglik` (`diff_LL`), bin by quantile of `column_name`
#'   ([bin_by_quantile()]), keep the top `max_per_bin` rows by `loglik` in
#'   each bin, optionally remove outliers in `column_name`, then call
#'   [mcap()] on the thinned `(loglik, column_name)` pairs and plot the
#'   result with [mcap_plot()].
#' @export
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

#' Build a classical (non-MC-adjusted) profile-likelihood plot from binned output
#'
#' Same thinning-and-plotting idea as [plot_mcap_binned()], but fits
#' [simple_mcap()] instead of [mcap()], and removes outliers in
#' `column_name` before binning rather than after.
#'
#' @param dataset A data frame with (at least) columns `loglik` and
#'   `column_name`.
#' @param column_name Character scalar naming the parameter column to
#'   profile over.
#' @param n_quantiles Number of quantile bins passed to
#'   [bin_by_quantile()].
#' @param max_per_bin Maximum number of highest-`loglik` rows kept per
#'   quantile bin.
#' @param diff_LL Optional numeric; if supplied, rows with `loglik` more
#'   than `diff_LL` below the maximum `loglik` are dropped before outlier
#'   removal and binning.
#' @param remove_outlier Logical; whether to filter outlying values of
#'   `column_name` before binning.
#' @param zscore Logical; if `TRUE` outliers are removed using an absolute
#'   z-score threshold of 2, otherwise using 1.5*IQR fences.
#'
#' @return A ggplot object (from [mcap_plot()]) with `column_name` as the
#'   x-axis label and `"loglikelihood"` as the y-axis label.
#'
#' @details Pipeline: drop incomplete rows, optionally drop rows far below
#'   the maximum `loglik` (`diff_LL`), optionally remove outliers in
#'   `column_name`, bin by quantile of `column_name`
#'   ([bin_by_quantile()]) and keep the top `max_per_bin` rows by `loglik`
#'   in each bin, then call [simple_mcap()] on the thinned pairs. The
#'   resulting `delta` is then explicitly overwritten with
#'   `qchisq(0.95, df = 1) / 2` (the classical 95% cutoff regardless of
#'   any other confidence level) before plotting with [mcap_plot()].
#' @export
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

#' Compute an MCAP profile with diagnostics and safeguards against degenerate profiles
#'
#' A defensive variant of [mcap()] that sanitizes its inputs and adds
#' diagnostics so that pathological profiles (a non-concave local quadratic
#' fit, a maximum at the edge of the explored range, too few distinct
#' parameter values near the maximum, or a confidence cutoff never crossed
#' within the explored grid) are reported explicitly via a `status` flag
#' vector instead of silently producing a nonsensical or infinite-bounded
#' confidence interval.
#'
#' @param lp Numeric vector of profile log-likelihood values.
#' @param parameter Numeric vector of parameter values corresponding to
#'   `lp` (same length).
#' @param confidence Confidence level for the interval.
#' @param lambda Loess span, and (for the local quadratic fit) target
#'   fraction of points nearest the smoothed maximum to include.
#' @param Ngrid Number of points in the grid used to evaluate the loess
#'   smooth and locate the confidence interval.
#' @param external_mle Optional externally-known MLE for `parameter`. If
#'   supplied, its smoothed log-likelihood is compared to the profile's
#'   smoothed maximum and the gap is reported as `mle_gap`, as a
#'   diagnostic of agreement between the profile and an independently
#'   obtained MLE.
#' @param edge_fraction Fraction of the observed `parameter` range treated
#'   as "near the boundary"; used to flag (`peak_at_edge`) when the
#'   smoothed maximum falls within this margin of either end.
#' @param min_local_unique Minimum number of distinct `parameter` values
#'   required within the local quadratic-fit neighborhood; below this,
#'   `status` includes `"few_local_x"`.
#'
#' @return A list with the same core elements as [mcap()] (`lp`,
#'   `parameter`, `confidence`, `quadratic_fit`, `quadratic_max`,
#'   `smooth_fit`, `fit`, `mle`, `ci`, `delta`, `se_stat`, `se_mc`, `se`),
#'   plus diagnostic fields: `quadratic_stationary` (the quadratic fit's
#'   vertex regardless of concavity; `quadratic_max` is `NA` unless the fit
#'   is concave), `a`, `b` (quadratic fit coefficients), `concave`
#'   (logical), `ci_valid` (`FALSE` when the quadratic fit is not concave,
#'   in which case `ci` is only the explored `parameter` range rather than
#'   a true MCAP interval), `ci_truncated` (named logical vector with
#'   elements `lower`/`upper`, flagging sides where the log-likelihood
#'   cutoff was never crossed within the explored grid, so the edge of the
#'   explored range was used instead), `ci_type` (one of `"mcap"`,
#'   `"mcap_beyond_lower_range"`, `"mcap_beyond_upper_range"`,
#'   `"mcap_beyond_both_ranges"`, or `"explored_range"`), `explored_range`
#'   (range of the valid grid points), `peak_at_edge`, `n_unique_local`,
#'   `mle_gap`, and `status` (character vector of diagnostic flags such as
#'   `"non_concave"`, `"peak_at_edge"`, `"few_local_x"`,
#'   `"left_unbounded"`, `"right_unbounded"`, or `"ok"`).
#'
#' @details Non-finite `(lp, parameter)` pairs are dropped up front. The
#'   local quadratic is parameterized as
#'   `lp = c - a * parameter^2 + b * parameter`, so `a > 0` is required for
#'   it to be concave (a genuine local maximum) and for the MCAP variance
#'   and cutoff formulas to be valid. If the fit is not concave, the
#'   function returns early: `ci` is set to the explored `parameter` range,
#'   `ci_valid` is `FALSE`, and `status` includes `"non_concave"`. If the
#'   fit is concave, `delta` is computed the same way as in [mcap()]
#'   (MC-adjusted cutoff). The function then checks, separately on each
#'   side of the maximum, whether the log-likelihood actually drops below
#'   `delta` before the edge of the evaluated grid is reached; if not, that
#'   side of `ci` is set to the edge of the explored range instead of being
#'   left unbounded, and is flagged via `ci_truncated`/`status`
#'   (`"left_unbounded"`/`"right_unbounded"`) and reflected in `ci_type`.
#' @export
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
