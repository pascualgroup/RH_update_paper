#' Generate a logarithmically spaced sequence
#'
#' Small generic helper that returns a sequence of `length.out` values spaced
#' evenly on a log scale between `from` and `to`. Adapted from the `emdbook`
#' package (see inline comment) rather than depending on it for this single
#' function. Used here to place color-scale breaks on wavelet power plots.
#'
#' @param from Starting value of the sequence (must be positive).
#' @param to Ending value of the sequence (must be positive).
#' @param length.out Number of points to generate.
#' @return A numeric vector of length `length.out`, log-spaced between `from`
#'   and `to`.
#' @keywords internal
lseq <- function(from = 1, to = 100000, length.out = 6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

#' Compute the hour of the year for a date-time
#'
#' Convert a date-time into an integer hour-of-year count, 0-indexed from
#' midnight on January 1st of that year (i.e. Jan 1st 00:00 is hour 0).
#'
#' @param date_hour A date-time value (or vector of date-times).
#' @return A numeric vector of hour-of-year values.
#' @keywords internal
yhour <- function(date_hour) {
  day_of_year <- yday(date_hour)
  # Get the hour of the day
  hour_of_day <- hour(date_hour)
  hour_of_year <- (day_of_year - 1) * 24 + hour_of_day
  return(hour_of_year)
}

#' Number of days in a year
#'
#' Return 366 for leap years and 365 otherwise.
#'
#' @param year A numeric year (or vector of years).
#' @return A numeric vector: 365 or 366 for each element of `year`.
#' @keywords internal
days_in_year <- function(year) {
  result <- if_else(leap_year(year), 366, 365)
  return(result)
}

#' Number of hours in a year
#'
#' Return the number of hours in `year`, accounting for leap years via
#' [days_in_year()].
#'
#' @param year A numeric year (or vector of years).
#' @return A numeric vector: `days_in_year(year) * 24`.
#' @keywords internal
hours_in_year <- function(year) {
  result <- days_in_year(year) * 24
  return(result)
}

#' Number of minutes in a year
#'
#' Return the number of minutes in `year`, accounting for leap years via
#' [hours_in_year()].
#'
#' @param year A numeric year (or vector of years).
#' @return A numeric vector: `hours_in_year(year) * 60`.
#' @keywords internal
mins_in_year <- function(year) {
  result <- hours_in_year(year) * 60
  return(result)
}

#' Plot a wavelet power spectrum
#'
#' Build a ggplot raster of wavelet power over time and period from one or
#' more `WaveletComp::analyze.wavelet()` results, overlaid with a ridge
#' contour (0.95), a significance contour (p = 0.1), and cone-of-influence
#' shading. When several wavelet objects are supplied they are combined into
#' a single faceted plot rather than a list of separate plots.
#'
#' @param wc_list Either a single wavelet object as returned by
#'   `WaveletComp::analyze.wavelet()`, an unnamed list of such objects, or a
#'   named list whose elements are each a list of such objects. See Details
#'   for how the shape of `wc_list` is detected and what it controls.
#'
#' @details
#' The accepted shape of `wc_list` is detected via `vctrs::vec_depth()`:
#' \itemize{
#'   \item depth 3: `wc_list` is itself a single wavelet object. The result is
#'     one unfaceted power/period panel.
#'   \item depth 4: `wc_list` is an unnamed list of wavelet objects (e.g. one
#'     per series). Each object's power, ridge, p-value and cone-of-influence
#'     data are stacked and combined into one plot with
#'     `facet_grid(series ~ .)`.
#'   \item depth 5: `wc_list` is a named list whose elements are each a list
#'     of wavelet objects (e.g. one list per group, each holding one wavelet
#'     object per series). The outer names become the `group` facet variable
#'     and the plot uses `facet_grid(group ~ series)`.
#'   \item any other shape triggers
#'     `stop("wc_list is not list(wavelet) or list(list(wavelet)))")`.
#' }
#' In every case the period axis is log2-scaled and the fill scale uses
#' `scico::scale_fill_scico(palette = "vik")` on a log1p-transformed power
#' value.
#'
#' @return A single ggplot object.
#' @export
plot_wavelet <- function(wc_list) {
  all_power_vals <- c()
  all_ridge_vals <- c()
  all_pval_vals <- c()
  all_coi <- c()
  all_result <- c()
  if (wc_list |> vec_depth() == 3) {
    p <- with(wc_list, {
      pval_vals <- Power.pval |> as_tibble()
      ridge_vals <- Ridge |> as_tibble()
      power_vals <- Power |> as_tibble()
      colnames(power_vals) <- series |> pull("date")
      colnames(ridge_vals) <- series |> pull("date")
      colnames(pval_vals) <- series |> pull("date")

      series_name <- series |>
        select(-date) |>
        select(-contains("trend")) |>
        colnames()
      series_name <- series_name[1]
      time_series <- series |> pull("date")


      power_vals <- power_vals |>
        mutate(Period = Period) |>
        pivot_longer(-Period, names_to = "date", values_to = "Power") |>
        mutate(date = as.Date(date)) |>
        mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date))

      ridge_vals <- ridge_vals |>
        mutate(Period = Period) |>
        pivot_longer(-Period, names_to = "date", values_to = "Ridge") |>
        mutate(date = as.Date(date)) |>
        mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date))

      pval_vals <- pval_vals |>
        mutate(Period = Period) |>
        pivot_longer(-Period, names_to = "date", values_to = "pval") |>
        mutate(date = as.Date(date)) |>
        mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date))

      coi <- tibble(x = coi.1, y = 2**coi.2) |>
        mutate(y = if_else(y < min(power_vals$Period), min(power_vals$Period), y)) |>
        mutate(y = if_else(y > max(power_vals$Period), max(power_vals$Period), y))

      min_date <- min(time_series) |> as.Date()
      min_date <- year(min_date) + (month(min_date) - 1) / 12 + (day(min_date) - 1) / days_in_month(min_date)

      coi <- coi |> mutate(date = min_date + x - 1)

      power_vals |> ggplot(aes(x = date, y = Period, fill = Power)) +
        geom_raster(width = 1, height = 1) +
        scale_y_continuous(
          trans = "log2", breaks = 2**(-2:4), expand = c(0, 0),
          labels = scales::label_number(drop0trailing = TRUE)
        ) +
        scale_x_continuous(
          expand = c(0, 0),
          labels = ~ format(make_date(floor(.x), round(12 * (.x - floor(.x)) + 1), 1), "%b %Y")
        ) +
        scico::scale_fill_scico(
          palette = "vik", trans = "log1p",
          breaks = c(0, round(lseq(min(power_vals$Power) + 1, max(power_vals$Power) + 1, 5) - 1, digits = 1)),
          limits = c(0, NA), midpoint = log(2)
        ) +
        geom_contour(aes(x = date, y = Period, z = Ridge),
          inherit.aes = FALSE, data = ridge_vals,
          color = "black", linewidth = 1, fill = "black", breaks = 0.95
        ) +
        geom_contour(aes(x = date, y = Period, z = pval),
          inherit.aes = FALSE, data = pval_vals,
          color = "white", linewidth = 0.5, fill = NA, breaks = 0.1
        ) +
        geom_polygon(aes(x = date, y = y), data = coi, inherit.aes = F, fill = "white", alpha = 0.5) +
        theme_minimal(22) +
        labs(x = "Date", y = "Period") +
        theme(
          legend.key.height = unit(4, "lines"),
          text = element_text(family = "Times")
        )
    })
    return(p)
  } else if (wc_list |> vec_depth() == 4) {
    for (wc in wc_list) {
      p <- with(wc, {
        pval_vals <- Power.pval |> as_tibble()
        ridge_vals <- Ridge |> as_tibble()
        power_vals <- Power |> as_tibble()
        colnames(power_vals) <- series |> pull("date")
        colnames(ridge_vals) <- series |> pull("date")
        colnames(pval_vals) <- series |> pull("date")

        series_name <- series |>
          select(-date) |>
          select(-contains("trend")) |>
          colnames()
        series_name <- series_name[1]
        time_series <- series |> pull("date")


        power_vals <- power_vals |>
          mutate(Period = Period) |>
          pivot_longer(-Period, names_to = "date", values_to = "Power") |>
          mutate(date = as.Date(date)) |>
          mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
          mutate(series = series_name)

        ridge_vals <- ridge_vals |>
          mutate(Period = Period) |>
          pivot_longer(-Period, names_to = "date", values_to = "Ridge") |>
          mutate(date = as.Date(date)) |>
          mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
          mutate(series = series_name)

        pval_vals <- pval_vals |>
          mutate(Period = Period) |>
          pivot_longer(-Period, names_to = "date", values_to = "pval") |>
          mutate(date = as.Date(date)) |>
          mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
          mutate(series = series_name)

        coi <- tibble(x = coi.1, y = 2**coi.2)

        min_date <- min(time_series) |> as.Date()
        min_date <- year(min_date) + (month(min_date) - 1) / 12 + (day(min_date) - 1) / days_in_month(min_date)

        coi <- coi |>
          mutate(date = min_date + x - 1) |>
          mutate(series = series_name)


        result <- list(power = power_vals, ridge = ridge_vals, pval = pval_vals, coi = coi)
      })
      all_power_vals <- all_power_vals |> rbind(p$power)
      all_ridge_vals <- all_ridge_vals |> rbind(p$ridge)
      all_pval_vals <- all_pval_vals |> rbind(p$pval)
      all_coi <- all_coi |> rbind(p$coi)
    }
    all_result <- list(power = all_power_vals, ridge = all_ridge_vals, pval = all_pval_vals, coi = all_coi)
    p <- with(all_result, {
      coi <- coi |>
        mutate(y = if_else(y < min(power$Period), min(power$Period), y)) |>
        mutate(y = if_else(y > max(power$Period), max(power$Period), y))

      power |> ggplot(aes(x = date, y = Period, fill = Power)) +
        facet_grid(series ~ .) +
        geom_raster(width = 1, height = 1) +
        scale_y_continuous(
          trans = "log2", breaks = 2**(-2:4), expand = c(0, 0),
          labels = scales::label_number(drop0trailing = TRUE)
        ) +
        scale_x_continuous(
          expand = c(0, 0),
          labels = ~ format(make_date(floor(.x), round(12 * (.x - floor(.x)) + 1), 1), "%b %Y")
        ) +
        scico::scale_fill_scico(
          palette = "vik", trans = "log1p",
          breaks = c(0, round(lseq(min(power$Power) + 1, max(power$Power) + 1, 5) - 1, digits = 1)),
          limits = c(0, NA), midpoint = log(2)
        ) +
        geom_contour(aes(x = date, y = Period, z = Ridge),
          inherit.aes = FALSE, data = ridge,
          color = "black", linewidth = 1, fill = "black", breaks = 0.95
        ) +
        geom_contour(aes(x = date, y = Period, z = pval),
          inherit.aes = FALSE, data = pval,
          color = "white", linewidth = 0.5, fill = NA, breaks = 0.1
        ) +
        geom_polygon(aes(x = date, y = y), data = coi, inherit.aes = F, fill = "white", alpha = 0.5) +
        theme_minimal(22) +
        labs(x = "Date", y = "Period") +
        theme(
          legend.key.height = unit(4, "lines"),
          text = element_text(family = "Times")
        )
    })
    return(p)
  } else if (wc_list |> vec_depth() == 5) {
    for (i in 1:length(wc_list)) {
      group_name <- names(wc_list)[i]
      wc_list2 <- wc_list[[i]]
      for (wc in wc_list2) {
        p <- with(wc, {
          pval_vals <- Power.pval |> as_tibble()
          ridge_vals <- Ridge |> as_tibble()
          power_vals <- Power |> as_tibble()
          colnames(power_vals) <- series |> pull("date")
          colnames(ridge_vals) <- series |> pull("date")
          colnames(pval_vals) <- series |> pull("date")

          series_name <- series |>
            select(-date) |>
            select(-contains("trend")) |>
            colnames()
          series_name <- series_name[1]
          time_series <- series |> pull("date")


          power_vals <- power_vals |>
            mutate(Period = Period) |>
            pivot_longer(-Period, names_to = "date", values_to = "Power") |>
            mutate(date = as.Date(date)) |>
            mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
            mutate(series = series_name, group = group_name)

          ridge_vals <- ridge_vals |>
            mutate(Period = Period) |>
            pivot_longer(-Period, names_to = "date", values_to = "Ridge") |>
            mutate(date = as.Date(date)) |>
            mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
            mutate(series = series_name, group = group_name)

          pval_vals <- pval_vals |>
            mutate(Period = Period) |>
            pivot_longer(-Period, names_to = "date", values_to = "pval") |>
            mutate(date = as.Date(date)) |>
            mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
            mutate(series = series_name, group = group_name)

          coi <- tibble(x = coi.1, y = 2**coi.2)

          min_date <- min(time_series) |> as.Date()
          min_date <- year(min_date) + (month(min_date) - 1) / 12 + (day(min_date) - 1) / days_in_month(min_date)

          coi <- coi |>
            mutate(date = min_date + x - 1) |>
            mutate(series = series_name, group = group_name)


          result <- list(power = power_vals, ridge = ridge_vals, pval = pval_vals, coi = coi)
        })
        all_power_vals <- all_power_vals |> rbind(p$power)
        all_ridge_vals <- all_ridge_vals |> rbind(p$ridge)
        all_pval_vals <- all_pval_vals |> rbind(p$pval)
        all_coi <- all_coi |> rbind(p$coi)
      }
    }
    all_result <- list(power = all_power_vals, ridge = all_ridge_vals, pval = all_pval_vals, coi = all_coi)

    p <- with(all_result, {
      coi <- coi |>
        mutate(y = if_else(y < min(power$Period), min(power$Period), y)) |>
        mutate(y = if_else(y > max(power$Period), max(power$Period), y))

      power |> ggplot(aes(x = date, y = Period, fill = Power)) +
        # facet_grid(interaction(group, series, sep = ": ") ~ .) +
        facet_grid(group ~ series) +
        geom_raster(width = 1, height = 1) +
        scale_y_continuous(
          trans = "log2", breaks = 2**(-2:4), expand = c(0, 0),
          labels = scales::label_number(drop0trailing = TRUE)
        ) +
        scale_x_continuous(
          expand = c(0, 0),
          labels = ~ format(make_date(floor(.x), round(12 * (.x - floor(.x)) + 1), 1), "%b %Y")
        ) +
        scico::scale_fill_scico(
          palette = "vik", trans = "log1p",
          breaks = c(0, round(lseq(min(power$Power) + 1, max(power$Power) + 1, 5) - 1, digits = 1)),
          limits = c(0, NA), midpoint = log(2)
        ) +
        geom_contour(aes(x = date, y = Period, z = Ridge),
          inherit.aes = FALSE, data = ridge,
          color = "black", linewidth = 1, fill = "black", breaks = 0.95
        ) +
        geom_contour(aes(x = date, y = Period, z = pval),
          inherit.aes = FALSE, data = pval,
          color = "white", linewidth = 0.5, fill = NA, breaks = 0.1
        ) +
        geom_polygon(aes(x = date, y = y), data = coi, inherit.aes = F, fill = "white", alpha = 0.5) +
        theme_minimal(22) +
        labs(x = "Date", y = "Period") +
        theme(
          legend.key.height = unit(4, "lines"),
          text = element_text(family = "Times")
        )
    })


    return(p)
  } else {
    stop("wc_list is not list(wavelet) or list(list(wavelet))")
  }
}


#' Plot a cross-wavelet power or wavelet coherence spectrum
#'
#' Build a ggplot raster of cross-wavelet power or wavelet coherence over
#' time and period from one or more `WaveletComp::analyze.coherency()`
#' results, overlaid with phase-difference arrows, a significance contour
#' (p = 0.1), and cone-of-influence shading. As in [plot_wavelet()], several
#' coherency objects are combined into a single faceted plot rather than a
#' list of separate plots.
#'
#' @param wc_list Either a single coherency object as returned by
#'   `WaveletComp::analyze.coherency()`, an unnamed list of such objects, or a
#'   named list whose elements are each a list of such objects. The shape is
#'   detected the same way as in [plot_wavelet()] (via `vctrs::vec_depth()`):
#'   depth 3 is a single object (plotted unfaceted), depth 4 is an unnamed
#'   list of objects (combined with `facet_grid(. ~ series)`), and depth 5 is
#'   a named list of lists of objects (combined with
#'   `facet_grid(group ~ series)`). Unlike [plot_wavelet()], other shapes are
#'   not explicitly rejected; they will fail when the plot is built from the
#'   (empty) combined result.
#' @param which.image Character scalar selecting which quantity to plot when
#'   `wc_list` has depth 4 or 5 (i.e. holds more than one coherency object):
#'   `"wp"` (default) uses cross-wavelet power (`Power.xy`/`Power.xy.pval`);
#'   `"wc"` uses wavelet coherence (`Coherence`/`Coherence.pval`). It also
#'   sets the fill legend title (`"Power"` vs `"Coherence"`). It has no effect
#'   when `wc_list` is a single coherency object (depth 3), which is always
#'   plotted as cross-wavelet power.
#'
#' @details
#' Phase-difference arrows are drawn from the coherency object's `sAngle`
#' component via `metR::geom_arrow()`. Before plotting, the angle values are
#' thinned (every 3rd date, every 6th period) and filtered to points where
#' power is positive and the p-value is at most 0.05, to keep the arrow
#' overlay legible.
#'
#' @return A single ggplot object.
#' @export
plot_wavelet_coherence <- function(wc_list, which.image = "wp") {
  all_power_vals <- c()
  all_arrow_vals <- c()
  all_pval_vals <- c()
  all_coi <- c()
  all_result <- c()
  if (wc_list |> vec_depth() == 3) {
    p <- with(wc_list, {
      pval_vals <- Power.xy.pval |> as_tibble()
      arrow_vals <- sAngle |> as_tibble()
      power_vals <- Power.xy |> as_tibble()
      colnames(power_vals) <- series |> pull("date")
      colnames(arrow_vals) <- series |> pull("date")
      colnames(pval_vals) <- series |> pull("date")

      series_name <- series |>
        select(-date) |>
        select(-contains("trend")) |>
        colnames()
      series_name <- paste(series_name[1], series_name[2], sep = ":")
      time_series <- series |> pull("date")


      power_vals <- power_vals |>
        mutate(Period = Period) |>
        pivot_longer(-Period, names_to = "date", values_to = "Power") |>
        mutate(date = as.Date(date)) |>
        mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date))
      pval_vals <- pval_vals |>
        mutate(Period = Period) |>
        pivot_longer(-Period, names_to = "date", values_to = "pval") |>
        mutate(date = as.Date(date)) |>
        mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date))


      arrow_vals <- arrow_vals |>
        mutate(Period = Period) |>
        pivot_longer(-Period, names_to = "date", values_to = "Angle") |>
        mutate(date = as.Date(date)) |>
        mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
        left_join(power_vals) |>
        left_join(pval_vals) |>
        group_by(date) |>
        filter(row_number() %% 3 == 1) |>
        ungroup() |>
        group_by(Period) |>
        filter(row_number() %% 6 == 1) |>
        ungroup() |>
        filter(Power > 0, pval <= 0.05)


      coi <- tibble(x = coi.1, y = 2**coi.2) |>
        mutate(y = if_else(y < min(power_vals$Period), min(power_vals$Period), y)) |>
        mutate(y = if_else(y > max(power_vals$Period), max(power_vals$Period), y))

      min_date <- min(time_series) |> as.Date()
      min_date <- year(min_date) + (month(min_date) - 1) / 12 + (day(min_date) - 1) / days_in_month(min_date)

      coi <- coi |> mutate(date = min_date + x - 1)

      result <- list(power = power_vals, arrow = arrow_vals, pval = pval_vals, coi = coi)
    })
    all_power_vals <- all_power_vals |> rbind(p$power)
    all_arrow_vals <- all_arrow_vals |> rbind(p$arrow)
    all_pval_vals <- all_pval_vals |> rbind(p$pval)
    all_coi <- all_coi |> rbind(p$coi)
    all_result <- list(power = all_power_vals, arrow = all_arrow_vals, pval = all_pval_vals, coi = all_coi)
  } else if (wc_list |> vec_depth() == 4) {
    for (wc in wc_list) {
      p <- with(wc, {
        if (which.image == "wp") {
          pval_vals <- Power.xy.pval |> as_tibble()
          power_vals <- Power.xy |> as_tibble()
        } else if (which.image == "wc") {
          pval_vals <- Coherence.pval |> as_tibble()
          power_vals <- Coherence |> as_tibble()
        }
        arrow_vals <- sAngle |> as_tibble()
        colnames(power_vals) <- series |> pull("date")
        colnames(arrow_vals) <- series |> pull("date")
        colnames(pval_vals) <- series |> pull("date")

        series_name <- series |>
          select(-date) |>
          select(-contains("trend")) |>
          colnames()
        series_name <- paste(series_name[1], series_name[2], sep = ":")
        time_series <- series |> pull("date")


        power_vals <- power_vals |>
          mutate(Period = Period) |>
          pivot_longer(-Period, names_to = "date", values_to = "Power") |>
          mutate(date = as.Date(date)) |>
          mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
          mutate(series = series_name)


        pval_vals <- pval_vals |>
          mutate(Period = Period) |>
          pivot_longer(-Period, names_to = "date", values_to = "pval") |>
          mutate(date = as.Date(date)) |>
          mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
          mutate(series = series_name)


        arrow_vals <- arrow_vals |>
          mutate(Period = Period) |>
          pivot_longer(-Period, names_to = "date", values_to = "Angle") |>
          mutate(date = as.Date(date)) |>
          mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
          left_join(power_vals) |>
          left_join(pval_vals) |>
          group_by(date) |>
          filter(row_number() %% 3 == 1) |>
          ungroup() |>
          group_by(Period) |>
          filter(row_number() %% 6 == 1) |>
          ungroup() |>
          filter(Power > 0, pval <= 0.05) |>
          mutate(series = series_name)

        coi <- tibble(x = coi.1, y = 2**coi.2)

        min_date <- min(time_series) |> as.Date()
        min_date <- year(min_date) + (month(min_date) - 1) / 12 + (day(min_date) - 1) / days_in_month(min_date)

        coi <- coi |>
          mutate(date = min_date + x - 1) |>
          mutate(series = series_name)


        result <- list(power = power_vals, arrow = arrow_vals, pval = pval_vals, coi = coi)
      })
      all_power_vals <- all_power_vals |> rbind(p$power)
      all_arrow_vals <- all_arrow_vals |> rbind(p$arrow)
      all_pval_vals <- all_pval_vals |> rbind(p$pval)
      all_coi <- all_coi |> rbind(p$coi)
    }
    all_result <- list(power = all_power_vals, arrow = all_arrow_vals, pval = all_pval_vals, coi = all_coi)
  } else if (wc_list |> vec_depth() == 5) {
    for (i in 1:length(wc_list)) {
      group_name <- names(wc_list)[i]
      wc_list2 <- wc_list[[i]]
      for (wc in wc_list2) {
        p <- with(wc, {
          if (which.image == "wp") {
            pval_vals <- Power.xy.pval |> as_tibble()
            power_vals <- Power.xy |> as_tibble()
          } else if (which.image == "wc") {
            pval_vals <- Coherence.pval |> as_tibble()
            power_vals <- Coherence |> as_tibble()
          }
          arrow_vals <- sAngle |> as_tibble()
          colnames(power_vals) <- series |> pull("date")
          colnames(arrow_vals) <- series |> pull("date")
          colnames(pval_vals) <- series |> pull("date")

          series_name <- series |>
            select(-date) |>
            select(-contains("trend")) |>
            colnames()
          series_name <- paste(series_name[1], series_name[2], sep = ":")
          time_series <- series |> pull("date")


          power_vals <- power_vals |>
            mutate(Period = Period) |>
            pivot_longer(-Period, names_to = "date", values_to = "Power") |>
            mutate(date = as.Date(date)) |>
            mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
            mutate(series = series_name, group = group_name)


          pval_vals <- pval_vals |>
            mutate(Period = Period) |>
            pivot_longer(-Period, names_to = "date", values_to = "pval") |>
            mutate(date = as.Date(date)) |>
            mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
            mutate(series = series_name, group = group_name)

          arrow_vals <- arrow_vals |>
            mutate(Period = Period) |>
            pivot_longer(-Period, names_to = "date", values_to = "Angle") |>
            mutate(date = as.Date(date)) |>
            mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
            left_join(power_vals) |>
            left_join(pval_vals) |>
            group_by(date) |>
            filter(row_number() %% 3 == 1) |>
            ungroup() |>
            group_by(Period) |>
            filter(row_number() %% 6 == 1) |>
            ungroup() |>
            filter(Power > 0, pval <= 0.05) |>
            mutate(series = series_name, group = group_name)

          coi <- tibble(x = coi.1, y = 2**coi.2)

          min_date <- min(time_series) |> as.Date()
          min_date <- year(min_date) + (month(min_date) - 1) / 12 + (day(min_date) - 1) / days_in_month(min_date)

          coi <- coi |>
            mutate(date = min_date + x - 1) |>
            mutate(series = series_name, group = group_name)


          result <- list(power = power_vals, arrow = arrow_vals, pval = pval_vals, coi = coi)
        })
        all_power_vals <- all_power_vals |> rbind(p$power)
        all_arrow_vals <- all_arrow_vals |> rbind(p$arrow)
        all_pval_vals <- all_pval_vals |> rbind(p$pval)
        all_coi <- all_coi |> rbind(p$coi)
      }
    }
    all_result <- list(power = all_power_vals, arrow = all_arrow_vals, pval = all_pval_vals, coi = all_coi)
  }

  p <- with(all_result, {
    if (which.image == "wp") {
      break_vals <- c(0, round(lseq(min(power$Power) + 1, max(power$Power) + 1, 5) - 1, digits = 1))
    } else {
      break_vals <- round(c(seq(min(power$Power), max(power$Power), length.out = 6)), digits = 1)
      break_vals <- waiver()
    }
    coi <- coi |>
      mutate(y = if_else(y < min(power$Period), min(power$Period), y)) |>
      mutate(y = if_else(y > max(power$Period), max(power$Period), y))

    power |> ggplot(aes(x = date, y = Period, fill = Power)) +
      geom_raster(width = 1, height = 1) +
      scale_y_continuous(
        trans = "log2", breaks = 2**(-2:4), expand = c(0, 0),
        labels = scales::label_number(drop0trailing = TRUE)
      ) +
      scale_x_continuous(
        expand = c(0, 0),
        labels = ~ format(make_date(floor(.x), round(12 * (.x - floor(.x)) + 1), 1), "%Y"),
        n.breaks = 3
      ) +
      metR::geom_arrow(
        aes(
          x = date, y = Period, mag = 0.4, angle = 180 * Angle / pi
        ),
        inherit.aes = FALSE,
        data = arrow, size = 0.4, arrow.length = 0.3
      ) +
      geom_contour(aes(x = date, y = Period, z = pval),
        inherit.aes = FALSE, data = pval,
        color = "white", linewidth = 0.5, fill = NA, breaks = 0.1
      ) +
      geom_polygon(aes(x = date, y = y), data = coi, inherit.aes = F, fill = "white", alpha = 0.5) +
      theme_minimal(22) +
      labs(x = "Time", y = "Period", fill = if (which.image == "wc") "Coherence" else "Power") +
      theme(
        legend.key.height = unit(4, "lines"),
        text = element_text(family = "Times")
      )
  })
  vals <- quantile(all_result$power$Power, probs = seq(0, 1, length.out = 256)) |> as.numeric()
  if (wc_list |> vec_depth() == 3) {} else if (wc_list |> vec_depth() == 4) {
    p <- p + facet_grid(. ~ series)
  } else if (wc_list |> vec_depth() == 5) {
    p <- p + facet_grid(group ~ series)
  }

  if (which.image == "wp") {
    p <- p + scale_fill_gradientn(
      colours = rev(colorRampPalette(brewer.pal(9, "RdBu"))(256)),
      # values = scales::rescale(log(vals + 1)),
      values = scales::rescale(vals),
      limits = c(0, NA)
    )
  } else {
    p <- p + scale_fill_gradientn(
      colours = rev(colorRampPalette(brewer.pal(9, "RdBu"))(256)),
      values = scales::rescale(vals),
      limits = c(0, NA)
    )
  }

  return(p)
}
