library(tidyverse)
library(WaveletComp)

lseq <- function(from = 1, to = 100000, length.out = 6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
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

# plot_wavelet_coherence <- function(wc_list, which.image = "wp") {
#   all_power_vals <- c()
#   all_arrow_vals <- c()
#   all_pval_vals <- c()
#   all_coi <- c()
#   all_result <- c()

#   if (which.image == "wp") {
#     transf <- "log1p"
#     midpoint_val <- log(2)
#   } else {
#     transf <- "identity"
#     midpoint_val <- log(1.5)
#   }
#   if (wc_list |> vec_depth() == 3) {
#     p <- with(wc_list, {
#       if (which.image == "wp") {
#         pval_vals <- Power.xy.pval |> as_tibble()
#         power_vals <- Power.xy |> as_tibble()
#       } else if (which.image == "wc") {
#         pval_vals <- Coherence.pval |> as_tibble()
#         power_vals <- Coherence |> as_tibble()
#       }
#       arrow_vals <- sAngle |> as_tibble()
#       colnames(power_vals) <- series |> pull("date")
#       colnames(arrow_vals) <- series |> pull("date")
#       colnames(pval_vals) <- series |> pull("date")

#       series_name <- series |>
#         select(-date) |>
#         select(-contains("trend")) |>
#         colnames()
#       series_name <- paste(series_name[1], series_name[2], sep = ":")
#       time_series <- series |> pull("date")


#       power_vals <- power_vals |>
#         mutate(Period = Period) |>
#         pivot_longer(-Period, names_to = "date", values_to = "Power") |>
#         mutate(date = as.Date(date)) |>
#         mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date))



#       pval_vals <- pval_vals |>
#         mutate(Period = Period) |>
#         pivot_longer(-Period, names_to = "date", values_to = "pval") |>
#         mutate(date = as.Date(date)) |>
#         mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date))

#       arrow_vals <- arrow_vals |>
#         mutate(Period = Period) |>
#         pivot_longer(-Period, names_to = "date", values_to = "Angle") |>
#         mutate(date = as.Date(date)) |>
#         mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
#         left_join(power_vals) |>
#         left_join(pval_vals) |>
#         group_by(date) |>
#         filter(row_number() %% 3 == 1) |>
#         ungroup() |>
#         group_by(Period) |>
#         filter(row_number() %% 6 == 1) |>
#         ungroup() |>
#         filter(Power > 0, pval <= 0.05)

#       if (which.image == "wp") {
#         break_vals <- c(0, round(lseq(min(power_vals$Power) + 1, max(power_vals$Power) + 1, 5) - 1, digits = 1))
#       } else {
#         break_vals <- round(seq(min(power_vals$Power), max(power_vals$Power), length.out = 6), digits = 1)
#         midpoint_val <- quantile(power$Power, 0.7) |> as.numeric()
#       }

#       coi <- tibble(x = coi.1, y = 2**coi.2) |>
#         mutate(y = if_else(y < min(power_vals$Period), min(power_vals$Period), y)) |>
#         mutate(y = if_else(y > max(power_vals$Period), max(power_vals$Period), y))

#       min_date <- min(time_series) |> as.Date()
#       min_date <- year(min_date) + (month(min_date) - 1) / 12 + (day(min_date) - 1) / days_in_month(min_date)

#       coi <- coi |> mutate(date = min_date + x - 1)


#   } else if (wc_list |> vec_depth() == 4) {
#     for (wc in wc_list) {
#       p <- with(wc, {
#         if (which.image == "wp") {
#           pval_vals <- Power.xy.pval |> as_tibble()
#           power_vals <- Power.xy |> as_tibble()
#         } else if (which.image == "wc") {
#           pval_vals <- Coherence.pval |> as_tibble()
#           power_vals <- Coherence |> as_tibble()
#         }
#         arrow_vals <- sAngle |> as_tibble()
#         colnames(power_vals) <- series |> pull("date")
#         colnames(arrow_vals) <- series |> pull("date")
#         colnames(pval_vals) <- series |> pull("date")

#         series_name <- series |>
#           select(-date) |>
#           select(-contains("trend")) |>
#           colnames()
#         series_name <- paste(series_name[1], series_name[2], sep = ":")
#         time_series <- series |> pull("date")


#         power_vals <- power_vals |>
#           mutate(Period = Period) |>
#           pivot_longer(-Period, names_to = "date", values_to = "Power") |>
#           mutate(date = as.Date(date)) |>
#           mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
#           mutate(series = series_name)


#         pval_vals <- pval_vals |>
#           mutate(Period = Period) |>
#           pivot_longer(-Period, names_to = "date", values_to = "pval") |>
#           mutate(date = as.Date(date)) |>
#           mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
#           mutate(series = series_name)


#         arrow_vals <- arrow_vals |>
#           mutate(Period = Period) |>
#           pivot_longer(-Period, names_to = "date", values_to = "Angle") |>
#           mutate(date = as.Date(date)) |>
#           mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
#           left_join(power_vals) |>
#           left_join(pval_vals) |>
#           group_by(date) |>
#           filter(row_number() %% 3 == 1) |>
#           ungroup() |>
#           group_by(Period) |>
#           filter(row_number() %% 6 == 1) |>
#           ungroup() |>
#           filter(Power > 0, pval <= 0.05) |>
#           mutate(series = series_name)

#         coi <- tibble(x = coi.1, y = 2**coi.2)

#         min_date <- min(time_series) |> as.Date()
#         min_date <- year(min_date) + (month(min_date) - 1) / 12 + (day(min_date) - 1) / days_in_month(min_date)

#         coi <- coi |>
#           mutate(date = min_date + x - 1) |>
#           mutate(series = series_name)


#         result <- list(power = power_vals, arrow = arrow_vals, pval = pval_vals, coi = coi)
#       })
#       all_power_vals <- all_power_vals |> rbind(p$power)
#       all_arrow_vals <- all_arrow_vals |> rbind(p$arrow)
#       all_pval_vals <- all_pval_vals |> rbind(p$pval)
#       all_coi <- all_coi |> rbind(p$coi)
#     }
#     all_result <- list(power = all_power_vals, arrow = all_arrow_vals, pval = all_pval_vals, coi = all_coi)
#     p <- with(all_result, {
#       if (which.image == "wp") {
#         break_vals <- c(0, round(lseq(min(power$Power) + 1, max(power$Power) + 1, 5) - 1, digits = 1))
#       } else {
#         break_vals <- round(seq(min(power$Power), max(power$Power), length.out = 6), digits = 1)
#       }
#       coi <- coi |>
#         mutate(y = if_else(y < min(power$Period), min(power$Period), y)) |>
#         mutate(y = if_else(y > max(power$Period), max(power$Period), y))

#       power |> ggplot(aes(x = date, y = Period, fill = Power)) +
#         facet_grid(series ~ .) +
#         geom_raster(width = 1, height = 1) +
#         scale_y_continuous(
#           trans = "log2", breaks = 2**(-2:4), expand = c(0, 0),
#           labels = scales::label_number(drop0trailing = TRUE)
#         ) +
#         scale_x_continuous(
#           expand = c(0, 0),
#           labels = ~ format(make_date(floor(.x), round(12 * (.x - floor(.x)) + 1), 1), "%b %Y")
#         ) +
#         scale_fill_distiller(
#           palette = "RdBu", trans = transf,
#           # breaks = break_vals,
#           limits = c(0, NA)
#         ) +
#         metR::geom_arrow(
#           aes(
#             x = date, y = Period, mag = 0.4, angle = 180 * Angle / pi
#           ),
#           inherit.aes = FALSE,
#           data = arrow, size = 0.4, arrow.length = 0.3
#         ) +
#         geom_contour(aes(x = date, y = Period, z = pval),
#           inherit.aes = FALSE, data = pval,
#           color = "white", linewidth = 0.5, fill = NA, breaks = 0.1
#         ) +
#         geom_polygon(aes(x = date, y = y), data = coi, inherit.aes = F, fill = "white", alpha = 0.5) +
#         theme_bw(16) +
#         labs(x = "Date", y = "Period") +
#         theme(legend.key.height = unit(4, "lines"))
#     })
#     return(p)
#   } else if (wc_list |> vec_depth() == 5) {
#     for (i in 1:length(wc_list)) {
#       group_name <- names(wc_list)[i]
#       wc_list2 <- wc_list[[i]]
#       for (wc in wc_list2) {
#         p <- with(wc, {
#           if (which.image == "wp") {
#             pval_vals <- Power.xy.pval |> as_tibble()
#             power_vals <- Power.xy |> as_tibble()
#           } else if (which.image == "wc") {
#             pval_vals <- Coherence.pval |> as_tibble()
#             power_vals <- Coherence |> as_tibble()
#           }
#           arrow_vals <- sAngle |> as_tibble()
#           colnames(power_vals) <- series |> pull("date")
#           colnames(arrow_vals) <- series |> pull("date")
#           colnames(pval_vals) <- series |> pull("date")

#           series_name <- series |>
#             select(-date) |>
#             select(-contains("trend")) |>
#             colnames()
#           series_name <- paste(series_name[1], series_name[2], sep = ":")
#           time_series <- series |> pull("date")


#           power_vals <- power_vals |>
#             mutate(Period = Period) |>
#             pivot_longer(-Period, names_to = "date", values_to = "Power") |>
#             mutate(date = as.Date(date)) |>
#             mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
#             mutate(series = series_name, group = group_name)



#           pval_vals <- pval_vals |>
#             mutate(Period = Period) |>
#             pivot_longer(-Period, names_to = "date", values_to = "pval") |>
#             mutate(date = as.Date(date)) |>
#             mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
#             mutate(series = series_name, group = group_name)

#           arrow_vals <- arrow_vals |>
#             mutate(Period = Period) |>
#             pivot_longer(-Period, names_to = "date", values_to = "Angle") |>
#             mutate(date = as.Date(date)) |>
#             mutate(date = year(date) + (month(date) - 1) / 12 + (day(date) - 1) / days_in_month(date)) |>
#             left_join(power_vals) |>
#             left_join(pval_vals) |>
#             group_by(date) |>
#             filter(row_number() %% 3 == 1) |>
#             ungroup() |>
#             group_by(Period) |>
#             filter(row_number() %% 6 == 1) |>
#             ungroup() |>
#             filter(Power > 0, pval <= 0.05) |>
#             mutate(series = series_name, group = group_name)

#           coi <- tibble(x = coi.1, y = 2**coi.2)

#           min_date <- min(time_series) |> as.Date()
#           min_date <- year(min_date) + (month(min_date) - 1) / 12 + (day(min_date) - 1) / days_in_month(min_date)

#           coi <- coi |>
#             mutate(date = min_date + x - 1) |>
#             mutate(series = series_name, group = group_name)


#           result <- list(power = power_vals, arrow = arrow_vals, pval = pval_vals, coi = coi)
#         })
#         all_power_vals <- all_power_vals |> rbind(p$power)
#         all_arrow_vals <- all_arrow_vals |> rbind(p$arrow)
#         all_pval_vals <- all_pval_vals |> rbind(p$pval)
#         all_coi <- all_coi |> rbind(p$coi)
#       }
#     }
#     all_result <- list(power = all_power_vals, arrow = all_arrow_vals, pval = all_pval_vals, coi = all_coi)

#   } else {
#     stop("wc_list is not list(wavelet) or list(list(wavelet))")
#   }

# p <- with(all_result, {
#   if (which.image == "wp") {
#     break_vals <- c(0, round(lseq(min(power$Power) + 1, max(power$Power) + 1, 5) - 1, digits = 1))
#   } else {
#     break_vals <- round(c(seq(min(power$Power), max(power$Power), length.out = 6)), digits = 1)
#     break_vals <- waiver()
#   }
#   coi <- coi |>
#     mutate(y = if_else(y < min(power$Period), min(power$Period), y)) |>
#     mutate(y = if_else(y > max(power$Period), max(power$Period), y))

#   power |> ggplot(aes(x = date, y = Period, fill = Power)) +
#     # facet_grid(interaction(group, series, sep = ": ") ~ .) +
#     facet_grid(group ~ series) +
#     geom_raster(width = 1, height = 1) +
#     scale_y_continuous(
#       trans = "log2", breaks = 2**(-2:4), expand = c(0, 0),
#       labels = scales::label_number(drop0trailing = TRUE)
#     ) +
#     scale_x_continuous(
#       expand = c(0, 0),
#       labels = ~ format(make_date(floor(.x), round(12 * (.x - floor(.x)) + 1), 1), "%b %Y")
#     ) +
#     # scale_fill_distiller(
#     #   palette = "RdBu", trans = transf,
#     #   # breaks = break_vals,
#     #   limits = c(0, NA)
#     # ) +
#     scale_fill_gradientn(
#       colours = rev(colorRampPalette(brewer.pal(9, "RdBu"))(256)),
#       values = scales::rescale(log(vals + 1)),
#       # trans = transf,
#       # breaks = scales::trans_breaks("log10", function(x) 10^x),
#       # labels = scales::trans_format("log10", scales::math_format(10^.x))
#       limits = c(0, NA)
#     ) +
#     metR::geom_arrow(
#       aes(
#         x = date, y = Period, mag = 0.4, angle = 180 * Angle / pi
#       ),
#       inherit.aes = FALSE,
#       data = arrow, size = 0.4, arrow.length = 0.3
#     ) +
#     geom_contour(aes(x = date, y = Period, z = pval),
#       inherit.aes = FALSE, data = pval,
#       color = "white", linewidth = 0.5, fill = NA, breaks = 0.1
#     ) +
#     geom_polygon(aes(x = date, y = y), data = coi, inherit.aes = F, fill = "white", alpha = 0.5) +
#     theme_bw(16) +
#     labs(x = "Date", y = "Period") +
#     theme(legend.key.height = unit(4, "lines"))
# })



# return(p)
# }


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
      labs(x = "Time", y = "Period") +
      theme(
        legend.key.height = unit(4, "lines"),
        text = element_text(family = "Times")
      )
  })
  vals <- quantile(all_result$power$Power, probs = seq(0, 1, length.out = 256)) |> as.numeric()
  if (wc_list |> vec_depth() == 3) {

  } else if (wc_list |> vec_depth() == 4) {
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
