#' rhmalaria: Relative Humidity and Malaria Transmission Modeling Tools
#'
#' Simulation, model-fitting, and diagnostic utilities for a mechanistic
#' *Plasmodium falciparum* transmission model relating malaria transmission
#' to relative humidity and temperature in Ahmedabad and Surat, India. See
#' the package vignettes for the full analysis pipeline.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import dplyr
#' @import ggplot2
#' @import pomp
#' @import tidyr
#' @importFrom doFuture %dofuture%
#' @importFrom foreach foreach %do%
#' @importFrom future plan multicore sequential
#' @importFrom ggtext element_markdown
#' @importFrom lubridate yday hour leap_year
#' @importFrom readr read_csv problems
#' @importFrom rlang sym
#' @importFrom tibble tibble is_tibble column_to_rownames
## usethis namespace: end
NULL
