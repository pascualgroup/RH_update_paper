# Installation helper for all R packages referenced in this repository.
# Linux (Debian/Ubuntu) system packages typically needed to build these from source:
#   sudo apt-get update && sudo apt-get install -y \
#     build-essential libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
#     libfreetype6-dev libharfbuzz-dev libfribidi-dev libjpeg-dev libpng-dev \
#     libtiff5-dev libicu-dev
# If the 'circumstance' package is not on CRAN for you, set the correct GitHub
# repo in `github_pkgs` below (example uses kingaa/circumstance).

cran_pkgs <- c(
  "jsonlite",
  "readr",
  "dplyr",
  "pomp",
  "tidyverse",
  "lubridate",
  "doParallel",
  "doRNG",
  "ggtext",
  "doFuture",
  "tsibble",
  "feasts",
  "here",
  "ggplot2",
  "rstatix",
  "forecast",
  "caret",
  "future",
  "tibble",
  "rlang",
  "tidyr",
  "ggdist",
  "purrr"
)

# Name -> repo path format expected by remotes::install_github()
github_pkgs <- c(
  circumstance = "kingaa/circumstance"
)

cran_pkgs <- unique(cran_pkgs)

install_if_missing <- function(pkgs, install_fn) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install_fn(missing)
}

install_if_missing(
  cran_pkgs,
  function(missing) install.packages(missing, repos = "https://cloud.r-project.org")
)

if (length(github_pkgs)) {
  install_if_missing(
    names(github_pkgs),
    function(missing) {
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes", repos = "https://cloud.r-project.org")
      }
      remotes::install_github(unname(github_pkgs[missing]))
    }
  )
}
