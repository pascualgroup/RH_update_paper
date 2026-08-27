# rhmalaria

Research compendium for a mechanistic *Plasmodium falciparum* transmission
model relating malaria transmission to relative humidity and maximum
temperature in Ahmedabad and Surat, India. The transmission model is a
`pomp`-based compartmental model (susceptible/exposed/infectious stages,
climate-driven transmission rate) fit to monthly case-count and climate data
via iterated filtering and profile likelihood.

This repository is organized as a "semi-package": the reusable modeling code
lives in `R/` and is documented with roxygen2, while the analysis pipeline
itself lives in `vignettes/` as Quarto documents that call into `R/`.

## Installation (local development)

```r
# CRAN + GitHub dependencies (reads them straight from DESCRIPTION):
install.packages("remotes")
remotes::install_deps(dependencies = TRUE)

# Then, from the repository root:
devtools::load_all()
```

The vignettes themselves locate the repository root with `rprojroot`
(`rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))`) and call
`pkgload::load_all()` on it in their first chunk, so they don't need the
package pre-installed — running `quarto render vignettes/<name>.qmd` (or
`altdoc::render_docs()`) from anywhere inside the repository is enough.

## Cluster usage (SLURM)

The fitting/profiling pipeline is designed to run as SLURM array jobs from
`scripts/`. Unlike the vignettes, the job scripts (`scripts/run_fitting.R`,
`scripts/run_profile.R`) expect the package to be properly **installed**
(`library(rhmalaria)`), not `load_all()`-ed, since each array task starts a
fresh `R CMD BATCH` process.

1. **Install once per cluster environment** (e.g. in your `$HOME`, before
   syncing to scratch — see the sync logic in `scripts/run_fitting.sh`):

   ```bash
   module load r/4.5.1   # or whatever R module the cluster provides
   Rscript -e 'install.packages("remotes"); remotes::install_deps(dependencies = TRUE)'
   Rscript -e 'remotes::install_local(".", dependencies = FALSE, force = TRUE)'
   ```

   (`scripts/install_packages.R` is an older, hand-maintained equivalent of
   the `remotes::install_deps()` step above, kept for reference; DESCRIPTION
   is now the source of truth for dependencies.)

2. **Submit a fitting run**: `scripts/submit_run_and_combine.sh` submits the
   `run_fitting.sh` array job and, with a dependency, the
   `combine_array_outputs.sh` aggregator that runs once the array completes:

   ```bash
   ./scripts/submit_run_and_combine.sh --array 1-100
   ```

   `run_fitting.sh` itself syncs the repo to `$SCRATCH`, loads the `r`
   module, and runs `R CMD BATCH --vanilla --args <city> <model> <covariate>
   <window_start> <window_end> <VC_file> scripts/run_fitting.R` per array
   task — edit the hard-coded `--args` line (and `#SBATCH` resources) for
   the city/model/window you want to fit.

3. **Submit a profile-likelihood sweep**: `scripts/profile_loop.sh` loops
   over every model parameter and calls `scripts/submit_profile.sh` (which
   runs `scripts/run_profile.R`) for each one:

   ```bash
   RUN_ID=my_run bash scripts/profile_loop.sh
   ```

Both `run_fitting.R` and `run_profile.R` read their remaining configuration
from SLURM environment variables (`SLURM_ARRAY_TASK_ID`, `RUN_ID`,
`END_YEAR`, ...) and command-line arguments; see the comment block at the
top of each script for the full list.

## Data organization

Case/climate data lives under `data/<City>/`, one CSV per city (e.g.
`data/Ahmedabad/dataset_Ahmedabad_2023.csv`). This is the file you pass as
`fit_model(dataset = ...)`. It's a monthly time series, one row per month,
with (at least) these columns:

| column | meaning |
| --- | --- |
| `time` | fractional-year time index (`year + (month - 1) / 12`); the `pomp` time axis |
| `year`, `month` | calendar year and month |
| `pop` | population size that month |
| `dpopdt` | population growth rate (re-derived/smoothed internally; can be `0` in the raw file) |
| `PF` | observed *P. falciparum* case count — the value the model is fit to |
| `season1` ... `season6` | precomputed seasonal basis: the fraction of that month falling in each of 6 two-month blocks of the calendar year (the regressors for the seasonal transmission-rate coefficients `b1`...`b6`) |
| climate covariate columns | one column per candidate climate driver, e.g. `hadISD` (relative humidity), `Max_Temp`, `Min_Temp`, `Mean_Temp` |

`fit_model()`'s `covariate` argument names one of those climate columns;
`window_start`/`window_end` pick which calendar months of that covariate get
averaged into a single yearly value (e.g. `window_start = 8, window_end = 9`
averages August-September), which is then min-max normalized over the
fitting period and used as the `covariate` term in the transmission-rate
equation below.

## Implemented models

Three model variants live in `dynamical_model/fitting/pomp_object_<model>.R`
and are selected with `fit_model(model = ...)`. All three share the same
compartmental structure (S1 -> E -> I1 -> S2 -> I2 -> S1, i.e. susceptible,
exposed, first infection, partially-immune susceptible, reinfection, with
waning immunity); they differ only in the transmission-rate term `betaIN`:

- **`baseline`** — seasonal transmission only, no climate covariate:
  `betaIN = exp(b1*season1 + ... + b6*season6)`.
- **`inf_exponent`** — the main climate-driven model: adds a covariate term
  on top of the August-September (`season4`) seasonal block:
  `betaIN = exp(b1*season1 + ... + b6*season6 + season4*bH*covariate)`.
- **`inf_exponent_scaled`** — `inf_exponent` plus one extra multiplicative
  parameter `scaled_b`: `betaIN = scaled_b * exp(...)`. Used by the
  sliding-window sensitivity analysis (`vignettes/refitting_window.qmd`,
  `refitting_window_plot.qmd`, `results_reffiting.qmd`) to ask how much the
  already-fitted seasonal/climate shape would need to be rescaled to match
  observed cases in a given out-of-sample year.

## Fitting a model with `fit_model()`

`fit_model()` is the main entry point: it preprocesses `dataset`, builds the
`pomp` object for `model`, and runs `pomp::mif2()` (iterated filtering) from
one or more starting parameter vectors.

```r
library(rhmalaria)

# `parameters` must be a data.frame/tibble, one row per starting point —
# fit_model()/mif2() multi-start from every row (in parallel across
# `n_cores` when allow_parallel = TRUE). A bare named vector is *not*
# accepted here even though some lower-level helpers allow one.
start <- data.frame(
  sigOBS = 0.138, sigPRO = 0.14, muS2S1 = 2.06, muEI1 = 23.7, muI1S2 = 10.6,
  muI2S2 = 1.99, betaOUT = 5.77e-05, delta = 0.02, rho = 0.012, tau = 0.00619,
  q0 = 0.987, alpha = 0.502,
  b1 = -2.91, b2 = -1.15, b3 = -3.51, b4 = 1.09, b5 = -0.756, b6 = -0.0211,
  bH = 1.51,
  S1_0 = 0.302, E_0 = 0.000621, I1_0 = 0.00492, S2_0 = 0.606, I2_0 = 0.0868,
  K_0 = 0.0499, F_0 = 0.715
)

fit <- fit_model(
  parameters = start,
  city = "Ahmedabad",
  model = "inf_exponent",
  covariate = "hadISD",
  start_year = 1997, start_month = 1,
  end_year = 2011, end_month = 12,
  window_start = 8, window_end = 9,
  dataset = "data/Ahmedabad/dataset_Ahmedabad_2023.csv",
  mode = "refined_second", # `start` above is already a good fit; polish it
  output_to_file = FALSE
)

fit$result # one row per starting point: fitted params + loglik + flag (1 = converged)
```

`mode` controls how far `mif2` is allowed to move the parameters per
iteration: `"fit"` (large random-walk steps, for a rough/unfitted starting
guess) -> `"refined"` -> `"refined_second"` (tiny steps, for polishing an
already-good fit, as in the example above). The standard pipeline (see
`scripts/run_profile.R` / `scripts/run_fitting.R`) starts an unfitted grid at
`"fit"` and chains all three stages, each one re-starting `mif2` from the
best rows (`arrange(desc(loglik))`) of the previous stage's `fit$result`.
Iterated filtering is stochastic, so with a real (not artificially small)
`Np`/`Nmif` some starting rows can still fail to converge — that's what the
`flag` column and the `filter(flag == 1)` step are for; production runs
multi-start from many rows precisely because a few are expected to fail.

## Profiling one parameter

Profile likelihood for a parameter (say `b6`) reuses the exact same
`fit_model()` call, with two additions:

- `parameters` becomes a *grid*: a set of starting rows that all fix `b6` at
  the same value but vary every other parameter (one CSV per grid point,
  e.g. `param_grids/Surat/Max_Temp/param_grid_b6.csv`).
  `vignettes/create_profile_grids.qmd` builds these grids from an initial
  unconstrained fit.
- `profile_var = "b6"` tells `generate_rw_sd()` to give `b6` a random-walk
  standard deviation of `0`, so `mif2` never moves it while every other
  parameter is still refined:

```r
fit_model(
  parameters = read.csv("param_grids/Surat/Max_Temp/param_grid_b6.csv"),
  city = "Surat", model = "inf_exponent", covariate = "Max_Temp",
  window_start = 5, window_end = 9,
  dataset = "data/Surat/dataset_Surat_2023.csv",
  profile_var = "b6",
  mode = "fit", output_to_file = FALSE
)
```

Running this (through the `"fit"`/`"refined"`/`"refined_second"` stages)
across every grid point for `b6` and collecting the log-likelihoods traces
out the profile. `clean_profile()` (`R/cleaning_profile_results.R`) combines
the per-parameter result CSVs into one table, and `mcap()`/`mcap_checked()`
(`R/mcap.R`) turn the noisy Monte Carlo profile into a smoothed profile and
confidence interval — see `vignettes/parameter_profiles.qmd`. On the
cluster, this whole sweep (one array job per parameter) is what
`scripts/profile_loop.sh` / `scripts/submit_profile.sh` automate.

## Package contents (`R/`)

- **simulation_functions.R** — preprocess climate/case data, build and fit
  the `pomp` transmission model, simulate/forecast from fitted parameters
  (training-period simulation plus out-of-sample forecasting), and produce
  the associated diagnostic plots.
- **mcap.R** — the Monte Carlo Adjusted Profile (MCAP) method for turning a
  noisy profile log-likelihood into a smoothed profile and confidence
  interval, plus supporting binned-diagnostic plots.
- **cleaning_profile_results.R** — combine and clean per-parameter
  profile-likelihood result CSVs produced by the fitting pipeline.
- **plot_wavelet.R** — wavelet power/coherence diagnostic plots from
  `WaveletComp`-style wavelet-analysis objects.
- **run_manifest.R** / **save_safe.R** — atomic JSON/CSV write helpers and
  run-manifest bookkeeping used by the cluster (SLURM) fitting jobs in
  `scripts/`.

## Analysis pipeline (`vignettes/`)

- **create_profile_grids.qmd** — build the parameter grids used to launch
  the profile-likelihood fitting jobs.
- **parameter_profiles.qmd** — profile-likelihood curves and MCAP confidence
  intervals for the fitted model parameters.
- **quantitative_estimates_fits.qmd** — fit each city/covariate model,
  simulate from the fitted posterior, and score against observed cases with
  CRPS over training and out-of-sample (pre-/post-2019) periods.
- **crps_code.qmd** — CRPS scoring exploration used for the sliding-window
  refitting analysis (window-length sensitivity).
- **refitting_window.qmd** / **refitting_window_plot.qmd** — repeated
  re-fitting over a sliding window of training-period lengths.
- **results_reffiting.qmd** — results of the yearly re-fitting/prediction
  analysis.
- **plausibility_of_rates.qmd** — biologically-plausibility checks on the
  fitted model: susceptible-fraction trajectory and seasonal transmission
  rate over the full study period.
- **new_figures.qmd** — main manuscript figures, including model
  predictions, covariate correlation/wavelet-coherence, and covariate
  scaling-window diagnostics.
- **climate_correlations.qmd** — correlation between transmission-season
  (Aug-Nov) case counts and every possible start/end-month covariate window,
  for relative humidity and temperature in both cities; visualized as a
  start-vs-end-month heatmap (the basis for picking `window_start`/
  `window_end` elsewhere).
- **malaria_suitability.qmd** — a temperature-driven vectorial-capacity
  suitability index (from mosquito trait curves) for comparison stations in
  the endemic states of Odisha, Jharkhand, and Chhattisgarh, contrasting a
  baseline period against later years via permutation-style histograms.

## Documentation site

Function reference (from roxygen2) and the vignettes above are built into a
static site with [altdoc](https://altdoc.etiennebacher.com/), published at
<https://pascualgroup.github.io/RH_update_paper/> (served from `docs/` on
`main` via GitHub Pages, so `docs/` is committed rather than gitignored).

```r
altdoc::render_docs(freeze = TRUE)
```

Always pass `freeze = TRUE`: it skips re-knitting any vignette or man page
that hasn't changed since the last render (hashes are cached in
`altdoc/freeze.rds`, itself gitignored/machine-local), so adding one new
vignette only executes that vignette instead of re-running the whole
pipeline — several of these vignettes are expensive (real `pomp` fitting)
or memory-heavy enough to be impractical to rebuild on every change. Without
`freeze = TRUE` (the plain `altdoc::render_docs()` default), everything is
re-rendered from scratch every time.

If you change a function in `R/` that a vignette depends on, the freeze
cache won't know to invalidate that vignette (it only hashes the vignette's
own source) — delete `altdoc/freeze.rds`, or render that vignette on its own
with `quarto render vignettes/<name>.qmd --execute`, to force a refresh.
Remember to run `devtools::document()` first if you changed any roxygen
comments, so the reference pages pick up the change.

After rendering, push `docs/` along with your other changes to update the
published site.
