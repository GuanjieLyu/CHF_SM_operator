# Szasz–Mirakyan cumulative hazard estimation

R code for *Smooth cumulative hazard function using Szász–Mirakyan operators for right-censored data*, by Guanjie Lyu and Cindy Feng.

This repository contains only R source files and this README. It implements the estimators, simulations, figures, and tables. The clinical examples read public data directly from the `KMsurv` and `survival` packages. No datasets or previously generated results are bundled.

## Getting started

Use **R 4.5.1**. From the repository folder, run:

```sh
Rscript scripts/setup.R
Rscript scripts/check.R
Rscript scripts/reproduce.R --profile smoke --stage all
```

Setup explicitly installs the pinned dependencies into a local `.library` folder. It needs Internet access and may require Rtools on Windows or a compatible compiler on other systems because an archived package can be installed from source. The analysis scripts do not install packages automatically.

The smoke run uses three simulated datasets per setting and reduced numerical resolution. It checks the complete workflow. Use its results only for testing.

For the full study:

```sh
Rscript scripts/reproduce.R --profile paper --stage all
```

This generates 500 datasets in each of 20 settings, with 2001 integration points and 1000 paired bootstrap resamples. Outputs are saved under `outputs/paper/`, including eight PDF figures, three LaTeX tables, replicate losses, and numerical summaries. The smoke profile uses `outputs/smoke/`.

To generate only the introductory illustration and clinical figures:

```sh
Rscript scripts/reproduce.R --stage applications
```

Other stages are `simulate`, `summarize`, and `figures`. Use `--output PATH` to choose a different destination, or `--help` for all options. Completed simulation cells can be resumed with the same command. Use a new output directory after changing code, package versions, or the study configuration. Separate processes must use separate output directories.

## Files

| File | Purpose |
| --- | --- |
| `R/bootstrap.R` | Load the modules and local dependency library |
| `R/validation.R` | Input checks and reproducible random-number handling |
| `R/estimators.R` | Grouped Nelson–Aalen, Szasz–Mirakyan, and integrated kernel estimators |
| `R/design.R` | Event-time models, censoring calibration, and fixed seeds |
| `R/simulation.R` | Replicate losses, Monte Carlo standard errors, and paired bootstrap summaries |
| `R/pipeline.R` | Study settings, resumable execution, output checks, and summaries |
| `R/plot_simulations.R` | Simulation figures and tables |
| `R/plot_applications.R` | Introductory illustration and clinical figures |
| `scripts/setup.R` | Install the pinned dependencies |
| `scripts/reproduce.R` | Run the selected analysis stages |
| `scripts/check.R` | Run the mathematical and design checks |
| `scripts/example.R` | Minimal estimator example |
| `tests/test_estimators.R` | Mathematical identities, ties, boundaries, and invalid inputs |
| `tests/test_design.R` | Censoring calibration, summaries, and random-number checks |

## Use the estimators

```r
source("R/bootstrap.R")
load_project(".")

fit <- grouped_na(
  x = c(0.5, 1, 1, 2, 3),
  delta = c(1, 1, 0, 1, 0)
)
time <- seq(0, 3, length.out = 301)
na <- na_curve(time, fit)
sm <- sm_curve(time, fit, m = 2.5)
```

`delta = 1` denotes an observed event. Times must be finite and nonnegative. Both smooth estimators require positive observed failure times. The smoothing parameter can be any finite real `m > 0`. The grouped estimator retains subjects censored at an event time in that time's risk set.

## Statistical protocol

The main simulations use exponential, Weibull, and log-logistic event times, sample sizes 100 and 500, and independent censoring calibrated to 10%, 30%, or 50%. Two additional settings use a bounded piecewise-constant hazard with 30% censoring. Every generated sample is retained.

All estimators are compared on the same samples over `[0, 3.7]`. Smoothing is evaluated on `m = 1, 2, 4, ..., 100`. The kernel comparator uses local Epanechnikov smoothing with both boundary corrections, integration on its native grid, and constant cumulative-hazard continuation after the fitted endpoint. Its native sequential tie weights differ from the grouped increments used by Nelson–Aalen and Szasz–Mirakyan.

The minimum over the smoothing grid is an **oracle benchmark** based on the known simulation truth. It is not a data-driven selector. The paired bootstrap reselects that minimum and describes Monte Carlo variability. The theoretical claims concern cumulative hazards and do not establish a hazard-derivative consistency theorem.

## Dependencies and public data

The pinned versions are `muhaz` 1.2.6.4, `survival` 3.8-3, `Matrix` 1.7-3, `lattice` 0.22-7, and `KMsurv` 0.1-6. Their versions are specified in the R source, so no separate lock file is required in this repository.

- **Laryngeal cancer:** `KMsurv::larynx`, 90 observations, with 50 deaths and 40 censored observations. Follow-up values are interpreted in years, following Kardaun's primary report. The package documentation's label of months conflicts with that report. See the [package](https://cran.r-project.org/package=KMsurv) and [primary report](https://www.vvsor.nl/wp-content/uploads/2020/06/KM1983012001.pdf).
- **Primary biliary cholangitis:** `survival::pbc`, all 418 observations. Time is measured in days. The endpoint is death or transplantation, using `status != 0`, giving 186 composite events. See the [dataset documentation](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/pbc.html).

These data remain under their upstream package terms, GPL-3.0-or-later for `KMsurv` and LGPL-2.0-or-later for `survival`. They are loaded from installed packages and are not redistributed in this source folder. The MIT license below covers the original code and this README.

## Upload to GitHub

Upload this README together with the `R`, `scripts`, and `tests` folders. After running the code locally, leave the generated `outputs` and `.library` folders out of your upload.

## MIT License

Copyright (c) 2026 Guanjie Lyu, Cindy Feng, and contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
