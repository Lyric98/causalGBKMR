# causalBKMR: Causal Inference for Time-Varying Environmental Mixtures via g-BKMR

## About the Method

`causalBKMR` implements **g-BKMR**, a method for estimating the causal
effect of time-varying environmental mixtures on a health outcome.
The method combines two ideas:

- **The g-formula** (Robins, 1986): a standard tool from causal inference for
  handling time-varying confounding---the setting where intermediate
  covariates (e.g., body mass index, blood pressure) are affected by past
  exposures and simultaneously confound future exposure–outcome relationships.
- **Bayesian Kernel Machine Regression (BKMR)** (Bobb et al., 2015): a
  flexible Bayesian nonparametric model that captures nonlinear and
  non-additive effects of correlated environmental mixtures while
  simultaneously performing variable selection.

Together, g-BKMR provides valid causal estimates of mixture effects when
exposures and confounders evolve over multiple visits, while retaining
BKMR's ability to flexibly model complex exposure–response surfaces.

The estimand is the **average causal effect (ACE)** of a contrast between
two exposure trajectories $\bar{a}$ and $\bar{a}^*$ (by default, all mixture
components at the 25th vs. 75th percentile):

$$
\text{ACE} = E[Y^{\bar{a}^*}] - E[Y^{\bar{a}}]
$$

A full description of the methodology is available in the
[Methodology Overview](https://lyric98.github.io/causalBKMR/articles/gbkmr_method_overview.html).

## What this package provides

- A single-function interface (`gbkmr_run()`) that handles the full
  g-computation pipeline: automated data preparation, sequential BKMR
  fitting, Monte Carlo g-computation, and posterior summarisation.
- Support for **arbitrary numbers of time points, mixture components, and
  time-varying confounders**, with user-configurable intervention contrasts.
- A **fast-BKMR engine** (via the `fbkmr` package) that is automatically
  activated for large samples (n > 2,000), extending applicability to
  cohorts of MESA / ARIC scale.
- **Continuous and binary outcomes** (probit BKMR for the latter).
- An **Input Audit** block that prints the package's interpretation of
  your data before launching MCMC, plus **convergence diagnostics**
  (effective sample size, Geweke z-scores) for every fitted model.

We welcome your feedback (email <yl5465@cumc.columbia.edu>).

## Installation

Install required CRAN dependencies:

``` r
install.packages(c("bkmr", "dplyr", "fields"))
```

Install `causalBKMR` from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("Lyric98/causalBKMR")

library(causalBKMR)
```

Optional: for large-sample analyses, also install `fbkmr`:

``` r
remotes::install_github("junwei-lu/fbkmr")
```

## Quick Start

```r
library(causalBKMR)

# 1. Prepare your data (Y outcome, Z mixture, X covariates)
dat <- prepare_gbkmr_data(
  Y = outcome, Z = metals_matrix, X = covariates_matrix,
  time_points = 3, mixture_components = 2,
  td_covariates = 1, baseline_covariates = 1,
  td_covariate_names = "waist"
)

# 2. Run the analysis (engine is auto-selected based on sample size)
results <- gbkmr_run(
  data = dat, time_points = 3,
  iter = 15000, K = 1000
)

# 3. View results
print(results)    # ACE, 95% credible interval, counterfactual means
summary(results)  # + settings, diagnostics
```

See the [Quick Start guide](https://lyric98.github.io/causalBKMR/articles/gbkmr_quickstart.html)
for a complete walkthrough.

## Documentation

| Topic | Link |
|-------|------|
| Methodology overview | [Method](https://lyric98.github.io/causalBKMR/articles/gbkmr_method_overview.html) |
| Quick start | [Quick Start](https://lyric98.github.io/causalBKMR/articles/gbkmr_quickstart.html) |
| Custom contrasts | [Example 1](https://lyric98.github.io/causalBKMR/articles/gbkmr_example_contrasts.html) |
| Binary outcome | [Example 2](https://lyric98.github.io/causalBKMR/articles/gbkmr_example_binary.html) |
| Large-sample fastBKMR | [Example 3](https://lyric98.github.io/causalBKMR/articles/gbkmr_example_fastbkmr.html) |
| Diagnostics & troubleshooting | [Diagnostics](https://lyric98.github.io/causalBKMR/articles/gbkmr_diagnostics.html) |
| Function reference | [Functions](https://lyric98.github.io/causalBKMR/reference/) |

## Citation

If you use `causalBKMR` in your research, please cite the method paper
(Chai et al., g-BKMR) and the software paper (forthcoming).

## Related Packages

`causalBKMR` is part of a broader effort to make causal inference
methods for environmental mixtures accessible:

- [`bkmr`](https://github.com/jenfb/bkmr) — single-time-point BKMR (Bobb et al.)
- [`fbkmr`](https://github.com/junwei-lu/fbkmr) — scalable BKMR via divide-and-conquer
- [`causalmixtures`](https://arcedr.github.io/causalmixtures/) — BKMR-CMA and BKMR-MI
- `causalBKMR` *(this package)* — g-formula + BKMR for time-varying mixtures
