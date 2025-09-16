# causalGBKMR

> Causal Inference for Time-Varying Environmental Mixtures using g-BKMR

## Overview

`causalGBKMR` implements g-BKMR (g-computation with Bayesian Kernel Machine Regression) for causal inference with time-varying environmental mixtures and time-varying confounders. This package addresses the challenge of estimating causal effects when multiple correlated exposures and confounders change over time.

## Key Features

- **Time-varying exposure analysis**: Handle multiple correlated exposures across time points
- **Time-varying confounding adjustment**: Properly account for confounders that change over time
- **Flexible variable detection**: Automatically detect variable patterns in your data
- **User-friendly interface**: Simple API with sensible defaults
- **Comprehensive validation**: Built-in data validation and error checking
- **Rich output**: Detailed results with confidence intervals and diagnostics

## Installation

### From GitHub (Development Version)

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install causalGBKMR
devtools::install_github("Lyric98/causalGBKMR")
```

### Dependencies

The package requires the following R packages:

```r
install.packages(c("bkmr", "fields", "dplyr", "ggplot2", "parallel"))

# Note: bkmr package may need to be installed from GitHub
devtools::install_github("jenfb/bkmr")
```

## Quick Start

### Basic Example

```r
library(causalGBKMR)

# Step 1: Prepare your data
# Your raw data: Y (outcome), Z (exposures), X (covariates)
n <- 500
Y <- rnorm(n)  # Outcome
Z <- matrix(rlnorm(n * 6), nrow = n, ncol = 6)  # 2 metals × 3 time points
X <- matrix(rnorm(n * 5), nrow = n, ncol = 5)   # 1 TD covariate × 3 time + 2 baseline

# Step 2: Convert to g-BKMR format
prepared_data <- prepare_gbkmr_data(
  Y = Y,
  Z = Z, 
  X = X,
  time_points = 3,
  mixture_components = 2,
  td_covariates = 1,
  baseline_covariates = 2,
  td_covariate_names = "bmi"
)

# Step 3: Run g-BKMR analysis
results <- gbkmr_run(
  data = prepared_data,
  time_points = 3
)

# Step 4: View results
print(results)
summary(results)
```

## Data Requirements

### Input Data Format

Your data should be organized as:

1. **Y**: Outcome vector (length n)
2. **Z**: Exposure matrix (n × [number_of_exposures × time_points])
   - Organized as: [Exp1_T0, Exp2_T0, ..., Exp1_T1, Exp2_T1, ...]
3. **X**: Covariate matrix (n × [td_covariates × time_points + baseline_covariates])
   - Time-dependent covariates for T-1 time points, followed by baseline covariates

### Example Data Structure

For 3 time points, 2 exposures, 1 time-dependent covariate, 2 baseline covariates:

- **Z matrix**: 6 columns [Metal1_T0, Metal2_T0, Metal1_T1, Metal2_T1, Metal1_T2, Metal2_T2]
- **X matrix**: 5 columns [BMI_T1, BMI_T2, Sex, Age]

## Main Functions

### Core Functions

- `gbkmr_run()`: Main user interface for g-BKMR analysis
- `prepare_gbkmr_data()`: Convert raw matrices to g-BKMR format
- `detect_variable_patterns()`: Automatically detect variable structure

### Utility Functions

- `validate_user_matrices()`: Validate input data dimensions
- `run_gbkmr_panel()`: Core analysis engine (internal)

## Advanced Usage

### Custom Variable Names

```r
prepared_data <- prepare_gbkmr_data(
  Y = Y, Z = Z, X = X,
  time_points = 3,
  mixture_components = 2,
  td_covariates = 2,
  baseline_covariates = 1,
  td_covariate_names = c("bmi", "blood_pressure")
)
```

### Analysis Parameters

```r
results <- gbkmr_run(
  data = prepared_data,
  time_points = 3,
  iter = 20000,           # MCMC iterations
  n = 500,                # Sample size for analysis
  use_knots = TRUE,       # Use kernel approximation
  n_knots = 50,          # Number of knots
  verbose = TRUE          # Print progress
)
```

## Interpretation

The main result is the **causal effect estimate**: the expected change in outcome when all exposures change from their 25th percentile to their 75th percentile across all time points.

### Results Structure

```r
# Main causal effect
results$causal_effect$estimate
results$causal_effect$lower    # 95% CI lower bound
results$causal_effect$upper    # 95% CI upper bound

# Counterfactual means
results$counterfactual_means$low   # Mean under low exposure
results$counterfactual_means$high  # Mean under high exposure

# Variable importance
results$variable_importance

# Detection information
results$detection_info
```

## Methodological Background

g-BKMR combines:

1. **g-computation**: A causal inference method for time-varying treatments
2. **BKMR**: Flexible modeling of exposure mixtures using kernel machine regression
3. **Time-varying confounding adjustment**: Proper handling of confounders that change over time


## Troubleshooting

### Common Issues

1. **"number of design points >= the number of candidates"**
   - Solution: Use `use_knots = FALSE` for small datasets (n < 100)

2. **Dimension errors in data preparation**
   - Check that your Z and X matrices have the correct number of columns
   - Use `validate_user_matrices()` to verify dimensions

3. **Variable detection failures**
   - Ensure your prepared data follows the expected naming conventions
   - Check column names match the expected patterns


## Contact

- **Maintainer**: Yanran Li
- **Email**: yl5465@cumc.columbia.edu