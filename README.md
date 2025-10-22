# causalGBKMR: Causal Inference for Time-Varying Environmental Mixtures

> **g-BKMR**: Combining g-formula with Bayesian Kernel Machine Regression for causal inference with time-varying correlated environmental mixtures

## Overview

`causalGBKMR` implements a novel causal inference method that combines the g-formula for handling time-varying confounding with Bayesian Kernel Machine Regression (BKMR) for flexible modeling of environmental mixture effects.

**Key Features:**
- ✅ Estimate causal effects of **multiple, correlated, time-varying exposures**
- ✅ Account for **time-varying confounding** using g-computation
- ✅ Model **non-linear and non-additive** exposure-response relationships
- ✅ Perform **variable selection** in high-dimensional settings
- ✅ Handle **multiple time-dependent confounders** simultaneously

**Ideal for:**
- Environmental health studies with longitudinal exposure data
- Assessing health effects of pollutant mixtures over time
- Studies where confounders are affected by past exposures

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Data Structure](#data-structure)
- [Detailed Examples](#detailed-examples)
- [Function Reference](#function-reference)
- [Methodological Details](#methodological-details)
- [Best Practices](#best-practices)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

---

## Installation

### Dependencies

```r
# Install required packages
install.packages("bkmr")
install.packages("dplyr")
install.packages("fields")
```

### Install causalGBKMR

```r
# From GitHub (when available)
# devtools::install_github("username/causalGBKMR")

# For development version, source files directly:
source("R/01-validation.R")
source("R/02-data-preparation.R")
source("R/03-variable-detection.R")
source("R/04-core-analysis.R")
```

---

## Quick Start

### Minimal Working Example

```r
library(causalGBKMR)

# 1. Simulate example data
set.seed(123)
n <- 200
Y <- rnorm(n)
Z <- matrix(rlnorm(n * 6), nrow = n, ncol = 6)  # 2 exposures x 3 time points
X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)   # 1 TD covariate x 2 times + 1 baseline

# 2. Prepare data
prepared_data <- prepare_gbkmr_data(
  Y = Y, Z = Z, X = X,
  time_points = 3,
  mixture_components = 2,
  td_covariates = 1,
  baseline_covariates = 1,
  td_covariate_names = "waist"
)

# 3. Run analysis
results <- run_gbkmr_panel(
  sim_popn = prepared_data,
  T = 3, p = 2,
  mediator_basenames = c("waist"),
  common_covariates = c("sex", "waist0"),
  iter = 10000,  # Use more iterations in practice (e.g., 24000)
  sel = seq(8000, 10000, by = 50),
  K = 500
)

# 4. View results
cat("ATE:", results$diff_gBKMR, "\n")
ci <- quantile(rowMeans(results$Yastar_mat) - rowMeans(results$Ya_mat), c(0.025, 0.975))
cat("95% CI: [", ci[1], ",", ci[2], "]\n")
```

**Expected output:**
```
ATE: 0.127
95% CI: [ -0.089 , 0.343 ]
```

---

## Data Structure

### Input Formats

The package supports two input formats:

#### Format 1: List Structure (Recommended)

More intuitive and self-documenting:

```r
# Time-dependent covariates organized by variable name
TD_covariates <- list(
  parity = list(t1 = C_parity_t1, t2 = C_parity_t2, t3 = C_parity_t3),
  breastfeeding = list(t1 = C_bf_t1, t2 = C_bf_t2, t3 = C_bf_t3),
  BMI = list(t1 = C_bmi_t1, t2 = C_bmi_t2, t3 = C_bmi_t3)
)

# Exposures organized by time point (each element is n x p matrix)
exposures <- list(
  t0 = X_exp_t0,  # Baseline exposures
  t1 = X_exp_t1,  # Time 1 exposures
  t2 = X_exp_t2,
  t3 = X_exp_t3
)

# Baseline confounders (single matrix: n x k)
C0 <- C_baseline

# Outcome (vector: length n)
Y <- outcome
```

**To use this format**, convert to matrices using:

```r
convert_lists_to_matrices <- function(exposures, TD_covariates, C0, Y) {
  # Flatten exposures
  Z <- do.call(cbind, exposures)
  
  # Flatten TD covariates (time-first order)
  n_times <- length(TD_covariates[[1]])
  X_td_list <- list()
  for (t_idx in 1:n_times) {
    for (cov_name in names(TD_covariates)) {
      X_td_list[[length(X_td_list) + 1]] <- TD_covariates[[cov_name]][[t_idx]]
    }
  }
  X_td <- do.call(cbind, X_td_list)
  X <- cbind(X_td, C0)
  
  list(Y = Y, Z = Z, X = X)
}

# Convert and use
data_matrices <- convert_lists_to_matrices(exposures, TD_covariates, C0, Y)
prepared_data <- prepare_gbkmr_data(data_matrices$Y, data_matrices$Z, data_matrices$X, ...)
```

#### Format 2: Matrix Structure (Direct)

Traditional wide-format matrices:

```r
# Matrix Z: Exposures (n x Adim*T)
# Column order: [M1_T0, M2_T0, ..., Mp_T0, M1_T1, M2_T1, ..., Mp_T(T-1)]
Z <- cbind(As_T0, Cd_T0, As_T1, Cd_T1, As_T2, Cd_T2)

# Matrix X: Covariates (n x Ldim*T + n_baseline)
# Column order: [TD1_T1, TD2_T1, ..., TDLdim_T(T-1), Baseline1, Baseline2, ...]
X <- cbind(waist_T1, waist_T2, sex, age)

# Vector Y: Outcome (length n)
Y <- blood_glucose
```

**Important Notes:**
- Exposures (Z) include **baseline (T0)**: indices 0 to T-1
- Time-dependent covariates (X) start at **T1**: indices 1 to T-1
- This reflects that TD covariates are measured *after* baseline exposure

### Output Format

After `prepare_gbkmr_data()`, variables are ordered as:

```
[sex, baseline_2, ..., waist_0, logM1_0, logM2_0, logM1_1, logM2_1, waist_1, ..., Y, id]
```

Example for T=3, p=2, Ldim=1:
```
[sex, waist_0, logM1_0, logM2_0, logM1_1, logM2_1, waist_1, logM1_2, logM2_2, waist_2, Y, id]
```

---

## Detailed Examples

### Example 1: Environmental Health Study (List Format)

Study of metal mixtures and cardiovascular outcomes with multiple time-dependent confounders.

```r
# Study design:
# - 3 metals (As, Cd, Pb) at 4 time points
# - 2 time-dependent confounders (waist, BMI)
# - Outcome: blood glucose
# - Baseline: sex, age

# Step 1: Organize data into lists
TD_covariates <- list(
  waist = list(t1 = waist_v1, t2 = waist_v2, t3 = waist_v3),
  bmi = list(t1 = bmi_v1, t2 = bmi_v2, t3 = bmi_v3)
)

exposures <- list(
  t0 = cbind(As_baseline, Cd_baseline, Pb_baseline),
  t1 = cbind(As_v1, Cd_v1, Pb_v1),
  t2 = cbind(As_v2, Cd_v2, Pb_v2),
  t3 = cbind(As_v3, Cd_v3, Pb_v3)
)

C0 <- cbind(sex, age)
Y <- blood_glucose

# Step 2: Convert to matrix format
data_matrices <- convert_lists_to_matrices(exposures, TD_covariates, C0, Y)

# Step 3: Prepare for g-BKMR
prepared_data <- prepare_gbkmr_data(
  Y = data_matrices$Y,
  Z = data_matrices$Z,
  X = data_matrices$X,
  time_points = 4,
  mixture_components = 3,
  td_covariates = 2,
  baseline_covariates = 2,
  td_covariate_names = c("waist", "bmi"),
  log_transform_mixtures = TRUE
)

# Step 4: Run analysis
results <- run_gbkmr_panel(
  sim_popn = prepared_data,
  T = 4, p = 3,
  mediator_basenames = c("waist", "bmi"),
  common_covariates = c("sex", "age", "waist0", "bmi0"),
  n = 500,
  K = 1000,
  iter = 24000,
  sel = seq(22000, 24000, by = 25),
  n_knots = 50
)

# Step 5: Interpret results
cat("=== CAUSAL EFFECT ESTIMATES ===\n")
cat("ATE:", round(results$diff_gBKMR, 4), "\n")
cat("E[Y | low exposure]:", round(mean(results$Ya), 4), "\n")
cat("E[Y | high exposure]:", round(mean(results$Yastar), 4), "\n")

diff_posterior <- rowMeans(results$Yastar_mat) - rowMeans(results$Ya_mat)
ci <- quantile(diff_posterior, c(0.025, 0.975))
cat("95% CI: [", round(ci[1], 4), ",", round(ci[2], 4), "]\n")
cat("Significant:", sign(ci[1]) == sign(ci[2]), "\n")
cat("Runtime:", results$meta$total_time_minutes, "min\n")
```

### Example 2: Air Pollution Study (Matrix Format)

Simpler study with single time-dependent confounder.

```r
# Study design:
# - 2 pollutants (PM2.5, NO2) at 3 time points
# - 1 time-dependent confounder (temperature)
# - Outcome: respiratory function
# - Baseline: sex, age, baseline temperature

# Prepare matrices directly
Z <- cbind(
  PM25_T0, NO2_T0,
  PM25_T1, NO2_T1,
  PM25_T2, NO2_T2
)

X <- cbind(
  temp_T1, temp_T2,
  sex, age
)

Y <- lung_function

# Prepare and run
prepared_data <- prepare_gbkmr_data(
  Y, Z, X,
  time_points = 3,
  mixture_components = 2,
  td_covariates = 1,
  baseline_covariates = 2,
  td_covariate_names = "temperature"
)

results <- run_gbkmr_panel(
  prepared_data,
  T = 3, p = 2,
  mediator_basenames = "temperature",
  common_covariates = c("sex", "age", "temperature0"),
  iter = 20000,
  sel = seq(18000, 20000, by = 25)
)

# Results
print(results$diff_gBKMR)
```

---

## Function Reference

### Core Functions

#### `prepare_gbkmr_data()`

Converts user matrices to g-BKMR format.

```r
prepare_gbkmr_data(
  Y,                          # Outcome vector
  Z,                          # Exposure matrix (n x Adim*T)
  X,                          # Covariate matrix (n x Ldim*T + baseline)
  time_points,                # Number of time points (T)
  mixture_components,         # Exposures per time (Adim)
  td_covariates,              # TD covariates per time (Ldim)
  baseline_covariates = 1,    # Number of baseline covariates
  td_covariate_names = NULL,  # Names for TD covariates
  log_transform_mixtures = TRUE,
  validate_input = TRUE
)
```

#### `run_gbkmr_panel()`

Main analysis function implementing g-BKMR.

```r
run_gbkmr_panel(
  sim_popn,                   # Prepared data frame
  T = 5,                      # Time points
  p = 3,                      # Exposures per time
  mediator_basenames = c("waist"),
  common_covariates = c("sex", "waist0"),
  currind = 1,                # Random seed
  n = 500,                    # Sample size
  K = 1000,                   # Monte Carlo samples
  sel = seq(22000, 24000, by = 25),  # MCMC samples to use
  iter = 24000,               # MCMC iterations (mediators)
  n_iter = NULL,              # MCMC iterations (outcome)
  n_knots = 50,               # Kernel approximation knots
  verbose_every = 50          # Progress frequency
)
```

**Returns:** List with:
- `diff_gBKMR`: Average treatment effect
- `Ya`, `Yastar`: Counterfactual outcomes
- `Ya_mat`, `Yastar_mat`: Full posterior samples
- `fit_mediators`, `fit_y`: BKMR model objects
- `meta`: Analysis metadata

#### `detect_variable_patterns()`

Automatically detects variable structure in prepared data.

```r
detect_variable_patterns(data, T)
```

### Parameter Guide

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `T` | Integer | 5 | Total time points (including baseline) |
| `p` | Integer | 3 | Exposures per time point |
| `mediator_basenames` | Character | c("waist") | TD confounder names |
| `iter` | Integer | 24000 | MCMC iterations |
| `sel` | Numeric vector | seq(22000, 24000, by=25) | Post-burnin samples |
| `K` | Integer | 1000 | Monte Carlo samples |
| `n` | Integer | 500 | Sample size for analysis |
| `n_knots` | Integer | 50 | Kernel approximation knots |

---

## Methodological Details

### Causal Framework

**Target estimand (g-formula):**

$$E[Y^{\bar{a}}] = \sum_{\bar{l}} E[Y | \bar{A}=\bar{a}, \bar{L}=\bar{l}, C_0] \prod_{t=1}^{T-1} f(l_t | \bar{a}_{0:(t-1)}, \bar{l}_{1:(t-1)}, C_0) f(c_0)$$

Where:
- $Y^{\bar{a}}$ = potential outcome under exposure history $\bar{a}$
- $\bar{L}$ = time-dependent confounders
- $C_0$ = baseline covariates

### BKMR Models

**Mediator models** (for each $t = 1, ..., T-1$):
$$L_t = h_L(\bar{A}_{0:(t-1)}, \bar{L}_{1:(t-1)}) + C_0^T\beta + \epsilon_t$$

**Outcome model:**
$$Y = h_Y(\bar{A}_{0:(T-1)}, \bar{L}_{1:(T-1)}) + C_0^T\theta + \epsilon_Y$$

Where $h(\cdot)$ uses Gaussian kernel:
$$\text{cor}(h_i, h_j) = \exp\left\{-\frac{1}{\rho}\sum_{m=1}^M (z_{im} - z_{jm})^2\right\}$$

### G-Computation Algorithm

1. **Fit BKMR models** for all mediators and outcome
2. **For each MCMC iteration** $j$:
   - **For each Monte Carlo sample** $k = 1, ..., K$:
     - Sequentially sample mediators: $L_t^{(j,k)} \sim P(L_t | \bar{A}, \bar{L}, C_0)$
     - Sample outcome: $Y^{(j,k)} \sim P(Y | \bar{A}, \bar{L}, C_0)$
   - Average: $Y^{(j)} = \frac{1}{K}\sum_k Y^{(j,k)}$
3. **Estimate ATE**: $\widehat{ATE} = \frac{1}{J}\sum_j (Y^{(j)}_{a^*} - Y^{(j)}_a)$

### Default Interventions

- **Low exposure (a)**: All exposures at 25th percentile
- **High exposure (a*)**: All exposures at 75th percentile

---

## Best Practices

### 1. Data Quality Checks

```r
# Check distributions
summary(Z)
summary(X)
hist(Y)

# Check for missing data
sum(is.na(Y))
sum(is.na(Z))
sum(is.na(X))

# Check correlations
cor(Z)
pairs(Z[, 1:min(6, ncol(Z))])
```

### 2. Start with Moderate Settings

```r
# Initial run: balance speed and precision
results_initial <- run_gbkmr_panel(
  ...,
  iter = 15000,
  sel = seq(13000, 15000, by = 50),
  K = 500,
  n_knots = 40
)
```

### 3. Refine Based on Results

```r
# If estimates unstable or need more precision
results_refined <- run_gbkmr_panel(
  ...,
  iter = 30000,
  sel = seq(28000, 30000, by = 20),
  K = 2000,
  n_knots = 60
)
```

### 4. Parameter Tuning Guide

| Situation | Adjustment | Reason |
|-----------|------------|--------|
| Analysis too slow | ↓ `iter`, `K`, or `n` | Reduce computation |
| Unstable estimates | ↑ `K` or `length(sel)` | More samples |
| Memory errors | ↓ `n_knots` or `n` | Lower memory |
| Need precision | ↑ `iter`, narrow `sel` | Better convergence |
| Large dataset (n>1000) | `n=500-800`, `n_knots=40-50` | Efficiency |

### 5. Diagnostic Checks

```r
# Posterior distribution
diff_posterior <- rowMeans(results$Yastar_mat) - rowMeans(results$Ya_mat)
hist(diff_posterior, main = "Posterior Distribution of ATE")

# Trace plots (visual convergence check)
plot(rowMeans(results$Ya_mat), type = "l", main = "E[Y|a] trace")
plot(rowMeans(results$Yastar_mat), type = "l", main = "E[Y|a*] trace")

# Credible interval
ci <- quantile(diff_posterior, c(0.025, 0.975))
cat("95% CI:", ci, "\n")
cat("Contains zero:", ci[1] <= 0 & ci[2] >= 0, "\n")
```

---

## Troubleshooting

### Common Issues

#### 1. "number of design points >= number of candidates"

**Cause:** Too few observations for requested knots.

**Solutions:**
```r
# Reduce knots
run_gbkmr_panel(..., n_knots = 30)

# Or increase sample size
run_gbkmr_panel(..., n = 300)
```

#### 2. "Missing columns: logM1_0, ..."

**Cause:** Variable names don't match expected format.

**Solutions:**
```r
# Check prepared data
names(prepared_data)

# Verify detection
detection <- detect_variable_patterns(prepared_data, T = your_T)
print(detection)

# Re-prepare if needed
```

#### 3. Analysis takes too long

**Solutions:**
```r
# Fast prototyping settings
run_gbkmr_panel(
  ...,
  iter = 5000,
  sel = seq(4000, 5000, by = 100),
  K = 200,
  n = 200,
  n_knots = 25
)
```

#### 4. Memory issues

**Solutions:**
```r
# Reduce memory footprint
run_gbkmr_panel(..., n = 300, K = 500, n_knots = 30)

# Clear memory between runs
rm(list = ls())
gc()
```

#### 5. Dimension mismatch

**Check dimensions:**
```r
T <- 3; p <- 2; Ldim <- 1; n_baseline <- 1
expected_Z_cols <- p * T                      # Should be 6
expected_X_cols <- Ldim * T + n_baseline      # Should be 4

cat("Z:", ncol(Z), "expected", expected_Z_cols, "\n")
cat("X:", ncol(X), "expected", expected_X_cols, "\n")
```

### Getting Help

1. **Check documentation:** `?prepare_gbkmr_data`, `?run_gbkmr_panel`
2. **Verify data structure:** `str(prepared_data)`
3. **Create minimal example:**

```r
set.seed(123)
n <- 100
Y <- rnorm(n)
Z <- matrix(rlnorm(n * 4), n, 4)
X <- matrix(rnorm(n * 2), n, 2)

data_test <- prepare_gbkmr_data(Y, Z, X, 
  time_points = 2, mixture_components = 2,
  td_covariates = 1, baseline_covariates = 1)
```

## Package Structure

```
causalGBKMR/
├── R/
│   ├── 01-validation.R           # Input validation
│   ├── 02-data-preparation.R     # Data formatting
│   ├── 03-variable-detection.R   # Pattern detection
│   └── 04-core-analysis.R        # Main g-BKMR function
├── man/                          # Documentation
├── tests/                        # Unit tests
├── vignettes/                    # Tutorials
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

---

**Version:** 0.1.0 (Development)

**Status:** Under active development. Please report issues or suggestions via GitHub.

**Planned Features:**
- Visualization functions
- Sensitivity analysis tools
- Binary outcome support
- Parallel computing
- Additional intervention types

---

*Last updated: 2024*