## Troubleshooting

### Common Issues

#### 1. "number of design points >= number of candidates" (knots error)

**Cause**: Too few data points for the requested number of knots

**Solutions**:
```r
# Option 1: Reduce number of knots
results <- run_gbkmr_panel(..., n_knots = 30)

# Option 2: Increase sample size
results <- run_gbkmr_panel(..., n = 300)

# Option 3: Both
results <- run_gbkmr_panel(..., n = 400, n_knots = 40)
```

#### 2. "Missing columns: logM1_0, logM2_0, ..."

**Cause**: Variable names don't match expected g-BKMR format

**Solutions**:
```r
# Check your prepared data
names(prepared_data)

# Verify column detection
detection <- detect_variable_patterns(prepared_data, T = your_T)
print(detection)

# Re-prepare data if needed with correct parameters
prepared_data <- prepare_gbkmr_data(
  Y, Z, X,
  time_points = correct_T,
  mixture_components = correct_p,
  td_covariates = correct_Ldim,
  ...
)
```

#### 3. Analysis takes too long

**Solutions**:
```r
# Reduce MCMC iterations
results <- run_gbkmr_panel(..., iter = 10000, sel = seq(8000, 10000, by = 50))

# Reduce Monte Carlo samples
results <- run_gbkmr_panel(..., K = 500)

# Use smaller sample
results <- run_gbkmr_panel(..., n = 300)

# Reduce number of knots
results <- run_gbkmr_panel(..., n_knots = 30)

# Combine strategies for fast prototyping
results_fast <- run_gbkmr_panel(
  ...,
  iter = 5000,
  sel = seq(4000, 5000, by = 100),
  K = 200,
  n = 200,
  n_knots = 25
)
```

#### 4. Memory issues / R crashes

**Solutions**:
```r
# Reduce sample size
results <- run_gbkmr_panel(..., n = 300)

# Reduce Monte Carlo samples  
results <- run_gbkmr_panel(..., K = 500)

# Reduce posterior samples
results <- run_gbkmr_panel(..., sel = seq(22000, 24000, by = 50))

# Clear memory between runs
rm(list = ls())
gc()
```

#### 5. "sel contains indices beyond total MCMC iterations"

**Cause**: Selection indices exceed the number of iterations

**Solution**:
```r
# Make sure max(sel) <= iter
results <- run_gbkmr_panel(
  ...,
  iter = 24000,
  sel = seq(22000, 24000, by = 25)  # max is 24000
)
```

#### 6. Dimension mismatch errors in prepare_gbkmr_data()

**Cause**: Matrix dimensions don't match specified parameters

**Check your dimensions**:
```r
# For T time points, p exposures, Ldim TD covariates, n_baseline baseline covs:
# Z should be: n × (p * T)
# X should be: n × (Ldim * T + n_baseline)

# Example check:
T <- 3
p <- 2  
Ldim <- 1
n_baseline <- 1

expected_Z_cols <- p * T          # Should be 6
expected_X_cols <- Ldim * T + n_baseline  # Should be 4

cat("Z has", ncol(Z), "columns, expected", expected_Z_cols, "\n")
cat("X has", ncol(X), "columns, expected", expected_X_cols, "\n")
```

### Diagnostic Checks

After running analysis:

```r
# 1. Check if results are reasonable
cat("ATE:", results$diff_gBKMR, "\n")
cat("E[Y|a]:", mean(results$Ya), "\n")
cat("E[Y|a*]:", mean(results$Yastar), "\n")

# 2. Check posterior distribution
diff_posterior <- rowMeans(results$Yastar_mat) - rowMeans(results$Ya_mat)
hist(diff_posterior, main = "Posterior Distribution of ATE", xlab = "ATE")
abline(v = results$diff_gBKMR, col = "red", lwd = 2)

# 3. Check convergence (visual inspection of samples)
plot(rowMeans(results$Ya_mat), type = "l", 
     main = "E[Y|a] across MCMC iterations",
     ylab = "E[Y|a]", xlab = "Iteration")
plot(rowMeans(results$Yastar_mat), type = "l",
     main = "E[Y|a*] across MCMC iterations", 
     ylab = "E[Y|a*]", xlab = "Iteration")

# 4. Check credible interval
ci <- quantile(diff_posterior, c(0.025, 0.975))
cat("95% CI: [", ci[1], ",", ci[2], "]\n")
cat("Contains zero:", ci[1] <= 0 & ci[2] >= 0, "\n")

# 5. Examine posterior mean coefficients
cat("\nOutcome model coefficients (beta_y):\n")
print(results$beta_y)
```

### Getting Help

If you encounter issues not covered here:

1. **Check function documentation**:
   ```r
   ?prepare_gbkmr_data
   ?run_gbkmr_panel
   ?detect_variable_patterns
   ```

2. **Verify data structure**:
   ```r
   str(prepared_data)
   head(prepared_data)
   detection <- detect_variable_patterns(prepared_data, T = your_T)
   print(detection)
   ```

3. **Create minimal reproducible example**:
   ```r
   # Minimal example
   set.seed(123)
   n <- 100
   Y <- rnorm(n)
   Z <- matrix(rlnorm(n * 4), n, 4)  # 2 exposures × 2 time points
   X <- matrix(rnorm(n * 2), n, 2)   # 1 TD cov × 1 time point + 1 baseline
   
   data_test <- prepare_gbkmr_data(
     Y, Z, X,
     time_points = 2,
     mixture_components = 2,
     td_covariates = 1,
     baseline_covariates = 1
   )
   
   # Try your problematic code here
   ```

4. **Search GitHub issues** for similar problems

5. **Contact package maintainers** with reproducible example
- Account for time-varying confounding
- Model non-linear and non-additive exposure-response relationships
- Perform variable selection in high-dimensional mixture settings

# causalGBKMR: g-BKMR for Causal Inference with Time-Varying Environmental Mixtures

An R package implementing g-BKMR (g-formula with Bayesian Kernel Machine Regression) for estimating causal effects of time-varying correlated environmental mixtures on health outcomes.

## Overview

`causalGBKMR` combines the g-formula for handling time-varying confounding with Bayesian Kernel Machine Regression (BKMR) for flexible modeling of mixture effects. This approach allows researchers to:

- Estimate joint effects of multiple, correlated, time-varying exposures
- Account for time-varying confounding
- Model non-linear and non-additive exposure-response relationships
- Perform variable selection in high-dimensional mixture settings

## Installation

```r
# Install dependencies
install.packages("bkmr")
install.packages("dplyr")
install.packages("fields")

# Install from GitHub (when available)
# devtools::install_github("username/causalGBKMR")

# For now, source the R files directly
source("R/01-validation.R")
source("R/02-data-preparation.R")
source("R/03-variable-detection.R")
source("R/04-core-analysis.R")
```

## Quick Start

### 1. Prepare Your Data

#### Recommended Input Structure (List Format)

The package supports flexible input via lists for better organization:

```r
library(causalGBKMR)

# Example with 3 time points, 3 exposures per time, 2 time-dependent covariates

# Time-dependent covariates organized by covariate name
TD_covariates <- list(
  waist = list(
    t1 = waist_measurements_t1,
    t2 = waist_measurements_t2,
    t3 = waist_measurements_t3
  ),
  bmi = list(
    t1 = bmi_measurements_t1,
    t2 = bmi_measurements_t2,
    t3 = bmi_measurements_t3
  )
)

# Exposures organized by time point
exposures <- list(
  t0 = exposure_matrix_t0,  # n x p matrix (p exposures)
  t1 = exposure_matrix_t1,
  t2 = exposure_matrix_t2,
  t3 = exposure_matrix_t3
)

# Baseline confounders
C0 <- baseline_covariates  # n x k matrix (k baseline covariates)

# Outcome
Y <- outcome_vector  # length n
```

#### Alternative: Matrix Format (Legacy)

For backward compatibility, you can also provide data as matrices:

```r
# Y: Outcome vector (n x 1)
Y <- rnorm(200)

# Z: Exposure matrix (n x Adim*T)
# Organized as: [M1_T0, M2_T0, ..., Mp_T0, M1_T1, ..., Mp_T(T-1)]
Z <- matrix(rlnorm(200 * 6), nrow = 200, ncol = 6)  # 2 exposures x 3 time points

# X: Covariate matrix (n x Ldim*T + baseline)
# Organized as: [TD1_T1, TD2_T1, ..., TDLdim_T(T-1), Baseline1, Baseline2, ...]
X <- matrix(rnorm(200 * 3), nrow = 200, ncol = 3)   # 1 TD cov x 2 times + 1 baseline

# Prepare data for g-BKMR
prepared_data <- prepare_gbkmr_data(
  Y = Y,
  Z = Z,
  X = X,
  time_points = 3,
  mixture_components = 2,
  td_covariates = 1,
  baseline_covariates = 1,
  td_covariate_names = "waist"
)
```

### 2. Run g-BKMR Analysis

```r
# Run the analysis
results <- run_gbkmr_panel(
  sim_popn = prepared_data,
  T = 3,                                    # Time points
  p = 2,                                    # Exposures per time point
  mediator_basenames = c("waist"),          # Time-dependent confounders
  common_covariates = c("sex", "waist0"),   # Baseline covariates
  iter = 24000,                             # MCMC iterations
  sel = seq(22000, 24000, by = 25),        # Post burn-in samples
  K = 1000,                                 # Monte Carlo samples
  n_knots = 50,                            # Kernel approximation knots
  verbose_every = 50                        # Progress frequency
)

# View results
print(results$diff_gBKMR)  # Average treatment effect
```

### 3. Examine Results

```r
# Causal effect estimate
cat("ATE:", results$diff_gBKMR, "\n")

# Counterfactual means
cat("E[Y under low exposure]:", mean(results$Ya), "\n")
cat("E[Y under high exposure]:", mean(results$Yastar), "\n")

# 95% Credible interval (approximate)
diff_samples <- rowMeans(results$Yastar_mat) - rowMeans(results$Ya_mat)
ci <- quantile(diff_samples, c(0.025, 0.975))
cat("95% CI:", ci[1], "to", ci[2], "\n")
```

## Data Structure Details

### Understanding the Input Formats

#### List Format (Recommended)

This format provides better organization and readability:

```r
# Time-dependent covariates: nested list
# Structure: TD_covariates[[covariate_name]][[time_point]]
TD_covariates <- list(
  parity = list(t1 = C_parity_t1, t2 = C_parity_t2, t3 = C_parity_t3),
  breastfeeding = list(t1 = C_bf_t1, t2 = C_bf_t2, t3 = C_bf_t3),
  BMI = list(t1 = C_bmi_t1, t2 = C_bmi_t2, t3 = C_bmi_t3)
)

# Exposures: list by time point, each element is n x p matrix
exposures <- list(
  t0 = X_exp_t0,  # Matrix: n x p (p = number of mixture components)
  t1 = X_exp_t1,
  t2 = X_exp_t2,
  t3 = X_exp_t3
)

# Baseline confounders: single matrix
C0 <- C_baseline  # Matrix: n x k (k = number of baseline covariates)

# Outcome: vector
Y <- outcome  # Vector: length n
```

**Advantages of list format:**
- Clear separation of different covariates
- Easy to add/remove variables
- Natural organization by time
- Self-documenting variable names

#### Matrix Format (Wide Format)

Traditional format for backward compatibility:

**Matrix Z (Exposures)**
```
Columns: [M1_T0, M2_T0, ..., Mp_T0, M1_T1, M2_T1, ..., Mp_T1, ..., M1_T(T-1), ..., Mp_T(T-1)]
Dimensions: n x (p * T)
```

Example for 2 metals x 3 time points:
```
[As_T0, Cd_T0, As_T1, Cd_T1, As_T2, Cd_T2]
```

**Matrix X (Covariates)**
```
Columns: [TD1_T1, TD2_T1, ..., TDLdim_T1, ..., TDLdim_T(T-1), Baseline1, Baseline2, ...]
Dimensions: n x (Ldim * (T-1) + n_baseline)
```

Example for 1 TD covariate x 2 time points + 1 baseline:
```
[waist_T1, waist_T2, sex]
```

**Note on time indexing:**
- Exposures (Z) include baseline (T0): indices 0 to T-1
- Time-dependent covariates (X) start at T1: indices 1 to T-1
- This reflects that TD covariates are measured AFTER baseline exposure

### Output Format from prepare_gbkmr_data()

The prepared data frame has variables in this order:

1. **sex**: First baseline covariate (required)
2. **baseline_2, baseline_3, ...**: Additional baseline covariates
3. **waist_0** (or specified TD covariate names with _0): Baseline TD covariates
4. **logM1_0, logM2_0, ...**: Log-transformed exposures at time 0
5. **logM1_1, logM2_1, ...**: Log-transformed exposures at time 1
6. **waist_1, waist_2, ...**: Time-dependent covariates at times 1, 2, ...
7. **Y**: Outcome variable
8. **id**: Subject identifier

Example column structure for T=3, p=2, Ldim=1:
```
[sex, waist_0, logM1_0, logM2_0, logM1_1, logM2_1, waist_1, logM1_2, logM2_2, waist_2, Y, id]
```

## Key Functions

### Data Preparation
- **`prepare_gbkmr_data()`**: Convert user matrices (Y, Z, X) to g-BKMR format
  - Handles log transformation of exposures
  - Creates proper variable naming (logM1_0, waist_1, etc.)
  - Adds metadata for downstream analysis
  - Validates input dimensions

- **`validate_user_matrices()`**: Validate input data structure
  - Checks dimensions match specified parameters
  - Warns about missing values
  - Ensures data types are correct

- **`detect_variable_patterns()`**: Automatically detect variable structure
  - Identifies exposure variables (logM patterns)
  - Detects time-dependent covariate names
  - Recognizes multiple naming conventions
  - Returns structure information for analysis

### Core Analysis
- **`run_gbkmr_panel()`**: Main analysis function
  - Fits BKMR models for mediators and outcome
  - Implements g-computation with sequential sampling
  - Estimates counterfactual outcomes under interventions
  - Returns comprehensive results with timing information

### Analysis Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `T` | Integer | 5 | Total time points (including t=0) |
| `p` | Integer | 3 | Number of exposures per time point |
| `mediator_basenames` | Character vector | c("waist") | Names of time-dependent confounders |
| `common_covariates` | Character vector | c("sex", "waist0") | Baseline covariate names |
| `n` | Integer | 500 | Sample size for analysis |
| `K` | Integer | 1000 | Monte Carlo samples for g-computation |
| `iter` | Integer | 24000 | MCMC iterations for mediator models |
| `n_iter` | Integer | NULL | MCMC iterations for outcome (uses `iter` if NULL) |
| `sel` | Numeric vector | seq(22000, 24000, by=25) | Post burn-in MCMC samples |
| `n_knots` | Integer | 50 | Knots for kernel approximation |
| `verbose_every` | Integer | 50 | Progress report frequency |

## Real Data Example

### Environmental Health Study: Metal Mixtures and Cardiovascular Outcomes

#### Example 1: Using List Format (Recommended)

```r
# Study design:
# - 3 metals (As, Cd, Pb) measured at 4 time points (baseline + 3 follow-ups)
# - Waist circumference and BMI as time-dependent confounders
# - Blood glucose as outcome
# - Baseline covariates: sex, age

# Step 1: Organize your data into lists

# Time-dependent covariates
TD_covariates <- list(
  waist = list(
    t1 = waist_visit1,  # Vector of length n
    t2 = waist_visit2,
    t3 = waist_visit3
  ),
  bmi = list(
    t1 = bmi_visit1,
    t2 = bmi_visit2,
    t3 = bmi_visit3
  )
)

# Exposures by time point (each is n x 3 matrix for 3 metals)
exposures <- list(
  t0 = cbind(As_baseline, Cd_baseline, Pb_baseline),
  t1 = cbind(As_visit1, Cd_visit1, Pb_visit1),
  t2 = cbind(As_visit2, Cd_visit2, Pb_visit2),
  t3 = cbind(As_visit3, Cd_visit3, Pb_visit3)
)

# Baseline confounders
C0 <- cbind(sex, age)
colnames(C0) <- c("sex", "age")

# Outcome
Y <- blood_glucose

# Step 2: Convert lists to matrix format for prepare_gbkmr_data()
# Helper function to convert list format to matrix format
convert_lists_to_matrices <- function(exposures, TD_covariates, C0, Y) {
  # Flatten exposure list into matrix
  Z <- do.call(cbind, exposures)  # Combines all time points
  
  # Flatten TD covariates
  # Extract all time points for each covariate
  n_covariates <- length(TD_covariates)
  n_times <- length(TD_covariates[[1]])
  
  X_td_list <- list()
  for (t_idx in 1:n_times) {
    for (cov_name in names(TD_covariates)) {
      X_td_list[[length(X_td_list) + 1]] <- TD_covariates[[cov_name]][[t_idx]]
    }
  }
  X_td <- do.call(cbind, X_td_list)
  
  # Combine TD covariates with baseline
  X <- cbind(X_td, C0)
  
  list(Y = Y, Z = Z, X = X)
}

# Convert
data_matrices <- convert_lists_to_matrices(exposures, TD_covariates, C0, Y)

# Step 3: Prepare data for g-BKMR
prepared_data <- prepare_gbkmr_data(
  Y = data_matrices$Y,
  Z = data_matrices$Z,
  X = data_matrices$X,
  time_points = 4,
  mixture_components = 3,
  td_covariates = 2,
  baseline_covariates = 2,
  td_covariate_names = c("waist", "bmi"),
  log_transform_mixtures = TRUE,
  validate_input = TRUE
)

# Step 4: Run g-BKMR analysis
results <- run_gbkmr_panel(
  sim_popn = prepared_data,
  T = 4,
  p = 3,
  mediator_basenames = c("waist", "bmi"),
  common_covariates = c("sex", "age", "waist0", "bmi0"),
  n = min(500, nrow(prepared_data)),
  K = 1000,
  iter = 24000,
  sel = seq(22000, 24000, by = 25),
  n_knots = 50,
  verbose_every = 50
)

# Step 5: Extract and interpret results
cat("\n=== RESULTS ===\n")
cat("Average Treatment Effect:", round(results$diff_gBKMR, 4), "\n")
cat("E[Y | low exposure]:", round(mean(results$Ya), 4), "\n")
cat("E[Y | high exposure]:", round(mean(results$Yastar), 4), "\n")

# Calculate 95% credible interval
diff_posterior <- rowMeans(results$Yastar_mat) - rowMeans(results$Ya_mat)
ci <- quantile(diff_posterior, c(0.025, 0.975))
cat("95% CI: [", round(ci[1], 4), ",", round(ci[2], 4), "]\n")

# Check if effect is significant
is_significant <- sign(ci[1]) == sign(ci[2])
cat("Statistically significant:", is_significant, "\n")

# Computation time
cat("Total time:", results$meta$total_time_minutes, "minutes\n")
```

#### Example 2: Using Matrix Format (Legacy)

```r
# Study design:
# - 3 metals (As, Cd, Pb) measured at 3 time points
# - Waist circumference as time-dependent confounder
# - Blood glucose as outcome
# - Baseline covariates: sex, age, baseline waist

# Prepare exposure matrix Z (n x 9: 3 metals x 3 time points)
# Order: As_T0, Cd_T0, Pb_T0, As_T1, Cd_T1, Pb_T1, As_T2, Cd_T2, Pb_T2
Z <- cbind(
  arsenic_t0, cadmium_t0, lead_t0,
  arsenic_t1, cadmium_t1, lead_t1,
  arsenic_t2, cadmium_t2, lead_t2
)

# Prepare covariate matrix X (n x 4: waist at T1, T2 + sex + age)
# Order: TD covariates first, then baseline
X <- cbind(
  waist_t1, waist_t2,
  sex, age
)

# Outcome
Y <- blood_glucose

# Prepare data
prepared_data <- prepare_gbkmr_data(
  Y = Y,
  Z = Z,
  X = X,
  time_points = 3,
  mixture_components = 3,
  td_covariates = 1,
  baseline_covariates = 2,
  td_covariate_names = "waist",
  log_transform_mixtures = TRUE,
  validate_input = TRUE
)

# Detect variable patterns (optional check)
detection <- detect_variable_patterns(prepared_data, T = 3)
print(detection)

# Run g-BKMR analysis
results <- run_gbkmr_panel(
  sim_popn = prepared_data,
  T = 3,
  p = 3,
  mediator_basenames = c("waist"),
  common_covariates = c("sex", "age", "waist0"),
  n = min(500, nrow(prepared_data)),
  K = 1000,
  iter = 24000,
  sel = seq(22000, 24000, by = 25),
  n_knots = 50,
  verbose_every = 50
)

# Extract results
cat("\n=== RESULTS ===\n")
cat("Average Treatment Effect:", round(results$diff_gBKMR, 4), "\n")
cat("E[Y | low exposure]:", round(mean(results$Ya), 4), "\n")
cat("E[Y | high exposure]:", round(mean(results$Yastar), 4), "\n")

diff_posterior <- rowMeans(results$Yastar_mat) - rowMeans(results$Ya_mat)
ci <- quantile(diff_posterior, c(0.025, 0.975))
cat("95% CI: [", round(ci[1], 4), ",", round(ci[2], 4), "]\n")
cat("Significant:", sign(ci[1]) == sign(ci[2]), "\n")
cat("Total time:", results$meta$total_time_minutes, "minutes\n")
```

## Methodological Details

### Causal Framework

The package implements g-computation to estimate causal effects under time-varying confounding:

**Target estimand:**
$E[Y^{\bar{a}}] = \sum_{\bar{l}} E[Y | \bar{A}=\bar{a}, \bar{L}=\bar{l}, C_0] \prod_{t=1}^{T-1} f(l_t | \bar{a}_{0:(t-1)}, \bar{l}_{1:(t-1)}, C_0) f(c_0)$

Where:
- $Y^{\bar{a}}$ is the potential outcome under exposure history $\bar{a}$
- $\bar{L} = (L_1, ..., L_{T-1})$ are time-dependent confounders
- $C_0$ are baseline covariates

### BKMR Models

**Mediator (confounder) models** at each time point $t = 1, ..., T-1$:
$L_t = h_L(\bar{A}_{0:(t-1)}, \bar{L}_{1:(t-1)}) + C_0^T\beta + \epsilon_t$

**Outcome model:**
$Y = h_Y(\bar{A}_{0:(T-1)}, \bar{L}_{1:(T-1)}) + C_0^T\theta + \epsilon_Y$

Where $h(\cdot)$ is flexibly modeled using Gaussian kernel machine regression:
$\text{cor}(h_i, h_j) = \exp\left\{-\frac{1}{\rho}\sum_{m=1}^M (z_{im} - z_{jm})^2\right\}$

### G-Computation Algorithm

1. **Fit BKMR models** for all mediators and outcome
2. **For each MCMC iteration** $j \in \text{sel}$:
   - **For each Monte Carlo sample** $k = 1, ..., K$:
     - Sequentially sample mediators forward through time:
       $L_t^{(j,k)} \sim P(L_t | \bar{A}_{0:(t-1)}=\bar{a}, \bar{L}_{1:(t-1)}^{(j,k)}, C_0)$
     - Sample final outcome:
       $Y^{(j,k)} \sim P(Y | \bar{A}_{0:(T-1)}=\bar{a}, \bar{L}_{1:(T-1)}^{(j,k)}, C_0)$
   - Average over Monte Carlo samples: $Y^{(j)} = \frac{1}{K}\sum_k Y^{(j,k)}$
3. **Estimate ATE**: $\widehat{ATE} = \frac{1}{|\text{sel}|}\sum_j (Y^{(j)}_{a^*} - Y^{(j)}_a)$

### Interventions

By default, the package compares:
- **Low exposure (a)**: All exposures at 25th percentile
- **High exposure (a*)**: All exposures at 75th percentile

These interventions represent realistic policy-relevant scenarios.

## Citation

If you use this package in your research, please cite:

```bibtex
@article{chai2024gbkmr,
  title={g-BKMR: Causal Inference for health effects of time-varying correlated environmental mixtures},
  author={Chai, Zilan and Navas-Acien, Ana and Coull, Brent and Valeri, Linda},
  journal={Under Review},
  year={2024}
}
```

## References

### Key Papers

1. **g-BKMR methodology**:
   - Chai, Z., Navas-Acien, A., Coull, B., & Valeri, L. (2024). g-BKMR: Causal Inference for health effects of time-varying correlated environmental mixtures. *Under Review*.

2. **BKMR foundation**:
   - Bobb, J.F., Valeri, L., Claus Henn, B., et al. (2015). Bayesian kernel machine regression for estimating the health effects of multi-pollutant mixtures. *Biostatistics*, 16(3), 493-508.
   - Bobb, J.F., Claus Henn, B., Valeri, L., & Coull, B.A. (2018). Statistical software for analyzing the health effects of multiple concurrent exposures via Bayesian kernel machine regression. *Environmental Health*, 17(1), 1-10.

3. **G-formula and causal inference**:
   - Robins, J. (1986). A new approach to causal inference in mortality studies with sustained exposure periods. *Mathematical Modelling*, 7(9-12), 1393-1512.
   - Hernán, M.A., & Robins, J.M. (2020). *Causal Inference: What If*. Chapman & Hall/CRC.
   - Naimi, A.I., Cole, S.R., & Kennedy, E.H. (2017). An introduction to g methods. *International Journal of Epidemiology*, 46(2), 756-762.

4. **Time-varying confounding**:
   - Hernán, M.A., Brumback, B., & Robins, J.M. (2000). Marginal structural models to estimate the causal effect of zidovudine on the survival of HIV-positive men. *Epidemiology*, 11(5), 561-570.

5. **Environmental mixtures**:
   - Billionnet, C., Sherrill, D., & Annesi-Maesano, I. (2012). Estimating the health effects of exposure to multi-pollutant mixture. *Annals of Epidemiology*, 22(2), 126-141.
   - Carrico, C., Gennings, C., Wheeler, D.C., & Factor-Litvak, P. (2015). Characterization of weighted quantile sum regression for highly correlated data in a risk analysis setting. *Journal of Agricultural, Biological, and Environmental Statistics*, 20, 100-120.

## Package Structure

```
causalGBKMR/
├── R/
│   ├── 01-validation.R           # Input validation functions
│   ├── 02-data-preparation.R     # Data formatting and preparation
│   ├── 03-variable-detection.R   # Automatic variable pattern detection
│   └── 04-core-analysis.R        # Main g-BKMR analysis function
├── man/                           # Documentation files
├── tests/                         # Unit tests
├── vignettes/                     # Tutorial vignettes
├── DESCRIPTION                    # Package metadata
├── NAMESPACE                      # Package exports
└── README.md                      # This file
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/your-feature`)
5. Create a Pull Request

Please include:
- Clear description of changes
- Tests for new functionality
- Updated documentation
- Examples demonstrating new features

## License

MIT License

Copyright (c) 2024 Zilan Chai, Linda Valeri

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

## Contact

**Maintainer**: Zilan Chai, Linda Valeri

**Institution**: Columbia University Mailman School of Public Health

**Email**: lv2424@cumc.columbia.edu

For bug reports and feature requests, please open an issue on GitHub.

## Acknowledgments

This work was supported by [funding information]. We thank the Strong Heart Study participants and research team for their contributions to the application example.

## Version History

### Version 0.1.0 (Current)
- Initial release
- Core g-BKMR functionality
- Data preparation utilities
- Automatic variable detection
- Comprehensive documentation

### Planned Features
- Visualization functions for results
- Sensitivity analysis tools
- Support for binary outcomes
- Parallel computing options
- Additional intervention types

---

**Note**: This package is under active development. Please report any issues or suggestions via GitHub Issues.