# User Guide for causalGBKMR

This comprehensive guide covers how to use the `causalGBKMR` package for causal inference with time-varying environmental mixtures.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Data Preparation](#data-preparation)
3. [Running Analysis](#running-analysis)
4. [Interpreting Results](#interpreting-results)
5. [Advanced Usage](#advanced-usage)
6. [Best Practices](#best-practices)
7. [Troubleshooting](#troubleshooting)

## Getting Started

### Package Overview

`causalGBKMR` implements g-computation with Bayesian Kernel Machine Regression to estimate causal effects of time-varying exposure mixtures while properly accounting for time-varying confounding.

### Basic Workflow

```r
library(causalGBKMR)

# 1. Prepare data
data_gbkmr <- prepare_gbkmr_data(Y, Z, X, ...)

# 2. Run analysis  
results <- gbkmr_run(data_gbkmr, time_points = T)

# 3. Examine results
print(results)
summary(results)
```

## Data Preparation

### Understanding Your Data Structure

Before using the package, you need to organize your data into three components:

- **Y**: Outcome vector (continuous or binary)
- **Z**: Exposure matrix (mixtures × time points)  
- **X**: Covariate matrix (time-dependent + baseline covariates)

### Example Data Organization

For a study with:
- 3 time points (T₀, T₁, T₂)
- 2 metals (As, Cd) measured at each time point
- 1 time-dependent covariate (BMI) measured at T₁, T₂
- 2 baseline covariates (sex, age)

Your data should be:

```r
# Y: Outcome (length n)
Y <- c(outcome_value_1, outcome_value_2, ..., outcome_value_n)

# Z: Exposure matrix (n × 6)
# Columns: As_T0, Cd_T0, As_T1, Cd_T1, As_T2, Cd_T2
Z <- matrix(c(
  # Subject 1: As_T0, Cd_T0, As_T1, Cd_T1, As_T2, Cd_T2
  c(1.2, 0.8, 1.1, 0.9, 1.3, 0.7),
  # Subject 2: ...
  c(0.9, 1.1, 1.0, 1.2, 0.8, 1.0),
  # ... more subjects
), nrow = n, ncol = 6, byrow = TRUE)

# X: Covariate matrix (n × 4)  
# Columns: BMI_T1, BMI_T2, sex, age
X <- matrix(c(
  # Subject 1: BMI_T1, BMI_T2, sex, age
  c(25.2, 25.8, 1, 45),
  # Subject 2: ...
  c(28.1, 28.5, 0, 52),
  # ... more subjects
), nrow = n, ncol = 4, byrow = TRUE)
```

### Preparing Data for Analysis

```r
prepared_data <- prepare_gbkmr_data(
  Y = Y,                          # Outcome vector
  Z = Z,                          # Exposure matrix
  X = X,                          # Covariate matrix
  time_points = 3,                # Number of time points
  mixture_components = 2,         # Number of exposures per time point
  td_covariates = 1,             # Number of time-dependent covariates
  baseline_covariates = 2,        # Number of baseline covariates
  td_covariate_names = "bmi"     # Names for time-dependent covariates
)
```

### Data Validation

The package automatically validates your data:

```r
# Check if your data passes validation
validate_user_matrices(Y, Z, X, T = 3, Adim = 2, Ldim = 1, n_baseline = 2)
```

Common validation errors:
- **Dimension mismatch**: Check that matrix dimensions match parameters
- **Missing values**: The function will warn about NA values
- **Data types**: Ensure Y, Z, X are numeric

## Running Analysis

### Basic Analysis

```r
results <- gbkmr_run(
  data = prepared_data,
  time_points = 3
)
```

### Analysis with Custom Parameters

```r
results <- gbkmr_run(
  data = prepared_data,
  time_points = 3,
  iter = 20000,              # MCMC iterations
  n = 500,                   # Sample size for analysis
  use_knots = TRUE,          # Use kernel approximation
  n_knots = 50,             # Number of knots
  parallel = TRUE,           # Use parallel processing
  verbose = TRUE             # Print progress
)
```

### Parameter Guidance

| Parameter | Description | Default | Guidance |
|-----------|-------------|---------|----------|
| `iter` | MCMC iterations | 15000 | Increase for more precision, decrease for speed |
| `n` | Sample size | min(500, nrow(data)) | Use smaller values for faster analysis |
| `use_knots` | Kernel approximation | TRUE | Set FALSE for small datasets (n < 100) |
| `n_knots` | Number of knots | 50 | Reduce if getting knots-related errors |
| `parallel` | Parallel processing | TRUE | Set FALSE if experiencing memory issues |

## Interpreting Results

### Main Results

```r
print(results)
```

Output interpretation:
- **Causal Effect Estimate**: Expected change in outcome when all exposures change from 25th to 75th percentile
- **95% CI**: Confidence interval for the causal effect
- **Counterfactual Means**: Expected outcomes under low vs high exposure scenarios

### Detailed Results

```r
# Causal effect and confidence interval
effect_size <- results$causal_effect$estimate
ci_lower <- results$causal_effect$lower  
ci_upper <- results$causal_effect$upper

# Statistical significance (approximate)
is_significant <- sign(ci_lower) == sign(ci_upper)

# Effect size interpretation
if (abs(effect_size) < 0.1) {
  cat("Small effect size")
} else if (abs(effect_size) < 0.5) {
  cat("Medium effect size") 
} else {
  cat("Large effect size")
}
```

### Variable Importance

```r
# Get variable importance scores
importance_scores <- results$variable_importance

# Sort by absolute importance
sorted_importance <- sort(abs(importance_scores), decreasing = TRUE)
print(head(sorted_importance, 10))  # Top 10 most important variables
```

### Summary Statistics

```r
summary(results)
```

This provides:
- Analysis details (iterations, sample size, etc.)
- Top important variables
- Effect size interpretation
- Statistical significance assessment

## Advanced Usage

### Multiple Analyses

Compare different exposure scenarios:

```r
# Analysis 1: Standard percentiles (25th vs 75th)
results_standard <- gbkmr_run(data = prepared_data, time_points = 3)

# Analysis 2: More extreme comparison (need custom implementation)
# Note: Package currently supports 25th vs 75th percentile comparison
```

### Working with Different Data Types

#### Binary Outcomes

```r
# For binary outcomes, prepare data normally
results_binary <- gbkmr_run(
  data = prepared_data,
  outcome_type = "binary",
  time_points = 3
)
```

#### Multiple Time-Dependent Covariates

```r
# Prepare data with multiple TD covariates
prepared_data_multi <- prepare_gbkmr_data(
  Y = Y, Z = Z, X = X,
  time_points = 3,
  mixture_components = 2, 
  td_covariates = 3,                           # 3 TD covariates
  baseline_covariates = 2,
  td_covariate_names = c("bmi", "bp", "glucose")
)
```

### Memory and Performance Optimization

#### For Large Datasets

```r
# Use smaller sample size
results <- gbkmr_run(
  data = prepared_data,
  time_points = 3,
  n = 300,                  # Smaller sample
  iter = 10000,            # Fewer iterations
  use_knots = TRUE,        # Use approximation
  n_knots = 30            # Fewer knots
)
```

#### For Fast Prototyping

```r
# Quick test run
results_quick <- gbkmr_run(
  data = prepared_data,
  time_points = 3,
  iter = 1000,             # Very few iterations
  n = 100,                 # Small sample
  use_knots = FALSE,       # No approximation
  verbose = FALSE          # No progress output
)
```

## Best Practices

### Data Quality

1. **Check for outliers** in exposure and covariate data
2. **Handle missing values** before analysis
3. **Validate temporal ordering** of variables
4. **Use appropriate transformations** (e.g., log-transform exposures)

```r
# Example data quality checks
summary(Z)  # Check exposure distributions
sum(is.na(X))  # Count missing values
cor(Z)  # Check exposure correlations
```

### Analysis Strategy

1. **Start with exploratory analysis**
   ```r
   # Check variable detection
   detection_info <- detect_variable_patterns(prepared_data, T = 3)
   print(detection_info)
   ```

2. **Run initial analysis with default settings**
   ```r
   initial_results <- gbkmr_run(data = prepared_data, time_points = 3)
   ```

3. **Refine parameters based on initial results**
   ```r
   # If initial analysis suggests need for more precision
   refined_results <- gbkmr_run(
     data = prepared_data, 
     time_points = 3,
     iter = 25000,
     n = 600
   )
   ```

### Reporting Results

Include in your analysis report:

1. **Data description**: Sample size, time points, variables
2. **Analysis parameters**: MCMC iterations, sample size used
3. **Main results**: Causal effect estimate with CI
4. **Variable importance**: Top contributing variables
5. **Sensitivity analysis**: Results with different parameters

Example results summary:
```r
cat("g-BKMR Analysis Results\n")
cat("Sample size:", nrow(prepared_data), "\n")
cat("Time points:", 3, "\n") 
cat("Causal effect (95% CI):", 
    round(results$causal_effect$estimate, 3), " (",
    round(results$causal_effect$lower, 3), ", ",
    round(results$causal_effect$upper, 3), ")\n")
```

## Troubleshooting

### Common Issues and Solutions

#### Issue 1: "number of design points >= the number of candidates"

**Cause**: Too few data points for the number of knots requested

**Solutions**:
```r
# Option 1: Disable knots
results <- gbkmr_run(data = prepared_data, time_points = 3, use_knots = FALSE)

# Option 2: Reduce number of knots  
results <- gbkmr_run(data = prepared_data, time_points = 3, n_knots = 20)

# Option 3: Use larger sample size
results <- gbkmr_run(data = prepared_data, time_points = 3, n = 200)
```

#### Issue 2: Analysis takes too long

**Solutions**:
```r
# Reduce iterations
results <- gbkmr_run(data = prepared_data, time_points = 3, iter = 10000)

# Use smaller sample
results <- gbkmr_run(data = prepared_data, time_points = 3, n = 300)

# Enable knots for approximation
results <- gbkmr_run(data = prepared_data, time_points = 3, use_knots = TRUE)
```

#### Issue 3: Memory issues

**Solutions**:
```r
# Disable parallel processing
results <- gbkmr_run(data = prepared_data, time_points = 3, parallel = FALSE)

# Use smaller sample size
results <- gbkmr_run(data = prepared_data, time_points = 3, n = 200)

# Clear memory between analyses
gc()
```

#### Issue 4: Variable detection errors

**Check your data structure**:
```r
# Verify column names
names(prepared_data)

# Check for expected patterns
grep("logM\\d+_\\d+", names(prepared_data), value = TRUE)  # Exposure variables
grep("_\\d+$", names(prepared_data), value = TRUE)        # Time-varying variables
```

#### Issue 5: Convergence issues

**Solutions**:
```r
# Increase iterations
results <- gbkmr_run(data = prepared_data, time_points = 3, iter = 30000)

# Check for data quality issues
summary(prepared_data)
```

### Diagnostic Checks

After analysis, perform these checks:

```r
# 1. Check convergence (if possible with more iterations)
# This requires examining raw BKMR output

# 2. Check effect size reasonableness  
effect_size <- results$causal_effect$estimate
cat("Effect size:", effect_size, "\n")

# 3. Check confidence interval width
ci_width <- results$causal_effect$upper - results$causal_effect$lower
cat("CI width:", ci_width, "\n")

# 4. Examine variable importance
importance <- results$variable_importance
print(head(sort(abs(importance), decreasing = TRUE), 10))
```

### Getting Additional Help

If you encounter issues not covered here:

1. **Check function documentation**:
   ```r
   ?gbkmr_run
   ?prepare_gbkmr_data
   ?detect_variable_patterns
   ```

2. **Run built-in examples**:
   ```r
   example(gbkmr_run)
   example(prepare_gbkmr_data)
   ```

3. **Search existing issues** on GitHub

4. **Create a reproducible example**:
   ```r
   # Minimal example that reproduces your issue
   set.seed(123)
   n <- 50
   Y <- rnorm(n)
   Z <- matrix(rlnorm(n * 4), nrow = n, ncol = 4)
   X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
   
   data_test <- prepare_gbkmr_data(Y, Z, X, 
                                   time_points = 2, 
                                   mixture_components = 2,
                                   td_covariates = 1, 
                                   baseline_covariates = 1)
   
   # Your problematic code here
   ```

## Example Workflows

### Workflow 1: Environmental Health Study

```r
# Study: Effect of metal mixtures on cardiovascular outcomes
# 3 metals (As, Cd, Pb), 4 time points, BMI as time-dependent confounder

library(causalGBKMR)

# Load your data
# Y: cardiovascular risk score
# Z: metal concentrations (n × 12 matrix: 3 metals × 4 time points)  
# X: BMI measurements + baseline covariates (n × 6 matrix: 3 BMI + 3 baseline)

# Step 1: Prepare data
cv_data <- prepare_gbkmr_data(
  Y = cardiovascular_score,
  Z = metal_concentrations,
  X = covariates_matrix,
  time_points = 4,
  mixture_components = 3,
  td_covariates = 1,
  baseline_covariates = 3,
  td_covariate_names = "bmi"
)

# Step 2: Run analysis
cv_results <- gbkmr_run(
  data = cv_data,
  time_points = 4,
  iter = 20000,
  verbose = TRUE
)

# Step 3: Interpret results
print(cv_results)
summary(cv_results)

# Step 4: Extract key findings
effect_estimate <- cv_results$causal_effect$estimate
ci_bounds <- c(cv_results$causal_effect$lower, cv_results$causal_effect$upper)
top_variables <- head(sort(abs(cv_results$variable_importance), decreasing = TRUE), 5)
```

### Workflow 2: Air Pollution Study

```r
# Study: Effect of air pollutants on respiratory function
# 4 pollutants (PM2.5, NO2, O3, SO2), 3 time points, multiple confounders

# Step 1: Prepare data with multiple time-dependent confounders
resp_data <- prepare_gbkmr_data(
  Y = lung_function,
  Z = pollutant_matrix,        # n × 12: 4 pollutants × 3 time points
  X = confounder_matrix,       # n × 8: 2 TD covariates × 3 time points + 2 baseline
  time_points = 3,
  mixture_components = 4,
  td_covariates = 2,
  baseline_covariates = 2,
  td_covariate_names = c("temperature", "humidity")
)

# Step 2: Check variable detection
detection <- detect_variable_patterns(resp_data, T = 3)
print(detection)

# Step 3: Run analysis with optimization for larger dataset
resp_results <- gbkmr_run(
  data = resp_data,
  time_points = 3,
  iter = 15000,
  n = 800,                    # Use larger sample
  use_knots = TRUE,
  n_knots = 60,
  parallel = TRUE
)

# Step 4: Detailed examination
summary(resp_results)

# Check effect size categories
effect <- resp_results$causal_effect$estimate
if (abs(effect) < 0.05) {
  cat("Small clinical effect\n")
} else if (abs(effect) < 0.2) {
  cat("Moderate clinical effect\n") 
} else {
  cat("Large clinical effect\n")
}
```

### Workflow 3: Sensitivity Analysis

```r
# Compare results with different analysis parameters

# Standard analysis
results_standard <- gbkmr_run(data = prepared_data, time_points = 3)

# Higher precision analysis
results_precise <- gbkmr_run(
  data = prepared_data, 
  time_points = 3,
  iter = 30000,
  n = 700
)

# Different knot settings
results_no_knots <- gbkmr_run(
  data = prepared_data,
  time_points = 3, 
  use_knots = FALSE
)

# Compare results
comparison <- data.frame(
  Analysis = c("Standard", "High Precision", "No Knots"),
  Estimate = c(results_standard$causal_effect$estimate,
               results_precise$causal_effect$estimate,
               results_no_knots$causal_effect$estimate),
  CI_Lower = c(results_standard$causal_effect$lower,
               results_precise$causal_effect$lower, 
               results_no_knots$causal_effect$lower),
  CI_Upper = c(results_standard$causal_effect$upper,
               results_precise$causal_effect$upper,
               results_no_knots$causal_effect$upper)
)

print(comparison)
```

## Conclusion

This guide covers the essential aspects of using `causalGBKMR` for your research. The package provides a powerful framework for causal inference with time-varying environmental mixtures, but requires careful attention to data preparation and parameter selection.

Key takeaways:
- **Data preparation is crucial** - ensure your matrices are correctly organized
- **Start with default parameters** and refine based on your specific needs  
- **Validate your results** through sensitivity analyses
- **Interpret effects in context** of your specific research question

For additional help, consult the function documentation, GitHub issues, or contact the package maintainers.