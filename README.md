### Advanced Usage (Statistician Mode)

For users who need more control over the modeling process:

```r
# Advanced analysis with custom parameters
results_advanced <- gbkmr_run(
  data = prepared_data,
  outcome = "Y",
  outcome_type = "continuous",
  time_points = 4,
  
  #### Performance Optimization

```r
# For large datasets and complex models
optimize_gbkmr <- function(sim_popn, T, Adim, Ldim) {
  
  # Recommended settings based on data complexity
  if (nrow(sim_popn) > 10000 && Adim * T > 10) {
    # Large, complex data
    settings <- list(
      n = min(1000, nrow(sim_popn) * 0.1),    # Sample 10% or max 1000
      iter = 20000,                            # More iterations for stability
      sel = seq(15000, 20000, by = 25),       # Conservative burn-in
      use_knots = TRUE,                        # Essential for efficiency
      n_knots = 100,                          # More knots for complex relationships
      K = 2000                                # More prediction samples
    )
  } else if (Adim * Ldim * T > 20) {
    # High-dimensional but smaller sample
    settings <- list(
      n = min(500, nrow(sim_popn) * 0.2),
      iter = 15000,
      sel = seq(10000, 15000, by = 20),
      use_knots = TRUE,
      n_knots = 75,
      K = 1500
    )
  } else {
    # Standard settings
    settings <- list(
      n = min(500, nrow(sim_popn)),
      iter = 12000,
      sel = seq(8000, 12# g-Bayesian Kernel Machine Regression (g-BKMR) Package
```
In this document, we illustrate the main features of the `causalGBKMR` R package through examples. This approach enables causal inference for health effects of time-varying correlated environmental mixtures while accounting for time-varying confounding.

## Cite the method

The g-BKMR method is introduced in the paper:

Chai, Z., Navas-Acien, A., Coull, B., & Valeri, L. (2024). g-BKMR: Causal Inference for health effects of time-varying correlated environmental mixtures. *Journal of Statistical Methods*.

## g-BKMR Method Overview

g-BKMR extends Bayesian Kernel Machine Regression (BKMR) to longitudinal settings with time-varying exposures and confounders. The method combines:

- **Flexible modeling**: Non-linear, non-additive effects of environmental mixtures
- **Causal inference**: G-formula approach for time-varying confounding
- **Variable selection**: Bayesian variable selection for high-dimensional exposures

The general framework models:
- Time-varying confounders: $L_t = h_L(A_{t-1}, L_{t-1}, C_0) + \epsilon_{L_t}$
- Outcome: $Y = h_Y(A_0, ..., A_{T-1}, L_1, ..., L_{T-1}, C_0) + \epsilon_Y$

where $A_t$ represents exposures at time $t$, $L_t$ are time-varying confounders, and $C_0$ are baseline covariates.

## Installation

```r
# Install required packages
install.packages(c("bkmr", "dplyr", "ggplot2", "parallel", 
                   "fields", "lubridate", "data.table", "corrplot"))

# Install causalGBKMR package (development version)
devtools::install_github("your-username/causalGBKMR")
library(causalGBKMR)
```

## Data Preparation from User Matrices

Most users have data in matrix format. The package provides `prepare_gbkmr_data()` to convert user matrices into the required g-BKMR format.

### Required Input

**Three Matrices:**
1. **`Y`**: Outcome vector (length n)
2. **`Z`**: Mixture exposure matrix (n × (Adim × T)) 
3. **`X`**: Covariate matrix (n × (Ldim × T + baseline_covs))

**Three Parameters:**
1. **`time_points`**: Number of time points (T)
2. **`mixture_components`**: Number of mixture components per time point (Adim)
3. **`td_covariates`**: Number of time-dependent covariates per time point (Ldim)

### Critical Matrix Organization

**Z Matrix (Mixtures)**: Chronological order by time point
```r
# Example: 3 metals × 4 time points = 12 columns
# [As_T0, Mn_T0, Pb_T0, As_T1, Mn_T1, Pb_T1, As_T2, Mn_T2, Pb_T2, As_T3, Mn_T3, Pb_T3]
```

**X Matrix (Covariates)**: Time-dependent covariates + baseline covariates  
```r
# Example: 2 TD covariates × 4 time points + 2 baseline = 10 columns
# [BMI_T1, BP_T1, BMI_T2, BP_T2, BMI_T3, BP_T3, BMI_T4, BP_T4, Sex, Age]
# Note: TD covariates start from T1 (not T0), baseline covariates at end
```

### Smart Variable Naming

The package now uses **intelligent naming** for time-dependent covariates:

**When you provide covariate names:**
```r
prepared_data <- prepare_gbkmr_data(
  Y = Y, Z = Z, X = X,
  time_points = 4, mixture_components = 3, td_covariates = 2,
  td_covariate_names = c("bmi", "blood_pressure")  # Your actual names
)

# Creates variables: bmi_0, blood_pressure_0, bmi_1, blood_pressure_1, bmi_2, blood_pressure_2, etc.
```

**When you don't provide names (generic):**
```r
prepared_data <- prepare_gbkmr_data(
  Y = Y, Z = Z, X = X,
  time_points = 4, mixture_components = 3, td_covariates = 2
  # No td_covariate_names provided
)

# Creates variables: td_covariate1_0, td_covariate2_0, td_covariate1_1, td_covariate2_1, etc.
```

### Usage Example

```r
library(causalGBKMR)

# Example 1: With custom covariate names (recommended)
prepared_data <- prepare_gbkmr_data(
  Y = your_outcome_vector,              # Your outcome data
  Z = your_mixture_matrix,              # Your mixture exposure data  
  X = your_covariate_matrix,            # Your covariate data
  time_points = 4,                      # 4 time points
  mixture_components = 3,               # 3 metals per time point  
  td_covariates = 2,                    # 2 time-dependent covariates per time point
  baseline_covariates = 2,              # 2 baseline covariates
  td_covariate_names = c("bmi", "bp"),  # Meaningful names for your covariates
  log_transform_mixtures = TRUE         # Log-transform exposures (recommended)
)

# Example 2: Without custom names (uses generic names)
prepared_data_generic <- prepare_gbkmr_data(
  Y = your_outcome_vector,
  Z = your_mixture_matrix,
  X = your_covariate_matrix,
  time_points = 4,
  mixture_components = 3,
  td_covariates = 2,
  baseline_covariates = 2
  # No td_covariate_names - will use td_covariate1, td_covariate2
)

# Validate the converted data
check_gbkmr_data(prepared_data)

# Run g-BKMR analysis
results <- gbkmr_run(
  data = prepared_data,
  outcome = "Y",
  outcome_type = "continuous",
  time_points = 4
)
```

### Key Constraints

- **Time-dependent covariates MUST be measured at ALL time points**
- **Matrix dimensions must match exactly** (function validates this)
- **Baseline covariates go at the END of X matrix**
- **Missing time points not supported** in basic version
```

## Tutorial

### Option 1: Using Your Own Data (Recommended)

For users with real data in matrix format:

### Basic Usage (Novice Mode)

For users who want a simple analysis, the package now **automatically detects** variable patterns in your data:

```r
library(causalGBKMR)

# Step 1: Prepare your data (if from matrices)
prepared_data <- prepare_gbkmr_data(
  Y = your_outcome_vector,
  Z = your_mixture_matrix,
  X = your_covariate_matrix,
  time_points = 4,
  mixture_components = 3,
  td_covariates = 2,
  baseline_covariates = 2
)

# Step 2: Run g-BKMR analysis (simplified interface)
results <- gbkmr_run(
  data = prepared_data,           # Your prepared dataset
  outcome = "Y",                  # Outcome variable name
  outcome_type = "continuous",    # "continuous" or "binary"
  time_points = 4                 # Number of time points
)

# The function automatically detects:
# - Number of exposures per time point
# - Number of time-dependent covariates per time point  
# - Time-dependent covariate naming patterns
# - All variable relationships across time points

print(results)
```

**What you need to provide (minimal input):**
1. `data`: Your dataset (prepared with `prepare_gbkmr_data()` or in sim_popn format)
2. `outcome`: Name of your outcome variable
3. `outcome_type`: Whether outcome is "continuous" or "binary"
4. `time_points`: Number of time points in your study

**What the package does automatically:**
- **Auto-detects exposure variables** (logM1_0, logM2_0, etc.) across all time points
- **Auto-detects time-dependent covariates** with flexible naming patterns
- **Auto-detects baseline covariates** and their patterns
- **Fits sequential BKMR models** for each time-dependent confounder
- **Handles g-formula computation** automatically
- **Returns causal effect estimates** with credible intervals

### Advanced Usage (Statistician Mode)

For users who need more control over the modeling process:

```r
# Generate data with different dimensions for advanced example
sim_popn_advanced <- generate_panel_data(
  popN = 1e5,
  T = 4,                      # 4 time points
  Adim = 4,                   # 4 exposures per time point
  Ldim = 2,                   # 2 confounders per time point
  outcome_type = "continuous",
  relationship_type = "quadratic+interaction",
  confounding = "high"
)

# Advanced analysis with custom parameters
results_advanced <- run_gbkmr_panel(
  sim_popn = sim_popn_advanced,
  T = 4,
  currind = 1,
  sel = seq(20000, 25000, by = 30),   # Different sampling scheme
  n = 800,                            # Larger sample size
  K = 2000,                           # More prediction samples
  iter = 30000,                       # More MCMC iterations
  parallel = TRUE,
  save_exposure_preds = TRUE,
  return_ci = TRUE,
  make_plots = TRUE,
  use_knots = TRUE,
  n_knots = 100                       # More knots for complex relationships
)
```

**Advanced Parameters:**
- `T`: Number of time points (must match data generation)
- `sel`: Custom MCMC sampling scheme for inference
- `n`: Analysis sample size (allows subsampling from large populations)
- `K`: Number of prediction samples for g-formula computation
- `iter`: Total MCMC iterations (adjust based on convergence needs)
- `n_knots`: Number of knots (increase for complex exposure-response relationships)

### Flexible Variable Detection

The updated package now **automatically detects** various time-dependent covariate naming patterns:

```r
# Supported naming patterns for time-dependent covariates:

# Pattern 1: User-provided meaningful names
# bmi_0, bp_0 (baseline)
# bmi_1, bp_1 (time 1) 
# bmi_2, bp_2 (time 2)

# Pattern 2: Alternative user format
# bmi0, bp0 (baseline)
# bmi1, bp1 (time 1)
# bmi2, bp2 (time 2)

# Pattern 3: Generated data format (legacy)
# waist0_1, waist0_2 (baseline)
# waist1_1, waist1_2 (time 1)  
# waist2_1, waist2_2 (time 2)

# Pattern 4: Generic format (when names not provided)
# td_covariate1_0, td_covariate2_0 (baseline)
# td_covariate1_1, td_covariate2_1 (time 1)
# td_covariate1_2, td_covariate2_2 (time 2)

# The detection function automatically identifies the pattern used
results <- gbkmr_run(
  data = your_data_with_any_naming,
  outcome = "Y", 
  outcome_type = "continuous",
  time_points = 4
)

# Check what was detected
cat("Auto-detected TD covariate names:\n")
print(results$detection_info$td_covariate_names)
cat("Detection pattern:", results$detection_info$detected_pattern, "\n")
```

## Results Interpretation

### Main Causal Effect

```r
# Extract main causal effect estimate
causal_effect <- results_basic$causal_effect
cat("Causal Effect (75th vs 25th percentile):", causal_effect$estimate, "\n")
cat("95% Credible Interval:", causal_effect$lower, "to", causal_effect$upper, "\n")

# Access detailed results
summary(results_basic)
```

### Visualization

```r
library(ggplot2)

# Plot main results
plot(results_basic)

# Or create custom plots
# Plot counterfactual means comparison
df_plot <- data.frame(
  Scenario = c("25th percentile", "75th percentile"), 
  Mean = c(results_basic$counterfactual_means$low, 
           results_basic$counterfactual_means$high)
)

ggplot(df_plot, aes(x = Scenario, y = Mean, fill = Scenario)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Counterfactual Outcome Comparison", 
       y = "Mean Outcome")

# Plot variable importance if available
if (!is.null(results_basic$variable_importance)) {
  barplot(results_basic$variable_importance, 
          main = "Variable Importance", 
          las = 2, cex.names = 0.8)
}
```

### Accessing Detailed Results

```r
# Examine results structure
str(results_basic, max.level = 1)

# Access fitted models (if saved)
if (!is.null(results_basic$fitted_models)) {
  outcome_model <- results_basic$fitted_models$outcome
  confounder_models <- results_basic$fitted_models$confounders
}

# Access variable selection results (if enabled)
if (!is.null(results_basic$variable_selection)) {
  inclusion_probs <- results_basic$variable_selection$inclusion_probabilities
  print("Variable inclusion probabilities:")
  print(inclusion_probs)
}

# Access convergence diagnostics (if enabled)
if (!is.null(results_basic$diagnostics)) {
  print("MCMC Convergence Summary:")
  print(results_basic$diagnostics$convergence_summary)
}
```

## Model Comparison and Sensitivity Analysis

### Comparing Different Data Generation Scenarios

```r
# Function to compare different relationship types
compare_relationships <- function() {
  relationships <- c("linear", "quadratic", "quadratic+interaction")
  results <- list()
  
  for (rel in relationships) {
    # Generate data
    sim_data <- generate_panel_data(
      popN = 5e4, T = 3, Adim = 3, Ldim = 2,
      outcome_type = "continuous",
      relationship_type = rel,
      confounding = "high"
    )
    
    # Fit model
    fit <- run_gbkmr_panel(
      sim_popn = sim_data, T = 3, currind = 1,
      sel = seq(8000, 10000, by = 25),
      n = 300, iter = 12000, parallel = TRUE,
      make_plots = FALSE
    )
    
    results[[rel]] <- fit$diff_gBKMR
  }
  
  return(results)
}

# Run comparison
relationship_results <- compare_relationships()

# Plot comparison
df_comparison <- data.frame(
  Relationship = names(relationship_results),
  Effect = unlist(relationship_results)
)

ggplot(df_comparison, aes(x = Relationship, y = Effect, fill = Relationship)) +
  geom_col() +
  theme_minimal() +
  labs(title = "g-BKMR Effects Across Different Functional Forms",
       y = "Causal Effect Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### Sensitivity to Sample Size and Dimensions

```r
# Function for comprehensive sensitivity analysis
sensitivity_analysis <- function() {
  # Test different sample sizes
  sample_sizes <- c(200, 500, 1000)
  
  # Test different exposure dimensions
  exposure_dims <- c(2, 3, 5)
  
  # Test different confounder dimensions
  confounder_dims <- c(1, 2, 3)
  
  results <- expand.grid(
    n = sample_sizes,
    Adim = exposure_dims,
    Ldim = confounder_dims
  )
  results$effect <- NA
  
  for (i in 1:nrow(results)) {
    # Generate data
    sim_data <- generate_panel_data(
      popN = max(results$n) * 2,  # Ensure enough population
      T = 3,
      Adim = results$Adim[i],
      Ldim = results$Ldim[i],
      outcome_type = "continuous",
      relationship_type = "quadratic",
      confounding = "high"
    )
    
    # Fit model
    fit <- run_gbkmr_panel(
      sim_popn = sim_data, T = 3, currind = i,
      sel = seq(5000, 7000, by = 25),
      n = results$n[i], iter = 8000, parallel = TRUE,
      make_plots = FALSE
    )
    
    results$effect[i] <- fit$diff_gBKMR
  }
  
  return(results)
}

# Run sensitivity analysis (this may take a while)
# sens_results <- sensitivity_analysis()

# Example of plotting sensitivity results
# ggplot(sens_results, aes(x = n, y = effect, color = factor(Adim))) +
#   geom_line() + geom_point() +
#   facet_wrap(~Ldim, labeller = label_both) +
#   theme_minimal() +
#   labs(title = "Sensitivity Analysis",
#        x = "Sample Size", y = "Causal Effect",
#        color = "Exposure Dimension")
```

## Data Preparation for Real Studies

When working with your own data, you'll need to format it to match the expected structure. The package expects wide-format data with specific naming conventions.

### Data Formatting Guidelines

```r
# Function to help format real data
format_real_data <- function(your_data, T, Adim, Ldim) {
  # Example transformation for longitudinal environmental health data
  
  # Ensure exposure variables are log-transformed and named correctly
  # Expected names: logM1_0, logM2_0, ..., logMAdim_0, logM1_1, ..., logMAdim_T-1
  
  # Ensure confounder variables are named correctly  
  # Expected names: waist1_1, ..., waistLdim_1, ..., waist1_T-1, ..., waistLdim_T-1
  
  # Ensure baseline variables are named correctly
  # Expected names: sex, waist0_1, ..., waist0_Ldim
  
  # Ensure outcome variable is named Y
  
  return(formatted_data)
}

# Validate data structure
validate_gbkmr_data <- function(data, T, Adim, Ldim) {
  required_cols <- c(
    "sex", "Y", "id",
    paste0("waist0_", 1:Ldim),
    unlist(lapply(0:(T-1), function(t) paste0("logM", 1:Adim, "_", t))),
    if (T > 1) unlist(lapply(1:(T-1), function(t) paste0("waist", t, "_", 1:Ldim)))
  )
  
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    cat("Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    return(FALSE)
  }
  
  cat("Data validation passed!\n")
  return(TRUE)
}
```

## Performance Optimization

### For Large Datasets

```r
# Optimize for computational efficiency
results_optimized <- gbkmr_run(
  data = large_dataset,
  outcome = "y_binary",
  outcome_type = "binary",
  exposure = "arsenic",
  t_covariates = c("bp", "bmi"),
  base_covariates = c("age", "sex"),
  time_points = 4,
  
  computation = list(
    parallel = TRUE,
    n_cores = 8,              # Use more cores
    use_knots = TRUE,         # Essential for large datasets
    n_knots = 100,            # Increase knots for better approximation
    chunk_size = 1000         # Process data in chunks
  ),
  
  mcmc_control = list(
    iter = 10000,             # Reduce iterations for speed
    burn = 2000,
    thin = 2                  # Use thinning to reduce memory
  )
)
```

### Memory Management

```r
# For memory-constrained environments
results_memory_efficient <- gbkmr_run(
  data = example_data,
  outcome = "y_binary",
  outcome_type = "binary",
  exposure = "arsenic",
  t_covariates = c("bp", "bmi"),
  base_covariates = c("age", "sex"),
  time_points = 4,
  
  memory_control = list(
    save_posterior = FALSE,    # Don't save all posterior samples
    save_models = FALSE,       # Don't save intermediate models
    stream_output = TRUE       # Stream results to disk
  )
)
```

## Real Data Example

Here's an example using a realistic environmental health dataset:

```r
# Load example environmental health data
data("environmental_cohort")  # Example dataset included in package

# Examine data structure
str(environmental_cohort)
head(environmental_cohort)

# Check for missing data
check_missing_data(environmental_cohort)

# Fit g-BKMR model
environmental_results <- gbkmr_run(
  data = environmental_cohort,
  outcome = "cardiovascular_risk",
  outcome_type = "binary",
  exposure = "pm25",                    # Air pollution exposure
  t_covariates = c("blood_pressure", "bmi", "stress_score"),
  base_covariates = c("age", "sex", "education", "smoking_status"),
  time_points = 5,
  id = "participant_id",
  time = "visit_number",
  
  mcmc_control = list(
    iter = 25000,
    burn = 5000
  ),
  
  variable_selection = TRUE,
  
  diagnostics = list(
    trace_plots = TRUE,
    convergence_check = TRUE
  )
)

# Interpret results
summary(environmental_results)
plot(environmental_results, type = "all")

# Export results for publication
export_results(environmental_results, 
               format = "table", 
               file = "gbkmr_results.csv")
```

## Best Practices

### Data Quality

1. **Missing Data**: Handle missing values before analysis
2. **Outliers**: Check for and address extreme values
3. **Temporal Consistency**: Ensure time points are properly ordered
4. **Variable Scaling**: Consider standardizing exposures and covariates

### Model Selection

1. **Start Simple**: Begin with vanilla kernel configuration
2. **Progressive Complexity**: Add variables to kernel as needed
3. **Cross-Validation**: Use cross-validation for model comparison
4. **Domain Knowledge**: Incorporate subject matter expertise

### Computational Considerations

1. **Parallel Processing**: Always enable for datasets with >100 subjects
2. **Knot Selection**: Use more knots for larger datasets
3. **MCMC Tuning**: Monitor convergence and adjust iterations accordingly
4. **Memory Planning**: Consider memory requirements for large datasets

## Troubleshooting

### Common Issues

1. **Convergence Problems**:
   ```r
   # Increase iterations and check trace plots
   results <- gbkmr_run(..., 
                       mcmc_control = list(iter = 50000, burn = 10000),
                       diagnostics = list(trace_plots = TRUE))
   ```

2. **Memory Issues**:
   ```r
   # Use memory-efficient settings
   results <- gbkmr_run(...,
                       computation = list(use_knots = TRUE, n_knots = 50),
                       memory_control = list(save_posterior = FALSE))
   ```

3. **Long Computing Times**:
   ```r
   # Optimize computational settings
   results <- gbkmr_run(...,
                       computation = list(parallel = TRUE, n_cores = 8),
                       mcmc_control = list(iter = 10000))
   ```

## Limitations and Considerations

- **Computational Complexity**: g-BKMR is computationally intensive for large numbers of time points and subjects
- **Model Assumptions**: Relies on standard causal inference assumptions (consistency, exchangeability, positivity)
- **Data Requirements**: Requires complete time series for all subjects
- **Variable Selection**: May be conservative in high-dimensional settings

## Conclusion

The g-BKMR package provides a user-friendly interface for causal inference with time-varying environmental mixtures. The progressive complexity design accommodates both novice users seeking straightforward analysis and expert statisticians requiring detailed model customization. The package's automation of complex longitudinal modeling while maintaining flexibility makes it valuable for environmental health research.
