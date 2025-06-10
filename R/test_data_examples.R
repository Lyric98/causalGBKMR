# ===============================================================================
# Test Data Generation for g-BKMR Package
# ===============================================================================

# ===============================================================================
# Example 1: Small Test Dataset (Quick Testing)
# ===============================================================================

create_test_data_small <- function() {
  cat("Creating small test dataset for quick testing...\n")

  set.seed(42)
  n <- 100          # Small sample for quick testing
  T <- 3            # 3 time points
  Adim <- 2         # 2 metals (Arsenic, Lead)
  Ldim <- 2         # 2 time-dependent covariates (BMI, Blood Pressure)
  n_baseline <- 2   # 2 baseline covariates (Sex, Age)

  # Generate outcome (cardiovascular risk score)
  Y <- rnorm(n, mean = 50, sd = 10)

  # Generate mixture matrix Z (metals with realistic concentrations)
  # Columns: [As_T0, Pb_T0, As_T1, Pb_T1, As_T2, Pb_T2]
  Z <- matrix(
    c(
      rlnorm(n, meanlog = log(5), sdlog = 0.5),    # Arsenic T0 (μg/L)
      rlnorm(n, meanlog = log(10), sdlog = 0.6),   # Lead T0 (μg/L)
      rlnorm(n, meanlog = log(5.5), sdlog = 0.5),  # Arsenic T1
      rlnorm(n, meanlog = log(11), sdlog = 0.6),   # Lead T1
      rlnorm(n, meanlog = log(6), sdlog = 0.5),    # Arsenic T2
      rlnorm(n, meanlog = log(12), sdlog = 0.6)    # Lead T2
    ),
    nrow = n, ncol = Adim * T
  )

  # Generate covariate matrix X
  # Structure: [BMI_T1, BP_T1, BMI_T2, BP_T2, BMI_T3, BP_T3, Sex, Age]
  # Note: TD covariates at all T time points + baseline covariates
  bmi_base <- rnorm(n, mean = 25, sd = 4)
  bp_base <- rnorm(n, mean = 120, sd = 15)

  X <- cbind(
    # BMI at all time points (with trend)
    bmi_base + rnorm(n, 0, 1),                    # BMI T1
    bp_base + rnorm(n, 0, 5),                     # BP T1
    bmi_base + rnorm(n, 0.5, 1),                  # BMI T2
    bp_base + rnorm(n, 2, 5),                     # BP T2
    bmi_base + rnorm(n, 1, 1),                    # BMI T3
    bp_base + rnorm(n, 3, 5),                     # BP T3
    # Baseline covariates
    rbinom(n, 1, 0.5),                            # Sex (0/1)
    rnorm(n, mean = 45, sd = 12)                  # Age
  )

  cat("Small test data dimensions:\n")
  cat("Y: length", length(Y), "\n")
  cat("Z:", nrow(Z), "×", ncol(Z), "(", Adim, "metals ×", T, "time points)\n")
  cat("X:", nrow(X), "×", ncol(X), "(", Ldim, "TD covs ×", T, "time points +", n_baseline, "baseline)\n\n")

  return(list(Y = Y, Z = Z, X = X,
              n = n, T = T, Adim = Adim, Ldim = Ldim, n_baseline = n_baseline,
              description = "Small dataset: 100 subjects, 3 time points, 2 metals, 2 TD covariates"))
}

# ===============================================================================
# Example 2: Medium Test Dataset (Realistic Size)
# ===============================================================================

create_test_data_medium <- function() {
  cat("Creating medium test dataset for realistic testing...\n")

  set.seed(123)
  n <- 300          # Medium sample size
  T <- 4            # 4 time points
  Adim <- 3         # 3 metals (Arsenic, Lead, Cadmium)
  Ldim <- 2         # 2 time-dependent covariates (BMI, SBP)
  n_baseline <- 3   # 3 baseline covariates (Sex, Age, Education)

  # Generate outcome with some realistic relationships
  metal_effect <- rnorm(n, 0, 2)  # Metal mixture effect
  covariate_effect <- rnorm(n, 0, 1)  # Covariate effect
  Y <- 50 + metal_effect + covariate_effect + rnorm(n, 0, 5)

  # Generate mixture matrix Z with temporal correlation
  # Metals: Arsenic, Lead, Cadmium at 4 time points
  baseline_metals <- matrix(
    c(
      rlnorm(n, meanlog = log(3), sdlog = 0.8),   # Arsenic baseline
      rlnorm(n, meanlog = log(8), sdlog = 0.7),   # Lead baseline
      rlnorm(n, meanlog = log(0.5), sdlog = 0.9)  # Cadmium baseline
    ),
    nrow = n, ncol = Adim
  )

  Z <- matrix(0, nrow = n, ncol = Adim * T)
  for (t in 0:(T-1)) {
    for (a in 1:Adim) {
      col_idx <- t * Adim + a
      # Add temporal correlation and trend
      Z[, col_idx] <- baseline_metals[, a] * exp(rnorm(n, mean = t * 0.1, sd = 0.3))
    }
  }

  # Generate covariate matrix X with realistic health trajectories
  # BMI and SBP with age-related trends
  age <- rnorm(n, mean = 50, sd = 15)
  education <- sample(1:4, n, replace = TRUE)  # 1=< HS, 2=HS, 3=College, 4=Grad
  sex <- rbinom(n, 1, 0.6)  # 60% female

  # BMI trajectory (tends to increase with age and time)
  bmi_baseline <- 22 + 0.1 * age + rnorm(n, 0, 3)
  bmi_trajectory <- matrix(0, nrow = n, ncol = T)
  for (t in 1:T) {
    bmi_trajectory[, t] <- bmi_baseline + t * 0.2 + rnorm(n, 0, 0.5)
  }

  # SBP trajectory (increases with age and BMI)
  sbp_baseline <- 100 + 0.8 * age + 0.5 * bmi_baseline + rnorm(n, 0, 10)
  sbp_trajectory <- matrix(0, nrow = n, ncol = T)
  for (t in 1:T) {
    sbp_trajectory[, t] <- sbp_baseline + t * 1.5 + 0.3 * bmi_trajectory[, t] + rnorm(n, 0, 3)
  }

  # Construct X matrix: [BMI_T1, SBP_T1, BMI_T2, SBP_T2, ..., BMI_T4, SBP_T4, Sex, Age, Education]
  X <- cbind(
    matrix(c(bmi_trajectory, sbp_trajectory), nrow = n, ncol = Ldim * T),
    sex, age, education
  )

  cat("Medium test data dimensions:\n")
  cat("Y: length", length(Y), "\n")
  cat("Z:", nrow(Z), "×", ncol(Z), "(", Adim, "metals ×", T, "time points)\n")
  cat("X:", nrow(X), "×", ncol(X), "(", Ldim, "TD covs ×", T, "time points +", n_baseline, "baseline)\n\n")

  return(list(Y = Y, Z = Z, X = X,
              n = n, T = T, Adim = Adim, Ldim = Ldim, n_baseline = n_baseline,
              description = "Medium dataset: 300 subjects, 4 time points, 3 metals, 2 TD covariates"))
}

# ===============================================================================
# Example 3: Complex Test Dataset (High Dimensional)
# ===============================================================================

create_test_data_complex <- function() {
  cat("Creating complex test dataset for comprehensive testing...\n")

  set.seed(456)
  n <- 500          # Larger sample
  T <- 5            # 5 time points
  Adim <- 4         # 4 metals (As, Pb, Cd, Hg)
  Ldim <- 3         # 3 time-dependent covariates (BMI, SBP, Glucose)
  n_baseline <- 4   # 4 baseline covariates (Sex, Age, Education, Smoking)

  # Generate complex outcome with interactions
  Y <- rnorm(n, mean = 100, sd = 20)  # More complex health outcome

  # Generate highly correlated metal mixtures
  # Simulate environmental co-contamination
  contamination_level <- rnorm(n, 0, 1)  # Individual contamination exposure

  Z <- matrix(0, nrow = n, ncol = Adim * T)
  for (t in 0:(T-1)) {
    for (a in 1:Adim) {
      col_idx <- t * Adim + a
      base_level <- c(2, 5, 0.3, 0.1)[a]  # Typical levels for As, Pb, Cd, Hg

      # Correlated contamination with temporal trends
      Z[, col_idx] <- exp(
        log(base_level) +
          0.8 * contamination_level +  # Common contamination factor
          rnorm(n, mean = t * 0.05, sd = 0.4)  # Temporal variation
      )
    }
  }

  # Generate complex time-dependent covariates
  # Baseline characteristics
  age <- pmax(18, pmin(80, rnorm(n, mean = 55, sd = 15)))
  sex <- rbinom(n, 1, 0.55)
  education <- sample(1:5, n, replace = TRUE, prob = c(0.1, 0.2, 0.3, 0.3, 0.1))
  smoking <- rbinom(n, 1, 0.25)

  # Complex health trajectories
  X_td <- array(0, dim = c(n, Ldim, T))

  # BMI trajectory
  bmi_0 <- 20 + 0.15 * age + rnorm(n, 0, 4)
  for (t in 1:T) {
    X_td[, 1, t] <- bmi_0 + (t-1) * 0.3 + rnorm(n, 0, 0.8)
  }

  # SBP trajectory (depends on BMI and age)
  sbp_0 <- 90 + 0.9 * age + rnorm(n, 0, 12)
  for (t in 1:T) {
    X_td[, 2, t] <- sbp_0 + (t-1) * 2 + 0.8 * X_td[, 1, t] + rnorm(n, 0, 5)
  }

  # Glucose trajectory (metabolic changes)
  glucose_0 <- 85 + 0.2 * age + 0.5 * bmi_0 + rnorm(n, 0, 8)
  for (t in 1:T) {
    X_td[, 3, t] <- glucose_0 + (t-1) * 1.5 + 0.3 * X_td[, 1, t] + rnorm(n, 0, 4)
  }

  # Reshape to matrix format
  X <- cbind(
    matrix(X_td, nrow = n, ncol = Ldim * T),
    sex, age, education, smoking
  )

  cat("Complex test data dimensions:\n")
  cat("Y: length", length(Y), "\n")
  cat("Z:", nrow(Z), "×", ncol(Z), "(", Adim, "metals ×", T, "time points)\n")
  cat("X:", nrow(X), "×", ncol(X), "(", Ldim, "TD covs ×", T, "time points +", n_baseline, "baseline)\n\n")

  return(list(Y = Y, Z = Z, X = X,
              n = n, T = T, Adim = Adim, Ldim = Ldim, n_baseline = n_baseline,
              description = "Complex dataset: 500 subjects, 5 time points, 4 metals, 3 TD covariates"))
}

# ===============================================================================
# Complete Test Function
# ===============================================================================

test_gbkmr_workflow <- function(dataset_type = "small") {
  cat("=== Testing g-BKMR Workflow ===\n\n")

  # Generate test data
  if (dataset_type == "small") {
    test_data <- create_test_data_small()
  } else if (dataset_type == "medium") {
    test_data <- create_test_data_medium()
  } else if (dataset_type == "complex") {
    test_data <- create_test_data_complex()
  } else {
    stop("dataset_type must be 'small', 'medium', or 'complex'")
  }

  cat("Generated test data:", test_data$description, "\n\n")

  # Step 1: Prepare data for g-BKMR
  cat("Step 1: Preparing data with meaningful covariate names...\n")

  covariate_names <- switch(test_data$Ldim,
                            "1" = "bmi",
                            "2" = c("bmi", "systolic_bp"),
                            "3" = c("bmi", "systolic_bp", "glucose")
  )

  tryCatch({
    prepared_data <- prepare_gbkmr_data(
      Y = test_data$Y,
      Z = test_data$Z,
      X = test_data$X,
      time_points = test_data$T,
      mixture_components = test_data$Adim,
      td_covariates = test_data$Ldim,
      baseline_covariates = test_data$n_baseline,
      td_covariate_names = covariate_names,
      log_transform_mixtures = TRUE
    )

    cat("✓ Data preparation successful!\n")
    cat("✓ Generated data frame:", nrow(prepared_data), "×", ncol(prepared_data), "\n")
    cat("✓ Column names:", paste(colnames(prepared_data)[1:min(10, ncol(prepared_data))], collapse = ", "))
    if (ncol(prepared_data) > 10) cat(", ...")
    cat("\n\n")

  }, error = function(e) {
    cat("✗ Data preparation failed:", e$message, "\n")
    return(NULL)
  })

  # Step 2: Run g-BKMR analysis (with reduced iterations for testing)
  cat("Step 2: Running g-BKMR analysis...\n")

  tryCatch({
    results <- gbkmr_run(
      data = prepared_data,
      outcome = "Y",
      outcome_type = "continuous",
      time_points = test_data$T,
      # Reduced parameters for faster testing
      iter = 1000,
      n = min(100, nrow(prepared_data)),
      use_knots = TRUE,
      n_knots = 20,
      verbose = TRUE
    )

    cat("\n✓ g-BKMR analysis completed successfully!\n")
    cat("✓ Causal effect estimate:", round(results$causal_effect$estimate, 4), "\n")
    cat("✓ 95% CI: (", round(results$causal_effect$lower, 4), ",",
        round(results$causal_effect$upper, 4), ")\n")
    cat("✓ Detection info available: p =", results$detection_info$p,
        ", Ldim =", results$detection_info$Ldim, "\n")
    cat("✓ TD covariate names detected:", paste(results$detection_info$td_covariate_names, collapse = ", "), "\n\n")

    return(list(test_data = test_data, prepared_data = prepared_data, results = results))

  }, error = function(e) {
    cat("✗ g-BKMR analysis failed:", e$message, "\n")
    return(list(test_data = test_data, prepared_data = prepared_data, results = NULL))
  })
}

# ===============================================================================
# Quick Test Functions
# ===============================================================================

# Quick test for small dataset
test_small <- function() {
  cat("Running quick test with small dataset...\n")
  return(test_gbkmr_workflow("small"))
}

# Quick test for medium dataset
test_medium <- function() {
  cat("Running test with medium dataset...\n")
  return(test_gbkmr_workflow("medium"))
}

# Quick test for complex dataset
test_complex <- function() {
  cat("Running test with complex dataset...\n")
  return(test_gbkmr_workflow("complex"))
}

# ===============================================================================
# Example Usage Instructions
# ===============================================================================

print_usage_instructions <- function() {
  cat("=== How to Test g-BKMR Package ===\n\n")
  cat("1. Quick Test (recommended to start):\n")
  cat("   result_small <- test_small()\n\n")

  cat("2. Medium Test (more realistic):\n")
  cat("   result_medium <- test_medium()\n\n")

  cat("3. Complex Test (full functionality):\n")
  cat("   result_complex <- test_complex()\n\n")

  cat("4. Manual Testing:\n")
  cat("   # Generate data\n")
  cat("   test_data <- create_test_data_small()  # or medium, complex\n")
  cat("   \n")
  cat("   # Prepare for g-BKMR\n")
  cat("   prepared <- prepare_gbkmr_data(\n")
  cat("     Y = test_data$Y, Z = test_data$Z, X = test_data$X,\n")
  cat("     time_points = test_data$T, mixture_components = test_data$Adim,\n")
  cat("     td_covariates = test_data$Ldim, baseline_covariates = test_data$n_baseline,\n")
  cat("     td_covariate_names = c('bmi', 'bp')  # customize as needed\n")
  cat("   )\n")
  cat("   \n")
  cat("   # Run analysis\n")
  cat("   results <- gbkmr_run(prepared, outcome_type = 'continuous', time_points = test_data$T)\n\n")

  cat("5. Inspect Results:\n")
  cat("   print(results)\n")
  cat("   summary(results)\n")
  cat("   results$detection_info  # See what was auto-detected\n\n")
}

# Print usage instructions when this file is sourced
print_usage_instructions()
