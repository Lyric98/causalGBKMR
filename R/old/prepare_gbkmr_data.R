prepare_gbkmr_data <- function(
    Y,                     # Outcome vector (length n)
    Z,                     # Mixture exposure matrix (n × (Adim * T))
    X,                     # Covariate matrix (n × (Ldim * T + baseline_covs))
    time_points,           # Number of time points (T)
    mixture_components,    # Number of mixture components per time point (Adim)
    td_covariates,         # Number of time-dependent covariates per time point (Ldim)
    baseline_covariates = 1,  # Number of baseline covariates (default: sex only)
    log_transform_mixtures = TRUE,   # Whether to log-transform mixtures
    validate_input = TRUE            # Whether to validate input dimensions
) {

  # Get dimensions
  n <- length(Y)
  T <- time_points
  Adim <- mixture_components
  Ldim <- td_covariates
  n_baseline <- baseline_covariates

  # Input validation
  if (validate_input) {
    validate_user_matrices(Y, Z, X, T, Adim, Ldim, n_baseline)
  }

  cat("Converting user matrices to g-BKMR format...\n")
  cat("Data dimensions: n =", n, ", T =", T, ", Adim =", Adim, ", Ldim =", Ldim, "\n")

  # Initialize the data frame with required columns in sim_popn order
  df <- data.frame(id = 1:n)

  # Add baseline covariates
  # Extract baseline covariates from the END of X matrix
  if (n_baseline > 0) {
    baseline_start_col <- ncol(X) - n_baseline + 1
    baseline_data <- X[, baseline_start_col:ncol(X), drop = FALSE]

    # First baseline covariate is always 'sex' (required by sim_popn format)
    df$sex <- baseline_data[, 1]

    # Add additional baseline covariates if any
    if (n_baseline > 1) {
      for (i in 2:n_baseline) {
        df[[paste0("baseline_", i)]] <- baseline_data[, i]
      }
    }
  } else {
    # Create default sex variable if no baseline covariates
    df$sex <- rbinom(n, 1, 0.5)
  }

  # Add baseline time-dependent covariates (waist0_1, waist0_2, ..., waist0_Ldim)
  # These represent initial values - use first time point or generate
  for (l in 1:Ldim) {
    var_name <- paste0("waist0_", l)

    # Extract first time point of this time-dependent covariate from X
    first_timepoint_col <- l  # First Ldim columns are T=1 values
    if (first_timepoint_col <= (ncol(X) - n_baseline)) {
      # Use first time point value with small noise as baseline
      baseline_val <- X[, first_timepoint_col] + rnorm(n, 0, sd = 0.1 * sd(X[, first_timepoint_col]))
    } else {
      # Generate if not available
      baseline_val <- rnorm(n, mean = 0, sd = 1)
    }

    df[[var_name]] <- baseline_val
  }

  # Add mixture exposures in chronological order (logM1_0, logM2_0, ..., logMAdim_T-1)
  for (t in 0:(T-1)) {
    for (a in 1:Adim) {
      # Calculate column index in Z matrix
      col_idx <- t * Adim + a
      var_name <- paste0("logM", a, "_", t)

      mixture_data <- Z[, col_idx]

      # Log-transform if requested
      if (log_transform_mixtures) {
        if (any(mixture_data <= 0, na.rm = TRUE)) {
          # Handle non-positive values
          min_positive <- min(mixture_data[mixture_data > 0], na.rm = TRUE)
          shift_value <- abs(min(mixture_data, na.rm = TRUE)) + min_positive * 0.001
          mixture_data <- log(mixture_data + shift_value)
          warning("Some mixture values were <= 0. Added constant ", round(shift_value, 6), " before log transformation.")
        } else {
          mixture_data <- log(mixture_data)
        }
      }

      df[[var_name]] <- mixture_data
    }
  }

  # Add time-dependent covariates (waist1_1, waist2_1, ..., waistT-1_Ldim)
  # Note: No time-dependent covariates at T=0 (that's baseline)
  for (t in 1:(T-1)) {
    for (l in 1:Ldim) {
      var_name <- paste0("waist", t, "_", l)

      # Calculate column index in X matrix
      # X matrix structure: [T1_cov1, T1_cov2, ..., T1_covLdim, T2_cov1, T2_cov2, ..., T(T-1)_covLdim, baseline1, baseline2, ...]
      col_idx <- (t - 1) * Ldim + l

      if (col_idx <= (ncol(X) - n_baseline)) {
        df[[var_name]] <- X[, col_idx]
      } else {
        stop("Error in indexing time-dependent covariates. Check matrix dimensions.")
      }
    }
  }

  # Add outcome
  df$Y <- Y

  # Reorder columns to match sim_popn format exactly
  # Expected order: sex, [additional_baseline], waist0_1, ..., waist0_Ldim, logM1_0, logM2_0, ..., logMAdim_T-1, waist1_1, ..., waistT-1_Ldim, Y, id

  # Build the correct column order
  ordered_cols <- "sex"

  # Add additional baseline covariates
  if (n_baseline > 1) {
    ordered_cols <- c(ordered_cols, paste0("baseline_", 2:n_baseline))
  }

  # Add baseline time-dependent covariates
  if (Ldim > 0) {
    ordered_cols <- c(ordered_cols, paste0("waist0_", 1:Ldim))
  }

  # Add mixture exposures in chronological order
  for (t in 0:(T-1)) {
    for (a in 1:Adim) {
      ordered_cols <- c(ordered_cols, paste0("logM", a, "_", t))
    }
  }

  # Add time-dependent covariates in chronological order
  if (T > 1 && Ldim > 0) {
    for (t in 1:(T-1)) {
      for (l in 1:Ldim) {
        ordered_cols <- c(ordered_cols, paste0("waist", t, "_", l))
      }
    }
  }

  # Add outcome and id
  ordered_cols <- c(ordered_cols, "Y", "id")

  # Reorder the dataframe
  df <- df[, ordered_cols]

  # Add metadata
  attr(df, "data_info") <- list(
    n = n,
    time_points = T,
    mixture_components = Adim,
    td_covariates = Ldim,
    baseline_covariates = n_baseline,
    log_transformed = log_transform_mixtures,
    source = "user_matrices"
  )

  # Print summary
  cat("✓ Data conversion successful!\n")
  cat("✓ Generated", nrow(df), "×", ncol(df), "data frame\n")
  cat("✓ Column structure matches sim_popn format\n")
  cat("✓ Ready for g-BKMR analysis\n\n")

  cat("Variable summary:\n")
  cat("- Baseline covariates:", n_baseline, "(including sex)\n")
  cat("- Mixture exposures:", Adim, "components ×", T, "time points =", Adim * T, "variables\n")
  cat("- Time-dependent covariates:", Ldim, "variables ×", (T-1), "time points =", Ldim * (T-1), "variables\n")
  cat("- Total variables:", ncol(df), "\n\n")

  return(df)
}

# Validation function for user matrices
validate_user_matrices <- function(Y, Z, X, T, Adim, Ldim, n_baseline) {

  n <- length(Y)

  # Basic input checks
  if (!is.numeric(Y)) stop("Y must be numeric")
  if (any(is.na(Y))) warning("Y contains missing values")

  if (!is.matrix(Z) && !is.data.frame(Z)) stop("Z must be a matrix or data frame")
  if (nrow(Z) != n) stop("Z must have the same number of rows as length of Y")

  if (!is.matrix(X) && !is.data.frame(X)) stop("X must be a matrix or data frame")
  if (nrow(X) != n) stop("X must have the same number of rows as length of Y")

  # Parameter checks
  if (T < 1) stop("time_points must be >= 1")
  if (Adim < 1) stop("mixture_components must be >= 1")
  if (Ldim < 0) stop("td_covariates must be >= 0")
  if (n_baseline < 0) stop("baseline_covariates must be >= 0")

  # Matrix dimension checks
  expected_Z_cols <- Adim * T
  expected_X_cols <- Ldim * T + n_baseline  # All T time points for TD covariates + baseline

  if (ncol(Z) != expected_Z_cols) {
    stop("Z matrix dimension error!\n",
         "Expected: ", expected_Z_cols, " columns (", Adim, " mixtures × ", T, " time points)\n",
         "Actual: ", ncol(Z), " columns\n",
         "Z should contain mixtures in chronological order: [Mix1_T0, Mix2_T0, ..., MixAdim_T0, Mix1_T1, Mix2_T1, ..., MixAdim_T(T-1)]")
  }

  if (ncol(X) != expected_X_cols) {
    stop("X matrix dimension error!\n",
         "Expected: ", expected_X_cols, " columns (", Ldim, " TD covariates × ", T, " time points + ", n_baseline, " baseline covariates)\n",
         "Actual: ", ncol(X), " columns\n",
         "X should contain: [TD_Cov1_T1, TD_Cov2_T1, ..., TD_CovLdim_T1, TD_Cov1_T2, ..., TD_CovLdim_T(T-1), Baseline1, Baseline2, ...]")
  }

  # Check for missing values
  if (any(is.na(Z))) warning("Z contains missing values - this may cause issues in g-BKMR analysis")
  if (any(is.na(X))) warning("X contains missing values - this may cause issues in g-BKMR analysis")

  cat("✓ Input validation passed\n")
}

# Example usage and demonstration
example_matrix_conversion <- function() {
  cat("=== Example: Converting User Matrices to g-BKMR Format ===\n\n")

  # Simulate realistic user data
  set.seed(42)
  n <- 200          # 200 subjects
  T <- 4            # 4 time points
  Adim <- 3         # 3 metals (Arsenic, Manganese, Lead)
  Ldim <- 2         # 2 time-dependent covariates (BMI, Blood Pressure)
  n_baseline <- 2   # 2 baseline covariates (Sex, Age)

  cat("User scenario:\n")
  cat("- Analyzing the effect of metal mixtures on cardiovascular outcomes\n")
  cat("- 3 metals measured at 4 time points each\n")
  cat("- 2 time-dependent confounders (BMI, Blood Pressure) measured at all time points\n")
  cat("- 2 baseline covariates (Sex, Age)\n\n")

  # Generate user matrices

  # Y: Cardiovascular risk score (continuous outcome)
  Y <- rnorm(n, mean = 50, sd = 15)

  # Z: Metal concentrations (positive values, log-normal distributed)
  # Columns: [As_T0, Mn_T0, Pb_T0, As_T1, Mn_T1, Pb_T1, As_T2, Mn_T2, Pb_T2, As_T3, Mn_T3, Pb_T3]
  Z <- matrix(rlnorm(n * Adim * T, meanlog = 1, sdlog = 0.5), nrow = n, ncol = Adim * T)
  colnames(Z) <- rep(c("Arsenic", "Manganese", "Lead"), T)

  # X: Time-dependent covariates + baseline covariates
  # Structure: [BMI_T1, BP_T1, BMI_T2, BP_T2, BMI_T3, BP_T3, BMI_T4, BP_T4, Sex, Age]
  # Note: No TD covariates at T0 (that becomes baseline)
  td_data <- matrix(rnorm(n * Ldim * T, mean = c(25, 120), sd = c(5, 15)), nrow = n, ncol = Ldim * T)
  baseline_data <- cbind(
    sex = rbinom(n, 1, 0.5),           # Sex (0/1)
    age = rnorm(n, mean = 45, sd = 12) # Age
  )
  X <- cbind(td_data, baseline_data)

  cat("User input matrices:\n")
  cat("Y: ", class(Y), " vector of length ", length(Y), "\n")
  cat("Z: ", nrow(Z), "×", ncol(Z), " matrix (metal concentrations)\n")
  cat("X: ", nrow(X), "×", ncol(X), " matrix (time-dependent + baseline covariates)\n\n")

  cat("Expected matrix organization:\n")
  cat("Z columns: As_T0, Mn_T0, Pb_T0, As_T1, Mn_T1, Pb_T1, As_T2, Mn_T2, Pb_T2, As_T3, Mn_T3, Pb_T3\n")
  cat("X columns: BMI_T1, BP_T1, BMI_T2, BP_T2, BMI_T3, BP_T3, BMI_T4, BP_T4, Sex, Age\n\n")

  # Convert to g-BKMR format
  gbkmr_data <- prepare_gbkmr_data(
    Y = Y,
    Z = Z,
    X = X,
    time_points = T,
    mixture_components = Adim,
    td_covariates = Ldim,
    baseline_covariates = n_baseline,
    log_transform_mixtures = TRUE
  )

  cat("Resulting g-BKMR data structure:\n")
  cat("Columns (", ncol(gbkmr_data), " total): ", paste(colnames(gbkmr_data), collapse = ", "), "\n\n")

  cat("First 3 rows of converted data:\n")
  print(head(gbkmr_data, 3))

  cat("\n=== Ready for g-BKMR analysis! ===\n")
  cat("Next step: Use gbkmr_run() with this prepared data\n\n")

  return(gbkmr_data)
}

# Quick data validation function
check_gbkmr_data <- function(data) {
  info <- attr(data, "data_info")
  if (is.null(info)) {
    cat("⚠ This data was not prepared with prepare_gbkmr_data()\n")
    return(FALSE)
  }

  cat("✓ Data validation summary:\n")
  cat("  Sample size:", info$n, "\n")
  cat("  Time points:", info$time_points, "\n")
  cat("  Mixture components:", info$mixture_components, "\n")
  cat("  Time-dependent covariates:", info$td_covariates, "\n")
  cat("  Baseline covariates:", info$baseline_covariates, "\n")
  cat("  Log-transformed:", info$log_transformed, "\n")
  cat("  Data source:", info$source, "\n")

  # Check for required columns
  required_cols <- c("sex", "Y", "id")
  missing_req <- setdiff(required_cols, colnames(data))
  if (length(missing_req) > 0) {
    cat("⚠ Missing required columns:", paste(missing_req, collapse = ", "), "\n")
    return(FALSE)
  }

  cat("✓ Data is ready for g-BKMR analysis\n")
  return(TRUE)
}
