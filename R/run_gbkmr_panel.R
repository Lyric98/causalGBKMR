run_gbkmr_panel <- function(
    sim_popn,
    T = 5,
    currind = 1,
    sel = seq(22000, 24000, by = 25),
    n = 500,
    K = 1000,
    iter = 24000,
    parallel = TRUE,
    save_exposure_preds = TRUE,
    return_ci = TRUE,
    make_plots = TRUE,
    use_knots = TRUE,
    n_knots = 50
) {
  library(bkmr)
  library(fields)
  library(dplyr)
  library(ggplot2)
  library(parallel)

  if (max(sel) > iter) stop("sel contains indices beyond total MCMC iterations!")

  message("Subsampling population...")
  set.seed(currind)
  dat_sim <- sim_popn[sample(sim_popn$id, n, replace = FALSE), ]

  # Auto-detect dimensions and variable patterns
  detection_results <- detect_variable_patterns(dat_sim, T)
  p <- detection_results$p
  Ldim <- detection_results$Ldim
  td_covariate_names <- detection_results$td_covariate_names
  baseline_td_vars <- detection_results$baseline_td_vars

  message("Detected data structure:")
  message("  - Number of exposures per time point (p): ", p)
  message("  - Number of time-dependent covariates per time point (Ldim): ", Ldim)
  message("  - Time-dependent covariate names: ", paste(td_covariate_names, collapse = ", "))
  message("  - Baseline TD variables: ", paste(baseline_td_vars, collapse = ", "))

  # Common covariates - use detected baseline TD variables
  cov_names_common <- c("sex", baseline_td_vars)
  X_common <- dat_sim[, cov_names_common]

  fitkm_list <- list()
  L_values_a <- list()
  L_values_astar <- list()

  for (t in 1:(T - 1)) {
    message(paste("Fitting mediator model L", t))

    # Get time-dependent covariate names for this time point
    td_vars_t <- detection_results$td_vars_by_time[[t]]

    if (Ldim == 1) {
      y_L <- dat_sim[, td_vars_t]
    } else {
      y_L <- as.matrix(dat_sim[, td_vars_t])
    }

    # Build exposure variable names for times 0 to t-1
    exp_names <- unlist(lapply(0:(t - 1), function(s) paste0("logM", 1:p, "_", s)))

    # Add previous time-dependent covariates
    if (t > 1) {
      for (j in 1:(t - 1)) {
        prev_td_vars <- detection_results$td_vars_by_time[[j]]
        exp_names <- c(exp_names, prev_td_vars)
      }
    }

    Z <- dat_sim[, exp_names]
    X <- X_common

    scale_Z <- scale(Z)
    attr_list <- list(center = attr(scale_Z, "scaled:center"),
                      scale = attr(scale_Z, "scaled:scale"))
    knots <- if (use_knots) {
      set.seed(1000)
      fields::cover.design(scale_Z, nd = n_knots)$design
    } else NULL

    fit <- kmbayes(
      y = if (is.matrix(y_L)) y_L[, 1] else y_L,
      Z = scale_Z,
      X = X,
      iter = iter,
      varsel = TRUE,
      verbose = FALSE,
      knots = knots
    )
    fitkm_list[[paste0("L", t)]] <- list(fit = fit, scale_info = attr_list)
  }

  message("Fitting outcome model Y")
  Y <- dat_sim$Y

  # Auto-detect all exposure and time-dependent covariate variables
  all_exp_vars <- grep("^logM\\d+_\\d+$", names(dat_sim), value = TRUE)
  all_td_vars <- unlist(detection_results$td_vars_by_time)
  exp_names_y <- c(all_exp_vars, all_td_vars)

  Z_y <- dat_sim[, exp_names_y]
  scale_Zy <- scale(Z_y)
  scale_info_y <- list(center = attr(scale_Zy, "scaled:center"),
                       scale = attr(scale_Zy, "scaled:scale"))
  knots_y <- if (use_knots) {
    set.seed(1000)
    fields::cover.design(scale_Zy, nd = n_knots * 2)$design
  } else NULL

  fit_y <- kmbayes(y = Y, Z = scale_Zy, X = X_common, iter = iter, varsel = TRUE,
                   verbose = FALSE, knots = knots_y)
  fitkm_y <- list(fit = fit_y, scale_info = scale_info_y)

  message("Computing counterfactual exposure vectors...")
  X.predict <- matrix(colMeans(X_common), nrow = 1)

  for (t in 1:(T - 1)) {
    message(paste("Predicting mediator L", t))
    scale_info <- fitkm_list[[paste0("L", t)]]$scale_info
    fit <- fitkm_list[[paste0("L", t)]]$fit

    cols_exp <- unlist(lapply(0:(t - 1), function(s) paste0("logM", 1:p, "_", s)))
    cols_med <- if (t > 1) {
      unlist(detection_results$td_vars_by_time[1:(t-1)])
    } else character(0)
    cols_all <- c(cols_exp, cols_med)

    A_full <- dat_sim[, cols_all, drop = FALSE]
    a_vec <- apply(A_full, 2, quantile, probs = 0.25)
    astar_vec <- apply(A_full, 2, quantile, probs = 0.75)
    newz <- rbind(a_vec[cols_all], astar_vec[cols_all])

    Znew_scaled <- scale(newz, center = scale_info$center, scale = scale_info$scale)
    L_pred <- SamplePred(fit, Znew = Znew_scaled, Xnew = X.predict, sel = sel)

    L_values_a[[t]] <- as.vector(L_pred[, "znew1"])
    L_values_astar[[t]] <- as.vector(L_pred[, "znew2"])
  }

  message("Predicting outcome Y")
  A_full_y <- dat_sim[, exp_names_y, drop = FALSE]
  a_vec <- apply(A_full_y, 2, quantile, probs = 0.25)
  astar_vec <- apply(A_full_y, 2, quantile, probs = 0.75)
  Z_final <- rbind(a_vec[exp_names_y], astar_vec[exp_names_y])
  Z_final_scaled <- scale(Z_final, center = fitkm_y$scale_info$center,
                          scale = fitkm_y$scale_info$scale)
  Y_pred <- SamplePred(fitkm_y$fit, Znew = Z_final_scaled, Xnew = X.predict, sel = sel)

  Ya <- Y_pred[, "znew1"]
  Yastar <- Y_pred[, "znew2"]
  diff_gBKMR <- mean(Yastar) - mean(Ya)

  message(sprintf("g-BKMR effect estimate: %.4f", diff_gBKMR))

  if (make_plots) {
    df_plot <- data.frame(Scenario = c("a", "astar"), Mean = c(mean(Ya), mean(Yastar)))
    print(ggplot(df_plot, aes(x = Scenario, y = Mean)) +
            geom_col(fill = "skyblue") +
            theme_minimal() +
            ggtitle("Counterfactual Means"))
  }

  # Add detection results to output for transparency
  results <- list(
    diff_gBKMR = diff_gBKMR,
    Ya = Ya,
    Yastar = Yastar,
    beta_all = c(
      unlist(lapply(fitkm_list, function(l) colMeans(l$fit$beta[sel, , drop = FALSE]))),
      colMeans(fitkm_y$fit$beta[sel, , drop = FALSE])
    ),
    L_values_a = L_values_a,
    L_values_astar = L_values_astar,
    detection_info = detection_results,  # Add detection info for user reference
    fitted_models = if (save_exposure_preds) list(
      confounders = fitkm_list,
      outcome = fitkm_y
    ) else NULL
  )

  return(results)
}

# Updated helper function to detect variable patterns with better naming
detect_variable_patterns <- function(data, T) {

  # Detect number of exposures per time point (p)
  p <- length(grep("^logM\\d+_0$", names(data)))

  # Try to detect actual time-dependent covariate names from data
  # Look for baseline TD covariates with various patterns
  baseline_td_patterns <- list(
    # Pattern 1: Known covariate names with 0 (e.g., bmi_0, bp_0)
    "known_with_underscore" = "^[a-zA-Z_]+_0$",
    # Pattern 2: Known covariate names ending with 0 (e.g., bmi0, bp0)
    "known_ending_zero" = "^[a-zA-Z_]+0$",
    # Pattern 3: Generated format (waist0_1, waist0_2, etc.)
    "generated_format" = "^[a-zA-Z_]+0_\\d+$",
    # Pattern 4: Generic td_covariate format
    "generic_format" = "^td_covariate\\d+_0$"
  )

  baseline_td_vars <- character(0)
  detected_pattern <- ""
  td_covariate_names <- character(0)

  # Try each pattern to find baseline TD covariates
  for (pattern_name in names(baseline_td_patterns)) {
    pattern <- baseline_td_patterns[[pattern_name]]
    matches <- grep(pattern, names(data), value = TRUE)
    if (length(matches) > 0) {
      baseline_td_vars <- matches
      detected_pattern <- pattern_name

      # Extract base names from the matches
      if (pattern_name == "known_with_underscore") {
        # Extract names like "bmi", "bp" from "bmi_0", "bp_0"
        td_covariate_names <- gsub("_0$", "", matches)
      } else if (pattern_name == "known_ending_zero") {
        # Extract names like "bmi", "bp" from "bmi0", "bp0"
        td_covariate_names <- gsub("0$", "", matches)
      } else if (pattern_name == "generated_format") {
        # For waist0_1, waist0_2, extract "waist" and number the covariates
        base_name <- gsub("0_\\d+$", "", matches[1])
        n_covariates <- length(matches)
        td_covariate_names <- paste0(base_name, 1:n_covariates)
      } else if (pattern_name == "generic_format") {
        # For td_covariate1_0, td_covariate2_0, extract the names
        td_covariate_names <- gsub("_0$", "", matches)
      }
      break
    }
  }

  # If no baseline TD covariates found, set defaults
  if (length(baseline_td_vars) == 0) {
    baseline_td_vars <- character(0)
    Ldim <- 0
    td_covariate_names <- character(0)
    warning("No baseline time-dependent covariates detected. Setting Ldim = 0.")
  } else {
    Ldim <- length(baseline_td_vars)
  }

  # Detect time-dependent covariates for each time point t = 1, 2, ..., T-1
  td_vars_by_time <- list()

  if (Ldim > 0 && T > 1) {
    for (t in 1:(T-1)) {
      td_vars_t <- character(0)

      # Generate possible variable names based on detected pattern
      if (detected_pattern == "known_with_underscore") {
        # Look for bmi_1, bp_1, bmi_2, bp_2, etc.
        potential_vars <- paste0(td_covariate_names, "_", t)
      } else if (detected_pattern == "known_ending_zero") {
        # Look for bmi1, bp1, bmi2, bp2, etc.
        potential_vars <- paste0(td_covariate_names, t)
      } else if (detected_pattern == "generated_format") {
        # Look for waist1_1, waist1_2, waist2_1, waist2_2, etc.
        base_name <- gsub("0_\\d+$", "", baseline_td_vars[1])
        potential_vars <- paste0(base_name, t, "_", 1:Ldim)
      } else if (detected_pattern == "generic_format") {
        # Look for td_covariate1_1, td_covariate2_1, etc.
        potential_vars <- paste0(td_covariate_names, "_", t)
      }

      # Check which of these potential variables exist in the data
      existing_vars <- intersect(potential_vars, names(data))

      if (length(existing_vars) == Ldim) {
        td_vars_t <- existing_vars
      } else {
        # Fallback: try pattern matching
        time_patterns <- c(
          paste0("^.*", t, "_\\d+$"),    # any_name1_1, any_name1_2, etc.
          paste0("^.*_", t, "_\\d+$"),   # any_name_1_1, any_name_1_2, etc.
          paste0("^.*", t, "$")          # any_name1 (single TD covariate)
        )

        for (pattern in time_patterns) {
          matches <- grep(pattern, names(data), value = TRUE)
          # Exclude exposure variables and outcome/id
          matches <- matches[!grepl("^logM\\d+_", matches)]
          matches <- matches[!matches %in% c("Y", "id")]

          if (length(matches) >= Ldim) {
            td_vars_t <- matches[1:Ldim]
            break
          }
        }
      }

      if (length(td_vars_t) != Ldim) {
        stop("Cannot detect time-dependent covariates for time ", t,
             ". Expected ", Ldim, " variables but found ", length(td_vars_t),
             ". Available variables: ", paste(names(data), collapse = ", "))
      }

      td_vars_by_time[[t]] <- td_vars_t
    }
  }

  return(list(
    p = p,
    Ldim = Ldim,
    td_covariate_names = td_covariate_names,
    detected_pattern = detected_pattern,
    baseline_td_vars = baseline_td_vars,
    td_vars_by_time = td_vars_by_time
  ))
}

# Updated prepare_gbkmr_data function with proper naming
prepare_gbkmr_data <- function(
    Y,                     # Outcome vector (length n)
    Z,                     # Mixture exposure matrix (n × (Adim * T))
    X,                     # Covariate matrix (n × (Ldim * T + baseline_covs))
    time_points,           # Number of time points (T)
    mixture_components,    # Number of mixture components per time point (Adim)
    td_covariates,         # Number of time-dependent covariates per time point (Ldim)
    baseline_covariates = 1,  # Number of baseline covariates (default: sex only)
    td_covariate_names = NULL,  # User-provided names for time-dependent covariates
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

  # Determine naming strategy for time-dependent covariates
  use_user_names <- !is.null(td_covariate_names) && length(td_covariate_names) == Ldim

  if (use_user_names) {
    cat("Using user-provided time-dependent covariate names: ", paste(td_covariate_names, collapse = ", "), "\n")
  } else {
    # Use generic names when user doesn't provide them
    td_covariate_names <- paste0("td_covariate", 1:Ldim)
    cat("Using generic time-dependent covariate names: ", paste(td_covariate_names, collapse = ", "), "\n")
  }

  cat("Converting user matrices to g-BKMR format...\n")
  cat("Data dimensions: n =", n, ", T =", T, ", Adim =", Adim, ", Ldim =", Ldim, "\n")

  # Initialize the data frame
  df <- data.frame(id = 1:n)

  # Add baseline covariates
  if (n_baseline > 0) {
    baseline_start_col <- ncol(X) - n_baseline + 1
    baseline_data <- X[, baseline_start_col:ncol(X), drop = FALSE]

    # First baseline covariate is always 'sex' (required by format)
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

  # Add baseline time-dependent covariates with proper naming
  for (l in 1:Ldim) {
    if (use_user_names) {
      # Use user's actual covariate name with _0 suffix
      var_name <- paste0(td_covariate_names[l], "_0")
    } else {
      # Use generic name
      var_name <- paste0("td_covariate", l, "_0")
    }

    # Extract first time point of this time-dependent covariate from X
    first_timepoint_col <- l
    if (first_timepoint_col <= (ncol(X) - n_baseline)) {
      baseline_val <- X[, first_timepoint_col] + rnorm(n, 0, sd = 0.1 * sd(X[, first_timepoint_col]))
    } else {
      baseline_val <- rnorm(n, mean = 0, sd = 1)
    }

    df[[var_name]] <- baseline_val
  }

  # Add mixture exposures in chronological order (logM1_0, logM2_0, ..., logMAdim_T-1)
  for (t in 0:(T-1)) {
    for (a in 1:Adim) {
      col_idx <- t * Adim + a
      var_name <- paste0("logM", a, "_", t)

      mixture_data <- Z[, col_idx]

      # Log-transform if requested
      if (log_transform_mixtures) {
        if (any(mixture_data <= 0, na.rm = TRUE)) {
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

  # Add time-dependent covariates with proper naming
  for (t in 1:(T-1)) {
    for (l in 1:Ldim) {
      if (use_user_names) {
        # Use user's actual covariate name with _t suffix
        var_name <- paste0(td_covariate_names[l], "_", t)
      } else {
        # Use generic name
        var_name <- paste0("td_covariate", l, "_", t)
      }

      # Calculate column index in X matrix
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

  # Build the correct column order
  ordered_cols <- "sex"

  # Add additional baseline covariates
  if (n_baseline > 1) {
    ordered_cols <- c(ordered_cols, paste0("baseline_", 2:n_baseline))
  }

  # Add baseline time-dependent covariates
  if (Ldim > 0) {
    if (use_user_names) {
      ordered_cols <- c(ordered_cols, paste0(td_covariate_names, "_0"))
    } else {
      ordered_cols <- c(ordered_cols, paste0("td_covariate", 1:Ldim, "_0"))
    }
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
      if (use_user_names) {
        ordered_cols <- c(ordered_cols, paste0(td_covariate_names, "_", t))
      } else {
        ordered_cols <- c(ordered_cols, paste0("td_covariate", 1:Ldim, "_", t))
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
    td_covariate_names = td_covariate_names,
    use_user_names = use_user_names,
    log_transformed = log_transform_mixtures,
    source = "user_matrices"
  )

  # Print summary
  cat("✓ Data conversion successful!\n")
  cat("✓ Generated", nrow(df), "×", ncol(df), "data frame\n")
  if (use_user_names) {
    cat("✓ Using user covariate names:", paste(td_covariate_names, collapse = ", "), "\n")
  } else {
    cat("✓ Using generic covariate names:", paste(td_covariate_names, collapse = ", "), "\n")
  }
  cat("✓ Ready for g-BKMR analysis\n\n")

  return(df)
}

# Validation function (unchanged)
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
  expected_X_cols <- Ldim * T + n_baseline

  if (ncol(Z) != expected_Z_cols) {
    stop("Z matrix dimension error!\n",
         "Expected: ", expected_Z_cols, " columns (", Adim, " mixtures × ", T, " time points)\n",
         "Actual: ", ncol(Z), " columns")
  }

  if (ncol(X) != expected_X_cols) {
    stop("X matrix dimension error!\n",
         "Expected: ", expected_X_cols, " columns (", Ldim, " TD covariates × ", T, " time points + ", n_baseline, " baseline covariates)\n",
         "Actual: ", ncol(X), " columns")
  }

  if (any(is.na(Z))) warning("Z contains missing values")
  if (any(is.na(X))) warning("X contains missing values")

  cat("✓ Input validation passed\n")
}

# Updated gbkmr_run function (unchanged core logic)
gbkmr_run <- function(
    data,
    outcome = "Y",
    outcome_type = c("continuous", "binary"),
    time_points,

    # Optional advanced parameters
    currind = 1,
    sel = NULL,
    n = NULL,
    K = 1000,
    iter = 15000,
    parallel = TRUE,
    use_knots = TRUE,
    n_knots = 50,
    make_plots = FALSE,
    verbose = TRUE
) {

  # Validate inputs
  outcome_type <- match.arg(outcome_type)

  if (missing(data)) stop("Data must be provided")
  if (missing(time_points)) stop("Number of time points must be provided")

  # Set defaults
  if (is.null(sel)) sel <- seq(iter * 0.6, iter, by = 25)
  if (is.null(n)) n <- min(500, nrow(data))

  if (verbose) {
    cat("Starting g-BKMR analysis...\n")
    cat("Data dimensions:", nrow(data), "subjects,", ncol(data), "variables\n")
    cat("Time points:", time_points, "\n")
    cat("Sample size for analysis:", n, "\n")
  }

  # Run the core g-BKMR function with auto-detection
  results <- run_gbkmr_panel(
    sim_popn = data,
    T = time_points,
    currind = currind,
    sel = sel,
    n = n,
    K = K,
    iter = iter,
    parallel = parallel,
    save_exposure_preds = TRUE,
    return_ci = TRUE,
    make_plots = make_plots,
    use_knots = use_knots,
    n_knots = n_knots
  )

  # Format results for user-friendly output
  causal_effect <- list(
    estimate = results$diff_gBKMR,
    lower = quantile(results$Yastar - results$Ya, 0.025, na.rm = TRUE),
    upper = quantile(results$Yastar - results$Ya, 0.975, na.rm = TRUE)
  )

  counterfactual_means <- list(
    low = mean(results$Ya, na.rm = TRUE),
    high = mean(results$Yastar, na.rm = TRUE)
  )

  formatted_results <- list(
    causal_effect = causal_effect,
    counterfactual_means = counterfactual_means,
    variable_importance = results$beta_all,
    detection_info = results$detection_info,
    raw_results = results,
    call_info = list(
      outcome = outcome,
      outcome_type = outcome_type,
      time_points = time_points,
      sample_size = n,
      mcmc_iterations = iter
    )
  )

  class(formatted_results) <- "gbkmr_results"

  if (verbose) {
    cat("Analysis complete!\n")
    cat("Detected", results$detection_info$p, "exposures per time point\n")
    cat("Detected", results$detection_info$Ldim, "time-dependent covariates per time point\n")
    if (length(results$detection_info$td_covariate_names) > 0) {
      cat("TD covariate names:", paste(results$detection_info$td_covariate_names, collapse = ", "), "\n")
    }
    cat("Causal effect estimate:", round(causal_effect$estimate, 4), "\n")
    cat("95% CI: (", round(causal_effect$lower, 4), ",", round(causal_effect$upper, 4), ")\n")
  }

  return(formatted_results)
}
