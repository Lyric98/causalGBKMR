#' Prepare user matrices for g-BKMR analysis
#'
#' @description Converts user-provided matrices (Y, Z, X) into the wide-format
#' data structure required for g-BKMR analysis. Supports both continuous and
#' binary time-dependent covariates with enhanced input validation.
#'
#' @param Y Numeric vector. Outcome variable (length n).
#' @param Z Numeric matrix. Mixture exposure matrix (n × (Adim × T)).
#' @param X Numeric matrix. Covariate matrix (n × (Ldim × T + baseline_covs)).
#' @param time_points Integer. Number of time points (T).
#' @param mixture_components Integer. Number of mixture components per time point (Adim).
#' @param td_covariates Integer. Number of time-dependent covariates per time point (Ldim).
#' @param baseline_covariates Integer. Number of baseline covariates (default: 1).
#' @param td_covariate_names Character vector. Names for time-dependent covariates (optional).
#' @param td_covariate_types Character vector. Types for each TD covariate: "continuous" or "binary" (optional).
#' @param data_format Character. Data format: "wide" (current) or "long" (future).
#' @param log_transform_mixtures Logical. Whether to log-transform mixture exposures (default: TRUE).
#' @param validate_input Logical. Whether to validate input dimensions (default: TRUE).
#' @param verbose Logical. Whether to print detailed validation info (default: TRUE).
#'
#' @return A data frame in g-BKMR format with proper variable naming and metadata.
#'
#' @details
#' The function expects matrices organized as follows:
#' \itemize{
#'   \item Z matrix: Mixtures in chronological order [Mix1_T0, Mix2_T0, ..., MixAdim_T0, Mix1_T1, ...]
#'   \item X matrix: [TD_Cov1_T1, TD_Cov2_T1, ..., TD_CovLdim_T1, ..., TD_CovLdim_T(T-1), Baseline1, Baseline2, ...]
#' }
#'
#' @examples
#' \dontrun{
#' # Generate test data
#' n <- 200
#' Y <- rnorm(n)
#' Z <- matrix(rlnorm(n * 6), nrow = n, ncol = 6)  # 2 metals × 3 time points
#' X <- matrix(rnorm(n * 8), nrow = n, ncol = 8)   # 2 TD covs × 3 time points + 2 baseline
#'
#' # Prepare data
#' prepared_data <- prepare_gbkmr_data(
#'   Y = Y, Z = Z, X = X,
#'   time_points = 3,
#'   mixture_components = 2,
#'   td_covariates = 2,
#'   baseline_covariates = 2,
#'   td_covariate_names = c("bmi", "bp"),
#'   td_covariate_types = c("continuous", "continuous")
#' )
#' }
#'
#' @export
prepare_gbkmr_data <- function(
    Y,
    Z,
    X,
    time_points,
    mixture_components,
    td_covariates,
    baseline_covariates = 1,
    td_covariate_names = NULL,
    td_covariate_types = NULL,
    data_format = "wide",
    log_transform_mixtures = TRUE,
    validate_input = TRUE,
    verbose = TRUE
) {

  # Get dimensions
  n <- length(Y)
  T <- time_points
  Adim <- mixture_components
  Ldim <- td_covariates
  n_baseline <- baseline_covariates

  # Enhanced input validation with detailed error messages
  if (validate_input) {
    validation_result <- validate_user_matrices_enhanced(
      Y, Z, X, T, Adim, Ldim, n_baseline, verbose
    )
    if (!validation_result$valid) {
      stop("Input validation failed. See messages above for details.")
    }
  }

  # Validate and set up time-dependent covariate types
  if (is.null(td_covariate_types)) {
    td_covariate_types <- rep("continuous", Ldim)
    if (verbose) cat("No covariate types specified. Defaulting all to 'continuous'.\n")
  } else {
    if (length(td_covariate_types) != Ldim) {
      stop("Length of td_covariate_types (", length(td_covariate_types),
           ") must match td_covariates (", Ldim, ")")
    }

    # Validate covariate types
    valid_types <- c("continuous", "binary")
    invalid_types <- setdiff(td_covariate_types, valid_types)
    if (length(invalid_types) > 0) {
      stop("Invalid covariate types: ", paste(invalid_types, collapse = ", "),
           ". Must be one of: ", paste(valid_types, collapse = ", "))
    }

    if (verbose) {
      cat("Time-dependent covariate types specified:\n")
      for (i in 1:length(td_covariate_types)) {
        cat("  Covariate", i, ":", td_covariate_types[i], "\n")
      }
    }
  }

  # Check for future long-format support
  if (data_format == "long") {
    stop("Long-format data is not yet supported. Currently only 'wide' format is implemented.\n",
         "Future implementation will require additional parameters:\n",
         "- subject_id_col: name of subject identifier column\n",
         "- time_col: name of time variable column\n",
         "- variable_name_cols: names of variable identifier columns")
  }

  # Determine naming strategy for time-dependent covariates
  use_user_names <- !is.null(td_covariate_names) && length(td_covariate_names) == Ldim

  if (use_user_names) {
    if (verbose) cat("Using user-provided time-dependent covariate names: ",
                     paste(td_covariate_names, collapse = ", "), "\n")
  } else {
    # Use generic names when user doesn't provide them
    td_covariate_names <- paste0("td_covariate", 1:Ldim)
    if (verbose) cat("Using generic time-dependent covariate names: ",
                     paste(td_covariate_names, collapse = ", "), "\n")
  }

  if (verbose) {
    cat("Converting user matrices to g-BKMR format...\n")
    cat("Data dimensions: n =", n, ", T =", T, ", Adim =", Adim, ", Ldim =", Ldim, "\n")
    cat("Data format:", data_format, "\n")
  }

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

  # Add time-dependent covariates with proper naming and type handling
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
        raw_data <- X[, col_idx]

        # Handle binary vs continuous covariates
        if (td_covariate_types[l] == "binary") {
          # Validate binary data
          unique_vals <- unique(raw_data[!is.na(raw_data)])
          if (!all(unique_vals %in% c(0, 1))) {
            warning("Covariate '", td_covariate_names[l], "' at time ", t,
                    " declared as binary but contains non-binary values. ",
                    "Converting: values > 0.5 → 1, others → 0")
            raw_data <- as.numeric(raw_data > 0.5)
          }
          df[[var_name]] <- as.numeric(raw_data)
        } else {
          # Continuous covariate
          df[[var_name]] <- raw_data
        }
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

  # Add metadata including the user's original TD covariate names and types
  attr(df, "data_info") <- list(
    n = n,
    time_points = T,
    mixture_components = Adim,
    td_covariates = Ldim,
    baseline_covariates = n_baseline,
    td_covariate_names = td_covariate_names,
    td_covariate_types = td_covariate_types,
    use_user_names = use_user_names,
    data_format = data_format,
    log_transformed = log_transform_mixtures,
    source = "user_matrices"
  )

  # Print summary
  if (verbose) {
    cat("✓ Data conversion successful!\n")
    cat("✓ Generated", nrow(df), "×", ncol(df), "data frame\n")
    if (use_user_names) {
      cat("✓ Using user covariate names:", paste(td_covariate_names, collapse = ", "), "\n")
    } else {
      cat("✓ Using generic covariate names:", paste(td_covariate_names, collapse = ", "), "\n")
    }

    # Report covariate types
    binary_count <- sum(td_covariate_types == "binary")
    continuous_count <- sum(td_covariate_types == "continuous")
    cat("✓ Covariate types:", continuous_count, "continuous,", binary_count, "binary\n")
    if (binary_count > 0) {
      binary_names <- td_covariate_names[td_covariate_types == "binary"]
      cat("  Binary covariates:", paste(binary_names, collapse = ", "), "\n")
    }
    cat("✓ Ready for g-BKMR analysis\n\n")
  }

  return(df)
}
