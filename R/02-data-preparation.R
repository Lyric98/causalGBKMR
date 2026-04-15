#' @file 02-data-preparation.R
#' @title Data preparation and formatting for g-BKMR
#' @importFrom stats rbinom rnorm sd
#' @description Functions to convert user data into the wide-format data structure
#' required for g-BKMR analysis. Handles variable naming, transformations, and metadata.

#' Prepare user matrices for g-BKMR analysis
#'
#' @description Converts user-provided matrices (Y, Z, X) into the wide-format
#' data structure required for g-BKMR analysis. Supports both continuous and
#' binary time-dependent covariates with enhanced input validation.
#'
#' @param Y Numeric vector. Outcome variable (length n).
#' @param Z Numeric matrix. Mixture exposure matrix (n x (Adim x T)).
#' @param X Numeric matrix. Covariate matrix (n x (Ldim x T + baseline_covs)).
#' @param time_points Integer. Number of time points (T).
#' @param mixture_components Integer. Number of mixture components per time point (Adim).
#' @param td_covariates Integer. Number of time-dependent covariates per time point (Ldim).
#' @param baseline_covariates Integer. Number of baseline covariates (default: 1).
#' @param td_covariate_names Character vector. Names for time-dependent covariates (optional).
#' @param log_transform_mixtures Logical. Whether to log-transform mixture exposures (default: TRUE).
#' @param validate_input Logical. Whether to validate input dimensions (default: TRUE).
#'
#' @return A data frame in g-BKMR format with proper variable naming and metadata.
#'
#' @details
#' The function expects matrices organized as follows:
#' \itemize{
#'   \item Z matrix: Mixtures in chronological order (Mix1_T0, Mix2_T0, ..., MixAdim_T0, Mix1_T1, ...)
#'   \item X matrix: (TD_Cov1_T1, TD_Cov2_T1, ..., TD_CovLdim_T1, ..., TD_CovLdim_T(T-1), Baseline1, Baseline2, ...)
#' }
#'
#' The output data frame has the following structure:
#' \itemize{
#'   \item sex: First baseline covariate (required by g-BKMR format)
#'   \item baseline_2, baseline_3, ...: Additional baseline covariates
#'   \item td_covariate1_0, td_covariate2_0, ...: Baseline time-dependent covariates
#'   \item logM1_0, logM2_0, ...: Mixture exposures at time 0
#'   \item logM1_1, logM2_1, ...: Mixture exposures at time 1
#'   \item td_covariate1_1, td_covariate2_1, ...: Time-dependent covariates at time 1
#'   \item Y: Outcome variable
#'   \item id: Subject identifier
#' }
#'
#' @examples
#' \dontrun{
#' # Generate test data
#' n <- 200
#' Y <- rnorm(n)
#' Z <- matrix(rlnorm(n * 6), nrow = n, ncol = 6)  # 2 metals x 3 time points
#' X <- matrix(rnorm(n * 8), nrow = n, ncol = 8)   # 2 TD covs x 3 time points + 2 baseline
#'
#' # Prepare data
#' prepared_data <- prepare_gbkmr_data(
#'   Y = Y, Z = Z, X = X,
#'   time_points = 3,
#'   mixture_components = 2,
#'   td_covariates = 2,
#'   baseline_covariates = 2,
#'   td_covariate_names = c("bmi", "bp")
#' )
#'
#' # Check the structure
#' str(prepared_data)
#' head(prepared_data)
#' }
#'
#' @export
prepare_gbkmr_data <- function(
    Y,                     # Outcome vector (length n)
    Z,                     # Mixture exposure matrix (n x (Adim * T))
    X,                     # Covariate matrix (n x (Ldim * T + baseline_covs))
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

  # ========== MODIFICATION 1: Handle non-numeric matrices ==========
  # Convert X matrix to numeric if it contains non-numeric columns
  if (!is.numeric(X)) {
    # X is a character matrix (common when as.matrix() is applied to mixed types)
    # Convert each column appropriately
    X_numeric <- matrix(nrow = nrow(X), ncol = ncol(X))
    for (j in 1:ncol(X)) {
      col_data <- X[, j]
      # Try to convert to numeric
      col_numeric <- suppressWarnings(as.numeric(col_data))
      # If conversion fails (produces NAs where there weren't any), treat as factor
      if (sum(is.na(col_numeric)) > sum(is.na(col_data))) {
        # This is a categorical variable
        col_numeric <- as.numeric(factor(col_data))
      }
      X_numeric[, j] <- col_numeric
    }
    X <- X_numeric
  }

  # Ensure Z is numeric
  if (!is.numeric(Z)) {
    Z <- apply(Z, 2, as.numeric)
  }
  # ================================================================

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
      base_values <- as.numeric(X[, first_timepoint_col])
      baseline_val <- base_values + rnorm(n, 0, sd = 0.1 * sd(base_values, na.rm = TRUE))
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
  cat("[OK] Data conversion successful!\n")
  cat("[OK] Generated", nrow(df), "x", ncol(df), "data frame\n")
  if (use_user_names) {
    cat("[OK] Using user covariate names:", paste(td_covariate_names, collapse = ", "), "\n")
  } else {
    cat("[OK] Using generic covariate names:", paste(td_covariate_names, collapse = ", "), "\n")
  }
  cat("[OK] Ready for g-BKMR analysis\n\n")

  return(df)
}
