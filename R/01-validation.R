#' @file 01-validation.R
#' @title Input validation utilities for g-BKMR
#' @description Contains functions to validate user inputs before g-BKMR analysis.
#' These functions ensure data quality and parameter consistency.

#' Validate user-provided matrices for g-BKMR analysis
#'
#' @description Performs comprehensive validation of input matrices Y, Z, X
#' to ensure they meet the requirements for g-BKMR analysis.
#'
#' @param Y Numeric vector. Outcome variable (length n).
#' @param Z Numeric matrix. Mixture exposure matrix (n x (Adim x T)).
#' @param X Numeric matrix. Covariate matrix (n x (Ldim x T + baseline_covs)).
#' @param T Integer. Number of time points.
#' @param Adim Integer. Number of mixture components per time point.
#' @param Ldim Integer. Number of time-dependent covariates per time point.
#' @param n_baseline Integer. Number of baseline covariates.
#'
#' @return Invisible TRUE if validation passes, stops with error if validation fails.
#'
#' @details
#' The function validates:
#' \itemize{
#'   \item Data types and dimensions
#'   \item Missing values (warnings)
#'   \item Parameter consistency
#'   \item Matrix dimension compatibility
#' }
#'
#' @examples
#' \dontrun{
#' n <- 100
#' Y <- rnorm(n)
#' Z <- matrix(rnorm(n * 6), nrow = n, ncol = 6)  # 2 metals x 3 time points
#' X <- matrix(rnorm(n * 5), nrow = n, ncol = 5)  # 1 TD cov x 3 time points + 2 baseline
#'
#' validate_user_matrices(Y, Z, X, T = 3, Adim = 2, Ldim = 1, n_baseline = 2)
#' }
#'
#' @keywords internal
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
         "Expected: ", expected_Z_cols, " columns (", Adim, " mixtures x ", T, " time points)\n",
         "Actual: ", ncol(Z), " columns")
  }

  if (ncol(X) != expected_X_cols) {
    stop("X matrix dimension error!\n",
         "Expected: ", expected_X_cols, " columns (", Ldim, " TD covariates x ", T, " time points + ", n_baseline, " baseline covariates)\n",
         "Actual: ", ncol(X), " columns")
  }

  if (any(is.na(Z))) warning("Z contains missing values")
  if (any(is.na(X))) warning("X contains missing values")

  cat("[OK] Input validation passed\n")
  invisible(TRUE)
}
