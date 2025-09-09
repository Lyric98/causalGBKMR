#' @file 05-user-interface.R
#' @title User-friendly interface functions for g-BKMR
#' @description High-level functions that provide a clean, user-friendly API
#' for g-BKMR analysis. These functions handle parameter defaults and result formatting.

#' Run g-BKMR analysis with user-friendly interface
#'
#' @description Main user interface for g-BKMR analysis. This function provides
#' a simplified API for conducting g-BKMR analysis with sensible defaults and
#' user-friendly output formatting.
#'
#' @param data Data frame. Input data in g-BKMR format. Can be prepared using
#'   \code{\link{prepare_gbkmr_data}}.
#' @param outcome Character. Name of the outcome variable (default: "Y").
#' @param outcome_type Character. Type of outcome: "continuous" or "binary" (default: "continuous").
#' @param time_points Integer. Number of time points in the study.
#' @param currind Integer. Random seed for reproducibility (default: 1).
#' @param sel Numeric vector. MCMC iterations to use for inference (default: auto-calculated).
#' @param n Integer. Sample size for analysis (default: min(500, nrow(data))).
#' @param K Integer. Number of Monte Carlo samples (default: 1000).
#' @param iter Integer. Total MCMC iterations (default: 15000).
#' @param parallel Logical. Whether to use parallel processing (default: TRUE).
#' @param use_knots Logical. Whether to use kernel approximation (default: TRUE).
#' @param n_knots Integer. Number of knots for approximation (default: 50).
#' @param make_plots Logical. Whether to generate diagnostic plots (default: FALSE).
#' @param verbose Logical. Whether to print detailed progress (default: TRUE).
#'
#' @return An object of class "gbkmr_results" containing:
#' \describe{
#'   \item{causal_effect}{List with estimate, lower CI, and upper CI}
#'   \item{counterfactual_means}{List with low and high exposure means}
#'   \item{variable_importance}{Numeric vector of regression coefficients}
#'   \item{detection_info}{Information about detected variable structure}
#'   \item{raw_results}{Complete results from the core analysis}
#'   \item{call_info}{Information about the function call}
#' }
#'
#' @details
#' This function is the main entry point for users. It:
#' \enumerate{
#'   \item Validates inputs and sets sensible defaults
#'   \item Calls the core g-BKMR analysis function
#'   \item Formats results for user-friendly interpretation
#'   \item Provides confidence intervals and summary statistics
#' }
#'
#' The analysis estimates the causal effect of changing all exposures from their
#' 25th percentile to their 75th percentile, while properly accounting for
#' time-varying confounding.
#'
#' @examples
#' \dontrun{
#' # Prepare your data first
#' prepared_data <- prepare_gbkmr_data(
#'   Y = Y, Z = Z, X = X,
#'   time_points = 3,
#'   mixture_components = 2,
#'   td_covariates = 1,
#'   baseline_covariates = 1
#' )
#'
#' # Run g-BKMR analysis
#' results <- gbkmr_run(
#'   data = prepared_data,
#'   time_points = 3,
#'   verbose = TRUE
#' )
#'
#' # View results
#' print(results)
#' summary(results)
#'
#' # Extract specific components
#' causal_effect <- results$causal_effect$estimate
#' confidence_interval <- c(results$causal_effect$lower, results$causal_effect$upper)
#' }
#'
#' @seealso \code{\link{prepare_gbkmr_data}} for data preparation,
#'   \code{\link{detect_variable_patterns}} for variable detection
#'
#' @export
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
  if (!outcome %in% names(data)) stop("Outcome variable '", outcome, "' not found in data")

  # Set intelligent defaults
  if (is.null(sel)) sel <- seq(iter * 0.6, iter, by = 25)
  if (is.null(n)) n <- min(500, nrow(data))

  if (verbose) {
    cat("Starting g-BKMR analysis...\n")
    cat("Data dimensions:", nrow(data), "subjects,", ncol(data), "variables\n")
    cat("Time points:", time_points, "\n")
    cat("Sample size for analysis:", n, "\n")
    cat("MCMC iterations:", iter, "\n")
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

#' Print method for gbkmr_results objects
#'
#' @param x A gbkmr_results object from \code{\link{gbkmr_run}}
#' @param ... Additional arguments (not used)
#'
#' @return Invisible x
#' @export
print.gbkmr_results <- function(x, ...) {
  cat("g-BKMR Analysis Results\n")
  cat("======================\n\n")

  cat("Causal Effect Estimate:\n")
  cat("  Estimate:", round(x$causal_effect$estimate, 4), "\n")
  cat("  95% CI: (", round(x$causal_effect$lower, 4), ",",
      round(x$causal_effect$upper, 4), ")\n\n")

  cat("Counterfactual Means:\n")
  cat("  Low exposure (25th percentile):", round(x$counterfactual_means$low, 4), "\n")
  cat("  High exposure (75th percentile):", round(x$counterfactual_means$high, 4), "\n\n")

  cat("Data Structure:\n")
  cat("  Exposures per time point:", x$detection_info$p, "\n")
  cat("  Time-dependent covariates per time point:", x$detection_info$Ldim, "\n")
  if (length(x$detection_info$td_covariate_names) > 0) {
    cat("  TD covariate names:", paste(x$detection_info$td_covariate_names, collapse = ", "), "\n")
  }
  cat("  Time points:", x$call_info$time_points, "\n")
  cat("  Sample size:", x$call_info$sample_size, "\n")

  invisible(x)
}

#' Summary method for gbkmr_results objects
#'
#' @param object A gbkmr_results object from \code{\link{gbkmr_run}}
#' @param ... Additional arguments (not used)
#'
#' @return Invisible object
#' @export
summary.gbkmr_results <- function(object, ...) {
  cat("g-BKMR Analysis Summary\n")
  cat("=======================\n\n")

  # Main results
  print(object)

  # Additional details
  cat("\nAnalysis Details:\n")
  cat("  Outcome type:", object$call_info$outcome_type, "\n")
  cat("  MCMC iterations:", object$call_info$mcmc_iterations, "\n")
  cat("  Detection pattern:", object$detection_info$detected_pattern, "\n")

  if (length(object$variable_importance) > 0) {
    cat("\nVariable Importance (Top 5):\n")
    top_vars <- head(sort(abs(object$variable_importance), decreasing = TRUE), 5)
    for (i in 1:length(top_vars)) {
      cat("  ", names(top_vars)[i], ":", round(top_vars[i], 4), "\n")
    }
  }

  # Interpretation
  cat("\nInterpretation:\n")
  effect_size <- abs(object$causal_effect$estimate)
  if (effect_size < 0.1) {
    cat("  Effect size: Small\n")
  } else if (effect_size < 0.5) {
    cat("  Effect size: Medium\n")
  } else {
    cat("  Effect size: Large\n")
  }

  ci_width <- object$causal_effect$upper - object$causal_effect$lower
  cat("  95% CI width:", round(ci_width, 4), "\n")

  # Statistical significance (rough approximation)
  if (sign(object$causal_effect$lower) == sign(object$causal_effect$upper)) {
    cat("  95% CI excludes zero: Yes (statistically significant)\n")
  } else {
    cat("  95% CI excludes zero: No (not statistically significant)\n")
  }

  invisible(object)
}
