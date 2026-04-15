#' @file 05-user-interface.R
#' @title User-friendly interface for g-BKMR
#' @importFrom stats quantile

#' Run g-BKMR analysis
#'
#' @param data Data frame in g-BKMR format (see \code{\link{prepare_gbkmr_data}}).
#' @param outcome Character. Outcome variable name (default: "Y").
#' @param outcome_type Character. "continuous" or "binary".
#' @param time_points Integer. Number of time points.
#' @param currind Integer. Random seed.
#' @param sel Numeric vector. Post-burn-in MCMC indices (auto-calculated if NULL).
#' @param n Integer. Sample size (default: min(500, nrow(data))).
#' @param K Integer. Monte Carlo samples.
#' @param iter Integer. Total MCMC iterations.
#' @param n_knots Integer. Knots for kernel approximation.
#' @param engine Character. "bkmr" or "fastbkmr".
#' @param n_subset Integer. Number of subsets for fastBKMR.
#' @param n_cores Integer. Number of cores for fastBKMR.
#' @param a_probs Numeric vector of length 2. Quantile probabilities for
#'   intervention levels (default: c(0.25, 0.75)).
#' @param a_vals Named numeric vector or NULL. Custom low-exposure values.
#' @param astar_vals Named numeric vector or NULL. Custom high-exposure values.
#' @param verbose Logical. Print progress.
#'
#' @return Object of class "gbkmr_results".
#' @export
gbkmr_run <- function(
    data,
    outcome = "Y",
    outcome_type = c("continuous", "binary"),
    time_points,
    currind = 1,
    sel = NULL,
    n = NULL,
    K = 1000,
    iter = 15000,
    n_knots = 50,
    engine = c("bkmr", "fastbkmr"),
    n_subset = 10,
    n_cores = 10,
    a_probs = c(0.25, 0.75),
    a_vals = NULL,
    astar_vals = NULL,
    verbose = TRUE
) {

  outcome_type <- match.arg(outcome_type)
  engine <- match.arg(engine)
  if (missing(data)) stop("Data must be provided")
  if (missing(time_points)) stop("Number of time points must be provided")
  if (!outcome %in% names(data)) stop("Outcome variable '", outcome, "' not found in data")

  if (is.null(sel)) sel <- seq(floor(iter * 0.6), iter, by = 25)
  if (is.null(n)) n <- min(500, nrow(data))

  if (verbose) {
    cat("Starting g-BKMR analysis...\n")
    cat("  Engine:", engine, "| n:", n, "| iter:", iter,
        "| T:", time_points, "\n")
    if (engine == "fastbkmr")
      cat("  fastBKMR: n_subset=", n_subset, ", n_cores=", n_cores, "\n")
    cat("  Intervention: a_probs =", a_probs[1], "/", a_probs[2], "\n")
  }

  # Detect variable structure
  detection <- detect_variable_patterns(data, time_points)

  results <- run_gbkmr_panel(
    sim_popn = data,
    T = time_points,
    p = detection$p,
    mediator_basenames = detection$td_covariate_names,
    common_covariates = c("sex", detection$baseline_td_vars),
    currind = currind,
    sel = sel,
    n = n,
    K = K,
    iter = iter,
    n_knots = n_knots,
    engine = engine,
    n_subset = n_subset,
    n_cores = n_cores,
    a_probs = a_probs,
    a_vals = a_vals,
    astar_vals = astar_vals
  )

  # Format output
  causal_effect <- list(
    estimate = results$diff_gBKMR,
    lower = quantile(results$Yastar - results$Ya, 0.025, na.rm = TRUE),
    upper = quantile(results$Yastar - results$Ya, 0.975, na.rm = TRUE)
  )

  formatted_results <- list(
    causal_effect = causal_effect,
    counterfactual_means = list(
      low  = mean(results$Ya, na.rm = TRUE),
      high = mean(results$Yastar, na.rm = TRUE)
    ),
    variable_importance = results$beta_y,
    detection_info = detection,
    raw_results = results,
    call_info = list(
      outcome = outcome,
      outcome_type = outcome_type,
      time_points = time_points,
      sample_size = n,
      mcmc_iterations = iter,
      engine = engine,
      a_probs = a_probs
    )
  )

  class(formatted_results) <- "gbkmr_results"

  if (verbose) {
    cat("Analysis complete!\n")
    cat("  ATE:", round(causal_effect$estimate, 4),
        "  95% CI: (", round(causal_effect$lower, 4), ",",
        round(causal_effect$upper, 4), ")\n")
  }

  formatted_results
}

#' @export
print.gbkmr_results <- function(x, ...) {
  cat("g-BKMR Analysis Results\n")
  cat("======================\n\n")
  cat("Engine:", x$call_info$engine, "\n")
  cat("Causal Effect (ATE):", round(x$causal_effect$estimate, 4), "\n")
  cat("95% CI: (", round(x$causal_effect$lower, 4), ",",
      round(x$causal_effect$upper, 4), ")\n\n")
  cat("Counterfactual Means:\n")
  cat("  E[Y^a]  (low):", round(x$counterfactual_means$low, 4), "\n")
  cat("  E[Y^a*] (high):", round(x$counterfactual_means$high, 4), "\n")
  invisible(x)
}

#' @export
summary.gbkmr_results <- function(object, ...) {
  print(object)
  cat("\nSettings:\n")
  cat("  Time points:", object$call_info$time_points, "\n")
  cat("  Sample size:", object$call_info$sample_size, "\n")
  cat("  MCMC iterations:", object$call_info$mcmc_iterations, "\n")
  cat("  Intervention quantiles:", object$call_info$a_probs[1],
      "vs", object$call_info$a_probs[2], "\n")

  if (sign(object$causal_effect$lower) == sign(object$causal_effect$upper)) {
    cat("\n  95% CI excludes zero: Yes\n")
  } else {
    cat("\n  95% CI excludes zero: No\n")
  }
  invisible(object)
}
