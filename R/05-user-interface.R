#' @file 05-user-interface.R
#' @title User-friendly interface for g-BKMR
#' @importFrom stats quantile

#' Detect whether a vector is binary (0/1) or continuous
#'
#' @description Returns "binary" if the variable takes only 0/1 (or TRUE/FALSE)
#'   values, "continuous" otherwise.
#' @param x A numeric or logical vector.
#' @return Character string: "binary" or "continuous".
#' @keywords internal
detect_variable_type <- function(x) {
  if (is.logical(x)) return("binary")
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return("continuous")
  unique_vals <- unique(x_clean)
  if (length(unique_vals) == 2L && all(unique_vals %in% c(0, 1))) {
    return("binary")
  }
  "continuous"
}

#' Quick convergence diagnostics for a fitted BKMR model
#'
#' @description Computes two standard MCMC diagnostics for the beta parameters
#'   of a BKMR fit: effective sample size (ESS) and Geweke z-score. Returns
#'   a data.frame of per-parameter diagnostics plus an overall flag.
#' @param fit A bkmrfit object, or a list of bkmrfit objects (fastBKMR).
#' @param sel_idx Integer vector of post-burn-in indices to use.
#' @return A list with: `ess` (min effective sample size across betas),
#'   `geweke_max_abs` (max absolute Geweke z-score across betas),
#'   `warning_flags` (character vector of issues found, or NULL if OK).
#' @keywords internal
check_convergence <- function(fit, sel_idx) {
  if (is.list(fit) && !inherits(fit, "bkmrfit")) {
    # fastBKMR: check each subset, return worst-case
    results_list <- lapply(fit, check_convergence, sel_idx = sel_idx)
    min_ess <- min(sapply(results_list, function(r) r$ess), na.rm = TRUE)
    max_gew <- max(sapply(results_list, function(r) r$geweke_max_abs), na.rm = TRUE)
    all_warnings <- unique(unlist(lapply(results_list, function(r) r$warning_flags)))
    return(list(ess = min_ess, geweke_max_abs = max_gew,
                warning_flags = all_warnings))
  }

  beta_post <- fit$beta[sel_idx, , drop = FALSE]
  n_samp <- nrow(beta_post)
  flags <- character(0)

  # Effective sample size via simple autocorrelation
  ess_per_col <- apply(beta_post, 2, function(x) {
    if (stats::sd(x) < 1e-10) return(n_samp)  # constant column
    acf_vals <- tryCatch(
      stats::acf(x, lag.max = min(50, n_samp - 1), plot = FALSE)$acf[-1],
      error = function(e) 0
    )
    # Truncate at first negative autocorr
    first_neg <- which(acf_vals <= 0)[1]
    if (is.na(first_neg)) first_neg <- length(acf_vals)
    rho_sum <- if (first_neg > 1) sum(acf_vals[1:(first_neg - 1)]) else 0
    n_samp / (1 + 2 * rho_sum)
  })
  min_ess <- min(ess_per_col, na.rm = TRUE)

  # Geweke: compare first 10% to last 50% of chain
  geweke_per_col <- apply(beta_post, 2, function(x) {
    n1 <- max(2, floor(n_samp * 0.1))
    n2 <- max(2, floor(n_samp * 0.5))
    if (n_samp < n1 + n2) return(0)
    x1 <- x[1:n1]
    x2 <- x[(n_samp - n2 + 1):n_samp]
    s1 <- stats::var(x1); s2 <- stats::var(x2)
    if (!is.finite(s1) || !is.finite(s2)) return(0)
    denom <- sqrt(s1 / n1 + s2 / n2)
    if (!is.finite(denom) || denom < 1e-10) return(0)
    (mean(x1) - mean(x2)) / denom
  })
  max_abs_geweke <- suppressWarnings(max(abs(geweke_per_col), na.rm = TRUE))
  if (!is.finite(max_abs_geweke)) max_abs_geweke <- 0

  if (min_ess < 100) flags <- c(flags, "low_ess")
  if (max_abs_geweke > 2.5) flags <- c(flags, "geweke_nonstationary")

  list(ess = min_ess, geweke_max_abs = max_abs_geweke,
       warning_flags = flags)
}

#' Print a summary of what the package detected from user input
#'
#' @description Shows the user exactly which variables the package identified
#'   and how it plans to analyze them. Helps users catch misspecifications
#'   before launching a long MCMC run.
#' @keywords internal
print_input_audit <- function(data, outcome, outcome_type, time_points,
                               detection, engine, n, iter, K, n_knots,
                               a_probs, n_subset = NULL, n_cores = NULL) {
  cat("\n")
  cat("==============================================================\n")
  cat("  causalBKMR: Input Audit\n")
  cat("==============================================================\n")

  # --- Outcome ---
  y_type_detected <- detect_variable_type(data[[outcome]])
  cat("\n[Outcome]\n")
  cat(sprintf("  Variable:        %s\n", outcome))
  cat(sprintf("  Detected type:   %s\n", y_type_detected))
  cat(sprintf("  User specified:  %s\n", outcome_type))
  if (y_type_detected != outcome_type) {
    cat("  WARNING: detected type differs from user-specified type.\n")
    cat("           Proceeding with user-specified type. Check your data!\n")
  }
  if (outcome_type == "binary") {
    cat("  Link function:   probit (via bkmr::kmbayes family='binomial')\n")
    cat("  NOTE: Only engine='bkmr' supports binary outcomes.\n")
    cat("        Output is returned on the probability scale.\n")
  } else {
    cat("  Link function:   identity (Gaussian BKMR)\n")
  }

  # --- Exposures (mixture) ---
  cat("\n[Mixture exposures]\n")
  cat(sprintf("  Components per time (p): %d\n", detection$p))
  cat(sprintf("  Time points (T):         %d\n", time_points))
  cat(sprintf("  Total exposure columns:  %d\n", detection$p * time_points))

  # --- Time-varying confounders ---
  cat("\n[Time-varying confounders]\n")
  if (detection$Ldim == 0L) {
    cat("  None detected.\n")
  } else {
    cat(sprintf("  Number per time (Ldim):  %d\n", detection$Ldim))
    cat(sprintf("  Names:                   %s\n",
                paste(detection$td_covariate_names, collapse = ", ")))
    for (nm in detection$td_covariate_names) {
      # Detect type using all time-indexed instances (e.g., bmi_0, bmi_1, ...)
      cols <- grep(paste0("^", nm, "_\\d+$"), names(data), value = TRUE)
      if (length(cols) > 0L) {
        t <- detect_variable_type(unlist(data[, cols]))
        cat(sprintf("    - %-20s type: %s\n", nm, t))
        if (t == "binary") {
          cat("      NOTE: binary TD confounder -- BKMR mediator model\n")
          cat("            will treat it as continuous in current version.\n")
        }
      }
    }
  }

  # --- Baseline covariates ---
  cat("\n[Baseline covariates]\n")
  baseline_vars <- c("sex", detection$baseline_td_vars)
  cat(sprintf("  Variables: %s\n", paste(baseline_vars, collapse = ", ")))
  cat("  Role:      linear term X*beta (NOT in kernel h())\n")

  # --- Model configuration ---
  cat("\n[Model configuration]\n")
  cat(sprintf("  Engine:          %s\n", engine))
  if (engine == "fastbkmr" && !is.null(n_subset)) {
    cat(sprintf("  fastBKMR setup:  %d subsets x ~%d obs, %d cores\n",
                n_subset, as.integer(n / n_subset), n_cores))
  }
  if (engine == "bkmr") {
    cat(sprintf("  Knots:           %d (Gaussian process approximation)\n", n_knots))
  }
  cat(sprintf("  Sample size:     %d\n", n))
  cat(sprintf("  MCMC iterations: %d\n", iter))
  cat(sprintf("  MC samples (K):  %d\n", K))

  # --- Intervention contrast ---
  cat("\n[Causal contrast]\n")
  cat(sprintf("  Low exposure (a):   all mixture components at %dth percentile\n",
              round(100 * a_probs[1])))
  cat(sprintf("  High exposure (a*): all mixture components at %dth percentile\n",
              round(100 * a_probs[2])))
  cat("  Target estimand:    ACE = E[Y(a*)] - E[Y(a)]\n")
  cat("  IMPORTANT: Only mixture columns vary between a and a*.\n")
  cat("             Time-varying confounders are sampled via g-computation.\n")
  cat("             Baseline covariates are held at their mean.\n")

  cat("\n==============================================================\n")
  cat("  If anything above is unexpected, stop and check your input.\n")
  cat("==============================================================\n\n")
}

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
#' @param engine Character. "auto" (default), "bkmr", or "fastbkmr". When
#'   "auto", the engine is selected based on sample size: standard BKMR for
#'   n <= 2000, fast BKMR for n > 2000 (if fbkmr is installed).
#' @param n_subset Integer or NULL. Number of subsets for fastBKMR. If NULL
#'   (default), auto-calculated as max(5, floor(n / 1000)).
#' @param n_cores Integer or NULL. Number of parallel cores for fastBKMR.
#'   If NULL (default), auto-calculated as min(n_subset, available cores, 10).
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
    engine = c("auto", "bkmr", "fastbkmr"),
    n_subset = NULL,
    n_cores = NULL,
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

  # --- Auto engine selection ---
  # Threshold: n > 2000 triggers fastBKMR (Sonabend et al. 2024 recommend
  # ~1000 per subset; 2 subsets x 1000 = 2000 minimum for fastBKMR to be
  # meaningful; standard BKMR with knots is still tractable below this).
  if (engine == "auto") {
    if (n > 2000 && requireNamespace("fbkmr", quietly = TRUE)) {
      engine <- "fastbkmr"
      if (verbose) message("Auto-selected engine: fastbkmr (n = ", n, " > 2000)")
    } else {
      if (n > 2000 && verbose) {
        message("Note: n = ", n, " > 2000. Consider installing 'fbkmr' for faster fitting:\n",
                "  remotes::install_github('junwei-lu/fbkmr')")
      }
      engine <- "bkmr"
      if (verbose && n <= 2000) message("Auto-selected engine: bkmr (n = ", n, " <= 2000)")
    }
  }

  # Default n_subset / n_cores when using fastBKMR
  if (engine == "fastbkmr") {
    if (is.null(n_subset)) n_subset <- max(5L, as.integer(n / 1000))
    if (is.null(n_cores))  n_cores  <- min(n_subset, parallel::detectCores() - 1L, 10L)
  } else {
    if (is.null(n_subset)) n_subset <- 10L
    if (is.null(n_cores))  n_cores  <- 10L
  }

  # Detect variable structure
  detection <- detect_variable_patterns(data, time_points)

  # Print full input audit so user knows what the package understood
  if (verbose) {
    print_input_audit(
      data = data, outcome = outcome, outcome_type = outcome_type,
      time_points = time_points, detection = detection,
      engine = engine, n = n, iter = iter, K = K, n_knots = n_knots,
      a_probs = a_probs, n_subset = n_subset, n_cores = n_cores
    )
  }

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
    outcome_type = outcome_type,
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

  # --- Layer 1 convergence diagnostics ---
  # Check each fitted BKMR model for ESS and Geweke stationarity
  diagnostics <- list()
  for (t in seq_along(results$fit_mediators)) {
    for (li in seq_along(results$fit_mediators[[t]])) {
      nm <- paste0("mediator_t", t, "_", li)
      diagnostics[[nm]] <- check_convergence(
        results$fit_mediators[[t]][[li]], sel_idx = sel)
    }
  }
  diagnostics[["outcome_Y"]] <- check_convergence(results$fit_y, sel_idx = sel)

  formatted_results <- list(
    causal_effect = causal_effect,
    counterfactual_means = list(
      low  = mean(results$Ya, na.rm = TRUE),
      high = mean(results$Yastar, na.rm = TRUE)
    ),
    variable_importance = results$beta_y,
    detection_info = detection,
    diagnostics = diagnostics,
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
    print_output_summary(formatted_results, detection)
  }

  formatted_results
}

#' Print a post-analysis summary describing the result
#'
#' @description After MCMC completes, prints what the result means: which
#'   models were run, what scale the ACE is on, what assumptions are
#'   required for the causal interpretation to be valid, plus convergence
#'   diagnostics. Helps users catch mis-specifications they made earlier.
#' @keywords internal
print_output_summary <- function(results, detection) {
  ci <- results$causal_effect
  cf <- results$counterfactual_means
  ct <- results$call_info

  ace_scale <- if (ct$outcome_type == "binary") {
    "risk difference (probability scale, in [-1, 1])"
  } else {
    "same scale as the outcome variable Y"
  }

  cat("\n")
  cat("==============================================================\n")
  cat("  causalBKMR: Results Summary\n")
  cat("==============================================================\n")

  # --- Headline result ---
  cat("\n[Causal effect estimate]\n")
  cat(sprintf("  ACE       = %.4f\n", ci$estimate))
  cat(sprintf("  95%% CrI  = (%.4f, %.4f)\n", ci$lower, ci$upper))
  cat(sprintf("  Scale     : %s\n", ace_scale))
  excl_zero <- sign(ci$lower) == sign(ci$upper) && ci$lower != 0
  cat(sprintf("  95%% CrI excludes zero: %s\n",
              if (excl_zero) "YES" else "NO"))

  # --- Counterfactual means ---
  cat("\n[Counterfactual means]\n")
  cat(sprintf("  E[Y under low exposure  (a) ] = %.4f\n", cf$low))
  cat(sprintf("  E[Y under high exposure (a*)] = %.4f\n", cf$high))
  cat("  (ACE = E[Y(a*)] - E[Y(a)])\n")

  # --- What was actually fit ---
  cat("\n[Models fit]\n")
  T_total <- ct$time_points
  Ldim <- detection$Ldim
  n_med <- Ldim * (T_total - 1)
  cat(sprintf("  %d mediator BKMR model(s):", n_med))
  if (n_med == 0) {
    cat(" none (no time-varying confounders)\n")
  } else {
    cat("\n")
    for (t in 1:(T_total - 1)) {
      for (l in seq_along(detection$td_covariate_names)) {
        cat(sprintf("    - %s at visit t=%d  (Gaussian BKMR)\n",
                    detection$td_covariate_names[l], t))
      }
    }
  }
  outcome_fam <- if (ct$outcome_type == "binary") "probit BKMR (family='binomial')" else "Gaussian BKMR"
  cat(sprintf("  1 outcome BKMR model:\n"))
  cat(sprintf("    - %s (%s)\n", ct$outcome, outcome_fam))
  cat(sprintf("  Engine: %s  |  Sample size: %d  |  MCMC iter: %d\n",
              ct$engine, ct$sample_size, ct$mcmc_iterations))

  # --- Causal assumptions ---
  cat("\n[Required causal assumptions]\n")
  cat("  This ACE has a causal interpretation only if:\n")
  cat("    1. CONSISTENCY     -- observed Y under the assigned exposure\n")
  cat("                          equals the counterfactual Y under that exposure.\n")
  cat("    2. SEQUENTIAL EXCH -- no unmeasured time-varying confounding\n")
  cat("                          given the observed history at each visit.\n")
  cat("    3. POSITIVITY      -- every covariate stratum has a non-zero\n")
  cat("                          probability of each exposure level.\n")
  cat("  Violations of these assumptions yield biased estimates.\n")

  # --- Convergence ---
  cat("\n[Convergence check]\n")
  any_warn <- FALSE
  for (nm in names(results$diagnostics)) {
    d <- results$diagnostics[[nm]]
    flags <- d$warning_flags
    if (length(flags) == 0L) {
      cat(sprintf("  [OK]      %-20s ESS=%.0f, Geweke |z|=%.2f\n",
                  nm, d$ess, d$geweke_max_abs))
    } else {
      any_warn <- TRUE
      cat(sprintf("  [WARNING] %-20s ESS=%.0f, Geweke |z|=%.2f  [%s]\n",
                  nm, d$ess, d$geweke_max_abs, paste(flags, collapse = ", ")))
    }
  }
  if (any_warn) {
    cat("\n  Some models show signs of poor convergence. Consider:\n")
    cat("    - Increasing iter (try 30000 or 50000)\n")
    cat("    - Increasing burn-in: sel = seq(floor(iter*0.8), iter, by=25)\n")
    cat("    - Inspecting trace plots:\n")
    cat("        bkmr::TracePlot(results$raw_results$fit_y)\n")
  }

  cat("\n==============================================================\n\n")

  invisible(NULL)
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
