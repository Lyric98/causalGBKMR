#' @file 04-core-analysis.R
#' @title Core g-BKMR analysis functions
#' @description Main computational engines for g-BKMR analysis. Contains the core
#' algorithm that performs Bayesian kernel machine regression with time-varying confounding.

#' Run g-BKMR panel analysis
#'
#' @description Performs the core g-BKMR analysis on a prepared dataset. This function
#' implements the g-formula with Bayesian kernel machine regression to estimate
#' causal effects of time-varying exposure mixtures while accounting for time-varying
#' confounding.
#'
#' @param sim_popn Data frame. Prepared data in g-BKMR format.
#' @param T Integer. Number of time points (default: 5).
#' @param currind Integer. Current iteration index for reproducibility (default: 1).
#' @param sel Numeric vector. MCMC iterations to use for inference (default: seq(22000, 24000, by = 25)).
#' @param n Integer. Sample size to use for analysis (default: 500).
#' @param K Integer. Number of samples for Monte Carlo integration (default: 1000).
#' @param iter Integer. Total MCMC iterations (default: 24000).
#' @param parallel Logical. Whether to use parallel processing (default: TRUE).
#' @param save_exposure_preds Logical. Whether to save exposure predictions (default: TRUE).
#' @param return_ci Logical. Whether to return confidence intervals (default: TRUE).
#' @param make_plots Logical. Whether to generate diagnostic plots (default: TRUE).
#' @param use_knots Logical. Whether to use knots for kernel approximation (default: TRUE).
#' @param n_knots Integer. Number of knots to use (default: 50).
#'
#' @return A list containing:
#' \describe{
#'   \item{diff_gBKMR}{Numeric. The main causal effect estimate}
#'   \item{Ya}{Numeric vector. Potential outcomes under low exposure}
#'   \item{Yastar}{Numeric vector. Potential outcomes under high exposure}
#'   \item{beta_all}{Numeric vector. All regression coefficients}
#'   \item{L_values_a}{List. Mediator values under low exposure}
#'   \item{L_values_astar}{List. Mediator values under high exposure}
#'   \item{detection_info}{List. Variable detection information}
#'   \item{fitted_models}{List. Fitted BKMR models (if save_exposure_preds = TRUE)}
#' }
#'
#' @details
#' The analysis proceeds in several steps:
#' \enumerate{
#'   \item Auto-detect variable patterns using \code{\link{detect_variable_patterns}}
#'   \item Fit mediator models for each time point using BKMR
#'   \item Fit outcome model using BKMR
#'   \item Compute counterfactual predictions using the g-formula
#'   \item Calculate causal effect estimates
#' }
#'
#' The function compares outcomes when all exposures are set to their 25th percentile
#' versus when all exposures are set to their 75th percentile.
#'
#' @examples
#' \dontrun{
#' # Assume you have prepared data
#' library(bkmr)
#'
#' # Run analysis with default settings
#' results <- run_gbkmr_panel(prepared_data, T = 3)
#'
#' # Extract main results
#' causal_effect <- results$diff_gBKMR
#' detection_info <- results$detection_info
#'
#' # Run with custom settings
#' results <- run_gbkmr_panel(
#'   sim_popn = prepared_data,
#'   T = 3,
#'   n = 300,
#'   iter = 15000,
#'   make_plots = FALSE
#' )
#' }
#'
#' @importFrom bkmr kmbayes SamplePred
#' @importFrom fields cover.design
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot geom_col theme_minimal ggtitle aes
#' @importFrom parallel mclapply
#'
#' @keywords internal
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
  # Libraries are now imported via NAMESPACE (see @importFrom above)

  # Input validation
  if (max(sel) > iter) stop("sel contains indices beyond total MCMC iterations!")
  if (!"Y" %in% names(sim_popn)) stop("Data must contain outcome variable 'Y'")
  if (n > nrow(sim_popn)) {
    warning("n is larger than available data. Using all available data.")
    n <- nrow(sim_popn)
  }

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

  # Fit mediator models for each time point
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

    # Scale predictors
    scale_Z <- scale(Z)
    attr_list <- list(center = attr(scale_Z, "scaled:center"),
                      scale = attr(scale_Z, "scaled:scale"))

    # Set up knots if requested
    knots <- if (use_knots) {
      set.seed(1000)
      fields::cover.design(scale_Z, nd = n_knots)$design
    } else NULL

    # Fit BKMR model
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

  # Fit outcome model
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

  # Compute counterfactual exposure vectors
  message("Computing counterfactual exposure vectors...")
  X.predict <- matrix(colMeans(X_common), nrow = 1)

  # Predict mediators under different exposure scenarios
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

  # Predict final outcome
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

  # Generate plots if requested
  if (make_plots) {
    df_plot <- data.frame(Scenario = c("a", "astar"), Mean = c(mean(Ya), mean(Yastar)))
    print(ggplot(df_plot, aes(x = Scenario, y = Mean)) +
            geom_col(fill = "skyblue") +
            theme_minimal() +
            ggtitle("Counterfactual Means"))
  }

  # Compile results
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
    detection_info = detection_results,
    fitted_models = if (save_exposure_preds) list(
      confounders = fitkm_list,
      outcome = fitkm_y
    ) else NULL
  )

  return(results)
}
