#' @file 04-core-analysis.R
#' @title Core g-BKMR analysis implementation
#' @importFrom stats complete.cases quantile rnorm

#' Run g-BKMR panel analysis
#'
#' @param sim_popn Data frame in g-BKMR format (see \code{\link{prepare_gbkmr_data}}).
#' @param T Integer. Number of time points (including t=0).
#' @param p Integer. Number of exposures per time point.
#' @param mediator_basenames Character vector. Base names for time-dependent confounders.
#' @param common_covariates Character vector. Baseline covariate names.
#' @param currind Integer. Random seed.
#' @param n Integer. Sample size for analysis.
#' @param K Integer. Monte Carlo samples for g-computation.
#' @param sel Numeric vector. Post-burn-in MCMC indices for inference.
#' @param iter Integer. MCMC iterations for mediator models.
#' @param n_iter Integer or NULL. MCMC iterations for outcome model (default: iter).
#' @param n_knots Integer. Number of knots for BKMR kernel approximation.
#' @param engine Character. Fitting engine: "bkmr" or "fastbkmr".
#' @param n_subset Integer. Number of subsets for fastBKMR.
#' @param n_cores Integer. Number of cores for fastBKMR parallel.
#' @param outcome_type Character. "continuous" (Gaussian BKMR, default) or
#'   "binary" (probit BKMR via family="binomial"). Binary outcome requires
#'   engine="bkmr" (fastBKMR does not yet support non-Gaussian outcomes).
#' @param a_probs Numeric vector of length 2. Quantile probabilities for
#'   intervention levels (default: c(0.25, 0.75)).
#' @param a_vals Named numeric vector or NULL. Custom intervention values for
#'   low-exposure scenario. Overrides a_probs if provided.
#' @param astar_vals Named numeric vector or NULL. Custom intervention values for
#'   high-exposure scenario. Overrides a_probs if provided.
#' @param verbose_every Integer. Print progress every N iterations.
#'
#' @return A list with causal effect estimate and model fits.
#' @export
run_gbkmr_panel <- function(
    sim_popn,
    T = 5,
    p = 3,
    mediator_basenames = c("waist"),
    common_covariates = c("sex", "waist0"),
    currind = 1,
    n = 500,
    K = 1000,
    sel = seq(22000, 24000, by = 25),
    iter = 24000,
    n_iter = NULL,
    n_knots = 50,
    engine = c("bkmr", "fastbkmr"),
    n_subset = 10,
    n_cores = 10,
    outcome_type = c("continuous", "binary"),
    a_probs = c(0.25, 0.75),
    a_vals = NULL,
    astar_vals = NULL,
    verbose_every = 50) {

  engine <- match.arg(engine)
  outcome_type <- match.arg(outcome_type)

  # Binary outcome is only supported for standard BKMR; fbkmr::skmbayes()
  # internally calls kmbayes() without a family argument, so it cannot
  # fit probit models in the current version.
  if (outcome_type == "binary" && engine == "fastbkmr") {
    stop("Binary outcome is not supported with engine='fastbkmr'.\n",
         "  fbkmr::skmbayes() does not expose the family argument.\n",
         "  Use engine='bkmr' for binary outcomes, or subsample your data.")
  }

  if (is.null(n_iter)) n_iter <- iter
  if (!"Y" %in% names(sim_popn)) stop("Data must contain outcome variable 'Y'")
  if (!"id" %in% names(sim_popn)) stop("Data must contain 'id' column")
  if (n > nrow(sim_popn)) {
    warning("n larger than data; using all rows.")
    n <- nrow(sim_popn)
  }
  if (max(sel) > max(iter, n_iter)) stop("sel contains indices beyond total MCMC iterations!")

  # --- Internal helpers ---
  .fit_model <- function(y, Z_sc, X, it, knots = NULL, family = "gaussian") {
    if (engine == "fastbkmr") {
      if (!requireNamespace("fbkmr", quietly = TRUE))
        stop("Package 'fbkmr' is required for engine='fastbkmr'.\n",
             "Install with: remotes::install_github('junwei-lu/fbkmr')")
      nc <- min(n_subset, parallel::detectCores() - 1, n_cores)
      use_parallel <- nc > 1
      tryCatch(
        fbkmr::skmbayes(Z = Z_sc, X = X, y = y,
                         n_subset = n_subset, n_samp = 200,
                         iter = it, varsel = TRUE, est.h = FALSE,
                         parallel = use_parallel, n_cores = nc),
        error = function(e) {
          if (use_parallel && grepl("doSnowGlobals|parallel|snow", e$message, ignore.case = TRUE)) {
            warning("Parallel execution failed, falling back to sequential mode.\n",
                    "  Original error: ", e$message, call. = FALSE)
            fbkmr::skmbayes(Z = Z_sc, X = X, y = y,
                             n_subset = n_subset, n_samp = 200,
                             iter = it, varsel = TRUE, est.h = FALSE,
                             parallel = FALSE, n_cores = 1)
          } else {
            stop(e)
          }
        }
      )
    } else {
      bkmr::kmbayes(y = y, Z = Z_sc, X = X, iter = it,
                     family = family,
                     varsel = TRUE, verbose = FALSE, knots = knots)
    }
  }

  .sample_pred <- function(fit, Znew, Xnew, sel_j) {
    if (is.list(fit) && !inherits(fit, "bkmrfit")) {
      bkmr::SamplePred(fit[[sample(length(fit), 1)]],
                        Znew = Znew, Xnew = Xnew, sel = sel_j)
    } else {
      bkmr::SamplePred(fit, Znew = Znew, Xnew = Xnew, sel = sel_j)
    }
  }

  .extract_beta <- function(fit, sel_idx) {
    if (is.list(fit) && !inherits(fit, "bkmrfit")) {
      betas <- lapply(fit, function(f) colMeans(f$beta[sel_idx, , drop = FALSE]))
      colMeans(do.call(rbind, betas))
    } else {
      colMeans(fit$beta[sel_idx, , drop = FALSE])
    }
  }

  # --- Sampling and naming ---
  set.seed(currind)
  dat_sim <- sim_popn[sample(seq_len(nrow(sim_popn)), n, replace = FALSE), ]

  exposure_names_at_t <- function(t) paste0("logM", 1:p, "_", t)
  all_exposure_names <- unlist(lapply(0:(T - 1), exposure_names_at_t))
  mediator_names_at_t <- function(t) paste0(mediator_basenames, "_", t)
  all_mediator_names  <- unlist(lapply(1:(T - 1), mediator_names_at_t))

  needed_cols <- c("Y", "id", common_covariates, all_exposure_names, all_mediator_names)
  miss <- setdiff(needed_cols, names(dat_sim))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  X_common <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(common_covariates)))
  X_predict_common <- matrix(colMeans(X_common), nrow = 1)

  fitkm_list <- vector("list", T - 1)
  names(fitkm_list) <- paste0("L", 1:(T - 1))
  scaleinfo_list <- vector("list", T - 1)
  names(scaleinfo_list) <- names(fitkm_list)

  # --- Knot helper (only used for engine == "bkmr") ---
  .compute_knots <- function(Z_sc, n_knots) {
    n_unique <- nrow(unique(round(Z_sc, 10)))
    if (n_unique < nrow(Z_sc)) {
      Z_sc <- Z_sc + matrix(rnorm(length(Z_sc), 0, 1e-6), nrow = nrow(Z_sc))
      n_unique <- nrow(unique(round(Z_sc, 10)))
    }
    nd <- min(n_knots, n_unique - 1, floor(n_unique * 0.9))
    if (nd < 2) stop("Not enough unique rows to place knots")
    tryCatch(
      fields::cover.design(Z_sc, nd = nd)$design,
      error = function(e) Z_sc[sample(nrow(Z_sc), nd), , drop = FALSE]
    )
  }

  # =========================================================================
  # 1) Fit mediator models
  # =========================================================================
  message("Fitting mediator models ...")

  for (t in 1:(T - 1)) {
    y_cols <- mediator_names_at_t(t)
    y_mat  <- as.matrix(dat_sim[, y_cols, drop = FALSE])
    colnames(y_mat) <- y_cols

    # Z: exposures 0..t-1 + mediators 1..t-1
    Z_names <- unlist(lapply(0:(t - 1), exposure_names_at_t))
    if (t > 1) Z_names <- c(Z_names, unlist(lapply(1:(t - 1), mediator_names_at_t)))

    Z_raw <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(Z_names)))
    rows_ok_ZX <- complete.cases(Z_raw, X_common)
    if (sum(rows_ok_ZX) < 3) stop("Not enough complete rows for mediator at t=", t)

    Z_sc <- scale(Z_raw[rows_ok_ZX, , drop = FALSE])
    sc_center <- attr(Z_sc, "scaled:center")
    sc_scale  <- attr(Z_sc, "scaled:scale")
    scaleinfo_list[[t]] <- list(center = sc_center, scale = sc_scale)

    knots_t <- if (engine == "bkmr") .compute_knots(Z_sc, n_knots) else NULL

    fit_list_t <- vector("list", ncol(y_mat))
    for (li in seq_len(ncol(y_mat))) {
      y_vec <- y_mat[, li]
      y_ok  <- y_vec[rows_ok_ZX]
      mask_y <- !is.na(y_ok)
      if (sum(mask_y) < 3) {
        stop(sprintf("Not enough complete rows for mediator %s at t=%d",
                     colnames(y_mat)[li], t))
      }
      valid_idx <- which(rows_ok_ZX)[mask_y]
      Z_sc_fit <- scale(Z_raw[valid_idx, , drop = FALSE],
                        center = sc_center, scale = sc_scale)
      X_common_fit <- X_common[valid_idx, , drop = FALSE]
      y_vec_fit <- y_vec[valid_idx]

      message(sprintf("  L%d: fitting %s [engine=%s, Z=%d cols, n=%d]",
                      t, colnames(y_mat)[li], engine, ncol(Z_sc_fit), length(y_vec_fit)))

      fit_list_t[[li]] <- .fit_model(y_vec_fit, Z_sc_fit, X_common_fit,
                                      iter, knots_t)
    }
    fitkm_list[[t]] <- fit_list_t
  }

  # =========================================================================
  # 2) Fit outcome model
  # =========================================================================
  message("Fitting outcome model Y ...")
  Y <- dat_sim$Y
  valid_y <- !is.na(Y)
  if (any(!valid_y)) {
    message(sprintf("  Removing %d NA outcomes", sum(!valid_y)))
    dat_sim <- dat_sim[valid_y, ]
    Y <- Y[valid_y]
    X_common <- X_common[valid_y, , drop = FALSE]
    X_predict_common <- matrix(colMeans(X_common), nrow = 1)
  }

  Zy_names <- c(all_exposure_names, all_mediator_names)
  Zy_raw   <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(Zy_names)))
  Zy_sc    <- scale(Zy_raw)
  scale_info_y <- list(center = attr(Zy_sc, "scaled:center"),
                       scale  = attr(Zy_sc, "scaled:scale"))

  knots_y <- if (engine == "bkmr") .compute_knots(Zy_sc, n_knots) else NULL

  y_family <- if (outcome_type == "binary") "binomial" else "gaussian"
  message(sprintf("  [engine=%s, Z=%d cols, n=%d, family=%s]",
                  engine, ncol(Zy_sc), length(Y), y_family))
  fit_y <- .fit_model(Y, Zy_sc, X_common, n_iter, knots_y, family = y_family)

  # =========================================================================
  # 3) Intervention levels (a / a*)
  # =========================================================================
  A_all <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(all_exposure_names)))
  if (!is.null(a_vals) && !is.null(astar_vals)) {
    a_vec     <- a_vals[all_exposure_names]
    astar_vec <- astar_vals[all_exposure_names]
  } else {
    a_vec     <- apply(A_all, 2, quantile, probs = a_probs[1])
    astar_vec <- apply(A_all, 2, quantile, probs = a_probs[2])
  }

  # Containers
  L_samp_a     <- vector("list", T - 1)
  L_samp_astar <- vector("list", T - 1)

  scale_like <- function(newZ, center, sc) scale(newZ, center = center, scale = sc)

  # =========================================================================
  # 4) Sequential mediator sampling
  # =========================================================================
  message("\n=== Sampling mediators sequentially ===")
  start_time_global <- proc.time()

  for (t in 1:(T - 1)) {
    message(sprintf("\n--- Time point t=%d ---", t))
    start_time_t <- proc.time()

    Z_names_t <- unlist(lapply(0:(t - 1), exposure_names_at_t))
    if (t > 1) Z_names_t <- c(Z_names_t, unlist(lapply(1:(t - 1), mediator_names_at_t)))

    a_exp_t     <- unlist(lapply(0:(t - 1), function(s) a_vec[exposure_names_at_t(s)]))
    astar_exp_t <- unlist(lapply(0:(t - 1), function(s) astar_vec[exposure_names_at_t(s)]))

    L_samp_a_t     <- vector("list", length(mediator_basenames))
    L_samp_astar_t <- vector("list", length(mediator_basenames))

    for (li in seq_along(mediator_basenames)) {
      message(sprintf("  Sampling %s at t=%d", mediator_basenames[li], t))

      fit_li   <- fitkm_list[[t]][[li]]
      scinfo_t <- scaleinfo_list[[t]]

      L_a_mat     <- matrix(NA, nrow = length(sel), ncol = K)
      L_astar_mat <- matrix(NA, nrow = length(sel), ncol = K)

      for (j in seq_along(sel)) {
        # Build historical mediator block
        if (t == 1) {
          L_hist_a_j <- NULL
          L_hist_astar_j <- NULL
        } else {
          L_hist_a_blocks <- L_hist_astar_blocks <- list()
          for (tt in 1:(t - 1)) {
            for (lj in seq_along(mediator_basenames)) {
              L_hist_a_blocks[[length(L_hist_a_blocks) + 1]] <- L_samp_a[[tt]][[lj]][j, ]
              L_hist_astar_blocks[[length(L_hist_astar_blocks) + 1]] <- L_samp_astar[[tt]][[lj]][j, ]
            }
          }
          L_hist_a_j     <- do.call(cbind, L_hist_a_blocks)
          L_hist_astar_j <- do.call(cbind, L_hist_astar_blocks)
        }

        # Build full Z: [exposure block] + [historical mediator block]
        Za_exp_mat     <- matrix(a_exp_t,     nrow = K, ncol = length(a_exp_t),     byrow = TRUE)
        Zastar_exp_mat <- matrix(astar_exp_t, nrow = K, ncol = length(astar_exp_t), byrow = TRUE)

        if (is.null(L_hist_a_j)) {
          aL_a_j         <- Za_exp_mat
          astarL_astar_j <- Zastar_exp_mat
        } else {
          aL_a_j         <- cbind(Za_exp_mat,     L_hist_a_j)
          astarL_astar_j <- cbind(Zastar_exp_mat, L_hist_astar_j)
        }

        for (k in 1:K) {
          row_a     <- matrix(aL_a_j[k, ],         nrow = 1)
          row_astar <- matrix(astarL_astar_j[k, ], nrow = 1)
          newz    <- rbind(row_a, row_astar)
          newz_sc <- scale_like(newz, scinfo_t$center, scinfo_t$scale)

          set.seed(j + 10000 + li)
          L_pred <- .sample_pred(fit_li, Znew = newz_sc,
                                 Xnew = X_predict_common, sel_j = sel[j])
          L_a_mat[j, k]     <- L_pred[, "znew1"]
          L_astar_mat[j, k] <- L_pred[, "znew2"]
        }

        if (j %% verbose_every == 0) {
          elapsed <- round((proc.time() - start_time_t)["elapsed"] / 60, 2)
          message(sprintf("    iter %d/%d | %.2f min", j, length(sel), elapsed))
        }
      }

      L_samp_a_t[[li]]     <- L_a_mat
      L_samp_astar_t[[li]] <- L_astar_mat
    }

    L_samp_a[[t]]     <- L_samp_a_t
    L_samp_astar[[t]] <- L_samp_astar_t
  }

  # =========================================================================
  # 5) Sample outcome Y
  # =========================================================================
  message("\n=== Sampling outcome Y ===")
  start_time_y <- proc.time()

  Ya_mat     <- matrix(NA, nrow = length(sel), ncol = K)
  Yastar_mat <- matrix(NA, nrow = length(sel), ncol = K)

  pT     <- length(all_exposure_names)
  Ltotal <- length(all_mediator_names)

  for (j in seq_along(sel)) {
    exp_a_block     <- matrix(a_vec[all_exposure_names],     nrow = K, ncol = pT, byrow = TRUE)
    exp_astar_block <- matrix(astar_vec[all_exposure_names], nrow = K, ncol = pT, byrow = TRUE)

    L_a_blocks <- L_astar_blocks <- list()
    for (t in 1:(T - 1)) {
      for (li in seq_along(mediator_basenames)) {
        L_a_blocks[[length(L_a_blocks) + 1]]         <- L_samp_a[[t]][[li]][j, ]
        L_astar_blocks[[length(L_astar_blocks) + 1]] <- L_samp_astar[[t]][[li]][j, ]
      }
    }
    L_a_mat2     <- do.call(cbind, L_a_blocks)
    L_astar_mat2 <- do.call(cbind, L_astar_blocks)

    aL_a_j         <- cbind(exp_a_block,     L_a_mat2)
    astarL_astar_j <- cbind(exp_astar_block, L_astar_mat2)

    for (k in 1:K) {
      newz    <- rbind(aL_a_j[k, ], astarL_astar_j[k, ])
      newz_sc <- scale_like(newz, scale_info_y$center, scale_info_y$scale)

      set.seed(j + 10000)
      Y_jk <- .sample_pred(fit_y, Znew = newz_sc,
                            Xnew = X_predict_common, sel_j = sel[j])
      # For binary outcome (probit BKMR), SamplePred returns the linear
      # predictor h(Z) + X*beta. Convert to probability scale via Phi().
      if (outcome_type == "binary") {
        Ya_mat[j, k]     <- stats::pnorm(Y_jk[, "znew1"])
        Yastar_mat[j, k] <- stats::pnorm(Y_jk[, "znew2"])
      } else {
        Ya_mat[j, k]     <- Y_jk[, "znew1"]
        Yastar_mat[j, k] <- Y_jk[, "znew2"]
      }
    }

    if (j %% verbose_every == 0) {
      elapsed <- round((proc.time() - start_time_y)["elapsed"] / 60, 2)
      message(sprintf("  iter %d/%d | %.2f min", j, length(sel), elapsed))
    }
  }

  # =========================================================================
  # 6) Aggregate results
  # =========================================================================
  Ya     <- rowMeans(Ya_mat)
  Yastar <- rowMeans(Yastar_mat)
  diff_gBKMR <- mean(Yastar) - mean(Ya)

  beta_L <- lapply(seq_along(fitkm_list), function(t) {
    lapply(seq_along(fitkm_list[[t]]), function(li) {
      .extract_beta(fitkm_list[[t]][[li]], sel)
    })
  })
  beta_y <- .extract_beta(fit_y, sel)

  end_time_global <- proc.time()
  total_time_min <- round((end_time_global - start_time_global)["elapsed"] / 60, 2)
  message(sprintf("\n=== Total time: %s minutes ===", total_time_min))

  list(
    diff_gBKMR = diff_gBKMR,
    Ya = Ya,
    Yastar = Yastar,
    Ya_mat = Ya_mat,
    Yastar_mat = Yastar_mat,
    L_samp_a = L_samp_a,
    L_samp_astar = L_samp_astar,
    fit_mediators = fitkm_list,
    fit_y = fit_y,
    beta_L = beta_L,
    beta_y = beta_y,
    meta = list(
      T = T, p = p,
      mediator_basenames = mediator_basenames,
      common_covariates = common_covariates,
      K = K, sel = sel, iter = iter, n_iter = n_iter,
      n_knots = n_knots, n = n,
      engine = engine, n_subset = n_subset, n_cores = n_cores,
      outcome_type = outcome_type,
      a_probs = a_probs, a_vals = a_vals, astar_vals = astar_vals,
      a_vec = a_vec, astar_vec = astar_vec,
      exposure_names = all_exposure_names,
      mediator_names = all_mediator_names,
      total_time_minutes = total_time_min
    )
  )
}
