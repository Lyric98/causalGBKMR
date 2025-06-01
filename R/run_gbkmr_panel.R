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

  if (max(sel) > iter) stop("sel contains indices beyond total MCMC iterations!")

  message("Subsampling population...")
  set.seed(currind)
  dat_sim <- sim_popn[sample(sim_popn$id, n, replace = FALSE), ]

  # Detect p and Ldim
  p <- length(grep("^logM\\d+_0$", names(dat_sim)))
  Ldim <- if (any(grepl("^waist0_\\d+$", names(dat_sim)))) {
    length(grep("^waist0_\\d+$", names(dat_sim)))
  } else 1

  # Common covariates
  cov_names_common <- if (Ldim == 1) c("sex", "waist0") else c("sex", paste0("waist0_", 1:Ldim))
  X_common <- dat_sim[, cov_names_common]

  fitkm_list <- list()
  L_values_a <- list()
  L_values_astar <- list()

  for (t in 1:(T - 1)) {
    message(paste("Fitting mediator model L", t))
    y_L <- if (Ldim == 1) {
      dat_sim[, paste0("waist", t)]
    } else {
      as.matrix(dat_sim[, paste0("waist", t, "_", 1:Ldim)])
    }

    exp_names <- unlist(lapply(0:(t - 1), function(s) paste0("logM", 1:p, "_", s)))
    if (t > 1) {
      for (j in 1:(t - 1)) {
        waist_cols <- if (Ldim == 1) paste0("waist", j) else paste0("waist", j, "_", 1:Ldim)
        exp_names <- c(exp_names, waist_cols)
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
  exp_names_y <- grep("^logM\\d+_\\d+$|^waist\\d+(_\\d+)?$", names(dat_sim), value = TRUE)
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
      unlist(lapply(1:(t - 1), function(s) {
        if (Ldim == 1) paste0("waist", s) else paste0("waist", s, "_", 1:Ldim)
      }))
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

  return(list(
    diff_gBKMR = diff_gBKMR,
    Ya = Ya,
    Yastar = Yastar,
    beta_all = c(
      unlist(lapply(fitkm_list, function(l) colMeans(l$fit$beta[sel, , drop = FALSE]))),
      colMeans(fitkm_y$fit$beta[sel, , drop = FALSE])
    ),
    L_values_a = L_values_a,
    L_values_astar = L_values_astar
  ))
}
