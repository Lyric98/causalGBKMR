#' @file 04-core-analysis.R
#' @title Core g-BKMR analysis implementation
#' @description Main function to run g-BKMR analysis with time-varying mixtures
#' and time-dependent confounders. Implements g-formula with BKMR using sequential
#' mediator sampling and outcome prediction.

#' Run g-BKMR panel analysis
#'
#' @description Performs g-BKMR analysis on longitudinal data with time-varying
#' exposures and confounders. Uses Bayesian Kernel Machine Regression (BKMR) to
#' model both mediator (confounder) and outcome processes, then implements
#' g-computation to estimate causal effects.
#'
#' @param sim_popn Data frame. Must contain prepared g-BKMR data with columns:
#'   \itemize{
#'     \item \strong{Y}: Outcome variable
#'     \item \strong{id}: Subject identifier
#'     \item \strong{Exposure variables}: Named as logM1_0, logM2_0, ..., logMp_T-1
#'     \item \strong{Time-dependent confounders}: Named according to mediator_basenames
#'     \item \strong{Baseline covariates}: As specified in common_covariates
#'   }
#' @param T Integer. Total number of time points (including t=0). Default: 5.
#' @param p Integer. Number of exposure variables per time point. Default: 3.
#' @param mediator_basenames Character vector. Base names for time-dependent
#'   confounders (e.g., c("waist", "bmi")). Variables should be named as
#'   basename0, basename1, etc. Default: c("waist").
#' @param common_covariates Character vector. Names of baseline covariates to
#'   adjust for in all models. Default: c("sex", "waist0").
#' @param currind Integer. Random seed for reproducibility. Default: 1.
#' @param n Integer. Sample size to use for analysis (sampled from sim_popn).
#'   Default: 500.
#' @param K Integer. Number of Monte Carlo samples for g-computation. Higher
#'   values provide more stable estimates. Default: 1000.
#' @param sel Numeric vector. MCMC iteration indices to use for inference
#'   (after burn-in). Default: seq(22000, 24000, by = 25).
#' @param iter Integer. Number of MCMC iterations for mediator models.
#'   Default: 24000.
#' @param n_iter Integer or NULL. Number of MCMC iterations for outcome model.
#'   If NULL, uses same as iter. Default: NULL.
#' @param n_knots Integer. Number of knots for kernel approximation via
#'   cover.design(). Reduce if memory issues occur. Default: 50.
#' @param verbose_every Integer. Print progress every N iterations. Set to
#'   large number to reduce output. Default: 50.
#'
#' @return A list containing:
#' \describe{
#'   \item{diff_gBKMR}{Numeric. Estimated average treatment effect (ATE)}
#'   \item{Ya}{Numeric vector. Posterior means of counterfactual outcome under intervention a (25th percentile)}
#'   \item{Yastar}{Numeric vector. Posterior means of counterfactual outcome under intervention a* (75th percentile)}
#'   \item{Ya_mat}{Matrix (length(sel) × K). Full posterior samples for Y under a}
#'   \item{Yastar_mat}{Matrix (length(sel) × K). Full posterior samples for Y under a*}
#'   \item{L_samp_a}{List. Sampled mediator values under intervention a}
#'   \item{L_samp_astar}{List. Sampled mediator values under intervention a*}
#'   \item{fit_mediators}{List. Fitted BKMR models for each mediator at each time point}
#'   \item{fit_y}{BKMR model object. Fitted outcome model}
#'   \item{beta_L}{List. Posterior mean regression coefficients for mediator models}
#'   \item{beta_y}{Numeric vector. Posterior mean regression coefficients for outcome model}
#'   \item{meta}{List. Metadata including all analysis parameters and timing information}
#' }
#'
#' @details
#' The function implements the following workflow:
#'
#' \strong{1. Model Fitting Phase}
#' \itemize{
#'   \item Fits BKMR models for each time-dependent confounder at each time point
#'   \item Each mediator model conditions on: past exposures, past mediators, baseline covariates
#'   \item Fits outcome model conditioning on: all exposures, all mediators, baseline covariates
#'   \item Uses kernel approximation (cover.design) for computational efficiency
#' }
#'
#' \strong{2. G-Computation Phase}
#' \itemize{
#'   \item Sequentially samples mediators forward through time under two intervention scenarios:
#'     \itemize{
#'       \item \strong{a}: All exposures set to 25th percentile
#'       \item \strong{a*}: All exposures set to 75th percentile
#'     }
#'   \item For each MCMC iteration j and Monte Carlo sample k:
#'     \itemize{
#'       \item Sample L_t from P(L_t | A_{0:t-1}, L_{1:t-1}, C_0) under interventions
#'       \item Continue sequentially through all time points
#'       \item Sample final outcome from P(Y | A_{0:T-1}, L_{1:T-1}, C_0)
#'     }
#'   \item Averages over K samples and posterior draws to get E[Y^a] and E[Y^{a*}]
#' }
#'
#' \strong{3. Causal Effect Estimation}
#' \itemize{
#'   \item Computes ATE = E[Y^{a*}] - E[Y^a]
#'   \item Provides full posterior distribution for uncertainty quantification
#' }
#'
#' @section Computational Considerations:
#' \itemize{
#'   \item Total computation time scales with: T, p, |mediator_basenames|, length(sel), K, iter
#'   \item Memory usage scales with: n, n_knots, length(sel) × K
#'   \item For large datasets: reduce n, n_knots, or length(sel)
#'   \item For faster prototyping: reduce iter, length(sel), or K
#' }
#'
#' @section Data Requirements:
#' The input data frame \code{sim_popn} must follow g-BKMR naming conventions:
#' \itemize{
#'   \item Exposures: logM1_0, logM2_0, ..., logMp_0, logM1_1, ..., logMp_{T-1}
#'   \item Mediators: waist0, waist1, ..., waist_{T-1} (if mediator_basenames = "waist")
#'   \item For multiple mediators: waist0, bmi0, waist1, bmi1, etc.
#'   \item Baseline covariates: sex, age, or as specified in common_covariates
#'   \item Outcome: Y (continuous or binary)
#'   \item ID: id (subject identifier)
#' }
#'
#' @examples
#' \dontrun{
#' # Prepare your data first using prepare_gbkmr_data()
#' # Assume data is already in correct format
#'
#' # Basic usage with defaults
#' results <- run_gbkmr_panel(
#'   sim_popn = my_prepared_data,
#'   T = 3,
#'   p = 2,
#'   mediator_basenames = c("waist")
#' )
#'
#' # Access results
#' cat("Estimated ATE:", results$diff_gBKMR, "\n")
#' cat("95% CI:", quantile(results$Yastar_mat - results$Ya_mat, c(0.025, 0.975)), "\n")
#'
#' # Advanced usage with custom parameters
#' results_advanced <- run_gbkmr_panel(
#'   sim_popn = my_prepared_data,
#'   T = 4,
#'   p = 3,
#'   mediator_basenames = c("waist", "bmi"),
#'   common_covariates = c("sex", "age", "waist0", "bmi0"),
#'   n = 600,
#'   K = 2000,
#'   sel = seq(30000, 40000, by = 50),
#'   iter = 40000,
#'   n_knots = 60,
#'   verbose_every = 100
#' )
#' }
#'
#' @references
#' Bobb, J.F., Valeri, L., Claus Henn, B., et al. (2015). Bayesian kernel machine
#' regression for estimating the health effects of multi-pollutant mixtures.
#' \emph{Biostatistics}, 16(3), 493-508.
#'
#' Robins, J. (1986). A new approach to causal inference in mortality studies with
#' sustained exposure periods. \emph{Mathematical Modelling}, 7(9-12), 1393-1512.
#'
#' @seealso
#' \code{\link{prepare_gbkmr_data}} for data preparation
#' \code{\link{detect_variable_patterns}} for automatic variable detection
#' \code{\link[bkmr]{kmbayes}} for BKMR model fitting
#'
#' @export

run_gbkmr_panel <- function(
    sim_popn,
    T = 5,                               #
    p = 3,                               #
    mediator_basenames = c("waist"),     #
    common_covariates = c("sex", "waist0"),
    currind = 1,
    n = 500,
    K = 1000,
    sel = seq(22000, 24000, by = 25),
    iter = 24000,                        # 中介模型迭代
    n_iter = NULL,                       # 结局模型迭代（默认用 iter）
    n_knots = 50,                        # cover.design 的 nd
    verbose_every = 50) {

  if (is.null(n_iter)) n_iter <- iter
  if (!"Y" %in% names(sim_popn)) stop("Data must contain outcome variable 'Y'")
  if (!"id" %in% names(sim_popn)) stop("Data must contain 'id' column")
  if (n > nrow(sim_popn)) {
    warning("n larger than data; using all rows.")
    n <- nrow(sim_popn)
  }
  if (max(sel) > max(iter, n_iter)) stop("sel contains indices beyond total MCMC iterations!")

  set.seed(currind)
  #dat_sim <- sim_popn[sample(sim_popn$id, n, replace = FALSE), ]
  dat_sim <- sim_popn[sample(seq_len(nrow(sim_popn)), n, replace = FALSE), ]

  # ===== MODIFICATION 2: Add underscore to mediator names =====
  # 名称工具
  exposure_names_at_t <- function(t) paste0("logM", 1:p, "_", t)
  all_exposure_names <- unlist(lapply(0:(T - 1), exposure_names_at_t))
  mediator_names_at_t <- function(t) paste0(mediator_basenames, "_", t)  # ADDED UNDERSCORE
  all_mediator_names  <- unlist(lapply(1:(T - 1), mediator_names_at_t))
  # =======================================================

  # # 名称工具
  # exposure_names_at_t <- function(t) paste0("logM", 1:p, "_", t)
  # all_exposure_names <- unlist(lapply(0:(T - 1), exposure_names_at_t))
  # mediator_names_at_t <- function(t) paste0(mediator_basenames, t)
  # all_mediator_names  <- unlist(lapply(1:(T - 1), mediator_names_at_t))

  # 校验列
  needed_cols <- c("Y", "id", common_covariates, all_exposure_names, all_mediator_names)
  miss <- setdiff(needed_cols, names(dat_sim))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  X_common <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(common_covariates)))
  X_predict_common <- matrix(colMeans(X_common), nrow = 1)

  fitkm_list <- vector("list", T - 1)
  names(fitkm_list) <- paste0("L", 1:(T - 1))
  scaleinfo_list <- vector("list", T - 1)
  names(scaleinfo_list) <- names(fitkm_list)

  message("Fitting mediator models ...")
  # for (t in 1:(T - 1)) {
  #   y_cols <- mediator_names_at_t(t)
  #   y_mat  <- as.matrix(dat_sim[, y_cols, drop = FALSE])
  #   colnames(y_mat) <- y_cols
  #
  #   # Z: 0..t-1 的暴露 + 1..t-1 的所有中介
  #   Z_names <- unlist(lapply(0:(t - 1), exposure_names_at_t))
  #   if (t > 1) Z_names <- c(Z_names, unlist(lapply(1:(t - 1), mediator_names_at_t)))
  #
  #   Z_raw <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(Z_names)))
  #   Z_sc  <- scale(Z_raw)
  #   sc_center <- attr(Z_sc, "scaled:center")
  #   sc_scale  <- attr(Z_sc, "scaled:scale")
  #   scaleinfo_list[[t]] <- list(center = sc_center, scale = sc_scale)
  #
  #   # knots
  #   nd_t <- min(n_knots, nrow(Z_sc) - 1)
  #   if (nd_t < 2) stop("Not enough rows to place knots for mediator at t=", t)
  #   knots_t <- fields::cover.design(Z_sc, nd = nd_t)$design
  #
  #   fit_list_t <- vector("list", ncol(y_mat))
  #   for (li in seq_len(ncol(y_mat))) {
  #     y_vec <- y_mat[, li]
  #     message(sprintf("  - L%d: fitting %s (Z %d cols, knots %d)",
  #                     t, colnames(y_mat)[li], ncol(Z_sc), nrow(knots_t)))
  #     fit_list_t[[li]] <- bkmr::kmbayes(
  #       y = y_vec,
  #       Z = Z_sc,
  #       X = X_common,
  #       iter = iter,
  #       varsel = TRUE,
  #       verbose = FALSE,
  #       knots = knots_t
  #     )
  #   }
  #   fitkm_list[[t]] <- fit_list_t
  # }

  for (t in 1:(T - 1)) {
    y_cols <- mediator_names_at_t(t)
    y_mat  <- as.matrix(dat_sim[, y_cols, drop = FALSE])
    colnames(y_mat) <- y_cols

    # Z: 0..t-1 的暴露 + 1..t-1 的所有中介
    Z_names <- unlist(lapply(0:(t - 1), exposure_names_at_t))
    if (t > 1) Z_names <- c(Z_names, unlist(lapply(1:(t - 1), mediator_names_at_t)))

    Z_raw <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(Z_names)))

    rows_ok_ZX <- complete.cases(Z_raw, X_common)
    if (sum(rows_ok_ZX) < 3) {
      stop("Not enough complete Z/X rows to fit mediator model at t=", t)
    }

    Z_sc  <- scale(Z_raw[rows_ok_ZX, , drop = FALSE])

    # ===== MODIFICATION 1-A: Check duplicates after scaling =====
    n_unique <- nrow(unique(round(Z_sc, 10)))
    if (n_unique < nrow(Z_sc)) {
      Z_sc <- Z_sc + matrix(rnorm(length(Z_sc), 0, 1e-6), nrow = nrow(Z_sc))
      n_unique <- nrow(unique(round(Z_sc, 10)))
    }
    # ===========================================================

    sc_center <- attr(Z_sc, "scaled:center")
    sc_scale  <- attr(Z_sc, "scaled:scale")
    scaleinfo_list[[t]] <- list(center = sc_center, scale = sc_scale)

    # ===== MODIFICATION 1-B: Safe knots calculation =====
    nd_t <- min(n_knots, n_unique - 1, floor(n_unique * 0.9))
    if (nd_t < 2) stop("Not enough rows to place knots for mediator at t=", t)

    knots_t <- tryCatch({
      fields::cover.design(Z_sc, nd = nd_t)$design
    }, error = function(e) {
      message("  cover.design failed, using random subset")
      Z_sc[sample(nrow(Z_sc), nd_t), , drop = FALSE]
    })

    # ===========================================================

    # 真正拟合每个中介并存到 fitkm_list
    fit_list_t <- vector("list", ncol(y_mat))
    for (li in seq_len(ncol(y_mat))) {
      y_vec <- y_mat[, li]

      # ==== FIX-B: 针对该中介，再把 y 的 NA 过滤掉 ====
      # 先限制到 Z/X 可用行，再去掉 y 的 NA
      y_ok    <- y_vec[rows_ok_ZX]
      mask_y  <- !is.na(y_ok)
      if (sum(mask_y) < 3) {
        stop(sprintf("Not enough complete rows for mediator %s at t=%d",
                     colnames(y_mat)[li], t))
      }
      # 训练用的行索引（相对于原始 dat_sim）
      valid_idx <- which(rows_ok_ZX)[mask_y]

      # 用同一套 center/scale 标准化训练用 Z
      Z_sc_fit <- scale(Z_raw[valid_idx, , drop = FALSE],
                        center = sc_center, scale = sc_scale)
      X_common_fit <- X_common[valid_idx, , drop = FALSE]
      y_vec_fit    <- y_vec[valid_idx]
      # ===============================================

      message(sprintf("  - L%d: fitting %s (Z %d cols, knots %d, n=%d)",
                      t, colnames(y_mat)[li], ncol(Z_sc_fit), nrow(knots_t), length(y_vec_fit)))

      fit_list_t[[li]] <- bkmr::kmbayes(
        y = y_vec_fit,
        Z = Z_sc_fit,
        X = X_common_fit,
        iter = iter,
        varsel = TRUE,
        verbose = FALSE,
        knots = knots_t
      )
    }
    fitkm_list[[t]] <- fit_list_t
    # ===========================================================
    }
    # ===========================================================

  # =========================
  # 2) 结局模型 Y（缩放 + cover.design）
  # =========================
  message("Fitting outcome model Y ...")
  Y <- dat_sim$Y
  # ===== MODIFICATION 3: Remove NA in outcome =====
  valid_idx <- !is.na(Y)
  if (sum(!valid_idx) > 0) {
    message(sprintf("  Removing %d observations with NA in outcome Y", sum(!valid_idx)))
    dat_sim <- dat_sim[valid_idx, ]
    Y <- Y[valid_idx]
    X_common <- X_common[valid_idx, , drop = FALSE]
    X_predict_common <- matrix(colMeans(X_common), nrow = 1)
  }
  # =======================================================


  Zy_names <- c(all_exposure_names, all_mediator_names)
  Zy_raw   <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(Zy_names)))
  Zy_sc    <- scale(Zy_raw)

########################
  n_unique_y <- nrow(unique(round(Zy_sc, 10)))
  if (n_unique_y < nrow(Zy_sc)) {
    Zy_sc <- Zy_sc + matrix(rnorm(length(Zy_sc), 0, 1e-6), nrow = nrow(Zy_sc))
    n_unique_y <- nrow(unique(round(Zy_sc, 10)))
  }
########################
  scale_info_y <- list(
    center = attr(Zy_sc, "scaled:center"),
    scale  = attr(Zy_sc, "scaled:scale")
  )

  # nd_y <- min(n_knots, nrow(Zy_sc) - 1)
  # if (nd_y < 2) stop("Not enough rows to place knots for outcome")
  # knots_y <- fields::cover.design(Zy_sc, nd = nd_y)$design

  # ===== MODIFICATION 2-B: Safe knots calculation =====
  nd_y <- min(n_knots, n_unique_y - 1, floor(n_unique_y * 0.9))
  if (nd_y < 2) stop("Not enough rows to place knots for outcome")

  knots_y <- tryCatch({
    fields::cover.design(Zy_sc, nd = nd_y)$design
  }, error = function(e) {
    message("cover.design failed for outcome, using random subset")
    Zy_sc[sample(nrow(Zy_sc), nd_y), , drop = FALSE]
  })
  # =======================================================

  fit_y <- bkmr::kmbayes(
    y = Y,
    Z = Zy_sc,
    X = X_common,
    iter = n_iter,
    varsel = TRUE,
    verbose = FALSE,
    knots = knots_y
  )

  # =========================
  # 3) a / astar（分位数）
  # =========================
  A_all <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(all_exposure_names)))
  a_vec     <- apply(A_all, 2, quantile, probs = 0.25)
  astar_vec <- apply(A_all, 2, quantile, probs = 0.75)

  # 容器：每期每个 L 的 (length(sel) × K) 采样矩阵
  L_samp_a     <- vector("list", T - 1)
  L_samp_astar <- vector("list", T - 1)

  scale_like <- function(newZ, center, scale) {
    scale(newZ, center = center, scale = scale)
  }

  # =========================
  # 4) 顺序采样各期中介（逐期累积历史）
  # =========================
  message("\n=== Sampling mediators sequentially ===")
  start_time_global <- proc.time()

  for (t in 1:(T - 1)) {
    message(sprintf("\n--- Time point t=%d ---", t))
    start_time_t <- proc.time()

    # 训练时 Z 的列顺序：先暴露（至 t-1），再历史所有中介（至 t-1）
    Z_names_t <- unlist(lapply(0:(t - 1), exposure_names_at_t))
    if (t > 1) Z_names_t <- c(Z_names_t, unlist(lapply(1:(t - 1), mediator_names_at_t)))

    # 对应 a / astar 的暴露块（长度 t*p）
    a_exp_t     <- unlist(lapply(0:(t - 1), function(s) a_vec[exposure_names_at_t(s)]))
    astar_exp_t <- unlist(lapply(0:(t - 1), function(s) astar_vec[exposure_names_at_t(s)]))

    # 该期每个 L 的容器
    L_samp_a_t     <- vector("list", length(mediator_basenames))
    L_samp_astar_t <- vector("list", length(mediator_basenames))

    for (li in seq_along(mediator_basenames)) {
      message(sprintf("  Sampling %s at t=%d", mediator_basenames[li], t))

      fit_li   <- fitkm_list[[t]][[li]]
      scinfo_t <- scaleinfo_list[[t]]

      # 初始化本 L 的采样矩阵
      L_a_mat     <- matrix(NA, nrow = length(sel), ncol = K)
      L_astar_mat <- matrix(NA, nrow = length(sel), ncol = K)

      # MCMC 循环：遍历 sel 中的每个迭代
      for (j in seq_along(sel)) {
        # 提取历史中介在第 j 次迭代的值
        if (t == 1) {
          # t=1 时无历史中介
          L_hist_a_j     <- NULL
          L_hist_astar_j <- NULL
        } else {
          # t>1: 需要从 L_samp_a[[1:(t-1)]] 中提取第 j 行
          L_hist_a_blocks     <- list()
          L_hist_astar_blocks <- list()
          for (tt in 1:(t - 1)) {
            for (lj in seq_along(mediator_basenames)) {
              # L_samp_a[[tt]][[lj]] 是 (length(sel) × K) 矩阵
              L_hist_a_blocks[[length(L_hist_a_blocks) + 1]]         <- L_samp_a[[tt]][[lj]][j, ]
              L_hist_astar_blocks[[length(L_hist_astar_blocks) + 1]] <- L_samp_astar[[tt]][[lj]][j, ]
            }
          }
          # 转为 K × n_hist_meds 矩阵
          L_hist_a_j     <- do.call(cbind, L_hist_a_blocks)
          L_hist_astar_j <- do.call(cbind, L_hist_astar_blocks)
        }

        # 构建完整的 Z：[暴露块] + [历史中介块]
        # 暴露块：K × (t*p)
        Za_exp_mat     <- matrix(a_exp_t,     nrow = K, ncol = length(a_exp_t),     byrow = TRUE)
        Zastar_exp_mat <- matrix(astar_exp_t, nrow = K, ncol = length(astar_exp_t), byrow = TRUE)

        if (is.null(L_hist_a_j)) {
          aL_a_j         <- Za_exp_mat
          astarL_astar_j <- Zastar_exp_mat
        } else {
          aL_a_j         <- cbind(Za_exp_mat,     L_hist_a_j)
          astarL_astar_j <- cbind(Zastar_exp_mat, L_hist_astar_j)
        }

        # 采样 K 个新值
        for (k in 1:K) {
          #newz <- rbind(aL_a_j[k, ], astarL_astar_j[k, ])
          #=====================Nov9===============
          # Explicitly maintain matrix structure
          row_a_vec <- aL_a_j[k, ]
          row_astar_vec <- astarL_astar_j[k, ]

          row_a <- matrix(row_a_vec, nrow = 1, ncol = length(row_a_vec))
          row_astar <- matrix(row_astar_vec, nrow = 1, ncol = length(row_astar_vec))

          newz <- rbind(row_a, row_astar)
          # newz is guaranteed to be a 2×n matrix
          #============================================
          newz_sc <- scale_like(newz, scinfo_t$center, scinfo_t$scale)

          set.seed(j + 10000 + li)
          L_pred <- bkmr::SamplePred(
            fit_li,
            Znew = newz_sc,
            Xnew = X_predict_common,
            sel = sel[j]
          )
          L_a_mat[j, k]     <- L_pred[, "znew1"]
          L_astar_mat[j, k] <- L_pred[, "znew2"]
        }

        # 进度提示
        if (j %% verbose_every == 0) {
          end_time_temp <- proc.time()
          elapsed_min <- round((end_time_temp - start_time_t)["elapsed"] / 60, 2)
          message(sprintf("    iter %d/%d | time: %s min", j, length(sel), elapsed_min))
        }
      }

      L_samp_a_t[[li]]     <- L_a_mat
      L_samp_astar_t[[li]] <- L_astar_mat
    }

    L_samp_a[[t]]     <- L_samp_a_t
    L_samp_astar[[t]] <- L_samp_astar_t
  }

  # =========================
  # 5) 采样结局 Y
  # =========================
  message("\n=== Sampling outcome Y ===")
  start_time_y <- proc.time()

  Ya_mat     <- matrix(NA, nrow = length(sel), ncol = K)
  Yastar_mat <- matrix(NA, nrow = length(sel), ncol = K)

  # 暴露块维度
  pT <- length(all_exposure_names)              # p * T
  Ltotal <- length(all_mediator_names)          # (T-1) * |mediator_basenames|

  for (j in seq_along(sel)) {
    # 构建暴露块：K × pT
    exp_a_block     <- matrix(a_vec[all_exposure_names],     nrow = K, ncol = pT, byrow = TRUE)
    exp_astar_block <- matrix(astar_vec[all_exposure_names], nrow = K, ncol = pT, byrow = TRUE)

    # 构建中介块：顺序为 t=1..T-1，每期按 mediator_basenames 顺序
    L_a_blocks     <- list()
    L_astar_blocks <- list()
    for (t in 1:(T - 1)) {
      for (li in seq_along(mediator_basenames)) {
        L_a_blocks[[length(L_a_blocks) + 1]]         <- L_samp_a[[t]][[li]][j, ]
        L_astar_blocks[[length(L_astar_blocks) + 1]] <- L_samp_astar[[t]][[li]][j, ]
      }
    }
    L_a_mat     <- do.call(cbind, L_a_blocks)         # K × Ltotal
    L_astar_mat <- do.call(cbind, L_astar_blocks)

    # 组合：[暴露] + [中介]
    aL_a_j         <- cbind(exp_a_block,     L_a_mat)
    astarL_astar_j <- cbind(exp_astar_block, L_astar_mat)

    # 采样 K 个 Y
    for (k in 1:K) {
      newz <- rbind(aL_a_j[k, ], astarL_astar_j[k, ])
      newz_sc <- scale_like(newz, scale_info_y$center, scale_info_y$scale)

      set.seed(j + 10000)
      Y_jk <- bkmr::SamplePred(
        fit_y,
        Znew = newz_sc,
        Xnew = X_predict_common,
        sel = sel[j]
      )
      Ya_mat[j, k]     <- Y_jk[, "znew1"]
      Yastar_mat[j, k] <- Y_jk[, "znew2"]
    }

    # 进度提示
    if (j %% verbose_every == 0) {
      end_time_temp <- proc.time()
      elapsed_min <- round((end_time_temp - start_time_y)["elapsed"] / 60, 2)
      message(sprintf("  iter %d/%d | time: %s min", j, length(sel), elapsed_min))
    }
  }

  # =========================
  # 6) 汇总结果
  # =========================
  Ya     <- rowMeans(Ya_mat)
  Yastar <- rowMeans(Yastar_mat)
  diff_gBKMR <- mean(Yastar) - mean(Ya)

  # 提取 beta 系数
  beta_L <- lapply(seq_along(fitkm_list), function(t) {
    lapply(seq_along(fitkm_list[[t]]), function(li) {
      apply(fitkm_list[[t]][[li]]$beta[sel, , drop = FALSE], 2, mean)
    })
  })
  beta_y <- apply(fit_y$beta[sel, , drop = FALSE], 2, mean)

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
      T = T,
      p = p,
      mediator_basenames = mediator_basenames,
      common_covariates = common_covariates,
      K = K,
      sel = sel,
      iter = iter,
      n_iter = n_iter,
      n_knots = n_knots,
      n = n,
      exposure_names = all_exposure_names,
      mediator_names = all_mediator_names,
      total_time_minutes = total_time_min
    )
  )
}
