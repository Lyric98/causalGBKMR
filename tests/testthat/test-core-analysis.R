test_that("gbkmr_run completes with bkmr engine", {
  skip_if_not_installed("bkmr")
  set.seed(42)
  n <- 60
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * 6)) + 0.01, n, 6)
  X <- matrix(rnorm(n * 4), n, 4)

  dat <- prepare_gbkmr_data(Y, Z, X,
    time_points = 3, mixture_components = 2,
    td_covariates = 1, baseline_covariates = 1,
    td_covariate_names = "waist", log_transform_mixtures = TRUE)

  res <- gbkmr_run(
    data = dat, outcome = "Y", outcome_type = "continuous",
    time_points = 3, iter = 200, n = 30, K = 3,
    n_knots = 10, engine = "bkmr", verbose = FALSE)

  expect_s3_class(res, "gbkmr_results")
  expect_true(is.numeric(res$causal_effect$estimate))
  expect_true(is.numeric(res$causal_effect$lower))
  expect_true(is.numeric(res$causal_effect$upper))
  expect_equal(res$call_info$engine, "bkmr")
})

test_that("gbkmr_run completes with fastbkmr engine", {
  skip_if_not_installed("bkmr")
  skip_if_not_installed("fbkmr")
  set.seed(42)
  n <- 80
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * 6)) + 0.01, n, 6)
  X <- matrix(rnorm(n * 4), n, 4)

  dat <- prepare_gbkmr_data(Y, Z, X,
    time_points = 3, mixture_components = 2,
    td_covariates = 1, baseline_covariates = 1,
    td_covariate_names = "waist", log_transform_mixtures = TRUE)

  res <- gbkmr_run(
    data = dat, outcome = "Y", outcome_type = "continuous",
    time_points = 3, iter = 200, n = 60, K = 3,
    engine = "fastbkmr", n_subset = 3, n_cores = 1, verbose = FALSE)

  expect_s3_class(res, "gbkmr_results")
  expect_true(is.numeric(res$causal_effect$estimate))
  expect_equal(res$call_info$engine, "fastbkmr")
})

test_that("run_gbkmr_panel rejects sel beyond iter", {
  skip_if_not_installed("bkmr")
  set.seed(42)
  n <- 40
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * 4)) + 0.01, n, 4)
  X <- matrix(rnorm(n * 3), n, 3)

  dat <- prepare_gbkmr_data(Y, Z, X,
    time_points = 2, mixture_components = 2,
    td_covariates = 1, baseline_covariates = 1,
    td_covariate_names = "waist", log_transform_mixtures = TRUE)

  expect_error(
    run_gbkmr_panel(sim_popn = dat, T = 2, p = 2,
      mediator_basenames = "waist", common_covariates = c("sex", "waist_0"),
      iter = 200, sel = seq(500, 1000, by = 25), n = 30, K = 3),
    "sel contains indices beyond"
  )
})

test_that("binary outcome produces probability-scale ATE", {
  skip_if_not_installed("bkmr")
  set.seed(42)
  n <- 80
  Y <- rbinom(n, 1, 0.4)
  Z <- matrix(abs(rnorm(n * 6)) + 0.01, n, 6)
  X <- matrix(rnorm(n * 4), n, 4)

  dat <- prepare_gbkmr_data(Y, Z, X,
    time_points = 3, mixture_components = 2,
    td_covariates = 1, baseline_covariates = 1,
    td_covariate_names = "waist", log_transform_mixtures = TRUE,
    validate_input = FALSE)

  res <- gbkmr_run(
    data = dat, outcome_type = "binary", time_points = 3,
    iter = 300, n = 40, K = 3, n_knots = 10,
    engine = "bkmr", verbose = FALSE)

  # Counterfactual means must be probabilities in [0, 1]
  expect_gte(res$counterfactual_means$low, 0)
  expect_lte(res$counterfactual_means$low, 1)
  expect_gte(res$counterfactual_means$high, 0)
  expect_lte(res$counterfactual_means$high, 1)
  expect_equal(res$call_info$outcome_type, "binary")
})

test_that("binary outcome with fastbkmr engine errors clearly", {
  skip_if_not_installed("bkmr")
  set.seed(42)
  n <- 60
  Y <- rbinom(n, 1, 0.4)
  Z <- matrix(abs(rnorm(n * 6)) + 0.01, n, 6)
  X <- matrix(rnorm(n * 4), n, 4)

  dat <- prepare_gbkmr_data(Y, Z, X,
    time_points = 3, mixture_components = 2,
    td_covariates = 1, baseline_covariates = 1,
    td_covariate_names = "waist", log_transform_mixtures = TRUE,
    validate_input = FALSE)

  expect_error(
    gbkmr_run(data = dat, outcome_type = "binary",
              time_points = 3, iter = 200, n = 40, K = 3,
              engine = "fastbkmr", verbose = FALSE),
    "Binary outcome is not supported"
  )
})

test_that("convergence diagnostics are returned", {
  skip_if_not_installed("bkmr")
  set.seed(42)
  n <- 60
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * 6)) + 0.01, n, 6)
  X <- matrix(rnorm(n * 4), n, 4)

  dat <- prepare_gbkmr_data(Y, Z, X,
    time_points = 3, mixture_components = 2,
    td_covariates = 1, baseline_covariates = 1,
    td_covariate_names = "waist", log_transform_mixtures = TRUE,
    validate_input = FALSE)

  res <- gbkmr_run(
    data = dat, time_points = 3, iter = 300, n = 40, K = 3,
    n_knots = 10, engine = "bkmr", verbose = FALSE)

  expect_true("diagnostics" %in% names(res))
  expect_true("outcome_Y" %in% names(res$diagnostics))
  expect_true(is.numeric(res$diagnostics$outcome_Y$ess))
  expect_true(is.numeric(res$diagnostics$outcome_Y$geweke_max_abs))
})

test_that("print and summary methods work", {
  skip_if_not_installed("bkmr")
  set.seed(42)
  n <- 60
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * 6)) + 0.01, n, 6)
  X <- matrix(rnorm(n * 4), n, 4)

  dat <- prepare_gbkmr_data(Y, Z, X,
    time_points = 3, mixture_components = 2,
    td_covariates = 1, baseline_covariates = 1,
    td_covariate_names = "waist", log_transform_mixtures = TRUE)

  res <- gbkmr_run(
    data = dat, outcome = "Y", outcome_type = "continuous",
    time_points = 3, iter = 200, n = 30, K = 3,
    n_knots = 10, engine = "bkmr", verbose = FALSE)

  expect_output(print(res), "g-BKMR Analysis Results")
  expect_output(summary(res), "Settings")
})
