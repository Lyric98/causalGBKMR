test_that("detect_variable_patterns finds user-named covariates", {
  set.seed(1)
  n <- 50
  dat <- data.frame(
    id = 1:n, sex = rbinom(n, 1, 0.5),
    waist_0 = rnorm(n), logM1_0 = rnorm(n), logM2_0 = rnorm(n),
    logM1_1 = rnorm(n), logM2_1 = rnorm(n), waist_1 = rnorm(n),
    Y = rnorm(n)
  )

  det <- detect_variable_patterns(dat, T = 2)
  expect_equal(det$p, 2)
  expect_equal(det$Ldim, 1)
  expect_equal(det$td_covariate_names, "waist")
  expect_s3_class(det, "gbkmr_detection")
})

test_that("detect_variable_patterns finds generic covariates", {
  set.seed(1)
  n <- 30
  dat <- data.frame(
    id = 1:n, sex = rbinom(n, 1, 0.5),
    td_covariate1_0 = rnorm(n),
    logM1_0 = rnorm(n), logM2_0 = rnorm(n),
    logM1_1 = rnorm(n), logM2_1 = rnorm(n),
    td_covariate1_1 = rnorm(n),
    Y = rnorm(n)
  )

  det <- detect_variable_patterns(dat, T = 2)
  expect_equal(det$p, 2)
  expect_equal(det$Ldim, 1)
  expect_equal(det$detected_pattern, "generic_format")
})

test_that("detect_variable_patterns counts exposures correctly", {
  set.seed(1)
  n <- 20
  dat <- data.frame(
    id = 1:n, sex = rbinom(n, 1, 0.5),
    bmi_0 = rnorm(n), bp_0 = rnorm(n),
    logM1_0 = rnorm(n), logM2_0 = rnorm(n), logM3_0 = rnorm(n),
    logM1_1 = rnorm(n), logM2_1 = rnorm(n), logM3_1 = rnorm(n),
    bmi_1 = rnorm(n), bp_1 = rnorm(n),
    Y = rnorm(n)
  )

  det <- detect_variable_patterns(dat, T = 2)
  expect_equal(det$p, 3)
  expect_equal(det$Ldim, 2)
})
