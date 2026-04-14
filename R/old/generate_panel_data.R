generate_panel_data <- function(
    popN = 1e6,
    T = 5,
    Adim = 3,
    Ldim = 1,
    outcome_type = "binary",
    relationship_type = "quadratic",
    confounding = "high"
) {

  k <- popN
  p <- Adim

  covmat <- switch(confounding,
                   "low" = diag(p) * 1 + 0.05 * (1 - diag(p)),  # nearly independent

                   "high" = if (p == 3) {
                     matrix(c(
                       1,   0.05, 0.2,
                       0.05, 1,   0.1,
                       0.2,  0.1, 1
                     ), nrow = 3, byrow = TRUE)
                   } else {
                     # default high correlation for p > 3
                     toeplitz(0.5 ^ (0:(p - 1)))
                   }
  )


  sigsq.trueL <- 1
  sigsq.trueY <- 1
  sigsq.trueA <- rep(0.05, T - 1)

  set.seed(1)
  dat <- list()
  dat$sex <- rbinom(k, 1, 0.5)
  dat$C0 <- matrix(rnorm(k * Ldim, mean = rep(0.25 * (1 - dat$sex), each = Ldim), sd = 1), ncol = Ldim)

  dat$At0 <- mvtnorm::rmvnorm(k, rep(0, p), sigma = covmat)
  colnames(dat$At0) <- paste0("At0_", 1:p)

  make_hfunL <- function(type) {
    function(z) {
      out <- 0
      if (type == "linear") {
        out <- rowSums(z[, 1:min(2, ncol(z)), drop = FALSE] * c(0.25, 0.5)[1:min(2, ncol(z))])
      } else if (type == "quadratic") {
        out <- rowSums(z[, 1:min(2, ncol(z)), drop = FALSE] * c(0.25, 0.5)[1:min(2, ncol(z))])
        if (ncol(z) >= 2) out <- out + 0.25 * z[, 2]^2
      } else if (type == "quadratic+interaction") {
        out <- rowSums(z[, 1:min(2, ncol(z)), drop = FALSE] * c(0.25, 0.5)[1:min(2, ncol(z))])
        if (ncol(z) >= 2) out <- out + 0.25 * z[, 2]^2 + 0.25 * z[, 1] * z[, 2]
      }
      return(out)
    }
  }

  for (t in 1:(T - 1)) {
    prev_At <- dat[[paste0("At", t - 1)]]
    prev_L <- if (t == 1) dat$C0 else dat[[paste0("L", t - 1)]]
    hAt_t <- matrix(NA, k, p)
    for (j in 1:p) {
      hAt_t[, j] <- 1/3 * prev_At[, j]^2 + 0.1 * prev_At[, j] +
        rowSums(0.1 * prev_L)
    }
    dat[[paste0("hAt", t)]] <- hAt_t
    dat[[paste0("epsAt", t)]] <- matrix(rnorm(k * p, sd = sqrt(sigsq.trueA[t])), k, p)
    dat[[paste0("At", t)]] <- hAt_t + dat[[paste0("epsAt", t)]]
    colnames(dat[[paste0("At", t)]]) <- paste0("At", t, "_", 1:p)

    hfunL <- make_hfunL(relationship_type)
    At_curr <- dat[[paste0("At", t)]]
    if (ncol(At_curr) >= 2) {
      A_t2 <- At_curr[, 2]
    } else {
      A_t2 <- At_curr[, 1]  # fall back to the only available exposure
    }
    L_input <- cbind(A_t2, prev_L[, 1])

    dat[[paste0("hL", t)]] <- hfunL(L_input)

    # All L are noisy around the same mean
    epsL_t <- matrix(rnorm(k * Ldim, sd = sqrt(sigsq.trueL)), ncol = Ldim)
    hL_t_matrix <- matrix(dat[[paste0("hL", t)]], ncol = 1)
    dat[[paste0("L", t)]] <- matrix(rep(hL_t_matrix, Ldim), ncol = Ldim) + epsL_t

  }

  At_list <- lapply(0:(T - 1), function(t) dat[[paste0("At", t)]])
  names(At_list) <- paste0("At", 0:(T - 1))
  L_list <- c(list(L0 = dat$C0), dat[paste0("L", 1:(T - 1))])
  names(L_list) <- paste0("L", 0:(T - 1))

  dat$ALL <- do.call(cbind, c(At_list, L_list, list(sex = matrix(dat$sex, ncol = 1))))

  colnames(dat$ALL) <- c(
    unlist(lapply(0:(T - 1), function(t) paste0("At", t, "_", 1:p))),
    unlist(lapply(0:(T - 1), function(t) paste0("L", t, "_", 1:Ldim))),
    "sex"
  )

  hfunY <- function(z) {
    out <- 0
    for (t in 0:(T - 1)) {
      L <- z[[paste0("L", t, "_1")]]
      M2 <- if (paste0("At", t, "_2") %in% names(z)) {
        z[[paste0("At", t, "_2")]]
      } else {
        z[[paste0("At", t, "_1")]]  # fallback if only 1 exposure
      }

      term <- 0.25 * M2 + 0.25 * L + 0.25 * L^2
      if (relationship_type == "quadratic+interaction") {
        term <- term + 0.1 * M2 * L
      }
      out <- out + term
    }
    out
  }


  dat$hY <- apply(dat$ALL, 1, function(row) {
    row_list <- as.list(row)
    names(row_list) <- colnames(dat$ALL)
    hfunY(row_list)
  })
  dat$epsY <- rnorm(k, sd = sqrt(sigsq.trueY))
  dat$z <- dat$hY + dat$epsY
  dat$y <- if (outcome_type == "binary") rbinom(k, 1, pnorm(dat$z)) else dat$z

  for (t in 0:(T - 1)) {
    colnames(dat[[paste0("At", t)]]) <- paste0("logM", 1:p, "_", t)
  }

  df <- data.frame(
    sex = dat$sex,
    dat$C0,
    do.call(cbind, dat[paste0("At", 0:(T - 1))]),
    do.call(cbind, L_list[-1]),
    Y = dat$y
  )
  df$id <- 1:k

  waist0_names <- if (Ldim == 1) "waist0" else paste0("waist0_", 1:Ldim)
  waist_t_names <- if (Ldim == 1) {
    paste0("waist", 1:(T - 1))
  } else {
    unlist(lapply(1:(T - 1), function(t) paste0("waist", t, "_", 1:Ldim)))
  }
  logM_names <- if (p == 1) {
    paste0("logM", 0:(T - 1))
  } else {
    unlist(lapply(0:(T - 1), function(t) paste0("logM", 1:p, "_", t)))
  }

  colnames(df) <- c("sex", waist0_names, logM_names, waist_t_names, "Y", "id")


  corrplot(cor(df[setdiff(names(df), c("id"))][sapply(df[setdiff(names(df), c("id"))], is.numeric)]), type = "upper")
  saveRDS(df, file = paste0("popn_", T, "t_", relationship_type, "_", confounding, "_", outcome_type,
                            "_Adim", Adim, "_Ldim", Ldim, ".rds"))

  return(df)
}











