#' Coefficient Standard Errors Conditional on Break Dates
#'
#' Computes variance-covariance matrices and standard errors for the stacked
#' regime coefficient vector, treating estimated break dates as fixed.
coefse <- function(maty, matz, n, m, bigt, br, beta, vv,
                                        R = NULL, hac = FALSE,
                                        brv = FALSE,
                                        prewhit = FALSE) {
  q <- ncol(matz)
  xbar <- pzbar(matz, m, br, bigt)
  beta <- matrix(beta, ncol = 1)

  R <- valR(R, m = m, q = q)
  xtheta <- xbar %*% R

  residual <- maty - xbar %*% beta
  weighted <- wsys(
    x = xtheta,
    residual = residual,
    vv = vv,
    n = n,
    m = m,
    bigt = bigt,
    br = br,
    brv = brv
  )

  vcov_theta <- if (hac) {
    vcovhac(
      xstar = weighted$x,
      ustar = weighted$residuals,
      n = n,
      m = m,
      bigt = bigt,
      br = br,
      brv = brv,
      prewhit = prewhit
    )
  } else {
    pinv(crossprod(weighted$x))
  }

  vcov_beta <- R %*% vcov_theta %*% t(R)
  vcov_beta <- (vcov_beta + t(vcov_beta)) / 2
  vcov_theta <- (vcov_theta + t(vcov_theta)) / 2

  beta_names <- bnames(m = m, q = q)
  dimnames(vcov_beta) <- list(beta_names, beta_names)
  se_beta <- sqdiag(vcov_beta)
  names(se_beta) <- beta_names

  list(
    vcov = vcov_beta,
    se = se_beta,
    se_by_regime = byregime(se_beta, m, q),
    vcov_theta = vcov_theta,
    se_theta = sqdiag(vcov_theta)
  )
}

#' Whiten Coefficient Design and Residuals
#'
#' Applies the covariance weighting used by the GLS objective.
wsys <- function(x, residual, vv, n, m, bigt, br,
                                        brv) {
  xstar <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  ustar <- matrix(0, nrow = bigt, ncol = n)
  endpoints <- c(0, br, bigt)

  for (seg in seq_len(m + 1)) {
    sigma_idx <- if (brv) seg else 1
    sigma <- vvblock(vv, sigma_idx, n)
    sigma_inv_sqrt <- pinv(sqrm(sigma))

    t_start <- endpoints[[seg]] + 1
    t_end <- endpoints[[seg + 1]]
    for (tt in t_start:t_end) {
      rows <- stacked_rows(tt, tt, n)
      xstar[rows, ] <- sigma_inv_sqrt %*% x[rows, , drop = FALSE]
      ustar[tt, ] <- as.vector(sigma_inv_sqrt %*% residual[rows, , drop = FALSE])
    }
  }

  list(x = xstar, residuals = ustar)
}

#' HAC Variance for a Weighted Coefficient Design
vcovhac <- function(xstar, ustar, n, m, bigt, br,
                                 brv, prewhit) {
  scores <- scores(xstar, ustar, n = n, bigt = bigt)
  xx <- crossprod(xstar) / bigt
  bread <- pinv(xx)

  if (brv) {
    endpoints <- c(0, br, bigt)
    meat <- matrix(0, nrow = ncol(scores), ncol = ncol(scores))
    for (seg in seq_len(m + 1)) {
      t_start <- endpoints[[seg]] + 1
      t_end <- endpoints[[seg + 1]]
      seg_sc <- scores[t_start:t_end, , drop = FALSE]
      weight <- (t_end - t_start + 1) / bigt
      meat <- meat + weight * correct(seg_sc, prewhit = prewhit)
    }
  } else {
    meat <- correct(scores, prewhit = prewhit)
  }

  bread %*% (meat / bigt) %*% bread
}

#' Coefficient Score Contributions
scores <- function(xstar, ustar, n, bigt) {
  scores <- matrix(0, nrow = bigt, ncol = ncol(xstar))
  for (tt in seq_len(bigt)) {
    rows <- stacked_rows(tt, tt, n)
    scores[tt, ] <- as.vector(
      t(xstar[rows, , drop = FALSE]) %*% matrix(ustar[tt, ], ncol = 1)
    )
  }
  scores
}

#' Extract a Regime Covariance Block
vvblock <- function(vv, seg, n) {
  rows <- ((seg - 1) * n + 1):(seg * n)
  vv[rows, , drop = FALSE]
}

#' Stable Square Root of a Variance Diagonal
sqdiag <- function(x) {
  variances <- diag(x)
  variances[variances < 0 & variances > -sqrt(.Machine$double.eps)] <- 0
  sqrt(variances)
}

#' Default Names for Stacked Regime Coefficients
bnames <- function(m, q) {
  paste0(
    "regime", rep(seq_len(m + 1), each = q),
    ".coef", rep(seq_len(q), times = m + 1)
  )
}
