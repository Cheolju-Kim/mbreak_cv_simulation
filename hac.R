#' HAC Correction Matrix
#'
#' Computes the long-run covariance matrix for score-like residual products,
#' optionally applying VAR(1)-style prewhitening and recoloring.
correct <- function(res, prewhit = FALSE) {
  d <- ncol(res)
  nt <- nrow(res)

  vmat <- res
  if (prewhit) {
    bmat <- matrix(0, nrow = d, ncol = d)
    vstar <- matrix(0, nrow = nt - 1, ncol = d)

    for (i in seq_len(d)) {
      x_lag <- vmat[1:(nt - 1), , drop = FALSE]
      y_now <- vmat[2:nt, i, drop = FALSE]
      b <- pinv(x_lag) %*% y_now
      bmat[, i] <- as.vector(b)
      vstar[, i] <- vmat[2:nt, i] - as.vector(vmat[1:(nt - 1), , drop = FALSE] %*% b)
    }

    jh <- jhatpr(vstar)
    recolor <- pinv(diag(d) - bmat)
    return(recolor %*% jh %*% t(recolor))
  }

  jhatpr(vmat)
}

#' Andrews HAC Long-Run Covariance
#'
#' Computes the HAC covariance estimate using the automatic bandwidth and
#' quadratic spectral kernel used by the original mbreak routines.
jhatpr <- function(vmat) {
  nt <- nrow(vmat)
  d <- ncol(vmat)

  st <- bandw(vmat)
  jhat <- t(vmat) %*% vmat

  for (j in seq_len(nt - 1)) {
    weight <- kern(j / st)
    jhat <- jhat + weight * t(vmat[(j + 1):nt, , drop = FALSE]) %*%
      vmat[1:(nt - j), , drop = FALSE]
  }

  for (j in seq_len(nt - 1)) {
    weight <- kern(j / st)
    jhat <- jhat + weight * t(vmat[1:(nt - j), , drop = FALSE]) %*%
      vmat[(j + 1):nt, , drop = FALSE]
  }

  jhat / (nt - d)
}

#' Automatic HAC Bandwidth
#'
#' Computes the plug-in bandwidth used by the HAC covariance estimator.
bandw <- function(vhat) {
  nt <- nrow(vhat)
  d <- ncol(vhat)

  a2n <- 0
  a2d <- 0
  for (i in seq_len(d)) {
    x_lag <- vhat[1:(nt - 1), i, drop = FALSE]
    y_now <- vhat[2:nt, i, drop = FALSE]
    b <- as.numeric(pinv(x_lag) %*% y_now)
    sig <- sum((vhat[2:nt, i] - b * vhat[1:(nt - 1), i])^2) / (nt - 1)
    a2n <- a2n + 4 * b * b * sig * sig / (1 - b)^8
    a2d <- a2d + sig * sig / (1 - b)^4
  }

  a2 <- a2n / a2d
  1.3221 * (a2 * nt)^0.2
}

#' Quadratic Spectral Kernel
#'
#' Evaluates the quadratic spectral kernel at a scalar argument.
kern <- function(x) {
  if (abs(x) < .Machine$double.eps) {
    return(1)
  }

  del <- 6 * pi * x / 5
  3 * (sin(del) / del - cos(del)) / (del * del)
}
