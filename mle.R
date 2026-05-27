# Unrestricted segment MLE and dynamic programming routines.

#' Segment MLE Objective Values
#'
#' Computes Gaussian negative log-likelihood values for all segments beginning
#' at a fixed start date and ending at admissible dates.
mlef <- function(start, maty, matz, h, n, last, tol = 1e-6, max_iter = 1000) {
  vecmle <- matrix(0, nrow = last, ncol = 1)
  i <- start
  j <- start + h - 1

  idx <- stacked_rows(i, j, n)
  segx <- matz[idx, , drop = FALSE]
  segy <- maty[idx, , drop = FALSE]
  fit <- fgls(segx, segy, n = n, tol = tol, max_iter = max_iter)

  bstar <- fit$beta
  vvar <- fit$sigma
  ibigv <- kronecker(diag(h), pinv(vvar))
  ihstar <- pinv(t(segx) %*% ibigv %*% segx)
  ystar <- segy
  zstar <- segx
  warm_sigma <- vvar

  vecmle[j, 1] <- (h * n / 2) * (log(2 * pi) + 1) +
    (h / 2) * logdet(vvar)
  j <- j + 1

  while (j <= last) {
    idx_j <- stacked_rows(j, j, n)
    x_j <- matz[idx_j, , drop = FALSE]
    y_j <- maty[idx_j, , drop = FALSE]

    gstar <- ihstar %*% t(x_j)
    tempz <- rbind(zstar, x_j)
    tempy <- rbind(ystar, y_j)
    b <- NULL
    icstar <- NULL

    upd <- function(prev_sigma) {
      icstar <<- pinv(prev_sigma + x_j %*% ihstar %*% t(x_j))
      b <<- bstar + gstar %*% icstar %*% (y_j - x_j %*% bstar)
      res <- tempy - tempz %*% b
      umat <- matrix(as.vector(res), ncol = n, byrow = TRUE)
      t(umat) %*% umat / (j - i + 1)
    }
    fp <- fixed_point(warm_sigma, upd, tol, max_iter,
                      "mlef() recursive FGLS update")
    vvar <- fp$value

    vecmle[j, 1] <- ((j - i + 1) * n / 2) * (log(2 * pi) + 1) +
      ((j - i + 1) / 2) * logdet(vvar)

    bstar <- b
    zstar <- tempz
    ystar <- tempy
    warm_sigma <- vvar
    ihstar <- ihstar - gstar %*% icstar %*% t(gstar)
    j <- j + 1
  }

  vecmle
}

#' Two-Segment Partition Search
#'
#' Finds the best single split between `b1` and `b2` for a segment from
#' `start` to `last`, using precomputed segment objective values.
parti <- function(start, b1, b2, last, bigvec, bigt) {
  dvec <- matrix(0, nrow = bigt, ncol = 1)
  for (j in b1:b2) {
    dvec[j, 1] <- bigvec[bigvec_index(start, j, bigt), 1] +
      bigvec[bigvec_index(j + 1, last, bigt), 1]
  }

  candidates <- dvec[b1:b2, 1]
  list(
    ssrmin = min(candidates),
    dx = as.integer((b1 - 1) + which.min(candidates))
  )
}

#' Dynamic Programming Break-Date Search
#'
#' Computes optimal unrestricted break dates and objective values for
#' `1, ..., m` breaks using dynamic programming.
dating_MLE <- function(maty, matz, n, h, m, bigt, tol = 1e-6, max_iter = 1000) {
  datevec <- matrix(0, nrow = m, ncol = m)
  optdat <- matrix(0, nrow = bigt, ncol = m)
  optmle <- matrix(0, nrow = bigt, ncol = m)
  dvec <- matrix(0, nrow = bigt, ncol = 1)
  global <- matrix(0, nrow = m, ncol = 1)
  bigvec <- matrix(0, nrow = bigt * (bigt + 1) / 2, ncol = 1)

  for (i in seq_len(bigt - h + 1)) {
    vecmle <- mlef(i, maty, matz, h, n, bigt, tol = tol, max_iter = max_iter)
    idx <- bigvec_index(i, i, bigt):bigvec_index(i, bigt, bigt)
    bigvec[idx, 1] <- vecmle[i:bigt, 1]
  }

  if (m == 1) {
    one <- parti(1, h, bigt - h, bigt, bigvec, bigt)
    datevec[1, 1] <- one$dx
    global[1, 1] <- one$ssrmin
  } else {
    for (j1 in (2 * h):bigt) {
      one <- parti(1, h, j1 - h, j1, bigvec, bigt)
      optmle[j1, 1] <- one$ssrmin
      optdat[j1, 1] <- one$dx
    }
    global[1, 1] <- optmle[bigt, 1]
    datevec[1, 1] <- optdat[bigt, 1]

    for (ib in 2:m) {
      if (ib == m) {
        jlast <- bigt
        for (jb in (ib * h):(jlast - h)) {
          dvec[jb, 1] <- optmle[jb, ib - 1] +
            bigvec[bigvec_index(jb + 1, jlast, bigt), 1]
        }
        candidates <- dvec[(ib * h):(jlast - h), 1]
        optmle[jlast, ib] <- min(candidates)
        optdat[jlast, ib] <- as.integer((ib * h - 1) + which.min(candidates))
      } else {
        for (jlast in ((ib + 1) * h):bigt) {
          for (jb in (ib * h):(jlast - h)) {
            dvec[jb, 1] <- optmle[jb, ib - 1] +
              bigvec[bigvec_index(jb + 1, jlast, bigt), 1]
          }
          candidates <- dvec[(ib * h):(jlast - h), 1]
          optmle[jlast, ib] <- min(candidates)
          optdat[jlast, ib] <- as.integer((ib * h - 1) + which.min(candidates))
        }
      }

      datevec[ib, ib] <- optdat[bigt, ib]
      for (i in seq_len(ib - 1)) {
        xx <- ib - i
        datevec[xx, ib] <- optdat[datevec[xx + 1, ib], xx]
      }
      global[ib, 1] <- optmle[bigt, ib]
    }
  }

  list(global = -global, datevec = datevec, bigvec = bigvec)
}

#' Unrestricted Coefficients and Covariances
#'
#' Estimates regime-specific coefficients and covariance matrices for fixed
#' unrestricted break dates.
estim <- function(maty, matz, n, m, br, tol = 1e-6, max_iter = 200) {
  bigt <- nrow(maty) / n
  br <- if (m == 0) integer() else br
  endpoints <- c(0, br, bigt)
  beta <- matrix(0, nrow = m + 1, ncol = ncol(matz))
  vv <- matrix(0, nrow = (m + 1) * n, ncol = n)

  for (k in seq_len(m + 1)) {
    i <- endpoints[[k]] + 1
    j <- endpoints[[k + 1]]
    idx <- stacked_rows(i, j, n)
    segx <- matz[idx, , drop = FALSE]
    segy <- maty[idx, , drop = FALSE]
    fit <- fgls(segx, segy, n = n, tol = tol, max_iter = max_iter)
    beta[k, ] <- as.vector(fit$beta)
    vv[((k - 1) * n + 1):(k * n), ] <- fit$sigma
  }

  list(beta = beta, vv = vv)
}

#' Feasible GLS Fit
#'
#' Fits a stacked multivariate regression segment by iterating between GLS
#' coefficients and the residual covariance matrix.
fgls <- function(segx, segy, n, tol = 1e-6, max_iter = 200) {
  seg_len <- nrow(segy) / n
  beta <- pinv(t(segx) %*% segx) %*% t(segx) %*% segy
  residual <- segy - segx %*% beta
  umat <- matrix(as.vector(residual), ncol = n, byrow = TRUE)
  sigma <- t(umat) %*% umat / seg_len

  upd <- function(prev_sigma) {
    weight <- kronecker(diag(seg_len), pinv(prev_sigma))
    beta <<- pinv(t(segx) %*% weight %*% segx) %*% t(segx) %*% weight %*% segy
    residual <- segy - segx %*% beta
    umat <- matrix(as.vector(residual), ncol = n, byrow = TRUE)
    t(umat) %*% umat / seg_len
  }
  fp <- fixed_point(sigma, upd, tol, max_iter, "fgls()")

  list(beta = beta, sigma = fp$value,
       iterations = fp$iterations, converged = fp$converged)
}
