# Confidence intervals for break dates under shrinking-shift asymptotics.

#' Break-Date Confidence Intervals
#'
#' Computes confidence intervals for estimated break dates using the
#' shrinking-shift asymptotic approximation.
interval2 <- function(maty, matx, m, n, bigt, br, beta, vv,
                      hetq = FALSE,
                      hac = FALSE, brv = FALSE,
                      brbeta = TRUE, prewhit = FALSE) {
  q <- ncol(matx)
  bound <- matrix(0, nrow = m, ncol = 4)
  diagx <- pzbar(matx, m, br, bigt)
  br_all <- c(br, nrow(maty) / n)
  res <- maty - diagx %*% beta
  res <- matrix(as.vector(res), ncol = n, byrow = TRUE)

  for (j in seq_len(m + 1)) {
    if (j == 1) {
      idx <- 1:br_all[[j]]
    } else {
      idx <- (br_all[[j - 1]] + 1):br_all[[j]]
    }
    v_j <- vv[((j - 1) * n + 1):(j * n), , drop = FALSE]
    res[idx, ] <- res[idx, , drop = FALSE] %*% pinv(sqrm(v_j))
  }

  for (j in seq_len(m)) {
    vvar1 <- vv[((j - 1) * n + 1):(j * n), , drop = FALSE]
    vvar2 <- vv[(j * n + 1):((j + 1) * n), , drop = FALSE]
    A1 <- sqrm(vvar1) %*% pinv(vvar2) %*% (vvar2 - vvar1) %*% pinv(sqrm(vvar1))
    A2 <- sqrm(vvar2) %*% pinv(vvar1) %*% (vvar2 - vvar1) %*% pinv(sqrm(vvar2))
    start_left <- if (j == 1) 1 else br_all[[j - 1]] + 1
    end_left <- br_all[[j]]
    start_right <- br_all[[j]] + 1
    end_right <- br_all[[j + 1]]
    diffb <- beta[(j * q + 1):((j + 1) * q), , drop = FALSE] -
      beta[((j - 1) * q + 1):(j * q), , drop = FALSE]

    pieces <- intpc(
      matx = matx,
      res = res,
      n = n,
      q = q,
      vvar1 = vvar1,
      vvar2 = vvar2,
      A1 = A1,
      A2 = A2,
      start_left = start_left,
      end_left = end_left,
      start_right = start_right,
      end_right = end_right,
      bigt = bigt,
      hetq = hetq,
      hac = hac,
      brv = brv,
      brbeta = brbeta,
      prewhit = prewhit,
      diffb = diffb
    )

    Delt1 <- pieces$Delt1
    Delt2 <- pieces$Delt2
    Gam1 <- pieces$Gam1
    Gam2 <- pieces$Gam2
    eta <- Delt2 / Delt1
    phi1s <- Gam1 * Gam1 / Delt1
    phi2s <- Gam2 * Gam2 / Delt2
    cvf <- cvg(eta, phi1s, phi2s)
    a <- (Delt1 / Gam1)^2

    bound[j, 1] <- br_all[[j]] - cvf[4, 1] / a
    bound[j, 2] <- br_all[[j]] - cvf[1, 1] / a
    bound[j, 3] <- br_all[[j]] - cvf[3, 1] / a
    bound[j, 4] <- br_all[[j]] - cvf[2, 1] / a
  }

  bound[, 1] <- floor(bound[, 1])
  bound[, 2] <- floor(bound[, 2]) + 1
  bound[, 3] <- floor(bound[, 3])
  bound[, 4] <- floor(bound[, 4]) + 1
  bound
}

#' Critical Values for Break-Date Limits
#'
#' Computes critical values for the limiting distribution used in break-date
#' confidence intervals.
cvg <- function(eta, phi1s, phi2s) {
  if (!all(is.finite(c(eta, phi1s, phi2s))) || eta <= 0 || phi1s <= 0 || phi2s <= 0) {
    stop("`eta`, `phi1s`, and `phi2s` must be positive finite values.")
  }

  cvec <- matrix(0, nrow = 4, ncol = 1)
  a <- phi1s / phi2s
  gam <- ((phi2s / phi1s) + 1) * eta / 2
  b <- sqrt(phi1s * eta / phi2s)
  deld <- sqrt(phi2s * eta / phi1s) + b / 2
  alph <- a * (1 + a) / 2
  bet <- (1 + 2 * a) / 2
  sig <- c(0.025, 0.05, 0.95, 0.975)

  for (isig in seq_along(sig)) {
    upb <- max(c(200, 100 / gam))
    lwb <- -200
    crit <- 999999
    cct <- 1
    xx <- NA_real_

    while (abs(crit) >= 0.00001) {
      cct <- cct + 1
      if (cct > 1000) {
        warning("`cvg()` reached the iteration limit; CI may be inaccurate.")
        return(cvec)
      }

      xx <- lwb + (upb - lwb) / 20
      pval <- funcg(xx, bet, alph, b, deld, gam)
      crit <- pval - sig[[isig]]
      if (crit <= 0) {
        lwb <- xx
      } else {
        upb <- xx
      }
    }
    cvec[isig, 1] <- xx
  }

  cvec
}

#' Distribution Function for Break-Date Limits
#'
#' Evaluates the distribution function used by `cvg()` for break-date
#' confidence interval critical values.
funcg <- function(x, bet, alph, b, deld, gam) {
  if (x <= 0) {
    xb <- bet * sqrt(abs(x))
    if (abs(xb) <= 30) {
      return(
        -sqrt(-x / (2 * pi)) * exp(x / 8) -
          (bet / alph) * exp(-alph * x) * stats::pnorm(-bet * sqrt(abs(x))) +
          ((2 * bet * bet / alph) - 2 - x / 2) * stats::pnorm(-sqrt(abs(x)) / 2)
      )
    }

    aa <- log(bet / alph) - alph * x - xb^2 / 2 - log(sqrt(2 * pi)) - log(xb)
    return(
      -sqrt(-x / (2 * pi)) * exp(x / 8) -
        exp(aa) * stats::pnorm(-sqrt(abs(x)) / 2) +
        ((2 * bet * bet / alph) - 2 - x / 2) * stats::pnorm(-sqrt(abs(x)) / 2)
    )
  }

  xb <- deld * sqrt(x)
  if (abs(xb) <= 30) {
    return(
      1 + (b / sqrt(2 * pi)) * sqrt(x) * exp(-b * b * x / 8) +
        (b * deld / gam) * exp(gam * x) * stats::pnorm(-deld * sqrt(x)) +
        (2 - b * b * x / 2 - 2 * deld * deld / gam) * stats::pnorm(-b * sqrt(x) / 2)
    )
  }

  aa <- log(b * deld / gam) + gam * x - xb^2 / 2 - log(sqrt(2 * pi)) - log(xb)
  1 + (b / sqrt(2 * pi)) * sqrt(x) * exp(-b * b * x / 8) +
    exp(aa) +
    (2 - b * b * x / 2 - 2 * deld * deld / gam) * stats::pnorm(-b * sqrt(x) / 2)
}

#' Confidence-Interval Drift and Variance Pieces
#'
#' Computes the drift and variance terms that enter the limiting distribution
#' for a single break-date confidence interval.
intpc <- function(matx, res, n, q, vvar1, vvar2, A1, A2,
                                start_left, end_left, start_right, end_right,
                                bigt, hetq, hac, brv,
                                brbeta, prewhit, diffb) {
  Q1 <- matrix(0, nrow = q, ncol = q)
  Q2 <- matrix(0, nrow = q, ncol = q)
  Pi1 <- matrix(0, nrow = q, ncol = q)
  Pi2 <- matrix(0, nrow = q, ncol = q)

  if (!hetq && !brv && brbeta) {
    if (hac) {
      tempmat <- pitemp(matx, res, n, q, 1, bigt, vvar2, vvar1)
      pi1 <- correct(tempmat, prewhit = prewhit)
      Q1 <- qavg(matx, n, q, 1, bigt, vvar2)
      pi2 <- pi1
      Q2 <- Q1
      Delt1 <- as.numeric(t(diffb) %*% Q1 %*% diffb)
      Delt2 <- as.numeric(t(diffb) %*% Q2 %*% diffb)
      Gam1 <- sqrt(as.numeric(t(diffb) %*% pi1 %*% diffb))
      Gam2 <- sqrt(as.numeric(t(diffb) %*% pi2 %*% diffb))
      return(list(Delt1 = Delt1, Delt2 = Delt2, Gam1 = Gam1, Gam2 = Gam2))
    }

    for (k in seq_len(bigt)) {
      xk <- matx[stacked_rows(k, k, n), , drop = FALSE]
      Q1 <- Q1 + t(xk) %*% pinv(vvar2) %*% xk
      Pi1 <- Pi1 + t(xk) %*% pinv(vvar2) %*% vvar1 %*% pinv(vvar2) %*% xk
    }
    Q1 <- Q1 / bigt
    Pi1 <- Pi1 / bigt
    Q2 <- Q1
    Pi2 <- Pi1
    Delt1 <- as.numeric(t(diffb) %*% Q1 %*% diffb)
    Delt2 <- as.numeric(t(diffb) %*% Q2 %*% diffb)
    Gam1 <- sqrt(as.numeric(t(diffb) %*% Pi1 %*% diffb))
    Gam2 <- sqrt(as.numeric(t(diffb) %*% Pi2 %*% diffb))
    return(list(Delt1 = Delt1, Delt2 = Delt2, Gam1 = Gam1, Gam2 = Gam2))
  }

  if (hac) {
    pi1 <- correct(pitemp(matx, res, n, q, start_left, end_left, vvar2, vvar1),
                   prewhit = prewhit)
    omiga1 <- correct(omtemp(res, n, start_left, end_left),
                      prewhit = prewhit)
    Q1 <- qavg(matx, n, q, start_left, end_left, vvar2)

    pi2 <- correct(pitemp(matx, res, n, q, start_right, end_right, vvar1, vvar2),
                   prewhit = prewhit)
    omiga2 <- correct(omtemp(res, n, start_right, end_right),
                      prewhit = prewhit)
    Q2 <- qavg(matx, n, q, start_right, end_right, vvar1)

    Delt1 <- as.numeric(0.5 * sum(diag(A1 %*% A1)) + t(diffb) %*% Q1 %*% diffb)
    Delt2 <- as.numeric(0.5 * sum(diag(A2 %*% A2)) + t(diffb) %*% Q2 %*% diffb)
    Gam1 <- sqrt(as.numeric(0.25 * t(vec(A1)) %*% omiga1 %*% vec(A1) + t(diffb) %*% pi1 %*% diffb))
    Gam2 <- sqrt(as.numeric(0.25 * t(vec(A2)) %*% omiga2 %*% vec(A2) + t(diffb) %*% pi2 %*% diffb))
    return(list(Delt1 = Delt1, Delt2 = Delt2, Gam1 = Gam1, Gam2 = Gam2))
  }

  left <- intside(matx, res, n, q, start_left, end_left, vvar1, vvar2)
  right <- intside(matx, res, n, q, start_right, end_right, vvar2, vvar1)

  Delt1 <- as.numeric(0.5 * sum(diag(A1 %*% A1)) + t(diffb) %*% left$Q %*% diffb)
  Delt2 <- as.numeric(0.5 * sum(diag(A2 %*% A2)) + t(diffb) %*% right$Q %*% diffb)
  Gam1 <- sqrt(as.numeric(0.25 * t(vec(A1)) %*% left$omega %*% vec(A1) + t(diffb) %*% left$Pi %*% diffb))
  Gam2 <- sqrt(as.numeric(0.25 * t(vec(A2)) %*% right$omega %*% vec(A2) + t(diffb) %*% right$Pi %*% diffb))
  list(Delt1 = Delt1, Delt2 = Delt2, Gam1 = Gam1, Gam2 = Gam2)
}

#' Average Regressor Moment
#'
#' Computes an average transformed regressor moment over one segment.
qavg <- function(matx, n, q, start, end, v_other) {
  Q <- matrix(0, nrow = q, ncol = q)
  for (i in start:end) {
    x_i <- matx[stacked_rows(i, i, n), , drop = FALSE]
    Q <- Q + t(x_i) %*% pinv(v_other) %*% x_i
  }
  Q / (end - start + 1)
}

#' Non-HAC Confidence-Interval Side Pieces
#'
#' Computes covariance, regressor moment, and score variance pieces for one
#' side of a break under the non-HAC approximation.
intside <- function(matx, res, n, q, start, end, v_own, v_other) {
  omega <- matrix(0, nrow = n * n, ncol = n * n)
  Q <- matrix(0, nrow = q, ncol = q)
  Pi <- matrix(0, nrow = q, ncol = q)
  for (i in start:end) {
    ri <- matrix(res[i, ], nrow = 1)
    x_i <- matx[stacked_rows(i, i, n), , drop = FALSE]
    omega <- omega + vec(t(ri) %*% ri - diag(n)) %*% t(vec(t(ri) %*% ri - diag(n)))
    Q <- Q + t(x_i) %*% pinv(v_other) %*% x_i
    Pi <- Pi + t(x_i) %*% pinv(v_other) %*% v_own %*% pinv(v_other) %*% x_i
  }
  scale <- end - start + 1
  list(omega = omega / scale, Q = Q / scale, Pi = Pi / scale)
}

#' HAC Score Products for Coefficient Shifts
#'
#' Builds the observation-by-observation score-product matrix used in HAC
#' calculations for coefficient-break terms.
pitemp <- function(matx, res, n, q, start, end, v_other, v_own) {
  tempmat <- matrix(0, nrow = end - start + 1, ncol = q)
  for (k in seq_len(nrow(tempmat))) {
    i <- start + k - 1
    x_i <- matx[stacked_rows(i, i, n), , drop = FALSE]
    tempmat[k, ] <- as.vector(t(x_i) %*% pinv(v_other) %*% sqrm(v_own) %*% matrix(res[i, ], ncol = 1))
  }
  tempmat
}

#' HAC Score Products for Covariance Shifts
#'
#' Builds the observation-by-observation covariance score matrix used in HAC
#' calculations for covariance-break terms.
omtemp <- function(res, n, start, end) {
  tempmat <- matrix(0, nrow = end - start + 1, ncol = n * n)
  for (k in seq_len(nrow(tempmat))) {
    i <- start + k - 1
    ri <- matrix(res[i, ], nrow = 1)
    tempmat[k, ] <- as.vector(t(vec(t(ri) %*% ri - diag(n))))
  }
  tempmat
}
