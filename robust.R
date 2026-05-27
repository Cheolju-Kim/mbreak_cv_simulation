# Robust SupF-style tests for serially correlated errors.

#' HAC-Robust SupF Statistic
#'
#' Computes the robust SupF/Wald statistic for coefficient breaks at fixed break
#' dates in the unrestricted model.
supft <- function(maty, matx, m, n, br, rbeta, bigt,
                  prewhit = FALSE) {
  pmatx <- pzbar(matx, m, br, bigt)
  res <- maty - pmatx %*% rbeta
  umat <- matrix(as.vector(res), ncol = n, byrow = TRUE)
  sigma <- crossprod(umat) / bigt
  bigsig <- pinv(sqrm(kronecker(diag(bigt), sigma)))
  xbstar <- bigsig %*% pmatx
  ubstar <- matrix(as.vector(bigsig %*% res), ncol = n, byrow = TRUE)

  tempmat <- matrix(0, nrow = bigt, ncol = ncol(xbstar))
  for (k in seq_len(bigt)) {
    idx <- stacked_rows(k, k, n)
    tempmat[k, ] <- as.vector(t(xbstar[idx, , drop = FALSE]) %*% matrix(ubstar[k, ], ncol = 1))
  }

  vmat <- correct(tempmat, prewhit = prewhit)
  xx <- crossprod(xbstar) / bigt
  vbeta <- pinv(xx) %*% (vmat / bigt) %*% pinv(xx)
  rmat <- rdiff(m, ncol(matx))

  ftest <- t(rbeta) %*% t(rmat) %*% pinv(rmat %*% vbeta %*% t(rmat)) %*% rmat %*% rbeta
  as.numeric((bigt - (m + 1) * ncol(matx)) * ftest / bigt)
}

#' Restricted HAC-Robust SupF Statistic
#'
#' Computes the robust SupF/Wald statistic for coefficient breaks when the
#' coefficient vector is reparameterized by a restriction matrix.
suprft <- function(maty, matx, m, n, br, rbeta, R, bigt,
                   prewhit = FALSE) {
  pmatx <- pzbar(matx, m, br, bigt)
  res <- maty - pmatx %*% rbeta
  umat <- matrix(as.vector(res), ncol = n, byrow = TRUE)
  sigma <- crossprod(umat) / bigt
  bigsig <- pinv(sqrm(kronecker(diag(bigt), sigma)))
  pmatx_r <- pmatx %*% R
  xbstar <- bigsig %*% pmatx_r
  ubstar <- matrix(as.vector(bigsig %*% res), ncol = n, byrow = TRUE)

  tempmat <- matrix(0, nrow = bigt, ncol = ncol(xbstar))
  for (k in seq_len(bigt)) {
    idx <- stacked_rows(k, k, n)
    tempmat[k, ] <- as.vector(t(xbstar[idx, , drop = FALSE]) %*% matrix(ubstar[k, ], ncol = 1))
  }

  vmat <- correct(tempmat, prewhit = prewhit)
  xx <- crossprod(xbstar) / bigt
  vdelta <- pinv(xx) %*% (vmat / bigt) %*% pinv(xx)
  vbeta <- R %*% vdelta %*% t(R)
  rmat <- rdiff(m, ncol(matx))

  ftest <- t(rbeta) %*% t(rmat) %*% pinv(rmat %*% vbeta %*% t(rmat)) %*% rmat %*% rbeta
  as.numeric((bigt - (m + 1) * ncol(matx)) * ftest / bigt)
}

#' Sequential Robust F Statistic
#'
#' Computes the sequential `l + 1 | l` robust F statistic by inserting one
#' additional coefficient break into each admissible segment.
seqf <- function(maty, matz, n, br, bigt, h,
                 prewhit = FALSE, trm = NULL) {
  if (is.null(trm)) {
    trm <- h / bigt
  }

  nseg <- length(br) + 1
  dv <- c(0, br, bigt)
  lrtest <- numeric(nseg)
  no_room <- 0
  q <- ncol(matz)

  for (is in seq_len(nseg)) {
    seg_length <- dv[[is + 1]] - dv[[is]]
    if (seg_length >= 2 * h) {
      starti <- dv[[is]] + 1
      endi <- dv[[is + 1]]
      idx <- stacked_rows(starti, endi, n)
      segy <- maty[idx, , drop = FALSE]
      segz <- matz[idx, , drop = FALSE]
      seg_bigt <- endi - starti + 1
      fit <- est(
        segy, segz, n, 1, seg_bigt, trm, diag(q * 2),
        brbeta = TRUE,
        brv = FALSE
      )
      lrtest[is] <- supft(
        segy, segz, 1, n, fit$dx, fit$rbeta, seg_bigt,
        prewhit = prewhit
      )
    } else {
      no_room <- no_room + 1
      lrtest[is] <- 0
    }
  }

  if (no_room == nseg) {
    warning("No segment has enough observations to insert an additional break.")
  }

  max(lrtest)
}

#' Adjacent-Regime Difference Matrix
#'
#' Builds the matrix that takes adjacent differences of regime-specific
#' coefficient vectors.
rdiff <- function(m, q) {
  rsub <- matrix(0, nrow = m, ncol = m + 1)
  for (j in seq_len(m)) {
    rsub[j, j] <- -1
    rsub[j, j + 1] <- 1
  }
  kronecker(rsub, diag(q))
}
