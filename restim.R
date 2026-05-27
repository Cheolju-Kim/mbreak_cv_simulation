# Restricted and conditional estimation routines.

#' Restricted Coefficient Estimation
#'
#' Estimates coefficients and covariance matrices for fixed break dates under a
#' restriction or reparameterization matrix.
r_estim <- function(maty, matz, bigt, n, m, br, R = NULL,
                    brbeta = TRUE, brv = FALSE,
                    tol = 1e-6, max_iter = 1000) {
  q <- ncol(matz)
  br <- if (m == 0) integer() else br
  endpoints <- c(0, br, bigt)

  R <- valR(R, m = m, q = q)

  unrestricted <- estim(maty, matz, n = n, m = m, br = br)
  beta <- unrestricted$beta
  vv <- unrestricted$vv

  if (brv && brbeta) {
    xbar <- pzbar(matz, m, br, bigt)
    pmatz <- xbar %*% R
    ibigv <- bwt(endpoints, n = n, bigt = bigt, vv = vv)
    b <- pinv(t(pmatz) %*% ibigv %*% pmatz) %*% t(pmatz) %*% ibigv %*% maty

    upd <- function(prev) {
      vv    <<- segvv(pmatz, maty, prev, endpoints, n)
      ibigv <<- bwt(endpoints, n = n, bigt = bigt, vv = vv)
      pinv(t(pmatz) %*% ibigv %*% pmatz) %*% t(pmatz) %*% ibigv %*% maty
    }
    res <- fixed_point(b, upd, tol, max_iter,
                       "r_estim() coefficient-and-covariance FGLS")
    b <- res$value

    nbeta <- R %*% b
    nvv <- vv
    return(list(nbeta = nbeta, nvv = nvv, theta = b,
                iterations = res$iterations, converged = res$converged))
  }

  if (brv && !brbeta) {
    iter_limit <- min(max_iter, 200)
    ibigv <- bwt(endpoints, n = n, bigt = bigt, vv = vv)
    b <- pinv(t(matz) %*% ibigv %*% matz) %*% t(matz) %*% ibigv %*% maty

    upd <- function(prev) {
      vv    <<- segvv(matz, maty, prev, endpoints, n)
      ibigv <<- bwt(endpoints, n = n, bigt = bigt, vv = vv)
      pinv(t(matz) %*% ibigv %*% matz) %*% t(matz) %*% ibigv %*% maty
    }
    res <- fixed_point(b, upd, tol, iter_limit,
                       "r_estim() covariance-only FGLS")
    b <- res$value

    nbeta <- do.call(rbind, replicate(m + 1, b, simplify = FALSE))
    nvv <- vv
    return(list(nbeta = nbeta, nvv = nvv, theta = b,
                iterations = res$iterations, converged = res$converged))
  }

  if (brbeta && !brv) {
    vvar <- matrix(0, nrow = n, ncol = n)
    for (k in seq_len(m + 1)) {
      i <- endpoints[[k]] + 1
      j <- endpoints[[k + 1]]
      idx <- stacked_rows(i, j, n)
      tempx <- matz[idx, , drop = FALSE]
      tempy <- maty[idx, , drop = FALSE]
      resk <- tempy - tempx %*% matrix(beta[k, ], ncol = 1)
      umat <- matrix(as.vector(resk), ncol = n, byrow = TRUE)
      vvar <- vvar + t(umat) %*% umat
    }

    vvar <- vvar / bigt
    ibigv <- kronecker(diag(bigt), pinv(vvar))
    xbar <- pzbar(matz, m, br, bigt)
    pmatz <- xbar %*% R
    b <- pinv(t(pmatz) %*% ibigv %*% pmatz) %*% t(pmatz) %*% ibigv %*% maty

    upd <- function(prev) {
      resk <- maty - pmatz %*% prev
      umat <- matrix(as.vector(resk), ncol = n, byrow = TRUE)
      vvar  <<- t(umat) %*% umat / bigt
      ibigv <<- kronecker(diag(bigt), pinv(vvar))
      pinv(t(pmatz) %*% ibigv %*% pmatz) %*% t(pmatz) %*% ibigv %*% maty
    }
    res <- fixed_point(b, upd, tol, max_iter,
                       "r_estim() coefficient-break FGLS")
    b <- res$value

    nbeta <- R %*% b
    nvv <- do.call(rbind, replicate(m + 1, vvar, simplify = FALSE))
    return(list(nbeta = nbeta, nvv = nvv, theta = b,
                iterations = res$iterations, converged = res$converged))
  }

  iter_limit <- min(max_iter, 200)
  fit <- fgls(matz, maty, n = n, tol = tol, max_iter = iter_limit)
  nbeta <- do.call(rbind, replicate(m + 1, fit$beta, simplify = FALSE))
  nvv <- do.call(rbind, replicate(m + 1, fit$sigma, simplify = FALSE))
  list(nbeta = nbeta, nvv = nvv, theta = fit$beta, iterations = fit$iterations, converged = fit$converged)
}

#' Validate a beta = R theta Reparameterization Matrix
valR <- function(R, m, q, arg = "R") {
  exp_rows <- (m + 1) * q
  if (is.null(R)) {
    return(diag(exp_rows))
  }

  R <- as.matrix(R)
  if (!is.numeric(R)) {
    stop(sprintf("`%s` must be a numeric matrix.", arg))
  }
  if (nrow(R) != exp_rows) {
    stop(sprintf("`%s` must have (m + 1) * q rows.", arg))
  }
  if (ncol(R) < 1) {
    stop(sprintf("`%s` must have at least one column.", arg))
  }

  R
}

#' Restricted Residual Blocks
#'
#' Computes residuals from each candidate regime coefficient block for use in
#' restricted break-date updating.
rresid <- function(maty, matz, nbeta, q, n, m) {
  bigt <- nrow(maty) / n
  bigvec2 <- matrix(0, nrow = bigt * (m + 1), ncol = n)
  for (i in seq_len(m + 1)) {
    beta_i <- nbeta[((i - 1) * q + 1):(i * q), , drop = FALSE]
    bigvec2[((i - 1) * bigt + 1):(i * bigt), ] <- matrix(as.vector(maty - matz %*% beta_i), ncol = n, byrow = TRUE)
  }

  bigvec2
}

#' Restricted Two-Segment Partition Search
#'
#' Finds the best split for a segment using residual blocks from restricted
#' coefficient estimates.
parti2 <- function(start, b1, b2, last, bigvec2, bigt, n) {
  dvec <- matrix(0, nrow = bigt, ncol = 1)
  for (j in b1:b2) {
    left_res <- bigvec2[start:j, , drop = FALSE]
    right_res <- bigvec2[(bigt + j + 1):(bigt + last), , drop = FALSE]
    m1 <- crossprod(left_res) / (j - start + 1)
    m2 <- crossprod(right_res) / (last - j)
    dvec[j, 1] <- ((last - start + 1) * n / 2) * (log(2 * pi) + 1) +
      ((j - start + 1) / 2) * slogdet(m1) +
      ((last - j) / 2) * slogdet(m2)
  }

  candidates <- dvec[b1:b2, 1]
  list(
    optmle = min(candidates),
    dx = as.integer((b1 - 1) + which.min(candidates))
  )
}

#' Restricted Dynamic Programming Search
#'
#' Updates break dates conditional on restricted coefficient estimates using
#' dynamic programming over residual covariance objective values.
dating_M2 <- function(bigvec2, h, m, n, bigt) {
  datevec <- matrix(0, nrow = m, ncol = m)
  optdat <- matrix(0, nrow = bigt, ncol = m)
  optmle <- matrix(0, nrow = bigt, ncol = m)
  dvec <- matrix(0, nrow = bigt, ncol = 1)
  global <- matrix(0, nrow = m, ncol = 1)
  bigvec <- matrix(0, nrow = bigt * (bigt + 1) / 2, ncol = 1)

  if (m == 1) {
    one <- parti2(1, h, bigt - h, bigt, bigvec2, bigt, n)
    datevec[1, 1] <- one$dx
    global[1, 1] <- one$optmle
  } else {
    for (j1 in (2 * h):bigt) {
      one <- parti2(1, h, j1 - h, j1, bigvec2[1:(2 * bigt), , drop = FALSE], bigt, n)
      optmle[j1, 1] <- one$optmle
      optdat[j1, 1] <- one$dx
    }
    global[1, 1] <- optmle[bigt, 1]
    datevec[1, 1] <- optdat[bigt, 1]

    for (ib in 2:m) {
      if (ib == m) {
        jlast <- bigt
        for (jb in (ib * h):(jlast - h)) {
          rows_j <- (bigt * m + jb + 1):(bigt * (m + 1))
          res_j <- bigvec2[rows_j, , drop = FALSE]
          pmle <- crossprod(res_j) / (bigt - jb)
          dvec[jb, 1] <- optmle[jb, ib - 1] +
            (((bigt - jb) * n / 2) * (log(2 * pi) + 1)) +
            (((bigt - jb) / 2) * slogdet(pmle))
        }
        candidates <- dvec[(ib * h):(jlast - h), 1]
        optmle[jlast, ib] <- min(candidates)
        optdat[jlast, ib] <- as.integer((ib * h - 1) + which.min(candidates))
      } else {
        for (jlast in ((ib + 1) * h):bigt) {
          for (jb in (ib * h):(jlast - h)) {
            rows_j <- (bigt * ib + jb + 1):(bigt * ib + jlast)
            res_j <- bigvec2[rows_j, , drop = FALSE]
            pmle <- crossprod(res_j) / (jlast - jb)
            dvec[jb, 1] <- optmle[jb, ib - 1] +
              (((jlast - jb) * n / 2) * (log(2 * pi) + 1)) +
              (((jlast - jb) / 2) * slogdet(pmle))
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

#' Iterative Restricted Break Estimation
#'
#' Alternates between restricted coefficient estimation and break-date updating
#' until the break dates converge.
est <- function(maty, matz, n, m, bigt, trm, R = NULL,
                brbeta = TRUE, brv = FALSE,
                tol = 1e-6, max_iter = 1000, max_outer = 1000,
                initial_breaks = NULL) {
  q <- ncol(matz)
  R <- valR(R, m = m, q = q)
  h <- as.integer(round(trm * bigt))

  if (is.null(initial_breaks)) {
    unrestricted <- dating_MLE(maty, matz, n = n, h = h, m = m, bigt = bigt, tol = tol, max_iter = max_iter)
    br <- as.integer(unrestricted$datevec[, m])
  } else {
    unrestricted <- NULL
    br <- valbrk(initial_breaks, m = m, bigt = bigt, h = h)
  }
  start_br <- br

  restricted <- r_estim(
    maty, matz, bigt, n, m, br, R,
    brbeta = brbeta,
    brv = brv,
    tol = tol,
    max_iter = max_iter
  )

  outer_iter <- 0
  repeat {
    bigvec2 <- rresid(maty, matz, restricted$nbeta, q, n, m)
    redated <- dating_M2(bigvec2, h, m, n, bigt)
    new_br <- as.integer(redated$datevec[, m])

    if (identical(new_br, br)) {
      return(list(
        dx = br,
        rbeta = restricted$nbeta,
        rvv = restricted$nvv,
        initial_breaks = start_br,
        unrestricted = unrestricted,
        redated = redated,
        iterations = outer_iter,
        r_estim = restricted
      ))
    }

    br <- new_br
    restricted <- r_estim(
      maty, matz, bigt, n, m, br, R,
      brbeta = brbeta,
      brv = brv,
      tol = tol,
      max_iter = max_iter
    )
    outer_iter <- outer_iter + 1
    if (outer_iter >= max_outer) {
      warning("`est()` reached `max_outer` before break dates converged.")
      return(list(
        dx = br,
        rbeta = restricted$nbeta,
        rvv = restricted$nvv,
        initial_breaks = start_br,
        unrestricted = unrestricted,
        redated = redated,
        iterations = outer_iter,
        r_estim = restricted
      ))
    }
  }
}

#' Block-Diagonal Weight Matrix
#'
#' Constructs the full stacked inverse-covariance weight matrix from
#' segment-specific covariance estimates.
bwt <- function(endpoints, n, bigt, vv) {
  out <- matrix(0, nrow = bigt * n, ncol = bigt * n)
  for (k in seq_len(length(endpoints) - 1)) {
    i <- endpoints[[k]] + 1
    j <- endpoints[[k + 1]]
    idx <- stacked_rows(i, j, n)
    v_k <- vv[((k - 1) * n + 1):(k * n), , drop = FALSE]
    out[idx, idx] <- kronecker(diag(j - i + 1), pinv(v_k))
  }
  out
}

#' Segment Covariance Updates
#'
#' Computes covariance matrices for each segment from residuals under a common
#' coefficient vector in the current restricted parameterization.
segvv <- function(x, y, b, endpoints, n) {
  vv <- matrix(0, nrow = (length(endpoints) - 1) * n, ncol = n)
  for (k in seq_len(length(endpoints) - 1)) {
    i <- endpoints[[k]] + 1
    j <- endpoints[[k + 1]]
    idx <- stacked_rows(i, j, n)
    res <- y[idx, , drop = FALSE] - x[idx, , drop = FALSE] %*% b
    umat <- matrix(as.vector(res), ncol = n, byrow = TRUE)
    vv[((k - 1) * n + 1):(k * n), ] <- t(umat) %*% umat / (j - i + 1)
  }
  vv
}
