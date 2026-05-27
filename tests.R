# Likelihood-ratio tests and asymptotic critical-value approximations.

#' SupLR Statistic for Fixed Break Dates
#'
#' Computes the likelihood-ratio statistic comparing a fitted break model at
#' fixed break dates with the no-break null model.
suplr <- function(maty, matx, m, n, br, rbeta, bigt,
                  brbeta = TRUE, brv = FALSE,
                  tol = 1e-6, max_iter = 1000) {
  br <- if (m == 0) integer() else br
  pmatx <- pzbar(matx, m, br, bigt)

  res <- maty - pmatx %*% rbeta
  res <- matrix(as.vector(res), ncol = n, byrow = TRUE)
  endpoints <- c(0, br, bigt)

  if (brv) {
    ulr <- 0
    for (k in seq_len(m + 1)) {
      i <- endpoints[[k]] + 1
      j <- endpoints[[k + 1]]
      res_k <- res[i:j, , drop = FALSE]
      vvar <- crossprod(res_k) / (j - i + 1)
      ulr <- ulr + nllk(j - i + 1, n, vvar)
    }
  } else if (brbeta) {
    vvar <- crossprod(res[seq_len(bigt), , drop = FALSE]) / bigt
    ulr <- nllk(bigt, n, vvar)
  } else {
    stop("At least one of `brbeta` or `brv` must be `TRUE` for `suplr()`.")
  }

  null_fit <- fgls(matx, maty, n = n, tol = tol, max_iter = max_iter)
  rlr <- nllk(bigt, n, null_fit$sigma)

  as.numeric(2 * (ulr - rlr))
}

#' Sequential LR Test Statistic
#'
#' Computes the sequential `l + 1 | l` likelihood-ratio statistic by inserting
#' one additional break into each admissible segment of the current partition.
seqlr <- function(maty, matz, n, br, bigt, h,
                  brbeta = TRUE, brv = FALSE,
                  trm = NULL, tol = 1e-6, max_iter = 1000) {
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

      if (brbeta && !brv) {
        fit <- est(
          segy, segz, n, 1, seg_bigt, trm, diag(q * 2),
          brbeta = TRUE,
          brv = FALSE,
          tol = tol,
          max_iter = max_iter
        )
        lrtest[is] <- suplr(
          segy, segz, 1, n, fit$dx, fit$rbeta, seg_bigt,
          brbeta = TRUE,
          brv = FALSE,
          tol = tol,
          max_iter = max_iter
        )
      } else if (!brbeta && brv) {
        fit <- est(
          segy, segz, n, 1, seg_bigt, trm,
          kronecker(matrix(1, nrow = 2), diag(q)),
          brbeta = FALSE,
          brv = TRUE,
          tol = tol,
          max_iter = max_iter
        )
        lrtest[is] <- suplr(
          segy, segz, 1, n, fit$dx, fit$rbeta, seg_bigt,
          brbeta = FALSE,
          brv = TRUE,
          tol = tol,
          max_iter = max_iter
        )
      } else if (brbeta && brv) {
        fit <- est(
          segy, segz, n, 1, seg_bigt, trm, diag(q * 2),
          brbeta = TRUE,
          brv = TRUE,
          tol = tol,
          max_iter = max_iter
        )
        lrtest[is] <- suplr(
          segy, segz, 1, n, fit$dx, fit$rbeta, seg_bigt,
          brbeta = TRUE,
          brv = TRUE,
          tol = tol,
          max_iter = max_iter
        )
      }
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

#' UDmax Test Statistic
#'
#' Computes the unweighted double-maximum statistic from a vector of SupLR or
#' SupF statistics over candidate break counts.
#'
#' @param lrtest Numeric vector of test statistics for \eqn{m=1,\dots,M}.
#'
#' @return The maximum value in \code{lrtest}.
#'
#' @export
udmax <- function(lrtest) {
  max(lrtest)
}

#' WDmax Test Statistic
#'
#' Computes the weighted double-maximum statistic using Bai-Perron
#' response-surface critical values.
#'
#' @param lrtest Numeric vector of test statistics for \eqn{m=1,\dots,M}.
#' @param q Integer number of changing parameters used for the critical values.
#' @param trm Numeric trimming parameter.
#'
#' @return Numeric vector of weighted double-maximum statistics at the 90%,
#'   95%, 97.5%, and 99% levels.
#'
#' @export
wdmax <- function(lrtest, q, trm) {
  cvs <- t(vapply(seq_along(lrtest), function(k) cvlr(k, q, trm), numeric(4)))
  out <- numeric(4)
  for (i in seq_along(lrtest)) {
    out <- pmax(out, lrtest[[i]] * (cvs[1, ] / cvs[i, ]))
  }
  out
}

#' Critical Values for SupLR/SupF Tests
#'
#' Computes Bai-Perron response-surface critical values for the sup test with a
#' fixed number of breaks.
#'
#' @param k Integer number of breaks under the alternative.
#' @param q Integer number of changing parameters.
#' @param trm Numeric trimming parameter.
#'
#' @return Numeric vector of critical values at the 90%, 95%, 97.5%, and 99%
#'   quantiles.
#'
#' @references
#' Bai, J., and Perron, P. (2003). Critical values for multiple structural
#' change tests. \emph{The Econometrics Journal}, 6(1), 72--78.
#'
#' @examples
#' cvlr(k = 2, q = 4, trm = 0.15)
#'
#' @export
cvlr <- function(k, q, trm) {
  cv90 <- (7.551 + 1.718 * q - 0.041 * q * q - 0.610 * k - 15.846 * trm + 0.025 * (q / trm)) *
    exp(0.338 / k - 0.014 / (trm * k))
  cv95 <- (8.238 + 1.756 * q - 0.043 * q * q - 0.659 * k - 15.436 * trm + 0.025 * (q / trm)) *
    exp(0.389 / k - 0.013 / (trm * k))
  cv975 <- (8.968 + 1.788 * q - 0.045 * q * q - 0.715 * k - 15.255 * trm + 0.025 * (q / trm)) *
    exp(0.426 / k - 0.013 / (trm * k))
  cv99 <- (9.879 + 1.771 * q - 0.042 * q * q - 0.777 * k - 14.551 * trm + 0.025 * (q / trm)) *
    exp(0.471 / k - 0.012 / (trm * k))

  k * c(cv90, cv95, cv975, cv99)
}

#' Critical Values for the UDmax Test
#'
#' Computes Bai-Perron response-surface critical values for the unweighted
#' double-maximum test.
#'
#' @param q Integer number of changing parameters.
#' @param trm Numeric trimming parameter.
#'
#' @return Numeric vector of critical values at the 90%, 95%, 97.5%, and 99%
#'   quantiles.
#'
#' @references
#' Bai, J., and Perron, P. (2003). Critical values for multiple structural
#' change tests. \emph{The Econometrics Journal}, 6(1), 72--78.
#'
#' @examples
#' cvud(q = 4, trm = 0.15)
#'
#' @export
cvud <- function(q, trm) {
  cv90 <- (6.917 + 2.930 * q - 9.275 * trm) * exp(-0.028 * q - 0.406 * trm)
  cv95 <- (8.228 + 3.095 * q - 9.644 * trm) * exp(-0.029 * q - 0.291 * trm)
  cv975 <- (9.436 + 3.304 * q - 9.301 * trm) * exp(-0.030 * q - 0.259 * trm)
  cv99 <- (11.211 + 3.366 * q - 7.279 * trm) * exp(-0.027 * q - 0.268 * trm)

  c(cv90, cv95, cv975, cv99)
}

#' Critical Values for the WDmax Test
#'
#' Computes Bai-Perron response-surface critical values for the weighted
#' double-maximum test.
#'
#' @param q Integer number of changing parameters.
#' @param trm Numeric trimming parameter.
#'
#' @return Numeric vector of critical values at the 90%, 95%, 97.5%, and 99%
#'   quantiles.
#'
#' @references
#' Bai, J., and Perron, P. (2003). Critical values for multiple structural
#' change tests. \emph{The Econometrics Journal}, 6(1), 72--78.
#'
#' @examples
#' cvwd(q = 4, trm = 0.15)
#'
#' @export
cvwd <- function(q, trm) {
  cv90 <- (7.316 + 3.128 * q - 8.624 * trm) * exp(-0.029 * q - 0.412 * trm)
  cv95 <- (9.039 + 3.318 * q - 9.969 * trm) * exp(-0.030 * q - 0.327 * trm)
  cv975 <- (10.703 + 3.465 * q - 11.119 * trm) * exp(-0.031 * q - 0.250 * trm)
  cv99 <- (13.189 + 3.346 * q - 11.870 * trm) * exp(-0.026 * q - 0.173 * trm)

  c(cv90, cv95, cv975, cv99)
}

#' Critical Values for Sequential Break Tests
#'
#' Computes Bai-Perron response-surface critical values for the sequential
#' \eqn{l+1 | l} test.
#'
#' @param l Integer number of breaks under the null.
#' @param q Integer number of changing parameters.
#' @param trm Numeric trimming parameter.
#'
#' @return Numeric vector of critical values at the 90%, 95%, 97.5%, and 99%
#'   quantiles.
#'
#' @references
#' Bai, J., and Perron, P. (2003). Critical values for multiple structural
#' change tests. \emph{The Econometrics Journal}, 6(1), 72--78.
#'
#' @examples
#' supseq(l = 1, q = 4, trm = 0.15)
#'
#' @export
supseq <- function(l, q, trm) {
  cv90 <- (8.397 + 3.702 * q - 0.209 * q * q + 0.317 * (l + 1) - 3.736 / (l + 1) - 11.596 * trm) *
    exp(-0.027 * q + 0.005 * q * q)
  cv95 <- (9.879 + 4.086 * q - 0.232 * q * q + 0.322 * (l + 1) - 3.687 / (l + 1) - 11.931 * trm) *
    exp(-0.039 * q + 0.006 * q * q)
  cv975 <- (11.424 + 4.435 * q - 0.261 * q * q + 0.300 * (l + 1) - 3.700 / (l + 1) - 12.292 * trm) *
    exp(-0.047 * q + 0.006 * q * q)
  cv99 <- (13.073 + 4.954 * q - 0.317 * q * q + 0.285 * (l + 1) - 3.419 / (l + 1) - 12.452 * trm) *
    exp(-0.058 * q + 0.008 * q * q)

  c(cv90, cv95, cv975, cv99)
}

#' Gaussian Log Likelihood at the Segment MLE
#'
#' Computes the normal log likelihood contribution for a segment using the
#' segment covariance estimate.
nllk <- function(seg_len, n, sigma) {
  -((seg_len * n) / 2) * (log(2 * pi) + 1) -
    (seg_len / 2) * logdet(sigma)
}
