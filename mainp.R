# Main procedure for the general multivariate structural break model.

#' Main Procedure for Multivariate Structural Break Models
#'
#' Estimates multivariate structural break models with a fixed maximum number of
#' breaks. The procedure computes unrestricted break estimates, optional
#' restricted estimates, likelihood-ratio or HAC-robust test statistics,
#' Bai-Perron critical-value approximations, sequential tests, and break-date
#' confidence intervals.
#'
#' @usage
#' mainp(y, z, m, S, R = NULL, trim = 0.15, verbose = FALSE,
#'   brbeta = TRUE, brv = FALSE, hac = FALSE,
#'   hetq = FALSE, prewhit = FALSE,
#'   sumq = NULL,
#'   restart = c("none", "brbp", "brbp_residual"),
#'   restart_B = 20, restart_init = NULL,
#'   restart_verbose = FALSE)
#'
#' @param y Numeric matrix, vector, or data frame of dependent variables with
#'   \eqn{T} rows and \eqn{n} columns. A vector is treated as a one-equation
#'   system.
#' @param z Numeric matrix, vector, or data frame of regressors with \eqn{T}
#'   rows. Include an intercept column if desired.
#' @param m Integer number of breaks to estimate.
#' @param S Selection matrix mapping equation-stacked raw regressors into the
#'   transformed system design. For observation \eqn{t}, the transformed design
#'   row block is \eqn{(I_n \otimes z_t')S}.
#' @param R Optional restriction/reparameterization matrix for the stacked
#'   regime coefficients. If supplied, the restricted coefficient vector is
#'   estimated under \eqn{\beta = Rb}.
#' @param trim Trimming parameter. The minimum segment length is
#'   \code{round(trim * nrow(y))}.
#' @param verbose Logical. If \code{TRUE}, print a short estimation summary.
#' @param brbeta Logical. If \code{TRUE}, allow regression
#'   coefficients to change across regimes.
#' @param brv Logical. If \code{TRUE}, allow the disturbance
#'   covariance matrix to change across regimes.
#' @param hac Logical. If \code{TRUE}, use HAC-robust F/Wald statistics for
#'   coefficient-break tests. If \code{FALSE}, use likelihood-ratio statistics.
#' @param hetq Logical. If \code{TRUE}, the regressor distribution is allowed
#'   to change across regimes (Gauss's \code{hetq}); if \code{FALSE}, pooled
#'   regressor moments are used in the break-date confidence intervals.
#' @param prewhit Logical. If \code{TRUE}, apply VAR(1)-style prewhitening
#'   in HAC covariance estimation.
#' @param sumq Optional integer giving the number of changing
#'   coefficient restrictions used for critical values. If \code{NULL}, all
#'   transformed regressors are treated as changing.
#' @param restart Character string, one of \code{"none"},
#'   \code{"brbp"}, or \code{"brbp_residual"}, controlling whether bootstrap
#'   restarting is used in the restricted branch.
#' @param restart_B Integer number of bootstrap restart iterations when
#'   \code{restart} is not \code{"none"}.
#' @param restart_init Optional integer vector of starting break dates
#'   for bootstrap restarting. If \code{NULL}, the usual restricted estimate is
#'   used as the first restart point.
#' @param restart_verbose Logical. If \code{TRUE}, print bootstrap restart
#'   progress.
#'
#' @return An object of class \code{"mbreak_mainp"}, a list with components:
#' \describe{
#'   \item{\code{setup}}{Input settings and derived dimensions.}
#'   \item{\code{transformed}}{The stacked dependent variable and transformed
#'     design matrix used internally.}
#'   \item{\code{unrestricted}}{Unrestricted testing and estimation results,
#'     including break dates, coefficients, coefficient standard errors,
#'     covariance estimates, confidence intervals, SupLR/SupF statistics, UDmax,
#'     WDmax, and sequential tests.}
#'   \item{\code{restricted}}{Restricted testing and estimation results when
#'     \code{R} is supplied; otherwise \code{NULL}.}
#' }
#'
#' @details
#' The unrestricted branch first estimates the optimal break dates for
#' \eqn{1,\dots,m} breaks, computes test statistics and critical values, and then
#' returns final estimates for \eqn{m} breaks. The restricted branch, when
#' requested through \code{R}, iterates between restricted coefficient estimation
#' and break-date updating. With \code{restart = "brbp"}, the
#' restricted branch applies observation-resampling bootstrap restarting. With
#' \code{restart = "brbp_residual"}, it applies residual-bootstrap
#' restarting. Restarting is applied before final estimation, tests, standard
#' errors, and break-date confidence intervals are computed.
#'
#' With \code{hac = FALSE}, tests are likelihood-ratio tests and may include
#' coefficient breaks, covariance breaks, or both. With \code{hac = TRUE}, the
#' reported robust statistic is a coefficient-break Wald/F statistic; robust
#' covariance-only tests are not implemented.
#'
#' @references
#' Qu, Z., and Perron, P. (2007). Estimating and testing structural changes
#' in multivariate regressions. \emph{Econometrica}, 75(2), 459-502.
#'
#' @examples
#' set.seed(1)
#' bigt <- 60
#' z <- cbind(1, rnorm(bigt))
#' beta <- rbind(c(1, 0.5), c(2, -0.2), c(0.5, 0.8))
#' y <- matrix(NA_real_, nrow = bigt, ncol = 1)
#' y[1:20, ] <- z[1:20, ] %*% beta[1, ] + rnorm(20, sd = 0.4)
#' y[21:40, ] <- z[21:40, ] %*% beta[2, ] + rnorm(20, sd = 0.4)
#' y[41:60, ] <- z[41:60, ] %*% beta[3, ] + rnorm(20, sd = 0.4)
#' fit <- mainp(
#'   y = y,
#'   z = z,
#'   m = 2,
#'   S = diag(ncol(z)),
#'   trim = 0.15
#' )
#' fit$unrestricted$estimation$breaks
#'
#' @export
mainp <- function(y, z, m, S, R = NULL, trim = 0.15, verbose = FALSE,
                  brbeta = TRUE, brv = FALSE,
                  hac = FALSE, hetq = FALSE,
                  prewhit = FALSE, sumq = NULL,
                  restart = c("none", "brbp", "brbp_residual"),
                  restart_B = 20, restart_init = NULL,
                  restart_verbose = FALSE) {
  y <- asmat(y, "y")
  z <- asmat(z, "z")
  S <- asmat(S, "S", allow_vector = FALSE)
  n <- ncol(y)
  bigt <- nrow(y)
  m <- as.integer(m)
  trm <- as.numeric(trim)
  restart <- match.arg(restart)
  restart_B <- as.integer(restart_B)
  if (length(restart_B) != 1 || is.na(restart_B) || restart_B < 0) {
    stop("`restart_B` must be a non-negative integer.", call. = FALSE)
  }
  if (restart != "none" && is.null(R)) {
    stop("`restart` requires a restriction matrix `R`.", call. = FALSE)
  }

  transformed <- transf(y, z, S)
  maty <- transformed$maty
  matz <- transformed$matz
  q <- ncol(matz)
  if (!is.null(R)) {
    R <- valR(R, m = m, q = q)
  }
  if (is.null(sumq)) {
    sumq <- q
  }
  sumq <- as.integer(sumq)
  h <- as.integer(round(trm * bigt))

  unrestricted <- unrm(
    maty = maty,
    matz = matz,
    n = n,
    m = m,
    bigt = bigt,
    trm = trm,
    h = h,
    sumq = sumq,
    brbeta = brbeta,
    brv = brv,
    hac = hac,
    hetq = hetq,
    prewhit = prewhit
  )

  restricted <- NULL
  if (!is.null(R)) {
    restricted <- resm(
      maty = maty,
      matz = matz,
      n = n,
      m = m,
      bigt = bigt,
      trm = trm,
      R = R,
      brbeta = brbeta,
      brv = brv,
      hac = hac,
      hetq = hetq,
      prewhit = prewhit,
      restart = restart,
      restart_B = restart_B,
      restart_init = restart_init,
      restart_verbose = restart_verbose
    )
  }

  result <- list(
    setup = list(
      m = m,
      trim = trm,
      h = h,
      bigt = bigt,
      n = n,
      sumq = sumq,
      brbeta = brbeta,
      brv = brv,
      hac = hac,
      hetq = hetq,
      prewhit = prewhit,
      restart = restart,
      restart_B = restart_B,
      restart_init = restart_init,
      restart_verbose = restart_verbose,
      S = S,
      R = R
    ),
    unrestricted = unrestricted,
    restricted = restricted
  )
  class(result) <- c("mbreak_mainp", class(result))

  if (verbose) {
    print(result)
  }

  result
}

#' Print an mbreak Main-Procedure Result
#'
#' Prints the main test statistics, break-date estimates, confidence
#' intervals, coefficient estimates, and covariance estimates from an object
#' returned by \code{\link{mainp}}.
print.mbreak_mainp <- function(x, ...) {
  pheader(x$setup)

  if (!is.null(x$unrestricted)) {
    pbranch(
      branch = x$unrestricted,
      title = "Output from the unrestricted structural change procedures",
      n = x$setup$n,
      restricted = FALSE
    )
  }

  if (!is.null(x$restricted)) {
    pbranch(
      branch = x$restricted,
      title = "Output from the restricted structural change procedures",
      n = x$setup$n,
      restricted = TRUE
    )
  }

  invisible(x)
}

#' Print Main-Procedure Result Header
#'
#' Prints the compact header used by the S3 print method for `mainp()`.
pheader <- function(setup) {
  cat("mbreak mainp result\n")
  cat("------------------------------------------------------------------------\n")
  cat("The basic specifications for testing and estimation:\n")
  cat("-------------------------------------------\n")
  cat(sprintf("(1) M = %d\n", setup$m))
  cat(sprintf("(2) Trimming = %.3f\n", setup$trim))
  cat(sprintf("(3) T = %d\n", setup$bigt))
  cat(sprintf("(4) brbeta = %s, brv = %s\n",
              setup$brbeta, setup$brv))
  cat(sprintf("(5) The number of coefficients in each regime is: %d\n",
              setup$sumq))
  cat("-------------------------------------------\n")
  cat(sprintf("hac = %s, hetq = %s, prewhit = %s\n",
              setup$hac, setup$hetq, setup$prewhit))
  cat("------------------------------------------------------------------------\n")
  invisible(NULL)
}

#' Print One Main-Procedure Branch
#'
#' Prints the testing and estimation output for either the unrestricted or
#' restricted branch of a `mainp()` result.
pbranch <- function(branch, title, n, restricted) {
  cat("\n")
  cat(title, "\n")
  cat("***************************************************\n")
  cat("\n")

  ptest(branch$testing, restricted = restricted)
  pestim(branch$estimation, n = n)

  invisible(NULL)
}

#' Print Main-Procedure Testing Output
#'
#' Prints SupLR, sequential, UDmax, WDmax, or restricted test statistics.
ptest <- function(testing, restricted) {
  if (is.null(testing)) {
    return(invisible(NULL))
  }

  cat("Output from the testing procedure\n")
  cat("__________________________________\n")
  cat("\n")

  if (restricted) {
    if (!is.null(testing$statistic)) {
      cat(sprintf("The value of the restricted test is: %.3f\n",
                  testing$statistic))
    }
    pcv(testing$critical_values)
    cat("\n")
    return(invisible(NULL))
  }

  if (!is.null(testing$sup_lr)) {
    for (i in seq_along(testing$sup_lr)) {
      cat(sprintf("The SupLR test for 0 versus %d break(s) is: %.3f\n",
                  i, testing$sup_lr[i]))
    }
  }

  pcv(testing$critical_values)

  if (!is.null(testing$udmax)) {
    cat(sprintf("The UDmax statistic is: %.3f\n", testing$udmax))
    pcv(testing$udmax_critical_values,
                                label = "UDmax critical values")
  }

  if (!is.null(testing$wdmax)) {
    if (length(testing$wdmax) == 1) {
      cat(sprintf("The WDmax statistic is: %.3f\n", testing$wdmax))
    } else {
      pcv(testing$wdmax, label = "WDmax statistics")
    }
    pcv(testing$wdmax_critical_values,
                                label = "WDmax critical values")
  }

  if (!is.null(testing$sequential)) {
    cat("Sequential tests:\n")
    for (item in testing$sequential) {
      cat(sprintf("The sequential test for %d versus %d break(s) is: %.3f\n",
                  item$l, item$l + 1, item$statistic))
      pcv(item$critical_values)
    }
  }

  cat("\n")
  invisible(NULL)
}

#' Print Critical Values
#'
#' Prints a vector or matrix of critical values when available.
pcv <- function(values, label = "Critical values") {
  if (is.null(values)) {
    return(invisible(NULL))
  }

  cat(label, ":\n", sep = "")
  print(round(values, 3))
  cat("\n")

  invisible(NULL)
}

#' Print Main-Procedure Estimation Output
#'
#' Prints estimated break dates, confidence intervals, coefficients, and
#' covariance matrices.
pestim <- function(estimation, n) {
  if (is.null(estimation)) {
    return(invisible(NULL))
  }

  cat("Output from the estimation procedure\n")
  cat("__________________________________\n")
  cat("\n")

  pbr(estimation$breaks, estimation$ci_breaks)
  pcoef(
    estimation$byregime,
    estimation$coefficient_se_by_regime
  )
  pcov(estimation$covariance, n = n)

  invisible(NULL)
}

#' Print Break-Date Estimates and Confidence Intervals
#'
#' Prints the estimated break dates and their 95 and 90 percent confidence
#' intervals in the style of rbreak.
pbr <- function(br, ci) {
  if (is.null(br) || length(br) == 0) {
    return(invisible(NULL))
  }

  cat("--------------------------------------------------------------\n")
  cat("The estimates and the confidence intervals for the break dates\n")
  cat("--------------------------------------------------------------\n")
  for (i in seq_along(br)) {
    cat(sprintf("The estimate of the %dth break is: %d\n", i, br[i]))
    if (!is.null(ci) && nrow(ci) >= i && ncol(ci) >= 4) {
      cat(sprintf("The 95%% C.I. for the %dth break is: [%d, %d]\n",
                  i, ci[i, 1], ci[i, 2]))
      cat(sprintf("The 90%% C.I. for the %dth break is: [%d, %d]\n",
                  i, ci[i, 3], ci[i, 4]))
    }
  }
  cat("********************************************************\n")
  cat("\n")

  invisible(NULL)
}

#' Print Coefficients by Regime
#'
#' Prints regime-specific coefficient estimates and, when available, standard
#' errors.
pcoef <- function(coefficients, standard_errors = NULL) {
  if (is.null(coefficients)) {
    return(invisible(NULL))
  }

  cat("The estimated coefficients are:\n")
  for (i in seq_len(nrow(coefficients))) {
    cat(sprintf("For the %dth regime:\n", i))
    if (is.null(standard_errors)) {
      print(round(coefficients[i, , drop = FALSE], 5))
    } else {
      out <- cbind(
        estimate = as.vector(coefficients[i, ]),
        std.error = as.vector(standard_errors[i, ])
      )
      rownames(out) <- paste0("coef", seq_len(ncol(coefficients)))
      print(round(out, 5))
    }
  }
  cat("\n")

  invisible(NULL)
}

#' Print Covariance Matrices by Regime
#'
#' Prints stacked covariance estimates as one block per regime.
pcov <- function(covariance, n) {
  if (is.null(covariance) || is.na(n) || n <= 0) {
    return(invisible(NULL))
  }

  nreg <- floor(nrow(covariance) / n)
  if (nreg == 0) {
    return(invisible(NULL))
  }

  cat("The estimated covariance matrices of the errors are:\n")
  for (i in seq_len(nreg)) {
    rows <- ((i - 1) * n + 1):(i * n)
    cat(sprintf("For the %dth regime:\n", i))
    print(round(covariance[rows, , drop = FALSE], 5))
  }
  cat("\n")

  invisible(NULL)
}

#' Unrestricted Main-Procedure Branch
#'
#' Runs unrestricted estimation, testing, critical-value calculations,
#' sequential tests, and confidence intervals for `mainp()`.
unrm <- function(maty, matz, n, m, bigt, trm, h, sumq,
                                   brbeta, brv, hac,
                                   hetq, prewhit) {
  q <- ncol(matz)
  q_test <- qtest(n, sumq, brbeta, brv, hac)
  lrtest <- matrix(NA_real_, nrow = m, ncol = 1)
  breaks_by_m <- vector("list", m)
  coefficients_by_m <- vector("list", m)

  for (i in seq_len(m)) {
    R_i <- defR(
      i, q,
      brbeta = brbeta,
      brv = brv
    )
    fit_i <- est(
      maty, matz, n, i, bigt, trm, R_i,
      brbeta = brbeta,
      brv = brv
    )
    final_i <- r_estim(
      maty, matz, bigt, n, i, fit_i$dx, R_i,
      brbeta = brbeta,
      brv = brv
    )
    breaks_by_m[[i]] <- fit_i$dx
    coefficients_by_m[[i]] <- final_i$nbeta

    lrtest[i, 1] <- tstat(
      maty = maty,
      matz = matz,
      m = i,
      n = n,
      br = fit_i$dx,
      beta = final_i$nbeta,
      bigt = bigt,
      brbeta = brbeta,
      brv = brv,
      hac = hac,
      prewhit = prewhit,
      R = R_i
    )
  }

  cvs <- cvlrmat(m, q_test, trm)
  colnames(cvs) <- c("90%", "95%", "97.5%", "99%")
  rownames(cvs) <- paste0("m=", seq_len(m))

  sequential <- NULL
  if (m > 1) {
    sequential <- vector("list", m - 1)
    for (i in seq_len(m - 1)) {
      stat_i <- if (is.na(q_test)) {
        NA_real_
      } else if (!hac) {
        seqlr(
          maty, matz, n, breaks_by_m[[i]], bigt, h,
          brbeta = brbeta,
          brv = brv,
          trm = trm
        )
      } else {
        seqf(
          maty, matz, n, breaks_by_m[[i]], bigt, h,
          prewhit = prewhit,
          trm = trm
        )
      }
      sequential[[i]] <- list(
        l = i,
        statistic = stat_i,
        critical_values = seq_na(i, q_test, trm)
      )
    }
  }

  R_m <- defR(
    m, q,
    brbeta = brbeta,
    brv = brv
  )
  final <- r_estim(
    maty, matz, bigt, n, m, breaks_by_m[[m]], R_m,
    brbeta = brbeta,
    brv = brv
  )
  coefficient_se <- coefse(
    maty = maty,
    matz = matz,
    n = n,
    m = m,
    bigt = bigt,
    br = breaks_by_m[[m]],
    beta = final$nbeta,
    vv = final$nvv,
    R = R_m,
    hac = hac,
    brv = brv,
    prewhit = prewhit
  )
  ci <- interval2(maty, matz, m, n, bigt, breaks_by_m[[m]], final$nbeta, final$nvv,
                  hetq = hetq,
                  hac = hac,
                  brv = brv,
                  brbeta = brbeta,
                  prewhit = prewhit)

  list(
    testing = list(
      sup_lr = as.numeric(lrtest[, 1]),
      critical_values = cvs,
      udmax = ud_na(lrtest[, 1]),
      udmax_critical_values = cvud_na(q_test, trm),
      wdmax = wd_na(lrtest[, 1], q_test, trm),
      wdmax_critical_values = cvwd_na(q_test, trm),
      sequential = sequential,
      changing_parameters = q_test
    ),
    estimation = list(
      breaks = breaks_by_m[[m]],
      coefficients = final$nbeta,
      byregime = byregime(final$nbeta, m, q),
      coefficient_vcov = coefficient_se$vcov,
      coefficient_se = coefficient_se$se,
      coefficient_se_by_regime = coefficient_se$se_by_regime,
      theta_vcov = coefficient_se$vcov_theta,
      theta_se = coefficient_se$se_theta,
      covariance = final$nvv,
      ci_breaks = ci,
      breaks_by_m = breaks_by_m,
      coefficients_by_m = coefficients_by_m
    )
  )
}

#' Restricted Main-Procedure Branch
#'
#' Runs restricted estimation, testing, and confidence intervals for `mainp()`
#' when a restriction matrix is supplied.
resm <- function(maty, matz, n, m, bigt, trm, R,
                                 brbeta, brv, hac,
                                 hetq, prewhit,
                                 restart = "none", restart_B = 20,
                                 restart_init = NULL,
                                 restart_verbose = FALSE) {
  q <- ncol(matz)
  R <- valR(R, m = m, q = q)
  info <- NULL
  if (restart %in% c("brbp", "brbp_residual")) {
    r_res <- brbp(
      maty = maty,
      matz = matz,
      n = n,
      m = m,
      bigt = bigt,
      trm = trm,
      R = R,
      method = restart,
      initial_breaks = restart_init,
      B = restart_B,
      brbeta = brbeta,
      brv = brv,
      verbose = restart_verbose
    )
    fit <- r_res$fit
    info <- r_res$diagnostics
  } else {
    fit <- est(
      maty, matz, n, m, bigt, trm, R,
      brbeta = brbeta,
      brv = brv
    )
  }

  final <- r_estim(
    maty, matz, bigt, n, m, fit$dx, R,
    brbeta = brbeta,
    brv = brv
  )
  coefficient_se <- coefse(
    maty = maty,
    matz = matz,
    n = n,
    m = m,
    bigt = bigt,
    br = fit$dx,
    beta = final$nbeta,
    vv = final$nvv,
    R = R,
    hac = hac,
    brv = brv,
    prewhit = prewhit
  )
  stat <- tstat(
    maty = maty,
    matz = matz,
    m = m,
    n = n,
    br = fit$dx,
    beta = final$nbeta,
    bigt = bigt,
    brbeta = brbeta,
    brv = brv,
    hac = hac,
    prewhit = prewhit,
    R = R
  )
  ci <- interval2(maty, matz, m, n, bigt, fit$dx, final$nbeta, final$nvv,
                  hetq = hetq,
                  hac = hac,
                  brv = brv,
                  brbeta = brbeta,
                  prewhit = prewhit)
  cq <- qres(
    n = n, m = m, R = R,
    brbeta = brbeta,
    brv = brv,
    hac = hac
  )

  list(
    testing = list(
      statistic = stat,
      critical_values = cvlr_na(m, cq, trm),
      changing_parameters = cq
    ),
    estimation = list(
      breaks = fit$dx,
      coefficients = final$nbeta,
      byregime = byregime(final$nbeta, m, q),
      coefficient_vcov = coefficient_se$vcov,
      coefficient_se = coefficient_se$se,
      coefficient_se_by_regime = coefficient_se$se_by_regime,
      theta_vcov = coefficient_se$vcov_theta,
      theta_se = coefficient_se$se_theta,
      covariance = final$nvv,
      ci_breaks = ci,
      initial_breaks = fit$initial_breaks,
      iterations = fit$iterations,
      restart = info
    )
  )
}

#' Main-Procedure Test Statistic Dispatcher
#'
#' Chooses and computes the appropriate LR or HAC-robust test statistic for a
#' fitted `mainp()` break model.
tstat <- function(maty, matz, m, n, br, beta, bigt,
                                brbeta, brv, hac,
                                prewhit, R = NULL) {
  if (!hac) {
    return(suplr(
      maty, matz, m, n, br, beta, bigt,
      brbeta = brbeta,
      brv = brv
    ))
  }

  if (brv && !brbeta) {
    warning(
      "HAC robust tests for covariance-only breaks are not available; returning `NA` for the test statistic.",
      call. = FALSE
    )
    return(NA_real_)
  }

  if (brv && brbeta) {
    warning(
      "With `hac = TRUE` and `brv = TRUE`, the reported robust test is for coefficient breaks only; covariance-break critical values are not available under serial correlation.",
      call. = FALSE
    )
  }

  if (is.null(R) || identical(R, diag(nrow(R)))) {
    return(supft(
      maty, matz, m, n, br, beta, bigt,
      prewhit = prewhit
    ))
  }

  suprft(
    maty, matz, m, n, br, beta, R, bigt,
    prewhit = prewhit
  )
}

#' Default Restriction Matrix for Break Type
#'
#' Constructs the default reparameterization matrix implied by the requested
#' coefficient-break and covariance-break specification.
defR <- function(m, q, brbeta, brv) {
  if (brbeta) {
    return(diag((m + 1) * q))
  }
  if (brv) {
    return(kronecker(matrix(1, nrow = m + 1, ncol = 1), diag(q)))
  }
  diag((m + 1) * q)
}

#' Number of Parameters for Unrestricted Critical Values
#'
#' Computes the number of changing parameters used to select Bai-Perron critical
#' values for the unrestricted branch.
qtest <- function(n, sumq, brbeta, brv, hac = FALSE) {
  if (hac) {
    if (brbeta) {
      return(as.integer(sumq))
    }
    return(NA_integer_)
  }

  if (brv && !brbeta) {
    return(as.integer(n * (n + 1) / 2))
  }
  if (brv && brbeta) {
    return(as.integer(sumq + n * (n + 1) / 2))
  }
  as.integer(sumq)
}

#' Number of Parameters for Restricted Critical Values
#'
#' Computes the implied number of changing parameters for a restricted model.
qres <- function(n, m, R, brbeta, brv, hac = FALSE) {
  cq <- nrow(R) / (m + 1) - (nrow(R) - ncol(R)) / m
  if (hac) {
    if (!brbeta) {
      return(NA_integer_)
    }
    if (abs(cq - round(cq)) > .Machine$double.eps^0.5) {
      stop("The implied number of changing coefficient parameters is not integer.")
    }
    return(as.integer(round(cq)))
  }

  if (brv && !brbeta) {
    return(as.integer(n * (n + 1) / 2))
  }

  cq <- n * (n + 1) / 2 * brv + cq
  if (abs(cq - round(cq)) > .Machine$double.eps^0.5) {
    stop("The implied number of changing parameters is not integer.")
  }
  as.integer(round(cq))
}

#' Critical Values or Missing Values
#'
#' Returns SupLR/SupF critical values when the parameter count is available,
#' otherwise returns missing values.
cvlr_na <- function(k, q, trm) {
  if (length(q) != 1 || is.na(q)) {
    return(rep(NA_real_, 4))
  }
  cvlr(k, q, trm)
}

#' Critical-Value Matrix or Missing Values
#'
#' Builds a matrix of SupLR/SupF critical values for break counts up to
#' `k_max`, or an `NA` matrix if unavailable.
cvlrmat <- function(k_max, q, trm) {
  if (length(q) != 1 || is.na(q)) {
    return(matrix(NA_real_, nrow = k_max, ncol = 4))
  }
  t(vapply(seq_len(k_max), function(i) cvlr(i, q, trm), numeric(4)))
}

#' Sequential Critical Values or Missing Values
#'
#' Returns sequential-test critical values when the parameter count is
#' available, otherwise returns missing values.
seq_na <- function(l, q, trm) {
  if (length(q) != 1 || is.na(q)) {
    return(rep(NA_real_, 4))
  }
  supseq(l, q, trm)
}

#' UDmax or Missing Value
#'
#' Computes UDmax when all input statistics are finite.
ud_na <- function(lrtest) {
  if (length(lrtest) == 0 || anyNA(lrtest) || !all(is.finite(lrtest))) {
    return(NA_real_)
  }
  udmax(lrtest)
}

#' UDmax Critical Values or Missing Values
#'
#' Returns UDmax critical values when the parameter count is available.
cvud_na <- function(q, trm) {
  if (length(q) != 1 || is.na(q) || q > 10) {
    return(rep(NA_real_, 4))
  }
  cvud(q, trm)
}

#' WDmax or Missing Values
#'
#' Computes WDmax when test statistics and critical-value parameters are
#' available.
wd_na <- function(lrtest, q, trm) {
  if (length(q) != 1 || is.na(q) ||
      length(lrtest) == 0 || anyNA(lrtest) || !all(is.finite(lrtest))) {
    return(rep(NA_real_, 4))
  }
  wdmax(lrtest, q, trm)
}

#' WDmax Critical Values or Missing Values
#'
#' Returns WDmax critical values when the parameter count is available.
cvwd_na <- function(q, trm) {
  if (length(q) != 1 || is.na(q)) {
    return(rep(NA_real_, 4))
  }
  cvwd(q, trm)
}

#' Coefficients by Regime
#'
#' Reshapes a stacked coefficient vector into a regime-by-coefficient matrix.
byregime <- function(beta, m, q) {
  matrix(as.vector(beta), nrow = m + 1, ncol = q, byrow = TRUE)
}

