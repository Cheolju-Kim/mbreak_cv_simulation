# Ordered-break model routines.

#' Ordered Break-Date Search
#'
#' Searches over ordered two-break candidates where one coefficient block
#' changes first and another block changes second.
brorder <- function(maty, matz, n, q1, q2, h, bigt, trun, tol = 1e-8, max_iter = 1000) {
  q <- ncol(matz)
  vecmle <- matrix(-1e12, nrow = bigt * bigt, ncol = 3 + 2 * q)
  zz <- cbind(matz, matz)
  prev_v <- NULL

  for (ii in (h + 1):(bigt - h)) {
    jj <- 0
    while ((ii + jj) <= bigt - h) {
      z1 <- zz
      z1[(ii * n + 1):nrow(zz), (q + 1):(q + q1)] <-
        zz[(ii * n + 1):nrow(zz), (q + 1):(q + q1)]
      z1[1:(ii * n), (q + 1):(q + q1)] <- 0

      z1[((ii + jj) * n + 1):nrow(zz), (q + q1 + 1):(2 * q)] <-
        zz[((ii + jj) * n + 1):nrow(zz), (q + q1 + 1):(2 * q)]
      z1[1:((ii + jj) * n), (q + q1 + 1):(2 * q)] <- 0

      fit <- ofgls(
        segx = z1,
        segy = maty,
        n = n,
        init_s = if (jj == 0) NULL else prev_v,
        tol = tol,
        max_iter = max_iter
      )
      prev_v <- fit$sigma

      row_i <- as.integer((ii - 1) * bigt + jj)
      vecmle[row_i, 1] <- ii
      vecmle[row_i, 2] <- ii + jj
      vecmle[row_i, 3] <- -bigt / 2 * logdet(fit$sigma)
      vecmle[row_i, 4:(3 + 2 * q)] <- as.vector(fit$beta)

      jj <- jj + 1
    }
  }

  vecmle2 <- vecmle
  valid_trunc <- (vecmle2[, 2] - vecmle2[, 1]) <= trun
  vecmle2[, 3] <- ifelse(valid_trunc, vecmle2[, 3], -1e12)

  best_unr <- vecmle[which.max(vecmle[, 3]), , drop = FALSE]
  best_tr <- vecmle2[which.max(vecmle2[, 3]), , drop = FALSE]
  H1like <- best_tr[1, 3]

  null_fit <- ofgls(matz, maty, n = n, init_s = NULL, tol = tol, max_iter = max_iter)
  H0like <- -bigt / 2 * logdet(null_fit$sigma)
  lr <- 2 * (H1like - H0like)

  beta1 <- matrix(best_unr[1, 4:(3 + q)], ncol = 1)
  beta2 <- matrix(best_unr[1, (4 + q):(3 + 2 * q)], ncol = 1) + beta1

  list(
    brdate = as.integer(best_unr[1, 1:2]),
    beta1 = beta1,
    beta2 = beta2,
    umat = null_fit$umat,
    lr = as.numeric(lr),
    truncated_brdate = as.integer(best_tr[1, 1:2]),
    H1like = as.numeric(H1like),
    H0like = as.numeric(H0like),
    candidates = vecmle
  )
}

#' Ordered-Break Confidence Bounds
#'
#' Simulates the limiting distribution used to construct confidence bounds for
#' the two ordered break dates.
interval3 <- function(matz, maty, n, start, last, dx1, dx2, delt1, delt2, res,
                      segT, bigM, rep, hetq = TRUE,
                      step = 1000, interv = 0.02, seed = NULL,
                      verbose = FALSE) {
  if (!is.numeric(interv) || length(interv) != 1 || !is.finite(interv) || interv <= 0) {
    stop("`interv` must be a positive finite number.")
  }
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  sigma <- t(res) %*% res / nrow(res)
  if (hetq) {
    left_idx <- stacked_rows(start, dx1, n)
    right_idx <- stacked_rows(dx2 + 1, last, n)
    Qj <- t(matz[left_idx, , drop = FALSE]) %*%
      kronecker(diag(dx1 - start + 1), pinv(sigma)) %*%
      matz[left_idx, , drop = FALSE] / (dx1 - start + 1)
    Qjp1 <- t(matz[right_idx, , drop = FALSE]) %*%
      kronecker(diag(last - dx2), pinv(sigma)) %*%
      matz[right_idx, , drop = FALSE] / (last - dx2)

    P1j <- as.numeric(t(delt1) %*% Qj %*% delt1)
    P2j <- as.numeric(t(delt2) %*% Qj %*% delt2)
    P12j <- as.numeric(t(delt1) %*% Qj %*% delt2)
    P1jp1 <- as.numeric(t(delt1) %*% Qjp1 %*% delt1)
    P2jp1 <- as.numeric(t(delt2) %*% Qjp1 %*% delt2)
    P12jp1 <- as.numeric(t(delt1) %*% Qjp1 %*% delt2)
  } else {
    idx <- stacked_rows(start, last, n)
    Qj <- t(matz[idx, , drop = FALSE]) %*%
      kronecker(diag(last - start + 1), pinv(sigma)) %*%
      matz[idx, , drop = FALSE] / (last - start + 1)

    P1j <- as.numeric(t(delt1) %*% Qj %*% delt1)
    P2j <- as.numeric(t(delt2) %*% Qj %*% delt2)
    P12j <- as.numeric(t(delt1) %*% Qj %*% delt2)
    P1jp1 <- P1j
    P2jp1 <- P2j
    P12jp1 <- P12j
  }

  bigN <- step * bigM
  grids <- as.integer(round(1 / interv))
  gridN <- grids * bigM
  grid_step <- as.integer(round(step * interv))
  sv <- -bigM + (1 / grids) * seq.int(0, 2 * grids * bigM)
  sv_len <- length(sv)
  V1 <- rep(sv, each = sv_len)
  V2 <- rep(sv, times = sv_len)
  obd <- matrix(0, nrow = rep, ncol = 3)

  for (simu in seq_len(rep)) {
    if (verbose && simu %% 10 == 0) {
      message("rep: ", simu)
    }

    bigW <- matrix(stats::rnorm(bigN * 4), nrow = bigN, ncol = 4)
    bigW <- apply(bigW, 2, cumsum) / sqrt(step)
    W <- bigW[seq.int(1, by = grid_step, length.out = gridN), , drop = FALSE]
    sbv1 <- c(base::rev(W[, 1]), 0, W[, 4])
    sbv2 <- c(base::rev(W[, 2]), 0, W[, 3])
    Bv1 <- rep(sbv1, each = sv_len)
    Bv2 <- rep(sbv2, times = sv_len)

    if (hetq) {
      Z1 <- (V1 < 0) * (V2 <= 0) * (V1 <= V2) *
        ((Bv1 + V1 / 2) + sqrt(P2j / P1j) * Bv2 + V2 / 2 * (P2j + 2 * P12j) / P1j)
      Z2 <- (V1 < 0) * (V2 > 0) *
        (Bv1 + V1 / 2 + sqrt(P2jp1 / P1j) * Bv2 - V2 / 2 * (P2jp1 / P1j))
      Z3 <- (V1 >= 0) * (V2 > 0) * (V1 <= V2) *
        (sqrt(P1jp1 / P1j) * Bv1 + sqrt(P2jp1 / P1j) * Bv2 -
           (P2jp1 / P1j) * V2 / 2 - ((P1jp1 + 2 * P12jp1) / P1j) * V1 / 2)
    } else {
      Z1 <- (V1 < 0) * (V2 <= 0) * (V1 <= V2) *
        ((Bv1 + V1 / 2) + sqrt(P2j / P1j) * Bv2 + V2 / 2 * (P2j + 2 * P12j) / P1j)
      Z2 <- (V1 < 0) * (V2 > 0) *
        (Bv1 + V1 / 2 + sqrt(P2j / P1j) * Bv2 - V2 / 2 * (P2j / P1j))
      Z3 <- (V1 >= 0) * (V2 > 0) * (V1 <= V2) *
        (Bv1 + sqrt(P2j / P1j) * Bv2 -
           (P2j / P1j) * V2 / 2 - ((P1j + 2 * P12j) / P1j) * V1 / 2)
    }
    Z4 <- (V1 > V2) * (-1)
    Z <- Z1 + Z2 + Z3 + Z4
    best <- which.max(Z)
    obd[simu, ] <- c(V1[[best]], V2[[best]], Z[[best]])
  }

  obd[, 1] <- sort(obd[, 1])
  obd[, 2] <- sort(obd[, 2])
  count <- sum(abs(abs(obd[, 1]) - bigM) < 1e-4) + sum(abs(abs(obd[, 2]) - bigM) < 1e-4)
  if (count >= rep * 0.05) {
    warning("Some simulated ordered-break bounds reached `bigM`.")
  }

  q025 <- max(1, floor(rep * 0.025))
  q05 <- max(1, floor(rep * 0.05))
  q95 <- max(1, floor(rep * 0.95))
  q975 <- max(1, floor(rep * 0.975))

  bound <- matrix(0, nrow = 4, ncol = 2)
  bound[1, ] <- floor(c(dx1, dx2) - obd[q025, 1:2] / P1j) + c(1, 1)
  bound[2, ] <- floor(c(dx1, dx2) - obd[q05, 1:2] / P1j) + c(1, 1)
  bound[3, ] <- floor(c(dx1, dx2) - obd[q95, 1:2] / P1j)
  bound[4, ] <- floor(c(dx1, dx2) - obd[q975, 1:2] / P1j)
  bound
}

#' Ordered Two-Break Multivariate Structural Change Model
#'
#' Estimates an ordered two-break model in which one subset of coefficients is
#' allowed to change at the first break and another subset is allowed to change
#' at the second break. The function searches over ordered break-date pairs,
#' computes a likelihood-ratio statistic, and optionally simulates confidence
#' bounds for the ordered break dates.
#'
#' @usage
#' mainop(y, z, S, q1 = NULL,
#'   q2 = NULL, trim = 0.2,
#'   hetq = TRUE, max_break_gap = NULL,
#'   confidence_intervals = FALSE, ci_grid_radius = 20,
#'   ci_replications = 5000, ci_step = 1000, ci_grid_step = 0.02,
#'   ci_seed = NULL, ci_verbose = FALSE, tol = 1e-8,
#'   max_iter = 1000)
#'
#' @param y Numeric matrix, vector, or data frame of dependent variables with
#'   \eqn{T} rows and \eqn{n} columns.
#' @param z Numeric matrix, vector, or data frame of regressors with \eqn{T}
#'   rows. Include an intercept column if desired.
#' @param S Selection matrix mapping equation-stacked raw regressors into the
#'   transformed system design.
#' @param q1 Integer number of transformed coefficients
#'   that may change at the first ordered break.
#' @param q2 Integer number of transformed coefficients
#'   that may change at the second ordered break.
#' @param trim Trimming parameter. The minimum segment length is
#'   \code{floor(trim * nrow(y))}.
#' @param hetq Logical. If \code{TRUE}, allow the regressor distribution to
#'   change across regimes (Gauss's \code{hetq}); if \code{FALSE}, pooled
#'   regressor moments are used when computing ordered-break confidence intervals.
#' @param max_break_gap Optional maximum allowed distance between the first and
#'   second ordered break. If \code{NULL}, defaults to \code{floor(sqrt(T))}.
#' @param confidence_intervals Logical. If \code{TRUE}, compute simulated
#'   confidence bounds for the ordered break dates.
#' @param ci_grid_radius Integer simulation grid radius used by the confidence
#'   interval routine.
#' @param ci_replications Integer number of simulation replications used by the
#'   confidence interval routine.
#' @param ci_step Integer grid scaling parameter used by the confidence interval
#'   routine.
#' @param ci_grid_step Positive numeric grid spacing used by the confidence
#'   interval routine.
#' @param ci_seed Optional random seed for confidence interval simulation.
#' @param ci_verbose Logical. If \code{TRUE}, print confidence interval warnings
#'   and progress messages.
#' @param tol Positive convergence tolerance for ordered FGLS estimation.
#' @param max_iter Integer maximum number of ordered FGLS iterations.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{setup}}{Input settings and derived dimensions.}
#'   \item{\code{transformed}}{Stacked dependent variable and transformed design
#'     matrix.}
#'   \item{\code{testing}}{Ordered SupLR statistic and critical values.}
#'   \item{\code{estimation}}{Estimated ordered break dates, truncated break
#'     dates, coefficients, residuals, and coefficient shifts.}
#'   \item{\code{confidence_intervals}}{Confidence bounds when requested;
#'     otherwise \code{NULL}.}
#'   \item{\code{brorder}}{Raw output from the ordered break-date search.}
#' }
#'
#' @references
#' Qu, Z., and Perron, P. (2007). Estimating and testing structural changes
#' in multivariate regressions. \emph{Econometrica}, 75(2), 459-502.
#'
#' @examples
#' set.seed(2)
#' bigt <- 50
#' z <- cbind(1, rnorm(bigt), rnorm(bigt))
#' y <- matrix(rnorm(bigt), ncol = 1)
#' fit <- mainop(
#'   y = y,
#'   z = z,
#'   S = diag(ncol(z)),
#'   q1 = 1,
#'   q2 = 1,
#'   trim = 0.2,
#'   max_break_gap = 10
#' )
#' fit$estimation$breaks
#'
#' @export
mainop <- function(y, z, S,
                   q1 = NULL,
                   q2 = NULL,
                   trim = 0.2,
                   hetq = TRUE,
                   max_break_gap = NULL,
                   confidence_intervals = FALSE,
                   ci_grid_radius = 20,
                   ci_replications = 5000,
                   ci_step = 1000,
                   ci_grid_step = 0.02,
                   ci_seed = NULL,
                   ci_verbose = FALSE,
                   tol = 1e-8,
                   max_iter = 1000) {
  y <- asmat(y, "y")
  z <- asmat(z, "z")
  S <- asmat(S, "S", allow_vector = FALSE)
  valyz(y, z)

  n <- ncol(y)
  bigt <- nrow(y)
  q1 <- as.integer(q1)[[1]]
  q2 <- as.integer(q2)[[1]]
  trm <- as.numeric(trim)[[1]]
  if (is.null(max_break_gap)) {
    max_break_gap <- floor(sqrt(bigt))
  }
  trun <- as.integer(max_break_gap)[[1]]

  ci_grid_radius <- as.integer(ci_grid_radius)[[1]]
  ci_replications <- as.integer(ci_replications)[[1]]
  ci_step <- as.integer(ci_step)[[1]]
  if (!is.numeric(ci_grid_step) || length(ci_grid_step) != 1 ||
      !is.finite(ci_grid_step) || ci_grid_step <= 0) {
    stop("`ci_grid_step` must be a positive finite number.", call. = FALSE)
  }
  if (!is.logical(ci_verbose) || length(ci_verbose) != 1 || is.na(ci_verbose)) {
    stop("`ci_verbose` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  h <- floor(trm * bigt)
  transformed <- transf(y, z, S)
  ordered <- brorder(
    maty = transformed$maty,
    matz = transformed$matz,
    n = n,
    q1 = q1,
    q2 = q2,
    h = h,
    bigt = bigt,
    trun = trun,
    tol = tol,
    max_iter = max_iter
  )
  delt1 <- rbind(ordered$beta2[seq_len(q1), , drop = FALSE] - ordered$beta1[seq_len(q1), , drop = FALSE],
                 matrix(0, nrow = q2, ncol = 1))
  delt2 <- rbind(matrix(0, nrow = q1, ncol = 1),
                 ordered$beta2[(q1 + 1):(q1 + q2), , drop = FALSE] -
                   ordered$beta1[(q1 + 1):(q1 + q2), , drop = FALSE])

  bound <- NULL
  if (confidence_intervals) {
    bound <- interval3(
      matz = transformed$matz,
      maty = transformed$maty,
      n = n,
      start = 1,
      last = bigt,
      dx1 = ordered$brdate[[1]],
      dx2 = ordered$brdate[[2]],
      delt1 = delt1,
      delt2 = delt2,
      res = ordered$umat,
      segT = bigt,
      bigM = ci_grid_radius,
      rep = ci_replications,
      hetq = hetq,
      step = ci_step,
      interv = ci_grid_step,
      seed = ci_seed,
      verbose = ci_verbose
    )
  }

  list(
    setup = list(
      trim = trm,
      bigt = bigt,
      n = n,
      h = h,
      q1 = q1,
      q2 = q2,
      hetq = hetq,
      trun = trun,
      bigM = ci_grid_radius,
      rep = ci_replications,
      confidence_intervals = confidence_intervals
    ),
    transformed = transformed,
    testing = list(
      sup_lr_ordered = ordered$lr,
      critical_values = cvlr(1, q1 + q2, h / bigt)
    ),
    estimation = list(
      breaks = ordered$brdate,
      truncated_breaks = ordered$truncated_brdate,
      beta1 = ordered$beta1,
      beta2 = ordered$beta2,
      residuals = ordered$umat,
      delta1 = delt1,
      delta2 = delt2
    ),
    confidence_intervals = bound,
    brorder = ordered
  )
}

#' Ordered Two-Break convenience interface
#'
#' Formula/data-frame interface for the ordered two-break model. If \code{S} is
#' omitted, the interface uses the common-regressor identity selection matrix and
#' calls \code{\link{mainop}} internally.
#'
#' @usage
#' obreak(y = NULL, z = NULL, S = NULL,
#'   q1 = NULL, q2 = NULL,
#'   trim = 0.2, hetq = TRUE,
#'   max_break_gap = NULL, confidence_intervals = FALSE,
#'   ci_grid_radius = 20, ci_replications = 5000, ci_step = 1000,
#'   ci_grid_step = 0.02, ci_seed = NULL, ci_verbose = FALSE,
#'   data = NULL, formula = NULL, tol = 1e-8, max_iter = 1000)
#'
#' @param y Numeric response matrix/data frame, or a formula such as
#'   \code{y ~ x1 + x2} or \code{cbind(y1, y2) ~ x1 + x2}.
#' @param z Numeric regressor matrix/data frame shared by all equations.
#'   Include an intercept column if desired. Ignored when formula input is used.
#' @param S Optional selection matrix mapping equation-stacked raw regressors
#'   into the transformed system design. If \code{NULL}, the common-regressor
#'   identity selection matrix is used.
#' @param q1 Integer number of transformed coefficients
#'   that may change at the first ordered break.
#' @param q2 Integer number of transformed coefficients
#'   that may change at the second ordered break.
#' @param trim Trimming parameter. The minimum segment length is
#'   \code{floor(trim * nrow(y))}.
#' @param hetq Logical. If \code{TRUE}, allow the regressor distribution to
#'   change across regimes (Gauss's \code{hetq}); if \code{FALSE}, pooled
#'   regressor moments are used when computing ordered-break confidence intervals.
#' @param max_break_gap Optional maximum allowed distance between the first and
#'   second ordered break. If \code{NULL}, defaults to \code{floor(sqrt(T))}.
#' @param confidence_intervals Logical. If \code{TRUE}, compute simulated
#'   confidence bounds for the ordered break dates.
#' @param ci_grid_radius Integer simulation grid radius used by the confidence
#'   interval routine.
#' @param ci_replications Integer number of simulation replications used by the
#'   confidence interval routine.
#' @param ci_step Integer grid scaling parameter used by the confidence interval
#'   routine.
#' @param ci_grid_step Positive numeric grid spacing used by the confidence
#'   interval routine.
#' @param ci_seed Optional random seed for confidence interval simulation.
#' @param ci_verbose Logical. If \code{TRUE}, print confidence interval warnings
#'   and progress messages.
#' @param data Optional data frame used when formula input is supplied.
#' @param formula Optional formula alternative to supplying \code{y} as a
#'   formula.
#' @param tol Positive convergence tolerance for ordered FGLS estimation.
#' @param max_iter Integer maximum number of ordered FGLS iterations.
#'
#' @return A list returned by \code{\link{mainop}}, with spec metadata added
#'   in component \code{spec}.
#'
#' @export
obreak <- function(y = NULL, z = NULL, S = NULL,
                           q1 = NULL,
                           q2 = NULL,
                           trim = 0.2,
                           hetq = TRUE,
                           max_break_gap = NULL,
                           confidence_intervals = FALSE,
                           ci_grid_radius = 20,
                           ci_replications = 5000,
                           ci_step = 1000,
                           ci_grid_step = 0.02,
                           ci_seed = NULL,
                           ci_verbose = FALSE,
                           data = NULL,
                           formula = NULL,
                           tol = 1e-8,
                           max_iter = 1000) {
  formula_mode <- !is.null(formula) || inherits(y, "formula")
  user_supplied_S <- !is.null(S)

  if (formula_mode) {
    if (is.null(formula)) {
      formula <- y
    }
    if (is.null(data) && is.data.frame(z)) {
      data <- z
    }
    pieces <- mb_mats(formula = formula, data = data)
    y <- pieces$y
    z <- pieces$z
  } else {
    y <- asmat(y, "y")
    z <- asmat(z, "z")
  }

  valyz(y, z)
  if (is.null(colnames(y))) {
    colnames(y) <- paste0("y", seq_len(ncol(y)))
  }
  if (is.null(colnames(z))) {
    colnames(z) <- paste0("z", seq_len(ncol(z)))
  }
  if (is.null(S)) {
    S <- diag(ncol(y) * ncol(z))
  } else {
    S <- asmat(S, "S", allow_vector = FALSE)
  }

  fit <- mainop(
    y = y,
    z = z,
    S = S,
    q1 = q1,
    q2 = q2,
    trim = trim,
    hetq = hetq,
    max_break_gap = max_break_gap,
    confidence_intervals = confidence_intervals,
    ci_grid_radius = ci_grid_radius,
    ci_replications = ci_replications,
    ci_step = ci_step,
    ci_grid_step = ci_grid_step,
    ci_seed = ci_seed,
    ci_verbose = ci_verbose,
    tol = tol,
    max_iter = max_iter
  )
  fit$spec <- list(
    model_type = if (user_supplied_S) "ordered_custom" else "ordered_common",
    formula = if (formula_mode) pieces$formula else NULL,
    y_names = colnames(y),
    z_names = colnames(z),
    eqz = rep(list(colnames(z)), ncol(y)),
    complete_cases = if (formula_mode) pieces$complete_cases else seq_len(nrow(y)),
    S = S
  )
  fit
}

#' Ordered Feasible GLS Fit
#'
#' Fits the ordered-break design by feasible GLS, optionally using a previous
#' covariance matrix as the initial weight.
ofgls <- function(segx, segy, n, init_s = NULL, tol = 1e-8, max_iter = 1000) {
  if (is.null(init_s)) {
    b <- pinv(t(segx) %*% segx) %*% t(segx) %*% segy
  } else {
    ibigv <- kronecker(diag(nrow(segy) / n), pinv(init_s))
    b <- pinv(t(segx) %*% ibigv %*% segx) %*% t(segx) %*% ibigv %*% segy
  }

  res <- segy - segx %*% b
  umat <- matrix(as.vector(res), ncol = n, byrow = TRUE)
  vvar <- t(umat) %*% umat / nrow(umat)

  upd <- function(prev_sigma) {
    ibigv <- kronecker(diag(nrow(segy) / n), pinv(prev_sigma))
    b   <<- pinv(t(segx) %*% ibigv %*% segx) %*% t(segx) %*% ibigv %*% segy
    res <<- segy - segx %*% b
    umat <<- matrix(as.vector(res), ncol = n, byrow = TRUE)
    t(umat) %*% umat / nrow(umat)
  }
  fp <- fixed_point(vvar, upd, tol, max_iter, "ofgls()")

  list(beta = b, sigma = fp$value, residuals = res, umat = umat,
       iterations = fp$iterations, converged = fp$converged)
}
