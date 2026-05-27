#' Symmetric Matrix Square Root
#'
#' Computes the symmetric square root of a positive semidefinite matrix after
#' truncating small negative eigenvalues caused by numerical error.
sqrm <- function(M, tol = sqrt(.Machine$double.eps)) {
  M <- (M + t(M)) / 2

  eig <- eigen(M, symmetric = TRUE)
  values <- eig$values

  if (any(values < -tol)) {
    warning("`M` has negative eigenvalues; setting them to zero in `sqrm()`.")
  }
  values[values < 0] <- 0

  eig$vectors %*% diag(sqrt(values), nrow = length(values)) %*% t(eig$vectors)
}

#' Column Vectorization
#'
#' Converts an object to a single-column matrix using R's column-major order.
vec <- function(x) {
  matrix(as.vector(x), ncol = 1)
}

#' Moore-Penrose Pseudoinverse
#'
#' Computes a generalized inverse using `MASS::ginv()` when available, with an
#' SVD fallback otherwise.
pinv <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (requireNamespace("MASS", quietly = TRUE)) {
    return(MASS::ginv(x, tol = tol))
  }

  sv <- svd(x)
  if (length(sv$d) == 0) {
    return(matrix(0, nrow = ncol(x), ncol = nrow(x)))
  }

  keep <- sv$d > tol * max(sv$d)
  if (!any(keep)) {
    return(matrix(0, nrow = ncol(x), ncol = nrow(x)))
  }

  sv$v[, keep, drop = FALSE] %*%
    (diag(1 / sv$d[keep], nrow = sum(keep)) %*% t(sv$u[, keep, drop = FALSE]))
}

#' Log Determinant
#'
#' Computes the log determinant of a positive definite matrix and stops if the
#' determinant is non-positive.
logdet <- function(x) {
  det_x <- determinant(x, logarithm = TRUE)

  if (det_x$sign <= 0) {
    stop("A covariance matrix has a non-positive determinant.")
  }

  as.numeric(det_x$modulus)
}

#' Safe Log Determinant
#'
#' Computes the log determinant and returns `Inf` when the determinant is
#' non-positive, optionally symmetrizing the input first.
slogdet <- function(x, symmetrize = FALSE) {
  if (symmetrize) {
    x <- (x + t(x)) / 2
  }

  det_x <- determinant(x, logarithm = TRUE)
  if (det_x$sign <= 0) {
    return(Inf)
  }
  as.numeric(det_x$modulus)
}

#' Iterate a Function to a Fixed Point
#'
#' Repeatedly applies `update()` until successive iterates differ by at most
#' `tol` in max-abs norm, or `max_iter` iterations have run.  Sibling state
#' (covariance matrices, weighting matrices, etc.) is normally captured in
#' the enclosing scope via `<<-` inside `update`.
fixed_point <- function(x, update, tol, max_iter, label) {
  for (itr in seq_len(max_iter)) {
    prev <- x
    x <- update(prev)
    if (max(abs(x - prev)) <= tol) {
      return(list(value = x, iterations = itr, converged = TRUE))
    }
  }
  warnconv(label, max_iter, tol)
  list(value = x, iterations = max_iter, converged = FALSE)
}

#' Warn About Iteration Non-Convergence
warnconv <- function(routine, max_iter, tol) {
  warning(
    sprintf(
      "%s did not converge within max_iter = %d (tol = %g).",
      routine, as.integer(max_iter), tol
    ),
    call. = FALSE
  )
  invisible(NULL)
}
