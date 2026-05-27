# Transform multivariate regression data into the stacked system used by mbreak.

#' Convert User Input to a Numeric Matrix
#'
#' Converts numeric vectors, matrices, and data frames to numeric matrices while
#' preserving column names where possible.
asmat <- function(x, arg, allow_vector = TRUE) {
  if (is.null(x)) {
    stop(sprintf("`%s` must not be NULL.", arg))
  }

  if (is.data.frame(x)) {
    numeric_cols <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      bad <- names(x)[!numeric_cols]
      stop(
        sprintf("`%s` must contain only numeric columns; nonnumeric column(s): %s.",
                arg, paste(bad, collapse = ", ")),
        call. = FALSE
      )
    }
    out <- as.matrix(x)
  } else if (is.matrix(x)) {
    out <- x
  } else if (allow_vector && is.numeric(x) && is.atomic(x)) {
    out <- matrix(x, ncol = 1)
    if (!is.null(names(x))) {
      rownames(out) <- names(x)
    }
  } else {
    stop(sprintf("`%s` must be a numeric matrix, vector, or data frame.", arg))
  }

  if (!is.numeric(out)) {
    stop(sprintf("`%s` must be numeric.", arg))
  }
  if (nrow(out) < 1 || ncol(out) < 1) {
    stop(sprintf("`%s` must have at least one row and one column.", arg))
  }

  out
}

#' Stack a Multivariate Regression System
#'
#' Converts multivariate responses and regressors into the equation-stacked
#' system used by the mbreak likelihood and test routines.
transf <- function(y, z, S) {
  bigt <- nrow(y)
  n <- ncol(y)
  maty <- matrix(as.vector(t(y)), ncol = 1)
  matz <- matrix(0, nrow = bigt * n, ncol = ncol(S))
  eye_n <- diag(n)

  for (i in seq_len(bigt)) {
    rows_i <- ((i - 1) * n + 1):(i * n)
    matz[rows_i, ] <- kronecker(eye_n, z[i, , drop = FALSE]) %*% S
  }

  list(maty = maty, matz = matz)
}

#' Regime-Block Design Matrix
#'
#' Expands a stacked design matrix into regime-specific blocks for a given set
#' of break dates.
pzbar <- function(matz, m, br, bigt = NULL) {
  if (m == 0) {
    return(matz)
  }

  n <- nrow(matz) / bigt
  q <- ncol(matz)
  endpoints <- c(0, br, bigt)

  matzb <- matrix(0, nrow = nrow(matz), ncol = (m + 1) * q)
  for (seg in seq_len(m + 1)) {
    t_start <- endpoints[[seg]] + 1
    t_end <- endpoints[[seg + 1]]
    row_idx <- stacked_rows(t_start, t_end, n)
    col_idx <- ((seg - 1) * q + 1):(seg * q)
    matzb[row_idx, col_idx] <- matz[row_idx, , drop = FALSE]
  }

  matzb
}
