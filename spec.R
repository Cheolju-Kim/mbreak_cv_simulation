#' Multivariate Break Model with Common Regressors
#'
#' Estimates an mbreak model where every equation uses the same regressor
#' matrix. This is the direct convenience interface for the case where the
#' selection matrix is the identity.
#'
#' @usage
#' mbreak(y = NULL, z = NULL, m, R = NULL, coef_breaks = NULL,
#'   data = NULL, formula = NULL, ...)
#'
#' @param y Numeric response matrix/data frame with one column per equation, or
#'   a formula such as \code{cbind(y1, y2) ~ x1 + x2}.
#' @param z Numeric regressor matrix/data frame shared by all equations. Include
#'   an intercept column if desired. Ignored when formula input is used.
#' @param m Integer number of breaks to estimate.
#' @param R Optional restriction/reparameterization matrix passed to
#'   \code{\link{mainp}}. Do not supply both \code{R} and
#'   \code{coef_breaks}.
#' @param coef_breaks Optional coefficient-break selector used to construct
#'   \code{R} automatically. Use \code{"all"} for all coefficients changing,
#'   \code{"none"} for no coefficient changes, or a list of formulas such as
#'   \code{list(y1 ~ x1, y2 ~ 0)} or \code{list(y1 ~ ., y2 ~ 0)}. These
#'   formulas are coefficient selectors: \code{y1 ~ x1} selects only
#'   \code{x1}; use \code{y1 ~ 1} or \code{y1 ~ .} to select the intercept.
#' @param data Optional data frame used when \code{y} is a formula.
#' @param formula Optional formula alternative to supplying \code{y} as a
#'   formula.
#' @param ... Additional arguments passed to \code{\link{mainp}}.
#' @return An object returned by \code{\link{mainp}}, with spec metadata.
#' @export
mbreak <- function(y = NULL, z = NULL, m, R = NULL,
                          coef_breaks = NULL, data = NULL,
                          formula = NULL, ...) {
  formula_mode <- !is.null(formula) || inherits(y, "formula")
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

  eqz <- rep(list(colnames(z)), ncol(y))
  names(eqz) <- colnames(y)
  S <- diag(ncol(y) * ncol(z))
  coef_names <- eqnames(colnames(y), eqz)
  rownames(S) <- coef_names
  colnames(S) <- coef_names
  fit_info <- call_mainp(
    y = y, z = z, m = m, S = S, R = R,
    coef_breaks = coef_breaks,
    eqz = eqz,
    ...
  )
  fit <- fit_info$fit
  fit$spec <- list(
    model_type = "common",
    formula = if (formula_mode) pieces$formula else NULL,
    y_names = colnames(y),
    z_names = colnames(z),
    eqz = eqz,
    complete_cases = if (formula_mode) pieces$complete_cases else seq_len(nrow(y)),
    coef_breaks = fit_info$coef_breaks,
    S = S
  )
  fit
}

#' Multivariate Break Model for a Seemingly Unrelated Regression System
#'
#' Estimates an mbreak model from a list of equation formulas. The interface
#' builds the response matrix, a union regressor matrix, and the selection
#' matrix that maps each equation to its own regressors.
#'
#' @usage
#' mbsur(formulas, data, m, R = NULL, coef_breaks = NULL, ...)
#'
#' @param formulas List of formulas, one per equation, such as
#'   \code{list(y1 ~ x1 + x2, y2 ~ x2 + x4)}.
#' @param data Data frame containing all variables used by the formulas.
#' @param m Integer number of breaks to estimate.
#' @param R Optional restriction/reparameterization matrix passed to
#'   \code{\link{mainp}}. Do not supply both \code{R} and
#'   \code{coef_breaks}.
#' @param coef_breaks Optional coefficient-break selector used to construct
#'   \code{R} automatically. Use \code{"all"} for all coefficients changing,
#'   \code{"none"} for no coefficient changes, or a list of formulas such as
#'   \code{list(y1 ~ x1, y2 ~ x3)} or \code{list(y1 ~ ., y2 ~ 0)}. These
#'   formulas are coefficient selectors: \code{y1 ~ x1} selects only
#'   \code{x1}; use \code{y1 ~ 1} or \code{y1 ~ .} to select the intercept.
#' @param ... Additional arguments passed to \code{\link{mainp}}.
#' @return An object returned by \code{\link{mainp}}, with spec metadata.
#' @export
mbsur <- function(formulas, data, m, R = NULL, coef_breaks = NULL, ...) {
  pieces <- sur_mats(formulas = formulas, data = data)
  fit_info <- call_mainp(
    y = pieces$y, z = pieces$z, m = m, S = pieces$S, R = R,
    coef_breaks = coef_breaks,
    eqz = pieces$eqz,
    ...
  )
  fit <- fit_info$fit
  fit$spec <- list(
    model_type = "sur",
    formulas = pieces$formulas,
    y_names = colnames(pieces$y),
    z_names = colnames(pieces$z),
    eqz = pieces$eqz,
    complete_cases = pieces$complete_cases,
    coef_breaks = fit_info$coef_breaks,
    S = pieces$S
  )
  fit
}

#' Multivariate Break Model for a VAR
#'
#' Estimates an mbreak model where every equation uses the same lagged
#' endogenous variables as regressors.
#'
#' @usage
#' mbvar(y, lags, m, intercept = TRUE, R = NULL,
#'   coef_breaks = NULL, ...)
#'
#' @param y Numeric response matrix with one column per variable.
#' @param lags Positive integer number of VAR lags.
#' @param m Integer number of breaks to estimate.
#' @param intercept Logical. If \code{TRUE}, include an intercept in the VAR
#'   regressor matrix.
#' @param R Optional restriction/reparameterization matrix passed to
#'   \code{\link{mainp}}. Do not supply both \code{R} and
#'   \code{coef_breaks}.
#' @param coef_breaks Optional coefficient-break selector used to construct
#'   \code{R} automatically after the VAR lag matrix is built. Use
#'   \code{"all"}, \code{"none"}, or a list such as
#'   \code{list(y1 ~ ., y2 ~ 0)} or
#'   \code{list(y1 ~ L1.y1 + L2.y2, y2 ~ 0)}. These formulas are coefficient
#'   selectors: use \code{y1 ~ 1} or \code{y1 ~ .} to select the intercept.
#' @param ... Additional arguments passed to \code{\link{mainp}}.
#' @return An object returned by \code{\link{mainp}}, with spec metadata.
#' @export
mbvar <- function(y, lags, m, intercept = TRUE, R = NULL,
                       coef_breaks = NULL, ...) {
  pieces <- var_mats(y = y, lags = lags, intercept = intercept)
  S <- diag(ncol(pieces$y) * ncol(pieces$z))
  eqz <- rep(list(colnames(pieces$z)), ncol(pieces$y))
  names(eqz) <- colnames(pieces$y)
  coef_names <- eqnames(colnames(pieces$y), eqz)
  rownames(S) <- coef_names
  colnames(S) <- coef_names
  fit_info <- call_mainp(
    y = pieces$y, z = pieces$z, m = m, S = S, R = R,
    coef_breaks = coef_breaks,
    eqz = eqz,
    ...
  )
  fit <- fit_info$fit
  fit$spec <- list(
    model_type = "var",
    lags = as.integer(lags),
    intercept = intercept,
    y_names = colnames(pieces$y),
    z_names = colnames(pieces$z),
    eqz = eqz,
    response_rows = pieces$response_rows,
    coef_breaks = fit_info$coef_breaks,
    S = S
  )
  fit
}

#' Call mainp with Optional Coefficient-Break Selectors
call_mainp <- function(y, z, m, S, R = NULL, coef_breaks = NULL,
                               eqz, ...) {
  dots <- list(...)
  coef_info <- NULL

  if (!is.null(coef_breaks)) {
    if (!is.null(R)) {
      stop("Do not supply both `R` and `coef_breaks`.")
    }
    coef_info <- breaks2R(
      coef_breaks = coef_breaks,
      m = m,
      y_names = colnames(y),
      eqz = eqz
    )
    R <- coef_info$R

    if (!is.null(dots$sumq)) {
      changing_user <- as.integer(dots$sumq)
      if (length(changing_user) != 1 || is.na(changing_user) ||
          changing_user != coef_info$sumq) {
        stop(
          "`sumq` conflicts with the number implied by `coef_breaks`.",
          call. = FALSE
        )
      }
    } else {
      dots$sumq <- coef_info$sumq
    }

    if (!is.null(dots$brbeta) &&
        !isTRUE(dots$brbeta) &&
        coef_info$sumq > 0) {
      stop(
        "`coef_breaks` selects changing coefficients, so `brbeta` cannot be FALSE.",
        call. = FALSE
      )
    }

    if (coef_info$sumq == 0L) {
      if (!isTRUE(dots$brv)) {
        stop(
          "`coef_breaks` selects no changing coefficients; use `brv = TRUE` for a pure covariance-break model.",
          call. = FALSE
        )
      }
      if (!is.null(dots$brbeta) && isTRUE(dots$brbeta)) {
        stop(
          "`coef_breaks` selects no changing coefficients, so `brbeta` cannot be TRUE.",
          call. = FALSE
        )
      }
      dots$brbeta <- FALSE
    }
  }

  fit <- do.call(
    mainp,
    c(list(y = y, z = z, m = m, S = S, R = R), dots)
  )
  list(fit = fit, coef_breaks = coef_info)
}

#' Build Reparameterization Matrix from Coefficient-Break Selectors
breaks2R <- function(coef_breaks, m, y_names,
                                     eqz) {
  m <- as.integer(m)
  if (length(m) != 1 || is.na(m) || m < 0) {
    stop("`m` must be a non-negative integer.")
  }

  cmap <- coefmap(y_names, eqz)
  q <- nrow(cmap)
  selected <- breaks_sel(
    coef_breaks = coef_breaks,
    cmap = cmap,
    y_names = y_names
  )
  sel_idx <- which(selected)
  nchg <- length(sel_idx)
  theta_count <- q + m * nchg

  R <- matrix(0, nrow = (m + 1) * q, ncol = theta_count)
  rownames(R) <- unlist(lapply(seq_len(m + 1), function(regime) {
    paste0("regime", regime, ":", cmap$coefficient)
  }), use.names = FALSE)
  extra_names <- character()
  if (nchg > 0L) {
    extra_names <- unlist(lapply(seq_len(m), function(break_i) {
      paste0("regime", break_i + 1, ":", cmap$coefficient[sel_idx])
    }), use.names = FALSE)
  }
  colnames(R) <- c(paste0("base:", cmap$coefficient), extra_names)

  sel_pos <- match(seq_len(q), sel_idx)
  for (regime in seq_len(m + 1)) {
    for (coef_i in seq_len(q)) {
      row_i <- (regime - 1) * q + coef_i
      if (regime == 1 || !selected[[coef_i]]) {
        R[row_i, coef_i] <- 1
      } else {
        col_i <- q + (regime - 2) * nchg + sel_pos[[coef_i]]
        R[row_i, col_i] <- 1
      }
    }
  }

  list(
    R = R,
    selector = breaks_lab(coef_breaks),
    sumq = nchg,
    selected = cmap[selected, , drop = FALSE],
    all_coefficients = cmap
  )
}

#' Select Changing Coefficients
breaks_sel <- function(coef_breaks, cmap,
                                         y_names) {
  if (is.character(coef_breaks) && length(coef_breaks) == 1L) {
    mode <- match.arg(tolower(coef_breaks), c("all", "none"))
    return(rep(mode == "all", nrow(cmap)))
  }

  if (!is.list(coef_breaks) || length(coef_breaks) < 1L) {
    stop(
      "`coef_breaks` must be \"all\", \"none\", or a list of formulas.",
      call. = FALSE
    )
  }

  selected <- rep(FALSE, nrow(cmap))
  list_names <- names(coef_breaks)
  if (is.null(list_names)) {
    list_names <- rep("", length(coef_breaks))
  }

  for (i in seq_along(coef_breaks)) {
    spec <- coef_breaks[[i]]
    if (!inherits(spec, "formula")) {
      stop("Every element of `coef_breaks` must be a formula.")
    }

    if (length(spec) == 3L) {
      equation <- rname(spec)
      rhs <- spec[[3]]
    } else if (length(spec) == 2L && nzchar(list_names[[i]])) {
      equation <- list_names[[i]]
      rhs <- spec[[2]]
    } else {
      stop(
        "Each `coef_breaks` formula must be two-sided, or named if one-sided.",
        call. = FALSE
      )
    }

    eq_i <- match(equation, y_names)
    if (is.na(eq_i)) {
      stop(
        sprintf("Equation `%s` in `coef_breaks` is not one of: %s.",
                equation, paste(y_names, collapse = ", ")),
        call. = FALSE
      )
    }

    available <- cmap$term[cmap$eq_index == eq_i]
    terms_i <- break_terms(rhs, available, equation)
    if (length(terms_i) == 0L) {
      next
    }

    idx <- which(cmap$eq_index == eq_i &
                   cmap$term %in% terms_i)
    selected[idx] <- TRUE
  }

  selected
}

#' Parse the Right-Hand Side of a Coefficient-Break Formula
break_terms <- function(rhs, available, equation) {
  rhs_text <- paste(deparse(rhs), collapse = " ")
  rhs_compact <- gsub("\\s+", "", rhs_text)

  if (rhs_compact %in% c("0", "-1")) {
    return(character())
  }
  if (rhs_compact == ".") {
    return(available)
  }

  terms_obj <- stats::terms(stats::as.formula(paste("~", rhs_text)))
  terms_i <- attr(terms_obj, "term.labels")
  if ("." %in% terms_i) {
    stop(
      "In `coef_breaks`, `.` is only supported by itself, as in `y1 ~ .`.",
      call. = FALSE
    )
  }
  rhs_pieces <- trimws(strsplit(rhs_text, "\\+", fixed = FALSE)[[1]])
  if (any(rhs_pieces == "1")) {
    terms_i <- c("(Intercept)", terms_i)
  }
  terms_i <- unique(terms_i)
  if (length(terms_i) == 0L) {
    return(character())
  }

  terms_i <- ifelse(tolower(terms_i) == "intercept", "(Intercept)", terms_i)
  missing <- setdiff(terms_i, available)
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "Coefficient(s) %s in `coef_breaks` are not available for equation `%s`. Available coefficients: %s.",
        paste(missing, collapse = ", "),
        equation,
        paste(available, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  terms_i
}

#' Build Coefficient Metadata
coefmap <- function(y_names, eqz) {
  if (length(y_names) != length(eqz)) {
    stop("`eqz` must have one element per equation.")
  }

  eq_index <- rep(seq_along(eqz), lengths(eqz))
  term <- unlist(eqz, use.names = FALSE)
  equation <- y_names[eq_index]
  data.frame(
    index = seq_along(term),
    eq_index = eq_index,
    equation = equation,
    term = term,
    coefficient = paste0(equation, ":", term),
    stringsAsFactors = FALSE
  )
}

#' Build Coefficient Names
eqnames <- function(y_names, eqz) {
  coefmap(y_names, eqz)$coefficient
}

#' Label Coefficient-Break Selectors
breaks_lab <- function(coef_breaks) {
  if (is.character(coef_breaks) && length(coef_breaks) == 1L) {
    return(coef_breaks)
  }
  if (is.list(coef_breaks)) {
    return(vapply(coef_breaks, flabel, character(1)))
  }
  NA_character_
}

#' Build SUR Matrices
sur_mats <- function(formulas, data) {
  formulas <- valfor(formulas)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame for `mbsur()`.")
  }

  neq <- length(formulas)
  frames <- vector("list", neq)
  x_list <- vector("list", neq)
  y_list <- vector("list", neq)
  y_names <- character(neq)

  for (i in seq_along(formulas)) {
    frames[[i]] <- stats::model.frame(
      formulas[[i]],
      data = data,
      na.action = stats::na.pass
    )
    terms_i <- stats::terms(formulas[[i]], data = data)
    y_i <- stats::model.response(frames[[i]])
    if (is.null(y_i) || !is.numeric(y_i)) {
      stop("Each SUR equation must have a numeric response.")
    }

    x_i <- stats::model.matrix(stats::delete.response(terms_i), data = frames[[i]])
    y_list[[i]] <- as.numeric(y_i)
    x_list[[i]] <- x_i
    y_names[[i]] <- rname(formulas[[i]])
  }

  keep <- Reduce(
    `&`,
    c(
      lapply(y_list, function(x) stats::complete.cases(x)),
      lapply(x_list, stats::complete.cases)
    )
  )
  if (!any(keep)) {
    stop("No complete observations remain after evaluating the SUR formulas.")
  }

  y <- do.call(cbind, lapply(y_list, function(x) x[keep]))
  colnames(y) <- make.unique(y_names)
  x_list <- lapply(x_list, function(x) x[keep, , drop = FALSE])

  z <- unionz(x_list)
  S <- surS(x_list, union_names = colnames(z))
  eqz <- lapply(x_list, colnames)
  names(eqz) <- colnames(y)

  list(
    y = y,
    z = z,
    S = S,
    formulas = vapply(formulas, flabel, character(1)),
    eqz = eqz,
    complete_cases = which(keep)
  )
}

#' Build VAR Matrices
var_mats <- function(y, lags, intercept = TRUE) {
  y <- asmat(y, "y")
  lags <- as.integer(lags)
  if (length(lags) != 1 || is.na(lags) || lags < 1) {
    stop("`lags` must be a positive integer.")
  }
  if (nrow(y) <= lags) {
    stop("`lags` must be smaller than the number of observations in `y`.")
  }
  if (is.null(colnames(y))) {
    colnames(y) <- paste0("y", seq_len(ncol(y)))
  }

  response_rows <- (lags + 1):nrow(y)
  y_response <- y[response_rows, , drop = FALSE]
  z_parts <- list()
  z_names <- character()

  if (intercept) {
    z_parts[[length(z_parts) + 1]] <- matrix(1, nrow = length(response_rows), ncol = 1)
    z_names <- c(z_names, "(Intercept)")
  }

  for (lag in seq_len(lags)) {
    lagged <- y[(lags + 1 - lag):(nrow(y) - lag), , drop = FALSE]
    z_parts[[length(z_parts) + 1]] <- lagged
    z_names <- c(z_names, paste0("L", lag, ".", colnames(y)))
  }

  z <- do.call(cbind, z_parts)
  colnames(z) <- z_names
  list(y = y_response, z = z, response_rows = response_rows)
}

#' Build Common-Regressor Matrices from a Formula
mb_mats <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula.")
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame when using formula input.")
  }

  frame <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  terms_i <- stats::terms(formula, data = data)
  y <- stats::model.response(frame)
  y <- asmat(y, "formula response")
  if (is.null(colnames(y))) {
    colnames(y) <- rnames(formula, ncol(y))
  }

  z <- stats::model.matrix(stats::delete.response(terms_i), data = frame)
  keep <- stats::complete.cases(y) & stats::complete.cases(z)
  if (!any(keep)) {
    stop("No complete observations remain after evaluating the formula.")
  }

  list(
    y = y[keep, , drop = FALSE],
    z = z[keep, , drop = FALSE],
    formula = flabel(formula),
    complete_cases = which(keep)
  )
}

#' Validate Response and Regressor Matrices
valyz <- function(y, z) {
  if (!is.numeric(y) || !is.numeric(z)) {
    stop("`y` and `z` must be numeric.")
  }
  if (nrow(y) != nrow(z)) {
    stop("`y` and `z` must have the same number of rows.")
  }
  invisible(NULL)
}

#' Validate SUR Formulas
valfor <- function(formulas) {
  if (!is.list(formulas)) {
    stop("`formulas` must be a list of formulas.")
  }
  if (length(formulas) < 1) {
    stop("`formulas` must contain at least one equation.")
  }
  if (!all(vapply(formulas, inherits, logical(1), what = "formula"))) {
    stop("Every element of `formulas` must be a formula.")
  }
  if (!all(vapply(formulas, length, integer(1)) == 3L)) {
    stop("Every SUR formula must include a response.")
  }
  formulas
}

#' Build Union Regressor Matrix
unionz <- function(x_list) {
  union_names <- unique(unlist(lapply(x_list, colnames), use.names = FALSE))
  z <- matrix(0, nrow = nrow(x_list[[1]]), ncol = length(union_names))
  colnames(z) <- union_names

  for (name in union_names) {
    values <- NULL
    for (x in x_list) {
      idx <- match(name, colnames(x))
      if (!is.na(idx)) {
        if (is.null(values)) {
          values <- x[, idx]
        } else if (!isTRUE(all.equal(values, x[, idx], check.attributes = FALSE))) {
          stop(
            sprintf("Regressor column `%s` has conflicting values across equations.", name),
            call. = FALSE
          )
        }
      }
    }
    z[, name] <- values
  }

  z
}

#' Build SUR Selection Matrix
surS <- function(x_list, union_names) {
  n <- length(x_list)
  p <- length(union_names)
  q_total <- sum(vapply(x_list, ncol, integer(1)))
  S <- matrix(0, nrow = n * p, ncol = q_total)
  rownames(S) <- paste0(
    "eq", rep(seq_len(n), each = p),
    ":", rep(union_names, times = n)
  )

  sel_names <- character(q_total)
  col_offset <- 0
  for (eq in seq_len(n)) {
    x_names <- colnames(x_list[[eq]])
    for (j in seq_along(x_names)) {
      raw_col <- (eq - 1) * p + match(x_names[[j]], union_names)
      selected_col <- col_offset + j
      S[raw_col, selected_col] <- 1
      sel_names[[selected_col]] <- paste0("eq", eq, ":", x_names[[j]])
    }
    col_offset <- col_offset + length(x_names)
  }
  colnames(S) <- sel_names
  S
}

#' Formula Label
flabel <- function(formula) {
  paste(deparse(formula), collapse = " ")
}

#' Formula Response Name
rname <- function(formula) {
  paste(deparse(formula[[2]]), collapse = " ")
}

#' Formula Response Names
rnames <- function(formula, n) {
  response <- formula[[2]]
  if (is.call(response) && identical(response[[1]], as.name("cbind"))) {
    names <- vapply(as.list(response[-1]), function(x) {
      paste(deparse(x), collapse = " ")
    }, character(1))
    return(make.unique(names))
  }
  if (n == 1L) {
    return(rname(formula))
  }
  paste0(rname(formula), seq_len(n))
}
