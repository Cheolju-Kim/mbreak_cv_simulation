#' Bootstrap Restarting for Restricted Break Estimation
#'
#' Applies the observation-resampling or residual-bootstrap restarting
#' algorithms used by rbreak to the restricted mbreak estimator.
brbp <- function(maty, matz, n, m, bigt, trm, R,
                            method = c("brbp", "brbp_residual"),
                            initial_breaks = NULL, B = 20,
                            brbeta = TRUE,
                            brv = FALSE,
                            tol = 1e-6, max_iter = 1000,
                            max_outer = 1000, verbose = FALSE) {
  method <- match.arg(method)
  B <- as.integer(B)
  if (length(B) != 1 || is.na(B) || B < 0) {
    stop("`restart_B` must be a non-negative integer.")
  }

  h <- as.integer(round(trm * bigt))
  q <- ncol(matz)
  R <- valR(R, m = m, q = q)
  current_fit <- est(
    maty, matz, n, m, bigt, trm, R,
    brbeta = brbeta,
    brv = brv,
    tol = tol,
    max_iter = max_iter,
    max_outer = max_outer,
    initial_breaks = initial_breaks
  )
  cur_obj <- bobj(
    maty = maty,
    matz = matz,
    beta = current_fit$r_estim$nbeta,
    br = current_fit$dx,
    n = n,
    m = m,
    bigt = bigt
  )

  break_history <- matrix(NA_integer_, nrow = B + 1, ncol = m)
  objective_history <- rep(NA_real_, B + 1)
  candidate_breaks <- matrix(NA_integer_, nrow = B, ncol = m)
  candidate_objectives <- rep(NA_real_, B)
  improved <- rep(FALSE, B)
  skipped <- rep(FALSE, B)
  skip_reasons <- rep(NA_character_, B)

  break_history[1, ] <- current_fit$dx
  objective_history[1] <- cur_obj

  if (verbose) {
    cat(sprintf("%s initial objective: %.6f\n", toupper(method), cur_obj))
    cat(sprintf("%s initial breaks: %s\n", toupper(method), paste(current_fit$dx, collapse = ", ")))
  }

  for (b in seq_len(B)) {
    if (method == "brbp") {
      boot_indices <- sort(sample.int(bigt, size = bigt, replace = TRUE))
      boot_rows <- orows(boot_indices, n)
      boot_maty <- maty[boot_rows, , drop = FALSE]
      boot_matz <- matz[boot_rows, , drop = FALSE]

      boot_init <- mapbs(
        current_fit$dx,
        orig_idx = boot_indices,
        bigt = bigt,
        h = h
      )
      if (!okbrk(boot_init, m = m, bigt = bigt, h = h)) {
        skipped[b] <- TRUE
        skip_reasons[b] <- "could not map current breaks to bootstrap sample"
        break_history[b + 1, ] <- current_fit$dx
        objective_history[b + 1] <- cur_obj
        next
      }
    } else {
      boot_indices <- NULL
      boot_matz <- matz
      fitted <- pzbar(matz, m, current_fit$dx, bigt) %*% current_fit$r_estim$nbeta
      residual <- maty - fitted
      umat <- matrix(as.vector(residual), ncol = n, byrow = TRUE)
      boot_res <- umat[sample.int(bigt, size = bigt, replace = TRUE), , drop = FALSE]
      boot_maty <- fitted + matrix(as.vector(t(boot_res)), ncol = 1)
      boot_init <- current_fit$dx
    }

    boot_fit <- tryCatch(
      est(
        boot_maty, boot_matz, n, m, bigt, trm, R,
        brbeta = brbeta,
        brv = brv,
        tol = tol,
        max_iter = max_iter,
        max_outer = max_outer,
        initial_breaks = boot_init
      ),
      error = function(e) e
    )
    if (inherits(boot_fit, "error")) {
      skipped[b] <- TRUE
      skip_reasons[b] <- conditionMessage(boot_fit)
      break_history[b + 1, ] <- current_fit$dx
      objective_history[b + 1] <- cur_obj
      next
    }

    if (method == "brbp") {
      reest_init <- mapbo(
        boot_fit$dx,
        orig_idx = boot_indices,
        bigt = bigt,
        h = h
      )
      if (!okbrk(reest_init, m = m, bigt = bigt, h = h)) {
        skipped[b] <- TRUE
        skip_reasons[b] <- "could not map bootstrap breaks to original sample"
        break_history[b + 1, ] <- current_fit$dx
        objective_history[b + 1] <- cur_obj
        next
      }
    } else {
      reest_init <- boot_fit$dx
    }

    reest_fit <- tryCatch(
      est(
        maty, matz, n, m, bigt, trm, R,
        brbeta = brbeta,
        brv = brv,
        tol = tol,
        max_iter = max_iter,
        max_outer = max_outer,
        initial_breaks = reest_init
      ),
      error = function(e) e
    )
    if (inherits(reest_fit, "error")) {
      skipped[b] <- TRUE
      skip_reasons[b] <- conditionMessage(reest_fit)
      break_history[b + 1, ] <- current_fit$dx
      objective_history[b + 1] <- cur_obj
      next
    }

    reest_obj <- bobj(
      maty = maty,
      matz = matz,
      beta = reest_fit$r_estim$nbeta,
      br = reest_fit$dx,
      n = n,
      m = m,
      bigt = bigt
    )
    candidate_breaks[b, ] <- reest_fit$dx
    candidate_objectives[b] <- reest_obj

    if (is.finite(reest_obj) && reest_obj < cur_obj) {
      current_fit <- reest_fit
      cur_obj <- reest_obj
      improved[b] <- TRUE
      if (verbose) {
        cat(sprintf(
          "%s iteration %d improved objective to %.6f with breaks: %s\n",
          toupper(method),
          b, cur_obj, paste(current_fit$dx, collapse = ", ")
        ))
      }
    }

    break_history[b + 1, ] <- current_fit$dx
    objective_history[b + 1] <- cur_obj
  }

  diagnostics <- list(
    method = method,
    B = B,
    final_objective = cur_obj,
    break_history = break_history,
    objective_history = objective_history,
    candidate_breaks = candidate_breaks,
    candidate_objectives = candidate_objectives,
    improved = improved,
    skipped = skipped,
    skip_reasons = skip_reasons
  )

  list(fit = current_fit, diagnostics = diagnostics)
}

#' Restricted Objective for a Fixed Restricted Fit
#'
#' Computes the Gaussian negative log-likelihood objective minimized by
#' `dating_M2()` for fixed break dates and restricted coefficients.
bobj <- function(maty, matz, beta, br, n, m, bigt) {
  q <- ncol(matz)
  bigvec2 <- rresid(maty, matz, beta, q, n, m)
  endpoints <- c(0, br, bigt)
  objective <- 0

  for (seg in seq_len(m + 1)) {
    t_start <- endpoints[[seg]] + 1
    t_end <- endpoints[[seg + 1]]
    seg_len <- t_end - t_start + 1
    rows <- ((seg - 1) * bigt + t_start):((seg - 1) * bigt + t_end)
    residuals <- bigvec2[rows, , drop = FALSE]
    sigma <- crossprod(residuals) / seg_len
    objective <- objective +
      (seg_len * n / 2) * (log(2 * pi) + 1) +
      (seg_len / 2) * slogdet(sigma)
  }

  objective
}

#' Observation Rows for a Bootstrap Sample
#'
#' Converts sampled observation indices to stacked row indices.
orows <- function(indices, n) {
  rows <- integer(length(indices) * n)
  for (i in seq_along(indices)) {
    rows[((i - 1) * n + 1):(i * n)] <- stacked_rows(indices[[i]], indices[[i]], n)
  }
  rows
}

#' Map Break Dates to a Bootstrap Sample
#'
#' Maps original-sample break dates to chronologically sorted bootstrap
#' positions.
mapbs <- function(br, orig_idx, bigt, h) {
  mapped <- integer(length(br))
  for (i in seq_along(br)) {
    mapped[[i]] <- sum(orig_idx <= br[[i]])
  }
  adjbrk(mapped, m = length(br), bigt = length(orig_idx), h = h)
}

#' Map Bootstrap Break Dates to Original Time
#'
#' Maps bootstrap-sample break dates back to original-sample observation
#' indices.
mapbo <- function(br, orig_idx, bigt, h) {
  positions <- pmin(pmax(as.integer(round(br)), 1L), length(orig_idx))
  mapped <- orig_idx[positions]
  adjbrk(mapped, m = length(br), bigt = bigt, h = h)
}

#' Validate Restricted Break Dates
valbrk <- function(br, m, bigt, h, arg = "initial_breaks") {
  if (!okbrk(br, m = m, bigt = bigt, h = h)) {
    stop(
      sprintf(
        "`%s` must contain %d finite break date(s) satisfying the trimming and spacing constraints.",
        arg, m
      ),
      call. = FALSE
    )
  }
  sort(as.integer(round(br)))
}

#' Check Restricted Break-Date Validity
okbrk <- function(br, m, bigt, h) {
  if (m == 0) {
    return(length(br) == 0)
  }
  if (length(br) != m || anyNA(br) || !all(is.finite(br))) {
    return(FALSE)
  }

  br <- sort(as.integer(round(br)))
  if (any(br < h) || any(br > bigt - h)) {
    return(FALSE)
  }
  if (m > 1 && any(diff(br) < h)) {
    return(FALSE)
  }

  TRUE
}

#' Adjust Break Dates to Trimming Constraints
#'
#' Projects mapped break dates back into the admissible trimmed region.
adjbrk <- function(br, m, bigt, h) {
  if (m == 0) {
    return(integer())
  }
  if (length(br) != m || anyNA(br) || !all(is.finite(br))) {
    return(integer())
  }
  if ((m + 1) * h > bigt) {
    return(integer())
  }

  br <- sort(as.integer(round(br)))
  for (i in seq_len(m)) {
    lower <- i * h
    upper <- bigt - (m - i + 1) * h
    br[[i]] <- min(max(br[[i]], lower), upper)
  }
  if (m > 1) {
    for (i in 2:m) {
      br[[i]] <- max(br[[i]], br[[i - 1]] + h)
    }
    for (i in (m - 1):1) {
      br[[i]] <- min(br[[i]], br[[i + 1]] - h)
    }
  }

  if (!okbrk(br, m = m, bigt = bigt, h = h)) {
    return(integer())
  }
  br
}
