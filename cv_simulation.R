# Critical-value simulation and response-surface fitting for Bai-Perron tests.
#
# This file is intentionally separate from tests.R.  The functions here are
# offline tools: simulate critical values on a grid, fit response-surface
# coefficients, and then use those coefficients to update cvlr(), cvud(),
# cvwd(), and supseq() if desired.

bp_cv_probs <- c("90" = 0.90, "95" = 0.95, "97.5" = 0.975, "99" = 0.99)

bp_feasible_max_breaks <- function(trm, n_grid = 1000L, cap = Inf) {
  if (!is.numeric(trm) || length(trm) != 1L || !is.finite(trm) || trm <= 0 || trm >= 0.5) {
    stop("`trm` must be one numeric value in (0, 0.5).")
  }
  h <- ceiling(trm * n_grid)
  max_k <- floor((n_grid - 1L) / h) - 1L
  cap_value <- if (is.infinite(cap)) max_k else as.integer(cap)
  as.integer(max(1L, min(cap_value, max_k)))
}

bp_default_doublemax_m <- function(trm, n_grid = 1000L, cap = 5L) {
  vapply(
    trm,
    function(eps) as.integer(min(cap, bp_feasible_max_breaks(eps, n_grid = n_grid, cap = Inf))),
    integer(1L)
  )
}

bp_response_grid <- function(q = 1:10,
                             trm = c(0.05, 0.10, 0.15, 0.20, 0.25),
                             k_cap = 9L,
                             l_cap = 9L,
                             double_m_cap = 5L,
                             n_grid = 1000L) {
  fixed <- do.call(rbind, lapply(trm, function(eps) {
    max_k <- bp_feasible_max_breaks(eps, n_grid = n_grid, cap = k_cap)
    expand.grid(
      q = q,
      trm = eps,
      k = seq_len(max_k),
      KEEP.OUT.ATTRS = FALSE
    )
  }))

  double <- expand.grid(
    q = q,
    trm = trm,
    KEEP.OUT.ATTRS = FALSE
  )
  double$M <- bp_default_doublemax_m(double$trm, n_grid = n_grid, cap = double_m_cap)

  sequential <- do.call(rbind, lapply(trm, function(eps) {
    expand.grid(
      q = q,
      trm = eps,
      l = 0:as.integer(l_cap),
      KEEP.OUT.ATTRS = FALSE
    )
  }))

  list(fixed = fixed, double = double, sequential = sequential)
}

bp_extended_response_grid <- function(q = sort(unique(c(1:20, seq(25, 60, by = 5), 70, 80, 90, 100))),
                                      trm = c(
                                        0.05, 0.075, 0.10, 0.125, 0.15,
                                        0.175, 0.20, 0.225, 0.25, 0.30
                                      ),
                                      k_cap = 10L,
                                      l_cap = 10L,
                                      double_m_cap = 5L,
                                      n_grid = 1000L) {
  bp_response_grid(
    q = q,
    trm = trm,
    k_cap = k_cap,
    l_cap = l_cap,
    double_m_cap = double_m_cap,
    n_grid = n_grid
  )
}

bp_quantiles <- function(x, probs = bp_cv_probs) {
  as.numeric(stats::quantile(x, probs = probs, names = FALSE, type = 7))
}

bp_bridge <- function(n_grid, q) {
  increments <- matrix(stats::rnorm(n_grid * q), nrow = n_grid, ncol = q) / sqrt(n_grid)
  wiener <- rbind(rep(0, q), apply(increments, 2L, cumsum))
  time <- (0:n_grid) / n_grid
  wiener - tcrossprod(time, wiener[n_grid + 1L, ])
}

bp_segment_score_matrix <- function(bridge, h) {
  n_grid <- nrow(bridge) - 1L
  score <- matrix(-Inf, nrow = n_grid + 1L, ncol = n_grid + 1L)

  for (start in 0:(n_grid - h)) {
    ends <- (start + h):n_grid
    diff <- sweep(bridge[ends + 1L, , drop = FALSE], 2L, bridge[start + 1L, ], "-")
    score[start + 1L, ends + 1L] <- rowSums(diff * diff) / ((ends - start) / n_grid)
  }

  score
}

bp_partition_stats <- function(q, trm, k_max, n_grid = 1000L,
                               scale = c("mbreak_lr", "bp_table")) {
  scale <- match.arg(scale)
  h <- ceiling(trm * n_grid)
  max_allowed <- bp_feasible_max_breaks(trm, n_grid = n_grid, cap = Inf)
  if (k_max > max_allowed) {
    stop("`k_max` is not feasible for this trimming parameter and grid size.")
  }

  bridge <- bp_bridge(n_grid, q)
  score <- bp_segment_score_matrix(bridge, h)
  stats <- numeric(k_max)

  dp_prev <- score[1L, ]
  for (segments in 2:(k_max + 1L)) {
    dp_curr <- rep(-Inf, n_grid + 1L)
    min_end <- segments * h
    if (min_end <= n_grid) {
      for (end in min_end:n_grid) {
        breaks <- ((segments - 1L) * h):(end - h)
        dp_curr[end + 1L] <- max(dp_prev[breaks + 1L] + score[breaks + 1L, end + 1L])
      }
    }
    stats[segments - 1L] <- dp_curr[n_grid + 1L]
    dp_prev <- dp_curr
  }

  if (scale == "bp_table") {
    stats <- stats / seq_len(k_max)
  }

  stats
}

bp_one_break_stat <- function(q, trm, n_grid = 1000L) {
  h <- ceiling(trm * n_grid)
  bridge <- bp_bridge(n_grid, q)
  idx <- h:(n_grid - h)
  mu <- idx / n_grid
  max(rowSums(bridge[idx + 1L, , drop = FALSE]^2) / (mu * (1 - mu)))
}

bp_progress <- function(verbose, label, i, n, start_time) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  step <- max(1L, floor(n / 10L))
  if (i == 1L || i == n || i %% step == 0L) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    message(sprintf("%s: %d/%d (%.1f%%), elapsed %.1fs", label, i, n, 100 * i / n, elapsed))
  }
  invisible(NULL)
}

bp_parallel_lapply <- function(X, FUN, n_cores = 1L, seed = NULL, export = character()) {
  n_cores <- max(1L, min(as.integer(n_cores), length(X)))
  if (n_cores <= 1L) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    return(lapply(X, FUN))
  }

  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The `parallel` package is required for `n_cores > 1`.")
  }

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, seed)
  }
  if (length(export)) {
    parallel::clusterExport(cl, export, envir = environment(bp_parallel_lapply))
  }

  parallel::parLapplyLB(cl, X, FUN)
}

bp_chunks <- function(n, n_cores = 1L, chunk_factor = 4L) {
  n_chunks <- max(1L, min(n, as.integer(n_cores) * as.integer(chunk_factor)))
  split(seq_len(n), rep(seq_len(n_chunks), length.out = n))
}

bp_simulate_partition_matrix <- function(q, trm, k_max, rep, n_grid = 1000L,
                                         n_cores = 1L, seed = NULL,
                                         chunk_factor = 4L) {
  if (n_cores <= 1L) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    mat <- matrix(NA_real_, nrow = rep, ncol = k_max)
    for (i in seq_len(rep)) {
      mat[i, ] <- bp_partition_stats(q, trm, k_max, n_grid = n_grid, scale = "bp_table")
    }
    return(mat)
  }

  chunks <- bp_chunks(rep, n_cores = n_cores, chunk_factor = chunk_factor)
  chunk_fun <- local({
    q0 <- q
    trm0 <- trm
    k_max0 <- k_max
    n_grid0 <- n_grid
    function(idx) {
      mat <- matrix(NA_real_, nrow = length(idx), ncol = k_max0)
      for (i in seq_along(idx)) {
        mat[i, ] <- bp_partition_stats(q0, trm0, k_max0, n_grid = n_grid0, scale = "bp_table")
      }
      mat
    }
  })

  draws <- bp_parallel_lapply(
    X = chunks,
    FUN = chunk_fun,
    n_cores = n_cores,
    seed = seed,
    export = c(
      "bp_chunks",
      "bp_partition_stats",
      "bp_feasible_max_breaks",
      "bp_bridge",
      "bp_segment_score_matrix"
    )
  )

  do.call(rbind, draws)
}

bp_simulate_one_break_draws <- function(q, trm, n_draws, n_grid = 1000L,
                                        n_cores = 1L, seed = NULL,
                                        chunk_factor = 4L) {
  if (n_cores <= 1L) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    v <- numeric(n_draws)
    for (i in seq_len(n_draws)) {
      v[i] <- bp_one_break_stat(q, trm, n_grid = n_grid)
    }
    return(v)
  }

  chunks <- bp_chunks(n_draws, n_cores = n_cores, chunk_factor = chunk_factor)
  chunk_fun <- local({
    q0 <- q
    trm0 <- trm
    n_grid0 <- n_grid
    function(idx) {
      v <- numeric(length(idx))
      for (i in seq_along(idx)) {
        v[i] <- bp_one_break_stat(q0, trm0, n_grid = n_grid0)
      }
      v
    }
  })

  draws <- bp_parallel_lapply(
    X = chunks,
    FUN = chunk_fun,
    n_cores = n_cores,
    seed = seed,
    export = c("bp_chunks", "bp_one_break_stat", "bp_bridge")
  )

  unlist(draws, use.names = FALSE)
}

bp_cv_rows <- function(stat, q, trm, values, probs, k = NA_integer_, l = NA_integer_, M = NA_integer_,
                       rep, n_grid, scale = NA_character_) {
  data.frame(
    stat = stat,
    q = q,
    trm = trm,
    k = k,
    l = l,
    M = M,
    prob = as.numeric(probs),
    prob_label = names(probs),
    cv = as.numeric(values),
    rep = rep,
    n_grid = n_grid,
    scale = scale,
    stringsAsFactors = FALSE
  )
}

simulate_bp_response_surface_data <- function(grid = bp_response_grid(),
                                              rep = 10000L,
                                              n_grid = 1000L,
                                              probs = bp_cv_probs,
                                              seed = NULL,
                                              fixed_scale = c("mbreak_lr", "bp_table"),
                                              n_cores = 1L,
                                              verbose = TRUE) {
  fixed_scale <- match.arg(fixed_scale)
  n_cores <- max(1L, as.integer(n_cores))
  if (!is.null(seed) && n_cores <= 1L) {
    set.seed(seed)
  }

  out <- list()
  out_i <- 1L

  fixed <- grid$fixed
  double <- grid$double
  sequential <- grid$sequential

  fixed_group_parts <- list()
  if (NROW(fixed)) {
    fixed_group_parts[[length(fixed_group_parts) + 1L]] <- fixed[, c("q", "trm"), drop = FALSE]
  }
  if (NROW(double)) {
    fixed_group_parts[[length(fixed_group_parts) + 1L]] <- double[, c("q", "trm"), drop = FALSE]
  }
  fixed_groups <- if (length(fixed_group_parts)) {
    unique(do.call(rbind, fixed_group_parts))
  } else {
    data.frame(q = numeric(), trm = numeric())
  }

  for (g in seq_len(NROW(fixed_groups))) {
    qg <- fixed_groups$q[g]
    trmg <- fixed_groups$trm[g]
    fixed_sub <- fixed[fixed$q == qg & fixed$trm == trmg, , drop = FALSE]
    double_sub <- double[double$q == qg & double$trm == trmg, , drop = FALSE]

    max_k <- 0L
    if (NROW(fixed_sub)) {
      max_k <- max(max_k, max(fixed_sub$k))
    }
    if (NROW(double_sub)) {
      max_k <- max(max_k, max(double_sub$M))
    }
    if (max_k == 0L) {
      next
    }

    label <- sprintf("partition q=%s trm=%s", qg, trmg)
    start_time <- Sys.time()
    if (isTRUE(verbose)) {
      message(sprintf("%s: starting %d replications on %d core(s)", label, rep, n_cores))
    }
    table_stats <- bp_simulate_partition_matrix(
      q = qg,
      trm = trmg,
      k_max = max_k,
      rep = rep,
      n_grid = n_grid,
      n_cores = n_cores,
      seed = if (is.null(seed)) NULL else as.integer(seed + g)
    )
    if (isTRUE(verbose)) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      message(sprintf("%s: finished in %.1fs", label, elapsed))
    }
    fixed_stats <- table_stats
    if (fixed_scale == "mbreak_lr") {
      fixed_stats <- t(t(table_stats) * seq_len(max_k))
    }

    if (NROW(fixed_sub)) {
      for (kg in sort(unique(fixed_sub$k))) {
        out[[out_i]] <- bp_cv_rows(
          stat = "cvlr",
          q = qg,
          trm = trmg,
          k = kg,
          values = bp_quantiles(fixed_stats[, kg], probs),
          probs = probs,
          rep = rep,
          n_grid = n_grid,
          scale = fixed_scale
        )
        out_i <- out_i + 1L
      }
    }

    if (NROW(double_sub)) {
      for (Mg in sort(unique(double_sub$M))) {
        sub_stats <- table_stats[, seq_len(Mg), drop = FALSE]
        ud_stats <- apply(sub_stats, 1L, max)
        out[[out_i]] <- bp_cv_rows(
          stat = "cvud",
          q = qg,
          trm = trmg,
          M = Mg,
          values = bp_quantiles(ud_stats, probs),
          probs = probs,
          rep = rep,
          n_grid = n_grid,
          scale = "bp_table"
        )
        out_i <- out_i + 1L

        fixed_cv <- apply(sub_stats, 2L, bp_quantiles, probs = probs)
        wd_cv <- numeric(length(probs))
        for (a in seq_along(probs)) {
          weights <- fixed_cv[a, 1L] / fixed_cv[a, ]
          wd_stats <- apply(t(t(sub_stats) * weights), 1L, max)
          wd_cv[a] <- bp_quantiles(wd_stats, probs[a])
        }
        out[[out_i]] <- bp_cv_rows(
          stat = "cvwd",
          q = qg,
          trm = trmg,
          M = Mg,
          values = wd_cv,
          probs = probs,
          rep = rep,
          n_grid = n_grid,
          scale = "bp_table"
        )
        out_i <- out_i + 1L
      }
    }
  }

  if (NROW(sequential)) {
    seq_groups <- unique(sequential[, c("q", "trm"), drop = FALSE])
    for (g in seq_len(NROW(seq_groups))) {
      qg <- seq_groups$q[g]
      trmg <- seq_groups$trm[g]
      seq_sub <- sequential[sequential$q == qg & sequential$trm == trmg, , drop = FALSE]
      max_l <- max(seq_sub$l)
      one_break <- matrix(NA_real_, nrow = rep, ncol = max_l + 1L)
      label <- sprintf("sequential q=%s trm=%s", qg, trmg)
      start_time <- Sys.time()
      if (isTRUE(verbose)) {
        message(sprintf(
          "%s: starting %d one-break draws on %d core(s)",
          label,
          rep * (max_l + 1L),
          n_cores
        ))
      }
      one_break <- matrix(
        bp_simulate_one_break_draws(
          q = qg,
          trm = trmg,
          n_draws = rep * (max_l + 1L),
          n_grid = n_grid,
          n_cores = n_cores,
          seed = if (is.null(seed)) NULL else as.integer(seed + 100000L + g)
        ),
        nrow = rep,
        ncol = max_l + 1L
      )
      if (isTRUE(verbose)) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        message(sprintf("%s: finished in %.1fs", label, elapsed))
      }

      for (lg in sort(unique(seq_sub$l))) {
        seq_stats <- apply(one_break[, seq_len(lg + 1L), drop = FALSE], 1L, max)
        out[[out_i]] <- bp_cv_rows(
          stat = "supseq",
          q = qg,
          trm = trmg,
          l = lg,
          values = bp_quantiles(seq_stats, probs),
          probs = probs,
          rep = rep,
          n_grid = n_grid,
          scale = "bp_table"
        )
        out_i <- out_i + 1L
      }
    }
  }

  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}

bp_design_matrices <- function(data, stat) {
  if (stat == "cvlr") {
    x1 <- cbind(
      intercept = 1,
      q = data$q,
      q2 = data$q^2,
      k = data$k,
      trm = data$trm,
      q_over_trm = data$q / data$trm
    )
    x2 <- cbind(
      inv_k = 1 / data$k,
      inv_trm_k = 1 / (data$trm * data$k)
    )
    is_bp_table <- "scale" %in% names(data) && all(stats::na.omit(unique(data$scale)) == "bp_table")
    y <- if (is_bp_table) data$cv else data$cv / data$k
  } else if (stat %in% c("cvud", "cvwd")) {
    x1 <- cbind(
      intercept = 1,
      q = data$q,
      trm = data$trm
    )
    x2 <- cbind(
      q = data$q,
      trm = data$trm
    )
    y <- data$cv
  } else if (stat == "supseq") {
    x1 <- cbind(
      intercept = 1,
      q = data$q,
      q2 = data$q^2,
      lp1 = data$l + 1,
      inv_lp1 = 1 / (data$l + 1),
      trm = data$trm
    )
    x2 <- cbind(
      q = data$q,
      q2 = data$q^2
    )
    y <- data$cv
  } else {
    stop("Unknown response-surface statistic: ", stat)
  }

  list(y = as.numeric(y), x1 = x1, x2 = x2)
}

bp_fit_response_one <- function(y, x1, x2, maxit = 5000L) {
  gamma0 <- qr.coef(qr(x1), y)
  gamma0[is.na(gamma0)] <- 0
  delta0 <- rep(0, ncol(x2))
  par0 <- c(gamma0, delta0)
  p1 <- ncol(x1)

  objective <- function(par) {
    gamma <- par[seq_len(p1)]
    delta <- par[-seq_len(p1)]
    pred <- as.vector(x1 %*% gamma) * exp(as.vector(x2 %*% delta))
    sum((y - pred)^2)
  }

  fit <- stats::optim(
    par = par0,
    fn = objective,
    method = "BFGS",
    control = list(maxit = maxit)
  )

  gamma <- fit$par[seq_len(p1)]
  delta <- fit$par[-seq_len(p1)]
  names(gamma) <- colnames(x1)
  names(delta) <- colnames(x2)

  list(
    gamma = gamma,
    delta = delta,
    sse = fit$value,
    convergence = fit$convergence,
    counts = fit$counts
  )
}

fit_bp_response_surfaces <- function(cv_data,
                                     stats = c("cvlr", "cvud", "cvwd", "supseq"),
                                     probs = bp_cv_probs,
                                     maxit = 5000L) {
  fits <- list()

  for (stat in stats) {
    stat_data <- cv_data[cv_data$stat == stat, , drop = FALSE]
    if (!NROW(stat_data)) {
      next
    }

    fits[[stat]] <- list()
    for (prob in as.numeric(probs)) {
      sub <- stat_data[abs(stat_data$prob - prob) < 1e-12, , drop = FALSE]
      if (!NROW(sub)) {
        next
      }
      design <- bp_design_matrices(sub, stat)
      fits[[stat]][[as.character(prob)]] <- bp_fit_response_one(
        y = design$y,
        x1 = design$x1,
        x2 = design$x2,
        maxit = maxit
      )
    }
  }

  class(fits) <- "bp_response_surface_fits"
  fits
}

predict_bp_response_surface <- function(fits, newdata, stat, prob) {
  fit <- fits[[stat]][[as.character(prob)]]
  if (is.null(fit)) {
    stop("No fitted response surface for this statistic/probability.")
  }

  design <- bp_design_matrices(newdata, stat)
  pred <- as.vector(design$x1 %*% fit$gamma) * exp(as.vector(design$x2 %*% fit$delta))
  is_bp_table <- "scale" %in% names(newdata) && all(stats::na.omit(unique(newdata$scale)) == "bp_table")
  if (stat == "cvlr" && !is_bp_table) {
    pred <- pred * newdata$k
  }
  pred
}

bp_surface_coefficients <- function(fits) {
  rows <- list()
  i <- 1L
  for (stat in names(fits)) {
    for (prob in names(fits[[stat]])) {
      fit <- fits[[stat]][[prob]]
      rows[[i]] <- data.frame(
        stat = stat,
        prob = as.numeric(prob),
        block = "gamma",
        name = names(fit$gamma),
        value = as.numeric(fit$gamma),
        sse = fit$sse,
        convergence = fit$convergence,
        stringsAsFactors = FALSE
      )
      i <- i + 1L
      rows[[i]] <- data.frame(
        stat = stat,
        prob = as.numeric(prob),
        block = "delta",
        name = names(fit$delta),
        value = as.numeric(fit$delta),
        sse = fit$sse,
        convergence = fit$convergence,
        stringsAsFactors = FALSE
      )
      i <- i + 1L
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
