#!/usr/bin/env Rscript

# Offline driver for Bai-Perron critical-value simulation.
#
# Examples:
#   Rscript run_bp_cv_simulation.R --mode=smoke --rep=500 --n-grid=300 --cores=4
#   Rscript run_bp_cv_simulation.R --mode=bp --rep=10000 --n-grid=1000 --cores=14
#   Rscript run_bp_cv_simulation.R --mode=final --rep=10000 --n-grid=1000 --cores=14
#   Rscript run_bp_cv_simulation.R --mode=final --q-min=70 --q-max=80

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "run_bp_cv_simulation.R"
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
if (!file.exists(file.path(script_dir, "cv_simulation.R"))) {
  script_dir <- getwd()
}

source(file.path(script_dir, "cv_simulation.R"))

parse_cli <- function(args) {
  opts <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    arg <- sub("^--", "", arg)
    kv <- strsplit(arg, "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    value <- if (length(kv) == 1L) TRUE else paste(kv[-1], collapse = "=")
    opts[[key]] <- value
  }
  opts
}

arg_value <- function(opts, key, default = NULL) {
  if (!is.null(opts[[key]])) opts[[key]] else default
}

arg_int <- function(opts, key, default) {
  as.integer(arg_value(opts, key, default))
}

arg_bool <- function(opts, key, default = TRUE) {
  value <- arg_value(opts, key, default)
  if (is.logical(value)) {
    return(value)
  }
  tolower(as.character(value)) %in% c("1", "true", "t", "yes", "y")
}

arg_optional_int <- function(opts, key) {
  value <- arg_value(opts, key, NULL)
  if (is.null(value) || identical(value, "") || is.na(value)) {
    return(NULL)
  }
  as.integer(value)
}

default_cores <- function() {
  cores <- parallel::detectCores(logical = TRUE)
  max(1L, min(14L, cores - 4L))
}

mode_defaults <- function(mode) {
  switch(
    mode,
    smoke = list(rep = 500L, n_grid = 300L, fixed_scale = "bp_table", seed = 101L),
    bp = list(rep = 10000L, n_grid = 1000L, fixed_scale = "bp_table", seed = 202L),
    final = list(rep = 10000L, n_grid = 1000L, fixed_scale = "mbreak_lr", seed = 303L),
    stop("Unknown mode: ", mode)
  )
}

make_mode_grid <- function(mode, n_grid) {
  switch(
    mode,
    smoke = bp_response_grid(
      q = 1:3,
      trm = c(0.05, 0.15),
      k_cap = 3L,
      l_cap = 3L,
      double_m_cap = 3L,
      n_grid = n_grid
    ),
    bp = bp_response_grid(
      q = 1:10,
      trm = c(0.05, 0.10, 0.15, 0.20, 0.25),
      k_cap = 9L,
      l_cap = 9L,
      double_m_cap = 5L,
      n_grid = n_grid
    ),
    final = bp_extended_response_grid(n_grid = n_grid),
    stop("Unknown mode: ", mode)
  )
}

filter_grid_q <- function(grid, q_min = NULL, q_max = NULL) {
  if (is.null(q_min) && is.null(q_max)) {
    return(grid)
  }

  keep <- function(data) {
    if (!NROW(data)) {
      return(data)
    }
    ok <- rep(TRUE, NROW(data))
    if (!is.null(q_min)) {
      ok <- ok & data$q >= q_min
    }
    if (!is.null(q_max)) {
      ok <- ok & data$q <= q_max
    }
    data[ok, , drop = FALSE]
  }

  list(
    fixed = keep(grid$fixed),
    double = keep(grid$double),
    sequential = keep(grid$sequential)
  )
}

write_config <- function(path, config) {
  lines <- sprintf("%s: %s", names(config), unlist(config, use.names = FALSE))
  writeLines(lines, path)
}

bp1998_table_i <- function() {
  txt <- "
q,prob,k1,k2,k3,k4,k5,k6,k7,k8,k9,cvud,cvwd
1,0.90,8.02,7.87,7.07,6.61,6.14,5.74,5.40,5.09,4.81,8.78,9.14
1,0.95,9.63,8.78,7.85,7.21,6.69,6.23,5.86,5.51,5.20,10.17,10.91
1,0.975,11.17,9.81,8.52,7.79,7.22,6.70,6.27,5.92,5.56,11.52,12.53
1,0.99,13.58,10.95,9.37,8.50,7.85,7.21,6.75,6.33,5.98,13.74,15.02
2,0.90,11.02,10.48,9.61,8.99,8.50,8.06,7.66,7.32,7.01,11.69,12.33
2,0.95,12.89,11.60,10.46,9.71,9.12,8.65,8.19,7.79,7.46,13.27,14.19
2,0.975,14.53,12.64,11.20,10.29,9.69,9.10,8.64,8.18,7.80,14.69,16.04
2,0.99,16.64,13.78,12.06,11.00,10.28,9.65,9.11,8.66,8.22,16.79,18.11
3,0.90,13.43,12.73,11.76,11.04,10.49,10.02,9.59,9.21,8.86,14.05,14.76
3,0.95,15.37,13.84,12.64,11.83,11.15,10.61,10.14,9.71,9.32,15.80,16.82
3,0.975,17.17,14.91,13.44,12.49,11.75,11.13,10.62,10.14,9.72,17.36,18.79
3,0.99,19.25,16.27,14.48,13.40,12.56,11.80,11.22,10.67,10.19,19.38,20.81
4,0.90,15.53,14.65,13.63,12.91,12.33,11.79,11.34,10.93,10.55,16.17,16.95
4,0.95,17.60,15.84,14.63,13.71,12.99,12.42,11.91,11.49,11.04,17.88,19.07
4,0.975,19.35,16.85,15.44,14.43,13.64,13.01,12.46,11.94,11.49,19.51,20.89
4,0.99,21.20,18.21,16.43,15.21,14.45,13.70,13.04,12.48,12.02,21.25,22.81
5,0.90,17.42,16.45,15.44,14.69,14.05,13.51,13.02,12.59,12.18,17.94,18.85
5,0.95,19.50,17.60,16.40,15.52,14.79,14.19,13.63,13.16,12.70,19.74,20.95
5,0.975,21.47,18.75,17.26,16.13,15.40,14.75,14.19,13.66,13.17,21.57,23.04
5,0.99,23.99,20.18,18.19,17.09,16.14,15.34,14.81,14.26,13.72,24.00,25.46
6,0.90,19.38,18.15,17.17,16.39,15.74,15.18,14.63,14.18,13.74,19.92,20.89
6,0.95,21.59,19.61,18.23,17.27,16.50,15.86,15.29,14.77,14.30,21.90,23.27
6,0.975,23.73,20.80,19.15,18.07,17.21,16.49,15.84,15.29,14.78,23.83,25.22
6,0.99,25.95,22.18,20.29,18.93,17.97,17.20,16.54,15.94,15.35,26.07,27.63
7,0.90,21.23,19.93,18.75,17.98,17.28,16.69,16.16,15.69,15.24,21.79,22.81
7,0.95,23.50,21.30,19.83,18.91,18.10,17.43,16.83,16.28,15.79,23.77,25.02
7,0.975,25.23,22.54,20.85,19.68,18.79,18.03,17.38,16.79,16.31,25.46,26.92
7,0.99,28.01,24.07,21.89,20.68,19.68,18.81,18.10,17.49,16.96,28.02,29.57
8,0.90,22.92,21.56,20.43,19.58,18.84,18.21,17.69,17.19,16.70,23.53,24.55
8,0.95,25.22,23.03,21.48,20.46,19.66,18.97,18.37,17.80,17.30,25.51,26.83
8,0.975,27.21,24.20,22.41,21.29,20.39,19.63,18.98,18.34,17.78,27.32,28.98
8,0.99,29.60,25.66,23.44,22.22,21.22,20.40,19.66,19.03,18.46,29.60,31.32
9,0.90,24.75,23.15,21.98,21.12,20.37,19.72,19.13,18.58,18.09,25.19,26.40
9,0.95,27.08,24.55,23.16,22.08,21.22,20.49,19.90,19.29,18.79,27.28,28.78
9,0.975,29.13,25.92,24.14,22.97,21.98,21.28,20.59,19.98,19.39,29.20,30.82
9,0.99,31.66,27.42,25.13,24.01,23.06,22.18,21.35,20.63,19.94,31.72,33.32
10,0.90,26.13,24.70,23.48,22.57,21.83,21.16,20.57,20.03,19.55,26.66,27.79
10,0.95,28.49,26.17,24.59,23.59,22.71,21.93,21.34,20.74,20.17,28.75,30.16
10,0.975,30.67,27.52,25.69,24.47,23.45,22.71,21.95,21.34,20.79,30.84,32.46
10,0.99,33.62,29.14,26.90,25.58,24.44,23.49,22.75,22.09,21.47,33.86,35.47
"
  wide <- read.csv(text = txt, stringsAsFactors = FALSE)
  fixed <- do.call(rbind, lapply(seq_len(9L), function(k) {
    data.frame(
      stat = "cvlr",
      q = wide$q,
      trm = 0.05,
      k = k,
      l = NA_integer_,
      M = NA_integer_,
      prob = wide$prob,
      bp_cv = wide[[paste0("k", k)]],
      stringsAsFactors = FALSE
    )
  }))
  cvud <- data.frame(
    stat = "cvud", q = wide$q, trm = 0.05, k = NA_integer_,
    l = NA_integer_, M = 5L, prob = wide$prob, bp_cv = wide$cvud
  )
  cvwd <- data.frame(
    stat = "cvwd", q = wide$q, trm = 0.05, k = NA_integer_,
    l = NA_integer_, M = 5L, prob = wide$prob, bp_cv = wide$cvwd
  )
  rbind(fixed, cvud, cvwd)
}

bp1998_table_ii <- function() {
  txt <- "
q,prob,l0,l1,l2,l3,l4,l5,l6,l7,l8,l9
1,0.90,8.02,9.56,10.45,11.07,11.65,12.07,12.47,12.70,13.07,13.34
1,0.95,9.63,11.14,12.16,12.83,13.45,14.05,14.29,14.50,14.69,14.88
1,0.975,11.17,12.88,14.05,14.50,15.03,15.37,15.56,15.73,16.02,16.39
1,0.99,13.58,15.03,15.62,16.39,16.60,16.90,17.04,17.27,17.32,17.61
2,0.90,11.02,12.79,13.72,14.45,14.90,15.35,15.81,16.12,16.44,16.58
2,0.95,12.89,14.50,15.42,16.16,16.61,17.02,17.27,17.55,17.76,17.97
2,0.975,14.53,16.19,17.02,17.55,17.98,18.15,18.46,18.74,18.98,19.22
2,0.99,16.64,17.98,18.66,19.22,20.03,20.87,20.97,21.19,21.43,21.74
3,0.90,13.43,15.26,16.38,17.07,17.52,17.91,18.35,18.61,18.92,19.19
3,0.95,15.37,17.15,17.97,18.72,19.23,19.59,19.94,20.31,21.05,21.20
3,0.975,17.17,18.75,19.61,20.31,21.33,21.59,21.78,22.07,22.41,22.73
3,0.99,19.25,21.33,22.01,22.73,23.13,23.48,23.70,23.79,23.84,24.59
4,0.90,15.53,17.54,18.55,19.30,19.80,20.15,20.48,20.73,20.94,21.10
4,0.95,17.60,19.33,20.22,20.75,21.15,21.55,21.90,22.27,22.63,22.83
4,0.975,19.35,20.76,21.60,22.27,22.84,23.44,23.74,24.14,24.36,24.54
4,0.99,21.20,22.84,24.04,24.54,24.96,25.36,25.51,25.58,25.63,25.88
5,0.90,17.42,19.38,20.46,21.37,21.96,22.47,22.77,23.23,23.56,23.81
5,0.95,19.50,21.43,22.57,23.33,23.90,24.34,24.62,25.14,25.34,25.51
5,0.975,21.47,23.34,24.37,25.14,25.58,25.79,25.96,26.39,26.60,26.84
5,0.99,23.99,25.58,26.32,26.84,27.39,27.86,27.90,28.32,28.38,28.39
6,0.90,19.38,21.51,22.81,23.64,24.19,24.59,24.86,25.27,25.53,25.87
6,0.95,21.59,23.72,24.66,25.29,25.89,26.36,26.84,27.10,27.26,27.40
6,0.975,23.73,25.41,26.37,27.10,27.42,28.02,28.39,28.75,29.13,29.44
6,0.99,25.95,27.42,28.60,29.44,30.18,30.52,30.64,30.99,31.25,31.33
7,0.90,21.23,23.41,24.51,25.07,25.75,26.30,26.74,27.06,27.46,27.70
7,0.95,23.50,25.17,26.34,27.19,27.96,28.25,28.64,28.84,28.97,29.14
7,0.975,25.23,27.24,28.25,28.84,29.14,29.72,30.41,30.76,31.09,31.43
7,0.99,28.01,29.14,30.61,31.43,32.56,32.75,32.90,33.25,33.25,33.85
8,0.90,22.92,25.15,26.38,27.09,27.77,28.15,28.61,28.90,29.19,29.49
8,0.95,25.22,27.18,28.21,28.99,29.54,30.05,30.45,30.79,31.29,31.75
8,0.975,27.21,29.01,30.09,30.79,31.80,32.50,32.81,32.86,33.20,33.60
8,0.99,29.60,31.80,32.84,33.60,34.23,34.57,34.75,35.01,35.50,35.65
9,0.90,24.75,26.99,28.11,29.03,29.69,30.18,30.61,30.93,31.14,31.46
9,0.95,27.08,29.10,30.24,30.99,31.48,32.46,32.71,32.89,33.15,33.43
9,0.975,29.13,31.04,32.48,32.89,33.47,33.98,34.25,34.74,34.88,35.07
9,0.99,31.66,33.47,34.60,35.07,35.49,37.08,37.12,37.23,37.47,37.68
10,0.90,26.13,28.40,29.68,30.62,31.25,31.81,32.37,32.78,33.09,33.53
10,0.95,28.49,30.65,31.90,32.83,33.57,34.27,34.53,35.01,35.33,35.65
10,0.975,30.67,32.87,34.27,35.01,35.86,36.32,36.65,36.90,37.15,37.41
10,0.99,33.62,35.86,36.68,37.41,38.20,38.70,38.91,39.09,39.11,39.12
"
  wide <- read.csv(text = txt, stringsAsFactors = FALSE)
  do.call(rbind, lapply(0:9, function(l) {
    data.frame(
      stat = "supseq",
      q = wide$q,
      trm = 0.05,
      k = NA_integer_,
      l = l,
      M = NA_integer_,
      prob = wide$prob,
      bp_cv = wide[[paste0("l", l)]],
      stringsAsFactors = FALSE
    )
  }))
}

bp1998_benchmark <- function() {
  rbind(bp1998_table_i(), bp1998_table_ii())
}

make_key <- function(data) {
  paste(
    data$stat,
    data$q,
    sprintf("%.6f", data$trm),
    ifelse(is.na(data$k), "", data$k),
    ifelse(is.na(data$l), "", data$l),
    ifelse(is.na(data$M), "", data$M),
    sprintf("%.3f", data$prob),
    sep = "|"
  )
}

compare_bp1998 <- function(cv_data) {
  bench <- bp1998_benchmark()
  sim <- cv_data[abs(cv_data$trm - 0.05) < 1e-12, , drop = FALSE]
  bench$key <- make_key(bench)
  sim$key <- make_key(sim)

  merged <- merge(
    bench,
    sim[, c("key", "cv", "rep", "n_grid", "scale")],
    by = "key",
    all.x = TRUE,
    suffixes = c("_bp", "_sim")
  )
  merged$abs_diff <- abs(merged$cv - merged$bp_cv)
  merged$rel_diff <- merged$abs_diff / merged$bp_cv
  merged[order(merged$stat, merged$q, merged$prob, merged$k, merged$l), ]
}

run_one_mode <- function(mode, opts) {
  defaults <- mode_defaults(mode)
  rep <- arg_int(opts, "rep", defaults$rep)
  n_grid <- arg_int(opts, "n-grid", defaults$n_grid)
  seed <- arg_int(opts, "seed", defaults$seed)
  n_cores <- arg_int(opts, "cores", default_cores())
  fixed_scale <- as.character(arg_value(opts, "fixed-scale", defaults$fixed_scale))
  do_fit <- arg_bool(opts, "fit", TRUE)
  q_min <- arg_optional_int(opts, "q-min")
  q_max <- arg_optional_int(opts, "q-max")

  out_root <- normalizePath(
    as.character(arg_value(opts, "out", file.path(script_dir, "cv_output"))),
    mustWork = FALSE
  )
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  job_id <- Sys.getenv("JOB_ID", unset = "")
  q_tag <- if (is.null(q_min) && is.null(q_max)) {
    NULL
  } else {
    paste0("q", if (is.null(q_min)) "min" else q_min, "_", if (is.null(q_max)) "max" else q_max)
  }
  job_tag <- if (nzchar(job_id)) paste0("job", job_id) else paste0("seed", seed)
  out_name <- paste(c(mode, q_tag, job_tag, stamp), collapse = "_")
  out_dir <- file.path(out_root, out_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  grid <- make_mode_grid(mode, n_grid)
  grid <- filter_grid_q(grid, q_min = q_min, q_max = q_max)
  if (!NROW(grid$fixed) && !NROW(grid$double) && !NROW(grid$sequential)) {
    stop("The selected q range has no grid points.")
  }
  config <- list(
    mode = mode,
    rep = rep,
    n_grid = n_grid,
    seed = seed,
    n_cores = n_cores,
    fixed_scale = fixed_scale,
    q_min = if (is.null(q_min)) "" else q_min,
    q_max = if (is.null(q_max)) "" else q_max,
    fixed_points = nrow(grid$fixed),
    double_points = nrow(grid$double),
    sequential_points = nrow(grid$sequential),
    output = out_dir
  )
  write_config(file.path(out_dir, "config.txt"), config)

  message("Output directory: ", out_dir)
  message("Running mode=", mode, ", rep=", rep, ", n_grid=", n_grid,
          ", cores=", n_cores, ", fixed_scale=", fixed_scale)

  start_time <- Sys.time()
  cv_data <- simulate_bp_response_surface_data(
    grid = grid,
    rep = rep,
    n_grid = n_grid,
    seed = seed,
    fixed_scale = fixed_scale,
    n_cores = n_cores,
    verbose = TRUE
  )
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  message(sprintf("Simulation finished in %.1f seconds.", elapsed))

  saveRDS(cv_data, file.path(out_dir, paste0("cv_", mode, ".rds")))
  write.csv(cv_data, file.path(out_dir, paste0("cv_", mode, ".csv")), row.names = FALSE)

  if (mode == "bp") {
    comparison <- compare_bp1998(cv_data)
    write.csv(comparison, file.path(out_dir, "bp1998_table_comparison.csv"), row.names = FALSE)

    summary <- stats::aggregate(
      abs_diff ~ stat + prob,
      data = comparison,
      FUN = function(x) c(mean = mean(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
    )
    write.csv(summary, file.path(out_dir, "bp1998_table_comparison_summary.csv"), row.names = FALSE)
  }

  if (do_fit) {
    fits <- fit_bp_response_surfaces(cv_data)
    coef_tab <- bp_surface_coefficients(fits)
    saveRDS(fits, file.path(out_dir, paste0("response_surface_fits_", mode, ".rds")))
    write.csv(coef_tab, file.path(out_dir, paste0("response_surface_coefficients_", mode, ".csv")),
              row.names = FALSE)
  }

  invisible(out_dir)
}

main <- function() {
  opts <- parse_cli(commandArgs(trailingOnly = TRUE))
  mode <- as.character(arg_value(opts, "mode", "smoke"))
  modes <- if (mode == "all") c("smoke", "bp", "final") else mode

  for (m in modes) {
    run_one_mode(m, opts)
  }
}

main()
