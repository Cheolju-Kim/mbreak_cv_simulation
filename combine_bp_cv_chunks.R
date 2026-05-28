#!/usr/bin/env Rscript

# Combine chunked BP critical-value simulations and fit response surfaces.
#
# Example:
#   Rscript combine_bp_cv_chunks.R --mode=final --rep=10000 --n-grid=1000

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else "combine_bp_cv_chunks.R"
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
    kv <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- kv[[1L]]
    value <- if (length(kv) == 1L) TRUE else paste(kv[-1L], collapse = "=")
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

read_config <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parts <- strsplit(lines, ":", fixed = TRUE)
  keys <- trimws(vapply(parts, `[`, character(1L), 1L))
  values <- trimws(vapply(parts, function(x) paste(x[-1L], collapse = ":"), character(1L)))
  stats::setNames(as.list(values), keys)
}

config_value <- function(config, key, default = NA_character_) {
  value <- config[[key]]
  if (is.null(value) || identical(value, "")) {
    default
  } else {
    value
  }
}

make_cv_key <- function(data) {
  paste(
    data$stat,
    data$q,
    sprintf("%.6f", data$trm),
    ifelse(is.na(data$k), "", data$k),
    ifelse(is.na(data$l), "", data$l),
    ifelse(is.na(data$M), "", data$M),
    sprintf("%.3f", data$prob),
    data$rep,
    data$n_grid,
    data$scale,
    sep = "|"
  )
}

opts <- parse_cli(commandArgs(trailingOnly = TRUE))
mode <- as.character(arg_value(opts, "mode", "final"))
rep <- arg_int(opts, "rep", 10000L)
n_grid <- arg_int(opts, "n-grid", 1000L)
do_fit <- arg_bool(opts, "fit", TRUE)

input_root <- normalizePath(
  as.character(arg_value(opts, "input", file.path(script_dir, "cv_output", "scc"))),
  mustWork = FALSE
)
output_root <- normalizePath(
  as.character(arg_value(opts, "out", file.path(script_dir, "cv_output", "combined"))),
  mustWork = FALSE
)

if (!dir.exists(input_root)) {
  stop("Input directory does not exist: ", input_root)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

chunk_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
chunk_dirs <- chunk_dirs[file.exists(file.path(chunk_dirs, "config.txt"))]

manifest <- do.call(rbind, lapply(chunk_dirs, function(dir) {
  config <- read_config(file.path(dir, "config.txt"))
  rds_path <- file.path(dir, paste0("cv_", config_value(config, "mode"), ".rds"))
  csv_path <- file.path(dir, paste0("cv_", config_value(config, "mode"), ".csv"))
  data.frame(
    dir = dir,
    mode = config_value(config, "mode"),
    rep = as.integer(config_value(config, "rep")),
    n_grid = as.integer(config_value(config, "n_grid")),
    seed = as.integer(config_value(config, "seed")),
    q_min = suppressWarnings(as.integer(config_value(config, "q_min"))),
    q_max = suppressWarnings(as.integer(config_value(config, "q_max"))),
    fixed_points = as.integer(config_value(config, "fixed_points")),
    double_points = as.integer(config_value(config, "double_points")),
    sequential_points = as.integer(config_value(config, "sequential_points")),
    has_output = file.exists(rds_path) || file.exists(csv_path),
    stringsAsFactors = FALSE
  )
}))

selected <- manifest[
  manifest$mode == mode & manifest$rep == rep & manifest$n_grid == n_grid & manifest$has_output,
  ,
  drop = FALSE
]
if (!NROW(selected)) {
  stop("No matching chunk directories found for mode=", mode, ", rep=", rep, ", n_grid=", n_grid)
}

cv_list <- vector("list", NROW(selected))
for (i in seq_len(NROW(selected))) {
  dir <- selected$dir[i]
  rds_path <- file.path(dir, paste0("cv_", mode, ".rds"))
  csv_path <- file.path(dir, paste0("cv_", mode, ".csv"))
  if (file.exists(rds_path)) {
    data <- readRDS(rds_path)
  } else if (file.exists(csv_path)) {
    data <- read.csv(csv_path, stringsAsFactors = FALSE)
  } else {
    stop("No cv output found in: ", dir)
  }
  data$source_dir <- basename(dir)
  cv_list[[i]] <- data
}

cv_data <- do.call(rbind, cv_list)
rownames(cv_data) <- NULL

dups <- duplicated(make_cv_key(cv_data))
if (any(dups)) {
  dup_keys <- unique(make_cv_key(cv_data)[dups])
  stop("Duplicate CV rows found. First duplicate key: ", dup_keys[[1L]])
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(output_root, paste0(mode, "_rep", rep, "_ngrid", n_grid, "_", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

q_seen <- sort(unique(cv_data$q))
expected_q <- sort(unique(bp_extended_response_grid(n_grid = n_grid)$fixed$q))
missing_q <- setdiff(expected_q, q_seen)

summary_by_stat <- aggregate(
  cv ~ stat,
  data = cv_data,
  FUN = length
)
names(summary_by_stat) <- c("stat", "rows")

write.csv(selected, file.path(out_dir, "combine_manifest.csv"), row.names = FALSE)
write.csv(summary_by_stat, file.path(out_dir, "combined_summary_by_stat.csv"), row.names = FALSE)
write.csv(data.frame(q = q_seen), file.path(out_dir, "combined_q_grid.csv"), row.names = FALSE)
write.csv(data.frame(q = missing_q), file.path(out_dir, "missing_q_grid.csv"), row.names = FALSE)
saveRDS(cv_data, file.path(out_dir, paste0("cv_", mode, "_combined.rds")))
write.csv(cv_data, file.path(out_dir, paste0("cv_", mode, "_combined.csv")), row.names = FALSE)

if (length(missing_q)) {
  warning("Missing expected q values: ", paste(missing_q, collapse = ", "))
}

if (do_fit) {
  fits <- fit_bp_response_surfaces(cv_data)
  coef_tab <- bp_surface_coefficients(fits)
  saveRDS(fits, file.path(out_dir, paste0("response_surface_fits_", mode, "_combined.rds")))
  write.csv(
    coef_tab,
    file.path(out_dir, paste0("response_surface_coefficients_", mode, "_combined.csv")),
    row.names = FALSE
  )
}

cat("Combined output directory:", out_dir, "\n")
cat("Chunks used:", NROW(selected), "\n")
cat("q grid:", paste(q_seen, collapse = ", "), "\n")
print(summary_by_stat, row.names = FALSE)
