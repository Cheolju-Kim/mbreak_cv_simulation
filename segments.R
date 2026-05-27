# Stacked-index helpers for break-date based calculations.

#' Stacked Row Indices
#'
#' Returns the row indices for observations from `start` to `end` in an
#' equation-stacked system with `n` equations.
stacked_rows <- function(start, end, n) {
  ((start - 1) * n + 1):(end * n)
}

#' Upper-Triangular Segment Index
#'
#' Maps a segment `(start, end)` into the compact vector index used to store
#' all possible segment objective values.
bigvec_index <- function(start, end, bigt) {
  (start - 1) * bigt - (start - 1) * (start - 2) / 2 + (end - start + 1)
}
