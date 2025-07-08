#' @exportS3Method print dmg
print.dmg <- function(x, ...) {
  cat("Directed adjacency matrix:\n\n")
  print(x$M1)
  cat("\nBidirected adjacency matrix:\n\n")
  print(x$M2)
}

.fill <- function(x) {
  ifelse(identical(x, numeric(0)), 0, x)
}
