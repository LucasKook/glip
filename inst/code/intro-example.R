### Example where oracle algorithms output wrong
### but global optimization approach works
### LK 2025

set.seed(42)

devtools::load_all()
library("comets")
library("pcalg")

dgp <- function(n, cor = 0.9) {
  U <- mvtnorm::rmvnorm(n, sigma = matrix(c(1, cor, cor, 1), nrow = 2, ncol = 2))
  X1 <- U[, 1]
  X2 <- U[, 2]
  X3 <- X1 + X2 + rnorm(n)
  X4 <- X3 + rnorm(n)
  data.frame(a = X1, b = X2, c = X3, d = X4)
}

cors <- c(0, 0.9, 0.9999)

tmp <- lapply(cors, \(tcor) {
  d <- dgp(1e4, cor = tcor)

  tmp <- capture.output(
    lP <- learn_graph(d, test_args = list(
      reg_YonZ = "lrm",
      reg_XonZ = "lrm"
    ), mode = "admg")
  )

  cat("\ncor =", tcor)
  cat("\nOptim:\n")
  print(lP$computed)

  cat("\nFCI:\n")
  out <- fci(list(C = cor(d), n = NROW(d)), gaussCItest, 0.05, colnames(d),
    selectionBias = FALSE
  )
  print(out@amat)
})
