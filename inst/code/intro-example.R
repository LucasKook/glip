### Example where oracle algorithms output wrong
### but global optimization approach works
### LK 2025

set.seed(42)

devtools::load_all()
library("tidyverse")
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
alpha <- 0.05
nsim <- 50
n <- 3e2

out <- lapply(cors, \(tcor) {
  cat("\nRunning cor =", tcor, "\n")
  pb <- txtProgressBar(0, nsim, style = 3, width = 60)
  lapply(1:nsim, \(iter) {
    setTxtProgressBar(pb, iter)
    ### Generate data
    d <- dgp(n, cor = tcor)
    ### Run GLIP
    tmp <- capture.output(
      lP <- learn_graph(d,
        alpha = alpha,
        comets = FALSE,
        test = "zf",
        trafo = \(x) 1 * (x <= alpha),
        mode = "admg"
      )
    )
    ### Run FCI
    out <- fci(list(C = cor(d), n = NROW(d)), gaussCItest, alpha, colnames(d),
      selectionBias = FALSE
    )
    ### Return
    tibble(cor = tcor, GLIP = list(lP$computed), FCI = list(out@amat), iter = iter)
  }) |> bind_rows()
}) |> bind_rows()

res <- out |>
  pivot_longer(c("GLIP", "FCI"), names_to = "method", values_to = "output") |>
  group_by(cor, method) |>
  count(output) |>
  filter(n == max(n))

res$output
