### Test discovery of long chains in DAGs
### LK 2025

set.seed(42)

### DEPs
devtools::load_all()
library("comets")
library("tidyverse")
save <- TRUE

### FUNs
dgp <- function(n, sd = 1) {
  A <- sd * rnorm(n)
  B <- A + sd * rnorm(n)
  C <- B + sd * rnorm(n)
  D <- C + sd * rnorm(n)
  E <- D + sd * rnorm(n)
  F <- E + sd * rnorm(n)
  data.frame(A = A, B = B, C = C, D = D, E = E, F = F)
}

gt <- matrix(
  c(
    0, 1, 0, 0, 0, 0,
    1, 0, 1, 0, 0, 0,
    0, 1, 0, 1, 0, 0,
    0, 0, 1, 0, 1, 0,
    0, 0, 0, 1, 0, 1,
    0, 0, 0, 0, 1, 0
  ),
  nrow = 6, ncol = 6, byrow = TRUE
)
dimnames(gt) <- list(LETTERS[1:6], LETTERS[1:6])

### PARs
nsim <- 20
ns <- c(50, 100, 300, 1000)
ds <- 3:4
levs <- c(0.001, 0.01, 0.1)
types <- c("hard", "sigmoid", "mixed", "inv-weight", "log-weight", "identity")

res <- lapply(ds, \(d) {
  lapply(ns, \(n) {
    lapply(seq_len(nsim), \(iter) {
      dd <- dgp(n)[, seq_len(d)]
      tsts <- learn_graph(dd,
        return_tests_only = TRUE,
        test_args = list(reg_YonZ = "lrm", reg_XonZ = "lrm")
      )
      lapply(levs, \(lev) {
        lapply(types, \(type) {
          cat("\nRunning d =", d, "n =", n, "lev =", lev, "type =", type, "iter =", iter, "\n")
          trf <- switch(type,
            "hard" = \(x) as.numeric(x <= lev),
            "sigmoid" = \(x) 1 - plogis(x, location = lev, scale = 0.001),
            "mixed" = \(x) (1 - x) * (x <= lev),
            "inv-weight" = \(x) as.numeric(x <= lev),
            "log-weight" = \(x) as.numeric(x <= lev),
            "identity" = \(x) 1 - x
          )
          wtype <- switch(type,
            "hard" = "const",
            "sigmoid" = "const",
            "mixed" = "const",
            "inv-weight" = "inv",
            "log-weight" = "log",
            "identity" = "const"
          )
          tmp <- capture.output(lG <- graph_optim(tsts,
            d = d, V = LETTERS[1:d],
            trafo = trf, weight_type = wtype, cache = TRUE,
            test_args = list(reg_YonZ = "lrm", reg_XonZ = "lrm")
          ))
          data.frame(
            iter = iter,
            d = d,
            lev = lev,
            n = n,
            type = type,
            SHD = mean(abs(1 * pcalg::dag2essgraph(lG$graph) - gt[1:d, ][, 1:d]))
          )
        }) |> do.call("rbind", args = _)
      }) |> do.call("rbind", args = _)
    }) |> do.call("rbind", args = _)
  }) |> do.call("rbind", args = _)
}) |> do.call("rbind", args = _)

res |>
  group_by(d, n, lev, type) |>
  summarize(SHD = mean(SHD), SD = sd(SHD))

ggplot(res, aes(x = n, y = SHD, color = ordered(lev))) +
  stat_summary(fun.data = mean_se) +
  stat_summary(fun.data = mean_se, geom = "line") +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Sample size", y = "SHD", color = "Level") +
  facet_grid(d ~ type, labeller = \(x) label_both(x, sep = " = ")) +
  theme(text = element_text(size = 13.5))

if (save) {
  ggsave("./inst/figures/sim-chain-new.pdf", height = 3 * length(ds), width = 3 * length(types))
}
