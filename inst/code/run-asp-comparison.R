### Timing comparison with Hyttinen et al 2014 ASP
### LK 2025

devtools::load_all()
setwd("./inst/asp/hyttinen2014uai_ver6/pkg/R")
source("./load.R")
loud()
library("tidyverse")

ncores <- max(7, parallel::detectCores(logical = TRUE) - 2)
cache <- TRUE

nsim <- 20
ds <- 3:4
N <- 1e3
n <- 0

out <- lapply(ds, \(d) {
  cat("\nRunning d =", d, "\n")
  pb <- txtProgressBar(0, nsim, style = 3, width = 60)
  lapply(1:nsim, \(iter) {
    setTxtProgressBar(pb, iter)
    tmp <- capture.output({
      n <<- d
      res <- pipeline(n = d, N = N, test = "classic", verbose = 0)
    })
    data.frame(d = d, n = N, iter = iter, res, row.names = NULL)
  }) |> do.call("rbind", args = _)
}) |> do.call("rbind", args = _)

mqr <- function(x, ...) {
  data.frame(y = median(x), ymin = quantile(x, 0.25), ymax = quantile(x, 0.75))
}

out |>
  pivot_longer(c("glip", "asp"), names_to = "method", values_to = "time") |>
  mutate(time = as.numeric(time)) |>
  ggplot(aes(x = ordered(d), y = time, color = method)) +
  geom_boxplot(width = 0.3) +
  theme_bw() +
  labs(x = "number of nodes", y = "median runtime time in seconds") +
  scale_y_log10()
