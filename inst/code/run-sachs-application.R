### Sachs application
### LK 2025

### Dependencies
devtools::load_all()
library("tidyverse")

### Parameters
ncores <- max(7, parallel::detectCores(logical = TRUE) - 2)
mode <- "dag"

### Read data
nms <- c(
  "Raf", "Mek", "PLCg", "PIP2", "PIP3", "Erk", "Akt", "PKA", "PKC",
  "p38", "JNK"
)
dat <- readxl::read_xls(file.path("./inst/data/sachs", "cd3cd28.xls")) |>
  mutate_all(log)
colnames(dat) <- nms

### Oracle graph
sachs <- matrix(c(
  0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0
), ncol = 11)
dimnames(sachs) <- list(nms, nms)

d <- NCOL(dat)

### Learn graph
lG <- learn_graph(
  dat,
  max_size = 1, mode = mode,
  trafo = \(x) as.numeric(x <= 0.01),
  test_args = list(reg_YonZ = "lrm", reg_XonZ = "lrm"),
  gurobi_args = list(Threads = ncores)
)

### Evaluate
G <- lG$graph$graph
learned <- .compute_graphical_representation(G, d - 2, mode)
groundtruth <- .compute_graphical_representation(sachs, d - 2, mode)

shd(learned, groundtruth)
sep(learned, groundtruth, mode, 1, precomputed_predicted = 1 - lG$graph$tests$dcon)
prf1(learned, groundtruth)
