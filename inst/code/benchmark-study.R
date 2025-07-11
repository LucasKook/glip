### Benchmark glip against existing algorithms
### LK 2025

set.seed(1)

### DEPs
devtools::load_all()

### PARs

# Parameters for generating random graph
d <- 5
pr <- 0.5
mode <- "dag"

# Parameters for simulating data from random graph
n <- 1e2

# Parameters for running the optimization
ms <- d - 2
ncores <- max(7, parallel::detectCores(logical = TRUE) - 2)

# Parameters for running the tests
targs <- list(reg_YonZ = "lrm", reg_XonZ = "lrm")

### RUN
graph <- random_graph(d = d, prob = pr, mode = mode)
data <- data.frame(rgraph(graph, n = n))
lG <- learn_graph(
  data = data, max_size = ms, mode = mode,
  gurobi_args = list(Threads = ncores), test_args = targs
)

lG$computed

gt <- switch(mode,
  "dag" = graph$DAG,
  "admg" = graph$ADMG
)
.compute_graphical_representation(gt, ms, mode)
