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
cache <- TRUE
alpha <- 0.05

# Parameters for running the tests
targs <- list(reg_YonZ = "lrm", reg_XonZ = "lrm")

### Generate random graph and data
graph <- random_graph(d = d, prob = pr, mode = mode)
data <- data.frame(rgraph(graph, n = n))
V <- colnames(data)

### Run CITs
tests <- learn_graph(
  data = data, max_size = ms, mode = mode, test_args = targs,
  return_tests_only = TRUE
)

### ORACLE
gt <- switch(mode,
  "dag" = graph$DAG,
  "admg" = graph$ADMG
)
ORACLE <- .compute_graphical_representation(gt, ms, mode)

### GLIP
capt <- capture.output(lG <- .get_opt(mode)(tests,
  d = d, max_size = ms,
  V = V, cache = cache,
  trafo = \(x) as.numeric(x <= alpha),
  gurobi_args = list(
    Threads = ncores
  ), mode = mode
))

GLIP <- .compute_graphical_representation(lG$graph, ms, mode)
runtime_GLIP <- lG$optim$runtime

### PC ALG
tstart <- Sys.time()
pcres <- pcalg::pc(list(tests = tests, V = V), lookup_ci, labels = V, alpha = alpha)
tstop <- Sys.time()
runtime_PC <- tstop - tstart
pcout <- as(pcres@graph, "matrix")
PC <- pcout

### FCI ALG
tstart <- Sys.time()
fcires <- pcalg::fciPlus(list(tests = tests, V = V), lookup_ci,
  labels = V, alpha = alpha, selectionBias = FALSE, verbose = FALSE
)
tstop <- Sys.time()
runtime_FCI <- tstop - tstart
fciout <- as(fcires@amat, "matrix")
FCI <- fciout
