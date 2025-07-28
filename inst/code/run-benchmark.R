### Benchmark glip against existing algorithms
### LK 2025

### DEPs
devtools::load_all()
library("pcalg")
library("reticulate")
use_condaenv("glip", required = TRUE)
utils <- import("dagma.utils", convert = TRUE)
dagma <- import("dagma.linear", convert = TRUE)
cd <- import("CausalDisco.baselines", convert = TRUE)

### PARs

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
mode <- darg(args[1], "dag")
d <- as.numeric(darg(args[2], 3))
ms <- as.numeric(darg(args[3], d - 2))
ms <- ifelse(ms == -1, d - 2, ms)
degree <- as.numeric(darg(args[4], 2))
n <- as.numeric(darg(args[5], 1e3))
seed <- as.numeric(darg(args[6], 12))
use_oracle_tests <- FALSE
save <- TRUE

# Parameters for running the optimization
ncores <- max(7, parallel::detectCores(logical = TRUE) - 2)
cache <- TRUE
alpha <- 0.01

# Parameters for running the tests
targs <- list(reg_YonZ = "lrm", reg_XonZ = "lrm")

# Output file
outdir <- file.path(
  "inst", "results", "benchmark", Sys.Date()
)
fout <- paste0(
  "res-mode_", mode, "-d_", d, "-ms_", ms, "-degree_",
  degree, "-n_", n, "-seed_", seed, ".rds"
)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

### Generate random graph and data
set.seed(tseed <- 1e4 + 3e4 * (mode == "dag") + n + seed)
graph <- random_graph(d = d, prob = pr, mode = mode)
data <- data.frame(py_data <- rgraph(graph, n = n))
py_data <- r_to_py(py_data)$copy()
V <- colnames(data)

### ORACLE
gt <- switch(mode,
  "dag" = graph$DAG,
  "admg" = graph$ADMG
)
ORACLE <- .compute_graphical_representation(gt, ms, mode)

### Run CITs
if (use_oracle_tests) {
  tests <- .compute_oracle_tests(gt, ms, mode)
} else {
  tests <- learn_graph(
    data = data, max_size = ms, mode = mode, test_args = targs,
    return_tests_only = TRUE
  )
}

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
runtime_GLIP <- as.difftime(lG$optim$runtime, units = "secs")

### PC ALG (only under causal sufficiency/DAG case)
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

### NOTEARS
model <- dagma$DagmaLinear(loss_type = "l2")
tstart <- Sys.time()
nto <- 1 * (model$fit(py_data, lambda1 = 0.02) != 0)
tstop <- Sys.time()
runtime_NOTEARS <- tstop - tstart
dimnames(nto) <- list(V, V)
NOTEARS <- .compute_graphical_representation(nto, ms, mode)

### R2sortability
tstart <- Sys.time()
r2s <- 1 * (cd$r2_sort_regress(py_data) != 0)
tstop <- Sys.time()
runtime_R2SORT <- tstop - tstart
dimnames(r2s) <- list(V, V)
R2SORT <- .compute_graphical_representation(r2s, ms, mode)

### Evaluate and summarize results
outputs <- list(
  GLIP = GLIP,
  PC = PC,
  FCI = FCI,
  NOTEARS = NOTEARS,
  R2SORT = R2SORT
)
timings <- list(
  GLIP = runtime_GLIP,
  PC = runtime_PC,
  FCI = runtime_FCI,
  NOTEARS = runtime_NOTEARS,
  R2SORT = runtime_R2SORT
)

if (mode == "admg") { # remove PC in case of ADMGs
  outputs <- outputs[-grep("PC", names(outputs))]
  timings <- timings[-grep("PC", names(timings))]
}

res <- lapply(seq_along(outputs), \(idx) {
  learned <- outputs[[idx]]
  SHD <- shd(learned, ORACLE)
  SEP <- sep(learned, ORACLE, ifelse(mode == "dag", "pdag", "mag"), ms)
  CM <- prf1(learned, ORACLE)
  data.frame(
    method = names(outputs)[[idx]],
    shd = SHD,
    sep = SEP,
    tail_prec = 1 - mean(CM$precision[CM$which == "tail"], na.rm = TRUE),
    tail_rec = 1 - mean(CM$recall[CM$which == "tail"], na.rm = TRUE),
    tail_f1 = 1 - mean(CM$f1[CM$which == "tail"], na.rm = TRUE),
    head_prec = 1 - mean(CM$precision[CM$which == "head"], na.rm = TRUE),
    head_rec = 1 - mean(CM$recall[CM$which == "head"], na.rm = TRUE),
    head_f1 = 1 - mean(CM$f1[CM$which == "head"], na.rm = TRUE),
    time = as.difftime(timings[[idx]], units = "secs"),
    d = d, ms = ms, n = n, degree = degree, mode = mode,
    use_oracle_tests = use_oracle_tests
  )
}) |> do.call("rbind", args = _)

if (save) {
  saveRDS(res, file.path(outdir, fout))
}
