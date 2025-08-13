### Real data benchmark comparison
### LK 2025

### Dependencies
devtools::load_all()
library("pcalg")
library("tidyverse")
library("dagitty")
library("igraph")
library("reticulate")
use_condaenv("glip", required = TRUE)
utils <- import("dagma.utils", convert = TRUE)
dagma <- import("dagma.linear", convert = TRUE)
cd <- import("CausalDisco.baselines", convert = TRUE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
mode <- darg(args[1], "dag")
dataset <- darg(args[2], "asia")
ms <- as.numeric(darg(args[3], 1))
alpha <- as.numeric(darg(args[4], 0.001))
use_oracle_tests <- as.numeric(darg(args[5], 0))
wtype <- darg(args[6], "const")
walltime <- as.numeric(darg(args[6], 30))
d_max <- as.numeric(darg(args[7], 11))
d_max <- min(d_max, ifelse(mode == "dag", 11, 8))
reg <- darg(args[8], "lrm")
test <- "gcm"
save <- TRUE

### Folders
inp <- "./inst/data/datasets"
inp_data <- file.path(inp, paste0(dataset, ".txt"))
inp_graph <- file.path(inp, paste0(dataset, ".graph.txt"))

### Read data
data <- read_table(inp_data, show_col_types = FALSE)
V <- colnames(data)
d <- NCOL(data)

### Ground truth graph
dagstr <- paste0(read_lines(inp_graph), collapse = ";")
dag <- dagitty(dagstr)
gt <- as.matrix(as_adjacency_matrix(graph_from_edgelist(
  as.matrix(dagitty::edges(dag)[, 1:2])
))[V, V])

if (d > d_max) {
  d <- d_max
  ms <- ifelse(ms > d - 2, d - 2, ms)
  V <- V[1:d]
  data <- data[, V]
  gt <- marginalize_dag_to_admg(gt, V)
}

### Parameters for running the optimization
ncores <- max(7, parallel::detectCores(logical = TRUE) - 2)
ms <- ifelse(ms == -1, d - 2, ms)
ms <- ifelse(ms > d - 2, d - 2, ms)
cache <- TRUE

### Parameters for running the tests
targs <- list(reg_YonZ = reg, reg_XonZ = reg)
alldiscr <- (dataset != "sachs")
fdata <- data
if (alldiscr) {
  test <- "mi"
  fdata <- data.frame(mutate_all(fdata, factor))
}

### Output file
outdir <- file.path("inst", "results", "datasets", dataset)
fout <- paste0(
  "res-mode_", mode, "-d_", d, "-ms_", ms,
  "-alpha_", alpha, "-oracle_", use_oracle_tests,
  "-dataset_", dataset, c("-all", "-texout", "-graphs", "-timings", "-sumtab"), ".rds"
)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

py_data <- as.matrix(data)
py_data <- r_to_py(py_data)$copy()

### ORACLE
ORACLE <- .compute_graphical_representation(gt, d - 2, mode)

### Run CITs
cat("\nRunning conditional independence tests\n")
use_ms <- ifelse(alldiscr, d - 2, ms)
tests <- otests <- tests_ms <- otests_ms <- .compute_oracle_tests(gt, use_ms, mode, TRUE)
if (!use_oracle_tests) {
  tests <- tests_ms <- learn_graph(
    data = fdata, max_size = use_ms, mode = mode, test_args = targs,
    return_tests_only = TRUE, all_discrete = alldiscr, test = test
  )
}
if (ms < d - 2) {
  tests_ms <- dplyr::filter(tests_ms, size <= ms)
  otests_ms <- dplyr::filter(otests_ms, size <= ms)
}
input_sep <- mean(otests$p.value != 1 * (tests$p.value > alpha))

### PC ALG (only under causal sufficiency/DAG case)
if (!alldiscr) {
  suff <- list(C = cov(data), n = nrow(data))
  cit <- gaussCItest
} else {
  suff <- list(tests = tests, V = V)
  cit <- lookup_ci
}
cat("\nRunning PC\n")
tstart <- Sys.time()
pcres <- pcalg::pc(suff, cit, labels = V, alpha = alpha)
tstop <- Sys.time()
runtime_PC <- tstop - tstart
pcout <- as(pcres@graph, "matrix")
PC <- pcout

### FCI ALG
cat("\nRunning FCI\n")
tstart <- Sys.time()
fcires <- pcalg::fciPlus(suff, cit,
  labels = V, alpha = alpha, selectionBias = FALSE, verbose = FALSE
)
tstop <- Sys.time()
runtime_FCI <- tstop - tstart
fciout <- as(fcires@amat, "matrix")
FCI <- fciout

### R2sortability
cat("\nRunning R2SORT\n")
tstart <- Sys.time()
r2s <- 1 * (cd$r2_sort_regress(py_data) != 0)
tstop <- Sys.time()
runtime_R2SORT <- tstop - tstart
dimnames(r2s) <- list(V, V)
R2SORT <- .compute_graphical_representation(r2s, d - 2, mode)

### GLIP
cat("\nRunning GLIP\n")
lG <- .get_opt(mode)(tests_ms,
  d = d, max_size = ms,
  V = V, cache = cache,
  trafo = \(x) as.numeric(x <= alpha),
  weight_type = wtype,
  warmstart = if (mode == "dag") gt else NULL,
  edgehints = if (mode == "dag") 1 * (PC != 0) else 1 * (FCI != 0),
  gurobi_args = list(
    Threads = ncores,
    TimeLimit = walltime
  ), mode = mode
)

GLIP <- .compute_graphical_representation(lG$graph, d - 2, mode)
runtime_GLIP <- as.difftime(lG$optim$runtime, units = "secs")

### NOTEARS
cat("\nRunning NOTEARS\n")
model <- dagma$DagmaLinear(loss_type = "l2")
tstart <- Sys.time()
nto <- 1 * (model$fit(py_data, lambda1 = 0.02) != 0)
tstop <- Sys.time()
runtime_NOTEARS <- tstop - tstart
dimnames(nto) <- list(V, V)
NOTEARS <- .compute_graphical_representation(nto, d - 2, mode)

### Evaluate and summarize results
class(PC) <- class(GLIP)
class(FCI) <- class(GLIP)
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
if (mode == "dag") { # remove FCI in case of DAGs
  outputs <- outputs[-grep("FCI", names(outputs))]
  timings <- timings[-grep("FCI", names(timings))]
}

cat("\nEvaluating and summarizing results\n")
res <- lapply(seq_along(outputs), \(idx) {
  method <- names(outputs)[[idx]]
  learned <- outputs[[idx]]
  SHD <- shd(learned, ORACLE)
  precomputed_predicted <- NULL
  if (method == "GLIP") {
    precomputed_predicted <- 1 - lG$tests$dcon
  }
  SEP <- sep(learned, ORACLE, ifelse(mode == "dag", "pdag", "mag"), ms,
    precomputed_predicted = precomputed_predicted,
    precomputed_groundtruth = otests_ms$p.value
  )
  CM <- prf1(learned, ORACLE)
  data.frame(
    method = method,
    shd = SHD,
    sep = 1 - SEP$acc,
    input_sep = input_sep,
    tail_prec = 1 - mean(CM$precision[CM$which == "tail"], na.rm = TRUE),
    tail_rec = 1 - mean(CM$recall[CM$which == "tail"], na.rm = TRUE),
    tail_f1 = 1 - mean(CM$f1[CM$which == "tail"], na.rm = TRUE),
    tail_fdr = mean(CM$fdr[CM$which == "tail"], na.rm = TRUE),
    head_prec = 1 - mean(CM$precision[CM$which == "head"], na.rm = TRUE),
    head_rec = 1 - mean(CM$recall[CM$which == "head"], na.rm = TRUE),
    head_f1 = 1 - mean(CM$f1[CM$which == "head"], na.rm = TRUE),
    head_fdr = mean(CM$fdr[CM$which == "head"], na.rm = TRUE),
    time = as.difftime(timings[[idx]], units = "secs"),
    d = d, ms = ms, mode = mode, wtype = wtype,
    use_oracle_tests = use_oracle_tests
  )
}) |> do.call("rbind", args = _)

sumtab <- res |>
  group_by(method) |>
  summarize(
    shd = shd * d^2,
    sep = sep,
    prec = mean(c(tail_prec, head_prec), na.rm = TRUE),
    rec = mean(c(tail_rec, head_rec), na.rm = TRUE),
    fdr = mean(c(tail_fdr, head_fdr), na.rm = TRUE),
  ) |>
  select(method, shd, sep, prec, rec, fdr)

texout <- knitr::kable(sumtab, format = "latex", booktabs = TRUE, digits = 2)

if (save) {
  cat("\nSaving results\n")
  saveRDS(res, file.path(outdir, fout[1]))
  saveRDS(texout, file.path(outdir, fout[2]))
  saveRDS(outputs, file.path(outdir, fout[3]))
  saveRDS(timings, file.path(outdir, fout[4]))
  saveRDS(sumtab, file.path(outdir, fout[5]))
}

res
sumtab
texout
