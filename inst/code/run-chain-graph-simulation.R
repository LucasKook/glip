### Benchmark glip for chain graphs
### LK 2025

### DEPs
devtools::load_all()
library("pcalg")
library("gRbase")
library("ggm")
library("lcd")

### PARs

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
d <- as.numeric(darg(args[1], 9))
ms <- as.numeric(darg(args[2], -1))
ms <- ifelse(ms == -1, d - 2, ms)
ms <- ifelse(ms > d - 2, d - 2, ms)
degree <- as.numeric(darg(args[3], 3))
n <- as.numeric(darg(args[4], 1e2))
nsim <- as.numeric(darg(args[5], 1))
alpha <- as.numeric(darg(args[6], 0.01))
use_oracle_tests <- as.numeric(darg(args[7], 0))
sim_name <- darg(args[8], "test-run")
wtype <- darg(args[9], "const")
reg <- darg(args[10], "lrm")
walltime <- as.numeric(darg(args[11], 30))
save <- TRUE
mode <- "chain"

use_comets <- FALSE # or TRUE and gcm/pcm test
test <- "gaussCItest"

# Parameters for running the optimization
ncores <- max(7, parallel::detectCores(logical = TRUE) - 2)
cache <- TRUE

# Parameters for running the tests
targs <- list(reg_YonZ = reg, reg_XonZ = reg)

out <- lapply(seq_len(nsim), \(seed) {
  # Output file
  outdir <- file.path(
    "inst", "results", "chain-graph-simulation", Sys.Date(), sim_name
  )
  fout <- paste0(
    "res-d_", d, "-ms_", ms, "-degree_",
    degree, "-n_", n, "-seed_", seed, "-alpha_", alpha,
    "-oracle_", use_oracle_tests, ".rds"
  )
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  ### Generate random graph and data
  cat("\nGenerating random graph and data\n")
  set.seed(tseed <- 1e4 + 3e4 * (mode == "dag") + n + seed)
  graph <- create_cg_ma(n = d, d = degree)
  data <- data.frame(rnorm.cg(n, graph, get.normal.dist(graph)))
  colnames(data) <- V <- letters[1:d]

  ### ORACLE
  gt <- graph
  ORACLE <<- .compute_graphical_representation(gt, d - 2, mode)

  ### Run CITs
  cat("\nRunning conditional independence tests\n")
  tests <- otests <- .compute_oracle_tests(gt, ms, mode)
  if (!use_oracle_tests) {
    tests <- learn_graph(
      data = data, max_size = ms, mode = mode, test_args = targs,
      return_tests_only = TRUE, comets = use_comets, test = test
    )
  }
  input_sep <- mean(otests$p.value != 1 * (tests$p.value > alpha))

  ### PC ALG (only under causal sufficiency/DAG case)
  cat("\nRunning PC\n")
  tstart <- Sys.time()
  pcres <- pcalg::pc(list(tests = tests, V = V), lookup_ci, labels = V, alpha = alpha)
  tstop <- Sys.time()
  runtime_PC <- tstop - tstart
  pcout <- as(pcres@graph, "matrix")
  PC <<- .compute_graphical_representation(pcout, d - 2, mode)

  ### GLIP
  cat("\nRunning GLIP\n")
  lG <- .get_opt(mode)(tests,
    d = d, max_size = ms,
    V = V, cache = cache,
    trafo = \(x) as.numeric(x <= alpha),
    weight_type = wtype,
    warmstart = pcout, # PC,
    edgehints = 1 * (pcout != 0), # 1 * (PC != 0),
    gurobi_args = list(
      Threads = ncores,
      TimeLimit = walltime
    ), mode = mode
  )

  GLIP <<- .compute_graphical_representation(lG$graph, d - 2, mode)
  runtime_GLIP <- as.difftime(lG$optim$runtime, units = "secs")

  outputs <- list(GLIP = GLIP, PC = PC)
  timings <- list(GLIP = runtime_GLIP, PC = runtime_PC)

  cat("\nEvaluating and summarizing results\n")
  res <- lapply(seq_along(outputs), \(idx) {
    method <- names(outputs)[[idx]]
    learned <- outputs[[idx]]
    SHD <- shd(learned, ORACLE)
    precomputed_predicted <- NULL
    if (method == "GLIP") {
      precomputed_predicted <- 1 - lG$tests$dcon
    }
    SEP <- sep(learned, ORACLE, mode, ms,
      precomputed_predicted = precomputed_predicted,
      precomputed_groundtruth = otests$p.value
    )
    CM <- prf1(learned, ORACLE)
    data.frame(
      method = method,
      shd = SHD,
      sep = 1 - SEP$acc,
      sep_prec = 1 - SEP$precision,
      sep_rec = 1 - SEP$recall,
      sep_f1 = 1 - SEP$f1,
      sep_mcc = SEP$mcc,
      sep_fdr = SEP$fdr,
      input_sep = input_sep,
      tail_prec = 1 - mean(CM$precision[CM$which == "tail"], na.rm = TRUE),
      tail_rec = 1 - mean(CM$recall[CM$which == "tail"], na.rm = TRUE),
      tail_f1 = 1 - mean(CM$f1[CM$which == "tail"], na.rm = TRUE),
      tail_fdr = mean(CM$fdr[CM$which == "tail"], na.rm = TRUE),
      tail_mcc = mean(CM$mcc[CM$which == "tail"], na.rm = TRUE),
      head_prec = 1 - mean(CM$precision[CM$which == "head"], na.rm = TRUE),
      head_rec = 1 - mean(CM$recall[CM$which == "head"], na.rm = TRUE),
      head_f1 = 1 - mean(CM$f1[CM$which == "head"], na.rm = TRUE),
      head_fdr = mean(CM$fdr[CM$which == "head"], na.rm = TRUE),
      head_mcc = mean(CM$mcc[CM$which == "head"], na.rm = TRUE),
      time = as.difftime(timings[[idx]], units = "secs"),
      d = d, ms = ms, n = n, degree = degree, mode = mode, wtype = wtype,
      use_oracle_tests = use_oracle_tests, iter = seed, seed = tseed,
      walltime = walltime, reg = reg
    )
  }) |> do.call("rbind", args = _)

  if (save) {
    cat("\nSaving results\n")
    saveRDS(res, file.path(outdir, fout))
  }
  res
})
out
