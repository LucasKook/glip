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
mode <- darg(args[1], "admg")
d <- as.numeric(darg(args[2], 3))
ms <- as.numeric(darg(args[3], d - 2))
ms <- ifelse(ms == -1, d - 2, ms)
ms <- ifelse(ms > d - 2, d - 2, ms)
degree <- as.numeric(darg(args[4], 3))
n <- as.numeric(darg(args[5], 1e3))
nsim <- as.numeric(darg(args[6], 1))
alpha <- as.numeric(darg(args[7], 0.01))
use_oracle_tests <- as.numeric(darg(args[8], 0))
sim_name <- darg(args[9], "test-run")
wtype <- darg(args[10], "log")
reg <- darg(args[11], "lrm")
walltime <- as.numeric(darg(args[12], 1800))
admg_add <- as.numeric(darg(args[13], 3))
save <- TRUE

# Parameters for running the optimization
ncores <- max(7, parallel::detectCores(logical = TRUE) - 2)
cache <- TRUE

# Parameters for running the tests
targs <- list(reg_YonZ = reg, reg_XonZ = reg)

out <- lapply(seq_len(nsim), \(seed) {
  # Output file
  outdir <- file.path(
    "inst", "results", "benchmark", Sys.Date(), sim_name
  )
  fout <- paste0(
    "res-mode_", mode, "-d_", d, "-ms_", ms, "-degree_",
    degree, "-n_", n, "-seed_", seed, "-alpha_", alpha,
    "-oracle_", use_oracle_tests, ".rds"
  )
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  ### Generate random graph and data
  cat("\nGenerating random graph and data\n")
  set.seed(tseed <- 1e4 + 3e4 * (mode == "dag") + n + seed)
  graph <- random_graph(d = d, prob = pr, mode = mode, admg_add = admg_add)
  data <- data.frame(py_data <- scale(rgraph(graph, n = n)))
  py_data <- r_to_py(py_data)$copy()
  V <- colnames(data)

  ### ORACLE
  gt <- switch(mode,
    "dag" = graph$DAG,
    "admg" = graph$ADMG
  )
  ORACLE <<- .compute_graphical_representation(gt, d - 2, mode)

  ### Run CITs
  cat("\nRunning conditional independence tests\n")
  tests <- otests <- tests_ms <- otests_ms <- .compute_oracle_tests(gt, d - 2, mode)
  if (!use_oracle_tests) {
    tests <- tests_ms <- learn_graph(
      data = data, max_size = d - 2, mode = mode, test_args = targs,
      return_tests_only = TRUE
    )
  }
  if (ms < d - 2) {
    tests_ms <- dplyr::filter(tests_ms, size <= ms)
    otests_ms <- dplyr::filter(otests_ms, size <= ms)
  }
  input_sep <- mean(otests$p.value != 1 * (tests$p.value > alpha))

  ### PC ALG (only under causal sufficiency/DAG case)
  cat("\nRunning PC\n")
  tstart <- Sys.time()
  pcres <- pcalg::pc(list(tests = tests, V = V), lookup_ci, labels = V, alpha = alpha)
  tstop <- Sys.time()
  runtime_PC <- tstop - tstart
  pcout <- as(pcres@graph, "matrix")
  PC <- pcout

  ### FCI ALG
  cat("\nRunning FCI\n")
  tstart <- Sys.time()
  fcires <- pcalg::fciPlus(list(tests = tests, V = V), lookup_ci,
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
    warmstart = if (mode == "dag") .ess_to_dag(PC) else .pag_to_admg(FCI),
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
