### Experiment B (paper Sec. 5.2, App. D.2): simulated DAGs and ADMGs.
### One task = one data set (mode, d, ms, n, rep), adapted from
### inst/code/run-simulation.R with the same graph/data generation and seed;
### changes: PC/FCI on gaussCItest from cor(data) (decision 2026-11-17),
### GLIP tests up to ms only, both GLIP encodings, one Gurobi thread.
### Files: <out>/<id>/input.rds, <out>/<id>/<method>.rds, <out>/<id>/result.rds
rr_task_simulation <- function(task, out) {
  mode <- task$mode
  d <- task$d
  ms <- min(task$ms, d - 2)
  n <- task$n
  seed <- task$rep
  alpha <- task$alpha
  degree <- min(task$degree, d - 1)
  dir <- file.path(out, task$id)
  input <- rr_cached(file.path(dir, "input.rds"), {
    set.seed(tseed <- 1e4 + 3e4 * (mode == "dag") + n + seed)
    graph <- random_graph(
      d = d, mode = mode, admg_add = task$admg_add, degree = degree
    )
    errMat <- mvtnorm::rmvnorm(n, sigma = diag(0.5 + stats::rchisq(ncol(graph$DAG), df = 2)))
    data <- data.frame(scale(rgraph(graph, n = n, errMat = errMat)))
    V <- colnames(data)
    gt <- switch(mode, dag = graph$DAG, admg = graph$ADMG)
    t0 <- Sys.time()
    tests <- learn_graph(
      data = data, max_size = ms, mode = mode,
      test_args = list(reg_YonZ = "lrm", reg_XonZ = "lrm"),
      return_tests_only = TRUE, comets = FALSE, test = "gaussCItest"
    )
    t_tests <- as.numeric(Sys.time() - t0, units = "secs")
    otests <- .compute_oracle_tests(gt, ms, mode)
    suff <- list(C = stats::cor(data), n = n)
    t0 <- Sys.time()
    pcm <- as(pcalg::pc(suff, pcalg::gaussCItest, labels = V, alpha = alpha)@graph, "matrix")
    t_pc <- as.numeric(Sys.time() - t0, units = "secs")
    fcim <- NULL
    t_fci <- NA
    if (mode == "admg") {
      t0 <- Sys.time()
      fcim <- as(pcalg::fciPlus(suff, pcalg::gaussCItest,
        labels = V, alpha = alpha, selectionBias = FALSE, verbose = FALSE
      )@amat, "matrix")
      t_fci <- as.numeric(Sys.time() - t0, units = "secs")
    }
    list(
      data = data, V = V, gt = gt, tests = tests, otests = otests,
      pc = pcm, fci = fcim, t_tests = t_tests, t_pc = t_pc, t_fci = t_fci,
      seed = tseed
    )
  })
  if (inherits(input, "rr_error")) stop(input$error)
  V <- input$V
  ORACLE <- .compute_graphical_representation(input$gt, d - 2, mode)
  if (mode == "dag") {
    start <- .ess_to_dag(input$pc)
    hints <- 1 * (input$pc != 0)
  } else {
    start <- tryCatch(.pag_to_admg(input$fci), error = function(e) {
      .to_admg(.ess_to_dag(input$pc)) # as in the paper driver: PC if FCI has cycles
    })
    hints <- 1 * (input$fci != 0)
  }
  input_sep <- mean(input$otests$p.value != 1 * (input$tests$p.value > alpha))

  ### Methods
  fits <- list()
  for (enc in task$encodings) {
    fits[[paste0("GLIP-", enc)]] <- rr_cached(file.path(dir, paste0("GLIP-", enc, ".rds")), {
      f <- rr_glip(enc, mode, input$tests, d, ms, V, alpha, task$walltime, start, hints)
      list(
        graph = f$graph, dcon = f$tests$dcon, info = rr_solver_info(f),
        time = f$optim$runtime
      )
    })
  }
  if (mode == "dag") {
    fits$PC <- list(graph = input$pc, time = input$t_pc)
    fits$HPC <- list(graph = .ess_to_dag(input$pc), time = input$t_pc)
  } else {
    fits$FCI <- list(graph = start, time = input$t_fci)
  }
  fits$NOTEARS <- rr_cached(file.path(dir, "NOTEARS.rds"), rr_notears(input$data, V))
  fits$R2SORT <- rr_cached(file.path(dir, "R2SORT.rds"), rr_r2sort(input$data, V))

  ### Evaluation as in the paper driver
  res <- lapply(names(fits), \(method) {
    f <- fits[[method]]
    err <- if (inherits(f, "rr_error")) f$error else NA
    learned <- NULL
    if (is.na(err) && !is.null(f$graph)) {
      if (method == "PC") {
        ### raw CPDAG with the class of the essential graph (paper driver)
        hpc <- .compute_graphical_representation(.ess_to_dag(f$graph), d - 2, "dag")
        learned <- f$graph
        class(learned) <- class(hpc)
      } else {
        learned <- .compute_graphical_representation(
          f$graph, d - 2, if (method == "HPC") "dag" else mode
        )
      }
    }
    pp <- if (grepl("^GLIP", method) && !is.null(f$dcon)) 1 - f$dcon else NULL
    met <- rr_metrics(learned, ORACLE, ifelse(mode == "dag", "pdag", "mag"), ms,
      input$otests,
      precomputed_predicted = pp
    )
    info <- if (!is.null(f$info)) f$info else rr_solver_info(list(optim = list()))
    data.frame(
      exp = "simulation", method = method, met, input_sep = input_sep,
      time = if (is.null(f$time)) NA else f$time, info, error = err,
      d = d, ms = ms, n = n, degree = degree, mode = mode, wtype = "const",
      iter = seed, seed = input$seed, walltime = task$walltime, alpha = alpha,
      t_tests = input$t_tests
    )
  }) |> do.call("rbind", args = _)
  res
}
