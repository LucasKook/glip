### Experiment D (paper App. D.4): LWF chain graphs (Ma et al., 2008).
### One task = one data set (d, n, rep), adapted from
### inst/code/run-chain-graph-simulation.R with the same graph/data generation
### and seed (full conditioning, so PC uses the GLIP tests via lookup_ci as
### in the paper); both GLIP encodings, one Gurobi thread. Needs lcd.
rr_task_chain <- function(task, out) {
  mode <- "chain"
  d <- task$d
  ms <- min(task$ms, d - 2)
  n <- task$n
  seed <- task$rep
  alpha <- task$alpha
  degree <- min(task$degree, d - 1)
  dir <- file.path(out, task$id)
  ### lcd calls igraph functions unqualified: attach the packages of the
  ### paper driver
  suppressPackageStartupMessages({
    library("igraph")
    library("gRbase")
    library("ggm")
    library("lcd")
  })
  input <- rr_cached(file.path(dir, "input.rds"), {
    set.seed(tseed <- 1e4 + n + seed + 1)
    graph <- create_cg_ma(n = d, d = degree, k = 2, prob = 0.8)
    data <- data.frame(lcd::rnorm.cg(n, graph, lcd::get.normal.dist(graph)))
    V <- sort(colnames(data))
    data <- data[, V]
    graph <- graph[V, V]
    t0 <- Sys.time()
    tests <- learn_graph(
      data = data, max_size = ms, mode = mode,
      test_args = list(reg_YonZ = "lrm", reg_XonZ = "lrm"),
      return_tests_only = TRUE, comets = FALSE, test = "gaussCItest"
    )
    t_tests <- as.numeric(Sys.time() - t0, units = "secs")
    otests <- .compute_oracle_tests(graph, ms, mode)
    t0 <- Sys.time()
    pcm <- as(pcalg::pc(list(tests = tests, V = V), lookup_ci,
      labels = V, alpha = alpha, skel.method = "stable"
    )@graph, "matrix")
    t_pc <- as.numeric(Sys.time() - t0, units = "secs")
    list(
      V = V, gt = graph, tests = tests, otests = otests, pc = pcm,
      t_tests = t_tests, t_pc = t_pc, seed = tseed
    )
  })
  if (inherits(input, "rr_error")) stop(input$error)
  V <- input$V
  ORACLE <- .compute_graphical_representation(input$gt, d - 2, mode)
  start <- .ess_to_dag(input$pc)
  hints <- 1 * (start != 0)
  input_sep <- mean(input$otests$p.value != 1 * (input$tests$p.value > alpha))

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
  fits$PC <- list(graph = start, time = input$t_pc)

  res <- lapply(names(fits), \(method) {
    f <- fits[[method]]
    err <- if (inherits(f, "rr_error")) f$error else NA
    learned <- NULL
    if (is.na(err) && !is.null(f$graph)) {
      learned <- .compute_graphical_representation(f$graph, d - 2, mode)
    }
    pp <- if (grepl("^GLIP", method) && !is.null(f$dcon)) 1 - f$dcon else NULL
    met <- rr_metrics(learned, ORACLE, mode, ms, input$otests,
      precomputed_predicted = pp
    )
    info <- if (!is.null(f$info)) f$info else rr_solver_info(list(optim = list()))
    data.frame(
      exp = "chain", method = method, met, input_sep = input_sep,
      time = if (is.null(f$time)) NA else f$time, info, error = err,
      d = d, ms = ms, n = n, degree = degree, mode = mode, wtype = "const",
      iter = seed, seed = input$seed, walltime = task$walltime, alpha = alpha,
      t_tests = input$t_tests
    )
  }) |> do.call("rbind", args = _)
  res
}
