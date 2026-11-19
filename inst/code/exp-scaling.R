### Scaling of the GLIP encodings with the number of nodes under weak
### conditioning (tests with conditioning sets of size at most k).
### One task = one data set (mode, d, k, rep) from the simulation designs of
### exp-simulation.R (DAG, ADMG) and exp-chain.R (LWF chain graphs) with
### their seed formulas; tests up to k by gaussCItest; warm start PC (DAG,
### chain) or FCI (ADMG) by gaussCItest on cor(data), since the weak set of
### tests lacks the tests PC/FCI need. Only the tested statements are
### evaluated (a full oracle is infeasible for large d): k-SEP = share of the
### tested statements on which a graph and the true graph disagree.
### Stopping rule: a fit is skipped if, for the same (mode, k, encoding),
### the next smaller d of the plan has all its fits and fewer than half of
### them are certified (status OPTIMAL).
### Files: <out>/<id>/input.rds, <out>/<id>/GLIP-<enc>.rds, <out>/<id>/result.rds
rr_task_scaling <- function(task, out, plan = NULL) {
  mode <- task$mode
  d <- task$d
  ms <- min(task$ms, d - 2)
  n <- task$n
  seed <- task$rep
  alpha <- task$alpha
  dir <- file.path(out, task$id)
  if (mode == "chain") {
    suppressPackageStartupMessages({
      library("igraph")
      library("gRbase")
      library("ggm")
      library("lcd")
    })
  }
  input <- rr_cached(file.path(dir, "input.rds"), {
    if (mode == "chain") {
      degree <- min(task$degree, d - 1)
      set.seed(tseed <- 1e4 + n + seed + 1)
      gt <- create_cg_ma(n = d, d = degree, k = 2, prob = 0.8)
      data <- data.frame(lcd::rnorm.cg(n, gt, lcd::get.normal.dist(gt)))
      V <- sort(colnames(data))
      data <- data[, V]
      gt <- gt[V, V]
    } else {
      degree <- min(task$degree, d - 1)
      set.seed(tseed <- 1e4 + 3e4 * (mode == "dag") + n + seed)
      ### node names as random_graph()'s default letters, extended beyond 26
      nm <- c(letters, paste0("a", letters))
      graph <- random_graph(
        d = d, mode = mode, admg_add = task$admg_add, degree = degree,
        V = nm[seq_len(d)], V2 = nm[d + seq_len(task$admg_add)]
      )
      errMat <- mvtnorm::rmvnorm(n, sigma = diag(0.5 + stats::rchisq(ncol(graph$DAG), df = 2)))
      data <- data.frame(scale(rgraph(graph, n = n, errMat = errMat)))
      V <- colnames(data)
      gt <- switch(mode, dag = graph$DAG, admg = graph$ADMG)
    }
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
      V = V, gt = gt, tests = tests, otests = otests, pc = pcm, fci = fcim,
      t_tests = t_tests, t_pc = t_pc, t_fci = t_fci, seed = tseed
    )
  })
  if (inherits(input, "rr_error")) stop(input$error)
  V <- input$V
  if (mode == "admg") {
    start <- tryCatch(.pag_to_admg(input$fci), error = function(e) {
      .to_admg(.ess_to_dag(input$pc))
    })
    hints <- 1 * (input$fci != 0)
    baseline <- "FCI"
  } else {
    start <- .ess_to_dag(input$pc)
    hints <- 1 * (start != 0)
    baseline <- "PC"
  }
  decision <- 1 * (input$tests$p.value > alpha) # 1 = independent
  truth <- input$otests$p.value # 1 = separated in the true graph
  ### Separations of a graph on the tested statements (same order as otests)
  seps <- function(g) .compute_oracle_tests(g, ms, mode)$p.value

  fits <- list()
  for (enc in task$encodings) {
    method <- paste0("GLIP-", enc)
    file <- file.path(dir, paste0(method, ".rds"))
    if (!file.exists(file) && rr_scaling_skip(task, enc, out, plan)) next
    fits[[method]] <- rr_cached(file, {
      f <- rr_glip(enc, mode, input$tests, d, ms, V, alpha, task$walltime, start, hints,
        extra = list(SoftMemLimit = task$memlimit)
      )
      list(
        graph = f$graph, info = rr_solver_info(f), time = f$optim$runtime,
        maxmemused = if (is.null(f$optim$maxmemused)) NA else f$optim$maxmemused,
        nodecount = if (is.null(f$optim$nodecount)) NA else f$optim$nodecount
      )
    })
  }
  fits[[baseline]] <- list(
    graph = start, time = if (mode == "admg") input$t_fci else input$t_pc
  )

  s_start <- seps(start)
  res <- lapply(names(fits), \(method) {
    f <- fits[[method]]
    err <- if (inherits(f, "rr_error")) f$error else NA
    s <- if (is.na(err) && !is.null(f$graph)) seps(f$graph) else NULL
    info <- if (!is.null(f$info)) f$info else rr_solver_info(list(optim = list()))
    data.frame(
      exp = "scaling", method = method,
      ksep = if (is.null(s)) NA else mean(s != truth),
      loss = if (is.null(s)) NA else sum(s != decision),
      start_loss = sum(s_start != decision),
      input_sep = mean(decision != truth), ntests = length(truth),
      nedges = if (is.null(f$graph)) NA else rr_nedges(f$graph),
      time = if (is.null(f$time)) NA else f$time, info,
      maxmemused = if (is.null(f$maxmemused)) NA else f$maxmemused,
      nodecount = if (is.null(f$nodecount)) NA else f$nodecount,
      error = err, d = d, ms = ms, n = n, mode = mode, iter = seed,
      seed = input$seed, walltime = task$walltime, alpha = alpha,
      t_tests = input$t_tests
    )
  }) |> do.call("rbind", args = _)
  res
}

### TRUE if the fit of `enc` should be skipped by the stopping rule
rr_scaling_skip <- function(task, enc, out, plan) {
  if (is.null(plan)) {
    return(FALSE)
  }
  same <- plan[plan$mode == task$mode & plan$ms == task$ms & plan$d < task$d, ]
  if (!nrow(same)) {
    return(FALSE)
  }
  prev <- same[same$d == max(same$d), ]
  files <- file.path(out, prev$id, paste0("GLIP-", enc, ".rds"))
  if (!all(file.exists(files))) {
    ### a skipped previous fit means the setting already stopped
    skipped <- file.exists(file.path(out, prev$id, paste0("skipped-", enc)))
    if (any(skipped)) {
      rr_scaling_mark(task, enc, out)
      return(TRUE)
    }
    return(FALSE)
  }
  cert <- vapply(files, \(f) {
    x <- readRDS(f)
    !inherits(x, "rr_error") && identical(x$info$status, "OPTIMAL")
  }, logical(1))
  if (mean(cert) < 0.5) {
    rr_scaling_mark(task, enc, out)
    return(TRUE)
  }
  FALSE
}

rr_scaling_mark <- function(task, enc, out) {
  dir.create(file.path(out, task$id), recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(out, task$id, paste0("skipped-", enc)))
}

### Number of edges of a graph (matrix, or list(M1 directed, M2 bidirected));
### an undirected chain-graph edge (symmetric entries) counts once
rr_nedges <- function(g) {
  if (is.list(g)) {
    return(sum(g$M1 != 0) + sum(g$M2[upper.tri(g$M2)] != 0))
  }
  g <- g != 0
  sum(g & !t(g)) + sum((g & t(g))[upper.tri(g)])
}
