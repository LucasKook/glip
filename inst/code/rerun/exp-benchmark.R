### Experiment C (paper Sec. 5.2.2, Table 4; App. D.3, Table 6): benchmark
### data sets (discrete, n 10,000), ADMG on random node subsets of size d_max.
### One task = one (data set, d_max, subset); GLIP for ms = 1, 2, d_max - 2
### on the same tests (chi-square/mi tests up to d - 2, as in the paper),
### competitors once. Adapted from inst/code/run-benchmark-datasets.R;
### subsets: all if choose(d, d_max) <= 50, else 50 distinct random ones
### (rr_bench_subsets(), paper text), instead of nsim = 30.

rr_bench_read <- function(dataset) {
  inp <- file.path(rr_root(), "inst", "data", "datasets")
  data <- utils::read.table(file.path(inp, paste0(dataset, ".txt")), header = TRUE)
  V <- colnames(data)
  dagstr <- paste0(readLines(file.path(inp, paste0(dataset, ".graph.txt"))), collapse = ";")
  dag <- dagitty::dagitty(dagstr)
  e <- as.matrix(dagitty::edges(dag)[, 1:2])
  gt <- matrix(0, length(V), length(V), dimnames = list(V, V))
  gt[cbind(e[, 1], e[, 2])] <- 1
  list(data = data, V = V, gt = gt)
}

### Node subsets of a data set (indices into its variables), with the code
### and seed of inst/code/run-benchmark-datasets.R (set.seed(42) at the start
### of each process) and nsim = 50 (paper text) instead of 30: unique() keeps
### the order of first occurrence, so the paper's subsets come first
rr_bench_subsets <- function(dataset, d, d_max, nsim = 50) {
  set.seed(42)
  d_max <- min(d_max, d)
  nsim <- ifelse(d < d_max, 1, nsim)
  perms <- unique(lapply(seq_len(2 * nsim), \(x) sort(sample.int(d, d_max))))[1:nsim]
  perms[!vapply(perms, is.null, TRUE)]
}

rr_task_benchmark <- function(task, out) {
  mode <- "admg"
  alpha <- task$alpha
  dir <- file.path(out, task$id)
  input <- rr_cached(file.path(dir, "input.rds"), {
    B <- rr_bench_read(task$dataset)
    sub <- rr_bench_subsets(task$dataset, length(B$V), task$d_max)[[task$rep]]
    newV <- B$V[sub]
    gt <- B$gt
    if (length(B$V) > task$d_max) gt <- marginalize_dag_to_admg(gt, newV)
    data <- B$data[, newV, drop = FALSE]
    d <- length(newV)
    fdata <- data.frame(lapply(data, factor))
    t0 <- Sys.time()
    tests <- learn_graph(
      data = fdata, max_size = d - 2, mode = mode, return_tests_only = TRUE,
      comets = FALSE, test = "mi"
    )
    t_tests <- as.numeric(Sys.time() - t0, units = "secs")
    otests <- .compute_oracle_tests(gt, d - 2, mode)
    t0 <- Sys.time()
    fcim <- as(pcalg::fciPlus(list(tests = tests, V = newV), lookup_ci,
      labels = newV, alpha = alpha, selectionBias = FALSE, verbose = FALSE
    )@amat, "matrix")
    t_fci <- as.numeric(Sys.time() - t0, units = "secs")
    list(
      data = data, V = newV, gt = gt, tests = tests, otests = otests,
      fci = fcim, t_tests = t_tests, t_fci = t_fci, subset = sub
    )
  })
  if (inherits(input, "rr_error")) stop(input$error)
  V <- input$V
  d <- length(V)
  ORACLE <- .compute_graphical_representation(input$gt, d - 2, mode)
  start <- .pag_to_admg(input$fci)
  hints <- 1 * (input$fci != 0)
  input_sep <- mean(input$otests$p.value != 1 * (input$tests$p.value > alpha))

  comp <- list(
    FCI = list(graph = start, time = input$t_fci),
    NOTEARS = rr_cached(file.path(dir, "NOTEARS.rds"), rr_notears(input$data, V)),
    R2SORT = rr_cached(file.path(dir, "R2SORT.rds"), rr_r2sort(input$data, V))
  )
  mss <- unique(pmin(c(1, 2, task$d_max - 2), d - 2)) # paper: k = 1, 2, d - 2
  lapply(mss, \(ms) {
    keep <- input$tests$size <= ms
    tests_ms <- input$tests[keep, ]
    otests_ms <- input$otests[input$otests$size <= ms, ]
    fits <- comp
    for (enc in task$encodings) {
      fits[[paste0("GLIP-", enc)]] <- rr_cached(
        file.path(dir, paste0("GLIP-", enc, "-ms", ms, ".rds")),
        {
          f <- rr_glip(enc, mode, tests_ms, d, ms, V, alpha, task$walltime, start, hints)
          list(
            graph = f$graph, dcon = f$tests$dcon, info = rr_solver_info(f),
            time = f$optim$runtime
          )
        }
      )
    }
    lapply(names(fits), \(method) {
      f <- fits[[method]]
      err <- if (inherits(f, "rr_error")) f$error else NA
      learned <- NULL
      if (is.na(err) && !is.null(f$graph)) {
        learned <- .compute_graphical_representation(f$graph, d - 2, mode)
      }
      pp <- if (grepl("^GLIP", method) && !is.null(f$dcon)) 1 - f$dcon else NULL
      met <- rr_metrics(learned, ORACLE, "mag", ms, otests_ms, precomputed_predicted = pp)
      ### k-SEP of the learned graph for k = 1, 2 and d - 2 (full SEP), the
      ### same for every method (Tables 4 and 6 report 1-SEP, 2-SEP, SEP)
      ksep <- function(k) {
        if (is.null(learned) || k > d - 2) return(NA)
        1 - sep(learned, ORACLE, "mag", k,
          precomputed_groundtruth = input$otests$p.value[input$otests$size <= k]
        )$acc
      }
      full_sep <- ksep(d - 2)
      info <- if (!is.null(f$info)) f$info else rr_solver_info(list(optim = list()))
      data.frame(
        exp = "benchmark", method = method, met, full_sep = full_sep,
        sep_k1 = ksep(1), sep_k2 = ksep(2),
        input_sep = input_sep, time = if (is.null(f$time)) NA else f$time,
        info, error = err, dataset = task$dataset, d_max = task$d_max,
        d = d, ms = ms, mode = mode, wtype = "const", iter = task$rep,
        subset = paste(input$subset, collapse = ","),
        walltime = task$walltime, alpha = alpha, t_tests = input$t_tests
      )
    }) |> do.call("rbind", args = _)
  }) |> do.call("rbind", args = _)
}
