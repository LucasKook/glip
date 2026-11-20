### Ablation of the cut families and solver logs of both encodings.
### One task = one existing data set (tests and PC/FCI warm start) from the
### simulation or scaling experiments; fits: the min-length encoding and the
### cut encoding with its families added step by step:
###   a: closure inequalities and dependence cuts (lazy at integer solutions,
###      one static cut per dependent decision) only
###   b: + minimum cuts at the root
###   c: + joint and lifted cuts
###   d: + subset consistency cuts (= default)
### Every fit writes a Gurobi log (<out>/<id>/<fit>.log) for the root bound
### and the incumbent and bound over time.
### Environment: ABL_RERUN (run folder of the simulation and chain
### experiments), ABL_SCALING (run folder of the scaling experiment).
ABL_CONFIGS <- list(
  a = list(nodecuts = "none", joint = FALSE, subset = "none"),
  b = list(nodecuts = "root", joint = FALSE, subset = "none"),
  c = list(nodecuts = "root", joint = TRUE, subset = "none"),
  d = list(nodecuts = "root", joint = TRUE, subset = "lp")
)

rr_ablation_input <- function(task) {
  rr <- Sys.getenv("ABL_RERUN", file.path(rr_root(), "inst/results/rerun-2026-11/raw/rerun-2026-11"))
  sc <- Sys.getenv("ABL_SCALING", file.path(rr_root(), "inst/results/scaling-2026-11/raw/scaling-2026-11"))
  m <- task$mode
  if (m == "chain") {
    if (task$k == "full") return(file.path(rr, "chain/full", sprintf("d%02d-n1000-r%02d", task$d, task$rep)))
    return(file.path(sc, "scaling", paste0("k", task$k), sprintf("chain-d%02d-r%02d", task$d, task$rep)))
  }
  if (task$k == "full") {
    return(file.path(rr, "simulation/full", sprintf("%s-d%02d-ms%d-n400-r%02d", m, task$d, task$d - 2, task$rep)))
  }
  file.path(rr, "simulation/weak", sprintf("%s-d%02d-ms%s-n400-r%02d", m, task$d, task$k, task$rep))
}

rr_task_ablation <- function(task, out) {
  dir <- file.path(out, task$id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  input <- readRDS(file.path(rr_ablation_input(task), "input.rds"))
  mode <- task$mode
  V <- input$V
  d <- length(V)
  ms <- if (task$k == "full") d - 2 else as.integer(task$k)
  alpha <- task$alpha
  if (mode == "admg") {
    start <- tryCatch(.pag_to_admg(input$fci), error = function(e) .to_admg(.ess_to_dag(input$pc)))
    hints <- 1 * (input$fci != 0)
  } else {
    start <- .ess_to_dag(input$pc)
    hints <- 1 * (start != 0)
  }
  decision <- 1 * (input$tests$p.value > alpha)
  loss_of <- function(g) sum(.compute_oracle_tests(g, ms, mode)$p.value != decision)
  start_loss <- loss_of(start)
  fits <- c("minlength", paste0("cuts-", names(ABL_CONFIGS)))
  res <- lapply(fits, \(fit) {
    log <- file.path(dir, paste0(fit, ".log"))
    f <- rr_cached(file.path(dir, paste0(fit, ".rds")), {
      if (file.exists(log)) file.remove(log)
      if (fit == "minlength") {
        g <- rr_glip("minlength", mode, input$tests, d, ms, V, alpha, task$walltime, start, hints,
          extra = list(LogFile = log, LogToConsole = 0)
        )
      } else {
        cfg <- ABL_CONFIGS[[sub("cuts-", "", fit)]]
        g <- rr_glip_cuts_opts(mode, input$tests, d, ms, V, alpha, task$walltime, start,
          c(cfg, list(gurobi_params = list(OutputFlag = 1, LogToConsole = 0, LogFile = log)))
        )
      }
      o <- g$optim
      list(graph = g$graph, info = rr_solver_info(g),
           nodecount = if (is.null(o$nodecount)) NA else o$nodecount,
           cuts = if (is.null(o$cuts)) NULL else unlist(o$cuts))
    })
    if (inherits(f, "rr_error")) {
      return(data.frame(fit = fit, error = f$error))
    }
    cuts <- f$cuts
    cnt <- function(nm) if (!is.null(cuts) && nm %in% names(cuts)) cuts[[nm]] else NA
    data.frame(
      fit = fit, f$info, nodecount = f$nodecount,
      loss = if (is.null(f$graph)) NA else loss_of(f$graph),
      lazy = cnt("lazy"), user = cnt("user"), subset_lp = cnt("subset_lp"), cycle = cnt("cycle"),
      error = NA
    )
  })
  res <- dplyr::bind_rows(res)
  cbind(data.frame(exp = "ablation", setting = task$setting, mode = mode, k = task$k, d = d,
                   ms = ms, rep = task$rep, ntests = length(decision), start_loss = start_loss,
                   walltime = task$walltime), res)
}

### cut encoding with explicit cut options (rr_glip passes none)
rr_glip_cuts_opts <- function(mode, tests, d, ms, V, alpha, walltime, warmstart, cut_options) {
  cl <- rr_class(mode)
  args <- list(tests,
    d = d, max_size = ms, V = V, class = cl$cuts,
    trafo = function(x) as.numeric(x <= alpha), weight_type = "const",
    warmstart = warmstart, gurobi_args = rr_gurobi_args(walltime), cut_options = cut_options
  )
  py <- Sys.getenv("GLIP_PYTHON")
  if (nzchar(py)) args$python <- py
  t0 <- Sys.time()
  fit <- do.call(cuts_optim, args)
  fit$wall <- as.numeric(Sys.time() - t0, units = "secs")
  fit
}

rr_plan_ablation <- function() {
  s <- rbind(
    expand.grid(mode = c("dag", "admg"), k = c("1", "2", "full"), stringsAsFactors = FALSE),
    data.frame(mode = c("dag", "admg"), k = "1"),
    data.frame(mode = "chain", k = c("full", "1", "2"))
  )
  s$d <- ifelse(s$mode == "chain", 6, ifelse(s$k == "full", 7, 10))
  s$d[nrow(s) - 3 - (1:0)] <- 12 # the added DAG/ADMG k = 1 rows: d = 12
  s$reps <- ifelse(s$mode == "chain" & s$k != "full", 10, 20)
  tasks <- do.call(rbind, lapply(seq_len(nrow(s)), \(i) data.frame(s[i, c("mode", "k", "d")], rep = seq_len(s$reps[i]))))
  tasks <- transform(tasks,
    exp = "ablation", setting = paste0(mode, "-k", k, "-d", d), alpha = .001,
    walltime = 300, nfits = 5, n = ifelse(mode == "chain" & k == "full", 1000, 400),
    ms = ifelse(k == "full", d - 2, suppressWarnings(as.integer(k))),
    id = sprintf("ablation/%s-k%s-d%02d-r%02d", mode, k, d, rep)
  )
  tasks$cost <- tasks$walltime * tasks$nfits
  tasks <- tasks[order(tasks$rep, -tasks$d), ]
  tasks$order <- seq_len(nrow(tasks))
  tasks
}
