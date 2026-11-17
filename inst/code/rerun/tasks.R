### Rerun of the paper experiments (arXiv:2601.20589 v2) as resumable tasks.
### From the repository root:
###   Rscript inst/code/rerun/tasks.R plan <out>
###       writes <out>/tasks.csv (all experiments of the paper)
###   Rscript inst/code/rerun/tasks.R run <out> <chunk> <n_chunks> [exp]
###       runs every n_chunks-th task (longest first), skipping finished ones
###       (<out>/<id>/result.rds); several chunks run in parallel, e.g.
###       seq 1 32 | xargs -P 32 -I{} Rscript .../tasks.R run <out> {} 32
###   Rscript inst/code/rerun/tasks.R prebuild <out>
###       builds the caches of the min-length encodings for all tasks, one
###       after another with the solver stubbed (before parallel runs)
###   Rscript inst/code/rerun/tasks.R collect <out>
###       one csv per experiment in <out>/collected/
###   Rscript inst/code/rerun/tasks.R smoke <out> [walltime]
###       one small task per experiment with a short time limit
### See common.R for the environment variables.
args <- commandArgs(trailingOnly = TRUE)
cmd <- args[1]
out <- args[2]
here <- file.path(Sys.getenv("GLIP_ROOT", getwd()), "inst", "code", "rerun")
for (f in c("common.R", "exp-simulation.R", "exp-chain.R", "exp-benchmark.R", "exp-asp.R")) {
  source(file.path(here, f))
}

ENCODINGS <- c("minlength", "cuts")

rr_plan <- function() {
  g <- function(...) expand.grid(..., stringsAsFactors = FALSE)
  pad <- function(r) sprintf("r%02d", r)
  ### A: ASP timing (paper: DAG and ADMG, full d 3 to 9, weak k = 1 d 3 to 14)
  asp <- rbind(
    transform(g(mode = c("dag", "admg"), d = 3:9, rep = 1:20), setting = "full", ms = d - 2),
    transform(g(mode = c("dag", "admg"), d = 3:14, rep = 1:20), setting = "weak", ms = 1)
  )
  asp <- transform(asp,
    exp = "asp", n = 1e4, alpha = .001, walltime = 600, nfits = 3,
    id = sprintf("asp/%s/%s-d%02d-ms%d-%s", setting, mode, d, ms, pad(rep))
  )
  ### B: simulations (full d 3 to 9; weak k = 1, 2 for d 3 to 12; n 400, 1e4)
  sim <- rbind(
    transform(g(mode = c("dag", "admg"), d = 3:9, n = c(400, 1e4), rep = 1:30),
      setting = "full", ms = d - 2
    ),
    transform(g(mode = c("dag", "admg"), d = 3:12, ms = 1:2, n = c(400, 1e4), rep = 1:30),
      setting = "weak"
    )
  )
  sim <- transform(sim,
    exp = "simulation", alpha = .001, walltime = 300, nfits = 2,
    degree = 2, admg_add = 3,
    id = sprintf("simulation/%s/%s-d%02d-ms%d-n%d-%s", setting, mode, d, ms, n, pad(rep))
  )
  ### C: benchmarks (ADMG, d_max 6 and 8, k = 1, 2, d - 2 in each task)
  bench <- do.call(rbind, lapply(c("alarm", "asia", "child", "hepar2", "sachs"), \(ds) {
    B <- rr_bench_read(ds)
    do.call(rbind, lapply(c(6, 8), \(dm) {
      ns <- length(rr_bench_subsets(ds, length(B$V), dm))
      data.frame(dataset = ds, d_max = dm, rep = seq_len(ns), d = min(dm, length(B$V)))
    }))
  }))
  bench <- transform(bench,
    exp = "benchmark", setting = paste0("d", d_max), mode = "admg", n = 1e4,
    alpha = .001, walltime = 300, nfits = 6, ms = NA,
    id = sprintf("benchmark/d%d/%s-%s", d_max, dataset, pad(rep))
  )
  ### D: chain graphs (full conditioning, d 3 to 10, n 1,000 and 1e5)
  chain <- transform(g(d = 3:10, n = c(1000, 1e5), rep = 1:30),
    exp = "chain", setting = "full", mode = "chain", ms = d - 2,
    alpha = .001, walltime = 300, nfits = 2, degree = 4,
    id = sprintf("chain/full/d%02d-n%d-%s", d, n, pad(rep))
  )
  cols <- c(
    "id", "exp", "setting", "mode", "d", "ms", "n", "rep", "alpha", "walltime",
    "nfits", "degree", "admg_add", "dataset", "d_max"
  )
  fill <- function(x) {
    for (k in setdiff(cols, names(x))) x[[k]] <- NA
    x[, cols]
  }
  tasks <- rbind(fill(asp), fill(sim), fill(bench), fill(chain))
  ### Upper bound of the solver time, used to start long tasks first
  tasks$cost <- tasks$walltime * tasks$nfits * (tasks$d >= 5)
  tasks <- tasks[order(-tasks$cost, -tasks$d, tasks$id), ]
  tasks$order <- seq_len(nrow(tasks))
  tasks
}

rr_run_task <- function(task, out) {
  task$encodings <- ENCODINGS
  f <- switch(task$exp,
    asp = rr_task_asp, simulation = rr_task_simulation,
    benchmark = rr_task_benchmark, chain = rr_task_chain
  )
  res <- f(task, out)
  attr(res, "versions") <- rr_versions()
  rr_save(res, file.path(out, task$id, "result.rds"))
  res
}

rr_run <- function(out, chunk, n_chunks, exp = NULL) {
  tasks <- utils::read.csv(file.path(out, "tasks.csv"), stringsAsFactors = FALSE)
  if (!is.null(exp)) tasks <- tasks[tasks$exp %in% exp, ]
  mine <- tasks[(seq_len(nrow(tasks)) - 1) %% n_chunks == chunk - 1, ]
  rr_setup()
  for (k in seq_len(nrow(mine))) {
    task <- as.list(mine[k, ])
    done <- file.path(out, task$id, "result.rds")
    if (file.exists(done)) next
    cat(format(Sys.time()), "start", task$id, "\n")
    t0 <- Sys.time()
    ok <- tryCatch(
      {
        rr_run_task(task, out)
        unlink(file.path(out, task$id, "error.txt"))
        TRUE
      },
      error = function(e) {
        dir.create(file.path(out, task$id), recursive = TRUE, showWarnings = FALSE)
        writeLines(conditionMessage(e), file.path(out, task$id, "error.txt"))
        FALSE
      }
    )
    cat(format(Sys.time()), if (ok) "done" else "ERROR", task$id,
      round(as.numeric(Sys.time() - t0, units = "secs")), "s\n")
  }
}

rr_collect <- function(out) {
  tasks <- utils::read.csv(file.path(out, "tasks.csv"), stringsAsFactors = FALSE)
  dir.create(file.path(out, "collected"), showWarnings = FALSE)
  status <- data.frame(
    id = tasks$id, exp = tasks$exp,
    done = file.exists(file.path(out, tasks$id, "result.rds")),
    error = file.exists(file.path(out, tasks$id, "error.txt"))
  )
  utils::write.csv(status, file.path(out, "collected", "status.csv"), row.names = FALSE)
  for (e in unique(tasks$exp)) {
    ids <- status$id[status$exp == e & status$done]
    if (!length(ids)) next
    rs <- lapply(ids, \(i) {
      r <- readRDS(file.path(out, i, "result.rds"))
      data.frame(task = i, r)
    })
    cols <- unique(unlist(lapply(rs, names)))
    rs <- lapply(rs, \(r) {
      for (k in setdiff(cols, names(r))) r[[k]] <- NA
      r[, cols]
    })
    utils::write.csv(do.call(rbind, rs), file.path(out, "collected", paste0(e, ".csv")),
      row.names = FALSE
    )
    v <- attr(readRDS(file.path(out, ids[1], "result.rds")), "versions")
    saveRDS(v, file.path(out, "collected", paste0(e, "-versions.rds")))
  }
  print(stats::aggregate(cbind(n = 1, done, error) ~ exp, status, sum))
}

### Min-length caches depend on (optimizer, mode, d, max_size) only
rr_prebuild <- function(out) {
  tasks <- utils::read.csv(file.path(out, "tasks.csv"), stringsAsFactors = FALSE)
  keys <- do.call(rbind, lapply(seq_len(nrow(tasks)), \(k) {
    t <- tasks[k, ]
    mss <- if (t$exp == "benchmark") unique(pmin(c(1, 2, t$d_max - 2), t$d - 2)) else min(t$ms, t$d - 2)
    data.frame(mode = t$mode, d = t$d, ms = mss)
  }))
  keys <- unique(keys)
  keys <- keys[order(keys$mode, keys$d, keys$ms), ]
  rr_setup()
  local_stub <- function(model, gurobi_args) list(status = "STUB")
  ns <- asNamespace("glip")
  unlockBinding(".gurobi_solve", ns)
  assign(".gurobi_solve", local_stub, envir = ns)
  for (k in seq_len(nrow(keys))) {
    x <- keys[k, ]
    V <- paste0("X", seq_len(x$d))
    sets <- .list_tests_graph(V, x$ms)$sets
    tests <- data.frame(
      X = vapply(sets, `[[`, "", "X"), Y = vapply(sets, `[[`, "", "Y"),
      Z = vapply(sets, \(q) paste(q$Z, collapse = ","), ""),
      size = lengths(lapply(sets, `[[`, "Z")), p.value = .5
    )
    t0 <- Sys.time()
    cl <- rr_class(x$mode)
    args <- list(tests,
      d = x$d, max_size = x$ms, V = V, cache = TRUE,
      cache_dir = file.path(Sys.getenv("GLIP_CACHE", file.path(rr_root(), ".cache-rerun")), cl$fun)
    )
    if (!is.null(cl$mode)) args$mode <- cl$mode
    invisible(utils::capture.output(do.call(cl$fun, args)))
    cat(sprintf("%s d=%d ms=%d %.1f s\n", x$mode, x$d, x$ms, as.numeric(Sys.time() - t0, units = "secs")))
  }
}

if (identical(cmd, "plan")) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  rr_setup()
  tasks <- rr_plan()
  utils::write.csv(tasks, file.path(out, "tasks.csv"), row.names = FALSE)
  print(table(tasks$exp, tasks$setting))
  cat("upper bound of the solver time:", round(sum(tasks$cost) / 3600), "h\n")
} else if (identical(cmd, "run")) {
  rr_run(out, as.integer(args[3]), as.integer(args[4]),
    if (length(args) > 4) strsplit(args[5], ",")[[1]] else NULL
  )
} else if (identical(cmd, "prebuild")) {
  rr_prebuild(out)
} else if (identical(cmd, "collect")) {
  rr_collect(out)
} else if (identical(cmd, "smoke")) {
  wt <- if (length(args) > 2) as.numeric(args[3]) else 10
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  rr_setup()
  tasks <- rr_plan()
  pick <- c(
    "asp/full/admg-d05-ms3-r01", "asp/weak/dag-d06-ms1-r01",
    "simulation/full/dag-d05-ms3-n400-r01", "simulation/weak/admg-d06-ms1-n10000-r01",
    "benchmark/d6/asia-r01", "chain/full/d05-n1000-r01"
  )
  only <- Sys.getenv("RR_SMOKE_EXP")
  tasks <- tasks[tasks$id %in% pick, ]
  if (nzchar(only)) tasks <- tasks[tasks$exp %in% strsplit(only, ",")[[1]], ]
  tasks$walltime <- wt
  utils::write.csv(tasks, file.path(out, "tasks.csv"), row.names = FALSE)
  rr_run(out, 1, 1)
  rr_collect(out)
}
