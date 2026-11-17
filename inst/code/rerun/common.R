### Shared helpers of the rerun of the paper experiments (inst/code/rerun)
###
### Environment (all optional, defaults for a local run from the repository
### root):
###   GLIP_ROOT        repository root (default: working directory)
###   GLIP_PYTHON      Python with gurobipy for the cut encoding
###   GLIP_COMP_VENV   virtualenv with dagma and CausalDisco (NOTEARS, R2SORT)
###   GLIP_CACHE       cache directory of the min-length encodings
###   GLIP_ASP_HOME    unzipped inst/asp/hyttinen2014uai_ver6.zip
###                    (default: <scratch>/asp-home, unzipped on first use)
### One Gurobi thread per solve (OMP/BLAS threads 1 as well); parallelism
### comes from running several tasks at once (tasks.R run).

rr_root <- function() normalizePath(Sys.getenv("GLIP_ROOT", getwd()))

rr_setup <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1,
    VECLIB_MAXIMUM_THREADS = 1, MKL_NUM_THREADS = 1
  )
  suppressPackageStartupMessages({
    pkgload::load_all(rr_root(), quiet = TRUE, export_all = TRUE)
    library("pcalg")
  })
  invisible(TRUE)
}

### Gurobi arguments of every GLIP fit
rr_gurobi_args <- function(walltime, extra = list()) {
  c(list(Threads = 1, TimeLimit = walltime), extra)
}

### Class of the cut encoding and min-length optimizer for a paper mode
rr_class <- function(mode) {
  switch(mode,
    dag = list(cuts = "dag", fun = "dcon_lean_optim", mode = "dag"),
    admg = list(cuts = "admg", fun = "admg_lean_optim", mode = "admg-lean"),
    chain = list(cuts = "chain", fun = "chain_lean_optim", mode = NULL),
    stop("unknown mode ", mode)
  )
}

### Fit GLIP with the min-length or the cut encoding (same tests, decisions,
### weights, warm start and edge hints). Returns the graphopt object (graph
### NULL without incumbent) with the wall time of the call
rr_glip <- function(encoding, mode, tests, d, ms, V, alpha, walltime,
                    warmstart, edgehints, wtype = "const", extra = list()) {
  cl <- rr_class(mode)
  args <- list(tests,
    d = d, max_size = ms, V = V,
    trafo = function(x) as.numeric(x <= alpha), weight_type = wtype,
    warmstart = warmstart, edgehints = edgehints,
    gurobi_args = rr_gurobi_args(walltime, extra)
  )
  t0 <- Sys.time()
  if (encoding == "minlength") {
    args$cache <- TRUE
    args$cache_dir <- file.path(
      Sys.getenv("GLIP_CACHE", file.path(rr_root(), ".cache-rerun")),
      cl$fun
    )
    if (!is.null(cl$mode)) args$mode <- cl$mode
    fit <- suppressWarnings(do.call(cl$fun, args))
  } else if (encoding == "cuts") {
    args$class <- cl$cuts
    args$edgehints <- NULL # hints are not used by the cut encoding
    py <- Sys.getenv("GLIP_PYTHON")
    if (nzchar(py)) args$python <- py
    fit <- do.call(cuts_optim, args)
  } else {
    stop("unknown encoding ", encoding)
  }
  fit$wall <- as.numeric(Sys.time() - t0, units = "secs")
  fit
}

### Solver summary of a GLIP fit
rr_solver_info <- function(fit) {
  o <- fit$optim
  g <- function(x) if (is.null(x) || length(x) == 0) NA else x[[1]]
  data.frame(
    status = g(o$status), objval = g(o$objval), objbound = g(o$objbound),
    mipgap = g(o$mipgap), runtime = g(o$runtime), wall = g(fit$wall)
  )
}

### Deterministic seed from task fields
rr_seed <- function(...) {
  s <- utf8ToInt(paste(..., sep = "|"))
  as.integer(sum(as.numeric(s) * seq_along(s) * 7919) %% 2147483647)
}

### Python competitors (dagma, CausalDisco) from GLIP_COMP_VENV
rr_python_competitors <- function() {
  if (!is.null(.rr$py)) {
    return(.rr$py)
  }
  venv <- Sys.getenv("GLIP_COMP_VENV")
  if (nzchar(venv)) reticulate::use_virtualenv(venv, required = TRUE)
  .rr$py <- list(
    dagma = reticulate::import("dagma.linear", convert = TRUE),
    cd = reticulate::import("CausalDisco.baselines", convert = TRUE)
  )
  .rr$py
}
.rr <- new.env()

rr_notears <- function(data, V) {
  py <- rr_python_competitors()
  X <- as.matrix(data)
  storage.mode(X) <- "double" # integer (discrete) data would become int32
  X <- reticulate::r_to_py(X)$copy()
  t0 <- Sys.time()
  model <- py$dagma$DagmaLinear(loss_type = "l2")
  G <- 1 * (model$fit(X, lambda1 = 0.02) != 0)
  dimnames(G) <- list(V, V)
  list(graph = G, time = as.numeric(Sys.time() - t0, units = "secs"))
}

rr_r2sort <- function(data, V) {
  py <- rr_python_competitors()
  X <- as.matrix(data)
  storage.mode(X) <- "double" # integer (discrete) data would become int32
  X <- reticulate::r_to_py(X)$copy()
  t0 <- Sys.time()
  G <- 1 * (py$cd$r2_sort_regress(X) != 0)
  dimnames(G) <- list(V, V)
  list(graph = G, time = as.numeric(Sys.time() - t0, units = "secs"))
}

### Metrics of the paper (as in inst/code/run-simulation.R) of a learned
### graph against the oracle; NA if the method returned no graph
rr_metrics <- function(learned, ORACLE, sepmode, ms, otests_ms,
                       precomputed_predicted = NULL) {
  na <- data.frame(
    shd = NA, sep = NA, sep_prec = NA, sep_rec = NA, sep_f1 = NA,
    sep_mcc = NA, sep_fdr = NA, tail_prec = NA, tail_rec = NA, tail_f1 = NA,
    tail_fdr = NA, tail_mcc = NA, head_prec = NA, head_rec = NA,
    head_f1 = NA, head_fdr = NA, head_mcc = NA
  )
  if (is.null(learned)) {
    return(na)
  }
  SEP <- sep(learned, ORACLE, sepmode, ms,
    precomputed_predicted = precomputed_predicted,
    precomputed_groundtruth = otests_ms$p.value
  )
  CM <- prf1(learned, ORACLE)
  m <- function(f, w) mean(CM[[f]][CM$which == w], na.rm = TRUE)
  data.frame(
    shd = shd(learned, ORACLE),
    sep = 1 - SEP$acc, sep_prec = 1 - SEP$precision,
    sep_rec = 1 - SEP$recall, sep_f1 = 1 - SEP$f1,
    sep_mcc = SEP$mcc, sep_fdr = SEP$fdr,
    tail_prec = 1 - m("precision", "tail"), tail_rec = 1 - m("recall", "tail"),
    tail_f1 = 1 - m("f1", "tail"), tail_fdr = m("fdr", "tail"),
    tail_mcc = m("mcc", "tail"),
    head_prec = 1 - m("precision", "head"), head_rec = 1 - m("recall", "head"),
    head_f1 = 1 - m("f1", "head"), head_fdr = m("fdr", "head"),
    head_mcc = m("mcc", "head")
  )
}

### Versions of R, packages, Gurobi, gurobipy, Python libraries and clingo
rr_versions <- function() {
  if (!is.null(.rr$versions)) {
    return(c(.rr$versions, clingo = .rr$clingo))
  }
  pv <- function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA)
  }
  pyv <- function(py, mods) {
    if (!nzchar(py) || !file.exists(py)) {
      return(NA)
    }
    code <- paste0(
      "import importlib.metadata as m\n",
      "for p in [", paste0("'", mods, "'", collapse = ", "), "]:\n",
      "    try: print(p, m.version(p))\n",
      "    except Exception: print(p, 'NA')"
    )
    out <- tryCatch(system2(py, c("-c", shQuote(code)), stdout = TRUE),
      error = function(e) NA, warning = function(w) NA
    )
    paste(out, collapse = "; ")
  }
  comp <- Sys.getenv("GLIP_COMP_VENV")
  git <- tryCatch(
    system2("git", c("-C", rr_root(), "rev-parse", "--short", "HEAD"),
      stdout = TRUE, stderr = FALSE
    ),
    error = function(e) NA, warning = function(w) NA
  )
  rev <- file.path(rr_root(), "REVISION")
  .rr$versions <- list(
    R = R.version.string, platform = R.version$platform,
    glip = pv("glip"), revision = if (file.exists(rev)) readLines(rev, 1) else git[1],
    gurobi_R = pv("gurobi"),
    pcalg = pv("pcalg"), bnlearn = pv("bnlearn"), lcd = pv("lcd"),
    reticulate = pv("reticulate"),
    glipcut_python = pyv(Sys.getenv("GLIP_PYTHON"), c("gurobipy", "networkx")),
    competitors_python = pyv(
      if (nzchar(comp)) file.path(comp, "bin", "python") else "",
      c("dagma", "CausalDisco", "numpy", "torch")
    ),
    time = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  c(.rr$versions, clingo = .rr$clingo)
}

### Atomic save (write to a temporary file in the same folder, then rename)
rr_save <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(file, ".tmp", Sys.getpid())
  saveRDS(x, tmp)
  file.rename(tmp, file)
}

### Run a fit unless its file exists; errors are saved, not raised
rr_cached <- function(file, expr) {
  if (file.exists(file)) {
    return(readRDS(file))
  }
  res <- tryCatch(expr, error = function(e) {
    structure(list(error = conditionMessage(e)), class = "rr_error")
  })
  rr_save(res, file)
  res
}
