#' Optimal Learning of DAGs, ADMGs and Chain Graphs with the Cut Encoding
#'
#' Solves the same optimization problem as the min-length encodings
#' (\code{\link{dcon_lean_optim}}, \code{\link{admg_lean_optim}},
#' \code{\link{chain_lean_optim}}): find a graph of the given class
#' minimizing the weighted number of disagreements with the supplied
#' conditional independence decisions (plus \code{edge_penalty} times the
#' number of edge indicators), subject to the background knowledge. Only the
#' mixed integer formulation differs: separation and connection are encoded
#' by reachability in a state graph of the separation criterion, with
#' closure rows for separated decisions, cuts separated on demand for
#' dependent decisions, joint cuts linking both, and consistency cuts on all
#' 3-node subsets. The model is solved by the Python module \code{glipcut}
#' (shipped in \code{inst/python}) with \code{gurobipy}, called as a
#' subprocess.
#'
#' @inheritParams learn_graph
#' @inheritParams dcon_lean_optim
#' @param class Graph class: \code{"dag"}, \code{"admg"} (acyclic directed
#'   mixed graphs, m-separation) or \code{"chain"} (LWF chain graphs,
#'   c-separation).
#' @param cut_options Named list of solver options of the cut encoding:
#'   \code{subset} (\code{"lp"}, default, \code{"static"} or \code{"none"}),
#'   \code{joint} (default \code{TRUE}), \code{prune} (default \code{FALSE}),
#'   \code{orient} (\code{"none"}, default, or \code{"cover"}),
#'   \code{nodecuts} (\code{"root"}, default, \code{"all"} or \code{"none"}).
#' @param python Python executable with \code{gurobipy} and \code{networkx}
#'   (default: option \code{glip.python}, environment variable
#'   \code{GLIP_PYTHON}, else \code{"python3"}).
#' @param keep_files Directory in which to keep the JSON input and output of
#'   the solver (default \code{NULL}: temporary files are deleted).
#' @param ... Not used (e.g., arguments of the min-length encodings such as
#'   \code{cache} or \code{mode} are ignored).
#'
#' @return An object of class \code{graphopt} as for the min-length
#'   encodings: \code{graph} (a matrix for DAGs and chain graphs, a list with
#'   \code{M1} and \code{M2} for ADMGs; \code{NULL} without incumbent),
#'   \code{tests} (the tests with \code{weight} and the connection indicator
#'   \code{dcon} of the learned graph) and \code{optim} (status, objective,
#'   bound, gap, runtime, node count, encoding and solver information).
#'
#' @details
#' Requires a Python installation with \code{gurobipy} and \code{networkx}
#' (see \code{inst/python/requirements.txt}) and a Gurobi licence usable by
#' \code{gurobipy}. Gurobi parameters \code{TimeLimit}, \code{Threads},
#' \code{Seed}, \code{MIPGap} and \code{MIPGapAbs} in \code{gurobi_args}
#' default to the Gurobi defaults as in the min-length encodings; other
#' parameters are passed on (except \code{OutputFlag}, \code{LogToConsole}
#' and callback related ones; the cut encoding sets \code{LazyConstraints}
#' and \code{PreCrush}). \code{edgehints} are not used. Only binary
#' decisions are supported: \code{trafo} must map p-values to 0 or 1.
#' A warm start that is not a graph of the class, or violates the
#' background knowledge, is ignored with a warning.
#'
#' @export
cuts_optim <- function(
    tests, d = 3, max_size = d - 2, V = letters[1:d],
    class = c("dag", "admg", "chain"),
    trafo = \(x) as.numeric(x <= 0.05),
    weight_type = c("const", "inv", "log", "size"),
    warmstart = NULL, edgehints = NULL, gurobi_args = list(),
    verbose = FALSE, forbidden = NULL, required = NULL,
    forbidden_bi = NULL, required_bi = NULL, forbidden_undirected = NULL,
    required_undirected = NULL, forbidden_anc = NULL, required_anc = NULL,
    edge_constraints = NULL, edge_penalty = 0, cut_options = list(),
    python = NULL, keep_files = NULL, ...) {
  class <- match.arg(class)
  weight_type <- match.arg(weight_type)
  bk <- .bk_validate(V,
    forbidden = forbidden, required = required,
    forbidden_bi = forbidden_bi, required_bi = required_bi,
    forbidden_undirected = forbidden_undirected,
    required_undirected = required_undirected,
    forbidden_anc = forbidden_anc, required_anc = required_anc,
    edge_penalty = edge_penalty
  )
  unsupported <- switch(class,
    dag = .bk_args_dag, admg = .bk_args_admg, chain = .bk_args_chain
  )
  given <- unsupported[!vapply(bk[unsupported], is.null, logical(1))]
  .bk_unsupported(stats::setNames(vector("list", length(given)), given),
    unsupported, paste0("the cut encoding for class \"", class, "\""),
    switch(class,
      dag = "DAGs have only directed edges",
      admg = "ADMGs have no undirected edges",
      chain = "chain graphs have no bidirected edges"
    )
  )
  if (class == "chain") {
    ### error for contradictory required edges, as in chain_lean_optim()
    .bk_chain(bk$forbidden, bk$required, bk$forbidden_undirected,
      bk$required_undirected, d)
  }
  ec <- .cuts_edge_constraints(edge_constraints, V,
    if (class == "admg") c("coef", "coef_bi") else "coef")

  ### Decisions and weights exactly as in the min-length encodings: one
  ### decision per unordered pair and conditioning set (the last of
  ### duplicated tests counts), tests outside the bank are ignored
  q <- .cuts_queries(tests, V, max_size)
  praw <- q$p.value
  w <- switch(weight_type,
    "const" = rep(1, nrow(q)),
    "inv" = 1 / pmax(praw, 0.001),
    "log" = -log2(pmax(praw, 2 * .Machine$double.eps)) + 0.1,
    "size" = 1 / (1 + q$size)
  )
  p <- trafo(praw)
  if (length(p) != length(praw) || anyNA(p) || !all(p %in% c(0, 1))) {
    stop("The cut encoding requires `trafo` to map p-values to 0 or 1.")
  }

  ### Warm start in glip's format, made consistent with the knowledge
  start <- NULL
  if (!is.null(warmstart)) {
    if (class == "admg") {
      ws <- if (is.list(warmstart)) warmstart else .to_admg(warmstart)
      M1 <- .bk_warmstart(1 * unname(as.matrix(ws$M1[, ])), bk$forbidden, bk$required)
      M2 <- .bk_warmstart(1 * unname(as.matrix(ws$M2[, ])), bk$forbidden_bi, bk$required_bi)
      diag(M1) <- diag(M2) <- 0
      start <- list(M1 = M1, M2 = 1 * ((M2 + t(M2)) > 0))
    } else {
      G0 <- 1 * unname(as.matrix(warmstart[, ]))
      if (class == "dag") G0 <- .bk_warmstart(G0, bk$forbidden, bk$required)
      diag(G0) <- 0
      start <- if (class == "dag") list(M1 = G0) else list(G = G0)
    }
  }

  ### Gurobi parameters: Gurobi defaults unless given
  ga <- gurobi_args
  num <- function(x, default) if (is.null(x)) default else as.numeric(x)
  opts <- c(list(
    time_limit = min(num(ga[["TimeLimit"]], 1e100), 1e100),
    threads = num(ga[["Threads"]], 0), seed = num(ga[["Seed"]], 0),
    mip_gap = num(ga[["MIPGap"]], 1e-4), mip_gap_abs = num(ga[["MIPGapAbs"]], 1e-10)
  ), cut_options)
  drop <- c("TimeLimit", "Threads", "Seed", "MIPGap", "MIPGapAbs", "OutputFlag",
    "LogToConsole", "LazyConstraints", "PreCrush")
  extra <- ga[setdiff(names(ga), drop)]
  if (length(extra)) opts$gurobi_params <- extra

  mat <- function(B) if (is.null(B)) NULL else unname(B)
  K <- list(
    forbidden = mat(bk$forbidden), required = mat(bk$required),
    forbidden_bi = mat(bk$forbidden_bi), required_bi = mat(bk$required_bi),
    forbidden_undirected = mat(bk$forbidden_undirected),
    required_undirected = mat(bk$required_undirected),
    forbidden_anc = mat(bk$forbidden_anc), required_anc = mat(bk$required_anc),
    edge_constraints = ec, edge_penalty = bk$edge_penalty
  )
  K <- K[!vapply(K, is.null, logical(1))]
  case <- list(
    class = class, d = d, max_size = max_size,
    queries = lapply(seq_len(nrow(q)), \(k) list(
      i = q$i[k] - 1L, j = q$j[k] - 1L, C = as.list(q$C[[k]] - 1L),
      p = p[k], w = w[k]
    )),
    start = start, knowledge = K, options = opts
  )
  if (verbose) {
    cat("\nSolving the cut encoding with Python (glipcut)...\n")
  }
  res <- .glipcut_run(case, python = python, keep_files = keep_files)
  for (msg in unlist(res[["warnings"]])) warning(msg, call. = FALSE)

  optim <- list(
    status = res[["status"]], objval = res[["objval"]], objbound = res[["objbound"]],
    mipgap = res[["gap"]], runtime = res[["runtime"]], nodecount = res[["nodes"]],
    objval_model = res[["objval_model"]], encoding = "cuts",
    solver = paste0("gurobipy ", paste(unlist(res[["version"]]), collapse = ".")),
    build = res[["build"]], wall = res[["wall"]], cuts = res[["cuts"]],
    nvars = res[["nvars"]], nconstrs = res[["nconstrs"]]
  )
  if (is.null(res[["M1"]])) {
    return(.no_incumbent(optim))
  }
  dn <- list(V, V)
  M1 <- matrix(unlist(res[["M1"]]), d, d, byrow = TRUE, dimnames = dn)
  graph <- switch(class,
    dag = M1,
    chain = matrix(unlist(res[["G"]]), d, d, byrow = TRUE, dimnames = dn),
    admg = list(M1 = M1, M2 = matrix(unlist(res[["M2"]]), d, d, byrow = TRUE, dimnames = dn))
  )
  merged <- q[, c("X", "Y", "Z", "p.value", "size")]
  merged$weight <- w
  merged$dcon <- unlist(res[["connected"]])
  structure(list(graph = graph, tests = merged, optim = optim),
    class = "graphopt")
}

### Graph class of the cut encoding for a mode of learn_graph()
.cuts_class <- function(mode) {
  cl <- switch(mode,
    "dag" = , "dag-dcon" = , "dag-dc" = "dag",
    "admg" = , "admg-lean" = , "admg-dc" = "admg",
    "chain" = , "chain-dcon" = "chain",
    NULL
  )
  if (is.null(cl)) {
    stop("`encoding = \"cuts\"` supports DAG, ADMG and chain graph modes, ",
      "not mode \"", mode, "\" (cyclic graphs are not supported).")
  }
  cl
}

### One decision per unordered pair and conditioning set, as in the min-length
### encodings (tests with an endpoint in Z or |Z| > max_size are ignored)
.cuts_queries <- function(tests, V, max_size) {
  Z <- strsplit(as.character(tests$Z), ",", fixed = TRUE)
  Z <- lapply(Z, \(z) z[nzchar(z)])
  i <- match(tests$X, V)
  j <- match(tests$Y, V)
  ok <- !is.na(i) & !is.na(j) & i != j &
    vapply(seq_along(Z), \(k) {
      zz <- match(Z[[k]], V)
      !anyNA(zz) && length(zz) <= max_size && !any(zz %in% c(i[k], j[k])) &&
        !anyDuplicated(zz)
    }, logical(1))
  C <- lapply(Z, \(z) sort(match(z, V)))
  key <- vapply(seq_along(Z), \(k) paste(min(i[k], j[k]), max(i[k], j[k]),
    paste(C[[k]], collapse = ","), sep = "|"), "")
  keep <- which(ok)
  keep <- keep[!duplicated(key[keep], fromLast = TRUE)]
  size <- if (!is.null(tests$size)) tests$size[keep] else lengths(C[keep])
  q <- data.frame(X = tests$X[keep], Y = tests$Y[keep], Z = tests$Z[keep],
    p.value = tests$p.value[keep], size = size, i = i[keep], j = j[keep])
  q$C <- C[keep]
  q
}

### Validate edge constraints (fields, sense, rhs, dimensions) and return
### them as lists of unnamed matrices for the solver
.cuts_edge_constraints <- function(cons, V, fields) {
  cons <- .ec_list(cons)
  lapply(seq_along(cons), \(k) {
    ck <- cons[[k]]
    if (!is.list(ck)) stop("Edge constraint ", k, " must be a list.")
    f <- setdiff(names(ck), c("sense", "rhs"))
    bad <- setdiff(f, fields)
    if (length(bad)) {
      stop("Edge constraint ", k, ": ", paste0("`", bad, "`", collapse = ", "),
        " not supported by this encoding (supported: ",
        paste0("`", fields, "`", collapse = ", "), ").")
    }
    if (!length(f)) stop("Edge constraint ", k, " has no coefficient matrix.")
    if (length(ck$sense) != 1 || !ck$sense %in% c("<=", ">=", "=")) {
      stop("Edge constraint ", k, ": `sense` must be one of \"<=\", \">=\", \"=\".")
    }
    if (!is.numeric(ck$rhs) || length(ck$rhs) != 1 || !is.finite(ck$rhs)) {
      stop("Edge constraint ", k, ": `rhs` must be a single finite number.")
    }
    out <- list(sense = ck$sense, rhs = as.numeric(ck$rhs))
    for (ff in f) {
      cf <- 1 * .bk_align(ck[[ff]], V, ff)
      if (anyNA(cf)) stop("Edge constraint ", k, ": `", ff, "` must not contain NAs.")
      if (ff == "coef_bi" && !isSymmetric(unname(cf))) {
        stop("Edge constraint ", k, ": `", ff, "` must be symmetric.")
      }
      out[[ff]] <- unname(cf)
    }
    out
  })
}

### Python executable for the cut encoding
.glipcut_python <- function(python = NULL) {
  py <- if (!is.null(python)) python else
    getOption("glip.python", Sys.getenv("GLIP_PYTHON", ""))
  ### system2() quotes the command, so expand "~" here
  if (nzchar(py)) path.expand(py) else "python3"
}

### Is the cut encoding available (Python with gurobipy and networkx)?
.has_glipcut <- function(python = NULL) {
  py <- .glipcut_python(python)
  out <- suppressWarnings(tryCatch(
    system2(py, c("-c", shQuote("import gurobipy, networkx")),
      stdout = FALSE, stderr = FALSE),
    error = function(e) 1L
  ))
  identical(as.integer(out), 0L)
}

### Run python -m glipcut on a case; returns the parsed output
.glipcut_run <- function(case, python = NULL, keep_files = NULL) {
  py <- .glipcut_python(python)
  pypath <- system.file("python", package = "glip")
  if (!nzchar(pypath) || !dir.exists(file.path(pypath, "glipcut"))) {
    stop("The Python module `glipcut` was not found in the glip installation.")
  }
  dir <- if (is.null(keep_files)) tempfile("glipcut-") else keep_files
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(keep_files)) on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  inp <- file.path(dir, "input.json")
  out <- file.path(dir, "output.json")
  log <- file.path(dir, "log.txt")
  unlink(out)
  jsonlite::write_json(case, inp, auto_unbox = TRUE, digits = NA, null = "null")
  old <- Sys.getenv("PYTHONPATH", NA)
  Sys.setenv(PYTHONPATH = if (is.na(old) || !nzchar(old)) pypath else
    paste(pypath, old, sep = .Platform$path.sep))
  on.exit(if (is.na(old)) Sys.unsetenv("PYTHONPATH") else
    Sys.setenv(PYTHONPATH = old), add = TRUE)
  status <- tryCatch(
    suppressWarnings(system2(py, c("-m", "glipcut", shQuote(inp), shQuote(out)),
      stdout = log, stderr = log)),
    error = function(e) {
      stop("Could not run Python (", py, ") for the cut encoding: ",
        conditionMessage(e), call. = FALSE)
    }
  )
  res <- if (file.exists(out)) jsonlite::read_json(out) else NULL
  if (!is.null(res[["error"]])) {
    stop("The cut encoding failed: ", res[["error"]], call. = FALSE)
  }
  if (!identical(as.integer(status), 0L) || is.null(res)) {
    msg <- if (file.exists(log)) paste(utils::tail(readLines(log, warn = FALSE), 5), collapse = "\n") else ""
    stop("The cut encoding requires Python with `gurobipy` and `networkx` ",
      "(set option `glip.python` or GLIP_PYTHON). Python (", py, ") failed:\n",
      msg, call. = FALSE)
  }
  res
}
