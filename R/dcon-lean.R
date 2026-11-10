#' Optimal Learning of DAGs via Mixed Integer Programming (Lean Encoding)
#'
#' Solves a global mixed integer optimization problem to learn an Directed
#' Acyclic Graph (DAG) or DG from supplied conditional independence
#' test results, using the `gurobi` solver. Supports warm-start, different
#' weight types, and optional caching for efficiency.
#'
#' @inheritParams admg_lean_optim
#' @inheritParams learn_graph
#'
#' @return An object of class \code{graphopt} with components:
#'   \describe{
#'     \item{graph}{Adjacency or edge matrices of the learned graph.}
#'     \item{tests}{Test results with optimization variables.}
#'     \item{optim}{Solver output.}
#'   }
#'
#' @details
#' Requires the `gurobi` package. Will warn if unavailable.
#'
#' @export
dcon_lean_optim <- function(
  tests, d = 3, max_size = d - 2, V = letters[1:d],
  trafo = \(x) as.numeric(x <= 0.05),
  weight_type = c("const", "inv", "log", "size"),
  warmstart = NULL,
  edgehints = NULL,
  gurobi_args = list(),
  verbose = FALSE,
  cache = TRUE,
  cache_dir = "./.cache-dcon-lean-new",
  mode = c("dag", "dg"),
  forbidden = NULL,
  required = NULL,
  forbidden_anc = NULL,
  required_anc = NULL,
  edge_constraints = NULL,
  edge_penalty = 0,
  ...
) {
  mode <- match.arg(mode)
  bk <- .bk_validate(V,
    forbidden = forbidden, required = required,
    forbidden_anc = forbidden_anc, required_anc = required_anc,
    edge_penalty = edge_penalty
  )
  .bk_unsupported(list(...), .bk_args_dag, "dcon_lean_optim()",
    "DAGs and DGs have only directed edges")
  warmstart <- .bk_warmstart(warmstart, bk$forbidden, bk$required)
  .require_gurobi()

  weight_type <- match.arg(weight_type)

  fn <- file.path(
    cache_dir,
    paste0(
      c("lhs", "rhs", "M1", "N1", "R1"), "dim", d, "max",
      max_size, mode, c(".mtx", ".rds", ".rds", ".rds", ".rds")
    )
  )

  ### COMPUTE DIMENSIONS

  # List of all conditioning sets of size at most max_size
  CC <- unlist(
    sapply(0:max_size, \(x)
    utils::combn(d, x, simplify = FALSE)),
    recursive = FALSE
  )
  if (max_size == 0) {
    CC <- list(integer(0))
  }
  n_C <- length(CC)

  ### Number of directed edges, excluding diagonal
  n_d <- d * (d - 1)

  ### tilde(n)
  nlc <- c(
    "L1" = d * (d - 1) * n_C,
    "L2" = d * (d - 1) * max(d - 2, 0) * n_C,
    "L3" = d * (d - 1) * max(d - 2, 0) * n_C,
    "L4" = d * (d - 1) * max(d - 2, 0) * n_C
  )

  indices <- .dcon_lean_indices_structured_r(d, CC)
  dic <- indices$dic
  zijc <- indices$zijc
  ndic <- NROW(dic)

  xij <- data.frame(expand.grid(i = seq_len(d), j = seq_len(d))) |>
    dplyr::filter(i != j)
  xij$idx <- 1:NROW(xij)

  ### Number of separation variables, excluding i,j in conditioning set
  n_z <- NROW(zijc)

  # Lookup table to place p-values in the correct positions
  # b/c self-edges and x -> y s.t. y in C are not tested
  l1 <- zijc
  l1$i <- V[zijc$i]
  l1$j <- V[zijc$j]
  l1$C <- sapply(CC[zijc$C], \(x) paste0(V[x], collapse = ","))
  l2 <- l1
  colnames(l1) <- c("X", "Y", "Z", "idx")
  colnames(l2) <- c("Y", "X", "Z", "idx")
  lookup <- rbind(l1, l2)

  # Merge p-values with lookup table for optimization
  # Positions that are not tested contain an edge by default
  # i -> j | C s.t. j not in C ignored => Set w to zero for those
  # i -> j | C s.t. i in C ignored => Set w to zero for those
  # i -> i | C ignored => Set w to zero for those
  w <- rep(0, n_z)
  s <- rep(0, n_z)
  praw <- rep(0, n_z)
  merged <- dplyr::left_join(tests, lookup, by = c("X", "Y", "Z"))
  praw[merged$idx] <- merged$p.value
  w[merged$idx] <- 1
  s[merged$idx] <- merged$size

  ### Weights
  w[w == 1] <- switch(weight_type,
    "const" = 1,
    "inv" = 1 / pmax(praw[w == 1], 0.001),
    "log" = -log2(pmax(praw[w == 1], 2 * .Machine$double.eps)) + 0.1,
    "size" = 1 / (1 + s[w == 1])
  )

  merged$weight <- w[merged$idx]
  p <- trafo(praw)

  # Linearize absolute value:
  # min w |z - b|
  # t >= z - b
  # t >= -(z - b)
  # min wt

  nN1 <- c(n_d, (d - 2) * n_d)
  if (max_size > 0) {
    nr1b <- d * sum(sapply(seq_len(max_size) - 1, \(x) {
      choose(d, x)
    }))
  } else {
    nr1b <- 1
  }
  bound <- c(M1 = sum(nlc), acyc = 2 * n_d, m1min = sum(nlc), r1b = nr1b, redundant = d * (d - 1) / 2)

  ### SETUP GUROBI MODEL
  # Model for given block sizes sz (bounds while building, exact sizes from
  # the compact cache otherwise; see .lean_compact()).
  make_model <- function(sz) {
  model <- list()
  model$obj <- c(
    w, # corresponds to t_{ij}^C
    "zijc" = rep(0, n_z), # corresponds to z_{ij}^C
    "xij" = rep(0, n_d), # corresponds to x^{->}_{ij}
    rep(0, n_z), # corresponds to l_{ij}^C
    rep(0, n_d), # corresponds to l_{ij}^{->}
    rep(0, sum(nN1)), # For min-constraint N1
    rep(0, n_d), # corresponds to d_ij^->
    rep(0, sz[["M1"]]), # aux M1
    rep(0, d * n_C), # diC
    0 # aux const for min R1
  )
  model$branchpriority <- rep(0, length(model$obj))
  model$branchpriority[grep("xij", names(model$obj))] <- 1
  guess <- rep(0, n_d)
  ws <- rep(0, n_d)
  for (i in seq_len(d)) {
    for (j in setdiff(seq_len(d), i)) {
      xidx <- xij$idx[xij$i == i & xij$j == j]
      zidx <- c(zijc$idx[zijc$i == i & zijc$j == j], zijc$idx[zijc$i == j & zijc$j == i])
      if (is.null(edgehints)) {
        guess[xidx] <- max(p[zidx])
      } else {
        guess[xidx] <- edgehints[i, j]
      }
      if (is.null(warmstart)) {
        ws[xidx] <- 0
      } else {
        ws[xidx] <- warmstart[i, j]
      }
    }
  }
  model$varhintval <- rep(NA, length(model$obj))
  model$varhintval[grep("zijc", names(model$obj))] <- p
  model$varhintval[grep("xij", names(model$obj))] <- guess
  model$varhintpri <- rep(0, length(model$obj))
  model$varhintpri[grep("zijc", names(model$obj))] <- 2
  model$varhintpri[grep("xij", names(model$obj))] <- 2 - guess
  model$varhintpri <- as.integer(model$varhintpri)
  model$start <- rep(NA, length(model$obj))
  model$start[grep("xij", names(model$obj))] <- ws
  model$modelsense <- "min"
  model$vtype <- rep(
    c("C", "B", "I", "B", "I", "B", "B"),
    c(n_z, n_z + n_d, n_z + n_d + sum(nN1), n_d, sz[["M1"]], d * n_C, 1)
  )
  model$lb <- rep(
    c(0, 1, 0, 1, 0, 1),
    c(2 * n_z + n_d, n_z + n_d + sum(nN1), n_d, sz[["M1"]], d * n_C, 1)
  )
  model$ub <- rep(
    c(Inf, 1, Inf, d, 2 * d, 1, 3 * d, 1, 1),
    c(n_z, n_z + n_d, n_z, n_d, sum(nN1), n_d, sz[["M1"]], d * n_C, 1)
  )
  ### Background knowledge as hard bounds on the edge indicators x_ij and
  ### the path indicators d_ij^-> (= 1 iff there is no directed path from i
  ### to j), edge penalty, and rows for linear edge constraints
  model <- .bk_bounds(model, list(
    list(
      pairs = xij, offset = 2 * n_z, zero = bk$forbidden, one = bk$required,
      penalty = bk$edge_penalty
    ),
    list(
      pairs = xij, offset = 3 * n_z + 2 * n_d + sum(nN1),
      zero = bk$required_anc, one = bk$forbidden_anc
    )
  ))
  ec <- .edge_constraints(edge_constraints,
    list(coef = list(pairs = xij, offset = 2 * n_z)), V, length(model$obj))
  model$rhs <- c(
    "labs" = -p, # linearize objective
    "labs" = p, # linearize objective
    "d1d2min" = rep(0, sum(nN1)), # D1, D2 auxiliary N1
    "acyc" = rep(0, sz[["acyc"]]), # CH1ab
    "indic" = rep(c(d - 1, -1), each = n_d), # dij->
    "m1min" = rep(0, sz[["m1min"]]), # aux M1
    "ZLcons" = rep(d, n_z), # (C4*) in the writeup
    "ZLcons" = rep(-d, n_z), # (C5*) in the writeup
    "r1b" = rep(0, sz[["r1b"]]),
    "redundant" = rep(0, sz[["redundant"]])
  )
  model$sense <- c(
    rep(">=", 2 * n_z), # linearize objective
    rep("=", sum(nN1)), # For min-constraint N1
    rep(">=", sz[["acyc"]]), # For acyclicity E1
    rep("<=", 2 * n_d), # dij->
    rep("=", sz[["m1min"]]), # aux M1
    rep("<=", 2 * n_z), # Z, L consistency
    rep("=", sz[["r1b"]]),
    rep("<=", sz[["redundant"]])
  )
  list(model = model, ec = ec)
  }

  if (.lean_cache_ok(cache, fn)) {
    cat("\nUsing cached constraint matrix and RHS...")
    A <- Matrix::readMM(fn[1])
    rhs <- readRDS(fn[2])
    sz <- .lean_cached_sizes(rhs, fn[2])
    M1 <- readRDS(fn[3])
    N1 <- readRDS(fn[4])
    R1 <- readRDS(fn[5])
  } else {
    # Exact layout from the closed form sizes; the reference builders build
    # the bounded layout (attr "bound_model") and are cut below.
    sz0 <- .lean_exact_sizes("dcon", mode, d, max_size)
    mb <- make_model(sz0)$model
    attr(mb, "sizes") <- sz0
    attr(mb, "bound_model") <- function() make_model(bound)$model
    built <- .dcon_lean_build_structured_r(d, max_size, mode, CC, dic, zijc, xij,
                                n_z, n_d, nN1, nlc, nr1b, mb, verbose)
    exact <- ncol(built$A) == length(mb$obj) && nrow(built$A) == length(mb$rhs)
    ls <- if (exact) sz0 else bound
    cmp <- .lean_compact(built$A, built$rhs, if (exact) mb$sense else make_model(bound)$model$sense, list(N1 = built$N1, M1 = built$M1, R1 = built$R1),
      cols = c(pre = 3 * n_z + 3 * n_d + sum(nN1), M1 = ls[["M1"]], post = d * n_C + 1), cut_cols = "M1",
      rows = c(pre = 2 * n_z + sum(nN1), acyc = ls[["acyc"]], indic = 2 * n_d, m1min = ls[["m1min"]], zl = 2 * n_z,
               r1b = ls[["r1b"]], redundant = ls[["redundant"]]),
      cut_rows = c("acyc", "m1min", "r1b", "redundant"))
    if (!identical(as.numeric(cmp$sizes[names(sz0)]), as.numeric(sz0))) {
      stop("Encoding block sizes differ from .lean_exact_sizes().")
    }
    A <- cmp$A
    rhs <- cmp$rhs
    sz <- cmp$sizes
    M1 <- cmp$gens$M1
    N1 <- cmp$gens$N1
    R1 <- cmp$gens$R1

    if (cache) {
      if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE)
      }
      Matrix::writeMM(A, fn[1])
      saveRDS(rhs, fn[2])
      saveRDS(M1, fn[3])
      saveRDS(N1, fn[4])
      saveRDS(R1, fn[5])
    }
  }

  ### Model with the exact block sizes
  mm <- make_model(sz)
  model <- mm$model
  ec <- mm$ec
  if (length(rhs) != length(model$sense) || ncol(A) != length(model$obj)) {
    stop("Cached model does not match the encoding; delete ", cache_dir, ".")
  }
  rhs[1:(2 * n_z)] <- c(-p, p)
  attr(rhs, "sizes") <- NULL
  model$rhs <- rhs

  ### INCLUDE CONSTRAINTS
  model$A <- A
  model$genconmin <- c(N1, M1, R1)
  ### Linear edge constraints, appended after the (cached) constraint matrix
  model <- .bk_rows(model, ec)

  if (verbose) {
    cat("\nModel setup done. Solving now...\n")
  }

  ### SOLVE
  if (requireNamespace("gurobi")) {
    # return(gurobi::gurobi_iis(model))
    sol <- gurobi::gurobi(model, gurobi_args)
  }

  ### No incumbent (infeasible constraints or solver stopped before finding a
  ### feasible solution): return the solver output only
  if (is.null(sol$x)) {
    return(.no_incumbent(sol))
  }

  ### Convert solution to DMG
  .to_dag <- function(x) {
    edge <- x[2 * n_z + 1:n_d]
    zijC <- x[n_z + 1:n_z]
    lijC <- x[2 * n_z + n_d + 1:n_z]
    leij <- x[3 * n_z + n_d + 1:n_d]
    nijk <- x[3 * n_z + 2 * n_d + 1:sum(nN1)]
    deij <- x[3 * n_z + 2 * n_d + sum(nN1) + 1:n_d]
    diC <- x[3 * n_z + 3 * n_d + sum(nN1) + sz[["M1"]] + 1:ndic]
    dag <- matrix(0, nrow = d, ncol = d)
    dimnames(dag) <- list(V, V)
    for (i in seq_len(d)) {
      for (j in seq_len(d)) {
        if (i != j) {
          dag[i, j] <- edge[xij$idx[xij$i == i & xij$j == j]]
        }
      }
    }
    structure(
      dag,
      class = class(dag), edge = edge, dcon = zijC,
      antlen = leij, minlen = lijC, pind = deij, diC = diC
    )
  }

  merged$dcon <- sol$x[n_z + 1:n_z][merged$idx]
  merged$minlen <- sol$x[2 * n_z + n_d + 1:n_z][merged$idx]

  ### RETURN
  structure(
    list(graph = .to_dag(sol$x), tests = merged, optim = sol),
    class = "graphopt"
  )
}
