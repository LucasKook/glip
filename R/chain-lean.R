#' Optimal Learning of LWF-Chain Graphs via Mixed Integer Programming (Lean Encoding)
#'
#' Solves a global mixed integer optimization problem to learn an LWF-Chain
#' graph from supplied conditional independence test results, using the
#' `gurobi` solver. Supports warm-start, different weight types, and optional
#' caching for efficiency.
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
chain_lean_optim <- function(
    tests, d = 3, max_size = d - 2, V = letters[1:d],
    trafo = \(x) as.numeric(x <= 0.05),
    weight_type = c("const", "inv", "log", "size"),
    warmstart = NULL,
    edgehints = NULL,
    gurobi_args = list(),
    verbose = FALSE,
    cache = TRUE,
    cache_dir = "./.cache-chain-lean-new",
    forbidden = NULL,
    required = NULL,
    forbidden_undirected = NULL,
    required_undirected = NULL,
    forbidden_anc = NULL,
    required_anc = NULL,
    edge_constraints = NULL,
    edge_penalty = 0,
    ...) {

  .require_gurobi()

  weight_type <- match.arg(weight_type)
  bk <- .bk_validate(V,
    forbidden = forbidden, required = required,
    forbidden_undirected = forbidden_undirected,
    required_undirected = required_undirected,
    forbidden_anc = forbidden_anc, required_anc = required_anc,
    edge_penalty = edge_penalty
  )
  .bk_unsupported(list(...), .bk_args_chain, "chain_lean_optim()",
    "chain graphs have no bidirected edges")
  bkc <- .bk_chain(bk$forbidden, bk$required, bk$forbidden_undirected,
    bk$required_undirected, d)
  warmstart <- .bk_warmstart(warmstart, bkc$zero, bkc$one)

  fn <- file.path(
    cache_dir,
    paste0(c("lhs", "rhs", "M1", "N1", "W1", "Z1", "R1"), "dim", d, "max",
      max_size, c(".mtx", ".rds", ".rds", ".rds", ".rds", ".rds", ".rds"))
  )

  ### COMPUTE DIMENSIONS

  # List of all conditioning sets of size at most max_size
  CC <- unlist(sapply(0:max_size, \(x)
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
  n_tilde <- ifelse(d < 4, d, 2 * d - 3) - 1
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

  # Array indices for bidirected edges
  luij <- data.frame(expand.grid(i = seq_len(d), j = seq_len(d))) |>
    dplyr::filter(i < j)
  luij$idx <- 1:NROW(luij)
  n_ul <- NROW(luij)
  tmp <- luij
  colnames(tmp) <- c("j", "i", "idx")
  luij <- rbind(luij, tmp)

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
  nW1 <- n_ul + n_ul * (d - 2)
  nZ1 <- n_d + d * (d - 1) * (d - 2)
  if (max_size > 0) {
    nr1b <- d * sum(sapply(seq_len(max_size) - 1, \(x) { choose(d, x) }))
  } else {
    nr1b <- 1
  }
  bound <- c(M1 = sum(nlc), acyc = 2 * n_d, m1min = sum(nlc), r1b = nr1b)

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
    rep(0, n_d), # *-> edges
    rep(0, n_ul), # undirected z variables
    rep(0, n_ul), # undirected l variables
    rep(0, nW1), # aux W1
    rep(0, nZ1), # aux Z1
    rep(0, d * n_C), # diC
    0, # aux const for min Z1
    0, # aux const for min R1
    0 # aux const for min W1
  )
  ### Br prio and warmstart
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
    c("C", "B", "I", "B", "I", "B", "I", "B", "I"),
    c(n_z, n_z + n_d, n_z + n_d + sum(nN1), n_d, sz[["M1"]],
      n_d + n_ul, n_ul + nW1 + nZ1, d * n_C, 3)
  )
  model$lb <- rep(
    c(0, 1, 0, 1, 0, 1, -2, 0, 0, 1, d),
    c(2 * n_z + n_d, n_z + n_d + sum(nN1), n_d, sz[["M1"]],
      n_d + n_ul, n_ul + nW1, nZ1, d * n_C, 1, 1, 1)
  )
  model$ub <- rep(
    c(Inf, 1, Inf, d, 2 * d, 1, 3 * d, 1, d, 3 * d - 1, 1, 1, 0, 1, d),
    c(n_z, n_z + n_d, n_z, n_d, sum(nN1), n_d, sz[["M1"]],
      n_d + n_ul, n_ul, nW1, nZ1, d * n_C, 1, 1, 1)
  )
  ### Background knowledge as hard bounds on the edge indicators x_ij (one if
  ### i -> j or i - j) and the path indicators d_ij^-> (= 1 iff there is no
  ### anterior path from i to j), edge penalty, and rows for linear edge
  ### constraints and forbidden directed and undirected edges
  model <- .bk_bounds(model, list(
    list(
      pairs = xij, offset = 2 * n_z, zero = bkc$zero, one = bkc$one,
      penalty = bk$edge_penalty
    ),
    list(
      pairs = xij, offset = 3 * n_z + 2 * n_d + sum(nN1),
      zero = bk$required_anc, one = bk$forbidden_anc
    )
  ))
  ec <- .edge_constraints(c(.ec_list(edge_constraints), bkc$cons),
    list(coef = list(pairs = xij, offset = 2 * n_z)), V, length(model$obj))
  model$rhs <- c(
    "labs" = -p, # linearize objective
    "labs" = p, # linearize objective
    "d1d2min" = rep(0, sum(nN1)), # D1, D2 auxiliary N1
    "acyc" = rep(0, sz[["acyc"]]), # CH1ab
    "indic" = rep(c(d - 1, -1), each = n_d), # dij->
    "m1min" = rep(0, sz[["m1min"]]), # aux M1
    "ZLcons" = rep(d, n_z), # (C4) in the writeup
    "ZLcons" = rep(- d, n_z), # (C5) in the writeup
    "w1min" = rep(0, nW1),
    "z1max" = rep(0, nZ1),
    "unZL" = rep(d, n_ul), # C12-13
    "unZL" = rep(-d, n_ul), # C12-13
    "r1b" = rep(0, sz[["r1b"]])
  )
  model$sense <- c(
    rep(">=", 2 * n_z), # linearize objective
    rep("=", sum(nN1)), # For min-constraint N1
    rep(">=", sz[["acyc"]]), # For acyclicity E1
    rep("<=", 2 * n_d), # dij->
    rep("=", sz[["m1min"]]), # aux M1
    rep("<=", 2 * n_z), # Z, L consistency
    rep("=", nW1),
    rep("=", nZ1),
    rep("<=", 2 * n_ul),
    rep("=", sz[["r1b"]])
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
    W1 <- readRDS(fn[5])
    Z1 <- readRDS(fn[6])
    R1 <- readRDS(fn[7])
  } else {
    # Exact layout from the closed form sizes; the reference builders build
    # the bounded layout (attr "bound_model") and are cut below.
    sz0 <- .lean_exact_sizes("chain", "chain", d, max_size)
    mb <- make_model(sz0)$model
    attr(mb, "sizes") <- sz0
    attr(mb, "bound_model") <- function() make_model(bound)$model
    built <- .chain_lean_build_structured_r(d, max_size, CC, n_C, n_d, n_z, n_ul, n_tilde, nlc, nN1, nW1, nZ1, nr1b, ndic, dic, zijc, xij, luij, mb, verbose)
    exact <- ncol(built$A) == length(mb$obj) && nrow(built$A) == length(mb$rhs)
    ls <- if (exact) sz0 else bound
    cmp <- .lean_compact(built$A, built$rhs, if (exact) mb$sense else make_model(bound)$model$sense,
      list(N1 = built$N1, M1 = built$M1, W1 = built$W1, Z1 = built$Z1, R1 = built$R1),
      cols = c(pre = 3 * n_z + 3 * n_d + sum(nN1), M1 = ls[["M1"]],
               post = n_d + 2 * n_ul + nW1 + nZ1 + d * n_C + 3),
      cut_cols = "M1",
      rows = c(pre = 2 * n_z + sum(nN1), acyc = ls[["acyc"]], indic = 2 * n_d, m1min = ls[["m1min"]], zl = 2 * n_z,
               w1 = nW1, z1 = nZ1, unzl = 2 * n_ul, r1b = ls[["r1b"]]),
      cut_rows = c("acyc", "m1min", "r1b"))
    if (!identical(as.numeric(cmp$sizes[names(sz0)]), as.numeric(sz0))) {
      stop("Encoding block sizes differ from .lean_exact_sizes().")
    }
    A <- cmp$A
    rhs <- cmp$rhs
    sz <- cmp$sizes
    M1 <- cmp$gens$M1
    N1 <- cmp$gens$N1
    W1 <- cmp$gens$W1
    Z1 <- cmp$gens$Z1
    R1 <- cmp$gens$R1

    if (cache) {
      if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE)
      }
      Matrix::writeMM(A, fn[1])
      saveRDS(rhs, fn[2])
      saveRDS(M1, fn[3])
      saveRDS(N1, fn[4])
      saveRDS(W1, fn[5])
      saveRDS(Z1, fn[6])
      saveRDS(R1, fn[7])
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
  model$genconmin <- c(N1, M1, W1, R1)
  model$genconmax <- Z1
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
    sedge <- x[3 * n_z + 3 * n_d + sum(nN1) + sz[["M1"]] + 1:n_d]
    uz <- x[3 * n_z + 4 * n_d + sum(nN1) + sz[["M1"]] + 1:n_ul]
    ul <- x[3 * n_z + 4 * n_d + sum(nN1) + sz[["M1"]] + n_ul + 1:n_ul]
    diC <- x[3 * n_z + 4 * n_d + sum(nN1) + sz[["M1"]] +
      nW1 + nZ1 + 2 * n_ul + 1:ndic]
    # z1l$learned <- x[3 * n_z + 4 * n_d + sum(nN1) + sz[["M1"]] + nW1 + 2 * n_ul + 1:nz1]
    # print(z1l)
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
      dag, class = class(dag), edge = edge, dcon = zijC,
      antlen = leij, minlen = lijC, pind = deij,
      sedge = sedge, zuij = uz, luij = ul, diC = diC
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

