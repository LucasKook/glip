#' Optimal Learning of ADMGs via Mixed Integer Programming with Lean Encoding
#'
#' Solves a global mixed integer optimization problem to learn an Acyclic
#' Directed Mixed Graph (ADMG) or DMG from supplied conditional independence
#' test results, using the `gurobi` solver. Supports warm-start, different
#' weight types, and optional caching for efficiency. This "lean" variant
#' focuses on a redundant but faster encoding.
#'
#' @param tests A data frame or list of conditional independence test results.
#' @param d Integer. Number of variables/nodes.
#' @param max_size Integer. Maximum size of conditioning set.
#' @param V Character vector of variable names, of length \code{d}.
#' @param trafo Transformation function applied to test p-values (default:
#'     \code{function(x) as.numeric(x <= 0.05)}).
#' @param weight_type Character. Type of weighting ("const", "inv", "log", "size").
#' @param warmstart Optional. Matrix, list, or NULL. Used to warm-start optimization.
#' @param edgehints Optional. Edge hints to speed up or constrain the solution.
#' @param gurobi_args List of arguments passed to the `gurobi` solver.
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#' @param cache Logical. If \code{TRUE}, cache constraint matrices and results.
#' @param cache_dir Character. Directory where cache files are saved.
#' @param mode Character. Learning mode: "admg" (default) or "dmg".
#' @param ... Additional arguments passed to internals.
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
admg_lean_optim <- function(
    tests, d = 3, max_size = d - 2, V = letters[1:d], trafo = \(x) as.numeric(x <= 0.05),
    weight_type = c("const", "inv", "log", "size"), warmstart = NULL, edgehints = NULL,
    gurobi_args = list(), verbose = FALSE, cache = TRUE, cache_dir = "./.cache-admg-lean-redundant",
    mode = c("admg-lean", "dmg-lean"), forbidden = NULL, required = NULL,
    forbidden_bi = NULL, required_bi = NULL, forbidden_anc = NULL,
    required_anc = NULL, edge_constraints = NULL, edge_penalty = 0,
    ...) {

  .require_gurobi()

  if (!is.null(warmstart) & !is.list(warmstart) & is.matrix(warmstart)) {
    warmstart <- .to_admg(warmstart)
  }

  mode <- match.arg(mode)
  weight_type <- match.arg(weight_type)
  bk <- .bk_validate(V,
    forbidden = forbidden, required = required,
    forbidden_bi = forbidden_bi, required_bi = required_bi,
    forbidden_anc = forbidden_anc, required_anc = required_anc,
    edge_penalty = edge_penalty
  )
  .bk_unsupported(list(...), .bk_args_admg, "admg_lean_optim()",
    "ADMGs and DMGs have no undirected edges")
  if (!is.null(warmstart)) {
    warmstart$M1 <- .bk_warmstart(warmstart$M1, bk$forbidden, bk$required)
    warmstart$M2 <- .bk_warmstart(warmstart$M2, bk$forbidden_bi, bk$required_bi)
  }

  fn <- file.path(
    cache_dir,
    paste0(c("lhs", "rhs", "P1", "N1", "O1", "R1"), "dim", d, "max",
      max_size, mode, c(".mtx", ".rds", ".rds", ".rds", ".rds", ".rds"))
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

  ### Number of bi-directed edges, excluding diagonal
  n_b <- (d * (d - 1)) / 2

  ### tilde(n)
  n_tilde <- ifelse(d < 4, d, 2 * d - 3) - 1
  nkc <- c(
    "K1" = n_d * n_C,
    "K2" = n_d * max(d - 2, 0) * n_C,
    "K3a" = n_d * max(d - 2, 0) * n_C,
    "K3b" = n_d * max(d - 3, 0) * max(d - 2, 0) * n_C,
    "K4" = n_d * max(d - 2, 0) * n_C
  )

  indices <- .dcon_lean_indices_structured_r(d, CC)
  dic <- indices$dic
  zijc <- indices$zijc
  ndic <- NROW(dic)
  lbijc <- .lean_pair_conditioning_indices(d, n_C)

  # Array indices for directed edges
  xij <- data.frame(expand.grid(i = seq_len(d), j = seq_len(d))) |>
    dplyr::filter(i != j)
  xij$idx <- 1:NROW(xij)

  # Array indices for bidirected edges
  bxij <- data.frame(expand.grid(i = seq_len(d), j = seq_len(d))) |>
    dplyr::filter(i < j)
  bxij$idx <- 1:NROW(bxij)
  tmp <- bxij
  colnames(tmp) <- c("j", "i", "idx")
  bxij <- rbind(bxij, tmp)

  ### Number of separation variables, excluding i,j in conditioning set
  n_z <- NROW(zijc)
  n_lb <- NROW(lbijc)

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

  ### SETUP GUROBI MODEL
  model <- list()
  model$obj <- c(
    w, # corresponds to t_{ij}^C
    "zijc" = rep(0, n_z), # corresponds to z_{ij}^C
    "dxij" = rep(0, n_d), # corresponds to x^{->}_{ij}
    rep(0, n_z), # corresponds to l_{ij}^C
    rep(0, n_d), # corresponds to l_{ij}^{->}
    rep(0, sum(nN1 <- c(n_d, (d - 2) * n_d))), # For min-constraint N1
    rep(0, n_d), # corresponds to d_ij^->
    rep(0, sum(nkc)), # aux P1
    "bxij" = rep(0, n_b), # bidirected edges
    "sxij" = rep(0, n_d), # *directed edges
    rep(0, n_lb), # corresponds to l_{ij}^{<->,C}
    rep(0, n_lb), # corresponds to z_{ij}^{<->,C}
    rep(0, n_lb + n_lb * (d - 2)), # corresponds to rhs of F1 and F2
    "diC" = rep(0, d * n_C), # diC
    0,
    0
  )
  model$branchpriority <- rep(0, length(model$obj))
  model$branchpriority[.multigrep(c("xij", "bxij", "sxij"), names(model$obj))] <- 1
  ### Warm start and edge hints
  guess <- ws_xij <- numeric(n_d)
  bguess <- ws_bxij <- numeric(n_b)
  sguess <- ws_sxij <- numeric(n_d)
  for (i in seq_len(d)) {
    for (j in setdiff(seq_len(d), i)) {
      bidx <- bxij$idx[bxij$i == i & bxij$j == j]
      xidx <- xij$idx[xij$i == i & xij$j == j]
      zidx <- c(zijc$idx[zijc$i == i & zijc$j == j], zijc$idx[zijc$i == j & zijc$j == i])
      if (is.null(edgehints)) {
        guess[xidx] <- as.integer(max(p[zidx]))
        bguess[bidx] <- as.integer(max(p[zidx]))
        sguess[xidx] <- as.integer(max(p[zidx]))
      } else {
        guess[xidx] <- edgehints[i, j]
        sguess[xidx] <- edgehints[i, j]
        bguess[bidx] <- 1 * (edgehints[i, j] | edgehints[j, i])
      }
      if (is.null(warmstart)) {
        ws_xij[xidx] <- 0
        ws_bxij[bidx] <- 0
        ws_sxij[xidx] <- 0
      } else {
        ws_xij[xidx] <- warmstart$M1[i, j]
        ws_bxij[bidx] <- warmstart$M2[i, j]
        ws_sxij[xidx] <- 1 * (warmstart$M1[i, j] | warmstart$M2[i, j])
      }
    }
  }
  model$varhintval <- rep(NA, length(model$obj))
  model$varhintval[grep("zijc", names(model$obj))] <- p
  model$varhintval[grep("dxij", names(model$obj))] <- guess
  model$varhintval[grep("bxij", names(model$obj))] <- bguess
  model$varhintval[grep("sxij", names(model$obj))] <- sguess
  model$varhintpri <- rep(0, length(model$obj))
  model$varhintpri[grep("zijc", names(model$obj))] <- 2
  model$varhintpri[grep("dxij", names(model$obj))] <- 2 - guess
  model$varhintpri[grep("bxij", names(model$obj))] <- 2 - bguess
  model$varhintpri[grep("sxij", names(model$obj))] <- 2 - sguess
  model$varhintpri <- as.integer(model$varhintpri)
  model$start <- rep(NA, length(model$obj))
  model$start[.multigrep(c("dxij", "bxij", "sxij"), names(model$obj))] <- c(ws_xij, ws_bxij, ws_sxij)
  model$modelsense <- "min"
  model$vtype <- rep(c("C", "B", "I", "B", "I", "B", "I", "B", "I", "B", "I", "B"),
    c(n_z, n_z + n_d, n_z + n_d + sum(nN1), n_d, sum(nkc), n_b + n_d,
      n_lb, n_lb, n_lb + n_lb * (d - 2), ndic, 1, 1))
  model$lb <- rep(c(0, 1, 0, 1, 0, 1, 0, 1, 0, d, 1),
    c(2 * n_z + n_d, n_z + n_d + sum(nN1), n_d, sum(nkc), n_b + n_d,
      n_lb, n_lb, n_lb + n_lb * (d - 2), ndic, 1, 1))
  model$ub <- rep(c(Inf, 1, Inf, d, 2 * d, 1, 4 * d, 1, d,
    1, 2 * (d - 1) + 1, 1, d, 1),
    c(n_z, n_z + n_d, n_z, n_d, sum(nN1), n_d, sum(nkc), n_b + n_d,
      n_lb, n_lb, n_lb + n_lb * (d - 2), ndic, 1, 1))
  ### Background knowledge as hard bounds on the directed and bidirected edge
  ### indicators and the path indicators d_ij^-> (= 1 iff there is no
  ### directed path from i to j), edge penalty, and rows for linear edge
  ### constraints
  bx_off <- 3 * n_z + 3 * n_d + sum(nN1) + sum(nkc)
  model <- .bk_bounds(model, list(
    list(
      pairs = xij, offset = 2 * n_z, zero = bk$forbidden, one = bk$required,
      penalty = bk$edge_penalty
    ),
    list(
      pairs = bxij, offset = bx_off, zero = bk$forbidden_bi,
      one = bk$required_bi, penalty = bk$edge_penalty
    ),
    list(
      pairs = xij, offset = 3 * n_z + 2 * n_d + sum(nN1),
      zero = bk$required_anc, one = bk$forbidden_anc
    )
  ))
  ec <- .edge_constraints(edge_constraints, list(
    coef = list(pairs = xij, offset = 2 * n_z),
    coef_bi = list(pairs = bxij, offset = bx_off, symmetric = TRUE)
  ), V, length(model$obj))
  if (max_size > 0) {
    nr1b <- d * sum(sapply(seq_len(max_size) - 1, \(x) {
      choose(d, x)
    }))
  } else {
    nr1b <- 1
  }
  model$rhs <- c(
    "labs" = -p, # linearize objective
    "labs" = p, # linearize objective
    "d1d2min" = rep(0, sum(nN1)), # D1, D2 auxiliary N1
    "acyc" = rep(0, n_d), # E1
    "indic" = rep(c(d - 1, -1), each = n_d), # dij->
    "p1min" = rep(0, sum(nkc)), # aux P1
    "ZLcons" = rep(d, n_z), # (C4) in the writeup
    "ZLcons" = rep(- d, n_z), # (C5) in the writeup
    "C6to9" = rep(0, 3 * d * (d - 1)), # C6to9
    "C10to11" = rep(c(d, -d), each = n_lb), # C10 to 11
    "f1f2min" = rep(0, n_lb + n_lb * (d - 2)), # Min O1 aux
    "r1b" = rep(0, nr1b),
    "redundant" = rep(0, d * (d - 1) / 2)
  )
  model$sense <- c(
    rep(">=", 2 * n_z), # linearize objective
    rep("=", sum(nN1)), # For min-constraint N1
    rep(">=", n_d), # For acyclicity E1
    rep("<=", 2 * n_d), # dij->
    rep("=", sum(nkc)), # aux P1
    rep("<=", 2 * n_z), # Z, L consistency
    rep("<=", 3 * d * (d - 1)), # C6to9
    rep("<=", 2 * n_lb), # C10to11
    rep("=", n_lb + n_lb * (d - 2)), # Min O1 aux
    rep("=", nr1b),
    rep("<=", d * (d - 1) / 2)
  )

  if (cache && all(file.exists(fn))) {
    cat("\nUsing cached constraint matrix and RHS...")
    A <- Matrix::readMM(fn[1])
    rhs <- readRDS(fn[2])
    rhs[1:(2*n_z)] <- c(-p, p)
    P1 <- readRDS(fn[3])
    N1 <- readRDS(fn[4])
    O1 <- readRDS(fn[5])
    R1 <- readRDS(fn[6])
    model$rhs <- rhs
  } else {
    built <- .admg_lean_build_structured_r(d, max_size, mode, CC, n_C, n_d, n_z, n_b, n_lb, n_tilde, nkc, nN1, nr1b, ndic, dic, zijc, xij, bxij, lbijc, model, verbose)
    A <- built$A
    model$rhs <- built$rhs
    P1 <- built$P1
    N1 <- built$N1
    O1 <- built$O1
    R1 <- built$R1

    if (cache) {
      if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE)
      }
      Matrix::writeMM(A, fn[1])
      saveRDS(model$rhs, fn[2])
      saveRDS(P1, fn[3])
      saveRDS(N1, fn[4])
      saveRDS(O1, fn[5])
      saveRDS(R1, fn[6])
    }
  }

  ### INCLUDE CONSTRAINTS
  model$A <- A
  model$genconmin <- c(N1, P1, O1, R1)
  ### Linear edge constraints, appended after the (cached) constraint matrix
  model <- .bk_rows(model, ec)

  if (verbose) {
    cat("\nModel setup done. Solving now...\n")
  }

  ### SOLVE
  if (requireNamespace("gurobi")) {
    sol <- gurobi::gurobi(model, gurobi_args)
  }

  ### No incumbent (infeasible constraints or solver stopped before finding a
  ### feasible solution): return the solver output only
  if (is.null(sol$x)) {
    return(.no_incumbent(sol))
  }

  ### Convert solution to DMG
  .to_dmg <- function(x) {
    edge <- x[2 * n_z + 1:n_d]
    zijC <- x[n_z + 1:n_z]
    lijC <- x[2 * n_z + n_d + 1:n_z]
    leij <- x[3 * n_z + n_d + 1:n_d]
    nijk <- x[3 * n_z + 2 * n_d + 1:sum(nN1)]
    deij <- x[3 * n_z + 2 * n_d + sum(nN1) + 1:n_d]
    bedge <- x[3 * n_z + 3 * n_d + sum(nN1) + sum(nkc) + 1:n_b]
    sedge <- x[3 * n_z + 3 * n_d + sum(nN1) + sum(nkc) + n_b + 1:n_d]
    bminlen <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + 1:n_lb]
    bminleni <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + n_lb + 1:n_lb]
    dic$dic <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + 3 * n_lb + n_lb * (d - 2) + 1:nrow(dic)]
    M1 <- M2 <- matrix(0, nrow = d, ncol = d)
    dimnames(M1) <- dimnames(M2) <- list(V, V)
    for (i in seq_len(d)) {
      for (j in seq_len(d)) {
        if (i != j) {
          M1[i, j] <- edge[xij$idx[xij$i == i & xij$j == j]]
          M2[i, j] <- bedge[bxij$idx[bxij$i == i & bxij$j == j]]
        }
      }
    }
    structure(
      list(M1 = M1, M2 = M2), edge = edge, bedge = bedge,
      dcon = zijC, antlen = leij, minlen = lijC,
      pind = deij, sedge = sedge, bminlen = bminlen,
      bminleni = bminleni, dic = dic
    )
  }

  merged$dcon <- sol$x[n_z + 1:n_z][merged$idx]
  merged$minlen <- sol$x[2 * n_z + n_d + 1:n_z][merged$idx]

  ### RETURN
  structure(
    list(graph = .to_dmg(sol$x), tests = merged, optim = sol),
    class = "graphopt"
  )
}

