.compute_oracle_tests <- function(
    G, max_size = NULL, mode = "dag", verbose = FALSE,
    restrict_to = NULL, parallel = FALSE, ncores = NULL) {
  V <- .get_node_set(G)
  if (!is.null(restrict_to)) {
    V <- restrict_to
  }
  d <- length(V)
  if (is.null(max_size)) {
    max_size <- d - 2
  }
  sets <- .list_tests_graph(V, max_size = max_size)$sets
  if (verbose) {
    pb <- txtProgressBar(0, length(sets), style = 3, width = 60)
  }

  if (parallel) {
    nc <- min(parallel::detectCores() - 1, 15, ncores)
    plan(multisession, workers = nc)
    my_apply <- \(...) future_lapply(..., future.seed = TRUE)
  } else {
    my_apply <- lapply
  }

  res <- my_apply(seq_along(sets), \(iter) {
    if (verbose) {
      setTxtProgressBar(pb, iter)
    }
    x <- sets[[iter]]
    data.frame(
      X = x$X,
      Y = x$Y,
      Z = paste0(x$Z, collapse = ","),
      size = length(x$Z),
      p.value = .check_separation(x$X, x$Y, x$Z, G, mode)
    )
  }) |> do.call("rbind", args = _)

  if (parallel) {
    plan(sequential)
  }

  res
}

.get_node_set <- function(G) {
  if (is.list(G)) {
    return(rownames(G$M1))
  }
  rownames(G)
}

.get_size <- function(G) {
  if (is.list(G)) {
    return(nrow(G$M1))
  }
  nrow(G)
}

### Modes (classes of graphs and encodings) supported by learn_graph()
.modes <- c(
  "dag", "dg", "dag-dcon", "dg-dcon", "dag-dc", "dg-dc", "admg", "dmg",
  "admg-lean", "dmg-lean", "admg-dc", "dmg-dc", "chain", "chain-dcon"
)

### Whether the solver package gurobi can be loaded (a function of its own so
### that tests can mock it)
.has_gurobi <- function() {
  requireNamespace("gurobi", quietly = TRUE)
}

.require_gurobi <- function() {
  if (!.has_gurobi()) {
    stop(
      "Package `gurobi` is required to solve the optimization problem. It is ",
      "not on CRAN but distributed with the Gurobi Optimizer ",
      "(https://www.gurobi.com), which also requires a license."
    )
  }
  invisible(TRUE)
}

.get_opt <- function(mode) {
  switch(mode,
    "dag-dc" = dag_optim,
    "dg-dc" = dag_optim,
    "admg" = admg_cubic_optim,
    "dmg" = admg_cubic_optim,
    "admg-dc" = admg_optim,
    "dmg-dc" = admg_optim,
    "chain" = chain_lean_optim,
    "chain-dcon" = chain_optim,
    "dag-dcon" = dcon_optim,
    "dg-dcon" = dcon_optim,
    "dag" = dcon_lean_optim,
    "dg" = dcon_lean_optim,
    "admg-lean" = admg_lean_optim,
    "dmg-lean" = admg_lean_optim
  )
}

.check_separation <- function(A, B, C, G, mode) {
  switch(mode,
    "dag" = .check_msep(A, B, C, G),
    "dg" = .check_msep(A, B, C, G),
    "admg" = .check_msep(A, B, C, G),
    "admg-lean" = .check_msep(A, B, C, G),
    "admg-dc" = .check_msep(A, B, C, G),
    "dmg-dc" = .check_msep(A, B, C, G),
    "dmg" = .check_msep(A, B, C, G),
    "dmg-lean" = .check_msep(A, B, C, G),
    "chain" = .check_csep(A, B, C, G),
    "chain-dcon" = .check_csep(A, B, C, G),
    "dag-dcon" = .check_msep(A, B, C, G),
    "dg-dcon" = .check_msep(A, B, C, G),
    "dag-dc" = .check_msep(A, B, C, G),
    "dg-dc" = .check_msep(A, B, C, G),
    "mag" = .check_dsepmag(A, B, C, G),
    "pdag" = .check_dseppdag(A, B, C, G)
  )
}

.check_equivalence <- function(G1, G2, max_size, mode = "dag") {
  V <- .get_node_set(G1)
  V2 <- .get_node_set(G2)
  stopifnot(isTRUE(all.equal(V, V2)))
  T1 <- .compute_oracle_tests(G1, max_size, mode)
  T2 <- .compute_oracle_tests(G2, max_size, mode)
  identical(
    which(T1$p.value == 1),
    which(T2$p.value == 1)
  )
}

.check_csep <- function(A, B, C, G) {
  1 * csep(G, A, B, C)
}

.check_msep <- function(A, B, C, G) {
  if (!is.list(G)) {
    tmp <- G
    tmp[] <- 0
    G <- list(M1 = G, M2 = tmp)
  }
  1 * is_m_separated(G, A, B, C)
}

.check_dsepmag <- function(A, B, C, G) {
  class(G) <- c("matrix", "array")
  V <- .get_node_set(G)
  MAG <- pcalg::pag2magAM(G, x = 1, max.chordal = nrow(G) + 1)
  1 * pcalg::dsepAM(match(A, V), match(B, V), match(C, V), MAG)
}

.check_dseppdag <- function(A, B, C, G) {
  class(G) <- c("matrix", "array")
  V <- .get_node_set(G)
  DAG <- pcalg::pdag2dag(as(G, "graphNEL"))$graph
  tmp <- as(DAG, "matrix")
  .check_msep(A, B, C, tmp)
}

.generate_random_graph <- function(
    d = 3, V = letters[1:d], mode, ...) {
  switch(mode,
    "dag" = .random_dmg(d, V, acyclic = TRUE, ...)$M1,
    "dg" = .random_dmg(d, V, acyclic = FALSE, ...)$M1,
    "admg" = .random_dmg(d, V, acyclic = TRUE, ...),
    "admg-lean" = .random_dmg(d, V, acyclic = TRUE, ...),
    "admg-dc" = .random_dmg(d, V, acyclic = TRUE, ...),
    "dmg-dc" = .random_dmg(d, V, acyclic = FALSE, ...),
    "dmg" = .random_dmg(d, V, acyclic = FALSE, ...),
    "dmg-lean" = .random_dmg(d, V, acyclic = FALSE, ...),
    "chain" = .random_cg(d, V, ...),
    "chain-dcon" = .random_cg(d, V, ...),
    "dag-dcon" = .random_dmg(d, V, acyclic = TRUE, ...)$M1,
    "dg-dcon" = .random_dmg(d, V, acyclic = FALSE, ...)$M1,
    "dag-dc" = .random_dmg(d, V, acyclic = TRUE, ...)$M1,
    "dg-dc" = .random_dmg(d, V, acyclic = FALSE, ...)$M1
  )
}

.random_cg <- function(d, V = letters[1:d], ...) {
  create_cg(d, V, ...)
}

.random_dmg <- function(d, V = letters[1:d], acyclic = FALSE, ...) {
  G <- create_dmg(d, nodeNames = V, diag = FALSE, ...)
  if (acyclic) {
    G$M1[upper.tri(G$M1)] <- 0
  }
  G
}

.compute_graphical_representation <- function(G, max_size, mode) {
  d <- .get_size(G)
  if (max_size < d - 2) {
    return(NULL)
  }
  switch(mode,
    "dag" = .dag2ess(G),
    "dg" = NULL,
    "admg" = .admg2pag(G),
    "admg-lean" = .admg2pag(G),
    "admg-dc" = .admg2pag(G),
    "dmg-dc" = NULL,
    "dmg" = NULL,
    "dmg-lean" = NULL,
    "chain" = .cg2lcg(G),
    "chain-dcon" = .cg2lcg(G),
    "dag-dcon" = .dag2ess(G),
    "dg-dcon" = NULL,
    "dag-dc" = .dag2ess(G),
    "dg-dc" = NULL
  )
}

.cg2lcg <- function(G) {
  compute_largest_cg(G)
}

.dag2ess <- function(G) {
  ret <- 1 * pcalg::dag2essgraph(G)
  class(ret) <- c("ess", class(ret))
  ret
}

.admg2pag <- function(G) {
  if (!is.list(G)) {
    G <- .to_admg(G)
  }
  ret <- compute_pag(G)
  class(ret) <- c("pag", class(ret))
  ret
}

.add <- function(x) paste0(x, collapse = "+")

.rm_int <- function(x) {
  if (all(x[, 1] == 1)) {
    return(x[, -1L, drop = FALSE])
  }
  x
}

### Helper for handling args
darg <- function(x, d) {
  if (is.na(x)) {
    return(d)
  } else {
    x
  }
}

.multigrep <- function(patterns, x, ...) {
  unlist(lapply(patterns, \(p) {
    grep(p, x, ...)
  }))
}

.to_admg <- function(G) {
  tmp <- G
  tmp[] <- 0
  list(M1 = G, M2 = tmp)
}

.pag_to_admg <- function(G) {
  if ("pag" %in% class(G)) {
    class(G) <- c("matrix", "array")
  }
  mag <- pcalg::pag2magAM(G, x = 1, max.chordal = nrow(G) + 1)
  V <- rownames(G)
  d <- length(V)
  tmp <- matrix(0, nrow = d, ncol = d, dimnames = list(V, V))
  admg <- list(M1 = tmp, M2 = tmp)
  for (i in V) {
    for (j in setdiff(V, i)) {
      admg$M1[i, j] <- 1 * (mag[i, j] == 2 & mag[j, i] == 3)
      admg$M2[i, j] <- 1 * (mag[i, j] == 2 & mag[j, i] == 2)
    }
  }
  admg
}

.ess_to_dag <- function(G) {
  if ("ess" %in% class(G) | "pag" %in% class(G)) {
    class(G) <- c("matrix", "array")
  }
  as(pcalg::pdag2dag(as(G, "graphNEL"))$graph, "matrix")
}

### Background knowledge helpers ------------------------------------------

### Check and normalize a d x d 0/1 background knowledge matrix; rows and
### columns are matched to V by name if B has dimnames
.bk_check <- function(B, V, name, symmetric = FALSE) {
  if (is.null(B)) {
    return(NULL)
  }
  B <- .bk_align(B, V, name)
  if (anyNA(B) || !all(B %in% c(0, 1))) {
    stop("`", name, "` must contain only 0/1 (or logical) entries.")
  }
  B <- 1 * B
  diag(B) <- 0
  if (symmetric) {
    B <- 1 * (B + t(B) > 0)
  }
  B
}

### Check dimensions of a d x d matrix and order rows and columns as in V
.bk_align <- function(B, V, name) {
  d <- length(V)
  B <- as.matrix(B)
  if (!identical(dim(B), c(d, d))) {
    stop("`", name, "` must be a ", d, " by ", d, " matrix.")
  }
  rn <- rownames(B)
  cn <- colnames(B)
  if (!is.null(rn) || !is.null(cn)) {
    if (!setequal(rn, V) || !setequal(cn, V)) {
      stop("Row and column names of `", name, "` must match `V`.")
    }
    B <- B[V, V, drop = FALSE]
  }
  if (!is.numeric(B) && !is.logical(B)) {
    stop("`", name, "` must be a numeric or logical matrix.")
  }
  B
}

### Error if a pair is both forbidden and required
.bk_conflict <- function(Fb, Rq, what) {
  if (!is.null(Fb) && !is.null(Rq) && any(Fb == 1 & Rq == 1)) {
    stop(what, " cannot be both forbidden and required.")
  }
  invisible(NULL)
}

### Validate and normalize all background knowledge arguments
.bk_validate <- function(V, forbidden = NULL, required = NULL,
                         forbidden_bi = NULL, required_bi = NULL,
                         forbidden_undirected = NULL,
                         required_undirected = NULL,
                         forbidden_anc = NULL, required_anc = NULL,
                         edge_penalty = 0) {
  bk <- list(
    forbidden = .bk_check(forbidden, V, "forbidden"),
    required = .bk_check(required, V, "required"),
    forbidden_bi = .bk_check(forbidden_bi, V, "forbidden_bi", symmetric = TRUE),
    required_bi = .bk_check(required_bi, V, "required_bi", symmetric = TRUE),
    forbidden_undirected = .bk_check(forbidden_undirected, V,
      "forbidden_undirected", symmetric = TRUE),
    required_undirected = .bk_check(required_undirected, V,
      "required_undirected", symmetric = TRUE),
    forbidden_anc = .bk_check(forbidden_anc, V, "forbidden_anc"),
    required_anc = .bk_check(required_anc, V, "required_anc"),
    edge_penalty = .bk_penalty(edge_penalty)
  )
  .bk_conflict(bk$forbidden, bk$required, "Edges")
  .bk_conflict(bk$forbidden_bi, bk$required_bi, "Bidirected edges")
  .bk_conflict(bk$forbidden_undirected, bk$required_undirected,
    "Undirected edges")
  .bk_conflict(bk$forbidden_anc, bk$required_anc, "Ancestral relations")
  bk
}

### Validate edge_penalty
.bk_penalty <- function(edge_penalty) {
  if (!is.numeric(edge_penalty) || length(edge_penalty) != 1 ||
    !is.finite(edge_penalty) || edge_penalty < 0) {
    stop("`edge_penalty` must be a single non-negative number.")
  }
  edge_penalty
}

### Background knowledge arguments not supported by all encodings
.bk_args_dag <- c(
  "forbidden_bi", "required_bi", "forbidden_undirected", "required_undirected"
)
.bk_args_admg <- c("forbidden_undirected", "required_undirected")
.bk_args_chain <- c("forbidden_bi", "required_bi")

### Error for background knowledge arguments (passed via ...) that the
### encoding does not support
.bk_unsupported <- function(dots, args, encoding, reason) {
  bad <- intersect(names(dots), args)
  if (length(bad)) {
    stop(
      paste0("`", bad, "`", collapse = ", "), " not supported by ",
      encoding, ": ", reason, "."
    )
  }
  invisible(NULL)
}

### Set bounds of the variables x[offset + idx] for the ordered pairs (i, j)
### with B[i, j] = 1; `idx` is the variable index of the pair in `pairs`
.bk_bound <- function(model, B, pairs, offset, which = c("ub", "lb"), value) {
  if (is.null(B)) {
    return(model)
  }
  which <- match.arg(which)
  sel <- pairs$idx[B[cbind(pairs$i, pairs$j)] == 1]
  model[[which]][offset + sel] <- value
  model
}

### Fix indicator variables by bounds and add a penalty to their objective
### coefficients; each element of `blocks` is a list with `pairs`, `offset`,
### and optionally `zero` and `one` (d x d 0/1 matrices, [i, j] = 1 fixes the
### variable of pair (i, j) to zero or one) and `penalty`
.bk_bounds <- function(model, blocks) {
  for (b in blocks) {
    model <- .bk_bound(model, b$zero, b$pairs, b$offset, "ub", 0)
    model <- .bk_bound(model, b$one, b$pairs, b$offset, "lb", 1)
    if (!is.null(b$penalty)) {
      model$obj[b$offset + unique(b$pairs$idx)] <- b$penalty
    }
  }
  model
}

### Make a warmstart adjacency matrix consistent with forbidden (zero) and
### required (one) edges
.bk_warmstart <- function(W, zero = NULL, one = NULL) {
  if (!is.matrix(W)) {
    return(W)
  }
  if (!is.null(zero)) W[zero == 1] <- 0
  if (!is.null(one)) W[one == 1] <- 1
  W
}

### Background knowledge for chain graphs, where the edge indicator x_ij is
### one if i -> j or i - j: required edges fix x_ij to one (and x_ji to zero
### for directed edges); forbidden directed (x_ij - x_ji <= 0) and undirected
### (x_ij + x_ji <= 1) edges become linear constraints
.bk_chain <- function(forbidden, required, forbidden_undirected,
                      required_undirected, d) {
  one <- zero <- matrix(0, d, d)
  if (!is.null(required)) {
    one[required == 1] <- 1
    zero[t(required) == 1] <- 1
  }
  if (!is.null(required_undirected)) {
    one[required_undirected == 1] <- 1
  }
  if (any(one == 1 & zero == 1)) {
    stop("Required edges are contradictory (e.g., `required[i, j]` and ",
      "`required[j, i]`, or `required[i, j]` and `required_undirected[i, j]`).")
  }
  .row <- function(i, j, vj) {
    cf <- matrix(0, d, d)
    cf[i, j] <- 1
    cf[j, i] <- vj
    list(coef = cf, sense = "<=", rhs = 1 * (vj == 1))
  }
  cons <- list()
  if (!is.null(forbidden)) {
    idx <- which(forbidden == 1, arr.ind = TRUE)
    cons <- c(cons, lapply(seq_len(NROW(idx)), \(k) .row(idx[k, 1], idx[k, 2], -1)))
  }
  if (!is.null(forbidden_undirected)) {
    idx <- which(forbidden_undirected == 1 & upper.tri(forbidden_undirected),
      arr.ind = TRUE)
    cons <- c(cons, lapply(seq_len(NROW(idx)), \(k) .row(idx[k, 1], idx[k, 2], 1)))
  }
  list(
    zero = if (any(zero == 1)) zero, one = if (any(one == 1)) one,
    cons = cons
  )
}

### Normalize edge_constraints to a list of constraints
.ec_list <- function(cons) {
  if (!length(cons)) {
    return(list())
  }
  if (!is.list(cons)) {
    stop("`edge_constraints` must be a list of constraints.")
  }
  if (!is.null(cons$sense)) {
    cons <- list(cons)
  }
  cons
}

### Rows for linear constraints on edge indicators; `blocks` is a named list
### whose names are the admissible coefficient fields (e.g., "coef",
### "coef_bi") and whose elements are lists with the `pairs` and column
### `offset` of the indicators (and `symmetric = TRUE` for unordered pairs);
### each constraint is a list with coefficient matrices, `sense` and `rhs`
.edge_constraints <- function(cons, blocks, V, ncol) {
  cons <- .ec_list(cons)
  if (!length(cons)) {
    return(NULL)
  }
  ii <- jj <- xx <- list()
  for (k in seq_along(cons)) {
    ck <- cons[[k]]
    if (!is.list(ck)) {
      stop("Edge constraint ", k, " must be a list.")
    }
    fields <- setdiff(names(ck), c("sense", "rhs"))
    bad <- setdiff(fields, names(blocks))
    if (length(bad)) {
      stop(
        "Edge constraint ", k, ": ", paste0("`", bad, "`", collapse = ", "),
        " not supported by this encoding (supported: ",
        paste0("`", names(blocks), "`", collapse = ", "), ")."
      )
    }
    if (!length(fields)) {
      stop("Edge constraint ", k, " has no coefficient matrix.")
    }
    if (length(ck$sense) != 1 || !ck$sense %in% c("<=", ">=", "=")) {
      stop("Edge constraint ", k, ": `sense` must be one of \"<=\", \">=\", \"=\".")
    }
    if (!is.numeric(ck$rhs) || length(ck$rhs) != 1 || !is.finite(ck$rhs)) {
      stop("Edge constraint ", k, ": `rhs` must be a single finite number.")
    }
    for (f in fields) {
      cf <- 1 * .bk_align(ck[[f]], V, f)
      if (anyNA(cf)) {
        stop("Edge constraint ", k, ": `", f, "` must not contain NAs.")
      }
      pairs <- blocks[[f]]$pairs
      if (isTRUE(blocks[[f]]$symmetric)) {
        if (!isSymmetric(unname(cf))) {
          stop("Edge constraint ", k, ": `", f, "` must be symmetric.")
        }
        pairs <- pairs[pairs$i < pairs$j, ]
      }
      v <- cf[cbind(pairs$i, pairs$j)]
      nz <- v != 0
      ii <- c(ii, list(rep(k, sum(nz))))
      jj <- c(jj, list(blocks[[f]]$offset + pairs$idx[nz]))
      xx <- c(xx, list(v[nz]))
    }
  }
  list(
    A = Matrix::sparseMatrix(
      i = unlist(ii), j = unlist(jj), x = unlist(xx),
      dims = c(length(cons), ncol)
    ),
    rhs = vapply(cons, \(x) as.numeric(x$rhs), numeric(1)),
    sense = vapply(cons, \(x) x$sense, character(1))
  )
}

### Append the rows of .edge_constraints() to the model
.bk_rows <- function(model, ec) {
  if (is.null(ec) || !length(ec$rhs)) {
    return(model)
  }
  model$A <- rbind(model$A, ec$A)
  model$rhs <- c(model$rhs, ec$rhs)
  model$sense <- c(model$sense, ec$sense)
  model
}

### Return value when the solver stops without an incumbent (e.g., the
### constraints are infeasible or the solver stopped early)
.no_incumbent <- function(sol) {
  structure(list(graph = NULL, tests = NULL, optim = sol), class = "graphopt")
}

### Background knowledge for the PC/FCI warm start in learn_graph(): pairs
### that cannot be adjacent (fixedGaps) or must be adjacent (fixedEdges), and
### orientations i -> j implied for adjacent pairs (orient[i, j] = TRUE)
.bk_pcalg <- function(bk, d, acyclic = TRUE, bidirected = FALSE) {
  nz <- function(B) if (is.null(B)) matrix(0, d, d) else B
  ### i -> j impossible if forbidden or i must not be an ancestor of j
  no_dir <- nz(bk$forbidden) == 1 | nz(bk$forbidden_anc) == 1
  adj <- nz(bk$required) == 1 | nz(bk$required_bi) == 1 |
    nz(bk$required_undirected) == 1
  adj <- adj | t(adj)
  gap <- no_dir & t(no_dir) & !adj
  if (bidirected) {
    gap <- gap & nz(bk$forbidden_bi) == 1
  }
  ### For adjacent pairs, i -> j if required or if j -> i is impossible; in
  ### acyclic graphs also if i must be an ancestor of j
  orient <- nz(bk$required) == 1 | t(no_dir)
  if (acyclic) {
    orient <- orient | nz(bk$required_anc) == 1
  }
  orient <- orient & !t(orient)
  list(
    fixedGaps = if (any(gap)) gap,
    fixedEdges = if (any(adj)) adj,
    orient = orient
  )
}

### Orient undirected edges of a PC result (graphNEL) according to `orient`,
### one at a time with Meek's rules (pcalg::addBgKnowledge()); orientations
### inconsistent with the graph are skipped. Returns an adjacency matrix with
### [i, j] = 1 for i -> j (both entries one for undirected edges)
.bk_orient_pdag <- function(g, orient) {
  M <- as(g, "matrix")
  V <- rownames(M)
  idx <- which(orient & M == 1 & t(M) == 1, arr.ind = TRUE)
  for (k in seq_len(NROW(idx))) {
    res <- tryCatch(
      pcalg::addBgKnowledge(g, x = V[idx[k, 1]], y = V[idx[k, 2]]),
      error = \(e) NULL
    )
    if (!is.null(res)) {
      g <- res
    }
  }
  as(g, "matrix")
}
