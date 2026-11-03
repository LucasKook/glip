### Background knowledge (hard constraints) for all encodings

### All tests in this file need the solver
skip_if_no_gurobi()

dag_modes <- c("dag", "dg", "dag-dcon", "dag-dc", "dg-dc")
admg_modes <- c("admg", "dmg", "admg-lean", "dmg-lean", "admg-dc", "dmg-dc")
chain_modes <- c("chain", "chain-dcon")
bk_modes <- c(dag_modes, admg_modes, chain_modes)

### Tests for simulated Gaussian data from E -> X -> Y ("chain") or
### E -> X <- Y ("collider")
.bk_tests <- function(type = c("chain", "collider"), n = 500, seed = 1) {
  type <- match.arg(type)
  set.seed(seed)
  E <- rnorm(n)
  if (type == "chain") {
    X <- E + rnorm(n)
    Y <- X + rnorm(n)
  } else {
    Y <- rnorm(n)
    X <- E + Y + rnorm(n)
  }
  dd <- data.frame(E = E, X = X, Y = Y)
  learn_graph(dd, comets = FALSE, test = "gaussCItest",
    return_tests_only = TRUE)
}

.bk_fit <- function(tests, mode, ..., cache = FALSE,
                    gurobi_args = list(OutputFlag = 0)) {
  V <- c("E", "X", "Y")
  tmp <- capture.output(opt <- .get_opt(mode)(tests,
    d = 3, V = V, mode = mode,
    cache = cache, gurobi_args = gurobi_args, ...
  ))
  opt
}

### Directed part (DAG, DG, chain graph: adjacency matrix; ADMG, DMG: M1)
.bk_dir <- function(G) {
  M <- if (is.list(G)) G$M1 else G
  matrix(as.numeric(M), nrow(M), dimnames = dimnames(M))
}

### Reachability: R[i, j] = 1 iff there is a path i -> ... -> j along M
.bk_reach <- function(M) {
  R <- M
  for (k in seq_len(NROW(M))) R <- 1 * ((R + R %*% M) > 0)
  R
}

### Number of edges (bidirected and undirected edges counted once)
.bk_nedges <- function(G, mode) {
  if (is.list(G)) {
    return(sum(G$M1) + sum(G$M2) / 2)
  }
  M <- .bk_dir(G)
  if (mode %in% chain_modes) {
    return(sum(M + t(M) > 0) / 2)
  }
  sum(M)
}

.bk_zero <- function() {
  V <- c("E", "X", "Y")
  matrix(0, 3, 3, dimnames = list(V, V))
}

test_that("Forbidden and required edges are respected in DAG and ADMG encodings", {
  tests <- .bk_tests("chain")
  Fd <- .bk_zero()
  Fd["E", "X"] <- Fd["X", "E"] <- 1
  Rq <- .bk_zero()
  Rq["Y", "X"] <- 1
  for (mode in c(dag_modes, admg_modes)) {
    opt <- .bk_fit(tests, mode, forbidden = Fd, required = Rq,
      edge_penalty = 1e-3)
    M <- .bk_dir(opt$graph)
    expect_true(all(M[Fd == 1] == 0), label = mode)
    expect_equal(M["Y", "X"], 1, label = mode)
  }
  ### Bidirected edges
  Fb <- .bk_zero()
  Fb["E", "X"] <- 1
  Rb <- .bk_zero()
  Rb["Y", "X"] <- 1
  for (mode in admg_modes) {
    opt <- .bk_fit(tests, mode, forbidden = Fd, forbidden_bi = Fb,
      required_bi = Rb, edge_penalty = 1e-3)
    expect_true(all(.bk_dir(opt$graph)[Fd == 1] == 0), label = mode)
    expect_equal(opt$graph$M2["E", "X"] + opt$graph$M2["X", "E"], 0,
      label = mode)
    expect_equal(opt$graph$M2["X", "Y"] + opt$graph$M2["Y", "X"], 2,
      label = mode)
  }
})

test_that("Forbidden and required edges are respected in chain graph encodings", {
  tests <- .bk_tests("chain")
  for (mode in chain_modes) {
    ### No directed edges between E and X and no undirected edge E - X:
    ### E and X are not adjacent
    Fd <- .bk_zero()
    Fd["E", "X"] <- Fd["X", "E"] <- 1
    Fu <- .bk_zero()
    Fu["X", "E"] <- 1
    G <- .bk_dir(.bk_fit(tests, mode, forbidden = Fd,
      forbidden_undirected = Fu, edge_penalty = 1e-3)$graph)
    expect_equal(G["E", "X"] + G["X", "E"], 0, label = mode)
    ### Only forbidding directed edges leaves the undirected edge E - X
    G <- .bk_dir(.bk_fit(tests, mode, forbidden = Fd,
      edge_penalty = 1e-3)$graph)
    expect_false(xor(G["E", "X"] == 1, G["X", "E"] == 1), label = mode)
    ### Required directed and undirected edges
    Rq <- .bk_zero()
    Rq["Y", "X"] <- 1
    Ru <- .bk_zero()
    Ru["E", "X"] <- 1
    G <- .bk_dir(.bk_fit(tests, mode, required = Rq,
      required_undirected = Ru, edge_penalty = 1e-3)$graph)
    expect_equal(c(G["Y", "X"], G["X", "Y"]), c(1, 0), label = mode)
    expect_equal(c(G["E", "X"], G["X", "E"]), c(1, 1), label = mode)
  }
})

test_that("Forbidden and required ancestral relations are respected", {
  tests <- .bk_tests("chain")
  ### Neither E nor Y is an ancestor (anterior for chain graphs) of the other
  Fa <- .bk_zero()
  Fa["E", "Y"] <- Fa["Y", "E"] <- 1
  ### Y is an ancestor (anterior for chain graphs) of E
  Ra <- .bk_zero()
  Ra["Y", "E"] <- 1
  for (mode in bk_modes) {
    opt <- .bk_fit(tests, mode, forbidden_anc = Fa, edge_penalty = 1e-3)
    R <- .bk_reach(.bk_dir(opt$graph))
    expect_equal(c(R["E", "Y"], R["Y", "E"]), c(0, 0), label = mode)
    opt <- .bk_fit(tests, mode, required_anc = Ra, edge_penalty = 1e-3)
    R <- .bk_reach(.bk_dir(opt$graph))
    expect_equal(R["Y", "E"], 1, label = mode)
    ### Combined with forbidden edges: Y -> E is excluded, so the path
    ### must go through X
    Fd <- .bk_zero()
    Fd["Y", "E"] <- 1
    opt <- .bk_fit(tests, mode, required_anc = Ra, forbidden = Fd,
      edge_penalty = 1e-3)
    M <- .bk_dir(opt$graph)
    expect_equal(.bk_reach(M)["Y", "E"], 1, label = mode)
    expect_equal(M["Y", "X"] * M["X", "E"], 1, label = mode)
  }
})

test_that("Linear edge constraints are respected and can worsen the fit", {
  tests <- .bk_tests("collider", seed = 2)
  Z <- .bk_zero()
  ### No common child X of E and Y (for ADMGs, no arrowheads at X from both
  ### E and Y)
  cf <- Z
  cf["E", "X"] <- cf["Y", "X"] <- 1
  cf_bi <- Z
  cf_bi["E", "X"] <- cf_bi["X", "E"] <- cf_bi["Y", "X"] <- cf_bi["X", "Y"] <- 1
  ### E and Y adjacent by a directed edge
  ca <- Z
  ca["E", "Y"] <- ca["Y", "E"] <- 1
  for (mode in bk_modes) {
    is_admg <- mode %in% admg_modes
    opt0 <- .bk_fit(tests, mode)
    con <- list(coef = cf, sense = "<=", rhs = 1)
    if (is_admg) con$coef_bi <- cf_bi
    opt1 <- .bk_fit(tests, mode, edge_constraints = list(con))
    M <- .bk_dir(opt1$graph)
    lhs <- M["E", "X"] + M["Y", "X"]
    if (is_admg) lhs <- lhs + opt1$graph$M2["E", "X"] + opt1$graph$M2["Y", "X"]
    expect_lte(lhs, 1, label = mode)
    expect_gt(opt1$optim$objval, opt0$optim$objval + 0.5, label = mode)
    ### A single constraint can be passed without the enclosing list
    opt2 <- .bk_fit(tests, mode,
      edge_constraints = list(coef = ca, sense = ">=", rhs = 1))
    M <- .bk_dir(opt2$graph)
    expect_gte(M["E", "Y"] + M["Y", "E"], 1, label = mode)
  }
  ### Constraints on bidirected edges only
  for (mode in admg_modes) {
    cb <- Z
    cb["E", "Y"] <- cb["Y", "E"] <- 1
    opt <- .bk_fit(tests, mode,
      edge_constraints = list(coef_bi = cb, sense = "=", rhs = 1))
    expect_equal(opt$graph$M2["E", "Y"], 1, label = mode)
  }
})

test_that("Edge penalty breaks ties in favour of sparse graphs", {
  ### Only marginal tests: all pairs dependent, so every connected graph
  ### without colliders has zero mismatches; the sparsest have two edges
  tests <- .bk_tests("chain")
  tests <- tests[tests$size == 0, ]
  for (mode in bk_modes) {
    opt0 <- .bk_fit(tests, mode)
    opt1 <- .bk_fit(tests, mode, edge_penalty = 1e-3)
    n0 <- .bk_nedges(opt0$graph, mode)
    n1 <- .bk_nedges(opt1$graph, mode)
    expect_equal(floor(opt1$optim$objval), round(opt0$optim$objval),
      label = mode)
    expect_equal(n1, 2, label = mode)
    expect_lte(n1, n0, label = mode)
  }
})

test_that("Constraints do not leak into the cache", {
  tests <- .bk_tests("chain")
  Fd <- .bk_zero()
  Fd["X", "E"] <- 1
  Ra <- .bk_zero()
  Ra["E", "Y"] <- 1
  cf <- .bk_zero()
  cf["E", "X"] <- cf["X", "Y"] <- 1
  con <- list(coef = cf, sense = "<=", rhs = 1)
  for (mode in bk_modes) {
    d1 <- file.path(tempdir(), paste0("bk-cache-1-", mode))
    d2 <- file.path(tempdir(), paste0("bk-cache-2-", mode))
    unlink(c(d1, d2), recursive = TRUE)
    fit <- function(...) {
      .bk_fit(tests, mode, forbidden = Fd, required_anc = Ra,
        edge_constraints = list(con), edge_penalty = 1e-3, ...)
    }
    opt0 <- fit(cache = FALSE)
    opt1 <- fit(cache = TRUE, cache_dir = d1) # writes cache
    opt2 <- fit(cache = TRUE, cache_dir = d1) # reads cache
    tmp <- .bk_fit(tests, mode, cache = TRUE, cache_dir = d2)
    expect_equal(opt1$optim$objval, opt0$optim$objval, label = mode)
    expect_equal(opt2$optim$objval, opt0$optim$objval, label = mode)
    expect_equal(.bk_dir(opt1$graph), .bk_dir(opt0$graph), label = mode)
    expect_equal(.bk_dir(opt2$graph), .bk_dir(opt0$graph), label = mode)
    M <- .bk_dir(opt2$graph)
    expect_equal(M["X", "E"], 0, label = mode)
    expect_equal(.bk_reach(M)["E", "Y"], 1, label = mode)
    f1 <- list.files(d1)
    expect_true(length(f1) > 0, label = mode)
    expect_equal(f1, list.files(d2), label = mode)
    expect_equal(unname(tools::md5sum(file.path(d1, f1))),
      unname(tools::md5sum(file.path(d2, f1))), label = mode)
    unlink(c(d1, d2), recursive = TRUE)
  }
})

test_that("Stopping without an incumbent returns the solver output", {
  tests <- .bk_tests("chain")
  ### Infeasible: more edges than there are pairs
  infeas <- list(coef = 1 - diag(3), sense = ">=", rhs = 100)
  for (mode in bk_modes) {
    opt <- .bk_fit(tests, mode, edge_constraints = infeas)
    expect_s3_class(opt, "graphopt")
    expect_null(opt$graph, label = mode)
    expect_null(opt$tests, label = mode)
    expect_equal(opt$optim$status, "INFEASIBLE", label = mode)
    ### Stopped at the root by BestBdStop: no error, a graph is returned
    ### only if an incumbent (e.g., from the warm start) exists
    opt <- .bk_fit(tests, mode,
      gurobi_args = list(OutputFlag = 0, BestBdStop = -1))
    expect_equal(opt$optim$status, "USER_OBJ_LIMIT", label = mode)
    expect_equal(is.null(opt$graph), is.null(opt$optim$x), label = mode)
  }
  ### Cyclic required ancestry is infeasible for acyclic graphs
  Ra <- .bk_zero()
  Ra["E", "Y"] <- Ra["Y", "E"] <- 1
  for (mode in c("dag", "dag-dcon", "dag-dc", "admg", "admg-lean", "admg-dc")) {
    expect_null(.bk_fit(tests, mode, required_anc = Ra)$graph, label = mode)
  }
})

test_that("learn_graph passes background knowledge to all encodings", {
  set.seed(1)
  n <- 500
  E <- rnorm(n)
  X <- E + rnorm(n)
  Y <- X + rnorm(n)
  dd <- data.frame(E = E, X = X, Y = Y)
  Fd <- .bk_zero()
  Fd[, "E"] <- 1
  Ra <- .bk_zero()
  Ra["E", "Y"] <- 1
  for (mode in c("dag", "admg", "chain")) {
    tmp <- capture.output(res <- suppressWarnings(learn_graph(dd,
      mode = mode, comets = FALSE, test = "gaussCItest", cache = FALSE,
      forbidden = Fd, required_anc = Ra, edge_penalty = 1e-3,
      gurobi_args = list(OutputFlag = 0)
    )))
    M <- .bk_dir(res$graph$graph)
    expect_true(all(M[, "E"] == 0), label = mode)
    expect_equal(.bk_reach(M)["E", "Y"], 1, label = mode)
  }
  ### Infeasible constraints: no graph and no graphical representation
  tmp <- capture.output(res <- learn_graph(dd,
    comets = FALSE, test = "gaussCItest", cache = FALSE,
    edge_constraints = list(coef = 1 - diag(3), sense = ">=", rhs = 100),
    gurobi_args = list(OutputFlag = 0)
  ))
  expect_null(res$graph$graph)
  expect_null(res$computed)
  ### Unsupported types of edges raise an error
  expect_error(learn_graph(dd, comets = FALSE, test = "gaussCItest",
    cache = FALSE, forbidden_bi = Fd), "not supported by dcon_lean_optim")
  expect_error(learn_graph(dd, mode = "admg", comets = FALSE,
    test = "gaussCItest", cache = FALSE, forbidden_undirected = Fd),
    "not supported by admg_cubic_optim")
  expect_error(learn_graph(dd, mode = "chain", comets = FALSE,
    test = "gaussCItest", cache = FALSE, required_bi = Fd),
    "not supported by chain_lean_optim")
})

test_that("Background knowledge is validated", {
  tests <- .bk_tests("chain")
  Z <- .bk_zero()
  Fd <- Z
  Fd["E", "X"] <- 1
  expect_error(.bk_fit(tests, "dag", forbidden = matrix(0, 2, 2)),
    "must be a 3 by 3 matrix")
  expect_error(.bk_fit(tests, "dag", forbidden = 2 * Fd), "only 0/1")
  expect_error(.bk_fit(tests, "dag", forbidden = Fd, required = Fd),
    "both forbidden and required")
  expect_error(.bk_fit(tests, "admg", forbidden_bi = Fd, required_bi = t(Fd)),
    "both forbidden and required")
  expect_error(.bk_fit(tests, "chain", required = Fd, required_undirected = Fd),
    "contradictory")
  expect_error(.bk_fit(tests, "chain", required = Fd + t(Fd)),
    "contradictory")
  expect_error(.bk_fit(tests, "dag", edge_penalty = -1), "non-negative")
  bad_names <- Fd
  dimnames(bad_names) <- list(c("a", "b", "c"), c("a", "b", "c"))
  expect_error(.bk_fit(tests, "dag", forbidden = bad_names), "must match `V`")
  expect_error(.bk_fit(tests, "dag",
    edge_constraints = list(coef = Fd, sense = "<", rhs = 1)), "`sense`")
  expect_error(.bk_fit(tests, "dag",
    edge_constraints = list(coef = Fd, sense = "<=")), "`rhs`")
  expect_error(.bk_fit(tests, "dag",
    edge_constraints = list(coef_bi = Fd, sense = "<=", rhs = 1)),
  "`coef_bi` not supported")
  expect_error(.bk_fit(tests, "admg",
    edge_constraints = list(coef_bi = Fd, sense = "<=", rhs = 1)),
  "must be symmetric")
  expect_error(.bk_fit(tests, "dag-dc", forbidden_undirected = Fd),
    "not supported by dag_optim")
  ### Matrices with dimnames are matched to V
  Fa <- Z
  Fa["E", "Y"] <- Fa["Y", "E"] <- 1
  perm <- c("Y", "E", "X")
  opt1 <- .bk_fit(tests, "dag", forbidden_anc = Fa, edge_penalty = 1e-3)
  opt2 <- .bk_fit(tests, "dag", forbidden_anc = Fa[perm, perm],
    edge_penalty = 1e-3)
  expect_equal(.bk_dir(opt2$graph), .bk_dir(opt1$graph))
})

test_that("DAG and DG modes do not share cached constraint matrices", {
  tests <- .bk_tests("chain")
  for (modes in list(c("dag", "dg"), c("dag-dcon", "dg-dcon"))) {
    cdir <- file.path(tempdir(), paste0("bk-cache-mode-", modes[1]))
    unlink(cdir, recursive = TRUE)
    opt <- .get_opt(modes[1])
    for (mode in modes) {
      tmp <- capture.output(opt(tests, d = 3, V = c("E", "X", "Y"),
        mode = mode, cache = TRUE, cache_dir = cdir,
        gurobi_args = list(OutputFlag = 0)))
    }
    fls <- list.files(cdir)
    for (mode in modes) {
      expect_true(any(grepl(paste0("max1", mode, "\\."), fls)), label = mode)
    }
    unlink(cdir, recursive = TRUE)
  }
})

test_that("Background knowledge is passed to the PC/FCI warm start", {
  V <- c("E", "X", "Y")
  Z <- .bk_zero()
  ### Gaps, fixed edges and orientations
  Fd <- Z
  Fd["E", "X"] <- Fd["X", "E"] <- 1
  Rq <- Z
  Rq["Y", "X"] <- 1
  bk <- .bk_validate(V, forbidden = Fd, required = Rq)
  out <- .bk_pcalg(bk, 3)
  expect_equal(which(out$fixedGaps), which(Fd == 1))
  expect_equal(which(out$fixedEdges), which(Rq + t(Rq) == 1))
  expect_true(out$orient["Y", "X"])
  expect_false(any(out$orient & t(out$orient)))
  ### Forbidden ancestry excludes one direction; required ancestry orients
  ### adjacent pairs only in acyclic graphs
  Fa <- Z
  Fa["E", "Y"] <- 1
  Ra <- Z
  Ra["X", "Y"] <- 1
  bk <- .bk_validate(V, forbidden_anc = Fa, required_anc = Ra)
  expect_true(.bk_pcalg(bk, 3)$orient["Y", "E"])
  expect_true(.bk_pcalg(bk, 3)$orient["X", "Y"])
  expect_false(.bk_pcalg(bk, 3, acyclic = FALSE)$orient["X", "Y"])
  ### ADMGs: no gap unless the bidirected edge is forbidden as well
  bk <- .bk_validate(V, forbidden = Fd)
  expect_null(.bk_pcalg(bk, 3, bidirected = TRUE)$fixedGaps)
  bk <- .bk_validate(V, forbidden = Fd, forbidden_bi = Fd)
  expect_equal(which(.bk_pcalg(bk, 3, bidirected = TRUE)$fixedGaps),
    which(Fd == 1))
  ### PC on data from E -> X -> Y: gaps are respected, orientations are
  ### propagated, inconsistent ones are skipped
  set.seed(1)
  n <- 500
  E <- rnorm(n)
  X <- E + rnorm(n)
  Y <- X + rnorm(n)
  suff <- list(C = cor(cbind(E, X, Y)), n = n)
  pcres <- pcalg::pc(suff, pcalg::gaussCItest, labels = V, alpha = 0.05,
    fixedGaps = .bk_pcalg(.bk_validate(V, forbidden = Fd), 3)$fixedGaps)
  PC <- as(pcres@graph, "matrix")
  expect_equal(PC["E", "X"] + PC["X", "E"], 0)
  pcres <- pcalg::pc(suff, pcalg::gaussCItest, labels = V, alpha = 0.05)
  orient <- Z == 1
  orient["E", "X"] <- orient["Y", "X"] <- TRUE
  PC <- .bk_orient_pdag(pcres@graph, orient)
  expect_equal(c(PC["E", "X"], PC["X", "E"]), c(1, 0))
  ### Y -> X contradicts E -> X -> Y (implied by Meek's rule 1) and is skipped
  expect_equal(c(PC["X", "Y"], PC["Y", "X"]), c(1, 0))
  ### learn_graph() with PC (all DAG modes) and FCI warm starts
  dd <- data.frame(E = E, X = X, Y = Y)
  for (mode in c("dag", "dg", "dag-dcon", "dag-dc", "dg-dc", "admg",
    "admg-dc")) {
    args <- list(dd,
      mode = mode, comets = FALSE, test = "gaussCItest", cache = FALSE,
      forbidden = Fd, required = Rq, gurobi_args = list(OutputFlag = 0)
    )
    if (mode %in% admg_modes) args$forbidden_bi <- Fd
    tmp <- capture.output(res <- suppressWarnings(do.call(learn_graph, args)))
    G <- res$graph$graph
    M <- .bk_dir(G)
    expect_equal(M["E", "X"] + M["X", "E"], 0, label = mode)
    expect_equal(M["Y", "X"], 1, label = mode)
    if (is.list(G)) expect_equal(G$M2["E", "X"], 0, label = mode)
  }
})
