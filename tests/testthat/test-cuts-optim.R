### cuts_optim(): same optimization problem as the min-length encodings

ga <- list(Threads = 1, Seed = 1, MIPGap = 0, MIPGapAbs = 0, OutputFlag = 0)

.std_fit <- function(class, tests, d, m, V, weight, ...) {
  switch(class,
    dag = dcon_lean_optim(tests, d = d, max_size = m, V = V, weight_type = weight,
      cache = FALSE, gurobi_args = ga, mode = "dag", ...),
    admg = admg_lean_optim(tests, d = d, max_size = m, V = V, weight_type = weight,
      cache = FALSE, gurobi_args = ga, mode = "admg-lean", ...),
    chain = chain_lean_optim(tests, d = d, max_size = m, V = V, weight_type = weight,
      cache = FALSE, gurobi_args = ga, ...)
  )
}

### Loss of a graph by the R oracle (plus edge penalty)
.oracle_objective <- function(class, G, tests, weight, pen = 0) {
  mode <- c(dag = "dag", admg = "admg-dc", chain = "chain")[[class]]
  z <- vapply(seq_len(nrow(tests)), \(k) {
    C <- strsplit(as.character(tests$Z[k]), ",", fixed = TRUE)[[1]]
    .check_separation(tests$X[k], tests$Y[k], C[nzchar(C)], G, mode) == 0
  }, logical(1))
  w <- if (weight == "const") 1 else -log2(pmax(tests$p.value, 2 * .Machine$double.eps)) + .1
  X <- if (class == "admg") G$M1 else G
  sum(w * (z != (tests$p.value <= 0.05))) +
    pen * (sum(X) + if (class == "admg") sum(G$M2[upper.tri(G$M2)]) else 0)
}

test_that("cut encoding finds the optimum of the min-length encodings", {
  skip_if_no_gurobi()
  skip_if_no_glipcut()
  set.seed(1)
  for (mode in c("dag", "admg-dc", "chain")) {
    class <- .cuts_class(mode)
    for (r in 1:4) {
      d <- 4
      V <- letters[seq_len(d)]
      m <- if (r == 4) 1 else d - 2
      weight <- if (r %% 2) "const" else "log"
      G <- .generate_random_graph(d, V = V, mode = mode, prob = 0.4)
      tests <- .cuts_bank(G, V, m, mode)
      lab <- paste(class, r)
      s <- .std_fit(class, tests, d, m, V, weight)
      f <- cuts_optim(tests, d = d, max_size = m, V = V, class = class,
        weight_type = weight, gurobi_args = ga)
      expect_s3_class(f, "graphopt")
      expect_equal(f$optim$status, "OPTIMAL", label = lab)
      expect_equal(f$optim$objval, s$optim$objval, tolerance = 1e-6, label = lab)
      Gf <- if (class == "admg") list(M1 = f$graph$M1, M2 = f$graph$M2) else f$graph
      expect_equal(.oracle_objective(class, Gf, tests, weight), f$optim$objval,
        tolerance = 1e-6, label = lab)
      ### connection indicators of the returned graph
      mode_o <- if (class == "admg") "admg-dc" else mode
      sep <- .compute_oracle_tests(Gf, max_size = m, mode = mode_o)
      key <- function(x) paste(pmin(x$X, x$Y), pmax(x$X, x$Y), x$Z)
      expect_equal(f$tests$dcon, 1 - sep$p.value[match(key(f$tests), key(sep))],
        label = lab)
      expect_identical(f$optim$encoding, "cuts")
    }
  }
})

test_that("all cut options give the same optimum", {
  skip_if_no_glipcut()
  set.seed(2)
  V <- letters[1:4]
  G <- .generate_random_graph(4, V = V, mode = "admg-dc", prob = 0.4)
  tests <- .cuts_bank(G, V, 2, "admg-dc", flip = 0.3)
  obj <- vapply(list(
    list(subset = "none", joint = FALSE), list(subset = "static", joint = TRUE),
    list(subset = "lp", joint = TRUE), list(orient = "cover", prune = TRUE),
    list(nodecuts = "all"), list(nodecuts = "none")
  ), \(o) cuts_optim(tests, d = 4, V = V, class = "admg", gurobi_args = ga,
    cut_options = o)$optim$objval, numeric(1))
  expect_equal(max(obj) - min(obj), 0, tolerance = 1e-8)
})

test_that("cut encoding respects background knowledge like the min-length encodings", {
  skip_if_no_gurobi()
  skip_if_no_glipcut()
  set.seed(3)
  d <- 4
  V <- letters[seq_len(d)]
  M <- function(...) {
    B <- matrix(0, d, d, dimnames = list(V, V))
    for (p in list(...)) B[p[1], p[2]] <- 1
    B
  }
  ec <- list(coef = M(c("a", "b"), c("b", "a"), c("c", "d"), c("d", "c")), sense = ">=", rhs = 2)
  cases <- list(
    list("dag", list(forbidden = M(c("a", "b"), c("b", "c")), required = M(c("c", "d")),
      edge_penalty = 0.3)),
    list("dag", list(forbidden_anc = M(c("d", "a")), required_anc = M(c("a", "d")))),
    list("dag", list(edge_constraints = ec)),
    list("admg", list(forbidden_bi = M(c("a", "b")), required_bi = M(c("c", "d")),
      required = M(c("a", "c")), edge_penalty = 0.5)),
    list("admg", list(forbidden_anc = M(c("b", "a"), c("c", "a")),
      edge_constraints = list(coef_bi = M(c("a", "d"), c("d", "a")), sense = "=", rhs = 1))),
    list("chain", list(required_undirected = M(c("a", "b")), forbidden = M(c("c", "d")),
      edge_penalty = 0.3)),
    list("chain", list(forbidden_undirected = M(c("b", "c"), c("c", "b")),
      required = M(c("b", "c")), forbidden_anc = M(c("d", "a")))),
    list("chain", list(required_anc = M(c("a", "d")), edge_constraints = ec))
  )
  for (k in seq_along(cases)) {
    class <- cases[[k]][[1]]
    K <- cases[[k]][[2]]
    mode <- c(dag = "dag", admg = "admg-dc", chain = "chain")[[class]]
    G <- .generate_random_graph(d, V = V, mode = mode, prob = 0.4)
    tests <- .cuts_bank(G, V, d - 2, mode)
    s <- do.call(.std_fit, c(list(class, tests, d, d - 2, V, "const"), K))
    f <- do.call(cuts_optim, c(list(tests, d = d, V = V, class = class,
      gurobi_args = ga), K))
    lab <- paste("case", k, class)
    expect_equal(f$optim$status, "OPTIMAL", label = lab)
    expect_equal(f$optim$objval, s$optim$objval, tolerance = 1e-6, label = lab)
    pen <- if (is.null(K$edge_penalty)) 0 else K$edge_penalty
    Gf <- if (class == "admg") list(M1 = f$graph$M1, M2 = f$graph$M2) else f$graph
    expect_equal(.oracle_objective(class, Gf, tests, "const", pen), f$optim$objval,
      tolerance = 1e-6, label = lab)
    ### knowledge holds in the learned graph
    X <- if (class == "admg") Gf$M1 else Gf
    dir <- if (class == "chain") X * (1 - t(X)) else X
    und <- if (class == "chain") X * t(X) else 0 * X
    An <- ancestor_matrix(X)
    chk <- function(B, A, v) if (!is.null(B)) expect_true(all(A[B == 1] == v), label = lab)
    chk(K[["forbidden"]], dir, 0)
    chk(K[["required"]], dir, 1)
    chk(K[["forbidden_undirected"]], und, 0)
    chk(K[["required_undirected"]], und, 1)
    if (class == "admg") {
      chk(K[["forbidden_bi"]], Gf$M2, 0)
      chk(K[["required_bi"]], Gf$M2, 1)
    }
    chk(K[["forbidden_anc"]], An > 0, 0)
    chk(K[["required_anc"]], An > 0, 1)
  }
})

test_that("infeasible knowledge returns no graph in both encodings", {
  skip_if_no_gurobi()
  skip_if_no_glipcut()
  V <- letters[1:3]
  tests <- .cuts_bank(matrix(0, 3, 3, dimnames = list(V, V)), V, 1, "dag", flip = 0)
  R <- matrix(0, 3, 3, dimnames = list(V, V))
  R["a", "b"] <- R["b", "c"] <- 1
  Fa <- matrix(0, 3, 3, dimnames = list(V, V))
  Fa["a", "c"] <- 1
  f <- cuts_optim(tests, d = 3, V = V, class = "dag", gurobi_args = ga,
    required = R, forbidden_anc = Fa)
  s <- dcon_lean_optim(tests, d = 3, V = V, cache = FALSE, gurobi_args = ga,
    required = R, forbidden_anc = Fa)
  expect_null(f$graph)
  expect_null(s$graph)
  expect_true(f$optim$status %in% c("INFEASIBLE", "INF_OR_UNBD"))
})

test_that("cuts_optim() validates its input", {
  V <- letters[1:3]
  tests <- data.frame(X = "a", Y = "b", Z = "", p.value = 0.01, size = 0)
  B <- matrix(0, 3, 3, dimnames = list(V, V))
  B["a", "b"] <- 1
  expect_error(cuts_optim(tests, d = 3, V = V, class = "dag", forbidden_bi = B),
    "not supported")
  expect_error(cuts_optim(tests, d = 3, V = V, class = "admg", forbidden_undirected = B),
    "not supported")
  expect_error(cuts_optim(tests, d = 3, V = V, class = "chain", required_bi = B),
    "not supported")
  expect_error(cuts_optim(tests, d = 3, V = V, class = "dag",
    edge_constraints = list(coef_bi = B + t(B), sense = "<=", rhs = 0)), "not supported")
  expect_error(cuts_optim(tests, d = 3, V = V, class = "dag",
    edge_constraints = list(coef = B, sense = "<", rhs = 0)), "sense")
  expect_error(cuts_optim(tests, d = 3, V = V, class = "dag", trafo = \(p) 1 - p),
    "0 or 1")
  expect_error(cuts_optim(tests, d = 3, V = V, class = "chain", required = B,
    required_undirected = B), "contradictory")
  expect_error(cuts_optim(tests, d = 3, V = V, class = "dag",
    python = file.path(tempdir(), "no-such-python")), "Python")
  expect_error(.cuts_class("dg"), "cyclic")
})

test_that("decisions are taken per unordered pair as in the min-length encodings", {
  V <- letters[1:3]
  CC <- .cond_sets(3, 1)
  zijc <- .dcon_lean_indices_structured_r(3, CC)$zijc
  tests <- data.frame(X = c("a", "b", "a"), Y = c("b", "a", "c"),
    Z = c("", "", "b"), p.value = c(0.5, 0.01, 0.2), size = c(0, 0, 1))
  dec <- .test_decisions(tests, zijc, CC, V, nrow(zijc), "log", \(x) as.numeric(x <= 0.05))
  key <- \(i, j, C) which(zijc$i == max(i, j) & zijc$j == min(i, j) &
    vapply(CC[zijc$C], \(z) identical(z, as.integer(C)), TRUE))
  ab <- key(1, 2, integer(0))
  ac <- key(1, 3, 2)
  ### the last of duplicated tests counts
  expect_equal(dec$p[ab], 1)
  expect_equal(dec$w[ab], -log2(0.01) + 0.1)
  expect_equal(dec$p[ac], 0)
  ### untested variables have weight zero
  expect_equal(sum(dec$w > 0), 2)
  expect_equal(dec$merged$weight, dec$w[dec$merged$idx])
})

test_that("invalid warm starts are ignored with a warning", {
  skip_if_no_glipcut()
  V <- letters[1:3]
  tests <- .cuts_bank(matrix(0, 3, 3, dimnames = list(V, V)), V, 1, "dag", flip = 0)
  W <- matrix(0, 3, 3, dimnames = list(V, V))
  W["a", "b"] <- W["b", "c"] <- W["c", "a"] <- 1
  expect_warning(f <- cuts_optim(tests, d = 3, V = V, class = "dag",
    warmstart = W, gurobi_args = ga), "warm start")
  expect_equal(f$optim$objval, 0)
  ### a valid warm start is accepted silently
  expect_no_warning(cuts_optim(tests, d = 3, V = V, class = "dag",
    warmstart = W * 0, gurobi_args = ga))
})

test_that("the Python executable is taken from the argument, option or environment", {
  expect_equal(.glipcut_python("~/venv/bin/python"), path.expand("~/venv/bin/python"))
  old <- options(glip.python = "/opt/py/bin/python")
  on.exit(options(old), add = TRUE)
  expect_equal(.glipcut_python(), "/opt/py/bin/python")
  options(glip.python = NULL)
  env <- Sys.getenv("GLIP_PYTHON", NA)
  on.exit(if (is.na(env)) Sys.unsetenv("GLIP_PYTHON") else Sys.setenv(GLIP_PYTHON = env),
    add = TRUE)
  Sys.setenv(GLIP_PYTHON = "/env/py")
  expect_equal(.glipcut_python(), "/env/py")
  Sys.unsetenv("GLIP_PYTHON")
  expect_equal(.glipcut_python(), "python3")
})
