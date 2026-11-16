### learn_graph(encoding = "cuts"): interface and agreement with the
### min-length encodings on data

.lg_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  a <- rnorm(n)
  b <- a + rnorm(n)
  c <- rnorm(n)
  d <- b + c + rnorm(n)
  data.frame(a = a, b = b, c = c, d = d)
}
ga <- list(Threads = 1, Seed = 1, MIPGap = 0, MIPGapAbs = 0, OutputFlag = 0)

test_that("learn_graph() with encoding = 'cuts' matches the min-length encodings", {
  skip_if_no_gurobi()
  skip_if_no_glipcut()
  dd <- .lg_data()
  for (mode in c("dag", "dg", "admg-dc", "dmg-lean", "chain")) {
    for (max_size in c(1, 2)) {
      args <- list(dd, mode = mode, max_size = max_size, comets = FALSE,
        test = "gaussCItest", cache = FALSE, gurobi_args = ga)
      ### (the shared PC warm start may warn that its graph is not
      ### connected; warnings of the cut encoding are tested in
      ### test-cuts-optim.R)
      s <- suppressWarnings(do.call(learn_graph, args))
      f <- suppressWarnings(do.call(learn_graph, c(args, encoding = "cuts")))
      lab <- paste(mode, max_size)
      expect_s3_class(f, "learned_graph")
      expect_identical(attr(f, "encoding"), "cuts")
      expect_identical(attr(s, "encoding"), "minlength")
      expect_identical(f$tests, s$tests)
      expect_equal(f$graph$optim$status, "OPTIMAL", label = lab)
      expect_equal(f$graph$optim$objval, s$graph$optim$objval, tolerance = 1e-6,
        label = lab)
      ### the graphical representation is computed for full banks
      if (max_size == 2 && mode %in% c("dag", "chain")) {
        expect_false(is.null(f$computed), label = lab)
      }
    }
  }
})

test_that("learn_graph() passes background knowledge to the cut encoding", {
  skip_if_no_gurobi()
  skip_if_no_glipcut()
  dd <- .lg_data()
  V <- names(dd)
  Fb <- matrix(0, 4, 4, dimnames = list(V, V))
  Fb["a", "b"] <- Fb["b", "a"] <- 1
  for (mode in c("dag", "chain")) {
    args <- list(dd, mode = mode, comets = FALSE, test = "gaussCItest",
      cache = FALSE, gurobi_args = ga, forbidden = Fb, edge_penalty = 0.01)
    s <- do.call(learn_graph, args)
    f <- do.call(learn_graph, c(args, encoding = "cuts"))
    expect_equal(f$graph$optim$objval, s$graph$optim$objval, tolerance = 1e-6,
      label = mode)
    G <- f$graph$graph
    expect_equal(G["a", "b"] * (1 - G["b", "a"]) + G["b", "a"] * (1 - G["a", "b"]), 0,
      label = mode)
  }
})

test_that("learn_graph() rejects unknown encodings and missing Python for the cuts", {
  dd <- .lg_data(50)
  expect_error(learn_graph(dd, mode = "dag", encoding = "cutz", comets = FALSE,
    test = "gaussCItest"), "should be one of")
  expect_error(learn_graph(dd, mode = "dag", encoding = "cuts", comets = FALSE,
    test = "gaussCItest", python = file.path(tempdir(), "no-such-python")),
    "requires Python")
  ### tests only: no Python and no solver needed
  res <- learn_graph(dd, mode = "dag", encoding = "cuts", comets = FALSE,
    test = "gaussCItest", return_tests_only = TRUE,
    python = file.path(tempdir(), "no-such-python"))
  expect_true(is.data.frame(res))
  expect_true(all(c("X", "Y", "Z", "p.value") %in% names(res)))
})
