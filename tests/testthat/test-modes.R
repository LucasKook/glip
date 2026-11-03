### Modes, conditioning set sizes, and solver availability

test_that("All encodings recover the separations of oracle tests for every max_size", {
  skip_if_no_gurobi()
  d <- 4
  V <- letters[1:d]
  for (mode in .modes) {
    for (max_size in 0:(d - 2)) {
      set.seed(100 * max_size + match(mode, .modes))
      G <- .generate_random_graph(d, V = V, mode = mode, prob = 0.4)
      tests <- .compute_oracle_tests(G, max_size = max_size, mode = mode)
      lab <- paste0(mode, ", max_size = ", max_size)
      tmp <- capture.output(opt <- .get_opt(mode)(tests,
        d = d, max_size = max_size, V = V, mode = mode, cache = FALSE,
        gurobi_args = list(OutputFlag = 0)
      ))
      expect_equal(opt$optim$objval, 0, tolerance = 1e-6, label = lab)
      Gl <- opt$graph
      Gl <- if (is.list(Gl)) list(M1 = Gl$M1[, ], M2 = Gl$M2[, ]) else Gl[, ]
      expect_identical(
        .compute_oracle_tests(Gl, max_size = max_size, mode = mode)$p.value,
        tests$p.value,
        label = lab
      )
    }
  }
})

test_that("max_size = 0 works with the cache and through learn_graph()", {
  skip_if_no_gurobi()
  d <- 4
  V <- letters[1:d]
  set.seed(3)
  dd <- as.data.frame(matrix(rnorm(400 * d), ncol = d, dimnames = list(NULL, V)))
  dd$b <- dd$a + rnorm(400)
  for (mode in .modes) {
    cdir <- file.path(tempdir(), paste0("max0-", mode))
    unlink(cdir, recursive = TRUE)
    G <- .generate_random_graph(d, V = V, mode = mode, prob = 0.4)
    tests <- .compute_oracle_tests(G, max_size = 0, mode = mode)
    fit <- function() {
      tmp <- capture.output(opt <- .get_opt(mode)(tests,
        d = d, max_size = 0, V = V, mode = mode, cache = TRUE,
        cache_dir = cdir, gurobi_args = list(OutputFlag = 0)
      ))
      opt
    }
    opt1 <- fit() # writes the cache
    opt2 <- fit() # reads the cache
    expect_equal(opt2$optim$objval, opt1$optim$objval, label = mode)
    expect_equal(opt1$optim$objval, 0, tolerance = 1e-6, label = mode)
    unlink(cdir, recursive = TRUE)
    tmp <- capture.output(res <- suppressWarnings(learn_graph(dd,
      max_size = 0, mode = mode, comets = FALSE, test = "gaussCItest",
      cache = FALSE, gurobi_args = list(OutputFlag = 0)
    )))
    expect_false(is.null(res$graph$graph), label = mode)
  }
})

test_that("learn_graph() supports mode = 'dg-dcon' and rejects unknown modes", {
  skip_if_no_gurobi()
  set.seed(1)
  n <- 300
  E <- rnorm(n)
  X <- E + rnorm(n)
  Y <- X + rnorm(n)
  dd <- data.frame(E = E, X = X, Y = Y)
  expect_warning(tmp <- capture.output(res <- learn_graph(dd,
    mode = "dg-dcon", comets = FALSE, test = "gaussCItest", cache = FALSE,
    gurobi_args = list(OutputFlag = 0)
  )), "d-separation")
  expect_equal(sum(res$graph$graph[, ]), 2)
  expect_null(res$computed)
  expect_error(learn_graph(dd, mode = "pdag", comets = FALSE,
    test = "gaussCItest"), "`mode` must be one of")
})

test_that("Optimizers fail early with an informative error without gurobi", {
  local_mocked_bindings(.has_gurobi = function() FALSE)
  V <- letters[1:3]
  for (mode in .modes) {
    G <- .generate_random_graph(3, V = V, mode = mode, prob = 0.4)
    tests <- .compute_oracle_tests(G, mode = mode)
    expect_error(.get_opt(mode)(tests, d = 3, V = V, mode = mode,
      cache = FALSE), "Package `gurobi` is required", label = mode)
  }
})
