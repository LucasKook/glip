test_that("Learning graph works", {
  set.seed(12)
  dd <- data.frame(X = rnorm(100), Y = rnorm(100), Z = rnorm(100))
  tmp <- capture.output(opt <- learn_graph(dd))
  expect_equal(sum(opt$computed), 0)
})

test_that("Falsifying graph works", {
  set.seed(12)
  dd <- data.frame(X = rnorm(100), Y = rnorm(100), Z = rnorm(100))
  tmp <- capture.output(opt <- learn_graph(dd))
  tmp <- capture.output(out <- falsify_graph(opt$graph$graph, dd))
  expect_true(all(p.adjust(out$p.value) >= 0.05))
})

test_that("Empty graph is feasible", {
  sapply(3:5, \(d) {
    V <- letters[1:d]
    sapply(c("dag", "chain", "admg", "dagdcon"), \(mode) {
      G <- .generate_random_graph(d, V = V, mode = mode, prob = 0)
      if (is.list(G)) {
        G$M1[] <- 0
        G$M2[] <- 0
      } else {
        G[] <- 0
      }
      tests <- .compute_oracle_tests(G, mode = mode)
      tmp <- capture.output(opt <- .get_opt(mode)(tests, d = d, max_size = d - 2, V = V,
        cache = TRUE, gurobi_args = list(Threads = 7), mode = mode))
      lgr <- .compute_graphical_representation(opt$graph, mode)
      expect_equal(sum(lgr), 0)
    })
  })
})
