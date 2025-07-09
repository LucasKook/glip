test_that("Check equivalence function works", {
  G1 <- G2 <- G3 <- create_dmg(4, prob = 0, M2prob = 0, diag = FALSE)
  G1$M1["a", "b"] <- G1$M1["b", "c"] <- 1
  G2$M1["c", "b"] <- G2$M1["b", "a"] <- 1
  G3$M1["c", "b"] <- G3$M1["a", "b"] <- 1
  sapply(seq_len(2), \(m) {
    sapply(c("dag", "chain", "admg", "dagdcon"), \(mode) {
      if (mode != "admg") {
        G1 <- G1$M1
        G2 <- G2$M1
        G3 <- G3$M1
      }
      expect_true(suppressWarnings(.check_equivalence(G1, G2, m, mode)))
      expect_false(suppressWarnings(.check_equivalence(G1, G3, m, mode)))
      expect_false(suppressWarnings(.check_equivalence(G2, G3, m, mode)))
    })
  })
})

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
      lgr <- .compute_graphical_representation(opt$graph, d - 2, mode)
      expect_equal(sum(lgr), 0)
    })
  })
})
