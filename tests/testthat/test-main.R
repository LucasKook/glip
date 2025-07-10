test_that("Check equivalence function works", {
  G1 <- G2 <- G3 <- create_dmg(4, prob = 0, M2prob = 0, diag = FALSE)
  G1$M1["a", "b"] <- G1$M1["b", "c"] <- 1
  G2$M1["c", "b"] <- G2$M1["b", "a"] <- 1
  G3$M1["c", "b"] <- G3$M1["a", "b"] <- 1
  sapply(seq_len(2), \(m) {
    sapply(c("dag", "chain", "admg", "dagdcon", "dg", "dmg"), \(mode) {
      if (!mode %in% c("admg", "dmg")) {
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
    sapply(c("dag", "chain", "admg", "dagdcon", "dg", "dmg"), \(mode) {
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
      lgr <- opt$graph
      out <- if (is.list(lgr)) {
        sum(lgr$M1, lgr$M2)
      } else {
        sum(lgr)
      }
      expect_equal(out, 0)
    })
  })
})

# test_that("msep agrees across implementations", {
#   library("ggm")
#   set.seed(12)
#   nsim <- 100
#   sapply(1:nsim, \(iter) {
#     sapply(5:10, \(d) {
#       V <- 1:d
#       sapply(c("dag", "admg", "dagdcon"), \(mode) {
#         print(d)
#         print(mode)
#         G <- .generate_random_graph(d, V = V, mode = mode, prob = 0.5)
#         print(G)
#         amat <- if (is.list(G)) {
#           1 * G$M1 + 100 * G$M2
#         } else {
#           G
#         }
#         C <- sample(setdiff(V, c(1, 2)), sample.int(d - 4, 1))
#         expect_equal(.check_msep(1, 2, 3, G), msep(amat, 1, 2, 3))
#       })
#     })
#   })
# })
