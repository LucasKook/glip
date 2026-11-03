all_modes <- c("dag", "chain", "admg", "dag-dcon", "dg", "dmg", "dag-dc", "admg-dc", "dg-dc", "dmg-dc", "chain-dcon")

test_that("Check equivalence function works", {
  G1 <- G2 <- G3 <- create_dmg(4, prob = 0, M2prob = 0, diag = FALSE)
  G1$M1["a", "b"] <- G1$M1["b", "c"] <- 1
  G2$M1["c", "b"] <- G2$M1["b", "a"] <- 1
  G3$M1["c", "b"] <- G3$M1["a", "b"] <- 1
  sapply(seq_len(2), \(m) {
    sapply(all_modes, \(mode) {
      if (!grepl("dmg", mode)) {
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
  skip_if_no_gurobi()
  set.seed(12)
  dd <- data.frame(X = rnorm(100), Y = rnorm(100), Z = rnorm(100))
  tmp <- capture.output(opt <- learn_graph(dd))
  expect_equal(sum(opt$computed), 0)
  skip_if_not_installed("bnlearn")
  tmp <- capture.output(opt <- learn_graph(dd, comets = FALSE, test = "zf"))
  expect_equal(sum(opt$computed), 0)
  tmp <- capture.output(opt <- learn_graph(dd, comets = FALSE, test = "cor"))
  expect_equal(sum(opt$computed), 0)
})

test_that("Falsifying graph works", {
  skip_if_no_gurobi()
  set.seed(12)
  dd <- data.frame(X = rnorm(100), Y = rnorm(100), Z = rnorm(100))
  tmp <- capture.output(opt <- learn_graph(dd))
  tmp <- capture.output(out <- falsify_graph(opt$graph$graph, dd))
  expect_true(all(p.adjust(out$p.value) >= 0.05))
})

test_that("Empty graph is feasible", {
  skip_if_no_gurobi()
  sapply(3:5, \(d) {
    V <- letters[1:d]
    sapply(all_modes, \(mode) {
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
#       sapply(c("dag", "admg", "dag-dcon"), \(mode) {
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

test_that("Background knowledge and edge penalty are enforced in dcon_lean_optim", {
  skip_if_no_gurobi()
  set.seed(1)
  n <- 500
  E <- rnorm(n)
  X <- E + rnorm(n)
  Y <- X + rnorm(n)
  dd <- data.frame(E = E, X = X, Y = Y)
  V <- colnames(dd)
  tests <- learn_graph(dd, comets = FALSE, test = "gaussCItest",
    return_tests_only = TRUE)
  forb <- matrix(0, 3, 3, dimnames = list(V, V))
  forb[, "E"] <- 1
  forb["E", "Y"] <- 1
  tmp <- capture.output(opt <- dcon_lean_optim(tests, d = 3, V = V,
    forbidden = forb, edge_penalty = 1e-3, cache = FALSE,
    gurobi_args = list(OutputFlag = 0)))
  G <- opt$graph[, ]
  expect_true(all(G[forb == 1] == 0))
  expect_equal(sum(G), 2)
  expect_equal(G["E", "X"], 1)
  expect_equal(G["X", "Y"], 1)
  req <- matrix(0, 3, 3, dimnames = list(V, V))
  req["Y", "X"] <- 1
  tmp <- capture.output(opt <- dcon_lean_optim(tests, d = 3, V = V,
    forbidden = forb, required = req, cache = FALSE,
    gurobi_args = list(OutputFlag = 0)))
  expect_equal(opt$graph["Y", "X"], 1)
})

test_that("Linear edge constraints are enforced in dcon_lean_optim", {
  skip_if_no_gurobi()
  set.seed(2)
  n <- 500
  E <- rnorm(n)
  Y <- rnorm(n)
  X <- E + Y + rnorm(n)
  dd <- data.frame(E = E, X = X, Y = Y)
  V <- colnames(dd)
  tests <- learn_graph(dd, comets = FALSE, test = "gaussCItest",
    return_tests_only = TRUE)
  fit <- function(cons) {
    tmp <- capture.output(opt <- dcon_lean_optim(tests, d = 3, V = V,
      edge_constraints = cons, edge_penalty = 1e-3, cache = FALSE,
      gurobi_args = list(OutputFlag = 0)))
    opt
  }
  Z <- matrix(0, 3, 3, dimnames = list(V, V))
  G0 <- fit(NULL)$graph[, ]
  expect_equal(G0["E", "X"] + G0["Y", "X"], 2)
  ### No common child of E and Y: the collider is excluded, fit gets worse
  cf <- Z
  cf["E", "X"] <- cf["Y", "X"] <- 1
  opt1 <- fit(list(list(coef = cf, sense = "<=", rhs = 1)))
  G1 <- opt1$graph[, ]
  expect_lte(G1["E", "X"] + G1["Y", "X"], 1)
  expect_gte(floor(opt1$optim$objval), 1)
  ### Require E and Y to be adjacent
  cf <- Z
  cf["E", "Y"] <- cf["Y", "E"] <- 1
  G2 <- fit(list(list(coef = cf, sense = ">=", rhs = 1)))$graph[, ]
  expect_equal(G2["E", "Y"] + G2["Y", "E"], 1)
})

test_that("Background knowledge is enforced in admg_cubic_optim", {
  skip_if_no_gurobi()
  set.seed(3)
  n <- 1000
  H <- rnorm(n)
  E <- rnorm(n)
  X <- E + H + rnorm(n)
  Y <- X + H + rnorm(n)
  dd <- data.frame(E = E, X = X, Y = Y)
  V <- colnames(dd)
  tests <- learn_graph(dd, comets = FALSE, test = "gaussCItest",
    return_tests_only = TRUE)
  Z <- matrix(0, 3, 3, dimnames = list(V, V))
  fit <- function(...) {
    tmp <- capture.output(opt <- admg_cubic_optim(tests, d = 3, V = V,
      mode = "admg", edge_penalty = 1e-3, cache = FALSE,
      gurobi_args = list(OutputFlag = 0), ...))
    opt
  }
  ### i is an ancestor of j in the directed part M1
  is_anc <- function(M1, i, j) {
    R <- M1
    for (k in seq_len(NROW(M1))) R <- 1 * ((R + R %*% M1) > 0)
    R[i, j] == 1
  }
  ### E and Y not adjacent, E not a descendant of Y, no bidirected edges at E
  Fd <- Z
  Fd["E", "Y"] <- Fd["Y", "E"] <- 1
  Fb <- Z
  Fb["E", ] <- Fb[, "E"] <- 1
  diag(Fb) <- 0
  Fa <- Z
  Fa["Y", "E"] <- 1
  opt <- fit(forbidden = Fd, forbidden_bi = Fb, forbidden_anc = Fa)
  G <- opt$graph
  expect_true(all(G$M1[Fd == 1] == 0))
  expect_true(all(G$M2[Fb == 1] == 0))
  expect_false(is_anc(G$M1, "Y", "E"))
  ### Forbidding X as an ancestor of Y is respected and costs mismatches
  Fa2 <- Fa
  Fa2["X", "Y"] <- 1
  opt2 <- fit(forbidden = Fd, forbidden_bi = Fb, forbidden_anc = Fa2)
  expect_false(is_anc(opt2$graph$M1, "X", "Y"))
  expect_gt(floor(opt2$optim$objval), floor(opt$optim$objval))
  ### Required edge and required ancestral relation
  Rq <- Z
  Rq["Y", "X"] <- 1
  opt3 <- fit(forbidden = Fd, required = Rq)
  expect_equal(opt3$graph$M1["Y", "X"], 1)
  Ra <- Z
  Ra["E", "Y"] <- 1
  opt4 <- fit(forbidden = Fd, forbidden_bi = Fb, forbidden_anc = Fa,
    required_anc = Ra)
  expect_true(is_anc(opt4$graph$M1, "E", "Y"))
})
