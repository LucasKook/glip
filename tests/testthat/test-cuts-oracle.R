### The separation automata of the cut encoding (Python module glipcut)
### against the R oracles: is_m_separated() for DAGs and ADMGs, csep() for
### LWF chain graphs; full banks on random graphs

test_that("cut encoding automata agree with the R separation oracles", {
  skip_if_no_glipcut_python()
  set.seed(20261116)
  for (mode in c("dag", "admg-dc", "chain")) {
    class <- .cuts_class(mode)
    for (d in 3:5) {
      for (r in 1:4) {
        V <- letters[seq_len(d)]
        G <- .generate_random_graph(d, V = V, mode = mode, prob = c(.2, .4, .6, .8)[r])
        tests <- .compute_oracle_tests(G, max_size = d - 2, mode = mode)
        res <- .glipcut_connected(class, G, tests, V)
        expect_true(res$valid, label = paste(mode, d, r))
        expect_identical(as.numeric(res$connected), 1 - tests$p.value,
          label = paste(mode, d, r))
      }
    }
  }
})

test_that("chain graph automaton handles long undirected sections", {
  skip_if_no_glipcut_python()
  V <- letters[1:6]
  G <- matrix(0, 6, 6, dimnames = list(V, V))
  ### a -> b - c - d <- e, d -> f: head-to-head section {b, c, d}
  G["a", "b"] <- G["e", "d"] <- G["d", "f"] <- 1
  G["b", "c"] <- G["c", "b"] <- G["c", "d"] <- G["d", "c"] <- 1
  tests <- .compute_oracle_tests(G, max_size = 4, mode = "chain")
  res <- .glipcut_connected("chain", G, tests, V)
  expect_true(res$valid)
  expect_identical(as.numeric(res$connected), 1 - tests$p.value)
  ### a and e are separated marginally and connected given c (in the section)
  expect_true(csep(G, "a", "e"))
  expect_false(csep(G, "a", "e", "c"))
})

test_that("semidirected cycles are not valid chain graphs", {
  skip_if_no_glipcut_python()
  V <- letters[1:3]
  G <- matrix(0, 3, 3, dimnames = list(V, V))
  G["a", "b"] <- 1
  G["b", "c"] <- G["c", "b"] <- G["c", "a"] <- G["a", "c"] <- 1
  tests <- .compute_oracle_tests(matrix(0, 3, 3, dimnames = list(V, V)), 1, "chain")
  expect_false(.glipcut_connected("chain", G, tests, V)$valid)
})
