test_that("evaluation metrics work", {
  set.seed(12)
  G <- random_graph(d = 10)$DAG

  ### Essential graph output
  G1 <- G2 <- .dag2ess(G)
  cm <- prf1(G1, G2)
  expect_equal(cm$precision, c(1, 1))
  expect_equal(cm$recall, c(1, 1))
  expect_equal(cm$f1, c(1, 1))
  expect_equal(cm$fdr, c(0, 0))
  expect_equal(cm$mcc, c(1, 1))
  expect_equal(shd(G1, G2), 0)
  expect_equal(sep(G1[1:5, 1:5], G2[1:5, 1:5], "pdag", 1), 0)

  ### PAG output
  G3 <- G4 <- .admg2pag(.to_admg(G))
  cm2 <- prf1(G3, G4)
  expect_equal(shd(G3, G4), 0)
  expect_true(all(cm2$precision == 1))
  expect_true(all(cm2$recall == 1))
  expect_true(all(cm2$f1 == 1))
  expect_true(all(cm2$fdr == 0))
  expect_true(all(cm2$mcc == 1))
  expect_equal(sep(G3[1:5, 1:5], G4[1:5, 1:5], "mag", 1), 0)
})
