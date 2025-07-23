test_that("evaluation metrics", {
  set.seed(12)
  G1 <- G2 <- random_graph(d = 10)$DAG
  expect_equal(shd(G1, G2), 0)
  sid <- sid(G1, G2)
  attr(sid, "full_output") <- NULL
  expect_equal(sid, 0)
  cm <- prf1(G1, G2)
  expect_equal(cm$precision, c(1, 1))
  expect_equal(cm$recall, c(1, 1))
  expect_equal(cm$f1, c(1, 1))
})
