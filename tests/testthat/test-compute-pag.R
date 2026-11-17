test_that("compute_pag works with more than 26 bidirected edges", {
  V <- paste0("X", 1:8)
  G <- list(
    M1 = matrix(0, 8, 8, dimnames = list(V, V)),
    M2 = matrix(1, 8, 8, dimnames = list(V, V))
  )
  diag(G$M2) <- 0
  G$M1["X8", c("X1", "X4", "X5")] <- 1
  G$M1["X1", "X4"] <- 1
  D <- invert_latent_projection(G)
  expect_equal(anyDuplicated(rownames(D$M1)), 0)
  expect_equal(NROW(D$M1), 8 + 28)
  pag <- compute_pag(G)
  expect_equal(dimnames(pag), list(V, V))
  ### Every pair is adjacent (bidirected edges everywhere)
  expect_true(all((pag != 0)[upper.tri(pag)]))
})

test_that("latent names stay distinct from observed names", {
  V <- c("n1", "n2", "X3")
  G <- list(
    M1 = matrix(0, 3, 3, dimnames = list(V, V)),
    M2 = matrix(0, 3, 3, dimnames = list(V, V))
  )
  G$M2["n1", "n2"] <- G$M2["n2", "n1"] <- 1
  G$M1["n2", "X3"] <- 1
  D <- invert_latent_projection(G)
  expect_equal(anyDuplicated(rownames(D$M1)), 0)
  expect_equal(rownames(D$M1)[1:3], V)
  expect_equal(dimnames(compute_pag(G)), list(V, V))
})

test_that("compute_pag is unchanged for at most 26 bidirected edges", {
  old_inv <- function(G) {
    V <- rownames(G$M1)
    d <- length(V)
    nb <- sum(G$M2) / 2
    if (nb == 0) {
      return(G)
    }
    VU <- c(V, paste0("n", letters[seq_len(nb)]))
    D0 <- create_dmg(length(VU),
      prob = 0, M2prob = 0,
      diag = FALSE, nodeNames = VU
    )
    D0$M1[V, V] <- G$M1
    bd <- which(G$M2 == 1, arr.ind = TRUE)
    bd <- bd[bd[, 1] < bd[, 2], , drop = FALSE]
    for (i in 1:nb) {
      D0$M1[d + i, bd[i, 1]] <- D0$M1[d + i, bd[i, 2]] <- 1
    }
    D0
  }
  set.seed(1)
  for (r in 1:20) {
    d <- sample(3:6, 1)
    V <- paste0("X", seq_len(d))
    M1 <- matrix(0, d, d, dimnames = list(V, V))
    M1[upper.tri(M1)] <- rbinom(d * (d - 1) / 2, 1, .4)
    M2 <- matrix(0, d, d, dimnames = list(V, V))
    M2[upper.tri(M2)] <- rbinom(d * (d - 1) / 2, 1, .3)
    M2 <- M2 + t(M2)
    G <- list(M1 = M1, M2 = M2)
    new <- suppressWarnings(compute_pag(G)) # bfs warnings for unconnected graphs
    old <- with_mocked_bindings(suppressWarnings(compute_pag(G)),
      invert_latent_projection = old_inv
    )
    expect_identical(old, new)
  }
})
