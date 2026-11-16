test_that("classic K4cd cannot connect an isolated endpoint", {
  skip_if_no_gurobi()
  V <- paste0("X", 1:5)
  M1 <- M2 <- matrix(0, 5, 5, dimnames = list(V, V))
  M1[2, 3] <- M1[5, 4] <- 1
  M2[2, 3] <- M2[3, 2] <- M2[3, 5] <- M2[5, 3] <- 1
  G <- list(M1 = M1, M2 = M2)
  expect_true(is_m_separated(G, "X4", "X1", c("X3", "X5")))
  tests <- .compute_oracle_tests(G, max_size = 3, mode = "admg-dc")
  F <- 1 - M1; FB <- 1 - M2; diag(F) <- diag(FB) <- 0
  for (mode in c("admg-dc", "dmg-dc")) {
    cache <- tempfile("classic-k4cd-")
    dir.create(cache)
    # All pre-v2 components exist but must not be read.
    stem <- paste0("dim5max3", mode)
    writeLines("stale matrix", file.path(cache, paste0("lhs", stem, ".mtx")))
    for (prefix in c("rhs", "P1", "N1", "O1")) {
      saveRDS(NULL, file.path(cache, paste0(prefix, stem, ".rds")))
    }
    run <- function(tests) {
      capture.output(fit <- admg_optim(tests, d = 5, V = V, max_size = 3,
        mode = mode, warmstart = G, required = M1, forbidden = F,
        required_bi = M2, forbidden_bi = FB, cache_dir = cache,
        gurobi_args = list(Threads = 1, Seed = 1, TimeLimit = 30,
                           OutputFlag = 0, MIPGap = 0, MIPGapAbs = 0)))
      fit
    }
    # Truth and opposite banks ensure the fixed graph has no adjustable,
    # spurious connection decisions, on both cache miss and cache hit.
    fit <- run(tests)
    expect_identical(fit$optim$status, "OPTIMAL")
    expect_equal(fit$optim$objval, 0, tolerance = 1e-7)
    expect_equal(fit$tests$dcon, 1 - fit$tests$p.value, tolerance = 1e-7)
    opposite <- tests; opposite$p.value <- 1 - opposite$p.value
    fit <- run(opposite)
    expect_identical(fit$optim$status, "OPTIMAL")
    expect_equal(fit$optim$objval, nrow(tests), tolerance = 1e-7)
    expect_equal(fit$tests$dcon, fit$tests$p.value, tolerance = 1e-7)
    expect_true(file.exists(file.path(cache, paste0("lhs", stem, "-v2.mtx"))))
  }
})
