### Complete warm starts of the min-length encodings (.lean_full_start())

### Build the model of a min-length optimizer without solving it
.capture_lean_model <- function(fun, ...) {
  captured <- NULL
  local_mocked_bindings(
    .has_gurobi = function() TRUE,
    .gurobi_solve = function(model, gurobi_args) {
      captured <<- model
      list(status = "STUB")
    }
  )
  utils::capture.output(fun(...))
  captured
}

.fs_cases <- function() {
  V <- paste0("X", 1:5)
  M <- matrix(0, 5, 5, dimnames = list(V, V))
  truth <- M
  truth[1, 2] <- truth[2, 3] <- truth[4, 3] <- truth[3, 5] <- 1
  start <- M
  start[1, 2] <- start[2, 4] <- start[5, 4] <- 1
  cyc <- start
  cyc[4, 2] <- 1 # 2-cycle for the cyclic modes
  B <- M
  B[1, 5] <- B[5, 1] <- 1
  ch <- start
  ch[3, 5] <- ch[5, 3] <- 1 # undirected edge
  tch <- truth
  tch[1, 4] <- tch[4, 1] <- 1
  list(
    list(fun = dcon_lean_optim, mode = "dag", truth = truth, start = start),
    list(fun = dcon_lean_optim, mode = "dg", truth = truth, start = cyc),
    list(fun = admg_lean_optim, mode = "admg-lean",
      truth = list(M1 = truth, M2 = B), start = list(M1 = start, M2 = B)),
    list(fun = admg_lean_optim, mode = "dmg-lean",
      truth = list(M1 = truth, M2 = B), start = list(M1 = cyc, M2 = M)),
    list(fun = chain_lean_optim, mode = "chain", truth = tch, start = ch)
  )
}

.fs_tests <- function(truth, ms, mode) {
  tests <- .compute_oracle_tests(truth, ms, mode)
  flip <- seq(1, nrow(tests), by = 7)
  tests$p.value[flip] <- 1 - tests$p.value[flip]
  tests$p.value <- ifelse(tests$p.value == 1, .5, .01)
  tests
}

.fs_model <- function(cs, tests, ms, ...) {
  args <- list(tests,
    d = 5, max_size = ms, V = paste0("X", 1:5),
    warmstart = cs$start, cache = FALSE, ...
  )
  if (cs$mode != "chain") args$mode <- cs$mode
  do.call(.capture_lean_model, c(list(cs$fun), args))
}

test_that("complete start satisfies all constraints and has the oracle loss", {
  for (cs in .fs_cases()) {
    for (ms in 0:3) {
      tests <- .fs_tests(cs$truth, ms, cs$mode)
      m <- .fs_model(cs, tests, ms)
      expect_false(anyNA(m$start), label = paste(cs$mode, ms))
      expect_true(.lean_check_start(m, m$start)$ok, label = paste(cs$mode, ms))
      o <- .compute_oracle_tests(cs$start, ms, cs$mode)
      mm <- merge(tests, o, by = c("X", "Y", "Z"))
      loss <- sum((mm$p.value.x <= .05) != (mm$p.value.y == 0))
      expect_equal(sum(m$obj * m$start), loss, label = paste(cs$mode, ms))
    }
  }
})

test_that("option glip.lean_full_start = FALSE keeps the edge values only", {
  withr_opt <- options(glip.lean_full_start = FALSE)
  on.exit(options(withr_opt))
  cs <- .fs_cases()[[1]]
  m <- .fs_model(cs, .fs_tests(cs$truth, 1, cs$mode), 1)
  expect_true(anyNA(m$start))
  full <- .lean_full_start(m)
  expect_false(anyNA(full))
  expect_true(.lean_check_start(m, full)$ok)
})

test_that("a start violating the knowledge keeps the partial start", {
  cs <- .fs_cases()[[1]]
  V <- paste0("X", 1:5)
  ra <- matrix(0, 5, 5, dimnames = list(V, V))
  ra[5, 1] <- 1 # required ancestry not in the start
  m <- .fs_model(cs, .fs_tests(cs$truth, 1, cs$mode), 1, required_anc = ra)
  expect_true(anyNA(m$start))
  expect_null(.lean_full_start(m))
  expect_type(.lean_full_start(m, why = TRUE), "character")
})

test_that("Gurobi accepts the complete start unchanged", {
  skip_if_no_gurobi()
  for (cs in .fs_cases()) {
    tests <- .fs_tests(cs$truth, 2, cs$mode)
    m <- .fs_model(cs, tests, 2)
    s <- gurobi::gurobi(m, list(
      Threads = 1, OutputFlag = 0, SolutionLimit = 1, StartNodeLimit = 0
    ))
    expect_equal(s$objval, sum(m$obj * m$start), label = cs$mode)
  }
})
