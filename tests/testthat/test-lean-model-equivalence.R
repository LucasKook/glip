# Compare complete in-memory models, including bounds, objective, hints and
# appended edge-constraint rows. No solver or Gurobi license is required.
test_that("compiled and reference encodings produce equivalent complete models", {
  ns <- asNamespace("glip")
  original_rows <- get(".bk_rows", ns)
  original_pairs <- function(d, n_C) {
    out <- data.frame(); counter <- 1
    for (i in seq_len(d)) for (j in seq_len(i - 1)) for (C in seq_len(n_C)) {
      out <- rbind(out, data.frame(i = i, j = j, C = C, idx = counter))
      counter <- counter + 1
    }
    out
  }
  no_solver <- function(...) invisible(NULL)
  done <- function() stop(structure(list(message = "model captured", call = NULL),
                                    class = c("model_captured", "error", "condition")))
  set.seed(20261105)
  for (d in 4:5) for (mode in c("dag", "dg", "admg-lean", "dmg-lean", "chain")) {
    max_size <- d - 2L
    V <- sample(LETTERS[seq_len(d)])
    cases <- list()
    for (i in seq_len(d)) for (j in seq_len(i - 1L)) {
      rest <- setdiff(seq_len(d), c(i, j))
      for (k in 0:max_size) for (C in combn(rest, k, simplify = FALSE)) {
        flip <- sample(c(FALSE, TRUE), 1)
        cases[[length(cases) + 1L]] <- data.frame(
          X = V[if (flip) j else i], Y = V[if (flip) i else j],
          Z = paste(V[C], collapse = ","), p.value = runif(1), size = k)
      }
    }
    tests <- do.call(rbind, cases)
    tests <- tests[sample(nrow(tests)), , drop = FALSE]
    forbidden <- matrix(0, d, d, dimnames = list(V, V))
    required <- forbidden; coef <- forbidden
    forbidden[1, 2] <- 1; required[2, 3] <- 1; coef[1, 3] <- 1
    kind <- if (mode %in% c("dag", "dg")) "dcon" else if (mode == "chain") "chain" else "admg"
    fit <- get(paste0(kind, "_lean_optim"), ns)
    args <- list(tests = tests, d = d, max_size = max_size, V = V, cache = FALSE,
                 weight_type = "log", forbidden = forbidden, required = required,
                 edge_penalty = 0.01,
                 edge_constraints = list(coef = coef, sense = "<=", rhs = 1))
    if (kind != "chain") args$mode <- mode
    captured <- new.env(parent = emptyenv())
    capture_rows <- function(model, ec) {
      captured$model <- original_rows(model, ec)
      done()
    }
    run <- function() {
      ok <- tryCatch({ capture.output(do.call(fit, args)); FALSE },
                     model_captured = function(e) TRUE)
      expect_true(ok, info = paste(mode, d))
      captured$model
    }
    actual <- testthat::with_mocked_bindings(run(),
      .bk_rows = capture_rows, .require_gurobi = no_solver, .package = "glip")
    reference <- if (kind == "dcon") {
      testthat::with_mocked_bindings(run(), .bk_rows = capture_rows,
        .require_gurobi = no_solver, .dcon_lean_indices_c = get(".dcon_lean_indices_r", ns),
        .dcon_lean_build_c = get(".dcon_lean_build_r", ns), .package = "glip")
    } else if (kind == "admg") {
      testthat::with_mocked_bindings(run(), .bk_rows = capture_rows,
        .require_gurobi = no_solver, .dcon_lean_indices_c = get(".dcon_lean_indices_r", ns),
        .lean_pair_conditioning_indices = original_pairs,
        .admg_lean_build_c = get(".admg_lean_build_r", ns), .package = "glip")
    } else {
      testthat::with_mocked_bindings(run(), .bk_rows = capture_rows,
        .require_gurobi = no_solver, .dcon_lean_indices_c = get(".dcon_lean_indices_r", ns),
        .chain_lean_build_c = get(".chain_lean_build_r", ns), .package = "glip")
    }
    expect_equal(actual, reference, tolerance = 0, info = paste(mode, d))
  }
})
