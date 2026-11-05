test_that("lean ADMG, DMG and chain caches match the R reference", {
  skip_if_no_gurobi()
  ns <- asNamespace("glip")
  set.seed(482)
  for (mode in c("admg-lean", "dmg-lean", "chain")) {
    for (d in 3:4) for (max_size in 0:(d - 2)) {
      V <- letters[seq_len(d)]
      graph <- get(".generate_random_graph", ns)(d, V = V, mode = mode, prob = 0.3)
      tests <- get(".compute_oracle_tests", ns)(graph, max_size = max_size, mode = mode)
      kind <- if (mode == "chain") "chain" else "admg"
      fit <- get(paste0(kind, "_lean_optim"), ns)
      args <- list(tests = tests, d = d, V = V, max_size = max_size,
                   cache = TRUE, gurobi_args = list(OutputFlag = 0, Threads = 1))
      if (kind == "admg") args$mode <- mode
      paths <- file.path(tempdir(), paste0("lean-", kind, "-", mode, "-", d,
                                           "-", max_size, c("-c", "-r")))
      unlink(paths, recursive = TRUE)
      run <- function(path) {
        capture.output(value <- do.call(fit, c(args, list(cache_dir = path))))
        value
      }
      fast <- run(paths[1])
      if (kind == "chain") {
        reference <- testthat::with_mocked_bindings(
          run(paths[2]), .chain_lean_build_c = get(".chain_lean_build_r", ns),
          .package = "glip")
      } else {
        reference <- testthat::with_mocked_bindings(
          run(paths[2]), .admg_lean_build_c = get(".admg_lean_build_r", ns),
          .package = "glip")
      }
      files <- list.files(paths[1])
      expect_identical(files, list.files(paths[2]))
      for (file in files) {
        read <- if (grepl("[.]mtx$", file)) Matrix::readMM else readRDS
        expect_identical(read(file.path(paths[1], file)), read(file.path(paths[2], file)),
                     info = paste(mode, d, max_size, file))
      }
      expect_equal(fast$optim$objval, reference$optim$objval)
      expect_equal(run(paths[2])$optim$objval, reference$optim$objval)
      unlink(paths, recursive = TRUE)
    }
  }
})

test_that("lean caches match the R reference without a solver", {
  ns <- asNamespace("glip")
  tests <- data.frame(X = character(), Y = character(), Z = character(),
                      p.value = numeric(), size = numeric())
  stop_after_cache <- function(...) {
    stop(structure(list(message = "cache built", call = NULL),
                   class = c("cache_built", "error", "condition")))
  }
  no_solver <- function(...) invisible(NULL)
  for (mode in c("admg-lean", "dmg-lean", "chain")) {
    for (spec in list(c(2L, 0L), c(3L, 0L), c(4L, 2L), c(5L, 3L))) {
      d <- spec[1]; max_size <- spec[2]
      kind <- if (mode == "chain") "chain" else "admg"
      fit <- get(paste0(kind, "_lean_optim"), ns)
      args <- list(tests = tests, d = d, V = letters[seq_len(d)],
                   max_size = max_size, cache = TRUE)
      if (kind == "admg") args$mode <- mode
      paths <- file.path(tempdir(), paste0("no-solver-", mode, "-", d, "-", max_size, c("-c", "-r")))
      unlink(paths, recursive = TRUE)
      run <- function(path) {
        ok <- tryCatch({
          capture.output(do.call(fit, c(args, list(cache_dir = path))))
          FALSE
        }, cache_built = function(e) TRUE)
        expect_true(ok)
      }
      testthat::with_mocked_bindings({
        run(paths[1])
        if (kind == "chain") {
          testthat::with_mocked_bindings(run(paths[2]),
            .chain_lean_build_c = get(".chain_lean_build_r", ns), .package = "glip")
        } else {
          testthat::with_mocked_bindings(run(paths[2]),
            .admg_lean_build_c = get(".admg_lean_build_r", ns), .package = "glip")
        }
      }, .bk_rows = stop_after_cache, .require_gurobi = no_solver, .package = "glip")
      files <- list.files(paths[1])
      expect_identical(files, list.files(paths[2]))
      for (file in files) {
        read <- if (grepl("[.]mtx$", file)) Matrix::readMM else readRDS
        expect_identical(read(file.path(paths[1], file)), read(file.path(paths[2], file)),
                     info = paste(mode, d, max_size, file))
      }
      unlink(paths, recursive = TRUE)
    }
  }
})
