#' Falsify a Graph by Testing Its Implied Conditional Independencies
#'
#' Performs conditional independence tests on all d-separations implied by a
#' (causal) graph, allowing users to empirically falsify a proposed model
#' using given data. Returns a data frame of each tested d-separation and test
#' result.
#'
#' @param G An adjacency matrix or graph object (e.g., representing a DAG,
#'     ADMG, etc.) whose implied d-separations will be tested.
#' @param data A data.frame or matrix; each column should correspond to a
#'     variable in \code{G}.
#' @param max_size Integer. Maximum size of conditioning set to consider
#'     (defaults to `d - 2`).
#' @param mode Character. Graph mode (e.g., "dag", "admg", etc.).
#' @param test Character. Conditional independence test to use.
#' @param test_args Optional list. Additional arguments for the chosen test.
#' @param parallel Logical. If \code{TRUE}, tests are run in parallel.
#' @param ncores Integer. Number of cores to use if running in parallel.
#' @param comets Logical. If \code{TRUE} (default), uses the "comets" test implementation.
#' @param ... Additional arguments passed to test functions.
#'
#' @return A data.frame of tested (X, Y | Z) triples, size of conditioning set,
#'     test p-value, selected weight, and formula.
#'
#' @details
#' Only d-separations present in the graph (and testable at the given max_size)
#' are tested. If no such testable d-separations exist, only a message is
#' returned.
#'
#' @export
falsify_graph <- function(
    G, data, max_size = NULL, mode = "dag", test = "gcm", test_args = NULL,
    parallel = FALSE, ncores = NULL, comets = TRUE, ...) {
  d <- .get_size(G)
  max_size <- min(max(1, d - 2), max_size)

  ### List all d separations with tests
  to_test <- list_separations(G, max_size = max_size, mode)

  if (length(to_test$sets) == 0) {
    message("No testable d-separations at the provided max_size.")
    return(invisible(NULL))
  }

  if (parallel) {
    nc <- min(parallel::detectCores() - 1, 15, ncores)
    plan(multisession, workers = nc)
    my_apply <- \(...) future_lapply(..., future.seed = TRUE)
  } else {
    my_apply <- lapply
  }

  res <- .run_tests(to_test, data, test, test_args, my_apply, comets, ...)

  if (parallel) {
    plan(sequential)
  }

  res
}

.run_tests <- function(to_test, data, test, test_args, my_apply, comets, ...) {
  sets <- to_test$sets
  fml <- to_test$formulas

  pb <- utils::txtProgressBar(min = 0, max = length(fml), style = 3, width = 60)
  my_apply(seq_along(fml), \(iter) {
    utils::setTxtProgressBar(pb, iter)
    if (!comets) {
      res <- tryCatch(
        {
          x <- sets[[iter]]$X
          y <- sets[[iter]]$Y
          z <- sets[[iter]]$Z
          if (test == "gaussCItest") {
            pv <- pcalg::gaussCItest(x, y, z, list(C = cor(data), n = NROW(data)))
          } else {
            pv <- if (identical(z, character(0))) {
              bnlearn::ci.test(x, y, data = data, test = test)$p.value
            } else {
              bnlearn::ci.test(x, y, z, data, test = test)$p.value
            }
          }
          list(p.value = pv)
        },
        error = \(e) {
          warning("Package `bnlearn`'s `ci.test()` failed. Consider using `comets = TRUE`.")
          list(p.value = 1 - 2 * .Machine$double.eps)
        }
      )
    } else {
      res <- do.call("comets", c(list(
        formula = fml[[iter]], data = data, test = test
      ), test_args))
    }
    data.frame(
      X = sets[[iter]][["X"]],
      Y = sets[[iter]][["Y"]],
      Z = paste0(sets[[iter]][["Z"]], collapse = ","),
      size = length(sets[[iter]][["Z"]]),
      formula = paste0(deparse(fml[[iter]]), collapse = ""),
      p.value = res$p.value, weight = 1
    )
  }) |> do.call("rbind", args = _)
}

list_separations <- function(G, max_size = NULL, mode, ...) {
  V <- .get_node_set(G)
  max_size <- min(max(1, length(V) - 2), max_size)
  tsts <- .list_tests_graph(V, max_size = max_size, ...)
  idx <- which(.compute_oracle_tests(G, max_size, mode)$p.value == 1)
  list(sets = tsts$sets[idx], formulas = tsts$formulas[idx])
}
