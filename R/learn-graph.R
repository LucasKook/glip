#' @export
learn_graph <- function(
    data, max_size = NULL, mode = "dag", test = "gcm", naive = FALSE,
    parallel = FALSE, ncores = NULL, alpha = 0.05, trafo = \(x) as.numeric(x <= alpha),
    weight_type = "const", warmstart = NULL, edgehints = NULL,
    gurobi_args = list(), test_args = NULL, return_tests_only = FALSE,
    verbose = FALSE, cache = TRUE, comets = TRUE, ...) {
  if (mode %in% c("dmg", "dg")) {
    warning("Using `mode = 'dmg'` or `mode = 'dg'` relies on d-separation
      and thus implicitly assumes a linear Gaussian SCM.")
  }

  ### Pre-process
  vars <- colnames(data)
  max_size <- min(max(1, length(vars) - 2), max_size)

  if (parallel) {
    nc <- min(parallel::detectCores() - 1, 15, ncores)
    plan(multisession, workers = nc)
    my_apply <- \(...) future_lapply(..., future.seed = TRUE)
  } else {
    my_apply <- lapply
  }

  ### Set naive = TRUE to run tests under causal sufficiency assumption
  all_tests <- .list_tests_graph(vars, max_size, naive = naive)
  res <- .run_tests(all_tests, data, test, test_args, my_apply, comets, ...)

  if (parallel) {
    plan(sequential)
  }

  if (return_tests_only) {
    return(res)
  }

  if (is.null(warmstart) & mode != "chain") {
    suff <- list(tests = res, V = vars)
    cit <- lookup_ci
    if (mode %in% c("dag", "dg")) {
      pcres <- pcalg::pc(suff, cit, labels = vars, alpha = alpha)
      PC <- as(pcres@graph, "matrix")
      warmstart <- .ess_to_dag(PC)
      edgehints <- 1 * (PC != 0)
    } else {
      fcires <- pcalg::fciPlus(suff, cit,
        labels = vars, alpha = alpha, selectionBias = FALSE, verbose = FALSE
      )
      FCI <- as(fcires@amat, "matrix")
      warmstart <- .pag_to_admg(FCI)
      edgehints <- 1 * (FCI != 0)
    }
  }

  opt <- .get_opt(mode)
  graph <- opt(res,
    d = length(vars), max_size = max_size,
    V = vars, trafo = trafo, weight_type = weight_type,
    warmstart = warmstart, edgehints = edgehints,
    gurobi_args = gurobi_args, verbose = verbose,
    cache = cache, mode = mode, ...
  )

  out <- .compute_graphical_representation(graph$graph, max_size, mode)

  structure(list(tests = res, graph = graph, computed = out),
    class = "learned_graph", vars = vars, max_size = max_size, test = test
  )
}

.list_tests_graph <- function(vars, max_size, naive = FALSE, ...) {
  pairs <- utils::combn(vars, 2)
  if (naive) {
    sets <- apply(pairs, 2, \(pair) {
      out <- utils::combn(setdiff(vars, pair), length(vars) - 2)
      apply(out, 2, \(x) {
        list(X = pair[1], Y = pair[2], Z = x)
      })
    }, simplify = FALSE) |> unlist(recursive = FALSE)
  } else {
    sets <- apply(pairs, 2, \(pair) {
      lapply(0:max_size, \(x) {
        out <- utils::combn(setdiff(vars, pair), x)
        apply(out, 2, \(x) {
          list(X = pair[1], Y = pair[2], Z = x)
        })
      })
    }, simplify = FALSE) |>
      unlist(recursive = FALSE) |>
      unlist(recursive = FALSE)
  }
  fml <- lapply(sets, \(set) {
    .to_formula_graph(set[["Y"]], set[["X"]], set[["Z"]], ...)
  })
  list(sets = sets, formulas = fml)
}

.to_formula_graph <- function(X, Y, Z) {
  if (identical(Z, character(0))) {
    Z <- "1"
  }
  as.formula(
    paste0(Y, "~", paste0(X, collapse = "+"), "|", paste0(Z, collapse = "+"))
  )
}

#' @exportS3Method print graphopt
print.graphopt <- function(x, ...) {
  print(x$graph)
}

#' @exportS3Method print learned_graph
print.learned_graph <- function(x, print_tests = FALSE, ...) {
  if (print_tests) {
    cat("\nTests:\n")
    print(x$tests[, c("X", "Y", "Z", "p.value")])
  }
  cat("\nLearned graph:\n")
  print(x$graph)
}
