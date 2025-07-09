list_separations <- function(G, max_size = NULL, mode, ...) {
  V <- .get_node_set(G)
  max_size <- min(max(1, length(V) - 2), max_size)
  tsts <- .list_tests_graph(V, max_size = max_size, ...)
  idx <- which(.compute_oracle_tests(G, max_size, mode)$p.value == 1)
  list(sets = tsts$sets[idx], formulas = tsts$formulas[idx])
}

falsify_graph <- function(
    G, data, max_size = NULL, mode = "dag", test = "gcm", test_args = NULL,
    parallel = FALSE, ncores = NULL, ...) {
  V <- colnames(G)
  max_size <- min(max(1, length(V) - 2), max_size)

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

  sets <- to_test$sets
  fml <- to_test$formulas

  pb <- utils::txtProgressBar(min = 0, max = length(fml), style = 3, width = 60)
  res <- my_apply(seq_along(fml), \(iter) {
    utils::setTxtProgressBar(pb, iter)
    res <- do.call("comets", c(list(
      formula = fml[[iter]], data = data, test = test
    ), test_args))
    data.frame(
      X = sets[[iter]][["X"]],
      Y = sets[[iter]][["Y"]],
      Z = paste0(sets[[iter]][["Z"]], collapse = ","),
      formula = paste0(deparse(fml[[iter]]), collapse = ""),
      p.value = res$p.value, weight = 1
    )
  }) |> do.call("rbind", args = _)

  if (parallel) {
    plan(sequential)
  }

  res
}
