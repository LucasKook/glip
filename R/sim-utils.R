random_graph <- function(
    d = 3, prob = 0.5, lb = -1, ub = 1, admg_add = 2,
    mode = c("dag", "admg"), V = letters[1:d], V2 = letters[(d + 1):(d + admg_add)]) {
  mode <- match.arg(mode)
  O <- V
  if (mode == "admg") {
    V <- c(V, V2)
    O <- V[sort(sample.int(d + admg_add, d, FALSE))]
    d <- d + admg_add
  }
  fg <- pcalg::randomDAG(n = d, prob = prob, lB = lb, uB = ub, V = V)
  G <- suppressWarnings(1 * (as(fg, "matrix") != 0))
  mg <- marginalize_dag_to_admg(G, O)
  list(graph = fg, DAG = G, ADMG = mg, O = O)
}

rgraph <- function(graphs, n, ...) {
  dat <- pcalg::rmvDAG(n = n, dag = graphs$graph, ...)
  dat[, graphs$O]
}

lookup_ci <- function(x, y, S, suffstat) {
  V <- suffstat$V
  if (identical(S, integer(0))) {
    pv <- suffstat$tests |>
      dplyr::filter((X == V[x] & Y == V[y]) | (X == V[y] & Y == V[x]), Z == "") |>
      dplyr::pull(p.value)
    return(pv)
  }
  pv <- suffstat$tests |>
    dplyr::filter((X == V[x] & Y == V[y]) | (X == V[y] & Y == V[x]), Z == paste0(V[sort(S)], collapse = ",")) |>
    dplyr::pull(p.value)
  pv
}

### TODO: Evaluation metrics
