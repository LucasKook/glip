### Generate random graph
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

### Sample from random graph
rgraph <- function(graphs, n, ...) {
  dat <- pcalg::rmvDAG(n = n, dag = graphs$graph, ...)
  dat[, graphs$O]
}

### Lookup-table based CI test for PC/FCI
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

### Evaluation metrics
shd <- function(G1, G2, ...) {
  .check_graphs(G1, G2)
  ### Hamming distance
  SID::hammingDist(G1, G2, ...)
}

sid <- function(G1, G2, ...) {
  .check_graphs(G1, G2)
  ### Structural intervention distance
  sid <- SID::structIntervDist(G1, G2, ...)
  structure(sid$sid, full_output = sid)
}

prf1 <- function(G1, G2) {
  .check_graphs(G1, G2)
  d <- nrow(G1)
  res <- lapply(seq_len(d), \(k) {
    rbind(
      data.frame(
        node = k,
        which = "head",
        .classification_metrics(G1[, k], G2[, k])
      ),
      data.frame(
        node = k,
        which = "tail",
        .classification_metrics(G1[k, ], G2[k, ])
      )
    )
  }) |> do.call("rbind", args = _)
  res
  res |>
    dplyr::group_by(which) |>
    dplyr::summarize_at(c("precision", "recall", "f1"), mean, na.rm = TRUE)
}

### Helpers
.check_graphs <- function(G1, G2) {
  stopifnot(is.matrix(G1) & is.matrix(G2))
  stopifnot(nrow(G1) == nrow(G2))
  stopifnot(rownames(G1) == rownames(G2))
}

.classification_metrics <- function(true_vec, pred_vec) {
  tp <- sum(true_vec == 1 & pred_vec == 1)
  fn <- sum(true_vec == 1 & pred_vec == 0)
  fp <- sum(true_vec == 0 & pred_vec == 1)

  precision <- ifelse(tp + fp == 0, NA, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, NA, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, NA, 2 * precision * recall / (precision + recall))

  data.frame(precision = precision, recall = recall, f1 = f1)
}
