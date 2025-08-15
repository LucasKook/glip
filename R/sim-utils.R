### Generate random graph
random_graph <- function(
    d = 3, prob = 0.5, lb = 0, ub = 1, admg_add = 2, degree = 2, method = "er",
    par1 = NULL, par2 = NULL, mode = c("dag", "admg"), V = letters[1:d],
    V2 = letters[(d + 1):(d + admg_add)], dag_gen = c("randDAG", "randomDAG")) {
  mode <- match.arg(mode)
  dag_gen <- match.arg(dag_gen)
  O <- V
  ### For marginalization of DAG to ADMG
  if (mode == "admg") {
    V <- c(V, V2) # expand node set
    O <- V[sort(sample.int(d + admg_add, d, FALSE))] # randomly choose unobservables
    d <- d + admg_add # Expand dimension
  }
  if (dag_gen == "randomDAG") {
    ### Already topologically sorted
    fg <- pcalg::randomDAG(n = d, prob = prob, lB = lb, uB = ub, V = V)
  } else {
    fg <- pcalg::randDAG(
      n = d, d = degree, method = method, par1 = par1, par2 = par2,
      weighted = TRUE, wFUN = list(runif, min = lb, max = ub)
    )
    ### Needs to be topologically sorted
    sorted <- igraph::topo_sort(
      igraph::graph_from_graphnel(fg),
      mode = "out"
    )
    wadj <- as(fg, "matrix")[sorted, ][, sorted]
    dimnames(wadj) <- list(V, V)
    ### Create new graphNEL and add edges manually with prev. simulated weight
    fg <- new("graphNEL", nodes = V, edgemode = "directed")
    for (i in seq_along(V)) {
      for (j in seq_along(V)) {
        w <- wadj[i, j]
        if (w != 0) {
          from <- V[i]
          to <- V[j]
          fg <- graph::addEdge(from, to, fg, w)
        }
      }
    }
  }
  ### Adjacency matrix and marginalization
  G <- suppressWarnings(1 * (as(fg, "matrix") != 0))
  mg <- marginalize_dag_to_admg(G, O)
  ### Return
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
  mean(1 * (G1 != G2))
}

prf1 <- function(predicted, groundtruth, ...) {
  UseMethod("prf1")
}

#' @exportS3Method prf1 default
prf1.default <- function(predicted, groundtruth, summarize = TRUE) {
  G1 <- predicted
  G2 <- groundtruth
  .check_graphs(G1, G2)
  d <- nrow(G1)
  V <- .get_node_set(G1)
  res <- lapply(seq_len(d), \(k) {
    rbind(
      data.frame(
        node = V[k],
        which = "head",
        .classification_metrics(G1[, k], G2[, k])
      ),
      data.frame(
        node = V[k],
        which = "tail",
        .classification_metrics(G1[k, ], G2[k, ])
      )
    )
  }) |> do.call("rbind", args = _)
  if (!summarize) {
    return(res)
  }
  res |>
    dplyr::group_by(which) |>
    dplyr::summarize_at(
      c("precision", "recall", "f1", "fdr", "mcc", "acc"), mean,
      na.rm = TRUE
    ) |>
    dplyr::mutate_all(~ ifelse(is.nan(.x), NA, .x))
}

#' @exportS3Method prf1 pag
prf1.pag <- function(predicted, groundtruth, summarize = TRUE) {
  G1 <- predicted
  G2 <- groundtruth
  .check_graphs(G1, G2)
  lapply(1:3, \(mark) {
    ret <- prf1.default(1 * (G1 == mark), 1 * (G2 == mark), summarize)
    ret$mark <- mark
    ret
  }) |> do.call("rbind", args = _)
}

### Separation agreement
sep <- function(
    predicted, groundtruth, mode, max_size = NULL,
    precomputed_predicted = NULL, precomputed_groundtruth = NULL, ...) {
  G1 <- predicted
  G2 <- groundtruth
  .check_graphs(G1, G2)
  if (isTRUE(all.equal(unclass(G1), unclass(G2)))) {
    return(data.frame(precision = 1, recall = 1, f1 = 1, fdr = 0, mcc = 1, acc = 1))
  }
  if (!is.null(precomputed_predicted)) {
    t1 <- list(p.value = precomputed_predicted)
  } else {
    t1 <- .compute_oracle_tests(G1, max_size = max_size, mode = mode)
  }
  if (!is.null(precomputed_groundtruth)) {
    t2 <- list(p.value = precomputed_groundtruth)
  } else {
    t2 <- .compute_oracle_tests(G2, max_size = max_size, mode = mode)
  }
  .classification_metrics(1 - t1$p.value, 1 - t2$p.value)
}

### Helpers
.check_graphs <- function(G1, G2) {
  stopifnot(is.matrix(G1) & is.matrix(G2))
  stopifnot(nrow(G1) == nrow(G2))
  stopifnot(rownames(G1) == rownames(G2))
}

.classification_metrics <- function(pred_vec, true_vec) {
  tp <- sum(true_vec == 1 & pred_vec == 1)
  fp <- sum(true_vec == 0 & pred_vec == 1)
  tn <- sum(true_vec == 0 & pred_vec == 0)
  fn <- sum(true_vec == 1 & pred_vec == 0)

  precision <- ifelse(tp + fp == 0, NA, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, NA, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, NA, 2 * precision * recall / (precision + recall))
  fdr <- ifelse(tp + fp == 0, NA, fp / (tp + fp))
  nf <- exp(sum(0.5 * log(c(tn + fn, fp + tp, tn + fp, fn + tp))))
  mcc <- ifelse(nf == 0, NA, (tp * tn - fp * fn) / nf)
  acc <- mean(pred_vec == true_vec)
  data.frame(precision = precision, recall = recall, f1 = f1, fdr = fdr, mcc = mcc, acc = acc)
}
