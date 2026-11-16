### Skip tests of the cut encoding if Python (option glip.python, GLIP_PYTHON
### or python3) lacks networkx (oracle tests) or gurobipy with a usable
### licence (solver tests)
skip_if_no_glipcut_python <- function() {
  if (is.null(.glipcut_ok$py)) {
    py <- .glipcut_python()
    .glipcut_ok$py <- identical(as.integer(suppressWarnings(tryCatch(
      system2(py, c("-c", shQuote("import networkx")), stdout = FALSE, stderr = FALSE),
      error = function(e) 1L))), 0L)
  }
  if (!.glipcut_ok$py) testthat::skip("Python with networkx is not available")
  invisible(TRUE)
}
skip_if_no_glipcut <- function() {
  skip_if_no_glipcut_python()
  if (is.null(.glipcut_ok$solver)) {
    .glipcut_ok$solver <- .has_glipcut() && tryCatch({
      tests <- data.frame(X = "a", Y = "b", Z = "", p.value = 0.01, size = 0)
      f <- cuts_optim(tests, d = 2, max_size = 0, V = c("a", "b"), class = "dag",
        gurobi_args = list(Threads = 1))
      identical(f$optim$status, "OPTIMAL")
    }, error = function(e) FALSE)
  }
  if (!.glipcut_ok$solver) {
    testthat::skip("Python with gurobipy and a Gurobi licence is not available")
  }
  invisible(TRUE)
}
.glipcut_ok <- new.env()

### Connection indicators of queries in a graph computed by the Python module
### (no solver), for comparison with the R oracles
.glipcut_connected <- function(class, G, tests, V) {
  C <- lapply(strsplit(as.character(tests$Z), ",", fixed = TRUE),
    \(z) match(z[nzchar(z)], V))
  graph <- switch(class,
    dag = , dg = list(M1 = unname(G[, ])),
    chain = list(G = unname(G[, ])),
    admg = , dmg = list(M1 = unname(G$M1[, ]), M2 = unname(G$M2[, ]))
  )
  case <- list(task = "connected", class = class, d = length(V), graph = graph,
    queries = lapply(seq_len(nrow(tests)), \(k) list(
      i = match(tests$X[k], V) - 1L, j = match(tests$Y[k], V) - 1L,
      C = as.list(C[[k]] - 1L))))
  res <- .glipcut_run(case)
  list(valid = isTRUE(res$valid), connected = unlist(res$connected))
}

### Full test bank (unordered pairs, all conditioning sets up to max_size)
### with p-values from the oracle of a graph, a share of decisions flipped
.cuts_bank <- function(G, V, max_size, mode, flip = 0.2) {
  tests <- .compute_oracle_tests(G, max_size = max_size, mode = mode)
  p <- 1 - tests$p.value
  if (flip > 0) {
    k <- sample.int(length(p), max(1, round(flip * length(p))))
    p[k] <- 1 - p[k]
  }
  tests$p.value <- ifelse(p == 1, stats::runif(length(p), 0.001, 0.04),
    stats::runif(length(p), 0.1, 0.99))
  tests
}
