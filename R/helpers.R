.compute_oracle_tests <- function(G, max_size = NULL, mode = "dag") {
  V <- .get_node_set(G)
  if (is.null(max_size)) {
    max_size <- length(V) - 2
  }
  sets <- .list_tests_graph(V, max_size = max_size)$sets
  tests <- lapply(sets, \(x) {
    data.frame(
      X = x$X,
      Y = x$Y,
      Z = paste0(x$Z, collapse = ","),
      p.value = .check_separation(x$X, x$Y, x$Z, G, mode)
    )
  }) |> do.call("rbind", args = _)
}

.get_node_set <- function(G) {
  if (is.list(G)) {
    return(rownames(G$M1))
  }
  rownames(G)
}

.get_size <- function(G) {
  if (is.list(G)) {
    return(nrow(G$M1))
  }
  nrow(G)
}

.get_opt <- function(mode = c("dag", "dg", "admg", "dmg", "chain", "dagdcon")) {
  mode <- match.arg(mode)
  switch(mode,
    "dag" = dag_optim,
    "dg" = dag_optim,
    "admg" = admg_optim,
    "dmg" = admg_optim,
    "chain" = chain_optim,
    "dagdcon" = dcon_optim
  )
}

.check_separation <- function(A, B, C, G, mode = c("dag", "dg", "admg", "dmg", "chain", "dagdcon")) {
  mode <- match.arg(mode)
  switch(mode,
    "dag" = .check_dsep(A, B, C, G),
    "dg" = .check_dsep(A, B, C, G), # TODO: Update
    "admg" = .check_dsep_dmg(A, B, C, G),
    "dmg" = .check_dsep_dmg(A, B, C, G), # TODO: Update
    "chain" = .check_csep(A, B, C, G),
    "dagdcon" = .check_dsep(A, B, C, G)
  )
}

.check_csep <- function(A, B, C, G) {
  1 * csep(G, A, B, C)
}

.check_dsep_dmg <- function(A, B, C, G) {
  D <- createD0(G)
  gD <- suppressWarnings(as(D$M1, "graphNEL"))
  1 * pcalg::dsep(A, B, C, g = gD)
}

.check_dsep <- function(A, B, C, G) {
  GG <- suppressWarnings(as(G, "graphNEL"))
  1 * pcalg::dsep(A, B, C, g = GG)
}

.generate_random_graph <- function(
    d = 3, V = letters[1:d], mode = c("dag", "dg", "admg", "dmg", "chain", "dagdcon"), ...) {
  mode <- match.arg(mode)
  switch(mode,
    "dag" = .random_dag(d, V, ...),
    "dg" = .random_dag(d, V, ...), # TODO: Update
    "admg" = .random_dmg(d, V, acyclic = TRUE, ...),
    "dmg" = .random_dmg(d, V, acyclic = FALSE, ...),
    "chain" = .random_cg(d, V, ...),
    "dagdcon" = .random_dag(d, V, ...)
  )
}

.random_cg <- function(d, V = letters[1:d], ...) {
  create_cg(d, V, ...)
}

.random_dmg <- function(d, V = letters[1:d], acyclic = FALSE, ...) {
  G <- create_dmg(d, nodeNames = V, diag = FALSE, ...)
  if (acyclic) {
    G$M1[upper.tri(G$M1)] <- 0
  }
  G
}

.random_dag <- function(d, V = letters[1:d], ...) {
  G <- 1 * (as(pcalg::randomDAG(d, ...), "matrix") != 0)
  dimnames(G) <- list(V, V)
  G
}

.compute_graphical_representation <- function(G, mode = c("dag", "dg", "admg", "dmg", "chain", "dagdcon")) {
  mode <- match.arg(mode)
  switch(mode,
    "dag" = .dag2ess(G),
    "dg" = .dag2ess(G), # TODO: Update
    "admg" = .admg2pag(G),
    "dmg" = .admg2pag(G), # TODO: Update
    "chain" = compute_largest_cg(G),
    "dagdcon" = .dag2ess(G)
  )
}

.dag2ess <- function(G) {
  ret <- 1 * pcalg::dag2essgraph(G)
  class(ret) <- c("ess", class(ret))
  ret
}

.admg2pag <- function(G) {
  ret <- compute_pag(G)
  class(ret) <- c("pag", class(ret))
  ret
}

.add <- function(x) paste0(x, collapse = "+")

.rm_int <- function(x) {
  if (all(x[, 1] == 1)) {
    return(x[, -1L, drop = FALSE])
  }
  x
}

### Helper for handling args
darg <- function(x, d) {
  if (is.na(x)) {
    return(d)
  } else {
    x
  }
}

.multigrep <- function(patterns, x) {
  unlist(lapply(patterns, \(p) {
    grep(p, x)
  }))
}
