.compute_oracle_tests <- function(
    G, max_size = NULL, mode = "dag", verbose = FALSE,
    restrict_to = NULL, parallel = FALSE, ncores = NULL) {
  V <- .get_node_set(G)
  if (!is.null(restrict_to)) {
    V <- restrict_to
  }
  d <- length(V)
  if (is.null(max_size)) {
    max_size <- d - 2
  }
  sets <- .list_tests_graph(V, max_size = max_size)$sets
  if (verbose) {
    pb <- txtProgressBar(0, length(sets), style = 3, width = 60)
  }

  if (parallel) {
    nc <- min(parallel::detectCores() - 1, 15, ncores)
    plan(multisession, workers = nc)
    my_apply <- \(...) future_lapply(..., future.seed = TRUE)
  } else {
    my_apply <- lapply
  }

  res <- my_apply(seq_along(sets), \(iter) {
    if (verbose) {
      setTxtProgressBar(pb, iter)
    }
    x <- sets[[iter]]
    data.frame(
      X = x$X,
      Y = x$Y,
      Z = paste0(x$Z, collapse = ","),
      size = length(x$Z),
      p.value = .check_separation(x$X, x$Y, x$Z, G, mode)
    )
  }) |> do.call("rbind", args = _)

  if (parallel) {
    plan(sequential)
  }

  res
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

.get_opt <- function(mode) {
  switch(mode,
    "dag-dc" = dag_optim,
    "dg-dc" = dag_optim,
    "admg" = admg_lean_optim,
    "dmg" = admg_lean_optim,
    "admg-dc" = admg_optim,
    "dmg-dc" = admg_optim,
    "chain" = chain_lean_optim,
    "chain-dcon" = chain_optim,
    "dag-dcon" = dcon_optim,
    "dag" = dcon_lean_optim,
    "dg" = dcon_lean_optim
  )
}

.check_separation <- function(A, B, C, G, mode) {
  switch(mode,
    "dag" = .check_msep(A, B, C, G),
    "dg" = .check_msep(A, B, C, G),
    "admg" = .check_msep(A, B, C, G),
    "admg-dc" = .check_msep(A, B, C, G),
    "dmg-dc" = .check_msep(A, B, C, G),
    "dmg" = .check_msep(A, B, C, G),
    "chain" = .check_csep(A, B, C, G),
    "chain-dcon" = .check_csep(A, B, C, G),
    "dag-dcon" = .check_msep(A, B, C, G),
    "dag-dc" = .check_msep(A, B, C, G),
    "dg-dc" = .check_msep(A, B, C, G),
    "mag" = .check_dsepmag(A, B, C, G),
    "pdag" = .check_dseppdag(A, B, C, G)
  )
}

.check_equivalence <- function(G1, G2, max_size, mode = "dag") {
  V <- .get_node_set(G1)
  V2 <- .get_node_set(G2)
  stopifnot(isTRUE(all.equal(V, V2)))
  T1 <- .compute_oracle_tests(G1, max_size, mode)
  T2 <- .compute_oracle_tests(G2, max_size, mode)
  identical(
    which(T1$p.value == 1),
    which(T2$p.value == 1)
  )
}

.check_csep <- function(A, B, C, G) {
  1 * csep(G, A, B, C)
}

.check_msep <- function(A, B, C, G) {
  if (!is.list(G)) {
    tmp <- G
    tmp[] <- 0
    G <- list(M1 = G, M2 = tmp)
  }
  1 * is_m_separated(G, A, B, C)
}

.check_dsepmag <- function(A, B, C, G) {
  class(G) <- c("matrix", "array")
  V <- .get_node_set(G)
  MAG <- pcalg::pag2magAM(G, x = 1, max.chordal = nrow(G) + 1)
  1 * pcalg::dsepAM(match(A, V), match(B, V), match(C, V), MAG)
}

.check_dseppdag <- function(A, B, C, G) {
  class(G) <- c("matrix", "array")
  V <- .get_node_set(G)
  DAG <- pcalg::pdag2dag(as(G, "graphNEL"))$graph
  tmp <- as(DAG, "matrix")
  .check_msep(A, B, C, tmp)
}

.generate_random_graph <- function(
    d = 3, V = letters[1:d], mode, ...) {
  switch(mode,
    "dag" = .random_dmg(d, V, acyclic = TRUE, ...)$M1,
    "dg" = .random_dmg(d, V, acyclic = FALSE, ...)$M1,
    "admg" = .random_dmg(d, V, acyclic = TRUE, ...),
    "admg-dc" = .random_dmg(d, V, acyclic = TRUE, ...),
    "dmg-dc" = .random_dmg(d, V, acyclic = FALSE, ...),
    "dmg" = .random_dmg(d, V, acyclic = FALSE, ...),
    "chain" = .random_cg(d, V, ...),
    "chain-dcon" = .random_cg(d, V, ...),
    "dag-dcon" = .random_dmg(d, V, acyclic = TRUE, ...)$M1,
    "dag-dc" = .random_dmg(d, V, acyclic = TRUE, ...)$M1,
    "dg-dc" = .random_dmg(d, V, acyclic = TRUE, ...)$M1
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

.compute_graphical_representation <- function(G, max_size, mode) {
  d <- .get_size(G)
  if (max_size < d - 2) {
    return(NULL)
  }
  switch(mode,
    "dag" = .dag2ess(G),
    "dg" = NULL,
    "admg" = .admg2pag(G),
    "admg-dc" = .admg2pag(G),
    "dmg-dc" = NULL,
    "dmg" = NULL,
    "chain" = .cg2lcg(G),
    "chain-dcon" = .cg2lcg(G),
    "dag-dcon" = .dag2ess(G),
    "dag-dc" = .dag2ess(G),
    "dg-dc" = NULL
  )
}

.cg2lcg <- function(G) {
  compute_largest_cg(G)
}

.dag2ess <- function(G) {
  ret <- 1 * pcalg::dag2essgraph(G)
  class(ret) <- c("ess", class(ret))
  ret
}

.admg2pag <- function(G) {
  if (!is.list(G)) {
    G <- .to_admg(G)
  }
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

.multigrep <- function(patterns, x, ...) {
  unlist(lapply(patterns, \(p) {
    grep(p, x, ...)
  }))
}

.to_admg <- function(G) {
  tmp <- G
  tmp[] <- 0
  list(M1 = G, M2 = tmp)
}

.pag_to_admg <- function(G) {
  if ("pag" %in% class(G)) {
    class(G) <- c("matrix", "array")
  }
  mag <- pcalg::pag2magAM(G, x = 1, max.chordal = nrow(G) + 1)
  V <- rownames(G)
  d <- length(V)
  tmp <- matrix(0, nrow = d, ncol = d, dimnames = list(V, V))
  admg <- list(M1 = tmp, M2 = tmp)
  for (i in V) {
    for (j in setdiff(V, i)) {
      admg$M1[i, j] <- 1 * (mag[i, j] == 2 & mag[j, i] == 3)
      admg$M2[i, j] <- 1 * (mag[i, j] == 2 & mag[j, i] == 2)
    }
  }
  admg
}

.ess_to_dag <- function(G) {
  if ("ess" %in% class(G) | "pag" %in% class(G)) {
    class(G) <- c("matrix", "array")
  }
  as(pcalg::pdag2dag(as(G, "graphNEL"))$graph, "matrix")
}
