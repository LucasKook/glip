create_dmg <- function(
    nNodes, prob = NULL, M2prob = NULL,
    nodeNames = NULL, mDG = FALSE, DG = FALSE,
    diag = TRUE) {
  if (!is.null(nodeNames)) {
    if (nNodes != length(nodeNames)) {
      stop("Number of node names is not equal to nNodes.")
    }
  }

  if (is.null(prob)) {
    prob <- stats::runif(1) / 2
  }
  if (is.null(M2prob)) {
    M2prob <- stats::runif(1) / 4
  }

  M1 <- matrix(stats::rbinom(nNodes^2, 1, prob),
    dimnames = list(
      letters[1:nNodes],
      letters[1:nNodes]
    ),
    nrow = nNodes, byrow = TRUE
  )

  M2 <- matrix(stats::rbinom(nNodes^2, 1, M2prob),
    dimnames = list(
      letters[1:nNodes],
      letters[1:nNodes]
    ),
    nrow = nNodes, ncol = nNodes, byrow = TRUE
  )
  M2[upper.tri(M2)] <- 0
  M2 <- M2 + t(M2)
  diag(M2) <- diag(M2) / 2

  if (mDG) {
    diag(M2)[rowSums(M2) > 0] <- 1
  }

  if (DG) {
    M2[] <- 0
  }

  if (is.null(nodeNames)) {
    nodeNames <- letters[1:nNodes]
  }

  if (diag) {
    diag(M1) <- diag(M2) <- 1
  } else {
    diag(M1) <- diag(M2) <- 0
  }

  rownames(M1) <- colnames(M1) <- rownames(M2) <- colnames(M2) <- nodeNames

  list(M1 = M1, M2 = M2)
}

ancestor_matrix <- function(M) {
  1 * (expm::expm(M) > 0)
}

createD0 <- function(G) {
  V <- rownames(G$M1)
  d <- length(V)
  nb <- sum(G$M2) / 2
  if (nb == 0) {
    return(G)
  }
  VU <- c(V, letters[d + seq_len(nb)])

  D0 <- create_dmg(length(VU),
    prob = 0, M2prob = 0,
    diag = FALSE, nodeNames = VU
  )

  D0$M1[V, V] <- G$M1
  bd <- which(G$M2 == 1, arr.ind = TRUE)
  bd <- bd[bd[, 1] < bd[, 2], , drop = FALSE]

  for (i in 1:nb) {
    j1 <- bd[i, 1]
    j2 <- bd[i, 2]
    D0$M1[d + i, j1] <- D0$M1[d + i, j2] <- 1
  }

  D0
}

compute_pag <- function(graph) {
  d <- NROW(graph$M1)
  orig <- rownames(graph$M1)
  sorted <- attr(igraph::topo_sort(
    igraph::graph_from_adjacency_matrix(graph$M1),
    mode = "out"
  ), "names")
  sgraph <- graph
  sgraph$M1 <- graph$M1[sorted, ][, sorted]
  sgraph$M2 <- graph$M2[sorted, ][, sorted]
  D <- createD0(sgraph)
  gD <- suppressWarnings(as(D$M1, "graphNEL"))
  lts <- setdiff(seq_len(NROW(D$M1)), seq_len(d))
  if (identical(lts, integer(0))) lts <- NULL
  out <- pcalg::dag2pag(
    list(g = gD, verbose = FALSE),
    pcalg::dsepTest,
    graph = gD,
    L = lts, alpha = 0.5
  )
  pag <- out@amat
  dimnames(pag) <- list(sorted, sorted)
  pag[orig, ][, orig]
}
