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

invert_latent_projection <- function(G) {
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
  D <- invert_latent_projection(sgraph)
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

moralize <- function(G, noWarning = FALSE) {
  M1 <- G$M1
  M2 <- G$M2

  if (is.null(rownames(M1))) {
    stop("Need named nodes.")
  }

  if (!noWarning && sum(abs(M2)) != 0) {
    warning("Should not be used on (non-LIG) MIGs")
  }

  # Compute (asymmetric) adjacency matrix
  M0 <- (M1 + M2) > 0

  # Compute (symmetric) adjacency matrix
  M <- M0 | t(M0)

  # Identify colliders and moralize the graph efficiently
  parent_pairs <- which(M0, arr.ind = TRUE)

  if (nrow(parent_pairs) > 0) {
    alpha_nodes <- unique(parent_pairs[, 2])

    for (alpha in alpha_nodes) {
      parents <- which(M0[, alpha])
      if (length(parents) > 1) {
        M[parents, parents] <- TRUE
      }
    }
  }

  storage.mode(M) <- "integer"
  M
}

augment <- function(G) {
  M1 <- G$M1
  M2 <- G$M2
  if (is.null(rownames(M1))) {
    stop("Need named nodes.")
  }
  V <- rownames(M1)

  # start from the moralGraph
  moralG <- moralize(G, noWarning = TRUE)

  # then add the collider connected nodes (collision path of length > 2)
  for (i in 1:length(V)) {
    for (j in i:length(V)) {
      alpha <- V[i]
      beta <- V[j]
      if (moralG[alpha, beta] > 0) next
      M1tmp <- M1
      # remove the other directed arrows, can never be in a collider path of
      # length > 2 (length = 2 is already in the moral graph)
      M1tmp[-which(V %in% c(alpha, beta)), ] <- 0
      M0tmp <- 1 * (M1tmp + t(M1tmp) + M2 > 0)
      Mconn <- ancestor_matrix(M0tmp)
      if (Mconn[alpha, beta] == 1) {
        moralG[alpha, beta] <- moralG[beta, alpha] <- 1
      }
    }
  }

  moralG
}

is_undirected_separated <- function(M, A, B, C) {
  if (length(C) > 0) {
    tmp <- which(rownames(M) %in% C)
    sepM <- M[-tmp, -tmp, drop = FALSE]
  } else {
    sepM <- M
  }
  !max(ancestor_matrix(sepM)[A, B])
}

is_m_separated <- function(G, A, B, C) {
  M1 <- G$M1
  M2 <- G$M2

  stopifnot(
    length(intersect(A, C)) == 0 &
      length(intersect(B, C)) == 0 &
      length(intersect(A, B)) == 0
  )

  # find G(B)_An(A\cup B \cup C)
  # determine ancestry of each node ([i,j] == 1 indicates that i is an ancestor
  # of j)
  An <- ancestor_matrix(G$M1)
  # find the relevant set (i.e., ancestors of A\cup B \cup C)
  isAn <- rowSums(An[, c(A, B, C), drop = FALSE]) > 0
  M1An <- G$M1[isAn, isAn, drop = FALSE]
  M2An <- G$M2[isAn, isAn, drop = FALSE]

  # find the augmented graph
  g <- augment(list(M1 = M1An, M2 = M2An))
  is_undirected_separated(g, A, B, C)
}
