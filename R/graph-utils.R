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

# G: DMG over V, O: observed subset of V
marginalize <- function(G, O) {
  V <- rownames(G$M1)
  M <- setdiff(V, O)

  Gold <- G
  Gold$M1[] <- Gold$M2[] <- 0

  while (!isTRUE(all.equal(G, Gold))) {
    Gold <- G

    for (m in M) {
      for (alpha in setdiff(V, m)) {
        for (beta in setdiff(V, m)) {
          if (G$M1[alpha, m] == 1 && G$M1[m, beta] == 1) {
            G$M1[alpha, beta] <- 1
          }
          if (G$M2[alpha, m] == 1 && G$M1[m, beta] == 1) {
            G$M2[alpha, beta] <- G$M2[beta, alpha] <- 1
          }
          if (G$M1[m, alpha] == 1 && G$M1[m, beta] == 1) {
            G$M2[alpha, beta] <- G$M2[beta, alpha] <- 1
          }
        }
      }
    }
  }

  list(
    M1 = (G$M1)[O, O, drop = FALSE],
    M2 = (G$M2)[O, O, drop = FALSE]
  )
}

is_max_dmg <- function(G, verbose = FALSE, max_size = NROW(G$M1) - 2) {
  M1 <- G$M1
  M2 <- G$M2
  nodes <- rownames(M1)

  for (alpha in nodes) {
    for (beta in nodes) {
      if (M1[alpha, beta] == 0) {
        Gtmp <- G
        Gtmp$M1[alpha, beta] <- 1
        aM <- is_markov_equivalent(G, Gtmp, max_size = max_size)
        if (aM) {
          if (verbose) {
            return(Gtmp)
          } else {
            return(FALSE)
          }
        }
      }
      if (M2[alpha, beta] == 0) {
        Gtmp <- G
        Gtmp$M2[alpha, beta] <- Gtmp$M2[beta, alpha] <- 1
        aM <- is_markov_equivalent(G, Gtmp, max_size = max_size)
        if (aM) {
          if (verbose) {
            return(Gtmp)
          } else {
            return(FALSE)
          }
        }
      }
    }
  }

  TRUE
}

is_markov_equivalent <- function(G1, G2, verbose = FALSE, O = NULL,
                                 max_size = NROW(G1$M1) - 2) {
  max_size <- max_size + 1
  if (!isTRUE(all.equal(names(G1), c("M1", "M2")))) {
    stop("G1 is not a DMG object.")
  }
  if (!isTRUE(all.equal(names(G2), c("M1", "M2")))) {
    stop("G2 is not a DMG object.")
  }
  if (!isTRUE(all.equal(rownames(G1$M1), rownames(G2$M1)))) {
    warning("Node names provided differ between the two graphs.")
  }

  nodes <- union(rownames(G1$M1), rownames(G2$M1))
  mu1 <- mu2 <- c()
  res.l <- list()

  if (is.null(O)) {
    checkNodes <- nodes
  } else {
    checkNodes <- O
  }


  for (alpha in checkNodes) {
    for (beta in checkNodes) {
      C <- setdiff(checkNodes, alpha)
      tmp1 <- is_mu_separated(G1, alpha, beta, c())
      tmp2 <- is_mu_separated(G2, alpha, beta, c())
      mu1 <- c(mu1, tmp1)
      mu2 <- c(mu2, tmp2)
      if (tmp1 != tmp2) {
        res.l[[length(res.l) + 1]] <- paste(alpha, "-", beta, "-", paste("",
          collapse = ""
        ), "-", tmp1, tmp2, collapse = " ")
      }
      if (length(C) > 0) {
        # for (k in 1:length(C)) {
        for (k in 1:max_size) {
          for (Z in utils::combn(C, k, simplify = FALSE)) {
            tmp1 <- is_mu_separated(G1, alpha, beta, Z)
            tmp2 <- is_mu_separated(G2, alpha, beta, Z)
            mu1 <- c(mu1, tmp1)
            mu2 <- c(mu2, tmp2)
            if (tmp1 != tmp2) {
              res.l[[length(res.l) + 1]] <- paste(alpha, "-", beta, "-",
                paste(Z, collapse = ""), "-", tmp1, tmp2,
                collapse = " "
              )
            }
          }
        }
      }
    }
  }

  if (verbose) {
    return(res.l)
  } else {
    isTRUE(all.equal(mu1, mu2))
  }
}

historize <- function(G, B) {
  nB <- length(B)
  if (nB == 0) stop("Empty B-set provided.")

  n <- nrow(G$M1)
  M1 <- M2 <- matrix(0, nrow = n + nB, ncol = n + nB)
  M1[1:n, 1:n] <- G$M1
  M2[1:n, 1:n] <- G$M2
  M1[1:n, (n + 1):(n + nB)] <- G$M1[, B]
  M2[1:n, (n + 1):(n + nB)] <- G$M2[, B]
  M2[(n + 1):(n + nB), 1:n] <- G$M2[B, ]
  nms <- c(rownames(G$M1), paste(B, "p", sep = ""))
  rownames(M1) <- colnames(M1) <- rownames(M2) <- colnames(M2) <- nms
  list(M1 = M1, M2 = M2)
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

is_mu_separated <- function(G, A, B, C, m_separation = FALSE) {
  M1 <- G$M1
  M2 <- G$M2
  # create extra nodes to be included in the B-history graph
  Bp <- paste(B, "p", sep = "")

  A <- setdiff(A, C)

  if (length(A) == 0) {
    return(TRUE)
  }

  if (!m_separation) {
    bH <- historize(G, B) # construct B-history graph
  } else { # this is for m-separation only
    bH <- G
    Bp <- B

    B <- setdiff(B, C)

    if (length(B) == 0) {
      return(TRUE)
    }
  }

  # find G(B)_An(A\cup Bp \cup C)
  # determine ancestry of each node ([i,j] == 1 indicates that i is an ancestor of j)
  An <- ancestor_matrix(bH$M1)
  # find the relevant set (i.e., ancestors of A\cup Bp \cup C)
  isAn <- rowSums(An[, c(A, Bp, C), drop = FALSE]) > 0
  M1An <- bH$M1[isAn, isAn, drop = FALSE]
  M2An <- bH$M2[isAn, isAn, drop = FALSE]

  # find the augmented graph
  g <- augment(list(M1 = M1An, M2 = M2An))

  is_undirected_separated(g, A, Bp, C)
}

ancestor_matrix <- function(M) {
  1 * (expm::expm(M) > 0)
}

# binary decision: \(x) as.numeric(x > 0.05)
objective <- function(
    graph, data, trafo = identity) {
  contribs <- apply(data, 1, \(x) {
    from <- x["from"]
    to <- x["to"]
    w <- as.numeric(x["weight"])
    given <- x["given"]
    given <- unlist(strsplit(given, ","))
    pv <- as.numeric(x["p.value"])
    z <- as.numeric(is_mu_separated(graph, from, to, given))
    w * abs(z - trafo(pv))
  })
  sum(contribs)
}

generate_all_DMGs <- function(n = 4, V) {
  # Number of possible directed edges
  num_directed_edges <- n * (n - 1)
  # Number of possible bidirected edges
  num_bidirected_edges <- (n * (n - 1)) / 2

  # Total number of graphs
  total_graphs <- 2^(num_directed_edges + num_bidirected_edges)

  # List to store all graphs
  all_graphs <- vector("list", total_graphs)

  # Generate all combinations
  for (i in 0:(total_graphs - 1)) {
    binary_string <- intToBits(i)[1:(num_directed_edges + num_bidirected_edges)]

    # Extract directed edges
    M1 <- matrix(0, n, n)
    directed_indices <- 1
    for (row in 1:n) {
      for (col in 1:n) {
        if (row != col) { # Avoid self-loops
          M1[row, col] <- as.integer(binary_string[directed_indices])
          directed_indices <- directed_indices + 1
        }
      }
    }
    diag(M1) <- 1
    colnames(M1) <- V
    rownames(M1) <- V

    # Extract bidirected edges
    M2 <- matrix(0, n, n)
    bidirected_indices <- num_directed_edges + 1
    for (row in 1:(n - 1)) {
      for (col in (row + 1):n) {
        value <- as.integer(binary_string[bidirected_indices])
        M2[row, col] <- value
        M2[col, row] <- value # Symmetric
        bidirected_indices <- bidirected_indices + 1
      }
    }
    diag(M2) <- 1
    colnames(M2) <- V
    rownames(M2) <- V

    # Store the graph as a list
    all_graphs[[i + 1]] <- list(M1 = M1, M2 = M2)
  }

  return(all_graphs)
}

find_max_dmg <- function(G, max_size = NROW(G$M1) - 2) {
  nodes <- rownames(G$M1)

  for (alpha in nodes) {
    for (beta in nodes) {
      if (G$M1[alpha, beta] == 0) {
        Gtmp <- G
        Gtmp$M1[alpha, beta] <- 1
        aM <- is_markov_equivalent(G, Gtmp, max_size = max_size)
        if (aM) {
          G$M1[alpha, beta] <- 1
        }
      }
      if (G$M2[alpha, beta] == 0) {
        Gtmp <- G
        Gtmp$M2[alpha, beta] <- Gtmp$M2[beta, alpha] <- 1
        aM <- is_markov_equivalent(G, Gtmp, max_size = max_size)
        if (aM) {
          G$M2[alpha, beta] <- G$M2[beta, alpha] <- 1
        }
      }
    }
  }

  G
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
