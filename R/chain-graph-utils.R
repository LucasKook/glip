### check if two CG are LWF-Markov equivalent (does not check that G1 and
### G2 do in fact correspond to chain graphs)
is_markov_equivalent_cg <- function(G1, G2) {
  d <- nrow(G1)
  nn <- colnames(G1)

  if (!isTRUE(all.equal(colnames(G1), colnames(G2)))) {
    stop("Column names do not match.")
  }

  for (i in nn) {
    for (j in setdiff(nn, i)) {
      for (k in seq(0, d - 2)) {
        for (C in combn(setdiff(nn, c(i, j)), k, simplify = FALSE)) {
          if (csep(G1, i, j, C) != csep(G2, i, j, C)) {
            return(FALSE)
          }
        }
      }
    }
  }

  TRUE
}

### generate random chain graph (Ma et al 2008 provides a different algorithm)
create_cg <- function(d, V = letters[1:d], prob = 0.5) {
  ss <- sample(c(0, 1), d^2, prob = c(1 - prob, prob), replace = TRUE)
  G <- matrix(ss, nrow = d, ncol = d)
  symG <- 1 * (G + t(G) > 1)
  com <- igraph::components(
    igraph::graph_from_adjacency_matrix(symG)
  )$membership

  diag(G) <- 0

  for (i in seq_len(d)) {
    for (j in setdiff(seq_len(d), i)) {
      if (com[i] > com[j]) {
        G[i, j] <- 0
      }
      if (com[i] == com[j]) {
        if (sum(G[i, j] + G[j, i]) == 1) {
          G[i, j] <- G[j, i] <- sample(c(0, 1), 1)
        }
      }
    }
  }

  oo <- sample(seq_len(d))
  G <- G[, oo][oo, ]
  rownames(G) <- colnames(G) <- V
  G
}

### compute largest LWF-Markov equivalent chain graph in the Markov equivalence
### class of G
compute_largest_cg <- function(G) {
  Gtmp <- G
  V <- colnames(G)
  G <- matrix(nrow = nrow(G), ncol = ncol(G))
  dimnames(G) <- list(V, V)
  G[] <- Gtmp
  class(G) <- c("lcg", class(G))


  notLarge <- TRUE
  while (notLarge) {
    gm <- compute_insub_metaarrow(G)

    if (gm$large) {
      return(G)
      notLarge <- FALSE
    } else {
      A <- gm$A
      B <- gm$B
      G[A, B][t(G[B, A]) == 1] <- 1
      G[B, A][t(G[A, B]) == 1] <- 1
    }
  }

  G
}

### auxiliary function to find non-empty meta-arrow (A \Rightarrow B) which is
### insubstantial, if it exists
compute_insub_metaarrow <- function(G) {
  symG <- 1 * (G + t(G) > 1)
  mem <- igraph::components(
    igraph::graph_from_adjacency_matrix(symG)
  )$membership
  com <- unique(mem)

  for (i in com) {
    A <- which(mem == i)
    for (j in setdiff(com, i)) {
      B <- which(mem == j)
      if (max(G[A, B]) > 0) { # nonempty meta-arrow A -> B
        aa <- which(seq_len(d) %in% A & rowSums(G[, B, drop = FALSE]) > 0) # pa(B) \cap A
        mm <- G[aa, aa, drop = FALSE]
        mm <- 1 * (mm + t(mm) > 0)
        diag(mm) <- 1

        r1 <- sum(mm) == nrow(mm)^2

        paB <- setdiff(which(rowSums(G[, B, drop = FALSE]) > 0), B)
        for (a in aa) {
          paAlpha <- setdiff(which(rowSums(G[, a, drop = FALSE]) > 0), a)
          r1 <- c(r1, length(setdiff(setdiff(paB, A), paAlpha)) == 0)
        }
        if (all(r1)) {
          return(list(large = FALSE, A = A, B = B))
        }
      }
    }
  }
  list(large = TRUE)
}

### decide separation in LWF-CG
csep <- function(G, i, j, C = c()) {
  an <- ancestor_matrix(G)
  anijC <- which(rowSums(an[, c(i, j, C)]) > 0)
  G <- G[anijC, anijC]
  d1 <- nrow(G)
  mor <- moralize_cg(G)
  if (length(C) > 0) {
    mor[, C] <- 0
    mor[C, ] <- 0
  }
  com <- igraph::components(
    igraph::graph_from_adjacency_matrix(mor)
  )$membership
  if (com[i] == com[j]) {
    FALSE
  } else {
    TRUE
  }
}

### compute moral graph
moralize_cg <- function(G) {
  symG <- 1 * (G + t(G) > 1)
  com <- igraph::components(
    igraph::graph_from_adjacency_matrix(symG)
  )$membership

  newG <- G
  for (i in unique(com)) {
    rs <- rowSums(G[, com == i, drop = FALSE])
    pa <- setdiff(which(rs > 0), which(com == i))
    if (length(pa) > 1) {
      tmp <- matrix(1, nrow = length(pa), ncol = length(pa))
      diag(tmp) <- 0
      newG[pa, pa] <- tmp
    }
  }
  1 * (newG + t(newG) > 0)
}

### Generate random chain graphs, Ma et al algorithm
### n:number of variables
### d:degree of nodes
create_cg_ma <- function(n, d) {
  order <- sample(1:n, n, replace = FALSE)

  amat <- matrix(0, n, n)
  prob <- d / (n - 1)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      amat[i, j] <- Rlab::rbern(1, prob = prob)
    }
  }
  amat <- amat + t(amat)

  k <- sample(1:n, 1)
  if (k != 1) {
    chain_cut <- cut(1:n, k, 1:k)
  } else {
    chain_cut <- 1:n
  }
  for (i in 1:n) {
    for (j in 1:n) {
      if (as.integer(chain_cut[which(order == i)]) > as.integer(chain_cut[which(order == j)])) {
        amat[i, j] <- 0
      }
    }
  }
  rownames(amat) <- colnames(amat) <- letters[1:n]
  amat
}
