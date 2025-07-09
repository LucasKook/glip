### check if two CG are LWF-Markov equivalent (does not check that amat1 and
### amat2 do in fact correspond to chain graphs)
is_markov_equivalent_cg <- function(amat1, amat2){

  d <- nrow(amat1)
  nn <- colnames(amat1)

  if (!isTRUE(all.equal(colnames(amat1), colnames(amat2)))){
    stop("Column names do not match.")
  }

  for (i in nn){
    for (j in setdiff(nn, i)){
      for (k in seq(0, d-2)){
        for (C in combn(setdiff(nn, c(i,j)),k,simplify=FALSE)){
          if (csep(amat1,i,j,C) != csep(amat2,i,j,C)){
            return(FALSE)
          }
        }
      }
    }
  }

  TRUE
}

### generate random chain graph (Ma et al 2008 provides a different algorithm)
create_cg <- function(d, V = letters[1:d], prob = 0.5){

  ss <- sample(c(0,1), d^2, prob = c(1-prob, prob), replace = TRUE)
  amat <- matrix(ss, nrow = d, ncol = d)
  symamat <- 1*(amat + t(amat) > 1)
  com <- igraph::components(
    igraph::graph_from_adjacency_matrix(symamat)
  )$membership

  diag(amat) <- 0

  for (i in seq_len(d)){
    for (j in setdiff(seq_len(d), i)){
      if (com[i] > com[j]){
        amat[i,j] <- 0
      }
      if (com[i] == com[j]){
        if (sum(amat[i,j] + amat[j,i]) == 1){
          amat[i,j] <- amat[j,i] <- sample(c(0,1), 1)
        }
      }
    }
  }

  oo <- sample(seq_len(d))
  amat <- amat[,oo][oo,]
  rownames(amat) <- colnames(amat) <- V
  amat
}

### compute largest LWF-Markov equivalent chain graph in the Markov equivalence
### class of amat
compute_largest_cg <- function(amat){ 

  G <- amat
  V <- colnames(amat)
  amat <- matrix(nrow = nrow(amat), ncol = ncol(amat))
  dimnames(amat) <- list(V, V)
  amat[] <- G
  class(amat) <- c("lcg", class(amat))


  notLarge <- TRUE
  while(notLarge){

    gm <- compute_insub_metaarrow(amat)

    if (gm$large){
      return(amat)
      notLarge <- FALSE
    } else {
      A <- gm$A
      B <- gm$B
      amat[A,B][t(amat[B,A]) == 1] <- 1
      amat[B,A][t(amat[A,B]) == 1] <- 1
    }
  }
  
  amat
}

### auxiliary function to find non-empty meta-arrow (A \Rightarrow B) which is
### insubstantial, if it exists
compute_insub_metaarrow <- function(amat){ 
  symamat <- 1*(amat + t(amat) > 1)
  mem <- igraph::components(
    igraph::graph_from_adjacency_matrix(symamat)
  )$membership
  com <- unique(mem)

  for (i in com){
    A <- which(mem == i)
    for (j in setdiff(com, i)){
      B <- which(mem == j)
      if (max(amat[A,B]) > 0){ # nonempty meta-arrow A -> B
        aa <- which(seq_len(d) %in% A & rowSums(amat[,B,drop=FALSE])>0) # pa(B) \cap A
        mm <- amat[aa,aa,drop=FALSE]
        mm <- 1*(mm + t(mm) > 0)
        diag(mm) <- 1

        r1 <- sum(mm) == nrow(mm)^2

        paB <- setdiff(which(rowSums(amat[,B,drop=FALSE]) > 0), B)
        for (a in aa){
          paAlpha <- setdiff(which(rowSums(amat[,a,drop=FALSE]) > 0), a)
          r1 <- c(r1, length(setdiff(setdiff(paB, A), paAlpha)) == 0)
        }
        if (all(r1)){
          return(list(large = FALSE, A = A, B = B))
        }
      }
    }
  }
  list(large = TRUE)
}

### decide separation in LWF-CG
csep <- function(amat,i,j,C = c()){ 
  an <- ancestor_matrix(amat)
  anijC <- which(rowSums(an[,c(i,j,C)]) > 0)
  amat <- amat[anijC,anijC]
  d1 <- nrow(amat)
  mor <- moralize_cg(amat)
  if (length(C) > 0){
    mor[,C] <- 0
    mor[C,] <- 0
  }
  com <- igraph::components(
    igraph::graph_from_adjacency_matrix(mor)
  )$membership
  if (com[i] == com[j]){
    FALSE
  } else {
    TRUE
  }
}

### compute moral graph
moralize_cg <- function(amat){ 
  symamat <- 1*(amat + t(amat) > 1)
  com <- igraph::components(
    igraph::graph_from_adjacency_matrix(symamat)
  )$membership

  newamat <- amat
  for (i in unique(com)){
    rs <- rowSums(amat[,com == i,drop = FALSE])
    pa <- setdiff(which(rs > 0), which(com == i))
    if (length(pa) > 1){
      tmp <- matrix(1, nrow = length(pa), ncol = length(pa))
      diag(tmp) <- 0
      newamat[pa, pa] <- tmp
    }
  }
  1*(newamat + t(newamat) > 0)
}
