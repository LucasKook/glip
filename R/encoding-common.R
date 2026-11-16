### Steps shared by all encodings (min-length, classic and cut encodings):
### background knowledge, conditioning sets, mapping of the tests to the
### decision variables with weights, and solving

### Validate the background knowledge of a class ("dag" for DAGs and DGs,
### "admg" for ADMGs and DMGs, "chain" for chain graphs), raise an error for
### knowledge (passed via `...` of the optimizer `fun`) that the class does
### not support, and make the warm start consistent with the knowledge.
### Returns the normalized knowledge `bk`, for chain graphs also `bkc` (see
### .bk_chain()), and the warm start
.bk_prepare <- function(class, fun, V, d, warmstart, dots, ...) {
  bk <- .bk_validate(V, ...)
  bkc <- NULL
  switch(class,
    dag = {
      .bk_unsupported(dots, .bk_args_dag, fun,
        "DAGs and DGs have only directed edges")
      warmstart <- .bk_warmstart(warmstart, bk$forbidden, bk$required)
    },
    admg = {
      .bk_unsupported(dots, .bk_args_admg, fun,
        "ADMGs and DMGs have no undirected edges")
      if (!is.null(warmstart)) {
        warmstart$M1 <- .bk_warmstart(warmstart$M1, bk$forbidden, bk$required)
        warmstart$M2 <- .bk_warmstart(warmstart$M2, bk$forbidden_bi, bk$required_bi)
      }
    },
    chain = {
      .bk_unsupported(dots, .bk_args_chain, fun,
        "chain graphs have no bidirected edges")
      bkc <- .bk_chain(bk$forbidden, bk$required, bk$forbidden_undirected,
        bk$required_undirected, d)
      warmstart <- .bk_warmstart(warmstart, bkc$zero, bkc$one)
    }
  )
  list(bk = bk, bkc = bkc, warmstart = warmstart)
}

### List of all conditioning sets of size at most max_size
.cond_sets <- function(d, max_size) {
  CC <- unlist(sapply(0:max_size, \(x)
    utils::combn(d, x, simplify = FALSE)),
    recursive = FALSE
  )
  if (max_size == 0) {
    CC <- list(integer(0))
  }
  CC
}

### Map the tests to the n_z separation variables indexed by `zijc` (one per
### unordered pair i, j and conditioning set CC[[C]]): p-values, weights and
### decisions p = trafo(p-value). Variables without a test keep weight zero;
### for duplicated tests (e.g., X, Y and Y, X) the last one counts. Returns
### the tests merged with the variable index (`merged`, with `weight`), the
### weights `w` and the decisions `p` of all variables
.test_decisions <- function(tests, zijc, CC, V, n_z, weight_type, trafo) {
  # Lookup table to place p-values in the correct positions
  # b/c self-edges and x -> y s.t. y in C are not tested
  l1 <- zijc
  l1$i <- V[zijc$i]
  l1$j <- V[zijc$j]
  l1$C <- sapply(CC[zijc$C], \(x) paste0(V[x], collapse = ","))
  l2 <- l1
  colnames(l1) <- c("X", "Y", "Z", "idx")
  colnames(l2) <- c("Y", "X", "Z", "idx")
  lookup <- rbind(l1, l2)

  # Merge p-values with lookup table for optimization
  # Positions that are not tested contain an edge by default
  # i -> j | C s.t. j not in C ignored => Set w to zero for those
  # i -> j | C s.t. i in C ignored => Set w to zero for those
  # i -> i | C ignored => Set w to zero for those
  w <- rep(0, n_z)
  s <- rep(0, n_z)
  praw <- rep(0, n_z)
  merged <- dplyr::left_join(tests, lookup, by = c("X", "Y", "Z"))
  praw[merged$idx] <- merged$p.value
  w[merged$idx] <- 1
  s[merged$idx] <- merged$size

  ### Weights
  w[w == 1] <- switch(weight_type,
    "const" = 1,
    "inv" = 1 / pmax(praw[w == 1], 0.001),
    "log" = -log2(pmax(praw[w == 1], 2 * .Machine$double.eps)) + 0.1,
    "size" = 1 / (1 + s[w == 1])
  )

  merged$weight <- w[merged$idx]
  list(merged = merged, w = w, p = trafo(praw))
}

### Solve a model with Gurobi
.gurobi_solve <- function(model, gurobi_args) {
  if (requireNamespace("gurobi")) {
    gurobi::gurobi(model, gurobi_args)
  }
}
