#' Optimal Learning of ADMGs via Mixed Integer Programming with Cubic Encoding
#'
#' Solves a global mixed integer optimization problem to learn an Acyclic
#' Directed Mixed Graph (ADMG) or DMG from supplied conditional independence
#' test results, using the `gurobi` solver. Supports warm-start, different
#' weight types, and optional caching for efficiency. This cubic variant
#' focuses on a redundant but faster encoding.
#'
#' @param tests A data frame or list of conditional independence test results.
#' @param d Integer. Number of variables/nodes.
#' @param max_size Integer. Maximum size of conditioning set.
#' @param V Character vector of variable names, of length \code{d}.
#' @param trafo Transformation function applied to test p-values (default:
#'     \code{function(x) as.numeric(x <= 0.05)}).
#' @param weight_type Character. Type of weighting ("const", "inv", "log", "size").
#' @param warmstart Optional. Matrix, list, or NULL. Used to warm-start optimization.
#' @param edgehints Optional. Edge hints to speed up or constrain the solution.
#' @param gurobi_args List of arguments passed to the `gurobi` solver.
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#' @param cache Logical. If \code{TRUE}, cache constraint matrices and results.
#' @param cache_dir Character. Directory where cache files are saved.
#' @param mode Character. Learning mode: "admg" (default) or "dmg".
#' @param ... Additional arguments passed to internals.
#'
#' @return An object of class \code{graphopt} with components:
#'   \describe{
#'     \item{graph}{Adjacency or edge matrices of the learned graph.}
#'     \item{tests}{Test results with optimization variables.}
#'     \item{optim}{Solver output.}
#'   }
#'
#' @details
#' Requires the `gurobi` package. Will warn if unavailable.
#'
#' @export
admg_cubic_optim <- function(
    tests, d = 3, max_size = d - 2, V = letters[1:d], trafo = \(x) as.numeric(x <= 0.05),
    weight_type = c("const", "inv", "log", "size"), warmstart = NULL, edgehints = NULL,
    gurobi_args = list(), verbose = FALSE, cache = TRUE, cache_dir = "./.cache-admg-cubic-redundant",
    mode = c("admg", "dmg"),
    ...) {

  if (!requireNamespace("gurobi")) {
    warning("Solver `gurobi` not available.")
  }

  if (!is.null(warmstart) & !is.list(warmstart) & is.matrix(warmstart)) {
    warmstart <- .to_admg(warmstart)
  }

  mode <- match.arg(mode)
  weight_type <- match.arg(weight_type)

  fn <- file.path(
    cache_dir,
    paste0(c("lhs", "rhs", "P1", "N1", "O1", "R1", "G1"), "dim", d, "max",
      max_size, mode, c(".mtx", ".rds", ".rds", ".rds", ".rds", ".rds", ".rds"))
  )

  ### COMPUTE DIMENSIONS

  # List of all conditioning sets of size at most max_size
  CC <- unlist(sapply(0:max_size, \(x)
    utils::combn(d, x, simplify = FALSE)),
    recursive = FALSE
  )
  n_C <- length(CC)

  ### Number of directed edges, excluding diagonal
  n_d <- d * (d - 1)

  ### Number of bi-directed edges, excluding diagonal
  n_b <- (d * (d - 1)) / 2

  ### tilde(n)
  n_tilde <- ifelse(d < 4, d, 2 * d - 3) - 1
  nkc <- c(
    "K1" = n_d * n_C,
    "K2" = n_d * max(d - 2, 0) * n_C,
    "K3a" = n_d * max(d - 2, 0) * n_C,
    # "K3b" = n_d * max(d - 3, 0) * max(d - 2, 0) * n_C,
    "K4" = n_d * max(d - 2, 0) * n_C
  )

  # Array indices for path indicators
  dic <- data.frame()
  counter <- 1
  for (i in seq_len(d)) {
    for (C in seq_along(CC)) {
      dic <- rbind(
        dic, data.frame(i = i, C = C, iinc = i %in% CC[[C]], idx = counter)
      )
      counter <- counter + 1
    }
  }
  ndic <- NROW(dic)

  # Array indices for d-separation statements
  zijc <- data.frame()
  counter <- 1
  for (i in seq_len(d)) {
    for (j in seq_len(i - 1)) {
      for (C in seq_along(CC)) {
        if (all(!(c(i, j) %in% CC[[C]]))) {
          zijc <- rbind(
            zijc,
            data.frame(i = i, j = j, C = C, idx = counter)
          )
          counter <- counter + 1
        }
      }
    }
  }

  # Array indices for bidirected length variables 
  lbijc <- data.frame()
  counter <- 1
  for (i in seq_len(d)) {
    for (j in seq_len(i - 1)) {
      for (C in seq_along(CC)) {
        lbijc <- rbind(
          lbijc,
          data.frame(i = i, j = j, C = C, idx = counter)
        )
        counter <- counter + 1
      }
    }
  }

  # Array indices for bidirected length variables 
  ldijc <- data.frame()
  counter <- 1
  for (i in seq_len(d)) {
    for (j in setdiff(seq_len(d), i)) {
      for (C in seq_along(CC)) {
        if (!(i %in% CC[[C]]) && (j %in% CC[[C]])) {
          ldijc <- rbind(
            ldijc,
            data.frame(i = i, j = j, C = C, idx = counter)
          )
          counter <- counter + 1
        }
      }
    }
  }

  # Array indices for directed edges
  xij <- data.frame(expand.grid(i = seq_len(d), j = seq_len(d))) |>
    dplyr::filter(i != j)
  xij$idx <- 1:NROW(xij)

  # Array indices for bidirected edges
  bxij <- data.frame(expand.grid(i = seq_len(d), j = seq_len(d))) |>
    dplyr::filter(i < j)
  bxij$idx <- 1:NROW(bxij)
  tmp <- bxij
  colnames(tmp) <- c("j", "i", "idx")
  bxij <- rbind(bxij, tmp)

  ### Number of separation variables, excluding i,j in conditioning set
  n_z <- NROW(zijc)
  n_lb <- NROW(lbijc)

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

  p <- trafo(praw)

  # Linearize absolute value:
  # min w |z - b|
  # t >= z - b
  # t >= -(z - b)
  # min wt

  ### SETUP GUROBI MODEL
  model <- list()
  model$obj <- c(
    w, # corresponds to t_{ij}^C
    "zijc" = rep(0, n_z), # corresponds to z_{ij}^C
    "dxij" = rep(0, n_d), # corresponds to x^{->}_{ij}
    rep(0, n_z), # corresponds to l_{ij}^C
    rep(0, n_d), # corresponds to l_{ij}^{->}
    rep(0, sum(nN1 <- c(n_d, (d - 2) * n_d))), # For min-constraint N1
    rep(0, n_d), # corresponds to d_ij^->
    rep(0, sum(nkc)), # aux P1
    "bxij" = rep(0, n_b), # bidirected edges
    "sxij" = rep(0, n_d), # *directed edges
    rep(0, n_lb), # corresponds to l_{ij}^{<->,C}
    rep(0, n_lb), # corresponds to z_{ij}^{<->,C}
    rep(0, n_lb + n_lb * (d - 2)), # corresponds to rhs of F1 and F2
    "diC" = rep(0, d * n_C), # diC
    "ldijc" = rep(0, nrow(ldijc)),
    rep(0, (d - 1) * nrow(ldijc)),
    0,
    0
  )
  model$branchpriority <- rep(0, length(model$obj))
  model$branchpriority[.multigrep(c("xij", "bxij", "sxij"), names(model$obj))] <- 1
  ### Warm start and edge hints
  guess <- ws_xij <- numeric(n_d)
  bguess <- ws_bxij <- numeric(n_b)
  sguess <- ws_sxij <- numeric(n_d)
  for (i in seq_len(d)) {
    for (j in setdiff(seq_len(d), i)) {
      bidx <- bxij$idx[bxij$i == i & bxij$j == j]
      xidx <- xij$idx[xij$i == i & xij$j == j]
      zidx <- c(zijc$idx[zijc$i == i & zijc$j == j], zijc$idx[zijc$i == j & zijc$j == i])
      if (is.null(edgehints)) {
        guess[xidx] <- as.integer(max(p[zidx]))
        bguess[bidx] <- as.integer(max(p[zidx]))
        sguess[xidx] <- as.integer(max(p[zidx]))
      } else {
        guess[xidx] <- edgehints[i, j]
        sguess[xidx] <- edgehints[i, j]
        bguess[bidx] <- 1 * (edgehints[i, j] | edgehints[j, i])
      }
      if (is.null(warmstart)) {
        ws_xij[xidx] <- 0
        ws_bxij[bidx] <- 0
        ws_sxij[xidx] <- 0
      } else {
        ws_xij[xidx] <- warmstart$M1[i, j]
        ws_bxij[bidx] <- warmstart$M2[i, j]
        ws_sxij[xidx] <- 1 * (warmstart$M1[i, j] | warmstart$M2[i, j])
      }
    }
  }
  model$varhintval <- rep(NA, length(model$obj))
  model$varhintval[grep("zijc", names(model$obj))] <- p
  model$varhintval[grep("dxij", names(model$obj))] <- guess
  model$varhintval[grep("bxij", names(model$obj))] <- bguess
  model$varhintval[grep("sxij", names(model$obj))] <- sguess
  model$varhintpri <- rep(0, length(model$obj))
  model$varhintpri[grep("zijc", names(model$obj))] <- 2
  model$varhintpri[grep("dxij", names(model$obj))] <- 2 - guess
  model$varhintpri[grep("bxij", names(model$obj))] <- 2 - bguess
  model$varhintpri[grep("sxij", names(model$obj))] <- 2 - sguess
  model$varhintpri <- as.integer(model$varhintpri)
  model$start <- rep(NA, length(model$obj))
  model$start[.multigrep(c("dxij", "bxij", "sxij"), names(model$obj))] <- c(ws_xij, ws_bxij, ws_sxij)
  model$modelsense <- "min"
  model$vtype <- rep(c("C", "B", "I", "B", "I", "B", "I", "B", "I", "B", "I", "I", "I", "B"),
    c(n_z, n_z + n_d, n_z + n_d + sum(nN1), n_d, sum(nkc), n_b + n_d,
      n_lb, n_lb, n_lb + n_lb * (d - 2), ndic, nrow(ldijc), (d - 1) * nrow(ldijc), 1, 1))
  model$lb <- rep(c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, d, 1),
    c(2 * n_z + n_d, n_z + n_d + sum(nN1), n_d, sum(nkc), n_b + n_d,
      n_lb, n_lb, n_lb + n_lb * (d - 2), ndic, nrow(ldijc), (d - 1) * nrow(ldijc), 1, 1))
  model$ub <- rep(c(Inf, 1, Inf, d, 2 * d, 1, 4 * d, 1, d, 1, 2 * (d - 1) + 1, 1, d, 2 * d - 1, d, 1),
    c(n_z, n_z + n_d, n_z, n_d, sum(nN1), n_d, sum(nkc), n_b + n_d,
      n_lb, n_lb, n_lb + n_lb * (d - 2), ndic, nrow(ldijc), (d - 1) * nrow(ldijc), 1, 1))
  model$rhs <- c(
    "labs" = -p, # linearize objective
    "labs" = p, # linearize objective
    "d1d2min" = rep(0, sum(nN1)), # D1, D2 auxiliary N1
    "acyc" = rep(0, n_d), # E1
    "indic" = rep(c(d - 1, -1), each = n_d), # dij->
    "p1min" = rep(0, sum(nkc)), # aux P1
    "ZLcons" = rep(d, n_z), # (C4) in the writeup
    "ZLcons" = rep(- d, n_z), # (C5) in the writeup
    "C6to9" = rep(0, 3 * d * (d - 1)), # C6to9
    "C10to11" = rep(c(d, -d), each = n_lb), # C10 to 11
    "f1f2min" = rep(0, n_lb + n_lb * (d - 2)), # Min O1 aux
    "r1b" = rep(0, nr1b <- d * sum(sapply(seq_len(max_size) - 1, \(x) {
      choose(d, x)
    }))),
    "redundant" = rep(0, d * (d - 1) / 2),
    "g1min" = rep(0, (d - 1) * nrow(ldijc)) # E1 + E2 u-variable constraints
  )
  model$sense <- c(
    rep(">=", 2 * n_z), # linearize objective
    rep("=", sum(nN1)), # For min-constraint N1
    rep(">=", n_d), # For acyclicity E1
    rep("<=", 2 * n_d), # dij->
    rep("=", sum(nkc)), # aux P1
    rep("<=", 2 * n_z), # Z, L consistency
    rep("<=", 3 * d * (d - 1)), # C6to9
    rep("<=", 2 * n_lb), # C10to11
    rep("=", n_lb + n_lb * (d - 2)), # Min O1 aux
    rep("=", nr1b),
    rep("<=", d * (d - 1) / 2),
    rep("=", (d - 1) * nrow(ldijc))
  )

  if (cache && all(file.exists(fn))) {
    cat("\nUsing cached constraint matrix and RHS...")
    A <- Matrix::readMM(fn[1])
    rhs <- readRDS(fn[2])
    rhs[1:(2*n_z)] <- c(-p, p)
    P1 <- readRDS(fn[3])
    N1 <- readRDS(fn[4])
    O1 <- readRDS(fn[5])
    R1 <- readRDS(fn[6])
    G1 <- readRDS(fn[7])
    model$rhs <- rhs
  } else {
    A <- MatrixExtra::emptySparse(
      nrow = length(model$rhs), ncol = length(model$obj), format = "C"
    )

    colnames(A) <- make.unique(rep(
      c("tijC", "zijC", "dxij", "lijc", "leij", "nijk", "deij", "uijklm", "bxij", "sxij", "lbijc", "zbijc", "fijc", "diC", "ldijC", "uldijkC", "constd", "const1"),
      c(n_z, n_z, n_d, n_z, n_d, sum(nN1), n_d, sum(nkc), n_b, n_d, n_lb, n_lb, n_lb + n_lb * (d - 2), ndic, nrow(ldijc), (d - 1) * nrow(ldijc), 1, 1)
    ))
    rownames(A) <- make.unique(c(
      rep("labs", 2 * n_z),
      rep("minN1", sum(nN1)),
      rep("acyc", n_d),
      rep("indic", 2 * n_d),
      rep("p1min", sum(nkc)),
      rep("ZLcons", 2 * n_z),
      rep("C6to9", 3 * d * (d - 1)),
      rep("C10to11", 2 * n_lb),
      rep("f1f2min", n_lb + n_lb * (d - 2)),
      rep("r1b", nr1b),
      rep("redundant", d * (d - 1) / 2),
      rep("g1min", (d - 1) * nrow(ldijc))
    ))

    ### Aux constraints for G1
    rn <- grep("g1min", rownames(A), value = TRUE)
    cn <- .multigrep(c("sxij", "lbijc", "uldijkC"), colnames(A), value = TRUE)
    re1e2 <- model$rhs[grep("g1min", names(model$rhs))]
    clist <- g1l <- list()[rep(1, length(rn))]
    skip <- n_d + n_lb
    cntr <- 1
    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        for (C in seq_len(n_C)) {
          if (!(i %in% CC[[C]]) && (j %in% CC[[C]])) {
            ij <- xij$idx[xij$i == i & xij$j == j]
            clist[[cntr]] <- data.frame(
              i = c(cntr, cntr),
              j = c(skip + cntr, ij),
              v = c(1, d - 1)
            )
            re1e2[cntr] <- d
            g1l[[cntr]] <- data.frame(i = i, j = j, k = NA, C = C, which = "E1", idx = cntr)
            cntr <- cntr + 1
            for (k in setdiff(seq_len(d), c(i, j))) {
              if (k %in% CC[[C]]) {
                kjC <- max(lbijc$idx[lbijc$i == k & lbijc$j == j & lbijc$C == C],
                          lbijc$idx[lbijc$i == j & lbijc$j == k & lbijc$C == C])
                ik <- xij$idx[xij$i == i & xij$j == k]
                clist[[cntr]] <- data.frame(
                  i = c(cntr, cntr, cntr),
                  j = c(skip + cntr, ik, n_d + kjC),
                  v = c(1, d - 2, -1)
                )
                re1e2[cntr] <- d - 1
                g1l[[cntr]] <- data.frame(i = i, j = j, k = k, C = C, which = "E2", idx = cntr)
                cntr <- cntr + 1
              }
            }
          }
        }
      }
    }
    clist <- do.call("rbind", clist)
    g1l <- do.call("rbind", g1l)
    A[grep("g1min", rownames(A)), .multigrep(c("sxij", "lbijc", "uldijkC"), colnames(A))] <- 
      Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
        dims = c(length(rn), length(cn)), dimnames = list(rn, cn))
    model$rhs[grep("g1min", names(model$rhs))] <- re1e2
    ng1 <- cntr - 1

    ### R1b
    if (verbose) {
      cat("\nWorking on constraint R1b")
    }

    rn <- grep("r1b", rownames(A), value = TRUE)
    cn <- grep("diC", colnames(A), value = TRUE)
    clist <- list()[rep(1, length(rn))]
    cntr <- 1
    for (i in seq_len(d)) {
      for (C in seq_len(n_C)) {
        if (i %in% CC[[C]]) {
          iC <- dic$idx[dic$i == i & dic$C == C]
          clist[[cntr]] <- data.frame(i = cntr, j = iC, v = 1)
          cntr <- cntr + 1
        }
      }
    }
    clist <- do.call("rbind", clist)
    A[grep("r1b", rownames(A)), grep("diC", colnames(A))] <-
      Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
        dims = c(length(rn), length(cn)), dimnames = list(rn, cn))

    ### RHS F1 F2
    if (max_size > 1) {
      rn <- grep("f1f2min", rownames(A), value = TRUE)
      cn <- .multigrep(c("bxij", "lbijc", "fijc"), colnames(A), value = TRUE)
      rf1f2 <- model$rhs[grep("f1f2min", names(model$rhs))]
      clist <- o1l <- list()[rep(1, length(rn))]
      skip <- n_b + n_lb
      cntr <- 1
      for (j in seq_len(d)) {
        for (i in seq_len(j - 1)) {
          for (C in seq_len(n_C)) {
            if (all(c(i, j) %in% CC[[C]])) {
              ij <- bxij$idx[bxij$i == i & bxij$j == j]
              ijC <- max(lbijc$idx[lbijc$i == i & lbijc$j == j & lbijc$C == C],
                        lbijc$idx[lbijc$i == j & lbijc$j == i & lbijc$C == C])
              clist[[cntr]] <- data.frame(
                i = c(cntr, cntr),
                j = c(skip + cntr, ij),
                v = c(1, d - 1)
              )
              rf1f2[cntr] <- d
              o1l[[cntr]] <- data.frame(i = i, j = j, k = NA, C = C, which = "F1", idx = cntr)
              cntr <- cntr + 1
              for (k in setdiff(seq_len(d), c(i, j))) {
                if (k %in% CC[[C]]) {
                  ikC <- max(lbijc$idx[lbijc$i == i & lbijc$j == k & lbijc$C == C],
                            lbijc$idx[lbijc$i == k & lbijc$j == i & lbijc$C == C])
                  jk <- bxij$idx[bxij$i == j & bxij$j == k]
                  clist[[cntr]] <- data.frame(
                    i = c(cntr, cntr, cntr),
                    j = c(skip + cntr, jk, n_b + ikC),
                    v = c(1, d - 2, -1)
                  )
                  rf1f2[cntr] <- d - 1
                  o1l[[cntr]] <- data.frame(i = i, j = j, k = k, C = C, which = "F2", idx = cntr)
                  cntr <- cntr + 1
                }
              }
            }
          }
        }
      }
      clist <- do.call("rbind", clist)
      o1l <- do.call("rbind", o1l)
      A[grep("f1f2min", rownames(A)), .multigrep(c("bxij", "lbijc", "fijc"), colnames(A))] <- 
        Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
          dims = c(length(rn), length(cn)), dimnames = list(rn, cn))
      model$rhs[grep("f1f2min", names(model$rhs))] <- rf1f2
    }

    ### CONSTRAINTS C10-11
    cm1011 <- rbind(
      cbind(diag(n_lb), diag(n_lb)), # C10
      cbind(-(d - 1) * diag(n_lb), -diag(n_lb)) # C11
    )
    A[grep("C10to11", rownames(A)), .multigrep(c("zbijc", "lbijc"), colnames(A))] <- cm1011

    ### CONSTRAINTS C6-9
    rn <- grep("C6to9", rownames(A), value = TRUE)
    cn <- .multigrep(c("dxij", "bxij", "sxij"), colnames(A), value = TRUE)
    clist <- list()[rep(1, length(rn))]
    cntr <- 1
    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        dij <- xij$idx[xij$i == i & xij$j == j]
        bij <- max(bxij$idx[bxij$i == i & bxij$j == j],
                   bxij$idx[bxij$i == j & bxij$j == i])
        ### C6
        clist[[cntr]] <-
          data.frame(
            i = c(cntr, cntr, cntr),
            j = c(dij, n_d + bij, n_d + n_b + dij),
            v = c(-1, -1, 1)
          )
        cntr <- cntr + 1
        ### C7
        clist[[cntr]] <-
          data.frame(
            i = c(cntr, cntr),
            j = c(dij, n_d + n_b + dij),
            v = c(1, -1)
          )
        cntr <- cntr + 1
        ### C8 + C9
        clist[[cntr]] <-
          data.frame(
            i = c(cntr, cntr),
            j = c(n_d + bij, n_d + n_b + dij),
            v = c(1, -1)
          )
        cntr <- cntr + 1
      }
    }
    clist <- do.call("rbind", clist)
    A[grep("C6to9", rownames(A)), .multigrep(c("dxij", "bxij", "sxij"), colnames(A))] <-
      Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
        dims = c(length(rn), length(cn)), dimnames = list(rn, cn))


    ### CONSTRAINTS FOR Z and L consistency (C4) + (C5)

    cmat_zl <- rbind(
      cbind(diag(n_z), diag(n_z)), # (C4)
      cbind(-diag(n_z) * (d - 1), -diag(n_z)) # (C5)
    )

    A[
      grep("ZLcons", rownames(A)),
      .multigrep(c("zijC", "lijc"), colnames(A))
    ] <- cmat_zl

    ### Auxiliary variables for P1 min constraint

    if (verbose) {
      cat("\nWorking on constraint K1-5")
    }

    rn <- grep("p1min", rownames(A), value = TRUE)
    cn <- .multigrep(c("dxij", "lijc", "sxij", "lbijc", "zbijc", "diC", "ldijC", "uijklm"), colnames(A), value = TRUE)
    skip <- n_d + n_z + n_d + n_lb + n_lb + ndic + nrow(ldijc)
    cntr <- 1
    rhs_m1 <- rep(0, sum(nkc))
    clist <- m1lup <- list()[rep(1, sum(nkc))]

    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        ### K1
        for (C in seq_len(n_C)) {
          if (all(!(c(i, j) %in% CC[[C]]))) {
            ij <- n_d + n_z + xij$idx[xij$i == i & xij$j == j]
            clist[[cntr]] <- data.frame(
              i = c(cntr, cntr),
              j = c(skip + cntr, ij),
              v = c(1, d - 1)
            )
            rhs_m1[cntr] <- d
            m1lup[[cntr]] <- data.frame(i = i, j = j, k = NA,
              l = NA, m = NA, C = C, idx = cntr, which = "K1")
            cntr <- cntr + 1
          }
        }
        for (k in setdiff(seq_len(d), c(i, j))) {
          ### K2
          for (C in seq_len(n_C)) {
            if (all(!c(i, j) %in% CC[[C]])) {
              if (!(k %in% CC[[C]])) {
                kj <- xij$idx[xij$i == k & xij$j == j]
                ikC <- max(zijc$idx[zijc$i == i & zijc$j == k & zijc$C == C],
                          zijc$idx[zijc$i == k & zijc$j == i & zijc$C == C])
                clist[[cntr]] <- data.frame(
                  i = c(cntr, cntr, cntr),
                  j = c(skip + cntr, kj, n_d + ikC),
                  v = c(1, d - 2, -1)
                )
                rhs_m1[cntr] <- d - 1
                m1lup[[cntr]] <- data.frame(i = i, j = j, k = k,
                  l = NA, m = NA, C = C, idx = cntr, which = "K2")
                cntr <- cntr + 1
                ### K4
                kC <- dic$idx[dic$i == k & dic$C == C]
                ikC <- max(zijc$idx[zijc$i == i & zijc$j == k & zijc$C == C],
                          zijc$idx[zijc$i == k & zijc$j == i & zijc$C == C])
                kjC <- max(zijc$idx[zijc$i == k & zijc$j == j & zijc$C == C],
                          zijc$idx[zijc$i == j & zijc$j == k & zijc$C == C])
                clist[[cntr]] <- data.frame(
                  i = c(cntr, cntr, cntr, cntr),
                  j = c(
                    skip + cntr, 
                    n_d + ikC, 
                    n_d + kjC, 
                    2 * n_d + n_z + 2 * n_lb + kC
                  ),
                  v = c(1, -1, -1, - d + 2)
                )
                rhs_m1[cntr] <- 0
                m1lup[[cntr]] <- data.frame(i = i, j = j, k = k,
                  l = NA, m = NA, C = C, idx = cntr, which = "K4")
                cntr <- cntr + 1
              } else {
                ### K3
                lskip <- n_d + n_z + n_d + n_lb + n_lb + ndic
                dikC <- lskip + ldijc$idx[ldijc$i == i & ldijc$j == k & ldijc$C == C]
                djkC <- lskip + ldijc$idx[ldijc$i == j & ldijc$j == k & ldijc$C == C]
                clist[[cntr]] <- data.frame(
                  i = c(cntr, cntr, cntr),
                  j = c(skip + cntr, dikC, djkC),
                  v = c(1, -1, -1)
                )
                rhs_m1[cntr] <- 0
                m1lup[[cntr]] <- data.frame(i = i, j = j, k = k,
                  l = NA, m = NA, C = C, idx = cntr, which = "K3a")
                cntr <- cntr + 1
                # for (k2 in setdiff(seq_len(d), c(i, j, k))) {
                #   if (k2 %in% CC[[C]]) {
                #     ### K3b
                #     ik <- n_d + n_z + xij$idx[xij$i == i & xij$j == k]
                #     jk2 <- n_d + n_z + xij$idx[xij$i == j & xij$j == k2]
                #     kk2c <- n_d + n_z + n_d
                #     kk2c <- kk2c + max(
                #       lbijc$idx[lbijc$i == k & lbijc$j == k2 & lbijc$C == C],
                #       lbijc$idx[lbijc$i == k2 & lbijc$j == k & lbijc$C == C]
                #     )
                #     clist[[cntr]] <- data.frame(
                #       i = c(cntr, cntr, cntr, cntr, cntr),
                #       j = c(skip + cntr, ik, jk2, kk2c, n_lb + kk2c),
                #       v = c(1, d - 3, d - 3, -1, d - 3)
                #     )
                #     rhs_m1[cntr] <- 3 * d - 7
                #     m1lup[[cntr]] <- data.frame(i = i, j = j, k = k,
                #       l = NA, m = NA, C = C, idx = cntr, which = "K3b")
                #     cntr <- cntr + 1
                #   }
                # }
              }
            }
          }
        }
      }
    }

    nP1 <- cntr - 1

    clist <- do.call("rbind", clist)
    m1lup <- do.call("rbind", m1lup)
    A[grep("p1min", rownames(A)),
      .multigrep(c("dxij", "lijc", "sxij", "lbijc", "zbijc", "diC", "ldijC", "uijklm"), colnames(A))] <- Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
        dims = c(length(rn), length(cn)), dimnames = list(rn, cn))
    model$rhs[grep("p1min", names(model$rhs))] <- rhs_m1

    ### Indicators dij-> (C2) + (C3)

    A[grep("indic", rownames(A)), .multigrep(c("leij", "deij"), colnames(A))] <-
    rbind(
      cbind(diag(n_d), -diag(n_d)), # (C2)
      cbind(-diag(n_d), (d - 1) * diag(n_d)) # (C3)
    )

    ### Swap constraints for ADMGs
    if (mode == "admg") {

      ### DAG1
      if (verbose) {
        cat("\nWorking on constraint DAG1")
      }

      rn <- grep("acyc", rownames(A), value = TRUE)
      cn <- grep("deij", colnames(A), value = TRUE)
      clist <- list()[rep(1, length(rn))]
      cntr <- 1
      for (i in seq_len(d)) {
        for (j in setdiff(seq_len(d), i)) {
          ij <- xij$idx[xij$i == i & xij$j == j]
          ji <- xij$idx[xij$i == j & xij$j == i]
          clist[[cntr]] <- data.frame(
            i = c(cntr, cntr),
            j = c(ij, ji),
            v = c(1, 1)
          )
          model$rhs[grep("acyc", names(model$rhs))[cntr]] <- 1
          cntr <- cntr + 1
        }
      }
      clist <- do.call("rbind", clist)
      A[grep("acyc", rownames(A)), grep("deij", colnames(A))] <-
        Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
          dims = c(length(rn), length(cn)), dimnames = list(rn, cn))

      ### Redundant constraint
      rn <- grep("redundant", rownames(A), value = TRUE)
      cn <- grep("dxij", colnames(A), value = TRUE)
      clist <- list()[rep(1, length(rn))]
      cntr <- 1
      for (i in seq_len(d)) {
        for (j in seq_len(i - 1)) {
          ij <- xij$idx[xij$i == i & xij$j == j]
          ji <- xij$idx[xij$i == j & xij$j == i]
          clist[[cntr]] <- data.frame(
            i = c(cntr, cntr),
            j = c(ij, ji),
            v = c(1, 1)
          )
          model$rhs[grep("redundant", names(model$rhs))[cntr]] <- 1
          cntr <- cntr + 1
        }
      }
      clist <- do.call("rbind", clist)
      A[grep("redundant", rownames(A)), grep("dxij", colnames(A))] <-
        Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
          dims = c(length(rn), length(cn)), dimnames = list(rn, cn))

    }

    ### Auxiliary variables for min constraint N1 (D1), (D2)

    if (verbose) {
      cat("\nWorking on constraints D1-D2")
    }

    rn <- grep("minN1", rownames(A), value = TRUE)
    cn <- .multigrep(c("dxij", "leij", "nijk"), colnames(A), value = TRUE)
    rN1 <- rep(0, sum(nN1))
    cntr <- 1
    skip <- 2 * n_d
    tab_N1 <- clist <- list()[rep(1, sum(nN1))]
    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        ij <- xij$idx[xij$i == i & xij$j == j]
        clist[[cntr]] <- data.frame(
          i = c(cntr, cntr),
          j = c(skip + cntr, ij),
          v = c(1, d - 1)
        )
        rN1[cntr] <- d
        tab_N1[[cntr]] <- data.frame(i = i, j = j, idx = cntr, name = "D1")
        cntr <- cntr + 1
        for (k in setdiff(seq_len(d), c(i, j))) {
          ik <- xij$idx[xij$i == i & xij$j == k]
          kj <- xij$idx[xij$i == k & xij$j == j]
          clist[[cntr]] <- data.frame(
            i = c(cntr, cntr, cntr),
            j = c(skip + cntr, kj, n_d + ik),
            v = c(1, d - 2, -1)
          )
          rN1[cntr] <- d - 1
          tab_N1[[cntr]] <- data.frame(i = i, j = j, idx = cntr, name = "D2")
          cntr <- cntr + 1
        }
      }
    }

    tab_N1 <- do.call("rbind", tab_N1)
    clist <- do.call("rbind", clist)
    model$rhs[grep("d1d2min", names(model$rhs))] <- rN1
    A[grep("minN1", rownames(A)), .multigrep(c("dxij", "leij", "nijk"), colnames(A))] <-
      Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
        dims = c(length(rn), length(cn)), dimnames = list(rn, cn))

    ### CONSTRAINTS FOR LINEARIZING THE OBJECTIVE

    cmat_t <- rbind(
      cbind(diag(n_z), -diag(n_z)),
      cbind(diag(n_z), diag(n_z))
    )

    A[
      grep("labs", rownames(A)),
      .multigrep(c("tijC", "zijC"), colnames(A))
    ] <- cmat_t

    ### Min constraint N1

    if (verbose) {
      cat("\nWorking on constraint N1")
    }

    N1 <- list()
    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        ij <- xij$idx[xij$i == i & xij$j == j]
        keep <- tab_N1$idx[tab_N1$i == i & tab_N1$j == j]
        N1 <- c(N1, list(list(
            resvar = 3 * n_z + n_d + ij,
            vars = sort(3 * n_z + 2 * n_d + keep)
          ))
        )
      }
    }

    ### Min constraint P1

    if (verbose) {
      cat("\nWorking on constraint P1")
    }

    P1 <- list()
    for (i in seq_len(d)) {
      for (j in seq_len(i - 1)) {
        for (C in seq_len(n_C)) {
          if (all(!c(i, j) %in% CC[[C]])) {
            ijC <- max(zijc$idx[zijc$i == i & zijc$j == j & zijc$C == C],
                      zijc$idx[zijc$i == j & zijc$j == i & zijc$C == C])
            keep <- c(
              m1lup$idx[m1lup$i == i & m1lup$j == j & m1lup$C == C],
              m1lup$idx[m1lup$i == j & m1lup$j == i & m1lup$C == C]
            )
            P1 <- c(P1, list(list(
                resvar = 2 * n_z + n_d + ijC,
                vars = sort(3 * n_z + 3 * n_d + sum(nN1) + keep)
              ))
            )
          }
        }
      }
    }

    ### Min constraint O1

    if (verbose) {
      cat("\nWorking on constraint O1")
    }

    O1 <- list()
    if (max_size > 1) {
      for (j in seq_len(d)) {
        for (i in seq_len(j - 1)) {
          for (C in seq_len(n_C)) {
            ijC <- max(lbijc$idx[lbijc$i == i & lbijc$j == j & lbijc$C == C],
                      lbijc$idx[lbijc$i == j & lbijc$j == i & lbijc$C == C])
            keep <- c(
              o1l$idx[o1l$i == i & o1l$j == j & o1l$C == C],
              o1l$idx[o1l$i == j & o1l$j == i & o1l$C == C]
            )
            O1 <- c(O1, list(list(
                resvar = 3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + ijC,
                vars = c(sort(3 * n_z + 2 * n_lb + 4 * n_d + sum(nN1) + sum(nkc) + n_b + keep), length(model$obj) - 1)
              ))
            )
          }
        }
      }
    }

    ### Min constraint G1

    if (verbose) {
      cat("\nWorking on constraint G1")
    }

    G1 <- list()
    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        for (C in seq_len(n_C)) {
          if (!(i %in% CC[[C]]) && (j %in% CC[[C]])) {
            ijC <- ldijc$idx[ldijc$i == i & ldijc$j == j & ldijc$C == C]
            keep <- g1l$idx[g1l$i == i & g1l$j == j & g1l$C == C]
            G1 <- c(G1, list(list(
                resvar = 3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + 
                3 * n_lb + n_lb * (d - 2) + d * n_C + ijC,
                vars = sort(3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + 
                3 * n_lb + n_lb * (d - 2) + d * n_C + nrow(ldijc) + keep)
              ))
            )
          }
        }
      }
    }

    ### Min constraint R1a

    if (verbose) {
      cat("\nWorking on constraint R1")
    }

    R1 <- list()
    for (i in seq_len(d)) {
      for (C in seq_len(n_C)) {
        if (!i %in% CC[[C]]) {
          ic <- dic$idx[dic$i == i & dic$C == C]
          keep <- c()
          for (k in CC[[C]]) {
            keep <- c(keep, xij$idx[xij$i == i & xij$j == k])
          }
          R1 <- c(R1, list(list(
              resvar = 3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + 
                3 * n_lb + n_lb * (d - 2) + ic,
              vars = c(sort(3 * n_z + 2 * n_d + sum(nN1) + keep), 
              length(model$obj))
            ))
          )
        }
      }
    }

    ### Cache

    if (cache) {
      if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE)
      }
      Matrix::writeMM(A, fn[1])
      saveRDS(model$rhs, fn[2])
      saveRDS(P1, fn[3])
      saveRDS(N1, fn[4])
      saveRDS(O1, fn[5])
      saveRDS(R1, fn[6])
      saveRDS(G1, fn[7])
    }
  }

  ### INCLUDE CONSTRAINTS
  model$A <- A
  model$genconmin <- c(N1, P1, O1, R1, G1)
  # print(gurobi::gurobi_iis(model))

  if (verbose) {
    cat("\nModel setup done. Solving now...\n")
  }

  ### SOLVE
  if (requireNamespace("gurobi")) {
    sol <- gurobi::gurobi(model, gurobi_args)
  }

  ### Convert solution to DMG
  .to_dmg <- function(x) {
    edge <- x[2 * n_z + 1:n_d]
    zijC <- x[n_z + 1:n_z]
    lijC <- x[2 * n_z + n_d + 1:n_z]
    leij <- x[3 * n_z + n_d + 1:n_d]
    nijk <- x[3 * n_z + 2 * n_d + 1:sum(nN1)]
    deij <- x[3 * n_z + 2 * n_d + sum(nN1) + 1:n_d]
    bedge <- x[3 * n_z + 3 * n_d + sum(nN1) + sum(nkc) + 1:n_b]
    sedge <- x[3 * n_z + 3 * n_d + sum(nN1) + sum(nkc) + n_b + 1:n_d]
    bminlen <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + 1:n_lb]
    bminleni <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + n_lb + 1:n_lb]
    dic$dic <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + 3 * n_lb + n_lb * (d - 2) + 1:nrow(dic)]
    ldijc$ldijc <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + 3 * n_lb + n_lb * (d - 2) + nrow(dic) + 1:nrow(ldijc)]
    ue1e2 <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nkc) + n_b + 3 * n_lb + n_lb * (d - 2) + nrow(dic) + nrow(ldijc) + 1:((d - 1) * nrow(ldijc))]
    M1 <- M2 <- matrix(0, nrow = d, ncol = d)
    dimnames(M1) <- dimnames(M2) <- list(V, V)
    for (i in seq_len(d)) {
      for (j in seq_len(d)) {
        if (i != j) {
          M1[i, j] <- edge[xij$idx[xij$i == i & xij$j == j]]
          M2[i, j] <- bedge[bxij$idx[bxij$i == i & bxij$j == j]]
        }
      }
    }
    structure(
      list(M1 = M1, M2 = M2), edge = edge, bedge = bedge,
      dcon = zijC, antlen = leij, minlen = lijC,
      pind = deij, sedge = sedge, bminlen = bminlen,
      bminleni = bminleni, dic = dic, ldijc = ldijc
    )
  }

  merged$dcon <- sol$x[n_z + 1:n_z][merged$idx]
  merged$minlen <- sol$x[2 * n_z + n_d + 1:n_z][merged$idx]

  ### RETURN
  structure(
    list(graph = .to_dmg(sol$x), tests = merged, optim = sol),
    class = "graphopt"
  )
}

