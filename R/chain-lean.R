chain_lean_optim <- function(
    tests, d = 3, max_size = d - 2, V = letters[1:d],
    trafo = \(x) as.numeric(x <= 0.05),
    weight_type = c("const", "inv", "log", "size"),
    warmstart = NULL,
    edgehints = NULL,
    gurobi_args = list(),
    verbose = FALSE,
    cache = TRUE,
    cache_dir = "./.cache-chain-lean-redundant",
    ...) {

  if (!requireNamespace("gurobi")) {
    warning("Solver `gurobi` not available.")
  }

  weight_type <- match.arg(weight_type)

  fn <- file.path(
    cache_dir,
    paste0(c("lhs", "rhs", "M1", "N1", "W1", "Z1", "R1"), "dim", d, "max",
      max_size, c(".mtx", ".rds", ".rds", ".rds", ".rds", ".rds", ".rds"))
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

  ### tilde(n)
  n_tilde <- ifelse(d < 4, d, 2 * d - 3) - 1
  nlc <- c(
    "L1" = d * (d - 1) * n_C,
    "L2" = d * (d - 1) * max(d - 2, 0) * n_C,
    "L3" = d * (d - 1) * max(d - 2, 0) * n_C,
    "L4" = d * (d - 1) * max(d - 2, 0) * n_C
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
  zijc

  xij <- data.frame(expand.grid(i = seq_len(d), j = seq_len(d))) |>
    dplyr::filter(i != j)
  xij$idx <- 1:NROW(xij)

  ### Number of separation variables, excluding i,j in conditioning set
  n_z <- NROW(zijc)

  # Array indices for bidirected edges
  luij <- data.frame(expand.grid(i = seq_len(d), j = seq_len(d))) |>
    dplyr::filter(i < j)
  luij$idx <- 1:NROW(luij)
  n_ul <- NROW(luij)
  tmp <- luij
  colnames(tmp) <- c("j", "i", "idx")
  luij <- rbind(luij, tmp)

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
    "xij" = rep(0, n_d), # corresponds to x^{->}_{ij}
    rep(0, n_z), # corresponds to l_{ij}^C
    rep(0, n_d), # corresponds to l_{ij}^{->}
    rep(0, sum(nN1 <- c(n_d, (d - 2) * n_d))), # For min-constraint N1
    rep(0, n_d), # corresponds to d_ij^->
    rep(0, sum(nlc)), # aux M1
    rep(0, n_d), # *-> edges
    rep(0, n_ul), # undirected z variables
    rep(0, n_ul), # undirected l variables
    rep(0, nW1 <- n_ul + n_ul * (d - 2)), # aux W1
    rep(0, nZ1 <- n_d + d * (d - 1) * (d - 2)), # aux Z1
    rep(0, d * n_C), # diC
    0, # aux const for min Z1
    0, # aux const for min R1
    0 # aux const for min W1
  )
  ### Br prio and warmstart
  model$branchpriority <- rep(0, length(model$obj))
  model$branchpriority[grep("xij", names(model$obj))] <- 1
  guess <- rep(0, n_d)
  ws <- rep(0, n_d)
  for (i in seq_len(d)) {
    for (j in setdiff(seq_len(d), i)) {
      xidx <- xij$idx[xij$i == i & xij$j == j]
      zidx <- c(zijc$idx[zijc$i == i & zijc$j == j], zijc$idx[zijc$i == j & zijc$j == i])
      if (is.null(edgehints)) {
        guess[xidx] <- max(p[zidx])
      } else {
        guess[xidx] <- edgehints[i, j]
      }
      if (is.null(warmstart)) {
        ws[xidx] <- 0
      } else {
        ws[xidx] <- warmstart[i, j]
      }
    }
  }
  model$varhintval <- rep(NA, length(model$obj))
  model$varhintval[grep("zijc", names(model$obj))] <- p
  model$varhintval[grep("xij", names(model$obj))] <- guess
  model$varhintpri <- rep(0, length(model$obj))
  model$varhintpri[grep("zijc", names(model$obj))] <- 2
  model$varhintpri[grep("xij", names(model$obj))] <- 2 - guess
  model$varhintpri <- as.integer(model$varhintpri)
  model$start <- rep(NA, length(model$obj))
  model$start[grep("xij", names(model$obj))] <- ws
  model$modelsense <- "min"
  model$vtype <- rep(
    c("C", "B", "I", "B", "I", "B", "I", "B", "I"),
    c(n_z, n_z + n_d, n_z + n_d + sum(nN1), n_d, sum(nlc),
      n_d + n_ul, n_ul + nW1 + nZ1, d * n_C, 3)
  )
  model$lb <- rep(
    c(0, 1, 0, 1, 0, 1, -2, 0, 0, 1, d),
    c(2 * n_z + n_d, n_z + n_d + sum(nN1), n_d, sum(nlc),
      n_d + n_ul, n_ul + nW1, nZ1, d * n_C, 1, 1, 1)
  )
  model$ub <- rep(
    c(Inf, 1, Inf, d, 2 * d, 1, 3 * d, 1, d, 3 * d - 1, 1, 1, 0, 1, d),
    c(n_z, n_z + n_d, n_z, n_d, sum(nN1), n_d, sum(nlc),
      n_d + n_ul, n_ul, nW1, nZ1, d * n_C, 1, 1, 1)
  )
  model$rhs <- c(
    "labs" = -p, # linearize objective
    "labs" = p, # linearize objective
    "d1d2min" = rep(0, sum(nN1)), # D1, D2 auxiliary N1
    "acyc" = rep(0, 2 * n_d), # CH1ab
    "indic" = rep(c(d - 1, -1), each = n_d), # dij->
    "m1min" = rep(0, sum(nlc)), # aux M1
    "ZLcons" = rep(d, n_z), # (C4) in the writeup
    "ZLcons" = rep(- d, n_z), # (C5) in the writeup
    "w1min" = rep(0, nW1),
    "z1max" = rep(0, nZ1),
    "unZL" = rep(d, n_ul), # C12-13
    "unZL" = rep(-d, n_ul), # C12-13
    "r1b" = rep(0, nr1b <- d * sum(sapply(seq_len(max_size) - 1, \(x) {
      choose(d, x)
    })))
  )
  model$sense <- c(
    rep(">=", 2 * n_z), # linearize objective
    rep("=", sum(nN1)), # For min-constraint N1
    rep(">=", 2 * n_d), # For acyclicity E1
    rep("<=", 2 * n_d), # dij->
    rep("=", sum(nlc)), # aux M1
    rep("<=", 2 * n_z), # Z, L consistency
    rep("=", nW1),
    rep("=", nZ1),
    rep("<=", 2 * n_ul),
    rep("=", nr1b)
  )

  if (cache && all(file.exists(fn))) {
    cat("\nUsing cached constraint matrix and RHS...")
    A <- Matrix::readMM(fn[1])
    rhs <- readRDS(fn[2])
    rhs[1:(2*n_z)] <- c(-p, p)
    M1 <- readRDS(fn[3])
    N1 <- readRDS(fn[4])
    W1 <- readRDS(fn[5])
    Z1 <- readRDS(fn[6])
    R1 <- readRDS(fn[7])
    model$rhs <- rhs
  } else {
    A <- MatrixExtra::emptySparse(
      nrow = length(model$rhs), ncol = length(model$obj), format = "C"
    )

    colnames(A) <- make.unique(rep(
      c("tijC", "zijC", "dxij", "lijc", "leij", "nijk", "deij", "uijklm",
        "sxij", "zuij", "luij", "wijk", "zijk", "diC", "const0", "const1",
        "constd"),
      c(n_z, n_z, n_d, n_z, n_d, sum(nN1), n_d, sum(nlc), n_d, n_ul,
        n_ul, nW1, nZ1, d * n_C, 1, 1, 1)
    ))
    rownames(A) <- make.unique(c(
      rep("labs", 2 * n_z),
      rep("minN1", sum(nN1)),
      rep("acyc", 2 * n_d),
      rep("indic", 2 * n_d),
      rep("m1min", sum(nlc)),
      rep("ZLcons", 2 * n_z),
      rep("w1min", nW1),
      rep("z1max", nZ1),
      rep("unZL", n_ul),
      rep("unZL", n_ul),
      rep("r1b", nr1b)
    ))

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

    ### Auxiliary Y1 and Y2 for Z1
    if (verbose) {
      cat("\nWorking on constraint Y1-2")
    }

    rn <- grep("z1max", rownames(A), value = TRUE)
    cn <- .multigrep(c("dxij", "zuij", "zijk"), colnames(A), value = TRUE)
    skip <- n_d + n_ul
    cntr <- 1
    rz1 <- rep(0, nZ1)
    z1l <- clist <- list()[rep(1, nZ1)]
    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        ### Y1
        ij <- xij$idx[xij$i == i & xij$j == j]
        ji <- xij$idx[xij$i == j & xij$j == i]
        clist[[cntr]] <- data.frame(
          i = c(cntr, cntr, cntr),
          j = c(skip + cntr, ij, ji),
          v = c(1, -1, 1)
        )
        rz1[cntr] <- 0
        z1l[[cntr]] <- data.frame(i = i, j = j, k = NA, idx = cntr, which = "Y1")
        cntr <- cntr + 1
        for (k in setdiff(seq_len(d), c(i, j))) {
          ### Y2
          ik <- xij$idx[xij$i == i & xij$j == k]
          ki <- xij$idx[xij$i == k & xij$j == i]
          kj <- n_d + luij$idx[luij$i == k & luij$j == j]
          clist[[cntr]] <- data.frame(
            i = c(cntr, cntr, cntr, cntr),
            j = c(skip + cntr, ik, ki, kj),
            v = c(1, -1, 1, -1)
          )
          rz1[cntr] <- -1
          z1l[[cntr]] <- data.frame(i = i, j = j, k = k, idx = cntr, which = "Y2")
          cntr <- cntr + 1
        }
      }
    }
    clist <- do.call("rbind", clist)
    z1l <- do.call("rbind", z1l)
    model$rhs[grep("z1max", names(model$rhs))] <- rz1
    A[grep("z1max", rownames(A)), .multigrep(c("dxij", "zuij", "zijk"), colnames(A))] <-
      Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
        dims = c(length(rn), length(cn)), dimnames = list(rn, cn))
    nz1 <- cntr - 1

    ### Auxiliary U1 U2 for W1
    if (verbose) {
      cat("\nWorking on constraints U1-2")
    }

    rn <- grep("w1min", rownames(A), value = TRUE)
    cn <- .multigrep(c("dxij", "luij", "wijk"), colnames(A), value = TRUE)
    skip <- n_d + n_ul
    cntr <- 1
    rw1 <- rep(0, nW1)
    w1l <- clist <- list()[rep(1, nW1)]
    for (i in seq_len(d)) {
      for (j in seq_len(i - 1)) {
        ### U1
        ij <- xij$idx[xij$i == i & xij$j == j]
        ji <- xij$idx[xij$i == j & xij$j == i]
        clist[[cntr]] <- data.frame(
          i = c(cntr, cntr, cntr),
          j = c(skip + cntr, ij, ji),
          v = c(1, d - 1, d - 1)
        )
        rw1[cntr] <- 2 * d - 1
        w1l[[cntr]] <- data.frame(i = i, j = j, k = NA, idx = cntr, which = "U1")
        cntr <- cntr + 1
        ### U2
        for (k in setdiff(seq_len(d), c(i, j))) {
          ik <- n_d + luij$idx[luij$i == i & luij$j == k]
          jk <- xij$idx[xij$i == j & xij$j == k]
          kj <- xij$idx[xij$i == k & xij$j == j]
          clist[[cntr]] <- data.frame(
            i = c(cntr, cntr, cntr, cntr),
            j = c(skip + cntr, ik, jk, kj),
            v = c(1, -1, d - 2, d - 2)
          )
          rw1[cntr] <- 2 * d - 3
          w1l[[cntr]] <- data.frame(i = i, j = j, k = k, idx = cntr, which = "U2")
          cntr <- cntr + 1
        }
      }
    }
    clist <- do.call("rbind", clist)
    w1l <- do.call("rbind", w1l)
    model$rhs[grep("w1min", names(model$rhs))] <- rw1
    A[grep("w1min", rownames(A)), .multigrep(c("dxij", "luij", "wijk"), colnames(A))] <-
      Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
        dims = c(length(rn), length(cn)), dimnames = list(rn, cn))

    ### CONSTRAINTS C12-13

    c1213m <- rbind(
      cbind(diag(n_ul), diag(n_ul)), # (C12)
      cbind(-diag(n_ul) * (d - 1), -diag(n_ul)) # (C13)
    )

    A[
      grep("unZL", rownames(A)),
      .multigrep(c("zuij", "luij"), colnames(A))
    ] <- c1213m

    ### CONSTRAINTS FOR Z and L consistency (C4) + (C5)

    cmat_zl <- rbind(
      cbind(diag(n_z), diag(n_z)), # (C4)
      cbind(-diag(n_z) * (d - 1), -diag(n_z)) # (C5)
    )

    A[
      grep("ZLcons", rownames(A)),
      .multigrep(c("zijC", "lijc"), colnames(A))
    ] <- cmat_zl

    ### Auxiliary variables for M1 min constraint

    if (verbose) {
      cat("\nWorking on constraint L1-5")
    }

    rn <- grep("m1min", rownames(A), value = TRUE)
    cn <- .multigrep(c("dxij", "lijc", "sxij", "diC", "uijklm"), colnames(A), value = TRUE)
    skip <- 2 * n_d + n_z + ndic
    cntr <- 1
    rm1 <- rep(0, sum(nlc))
    m1lup <- clist <- list()[rep(1, sum(nlc))]

    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        ### L1
        for (C in seq_len(n_C)) {
          if (all(!(c(i, j) %in% CC[[C]]))) {
            ij <- xij$idx[xij$i == i & xij$j == j]
            clist[[cntr]] <- data.frame(
              i = c(cntr, cntr),
              j = c(skip + cntr, ij),
              v = c(1, d - 1)
            )
            rm1[cntr] <- d
            m1lup[[cntr]] <- data.frame(i = i, j = j, k = NA,
              l = NA, m = NA, C = C, idx = cntr, which = "L1")
            cntr <- cntr + 1
          }
        }
        for (k in setdiff(seq_len(d), c(i, j))) {
          ### L2
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
                rm1[cntr] <- d - 1
                m1lup[[cntr]] <- data.frame(i = i, j = j, k = k,
                  l = NA, m = NA, C = C, idx = cntr, which = "L2")
                cntr <- cntr + 1
                ### L4
                ikC <- max(zijc$idx[zijc$i == i & zijc$j == k & zijc$C == C],
                          zijc$idx[zijc$i == k & zijc$j == i & zijc$C == C])
                kjC <- max(zijc$idx[zijc$i == k & zijc$j == j & zijc$C == C],
                          zijc$idx[zijc$i == j & zijc$j == k & zijc$C == C])
                kC <- 2 * n_d + n_z + dic$idx[dic$i == k & dic$C == C]
                clist[[cntr]] <- data.frame(
                  i = c(cntr, cntr, cntr, cntr),
                  j = c(skip + cntr, n_d + ikC, n_d + kjC, kC),
                  v = c(1, -1, -1, -(d - 2))
                )
                rm1[cntr] <- 0
                m1lup[[cntr]] <- data.frame(i = i, j = j, k = k,
                  l = NA, m = NA, C = C, idx = cntr, which = "L4")
                cntr <- cntr + 1
              }
              ### L3
              ik <- n_d + n_z + xij$idx[xij$i == i & xij$j == k]
              jk <- n_d + n_z + xij$idx[xij$i == j & xij$j == k]
              kC <- 2 * n_d + n_z + dic$idx[dic$i == k & dic$C == C]
              clist[[cntr]] <- data.frame(
                i = c(cntr, cntr, cntr, cntr),
                j = c(skip + cntr, ik, jk, kC),
                v = c(1, d - 1, d - 1, -(d - 1))
              )
              rm1[cntr] <- 2 * (d - 1) + 1
              m1lup[[cntr]] <- data.frame(i = i, j = j, k = k,
                l = NA, m = NA, C = C, idx = cntr, which = "L3")
              cntr <- cntr + 1
            }
          }
        }
      }
    }

    nM1 <- cntr - 1
    clist <- do.call("rbind", clist)
    m1lup <- do.call("rbind", m1lup)
    A[grep("m1min", rownames(A)), .multigrep(c("dxij", "lijc", "sxij", "diC", "uijklm"), colnames(A))] <-
      Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
        dims = c(length(rn), length(cn)), dimnames = list(rn, cn))
    model$rhs[grep("m1min", names(model$rhs))] <- rm1

    ### Indicators dij-> (C2) + (C3)

    A[grep("indic", rownames(A)), .multigrep(c("leij", "deij"), colnames(A))] <-
    rbind(
      cbind(diag(n_d), -diag(n_d)), # (C2)
      cbind(-diag(n_d), (d - 1) * diag(n_d)) # (C3)
    )

    ### Acyclicity constraint CH1ab
    if (verbose) {
      cat("\nWorking on constraint CH1ab")
    }

    rn <- grep("acyc", rownames(A), value = TRUE)
    cn <- .multigrep(c("leij", "deij"), colnames(A), value = TRUE)
    clist <- list()[rep(1, length(rn))]
    cntr <- 1
    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        ij <- xij$idx[xij$i == i & xij$j == j]
        ji <- xij$idx[xij$i == j & xij$j == i]
        clist[[cntr]] <- data.frame(
          i = c(cntr, cntr, cntr, cntr),
          j = c(ij, ji, n_d + ij, n_d + ji),
          v = c(1, -1, d - 1, d - 1)
        )
        cntr <- cntr + 1
        clist[[cntr]] <- data.frame(
          i = c(cntr, cntr, cntr, cntr),
          j = c(ij, ji, n_d + ij, n_d + ji),
          v = c(-1, 1, d - 1, d - 1)
        )
        cntr <- cntr + 1
      }
    }
    clist <- do.call("rbind", clist)
    A[grep("acyc", rownames(A)), .multigrep(c("leij", "deij"), colnames(A))] <-
      Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
        dims = c(length(rn), length(cn)), dimnames = list(rn, cn))

    ### Auxiliary variables for min constraint N1 (D1), (D2)

    if (verbose) {
      cat("\nWorking on constraints D1-D2")
    }

    rn <- grep("minN1", rownames(A), value = TRUE)
    cn <- .multigrep(c("dxij", "leij", "nijk"), colnames(A), value = TRUE)
    rN1 <- rep(0, sum(nN1))
    cntr <- 1
    skip <- 2 * n_d
    tab_N1 <- clist <- list()[rep(1, length(rn))]
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

    clist <- do.call("rbind", clist)
    tab_N1 <- do.call("rbind", tab_N1)
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

    ### Max constraint Z1

    if (verbose) {
      cat("\nWorking on constraint Z1")
    }

    Z1 <- list()
    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        ij <- xij$idx[xij$i == i & xij$j == j]
        keep <- z1l$idx[z1l$i == i & z1l$j == j]
        Z1 <- c(Z1, list(list(
            resvar = 3 * n_z + 3 * n_d + sum(nN1) + sum(nlc) + ij,
            vars = sort(c(
              3 * n_z + 4 * n_d + sum(nN1) + sum(nlc) + 2 * n_ul + nW1 + keep,
              length(model$obj) - 2
            ))
          ))
        )
      }
    }

    ### Min constraint W1

    if (verbose) {
      cat("\nWorking on constraint W1")
    }

    W1 <- list()
    for (i in seq_len(d)) {
      for (j in seq_len(i - 1)) {
        ij <- luij$idx[luij$i == i & luij$j == j]
        keep <- w1l$idx[w1l$i == i & w1l$j == j]
        W1 <- c(W1, list(list(
            resvar = 3 * n_z + 4 * n_d + sum(nN1) + sum(nlc) + n_ul + ij,
            vars = c(sort(3 * n_z + 4 * n_d + sum(nN1) + sum(nlc) + 2 * n_ul + keep), length(model$obj))
          ))
        )
      }
    }

    ### Min constraint R1

    if (verbose) {
      cat("\nWorking on constraint R1")
    }

    R1 <- list()
    for (i in seq_len(d)) {
      for (C in seq_len(n_C)) {
        if (!i %in% CC[[C]]) {
          ic <- dic$idx[dic$i == i & dic$C == C]
          keep <- c()
          for (k in seq_len(d)) {
            if (k %in% CC[[C]]) {
              keep <- c(keep, xij$idx[xij$i == i & xij$j == k])
            }
          }
          R1 <- c(R1, list(list(
              resvar = 3 * n_z + 4 * n_d + sum(nN1) + sum(nlc) + 2 * n_ul + nW1 + nZ1 + ic,
              vars = c(sort(3 * n_z + 2 * n_d + sum(nN1) + keep), length(model$obj) - 1)
            ))
          )
        }
      }
    }

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

    ### Min constraint M1

    if (verbose) {
      cat("\nWorking on constraint M1")
    }

    M1 <- list()
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
            M1 <- c(M1, list(list(
                resvar = 2 * n_z + n_d + ijC,
                vars = sort(3 * n_z + 3 * n_d + sum(nN1) + keep)
              ))
            )
          }
        }
      }
    }

    if (cache) {
      if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE)
      }
      Matrix::writeMM(A, fn[1])
      saveRDS(model$rhs, fn[2])
      saveRDS(M1, fn[3])
      saveRDS(N1, fn[4])
      saveRDS(W1, fn[5])
      saveRDS(Z1, fn[6])
      saveRDS(R1, fn[7])
    }
  }

  ### INCLUDE CONSTRAINTS
  model$A <- A
  model$genconmin <- c(N1, M1, W1, R1)
  model$genconmax <- Z1

  if (verbose) {
    cat("\nModel setup done. Solving now...\n")
  }

  ### SOLVE
  if (requireNamespace("gurobi")) {
    # return(gurobi::gurobi_iis(model))
    sol <- gurobi::gurobi(model, gurobi_args)
  }

  ### Convert solution to DMG
  .to_dag <- function(x) {
    edge <- x[2 * n_z + 1:n_d]
    zijC <- x[n_z + 1:n_z]
    lijC <- x[2 * n_z + n_d + 1:n_z]
    leij <- x[3 * n_z + n_d + 1:n_d]
    nijk <- x[3 * n_z + 2 * n_d + 1:sum(nN1)]
    deij <- x[3 * n_z + 2 * n_d + sum(nN1) + 1:n_d]
    sedge <- x[3 * n_z + 3 * n_d + sum(nN1) + sum(nlc) + 1:n_d]
    uz <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nlc) + 1:n_ul]
    ul <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nlc) + n_ul + 1:n_ul]
    diC <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nlc) +
      nW1 + nZ1 + 2 * n_ul + 1:ndic]
    # z1l$learned <- x[3 * n_z + 4 * n_d + sum(nN1) + sum(nlc) + nW1 + 2 * n_ul + 1:nz1]
    # print(z1l)
    dag <- matrix(0, nrow = d, ncol = d)
    dimnames(dag) <- list(V, V)
    for (i in seq_len(d)) {
      for (j in seq_len(d)) {
        if (i != j) {
          dag[i, j] <- edge[xij$idx[xij$i == i & xij$j == j]]
        }
      }
    }
    structure(
      dag, class = class(dag), edge = edge, dcon = zijC,
      antlen = leij, minlen = lijC, pind = deij,
      sedge = sedge, zuij = uz, luij = ul, diC = diC
    )
  }

  merged$dcon <- sol$x[n_z + 1:n_z][merged$idx]
  merged$minlen <- sol$x[2 * n_z + n_d + 1:n_z][merged$idx]

  ### RETURN
  structure(
    list(graph = .to_dag(sol$x), tests = merged, optim = sol),
    class = "graphopt"
  )
}

