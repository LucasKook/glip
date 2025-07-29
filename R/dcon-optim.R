dcon_optim <- function(
    tests, d = 3, max_size = d - 2, V = letters[1:d],
    trafo = \(x) as.numeric(x <= 0.05),
    weight_type = c("const", "inv", "log"),
    gurobi_args = list(),
    verbose = FALSE, 
    cache = TRUE,
    cache_dir = "./.cache-dcon",
    ...) {

  if (!requireNamespace("gurobi")) {
    warning("Solver `gurobi` not available.")
  }

  weight_type <- match.arg(weight_type)

  fn <- file.path(
    cache_dir,
    paste0(c("lhs", "rhs", "M1", "N1", "R1"), "dim", d, "max",
      max_size, c(".mtx", ".rds", ".rds", ".rds", ".rds"))
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
  nlc <- c(
    "L1" = d * (d - 1) * n_C,
    "L2" = d * (d - 1) * max(d - 2, 0) * n_C,
    "L3" = d * (d - 1) * max(d - 2, 0) * n_C,
    "L4" = d * (d - 1) * max(d - 2, 0) * max(d - 3, 0) * n_C,
    "L5" = d * (d - 1) * max(d - 2, 0) * max(d - 3, 0)^2 * n_C
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
  p <- rep(0, n_z)
  merged <- dplyr::left_join(tests, lookup, by = c("X", "Y", "Z"))
  p[merged$idx] <- merged$p.value
  w[merged$idx] <- 1
  merged$weight <- w[merged$idx]
  p <- trafo(p)

  w[w == 1] <- switch(weight_type,
    "const" = 1,
    "inv" = 1 / pmax(p[w == 1], 0.0001),
    "log" = -log2(pmax(p[w == 1], 2 * .Machine$double.eps))
  )

  # Linearize absolute value:
  # min w |z - b|
  # t >= z - b
  # t >= -(z - b)
  # min wt

  ### SETUP GUROBI MODEL
  model <- list()
  model$obj <- c(
    w, # corresponds to t_{ij}^C
    rep(0, n_z), # corresponds to z_{ij}^C
    rep(0, n_d), # corresponds to x^{->}_{ij}
    rep(0, n_z), # corresponds to l_{ij}^C
    rep(0, n_d), # corresponds to l_{ij}^{->}
    rep(0, sum(nN1 <- c(n_d, (d - 2) * n_d))), # For min-constraint N1
    rep(0, n_d), # corresponds to d_ij^->
    rep(0, sum(nlc)), # aux M1
    rep(0, d * n_C), # diC
    0 # aux const for min R1
  )
  model$modelsense <- "min"
  model$vtype <- rep(
    c("C", "B", "I", "B", "I", "B", "B"),
    c(n_z, n_z + n_d, n_z + n_d + sum(nN1), n_d, sum(nlc), d * n_C, 1)
  )
  model$lb <- rep(
    c(0, 1, 0, 1, 0, 1),
    c(2 * n_z + n_d, n_z + n_d + sum(nN1), n_d, sum(nlc), d * n_C, 1)
  )
  model$ub <- rep(
    c(Inf, 1, Inf, d, 2 * d, 1, 5 * d - 8, 1, 1),
    c(n_z, n_z + n_d, n_z, n_d, sum(nN1), n_d, sum(nlc), d * n_C, 1)
  )
  model$rhs <- c(
    "labs" = -p, # linearize objective
    "labs" = p, # linearize objective
    "d1d2min" = rep(0, sum(nN1)), # D1, D2 auxiliary N1
    "acyc" = rep(0, 2 * n_d), # CH1ab
    "indic" = rep(c(d - 1, -1), each = n_d), # dij->
    "m1min" = rep(0, sum(nlc)), # aux M1
    "ZLcons" = rep(d, n_z), # (C4*) in the writeup
    "ZLcons" = rep(- d, n_z), # (C5*) in the writeup
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
    rep("=", nr1b)
  )

  if (cache && all(file.exists(fn))) {
    cat("\nUsing cached constraint matrix and RHS...")
    A <- Matrix::readMM(fn[1])
    rhs <- readRDS(fn[2])
    rhs[1:(2*n_z)] <- c(-p, p)
    M1 <- readRDS(fn[3])
    N1 <- readRDS(fn[4])
    R1 <- readRDS(fn[5])
    model$rhs <- rhs
  } else {
    A <- MatrixExtra::emptySparse(
      nrow = length(model$rhs), ncol = length(model$obj)
    )

    colnames(A) <- make.unique(rep(
      c("tijC", "zijC", "dxij", "lijc", "leij", "nijk", "deij", "uijklm", "diC", "const1"),
      c(n_z, n_z, n_d, n_z, n_d, sum(nN1), n_d, sum(nlc), d * n_C, 1)
    ))
    rownames(A) <- make.unique(c(
      rep("labs", 2 * n_z),
      rep("minN1", sum(nN1)),
      rep("acyc", 2 * n_d),
      rep("indic", 2 * n_d),
      rep("m1min", sum(nlc)),
      rep("ZLcons", 2 * n_z),
      rep("r1b", nr1b)
    ))

    ### R1b
    if (verbose) {
      cat("\nWorking on constraint R1b")
    }

    tmp <- A[grep("r1b", rownames(A)), grep("diC", colnames(A))]
    cntr <- 1
    for (i in seq_len(d)) {
      for (C in seq_len(n_C)) {
        if (i %in% CC[[C]]) {
          iC <- dic$idx[dic$i == i & dic$C == C]
          tmp[cntr, iC] <- 1 + .fill(tmp[cntr, iC])
          cntr <- cntr + 1
        }
      }
    }
    A[grep("r1b", rownames(A)), grep("diC", colnames(A))] <- tmp

    ### CONSTRAINTS FOR Z and L consistency (C4*) + (C5*)
  
    cmat_zl <- rbind(
      cbind(diag(n_z), diag(n_z)), # (C4*)
      cbind(-diag(n_z) * (d - 1), -diag(n_z)) # (C5*)
    )

    A[
      grep("ZLcons", rownames(A)),
      .multigrep(c("zijC", "lijc"), colnames(A))
    ] <- cmat_zl

    ### Auxiliary variables for M1 min constraint

    if (verbose) {
      cat("\nWorking on constraint L1-5")
    }

    tmp <- A[
      grep("m1min", rownames(A)), 
      .multigrep(c("dxij", "lijc", "diC", "uijklm"), colnames(A))
    ]
    skip <- n_d + n_z + ndic
    cntr <- 1
    rm1 <- rep(0, sum(nlc))
    m1lup <- data.frame()

    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        ### L1
        for (C in seq_len(n_C)) {
          if (all(!(c(i, j) %in% CC[[C]]))) {
            ij <- xij$idx[xij$i == i & xij$j == j]
            tmp[cntr, skip + cntr] <- 1 + .fill(tmp[cntr, skip + cntr])
            tmp[cntr, ij] <- (d - 1) + .fill(tmp[cntr, ij])
            rm1[cntr] <- d
            m1lup <- rbind(m1lup, data.frame(i = i, j = j, k = NA, 
              l = NA, m = NA, C = C, idx = cntr, which = "L1"))
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
                tmp[cntr, skip + cntr] <- 1 + .fill(tmp[cntr, skip + cntr])
                tmp[cntr, kj] <- (d - 2) + .fill(tmp[cntr, kj])
                tmp[cntr, n_d + ikC] <- -1 + .fill(tmp[cntr, n_d + ikC])
                rm1[cntr] <- d - 1
                m1lup <- rbind(m1lup, data.frame(i = i, j = j, k = k, 
                  l = NA, m = NA, C = C, idx = cntr, which = "L2"))
                cntr <- cntr + 1
              }
              ### L3
              ik <- xij$idx[xij$i == i & xij$j == k]
              jk <- xij$idx[xij$i == j & xij$j == k]
              kC <- n_d + n_z + dic$idx[dic$i == k & dic$C == C]
              tmp[cntr, skip + cntr] <- 1 + .fill(tmp[cntr, skip + cntr])
              tmp[cntr, ik] <- (d - 2) + .fill(tmp[cntr, ik])
              tmp[cntr, jk] <- (d - 2) + .fill(tmp[cntr, jk])
              tmp[cntr, kC] <- -(d - 2) + .fill(tmp[cntr, kC])
              rm1[cntr] <- 2 * (d - 2) + 2
              m1lup <- rbind(m1lup, data.frame(i = i, j = j, k = k, 
                l = NA, m = NA, C = C, idx = cntr, which = "L3"))
              cntr <- cntr + 1
            }
          }
          for (l in setdiff(seq_len(d), c(i, j, k))) {
            ### L4
            for (C in seq_len(n_C)) {
              if (all(!c(i, j, l) %in% CC[[C]])) {
                kC <- n_d + n_z + dic$idx[dic$i == k & dic$C == C]
                ik <- xij$idx[xij$i == i & xij$j == k]
                lk <- xij$idx[xij$i == l & xij$j == k]
                jlC <- max(zijc$idx[zijc$i == j & zijc$j == l & zijc$C == C],
                          zijc$idx[zijc$i == l & zijc$j == j & zijc$C == C])
                tmp[cntr, skip + cntr] <- 1 + .fill(tmp[cntr, skip + cntr])
                tmp[cntr, n_d + jlC] <- -1 + .fill(tmp[cntr, n_d + jlC])
                tmp[cntr, ik] <- (d - 3) + .fill(tmp[cntr, ik])
                tmp[cntr, lk] <- (d - 3) + .fill(tmp[cntr, lk])
                tmp[cntr, kC] <- -(d - 3) + .fill(tmp[cntr, kC])
                rm1[cntr] <- 2 * (d - 3) + 2
                m1lup <- rbind(m1lup, data.frame(i = i, j = j, k = k, 
                  l = l, m = NA, C = C, idx = cntr, which = "L4"))
                cntr <- cntr + 1
              }
            }
            for (m in setdiff(seq_len(d), c(i, j, k, l))) {
              ### L5
              for (C in seq_len(n_C)) {
                if (all(!c(i, j, l, m) %in% CC[[C]])) {
                  lk <- xij$idx[xij$i == l & xij$j == k]
                  mk <- xij$idx[xij$i == m & xij$j == k]
                  kC <- n_d + n_z + dic$idx[dic$i == k & dic$C == C]
                  ilC <- max(zijc$idx[zijc$i == i & zijc$j == l & zijc$C == C],
                            zijc$idx[zijc$i == l & zijc$j == i & zijc$C == C])
                  mjC <- max(zijc$idx[zijc$i == m & zijc$j == j & zijc$C == C],
                            zijc$idx[zijc$i == j & zijc$j == m & zijc$C == C])
                  tmp[cntr, skip + cntr] <- 1 + .fill(tmp[cntr, skip + cntr])
                  tmp[cntr, n_d + ilC] <- -1 + .fill(tmp[cntr, n_d + ilC])
                  tmp[cntr, n_d + mjC] <- -1 + .fill(tmp[cntr, n_d + mjC])
                  tmp[cntr, lk] <- (d - 4) + .fill(tmp[cntr, lk])
                  tmp[cntr, mk] <- (d - 4) + .fill(tmp[cntr, mk])
                  tmp[cntr, kC] <- -(d - 4) + .fill(tmp[cntr, kC])
                  rm1[cntr] <- 2 * (d - 4) + 2
                  m1lup <- rbind(m1lup, data.frame(i = i, j = j, k = k, 
                    l = l, m = m, C = C, idx = cntr, which = "L5"))
                  cntr <- cntr + 1
                }
              }
            }
          }
        }
      }
    }

    nM1 <- cntr - 1

    A[grep("m1min", rownames(A)), .multigrep(c("dxij", "lijc", "diC", "uijklm"), colnames(A))] <- tmp
    model$rhs[grep("m1min", names(model$rhs))] <- rm1

    ### Indicators dij-> (C2) + (C3)

    A[grep("indic", rownames(A)), grep("leij|deij", colnames(A))] <-
    rbind(
      cbind(diag(n_d), -diag(n_d)), # (C2)
      cbind(-diag(n_d), (d - 1) * diag(n_d)) # (C3)
    )

    ### DAG1
    if (verbose) {
      cat("\nWorking on constraint DAG1")
    }

    tmp <- A[grep("acyc", rownames(A)), grep("deij", colnames(A))]
    cntr <- 1
    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        ij <- xij$idx[xij$i == i & xij$j == j]
        ji <- xij$idx[xij$i == j & xij$j == i]
        tmp[cntr, ij] <- 1 + .fill(tmp[cntr, ij])
        tmp[cntr, ji] <- 1 + .fill(tmp[cntr, ji])
        model$rhs[grep("acyc", names(model$rhs))[cntr]] <- 1
        cntr <- cntr + 1
      }
    }
    A[grep("acyc", rownames(A)), grep("deij", colnames(A))] <- tmp

    ### Auxiliary variables for min constraint N1 (D1), (D2)

    if (verbose) {
      cat("\nWorking on constraints D1-D2")
    }

    tmp <- A[grep("minN1", rownames(A)), grep("dxij|leij|nijk", colnames(A))]
    rN1 <- rep(0, sum(nN1))
    cntr <- 1
    skip <- 2 * n_d
    tab_N1 <- data.frame()
    for (i in seq_len(d)) {
      for (j in setdiff(seq_len(d), i)) {
        ij <- xij$idx[xij$i == i & xij$j == j]
        tmp[cntr, skip + cntr] <- 1 + .fill(tmp[cntr, skip + cntr])
        tmp[cntr, ij] <- d - 1 + .fill(tmp[cntr, ij])
        rN1[cntr] <- d
        tab_N1 <- rbind(tab_N1, data.frame(i = i, j = j, idx = cntr, name = "D1"))
        cntr <- cntr + 1
        for (k in setdiff(seq_len(d), c(i, j))) {
          ik <- xij$idx[xij$i == i & xij$j == k]
          kj <- xij$idx[xij$i == k & xij$j == j]
          tmp[cntr, skip + cntr] <- 1 + .fill(tmp[cntr, skip + cntr])
          tmp[cntr, kj] <- d - 2 + .fill(tmp[cntr, kj])
          tmp[cntr, n_d + ik] <- -1 + .fill(tmp[cntr, n_d + ik])
          rN1[cntr] <- d - 1
          tab_N1 <- rbind(tab_N1, data.frame(i = i, j = j, idx = cntr, name = "D2"))
          cntr <- cntr + 1
        }
      }
    }

    model$rhs[grep("d1d2min", names(model$rhs))] <- rN1
    A[grep("minN1", rownames(A)), grep("dxij|leij|nijk", colnames(A))] <- tmp

    ### CONSTRAINTS FOR LINEARIZING THE OBJECTIVE

    cmat_t <- rbind(
      cbind(diag(n_z), -diag(n_z)),
      cbind(diag(n_z), diag(n_z))
    )

    A[
      which(stringr::str_detect(rownames(A), "labs")),
      which(stringr::str_detect(colnames(A), "tijC|zijC"))
    ] <- cmat_t

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
              resvar = 3 * n_z + 3 * n_d + sum(nN1) + sum(nlc) + ic,
              vars = c(sort(3 * n_z + 2 * n_d + sum(nN1) + keep), length(model$obj))
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
      saveRDS(R1, fn[5])
    }
  }

  ### INCLUDE CONSTRAINTS
  model$A <- A
  model$genconmin <- c(N1, M1, R1)

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
    diC <- x[3 * n_z + 3 * n_d + sum(nN1) + sum(nlc) + 1:ndic]
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
      antlen = leij, minlen = lijC, pind = deij, diC = diC
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

