# Original admg cache builder, retained as the reference for compiled encoding.
.admg_lean_build_r <- function(d, max_size, mode, CC, n_C, n_d, n_z, n_b, n_lb, n_tilde, nkc, nN1, nr1b, ndic, dic, zijc, xij, bxij, lbijc, model, verbose) {
  A <- MatrixExtra::emptySparse(
    nrow = length(model$rhs), ncol = length(model$obj), format = "C"
  )

  colnames(A) <- make.unique(rep(
    c("tijC", "zijC", "dxij", "lijc", "leij", "nijk", "deij", "uijklm", "bxij", "sxij", "lbijc", "zbijc", "fijc", "diC", "constd", "const1"),
    c(n_z, n_z, n_d, n_z, n_d, sum(nN1), n_d, sum(nkc), n_b, n_d, n_lb, n_lb, n_lb + n_lb * (d - 2), ndic, 1, 1)
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
    rep("redundant", d * (d - 1) / 2)
  ))

  if (max_size > 0) {
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
  }

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
  cn <- .multigrep(c("dxij", "lijc", "sxij", "lbijc", "zbijc", "diC", "uijklm"), colnames(A), value = TRUE)
  skip <- n_d + n_z + n_d + n_lb + n_lb + ndic
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
              ik <- n_d + n_z + xij$idx[xij$i == i & xij$j == k]
              jk <- n_d + n_z + xij$idx[xij$i == j & xij$j == k]
              clist[[cntr]] <- data.frame(
                i = c(cntr, cntr, cntr),
                j = c(skip + cntr, ik, jk),
                v = c(1, d - 2, d - 2)
              )
              rhs_m1[cntr] <- 2 * (d - 1)
              m1lup[[cntr]] <- data.frame(i = i, j = j, k = k,
                l = NA, m = NA, C = C, idx = cntr, which = "K3a")
              cntr <- cntr + 1
              for (k2 in setdiff(seq_len(d), c(i, j, k))) {
                if (k2 %in% CC[[C]]) {
                  ### K3b
                  ik <- n_d + n_z + xij$idx[xij$i == i & xij$j == k]
                  jk2 <- n_d + n_z + xij$idx[xij$i == j & xij$j == k2]
                  kk2c <- n_d + n_z + n_d
                  kk2c <- kk2c + max(
                    lbijc$idx[lbijc$i == k & lbijc$j == k2 & lbijc$C == C],
                    lbijc$idx[lbijc$i == k2 & lbijc$j == k & lbijc$C == C]
                  )
                  clist[[cntr]] <- data.frame(
                    i = c(cntr, cntr, cntr, cntr, cntr),
                    j = c(skip + cntr, ik, jk2, kk2c, n_lb + kk2c),
                    v = c(1, d - 3, d - 3, -1, d - 3)
                  )
                  rhs_m1[cntr] <- 3 * d - 7
                  m1lup[[cntr]] <- data.frame(i = i, j = j, k = k,
                    l = NA, m = NA, C = C, idx = cntr, which = "K3b")
                  cntr <- cntr + 1
                }
              }
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
    .multigrep(c("dxij", "lijc", "sxij", "lbijc", "zbijc", "diC", "uijklm"), colnames(A))] <-
    Matrix::sparseMatrix(i = clist$i, j = clist$j, x = clist$v,
      dims = c(length(rn), length(cn)), dimnames = list(rn, cn))
  model$rhs[grep("p1min", names(model$rhs))] <- rhs_m1

  ### Indicators dij-> (C2) + (C3)

  A[grep("indic", rownames(A)), .multigrep(c("leij", "deij"), colnames(A))] <-
  rbind(
    cbind(diag(n_d), -diag(n_d)), # (C2)
    cbind(-diag(n_d), (d - 1) * diag(n_d)) # (C3)
  )

  ### Swap constraints for ADMGs
  if (mode == "admg-lean") {

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

  list(A = A, rhs = model$rhs, P1 = P1, N1 = N1, O1 = O1, R1 = R1)
}
