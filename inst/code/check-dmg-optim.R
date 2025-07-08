### Check dmg optim
### LK 2025

devtools::load_all()
library("pcalg")

d <- 5
max_size <- 3
V <- letters[1:d]

tries <- 1:100
glog <- list()

tmp <- lapply(tries, \(idx) {
  cat("\n", idx)
  set.seed(idx)

  G <- create_dmg(d, nodeNames = V, diag = FALSE, M2prob = 0.3)
  G$M1[upper.tri(G$M1)] <- 0 # Make ADMG

  if (any(sapply(glog, \(x) isTRUE(all.equal(x, G))))) {
    # cat("\nSkipping already tested graph")
    return(NULL)
  }
  glog <<- c(glog, list(G))

  D <- createD0(G)
  gD <- suppressWarnings(as(D$M1, "graphNEL"))

  sets <- .list_tests_graph(V, max_size = max_size)$sets
  tests <- lapply(sets, \(x) {
    data.frame(
      X = x$X,
      Y = x$Y,
      Z = paste0(x$Z, collapse = ","),
      p.value = 1 * pcalg::dsep(x$X, x$Y, x$Z, g = gD)
    )
  }) |> do.call("rbind", args = _)

  tmp <- capture.output(
    lG <- admg_optim(tests,
      d = d,
      max_size = max_size,
      V = V,
      verbose = TRUE,
      cache = TRUE,
      gurobi_args = list(
        # TimeLimit = 10,
        Threads = 8
      ),
      mode = "admg"
    )
  )

  lP <- compute_pag(lG$graph)
  P <- compute_pag(G)

  if (!isTRUE(all.equal(lP, P))) {
    cat("\nWrong output found. Writing...")

    tmp <- capture.output({
      cat("\nLearned ADMG:\n")
      print(lG$graph)

      cat("\nGround truth:\n")
      print(G)

      cat("\nLearned PAG:\n")
      print(lP)

      cat("\nGround truth PAG:\n")
      print(P)
    })

    writeLines(tmp, file.path(".wrong-output", paste0(
      idx, "-", d, "-",
      max_size, ".out"
    )))
  }
})
