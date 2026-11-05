# Construct the cache-miss encoding. C emits the constraint triplets and the
# Gurobi minimum constraints; Matrix assembles the compressed sparse matrix.
.dcon_lean_build_c <- function(d, max_size, mode, CC, dic, zijc, xij,
                               n_z, n_d, nN1, nlc, nr1b, model, verbose) {
  if (verbose) cat("\nBuilding DAG/DG constraint encoding")
  n_C <- length(CC)
  encoded <- .Call("glip_dcon_build", as.integer(d), CC, zijc, xij,
                   as.integer(sum(nlc)), as.integer(nr1b), model$rhs,
                   mode == "dag", PACKAGE = "glip")
  entries <- seq_len(encoded[[8]])
  col_names <- make.unique(rep(
    c("tijC", "zijC", "dxij", "lijc", "leij", "nijk", "deij", "uijklm", "diC", "const1"),
    c(n_z, n_z, n_d, n_z, n_d, sum(nN1), n_d, sum(nlc), d * n_C, 1)
  ))
  row_names <- make.unique(c(
    rep("labs", 2 * n_z), rep("minN1", sum(nN1)),
    rep("acyc", 2 * n_d), rep("indic", 2 * n_d),
    rep("m1min", sum(nlc)), rep("ZLcons", 2 * n_z),
    rep("r1b", nr1b), rep("redundant", d * (d - 1) / 2)
  ))
  A <- Matrix::sparseMatrix(i = encoded[[1]][entries], j = encoded[[2]][entries],
                            x = encoded[[3]][entries],
                            dims = c(length(model$rhs), length(model$obj)),
                            dimnames = list(row_names, col_names))
  list(A = A, rhs = encoded[[4]], M1 = encoded[[5]], N1 = encoded[[6]], R1 = encoded[[7]])
}
