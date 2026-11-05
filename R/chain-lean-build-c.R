# Construct the chain cache-miss encoding in C, then assemble sparse triplets.
.chain_lean_build_c <- function(d, max_size, CC, n_C, n_d, n_z, n_ul,
                                n_tilde, nlc, nN1, nW1, nZ1, nr1b,
                                ndic, dic, zijc, xij, luij, model, verbose) {
  if (verbose) cat("\nBuilding chain constraint encoding")
  path <- .Call("glip_lean_path_block", as.integer(d), CC, zijc, xij,
                0L, as.integer(sum(nlc)), FALSE, PACKAGE = "glip")
  encoded <- .Call("glip_chain_build", as.integer(d), CC, zijc, xij, luij,
                   as.integer(sum(nlc)), as.integer(nr1b), model$rhs, path,
                   PACKAGE = "glip")
  used <- seq_len(encoded[[10]])
  col_names <- make.unique(rep(
    c("tijC", "zijC", "dxij", "lijc", "leij", "nijk", "deij", "uijklm",
      "sxij", "zuij", "luij", "wijk", "zijk", "diC", "const0", "const1", "constd"),
    c(n_z, n_z, n_d, n_z, n_d, sum(nN1), n_d, sum(nlc), n_d, n_ul,
      n_ul, nW1, nZ1, d * n_C, 1, 1, 1)))
  row_names <- make.unique(c(rep("labs", 2 * n_z), rep("minN1", sum(nN1)),
    rep("acyc", 2 * n_d), rep("indic", 2 * n_d), rep("m1min", sum(nlc)),
    rep("ZLcons", 2 * n_z), rep("w1min", nW1), rep("z1max", nZ1),
    rep("unZL", 2 * n_ul), rep("r1b", nr1b)))
  A <- Matrix::sparseMatrix(i = encoded[[1]][used], j = encoded[[2]][used],
                           x = encoded[[3]][used],
                           dims = c(length(model$rhs), length(model$obj)),
                           dimnames = list(row_names, col_names))
  list(A = A, rhs = encoded[[4]], M1 = encoded[[5]], N1 = encoded[[6]],
       W1 = encoded[[7]], Z1 = encoded[[8]], R1 = encoded[[9]])
}
