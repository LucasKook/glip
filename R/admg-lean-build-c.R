# Construct the ADMG/DMG cache-miss encoding in C and assemble the triplets.
.admg_lean_build_c <- function(d, max_size, mode, CC, n_C, n_d, n_z, n_b,
                               n_lb, n_tilde, nkc, nN1, nr1b, ndic, dic,
                               zijc, xij, bxij, lbijc, model, verbose) {
  if (verbose) cat("\nBuilding ADMG/DMG constraint encoding")
  path <- .Call("glip_lean_path_block", as.integer(d), CC, zijc, xij,
                as.integer(n_lb), as.integer(sum(nkc)), TRUE, PACKAGE = "glip")
  encoded <- .Call("glip_admg_build", as.integer(d), CC, zijc, xij, bxij,
                   as.integer(sum(nkc)), as.integer(nr1b), model$rhs, path,
                   mode == "admg-lean", as.integer(max_size), PACKAGE = "glip")
  used <- seq_len(encoded[[9]])
  col_names <- make.unique(rep(
    c("tijC", "zijC", "dxij", "lijc", "leij", "nijk", "deij", "uijklm",
      "bxij", "sxij", "lbijc", "zbijc", "fijc", "diC", "constd", "const1"),
    c(n_z, n_z, n_d, n_z, n_d, sum(nN1), n_d, sum(nkc), n_b, n_d,
      n_lb, n_lb, n_lb + n_lb * (d - 2), ndic, 1, 1)))
  row_names <- make.unique(c(rep("labs", 2 * n_z), rep("minN1", sum(nN1)),
    rep("acyc", n_d), rep("indic", 2 * n_d), rep("p1min", sum(nkc)),
    rep("ZLcons", 2 * n_z), rep("C6to9", 3 * d * (d - 1)),
    rep("C10to11", 2 * n_lb), rep("f1f2min", n_lb + n_lb * (d - 2)),
    rep("r1b", nr1b), rep("redundant", d * (d - 1) / 2)))
  A <- Matrix::sparseMatrix(i = encoded[[1]][used], j = encoded[[2]][used],
                           x = encoded[[3]][used],
                           dims = c(length(model$rhs), length(model$obj)),
                           dimnames = list(row_names, col_names))
  list(A = A, rhs = encoded[[4]], P1 = encoded[[5]], N1 = encoded[[6]],
       O1 = encoded[[7]], R1 = encoded[[8]])
}
