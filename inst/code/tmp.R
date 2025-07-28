devtools::load_all()

g <- cbind(0, rbind(diag(4), 0))
V <- letters[1:5]
dimnames(g) <- list(V, V)
g

marginalize_dag_to_admg(g, c("a", "b", "e"))
