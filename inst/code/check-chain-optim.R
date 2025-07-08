### Check chain optim
### LK 2025

devtools::load_all()
library("comets")

d <- 5
max_size <- 3

V <- letters[1:d]
G <- matrix(0, nrow = d, ncol = d)
# G <- 1 * lower.tri(diag(d))
dimnames(G) <- list(V, V)
# G["a", "c"] <- 1
# G["b", "c"] <- 1
# G["d", "a"] <- 1

### FIXME: use c-separation
sets <- .list_tests_graph(V, max_size = max_size)$sets
tests <- lapply(sets, \(x) {
  data.frame(
    X = x$X,
    Y = x$Y,
    Z = paste0(x$Z, collapse = ","),
    p.value = 1 * pcalg::dsep(x$X, x$Y, x$Z, g = as(G, "graphNEL"))
  )
}) |> do.call("rbind", args = _)

# tmp <- capture.output(
lG <- chain_optim(tests,
  d = d,
  max_size = max_size,
  V = V,
  verbose = TRUE,
  cache = FALSE,
  gurobi_args = list(
    TimeLimit = 10,
    Threads = 8
  )
)

cat("\nGround truth graph:\n")
G
cat("\nLearned graph:\n")
lG$graph

# lG$tests

### Compute essential graph
cat("\nGround truth essential graph:\n")
1 * pcalg::dag2essgraph(G)
cat("\nLearned essential graph:\n")
1 * pcalg::dag2essgraph(lG$graph)
