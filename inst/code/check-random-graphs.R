### Small simulation / bug finder
### LK 2025

set.seed(1)
args <- commandArgs(trailingOnly = TRUE)

### Dependencies
devtools::load_all()
library("pcalg")

### Params
d <- as.numeric(darg(args[1], 3))
pp <- as.numeric(darg(args[2], 0.5))
mode <- darg(args[3], "admg")
seeds <- eval(parse(text = darg(args[4], "1:100")))
max_size <- as.numeric(darg(args[5], -1))
max_size <- ifelse(max_size == -1, d - 2, max_size)
walltime <- as.numeric(darg(args[6], 1))
V <- letters[1:d]
cache <- TRUE
ncores <- max(7, parallel::detectCores(logical = TRUE) - 2)

### Write
out <- ".wrong-output"
if (!dir.exists(out)) {
  dir.create(out)
}

### Run
glog <- list()
tmp <- sapply(seeds, \(idx) {
  cat("\n", idx)
  set.seed(idx)

  G <- .generate_random_graph(d, V, mode, prob = pp)

  if (any(sapply(glog, \(x) isTRUE(all.equal(x, G))))) {
    return(NULL)
  }
  glog <<- c(glog, list(G))

  tests <- .compute_oracle_tests(G, max_size, mode)
  capt <- capture.output(
    lG <<- .get_opt(mode)(tests,
      d = d, max_size = max_size,
      V = V, cache = cache,
      gurobi_args = list(
        Threads = ncores,
        TimeLimit = walltime
      ), mode = mode, verbose = TRUE
    )
  )

  learned <<- .compute_graphical_representation(lG$graph, max_size, mode)
  ground_truth <<- .compute_graphical_representation(G, max_size, mode)

  ### Compute output graph
  if (!isTRUE(all.equal(learned, ground_truth))) {
    message("\nWrong graph found, writing...")
    tmp <- capture.output({
      # print(capt)
      print(lG$tests)
      cat("\nTrue graph:\n")
      print(G)
      cat("\nLearned graph:\n")
      print(lG$graph)
      cat("\nTrue computed representation:\n")
      print(ground_truth)
      cat("\nLearned computed representation:\n")
      print(learned)
    })
    writeLines(tmp, file.path(out, paste0(
      mode, "-", idx, "-", d,
      "-", max_size, "-", pp, ".out"
    )))
  }
})
