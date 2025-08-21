### Example where oracle algorithms output wrong
### but global optimization approach works
### LK 2025

set.seed(42)

### DEPs
devtools::load_all()
library("tidyverse")
library("comets")
library("pcalg")

### Read data and groundtruth
dat <- read_csv("./inst/data/csuite/train.csv", col_names = FALSE)
gt <- as.matrix(read_csv("./inst/data/csuite/adj_matrix.csv", col_names = FALSE))
dimnames(gt) <- list(V <- colnames(dat), colnames(dat))

### Marginalize
which <- 5:9
gt <- marginalize_dag_to_admg(gt, V[which])
dat <- dat[, which]
V <- colnames(dat)
d <- ncol(dat)

###Params
wt <- 60
alpha <- 0.01
ms <- 2

### Run FCI
out <- fci(list(C = cor(dat), n = NROW(dat)), 
  gaussCItest, alpha, V, selectionBias = FALSE
)
FCI <- out@amat

### Run GLIP
lP <- learn_graph(dat,
  trafo = \(x) 1 * (x <= alpha),
  max_size = ms,
  test_args = list(
    reg_YonZ = "lrm",
    reg_XonZ = "lrm"
  ), mode = "admg",
  verbose = TRUE,
  warmstart = .pag_to_admg(FCI),
  edgehints = 1 * (FCI != 0),
  cache = TRUE,
  gurobi_args = list(TimeLimit = wt),
  weight = "const"
)

ORACLE <- .compute_graphical_representation(gt, d - 2, "admg")
GLIP <- .compute_graphical_representation(lP$graph$graph, d - 2, "admg")

ORACLE
GLIP
FCI

1 - sep(GLIP, ORACLE, mode = "mag", max_size = d - 2)$acc
1 - sep(FCI, ORACLE, mode = "mag", max_size = d - 2)$acc

shd(GLIP, ORACLE)
shd(FCI, ORACLE)
