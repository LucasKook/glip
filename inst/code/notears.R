### Code for testing R interface to NOTEARS via Python/cdt
### LK 2025

devtools::load_all()
library("reticulate")

# CLI: conda create --name glip python=3.12
use_condaenv("glip", required = TRUE)

# conda_install("glip", "cdt", pip = TRUE)
# cdt <- import("cdt", convert = TRUE)

# conda_install("glip", "torch", pip = TRUE)
# conda_install("glip", "dagma", pip = TRUE)
# conda_install("glip", "numpy", pip = TRUE)
# conda_install("glip", "CausalDisco", pip = TRUE)
utils <- import("dagma.utils", convert = TRUE)
dagma <- import("dagma.linear", convert = TRUE)
dnl <- import("dagma.nonlinear", convert = TRUE)

### Example usage of NOTEARS (linear)
data <- r_to_py(rd <- scale(iris[, 1:4])) # Normalize data
data <- data$copy()
model <- dagma$DagmaLinear(loss_type = "l2")
output <- model$fit(data, lambda1 = 0.02)
1 * (output != 0)

### Example usage of NOTEARS (nonlinear)
# eqm <- dnl$DagmaMLP(dims = list(4L, 10L, 1L), bias = TRUE)
# model <- dnl$DagmaNonlinear(eqm)
# output <- model$fit(data, lambda1 = 0.02, lambda2 = 0.005)
# 1 * (output != 0)

cd <- import("CausalDisco.baselines", convert = TRUE)
1 * (cd$r2_sort_regress(data) != 0)
1 * (cd$var_sort_regress(data) != 0)
1 * (cd$random_sort_regress(data) != 0)

lg <- learn_graph(data.frame(rd))
lg$computed
