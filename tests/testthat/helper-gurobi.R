### Skip tests that need the solver if the gurobi package or a Gurobi license
### is not available (e.g., on CI without license secrets)
skip_if_no_gurobi <- function() {
  if (is.null(.gurobi_ok$ok)) {
    .gurobi_ok$ok <- requireNamespace("gurobi", quietly = TRUE) &&
      tryCatch(
        {
          m <- list(
            A = matrix(1, 1, 1), obj = 1, rhs = 1, sense = ">=",
            vtype = "B"
          )
          !is.null(gurobi::gurobi(m, list(OutputFlag = 0))$x)
        },
        error = function(e) FALSE
      )
  }
  if (!.gurobi_ok$ok) {
    testthat::skip("gurobi or a Gurobi license is not available")
  }
  invisible(TRUE)
}
.gurobi_ok <- new.env()
