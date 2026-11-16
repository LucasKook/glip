#' @keywords internal
"_PACKAGE"

#' @importFrom methods as new
#' @importFrom stats as.formula cor rbinom runif
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom comets comets
#' @importFrom pcalg gaussCItest
NULL

### Column names used in dplyr verbs (lookup_ci() and the lean encodings)
utils::globalVariables(c("X", "Y", "Z", "p.value", "i", "j"))
