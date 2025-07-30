### Timing comparison with Hyttinen et al 2014 ASP
### LK 2025

# set.seed(1)

devtools::load_all()
odir <- "../../../../figures"
wdir <- "../../../../results/asp-comparison"
setwd("./inst/asp/hyttinen2014uai_ver6/pkg/R")
source("./load.R")
loud()
library("tidyverse")
save <- TRUE

if (!dir.exists(odir)) {
  dir.create(odir, recursive = TRUE)
}
if (!dir.exists(wdir)) {
  dir.create(wdir, recursive = TRUE)
}

ncores <- max(7, parallel::detectCores(logical = TRUE) - 2)

mode <- c("dag", "admg")[1]
test <- c("oracle", "classic")[1]
nsim <- 50
ds <- 3:8
N <- 1e3
n <- 0

restrict <- "acyclic"
if (mode == "dag") {
  restrict <- c(restrict, "sufficient")
}

out <- lapply(ds, \(d) {
  cat("\nRunning d =", d, "\n")
  pb <- txtProgressBar(0, nsim, style = 3, width = 60)
  lapply(1:nsim, \(iter) {
    setTxtProgressBar(pb, iter)
    n <<- d
    res <- pipeline(
      n = d, N = N, test = test, verbose = 0,
      clingoconf = "--configuration=crafty --time-limit=500 --quiet=1,0",
      restrict = restrict
    )
    res <- data.frame(d = d, n = N, iter = iter, res, row.names = NULL)
    if (save) {
      saveRDS(res, file.path(wdir, paste0(mode, "-iter_", iter, "-d_", d, "-", test, ".rds")))
    }
    res
  }) |> do.call("rbind", args = _)
}) |> do.call("rbind", args = _)
out
