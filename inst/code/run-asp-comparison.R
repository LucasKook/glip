### Timing comparison with Hyttinen et al 2014 ASP
### LK 2025

set.seed(12)

### Dependencies
devtools::load_all()
setwd("./inst/asp/hyttinen2014uai_ver6/pkg/R")
source("./load.R")
loud()
library("tidyverse")

### Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
mode <- darg(args[1], "dag")
d <- as.numeric(darg(args[2], 6))
ms <- as.numeric(darg(args[3], -1))
ms <- ifelse(ms == -1, d - 2, ms)
N <- as.numeric(darg(args[4], 3e2))
nsim <- as.numeric(darg(args[5], 1))
use_oracle_tests <- as.numeric(darg(args[6], 0))
sim_name <- darg(args[7], "test-run")
walltime <- as.numeric(darg(args[6], 30))
WTYPE <- darg(args[8], "log")
ncores <- max(7, parallel::detectCores(logical = TRUE) - 2)
GARGS <- list(Threads = ncores, TimeLimit = walltime, BestObjStop = 1e-4)
clstr <- paste0("--configuration=crafty --time-limit=", walltime, " --quiet=1,0")
save <- TRUE

### Output directory
wdir <- file.path("../../../../results/asp-comparison", Sys.Date(), sim_name)
if (!dir.exists(wdir)) {
  dir.create(wdir, recursive = TRUE)
}

### Prepare arguments
test <- c("oracle", "classic")[2 - use_oracle_tests]
n <- d
restrict <- "acyclic"
if (mode == "dag") {
  restrict <- c(restrict, "sufficient")
}

### Run
pb <- txtProgressBar(0, nsim, style = 3, width = 60)
out <- lapply(1:nsim, \(iter) {
  setTxtProgressBar(pb, iter)
  tmp <- pipeline(
    n = d, N = N, test = test, verbose = 0,
    clingoconf = clstr, restrict = restrict
  )
  res <- data.frame(
    mode = mode, d = d, n = N, ms = ms, nsim = nsim,
    use_oracle_tests = use_oracle_tests, iter = iter, tmp, row.names = NULL
  )
  if (save) {
    saveRDS(res, file.path(wdir, paste0(
      mode, "-iter_", iter, "-d_", d, "-ms_", ms, "-nsim_", nsim,
      "-n_", N, "-uot_", use_oracle_tests, "-wtype_", WTYPE,
      "-", test, ".rds"
    )))
  }
  res
}) |> do.call("rbind", args = _)
out
