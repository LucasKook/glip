### Experiment A (paper Sec. 5.1/6.1, Fig. 7; App. D.1): timing against the
### ASP method of Hyttinen et al. (2014). One task = one data set (mode, d,
### ms, rep). Uses the tracked inst/asp/hyttinen2014uai_ver6.zip (R code,
### ASP encodings, Linux clingo 4.3.0 as in the paper) with the pipeline of
### inst/asp-glip/pipeline.R, which calls GLIP_FIT() (below) for both
### encodings on the same data, tests, warm start and edge hints. Each task
### runs in its own copy of the package (the code writes to ../tmp).
### GLIP: Threads 1, TimeLimit walltime, BestObjStop 1e-4 (as in the paper);
### clingo: --configuration=crafty --time-limit=walltime (single thread).

rr_asp_home <- function() {
  home <- Sys.getenv("GLIP_ASP_HOME")
  if (nzchar(home)) {
    return(home)
  }
  home <- file.path(tempdir(), "asp-home")
  if (!dir.exists(file.path(home, "hyttinen2014uai_ver6"))) {
    dir.create(home, showWarnings = FALSE)
    utils::unzip(file.path(rr_root(), "inst", "asp", "hyttinen2014uai_ver6.zip"),
      exdir = home
    )
  }
  home
}

rr_task_asp <- function(task, out) {
  file <- file.path(out, task$id, "result.rds")
  if (file.exists(file)) {
    return(readRDS(file))
  }
  ### Private copy of the package
  src <- file.path(rr_asp_home(), "hyttinen2014uai_ver6", "pkg")
  work <- tempfile("asp-")
  dir.create(file.path(work, "pkg", "tmp"), recursive = TRUE)
  file.copy(file.path(src, c("R", "ASP")), file.path(work, "pkg"), recursive = TRUE)
  file.copy(file.path(rr_root(), "inst", "asp-glip", "pipeline.R"),
    file.path(work, "pkg", "R", "pipeline.R"),
    overwrite = TRUE
  )
  clingo <- file.path(work, "pkg", "ASP", "clingo430")
  if (nzchar(Sys.getenv("GLIP_CLINGO"))) {
    file.remove(clingo)
    file.symlink(normalizePath(Sys.getenv("GLIP_CLINGO")), clingo)
  }
  Sys.chmod(clingo, "755")
  if (is.null(.rr$clingo)) {
    .rr$clingo <- tryCatch(
      system2(clingo, "--version", stdout = TRUE)[1],
      error = function(e) NA
    )
  }
  owd <- setwd(file.path(work, "pkg", "R"))
  on.exit({
    setwd(owd)
    unlink(work, recursive = TRUE)
  })
  env <- new.env(parent = globalenv())
  sys.source("load.R", envir = env)
  suppressPackageStartupMessages(
    evalq({
      library("deal")
      library("combinat")
      library("graph")
      sourceDir("./", trace = FALSE, local = environment())
      sourceDir("./mods/", trace = FALSE, local = environment())
      others <- FALSE
    }, env)
  )
  ### Globals used by the pipeline
  assign("mode", task$mode, envir = env)
  assign("GLIP_FIT", function(data, ws, eh, schedule, p, md, true_graph) {
    V <- colnames(data)
    d <- length(V)
    t0 <- Sys.time()
    tests <- learn_graph(data,
      test = "gaussCItest", comets = FALSE, max_size = schedule, mode = md,
      return_tests_only = TRUE
    )
    t_tests <- as.numeric(Sys.time() - t0, units = "secs")
    do.call(rbind, lapply(task$encodings, \(enc) {
      f <- rr_glip(enc, md, tests, d, schedule, V, p, task$walltime, ws, eh,
        extra = list(BestObjStop = 1e-4)
      )
      data.frame(encoding = enc, rr_solver_info(f), t_tests = t_tests)
    }))
  }, envir = env)
  restrict <- switch(task$mode,
    dag = c("acyclic", "sufficient"),
    admg = c("acyclic")
  )
  set.seed(rr_seed("asp", task$mode, task$d, task$ms, task$rep))
  res <- evalq(pipeline, env)(
    n = task$d, N = task$n, schedule = task$ms, test = "classic", verbose = 1,
    clingoconf = paste0(
      "--configuration=crafty --time-limit=", task$walltime, " --quiet=1,0"
    ),
    restrict = restrict, weight = "constant", p = task$alpha
  )
  res <- data.frame(
    exp = "asp", mode = task$mode, d = task$d, ms = task$ms, n = task$n,
    iter = task$rep, walltime = task$walltime, alpha = task$alpha, res,
    row.names = NULL
  )
  rr_save(res, file)
  res
}
