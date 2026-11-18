### Example 1 (paper Sec. 1.1, Fig. 1): 4-node SCM where FCI fails for
### strongly correlated latent confounding; GLIP (ADMG) vs FCI with the
### settings of the paper (alpha .05, Fisher z test, n 300, 50 reps per
### correlation, seed 42); one Gurobi thread, default time limit; the GLIP
### encodings in GLIP_ENCODINGS (default: min-length).
### Usage (repository root): Rscript inst/code/intro-example.R [outdir]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[1] else "results/paper/intro-example"
source(file.path("inst", "code", "common.R"))
rr_setup()
set.seed(42)

dgp <- function(n, cor = 0.9) {
  U <- mvtnorm::rmvnorm(n, sigma = matrix(c(1, cor, cor, 1), nrow = 2, ncol = 2))
  X1 <- U[, 1]
  X2 <- U[, 2]
  X3 <- X1 + X2 + rnorm(n)
  X4 <- X3 + rnorm(n)
  data.frame(a = X1, b = X2, c = X3, d = X4)
}

cors <- c(0, 0.9, 0.9999)
alpha <- 0.05
nsim <- 50
n <- 3e2
encodings <- rr_encodings()
key <- function(M) paste(as.vector(M), collapse = "")

out <- lapply(cors, \(tcor) {
  lapply(seq_len(nsim), \(iter) {
    d <- dgp(n, cor = tcor)
    glip <- lapply(encodings, \(enc) {
      tmp <- utils::capture.output(
        lP <- learn_graph(d,
          alpha = alpha, comets = FALSE, test = "zf",
          trafo = \(x) 1 * (x <= alpha), mode = "admg-lean", encoding = enc,
          gurobi_args = list(Threads = 1)
        )
      )
      data.frame(
        cor = tcor, iter = iter, method = paste0("GLIP-", enc),
        output = key(lP$computed), status = lP$graph$optim$status,
        objval = lP$graph$optim$objval, runtime = lP$graph$optim$runtime
      )
    })
    fc <- pcalg::fci(list(C = stats::cor(d), n = NROW(d)), pcalg::gaussCItest, alpha,
      colnames(d),
      selectionBias = FALSE
    )
    rbind(do.call(rbind, glip), data.frame(
      cor = tcor, iter = iter, method = "FCI", output = key(fc@amat),
      status = NA, objval = NA, runtime = NA
    ))
  }) |> do.call("rbind", args = _)
}) |> do.call("rbind", args = _)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
attr(out, "versions") <- rr_versions()
saveRDS(out, file.path(outdir, "intro-example.rds"))
### Most frequent output per correlation and method
tab <- stats::aggregate(iter ~ cor + method + output, out, length)
names(tab)[4] <- "n"
tab <- do.call(rbind, lapply(split(tab, list(tab$cor, tab$method)), \(x) x[which.max(x$n), ]))
utils::write.csv(tab, file.path(outdir, "modal-output.csv"), row.names = FALSE)
print(tab[order(tab$cor, tab$method), ], row.names = FALSE)
cat("all GLIP fits optimal:", all(out$status[out$method != "FCI"] == "OPTIMAL"), "\n")
### With both encodings: agreement per replication (different graphs with
### equal objective are ties among optimal solutions)
if (all(c("minlength", "cuts") %in% encodings)) {
  g <- out[out$method != "FCI", ]
  w <- merge(g[g$method == "GLIP-minlength", ], g[g$method == "GLIP-cuts", ],
    by = c("cor", "iter"), suffixes = c(".ml", ".cuts"))
  agr <- stats::aggregate(cbind(
    same_graph = w$output.ml == w$output.cuts,
    same_objective = abs(w$objval.ml - w$objval.cuts) < 1e-6
  ) ~ cor, w, sum)
  utils::write.csv(agr, file.path(outdir, "encoding-agreement.csv"), row.names = FALSE)
  print(agr, row.names = FALSE)
}
