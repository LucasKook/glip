test_that('optimized indices and cache preserve the reference encoding', {
 ns <- asNamespace('glip')
 for (d in 2:5) for (max_size in 0:min(2,d)) for (mode in c('dag','dg')) {
 CC <- unlist(lapply(0:max_size, function(k) combn(d,k,simplify=FALSE)),recursive=FALSE)
 expect_identical(get('.dcon_lean_indices_structured_r', ns)(d, CC),
                  get('.dcon_lean_indices_r', ns)(d, CC))
 nc <- length(CC)
 dic <- do.call(rbind,lapply(seq_len(d),function(i) data.frame(i=i,C=seq_len(nc),iinc=vapply(CC,function(C)i %in% C,logical(1)),idx=(i-1)*nc+seq_len(nc))))
 zijc <- do.call(rbind,unlist(lapply(seq_len(d),function(i) lapply(seq_len(i-1),function(j) { C <- which(vapply(CC,function(z) !any(c(i,j) %in% z),logical(1))); data.frame(i=i,j=j,C=C) })),recursive=FALSE))
 zijc$idx <- seq_len(nrow(zijc))
 xij <- subset(expand.grid(i=seq_len(d),j=seq_len(d)),i!=j); xij$idx <- seq_len(nrow(xij))
 nz <- nrow(zijc); nd <- d*(d-1); nn <- c(nd,(d-2)*nd); nl <- c(nd*nc,rep(nd*(d-2)*nc,3)); nr <- if(max_size>0) d*sum(choose(d,0:(max_size-1))) else 1
 model <- list(obj=numeric(3*nz+3*nd+sum(nn)+sum(nl)+d*nc+1),rhs=c(labs=rep(0,2*nz),d1d2min=rep(0,sum(nn)),acyc=rep(0,2*nd),indic=rep(0,2*nd),m1min=rep(0,sum(nl)),ZLcons=rep(0,2*nz),r1b=rep(0,nr),redundant=rep(0,d*(d-1)/2)))
 model$rhs[seq_len(2 * nz)] <- seq_len(2 * nz) / 100
 args <- list(d,max_size,mode,CC,dic,zijc,xij,nz,nd,nn,nl,nr,model,FALSE)
 a <- do.call(get('.dcon_lean_build_r',ns),args); b <- do.call(get('.dcon_lean_build_structured_r',ns),args)
 expect_equal(b, a, info = paste(d, max_size, mode))
 structured <- do.call(get('.dcon_lean_build_structured_r', ns), args)
 expect_equal(structured, a, tolerance = 0, info = paste(d, max_size, mode))
 }
})

test_that("DAG and DG solver paths build compatible caches and reuse them", {
  skip_if_no_gurobi()
  ns <- asNamespace("glip")
  set.seed(91)
  d <- 4
  V <- letters[seq_len(d)]
  for (mode in c("dag", "dg")) for (max_size in 0:2) {
    graph <- get(".generate_random_graph", ns)(d, V = V, mode = mode, prob = 0.4)
    tests <- get(".compute_oracle_tests", ns)(graph, max_size = max_size, mode = mode)
    paths <- file.path(tempdir(), paste0("glip-encoding-", mode, "-", max_size, c("-c", "-r", "-structured")))
    unlink(paths, recursive = TRUE)
    fit <- function(path) {
      out <- capture.output(value <- dcon_lean_optim(
        tests, d = d, max_size = max_size, V = V, mode = mode,
        cache = TRUE, cache_dir = path, gurobi_args = list(OutputFlag = 0, Threads = 1)))
      value
    }
    fast <- fit(paths[1])
    reference <- testthat::with_mocked_bindings(
      fit(paths[2]),
      .dcon_lean_build_structured_r = get(".dcon_lean_build_r", ns), .package = "glip")
    structured <- fit(paths[3])
    suffix <- paste0("dim", d, "max", max_size, mode)
    for (name in c("lhs", "rhs", "M1", "N1", "R1")) {
      ext <- if (name == "lhs") ".mtx" else ".rds"
      file <- paste0(name, suffix, ext)
      read <- if (name == "lhs") Matrix::readMM else readRDS
      expect_identical(read(file.path(paths[1], file)), read(file.path(paths[2], file)),
                       info = paste(mode, max_size, name))
      expect_identical(read(file.path(paths[1], file)), read(file.path(paths[3], file)),
                       info = paste("structured", mode, max_size, name))
    }
    cached_fast <- fit(paths[1])
    cached_old <- fit(paths[2])
    expect_equal(fast$optim$objval, reference$optim$objval)
    expect_equal(fast$optim$objval, structured$optim$objval)
    expect_equal(fit(paths[3])$optim$objval, structured$optim$objval)
    expect_equal(cached_fast$optim$objval, fast$optim$objval)
    expect_equal(cached_old$optim$objval, reference$optim$objval)
    unlink(paths, recursive = TRUE)
  }
})
