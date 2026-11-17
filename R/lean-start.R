### Complete warm starts for the min-length encodings
###
### The min-length encodings fix only the edge variables in model$start; all
### other variables (separation and path indicators, path lengths, the
### auxiliaries of the min/max general constraints and the absolute values of
### the objective) are determined by the edges. Gurobi completes a partial
### start with a sub-MIP that can fail for large models. .lean_full_start()
### computes every variable from the model itself:
###   - resultants of general min/max constraints: min/max of their inputs;
###   - for each equality row with exactly one variable that is neither fixed
###     by the start nor a resultant: that variable, solved from the row;
###   - every other variable: the cheapest value (by its objective
###     coefficient) within the bounds implied by the rows in which it is the
###     only such variable; variables in no row take their lower bound.
### Resultants start at their upper bounds and the three rules are applied
### until no value changes. All rules are monotone, so this converges to the
### greatest fixed point (the shortest path lengths, capped at the bounds).
### The full start is returned only if it satisfies every bound, integrality
### requirement, row and general constraint (.lean_check_start()); otherwise
### NULL, and the caller keeps the partial start.
.lean_full_start <- function(model, max_iter = 1000L, tol = 1e-9, why = FALSE) {
  fail <- function(reason) if (why) reason else NULL
  start <- model$start
  if (is.null(start) || !any(!is.na(start))) {
    return(NULL)
  }
  A <- methods::as(methods::as(model$A, "CsparseMatrix"), "generalMatrix")
  A <- methods::as(A, "dMatrix")
  dimnames(A) <- NULL # rows appended by .bk_rows() leave partial dimnames
  n <- ncol(A)
  lb <- model$lb
  ub <- model$ub
  obj <- model$obj
  if (is.null(lb)) lb <- rep(0, n)
  if (is.null(ub)) ub <- rep(Inf, n)
  int <- if (is.null(model$vtype)) rep(FALSE, n) else model$vtype %in% c("B", "I")
  ub[model$vtype %in% "B"] <- pmin(ub[model$vtype %in% "B"], 1)
  sense <- model$sense
  rhs <- model$rhs
  if (length(sense) == 1) sense <- rep(sense, nrow(A))

  known <- !is.na(start)
  gens <- list(
    min = if (is.null(model$genconmin)) list() else model$genconmin,
    max = if (is.null(model$genconmax)) list() else model$genconmax
  )
  gres <- lapply(gens, \(g) vapply(g, \(x) as.integer(x$resvar), 1L))
  gvars <- lapply(gens, \(g) lapply(g, \(x) as.integer(x$vars)))
  gcon <- lapply(names(gens), \(k) vapply(gens[[k]], \(x) {
    if (is.null(x$con)) (if (k == "min") Inf else -Inf) else x$con
  }, 1))
  names(gcon) <- names(gens)
  resvar <- rep(FALSE, n)
  resvar[unlist(gres)] <- TRUE
  if (any(duplicated(unlist(gres))) || any(known & resvar)) {
    return(fail("resultant fixed by the start or of two general constraints"))
  }

  ### Triplets of A
  Tr <- Matrix::summary(A)
  Tr <- data.frame(i = Tr$i, j = Tr$j, x = Tr$x)
  Tr <- Tr[Tr$x != 0, ]
  eqrow <- sense == "="

  ### Evaluation order: layers of variables defined by an equality row whose
  ### other variables are available ("eq"), else by the bounds from the
  ### other rows whose other variables are available ("iv"); start values
  ### and resultants are available from the beginning
  ### If neither applies (variables linked only among themselves, e.g. the
  ### unused bidirected lengths of the ADMG encodings for max_size <= 1),
  ### the first remaining variable keeps its start value (upper bound if
  ### finite, else lower bound) and the layering continues
  avail <- known | resvar
  layers <- list()
  repeat {
    tn <- Tr[!avail[Tr$j], ]
    if (!nrow(tn)) break
    cnt <- tabulate(tn$i, nrow(A))
    e <- tn[eqrow[tn$i] & cnt[tn$i] == 1, ]
    e <- e[!duplicated(e$j), ]
    type <- "eq"
    if (!nrow(e)) {
      e <- tn[!eqrow[tn$i] & cnt[tn$i] == 1, ]
      type <- "iv"
    }
    if (!nrow(e)) {
      ### lowest-index remaining variable of every row (in one pass)
      first <- tn[order(tn$i, tn$j), ]
      avail[unique(first$j[!duplicated(first$i)])] <- TRUE
      next
    }
    rows <- sort(unique(e$i))
    e$r <- match(e$i, rows)
    e$v <- match(e$j, unique(e$j))
    layers[[length(layers) + 1L]] <- list(type = type, tr = e, A = A[rows, , drop = FALSE],
      vars = unique(e$j))
    avail[e$j] <- TRUE
  }
  free <- !avail
  if (any(free & !is.finite(lb))) {
    return(fail("free variable without a finite lower bound"))
  }

  v <- start
  v[!known] <- ifelse(is.finite(ub[!known]), ub[!known], lb[!known])
  v[free] <- lb[free]

  gmin_idx <- if (length(gres$min)) rep(seq_along(gres$min), lengths(gvars$min)) else integer(0)
  gmax_idx <- if (length(gres$max)) rep(seq_along(gres$max), lengths(gvars$max)) else integer(0)
  gmin_v <- unlist(gvars$min)
  gmax_v <- unlist(gvars$max)
  ### Grouped minimum of x by group g in 1..k (Inf for empty groups)
  gmin_f <- function(x, g, k) {
    r <- rep(Inf, k)
    if (length(x)) {
      o <- order(g, x)
      first <- !duplicated(g[o])
      r[g[o][first]] <- x[o][first]
    }
    r
  }

  for (it in seq_len(max_iter)) {
    old <- v
    for (L in layers) {
      tr <- L$tr
      act <- as.numeric(L$A %*% v)[tr$r]
      rest <- act - tr$x * v[tr$j]
      bnd <- (rhs[tr$i] - rest) / tr$x
      if (L$type == "eq") {
        v[tr$j] <- bnd
        next
      }
      s <- sense[tr$i]
      ### a x <= b: x <= b/a if a > 0, x >= b/a if a < 0
      is_lb <- (s == ">=" & tr$x > 0) | (s == "<=" & tr$x < 0)
      iv <- L$vars
      k <- length(iv)
      lo <- pmax(lb[iv], -gmin_f(-bnd[is_lb], tr$v[is_lb], k))
      up <- pmin(ub[iv], gmin_f(bnd[!is_lb], tr$v[!is_lb], k))
      lo[int[iv]] <- ceiling(lo[int[iv]] - 1e-6)
      up[int[iv]] <- floor(up[int[iv]] + 1e-6)
      pick <- ifelse(obj[iv] < 0, up, lo)
      alt <- ifelse(obj[iv] < 0, lo, up)
      pick[!is.finite(pick)] <- alt[!is.finite(pick)]
      v[iv] <- pick
    }
    ### Resultants of the general constraints
    if (length(gmin_v)) {
      v[gres$min] <- pmin(gmin_f(v[gmin_v], gmin_idx, length(gres$min)), gcon$min)
    }
    if (length(gmax_v)) {
      v[gres$max] <- pmax(-gmin_f(-v[gmax_v], gmax_idx, length(gres$max)), gcon$max)
    }
    if (isTRUE(all(abs(v - old) <= tol))) {
      break
    }
  }
  if (it == max_iter) {
    return(fail("no fixed point"))
  }
  chk <- .lean_check_start(model, v, A = A)
  if (!chk$ok) {
    return(fail(paste("violation:", paste(names(chk$violation), signif(chk$violation, 3), sep = "=", collapse = ", "))))
  }
  v
}

### Check a complete start: bounds, integrality, rows, general constraints
.lean_check_start <- function(model, v, A = NULL, tol = 1e-6) {
  if (is.null(A)) A <- model$A
  n <- length(v)
  lb <- if (is.null(model$lb)) rep(0, n) else model$lb
  ub <- if (is.null(model$ub)) rep(Inf, n) else model$ub
  vt <- if (is.null(model$vtype)) rep("C", n) else model$vtype
  ub[vt == "B"] <- pmin(ub[vt == "B"], 1)
  viol <- c(
    bounds = max(0, lb - v, v - ub),
    integer = max(0, abs(v - round(v))[vt %in% c("B", "I")]),
    rows = {
      act <- as.numeric(A %*% v)
      s <- model$sense
      d <- ifelse(s == "<=", act - model$rhs, ifelse(s == ">=", model$rhs - act, abs(act - model$rhs)))
      max(0, d)
    },
    genconmin = max(0, vapply(model$genconmin, \(g) {
      abs(v[g$resvar] - min(v[g$vars], if (is.null(g$con)) Inf else g$con))
    }, 1)),
    genconmax = max(0, vapply(model$genconmax, \(g) {
      abs(v[g$resvar] - max(v[g$vars], if (is.null(g$con)) -Inf else g$con))
    }, 1))
  )
  list(ok = all(viol <= tol), violation = viol)
}

### Replace the partial start of a min-length model by the complete start if
### it exists (see .lean_full_start()); otherwise keep the partial start
### (option glip.lean_full_start = FALSE: edge values only, as before)
.lean_complete_start <- function(model, verbose = FALSE) {
  if (!isTRUE(getOption("glip.lean_full_start", TRUE)) ||
    is.null(model$start) || !any(!is.na(model$start))) {
    return(model)
  }
  full <- .lean_full_start(model)
  if (is.null(full)) {
    if (verbose) {
      cat("\nNo complete warm start (the start violates a constraint);",
        "passing the edge values only.\n")
    }
  } else {
    model$start <- full
  }
  model
}
