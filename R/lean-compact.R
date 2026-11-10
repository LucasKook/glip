### Compact lean encodings (2026-11-10).
###
### The lean encodings size several auxiliary blocks by upper bounds (e.g. the
### M1 auxiliary variables of DAGs and chain graphs, d (d - 1) n_C
### (1 + 3 (d - 2)); P1 and F1/F2 of ADMGs and DMGs; the R1b rows; the
### acyclicity and redundant rows in cyclic modes). The builders fill every
### such block from its start, so only a prefix of each block is used (at full
### conditioning and d = 7 the model uses 14 to 18% of its rows and columns).
### The builders keep working in the bounded layout; right after building,
### .lean_compact() cuts each of these blocks to its used prefix, checks that
### the cut part is unused (no matrix entries, no general constraint, rows
### trivially satisfied), remaps the general constraints and returns the exact
### block sizes, which are cached with the RHS (attribute "sizes"). The model
### is then set up with the exact sizes; nothing is expanded back.
###
### Exact counting (2026-11-10): .lean_exact_sizes() gives the used size of
### every bounded block in closed form, so the structured builder works in
### the exact layout from the start (no bounded intermediate); it checks its
### final counters against these sizes. The reference builders keep the
### bounded layout (they obtain it from attr(model, "bound_model")) and their
### output is cut by .lean_compact(); tests compare both.

# A, rhs, sense: built model in the bounded layout; gens: named list of
# general constraint lists (resvar, vars); cols, rows: named block lengths in
# model order; cut_cols, cut_rows: names of the blocks to cut to their prefix.
.lean_compact <- function(A, rhs, sense, gens, cols, cut_cols, rows, cut_rows) {
  A <- methods::as(A, "CsparseMatrix")
  stopifnot(sum(cols) == ncol(A), sum(rows) == nrow(A), length(rhs) == nrow(A),
            length(sense) == nrow(A), all(c(cut_cols %in% names(cols), cut_rows %in% names(rows))))
  used_c <- diff(A@p) > 0
  for (g in unlist(gens, recursive = FALSE)) used_c[c(g$resvar, g$vars)] <- TRUE
  used_r <- tabulate(A@i + 1L, nrow(A)) > 0
  trivial <- ifelse(sense %in% c(">=", ">"), rhs <= 0,
                    ifelse(sense %in% c("<=", "<"), rhs >= 0, rhs == 0))
  cut <- function(used, len, which, kind) {
    keep <- rep(TRUE, length(used))
    sizes <- integer(0)
    start <- cumsum(c(0, utils::head(len, -1)))
    names(start) <- names(len)
    for (b in which) {
      idx <- start[[b]] + seq_len(len[[b]])
      k <- sum(used[idx])
      if (!all(used[idx][seq_len(k)])) {
        stop(".lean_compact(): the used ", kind, " of block ", b, " are not a prefix.")
      }
      keep[idx[seq_along(idx) > k]] <- FALSE # (idx[-seq_len(0)] would be empty)
      sizes[[b]] <- k
    }
    list(keep = keep, sizes = sizes)
  }
  cc <- cut(used_c, cols, cut_cols, "columns")
  cr <- cut(used_r, rows, cut_rows, "rows")
  if (any(!cr$keep & !trivial)) stop(".lean_compact(): a cut row is not trivially satisfied.")
  cmap <- cumsum(cc$keep)
  gens <- lapply(gens, \(gl) lapply(gl, \(g) {
    g$resvar <- cmap[g$resvar]
    g$vars <- cmap[g$vars]
    g
  }))
  rhs <- rhs[cr$keep]
  attr(rhs, "sizes") <- c(cc$sizes, cr$sizes)
  list(A = A[cr$keep, cc$keep, drop = FALSE], rhs = rhs, gens = gens, sizes = c(cc$sizes, cr$sizes))
}

# Block sizes of a cached compact model (stops for caches of the bounded layout).
.lean_cached_sizes <- function(rhs, file) {
  sz <- attr(rhs, "sizes")
  if (is.null(sz)) {
    stop("Cache ", file, " has the former (uncompacted) layout; delete it or use another cache_dir.")
  }
  sz
}

# Exact sizes of the bounded blocks (closed form; checked against the
# builders' counters). For an ordered pair (i, j) and a conditioning set C of
# size s with i, j not in C, the M1 (P1) rows are one L1 row, for each of the
# d - 2 - s nodes k outside C an L2 and an L4 row, and for each of the s nodes
# k in C one L3 row (ADMG/DMG: plus s - 1 further rows). F1/F2 (ADMG/DMG,
# max_size > 1): for an unordered pair inside C, one row plus one per other
# node of C. R1b: one row per node and set containing it.
.lean_exact_sizes <- function(kind, mode, d, max_size) {
  s <- 0:max_size
  n <- choose(d, s)
  nd <- d * (d - 1)
  m1 <- sum(n * (d - s) * (d - s - 1) *
    (1 + 2 * (d - 2 - s) + s + if (kind == "admg") s * (s - 1) else 0))
  r1b <- sum(n * s)
  acyclic <- mode %in% c("dag", "admg-lean")
  acyc <- if (kind == "chain") 2 * nd else if (acyclic) nd else 0
  red <- if (kind != "chain" && acyclic) nd / 2 else 0
  switch(kind,
    dcon = c(M1 = m1, acyc = acyc, m1min = m1, r1b = r1b, redundant = red),
    chain = c(M1 = m1, acyc = acyc, m1min = m1, r1b = r1b),
    admg = {
      f12 <- if (max_size > 1) sum(n * choose(s, 2) * (s - 1)) else 0
      c(P1 = m1, F12 = f12, acyc = acyc, p1min = m1, f1f2min = f12, r1b = r1b, redundant = red)
    })
}

# Whether a usable compact cache exists: all files present and the RHS has
# the block sizes. Caches of the former (bounded) layout are rebuilt and
# overwritten.
.lean_cache_ok <- function(cache, fn) {
  if (!cache || !all(file.exists(fn))) return(FALSE)
  if (is.null(attr(readRDS(fn[2]), "sizes"))) {
    message("Cache ", fn[2], " has the former (uncompacted) layout; rebuilding it.")
    return(FALSE)
  }
  TRUE
}
