# Original index enumeration, used as an independent reference in tests and
# benchmarks. Its ordering is part of the on-disk cache format.
.dcon_lean_indices_r <- function(d, CC) {
  # Array indices for path indicators
  dic <- data.frame()
  counter <- 1
  for (i in seq_len(d)) {
    for (C in seq_along(CC)) {
      dic <- rbind(
        dic, data.frame(i = i, C = C, iinc = i %in% CC[[C]], idx = counter)
      )
      counter <- counter + 1
    }
  }

  # Array indices for d-separation statements
  zijc <- data.frame()
  counter <- 1
  for (i in seq_len(d)) {
    for (j in seq_len(i - 1)) {
      for (C in seq_along(CC)) {
        if (all(!(c(i, j) %in% CC[[C]]))) {
          zijc <- rbind(
            zijc,
            data.frame(i = i, j = j, C = C, idx = counter)
          )
          counter <- counter + 1
        }
      }
    }
  }
  list(dic = dic, zijc = zijc)
}

.dcon_lean_indices_c <- function(d, CC) {
  tab <- .Call("glip_dcon_indices", as.integer(d), CC, PACKAGE = "glip")
  nc <- length(CC)
  dic <- data.frame(i = rep(seq_len(d), each = nc),
                    C = rep(seq_len(nc), times = d),
                    iinc = tab[[4]], idx = as.double(seq_len(d * nc)))
  zijc <- data.frame(i = tab[[1]], j = tab[[2]], C = tab[[3]],
                     idx = as.double(seq_along(tab[[1]])))
  list(dic = dic, zijc = zijc)
}

# All unordered pairs and conditioning sets in the ADMG cache order.
.lean_pair_conditioning_indices <- function(d, n_C) {
  i <- rep(seq_len(d), times = seq_len(d) - 1L)
  j <- sequence(seq_len(d) - 1L)
  n <- length(i) * n_C
  data.frame(i = rep(i, each = n_C), j = rep(j, each = n_C),
             C = rep(seq_len(n_C), times = length(i)),
             idx = as.double(seq_len(n)))
}
