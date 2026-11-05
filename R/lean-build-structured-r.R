# Structurally optimized R encoder for benchmarking and validation. It uses
# direct integer lookup arrays and assembles preallocated sparse rows once.
# The original R reference remains separate; production continues to use C.
.lean_build_structured_r <- function(kind, d, max_size, mode, CC, zijc, xij,
                                     n_z, n_d, nN1, nlc, nr1b, model,
                                     pair_edges = NULL, n_pair = 0L,
                                     nW1 = 0L, nZ1 = 0L, n_lb = 0L) {
  nc <- length(CC); nd <- n_d; nz <- n_z; nn <- sum(nN1); nl <- sum(nlc)
  member <- matrix(FALSE, d, nc)
  for (cc in seq_len(nc)) member[CC[[cc]], cc] <- TRUE
  z <- array(0L, c(d, d, nc))
  z[cbind(zijc$i, zijc$j, zijc$C)] <- seq_len(nz)
  z[cbind(zijc$j, zijc$i, zijc$C)] <- seq_len(nz)
  x <- matrix(0L, d, d)
  x[cbind(xij$i, xij$j)] <- xij$idx
  if (!is.null(pair_edges)) {
    b <- matrix(0L, d, d)
    b[cbind(pair_edges$i, pair_edges$j)] <- pair_edges$idx
  }
  nr <- length(model$rhs); ncol <- length(model$obj)
  cols <- vector("list", nr); vals <- vector("list", nr)
  rhs <- model$rhs
  ct <- 0; cz <- nz; cx <- 2*nz; cl <- cx+nd; cle <- cl+nz
  cn <- cle+nd; cde <- cn+nn; cu <- cde+nd
  rn <- 2*nz; ra <- rn+nn
  if (kind == "dag") {
    cs <- NA; cdi <- cu+nl
    rind <- ra+2*nd; rm <- rind+2*nd; rzl <- rm+nl
    rr <- rzl+2*nz; rred <- rr+nr1b
    col_names <- make.unique(rep(c("tijC","zijC","dxij","lijc","leij","nijk","deij","uijklm","diC","const1"),
                                  c(nz,nz,nd,nz,nd,nn,nd,nl,d*nc,1)))
    row_names <- make.unique(c(rep("labs",2*nz),rep("minN1",nn),rep("acyc",2*nd),
      rep("indic",2*nd),rep("m1min",nl),rep("ZLcons",2*nz),rep("r1b",nr1b),rep("redundant",nd/2)))
  } else if (kind == "chain") {
    cs <- cu+nl; czu <- cs+nd; clu <- czu+n_pair; cw <- clu+n_pair
    cy <- cw+nW1; cdi <- cy+nZ1
    rind <- ra+2*nd; rm <- rind+2*nd; rzl <- rm+nl
    rw <- rzl+2*nz; ry <- rw+nW1; ruz <- ry+nZ1; rr <- ruz+2*n_pair
    col_names <- make.unique(rep(c("tijC","zijC","dxij","lijc","leij","nijk","deij","uijklm",
      "sxij","zuij","luij","wijk","zijk","diC","const0","const1","constd"),
      c(nz,nz,nd,nz,nd,nn,nd,nl,nd,n_pair,n_pair,nW1,nZ1,d*nc,1,1,1)))
    row_names <- make.unique(c(rep("labs",2*nz),rep("minN1",nn),rep("acyc",2*nd),rep("indic",2*nd),
      rep("m1min",nl),rep("ZLcons",2*nz),rep("w1min",nW1),rep("z1max",nZ1),rep("unZL",2*n_pair),rep("r1b",nr1b)))
  } else {
    cb <- cu+nl; cs <- cb+n_pair; clb <- cs+nd; czb <- clb+n_lb
    nf <- n_lb*(d-1); cf <- czb+n_lb; cdi <- cf+nf
    rind <- ra+nd; rm <- rind+2*nd; rzl <- rm+nl
    rc6 <- rzl+2*nz; rc10 <- rc6+3*nd; rf <- rc10+2*n_lb; rr <- rf+nf; rred <- rr+nr1b
    col_names <- make.unique(rep(c("tijC","zijC","dxij","lijc","leij","nijk","deij","uijklm",
      "bxij","sxij","lbijc","zbijc","fijc","diC","constd","const1"),
      c(nz,nz,nd,nz,nd,nn,nd,nl,n_pair,nd,n_lb,n_lb,nf,d*nc,1,1)))
    row_names <- make.unique(c(rep("labs",2*nz),rep("minN1",nn),rep("acyc",nd),rep("indic",2*nd),
      rep("p1min",nl),rep("ZLcons",2*nz),rep("C6to9",3*nd),rep("C10to11",2*n_lb),
      rep("f1f2min",nf),rep("r1b",nr1b),rep("redundant",nd/2)))
  }
  for (q in seq_len(nz)) {
    cols[[q]] <- c(ct+q,cz+q); vals[[q]] <- c(1,-1)
    cols[[nz+q]] <- c(ct+q,cz+q); vals[[nz+q]] <- c(1,1)
    cols[[rzl+q]] <- c(cz+q,cl+q); vals[[rzl+q]] <- c(1,1)
    cols[[rzl+nz+q]] <- c(cz+q,cl+q); vals[[rzl+nz+q]] <- c(1-d,-1)
  }
  for (q in seq_len(nd)) {
    cols[[rind+q]] <- c(cle+q,cde+q); vals[[rind+q]] <- c(1,-1)
    cols[[rind+nd+q]] <- c(cle+q,cde+q); vals[[rind+nd+q]] <- c(-1,d-1)
  }
  if (kind == "chain") for (q in seq_len(n_pair)) {
    cols[[ruz+q]] <- c(czu+q,clu+q); vals[[ruz+q]] <- c(1,1)
    cols[[ruz+n_pair+q]] <- c(czu+q,clu+q); vals[[ruz+n_pair+q]] <- c(1-d,-1)
  }
  if (kind == "admg") for (q in seq_len(n_lb)) {
    cols[[rc10+q]] <- c(czb+q,clb+q); vals[[rc10+q]] <- c(1,1)
    cols[[rc10+n_lb+q]] <- c(czb+q,clb+q); vals[[rc10+n_lb+q]] <- c(1-d,-1)
  }
  N1 <- vector("list",nd); row <- pair <- 0L
  for (i in seq_len(d)) for (j in seq_len(d)) if(i!=j) {
    pair <- pair+1L; ij <- x[i,j]; ji <- x[j,i]; start <- row
    if (kind == "chain") {
      ar <- ra+2*(pair-1)
      cols[[ar+1]] <- cols[[ar+2]] <- c(cle+ij,cle+ji,cde+ij,cde+ji)
      vals[[ar+1]] <- c(1,-1,d-1,d-1); vals[[ar+2]] <- c(-1,1,d-1,d-1)
    } else if (mode %in% c("dag","admg-lean")) {
      cols[[ra+pair]] <- c(cde+ij,cde+ji); vals[[ra+pair]] <- c(1,1); rhs[ra+pair] <- 1
    }
    if (kind == "admg") {
      a <- rc6+3*(pair-1); bij <- b[i,j]
      cols[[a+1]] <- c(cx+ij,cb+bij,cs+ij); vals[[a+1]] <- c(-1,-1,1)
      cols[[a+2]] <- c(cx+ij,cs+ij); vals[[a+2]] <- c(1,-1)
      cols[[a+3]] <- c(cb+bij,cs+ij); vals[[a+3]] <- c(1,-1)
    }
    row <- row+1L; cols[[rn+row]] <- c(cn+row,cx+ij); vals[[rn+row]] <- c(1,d-1); rhs[rn+row] <- d
    for (k in seq_len(d)) if(k!=i && k!=j) {
      row <- row+1L; cols[[rn+row]] <- c(cn+row,cx+x[k,j],cle+x[i,k])
      vals[[rn+row]] <- c(1,d-2,-1); rhs[rn+row] <- d-1
    }
    N1[[pair]] <- list(resvar=cle+ij,vars=as.numeric(cn+seq.int(start+1L,row)))
  }
  if (kind != "chain" && mode %in% c("dag","admg-lean")) {
    row <- 0L
    for(i in seq_len(d)) for(j in seq_len(i-1L)) {
      row <- row+1L; cols[[rred+row]] <- c(cx+x[i,j],cx+x[j,i]); vals[[rred+row]] <- c(1,1); rhs[rred+row] <- 1
    }
  }
  row <- 0L; R1 <- vector("list", sum(!member)); r1 <- 0L
  for (i in seq_len(d)) for (cc in seq_len(nc)) {
    if (member[i,cc]) {
      row <- row+1L; cols[[rr+row]] <- cdi+(i-1L)*nc+cc; vals[[rr+row]] <- 1
    } else {
      r1 <- r1+1L
      keep <- which(member[,cc]); const <- if(kind=="chain") ncol-1 else ncol
      R1[[r1]] <- list(resvar=cdi+(i-1L)*nc+cc, vars=as.numeric(c(cde+x[i,keep],const)))
    }
  }
  # L1--L4 / K1--K4: each row is emitted once, and its separation
  # variable is recorded for the later grouped minimum constraints.
  owner <- integer(nl); row <- 0L
  for (i in seq_len(d)) for (j in seq_len(d)) if(i!=j) {
    for (cc in seq_len(nc)) if (!member[i,cc] && !member[j,cc]) {
      row <- row+1L; owner[row] <- z[i,j,cc]; rhs[rm+row] <- d
      cols[[rm+row]] <- c(cu+row, if(kind=="admg") cs+x[i,j] else cx+x[i,j]); vals[[rm+row]] <- c(1,d-1)
    }
    for (k in seq_len(d)) if(k!=i && k!=j) for (cc in seq_len(nc)) if (!member[i,cc] && !member[j,cc]) {
      if (!member[k,cc]) {
        row <- row+1L; owner[row] <- z[i,j,cc]; rhs[rm+row] <- d-1
        cols[[rm+row]] <- c(cu+row,cx+x[k,j],cl+z[i,k,cc]); vals[[rm+row]] <- c(1,d-2,-1)
        row <- row+1L; owner[row] <- z[i,j,cc]; rhs[rm+row] <- 0
        cols[[rm+row]] <- c(cu+row,cl+z[i,k,cc],cl+z[k,j,cc],cdi+(k-1L)*nc+cc)
        vals[[rm+row]] <- c(1,-1,-1,2-d)
      } else {
        row <- row+1L; owner[row] <- z[i,j,cc]
        rhs[rm+row] <- if(kind=="chain") 2*(d-1)+1 else 2*(d-1)
        base <- if(kind=="dag") cx else cs; coef <- if(kind=="chain") d-1 else d-2
        cols[[rm+row]] <- c(cu+row,base+x[i,k],base+x[j,k]); vals[[rm+row]] <- c(1,coef,coef)
        if(kind=="admg") for(k2 in seq_len(d)) if(k2!=i && k2!=j && k2!=k && member[k2,cc]) {
          row <- row+1L; owner[row] <- z[i,j,cc]; rhs[rm+row] <- 3*d-7
          hi <- max(k,k2); lo <- min(k,k2); lb <- (hi-1)*(hi-2)/2*nc+(lo-1)*nc+cc
          cols[[rm+row]] <- c(cu+row,cs+x[i,k],cs+x[j,k2],clb+lb,czb+lb)
          vals[[rm+row]] <- c(1,d-3,d-3,-1,d-3)
        }
      }
    }
  }
  groups <- split(seq_len(row), factor(owner[seq_len(row)],levels=seq_len(nz)))
  M1 <- lapply(seq_len(nz), function(q) list(resvar=cl+q,vars=as.numeric(cu+groups[[q]])))
  if (kind == "chain") {
    W1 <- vector("list",n_pair); row <- pair <- 0L
    for(i in seq_len(d)) for(j in seq_len(i-1L)) {
      start <- row; pair <- pair+1L; row <- row+1L
      cols[[rw+row]] <- c(cw+row,cx+x[i,j],cx+x[j,i]); vals[[rw+row]] <- c(1,d-1,d-1)
      rhs[rw+row] <- 2*d-1
      for(k in seq_len(d)) if(k!=i && k!=j) {
        row <- row+1L; cols[[rw+row]] <- c(cw+row,clu+b[i,k],cx+x[j,k],cx+x[k,j])
        vals[[rw+row]] <- c(1,-1,d-2,d-2); rhs[rw+row] <- 2*d-3
      }
      W1[[pair]] <- list(resvar=clu+b[i,j], vars=as.numeric(c(cw+seq.int(start+1L,row),ncol)))
    }
    Z1 <- vector("list",nd); row <- pair <- 0L
    for(i in seq_len(d)) for(j in seq_len(d)) if(i!=j) {
      start <- row; pair <- pair+1L; row <- row+1L
      cols[[ry+row]] <- c(cy+row,cx+x[i,j],cx+x[j,i]); vals[[ry+row]] <- c(1,-1,1)
      for(k in seq_len(d)) if(k!=i && k!=j) {
        row <- row+1L; cols[[ry+row]] <- c(cy+row,cx+x[i,k],cx+x[k,i],czu+b[k,j])
        vals[[ry+row]] <- c(1,-1,1,-1); rhs[ry+row] <- -1
      }
      Z1[[pair]] <- list(resvar=cs+x[i,j],vars=as.numeric(c(cy+seq.int(start+1L,row),ncol-2)))
    }
  }
  if (kind == "admg") {
    O1 <- list()
    if(max_size>1) {
      O1 <- vector("list",n_lb); row <- group <- 0L
      for(j in seq_len(d)) for(i in seq_len(j-1L)) for(cc in seq_len(nc)) {
        start <- row; group <- group+1L
        if(member[i,cc] && member[j,cc]) {
          row <- row+1L; cols[[rf+row]] <- c(cf+row,cb+b[i,j]); vals[[rf+row]] <- c(1,d-1); rhs[rf+row] <- d
          for(k in seq_len(d)) if(k!=i && k!=j && member[k,cc]) {
            hi <- max(i,k); lo <- min(i,k); lb <- (hi-1)*(hi-2)/2*nc+(lo-1)*nc+cc
            row <- row+1L; cols[[rf+row]] <- c(cf+row,cb+b[j,k],clb+lb)
            vals[[rf+row]] <- c(1,d-2,-1); rhs[rf+row] <- d-1
          }
        }
        keep <- if(row>start) cf+seq.int(start+1L,row) else numeric()
        O1[[group]] <- list(resvar=clb+group,vars=as.numeric(c(keep,ncol-1)))
      }
    }
  }
  count <- lengths(cols)
  ii <- rep.int(seq_len(nr),count); jj <- unlist(cols,use.names=FALSE)
  vv <- unlist(vals,use.names=FALSE); nonzero <- vv != 0
  A <- Matrix::sparseMatrix(i=ii[nonzero],j=jj[nonzero],x=vv[nonzero],
                            dims=c(nr,ncol),dimnames=list(row_names,col_names))
  if(kind=="dag") return(list(A=A,rhs=rhs,M1=M1,N1=N1,R1=R1))
  if(kind=="chain") return(list(A=A,rhs=rhs,M1=M1,N1=N1,W1=W1,Z1=Z1,R1=R1))
  list(A=A,rhs=rhs,P1=M1,N1=N1,O1=O1,R1=R1)
}

.dcon_lean_build_structured_r <- function(d, max_size, mode, CC, dic, zijc, xij,
                                         n_z, n_d, nN1, nlc, nr1b, model, verbose) {
  .lean_build_structured_r("dag", d, max_size, mode, CC, zijc, xij,
                           n_z, n_d, nN1, nlc, nr1b, model)
}
.chain_lean_build_structured_r <- function(d, max_size, CC, n_C, n_d, n_z, n_ul,
                                          n_tilde, nlc, nN1, nW1, nZ1, nr1b,
                                          ndic, dic, zijc, xij, luij, model, verbose) {
  .lean_build_structured_r("chain", d, max_size, "chain", CC, zijc, xij,
                           n_z, n_d, nN1, nlc, nr1b, model,
                           pair_edges=luij, n_pair=n_ul, nW1=nW1, nZ1=nZ1)
}
.admg_lean_build_structured_r <- function(d, max_size, mode, CC, n_C, n_d, n_z, n_b,
                                         n_lb, n_tilde, nkc, nN1, nr1b, ndic, dic,
                                         zijc, xij, bxij, lbijc, model, verbose) {
  .lean_build_structured_r("admg", d, max_size, mode, CC, zijc, xij,
                           n_z, n_d, nN1, nkc, nr1b, model,
                           pair_edges=bxij, n_pair=n_b, n_lb=n_lb)
}

.dcon_lean_indices_structured_r <- function(d, CC) {
  nc <- length(CC)
  member <- matrix(FALSE, d, nc)
  for (cc in seq_len(nc)) member[CC[[cc]], cc] <- TRUE
  dic <- data.frame(i=rep(seq_len(d),each=nc), C=rep(seq_len(nc),times=d),
                    iinc=as.vector(t(member)), idx=as.double(seq_len(d*nc)))
  pi <- rep(seq_len(d),times=seq_len(d)-1L)
  pj <- sequence(seq_len(d)-1L)
  i <- rep(pi,each=nc); j <- rep(pj,each=nc); cc <- rep(seq_len(nc),times=length(pi))
  keep <- !member[cbind(i,cc)] & !member[cbind(j,cc)]
  zijc <- data.frame(i=i[keep],j=j[keep],C=cc[keep],idx=as.double(seq_len(sum(keep))))
  list(dic=dic,zijc=zijc)
}
