#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <limits.h>
#include <string.h>

static SEXP constraint(int res, const int *vars, int n) {
    SEXP out=PROTECT(allocVector(VECSXP,2));
    SEXP names=PROTECT(allocVector(STRSXP,2));
    SEXP r=PROTECT(ScalarReal(res));
    SEXP v=PROTECT(allocVector(REALSXP,n));
    for(int i=0;i<n;++i) REAL(v)[i]=vars[i];
    SET_VECTOR_ELT(out,0,r); SET_VECTOR_ELT(out,1,v);
    SET_STRING_ELT(names,0,mkChar("resvar")); SET_STRING_ELT(names,1,mkChar("vars"));
    setAttrib(out,R_NamesSymbol,names);
    UNPROTECT(4); return out;
}

/* Assemble all non-path blocks and general constraints. Path triplets were
   emitted by glip_lean_path_block and are translated to global columns here. */
SEXP glip_chain_build(SEXP d_, SEXP cc, SEXP zijc, SEXP xij, SEXP luij,
                      SEXP nl_, SEXP nr_, SEXP rhs_, SEXP path) {
    int d=asInteger(d_), nc=LENGTH(cc), nl=asInteger(nl_), nr=asInteger(nr_);
    int nz=LENGTH(VECTOR_ELT(zijc,0));
    double nd_d=(double)d*(d-1), nn_d=nd_d*(d-1), ul_d=nd_d/2;
    double nw_d=ul_d*(d-1), ny_d=nd_d*(d-1);
    double rows_d=4.0*nz+nn_d+4*nd_d+nl+nw_d+ny_d+2*ul_d+nr;
    double cols_d=3.0*nz+4*nd_d+nn_d+nl+2*ul_d+nw_d+ny_d+(double)d*nc+3;
    if(d<2 || nl<0 || nr<0 || rows_d>INT_MAX/5 || cols_d>INT_MAX)
        error("chain encoding is too large for a sparse Matrix");
    int nd=(int)nd_d, nn=(int)nn_d, nu=(int)ul_d, nw=(int)nw_d, ny=(int)ny_d;
    int rows=(int)rows_d, cols=(int)cols_d;
    int *member=(int*)R_alloc((size_t)d*nc,sizeof(int));
    int *z=(int*)R_alloc((size_t)d*d*nc,sizeof(int));
    int *x=(int*)R_alloc((size_t)d*d,sizeof(int));
    int *u=(int*)R_alloc((size_t)d*d,sizeof(int));
    memset(member,0,(size_t)d*nc*sizeof(int)); memset(z,0,(size_t)d*d*nc*sizeof(int));
    memset(x,0,(size_t)d*d*sizeof(int)); memset(u,0,(size_t)d*d*sizeof(int));
    for(int c=0;c<nc;++c) {
        SEXP s=VECTOR_ELT(cc,c);
        for(int k=0;k<LENGTH(s);++k) member[(INTEGER(s)[k]-1)*nc+c]=1;
    }
    const int *zi=INTEGER(VECTOR_ELT(zijc,0)), *zj=INTEGER(VECTOR_ELT(zijc,1)), *zc=INTEGER(VECTOR_ELT(zijc,2));
    for(int q=0;q<nz;++q) { int i=zi[q]-1,j=zj[q]-1,c=zc[q]-1; z[(i*d+j)*nc+c]=z[(j*d+i)*nc+c]=q+1; }
    const int *xi=INTEGER(VECTOR_ELT(xij,0)), *xj=INTEGER(VECTOR_ELT(xij,1));
    for(int q=0;q<nd;++q) x[(xi[q]-1)*d+xj[q]-1]=q+1;
    const int *ui=INTEGER(VECTOR_ELT(luij,0)), *uj=INTEGER(VECTOR_ELT(luij,1));
    const int *uid=INTEGER(VECTOR_ELT(luij,2));
    for(int q=0;q<2*nu;++q) u[(ui[q]-1)*d+uj[q]-1]=uid[q];
    int num_r1=0;
    for(int i=0;i<d;++i) for(int c=0;c<nc;++c) if(!member[i*nc+c]) ++num_r1;
    SEXP out=PROTECT(allocVector(VECSXP,10));
    SET_VECTOR_ELT(out,0,allocVector(INTSXP,5*(R_xlen_t)rows));
    SET_VECTOR_ELT(out,1,allocVector(INTSXP,5*(R_xlen_t)rows));
    SET_VECTOR_ELT(out,2,allocVector(REALSXP,5*(R_xlen_t)rows));
    SET_VECTOR_ELT(out,3,duplicate(rhs_));
    SET_VECTOR_ELT(out,4,allocVector(VECSXP,nz));
    SET_VECTOR_ELT(out,5,allocVector(VECSXP,nd));
    SET_VECTOR_ELT(out,6,allocVector(VECSXP,nu));
    SET_VECTOR_ELT(out,7,allocVector(VECSXP,nd));
    SET_VECTOR_ELT(out,8,allocVector(VECSXP,num_r1));
    SET_VECTOR_ELT(out,9,allocVector(INTSXP,1));
    int *ri=INTEGER(VECTOR_ELT(out,0)), *cj=INTEGER(VECTOR_ELT(out,1));
    double *val=REAL(VECTOR_ELT(out,2)), *rhs=REAL(VECTOR_ELT(out,3));
    int *tmp=(int*)R_alloc((size_t)d+1,sizeof(int));
    int ct=0,cz=nz,cx=2*nz,cl=2*nz+nd,cle=3*nz+nd,cn=3*nz+2*nd;
    int cde=cn+nn,cu=cde+nd,cs=cu+nl,czu=cs+nd,clu=czu+nu;
    int cw=clu+nu,cy=cw+nw,cdi=cy+ny;
    int rn=2*nz,ra=rn+nn, rind=ra+2*nd,rm=rind+2*nd,rzl=rm+nl;
    int rw=rzl+2*nz,ry=rw+nw,ruz=ry+ny,rr=ruz+2*nu;
    int e=0;
#define ADD(r,c,v) do { double av=(v); if(av!=0) { ri[e]=(r)+1; cj[e]=(c)+1; val[e]=av; ++e; } } while(0)
    for(int q=0;q<nz;++q) {
        ADD(q,ct+q,1); ADD(q,cz+q,-1); ADD(nz+q,ct+q,1); ADD(nz+q,cz+q,1);
        ADD(rzl+q,cz+q,1); ADD(rzl+q,cl+q,1);
        ADD(rzl+nz+q,cz+q,1-d); ADD(rzl+nz+q,cl+q,-1);
    }
    for(int q=0;q<nd;++q) {
        ADD(rind+q,cle+q,1); ADD(rind+q,cde+q,-1);
        ADD(rind+nd+q,cle+q,-1); ADD(rind+nd+q,cde+q,d-1);
    }
    for(int q=0;q<nu;++q) {
        ADD(ruz+q,czu+q,1); ADD(ruz+q,clu+q,1);
        ADD(ruz+nu+q,czu+q,1-d); ADD(ruz+nu+q,clu+q,-1);
    }
    int row=0,pair=0;
    for(int i=0;i<d;++i) for(int j=0;j<d;++j) if(i!=j) {
        int ij=x[i*d+j],ji=x[j*d+i],start=row, ar=2*pair;
        ADD(ra+ar,cle+ij-1,1); ADD(ra+ar,cle+ji-1,-1);
        ADD(ra+ar,cde+ij-1,d-1); ADD(ra+ar,cde+ji-1,d-1);
        ADD(ra+ar+1,cle+ij-1,-1); ADD(ra+ar+1,cle+ji-1,1);
        ADD(ra+ar+1,cde+ij-1,d-1); ADD(ra+ar+1,cde+ji-1,d-1);
        ADD(rn+row,cn+row,1); ADD(rn+row,cx+ij-1,d-1); rhs[rn+row]=d; ++row;
        for(int k=0;k<d;++k) if(k!=i && k!=j) {
            ADD(rn+row,cn+row,1); ADD(rn+row,cx+x[k*d+j]-1,d-2);
            ADD(rn+row,cle+x[i*d+k]-1,-1); rhs[rn+row]=d-1; ++row;
        }
        for(int q=0;q<row-start;++q) tmp[q]=cn+start+q+1;
        SET_VECTOR_ELT(VECTOR_ELT(out,5),pair,constraint(cle+ij,tmp,row-start));
        ++pair;
    }
    row=0; pair=0;
    for(int i=0;i<d;++i) for(int j=0;j<i;++j) {
        int start=row,ij=x[i*d+j],ji=x[j*d+i];
        ADD(rw+row,cw+row,1); ADD(rw+row,cx+ij-1,d-1);
        ADD(rw+row,cx+ji-1,d-1); rhs[rw+row]=2*d-1; ++row;
        for(int k=0;k<d;++k) if(k!=i && k!=j) {
            ADD(rw+row,cw+row,1); ADD(rw+row,clu+u[i*d+k]-1,-1);
            ADD(rw+row,cx+x[j*d+k]-1,d-2); ADD(rw+row,cx+x[k*d+j]-1,d-2);
            rhs[rw+row]=2*d-3; ++row;
        }
        int count=0; for(int q=start;q<row;++q) tmp[count++]=cw+q+1;
        tmp[count++]=cols;
        SET_VECTOR_ELT(VECTOR_ELT(out,6),pair++,constraint(clu+u[i*d+j],tmp,count));
    }
    row=0; pair=0;
    for(int i=0;i<d;++i) for(int j=0;j<d;++j) if(i!=j) {
        int start=row;
        ADD(ry+row,cy+row,1); ADD(ry+row,cx+x[i*d+j]-1,-1);
        ADD(ry+row,cx+x[j*d+i]-1,1); ++row;
        for(int k=0;k<d;++k) if(k!=i && k!=j) {
            ADD(ry+row,cy+row,1); ADD(ry+row,cx+x[i*d+k]-1,-1);
            ADD(ry+row,cx+x[k*d+i]-1,1); ADD(ry+row,czu+u[k*d+j]-1,-1);
            rhs[ry+row]=-1; ++row;
        }
        int count=0; for(int q=start;q<row;++q) tmp[count++]=cy+q+1;
        tmp[count++]=cols-2;
        SET_VECTOR_ELT(VECTOR_ELT(out,7),pair++,constraint(cs+x[i*d+j],tmp,count));
    }
    row=0;
    for(int i=0;i<d;++i) for(int c=0;c<nc;++c) if(member[i*nc+c]) { ADD(rr+row,cdi+i*nc+c,1); ++row; }
    row=0;
    for(int i=0;i<d;++i) for(int c=0;c<nc;++c) if(!member[i*nc+c]) {
        int count=0;
        for(int k=0;k<d;++k) if(member[k*nc+c]) tmp[count++]=cde+x[i*d+k];
        tmp[count++]=cols-1;
        SET_VECTOR_ELT(VECTOR_ELT(out,8),row++,constraint(cdi+i*nc+c+1,tmp,count));
    }
    int m=INTEGER(VECTOR_ELT(path,7))[0], ne=INTEGER(VECTOR_ELT(path,7))[1];
    int *pr=INTEGER(VECTOR_ELT(path,0)), *pc=INTEGER(VECTOR_ELT(path,1));
    double *pv=REAL(VECTOR_ELT(path,2));
    for(int q=0;q<ne;++q) {
        int c=pc[q]-1, global;
        if(c<nd) global=cx+c;
        else if(c<nd+nz) global=cl+c-nd;
        else if(c<2*nd+nz) global=cs+c-nd-nz;
        else if(c<2*nd+nz+d*nc) global=cdi+c-2*nd-nz;
        else global=cu+c-2*nd-nz-d*nc;
        ADD(rm+pr[q]-1,global,pv[q]);
    }
    int *group=(int*)R_alloc((size_t)m,sizeof(int));
    int *counts=(int*)R_alloc((size_t)nz+1,sizeof(int));
    int *starts=(int*)R_alloc((size_t)nz+1,sizeof(int));
    int *pos=(int*)R_alloc((size_t)nz+1,sizeof(int));
    memset(counts,0,((size_t)nz+1)*sizeof(int));
    int *pi=INTEGER(VECTOR_ELT(path,3)),*pj=INTEGER(VECTOR_ELT(path,4));
    int *ppc=INTEGER(VECTOR_ELT(path,5)),*pv_rhs=INTEGER(VECTOR_ELT(path,6));
    for(int q=0;q<m;++q) {
        int g=z[((pi[q]-1)*d+pj[q]-1)*nc+ppc[q]-1];
        group[q]=g; ++counts[g]; rhs[rm+q]=pv_rhs[q];
    }
    starts[0]=0; for(int q=1;q<=nz;++q) starts[q]=starts[q-1]+counts[q-1];
    memcpy(pos,starts,((size_t)nz+1)*sizeof(int));
    int *group_rows=(int*)R_alloc((size_t)m+1,sizeof(int));
    for(int q=0;q<m;++q) group_rows[pos[group[q]]++]=cu+q+1;
    for(int q=1;q<=nz;++q) SET_VECTOR_ELT(VECTOR_ELT(out,4),q-1,
        constraint(cl+q,group_rows+starts[q],counts[q]));
#undef ADD
    INTEGER(VECTOR_ELT(out,9))[0]=e;
    UNPROTECT(1); return out;
}
