#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <limits.h>
#include <string.h>

/* Emit the L1--L4 chain or K1--K4 ADMG auxiliary block in the order of
   the original R loops. Column indices are relative to their selected block. */
SEXP glip_lean_path_block(SEXP d_, SEXP cc, SEXP zijc, SEXP xij,
                          SEXP nlb_, SEXP cap_, SEXP admg_) {
    int d=asInteger(d_), nc=LENGTH(cc), nlb=asInteger(nlb_);
    int cap=asInteger(cap_), admg=asLogical(admg_);
    int nz=LENGTH(VECTOR_ELT(zijc,0));
    double nd_d=(double)d*(d-1);
    if(d<2 || nlb<0 || cap<0 || (double)cap*5>INT_MAX ||
       2*nd_d+nz+2.0*nlb+(double)d*nc+cap>INT_MAX ||
       nd_d>INT_MAX || (double)d*d*nc>INT_MAX)
        error("lean path encoding is too large");
    int nd=(int)nd_d;
    int *member=(int*)R_alloc((size_t)d*nc,sizeof(int));
    int *z=(int*)R_alloc((size_t)d*d*nc,sizeof(int));
    int *x=(int*)R_alloc((size_t)d*d,sizeof(int));
    memset(member,0,(size_t)d*nc*sizeof(int));
    memset(z,0,(size_t)d*d*nc*sizeof(int));
    memset(x,0,(size_t)d*d*sizeof(int));
    for(int c=0;c<nc;++c) {
        SEXP s=VECTOR_ELT(cc,c);
        for(int q=0;q<LENGTH(s);++q) member[(INTEGER(s)[q]-1)*nc+c]=1;
    }
    const int *zi=INTEGER(VECTOR_ELT(zijc,0)), *zj=INTEGER(VECTOR_ELT(zijc,1)), *zc=INTEGER(VECTOR_ELT(zijc,2));
    for(int q=0;q<nz;++q) {
        int i=zi[q]-1,j=zj[q]-1,c=zc[q]-1;
        z[(i*d+j)*nc+c]=z[(j*d+i)*nc+c]=q+1;
    }
    const int *xi=INTEGER(VECTOR_ELT(xij,0)), *xj=INTEGER(VECTOR_ELT(xij,1));
    for(int q=0;q<nd;++q) x[(xi[q]-1)*d+xj[q]-1]=q+1;
    SEXP out=PROTECT(allocVector(VECSXP,8));
    SET_VECTOR_ELT(out,0,allocVector(INTSXP,5*(R_xlen_t)cap));
    SET_VECTOR_ELT(out,1,allocVector(INTSXP,5*(R_xlen_t)cap));
    SET_VECTOR_ELT(out,2,allocVector(REALSXP,5*(R_xlen_t)cap));
    for(int q=3;q<7;++q) SET_VECTOR_ELT(out,q,allocVector(INTSXP,cap));
    SET_VECTOR_ELT(out,7,allocVector(INTSXP,2));
    int *ri=INTEGER(VECTOR_ELT(out,0)),*cj=INTEGER(VECTOR_ELT(out,1));
    double *val=REAL(VECTOR_ELT(out,2));
    int *li=INTEGER(VECTOR_ELT(out,3)),*lj=INTEGER(VECTOR_ELT(out,4));
    int *lc=INTEGER(VECTOR_ELT(out,5)),*rhs=INTEGER(VECTOR_ELT(out,6));
    int skip=2*nd+nz+(admg ? 2*nlb : 0)+d*nc;
    int e=0,n=0;
#define ADD(col,value) do { ri[e]=n; cj[e]=(col); val[e]=(value); ++e; } while(0)
#define ROW(i,j,c,r) do { if(n>=cap) error("lean path block exceeds allocation"); ++n; li[n-1]=(i)+1; lj[n-1]=(j)+1; lc[n-1]=(c)+1; rhs[n-1]=(r); ADD(skip+n,1); } while(0)
    for(int i=0;i<d;++i) {
        R_CheckUserInterrupt();
        for(int j=0;j<d;++j) if(i!=j) {
            for(int c=0;c<nc;++c) if(!member[i*nc+c] && !member[j*nc+c]) {
                ROW(i,j,c,d); ADD((admg ? nd+nz : 0)+x[i*d+j],d-1);
            }
            for(int k=0;k<d;++k) if(k!=i && k!=j) for(int c=0;c<nc;++c)
                if(!member[i*nc+c] && !member[j*nc+c]) {
                    if(!member[k*nc+c]) {
                        ROW(i,j,c,d-1); ADD(x[k*d+j],d-2); ADD(nd+z[(i*d+k)*nc+c],-1);
                        ROW(i,j,c,0); ADD(nd+z[(i*d+k)*nc+c],-1);
                        ADD(nd+z[(k*d+j)*nc+c],-1);
                        ADD(2*nd+nz+(admg ? 2*nlb : 0)+k*nc+c+1,2-d);
                    } else {
                        ROW(i,j,c,admg ? 2*(d-1) : 2*(d-1)+1);
                        ADD(nd+nz+x[i*d+k],admg ? d-2 : d-1);
                        ADD(nd+nz+x[j*d+k],admg ? d-2 : d-1);
                        if(admg) for(int k2=0;k2<d;++k2)
                            if(k2!=i && k2!=j && k2!=k && member[k2*nc+c]) {
                                int hi=k>k2 ? k:k2, lo=k>k2 ? k2:k;
                                int pair=(hi*(hi-1)/2+lo)*nc+c+1;
                                int lb=2*nd+nz+pair;
                                ROW(i,j,c,3*d-7);
                                ADD(nd+nz+x[i*d+k],d-3); ADD(nd+nz+x[j*d+k2],d-3);
                                ADD(lb,-1); ADD(lb+nlb,d-3);
                            }
                    }
                }
        }
    }
#undef ROW
#undef ADD
    INTEGER(VECTOR_ELT(out,7))[0]=n;
    INTEGER(VECTOR_ELT(out,7))[1]=e;
    UNPROTECT(1);
    return out;
}
