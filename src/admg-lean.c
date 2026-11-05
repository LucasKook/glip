#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <limits.h>
#include <string.h>

static SEXP constraint(int res, const int *vars, int n) {
    SEXP out=PROTECT(allocVector(VECSXP,2)), names=PROTECT(allocVector(STRSXP,2));
    SEXP r=PROTECT(ScalarReal(res)), v=PROTECT(allocVector(REALSXP,n));
    for(int i=0;i<n;++i) REAL(v)[i]=vars[i];
    SET_VECTOR_ELT(out,0,r); SET_VECTOR_ELT(out,1,v);
    SET_STRING_ELT(names,0,mkChar("resvar")); SET_STRING_ELT(names,1,mkChar("vars"));
    setAttrib(out,R_NamesSymbol,names);
    UNPROTECT(4); return out;
}

SEXP glip_admg_build(SEXP d_, SEXP cc, SEXP zijc, SEXP xij, SEXP bxij,
                     SEXP nl_, SEXP nr_, SEXP rhs_, SEXP path,
                     SEXP dag_, SEXP max_size_) {
    int d=asInteger(d_), nc=LENGTH(cc), nl=asInteger(nl_), nr=asInteger(nr_);
    int dag=asLogical(dag_), max_size=asInteger(max_size_);
    int nz=LENGTH(VECTOR_ELT(zijc,0));
    double nd_d=(double)d*(d-1), nn_d=nd_d*(d-1), nb_d=nd_d/2, nlb_d=nb_d*nc;
    double nf_d=nlb_d*(d-1);
    double rows_d=4.0*nz+nn_d+6*nd_d+nl+2*nlb_d+nf_d+nr+nb_d;
    double cols_d=3.0*nz+4*nd_d+nn_d+nl+nb_d+2*nlb_d+nf_d+(double)d*nc+2;
    if(d<2 || nl<0 || nr<0 || rows_d>INT_MAX/5 || cols_d>INT_MAX)
        error("ADMG encoding is too large for a sparse Matrix");
    int nd=(int)nd_d, nn=(int)nn_d, nb=(int)nb_d, nlb=(int)nlb_d, nf=(int)nf_d;
    int rows=(int)rows_d, cols=(int)cols_d;
    int *member=(int*)R_alloc((size_t)d*nc,sizeof(int));
    int *z=(int*)R_alloc((size_t)d*d*nc,sizeof(int));
    int *x=(int*)R_alloc((size_t)d*d,sizeof(int));
    int *b=(int*)R_alloc((size_t)d*d,sizeof(int));
    memset(member,0,(size_t)d*nc*sizeof(int)); memset(z,0,(size_t)d*d*nc*sizeof(int));
    memset(x,0,(size_t)d*d*sizeof(int)); memset(b,0,(size_t)d*d*sizeof(int));
    for(int c=0;c<nc;++c) {
        SEXP s=VECTOR_ELT(cc,c);
        for(int k=0;k<LENGTH(s);++k) member[(INTEGER(s)[k]-1)*nc+c]=1;
    }
    const int *zi=INTEGER(VECTOR_ELT(zijc,0)), *zj=INTEGER(VECTOR_ELT(zijc,1)), *zc=INTEGER(VECTOR_ELT(zijc,2));
    for(int q=0;q<nz;++q) { int i=zi[q]-1,j=zj[q]-1,c=zc[q]-1; z[(i*d+j)*nc+c]=z[(j*d+i)*nc+c]=q+1; }
    const int *xi=INTEGER(VECTOR_ELT(xij,0)), *xj=INTEGER(VECTOR_ELT(xij,1));
    for(int q=0;q<nd;++q) x[(xi[q]-1)*d+xj[q]-1]=q+1;
    const int *bi=INTEGER(VECTOR_ELT(bxij,0)), *bj=INTEGER(VECTOR_ELT(bxij,1));
    const int *bid=INTEGER(VECTOR_ELT(bxij,2));
    for(int q=0;q<2*nb;++q) b[(bi[q]-1)*d+bj[q]-1]=bid[q];
    int num_r1=0;
    for(int i=0;i<d;++i) for(int c=0;c<nc;++c) if(!member[i*nc+c]) ++num_r1;
    SEXP out=PROTECT(allocVector(VECSXP,9));
    SET_VECTOR_ELT(out,0,allocVector(INTSXP,5*(R_xlen_t)rows));
    SET_VECTOR_ELT(out,1,allocVector(INTSXP,5*(R_xlen_t)rows));
    SET_VECTOR_ELT(out,2,allocVector(REALSXP,5*(R_xlen_t)rows));
    SET_VECTOR_ELT(out,3,duplicate(rhs_));
    SET_VECTOR_ELT(out,4,allocVector(VECSXP,nz));
    SET_VECTOR_ELT(out,5,allocVector(VECSXP,nd));
    SET_VECTOR_ELT(out,6,allocVector(VECSXP,max_size>1 ? nlb : 0));
    SET_VECTOR_ELT(out,7,allocVector(VECSXP,num_r1));
    SET_VECTOR_ELT(out,8,allocVector(INTSXP,1));
    int *ri=INTEGER(VECTOR_ELT(out,0)), *cj=INTEGER(VECTOR_ELT(out,1));
    double *val=REAL(VECTOR_ELT(out,2)), *rhs=REAL(VECTOR_ELT(out,3));
    int *tmp=(int*)R_alloc((size_t)d+1,sizeof(int));
    int ct=0,cz=nz,cx=2*nz,cl=2*nz+nd,cle=3*nz+nd,cn=3*nz+2*nd;
    int cde=cn+nn,cu=cde+nd,cb=cu+nl,cs=cb+nb,clb=cs+nd;
    int czb=clb+nlb,cf=czb+nlb,cdi=cf+nf;
    int rn=2*nz,ra=rn+nn,rind=ra+nd,rp=rind+2*nd,rzl=rp+nl;
    int rc6=rzl+2*nz,rc10=rc6+3*nd,rf=rc10+2*nlb,rr=rf+nf,rred=rr+nr;
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
    for(int q=0;q<nlb;++q) {
        ADD(rc10+q,czb+q,1); ADD(rc10+q,clb+q,1);
        ADD(rc10+nlb+q,czb+q,1-d); ADD(rc10+nlb+q,clb+q,-1);
    }
    int row=0,pair=0;
    for(int i=0;i<d;++i) for(int j=0;j<d;++j) if(i!=j) {
        int ij=x[i*d+j],ji=x[j*d+i],bij=b[i*d+j],start=row;
        if(dag) { ADD(ra+pair,cde+ij-1,1); ADD(ra+pair,cde+ji-1,1); rhs[ra+pair]=1; }
        int a=rc6+3*pair;
        ADD(a,cx+ij-1,-1); ADD(a,cb+bij-1,-1); ADD(a,cs+ij-1,1);
        ADD(a+1,cx+ij-1,1); ADD(a+1,cs+ij-1,-1);
        ADD(a+2,cb+bij-1,1); ADD(a+2,cs+ij-1,-1);
        ADD(rn+row,cn+row,1); ADD(rn+row,cx+ij-1,d-1); rhs[rn+row]=d; ++row;
        for(int k=0;k<d;++k) if(k!=i && k!=j) {
            ADD(rn+row,cn+row,1); ADD(rn+row,cx+x[k*d+j]-1,d-2);
            ADD(rn+row,cle+x[i*d+k]-1,-1); rhs[rn+row]=d-1; ++row;
        }
        for(int q=0;q<row-start;++q) tmp[q]=cn+start+q+1;
        SET_VECTOR_ELT(VECTOR_ELT(out,5),pair,constraint(cle+ij,tmp,row-start));
        ++pair;
    }
    if(dag) {
        row=0;
        for(int i=0;i<d;++i) for(int j=0;j<i;++j) {
            ADD(rred+row,cx+x[i*d+j]-1,1); ADD(rred+row,cx+x[j*d+i]-1,1);
            rhs[rred+row]=1; ++row;
        }
    }
    row=0;
    for(int i=0;i<d;++i) for(int c=0;c<nc;++c) if(member[i*nc+c]) { ADD(rr+row,cdi+i*nc+c,1); ++row; }
    row=0;
    for(int i=0;i<d;++i) for(int c=0;c<nc;++c) if(!member[i*nc+c]) {
        int count=0;
        for(int k=0;k<d;++k) if(member[k*nc+c]) tmp[count++]=cde+x[i*d+k];
        tmp[count++]=cols;
        SET_VECTOR_ELT(VECTOR_ELT(out,7),row++,constraint(cdi+i*nc+c+1,tmp,count));
    }
    if(max_size>1) {
        int frow=0, group=0;
        for(int j=0;j<d;++j) for(int i=0;i<j;++i) for(int c=0;c<nc;++c) {
            int start=frow;
            if(member[i*nc+c] && member[j*nc+c]) {
                ADD(rf+frow,cf+frow,1); ADD(rf+frow,cb+b[i*d+j]-1,d-1);
                rhs[rf+frow]=d; ++frow;
                for(int k=0;k<d;++k) if(k!=i && k!=j && member[k*nc+c]) {
                    int hi=i>k ? i:k, lo=i>k ? k:i;
                    int ikc=(hi*(hi-1)/2+lo)*nc+c;
                    ADD(rf+frow,cf+frow,1); ADD(rf+frow,cb+b[j*d+k]-1,d-2);
                    ADD(rf+frow,clb+ikc,-1); rhs[rf+frow]=d-1; ++frow;
                }
            }
            int count=0; for(int q=start;q<frow;++q) tmp[count++]=cf+q+1;
            tmp[count++]=cols-1;
            SET_VECTOR_ELT(VECTOR_ELT(out,6),group,constraint(clb+group+1,tmp,count));
            ++group;
        }
    }
    int m=INTEGER(VECTOR_ELT(path,7))[0], ne=INTEGER(VECTOR_ELT(path,7))[1];
    int *pr=INTEGER(VECTOR_ELT(path,0)), *pc=INTEGER(VECTOR_ELT(path,1));
    double *pv=REAL(VECTOR_ELT(path,2));
    for(int q=0;q<ne;++q) {
        int c=pc[q]-1, global;
        if(c<nd) global=cx+c;
        else if(c<nd+nz) global=cl+c-nd;
        else if(c<2*nd+nz) global=cs+c-nd-nz;
        else if(c<2*nd+nz+nlb) global=clb+c-2*nd-nz;
        else if(c<2*nd+nz+2*nlb) global=czb+c-2*nd-nz-nlb;
        else if(c<2*nd+nz+2*nlb+d*nc) global=cdi+c-2*nd-nz-2*nlb;
        else global=cu+c-2*nd-nz-2*nlb-d*nc;
        ADD(rp+pr[q]-1,global,pv[q]);
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
        group[q]=g; ++counts[g]; rhs[rp+q]=pv_rhs[q];
    }
    starts[0]=0; for(int q=1;q<=nz;++q) starts[q]=starts[q-1]+counts[q-1];
    memcpy(pos,starts,((size_t)nz+1)*sizeof(int));
    int *group_rows=(int*)R_alloc((size_t)m+1,sizeof(int));
    for(int q=0;q<m;++q) group_rows[pos[group[q]]++]=cu+q+1;
    for(int q=1;q<=nz;++q) SET_VECTOR_ELT(VECTOR_ELT(out,4),q-1,
        constraint(cl+q,group_rows+starts[q],counts[q]));
#undef ADD
    INTEGER(VECTOR_ELT(out,8))[0]=e;
    UNPROTECT(1); return out;
}
