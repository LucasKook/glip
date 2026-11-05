#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <limits.h>
#include <string.h>

/* Construct one Gurobi minimum constraint with numeric one-based indices. */
static SEXP min_constraint(int result, const int *vars, int n) {
    SEXP ans = PROTECT(allocVector(VECSXP, 2));
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SEXP result_sexp = PROTECT(ScalarReal(result));
    SEXP v = PROTECT(allocVector(REALSXP, n));
    for (int i = 0; i < n; ++i) REAL(v)[i] = vars[i];
    SET_VECTOR_ELT(ans, 0, result_sexp);
    SET_VECTOR_ELT(ans, 1, v);
    SET_STRING_ELT(names, 0, mkChar("resvar"));
    SET_STRING_ELT(names, 1, mkChar("vars"));
    setAttrib(ans, R_NamesSymbol, names);
    UNPROTECT(4);
    return ans;
}

/* All index offsets below are zero-based. The emitted triplets and Gurobi
   variable indices are one-based, matching the R implementation. */
static SEXP dcon_build(SEXP d_, SEXP cc, SEXP zijc, SEXP xij,
                       SEXP nl_, SEXP nr_, SEXP rhs_, SEXP dag_) {
    int d = asInteger(d_), nc = LENGTH(cc), nl = asInteger(nl_);
    int nr = asInteger(nr_), dag = asLogical(dag_);
    int nz = LENGTH(VECTOR_ELT(zijc, 0));
    double nd_d = (double)d * (d - 1), nn_d = nd_d * (d - 1);
    double rows_d = 4.0*nz + nn_d + 4.0*nd_d + nl + nr + nd_d/2;
    double cols_d = 3.0*nz + 3.0*nd_d + nn_d + nl + (double)d*nc + 1;
    if (d < 2 || nl < 0 || nr < 0 || rows_d > INT_MAX/4 || cols_d > INT_MAX)
        error("DAG encoding is too large for a sparse Matrix");
    int nd = (int)nd_d, nn = (int)nn_d, rows = (int)rows_d;
    int cols = (int)cols_d;
    size_t membership_size = (size_t)d * nc;
    size_t pair_size = (size_t)d * d * nc;
    int *member = (int *) R_alloc(membership_size, sizeof(int));
    int *z = (int *) R_alloc(pair_size, sizeof(int));
    int *x = (int *) R_alloc((size_t)d * d, sizeof(int));
    memset(member, 0, membership_size * sizeof(int));
    memset(z, 0, pair_size * sizeof(int));
    memset(x, 0, (size_t)d * d * sizeof(int));
    for (int c = 0; c < nc; ++c) {
        SEXP s = VECTOR_ELT(cc, c);
        for (int q = 0; q < LENGTH(s); ++q) member[(INTEGER(s)[q]-1)*nc+c] = 1;
    }
    int *zi=INTEGER(VECTOR_ELT(zijc,0)), *zj=INTEGER(VECTOR_ELT(zijc,1)), *zc=INTEGER(VECTOR_ELT(zijc,2));
    for (int q=0; q<nz; ++q) {
        int i=zi[q]-1, j=zj[q]-1, c=zc[q]-1;
        z[(i*d+j)*nc+c]=z[(j*d+i)*nc+c]=q+1;
    }
    int *xi=INTEGER(VECTOR_ELT(xij,0)), *xj=INTEGER(VECTOR_ELT(xij,1));
    for (int q=0; q<nd; ++q) x[(xi[q]-1)*d+xj[q]-1]=q+1;

    SEXP out=PROTECT(allocVector(VECSXP,8));
    SET_VECTOR_ELT(out,0,allocVector(INTSXP,4*(R_xlen_t)rows));
    SET_VECTOR_ELT(out,1,allocVector(INTSXP,4*(R_xlen_t)rows));
    SET_VECTOR_ELT(out,2,allocVector(REALSXP,4*(R_xlen_t)rows));
    SET_VECTOR_ELT(out,3,duplicate(rhs_));
    SET_VECTOR_ELT(out,4,allocVector(VECSXP,nz));
    SET_VECTOR_ELT(out,5,allocVector(VECSXP,nd));
    int count_r1 = 0;
    for (int i=0; i<d; ++i) for (int c=0; c<nc; ++c) if (!member[i*nc+c]) ++count_r1;
    SET_VECTOR_ELT(out,6,allocVector(VECSXP,count_r1));
    SET_VECTOR_ELT(out,7,allocVector(INTSXP,1));
    int *ri=INTEGER(VECTOR_ELT(out,0)), *cj=INTEGER(VECTOR_ELT(out,1));
    double *val=REAL(VECTOR_ELT(out,2)), *rhs=REAL(VECTOR_ELT(out,3));
    int *group_count=(int *)R_alloc((size_t)nz+1,sizeof(int));
    int *group_start=(int *)R_alloc((size_t)nz+1,sizeof(int));
    int *group_pos=(int *)R_alloc((size_t)nz+1,sizeof(int));
    int *m_group=(int *)R_alloc((size_t)nl+1,sizeof(int));
    memset(group_count,0,((size_t)nz+1)*sizeof(int));
    int *tmp=(int *)R_alloc((size_t)d+1,sizeof(int));
    int e=0;
    /* Each row has at most four coefficients; skip zero coefficients. */
#define ADD(r,c,v) do { double av=(v); if (av != 0) { ri[e]=(r)+1; cj[e]=(c)+1; val[e]=av; ++e; } } while(0)
    int ct=0, cz=nz, cx=2*nz, cl=2*nz+nd, cle=3*nz+nd;
    int cn=3*nz+2*nd, cde=cn+nn, cu=cde+nd, cdi=cu+nl;
    int rlabs=0, rn=2*nz, racyc=rn+nn, rind=racyc+2*nd;
    int rm=rind+2*nd, rzl=rm+nl, rr1=rzl+2*nz, rred=rr1+nr;
    for (int q=0; q<nz; ++q) {
        ADD(rlabs+q,ct+q,1); ADD(rlabs+q,cz+q,-1);
        ADD(rlabs+nz+q,ct+q,1); ADD(rlabs+nz+q,cz+q,1);
        ADD(rzl+q,cz+q,1); ADD(rzl+q,cl+q,1);
        ADD(rzl+nz+q,cz+q,-(d-1)); ADD(rzl+nz+q,cl+q,-1);
    }
    for (int q=0; q<nd; ++q) {
        ADD(rind+q,cle+q,1); ADD(rind+q,cde+q,-1);
        ADD(rind+nd+q,cle+q,-1); ADD(rind+nd+q,cde+q,d-1);
    }
    int row=0, pair=0;
    for (int i=0;i<d;++i) for (int j=0;j<d;++j) if(i!=j) {
        int ij=x[i*d+j], ji=x[j*d+i], start=row;
        if (dag) {
            ADD(racyc+pair,cde+ij-1,1); ADD(racyc+pair,cde+ji-1,1);
            rhs[racyc+pair]=1;
        }
        ADD(rn+row,cn+row,1); ADD(rn+row,cx+ij-1,d-1); rhs[rn+row]=d; ++row;
        for(int k=0;k<d;++k) if(k!=i && k!=j) {
            ADD(rn+row,cn+row,1); ADD(rn+row,cx+x[k*d+j]-1,d-2);
            ADD(rn+row,cle+x[i*d+k]-1,-1); rhs[rn+row]=d-1; ++row;
        }
        for(int q=0;q<row-start;++q) tmp[q]=cn+start+q+1;
        SET_VECTOR_ELT(VECTOR_ELT(out,5),pair,min_constraint(cle+ij,tmp,row-start));
        ++pair;
    }
    if (dag) {
        row=0;
        for(int i=0;i<d;++i) for(int j=0;j<i;++j) {
            ADD(rred+row,cx+x[i*d+j]-1,1); ADD(rred+row,cx+x[j*d+i]-1,1);
            rhs[rred+row]=1; ++row;
        }
    }
    row=0;
    for(int i=0;i<d;++i) for(int c=0;c<nc;++c) if(member[i*nc+c]) {
        ADD(rr1+row,cdi+i*nc+c,1); ++row;
    }
    row=0;
#define MROW(i,j,c,r) do { int group=z[((i)*d+(j))*nc+(c)]; m_group[row]=group; ++group_count[group]; rhs[rm+row]=(r); ADD(rm+row,cu+row,1); } while(0)
    for(int i=0;i<d;++i) {
        R_CheckUserInterrupt();
        for(int j=0;j<d;++j) if(i!=j) {
            for(int c=0;c<nc;++c) if(!member[i*nc+c] && !member[j*nc+c]) {
                MROW(i,j,c,d); ADD(rm+row,cx+x[i*d+j]-1,d-1); ++row;
            }
            for(int k=0;k<d;++k) if(k!=i && k!=j) for(int c=0;c<nc;++c)
                if(!member[i*nc+c] && !member[j*nc+c]) {
                    if(!member[k*nc+c]) {
                        MROW(i,j,c,d-1); ADD(rm+row,cx+x[k*d+j]-1,d-2);
                        ADD(rm+row,cl+z[(i*d+k)*nc+c]-1,-1); ++row;
                        MROW(i,j,c,0); ADD(rm+row,cl+z[(i*d+k)*nc+c]-1,-1);
                        ADD(rm+row,cl+z[(k*d+j)*nc+c]-1,-1);
                        ADD(rm+row,cdi+k*nc+c,-(d-2)); ++row;
                    } else {
                        MROW(i,j,c,2*(d-2)+2); ADD(rm+row,cx+x[i*d+k]-1,d-2);
                        ADD(rm+row,cx+x[j*d+k]-1,d-2); ++row;
                    }
                }
        }
    }
#undef MROW
    group_start[0]=0;
    for(int q=1;q<=nz;++q) group_start[q]=group_start[q-1]+group_count[q-1];
    memcpy(group_pos,group_start,((size_t)nz+1)*sizeof(int));
    int *group_rows=(int *)R_alloc((size_t)row+1,sizeof(int));
    for(int q=0;q<row;++q) group_rows[group_pos[m_group[q]]++]=cu+q+1;
    for(int q=1;q<=nz;++q) SET_VECTOR_ELT(VECTOR_ELT(out,4),q-1,
        min_constraint(cl+q,group_rows+group_start[q],group_count[q]));
    int r1=0;
    for(int i=0;i<d;++i) for(int c=0;c<nc;++c) if(!member[i*nc+c]) {
        int count=0;
        /* x indices increase with j for fixed i in expand.grid ordering. */
        for(int k=0;k<d;++k) if(member[k*nc+c]) tmp[count++]=cde+x[i*d+k];
        tmp[count++]=cols;
        SET_VECTOR_ELT(VECTOR_ELT(out,6),r1++,min_constraint(cdi+i*nc+c+1,tmp,count));
    }
#undef ADD
    INTEGER(VECTOR_ELT(out,7))[0]=e;
    UNPROTECT(1);
    return out;
}

/* Enumerate indices in the original i, j, conditioning-set order. */
static SEXP dcon_indices(SEXP d_, SEXP cc) {
    int d = asInteger(d_), nc = LENGTH(cc);
    if (d < 2) error("DAG encoding requires at least two vertices");
    size_t size = (size_t)d * nc;
    double max_pairs=(double)d*(d-1)*nc/2;
    if (size > INT_MAX || max_pairs > INT_MAX)
        error("DAG index table is too large");
    int *member = (int *)R_alloc(size, sizeof(int));
    memset(member, 0, size * sizeof(int));
    for (int c=0; c<nc; ++c) {
        SEXP s=VECTOR_ELT(cc,c);
        for (int k=0; k<LENGTH(s); ++k) member[(INTEGER(s)[k]-1)*nc+c]=1;
    }
    int nz=0;
    for(int i=0;i<d;++i) for(int j=0;j<i;++j) for(int c=0;c<nc;++c)
        if(!member[i*nc+c] && !member[j*nc+c]) ++nz;
    SEXP out=PROTECT(allocVector(VECSXP,4));
    for(int k=0;k<3;++k) SET_VECTOR_ELT(out,k,allocVector(INTSXP,nz));
    SET_VECTOR_ELT(out,3,allocVector(LGLSXP,(R_xlen_t)size));
    int *zi=INTEGER(VECTOR_ELT(out,0)), *zj=INTEGER(VECTOR_ELT(out,1));
    int *zc=INTEGER(VECTOR_ELT(out,2)), *inc=LOGICAL(VECTOR_ELT(out,3));
    int n=0;
    for(int i=0;i<d;++i) {
        R_CheckUserInterrupt();
        for(int c=0;c<nc;++c) inc[i*nc+c]=member[i*nc+c];
        for(int j=0;j<i;++j) for(int c=0;c<nc;++c)
            if(!member[i*nc+c] && !member[j*nc+c]) {
                zi[n]=i+1; zj[n]=j+1; zc[n]=c+1; ++n;
            }
    }
    UNPROTECT(1);
    return out;
}

extern SEXP glip_lean_path_block(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP glip_chain_build(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP glip_admg_build(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
    {"glip_dcon_build", (DL_FUNC) &dcon_build, 8},
    {"glip_dcon_indices", (DL_FUNC) &dcon_indices, 2},
    {"glip_lean_path_block", (DL_FUNC) &glip_lean_path_block, 7},
    {"glip_chain_build", (DL_FUNC) &glip_chain_build, 9},
    {"glip_admg_build", (DL_FUNC) &glip_admg_build, 11},
    {NULL, NULL, 0}
};
void R_init_glip(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
