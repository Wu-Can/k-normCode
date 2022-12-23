#include "mex.h"
#include "matrix.h"
#include "math.h" 

#ifndef mwSize
#define mwSize size_t
#endif
#ifndef mwIndex
#define mwIndex size_t 
#endif

#define A_in  prhs[0]
#define xt_in prhs[1]
#define ATDAvec_in prhs[2]
#define DA_in prhs[3]
#define Ax_in prhs[4]
#define lamATbeta_in prhs[5]
#define gradf_Hess_in prhs[6]
#define lamkk_in prhs[7]
#define Hessvec_in prhs[8]

#define xt_out plhs[0]
#define Ax_out plhs[1]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /* Inputs:
     *
     * A_in: data matrix
     * xt_in: initial point
     * ATDAvec_in: the vector of diagonal component of A'*D*A
     * DA_in: the matrix of D*A
     * Ax_in: the vector of A*x
     * lamATbeta_in: the vector of lambda*A'*beta
     * gradf_Hess_in: the vector of gradient and other information
     * lamkk_in: the parameter of lambda*k
     * Hessvec_in: the vector of diagonal component of Hession

     * Outputs:
     * xt_out: the output of a cyclic coordinate minimization
     * Ax_out: the vector of A*x
     */

    mwIndex *ir, *jc, i;
    mwSize nrow, dcol, dtmp;

    double *A, *xt, *ATDAvec, *DA, *Ax, *lamATbeta, *gradf_Hess, *Hessvec, *xt_new, *Ax_new;
    double lamkk;

    mwSize k, stop;
    double DAAx_ii, alpha = 0, tmp = 0, del_xt = 0;

    /* Check for proper number of arguments */
    if (nlhs != 2) {
        mexErrMsgTxt("OnceCD.c : requires 2 outputs");
    }
    if (nrhs != 9) {
        mexErrMsgTxt("OnceCD.c : requires 9 inputs");
    }

    /* Check for proper format of inputs */

    if (!mxIsSparse(A_in)) {
        mexErrMsgTxt("OnceCD.c : the input data matrix A must be sparse");
    }
    if (!mxIsSparse(DA_in)) {
        mexErrMsgTxt("OnceCD.c : the input matrix DA must be sparse");
    }

    nrow = mxGetM(A_in);
    dcol = mxGetN(A_in);
    dtmp = dcol + 1;

    /* Allocate output */
    xt_out = mxCreateDoubleMatrix(dtmp, 1, mxREAL);
    Ax_out = mxCreateDoubleMatrix(nrow, 1, mxREAL);

    /* I/O pointers */
    A = mxGetPr(A_in);
    xt = mxGetPr(xt_in);
    ATDAvec = mxGetPr(ATDAvec_in);
    DA = mxGetPr(DA_in);
    Ax = mxGetPr(Ax_in);
    lamATbeta = mxGetPr(lamATbeta_in);
    gradf_Hess = mxGetPr(gradf_Hess_in);
    lamkk = mxGetScalar(lamkk_in);
    Hessvec = mxGetPr(Hessvec_in);

    xt_new = mxGetPr(xt_out);
    Ax_new = mxGetPr(Ax_out);

    ir = mxGetIr(A_in);      /* Row indexing      */
    jc = mxGetJc(A_in);      /* Column count      */

    /* Initialization */
    for (i = 0; i < dtmp; i++) {
        xt_new[i] = xt[i];
    }

    for (i = 0; i < nrow; i++) {
        Ax_new[i] = Ax[i];
    }

    /* A cycle of Coordinate Minimization */
    for (i = 0; i < dcol; i++) {
        DAAx_ii = 0;
        stop = jc[i + 1];
        for (k = jc[i]; k < stop; k++) {
            DAAx_ii += DA[k] * Ax_new[ir[k]];
        }

        alpha = (xt_new[i]*ATDAvec[i] - DAAx_ii - xt[dcol]*lamATbeta[i] - gradf_Hess[i])/ Hessvec[i];
        tmp = lamkk/Hessvec[i];

        if (alpha > tmp) {
            xt_new[i] = alpha - tmp;
        }
        else if (alpha < -tmp) {
            xt_new[i] = alpha + tmp;
        }
        else {
            xt_new[i] = 0;
        }

        del_xt = xt_new[i] - xt[i];
        if (del_xt > 1e-12 || del_xt <-1e-12) {
            for (k = jc[i]; k < stop; k++) {
                Ax_new[ir[k]] += del_xt * A[k];
            }
        }

    }

}