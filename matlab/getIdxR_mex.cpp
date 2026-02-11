/*
 * getIdxR_mex.cpp  –  MATLAB MEX port of getIdxR() from src/getIdx.cpp
 *
 * Find the indices of the r smallest and r largest values in a vector,
 * excluding a set of already-selected indices.
 * Uses std::nth_element for O(n) selection (same algorithm as the R/C++ version).
 *
 * Usage (from MATLAB, after mex-compiling):
 *   idx = getIdxR_mex(r, z, del)
 *
 *   r   – scalar integer, number of extremes per tail
 *   z   – double column vector (n x 1)
 *   del – double vector of 1-based indices to exclude
 *
 * Returns:
 *   idx – (2*r x 1) double vector of 1-based indices
 *         [r smallest indices; r largest indices]
 */

#include "mex.h"
#include <algorithm>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* ---- input validation ---- */
    if (nrhs != 3)
        mexErrMsgIdAndTxt("getIdxR:nrhs",
                          "Three inputs required: r, z, del.");
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
        mexErrMsgIdAndTxt("getIdxR:notDouble", "z must be a real double vector.");

    int r = (int)mxGetScalar(prhs[0]);
    double *z = mxGetPr(prhs[1]);
    int n = (int)mxGetNumberOfElements(prhs[1]);
    double *del_dbl = mxGetPr(prhs[2]);
    int m = (int)mxGetNumberOfElements(prhs[2]);

    /* ---- convert del to int and sort ---- */
    int *del = new int[m];
    for (int i = 0; i < m; i++)
        del[i] = (int)del_dbl[i];
    std::sort(del, del + m);

    /* ---- build filtered values (exclude del indices) ---- */
    double *y = new double[n - m];
    int j = 0, k = 0;
    for (int i = 0; i < n; i++) {
        if (j >= m) {
            y[k++] = z[i];
        } else if (del[j] != i + 1) {
            y[k++] = z[i];
        } else {
            j++;
        }
    }

    /* ---- find r-th smallest value (lower threshold) ---- */
    std::nth_element(y, y + r - 1, y + n - m);
    double yrl = y[r - 1];

    /* ---- find r-th largest value (upper threshold) ---- */
    j = 0; k = 0;
    for (int i = 0; i < n; i++) {
        if (j >= m) {
            y[k++] = -z[i];
        } else if (del[j] != i + 1) {
            y[k++] = -z[i];
        } else {
            j++;
        }
    }
    std::nth_element(y, y + r - 1, y + n - m);
    double yru = -y[r - 1];

    delete[] y;

    /* ---- collect indices (excluding del) ---- */
    int *locl = new int[r];
    int *locu = new int[r];
    int jl = 0, ju = 0;
    j = 0;

    for (int i = 0; i < n; i++) {
        if (j >= m) {
            if (z[i] <= yrl && jl < r)
                locl[jl++] = i + 1;
            if (z[i] >= yru && ju < r)
                locu[ju++] = i + 1;
        } else if (del[j] != i + 1) {
            if (z[i] <= yrl && jl < r)
                locl[jl++] = i + 1;
            if (z[i] >= yru && ju < r)
                locu[ju++] = i + 1;
        } else {
            j++;
        }
        if (jl >= r && ju >= r)
            break;
    }

    /* ---- build output ---- */
    plhs[0] = mxCreateDoubleMatrix(2 * r, 1, mxREAL);
    double *idx = mxGetPr(plhs[0]);

    for (int i = 0; i < r; i++) {
        idx[i]     = (double)locl[i];
        idx[r + i] = (double)locu[i];
    }

    delete[] locl;
    delete[] locu;
    delete[] del;
}
