/*
 * getIdx_mex.cpp  –  MATLAB MEX port of getIdx() from src/getIdx.cpp
 *
 * Find the indices of the r smallest and r largest values in a vector.
 * Uses std::nth_element for O(n) selection (same algorithm as the R/C++ version).
 *
 * Usage (from MATLAB, after mex-compiling):
 *   idx = getIdx_mex(r, z)
 *
 *   r  – scalar integer, number of extremes per tail
 *   z  – double column vector (n x 1)
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
    if (nrhs != 2)
        mexErrMsgIdAndTxt("getIdx:nrhs",
                          "Two inputs required: r (integer scalar), z (double vector).");
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
        mexErrMsgIdAndTxt("getIdx:notDouble", "z must be a real double vector.");

    int r = (int)mxGetScalar(prhs[0]);
    double *z = mxGetPr(prhs[1]);
    int n = (int)mxGetNumberOfElements(prhs[1]);

    if (r < 1 || 2 * r > n)
        mexErrMsgIdAndTxt("getIdx:badR",
                          "r must satisfy 1 <= r <= n/2.");

    /* ---- find r-th smallest value (lower threshold) ---- */
    double *y = new double[n];
    for (int i = 0; i < n; i++)
        y[i] = z[i];

    std::nth_element(y, y + r - 1, y + n);
    double yrl = y[r - 1];

    /* ---- find r-th largest value (upper threshold) ---- */
    for (int i = 0; i < n; i++)
        y[i] = -z[i];

    std::nth_element(y, y + r - 1, y + n);
    double yru = -y[r - 1];

    delete[] y;

    /* ---- collect indices ---- */
    int *locl = new int[r];
    int *locu = new int[r];
    int jl = 0, ju = 0;

    for (int i = 0; i < n; i++) {
        if (z[i] <= yrl && jl < r)
            locl[jl++] = i + 1;          /* 1-based for MATLAB */
        if (z[i] >= yru && ju < r)
            locu[ju++] = i + 1;
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
}
