/**
 * @file fast_sweeping_mex.cpp
 * @brief MEX interface for solving the Eikonal equation using the Fast Sweeping Method in 2D.
 *
 * This MEX function takes a 2D speed function and initial mask as inputs and returns the
 * solution of the Eikonal equation.
 *
 * Input:
 *   - prhs[0]: 2D speed function (N x M) of type double
 *   - prhs[1]: 2D binary mask (N x M) of type double, serving as the initial conditions to the Eikonal equation. Mask should be = 1 at the source points and 0 elsewhere. 
 *
 * Output:
 *   - plhs[0]: 2D array (N x M) representing the solved Eikonal equation.
 *
 * Author: Liam Burrows
 * Date: 18/09/2024
 * Contact: lb2668@bath.ac.uk
 *
 * Implementation of [1]
 * [1] Zhao, Hongkai. "A fast sweeping method for eikonal equations." Mathematics of computation 74.250 (2005): 603-627.
 */

#include "mex.h"
#include <vector>
#include "fast_sweeping.h"

using namespace std;

//const double INF = numeric_limits<double>::max();




void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Inputs (f, mask): f - speed function. mask - binary mask indicating source (1 in source)
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:invalidNumInputs", "Two inputs are required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:invalidNumOutputs", "One output is required.");
    }
    
    double *speed = mxGetPr(prhs[0]);
    double *mask = mxGetPr(prhs[1]);
    
    int N = mxGetM(prhs[0]); //Num. rows
    int M = mxGetN(prhs[0]); //


    // Convert MATLAB arrays to 2D vectors
    vector<vector<double>> u(N, vector<double>(M, INF));
    vector<vector<double>> f(N, vector<double>(M));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            u[i][j] = INF*(1-mask[i + N * j]);  
            f[i][j] = speed[i + N * j];  
        }
    }

    fast_sweeping(u, f, N, M);
    
    plhs[0] = mxCreateDoubleMatrix(N, M, mxREAL);
    double* output = mxGetPr(plhs[0]);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            output[i + N * j] = u[i][j];
        }
    }
}