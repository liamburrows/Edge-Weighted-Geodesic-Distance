/**
 * @file fast_sweeping_mex.cpp
 * @brief MEX interface for solving the Eikonal equation using the Fast Sweeping Method in 3D.
 *
 * This MEX function takes a 3D speed function and initial mask as inputs and returns the
 * solution of the Eikonal equation.
 *
 * Input:
 *   - prhs[0]: 3D speed function (N x M x SL) of type double
 *   - prhs[1]: 3D binary mask (N x M X SL) of type double, serving as the initial conditions to the Eikonal equation. Mask should be = 1 at the source points and 0 elsewhere. 
 *
 * Output:
 *   - plhs[0]: 3D array (N x M x SL) representing the solved Eikonal equation.
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
#include <iostream>

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

    double* speed = mxGetPr(prhs[0]);
    double* mask = mxGetPr(prhs[1]);

    const mwSize *dims = mxGetDimensions(prhs[0]);
    //cout << dims;
    size_t N = dims[0]; //Num. rows
    size_t M = dims[1]; //
    size_t SL = dims[2];

    cout << "N = " << N << ", M = " << M << ", SL = " << SL << endl;


    // Convert MATLAB arrays to 2D vectors
    vector<vector<vector<double>>> u(N, vector<vector<double>>(M, vector<double>(SL, INF)));
    vector<vector<vector<double>>> f(N, vector<vector<double>>(M, vector<double>(SL)));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            for (int k = 0; k < SL; ++k) {
                size_t idx = i + N * (j + M * k);
                u[i][j][k] = INF * (1 - mask[idx]);
                f[i][j][k] = speed[idx];
            }
        }
    }

    cout << "u at point = " << u[199][499][19] << endl;
    cout << "u at point2 = " << u[149][249][19] << endl;

    fast_sweeping3D(u, f, N, M, SL);

    cout << "u at point = " << u[199][499][19] << endl;
    cout << "u at point2 = " << u[149][249][19] << endl;

    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double* output = mxGetPr(plhs[0]);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            for (int k = 0; k < SL; ++k) {
                size_t idx = i + N * (j + M * k);
                output[idx] = u[i][j][k];
            }
        }
    }
}