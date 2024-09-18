/**
 * @file fast_sweeping.cpp
 * @brief Implementation of the Fast Sweeping Method for the Eikonal equation in 2D and 3D.
 *
 * This file contains the main algorithm for solving the Eikonal equation using the
 * Fast Sweeping Method in 2D and 3D space. It iterates over all grid points in multiple
 * sweeping directions.
 *
 * Author: Liam Burrows
 * Date: 18/09/2024
 * Contact: lb2668@bath.ac.uk
 * 
 * Implementation of [1]
 * [1] Zhao, Hongkai. "A fast sweeping method for eikonal equations." Mathematics of computation 74.250 (2005): 603-627.
 */


#include "fast_sweeping.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

using namespace std;

double calculate_residual(const vector<vector<double>>& u, const vector<vector<double>>& oldu, int N, int M) {
    double res = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            res += pow(oldu[i][j] - u[i][j], 2);
        }
    }
    return res;
}

double calculate_residual3D(const vector<vector<vector<double>>>& u, const vector<vector<vector<double>>>& oldu, int N, int M, int SL) {
    double res = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            for (int k = 0; k < SL; ++k) {
                res += pow(oldu[i][j][k] - u[i][j][k], 2);
            }
        }
    }
    return res;
}

void fast_sweeping(vector<vector<double>>& u, const vector<vector<double>>& f, int N, int M) {

    // define sweep orders: 1:N, 1:M
    vector<int> xx(N), yy(M);

    // Initialise xx and yy
    for (int i = 0; i < N; ++i) xx[i] = i;
    for (int j = 0; j < M; ++j) yy[j] = j;

    // Order for sweeps
    vector<vector<int>> orderx = { xx, xx, vector<int>(xx.rbegin(), xx.rend()), vector<int>(xx.rbegin(), xx.rend()) };
    vector<vector<int>> ordery = { yy, vector<int>(yy.rbegin(), yy.rend()), vector<int>(yy.rbegin(), yy.rend()), yy };

    const double h = 0.005;
    const double stop = .0005;
    double R = 10 * stop;
    int order{};
    int count = 0;

    while (R > stop) {
    //for (int ii=0; ii<100; ++ii) {

        vector<vector<double>> oldu = u;

        order = count % 4;
        vector<int> x = orderx[order];
        vector<int> y = ordery[order];

        // Do sweep
        for (int i : x) {
            for (int j : y) {

                double u_x = INF, u_y = INF;

                if (i > 0) u_x = min(u_x, u[i - 1][j]);
                if (i < N - 1) u_x = min(u_x, u[i + 1][j]);

                if (j > 0) u_y = min(u_y, u[i][j - 1]);
                if (j < M - 1) u_y = min(u_y, u[i][j + 1]);

                double condition = abs(u_x - u_y);

                double ubar{};
                double fh = f[i][j] * h;
                if (condition >= f[i][j] * h) {
                    ubar = min(u_x, u_y) + fh;
                }
                else {
                    ubar = 0.5 * (u_x + u_y + sqrt(2 * pow(fh, 2) - pow(u_x - u_y, 2)));
                }

                u[i][j] = min(u[i][j], ubar);
            }
        }

        R = calculate_residual(u, oldu, N, M);
        count++;

    }


}


void fast_sweeping3D(vector<vector<vector<double>>>& u, const vector<vector<vector<double>>>& f, int N, int M, int SL) {

    // define sweep orders: 1:N, 1:M
    vector<int> xx(N), yy(M), zz(SL);

    // Initialise xx and yy
    for (int i = 0; i < N; ++i) xx[i] = i;
    for (int j = 0; j < M; ++j) yy[j] = j;
    for (int k = 0; k < SL; ++k) zz[k] = k;

        // Order for sweeps
    vector<vector<int>> orderx = { xx, xx, xx, xx, vector<int>(xx.rbegin(), xx.rend()), vector<int>(xx.rbegin(), xx.rend()), vector<int>(xx.rbegin(), xx.rend()), vector<int>(xx.rbegin(), xx.rend()) };
    vector<vector<int>> ordery = { yy, yy, vector<int>(yy.rbegin(), yy.rend()), vector<int>(yy.rbegin(), yy.rend()), yy, yy, vector<int>(yy.rbegin(), yy.rend()), vector<int>(yy.rbegin(), yy.rend()) };
    vector<vector<int>> orderz = { zz, vector<int>(zz.rbegin(), zz.rend()), zz, vector<int>(zz.rbegin(), zz.rend()), zz, vector<int>(zz.rbegin(), zz.rend()), zz, vector<int>(zz.rbegin(), zz.rend()) };


    const double h = 0.0005;
    const double stop = 0.1;
    double R = 10 * stop;
    int order{};
    int count = 0;


    while (R > stop) {
        vector<vector<vector<double>>> oldu = u;

        order = count % 8;
        vector<int> x = orderx[order];
        vector<int> y = ordery[order];
        vector<int> z = orderz[order];
        

        cout << "Iteration: " << count << endl;

        // Do sweep
        for (int i : x) {
            for (int j : y) {
                for (int k : z) {

                    double u_x = INF, u_y = INF, u_z = INF;

                    if (i > 0) u_x = min(u_x, u[i - 1][j][k]);
                    if (i < N - 1) u_x = min(u_x, u[i + 1][j][k]);

                    if (j > 0) u_y = min(u_y, u[i][j - 1][k]);
                    if (j < M - 1) u_y = min(u_y, u[i][j + 1][k]);

                    if (k > 0) u_z = min(u_z, u[i][j][k - 1]);
                    if (k < SL - 1) u_z = min(u_z, u[i][j][k + 1]);

                    vector<double> u_values = { u_x, u_y, u_z };
                    sort(u_values.begin(), u_values.end());

                    
                    double a1 = u_values[0];
                    double a2 = u_values[1];
                    double a3 = u_values[2];

                    double fh = h * f[i][j][k];

                    double xbar{};
                    xbar = a1 + fh;
                    
                    

                    if (xbar > a2) {
                        xbar = 0.5 * ((a1 + a2) + sqrt(2 * pow(fh, 2) - pow(a1 - a2, 2)));

                        if (xbar > a3) {
                            xbar = (1.0 / 3.0) * ((a1 + a2 + a3) + sqrt(3 * pow(fh, 2) - 2 * (pow(a1, 2) + pow(a2, 2) + pow(a3, 2) - a1 * a2 - a1 * a3 - a2 * a3)));
                        }

                    }


                    u[i][j][k] = min(u[i][j][k], xbar);
                    

                   
                }
            }
        }

        R = calculate_residual3D(u, oldu, N, M, SL);

        count++;
        //cout << R;
    }


}