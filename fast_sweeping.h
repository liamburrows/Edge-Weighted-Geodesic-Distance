#ifndef FAST_SWEEPING_H
#define FAST_SWEEPING_H

#include <vector>

const double INF = std::numeric_limits<double>::max();

double calculate_residual(const std::vector<std::vector<double>>& u, const std::vector<std::vector<double>>& oldu, int N, int M);
void fast_sweeping(std::vector<std::vector<double>>& u, const std::vector<std::vector<double>>& f, int N, int M);
void fast_sweeping3D(std::vector<std::vector<std::vector<double>>>& u, const std::vector<std::vector<std::vector<double>>>& f, int N, int M, int SL);

#endif