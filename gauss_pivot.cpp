#ifndef GAUSS_PIVOT_H
#define GAUSS_PIVOT_H

#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

// Метод Гаусса с частичным выбором ведущего элемента
bool gaussPivot(vector<vector<double>> A, vector<double> b, vector<double>& x) {
    int n = A.size();
    x.assign(n, 0.0);
    
    for (int k = 0; k < n; ++k) {
        int best = k;
        double best_val = fabs(A[k][k]);
        
        for (int i = k + 1; i < n; ++i) {
            if (fabs(A[i][k]) > best_val) {
                best_val = fabs(A[i][k]);
                best = i;
            }
        }
        
        if (best != k) {
            swap(A[k], A[best]);
            swap(b[k], b[best]);
        }
        
        if (fabs(A[k][k]) < 1e-12) {
            return false;
        }
        
        for (int i = k + 1; i < n; ++i) {
            double m = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) {
                A[i][j] -= m * A[k][j];
            }
            b[i] -= m * b[k];
        }
    }
    
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
    return true;
}

#endif
