#include "solvers.hpp"
#include <cmath>
#include <algorithm>

vec gaussNoPivoting(mat A, vec b) {
    int n = A.size();
    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }
    vec x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }
    return x;
}

vec gaussWithPivoting(mat A, vec b) {
    int n = A.size();
    for (int k = 0; k < n - 1; ++k) {
        int maxRow = k;
        double maxEl = std::abs(A[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(A[i][k]) > maxEl) {
                maxEl = std::abs(A[i][k]);
                maxRow = i;
            }
        }
        std::swap(A[k], A[maxRow]);
        std::swap(b[k], b[maxRow]);
        
        for (int i = k + 1; i < n; ++i) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }
    vec x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }
    return x;
}

void luDecomposition(const mat& A, mat& L, mat& U) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) sum += (L[i][j] * U[j][k]);
            U[i][k] = A[i][k] - sum;
        }
        for (int k = i; k < n; k++) {
            if (i == k) L[i][i] = 1.0;
            else {
                double sum = 0.0;
                for (int j = 0; j < i; j++) sum += (L[k][j] * U[j][i]);
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

vec luSolve(const mat& L, const mat& U, const vec& b) {
    int n = L.size();
    vec y(n), x(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) sum += L[i][j] * y[j];
        y[i] = b[i] - sum;
    }
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += U[i][j] * x[j];
        x[i] = (y[i] - sum) / U[i][i];
    }
    return x;
}
