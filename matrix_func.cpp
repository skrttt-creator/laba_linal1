#include "matrix_utils.hpp"
#include <random>
#include <cmath>

mat generateRandomMatrix(int n, int seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    mat A(n, vec(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = dis(gen);
        }
    }
    return A;
}

vec generateRandomVector(int n, int seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    vec b(n);
    for (int i = 0; i < n; ++i) {
        b[i] = dis(gen);
    }
    return b;
}

mat generateHilbertMatrix(int n) {
    mat H(n, vec(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            H[i][j] = 1.0 / (i + j + 1);
        }
    }
    return H;
}

vec multiplyMatrixVector(const mat& A, const vec& x) {
    int n = A.size();
    vec b(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            b[i] += A[i][j] * x[j];
        }
    }
    return b;
}

double vectorNorm(const vec& v) {
    double sum = 0.0;
    for (double val : v) sum += val * val;
    return std::sqrt(sum);
}

vec subtractVectors(const vec& a, const vec& b) {
    vec res(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] - b[i];
    }
    return res;
}
