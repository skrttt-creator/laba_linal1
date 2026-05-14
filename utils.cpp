#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <random>
#include <cmath>

using namespace std;

using Matr = vector<vector<double>>;
using Vec = vector<double>;

Matr copy_mat(const Matr& src) {
    int n = src.size();
    Matr dst(n, Vec(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dst[i][j] = src[i][j];
        }
    }
    return dst;
}

Matr rand_mat(int n, unsigned seed = 42) {
    mt19937 rng(seed);
    uniform_real_distribution<double> dist(-1.0, 1.0);
    
    Matr a(n, Vec(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i][j] = dist(rng);
        }
    }
    return a;
}

Vec rand_rhs(int n, unsigned seed = 123) {
    mt19937 rng(seed);
    uniform_real_distribution<double> dist(-1.0, 1.0);
    
    Vec b(n);
    for (int i = 0; i < n; ++i) {
        b[i] = dist(rng);
    }
    return b;
}

Matr hilbert_mat(int n) {
    Matr h(n, Vec(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            h[i][j] = 1.0 / (i + j + 1);
        }
    }
    return h;
}

double norm_resid(const Matr& a, const Vec& x, const Vec& b) {
    int n = a.size();
    double max_err = 0.0;
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) {
            s += a[i][j] * x[j];
        }
        double err = fabs(s - b[i]);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

#endif
