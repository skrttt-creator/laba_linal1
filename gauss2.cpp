#ifndef GAUSS2_H
#define GAUSS2_H

#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

bool gauss2(vector<vector<double>> a, vector<double> b, vector<double>& x) {
    int n = a.size();
    x.assign(n, 0.0);
    
    for (int k = 0; k < n; ++k) {
        int best = k;
        double best_val = fabs(a[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (fabs(a[i][k]) > best_val) {
                best_val = fabs(a[i][k]);
                best = i;
            }
        }
        
        if (best != k) {
            swap(a[k], a[best]);
            swap(b[k], b[best]);
        }
        
        if (fabs(a[k][k]) < 1e-12) {
            return false;
        }
        
        for (int i = k + 1; i < n; ++i) {
            double m = a[i][k] / a[k][k];
            for (int j = k; j < n; ++j) {
                a[i][j] -= m * a[k][j];
            }
            b[i] -= m * b[k];
        }
    }
    
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
    return true;
}

#endif
