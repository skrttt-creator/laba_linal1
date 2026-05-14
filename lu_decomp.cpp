#ifndef LU_DECOMP_H
#define LU_DECOMP_H

#include <vector>
#include <cmath>

using namespace std;

bool lu_decomp(const vector<vector<double>>& a, vector<vector<double>>& l, vector<vector<double>>& u) {
    int n = a.size();
    l.assign(n, vector<double>(n, 0.0));
    u.assign(n, vector<double>(n, 0.0));
    
    for (int i = 0; i < n; ++i) {
        l[i][i] = 1.0;
    }
    
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double s = 0.0;
            for (int k = 0; k < i; ++k) {
                s += l[i][k] * u[k][j];
            }
            u[i][j] = a[i][j] - s;
        }
        
        for (int j = i + 1; j < n; ++j) {
            double s = 0.0;
            for (int k = 0; k < i; ++k) {
                s += l[j][k] * u[k][i];
            }
            if (fabs(u[i][i]) < 1e-12) {
                return false;
            }
            l[j][i] = (a[j][i] - s) / u[i][i];
        }
    }
    return true;
}

#endif
