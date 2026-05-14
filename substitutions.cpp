  #ifndef SUBSTITUTIONS_H
#define SUBSTITUTIONS_H

#include <vector>
#include <cmath>

using namespace std;

vector<double> forward_solve(const vector<vector<double>>& l, const vector<double>& b) {
    int n = l.size();
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < i; ++j) {
            s += l[i][j] * y[j];
        }
        y[i] = b[i] - s;
    }
    return y;
}

vector<double> back_solve(const vector<vector<double>>& u, const vector<double>& y) {
    int n = u.size();
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        if (fabs(u[i][i]) < 1e-12) {
            return vector<double>();
        }
        double s = 0.0;
        for (int j = i + 1; j < n; ++j) {
            s += u[i][j] * x[j];
        }
        x[i] = (y[i] - s) / u[i][i];
    }
    return x;
}

#endif
