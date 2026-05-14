#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H

#include <vector>
#include <cmath>

using namespace std;

// LU-разложение (L - нижняя треугольная с единицами на диагонали, U - верхняя треугольная)
bool luDecomposition(const vector<vector<double>>& A, 
                      vector<vector<double>>& L, 
                      vector<vector<double>>& U) {
    int n = A.size();
    L.assign(n, vector<double>(n, 0.0));
    U.assign(n, vector<double>(n, 0.0));
    
    for (int i = 0; i < n; ++i) {
        L[i][i] = 1.0;
    }
    
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k) {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = A[i][j] - sum;
        }
        
        for (int j = i + 1; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k) {
                sum += L[j][k] * U[k][i];
            }
            if (fabs(U[i][i]) < 1e-12) {
                return false;
            }
            L[j][i] = (A[j][i] - sum) / U[i][i];
        }
    }
    return true;
}

// Прямая подстановка: L * y = b
vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& b) {
    int n = L.size();
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }
    return y;
}

// Обратная подстановка: U * x = y
vector<double> backwardSubstitution(const vector<vector<double>>& U, const vector<double>& y) {
    int n = U.size();
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        if (fabs(U[i][i]) < 1e-12) {
            return vector<double>();
        }
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }
    return x;
}

// Решение системы через LU-разложение
bool solveLU(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x) {
    vector<vector<double>> L, U;
    if (!luDecomposition(A, L, U)) {
        return false;
    }
    vector<double> y = forwardSubstitution(L, b);
    x = backwardSubstitution(U, y);
    return true;
}

#endif
