#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <iomanip>

using namespace std;
using namespace std::chrono;

typedef vector<double> vec;
typedef vector<vector<double>> mat;

mat randommatix(int n, int seed = 42) {
    mt19937 gen(seed);
    uniform_real_distribution<> dis(-1.0, 1.0);
    mat A(n, vec(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = dis(gen);
    return A;
}

vec randomvector(int n, int seed = 42 + 1) {
    mt19937 gen(seed);
    uniform_real_distribution<> dis(-1.0, 1.0);
    vec b(n);
    for (int i = 0; i < n; ++i)
        b[i] = dis(gen);
    return b;
}

mat genhilbermatrix(int n) {
    mat H(n, vec(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            H[i][j] = 1.0 / (i + j + 1); 
    return H;
}

vec matXvec(const mat& A, const vec& x) {
    int n = A.size();
    vec b(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            b[i] += A[i][j] * x[j];
    return b;
}

double vectorNorm(const vec& v) {
    double sum = 0.0;
    for (double val : v) sum += val * val;
    return sqrt(sum);
}

vec subtractVectors(const vec& a, const vec& b) {
    vec res(a.size());
    for (size_t i = 0; i < a.size(); ++i) res[i] = a[i] - b[i];
    return res;
}



vec gaussNoPivoting(mat A, vec b) {
    int n = A.size();
    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
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

vec gausswithpivot(mat A, vec b) {
    int n = A.size();
    for (int k = 0; k < n - 1; ++k) {
        int maxRow = k;
        double maxEl = abs(A[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (abs(A[i][k]) > maxEl) {
                maxEl = abs(A[i][k]);
                maxRow = i;
            }
        }
        swap(A[k], A[maxRow]);
        swap(b[k], b[maxRow]);
        for (int i = k + 1; i < n; ++i) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
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

void LUdecomp(const mat& A, mat& L, mat& U) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) sum += (L[i][j] * U[j][k]);
            U[i][k] = A[i][k] - sum;
        }
        for (int k = i; k < n; k++) {
            if (i == k)
                L[i][i] = 1.0; 
            else {
                double sum = 0.0;
                for (int j = 0; j < i; j++) sum += (L[k][j] * U[j][i]);
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

vec slau(const mat& L, const mat& U, const vec& b) {
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


//эксперементы

//сравнение времени
void experement1() {
    cout << "\n сравнение времени решения одной системы: " << endl;
    vector<int> sizes = {100, 200, 500, 1000};
    
    cout << left << setw(10) << "Размер(n)" 
         << setw(15) << "Гаусс(без)" 
         << setw(15) << "Гаусс(с)" 
         << setw(15) << "LU(общее)" 
         << setw(15) << "LU(разл)" 
         << setw(15) << "LU(подст)" << endl;

    for (int n : sizes) {
        mat A = randommatix(n);
        vec b = randomvector(n);
        
        auto start = high_resolution_clock::now();
        gaussNoPivoting(A, b);
        auto end = high_resolution_clock::now();
        double t_gauss_no = duration<double, milli>(end - start).count();

        start = high_resolution_clock::now();
        gausswithpivot(A, b);
        end = high_resolution_clock::now();
        double t_gauss_pivot = duration<double, milli>(end - start).count();

        mat L(n, vec(n, 0.0)), U(n, vec(n, 0.0));
        start = high_resolution_clock::now();
        LUdecomp(A, L, U);
        auto t_mid = high_resolution_clock::now();
        slau(L, U, b);
        end = high_resolution_clock::now();
        
        double t_lu_decomp = duration<double, milli>(t_mid - start).count();
        double t_lu_solve = duration<double, milli>(end - t_mid).count();
        double t_lu_total = t_lu_decomp + t_lu_solve;

        cout << left << setw(10) << n 
             << setw(15) << t_gauss_no 
             << setw(15) << t_gauss_pivot 
             << setw(15) << t_lu_total 
             << setw(15) << t_lu_decomp 
             << setw(15) << t_lu_solve << endl;
    }
}

//экономия времени
void experement2() {
    cout << "\n экономия времени при множественных правых частях:" << endl;
    int n = 500;
    mat A = randommatix(n);
    vector<int> k_values = {1, 10, 100};
    
    cout << "размер матрицы: " << n << "x" << n << endl;
    cout << left << setw(10) << "Кол-во b(k)" << setw(20) << "Гаусс с выбором (мс)" << setw(20) << "LU-разложение (мс)" << endl;

    for (int k : k_values) {
        vector<vec> B;
        for (int i = 0; i < k; ++i) B.push_back(randomvector(n, i));

        auto start = high_resolution_clock::now();
        for (int i = 0; i < k; ++i) {
            gausswithpivot(A, B[i]);
        }
        auto end = high_resolution_clock::now();
        double t_gauss = duration<double, milli>(end - start).count();

        mat L(n, vec(n, 0.0)), U(n, vec(n, 0.0));
        start = high_resolution_clock::now();
        LUdecomp(A, L, U);
        for (int i = 0; i < k; ++i) {
            slau(L, U, B[i]);
        }
        end = high_resolution_clock::now();
        double t_lu = duration<double, milli>(end - start).count();

        cout << left << setw(10) << k << setw(20) << t_gauss << setw(20) << t_lu << endl;
    }
}
//проверка точности
void experement3() {
    cout << "\n проверка точности на плохо обусловленных матрицах:" << endl;
    vector<int> sizes = {5, 10, 15};
    
    for (int n : sizes) {
        cout << "\n матрица Гильберта размерности n = " << n << endl;
        mat H = genhilbermatrix(n);
        vec exact_x(n, 1.0);
        vec b = matXvec(H, exact_x);

        double norm_exact_x = vectorNorm(exact_x);

        vec x_no = gaussNoPivoting(H, b);
        double err_no = vectorNorm(subtractVectors(x_no, exact_x)) / norm_exact_x;
        double res_no = vectorNorm(subtractVectors(matXvec(H, x_no), b));

        vec x_piv = gausswithpivot(H, b);
        double err_piv = vectorNorm(subtractVectors(x_piv, exact_x)) / norm_exact_x;
        double res_piv = vectorNorm(subtractVectors(matXvec(H, x_piv), b));

        cout << "Метод                   | Относит. погрешность | Невязка" << endl;
        cout << "Гаусс (без выбора)      | " << scientific << err_no << "         | " << res_no << endl;
        cout << "Гаусс (с выбором)       | " << err_piv << "         | " << res_piv << endl;
    }
}

int main() {
    cout << "таблицы эксперементов:" << endl;
    
    experement1();
    experement2();
    experement3();

    return 0;
}
