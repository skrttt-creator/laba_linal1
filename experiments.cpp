#include "experiments.hpp"
#include "matrix_utils.hpp"
#include "solvers.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace std::chrono;

void runExperiment1_SingleSystem() {
    vector<int> sizes = {100, 200, 500, 1000};
    cout << "\\n=== ЭКСПЕРИМЕНТ 1: Сравнение времени на одной системе ===\\n";
    cout << setw(10) << "Размер (n)" 
         << setw(18) << "Гаусс (без выб)" 
         << setw(18) << "Гаусс (с выб)" 
         << setw(18) << "LU (Общее)" 
         << setw(18) << "LU (Разлож)" 
         << setw(18) << "LU (Подст)" << "\\n";
    
    for (int n : sizes) {
        mat A = generateRandomMatrix(n);
        vec b = generateRandomVector(n);
        
        auto start = high_resolution_clock::now();
        gaussNoPivoting(A, b);
        double t1 = duration<double, milli>(high_resolution_clock::now() - start).count();

        start = high_resolution_clock::now();
        gaussWithPivoting(A, b);
        double t2 = duration<double, milli>(high_resolution_clock::now() - start).count();

        mat L(n, vec(n, 0.0)), U(n, vec(n, 0.0));
        start = high_resolution_clock::now();
        luDecomposition(A, L, U);
        auto mid = high_resolution_clock::now();
        luSolve(L, U, b);
        auto end = high_resolution_clock::now();
        
        double t_lu_d = duration<double, milli>(mid - start).count();
        double t_lu_s = duration<double, milli>(end - mid).count();
        
        cout << setw(10) << n 
             << setw(16) << fixed << setprecision(2) << t1 << " мс"
             << setw(16) << t2 << " мс"
             << setw(16) << (t_lu_d + t_lu_s) << " мс"
             << setw(16) << t_lu_d << " мс"
             << setw(16) << t_lu_s << " мс" << "\\n";
    }
}

void runExperiment2_MultipleRHS() {
    int n = 500;
    mat A = generateRandomMatrix(n);
    vector<int> k_values = {1, 10, 100};
    
    cout << "\\n=== ЭКСПЕРИМЕНТ 2: Несколько правых частей (n = 500) ===\\n";
    cout << setw(15) << "Кол-во систем (k)" 
         << setw(25) << "Метод Гаусса с выб." 
         << setw(30) << "LU-разложение + подст." << "\\n";
         
    for (int k : k_values) {
        vector<vec> B;
        for (int i = 0; i < k; ++i) B.push_back(generateRandomVector(n, 100 + i));

        auto start = high_resolution_clock::now();
        for (int i = 0; i < k; ++i) gaussWithPivoting(A, B[i]);
        double t_gauss = duration<double, milli>(high_resolution_clock::now() - start).count();

        mat L(n, vec(n, 0.0)), U(n, vec(n, 0.0));
        start = high_resolution_clock::now();
        luDecomposition(A, L, U);
        for (int i = 0; i < k; ++i) luSolve(L, U, B[i]);
        double t_lu = duration<double, milli>(high_resolution_clock::now() - start).count();

        cout << setw(15) << k 
             << setw(23) << fixed << setprecision(2) << t_gauss << " мс"
             << setw(28) << t_lu << " мс" << "\\n";
    }
}

void runExperiment3_HilbertAccuracy() {
    vector<int> sizes = {5, 10, 15};
    cout << "\\n=== ЭКСПЕРИМЕНТ 3: Точность на матрицах Гильберта ===\\n";
    
    for (int n : sizes) {
        mat H = generateHilbertMatrix(n);
        vec exact_x(n, 1.0); // Точное решение — вектор из единиц
        vec b = multiplyMatrixVector(H, exact_x);
        double norm_ex = vectorNorm(exact_x);

        vec x_no = gaussNoPivoting(H, b);
        double err_no = vectorNorm(subtractVectors(x_no, exact_x)) / norm_ex;
        double res_no = vectorNorm(subtractVectors(multiplyMatrixVector(H, x_no), b));

        vec x_piv = gaussWithPivoting(H, b);
        double err_piv = vectorNorm(subtractVectors(x_piv, exact_x)) / norm_ex;
        double res_piv = vectorNorm(subtractVectors(multiplyMatrixVector(H, x_piv), b));

        cout << "Размерность n = " << n << ":\\n";
        cout << "  Без выбора:  Погрешность = " << scientific << err_no << ", Невязка = " << res_no << "\\n";
        cout << "  С выбором:   Погрешность = " << scientific << err_piv << ", Невязка = " << res_piv << "\\n";
    }
}
