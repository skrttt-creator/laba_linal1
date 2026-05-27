#ifndef MATRIX_UTILS_HPP
#define MATRIX_UTILS_HPP

#include <vector>

// Удобные псевдонимы для типов данных
typedef std::vector<double> vec;
typedef std::vector<std::vector<double>> mat;

// Генерация случайной плотной матрицы с элементами из [-1.0, 1.0]
mat generateRandomMatrix(int n, int seed = 42);

// Генерация случайного вектора правой части
vec generateRandomVector(int n, int seed = 43);

// Генерация плохо обусловленной матрицы Гильберта
mat generateHilbertMatrix(int n);

// Умножение матрицы на вектор (для вычисления невязки)
vec multiplyMatrixVector(const mat& A, const vec& x);

// Вычисление евклидовой нормы вектора
double vectorNorm(const vec& v);

// Поэлементное вычитание векторов (a - b)
vec subtractVectors(const vec& a, const vec& b);

#endif // MATRIX_UTILS_HPP
