#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include "matrix_utils.hpp"

// Метод Гаусса без выбора ведущего элемента
vec gaussNoPivoting(mat A, vec b);

// Метод Гаусса с частичным выбором главного элемента по столбцу
vec gaussWithPivoting(mat A, vec b);

// Выполнение LU-разложения без перестановок (A = L * U)
void luDecomposition(const mat& A, mat& L, mat& U);

// Решение СЛАУ с использованием уже готовых матриц L и U (прямая и обратная подстановки)
vec luSolve(const mat& L, const mat& U, const vec& b);

#endif // SOLVERS_HPP
