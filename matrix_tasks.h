#ifndef MATRIX_TASKS_H
#define MATRIX_TASKS_H

#include "matrix.h"

// exp(A) по ряду: E + A + A^2/2! + ...
// Условие остановки: norm(ai) < eps, где ai = A^i/i!
// NULL при ошибке
matrix *matrix_exp(const matrix *A, double eps);

// Решение AX=B методом Гаусса с выбором главного элемента по столбцу.
// A: n×n, B: n×m. Выход X: n×m (должна быть заранее выделена правильного размера).
// 0 успех, -1 ошибка (размеры/память/вырожденность)
int matrix_solve_gauss(matrix *X, const matrix *A, const matrix *B, double eps);

#endif
