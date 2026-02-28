#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include "matrix.h"

// Транспонирование (квадратная — на месте), 0 успех, -1 ошибка
int matrix_transpose_inplace(matrix *m);

// Перестановки
int matrix_swap_rows(matrix *m, size_t r1, size_t r2); // 0/-1
int matrix_swap_cols(matrix *m, size_t c1, size_t c2); // 0/-1

// Строковые операции
int matrix_row_mul(matrix *m, size_t r, double k);      // r *= k
int matrix_row_div(matrix *m, size_t r, double k);      // r /= k
int matrix_row_add(matrix *m, size_t r_to, size_t r_from, double k); // r_to += k*r_from

// Норма: max по строкам sum |a_ij|
double matrix_norm(const matrix *m);

#endif
