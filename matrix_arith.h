#ifndef MATRIX_ARITH_H
#define MATRIX_ARITH_H

#include "matrix.h"

// m1 += m2, m1 -= m2, m *= d, m /= d
// 0 успех, -1 ошибка (разные размеры / NULL)
int  matrix_add(matrix *m1, const matrix *m2);
int  matrix_sub(matrix *m1, const matrix *m2);
void matrix_smul(matrix *m, double d);
void matrix_sdiv(matrix *m, double d); // через умножение на 1/d

// m = m1 + m2, m = m1 - m2, m = m1 * d, m = m1 / d
int matrix_add2(matrix *m, const matrix *m1, const matrix *m2);
int matrix_sub2(matrix *m, const matrix *m1, const matrix *m2);
int matrix_smul2(matrix *m, const matrix *m1, double d);
int matrix_sdiv2(matrix *m, const matrix *m1, double d);

// m1 *= m2  (через временную)
int matrix_mul_inplace(matrix *m1, const matrix *m2);

// m = m1 * m2  (если m совпадает с m1 или m2 — через временную)
int matrix_mul(matrix *m, const matrix *m1, const matrix *m2);

#endif
