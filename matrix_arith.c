#include "matrix_arith.h"   

#include <stdlib.h>        

// Проверка: обе матрицы существуют и имеют одинаковый размер (w и h)
static int same_size(const matrix *a, const matrix *b) {
    return a && b && matrix_w(a) == matrix_w(b) && matrix_h(a) == matrix_h(b);
}

// m1 += m2 (покомпонентное сложение), размеры должны совпадать
// 0 — успех, -1 — ошибка (NULL или разные размеры)
int matrix_add(matrix *m1, const matrix *m2) {
    if (!same_size(m1, m2)) return -1;
    for (size_t i = 0u; i < matrix_h(m1); ++i) {          // по строкам
        for (size_t j = 0u; j < matrix_w(m1); ++j) {      // по столбцам
            *matrix_ptr(m1, i, j) += *matrix_cptr(m2, i, j); // m1[i][j] += m2[i][j]
        }
    }
    return 0;
}

// m1 -= m2 (покомпонентное вычитание), размеры должны совпадать
int matrix_sub(matrix *m1, const matrix *m2) {
    if (!same_size(m1, m2)) return -1;
    for (size_t i = 0u; i < matrix_h(m1); ++i) {
        for (size_t j = 0u; j < matrix_w(m1); ++j) {
            *matrix_ptr(m1, i, j) -= *matrix_cptr(m2, i, j); // m1[i][j] -= m2[i][j]
        }
    }
    return 0;
}

// m *= d (умножение всех элементов на число)
void matrix_smul(matrix *m, double d) {
    if (!m) return;
    for (size_t i = 0u; i < matrix_h(m); ++i) {
        for (size_t j = 0u; j < matrix_w(m); ++j) {
            *matrix_ptr(m, i, j) *= d;                     // m[i][j] *= d
        }
    }
}

// m /= d (деление всех элементов на число)
// Реализовано через умножение на 1/d
void matrix_sdiv(matrix *m, double d) {
    if (!m) return;
    matrix_smul(m, 1.0 / d);
}

// m = m1 + m2 (создаёт сумму в уже существующую матрицу m)
// Требуется: размеры m, m1, m2 совпадают
int matrix_add2(matrix *m, const matrix *m1, const matrix *m2) {
    if (!same_size(m, m1) || !same_size(m, m2)) return -1;
    for (size_t i = 0u; i < matrix_h(m); ++i) {
        for (size_t j = 0u; j < matrix_w(m); ++j) {
            *matrix_ptr(m, i, j) = *matrix_cptr(m1, i, j) + *matrix_cptr(m2, i, j); // m[i][j] = m1[i][j] + m2[i][j]
        }
    }
    return 0;
}

// m = m1 - m2 (разность в матрицу m), размеры должны совпадать
int matrix_sub2(matrix *m, const matrix *m1, const matrix *m2) {
    if (!same_size(m, m1) || !same_size(m, m2)) return -1;
    for (size_t i = 0u; i < matrix_h(m); ++i) {
        for (size_t j = 0u; j < matrix_w(m); ++j) {
            *matrix_ptr(m, i, j) =
                *matrix_cptr(m1, i, j) - *matrix_cptr(m2, i, j); // m[i][j] = m1[i][j] - m2[i][j]
        }
    }
    return 0;
}

// m = m1 * d (скалярное умножение с записью результата в m)
// Требуется: размеры m и m1 совпадают
int matrix_smul2(matrix *m, const matrix *m1, double d) {
    if (!same_size(m, m1)) return -1;
    for (size_t i = 0u; i < matrix_h(m); ++i) {
        for (size_t j = 0u; j < matrix_w(m); ++j) {
            *matrix_ptr(m, i, j) = *matrix_cptr(m1, i, j) * d; // m[i][j] = m1[i][j] * d
        }
    }
    return 0;
}

// m = m1 / d (деление матрицы на число) — через умножение на 1/d
int matrix_sdiv2(matrix *m, const matrix *m1, double d) {
    return matrix_smul2(m, m1, 1.0 / d);
}

// m = m1 * m2 (матричное умножение)
// Условия размеров:
//   w(m1) == h(m2)
//   m имеет размер: h(m1) × w(m2)
int matrix_mul(matrix *m, const matrix *m1, const matrix *m2) {
    if (!m || !m1 || !m2) return -1;
    if (matrix_w(m1) != matrix_h(m2)) return -1; // проверка совместимости
    if (matrix_w(m) != matrix_w(m2) || matrix_h(m) != matrix_h(m1)) return -1; // размер результата

    // Если результат пишем в одну из входных матриц (m==m1 или m==m2),
    // нельзя затирать данные во время вычисления, значит нужен временный буфер
    const int need_tmp = (m == m1) || (m == m2);
    matrix *out = m;

    matrix *tmp = NULL;
    if (need_tmp) {
        tmp = matrix_alloc_zero(matrix_w(m2), matrix_h(m1)); // временная матрица результата
        if (!tmp) return -1;
        out = tmp;                                          // пишем в tmp
    } else {
        matrix_set_zero(out);                               // обнуляем m перед накоплением суммы
    }

    // Тройной цикл умножения:
    // out[i][j] += m1[i][k] * m2[k][j]
    for (size_t i = 0u; i < matrix_h(m1); ++i) {             // строки m1
        for (size_t k = 0u; k < matrix_w(m1); ++k) {         // общий индекс (столбцы m1 / строки m2)
            const double aik = *matrix_cptr(m1, i, k);       // кешируем m1[i][k]
            for (size_t j = 0u; j < matrix_w(m2); ++j) {     // столбцы m2
                *matrix_ptr(out, i, j) += aik * (*matrix_cptr(m2, k, j));
            }
        }
    }

    // Если был временный буфер — копируем tmp -> m и освобождаем tmp
    if (tmp) {
        const int rc = matrix_assign(m, tmp);
        matrix_free(tmp);
        return rc;
    }
    return 0;
}

// умножение: m1 = m1 * m2
// Делает через временную матрицу, затем копирует обратно.
int matrix_mul_inplace(matrix *m1, const matrix *m2) {
    if (!m1 || !m2) return -1;
    if (matrix_w(m1) != matrix_h(m2)) return -1; // совместимость

    // временная матрица результата (размер: h(m1) × w(m2))
    matrix *tmp = matrix_alloc_zero(matrix_w(m2), matrix_h(m1));
    if (!tmp) return -1;

    // считаем tmp = m1 * m2
    const int rc = matrix_mul(tmp, m1, m2);
    if (rc != 0) {
        matrix_free(tmp);
        return -1;
    }

    if (matrix_w(m1) != matrix_w(tmp) || matrix_h(m1) != matrix_h(tmp)) {
        matrix_free(tmp);
        return -1;
    }

    // m1 = tmp (копирование данных)
    const int rc2 = matrix_assign(m1, tmp);
    matrix_free(tmp);
    return rc2;
}
