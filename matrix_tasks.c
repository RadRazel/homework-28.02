#include "matrix_tasks.h"

#include "matrix_arith.h"
#include "matrix_ops.h"

#include <math.h>   // fabs

matrix *matrix_exp(const matrix *A, double eps) {
    if (!A) return NULL;
    if (matrix_w(A) != matrix_h(A)) return NULL;
    if (eps <= 0.0) return NULL;

    const size_t n = matrix_w(A);

    matrix *res = matrix_alloc_id(n);
    if (!res) return NULL;

    // term = a1 = A/1!
    matrix *term = matrix_copy(A);
    if (!term) { matrix_free(res); return NULL; }

    // res += term
    if (matrix_add(res, term) != 0) {
        matrix_free(term);
        matrix_free(res);
        return NULL;
    }

    // Итеративно: a_{k+1} = a_k * A / (k+1)
    // Старт: k=1 уже учтён.
    for (size_t k = 1u; k < 100000u; ++k) {
        matrix *tmp = matrix_alloc_zero(n, n);
        if (!tmp) { matrix_free(term); matrix_free(res); return NULL; }

        if (matrix_mul(tmp, term, A) != 0) {
            matrix_free(tmp); matrix_free(term); matrix_free(res);
            return NULL;
        }

        matrix_free(term);
        term = tmp;

        matrix_sdiv(term, (double)(k + 1u));

        if (matrix_add(res, term) != 0) {
            matrix_free(term); matrix_free(res);
            return NULL;
        }

        if (matrix_norm(term) < eps) break;
    }

    matrix_free(term);
    return res;
}

int matrix_solve_gauss(matrix *X, const matrix *A, const matrix *B, double eps) {
    if (!X || !A || !B) return -1;
    if (eps <= 0.0) return -1;

    const size_t n = matrix_h(A);
    if (matrix_w(A) != n) return -1;         // A квадратная
    if (matrix_h(B) != n) return -1;         // совместимость
    if (matrix_h(X) != n || matrix_w(X) != matrix_w(B)) return -1;

    const size_t m = matrix_w(B);            // число правых частей

    // Работаем на копиях, чтобы не портить вход
    matrix *M = matrix_copy(A);
    matrix *R = matrix_copy(B);
    if (!M || !R) {
        matrix_free(M);
        matrix_free(R);
        return -1;
    }

    // Прямой ход с частичным выбором главного элемента по столбцу
    for (size_t i = 0u; i < n; ++i) {
        // найти pivot row
        size_t piv = i;
        double best = fabs(*matrix_cptr(M, i, i));
        for (size_t r = i + 1u; r < n; ++r) {
            double v = fabs(*matrix_cptr(M, r, i));
            if (v > best) {
                best = v;
                piv = r;
            }
        }

        if (best < eps) { // вырожденная/почти
            matrix_free(M);
            matrix_free(R);
            return -1;
        }

        if (piv != i) {
            (void)matrix_swap_rows(M, piv, i);
            (void)matrix_swap_rows(R, piv, i);
        }

        // зануление ниже диагонали
        const double aii = *matrix_cptr(M, i, i);
        for (size_t r = i + 1u; r < n; ++r) {
            const double ari = *matrix_cptr(M, r, i);
            const double factor = ari / aii;

            // row_r -= factor * row_i  (для всех столбцов матрицы M)
            for (size_t c = i; c < n; ++c) {
                *matrix_ptr(M, r, c) -= factor * (*matrix_cptr(M, i, c));
            }

            // то же для правых частей R
            for (size_t c = 0u; c < m; ++c) {
                *matrix_ptr(R, r, c) -= factor * (*matrix_cptr(R, i, c));
            }
        }
    }

    // Обратный ход
    matrix_set_zero(X);

    for (size_t ii = 0u; ii < n; ++ii) {
        const size_t i = n - 1u - ii;

        const double aii = *matrix_cptr(M, i, i);
        if (fabs(aii) < eps) {
            matrix_free(M);
            matrix_free(R);
            return -1;
        }

        for (size_t c = 0u; c < m; ++c) {
            double s = *matrix_cptr(R, i, c);
            for (size_t j = i + 1u; j < n; ++j) {
                s -= (*matrix_cptr(M, i, j)) * (*matrix_cptr(X, j, c));
            }
            *matrix_ptr(X, i, c) = s / aii;
        }
    }

    matrix_free(M);
    matrix_free(R);
    return 0;
}
