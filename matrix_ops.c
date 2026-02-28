#include "matrix_ops.h"

#include <math.h> // fabs

int matrix_transpose_inplace(matrix *m) {
    if (!m) return -1;
    if (matrix_w(m) != matrix_h(m)) return -1;

    const size_t n = matrix_w(m);
    for (size_t i = 0u; i < n; ++i) {
        for (size_t j = i + 1u; j < n; ++j) {
            double tmp = *matrix_ptr(m, i, j);
            *matrix_ptr(m, i, j) = *matrix_ptr(m, j, i);
            *matrix_ptr(m, j, i) = tmp;
        }
    }
    return 0;
}

int matrix_swap_rows(matrix *m, size_t r1, size_t r2) {
    if (!m) return -1;
    if (r1 >= matrix_h(m) || r2 >= matrix_h(m)) return -1;
    if (r1 == r2) return 0;

    for (size_t j = 0u; j < matrix_w(m); ++j) {
        double tmp = *matrix_ptr(m, r1, j);
        *matrix_ptr(m, r1, j) = *matrix_ptr(m, r2, j);
        *matrix_ptr(m, r2, j) = tmp;
    }
    return 0;
}

int matrix_swap_cols(matrix *m, size_t c1, size_t c2) {
    if (!m) return -1;
    if (c1 >= matrix_w(m) || c2 >= matrix_w(m)) return -1;
    if (c1 == c2) return 0;

    for (size_t i = 0u; i < matrix_h(m); ++i) {
        double tmp = *matrix_ptr(m, i, c1);
        *matrix_ptr(m, i, c1) = *matrix_ptr(m, i, c2);
        *matrix_ptr(m, i, c2) = tmp;
    }
    return 0;
}

int matrix_row_mul(matrix *m, size_t r, double k) {
    if (!m) return -1;
    if (r >= matrix_h(m)) return -1;
    for (size_t j = 0u; j < matrix_w(m); ++j) {
        *matrix_ptr(m, r, j) *= k;
    }
    return 0;
}

int matrix_row_div(matrix *m, size_t r, double k) {
    return matrix_row_mul(m, r, 1.0 / k);
}

int matrix_row_add(matrix *m, size_t r_to, size_t r_from, double k) {
    if (!m) return -1;
    if (r_to >= matrix_h(m) || r_from >= matrix_h(m)) return -1;
    for (size_t j = 0u; j < matrix_w(m); ++j) {
        *matrix_ptr(m, r_to, j) += k * (*matrix_ptr(m, r_from, j));
    }
    return 0;
}

double matrix_norm(const matrix *m) {
    if (!m) return 0.0;
    double best = 0.0;
    for (size_t i = 0u; i < matrix_h(m); ++i) {
        double s = 0.0;
        for (size_t j = 0u; j < matrix_w(m); ++j) {
            s += fabs(*matrix_cptr(m, i, j));
        }
        if (s > best) best = s;
    }
    return best;
}
