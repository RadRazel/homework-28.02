#include "matrix_arith.h"

#include <stdlib.h>

static int same_size(const matrix *a, const matrix *b) {
    return a && b && matrix_w(a) == matrix_w(b) && matrix_h(a) == matrix_h(b);
}

int matrix_add(matrix *m1, const matrix *m2) {
    if (!same_size(m1, m2)) return -1;
    for (size_t i = 0u; i < matrix_h(m1); ++i) {
        for (size_t j = 0u; j < matrix_w(m1); ++j) {
            *matrix_ptr(m1, i, j) += *matrix_cptr(m2, i, j);
        }
    }
    return 0;
}

int matrix_sub(matrix *m1, const matrix *m2) {
    if (!same_size(m1, m2)) return -1;
    for (size_t i = 0u; i < matrix_h(m1); ++i) {
        for (size_t j = 0u; j < matrix_w(m1); ++j) {
            *matrix_ptr(m1, i, j) -= *matrix_cptr(m2, i, j);
        }
    }
    return 0;
}

void matrix_smul(matrix *m, double d) {
    if (!m) return;
    for (size_t i = 0u; i < matrix_h(m); ++i) {
        for (size_t j = 0u; j < matrix_w(m); ++j) {
            *matrix_ptr(m, i, j) *= d;
        }
    }
}

void matrix_sdiv(matrix *m, double d) {
    if (!m) return;
    matrix_smul(m, 1.0 / d);
}

int matrix_add2(matrix *m, const matrix *m1, const matrix *m2) {
    if (!same_size(m, m1) || !same_size(m, m2)) return -1;
    for (size_t i = 0u; i < matrix_h(m); ++i) {
        for (size_t j = 0u; j < matrix_w(m); ++j) {
            *matrix_ptr(m, i, j) = *matrix_cptr(m1, i, j) + *matrix_cptr(m2, i, j);
        }
    }
    return 0;
}

int matrix_sub2(matrix *m, const matrix *m1, const matrix *m2) {
    if (!same_size(m, m1) || !same_size(m, m2)) return -1;
    for (size_t i = 0u; i < matrix_h(m); ++i) {
        for (size_t j = 0u; j < matrix_w(m); ++j) {
            *matrix_ptr(m, i, j) = *matrix_cptr(m1, i, j) - *matrix_cptr(m2, i, j);
        }
    }
    return 0;
}

int matrix_smul2(matrix *m, const matrix *m1, double d) {
    if (!same_size(m, m1)) return -1;
    for (size_t i = 0u; i < matrix_h(m); ++i) {
        for (size_t j = 0u; j < matrix_w(m); ++j) {
            *matrix_ptr(m, i, j) = *matrix_cptr(m1, i, j) * d;
        }
    }
    return 0;
}

int matrix_sdiv2(matrix *m, const matrix *m1, double d) {
    return matrix_smul2(m, m1, 1.0 / d);
}

int matrix_mul(matrix *m, const matrix *m1, const matrix *m2) {
    if (!m || !m1 || !m2) return -1;
    if (matrix_w(m1) != matrix_h(m2)) return -1;
    if (matrix_w(m) != matrix_w(m2) || matrix_h(m) != matrix_h(m1)) return -1;

    const int need_tmp = (m == m1) || (m == m2);
    matrix *out = m;

    matrix *tmp = NULL;
    if (need_tmp) {
        tmp = matrix_alloc_zero(matrix_w(m2), matrix_h(m1));
        if (!tmp) return -1;
        out = tmp;
    } else {
        matrix_set_zero(out);
    }

    for (size_t i = 0u; i < matrix_h(m1); ++i) {
        for (size_t k = 0u; k < matrix_w(m1); ++k) {
            const double aik = *matrix_cptr(m1, i, k);
            for (size_t j = 0u; j < matrix_w(m2); ++j) {
                *matrix_ptr(out, i, j) += aik * (*matrix_cptr(m2, k, j));
            }
        }
    }

    if (tmp) {
        const int rc = matrix_assign(m, tmp);
        matrix_free(tmp);
        return rc;
    }
    return 0;
}

int matrix_mul_inplace(matrix *m1, const matrix *m2) {
    if (!m1 || !m2) return -1;
    if (matrix_w(m1) != matrix_h(m2)) return -1;

    matrix *tmp = matrix_alloc_zero(matrix_w(m2), matrix_h(m1));
    if (!tmp) return -1;

    const int rc = matrix_mul(tmp, m1, m2);
    if (rc != 0) {
        matrix_free(tmp);
        return -1;
    }

    // m1 должен иметь размер tmp? В inplace по смыслу меняется размер результата.
    // Методичка не даёт resize-assign, поэтому делаем строго: m1 должен быть нужного размера.
    // Чтобы не ломать API, требуем: h(m1)=h(tmp) и w(m1)=w(tmp) (то есть m2 квадратная и w(m1)=w(m2)).
    if (matrix_w(m1) != matrix_w(tmp) || matrix_h(m1) != matrix_h(tmp)) {
        matrix_free(tmp);
        return -1;
    }

    const int rc2 = matrix_assign(m1, tmp);
    matrix_free(tmp);
    return rc2;
}
