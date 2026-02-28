#include "matrix.h"

#include <stdlib.h>   // malloc, free
#include <string.h>   // memcpy
#include <stdint.h>   // SIZE_MAX

struct matrix {
    double *data;   // data[w*i + j]
    size_t w, h;    // w=cols, h=rows
};

size_t matrix_w(const matrix *m) { return m ? m->w : 0u; } //возвращает столбцы
size_t matrix_h(const matrix *m) { return m ? m->h : 0u; } //возвращает строки

// провёерка на переполнение size_t
static int mul_overflow_size_t(size_t a, size_t b, size_t *out) {
    if (a != 0u && b > SIZE_MAX / a) return 1;
    *out = a * b;
    return 0;
}

matrix *matrix_alloc(size_t w, size_t h) {
    if (w == 0u || h == 0u) return NULL;

    size_t n = 0u;
    if (mul_overflow_size_t(w, h, &n)) return NULL;
    if (n > SIZE_MAX / sizeof(double)) return NULL;

    matrix *m = (matrix *)malloc(sizeof(*m));
    if (!m) return NULL;

    m->w = w;
    m->h = h;
    m->data = (double *)malloc(n * sizeof(double));
    if (!m->data) {
        free(m);
        return NULL;
    }
    return m;
}

void matrix_free(matrix *m) {
    if (!m) return;
    free(m->data);
    free(m);
}

double *matrix_ptr(matrix *m, size_t i, size_t j) {
    return &m->data[m->w * i + j];
}

const double *matrix_cptr(const matrix *m, size_t i, size_t j) {
    return &m->data[m->w * i + j];
}

int matrix_get(const matrix *m, size_t i, size_t j, double *out) {
    if (!m || !out) return -1;
    if (i >= m->h || j >= m->w) return -1;
    *out = *matrix_cptr(m, i, j);
    return 0;
}

int matrix_set(matrix *m, size_t i, size_t j, double v) {
    if (!m) return -1;
    if (i >= m->h || j >= m->w) return -1;
    *matrix_ptr(m, i, j) = v;
    return 0;
}

void matrix_set_zero(matrix *m) {
    if (!m) return;
    const size_t n = m->w * m->h;
    for (size_t k = 0u; k < n; ++k) m->data[k] = 0.0;
}

matrix *matrix_set_id(matrix *m) {
    if (!m) return NULL;
    if (m->w != m->h) return NULL;

    matrix_set_zero(m);
    for (size_t i = 0u; i < m->w; ++i) {
        *matrix_ptr(m, i, i) = 1.0;
    }
    return m;
}

matrix *matrix_alloc_zero(size_t w, size_t h) {
    matrix *m = matrix_alloc(w, h);
    if (!m) return NULL;
    matrix_set_zero(m);
    return m;
}

matrix *matrix_alloc_id(size_t w) {
    matrix *m = matrix_alloc(w, w);
    if (!m) return NULL;
    if (!matrix_set_id(m)) {
        matrix_free(m);
        return NULL;
    }
    return m;
}

int matrix_assign(matrix *m1, const matrix *m2) {
    if (!m1 || !m2) return -1;
    if (m1->w != m2->w || m1->h != m2->h) return -1;

    const size_t n = m1->w * m1->h;
    memcpy(m1->data, m2->data, n * sizeof(double));
    return 0;
}

matrix *matrix_copy(const matrix *m) {
    if (!m) return NULL;
    matrix *c = matrix_alloc(m->w, m->h);
    if (!c) return NULL;
    if (matrix_assign(c, m) != 0) {
        matrix_free(c);
        return NULL;
    }
    return c;
}

int matrix_read(matrix *m, FILE *in) {
    if (!m || !in) return -1;
    for (size_t i = 0u; i < m->h; ++i) {
        for (size_t j = 0u; j < m->w; ++j) {
            double x = 0.0;
            if (fscanf(in, "%lf", &x) != 1) return -1;
            *matrix_ptr(m, i, j) = x;
        }
    }
    return 0;
}

void matrix_print(const matrix *m, FILE *out) {
    if (!out) return;
    if (!m) {
        fprintf(out, "(null)\n");
        return;
    }
    for (size_t i = 0u; i < m->h; ++i) {
        for (size_t j = 0u; j < m->w; ++j) {
            fprintf(out, "%12.6f ", *matrix_cptr(m, i, j));
        }
        fprintf(out, "\n");
    }
}
