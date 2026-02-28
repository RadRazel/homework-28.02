#include "matrix.h"     

#include <stdlib.h>     
#include <string.h>     // memcpy — копирование блоков памяти
#include <stdint.h>     // SIZE_MAX — максимум для size_t (для проверки переполнения)

// Внутренняя структура матрицы 
struct matrix {
    double *data;       // данные матрицы в одном массиве: data[w*i + j]
    size_t w, h;        // w = число столбцов, h = число строк
};

// Возвращает число столбцов (если m == NULL, возвращает 0)
size_t matrix_w(const matrix *m) { return m ? m->w : 0u; }

// Возвращает число строк (если m == NULL, возвращает 0)
size_t matrix_h(const matrix *m) { return m ? m->h : 0u; }

// Проверка переполнения при умножении a*b для size_t
// out = a*b, если всё ок; возвращает 1 при переполнении, 0 если ок
static int mul_overflow_size_t(size_t a, size_t b, size_t *out) {
    if (a != 0u && b > SIZE_MAX / a) return 1;
    *out = a * b;
    return 0;
}

// Выделяет матрицу w×h (без инициализации значений)
// Возвращает NULL при ошибке или если w/h = 0
matrix *matrix_alloc(size_t w, size_t h) {
    if (w == 0u || h == 0u) return NULL;

    // считаем количество элементов и проверяем переполнение
    size_t n = 0u;
    if (mul_overflow_size_t(w, h, &n)) return NULL;
    if (n > SIZE_MAX / sizeof(double)) return NULL;

    // выделяем сам объект matrix
    matrix *m = (matrix *)malloc(sizeof(*m));
    if (!m) return NULL;

    // запоминаем размеры и выделяем память под элементы
    m->w = w;
    m->h = h;
    m->data = (double *)malloc(n * sizeof(double));
    if (!m->data) {
        free(m);
        return NULL;
    }
    return m;
}

// Освобождение матрицы 
void matrix_free(matrix *m) {
    if (!m) return;
    free(m->data);  // сначала данные
    free(m);        // потом структура
}

// Быстрый доступ к элементу (i, j) без проверок границ
double *matrix_ptr(matrix *m, size_t i, size_t j) {
    return &m->data[m->w * i + j];
}

// То же самое, но для const-матрицы (только чтение)
const double *matrix_cptr(const matrix *m, size_t i, size_t j) {
    return &m->data[m->w * i + j];
}

// Безопасное чтение элемента: проверяет указатели и границы
// 0 — успех, -1 — ошибка
int matrix_get(const matrix *m, size_t i, size_t j, double *out) {
    if (!m || !out) return -1;
    if (i >= m->h || j >= m->w) return -1;
    *out = *matrix_cptr(m, i, j);
    return 0;
}

// Безопасная запись элемента: проверяет матрицу и границы
// 0 — успех, -1 — ошибка
int matrix_set(matrix *m, size_t i, size_t j, double v) {
    if (!m) return -1;
    if (i >= m->h || j >= m->w) return -1;
    *matrix_ptr(m, i, j) = v;
    return 0;
}

// Заполняет все элементы матрицы нулями
void matrix_set_zero(matrix *m) {
    if (!m) return;
    const size_t n = m->w * m->h;
    for (size_t k = 0u; k < n; ++k) m->data[k] = 0.0;
}

// Делает матрицу единичной (только если квадратная)
// Возвращает m при успехе, NULL если не квадратная/ошибка
matrix *matrix_set_id(matrix *m) {
    if (!m) return NULL;
    if (m->w != m->h) return NULL;

    matrix_set_zero(m);                 // сначала всё в ноль
    for (size_t i = 0u; i < m->w; ++i)  // по диагонали ставим 1
        *matrix_ptr(m, i, i) = 1.0;

    return m;
}

// matrix_alloc + заполнение нулями
matrix *matrix_alloc_zero(size_t w, size_t h) {
    matrix *m = matrix_alloc(w, h);
    if (!m) return NULL;
    matrix_set_zero(m);
    return m;
}

// matrix_alloc(w,w) + единичная
matrix *matrix_alloc_id(size_t w) {
    matrix *m = matrix_alloc(w, w);
    if (!m) return NULL;
    if (!matrix_set_id(m)) {    // если вдруг не получилось
        matrix_free(m);
        return NULL;
    }
    return m;
}

// Присваивает m1 = m2 (только если размеры совпадают)
// 0 — успех, -1 — ошибка
int matrix_assign(matrix *m1, const matrix *m2) {
    if (!m1 || !m2) return -1;
    if (m1->w != m2->w || m1->h != m2->h) return -1;

    // копируем все элементы как блок памяти
    const size_t n = m1->w * m1->h;
    memcpy(m1->data, m2->data, n * sizeof(double));
    return 0;
}

// Создаёт полную копию матрицы
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

// Читает элементы матрицы из потока FILE* 
// Ожидает m->h * m->w чисел типа double
// 0 — успех, -1 — ошибка чтения
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

// Печатает матрицу в поток FILE* (stdout/stderr)
void matrix_print(const matrix *m, FILE *out) {
    if (!out) return;
    if (!m) {
        fprintf(out, "(null)\n");
        return;
    }

    // печать построчно, с фиксированной шириной и точностью
    for (size_t i = 0u; i < m->h; ++i) {
        for (size_t j = 0u; j < m->w; ++j) {
            fprintf(out, "%12.6f ", *matrix_cptr(m, i, j));
        }
        fprintf(out, "\n");
    }
}
