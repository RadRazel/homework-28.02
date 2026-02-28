#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h> // size_t
#include <stdio.h>  // FILE

// Opaque type (инкапсуляция)
typedef struct matrix matrix;

// Размеры
size_t matrix_w(const matrix *m); // Возвращает число столбцов матрицы (m->w).
size_t matrix_h(const matrix *m); // Возвращает число строк матрицы (m->h).

// Выделение/копирование/освобождение
matrix *matrix_alloc(size_t w, size_t h);             // NULL при ошибке
matrix *matrix_copy(const matrix *m);                 // NULL при ошибке
void matrix_free(matrix *m);

// Быстрый доступ (без проверок диапазона)
double *matrix_ptr(matrix *m, size_t i, size_t j);
const double *matrix_cptr(const matrix *m, size_t i, size_t j);

// безопасные get/set (с проверками)
int matrix_get(const matrix *m, size_t i, size_t j, double *out);
int matrix_set(matrix *m, size_t i, size_t j, double v);

// Инициализация
void matrix_set_zero(matrix *m);
matrix *matrix_set_id(matrix *m);                     // NULL если не квадратная/ошибка

// alloc + set
matrix *matrix_alloc_zero(size_t w, size_t h);        // NULL при ошибке
matrix *matrix_alloc_id(size_t w);                    // NULL при ошибке

// Присваивание без изменения размеров: 0 успех, -1 ошибка
int matrix_assign(matrix *m1, const matrix *m2);

// Ввод/вывод
int matrix_read(matrix *m, FILE *in);                 // 0 успех, -1 ошибка
void matrix_print(const matrix *m, FILE *out);        // печать в читабельном виде

#endif
