#include "matrix_tasks.h"   

#include "matrix_arith.h"  
#include "matrix_ops.h"    

#include <math.h>         

// Экспонента матрицы exp(A) через ряд Тейлора:
// exp(A) = I + A/1! + A^2/2! + A^3/3! + ...
// eps — порог остановки: когда норма очередного члена < eps, прекращаем
matrix *matrix_exp(const matrix *A, double eps) {
    if (!A) return NULL;                         // нет входной матрицы
    if (matrix_w(A) != matrix_h(A)) return NULL;  // exp определяем для квадратной
    if (eps <= 0.0) return NULL;                 // порог должен быть положительным

    const size_t n = matrix_w(A);                // размер n×n

    // res = I (начинаем сумму ряда с единичной матрицы)
    matrix *res = matrix_alloc_id(n);
    if (!res) return NULL;

    // term = A/1! (первый добавляемый член ряда)
    matrix *term = matrix_copy(A);
    if (!term) { matrix_free(res); return NULL; }

    // res += term (добавляем A)
    if (matrix_add(res, term) != 0) {
        matrix_free(term);
        matrix_free(res);
        return NULL;
    }

    // Итеративно считаем следующие члены ряда:
    // a_{k+1} = a_k * A / (k+1)
    // Здесь term = a_k, стартуем с k=1 (т.к. A/1! уже добавили).
    for (size_t k = 1u; k < 100000u; ++k) {
        // tmp будет хранить term*A (то есть следующий числитель)
        matrix *tmp = matrix_alloc_zero(n, n);
        if (!tmp) { matrix_free(term); matrix_free(res); return NULL; }

        // tmp = term * A
        if (matrix_mul(tmp, term, A) != 0) {
            matrix_free(tmp); matrix_free(term); matrix_free(res);
            return NULL;
        }

        // term больше не нужен, заменяем его на tmp (новый член)
        matrix_free(term);
        term = tmp;

        // делим на (k+1): теперь term = A^(k+1)/(k+1)!
        matrix_sdiv(term, (double)(k + 1u));

        // прибавляем очередной член ряда к сумме
        if (matrix_add(res, term) != 0) {
            matrix_free(term); matrix_free(res);
            return NULL;
        }

        // критерий остановки: если вклад стал маленьким — выходим
        if (matrix_norm(term) < eps) break;
    }

    // term — временная матрица, освобождаем
    matrix_free(term);
    return res;   // res = exp(A)
}

// Решение системы A*X = B методом Гаусса с частичным выбором главного элемента
// A — n×n, B — n×m (m правых частей), X — n×m (результат)
// eps — порог для проверки "нулевого" pivot (вырожденность/почти ноль)
int matrix_solve_gauss(matrix *X, const matrix *A, const matrix *B, double eps) {
    if (!X || !A || !B) return -1;    // проверка указателей
    if (eps <= 0.0) return -1;  // eps должен быть > 0

    const size_t n = matrix_h(A);      // n = число строк A
    if (matrix_w(A) != n) return -1;     // A должна быть квадратной
    if (matrix_h(B) != n) return -1;     // B совместима по строкам
    if (matrix_h(X) != n || matrix_w(X) != matrix_w(B)) return -1; // X нужного размера

    const size_t m = matrix_w(B); // m = число столбцов B (число правых частей)

    // Работаем на копиях, чтобы не портить входные A и B
    // M ~ A (будет превращаться в верхнетреугольную)
    // R ~ B (будет меняться так же, как строки M)
    matrix *M = matrix_copy(A);
    matrix *R = matrix_copy(B);
    if (!M || !R) {
        matrix_free(M);
        matrix_free(R);
        return -1;
    }

    // Прямой ход (приведение к верхнетреугольному виду) 
    // Частичный выбор главного элемента по столбцу: ищем максимальный |M[r][i]| ниже/на диагонали.
    for (size_t i = 0u; i < n; ++i) {
        // 1) найти pivot row (строку с максимальным элементом в текущем столбце i)
        size_t piv = i;
        double best = fabs(*matrix_cptr(M, i, i));
        for (size_t r = i + 1u; r < n; ++r) {
            double v = fabs(*matrix_cptr(M, r, i));
            if (v > best) {
                best = v;
                piv = r;
            }
        }

        // если pivot слишком маленький — матрица вырождена/почти вырождена
        if (best < eps) {
            matrix_free(M);
            matrix_free(R);
            return -1;
        }

        // 2) если pivot не на диагонали — меняем строки местами (и в M, и в R)
        if (piv != i) {
            (void)matrix_swap_rows(M, piv, i);
            (void)matrix_swap_rows(R, piv, i);
        }

        // 3) зануление элементов ниже диагонали в столбце i
        const double aii = *matrix_cptr(M, i, i);     // ведущий элемент
        for (size_t r = i + 1u; r < n; ++r) {         // для каждой строки ниже
            const double ari = *matrix_cptr(M, r, i); // элемент под диагональю
            const double factor = ari / aii;          // коэффициент, чтобы занулить M[r][i]

            // row_r -= factor * row_i  (для M, начиная с текущего столбца i)
            for (size_t c = i; c < n; ++c) {
                *matrix_ptr(M, r, c) -= factor * (*matrix_cptr(M, i, c));
            }

            // то же самое делаем с правыми частями (R)
            for (size_t c = 0u; c < m; ++c) {
                *matrix_ptr(R, r, c) -= factor * (*matrix_cptr(R, i, c));
            }
        }
    }

    // Обратный ход (back substitution) 
    // Решаем верхнетреугольную систему: начиная с последней строки.
    matrix_set_zero(X);  // на всякий случай обнулим X

    for (size_t ii = 0u; ii < n; ++ii) {
        const size_t i = n - 1u - ii;  // идём i = n-1, n-2, ..., 0

        const double aii = *matrix_cptr(M, i, i);
        if (fabs(aii) < eps) {   // защита от деления на почти ноль
            matrix_free(M);
            matrix_free(R);
            return -1;
        }

        // решаем сразу для всех правых частей (для каждого столбца B)
        for (size_t c = 0u; c < m; ++c) {
            // s = R[i][c] - sum_{j=i+1..n-1} M[i][j]*X[j][c]
            double s = *matrix_cptr(R, i, c);
            for (size_t j = i + 1u; j < n; ++j) {
                s -= (*matrix_cptr(M, i, j)) * (*matrix_cptr(X, j, c));
            }
            // X[i][c] = s / M[i][i]
            *matrix_ptr(X, i, c) = s / aii;
        }
    }

    // очищаем временные копии
    matrix_free(M);
    matrix_free(R);
    return 0;                                     
}
