#include "matrix.h"
#include "matrix_tasks.h"
#include "matrix_arith.h"
#include "matrix_ops.h"

#include <stdio.h>

static void fill_diag(matrix *m, double a, double b, double c) {
    matrix_set_zero(m);
    *matrix_ptr(m, 0, 0) = a;
    *matrix_ptr(m, 1, 1) = b;
    *matrix_ptr(m, 2, 2) = c;
}

int main(void) {
    // === exp(A) test ===
    matrix *A = matrix_alloc(3, 3);
    if (!A) return 1;

    fill_diag(A, 1.0, 2.0, -1.0);

    matrix *E = matrix_exp(A, 1e-12);
    if (!E) {
        matrix_free(A);
        return 1;
    }

    printf("A:\n");
    matrix_print(A, stdout);
    printf("\nexp(A):\n");
    matrix_print(E, stdout);

    matrix_free(E);
    matrix_free(A);

    // === Gauss test: A x = b (n=3, m=1) ===
    matrix *G = matrix_alloc(3, 3);
    matrix *B = matrix_alloc(1, 3); // B: 3x1 => w=1, h=3
    matrix *X = matrix_alloc(1, 3);
    if (!G || !B || !X) {
        matrix_free(G); matrix_free(B); matrix_free(X);
        return 1;
    }

    // Пример системы:
    // [2 1 -1] [x]   [ 8]
    // [-3 -1 2] [y] =[ -11]
    // [-2 1 2] [z]   [ -3]
    *matrix_ptr(G,0,0)=2;  *matrix_ptr(G,0,1)=1;  *matrix_ptr(G,0,2)=-1;
    *matrix_ptr(G,1,0)=-3; *matrix_ptr(G,1,1)=-1; *matrix_ptr(G,1,2)=2;
    *matrix_ptr(G,2,0)=-2; *matrix_ptr(G,2,1)=1;  *matrix_ptr(G,2,2)=2;

    *matrix_ptr(B,0,0)=8;
    *matrix_ptr(B,1,0)=-11;
    *matrix_ptr(B,2,0)=-3;

    if (matrix_solve_gauss(X, G, B, 1e-12) != 0) {
        fprintf(stderr, "Gauss failed\n");
        matrix_free(G); matrix_free(B); matrix_free(X);
        return 1;
    }

    printf("\nSolution X:\n");
    matrix_print(X, stdout);

    // Проверка: ||G*X - B||
    matrix *GX = matrix_alloc_zero(1, 3);
    matrix_mul(GX, G, X);
    matrix_sub(GX, B);
    printf("\n||G*X - B|| = %.12g\n", matrix_norm(GX));

    matrix_free(GX);
    matrix_free(G); matrix_free(B); matrix_free(X);

    return 0;
}
