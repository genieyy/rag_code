#include <stdio.h>  
#include <stdlib.h>  
#include <string.h>
#include <math.h>
#include "lore.h"

#include <omp.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */ 
#define N 20
#define M 30
#define Q 40
#define P 6
/* end parameters define */

/* start kernel func */
void NPB_bt(double ****u, double ****rhs, int n, int m, int q, int p) {
    int i, j, k, l;
    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 1; i < n - 1; i++) {
        for (int j = 1; j < m - 1; j++) {
            for (int k = 1; k < q - 1; k++) {
                for (int l = 0; l < p - 1; l += 4) {
                    u[i][j][k][l] = u[i][j][k][l] + rhs[i][j][k][l];
                    u[i][j][k][l + 1] = u[i][j][k][l + 1] + rhs[i][j][k][l + 1];
                    u[i][j][k][l + 2] = u[i][j][k][l + 2] + rhs[i][j][k][l + 2];
                    if (l + 3 < p - 1) {
                        u[i][j][k][l + 3] = u[i][j][k][l + 3] + rhs[i][j][k][l + 3];
                    }
                }
            }
        }
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main() {
    ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, Q + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, Q + 1, P + 1);

    NPB_bt(array_0, array_1, N, M, Q, P);
    
    print_array_4d(array_0, N, M, Q, P);

    free_array_4d(array_0, N + 1, M + 1, Q + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, Q + 1, P + 1);

    return 0;
}