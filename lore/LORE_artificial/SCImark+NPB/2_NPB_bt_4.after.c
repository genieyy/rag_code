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
    for (int i = 1; i < n - 1; i += 4) {
        for (int j = 1; j < m - 1; j += 4) {
            for (int k = 1; k < q - 1; k += 4) {
                for (int l = 0; l < p - 1; l += 4) {
                    for (int ii = 0; ii < 4 && i + ii < n - 1; ii++) {
                        for (int jj = 0; jj < 4 && j + jj < m - 1; jj++) {
                            for (int kk = 0; kk < 4 && k + kk < q - 1; kk++) {
                                for (int ll = 0; ll < 4 && l + ll < p - 1; ll++) {
                                    u[i + ii][j + jj][k + kk][l + ll] = u[i + ii][j + jj][k + kk][l + ll] + rhs[i + ii][j + jj][k + kk][l + ll];
                                }
                            }
                        }
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