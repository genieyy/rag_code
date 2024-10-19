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
#define N 30
#define M 25
#define Q 20
#define P 15
#define R 10
/* end parameters define */

/* start kernel func */
void NPB_bt1(double *rms, double ****u, double *u_exact, int n, int m, int q, int p, int r)
{
    int i, j, k, l, d;
    double xi, eta, zeta, add;
    int dnxm1 = 1;
    int dnym1 = 2;
    int dnzm1 = 3;

    double time_start = omp_get_wtime();
#pragma scop
/**/

for (int iter = 0; iter < ITERATIONS; iter++){
    double *temp_rms = (double *)calloc(p, sizeof(double)); // Dynamic allocation for temporary array

    for (i = 0; i < n; i++) {
        xi = (double)i * dnxm1;
        for (j = 0; j < m; j++) {
            eta = (double)j * dnym1;
            for (k = 0; k < q; k++) {
                zeta = (double)k * dnzm1;

                for (l = 0; l < p; l++) {
                    add = u[i][j][k][l] - u_exact[l];
                    temp_rms[l] += add * add;
                }
            }
        }
    }

    for (i = 0; i < p; i++) {
        temp_rms[i] /= (double)(N - 2);
        rms[i] = sqrt(temp_rms[i]);
    }

    free(temp_rms); // Free the dynamically allocated memory
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main() {

    ARRAY_PREPARATION_1D(array_0, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);
    ARRAY_PREPARATION_1D(array_2, P+1);

    NPB_bt1(array_0, array_1, array_2, N, M, Q, P, R);

    print_array_1d(array_0, P+1);

    free_array_1d(array_0, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);
    free_array_1d(array_2, P+1);

    return 0;
}