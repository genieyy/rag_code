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
#define P 50
#define R 60
/* end parameters define */

/* start kernel func */
void NPB_bt2(double *rms, double ****rhs, int n_, int m_, int q_, int p_, int r_) {
    int i, j, k, d, m;
    double add;

    double time_start = omp_get_wtime();
#pragma scop
/**/

for (int iter = 0; iter < ITERATIONS; iter++){
    double temp_rms[p_];
    memset(temp_rms, 0, sizeof(temp_rms));

    for (i = 1; i < n_ - 1; i++) {
        for (j = 1; j < m_ - 1; j++) {
            for (k = 1; k < q_ - 1; k++) {
                for (m = 0; m < p_; m++) {
                    add = rhs[i][j][k][m];
                    temp_rms[m] += add * add;
                }
            }
        }
    }

    for (m = 0; m < p_; m++) {
        temp_rms[m] /= (double)(N - 2);
        rms[m] = sqrt(temp_rms[m]);
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main() {
    ARRAY_PREPARATION_1D(array_0, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt2(array_0, array_1, N, M, Q, P, R);

    print_array_1d(array_0, P+1);

    free_array_1d(array_0, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}