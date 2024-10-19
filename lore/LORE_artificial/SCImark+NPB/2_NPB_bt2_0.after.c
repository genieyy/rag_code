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

int t1, t2, t3, t4, t5;
register int lbv, ubv;

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (t1 = 0; t1 <= floord(n_ - 3, 32); t1++) {
        for (t2 = 0; t2 <= floord(m_ - 3, 32); t2++) {
            for (t3 = 0; t3 <= floord(q_ - 3, 32); t3++) {
                for (t4 = 0; t4 < p_; t4++) {
                    for (t5 = max(1, 32 * t1); t5 <= min(n_ - 2, 32 * t1 + 31); t5++) {
                        for (int j = max(1, 32 * t2); j <= min(m_ - 2, 32 * t2 + 31); j++) {
                            for (int k = max(1, 32 * t3); k <= min(q_ - 2, 32 * t3 + 31); k++) {
                                double add = rhs[t5][j][k][t4];
                                rms[t4] += add * add;
                            }
                        }
                    }
                }
            }
        }
    }

    for (t1 = 0; t1 < p_; t1++) {
        for (t2 = 0; t2 <= r_; t2++) {
            rms[t1] /= (double)(N - 2);
        }
        rms[t1] = sqrt(rms[t1]);
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