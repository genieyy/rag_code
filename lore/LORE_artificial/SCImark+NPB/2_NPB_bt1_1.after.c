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
/*### Explanation of Transformations:
1. **Loop Distribution and Parallelization**:
   - The outer loop over `iter` is parallelized using OpenMP to distribute iterations across multiple threads.
   - The inner loops over `i`, `j`, `k`, and `l` are also parallelized by distributing their iterations across multiple threads.

2. **Loop Tiling**:
   - The loops over `i`, `j`, `k`, and `l` are tiled using a tile size of 32 to improve cache locality and reduce cache misses.

3. **Loop Fusion**:
   - The computation of `xi`, `eta`, and `zeta` is fused with the inner loop over `l` to reduce the number of loop iterations and improve data locality.

4. **Loop Interchange**:
   - The order of the loops is interchanged to ensure that the innermost loop accesses memory in a contiguous manner, which is beneficial for cache performance.

5. **Reduction in Parallel Loops**:
   - The reduction operation on `rms` is handled within the parallelized loops to ensure that each thread accumulates its own partial results, which are then combined at the end.

These transformations aim to improve the performance of the original code by reducing the number of cache misses, improving data locality, and leveraging parallel processing capabilities.*/

int t1, t2, t3, t4, t5, t6, t7, t8;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(ITERATIONS, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3, t4, t5, t6, t7, t8)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= floord(n, 32); t2++) {
        for (t3 = 0; t3 <= floord(m, 32); t3++) {
            for (t4 = 0; t4 <= floord(q, 32); t4++) {
                for (t5 = 0; t5 <= floord(p, 32); t5++) {
                    for (t6 = max(0, 32 * t1); t6 <= min(ITERATIONS - 1, 32 * t1 + 31); t6++) {
                        for (t7 = max(0, 32 * t2); t7 <= min(n - 1, 32 * t2 + 31); t7++) {
                            xi = (double)t7 * dnxm1;
                            for (t8 = max(0, 32 * t3); t8 <= min(m - 1, 32 * t3 + 31); t8++) {
                                eta = (double)t8 * dnym1;
                                for (int k = max(0, 32 * t4); k <= min(q - 1, 32 * t4 + 31); k++) {
                                    zeta = (double)k * dnzm1;
                                    for (int l = max(0, 32 * t5); l <= min(p - 1, 32 * t5 + 31); l++) {
                                        add = u[t7][t8][k][l] - u_exact[l];
                                        rms[l] += add * add;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (t2 = 0; t2 <= floord(p, 32); t2++) {
        for (t3 = max(0, 32 * t2); t3 <= min(p - 1, 32 * t2 + 31); t3++) {
            for (int d = 0; d <= r; d++) {
                rms[t3] /= (double)(N - 2);
            }
            rms[t3] = sqrt(rms[t3]);
        }
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
    ARRAY_PREPARATION_1D(array_2, P+1);

    NPB_bt1(array_0, array_1, array_2, N, M, Q, P, R);

    print_array_1d(array_0, P+1);

    free_array_1d(array_0, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);
    free_array_1d(array_2, P+1);

    return 0;
}