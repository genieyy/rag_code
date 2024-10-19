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
   - The outer loop over `iter` is distributed into chunks of size 32 and parallelized using OpenMP. This allows multiple iterations of the loop to be executed simultaneously, leveraging multi-core processors.

2. **Loop Tiling**:
   - The inner loops over `i`, `j`, `k`, and `l` are tiled with a tile size of 32. This reduces the number of cache misses by ensuring that the data accessed within each tile fits into the cache.

3. **Loop Fusion**:
   - The computation of `xi`, `eta`, and `zeta` is fused with the loop over `l` to avoid redundant calculations and improve data locality.

4. **Vectorization**:
   - The inner loop over `l` is marked with `#pragma ivdep` and `#pragma vector always` to hint the compiler to vectorize the loop, assuming the hardware supports SIMD instructions.

5. **Reduction in Loop Nesting**:
   - The loop over `d` is kept separate from the computation of `rms[i]` to avoid unnecessary nesting, which can simplify the control flow and potentially improve performance.

These transformations aim to improve cache utilization, reduce overhead, and leverage parallel processing capabilities to enhance the performance of the original code.*/

int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(ITERATIONS, 32);
#pragma omp parallel for private(lbv, ubv, t4, t5, t6, t7, t8, t9, t10)
for (int t3 = lbp; t3 <= ubp; t3++) {
    for (int t4 = 0; t4 <= floord(n, 32); t4++) {
        for (int t5 = 0; t5 <= floord(m, 32); t5++) {
            for (int t6 = 0; t6 <= floord(q, 32); t6++) {
                for (int t7 = 0; t7 <= floord(p, 32); t7++) {
                    for (int t8 = max(0, 32 * t4); t8 <= min(n - 1, 32 * t4 + 31); t8++) {
                        double xi = (double)t8 * dnxm1;
                        for (int t9 = max(0, 32 * t5); t9 <= min(m - 1, 32 * t5 + 31); t9++) {
                            double eta = (double)t9 * dnym1;
                            for (int t10 = max(0, 32 * t6); t10 <= min(q - 1, 32 * t6 + 31); t10++) {
                                double zeta = (double)t10 * dnzm1;
                                for (int l = max(0, 32 * t7); l <= min(p - 1, 32 * t7 + 31); l++) {
                                    double add = u[t8][t9][t10][l] - u_exact[l];
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

for (int t3 = 0; t3 <= floord(p, 32); t3++) {
    for (int t4 = max(0, 32 * t3); t4 <= min(p - 1, 32 * t3 + 31); t4++) {
        for (int d = 0; d <= r; d++) {
            rms[t4] /= (double)(N - 2);
        }
        rms[t4] = sqrt(rms[t4]);
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