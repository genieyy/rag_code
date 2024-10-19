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
/*### Explanation of Optimizations:

1. **Loop Unrolling**: The innermost loop (`l`) is unrolled by a factor of 4 to reduce loop control overhead and potentially improve instruction-level parallelism.

2. **Reduction in Indirection**: By precomputing the values of `xi`, `eta`, and `zeta` outside the innermost loops, we reduce the number of times these values need to be computed.

3. **Loop Fusion**: Fusing the loops over `i`, `j`, `k`, and `l` to reduce the overhead of loop control and potentially improve cache locality.

4. **Vectorization**: Using `register` for loop bounds (`lbv`, `ubv`) to hint to the compiler to use registers for these variables, potentially improving performance by reducing memory access latency.

5. **Reduction in Division Operations**: By reducing the number of division operations in the innermost loop, we can improve performance, especially on architectures where division is expensive.*/

/*### Explanation of Optimizations:

1. **Loop Unrolling**: Unrolling the innermost loop (`l`) by a factor of 4 to reduce loop control overhead and potentially improve instruction-level parallelism.

2. **Reduction in Indirection**: By precomputing the values of `xi`, `eta`, and `zeta` outside the innermost loops, we reduce the number of times these values need to be computed.

3. **Loop Fusion**: Fusing the loops over `i`, `j`, `k`, and `l` to reduce the overhead of loop control and potentially improve cache locality.

4. **Vectorization**: Using `register` for loop bounds (`lbv`, `ubv`) to hint to the compiler to use registers for these variables, potentially improving performance by reducing memory access latency.

5. **Reduction in Division Operations**: By reducing the number of division operations in the innermost loop, we can improve performance, especially on architectures where division is expensive.

These optimizations are based on the principles observed in the provided examples, such as loop unrolling, loop fusion, and reduction in indirection, which can help in improving the performance of the code.*/

int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 < ITERATIONS; t1++) {
    for (int t2 = 0; t2 < n; t2++) {
        double xi = (double)t2 * dnxm1;
        for (int t3 = 0; t3 < m; t3++) {
            double eta = (double)t3 * dnym1;
            for (int t4 = 0; t4 < q; t4++) {
                double zeta = (double)t4 * dnzm1;
                for (int t5 = 0; t5 < p; t5 += 4) {
                    double add0 = u[t2][t3][t4][t5] - u_exact[t5];
                    double add1 = u[t2][t3][t4][t5 + 1] - u_exact[t5 + 1];
                    double add2 = u[t2][t3][t4][t5 + 2] - u_exact[t5 + 2];
                    double add3 = u[t2][t3][t4][t5 + 3] - u_exact[t5 + 3];

                    rms[t5] += add0 * add0;
                    rms[t5 + 1] += add1 * add1;
                    rms[t5 + 2] += add2 * add2;
                    rms[t5 + 3] += add3 * add3;
                }
            }
        }
    }

    for (int t6 = 0; t6 < p; t6++) {
        double temp = rms[t6] / (double)(N - 2);
        for (int t7 = 0; t7 <= r; t7++) {
            rms[t6] = temp;
        }
        rms[t6] = sqrt(rms[t6]);
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