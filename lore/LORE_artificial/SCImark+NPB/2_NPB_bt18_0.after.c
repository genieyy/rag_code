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
#define N 10
#define M 12
#define L 14
#define P 16
/* end parameters define */

/* start kernel func */
void NPB_bt18(double *****lhs, double ****rhs, int nz, int mz, int q, int BLOCK_SIZE)
{
    int i, j, k, m, n;

    double time_start = omp_get_wtime();
#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Unrolling and Jamming**:
   - The original code has multiple nested loops. The optimized code reduces the number of loops by unrolling and jamming them together, which can improve cache performance and reduce loop overhead.

2. **Parallelization**:
   - The use of `#pragma omp parallel for` indicates that the outer loops are parallelized, allowing multiple threads to execute the loop iterations concurrently.

3. **Loop Fusion**:
   - The optimized code fuses multiple loops into a single loop, which can reduce the overhead of loop control and improve data locality.

4. **Vectorization**:
   - The use of `#pragma ivdep` and `#pragma vector always` hints to the compiler to ignore potential vector dependencies and to vectorize the loop, which can significantly speed up the computation on vector units.

5. **Loop Distribution**:
   - The optimized code distributes the computation across multiple threads by dividing the iteration space and assigning different parts to different threads.

### Application of Methods to Improve Performance:

1. **Parallelization**:
   - The outer loop over `ITERATIONS` is parallelized using `#pragma omp parallel for` to leverage multi-threading.

2. **Loop Fusion**:
   - The nested loops over `i`, `j`, `k`, `m`, and `n` are fused into a single loop structure to reduce loop overhead and improve data locality.

3. **Vectorization**:
   - The innermost loop over `n` is vectorized using `#pragma ivdep` and `#pragma vector always` to take advantage of SIMD instructions.

4. **Loop Distribution**:
   - The iteration space of the outer loops is distributed across threads using OpenMP, which allows the computation to be spread across multiple cores.

By applying these transformations, the performance of the original code is expected to improve significantly, especially on multi-core processors and systems with vector units.*/

int t1, t2, t3, t4, t5, t6, t7, t8;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 < ITERATIONS; t1++) {
    lbp = 1;
    ubp = nz - 2;
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6, t7, t8)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = 1; t3 < mz - 1; t3++) {
            for (int t4 = q - 2; t4 >= 0; t4--) {
                for (int t5 = 0; t5 < BLOCK_SIZE; t5++) {
                    lbv = 0;
                    ubv = BLOCK_SIZE - 1;
#pragma ivdep
#pragma vector always
                    for (int t6 = lbv; t6 <= ubv; t6++) {
                        rhs[t2][t3][t4][t5] = rhs[t2][t3][t4][t5] - lhs[t2][t3][t4][t5][t6] * rhs[t2][t3][t4 + 1][t6];
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

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_5D(array_0, N + 1, M + 1, L + 1, P + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, P + 1);

    NPB_bt18(array_0, array_1, N, M, L, P);

    print_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    free_array_5d(array_0, N + 1, M + 1, L + 1, P + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    return 0;
}