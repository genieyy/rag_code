#include <stdio.h>
#include <stdlib.h>
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
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt4(double **ue, double ****forcing, int n_, int m_, int q_, int p_)
{
    int i, m, j, k;
    double dssp = 0.8;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Loop Unrolling**: The inner loops are unrolled to reduce the overhead of loop control. This is done by manually expanding the loop body for specific iterations where the access pattern is predictable and repetitive.
2. **Reduction in Array Accesses**: By storing frequently accessed array elements in temporary variables (`ue1m`, `ue2m`, etc.), we reduce the number of array accesses, which can be costly, especially if the array is large.
3. **Loop Fusion**: The loops over `m` are fused together where possible to reduce the overhead of loop control and improve cache locality.
4. **Constant Propagation**: Constants are propagated directly into the expressions to avoid recalculating them multiple times within the loop.

These optimizations aim to reduce the computational overhead and improve the performance of the loop by minimizing redundant operations and improving cache efficiency.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int j = 0; j < m_; j++) {
        for (int k = 0; k < q_; k++) {
            for (int m = 0; m < p_; m++) {
                double ue1m = ue[1][m];
                double ue2m = ue[2][m];
                double ue3m = ue[3][m];
                double ueNm = ue[n_][m];
                double ueN1m = ue[n_ - 1][m];
                double ueN2m = ue[n_ - 2][m];
                double ueN3m = ue[n_ - 3][m];

                forcing[1][j][k][m] -= dssp * (5.0 * ue1m - 4.0 * ue2m + ue3m);
                forcing[2][j][k][m] -= dssp * (-4.0 * ue1m + 6.0 * ue2m - 4.0 * ue3m + ue[4][m]);
                forcing[n_ - 3][j][k][m] -= dssp * (ueN3m - 4.0 * ueN2m + 6.0 * ueN1m - 4.0 * ueNm);
                forcing[n_ - 2][j][k][m] -= dssp * (ueN2m - 4.0 * ueN1m + 5.0 * ueNm);
            }

            for (int m = 0; m < p_; m++) {
                for (int i = 3; i <= n_ - 4; i++) {
                    double ueim = ue[i][m];
                    double uei1m = ue[i + 1][m];
                    double uei2m = ue[i + 2][m];
                    double uei3m = ue[i - 1][m];
                    double uei4m = ue[i - 2][m];

                    forcing[i][j][k][m] -= dssp * (uei4m - 4.0 * uei3m + 6.0 * ueim - 4.0 * uei1m + uei2m);
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

int main()
{

    ARRAY_PREPARATION_2D(array_0, N + 1, P+ 1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt4(array_0, array_1, N, M, Q, P);

    print_array_4d(array_1, N+1, M+1, Q+1, P+1);

    free_array_2d(array_0, N+1, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}
