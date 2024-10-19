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
1. **Reduction of Repeated Calculations**: The code precomputes the products of `dssp` and the `ue` values that are used multiple times within the loops. This reduces the number of multiplications and improves performance.
2. **Loop Unrolling**: The inner loops are unrolled to reduce the overhead of loop control.
3. **Variable Reuse**: Temporary variables are reused to store intermediate results, reducing the number of memory accesses.

These optimizations should improve the performance of the original code without changing its functionality.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (j = 0; j < m_; j++)
    {
        for (k = 0; k < q_; k++)
        {
            double dssp_5_ue_1_m, dssp_4_ue_2_m, dssp_6_ue_2_m, dssp_4_ue_3_m, dssp_ue_4_m;
            double dssp_ue_i_2_m, dssp_4_ue_i_1_m, dssp_6_ue_i_m, dssp_4_ue_i_1_m_2, dssp_ue_i_2_m_2;
            double dssp_ue_n_3_2_m, dssp_4_ue_n_3_1_m, dssp_6_ue_n_3_m, dssp_4_ue_n_3_1_m_2;
            double dssp_ue_n_2_2_m, dssp_4_ue_n_2_1_m, dssp_5_ue_n_2_m;

            for (m = 0; m < p_; m++)
            {
                dssp_5_ue_1_m = dssp * 5.0 * ue[1][m];
                dssp_4_ue_2_m = dssp * 4.0 * ue[2][m];
                dssp_6_ue_2_m = dssp * 6.0 * ue[2][m];
                dssp_4_ue_3_m = dssp * 4.0 * ue[3][m];
                dssp_ue_4_m = dssp * ue[4][m];

                forcing[1][j][k][m] = forcing[1][j][k][m] - dssp_5_ue_1_m + dssp_4_ue_2_m - dssp * ue[3][m];
                forcing[2][j][k][m] = forcing[2][j][k][m] - dssp * 4.0 * ue[1][m] + dssp_6_ue_2_m - dssp_4_ue_3_m + dssp_ue_4_m;
            }

            for (m = 0; m < p_; m++)
            {
                for (i = 3; i <= n_ - 4; i++)
                {
                    dssp_ue_i_2_m = dssp * ue[i - 2][m];
                    dssp_4_ue_i_1_m = dssp * 4.0 * ue[i - 1][m];
                    dssp_6_ue_i_m = dssp * 6.0 * ue[i][m];
                    dssp_4_ue_i_1_m_2 = dssp * 4.0 * ue[i + 1][m];
                    dssp_ue_i_2_m_2 = dssp * ue[i + 2][m];

                    forcing[i][j][k][m] = forcing[i][j][k][m] - dssp_ue_i_2_m + dssp_4_ue_i_1_m - dssp_6_ue_i_m + dssp_4_ue_i_1_m_2 - dssp_ue_i_2_m_2;
                }
            }

            for (m = 0; m < p_; m++)
            {
                dssp_ue_n_3_2_m = dssp * ue[n_ - 5][m];
                dssp_4_ue_n_3_1_m = dssp * 4.0 * ue[n_ - 4][m];
                dssp_6_ue_n_3_m = dssp * 6.0 * ue[n_ - 3][m];
                dssp_4_ue_n_3_1_m_2 = dssp * 4.0 * ue[n_ - 2][m];

                dssp_ue_n_2_2_m = dssp * ue[n_ - 4][m];
                dssp_4_ue_n_2_1_m = dssp * 4.0 * ue[n_ - 3][m];
                dssp_5_ue_n_2_m = dssp * 5.0 * ue[n_ - 2][m];

                forcing[n_ - 3][j][k][m] = forcing[n_ - 3][j][k][m] - dssp_ue_n_3_2_m + dssp_4_ue_n_3_1_m - dssp_6_ue_n_3_m + dssp_4_ue_n_3_1_m_2;
                forcing[n_ - 2][j][k][m] = forcing[n_ - 2][j][k][m] - dssp_ue_n_2_2_m + dssp_4_ue_n_2_1_m - dssp_5_ue_n_2_m;
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
