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
1. **Constant Multiplication**: Precompute the multiplication of `dssp` with constants (`5.0`, `4.0`, `6.0`) to avoid redundant multiplications inside the loops.
2. **Temporary Variables**: Store intermediate results of `ue` array accesses in temporary variables to reduce the number of array accesses, which can be costly.
3. **Loop Order**: The loop order remains the same, but the inner loops are optimized by reducing redundant calculations and array accesses.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (j = 0; j < m_; j++)
    {
        for (k = 0; k < q_; k++)
        {
            double dssp_times_5 = dssp * 5.0;
            double dssp_times_4 = dssp * 4.0;
            double dssp_times_6 = dssp * 6.0;

            for (m = 0; m < p_; m++)
            {
                double ue1_m = ue[1][m];
                double ue2_m = ue[1 + 1][m];
                double ue3_m = ue[1 + 2][m];

                forcing[1][j][k][m] -= dssp_times_5 * ue1_m - dssp_times_4 * ue2_m + dssp * ue3_m;

                double ue_m2 = ue[2 - 1][m];
                double ue_m3 = ue[2][m];
                double ue_m4 = ue[2 + 1][m];
                double ue_m5 = ue[2 + 2][m];

                forcing[2][j][k][m] -= -dssp_times_4 * ue_m2 + dssp_times_6 * ue_m3 - dssp_times_4 * ue_m4 + dssp * ue_m5;
            }

            for (m = 0; m < p_; m++)
            {
                for (i = 1 * 3; i <= n_ - 3 * 1 - 1; i++)
                {
                    double ue_im2 = ue[i - 2][m];
                    double ue_im1 = ue[i - 1][m];
                    double ue_i = ue[i][m];
                    double ue_ip1 = ue[i + 1][m];
                    double ue_ip2 = ue[i + 2][m];

                    forcing[i][j][k][m] -= ue_im2 - dssp_times_4 * ue_im1 + dssp_times_6 * ue_i - dssp_times_4 * ue_ip1 + dssp * ue_ip2;
                }
            }

            for (m = 0; m < p_; m++)
            {
                double ue_nm3m2 = ue[n_ - 3 - 2][m];
                double ue_nm3m1 = ue[n_ - 3 - 1][m];
                double ue_nm3 = ue[n_ - 3][m];
                double ue_nm3p1 = ue[n_ - 3 + 1][m];

                forcing[n_ - 3][j][k][m] -= ue_nm3m2 - dssp_times_4 * ue_nm3m1 + dssp_times_6 * ue_nm3 - dssp_times_4 * ue_nm3p1;

                double ue_nm2m2 = ue[n_ - 2 - 2][m];
                double ue_nm2m1 = ue[n_ - 2 - 1][m];
                double ue_nm2 = ue[n_ - 2][m];

                forcing[n_ - 2][j][k][m] -= ue_nm2m2 - dssp_times_4 * ue_nm2m1 + dssp_times_5 * ue_nm2;
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
