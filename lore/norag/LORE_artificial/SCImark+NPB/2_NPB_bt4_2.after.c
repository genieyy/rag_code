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
/**/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (j = 0; j < m_; j++) {
        for (k = 0; k < q_; k++) {
            double dssp_times_5 = dssp * 5.0;
            double dssp_times_4 = dssp * 4.0;
            double dssp_times_6 = dssp * 6.0;

            for (m = 0; m < p_; m++) {
                double ue_1_m = ue[1][m];
                double ue_2_m = ue[2][m];
                double ue_3_m = ue[3][m];
                double ue_4_m = ue[4][m];

                forcing[1][j][k][m] -= dssp_times_5 * ue_1_m - dssp_times_4 * ue_2_m + dssp * ue_3_m;
                forcing[2][j][k][m] -= dssp_times_4 * (ue_1_m - ue_3_m) + dssp_times_6 * ue_2_m - dssp * ue_4_m;
            }

            for (m = 0; m < p_; m++) {
                double ue_m_minus_2, ue_m_minus_1, ue_m, ue_m_plus_1, ue_m_plus_2;

                for (i = 3; i <= n_ - 4; i++) {
                    ue_m_minus_2 = ue[i - 2][m];
                    ue_m_minus_1 = ue[i - 1][m];
                    ue_m = ue[i][m];
                    ue_m_plus_1 = ue[i + 1][m];
                    ue_m_plus_2 = ue[i + 2][m];

                    forcing[i][j][k][m] -= dssp * (ue_m_minus_2 - dssp_times_4 * ue_m_minus_1 + dssp_times_6 * ue_m - dssp_times_4 * ue_m_plus_1 + ue_m_plus_2);
                }
            }

            for (m = 0; m < p_; m++) {
                double ue_n_minus_5 = ue[n_ - 5][m];
                double ue_n_minus_4 = ue[n_ - 4][m];
                double ue_n_minus_3 = ue[n_ - 3][m];
                double ue_n_minus_2 = ue[n_ - 2][m];

                forcing[n_ - 3][j][k][m] -= dssp * (ue_n_minus_5 - dssp_times_4 * ue_n_minus_4 + dssp_times_6 * ue_n_minus_3 - dssp_times_4 * ue_n_minus_2);
                forcing[n_ - 2][j][k][m] -= dssp * (ue_n_minus_4 - dssp_times_4 * ue_n_minus_3 + dssp_times_5 * ue_n_minus_2);
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
