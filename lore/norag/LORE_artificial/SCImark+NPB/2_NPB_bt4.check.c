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
for (int iter = 0; iter < ITERATIONS; iter++){
    for (j = 0; j < m_; j++)
    {
        for (k = 0; k < q_; k++)
        {
            for (m = 0; m < p_; m++)
            {
                forcing[1][j][k][m] = forcing[1][j][k][m] - dssp *
                                                                (5.0 * ue[1][m] - 4.0 * ue[1 + 1][m] + ue[1 + 2][m]);

                forcing[2][j][k][m] = forcing[2][j][k][m] - dssp *
                                                                (-4.0 * ue[2 - 1][m] + 6.0 * ue[2][m] -
                                                                 4.0 * ue[2 + 1][m] + ue[2 + 2][m]);
            }

            for (m = 0; m < p_; m++)
            {
                for (i = 1 * 3; i <= n_ - 3 * 1 - 1; i++)
                {
                    forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                                                    (ue[i - 2][m] - 4.0 * ue[i - 1][m] +
                                                                     6.0 * ue[i][m] - 4.0 * ue[i + 1][m] + ue[i + 2][m]);
                }
            }

            for (m = 0; m < p_; m++)
            {

                forcing[n_ - 3][j][k][m] = forcing[n_ - 3][j][k][m] - dssp *
                                                                          (ue[n_ - 3 - 2][m] - 4.0 * ue[n_ - 3 - 1][m] +
                                                                           6.0 * ue[N - 3][m] - 4.0 * ue[n_ - 3 + 1][m]);
                forcing[n_ - 2][j][k][m] = forcing[n_ - 2][j][k][m] - dssp *
                                                                          (ue[n_ - 2 - 2][m] - 4.0 * ue[n_ - 2 - 1][m] + 5.0 * ue[n_ - 2][m]);
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
